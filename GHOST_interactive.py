#WRITTEN BY DENE BOWDALO

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

#GHOST_interactive.py

#executable of the GHOST interactive tool

###------------------------------------------------------------------------------------###
###IMPORT CONFIGURATION FILE VARIABLES
###------------------------------------------------------------------------------------###

from configuration import n_CPUs, obs_root, exp_root, sequential_colourmap, sequential_colourmap_warm, diverging_colourmap, cartopy_data_dir

###------------------------------------------------------------------------------------###
###IMPORT BUILT-IN MODULES
###------------------------------------------------------------------------------------###
import bisect
from collections import OrderedDict
import copy
import datetime
from functools import partial
import gc
import glob
import multiprocessing
import os
import sys
from weakref import WeakKeyDictionary

###------------------------------------------------------------------------------------###
###IMPORT THIRD-PARTY MODULES
###------------------------------------------------------------------------------------###
import cartopy
#set cartopy data directory
cartopy.config['pre_existing_data_dir'] = cartopy_data_dir
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
# Make sure that we are using Qt5 backend with matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector
from matplotlib.gridspec import GridSpec
from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from PyQt5 import QtCore, QtWidgets, QtGui
import scipy.stats
import seaborn as sns

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

class NavigationToolbar(NavigationToolbar):

    '''define class that updates available buttons on matplotlib toolbar'''    

    #only display wanted buttons
    NavigationToolbar.toolitems = (
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        (None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        (None, None, None, None),
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
        (None, None, None, None)
    )

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

class QVLine(QtWidgets.QFrame):

    '''define class that generates vertical separator line'''

    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QtWidgets.QFrame.VLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)

class QHLine(QtWidgets.QFrame):

    '''define class that generates horizontal separator line'''

    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QtWidgets.QFrame.HLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)

#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#

class pop_up_window(QtWidgets.QWidget):

    '''define class that generates generalised pop-up window'''

    def __init__(self, window_type = '', window_titles=[], checkbox_labels=[], default_checkbox_selection=[], selected_indices={}):
        super(pop_up_window, self).__init__()
        
        #add input arguments to self
        self.window_type = window_type
        self.window_titles = window_titles
        self.checkbox_labels = checkbox_labels
        self.default_checkbox_selection = default_checkbox_selection
        self.selected_indices = selected_indices        

        #create UI
        self.initUI()

    def initUI(self):

        '''initialise user interface'''
        
        #set window title
        self.setWindowTitle(self.window_type)

        #get pop up window dimensions
        window_width = self.width()
        window_height = self.height()

        #create parent layout to hold N horizontally laid out windows
        parent_layout = QtWidgets.QHBoxLayout()
        parent_layout.setAlignment(QtCore.Qt.AlignTop) 
        
        #get N of nested windows from length of window_titles list
        n_nested_windows = len(self.window_titles)

        #create frame to hold checkboxes
        self.checkboxes = [[] for x in range(len(self.checkbox_labels))]

        #iterate through N nested windows, creating each one and placing in parent frame accordingly
        for nested_window_N in range(n_nested_windows):

            #------------------------------------------------------------------------#
            #create nested parent layout to pull together a title, button row, and grid of checkboxes
            nested_parent_layout = QtWidgets.QVBoxLayout() 
            nested_parent_layout.setAlignment(QtCore.Qt.AlignTop) 

            #define spacing/margin variables
            nested_parent_layout.setSpacing(10)
            nested_parent_layout.setContentsMargins(0,0,0,0)

            #------------------------------------------------------------------------#
            #create title label
            title_label = QtWidgets.QLabel(self, text = self.window_titles[nested_window_N])
            title_label.setAlignment(QtCore.Qt.AlignCenter)
            myFont=QtGui.QFont()
            myFont.setPointSize(22)
            title_label.setFont(myFont)

            #------------------------------------------------------------------------#
            #create row of buttons
            button_row = QtWidgets.QHBoxLayout()
            button_row.setAlignment(QtCore.Qt.AlignLeft)
            select_all_button = QtWidgets.QPushButton("Select All")
            select_all_button.setMinimumWidth(100)
            select_all_button.setMaximumWidth(100)
            clear_all_button = QtWidgets.QPushButton("Clear All")
            clear_all_button.setMinimumWidth(100)
            clear_all_button.setMaximumWidth(100)
            select_default_button = QtWidgets.QPushButton("Select Default")
            select_default_button.setMinimumWidth(120)
            select_default_button.setMaximumWidth(120)
            
            #order buttons in grid layout
            button_row.addWidget(select_all_button)
            button_row.addWidget(clear_all_button)
            button_row.addWidget(select_default_button)
            
            #add connectivity to buttons
            select_all_button.clicked.connect(partial(self.select_all, nested_window_N))
            clear_all_button.clicked.connect(partial(self.clear_all, nested_window_N))        
            select_default_button.clicked.connect(partial(self.select_all_default, nested_window_N))
        
            #------------------------------------------------------------------------#
            #create grid of checkboxes
            checkbox_grid = QtWidgets.QGridLayout()

            #define spacing/margin variables
            checkbox_grid.setHorizontalSpacing(1)
            checkbox_grid.setVerticalSpacing(1)
            checkbox_grid.setContentsMargins(0,0,0,0)
    
            #create checkboxes
            #force a new column to be started if the available vertical space for each row in grid goes below a critical value (18.2)
            #if checkbox has been previously selected (without updating network/resolution/species), then reselect it again
            row_n = 0
            column_n = 0
            current_selected_indices = self.selected_indices[self.window_type][nested_window_N]
            for ii, val in enumerate(self.checkbox_labels[nested_window_N]):
                row_available_space = window_height/(row_n+1)
                if row_available_space < 18.2:
                    column_n+=2
                    row_n = 0
                self.checkboxes[nested_window_N].append(QtWidgets.QCheckBox(val))
                if ii in current_selected_indices: 
                    self.checkboxes[nested_window_N][ii].setCheckState(QtCore.Qt.Checked)
                checkbox_grid.addWidget(self.checkboxes[nested_window_N][ii], row_n, column_n)
                row_n +=1

            #------------------------------------------------------------------------#
            #position title, button row and checkbox grid in nested parent layout 

            #add title to nested parent frame
            nested_parent_layout.addWidget(title_label)

            #add button row to nested parent frame
            nested_parent_layout.addLayout(button_row)
        
            #add checkbox grid to nested parent frame
            nested_parent_layout.addLayout(checkbox_grid)

            #------------------------------------------------------------------------#
            #add nested parent layout to parent frame
            parent_layout.addLayout(nested_parent_layout)

            #------------------------------------------------------------------------#
            #add vertical separation line (if not last nested window that is being iterated through)            
            if (nested_window_N+1) != n_nested_windows:
                parent_layout.addWidget(QVLine())
            
        #------------------------------------------------------------------------#
        #set finalised layout
        self.setLayout(parent_layout)

        #maximise window to fit screen
        self.showMaximized()  
        
        #------------------------------------------------------------------------#

        #setup event to get selected checkbox indices when closing window
        quit = QtWidgets.QAction("Quit", self)
        quit.triggered.connect(self.closeEvent)

    #------------------------------------------------------------------------#
    #------------------------------------------------------------------------#

    def select_all(self, nested_window_N):
        '''function to select all checkboxes'''
        for ii, val in enumerate(self.checkboxes[nested_window_N]):
            self.checkboxes[nested_window_N][ii].setCheckState(QtCore.Qt.Checked)
    
    def clear_all(self, nested_window_N):
        '''function to clear all checkboxes'''
        for ii, val in enumerate(self.checkboxes[nested_window_N]):
            self.checkboxes[nested_window_N][ii].setCheckState(QtCore.Qt.Unchecked)
                
    def select_all_default(self, nested_window_N):
        '''function to select all default selected checkboxes'''
        #unselect all checkboxes first 
        for ii, val in enumerate(self.checkboxes[nested_window_N]):
            self.checkboxes[nested_window_N][ii].setCheckState(QtCore.Qt.Unchecked)
        
        #now select only desired default checkboxes
        for ii in self.default_checkbox_selection[nested_window_N]:
            self.checkboxes[nested_window_N][ii].setCheckState(QtCore.Qt.Checked)
        
    #------------------------------------------------------------------------#
    #------------------------------------------------------------------------#

    def closeEvent(self, event):

        '''function to get indices of selected checkboxes upon closing of pop-up window'''

        selected_indices = []

        for ii in range(len(self.checkboxes)): 
            checked_indices = np.array([], dtype=np.uint8)
            for jj in range(len(self.checkboxes[ii])):
                if self.checkboxes[ii][jj].checkState() == QtCore.Qt.Checked:
                    checked_indices = np.append(checked_indices, jj)
            selected_indices.append(checked_indices)   

        self.selected_indices[self.window_type] = selected_indices

#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#
class generate_GHOST_interactive_dashboard(QtWidgets.QWidget):

    '''define class that generates GHOST interactive dashboard'''     
 
    def __init__(self, read_type):
        super(generate_GHOST_interactive_dashboard, self).__init__()
        
        #put read_type into self
        self.read_type = read_type

        #create UI
        self.initUI()

    #------------------------------------------------------------------------#

    def initUI(self):

        '''initialise user interface'''

        #set window title
        self.window_title = "GHOST Interactive"
        self.setWindowTitle(self.window_title)

        #create parent layout to pull together a configuration bar, a MPL navigation toolbar, and a MPL canvas of plots
        parent_layout = QtWidgets.QVBoxLayout()
        
        #define spacing/margin variables
        parent_layout.setSpacing(0)
        parent_layout.setContentsMargins(0,0,0,0)
        
        #------------------------------------------------------------------------#
        #setup configuration bar with combo boxes, input boxes and buttons
        #use a gridded layout to place objects
        config_bar = QtWidgets.QGridLayout()

        #define spacing/margin variables
        config_bar.setHorizontalSpacing(3)
        config_bar.setVerticalSpacing(1)
        config_bar.setContentsMargins(5,0,0,0)
        config_bar.setAlignment(QtCore.Qt.AlignLeft) 

        #define all configuration box objects (labels, comboboxes etc.)
        title_font=QtGui.QFont()
        title_font.setUnderline(True)
        self.lb_data_selection = QtWidgets.QLabel(self, text = "Data Selection")
        self.lb_data_selection.setFont(title_font)
        self.lb_data_selection.setToolTip('Setup configuration of data to read into memory')
        self.bu_read = QtWidgets.QPushButton('READ', self)
        self.bu_read.setMinimumWidth(40)
        self.bu_read.setMaximumWidth(40)
        self.bu_read.setStyleSheet("color: red;")
        self.bu_read.setToolTip('Read selected configuration of data into memory')
        self.ch_colocate = QtWidgets.QCheckBox("Colocate")
        self.cb_network = QtWidgets.QComboBox(self)
        self.cb_network.setMinimumWidth(95)
        self.cb_network.setMaximumWidth(95)
        self.cb_network.setToolTip('Select providing observational data network')
        self.cb_resolution = QtWidgets.QComboBox(self)
        self.cb_resolution.setMinimumWidth(95)
        self.cb_resolution.setMaximumWidth(95)
        self.cb_resolution.setToolTip('Select temporal resolution of data')
        self.cb_matrix = QtWidgets.QComboBox(self)
        self.cb_matrix.setMinimumWidth(95)
        self.cb_matrix.setMaximumWidth(95)
        self.cb_matrix.setToolTip('Select data matrix')
        self.cb_species = QtWidgets.QComboBox(self)
        self.cb_species.setMinimumWidth(95)
        self.cb_species.setMaximumWidth(95)
        self.cb_species.setToolTip('Select species')
        self.bu_experiments = QtWidgets.QPushButton('EXPS', self)
        self.bu_experiments.setMinimumWidth(44)
        self.bu_experiments.setMaximumWidth(44)
        self.bu_experiments.setToolTip('Select experiment/s data to read')
        self.le_start_date = QtWidgets.QLineEdit(self)
        self.le_start_date.setMinimumWidth(70)
        self.le_start_date.setMaximumWidth(70)
        self.le_start_date.setToolTip('Set data start date: YYYYMMDD')
        self.le_end_date = QtWidgets.QLineEdit(self)
        self.le_end_date.setMinimumWidth(70)
        self.le_end_date.setMaximumWidth(70)
        self.le_end_date.setToolTip('Set data end date: YYYYMMDD')
        self.bu_QA = QtWidgets.QPushButton('QA', self)
        self.bu_QA.setMinimumWidth(45)
        self.bu_QA.setMaximumWidth(45)
        self.bu_QA.setToolTip('Select standardised quality assurance flags to filter by')
        self.bu_flags = QtWidgets.QPushButton('FLAGS', self)
        self.bu_flags.setMinimumWidth(44)
        self.bu_flags.setMaximumWidth(44)
        self.bu_flags.setToolTip('Select standardised data reporter provided flags to filter by')
        self.bu_classifications = QtWidgets.QPushButton('CLASS', self)
        self.bu_classifications.setMinimumWidth(44)
        self.bu_classifications.setMaximumWidth(44)
        self.bu_classifications.setToolTip('Select standardised classifications to filter by')
        self.vertical_splitter_1 = QVLine()
        self.vertical_splitter_1.setMaximumWidth(20)
        self.lb_data_filter = QtWidgets.QLabel(self, text = "Data Filter")
        self.lb_data_filter.setMinimumWidth(65)
        self.lb_data_filter.setMaximumWidth(65)
        self.lb_data_filter.setFont(title_font)
        self.bu_screen = QtWidgets.QPushButton('FILTER', self)
        self.bu_screen.setMinimumWidth(46)
        self.bu_screen.setMaximumWidth(46)
        self.bu_screen.setStyleSheet("color: blue;")
        self.lb_data_bounds = QtWidgets.QLabel(self, text = "Bounds")
        self.lb_data_bounds.setMinimumWidth(47)
        self.lb_data_bounds.setMaximumWidth(47)
        self.le_minimum_value = QtWidgets.QLineEdit(self)
        self.le_minimum_value.setMinimumWidth(60)
        self.le_minimum_value.setMaximumWidth(60)
        self.le_maximum_value = QtWidgets.QLineEdit(self)
        self.le_maximum_value.setMinimumWidth(60)
        self.le_maximum_value.setMaximumWidth(60)
        self.lb_minimum_data_availability = QtWidgets.QLabel(self, text = "% Min.")
        self.lb_minimum_data_availability.setMinimumWidth(47)
        self.lb_minimum_data_availability.setMaximumWidth(47)
        self.le_minimum_data_availability = QtWidgets.QLineEdit(self)
        self.le_minimum_data_availability.setMinimumWidth(45)
        self.le_minimum_data_availability.setMaximumWidth(45)
        self.bu_methods = QtWidgets.QPushButton('METHOD', self)
        self.bu_methods.setMinimumWidth(60)
        self.bu_methods.setMaximumWidth(60)
        self.vertical_splitter_2 = QVLine()
        self.vertical_splitter_2.setMaximumWidth(20)
        self.lb_z = QtWidgets.QLabel(self, text = "Map Z")
        self.lb_z.setFont(title_font)
        self.cb_z_stat = QtWidgets.QComboBox(self)
        self.cb_z_stat.setMinimumWidth(80)
        self.cb_z_stat.setMaximumWidth(80)
        self.cb_z1 = QtWidgets.QComboBox(self)
        self.cb_z1.setMinimumWidth(125)
        self.cb_z1.setMaximumWidth(125)
        self.cb_z2 = QtWidgets.QComboBox(self)
        self.cb_z2.setMinimumWidth(125)
        self.cb_z2.setMaximumWidth(125)
        self.vertical_splitter_3 = QVLine()
        self.vertical_splitter_3.setMaximumWidth(20)
        self.lb_experiment_bias = QtWidgets.QLabel(self, text = "Exp. Bias")
        self.lb_experiment_bias.setFont(title_font)
        self.cb_experiment_bias_type = QtWidgets.QComboBox(self)
        self.cb_experiment_bias_type.setMinimumWidth(100)
        self.cb_experiment_bias_type.setMaximumWidth(100)
        self.cb_experiment_bias_stat = QtWidgets.QComboBox(self)
        self.cb_experiment_bias_stat.setMinimumWidth(100)
        self.cb_experiment_bias_stat.setMaximumWidth(100)
        self.vertical_splitter_4 = QVLine()
        self.vertical_splitter_4.setMaximumWidth(20)
        self.lb_station_selection = QtWidgets.QLabel(self, text = "Site Select")
        self.lb_station_selection.setFont(title_font)
        self.ch_select_all = QtWidgets.QCheckBox("All")
        self.ch_intersect = QtWidgets.QCheckBox("Intersect")

        #position objects on gridded configuration bar
        config_bar.addWidget(self.lb_data_selection, 0, 0, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_read, 0, 1, 1, 1, QtCore.Qt.AlignCenter)  
        config_bar.addWidget(self.ch_colocate, 0, 3, 1, 2, QtCore.Qt.AlignCenter)
        config_bar.addWidget(self.cb_network, 1, 0)
        config_bar.addWidget(self.cb_resolution, 2, 0)
        config_bar.addWidget(self.cb_matrix, 1, 1)
        config_bar.addWidget(self.cb_species, 2, 1)
        config_bar.addWidget(self.le_start_date, 1, 2)
        config_bar.addWidget(self.le_end_date, 2, 2) 
        config_bar.addWidget(self.bu_QA, 1, 3) 
        config_bar.addWidget(self.bu_flags, 2, 3) 
        config_bar.addWidget(self.bu_classifications, 1, 4)
        config_bar.addWidget(self.bu_experiments, 2, 4)
        config_bar.addWidget(self.vertical_splitter_1, 0, 5, 3, 1)
        config_bar.addWidget(self.lb_data_filter, 0, 6, 1, 2, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_screen, 0, 7, 1, 2, QtCore.Qt.AlignCenter) 
        config_bar.addWidget(self.lb_data_bounds, 1, 6)
        config_bar.addWidget(self.le_minimum_value, 1, 7)
        config_bar.addWidget(self.le_maximum_value, 1, 8)
        config_bar.addWidget(self.lb_minimum_data_availability, 2, 6)
        config_bar.addWidget(self.le_minimum_data_availability, 2, 7)
        config_bar.addWidget(self.bu_methods, 2, 7, 1, 2, QtCore.Qt.AlignRight)
        config_bar.addWidget(self.vertical_splitter_2, 0, 9, 3, 1)
        config_bar.addWidget(self.lb_z, 0, 10, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_z_stat, 0, 10, QtCore.Qt.AlignRight)
        config_bar.addWidget(self.cb_z1, 1, 10)
        config_bar.addWidget(self.cb_z2, 2, 10)
        config_bar.addWidget(self.vertical_splitter_3, 0, 11, 3, 1)        
        config_bar.addWidget(self.lb_experiment_bias, 0, 12, QtCore.Qt.AlignCenter)
        config_bar.addWidget(self.cb_experiment_bias_type, 1, 12)
        config_bar.addWidget(self.cb_experiment_bias_stat, 2, 12)
        config_bar.addWidget(self.vertical_splitter_4, 0, 13, 3, 1)
        config_bar.addWidget(self.lb_station_selection, 0, 14)
        config_bar.addWidget(self.ch_select_all, 1, 14)
        config_bar.addWidget(self.ch_intersect, 2, 14)

        #enable dynamic updating of configuration bar fields which filter data files
        self.cb_network.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.cb_resolution.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.cb_matrix.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.cb_species.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.le_start_date.textChanged.connect(self.config_bar_params_change_handler)
        self.le_end_date.textChanged.connect(self.config_bar_params_change_handler)

        #enable pop up configuration windows
        self.bu_experiments.clicked.connect(self.handle_pop_up_experiments_window)
        self.bu_flags.clicked.connect(self.handle_pop_up_flags_window)
        self.bu_QA.clicked.connect(self.handle_pop_up_qa_window)
        self.bu_classifications.clicked.connect(self.handle_pop_up_classifications_window)
        self.bu_methods.clicked.connect(self.handle_pop_up_methods_window)

        #define data provider flags
        self.flag_names = np.array(sorted(standard_data_flag_codes, key=standard_data_flag_codes.get))
        self.flag_codes = np.sort(list(standard_data_flag_codes.values())) 
        self.flag_default_codes = np.array([1, 2, 3, 10, 11, 12, 13, 14, 15, 16, 20, 21, 24, 25, 26, 29, 30, 31, 32, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 90, 150, 154, 155, 156, 157], dtype=np.uint8)
        self.flag_default_inds = np.array([np.where(self.flag_codes == code)[0][0] for code in self.flag_default_codes], dtype=np.uint8)

        #define qa flags
        self.qa_names = np.array(sorted(standard_qa_flag_codes, key=standard_qa_flag_codes.get))
        self.qa_codes = np.sort(list(standard_qa_flag_codes.values()))        
        self.qa_default_codes = np.array([0, 1, 2, 3, 4, 70, 80, 81, 82, 85, 87, 90, 94, 100, 106, 107, 108], dtype=np.uint8)
        self.qa_default_inds = np.array([np.where(self.qa_codes == code)[0][0] for code in self.qa_default_codes], dtype=np.uint8)

        #define classification flags        
        self.classification_names = np.array(sorted(standard_classification_flag_codes, key=standard_classification_flag_codes.get))
        self.classification_codes = np.sort(list(standard_classification_flag_codes.values()))
        self.classification_default_codes_to_retain = np.array([7], dtype=np.uint8)
        self.classification_default_codes_to_remove = np.array([0, 1, 5, 6], dtype=np.uint8)
        self.classification_default_inds_to_retain = np.array([np.where(self.classification_codes == code)[0][0] for code in self.classification_default_codes_to_retain], dtype=np.uint8)
        self.classification_default_inds_to_remove = np.array([np.where(self.classification_codes == code)[0][0] for code in self.classification_default_codes_to_remove], dtype=np.uint8)

        #create dictionary to hold indices of selected values in pop-up windows
        self.selected_indices = {'EXPERIMENTS':[[]], 'FLAGS':[[]], 'QA':[[]], 'CLASSIFICATIONS':[[],[]], 'METHODS':[[]]}
        
        #initialise configuration bar fields
        self.config_bar_initialisation = True
        self.update_configuration_bar_fields()  
        self.config_bar_initialisation = False

        #------------------------------------------------------------------------#
        #setup MPL canvas of plots

        #set variable that blocks updating of MPL canvas until some data has been read
        self.block_MPL_canvas_updates = True
        self.mpl_canvas = MPL_Canvas(self)

        #------------------------------------------------------------------------#
        #enable interactivity of functions which update MPL canvas

        #enable READ button
        self.bu_read.clicked.connect(self.handle_data_selection_update)

        #enable interactivity of colocation checkbox
        self.ch_colocate.stateChanged.connect(self.mpl_canvas.handle_colocate_update)

        #enable FILTER button
        self.bu_screen.clicked.connect(self.mpl_canvas.handle_data_filter_update)

        #enable updating of map z statistic
        self.cb_z_stat.currentTextChanged.connect(self.mpl_canvas.handle_map_z_statistic_update)
        self.cb_z1.currentTextChanged.connect(self.mpl_canvas.handle_map_z_statistic_update)
        self.cb_z2.currentTextChanged.connect(self.mpl_canvas.handle_map_z_statistic_update)

        #enable updating of experiment bias statistic
        self.cb_experiment_bias_type.currentTextChanged.connect(self.mpl_canvas.handle_experiment_bias_update) 
        self.cb_experiment_bias_stat.currentTextChanged.connect(self.mpl_canvas.handle_experiment_bias_update) 

        #enable interactivity of station selection checkboxes
        self.ch_select_all.stateChanged.connect(self.mpl_canvas.select_all_stations)
        self.ch_intersect.stateChanged.connect(self.mpl_canvas.select_intersect_stations)

        #------------------------------------------------------------------------#
        #generate MPL navigation toolbar
        self.navi_toolbar = NavigationToolbar(self.mpl_canvas, self)

        #------------------------------------------------------------------------#
        #position config bar, navigation toolbar and MPL canvas and elements in parent layout 
        
        #add config bar to parent frame
        parent_layout.addLayout(config_bar)

        #add MPL navigation toolbar to parent frame
        parent_layout.addWidget(self.navi_toolbar)

        #add MPL canvas of plots to parent frame
        parent_layout.addWidget(self.mpl_canvas)

        #------------------------------------------------------------------------#
        #set finalised layout
        self.setLayout(parent_layout)

        #------------------------------------------------------------------------#
        #plot whole dashboard
        self.show()

        #maximise window to fit screen
        self.showMaximized()  

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#
    
    def update_configuration_bar_fields(self):

        '''define function that initialises/updates configuration bar fields'''
        
        #set variable to block interactive handling while updating config bar parameters
        self.block_config_bar_handling_updates = True

        #set some default configuration values when initialising config bar
        if self.config_bar_initialisation == True:
            #set initially selected/active start-end date as default 201601-201701
            self.le_start_date.setText('20160101')
            self.le_end_date.setText('20170101')
            self.selected_start_date = int(self.le_start_date.text())
            self.selected_end_date = int(self.le_end_date.text())
            self.selected_start_date_firstdayofmonth = int(str(self.selected_start_date)[:6]+'01')
            self.active_start_date = int(self.le_start_date.text())
            self.active_end_date = int(self.le_end_date.text())
            self.date_range_has_changed = False

            #set initially selected minimum data availability % to 0.0
            self.le_minimum_data_availability.setText('0.0')

            #set selected/active values of other fields to be initially None
            self.selected_network = None
            self.active_network = None
            self.selected_resolution = None
            self.active_resolution = None
            self.selected_matrix = None
            self.active_matrix = None
            self.selected_species = None
            self.active_species = None

            #set selected/active values of variables associated with pop up windows to be empty lists
            self.active_experiment_grids = []  
            self.active_qa_inds = []
            self.active_flag_inds = []  
            self.active_classifications_to_retain_inds = []
            self.active_classifications_to_remove_inds = []      

            #set initial time array to be None
            self.time_array = None      

            #set initial station references to be empty list
            self.station_references = []
            #set initial unique station methods to be empty list
            self.station_unique_methods = [] 
        
            #--------------------------------------------------------------------------------------#
            #gather all observational data
            #create nested dictionary storing all observational species data by species matrix, by temporal resolution, by network, associated with list of start YYYYMM yearmonths of data files

            self.all_observation_data = {}

            #set all available networks
            available_networks = ['EBAS','EIONET']
        
            #set all available temporal resolutions
            available_resolutions = ['hourly','daily','monthly']

            #iterate through available networks
            for network in available_networks:

                self.all_observation_data[network] = {}

                #iterate through available resolutions
                for resolution in available_resolutions:

                    #write nested empty dictionary for resolution
                    self.all_observation_data[network][resolution] = {}

                    #get available species for network/resolution
                    available_species = os.listdir('%s/%s/%s'%(obs_root, network, resolution))

                    #iterate through available files per species
                    for species in available_species:   

                        #get all netCDF monthly files per species
                        species_files = os.listdir('%s/%s/%s/%s'%(obs_root, network, resolution, species))

                        #get monthly start date (YYYYMM) of all species files
                        species_files_yearmonths = [int(f.split('_')[-1][:6]+'01') for f in species_files if f != 'temporary'] 

                        #get matrix for current species
                        matrix = parameter_dictionary[species]['matrix']
            
                        if matrix not in list(self.all_observation_data[network][resolution].keys()):
                            #write nested empty dictionary for matrix
                            self.all_observation_data[network][resolution][matrix] = {}                            

                        #write nested dictionary for species, with associated file yearmonths
                        self.all_observation_data[network][resolution][matrix][species] = species_files_yearmonths

            #create dictionary of observational data inside date range
            self.get_valid_obs_files_in_date_range()

        #--------------------------------------------------------------------------------------#
        #if date range has changed then update available observational data dictionary 
        if self.date_range_has_changed == True:
             self.get_valid_obs_files_in_date_range()

        #--------------------------------------------------------------------------------------#
        #initialise/update fields - maintain previously selected values wherever possible

        #clear fields
        self.cb_network.clear()
        self.cb_resolution.clear()
        self.cb_matrix.clear()
        self.cb_species.clear()
        
        #if have no available observational data, return from function, updating variable informing that have no data 
        if len(self.available_observation_data) == 0:
            self.no_data_to_read = True
            #unset variable to allow interactive handling from now
            self.block_config_bar_handling_updates = False
            return
        else:
            self.no_data_to_read = False

        #update network field
        available_networks = sorted(list(self.available_observation_data.keys()))
        self.cb_network.addItems(available_networks)
        if self.selected_network in available_networks:
            self.cb_network.setCurrentText(self.selected_network)
        if self.config_bar_initialisation == True:
            self.selected_network = self.cb_network.currentText()
        
        #update resolution field
        available_resolutions = list(self.available_observation_data[self.cb_network.currentText()].keys())
        #manually force order of available resolutions
        resolution_order_dict = {'hourly':1, 'daily':2, 'monthly':3}
        available_resolutions = sorted(available_resolutions, key=resolution_order_dict.__getitem__) 
        self.cb_resolution.addItems(available_resolutions)
        if self.selected_resolution in available_resolutions:
            self.cb_resolution.setCurrentText(self.selected_resolution)
        if self.config_bar_initialisation == True:
            self.selected_resolution = self.cb_resolution.currentText()

        #update matrix field
        available_matrices = sorted(self.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()])
        self.cb_matrix.addItems(available_matrices)
        if self.selected_matrix in available_matrices:
            self.cb_matrix.setCurrentText(self.selected_matrix) 
        if self.config_bar_initialisation == True:
            self.selected_matrix = self.cb_matrix.currentText()

        #update species field
        available_species = sorted(self.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()][self.cb_matrix.currentText()])
        self.cb_species.addItems(available_species)
        if self.selected_species in available_species:
            self.cb_species.setCurrentText(self.selected_species) 
        if self.config_bar_initialisation == True:
            self.selected_species = self.cb_species.currentText()  

        #update selected indices for QA
        #if initialising config bar then check default selection
        if self.config_bar_initialisation == True:
            self.selected_indices['QA'] = [self.qa_default_inds]

        #unset variable to allow interactive handling from now
        self.block_config_bar_handling_updates = False

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#
    
    def get_valid_obs_files_in_date_range(self):
        
        '''define function that iterates through observational dictionary tree and returns a dictionary of available data in the selected date range'''

        #create dictionary to store available observational data
        self.available_observation_data = {}

        #iterate through networks
        for network in list(self.all_observation_data.keys()):
            #iterate through resolutions
            for resolution in list(self.all_observation_data[network].keys()):
                #iterate through matrices
                for matrix in list(self.all_observation_data[network][resolution].keys()):
                    #iterate through species
                    for species in list(self.all_observation_data[network][resolution][matrix].keys()):
                        #get all file yearmonths associated with species
                        species_file_yearmonths = self.all_observation_data[network][resolution][matrix][species]
                        #get file yearmonths within date range
                        valid_species_files_yearmonths = [ym for ym in species_file_yearmonths if (ym >= self.selected_start_date_firstdayofmonth) & (ym < self.selected_end_date)]
                        if len(valid_species_files_yearmonths) > 0:
                            #if network not in dictionary yet, add it
                            if network not in list(self.available_observation_data.keys()):
                                self.available_observation_data[network] = {}
                            #if resolution not in dictionary yet, add it
                            if resolution not in list(self.available_observation_data[network].keys()):
                                self.available_observation_data[network][resolution] = {}
                            #if matrix not in dictionary yet, add it
                            if matrix not in list(self.available_observation_data[network][resolution].keys()):
                                self.available_observation_data[network][resolution][matrix] = {}
                            #add species with associated list of file start yearmonths
                            self.available_observation_data[network][resolution][matrix][species] = valid_species_files_yearmonths
                            
    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#

    def config_bar_params_change_handler(self, changed_param):

        '''define function which handles interactive updates of combo box fields'''

        if (changed_param != '') & (self.block_config_bar_handling_updates == False):
        
            #get event origin source
            event_source = self.sender()
        
            #if network, resolution, matrix or species have changed then respective current selection for the changed param
            if event_source == self.cb_network:
                self.selected_network = changed_param
            elif event_source == self.cb_resolution:
                self.selected_resolution = changed_param
            elif event_source == self.cb_matrix:
                self.selected_matrix = changed_param
                self.selected_species = sorted(list(self.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()][self.cb_matrix.currentText()].keys()))[0]
            elif event_source == self.cb_species:
                self.selected_species = changed_param

            #set variable to check if date range changes
            self.date_range_has_changed = False

            #if have start date/end date have changed, make sure both have 8 characters (YYYYMMDD), are both numbers, and that end date is > start_date, before doing update of selection/fields
            if (event_source == self.le_start_date) or (event_source == self.le_end_date):
                valid_date = False 
                selected_start_date = self.le_start_date.text()
                selected_end_date = self.le_end_date.text()
                if (len(selected_start_date) == 8) & (len(selected_end_date) ==  8):
                    if (selected_start_date.isdigit() == True) & (selected_end_date.isdigit() == True):
                        self.date_range_has_changed = True
                        self.selected_start_date = int(selected_start_date)
                        self.selected_end_date = int(selected_end_date)
                        self.selected_start_date_firstdayofmonth = int(str(self.selected_start_date)[:6]+'01')
                    else:
                        return
                    
            #update configuration bar fields
            self.update_configuration_bar_fields()  

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#
    #define functions which generate pop up configuration windows for some fields

    def handle_pop_up_experiments_window(self):

        #gather available experiment data for selected network/resolution/species
        #create dictionary storing available experiment-grid names associated with valid files in set date range
        self.available_experiment_data = {}

        #get all different experiment names
        available_experiments = os.listdir('%s'%(exp_root))          

        #iterate through available experiments
        for experiment in available_experiments:      
 
            #get all available grid types by experiment 
            available_grids = os.listdir('%s/%s'%(exp_root,experiment))

            #iterate through all available grids
            for grid in available_grids:
            
                #test first if interpolated directory exists before trying to get files from it
                #if it does not exit, continue
                if not os.path.exists('%s/%s/%s/%s/%s/%s'%(exp_root,experiment,grid,self.selected_resolution,self.selected_species,self.selected_network)):       
                    continue
                else:
                    #get all experiment netCDF files by experiment/grid/selected resolution/selected species/selected network
                    network_files = os.listdir('%s/%s/%s/%s/%s/%s'%(exp_root,experiment,grid,self.selected_resolution,self.selected_species,self.selected_network))
                    #get start YYYYMM yearmonths of data files
                    network_files_yearmonths = [int(f.split('_')[-1][:6]+'01') for f in network_files] 
                    #limit data files to just those within date range
                    valid_network_files_yearmonths = [ym for ym in network_files_yearmonths if (ym >= self.selected_start_date_firstdayofmonth) & (ym < self.selected_end_date)]
                    #if have some valid data files for experiment-grid, add experiment grid (with associated yearmonths) to dictionary
                    if len(valid_network_files_yearmonths) > 0:
                        self.available_experiment_data['%s-%s'%(experiment,grid)] = valid_network_files_yearmonths

        #get list of available experiment-grid names
        self.available_experiment_grids = np.array(sorted(list(self.available_experiment_data.keys())))

        #update selected indices for experiments  
        #for experiments, keep previously selected values selected if still available
        #update previously selected experiments variable
        if len(self.selected_indices['EXPERIMENTS'][0]) > 0:
            previous_selected_experiments = self.previous_available_experiment_grids[self.selected_indices['EXPERIMENTS'][0]]
        else:
            previous_selected_experiments = []
        #set selected indices as previously selected indices in current available list of experiments
        selected_experiments = [previous_selected_experiment for previous_selected_experiment in previous_selected_experiments if previous_selected_experiment in self.available_experiment_grids]
        selected_experiment_inds = np.array([np.where(self.available_experiment_grids == selected_experiment)[0][0] for selected_experiment in selected_experiments], dtype=np.uint8)
        self.selected_indices['EXPERIMENTS'] = [selected_experiment_inds]   

        #set previous available experiments variable   
        self.previous_available_experiment_grids = np.array(self.available_experiment_grids)

        #setup pop up window
        self.experiments_window = pop_up_window(window_type='EXPERIMENTS', window_titles=['Select Experiment/s'], checkbox_labels=[self.available_experiment_grids], default_checkbox_selection=[[]], selected_indices=self.selected_indices)

    def handle_pop_up_flags_window(self):
        #setup pop up window
        self.qa_window = pop_up_window(window_type='FLAGS', window_titles=['Select standardised data reporter provided flags to filter by'], checkbox_labels=[self.flag_names], default_checkbox_selection=[self.flag_default_inds], selected_indices=self.selected_indices)

    def handle_pop_up_qa_window(self):
        #setup pop up window
        self.qa_window = pop_up_window(window_type='QA', window_titles=['Select standardised QA flags to filter by'], checkbox_labels=[self.qa_names], default_checkbox_selection=[self.qa_default_inds], selected_indices=self.selected_indices)

    def handle_pop_up_classifications_window(self):
        #setup pop up window
        self.qa_window = pop_up_window(window_type='CLASSIFICATIONS', window_titles=['Select standardised classifications to retain','Select standardised classifications to remove'], checkbox_labels=[self.classification_names,self.classification_names], default_checkbox_selection=[self.classification_default_inds_to_retain,self.classification_default_inds_to_remove], selected_indices=self.selected_indices)

    def handle_pop_up_methods_window(self):
        #only proceed if have some valid stations in memory
        if len(self.station_references) > 0:
            #setup pop up window
            self.qa_window = pop_up_window(window_type='METHODS', window_titles=['Select standardised measurement methodologies to retain'], checkbox_labels=[self.station_unique_methods], default_checkbox_selection=[np.arange(len(self.station_unique_methods), dtype=np.int)], selected_indices=self.selected_indices)

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#
    def handle_data_selection_update(self):

        '''define function which handles update of data selection and MPL canvas upon pressing of READ button'''

        #if have no data to read, then do not read any data
        if self.no_data_to_read == True:
            return

        #Update mouse cursor to a waiting cursor
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        #set variable that blocks updating of MPL canvas until all data has been updated
        self.block_MPL_canvas_updates = True

        #set all previously active variables as past_active variables
        self.previous_active_network = self.active_network
        self.previous_active_resolution = self.active_resolution
        self.previous_active_matrix = self.active_matrix
        self.previous_active_species = self.active_species  
        self.previous_active_start_date = self.active_start_date
        self.previous_active_end_date = self.active_end_date
        self.previous_active_experiment_grids = self.active_experiment_grids
        self.previous_active_qa_inds = self.active_qa_inds
        self.previous_active_flag_inds = self.active_flag_inds
        self.previous_active_classifications_to_retain_inds = self.active_classifications_to_retain_inds
        self.previous_active_classifications_to_remove_inds = self.active_classifications_to_remove_inds
        
        #set all currently selected variables as active variables
        self.active_network = self.selected_network
        self.active_resolution = self.selected_resolution
        self.active_matrix = self.selected_matrix
        self.active_species = self.selected_species
        self.active_start_date = self.selected_start_date
        self.active_end_date = self.selected_end_date
        if len(self.selected_indices['EXPERIMENTS'][0]) > 0:
            self.active_experiment_grids = self.available_experiment_grids[self.selected_indices['EXPERIMENTS'][0]]
        else:
            self.active_experiment_grids = []
        self.active_qa_inds = self.selected_indices['QA'][0]
        self.active_flag_inds = self.selected_indices['FLAGS'][0]
        self.active_classifications_to_retain_inds = self.selected_indices['CLASSIFICATIONS'][0]
        self.active_classifications_to_remove_inds = self.selected_indices['CLASSIFICATIONS'][1]
        
        #--------------------------------------------------------------------#
        #determine what data (if any) needs to be read

        #set variables that inform what data needs to be read (set all initially as False)
        read_all = False
        read_left = False
        read_right = False
        cut_left = False
        cut_right = False
        
        #determine if any of the key variables have changed 
        #(network, resolution, species, qa, flags, classifications_to_retain, classifications_to_remove)
        #if any have changed, observations and any selected experiments have to be re-read entirely
        if (self.active_network != self.previous_active_network) or (self.active_resolution != self.previous_active_resolution) or (self.active_species != self.previous_active_species) or (np.array_equal(self.active_qa_inds,self.previous_active_qa_inds) == False) or (np.array_equal(self.active_flag_inds,self.previous_active_flag_inds) == False) or (np.array_equal(self.active_classifications_to_retain_inds,self.previous_active_classifications_to_retain_inds) == False) or (np.array_equal(self.active_classifications_to_remove_inds,self.previous_active_classifications_to_remove_inds) == False):
            read_all = True
        #key variables have not changed, has start/end date?
        else:
            #determine if start date/end date have changed
            if (self.active_start_date != self.previous_active_start_date) or (self.active_end_date != self.previous_active_end_date):
                #if date range has changed then determine type of overlap with previous date range
                #no overlap (i.e. start date >= than previous end date, or end date <= than previous start date)?
                if (self.active_start_date >= self.previous_active_end_date) or (self.active_end_date <= self.previous_active_start_date):
                    read_all = True
                #data range fully inside previous data range (i.e. start date later and end date earlier)?
                elif (self.active_start_date > self.previous_active_start_date) & (self.active_end_date < self.previous_active_end_date):
                    cut_left = True
                    cut_right = True
                #need to read data on left edge and right edge of previous date range (i.e. start date earlier and end date later)?
                elif (self.active_start_date < self.previous_active_start_date) & (self.active_end_date > self.previous_active_end_date):
                    read_left = True
                    read_right = True
                #need to read data on left edge and cut on right edge of previous date range (i.e. start date earlier and end date earlier)?
                elif (self.active_start_date < self.previous_active_start_date) & (self.active_end_date < self.previous_active_end_date):
                    read_left =  True
                    cut_right = True
                #need to cut data on left edge and read data on right edge of previous date range (i.e. start date later and end date later)?
                elif (self.active_start_date > self.previous_active_start_date) & (self.active_end_date > self.previous_active_end_date): 
                    cut_left = True
                    read_right = True
                #need to read data on left edge of previous date range (i.e. start date earlier)?
                elif (self.active_start_date < self.previous_active_start_date):
                    read_left = True
                #need to read data on right edge of previous date range (i.e. end date later)?
                elif (self.active_end_date > self.previous_active_end_date):
                    read_right = True
                #need to cut data on left edge of previous date range (i.e. start date later)?
                elif (self.active_start_date > self.previous_active_start_date):
                    cut_left = True
                #need to cut data on right edge of previous date range (i.e. end date earlier)?
                elif (self.active_end_date < self.previous_active_end_date):
                    cut_right = True 

        #---------------------------------------#

        #determine if any of the active experiments have changed 
        #remove experiments that are no longer selected from data_in_memory dictionary
        experiments_to_remove = [experiment for experiment in self.previous_active_experiment_grids if experiment not in self.active_experiment_grids]
        for experiment in experiments_to_remove:
            del self.data_in_memory[experiment]
        #any new experiments will need completely re-reading
        experiments_to_read = [experiment for experiment in self.active_experiment_grids if experiment not in self.previous_active_experiment_grids]
        
        #has date range changed? 
        if (read_all == True) or (read_left == True) or (read_right == True) or (cut_left == True) or (cut_right == True):         

            #set new active time array/unique station references/longitudes/latitudes
            #adjust data arrays to account for potential changing number of stations 
            self.read_setup()

            #need to re-read all observations/experiments?
            if read_all == True:
                #reset data in memory dictionary
                self.data_in_memory = {}
                #read observations
                self.read_data('observations', self.active_start_date, self.active_end_date)    
                #read selected experiments (iterate through)
                for data_label in self.active_experiment_grids:
                    self.read_data(data_label, self.active_start_date, self.active_end_date)  
                    #if experiment in experiments_to_read list, remove it (as no longer need to read it)
                    if data_label in experiments_to_read:
                        experiments_to_read.remove(data_label)    
            else:
                #if station references array has changed then as cutting/appending to existing data need to rearrange existing data arrays accordingly  
                if np.array_equal(self.previous_station_references, self.station_references) == False:  
                    #get indices of stations in previous station references array in current station references array
                    old_station_inds = np.where(np.in1d(self.previous_station_references, self.station_references))[0]
                    #get indices of stations in current station references array that were in previous station references array
                    new_station_inds = np.where(np.in1d(self.station_references, self.previous_station_references))[0]
                    
                    #iterate through all keys in data in memory dictionary
                    for data_label in list(self.data_in_memory.keys()):
                        #create new data array in shape of current station references array, putting the old data into new array in the correct positions
                        new_data_array = np.full((len(self.previous_time_array),len(self.station_references)), np.NaN, dtype=np.float32)
                        new_data_array[:,new_station_inds] = self.data_in_memory[data_label]['data'][:,old_station_inds]
                        #overwrite data array with reshaped version
                        self.data_in_memory[data_label]['data'] = new_data_array
                
            #need to cut edges?
            if (cut_left == True) or (cut_right == True):

                #set default edge limits as current edges
                left_edge_ind = 0
                right_edge_ind = len(self.previous_time_array)

                #need to cut on left data edge?                                                                                      
                if cut_left == True:
                    left_edge_ind = np.where(self.previous_time_array == self.time_array[0])[0][0]
                
                #need to cut on right data edge?
                if cut_right == True:
                    right_edge_ind = np.where(self.previous_time_array == self.time_array[-1])[0][0]+1

                #iterate through all keys in data in memory dictionary and cut edges of the associated arrays appropriately 
                for data_label in list(self.data_in_memory.keys()):
                    self.data_in_memory[data_label]['data'] = self.data_in_memory[data_label]['data'][left_edge_ind:right_edge_ind,:]

            #need to read on left edge?
            if read_left == True:
                #get n number of new elements on left edge
                n_new_left_inds = np.where(self.time_array == self.previous_time_array[0])[0][0]
                #iterate through all keys in data in memory dictionary and insert read data on left edge of the associated arrays 
                for data_label in list(self.data_in_memory.keys()):
                    #add space on left edge to insert new read data
                    self.data_in_memory[data_label]['data'] = np.insert(self.data_in_memory[data_label]['data'], 0, np.full((n_new_left_inds,len(self.station_references)),np.NaN), axis=0)
                    self.read_data(data_label, self.active_start_date, self.previous_active_start_date)

            #need to read on right edge?
            if read_right == True:
                #get n number of new elements on right edge
                n_new_right_inds = (len(self.time_array) - 1) - np.where(self.time_array == self.previous_time_array[-1])[0][0]
                #iterate through all keys in data in memory dictionary and insert read data on right edge of the associated arrays
                for data_label in list(self.data_in_memory.keys()):
                    self.data_in_memory[data_label]['data'] = np.append(self.data_in_memory[data_label]['data'], np.full((n_new_right_inds,len(self.station_references)),np.NaN), axis=0)
                    self.read_data(data_label, self.previous_active_end_date, self.active_end_date)
        
        #if have new experiments to read, then read them now
        if len(experiments_to_read) > 0:
            for data_label in experiments_to_read:
                self.read_data(data_label, self.active_start_date, self.active_end_date)    

        #--------------------------------------------------------------------#
        #if species has changed, update default species specific lower/upper limits
        if (self.active_species != self.previous_active_species):
            #update default lower/upper species specific limits and filter data outside limits
            species_lower_limit = np.float32(parameter_dictionary[self.active_species]['extreme_lower_limit'])
            species_upper_limit = np.float32(parameter_dictionary[self.active_species]['extreme_upper_limit'])
            #set default limits
            self.le_minimum_value.setText(str(species_lower_limit))
            self.le_maximum_value.setText(str(species_upper_limit))  

        #--------------------------------------------------------------------#
        #update dictionary of plotting parameters (colour and zorder etc.) for each data array
        self.update_plotting_parameters()
        
        #--------------------------------------------------------------------#
        #run function to filter data outside lower/upper limits, not using desired measurement methods, and < desired minimum data availability
        self.mpl_canvas.handle_data_filter_update() 

        #--------------------------------------------------------------------#
        #update map z combobox fields based on data in memory

        #generate lists of basic and basis+bias statistics for using in the z statistic combobox
        self.basic_z_stats = np.array(list(OrderedDict(sorted(basic_stats_dict.items(), key=lambda x: x[1]['order'])).keys()))
        self.basic_and_bias_z_stats = np.append(self.basic_z_stats, list(OrderedDict(sorted(experiment_bias_stats_dict.items(), key=lambda x: x[1]['order'])).keys()))

        #generate list of sorted z1/z2 data arrays names in memory, putting observations before experiments, and empty string item as first element in z2 array list (for changing from 'difference' statistics to 'absolute')
        if len(list(self.data_in_memory.keys())) == 1:
            self.z1_arrays = np.array(['observations'])
        else:
            data_array_labels = np.array(list(self.data_in_memory.keys()))
            self.z1_arrays = np.append(['observations'], np.delete(data_array_labels, np.where(data_array_labels == 'observations')))
        self.z2_arrays = np.append([''], self.z1_arrays)

        #initialise map z statistic comboboxes
        self.mpl_canvas.handle_map_z_statistic_update()

        #--------------------------------------------------------------------#
        #update experiment bias combobox fields based on data in memory

        #if have no experiment data, all fields are empty
        if len(list(self.data_in_memory.keys())) == 1:
            self.experiment_bias_types = np.array([])
        #else, generate combobox lists
        else:
            #set all experiment bias types
            self.experiment_bias_types = np.array(['Aggregated'])
            
            #initialise experiment bias comboboxes
            self.mpl_canvas.handle_experiment_bias_update()

        #--------------------------------------------------------------------#
        #reset station select checkboxes to be unchecked
        self.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
        self.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
        #--------------------------------------------------------------------#

        #unset variable to allow updating of MPL canvas
        self.block_MPL_canvas_updates = False

        #update MPL canvas
        self.mpl_canvas.update_MPL_canvas()

        #Restore mouse cursor to normal 
        QtWidgets.QApplication.restoreOverrideCursor()

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#   
    
    def read_setup(self):

        '''function that setups key variables for new read of observational/experiment data
           a time array and arrays of unique station references/longitudes/latitudes are created.
        '''        

        #force garbage collection (to avoid memory issues)
        gc.collect()

        #set current time array, as previous time array
        self.previous_time_array = self.time_array    

        #set current station references, as previous station references 
        self.previous_station_references = self.station_references

        #get N time chunks between desired start date and end date to set time array
        if self.active_resolution == 'hourly':
            self.active_frequency_code = 'H'
        elif self.active_resolution == 'daily':
            self.active_frequency_code = 'D'
        elif self.active_resolution == 'monthly':
            self.active_frequency_code = 'MS'
        str_active_start_date = str(self.active_start_date)
        str_active_end_date = str(self.active_end_date)
        self.time_array = pd.date_range(start = datetime.datetime(int(str_active_start_date[:4]), int(str_active_start_date[4:6]), int(str_active_start_date[6:8])),end = datetime.datetime(int(str_active_end_date[:4]), int(str_active_end_date[4:6]), int(str_active_end_date[6:8])), freq = self.active_frequency_code)[:-1]

        #get all relevant observational files
        file_root = '%s/%s/%s/%s/%s_'%(obs_root, self.active_network, self.active_resolution, self.active_species, self.active_species)
        relevant_files = sorted([file_root+str(yyyymm)[:6]+'.nc' for yyyymm in self.available_observation_data[self.active_network][self.active_resolution][self.active_matrix][self.active_species]])       

        #redefine some key variables globally (for access by parallel netCDF reading functions)
        global time_array, active_species, selected_qa, selected_flags, selected_classifications_to_retain, selected_classifications_to_remove
        time_array = self.time_array
        active_species = self.active_species
        selected_qa = self.qa_codes[self.active_qa_inds]
        selected_flags = self.flag_codes[self.active_flag_inds]        
        selected_classifications_to_retain = self.classification_codes[self.active_classifications_to_retain_inds]    
        selected_classifications_to_remove = self.classification_codes[self.active_classifications_to_remove_inds]  

        #------------------------------------------------------------------------------------#
        #iterate through all relevant observational files and read station references/longitudes/latitudes (either in serial/parallel)

        #define dictionary to store all read metadata (i.e. per file)
        all_read_metadata = {}    

        #define dictionary with metadata variables to read, with associated data types
        metadata_dict = {'station_reference':np.object,
                         'station_name':np.object,
                         'latitude':np.float32, 
                         'longitude':np.float32,
                         'measurement_altitude':np.float32,
                         'country':np.object,
                         'network':np.object,
                         'standardised_network_provided_area_classification':np.object,
                         'standardised_network_provided_station_classification':np.object,
                         'standardised_network_provided_main_emission_source':np.object,
                         'standardised_network_provided_land_use':np.object,
                         'standardised_network_provided_terrain':np.object,
                         'standardised_network_provided_measurement_scale':np.object,
                         'representative_radius':np.float32,
                         'GSFC_coastline_proximity':np.float32,
                         'primary_sampling_type':np.object,
                         'sample_preparation_methodology_types':np.object,
                         'measurement_methodology':np.object,  
                         'measuring_instrument_name':np.object,    
                         'measuring_instrument_sampling_type':np.object
                        }

        #create list of metadata variables to read (make global)
        global metadata_vars_to_read
        metadata_vars_to_read = list(metadata_dict.keys())

        #add all metadata variables to read to station metadata dictionary with associated empty numpy arrays of the appropriate type
        for meta_var in metadata_vars_to_read:
            all_read_metadata[meta_var] = np.array([], dtype=metadata_dict[meta_var])

        #read serially
        if self.read_type == 'serial':
            #iterate through relevant files
            for relevant_file_ii, relevant_file in enumerate(relevant_files):
                file_metadata = read_netCDF_station_information(relevant_file)
                for meta_var in metadata_vars_to_read:
                    all_read_metadata[meta_var] = np.append(all_read_metadata[meta_var], file_metadata[meta_var])
                    
        #read in parallel
        elif self.read_type == 'parallel':

            #setup pool of N workers on N CPUs
            pool = multiprocessing.Pool(n_CPUs)

            #read netCDF files in parallel
            all_file_metadata = pool.map(read_netCDF_station_information, relevant_files)
            #will not submit more files to pool, so close access to it
            pool.close()
            #wait for worker processes to terminate before continuing
            pool.join()

            #iterate through relevant file data
            for file_metadata in all_file_metadata:
                for meta_var in metadata_vars_to_read:
                    all_read_metadata[meta_var] = np.append(all_read_metadata[meta_var], file_metadata[meta_var])

        #define dictionary to store read metadata, by station, by metadata variable
        self.station_metadata = {}

        #get unique sorted station references across all files (make global, and also add to self)
        global station_references
        station_references, unique_indices = np.unique(all_read_metadata['station_reference'], return_index=True)
        self.station_references = station_references  

        #get standard measurement methodology per station
        self.station_methods = np.array([ref.split('_')[-1] for ref in self.station_references]) 

        #for latitude/longitude/measurement altitude/GSFC coastline proximity variables, set static metadata variables (taking the first available value per station in the given time window)
        self.station_longitudes = all_read_metadata['longitude'][unique_indices]
        self.station_latitudes = all_read_metadata['latitude'][unique_indices]
        self.station_measurement_altitudes = all_read_metadata['measurement_altitude'][unique_indices] 
        self.station_GSFC_coastline_proximities = all_read_metadata['GSFC_coastline_proximity'][unique_indices] 

        #iterate through each unique station
        for station_reference in self.station_references:
            station_inds = np.where(all_read_metadata['station_reference'] == station_reference)[0]
            #create empty dictionary ro store station specific metadata
            self.station_metadata[station_reference] = {}
            #iterate through read metadata variables    
            for meta_var in metadata_vars_to_read:        
                #get all unique metadata for variable, by station
                unique_metadata = np.unique(all_read_metadata[meta_var][station_inds])
                #append all and unique metadata by station
                self.station_metadata[station_reference][meta_var] = {'all':all_read_metadata[meta_var][station_inds],'unique':unique_metadata}
                
        #update measurement units for species (take standard units from parameter dictionary)
        self.measurement_units = parameter_dictionary[self.active_species]['standard_units']

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def read_data(self, data_label, start_date_to_read, end_date_to_read):
 
        '''function that handles reading of observational/experiment data'''

        #force garbage collection (to avoid memory issues)
        gc.collect()
 
        #set global process_type variable (for access by parallel read function)
        #also get relevant file start dates 
        global process_type
        if data_label == 'observations':
            process_type = 'observations'
            file_root = '%s/%s/%s/%s/%s_'%(obs_root, self.active_network, self.active_resolution, self.active_species, self.active_species)
            relevant_file_start_dates = sorted(self.available_observation_data[self.active_network][self.active_resolution][self.active_matrix][self.active_species])  

        else:
            process_type = 'experiment'
            experiment_grid_split = data_label.split('-')
            active_experiment = experiment_grid_split[0]
            active_grid = experiment_grid_split[1]  
            file_root = '%s/%s/%s/%s/%s/%s/%s_'%(exp_root, active_experiment, active_grid, self.active_resolution, self.active_species, self.active_network, self.active_species)
            relevant_file_start_dates = sorted(self.available_experiment_data[data_label])

        #create list of relevant files to read
        relevant_files = [file_root+str(yyyymm)[:6]+'.nc' for yyyymm in relevant_file_start_dates]     

        #limit data files to required date range to read (i.e. taking care not to re-read what has already been read)
        first_valid_file_ind = bisect.bisect_right(relevant_file_start_dates, int(start_date_to_read))
        if first_valid_file_ind != 0:
            first_valid_file_ind = first_valid_file_ind - 1 
        last_valid_file_ind = bisect.bisect_left(relevant_file_start_dates, int(end_date_to_read))
        if first_valid_file_ind == last_valid_file_ind:
            relevant_files = [relevant_files[first_valid_file_ind]]   
        else:
            relevant_files = relevant_files[first_valid_file_ind:last_valid_file_ind]     
 
        #check if data label in data in memory dictionary
        if data_label not in list(self.data_in_memory.keys()):  
            #if not create empty array (filled with NaNs) to store species data and place it in the dictionary
            self.data_in_memory[data_label] = {'data': np.full((len(self.time_array),len(self.station_references)), np.NaN, dtype=np.float32)}
        
            #if process_type is experiment, get experiment specific grid edges from first relevant file, and save to data in memory dictionary
            if process_type == 'experiment':
                exp_nc_root = Dataset(relevant_files[0])
                self.data_in_memory[data_label]['grid_edge_longitude'] = exp_nc_root['grid_edge_longitude'][:]
                self.data_in_memory[data_label]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
                exp_nc_root.close()

        #------------------------------------------------------------------------------------#
        #iterate and read species data in all relevant netCDF files (either in serial/parallel)
 
        #read serially
        if self.read_type == 'serial':
 
            #iterate through relevant netCDF files
            for relevant_file in relevant_files: 
                #read file
                file_data, time_indices, full_array_station_indices = read_netCDF_data(relevant_file)
                #place read data into big array as appropriate
                self.data_in_memory[data_label]['data'][time_indices[:,np.newaxis], full_array_station_indices[np.newaxis,:]] = file_data
     
        #read in parallel 
        elif self.read_type == 'parallel':

            #setup pool of N workers on N CPUs
            pool = multiprocessing.Pool(n_CPUs)

            #read netCDF files in parallel
            all_file_data = pool.map(read_netCDF_data, relevant_files)
            #will not submit more files to pool, so close access to it
            pool.close()
            #wait for worker processes to terminate before continuing
            pool.join()

            #iterate through read file data and place data into data array as appropriate
            for file_data in all_file_data:
                self.data_in_memory[data_label]['data'][file_data[1][:,np.newaxis], file_data[2][np.newaxis,:]] = file_data[0]

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#    

    def update_plotting_parameters(self):

        '''function that updates plotting parameters (colour and zorder) for each selected data array'''
    
        #assign a colour/zorder to all selected data arrays 

        #define observations colour to be 'black'
        self.data_in_memory['observations']['colour'] = 'black'
        #define zorder of observations to be 5
        self.data_in_memory['observations']['zorder'] = 5
        
        #generate a list of RGB tuples for number of experiments there are
        sns.reset_orig()
        clrs = sns.color_palette('husl', n_colors=len(list(self.data_in_memory.keys()))-1)

        #iterate through sorted experiment names, assigning each experiment a new RGB colour tuple, and zorder
        experiment_ind = 1
        for experiment in sorted(list(self.data_in_memory.keys())):
            if experiment != 'observations':
                #define colour for experiment
                self.data_in_memory[experiment]['colour'] = clrs[experiment_ind-1]
                #define zorder for experiment (obs zorder + experiment_ind)
                self.data_in_memory[experiment]['zorder'] = self.data_in_memory['observations']['zorder'] + experiment_ind
                #update count of experiments
                experiment_ind +=1
#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#

def read_netCDF_station_information(relevant_file):

    '''function that handles reading of observational desired station metadata in a netCDF file, returning a dictionary with read metadata'''

    #read netCDF frame
    nCDF_root = Dataset(relevant_file) 

    #read all desired metadata, placing it within a dictionary by variable name
    read_metadata = {}
    for meta_var in metadata_vars_to_read:
        read_metadata[meta_var] = nCDF_root[meta_var][:]

    #close netCDF
    nCDF_root.close()

    return read_metadata

#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#

def read_netCDF_data(relevant_file):

    '''function that handles reading of observational/experiment netCDF data
       also handles filtering of observational data based on selected qa/flag/classification flags
    '''

    #read netCDF frame
    nCDF_root = Dataset(relevant_file) 

    #get time units
    time_units = nCDF_root['time'].units

    #get file time (handle monthly resolution data differently to hourly/daily as num2date does not support 'months since' units)
    if 'months' in time_units:
        monthly_start_date = time_units.split(' ')[2]
        file_time = pd.date_range(start=monthly_start_date, periods=1, freq='MS')
    else:
        file_time = num2date(nCDF_root['time'][:], time_units)
        #remove microseconds
        file_time = pd.to_datetime([t.replace(microsecond=0) for t in file_time])

    #get valid file time indices (i.e. those times in active full time array)
    valid_file_time_indices = np.array([i for i, val in enumerate(file_time) if (val >= time_array[0]) & (val <= time_array[-1])], dtype=np.int)
    #cut file time for valid indices 
    file_time = file_time[valid_file_time_indices]
    #get indices relative to active full time array
    full_array_time_indices = np.searchsorted(time_array, file_time)
    
    #get all station references in file
    file_station_references = nCDF_root['station_reference'][:]
    #get indices of all unique station references that are contained within file station references array
    full_array_station_indices = np.where(np.in1d(station_references, file_station_references))[0]
    #get indices of file station station references that are contained in all unique station references array
    current_file_station_indices = np.where(np.in1d(file_station_references, station_references))[0]

    #read in species data
    file_data = nCDF_root[active_species][valid_file_time_indices,:]
    #get masked data 
    data_mask = file_data.mask
    #set masked data as NaN
    file_data[data_mask] = np.NaN    

    #for observations, set species data based on selected qa flags/standard data provider flags/classifications to retain or remove as NaN
    if process_type == 'observations':
        #if some qa flags selected then screen
        if len(selected_qa) > 0:
            #screen out observations which are associated with any of the selected qa flags
            file_data[np.isin(nCDF_root['qa'][valid_file_time_indices,:,:], selected_qa).any(axis=2)] = np.NaN
        #if some data provider flags selected then screen
        if len(selected_flags) > 0:
            #screen out observations which are associated with any of the selected data provider flags
            file_data[np.isin(nCDF_root['flag'][valid_file_time_indices,:,:], selected_flags).any(axis=2)] = np.NaN
        #if some classification flags (retain  or remove) selected then screen
        if (len(selected_classifications_to_retain) > 0) or (len(selected_classifications_to_remove) > 0):
            file_classifications = nCDF_root['classification'][valid_file_time_indices,:,:]
            #screen out all observations that aren't associated with all of the selected classifications to retain
            if len(selected_classifications_to_retain) > 0:
                file_data[np.isin(file_classifications, selected_classifications_to_retain, invert=True).all(axis=2)] = np.NaN
            #screen out all observations that are associated with any of the selected classifications to remove 
            if len(selected_classifications_to_remove) > 0:
                file_data[np.isin(file_classifications, selected_classifications_to_remove).any(axis=2)] = np.NaN

    #close netCDF
    nCDF_root.close()

    #return valid species data, time indices relative to active full time array, file station indices relative to all unique station references array 
    return file_data[:,current_file_station_indices], full_array_time_indices, full_array_station_indices   

##-----------------------------------------------------------------------------------------------------------------------------------------------------##
##-----------------------------------------------------------------------------------------------------------------------------------------------------##

class MPL_Canvas(FigureCanvas):

    '''class that handles the creation and updates of a matplotlib canvas, and associated subplots'''

    def __init__(self, read_instance, parent=None):

        '''initialise the MPL canvas'''

        self.figure = Figure(dpi=100)
        self.canvas = FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)
        self.read_instance = read_instance
        
        #--------------------------------------------------#
        #setup gridding of canvas
        self.gridspec = GridSpec(100, 100)
        self.gridspec.update(left=0.01,right=0.99,top=0.96,bottom=0.04, wspace=0.00, hspace=0.00)
        
        #map_ax =              gridspec.new_subplotspec((0, 0),   rowspan=45, colspan=45)
        legend_ax =           self.gridspec.new_subplotspec((0, 46),  rowspan=8,  colspan=54)
        ts_ax =               self.gridspec.new_subplotspec((12, 54), rowspan=34, colspan=46)
        violin_hours_ax =     self.gridspec.new_subplotspec((57, 70), rowspan=19, colspan=30)    
        violin_months_ax =    self.gridspec.new_subplotspec((81, 70), rowspan=19, colspan=18)
        violin_days_ax =      self.gridspec.new_subplotspec((81, 91), rowspan=19, colspan=9)
        exp_bias_hours_ax =   self.gridspec.new_subplotspec((57, 35), rowspan=19, colspan=30)    
        exp_bias_months_ax =  self.gridspec.new_subplotspec((81, 35), rowspan=19, colspan=18)
        exp_bias_days_ax =    self.gridspec.new_subplotspec((81, 56), rowspan=19, colspan=9)
        station_metadata_ax = self.gridspec.new_subplotspec((55, 0),  rowspan=45, colspan=30)
        
        #create subplot axes on grid
        #self.map_ax =              self.figure.add_subplot(map_ax)
        self.legend_ax =           self.figure.add_subplot(legend_ax)
        self.ts_ax =               self.figure.add_subplot(ts_ax)
        self.violin_hours_ax =     self.figure.add_subplot(violin_hours_ax)
        self.violin_months_ax =    self.figure.add_subplot(violin_months_ax)
        self.violin_days_ax =      self.figure.add_subplot(violin_days_ax)
        self.exp_bias_hours_ax =   self.figure.add_subplot(exp_bias_hours_ax)
        self.exp_bias_months_ax =  self.figure.add_subplot(exp_bias_months_ax)
        self.exp_bias_days_ax =    self.figure.add_subplot(exp_bias_days_ax)
        self.station_metadata_ax = self.figure.add_subplot(station_metadata_ax)
        
        #map colorbar ax 
        self.cbar_ax = self.figure.add_axes([0.02, 0.525, 0.25, 0.0175])    

        #--------------------------------------------------# 
        #set map variables 
        #create variable to create map axis on first data read
        self.map_initialised = False

        #define projections for map plot and actual geographic coordinates 
        self.datacrs = ccrs.PlateCarree()
        self.plotcrs = ccrs.Robinson()
        
        #--------------------------------------------------# 
        #turning off specific spines of time series axis
        self.ts_ax.spines["top"].set_visible(False)      
        self.ts_ax.spines["right"].set_visible(False)   

        #--------------------------------------------------#
        #setup interactive picker/lasso on map
        #self.figure.canvas.mpl_connect('pick_event', self.on_click)
        #self.lasso = LassoSelector(self.map_ax, onselect=self.onlassoselect, useblit=True, lineprops=dict(alpha=0.5, color='hotpink', linewidth=1))
        #initialise variable that informs whether to use picker/lasso for updates
        #self.map_already_updated = False
        #initialise variable of valid station indices plotted on map as empty list
        #self.active_map_valid_station_inds = np.array([],dtype=np.int)

        #--------------------------------------------------#
        #hide all axes  
        self.cbar_ax.axis('off')
        self.legend_ax.axis('off')
        self.ts_ax.axis('off')
        self.violin_hours_ax.axis('off')
        self.violin_months_ax.axis('off')
        self.violin_days_ax.axis('off')
        self.exp_bias_hours_ax.axis('off')
        self.exp_bias_months_ax.axis('off')
        self.exp_bias_days_ax.axis('off')
        self.station_metadata_ax.axis('off')
            
    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#

    def update_MPL_canvas(self):

        '''function that updates MPL canvas upon clicking the 'READ' button, and when colocating data'''
            
        #clear all axes (except map)
        self.legend_ax.cla()
        self.ts_ax.cla()
        self.violin_hours_ax.cla()
        self.violin_months_ax.cla()
        self.violin_days_ax.cla()
        self.exp_bias_hours_ax.cla()
        self.exp_bias_months_ax.cla()
        self.exp_bias_days_ax.cla()
        self.station_metadata_ax.cla()

        #hide all axes (except map)
        self.legend_ax.axis('off')
        self.ts_ax.axis('off')
        self.violin_hours_ax.axis('off')
        self.violin_months_ax.axis('off')
        self.violin_days_ax.axis('off')
        self.exp_bias_hours_ax.axis('off')
        self.exp_bias_months_ax.axis('off')
        self.exp_bias_days_ax.axis('off')
        self.station_metadata_ax.axis('off')

        #add map axis if map not yet initialised (i.e. first time loading some data)
        if self.map_initialised == False:
    
            #create map axis
            map_ax = self.gridspec.new_subplotspec((0, 0),   rowspan=45, colspan=45)
            self.map_ax = self.figure.add_subplot(map_ax, projection=self.plotcrs)
            #set map extents
            self.map_ax.set_global() 

            #add coastlines and
            self.map_ax.add_feature(cfeature.LAND, facecolor='0.85')
            self.map_ax.gridlines(linestyle='-', alpha=0.4)
            
            #reset the navigation toolbar stack for the map axis with the current view limits
            self.reset_ax_navigation_toolbar_stack(self.map_ax)

            #setup interactive picker/lasso on map
            self.figure.canvas.mpl_connect('pick_event', self.on_click)
            self.lasso = LassoSelector(self.map_ax, onselect=self.onlassoselect, useblit=True, lineprops=dict(alpha=0.5, color='hotpink', linewidth=1))
            #initialise variable that informs whether to use picker/lasso for updates
            self.map_already_updated = False
            
            #initialise variable of valid station indices plotted on map as empty list
            self.active_map_valid_station_inds = np.array([],dtype=np.int)

            #update variable to indicate map is now initialised
            self.map_initialised = True
            
        #define all temporal aggregation resolutions that will be used to aggregate data (variable by temporal resolution of data in memory)
        if self.read_instance.active_resolution == 'hourly':
            self.temporal_aggregation_resolutions = ['hour','dayofweek','month']
        elif self.read_instance.active_resolution == 'daily':
            self.temporal_aggregation_resolutions = ['dayofweek','month']
        elif self.read_instance.active_resolution == 'monthly':
            self.temporal_aggregation_resolutions = ['month']

        #reset relative index lists of selected station on map as empty lists
        self.previous_relative_selected_station_inds = np.array([],dtype=np.int) 
        self.relative_selected_station_inds = np.array([],dtype=np.int) 
        self.absolute_selected_station_inds = np.array([],dtype=np.int)

        #calculate map z statistic (for selected z statistic) --> updating active map valid station indices
        self.calculate_z_statistic()

        #update plotted map z statistic
        self.update_map_z_statisitic()

        #plot experiment grid domain edges on map 
        self.update_experiment_grid_domain_edges()

        #update legend
        self.update_legend()

        #draw changes
        self.draw()

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def reset_ax_navigation_toolbar_stack(self, ax):

        '''function which resets the navigation toolbar stack for a given axis with the current view limits'''

        #check if have axes dictionaries in stack list
        if len(self.read_instance.navi_toolbar._nav_stack) == 0:
            #if don't have an axes dictionary in stack list, create one with current axis in dictionary, with current view limits
            self.read_instance.navi_toolbar._nav_stack.push(WeakKeyDictionary({ax: (ax._get_view(), (ax.get_position(True).frozen(),ax.get_position().frozen()))}))
        
        #if have existing axes dictionaries in stack list, iterate through stack list removing given axis from all stack list dictionaries
        else:   
            for axes_dict in self.read_instance.navi_toolbar._nav_stack:
                if ax in axes_dict.keyrefs():    
                    axes_dict.pop(ax)
        
            #now add axis to first dictionary in stack, with the current view limits
            self.read_instance.navi_toolbar._nav_stack[0][ax] = (ax._get_view(), (ax.get_position(True).frozen(),ax.get_position().frozen()))

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def handle_data_filter_update(self):

        '''function which handles updates data filtering by selected lower/upper limit bounds, selected measurement methods and selected minimum data availability %'''

        #get selected variables for minimum data availability, lower/upper limits and measurement methods
        selected_minimum_data_availability_percent = self.read_instance.le_minimum_data_availability.text()
        selected_lower_limit = self.read_instance.le_minimum_value.text()
        selected_upper_limit = self.read_instance.le_maximum_value.text()

        #check selected minimum data availability percent, selected lower limit, selected upper limits are numbers
        try:
            selected_minimum_data_availability_percent = np.float32(selected_minimum_data_availability_percent)
            selected_lower_limit = np.float32(selected_lower_limit)
            selected_upper_limit = np.float32(selected_upper_limit)
        #if any of the fields are not numbers, return from function
        except ValueError:
            return

        #if selected minimum data availability percent is < 0.0%, force it to be 0.0%
        if selected_minimum_data_availability_percent < 0.0:
            selected_minimum_data_availability_percent = 0.0
        #if selected minimum data availability percent is > 100.0%, force it to be 100.0%
        elif selected_minimum_data_availability_percent > 100.0:
            selected_minimum_data_availability_percent = 100.0 

        #reset data arrays to be un-filtered
        self.read_instance.data_in_memory_filtered = copy.deepcopy(self.read_instance.data_in_memory)

        #filter all observational data out of bounds of lower/upper limits
        inds_out_of_bounds = np.logical_or(self.read_instance.data_in_memory_filtered['observations']['data']<selected_lower_limit, self.read_instance.data_in_memory_filtered['observations']['data']>selected_upper_limit)
        self.read_instance.data_in_memory_filtered['observations']['data'][inds_out_of_bounds] = np.NaN  

        #colocate data (if necessary) 
        self.colocate_data()

        #get intersect of indices of stations with >= % minimum data availability percent, and with > 1 valid measurements ,in all observational data arrays (colocated and non-colocated)
        #then subset these indices with standard methods == selected methods, 
        #iterate through all data arrays
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()): 

            #check if data array is an observational data array
            if data_label.split('_')[0] == 'observations':
        
                #calculate data availability fraction per station in observational data array
                station_data_availability_percent = calculate_data_availability_fraction(self.read_instance.data_in_memory_filtered[data_label]['data'])
                #get indices of stations with >= selected_minimum_data_availability
                valid_station_indices_percent = np.arange(len(self.read_instance.station_references), dtype=np.int)[station_data_availability_percent >= selected_minimum_data_availability_percent]

                #get absolute data availability number per station in observational data array
                station_data_availability_number = calculate_data_availability_number(self.read_instance.data_in_memory_filtered[data_label]['data'])                
                #get indices of stations with > 1 available measurements
                valid_station_indices_absolute = np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]

                #get indices of valid stations in intersect of valid_station_indices_percent and valid_station_indices_absolute
                valid_station_indices_availability = np.intersect1d(valid_station_indices_percent, valid_station_indices_absolute)

                #get unique standard measurement methodologies across stations 
                self.read_instance.previous_station_unique_methods = copy.deepcopy(self.read_instance.station_unique_methods)
                self.read_instance.station_unique_methods = np.unique(self.read_instance.station_methods[valid_station_indices_availability])
                #if are reading new data into memory and unique methods have changed from previous, update all methods to be checked by default
                if (np.array_equal(self.read_instance.previous_station_unique_methods, self.read_instance.station_unique_methods) == False) & (self.read_instance.block_MPL_canvas_updates == True):
                    self.read_instance.selected_indices['METHODS'] = [np.arange(len(self.read_instance.station_unique_methods), dtype=np.int)]

                #get indices of subset stations which use checked standard methodologies
                checked_methods = self.read_instance.station_unique_methods[self.read_instance.selected_indices['METHODS'][0]]
                valid_station_indices = valid_station_indices_availability[np.isin(self.read_instance.station_methods[valid_station_indices_availability], checked_methods)]

                #save valid station indices with data array
                self.read_instance.data_in_memory_filtered[data_label]['valid_station_inds'] = valid_station_indices

        #write valid station indices calculated for observations across to associated experimental data arrays
        #iterate through all data arrays
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()): 
        
            #check if data array is not an observational data array
            if data_label.split('_')[0] != 'observations':
        
                #handle colocated experimental arrays
                if '_colocatedto_' in data_label:
                    exp_name = data_label.split('_colocatedto_')[0]
                    self.read_instance.data_in_memory_filtered[data_label]['valid_station_inds'] = copy.deepcopy(self.read_instance.data_in_memory_filtered['observations_colocatedto_%s'%(exp_name)]['valid_station_inds'])
                #handle non-located experimental arrays
                else:
                    self.read_instance.data_in_memory_filtered[data_label]['valid_station_inds'] = copy.deepcopy(self.read_instance.data_in_memory_filtered['observations']['valid_station_inds'])

        #after subsetting by pre-written associated observational valid stations, get indices of stations with > 1 valid measurements in all experiment data arrays (colocated and non-colocated)
        #iterate through all data arrays
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()): 
    
            #check if data array is not an observational data array
            if data_label.split('_')[0] != 'observations':
                #get indices of associated observational data array valid stations (pre-written to experiment data arrays)
                valid_station_inds = self.read_instance.data_in_memory_filtered[data_label]['valid_station_inds']
                #get absolute data availability number per station in experiment data array, after subsetting valid observational stations (i.e. number of non-NaN measurements)
                station_data_availability_number = calculate_data_availability_number(self.read_instance.data_in_memory_filtered[data_label]['data'][:,valid_station_inds])
                #get indices of stations with > 1 available measurements
                valid_station_inds = valid_station_inds[np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]]
                #overwrite previous written valid station indices (now at best a subset of those indices)
                self.read_instance.data_in_memory_filtered[data_label]['valid_station_inds'] = valid_station_inds       
            
        #update plotted map z statistic (if necessary)
        if self.read_instance.block_MPL_canvas_updates == False:
            #calculate map z statistic (for selected z statistic) --> updating active map valid station indices
            self.calculate_z_statistic()

            #update plotted map z statistic
            self.update_map_z_statisitic()

            #if selected stations have changed from previous selected, update associated plots
            if np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds) == False:
                self.update_associated_selected_station_plots()

            #draw changes
            self.draw()

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 
  
    def colocate_data(self):

        '''define function which colocates observational and experiment data'''

        #check if colocation is active or not (to inform subsequent plotting functions whether to use colocated data/or not)
        check_state = self.read_instance.ch_colocate.checkState()
        #update variable to inform plotting functions whether to use colocated data/or not 
        if check_state == QtCore.Qt.Checked: 
            self.colocate_active = True
        else:
            self.colocate_active = False

        #if do not have any experiment data loaded, no colocation is possible, therefore return from function (set colocation to be not active, regardless if the checkbox is ticked or not)
        if len(list(self.read_instance.data_in_memory.keys())) == 1:
            self.colocate_active = False
            return
        else:
            #otherwise, colocate observational and experiment data, creating new data arrays
        
            #colocate observational data array to every different experiment array in memory, and vice versa
            #wherever there is a NaN at one time in one of the observations/experiment arrays, the other array value is also made NaN             

            #get all instances observations are NaN
            nan_obs = np.isnan(self.read_instance.data_in_memory_filtered['observations']['data'])

            #create array for finding instances where have 0 valid values across all experiments
            #initialise as being all True, set as False on the occasion there is a valid value in an experiment 
            exps_all_nan = np.full(nan_obs.shape, True)

            #iterate through experiment data arrays in data in memory dictionary
            for data_label in list(self.read_instance.data_in_memory.keys()):
                if data_label != 'observations':
                    #get all instances experiment are NaN
                    nan_exp = np.isnan(self.read_instance.data_in_memory_filtered[data_label]['data'])
                    #get all instances where either the observational array or experiment array are NaN at a given time
                    nan_instances = np.any([nan_obs,nan_exp],axis=0)
                    #create new observational array colocated to experiment
                    obs_data = copy.deepcopy(self.read_instance.data_in_memory_filtered['observations']['data'])
                    obs_data[nan_instances] = np.NaN
                    self.read_instance.data_in_memory_filtered['observations_colocatedto_%s'%(data_label)] = {'data':obs_data, 'colour':self.read_instance.data_in_memory_filtered['observations']['colour'], 'zorder':self.read_instance.data_in_memory_filtered['observations']['zorder']} 
                    #create new experiment array colocated to observations
                    exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[data_label]['data'])
                    exp_data[nan_instances] = np.NaN
                    self.read_instance.data_in_memory_filtered['%s_colocatedto_observations'%(data_label)] = {'data':exp_data, 'colour':self.read_instance.data_in_memory_filtered[data_label]['colour'], 'zorder':self.read_instance.data_in_memory_filtered[data_label]['zorder']}
                    #update exps_all_nan array, making False all instances where have valid experiment data
                    exps_all_nan = np.all([exps_all_nan, nan_exp], axis=0) 
                    
            #create observational data array colocated to be non-NaN whenever there is a valid data in at least 1 experiment
            exps_all_nan = np.any([nan_obs,exps_all_nan],axis=0)
            obs_data = copy.deepcopy(self.read_instance.data_in_memory_filtered['observations']['data'])
            obs_data[exps_all_nan] = np.NaN
            self.read_instance.data_in_memory_filtered['observations_colocatedto_experiments'] = {'data':obs_data, 'colour':self.read_instance.data_in_memory_filtered['observations']['colour'], 'zorder':self.read_instance.data_in_memory_filtered['observations']['zorder']}   

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def handle_colocate_update(self):

        '''function that handles the update of the MPL canvas with colocated data upon checking of the colocate checkbox'''

        if self.read_instance.block_MPL_canvas_updates == False:

            #if only have 1 data array in memory (i.e. observations), no colocation is possible, therefore set colocation_active to be False, and return
            if len(list(self.read_instance.data_in_memory_filtered.keys())) == 1:
                self.colocate_active = False
                return

            #else, if have loaded experiment data, check if colocate checkbox is checked or unchecked  
            check_state = self.read_instance.ch_colocate.checkState()
            #update variable to inform plotting functions whether to use colocated data/or not 
            if check_state == QtCore.Qt.Checked: 
                self.colocate_active = True
            else:
                self.colocate_active = False

            #update z statistic/experiment bias comboboxes (without updating canvas) 
            self.read_instance.block_MPL_canvas_updates = True
            self.handle_map_z_statistic_update()
            self.handle_experiment_bias_update()
            self.read_instance.block_MPL_canvas_updates = False
            
            #update plotted map z statistic (if necessary)
            #calculate map z statistic (for selected z statistic) --> updating active map valid station indices
            self.calculate_z_statistic()

            #update plotted map z statistic
            self.update_map_z_statisitic()

            #update associated plots with selected stations
            self.update_associated_selected_station_plots()

            #draw changes
            self.draw()

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_map_z_statisitic(self):

        '''function that updates plotted z statistic on map, with colourbar'''

        #clear previously plotted station points
        try:
            self.map_points.remove()
        except:
            pass

        #clear/hide current colourbar axis
        self.cbar_ax.cla()
        self.cbar_ax.axis('off')

        #if have no valid active map indices, reset absolute/relative selected station indices to be empty lists
        #also uncheck select all/intersect checkboxes
        if len(self.active_map_valid_station_inds) == 0:
            #unselect all/intersect checkboxes
            self.read_instance.block_MPL_canvas_updates = True
            self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
            self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
            self.read_instance.block_MPL_canvas_updates = False
            #clear previously selected relative/absolute station indices
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)
            self.relative_selected_station_inds = np.array([],dtype=np.int)
            self.absolute_selected_station_inds = np.array([],dtype=np.int)

        #otherwise plot valid active map stations on map
        else:
            #if any of the currently selected stations are not in the current active map valid station indices --> unselect selected stations (and associated plots)
            #also uncheck select all/intersect checkboxes
            if np.all(np.in1d(self.relative_selected_station_inds, self.active_map_valid_station_inds)) == False:
                #unselect all/intersect checkboxes
                self.read_instance.block_MPL_canvas_updates = True
                self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
                self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
                self.read_instance.block_MPL_canvas_updates = False
                #reset relative/absolute selected station indices to be empty lists
                self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)
                self.relative_selected_station_inds = np.array([],dtype=np.int)
                self.absolute_selected_station_inds = np.array([],dtype=np.int)

            #plot new station points on map - coloured by currently active z statisitic, setting up plot picker
            self.map_points = self.map_ax.scatter(self.read_instance.station_longitudes[self.active_map_valid_station_inds],self.read_instance.station_latitudes[self.active_map_valid_station_inds], s=3.0, c=self.z_statistic, vmin=self.z_vmin, vmax=self.z_vmax, cmap=self.z_colourmap, picker = 1, zorder=2, transform=self.datacrs, linewidth=0.0, alpha=None)     
            #create 2D numpy array of plotted station coordinates
            self.map_points_coordinates = np.vstack((self.read_instance.station_longitudes[self.active_map_valid_station_inds],self.read_instance.station_latitudes[self.active_map_valid_station_inds])).T

            # create colour normalisation instance
            colour_norm = matplotlib.colors.Normalize(vmin=self.z_vmin,vmax=self.z_vmax)
            #normalise z statistic to colourmap range
            norm_z_statisitic =  colour_norm(self.z_statistic)
            #get list of rgba tuples to set map point colours
            self.rgba_tuples = matplotlib.cm.get_cmap(self.z_colourmap)(norm_z_statisitic)

            #colourbar
            #generate colourbar tick array 
            tick_array = np.linspace(self.z_vmin, self.z_vmax, 5, endpoint=True)
            #plot colourbar
            cb = self.figure.colorbar(self.map_points, orientation = 'horizontal', cax=self.cbar_ax, label='', ticks=tick_array, extend='max')
            self.cbar_ax.tick_params(labelsize=8.0)
            #plot colourbar label
            self.cbar_ax.yaxis.set_label_position("right")
            self.cbar_ax.set_ylabel(self.z_label, fontsize=8.0, rotation=0, ha='left', va='top')
            #turn colourbar axis on
            self.cbar_ax.axis('on')

            #call update of map drawing (this is a hack to force map plot object to be updated correctly --> only done when draw is called)
            self.draw()

        #update map selection appropriately for z statistic
        self.update_map_station_selection()

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_map_station_selection(self):

        '''function that updates the visual selection of stations on map'''

        #update map title
        if len(self.relative_selected_station_inds) == 1:
            self.map_ax.set_title('%s Selected'%(self.read_instance.station_references[self.relative_selected_station_inds][0]), fontsize=8.5, pad=3)
        else:
            self.map_ax.set_title('%s Stations Selected of %s Available'%(len(self.relative_selected_station_inds), len(self.active_map_valid_station_inds)), fontsize=8.5, pad=3)         

        #reset alphas of all plotted stations 
        self.rgba_tuples[:,-1] = 1.0
        marker_sizes = np.full(len(self.z_statistic), 3.0)
        self.map_points.set_facecolor(self.rgba_tuples)
        #reset marker sizes of all plotted stations
        self.map_points.set_sizes(marker_sizes)                

        #if have some selected stations, update map plot with station selection (reducing alpha of non-selected stations, and increasing marker size of selected stations)
        if len(self.relative_selected_station_inds) > 0:

            #decrease alpha of non-selected stations    
            absolute_non_selected_stations = np.nonzero(~np.in1d(range(self.z_statistic.shape[0]), self.absolute_selected_station_inds))[0]
            if len(absolute_non_selected_stations) > 0: 
                self.rgba_tuples[absolute_non_selected_stations,-1] = 0.25
                self.map_points.set_facecolor(self.rgba_tuples)

            #increase marker size of selected stations 
            marker_sizes[self.absolute_selected_station_inds] = 8.0    
            self.map_points.set_sizes(marker_sizes)
    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_associated_selected_station_plots(self):

        '''function that updates all plots associated with selected stations on map'''

        #clear/hide relevant axes
        #clear axes
        self.ts_ax.cla()
        self.violin_hours_ax.cla()
        self.violin_months_ax.cla()
        self.violin_days_ax.cla()
        self.exp_bias_hours_ax.cla()
        self.exp_bias_months_ax.cla()
        self.exp_bias_days_ax.cla()
        self.station_metadata_ax.cla()
            
        #hide axes
        self.ts_ax.axis('off')
        self.violin_hours_ax.axis('off')
        self.violin_months_ax.axis('off')
        self.violin_days_ax.axis('off')
        self.exp_bias_hours_ax.axis('off')
        self.exp_bias_months_ax.axis('off')
        self.exp_bias_days_ax.axis('off')
        self.station_metadata_ax.axis('off')

        #if have selected stations on map, then now remake plots
        if len(self.relative_selected_station_inds) > 0:

            #put selected data for each data array into pandas dataframes
            self.to_pandas_dataframe()

            #temporally aggregate selected data dataframes (by hour, day of week, month)
            self.pandas_temporal_aggregation()
        
            #if have some experiment data associated with selected stations, calculate temporally aggregated basic statistic differences and bias statistics between observations and experiment data arrays
            if len(list(self.selected_station_data.keys())) > 1:
                self.calculate_temporally_aggregated_experiment_bias_statistics()

            #update time series plot
            self.update_time_series_plot()

            #update violin plots for temporally aggregated data
            self.update_violin_plots()
        
            #if have some experiment data associated with selected stations, update experiment bias plots for temporally aggregated statistical differences/biases between observations and experiments
            if len(list(self.selected_station_data.keys())) > 1:
                self.update_experiment_bias_aggregated_plots()

            #update plotted station selected metadata
            self.update_selected_station_metadata()

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_experiment_grid_domain_edges(self):

        '''function that plots grid domain edges of experiments in memory'''

        #iterate through previously plotted experiment domain edge polygons, clearing them
        try:
            for grid_edge_polygon in self.grid_edge_polygons:
                grid_edge_polygon.remove()
        except:
            pass

        #create new array for storing plotted experiment domain edge polygons
        self.grid_edge_polygons = []

        #iterate through read experiments and plot grid domain edges on map
        for experiment in sorted(list(self.read_instance.data_in_memory.keys())):
            if experiment != 'observations':
                #compute map projection coordinates for each pair of longitude/latitude experiment grid edge coordinates
                #exp_x,exp_y = self.bm(self.read_instance.data_in_memory[experiment]['grid_edge_longitude'], self.read_instance.data_in_memory[experiment]['grid_edge_latitude'])
                #create matplotlib polygon object from experiment grid edge map projection coordinates 
                grid_edge_outline_poly = Polygon(np.vstack((self.read_instance.data_in_memory[experiment]['grid_edge_longitude'], self.read_instance.data_in_memory[experiment]['grid_edge_latitude'])).T, edgecolor=self.read_instance.data_in_memory[experiment]['colour'], linewidth=1, linestyle='--', fill=False, zorder=1, transform=self.datacrs)
                #plot grid edge polygon on map
                self.grid_edge_polygons.append(self.map_ax.add_patch(grid_edge_outline_poly))

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_legend(self):

        '''function that updates legend'''

        #create legend elements
        #add observations element
        legend_elements = [Line2D([0], [0], marker='o', color='white', markerfacecolor=self.read_instance.data_in_memory['observations']['colour'], markersize=11, label='observations')]
        #add element for each experiment
        for experiment_ind, experiment in enumerate(sorted(list(self.read_instance.data_in_memory.keys()))):
            if experiment != 'observations':
                #add experiment element
                legend_elements.append(Line2D([0], [0], marker='o', color='white', markerfacecolor=self.read_instance.data_in_memory[experiment]['colour'], markersize=11, label=experiment))
        
        #plot legend
        self.legend_ax.legend(handles=legend_elements, loc='best',mode='expand', ncol=4, fontsize=9.0)

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def to_pandas_dataframe(self):

        '''function that takes selected station data within data arrays and puts it into a pandas dataframe'''

        #create new dictionary to store selection station data by data array
        self.selected_station_data = {}

        #iterate through data arrays in data in memory filtered dictionary
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()):   
    
            #if colocation is not active, do not convert colocated data arrays to pandas data frames
            if self.colocate_active == False:
                if 'colocated' in data_label:
                    continue
            #else, if colocation is active, do not convert non-colocated data arrays to pandas data frames
            elif self.colocate_active == True:
                if 'colocated' not in data_label:
                    continue

            #observational arrays
            if data_label.split('_')[0] == 'observations':
                #get data for selected stations
                data_array = self.read_instance.data_in_memory_filtered[data_label]['data'][:,self.relative_selected_station_inds]
            #experiment arrays
            else:          
                #get intersect between selected station indices and valid available indices for experiment data array
                valid_selected_station_indices = np.intersect1d(self.relative_selected_station_inds, self.read_instance.data_in_memory_filtered[data_label]['valid_station_inds'])
                #get data for valid selected stations
                data_array = self.read_instance.data_in_memory_filtered[data_label]['data'][:,valid_selected_station_indices]

            #if data array has no valid data for selected stations, do not create a pandas dataframe
            #data array has valid data?
            if data_array.size:
                #add nested dictionary for data array name to selection station data dictionary
                self.selected_station_data[data_label] = {}
                #take cross station median of selected data for data array, and place it in a pandas dataframe -->  add to selected station data dictionary
                self.selected_station_data[data_label]['pandas_df'] = pd.DataFrame(np.nanmedian(data_array, axis=1), index=self.read_instance.time_array, columns=['data'])

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def pandas_temporal_aggregation(self):

        '''function that aggregates pandas dataframe data, for all data arrays, into desired temporal groupings
           also calculates all defined basic statistics for each individual temporal grouping
        '''

        #define statistics to calculate (all basic statistics)
        statistics_to_calculate = list(basic_stats_dict.keys())

        #iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:  

            #define all possible xticks for temporal resolution
            if temporal_aggregation_resolution == 'hour':  
                all_xticks = np.arange(24, dtype=np.int)
            elif temporal_aggregation_resolution == 'dayofweek':
                all_xticks = np.arange(7, dtype=np.int)
            elif temporal_aggregation_resolution == 'month':  
                all_xticks = np.arange(1, 13, dtype=np.int)

            #iterate through data arrays names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):  

                #create nested dictionary inside selected station data dictionary for storing aggregated data by data array label and temporal aggregation resolution
                self.selected_station_data[data_label][temporal_aggregation_resolution] = {}

                #aggregate data array into desired temporal groups (dropping NaNs)
                grouped_data = [g['data'].dropna() for n, g in self.selected_station_data[data_label]['pandas_df'].groupby(getattr(self.selected_station_data[data_label]['pandas_df'].index, temporal_aggregation_resolution))]
                #drop groups which have no data
                grouped_data = [group for group in grouped_data if len(group) > 0]              
                #get xticks for groups which have valid data (i.e. which hours/days/months have valid data)
                valid_xticks = [getattr(group.index, temporal_aggregation_resolution)[0] for group in grouped_data]
                #create array of the size of the full range of each aggregation periods, initialised with empty lists per element (i.e. 24 for hourly aggregation)
                full_grouped_data = [[] for _ in range(len(all_xticks))]
                #place valid grouped data in correct positions within full array
                full_group_indices_to_place = np.array([np.where(all_xticks==valid_xtick)[0][0] for valid_xtick in valid_xticks], dtype=np.int)
                for grouped_data_ii, full_group_index_to_place in enumerate(full_group_indices_to_place):
                    full_grouped_data[full_group_index_to_place] = grouped_data[grouped_data_ii]

                #add full grouped data to selected data dictionary
                self.selected_station_data[data_label][temporal_aggregation_resolution]['grouped_data'] = full_grouped_data
                #add valid xticks for group to selected data dictionary (i.e. the group xtick indexes which have valid data)
                self.selected_station_data[data_label][temporal_aggregation_resolution]['valid_xticks'] = valid_xticks              

                #calculate basic statistics in each group and add them to selected station data dictionary             
                for stat in statistics_to_calculate:
                    #get specific statistic dictionary (containing necessary information for calculation of selected statistic)
                    stats_dict = basic_stats_dict[stat]
                    #load default statistic arguments for passing to statistical function
                    function_arguments = stats_dict['arguments']
                    #create empty array for storing calculated statistic by group
                    stat_output_by_group = []
                    #iterate through grouped data
                    for group in full_grouped_data:
                        #only calculate statistic if have valid data in group
                        if len(group) > 0:
                            #add aggregated group data as argument to pass to statistical function 
                            function_arguments['data'] = group
                            #calculate statistic (appending to all group statistic output array) 
                            stat_output_by_group = np.append(stat_output_by_group, stats_dict['function'](**function_arguments))                      
                        #if no valid data in group, append NaN
                        else:
                            stat_output_by_group = np.append(stat_output_by_group, np.NaN)
                    #save statistical output by group to selected station data dictionary 
                    self.selected_station_data[data_label][temporal_aggregation_resolution][stat] = stat_output_by_group

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def calculate_temporally_aggregated_experiment_bias_statistics(self):

        '''function that calculates temporally aggregated basic statistic differences and bias statistics between observations and experiment data arrays'''

        #define all basic statistics that will be subtracted (each experiment - observations) for each temporal aggregation
        basic_statistics = list(basic_stats_dict.keys())
        #define all experiment bias statistics that will be calculated between each experiment and observations for each temporal aggregation
        bias_statistics = list(experiment_bias_stats_dict.keys())

        #iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            #iterate through data arrays names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):   

                #make sure the data array is an experimental one
                if data_label.split('_')[0] != 'observations':

                    #get relevant aggregated observational statistics dictionary (i.e. colocated or not)
                    if self.colocate_active == False:
                        relevant_aggregated_observations_dict = self.selected_station_data['observations'][temporal_aggregation_resolution]
                    else:
                        exp = data_label.split('_colocatedto_')[0]
                        relevant_aggregated_observations_dict = self.selected_station_data['observations_colocatedto_%s'%(exp)][temporal_aggregation_resolution]

                    #-----------------------------------------------------#                        
                    #calculate temporally aggregated basic statistic differences between experiment and observations

                    #iterate through basic statistics
                    for basic_stat in basic_statistics:

                        #create empty array for storing calculated experiment-observations difference statistic by group
                        stat_diff_by_group = []

                        #iterate through aggregated index times
                        for group_ii in range(len(relevant_aggregated_observations_dict[basic_stat])):

                            #get observational and experiment group aggregated statistic
                            group_obs_stat = relevant_aggregated_observations_dict[basic_stat][group_ii]
                            group_exp_stat = self.selected_station_data[data_label][temporal_aggregation_resolution][basic_stat][group_ii]

                            #take difference between observations and experiment statistics, if both values not NaN
                            if (np.isnan(group_obs_stat) == False) & (np.isnan(group_exp_stat) == False):
                                #calculate difference statistic (experiment - observations)    
                                stat_diff_by_group = np.append(stat_diff_by_group, group_exp_stat-group_obs_stat)
                            #else, if one (or both) of observations/experiment statistics are NaN, append NaN
                            else:
                                stat_diff_by_group = np.append(stat_diff_by_group, np.NaN)
                        #save statistical difference output by group to selected station data dictionary 
                        self.selected_station_data[data_label][temporal_aggregation_resolution]['%s_bias'%(basic_stat)] = stat_diff_by_group
                            
                    #-----------------------------------------------------#  
                    #if colocation is active, calculate temporally aggregated experiment bias statistical differences between experiment and observations
                    
                    if self.colocate_active == True:

                        #iterate through bias statistics
                        for bias_stat in bias_statistics:

                            #get specific statistic dictionary (containing necessary information for calculation of selected statistic)
                            stats_dict = experiment_bias_stats_dict[bias_stat]
                            #load default statistic arguments for passing to statistical function
                            function_arguments = stats_dict['arguments']
                            #create empty array for storing calculated statistic by group
                            stat_output_by_group = []
                            #iterate through experimental grouped data
                            for group_ii, exp_group in enumerate(self.selected_station_data[data_label][temporal_aggregation_resolution]['grouped_data']):
                                #add aggregated observational and experiment group data as arguments to pass to statistical function 
                                function_arguments['obs'] = relevant_aggregated_observations_dict['grouped_data'][group_ii]
                                function_arguments['exp'] = exp_group
                                
                                #calculate experiment bias statistic between observations and experiment group data (if have valid data in both groups)
                                if (len(function_arguments['obs']) > 0) & (len(function_arguments['exp']) > 0):
                                    #calculate experiment bias statistic
                                    stat_output_by_group = np.append(stat_output_by_group, stats_dict['function'](**function_arguments)) 
                                #if no valid data in one (or both) of observations/experiment groups, append NaN
                                else:
                                    stat_output_by_group = np.append(stat_output_by_group, np.NaN)
                            #save experiment bias statistic by group to selected station data dictionary 
                            self.selected_station_data[data_label][temporal_aggregation_resolution]['%s_bias'%(bias_stat)] = stat_output_by_group

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 
    
    def update_time_series_plot(self):

        '''function that updates time series plot upon selection of station/s'''

        #turn axis on
        self.ts_ax.axis('on')
    
        #iterate through data array names in selected station data dictionary
        for data_label in list(self.selected_station_data.keys()):   
    
            #if colocation is active, for observations use the 'observations_colocatedto_experiments', 
            if self.colocate_active == True:
                if data_label.split('_')[0] == 'observations':
                    if data_label != 'observations_colocatedto_experiments':
                        continue

            #plot time series data
            self.data_array_ts = self.ts_ax.plot(self.selected_station_data[data_label]['pandas_df'].dropna(), color=self.read_instance.data_in_memory_filtered[data_label]['colour'], marker = 'o', markeredgecolor = None, mew = 0, markersize = 1.1, linestyle = 'None', zorder=self.read_instance.data_in_memory_filtered[data_label]['zorder'])
              
        #set axes labels
        self.ts_ax.set_ylabel('Concentration (%s)'%(self.read_instance.measurement_units), fontsize=8.0)

        #plot grid
        self.ts_ax.grid(color='lightgrey',alpha=0.8)

        #set axis tick label sizes
        self.ts_ax.tick_params(labelsize=8.0)

        #as are re-plotting on time series axis, reset the navigation toolbar stack dictionaries entries associated with time series axis 
        self.reset_ax_navigation_toolbar_stack(self.ts_ax)

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_violin_plots(self):

        '''function that updates violin plots of temporally aggregated data upon selection of station/s'''

        #define dictionaries defining the relevant axis, axis titles, x axis ticks (and an empty nested dictionary for storing plot objects) for different temporal aggregation resolutions  
        hour_aggregation_dict =      {'ax':self.violin_hours_ax,  'title':'H',   'xticks':np.arange(24, dtype=np.int),    'plots':{}}
        dayofweek_aggregation_dict = {'ax':self.violin_days_ax,   'title':'DoW', 'xticks':np.arange(7, dtype=np.int),     'plots':{}}
        month_aggregation_dict =     {'ax':self.violin_months_ax, 'title':'M',   'xticks':np.arange(1, 13, dtype=np.int), 'plots':{}}

        #based on the temporal resolution of the data, combine the relevant temporal aggregation dictionaries  
        if self.read_instance.active_resolution == 'hourly':
            aggregation_dict = {'hour':hour_aggregation_dict, 'dayofweek':dayofweek_aggregation_dict, 'month':month_aggregation_dict}
        elif self.read_instance.active_resolution == 'daily':
            aggregation_dict = {'dayofweek':dayofweek_aggregation_dict, 'month':month_aggregation_dict}
        elif self.read_instance.active_resolution == 'monthly':
            aggregation_dict = {'month':month_aggregation_dict}

        #turn on all axes that will be plotted on, and add yaxis grid to each axis, and change axis label tick sizes
        for temporal_aggregation_resolution in list(aggregation_dict.keys()):
            #turn on axis
            aggregation_dict[temporal_aggregation_resolution]['ax'].axis('on')   
            #add yaxis grid
            aggregation_dict[temporal_aggregation_resolution]['ax'].yaxis.grid(color='lightgrey',alpha=0.8)
            #add axis aggregation resolution label
            aggregation_dict[temporal_aggregation_resolution]['ax'].annotate(aggregation_dict[temporal_aggregation_resolution]['title'], (0, 1), xytext=(2, -2), xycoords='axes fraction', textcoords='offset points', fontsize=9.0, ha='left', va='top')
            #change axis tick labels
            aggregation_dict[temporal_aggregation_resolution]['ax'].tick_params(labelsize=8.0)

        #------------------------------------------------------------------------------------------------# 
        #now, make violin plots for each temporally aggregated data array, for all relevant temporal aggregation resolutions       

        #iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            #create arrays for storing all calculated aggregated p5 and p95 (across all data arrays) for later limiting ylim 
            all_p5 = []
            all_p95 = []

            #iterate through data array names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):   

                #if colocation is active, for plotting observational aggregated data, use the 'observations_colocatedto_experiments' array 
                if self.colocate_active == True:
                    if data_label.split('_')[0] == 'observations':
                        if data_label != 'observations_colocatedto_experiments':
                            continue

                #get grouped data for current temporal aggregation resolution
                grouped_data = self.selected_station_data[data_label][temporal_aggregation_resolution]['grouped_data']
                #drop any groups which have no data
                grouped_data = [group for group in grouped_data if len(group) > 0]    

                #make violin plot --> add plotted object to aggregation dictionary
                aggregation_dict[temporal_aggregation_resolution]['plots'][data_label] = aggregation_dict[temporal_aggregation_resolution]['ax'].violinplot(grouped_data, positions=self.selected_station_data[data_label][temporal_aggregation_resolution]['valid_xticks'], points=100, widths=0.85, showmeans=False, showmedians=False, showextrema=False)
                #append aggregated p5/p95 for data array, to all_p5/all_p95 arrays
                all_p5 = np.append(all_p5, self.selected_station_data[data_label][temporal_aggregation_resolution]['p5']) 
                all_p95 = np.append(all_p95, self.selected_station_data[data_label][temporal_aggregation_resolution]['p95'])            

            #set x axis limits
            aggregation_dict[temporal_aggregation_resolution]['ax'].set_xlim(np.min(aggregation_dict[temporal_aggregation_resolution]['xticks'])-0.5, np.max(aggregation_dict[temporal_aggregation_resolution]['xticks'])+0.5)
            #set y axis limits (use minimum p5 and maximum p95 across aggregations, across all data arrays)
            aggregation_dict[temporal_aggregation_resolution]['ax'].set_ylim(np.nanmin(all_p5), np.nanmax(all_p95))
            #set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
            if temporal_aggregation_resolution == 'hour':
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(aggregation_dict[temporal_aggregation_resolution]['xticks'][::3])
            else:
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(aggregation_dict[temporal_aggregation_resolution]['xticks'])
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticklabels([temporal_axis_mapping_dict[temporal_aggregation_resolution][xtick] for xtick in aggregation_dict[temporal_aggregation_resolution]['xticks']])
                
        #------------------------------------------------------------------------------------------------#

        #format violin plots
            
        #iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            #iterate through data arrays have plotted data
            for data_label in list(aggregation_dict[temporal_aggregation_resolution]['plots'].keys()):

                #get plotted object per data array
                violin_plot = aggregation_dict[temporal_aggregation_resolution]['plots'][data_label]
                
                #update plotted objects with necessary colour, zorder and alpha
                for patch in violin_plot['bodies']:
                    patch.set_facecolor(self.read_instance.data_in_memory_filtered[data_label]['colour'])
                    patch.set_zorder(self.read_instance.data_in_memory_filtered[data_label]['zorder'])
                    if data_label.split('_')[0] == 'observations':                  
                        patch.set_alpha(0.7)
                    else:
                        patch.set_alpha(0.5)
                    #if have at least 1 experiment data array, split the violin plot across the horizontal (observations on left, experiments on right) 
                    if len(list(aggregation_dict[temporal_aggregation_resolution]['plots'].keys())) > 1:
                        m = np.mean(patch.get_paths()[0].vertices[:, 0])
                        #observations on left
                        if data_label.split('_')[0] == 'observations':
                            patch.get_paths()[0].vertices[:, 0] = np.clip(patch.get_paths()[0].vertices[:, 0], -np.inf, m)
                        #experiments on right
                        else:
                            patch.get_paths()[0].vertices[:, 0] = np.clip(patch.get_paths()[0].vertices[:, 0], m, np.inf)
                
                #-------------------------#
                #overplot time series of medians over boxes in necessary color

                #generate zorder to overplot medians in same order as violin plots are ordered, but on top of them
                median_zorder = (self.read_instance.data_in_memory_filtered['observations']['zorder']+len(list(aggregation_dict[temporal_aggregation_resolution]['plots'].keys())) - 1) + self.read_instance.data_in_memory_filtered[data_label]['zorder']

                #get xticks (all valid aggregated time indexes) and medians to plot
                xticks = aggregation_dict[temporal_aggregation_resolution]['xticks']
                medians = self.selected_station_data[data_label][temporal_aggregation_resolution]['p50']
                #split arrays if there are any temporal gaps to avoid line drawn being interpolated across missing values
                inds_to_split = np.where(np.diff(xticks) > 1)[0]
                if len(inds_to_split) == 0:
                    aggregation_dict[temporal_aggregation_resolution]['ax'].plot(xticks, medians, marker='o', color=self.read_instance.data_in_memory_filtered[data_label]['colour'], markersize=3, linewidth=0.5, zorder=median_zorder)
                else:
                    inds_to_split += 1
                    start_ind = 0
                    for end_ind in inds_to_split:
                        aggregation_dict[temporal_aggregation_resolution]['ax'].plot(xticks[start_ind:end_ind], medians[start_ind:end_ind], marker='o', color=self.read_instance.data_in_memory_filtered[data_label]['colour'], markersize=3, linewidth=0.5, zorder=median_zorder)
                        start_ind = end_ind   
                    aggregation_dict[temporal_aggregation_resolution]['ax'].plot(xticks[start_ind:], medians[start_ind:], marker='o', color=self.read_instance.data_in_memory_filtered[data_label]['colour'], markersize=3, linewidth=0.5, zorder=median_zorder) 

        #------------------------------------------------------------------------------------------------#
        #plot title (with units)
            
        #if selected data resolution is 'hourly', plot the title on off the hourly aggregation axis 
        if self.read_instance.active_resolution == 'hourly':
            self.violin_hours_ax.set_title('Temporal Distributions (%s)'%(self.read_instance.measurement_units), fontsize=8.0, loc='left') 
        #otherwise, plot the units on the monthly aggregation axis
        else:    
            self.violin_months_ax.set_title('Temporal Distributions (%s)'%(self.read_instance.measurement_units), fontsize=8.0, loc='left') 

        #------------------------------------------------------------------------------------------------#
        #as are re-plotting on violin plot axes, reset the navigation toolbar stack dictionaries entries associated with each of the axes 
        
        for temporal_aggregation_resolution in list(aggregation_dict.keys()):
            self.reset_ax_navigation_toolbar_stack(aggregation_dict[temporal_aggregation_resolution]['ax'])
        
    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_experiment_bias_aggregated_plots(self):

        '''function that updates the temporally aggregated experiment bias statistic plots'''
            
        #get currently selected experiment bias statistic 
        selected_stat = self.read_instance.cb_experiment_bias_stat.currentText()
        #get name of array to retrieve pre-calculated bias for selected statistic
        selected_experiment_bias_stat = '%s_bias'%(selected_stat)

        #define dictionaries defining the relevant axis, axis titles, x axis ticks (and an empty nested dictionary for storing plot objects) for different temporal aggregation resolutions  
        hour_aggregation_dict =      {'ax':self.exp_bias_hours_ax,  'title':'H',   'xticks':np.arange(24, dtype=np.int),    'plots':{}}
        dayofweek_aggregation_dict = {'ax':self.exp_bias_days_ax,   'title':'DoW', 'xticks':np.arange(7, dtype=np.int),     'plots':{}}
        month_aggregation_dict =     {'ax':self.exp_bias_months_ax, 'title':'M',   'xticks':np.arange(1, 13, dtype=np.int), 'plots':{}}

        #based on the temporal resolution of the data, combine the relevant temporal aggregation dictionaries  
        if self.read_instance.active_resolution == 'hourly':
            aggregation_dict = {'hour':hour_aggregation_dict, 'dayofweek':dayofweek_aggregation_dict, 'month':month_aggregation_dict}
        elif self.read_instance.active_resolution == 'daily':
            aggregation_dict = {'dayofweek':dayofweek_aggregation_dict, 'month':month_aggregation_dict}
        elif self.read_instance.active_resolution == 'monthly':
            aggregation_dict = {'month':month_aggregation_dict}

        #turn on all axes that will be plotted on, and add yaxis grid to each axis, and change axis label tick sizes
        for temporal_aggregation_resolution in list(aggregation_dict.keys()):
            #turn on axis
            aggregation_dict[temporal_aggregation_resolution]['ax'].axis('on')   
            #add yaxis grid
            aggregation_dict[temporal_aggregation_resolution]['ax'].yaxis.grid(color='lightgrey',alpha=0.8)
            #add axis aggregation resolution label
            aggregation_dict[temporal_aggregation_resolution]['ax'].annotate(aggregation_dict[temporal_aggregation_resolution]['title'], (0, 1), xytext=(2, -2), xycoords='axes fraction', textcoords='offset points', fontsize=9, ha='left', va='top')
            #change axis tick labels
            aggregation_dict[temporal_aggregation_resolution]['ax'].tick_params(labelsize=8.0)

        #------------------------------------------------------------------------------------------------#
        #plot title (with units)
        
        #create title string to plot (based on type of statistic plotting)
        if selected_stat not in self.read_instance.basic_z_stats:
            stats_dict = experiment_bias_stats_dict[selected_stat]
            plot_title = 'Experiment %s'%(stats_dict['label'])
        else:
            stats_dict = basic_stats_dict[selected_stat]
            if selected_stat != 'Data %':
                title_units = ' (%s)'%(self.read_instance.measurement_units)
            else:
                title_units = ''
            plot_title = 'Experiment %s bias%s'%(stats_dict['label'], title_units)

        #if selected data resolution is 'hourly', plot the title on off the hourly aggregation axis 
        if self.read_instance.active_resolution == 'hourly':
            self.exp_bias_hours_ax.set_title(plot_title, fontsize=8.0, loc='left') 
        #otherwise, plot the units on the monthly aggregation axis
        else:    
            self.exp_bias_months_ax.set_title(plot_title, fontsize=8.0, loc='left') 

        #get value/s of minimum bias for statistic
        minimum_bias = stats_dict['minimum_bias']

        #------------------------------------------------------------------------------------------------# 
        #now, make experiment bias plots for active bias statistic for all relevant temporal aggregation resolutions       

        #iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            #iterate through data array names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):   

                #if data array is observational, continue to next experiment data array
                if data_label.split('_')[0] == 'observations':
                    continue
                #else, make temporally aggregated plot for currently active experiment bias statistic 
                else:
                    aggregation_dict[temporal_aggregation_resolution]['plots'][data_label] = aggregation_dict[temporal_aggregation_resolution]['ax'].plot(aggregation_dict[temporal_aggregation_resolution]['xticks'], self.selected_station_data[data_label][temporal_aggregation_resolution][selected_experiment_bias_stat], color=self.read_instance.data_in_memory_filtered[data_label]['colour'], marker = 'o', zorder=self.read_instance.data_in_memory_filtered[data_label]['zorder'], markersize=3, linewidth=0.5)            

            #set x axis limits
            aggregation_dict[temporal_aggregation_resolution]['ax'].set_xlim(np.min(aggregation_dict[temporal_aggregation_resolution]['xticks'])-0.5, np.max(aggregation_dict[temporal_aggregation_resolution]['xticks'])+0.5)
            #set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
            if temporal_aggregation_resolution == 'hour':
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(aggregation_dict[temporal_aggregation_resolution]['xticks'][::3])
            else:
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(aggregation_dict[temporal_aggregation_resolution]['xticks'])
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticklabels([temporal_axis_mapping_dict[temporal_aggregation_resolution][xtick] for xtick in aggregation_dict[temporal_aggregation_resolution]['xticks']])
                
            #plot horizontal line/s across x axis at value/s of minimum experiment bias
            for mb in minimum_bias:
                aggregation_dict[temporal_aggregation_resolution]['ax'].axhline(y=mb, linestyle='--', linewidth=1.0, color='black', zorder=0)

        #------------------------------------------------------------------------------------------------#
        #as are re-plotting on experiment bias axes, reset the navigation toolbar stack dictionaries entries associated with each of the axes 
        
        for temporal_aggregation_resolution in list(aggregation_dict.keys()):
            self.reset_ax_navigation_toolbar_stack(aggregation_dict[temporal_aggregation_resolution]['ax'])

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 
    #function which creates a 

    #def update_experiment_bias_taylor_plot(self):

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_selected_station_metadata(self): 

        '''function which updates the plotted metadata detail of selected stations on the map'''

        #--------------------------------------------------# 
        #get some details of the station metadata axis --> to set limit for wrapping text
        #get axis bounding box 
        ax_bbox = self.station_metadata_ax.get_window_extent().transformed(self.figure.dpi_scale_trans.inverted())
        #get axis dimensions in inches
        ax_width_inches = ax_bbox.width
        #get axis dimensions in pixels
        ax_width_px = ax_width_inches * self.figure.dpi 

        #initialise string to plot on axis
        str_to_plot = ''

        #define dictionary with metadata variable names paired with metadata variable names to plot
        metadata_variable_naming = {'station_reference':'Station Reference',
                                    'station_name':'Station Name',
                                    'latitude':'Latitude', 
                                    'longitude':'Longitude',
                                    'measurement_altitude':'Measurement Altitude',
                                    'country':'Country',
                                    'network':'Network',
                                    'standardised_network_provided_area_classification':'Area Classification',
                                    'standardised_network_provided_station_classification':'Station Classification',
                                    'standardised_network_provided_main_emission_source':'Main Emission Source',
                                    'standardised_network_provided_land_use':'Land Use',
                                    'standardised_network_provided_terrain':'Terrain',
                                    'standardised_network_provided_measurement_scale':'Measurement Scale',
                                    'representative_radius':'Representative Radius',
                                    'GSFC_coastline_proximity':'To Coast',
                                    'primary_sampling_type':'Sampling Instrument Type',
                                    'sample_preparation_methodology_types':'Sample Preparation',
                                    'measurement_methodology':'Measurement Method',
                                    'measuring_instrument_name':'Measuring Instrument',   
                                    'measuring_instrument_sampling_type':'Measuring Instrument Sampling'
                                   }

        #is just 1 station selected?
        if len(self.relative_selected_station_inds) == 1: 

            #get station reference of selected station
            selected_station_reference = self.read_instance.station_references[self.relative_selected_station_inds][0]

            #add station reference, latitude, longitude and measurement altitude and GSFC coastline proximity to string (have 1 unique associated metadata value per station)
            str_to_plot += "%s   "%(selected_station_reference)
            str_to_plot += "Latitude: %s   "%(str(round(self.read_instance.station_latitudes[self.relative_selected_station_inds][0],4)))
            str_to_plot += "Longitude: %s\n"%(str(round(self.read_instance.station_longitudes[self.relative_selected_station_inds][0],4)))
            str_to_plot += "Measurement Altitude: %sm   "%(str(round(self.read_instance.station_measurement_altitudes[self.relative_selected_station_inds][0],2)))
            str_to_plot += "To Coast: %skm\n"%(str(round(self.read_instance.station_GSFC_coastline_proximities[self.relative_selected_station_inds][0],2)))
            
            #define other metadata variables to plot, in order to plot (plotting all unique associated metadata values)
            metadata_vars_to_plot = ['network', 'station_name', 'country', 'standardised_network_provided_area_classification', 
                                     'standardised_network_provided_station_classification', 'standardised_network_provided_terrain', 
                                     'standardised_network_provided_land_use', 'standardised_network_provided_main_emission_source', 
                                     'standardised_network_provided_measurement_scale', 'representative_radius',
                                     'measurement_methodology', 'measuring_instrument_name', 'measuring_instrument_sampling_type',
                                     'primary_sampling_type', 'sample_preparation_methodology_types']

            #iterate through metadata variables
            for meta_var in metadata_vars_to_plot:

                #get unique metadata values for selected station
                unique_station_meta = self.read_instance.station_metadata[selected_station_reference][meta_var]['unique']

                #is there just 1 unique value in the metadata array?   
                if len(unique_station_meta) == 1:
                    #set meta string as just the meta_var:unique value 
                    meta_string = '%s: %s\n'%(metadata_variable_naming[meta_var], unique_station_meta[0])

                #add meta string to str_to_plot
                str_to_plot += meta_string

        #more than 1 station selected?
        else:

            #get station references of all selected stations
            selected_station_references = self.read_instance.station_references[self.relative_selected_station_inds]

            #add N stations selected, in N countries
            str_to_plot += "%s Stations Selected\n"%(len(self.relative_selected_station_inds))
            #add median measurement altitude
            str_to_plot += "Median Measurement Altitude: %sm   "%(round(np.nanmedian(self.read_instance.station_measurement_altitudes[self.relative_selected_station_inds]),2))
            #add median GSFC coastline proximity
            str_to_plot += "Median To Coast: %skm\n"%(round(np.nanmedian(self.read_instance.station_GSFC_coastline_proximities[self.relative_selected_station_inds]),2))

            #get percentage of element occurrences across selected stations, for certain metadata variables
            metadata_vars_get_pc = ['network', 'country', 'standardised_network_provided_area_classification', 
                                    'standardised_network_provided_station_classification', 'standardised_network_provided_terrain', 
                                    'standardised_network_provided_land_use', 'standardised_network_provided_main_emission_source', 
                                    'standardised_network_provided_measurement_scale', 
                                    'measurement_methodology', 'measuring_instrument_name', 'measuring_instrument_sampling_type']
            #iterate through metadata variables
            for meta_var in metadata_vars_get_pc:

                #iterate through selected references, gathering all metadata across time        
                all_current_meta = []
                for selected_station_reference in selected_station_references:
                    all_current_meta = np.append(all_current_meta, self.read_instance.station_metadata[selected_station_reference][meta_var]['all'])

                #get counts of all unique metadata elements across selected stations 
                unique_meta, meta_counts = np.unique(all_current_meta, return_counts=True)
                #get number of unique metadata elements across selected stations
                n_unique_meta = len(unique_meta)
            
                #if have > 8 unique metadata elements, just return count of the elements across the selected stations
                if n_unique_meta > 8:
                    meta_string = '%s: %s unique elements\n'%(metadata_variable_naming[meta_var], n_unique_meta)
                #otherwise, get percentage of unique metadata elements across selected stations
                else:
                    meta_pc = (100./len(all_current_meta))*meta_counts
                    meta_pc = [str(round(meta,1))+'%' for meta in meta_pc]
                    #create string for variable to plot
                    meta_string = '%s: %s\n'%(metadata_variable_naming[meta_var], ', '.join([':'.join([str(var),pc]) for var, pc in zip(unique_meta, meta_pc)]))

                #add meta string to str_to_plot
                str_to_plot += meta_string

        #plot string to axis
        plot_txt = self.station_metadata_ax.text(0.0, 1.0, str_to_plot, ha='left', va='top', fontsize=8.0, transform=self.station_metadata_ax.transAxes, wrap=True, linespacing=1.5)
        #modify limit to wrap text as axis width in pixels  --> hack as matplotlib automatically sets limit as figure width
        plot_txt._get_wrap_line_width = lambda : ax_width_px

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def calculate_z_statistic(self):

        '''function that calculates selected z statistic for map'''
        
        #before doing anything check if have any valid station data for observations, if not update active map valid station indices to be empty list and return
        if len(self.read_instance.data_in_memory_filtered['observations']['valid_station_inds']) == 0:
            self.active_map_valid_station_inds = np.array([],dtype=np.int)
            return 

        #-------------------------------------------------#

        #get selected z1/z2 data arrays 
        z1_selected_name = self.read_instance.cb_z1.currentText()     
        z2_selected_name = self.read_instance.cb_z2.currentText()
        #check if have a selected z2 array (i.e. z2 name is not empty string)
        if z2_selected_name == '':
            have_z2 = False
        else:
            have_z2 = True
        
        #-------------------------------------------------#

        #get dictionary containing necessary information for calculation of selected statistic
        #check if the chosen statistic is a basic statistic
        z_statistic_name = self.read_instance.cb_z_stat.currentText()
        if z_statistic_name in list(basic_stats_dict.keys()):
            z_statistic_type = 'basic'
            stats_dict = basic_stats_dict[z_statistic_name]
            #set label units for statistic 
            if z_statistic_name != 'Data %':
                label_units = ' (%s)'%(self.read_instance.measurement_units)
            else:
                label_units = ''
        #if not a basic statistic, it must be an experiment bias statistic
        else:
            z_statistic_type = 'bias'
            stats_dict = experiment_bias_stats_dict[z_statistic_name]
            label_units = ''

        #-------------------------------------------------#

        #set colourbar for z statistic
        
        #first check if have defined colourbar for z statistic, if so use that
        if 'colourbar' in list(stats_dict.keys()):
            self.z_colourmap  = stats_dict['colourbar']
        #else, set appropriate colourmap for the type of statistic
        else:
            #if only have selected z1 array, the statistic is 'absolute', so use sequential colourbar
            if have_z2 == False:
                self.z_colourmap  = sequential_colourmap 
            #if have selected z1 and z2 arrays, the statistic is 'difference', so use diverging colourbar
            else:
                self.z_colourmap  = diverging_colourmap
    
        #-------------------------------------------------#

        #generate z colourbar label
        if have_z2 == False:
            self.z_label = '%s\n%s %s'%(z1_selected_name, stats_dict['label'], label_units)
        else:
            self.z_label = '%s - %s\n%s %s'%(z2_selected_name, z1_selected_name, stats_dict['label'], label_units)

        #-------------------------------------------------#
            
        #if colocation is active, set appropriate z1/z2 arrays to read to get colocated data arrays
        if (self.colocate_active == True):      
            #don't have z2 array?
            if have_z2 == False:            
                if z1_selected_name == 'observations':
                    z1_array_to_read = 'observations_colocatedto_experiments'
                else:
                    z1_array_to_read = '%s_colocatedto_observations'%(z1_selected_name)
            #have z2 array? 
            elif have_z2 == True:
                if z1_selected_name == 'observations':
                    z1_array_to_read = 'observations_colocatedto_%s'%(z2_selected_name)
                else:
                    z1_array_to_read = '%s_colocatedto_observations'%(z1_selected_name)       

                if z2_selected_name == 'observations':
                    z2_array_to_read = 'observations_colocatedto_%s'%(z1_selected_name)
                else:
                    z2_array_to_read = '%s_colocatedto_observations'%(z2_selected_name)   
        #else, simply use selected z1/z2 array names to read uncolocated data arrays  
        else:
            z1_array_to_read = copy.deepcopy(z1_selected_name)
            z2_array_to_read = copy.deepcopy(z2_selected_name)
            
        #-------------------------------------------------# 

        #read selected data arrays (after subsetting arrays by intersection of z1/z2 valid station indices) and calculate desired Z statistic (after removing NaNs from arrays)
        
        #get active map valid station indices (i.e. the indices of the stations data to plot on the map)
        #if only have z1, valid map indices are those simply for the z1 array
        if have_z2 == False:
            self.active_map_valid_station_inds = self.read_instance.data_in_memory_filtered[z1_array_to_read]['valid_station_inds']
        else:
            #if have z2 array, get intersection of z1 and z2 valid station indices
            self.active_map_valid_station_inds = np.intersect1d(self.read_instance.data_in_memory_filtered[z1_array_to_read]['valid_station_inds'], self.read_instance.data_in_memory_filtered[z2_array_to_read]['valid_station_inds'])

        #update absolute selected plotted station indices with respect to new active map valid station indices
        self.absolute_selected_station_inds = np.array([np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind in self.relative_selected_station_inds if selected_ind in self.active_map_valid_station_inds], dtype=np.int)

        #-------------------------------------------------#

        #read z1 data
        z1_array_data = self.read_instance.data_in_memory_filtered[z1_array_to_read]['data'][:,self.active_map_valid_station_inds]
        #drop NaNs and reshape to object list of station data arrays (if not checking data %)
        if z_statistic_name != 'Data %':
            z1_array_data = drop_NaNs(z1_array_data)
        else:
            z1_array_data = z1_array_data.transpose(1,0).tolist()
        
        #create empty array to store z statistic
        self.z_statistic = np.empty(len(z1_array_data))
        
        #if have no z2 data, calculate 'absolute' basic statistic
        if have_z2 == False:

            #load default selected z statistic arguments for passing to statistical function
            function_arguments = stats_dict['arguments']

            #iterate through stations calculating statistic
            for z_ii in range(len(self.z_statistic)):

                #set station z1 array as argument in argument dictionary
                function_arguments['data'] = z1_array_data[z_ii] 

                #calculate statistic
                self.z_statistic[z_ii] = stats_dict['function'](**function_arguments)
            
        #else, read z2 data then calculate 'difference' statistic
        else:
            #read z2 data 
            z2_array_data = self.read_instance.data_in_memory_filtered[z2_array_to_read]['data'][:,self.active_map_valid_station_inds]
            #drop NaNs and reshape to object list of station data arrays (if not checking data %)
            if z_statistic_name != 'Data %':
                z2_array_data = drop_NaNs(z2_array_data)
            else:
                z2_array_data = z2_array_data.transpose(1,0).tolist()
            
            #is the difference statistic basic (i.e. mean)?
            if z_statistic_type == 'basic':

                #load default selected z statistic arguments and make separate arguments dictionaries for z1/z2 calculations (as doing 2 separate calculations for z1/z2 and subtracting) 
                function_arguments_z1 = stats_dict['arguments']
                function_arguments_z2 = copy.deepcopy(function_arguments_z1)

                #iterate through stations calculating statistic
                for z_ii in range(len(self.z_statistic)):

                    #set station z1/z2 arrays as arguments in argument dictionaries
                    function_arguments_z1['data'] = z1_array_data[z_ii] 
                    function_arguments_z2['data'] = z2_array_data[z_ii]

                    #calculate statistics for z1/z2 arrays and subtract z2-z1
                    self.z_statistic[z_ii] = stats_dict['function'](**function_arguments_z2) - stats_dict['function'](**function_arguments_z1)

            #else, is the difference statistic an experiment bias statistic (i.e. r)?
            elif z_statistic_type == 'bias': 

                #load default selected z statistic arguments for passing to statistical function
                function_arguments = stats_dict['arguments']

                #iterate through stations calculating statistic
                for z_ii in range(len(self.z_statistic)):

                    #set station z1/z2 arrays as arguments in argument dictionary
                    function_arguments['obs'] = z1_array_data[z_ii] 
                    function_arguments['exp'] = z2_array_data[z_ii] 

                    #calculate statistic
                    self.z_statistic[z_ii] = stats_dict['function'](**function_arguments)

        #-------------------------------------------------#
        #if any station z statistics come out as NaN, remove respective stations from active map valid station indices
        #also cut z_statistic to remove invalid NaNs

        valid_z_statistic_boolean = ~np.isnan(self.z_statistic)
        self.active_map_valid_station_inds = self.active_map_valid_station_inds[valid_z_statistic_boolean]
        self.z_statistic = self.z_statistic[valid_z_statistic_boolean]

        #-------------------------------------------------#

        #calculate z vmin/vmax 
        #if have defined z vmin OR vmax for chosen statistic, use them by priority (if not taking a difference with a 'basic' statistic) 
        
        #check if have defined vmin
        have_defined_vmin = False
        if 'vmin' in list(stats_dict.keys()):
            #if z statistic type is 'basic' and taking a difference, then DO NOT use defined vmin 
            if (z_statistic_type == 'basic') & (have_z2 == False):
                have_defined_vmin = True
            elif z_statistic_type == 'bias':
                have_defined_vmin = True
        #have defined zmin?
        if have_defined_vmin  == True:
            self.z_vmin = stats_dict['vmin']
        #else, take vmin as minimum range value of calculated statistic
        else:
            self.z_vmin = np.nanmin(self.z_statistic)

        #check if have defined vmax
        have_defined_vmax = False
        if 'vmax' in list(stats_dict.keys()):
            #if z statistic type is 'basic' and taking a difference, then DO NOT use defined vmax
            if (z_statistic_type == 'basic') & (have_z2 == False):
                have_defined_vmax = True
            elif z_statistic_type == 'bias':
                have_defined_vmax = True
        #have defined zmax?
        if have_defined_vmax  == True:
            self.z_vmax = stats_dict['vmax']
        #else, take vmax as maximum range value of calculated statistic
        else:
            self.z_vmax = np.nanmax(self.z_statistic)
                
        #if z statistic is a 'difference', and do not have one of vmin/vmax pre-defined, force vmin/vmax to be symmetrical across 0
        if (have_z2 == True) & (have_defined_vmin  == False) & (have_defined_vmax  == False):
            limit_stat = np.max(np.abs([self.z_vmin,self.z_vmax]))
            self.z_vmin = -limit_stat
            self.z_vmax = limit_stat         

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def handle_map_z_statistic_update(self):

        '''define function which handles update of map z statistic'''
    
        if self.read_instance.block_config_bar_handling_updates == False:

            #-------------------------------------------#
            #update map z statistic comboboxes

            #set variable that blocks configuration bar handling updates until all changes to the z statistic comboboxes are made
            self.read_instance.block_config_bar_handling_updates = True

            #get currently selected items
            selected_z_stat = self.read_instance.cb_z_stat.currentText()
            selected_z1_array = self.read_instance.cb_z1.currentText()     
            selected_z2_array = self.read_instance.cb_z2.currentText()

            #if selected_z_stat and selected_z1_array are empty strings it is because they being initialised for the first time
            #force them to be 'observations' and first basic z statistic respectively
            if selected_z_stat == '':
                selected_z_stat = self.read_instance.basic_z_stats[0]
            if selected_z1_array == '':
                selected_z1_array  = 'observations'

            #update z statistic field to all basic stats if colocation not-active OR z2 array not selected, else select basic+bias stats
            if (self.colocate_active == False) or (selected_z2_array == ''):
                z_stat_items = copy.deepcopy(self.read_instance.basic_z_stats)
            else:
                z_stat_items = copy.deepcopy(self.read_instance.basic_and_bias_z_stats)

            #remove selected z1/z2 items from opposite z2/z1 comboboxes (if have value selected, i.e. z2 array not empty string)
            if selected_z2_array != '':
                z1_items = np.delete(self.read_instance.z1_arrays, np.where(self.read_instance.z1_arrays==selected_z2_array)[0])
            else:
                z1_items = self.read_instance.z1_arrays
            z2_items = np.delete(self.read_instance.z2_arrays, np.where(self.read_instance.z2_arrays==selected_z1_array)[0])
             
            #update all comboboxes (clear, then add items)
            self.read_instance.cb_z_stat.clear()
            self.read_instance.cb_z1.clear()
            self.read_instance.cb_z2.clear()
            self.read_instance.cb_z_stat.addItems(z_stat_items)
            self.read_instance.cb_z1.addItems(z1_items)
            self.read_instance.cb_z2.addItems(z2_items)
            #maintain currently selected z statistic (if exists in new item list)
            if selected_z_stat in z_stat_items:
                self.read_instance.cb_z_stat.setCurrentText(selected_z_stat)
            #maintain currently selected z1/z2 arrays
            self.read_instance.cb_z1.setCurrentText(selected_z1_array)
            self.read_instance.cb_z2.setCurrentText(selected_z2_array)

            #allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            #-------------------------------------------#
            if self.read_instance.block_MPL_canvas_updates == False:

                #calculate map z statistic (for selected z statistic) --> updating active map valid station indices
                self.calculate_z_statistic()

                #update plotted map z statistic
                self.update_map_z_statisitic()

                #draw changes
                self.draw()

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def handle_experiment_bias_update(self):

        '''define function that handles update of plotted experiment bias statisitics'''
        
        if self.read_instance.block_config_bar_handling_updates == False:

            #if no experiment data loaded, do not update
            if len(self.read_instance.experiment_bias_types) > 0:

                #-------------------------------------------#
                #update experiment bias comboboxes

                #set variable that blocks configuration bar handling updates until all changes to the experiment bias comboboxes are made
                self.read_instance.block_config_bar_handling_updates = True

                #get currently selected items
                selected_experiment_bias_type = self.read_instance.cb_experiment_bias_type.currentText()
                selected_experiment_bias_stat = self.read_instance.cb_experiment_bias_stat.currentText()

                #update experiment bias statistics (used for Aggregated field), to all basic stats if colocation not-active, and basic+bias stats if colocation active
                if self.colocate_active == False:
                    available_experiment_bias_stats = copy.deepcopy(self.read_instance.basic_z_stats)
                else:
                    available_experiment_bias_stats = copy.deepcopy(self.read_instance.basic_and_bias_z_stats)

                #if selected bias type is empty string, it is because fields are being initialised for the first time
                if selected_experiment_bias_type == '':
                    #set experiment bias type to be first available type  
                    selected_experiment_bias_type = self.read_instance.experiment_bias_types[0]
                    #set experiment bias stat to be first available stat    
                    selected_experiment_bias_stat = available_experiment_bias_stats[0]  

                #if selected bias type is 'Rank', then there are no stat options so force the available items and selected stat to be empty
                if selected_experiment_bias_type == 'Rank':            
                    available_experiment_bias_stats = []
                    selected_experiment_bias_stat = ''

                #update all comboboxes (clear, then add items)
                self.read_instance.cb_experiment_bias_type.clear()
                self.read_instance.cb_experiment_bias_stat.clear()
                self.read_instance.cb_experiment_bias_type.addItems(self.read_instance.experiment_bias_types)
                self.read_instance.cb_experiment_bias_stat.addItems(available_experiment_bias_stats)

                #update selected values
                self.read_instance.cb_experiment_bias_type.setCurrentText(selected_experiment_bias_type)
                #maintain currently selected bias statistic (if exists in new item list)
                if selected_experiment_bias_stat in available_experiment_bias_stats:
                    self.read_instance.cb_experiment_bias_stat.setCurrentText(selected_experiment_bias_stat)

                #allow handling updates to the configuration bar again
                self.read_instance.block_config_bar_handling_updates = False

                #-------------------------------------------#
                if self.read_instance.block_MPL_canvas_updates == False:

                    #update experiment bias plot/s if have some stations selected on map
                    if len(self.relative_selected_station_inds) > 0:

                        #clear and turn off all relevant axes before updating
                        self.exp_bias_hours_ax.cla()
                        self.exp_bias_months_ax.cla()
                        self.exp_bias_days_ax.cla()
                        self.exp_bias_hours_ax.axis('off')
                        self.exp_bias_months_ax.axis('off')
                        self.exp_bias_days_ax.axis('off')

                        #if experiment bias type == 'Aggregated' --> update plotted experiment bias plots
                        if selected_experiment_bias_type == 'Aggregated':
                            self.update_experiment_bias_aggregated_plots()

                        #draw changes
                        self.draw()

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 
    
    def select_all_stations(self):

        '''define function that selects/unselects all plotted stations (and associated plots) upon ticking of checkbox'''
        
        if self.read_instance.block_MPL_canvas_updates == False:

            #make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)
    
            #check if checkbox to select all stations is checked or unchecked   
            check_state = self.read_instance.ch_select_all.checkState()
        
            #if checkbox is checked, select all plotted stations
            if check_state == QtCore.Qt.Checked:
                self.relative_selected_station_inds = copy.deepcopy(self.active_map_valid_station_inds) 

                #if select intersect stations checkbox is checked then uncheck it (without updating canvas)
                if self.read_instance.ch_intersect.checkState() == QtCore.Qt.Checked:
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False  

            #if checkbox is unchecked then unselect all plotted stations
            elif check_state == QtCore.Qt.Unchecked:
                self.relative_selected_station_inds = np.array([], dtype=np.int) 

            #update absolute selected station indices (indices relative to plotted stations on map)
            self.absolute_selected_station_inds = np.arange(len(self.relative_selected_station_inds), dtype=np.int)

            #update map station selection
            self.update_map_station_selection() 

            #if selected stations have changed from previous selected, update associated plots
            if np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds) == False:
                self.update_associated_selected_station_plots()

            #draw changes
            self.draw()

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 
    
    def select_intersect_stations(self):

        '''define function that selects/unselects intersection of stations and all experiment domains (and associated plots) upon ticking of checkbox'''

        if self.read_instance.block_MPL_canvas_updates == False:

            #make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

            #check if checkbox to select intersection of stations is checked or unchecked   
            check_state = self.read_instance.ch_intersect.checkState()

            #if checkbox is unchecked then unselect all plotted stations
            if check_state == QtCore.Qt.Unchecked:
                self.relative_selected_station_inds = np.array([], dtype=np.int) 
                self.absolute_selected_station_inds = np.array([], dtype=np.int) 
        
            #else, if checkbox is checked then select all stations which intersect with all loaded experiment domains
            elif check_state == QtCore.Qt.Checked:

                #if have only observations loaded into memory, select all plotted stations
                if len(list(self.read_instance.data_in_memory.keys())) == 1:
                    self.relative_selected_station_inds = copy.deepcopy(self.active_map_valid_station_inds)
                    self.absolute_selected_station_inds = np.arange(len(self.relative_selected_station_inds), dtype=np.int)
                #else, define list of lists to get intersection between (active_map_valid_station_inds, and valid station indices associated with each loaded experiment array)
                else:
                    intersect_lists = [self.active_map_valid_station_inds]
                    for data_label in list(self.read_instance.data_in_memory.keys()):
                        if data_label != 'observations':
                            intersect_lists.append(self.read_instance.data_in_memory_filtered[data_label]['valid_station_inds'])
                    #get intersect between active map valid station indices and valid station indices associated with each loaded experiment array --> relative selected station indcies
                    self.relative_selected_station_inds = np.sort(list(set.intersection(*map(set,intersect_lists))))
                    #get absolute selected station indices (indices relative to plotted stations on map)
                    self.absolute_selected_station_inds = np.array([np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind in self.relative_selected_station_inds], dtype=np.int)

                #if select all stations checkbox is checked then uncheck it (without updating canvas)
                if self.read_instance.ch_select_all.checkState() == QtCore.Qt.Checked:
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False

            #update map station selection
            self.update_map_station_selection() 

            #if selected stations have changed from previous selected, update associated plots
            if np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds) == False:
                self.update_associated_selected_station_plots()

            #draw changes
            self.draw()
                
    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 
    #define functions that handle interactive station selection on map
    #the selection methods are individual station selection, or multiple selection with lasso

    def on_click(self, event):  

        '''function that handles single station selection upon mouse click'''
            
        #update variable to inform lasso handler that map as already been updated (to not redraw map)
        #the on_click function is only called when a station index has been selected  
        #the variable will be reset by lasso handler (which is always called after on_click)
        self.map_already_updated = True

        #make copy of current full array relative selected stations indices, before selecting new ones
        self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

        #get absolute selected index of station on map 
        self.absolute_selected_station_inds = np.array([event.ind[0]], dtype=np.int)

        #get selected station indices with respect to all available stations
        self.relative_selected_station_inds = self.map_selected_station_inds_to_all_available_inds(self.absolute_selected_station_inds)

        #update map station selection
        self.update_map_station_selection()

        #if selected stations have changed from previous selected, update associated plots
        if np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds) == False:
            self.update_associated_selected_station_plots()

        #draw changes
        self.draw()

    def onlassoselect(self, verts): 

        '''function that handles multiple station selection upon lasso drawing'''

        #unselect all/intersect checkboxes
        self.read_instance.block_MPL_canvas_updates = True
        self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.block_MPL_canvas_updates = False

        #check if have any plotted stations on map, if not, return
        if len(self.active_map_valid_station_inds) == 0:
            return

        #check if map as already been processed by on_click mouse click handling function, if so, return
        #the on_click function will always be called before lasso handler
        if self.map_already_updated == True: 
            self.map_already_updated = False
            return

        #make copy of current full array relative selected stations indices, before selecting new ones
        self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

        #get coordinates of drawn lasso
        lasso_path = Path(verts)
        lasso_path_vertices = lasso_path.vertices
        #transform lasso coordinates from projected to standard geographic coordinates
        lasso_path.vertices = self.datacrs.transform_points(self.plotcrs, lasso_path_vertices[:,0], lasso_path_vertices[:,1])[:,:2]
        #xys = self.map_points.get_offsets()
        #get absolute selected indices of stations on map (the station coordinates contained within lasso)
        #self.absolute_selected_station_inds = np.nonzero(lasso_path.contains_points(xys))[0]
        self.absolute_selected_station_inds = np.nonzero(lasso_path.contains_points(self.map_points_coordinates))[0]

        #get selected station indices with respect to all available stations
        self.relative_selected_station_inds = self.map_selected_station_inds_to_all_available_inds(self.absolute_selected_station_inds)

        #update map station selection
        self.update_map_station_selection()

        #hide lasso after selection
        self.lasso.set_visible(False)

        #if selected stations have changed from previous selected, update associated plots
        if np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds) == False:
            self.update_associated_selected_station_plots()

        #draw changes
        self.draw()

    def map_selected_station_inds_to_all_available_inds(self, selected_map_inds):

        '''function that takes the indices of selected stations on the map (potentially a subset of all available stations), and returns the indices of the stations inside the full loaded data arrays'''

        #index the array of indices of stations plotted on the map (indexed with respect to all available stations), with the absolute indices of the subset of plotted selected stations         
        return self.active_map_valid_station_inds[selected_map_inds]
        
#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

#define statistical calculation functions 

#---------------------------------------------------------#
#---------------------------------------------------------#

def drop_NaNs(data):

    '''function that returns numpy object of lists of station data with NaNs removed'''

    #reshape numpy array to have lists of data per station
    data = data.transpose(1,0).tolist()
    #iterate through each list of station data and remove NaNs
    for station_ii, station_data in enumerate(data):
        data[station_ii] = np.array(station_data)[~np.isnan(station_data)]
    #return numpy object of lists of station data with NaNs removed
    return np.array(data)

#---------------------------------------------------------#
#---------------------------------------------------------#
#Define basic statistic calculation functions:

#Mean
#Percentile
#Standard Deviation
#Variance
#Data Availability Fraction

def calculate_mean(data):
    '''calculate mean in a dataset'''
    return np.mean(data)

def calculate_percentile(data, percentile=50.0):
    '''calculate specific percentile in a dataset'''
    return np.percentile(data, percentile)

def calculate_standard_deviation(data):
    '''calculate standard deviation in a dataset'''
    return np.std(data)

def calculate_variance(data):
    '''calculate variance in a dataset'''
    return np.var(data)

def calculate_data_availability_fraction(data):
    '''calculate data availability fraction (i.e. fraction of total data array not equal to NaN)'''
    return (100./len(data)) * (np.count_nonzero(~np.isnan(data),axis=0))

def calculate_data_availability_number(data):
    '''calculate data availability absolute number (i.e. number of total data measurements not equal to NaN)'''
    return np.count_nonzero(~np.isnan(data),axis=0)

#define dictionary storing basic statistics that can be plotted
basic_stats_dict = {'Mean':  {'function':calculate_mean,                       'order':0,  'label':'Mean',                'arguments':{},                  'minimum_bias':[0.0]},
                    'StdDev':{'function':calculate_standard_deviation,         'order':1,  'label':'StdDev',              'arguments':{},                  'minimum_bias':[0.0]},
                    'Var':   {'function':calculate_variance,                   'order':2,  'label':'Variance',            'arguments':{},                  'minimum_bias':[0.0]},
                    'Data %':{'function':calculate_data_availability_fraction, 'order':3,  'label':'Data Availability %', 'arguments':{},                  'minimum_bias':[0.0],  'vmin':0.0, 'vmax':100.0},
                    'p1'  :  {'function':calculate_percentile,                 'order':4,  'label':'p1',                  'arguments':{'percentile':1.0},  'minimum_bias':[0.0]},
                    'p5'  :  {'function':calculate_percentile,                 'order':5,  'label':'p5',                  'arguments':{'percentile':5.0},  'minimum_bias':[0.0]},
                    'p10' :  {'function':calculate_percentile,                 'order':6,  'label':'p10',                 'arguments':{'percentile':10.0}, 'minimum_bias':[0.0]},
                    'p25' :  {'function':calculate_percentile,                 'order':7,  'label':'p25',                 'arguments':{'percentile':25.0}, 'minimum_bias':[0.0]},
                    'p50' :  {'function':calculate_percentile,                 'order':8,  'label':'p50',                 'arguments':{'percentile':50.0}, 'minimum_bias':[0.0]},
                    'p75' :  {'function':calculate_percentile,                 'order':9,  'label':'p75',                 'arguments':{'percentile':75.0}, 'minimum_bias':[0.0]},
                    'p90' :  {'function':calculate_percentile,                 'order':10, 'label':'p90',                 'arguments':{'percentile':90.0}, 'minimum_bias':[0.0]},
                    'p95' :  {'function':calculate_percentile,                 'order':11, 'label':'p95',                 'arguments':{'percentile':95.0}, 'minimum_bias':[0.0]},
                    'p99' :  {'function':calculate_percentile,                 'order':12, 'label':'p99',                 'arguments':{'percentile':99.0}, 'minimum_bias':[0.0]}}
    
#---------------------------------------------------------#
#---------------------------------------------------------#
#Define experiment bias evaluation statistics:

#Mean Absolute Error (MAE)
#Normalised Mean Absolute Error (NMAE)
#Mean Bias Error (MBE)
#Normalised Mean Bias Error (NMBE)
#Root Mean Squared Error (RMSE)
#Normalised Root Mean Squared Error (NRMSE)
#Absolute Percent Bias Error (APBE)
#Percent Bias Error (PBE)
#Coefficient of Efficiency (COE)
#Fraction of Predictions Within a Factor of Two (FAC2)
#Index of Agreement (IOA)
#Pearson correlation coefficient (r)
#Unpaired Peak Accuracy (UPA)

def calculate_APBE(obs, exp):
    '''Calculate absolute percent bias error (APBE) between observations and experiment'''
    return 100.0*np.sum(np.abs(exp-obs))/np.sum(obs)

def calculate_PBE(obs, exp):
    '''calculate percent bias error (PBE) between observations and experiment'''
    return 100.0*np.sum(exp-obs)/np.sum(obs)


def calculate_COE(obs, exp):
    '''Calculate coefficient of efficiency (COE) between observations and experiment, based on Legates and McCabe (1999, 2012).
       There have been many suggestions for measuring model performance over the years, but the COE is a simple formulation which is easy to interpret.
       A perfect model has a COE = 1. As noted by Legates and McCabe although the COE has no lower bound, a value of COE = 0.0 has a fundamental meaning. 
       It implies that the model is no more able to predict the observed values than does the observed mean. 
       Therefore, since the model can explain no more of the variation in the observed values than can the observed mean, such a model can have no predictive advantage.
       For negative values of COE, the model is less effective than the observed mean in predicting the variation in the observations.
       References:
       Legates DR, McCabe GJ. (1999). Evaluating the use of goodness-of-fit measures in hydrologic and hydroclimatic model validation. Water Resources Research 35(1): 233-241.
       Legates DR, McCabe GJ. (2012). A refined index of model performance: a rejoinder, International Journal of Climatology.
    å'''
    return 1.0 - (np.mean(np.abs(exp-obs)) / np.mean(np.abs(obs-np.mean(obs)))) 


def calculate_IOA(obs, exp):
    '''Calculate the Index of Agreement (IOA) between observations and experiment, based on Willmott et al. (2011)
       The metric spans between -1 and +1 with values approaching +1 representing better model performance.
       An IOA of 0.5, for example, indicates that the sum of the error-magnitudes is one half of the sum of the observed-deviation magnitudes. 
       When IOA = 0.0, it signifies that the sum of the magnitudes of the errors and the sum of the observed-deviation magnitudes are equivalent. 
       When IOA = -0.5, it indicates that the sum of the error-magnitudes is twice the sum of the perfect model-deviation and observed-deviation magnitudes. 
       Values of IOA near -1.0 can mean that the model-estimated deviations about O are poor estimates of the observed deviations; but, they also can mean that there simply is little observed variability - so some caution is needed when the IOA approaches -1.
       References;
       Willmott, C.J., Robeson, S.M., Matsuura, K., 2011. A refined index of model performance. International Journal of Climatology.
    ''' 
    return 1.0 -(np.sum((obs-exp)**2))/(np.sum((np.abs(exp-np.mean(obs))+np.abs(obs-np.mean(obs)))**2))

def calculate_MAE(obs, exp, normalisation_type='none'):
    '''Calculate mean absolute error (MAE)/ normalised mean absolute error (NMAE) between observations and experiment'''
    mae = np.mean(np.abs(exp-obs))
    #handle normalisation if desired
    if normalisation_type == 'max_min':
        mae = mae / (np.max(obs) - np.min(obs))  
    elif normalisation_type == 'mean':
        mae = mae / np.mean(obs) 
    elif normalisation_type == 'iq':
        mae = mae / (np.percentile(obs, 75) - np.percentile(obs, 25)) 
    elif normalisation_type == 'stdev':
        mae = mae / np.std(actual)
    return mae

def calculate_MBE(obs, exp, normalisation_type='none'):
    '''Calculate mean bias error (MBE)/ normalised mean bias error (NMBE) between observations and experiment'''
    mbe = np.mean(exp-obs)
    #handle normalisation if desired
    if normalisation_type == 'max_min':
        mbe = mbe / (np.max(obs) - np.min(obs))  
    elif normalisation_type == 'mean':
        mbe = mbe / np.mean(obs) 
    elif normalisation_type == 'iq':
        mbe = mae / (np.percentile(obs, 75) - np.percentile(obs, 25)) 
    elif normalisation_type == 'stdev':
        mbe = mbe / np.std(actual)
    return mbe

def calculate_RMSE(obs, exp, normalisation_type='none'):
    '''Calculate root mean squared error (RMSE) / normalised root mean squared error (NRMSE) between observations and experiment'''
    rmse = np.sqrt(np.mean((exp-obs)**2))
    #handle normalisation if desired
    if normalisation_type == 'max_min':
        rmse = rmse / (np.max(obs) - np.min(obs))  
    elif normalisation_type == 'mean':
        rmse = rmse / np.mean(obs) 
    elif normalisation_type == 'iq':
        rmse = rmse / (np.percentile(obs, 75) - np.percentile(obs, 25)) 
    elif normalisation_type == 'stdev':
        rmse = rmse / np.std(actual)
    return rmse


def calculate_r(obs, exp):
    '''Calculate the Pearson correlation coefficient (r) between observations and experiment
       The Pearson correlation coefficient measures the linear relationship between two datasets. 
       Strictly speaking, Pearson’s correlation requires that each dataset be normally distributed. 
       Like other correlation coefficients, this one varies between -1 and +1 with 0 implying no correlation. 
       Correlations of -1 or +1 imply an exact linear relationship. 
       Positive correlations imply that as x increases, so does y. 
       Negative correlations imply that as x increases, y decreases.
    '''
    return scipy.stats.pearsonr(obs, exp)[0]

def calculate_r_squared(obs, exp):
    '''Calculate the coefficient of determination, r squared, between observations and experiment
       It is the proportion of the variance in the dependent variable that is predictable from the independent variable(s).
       In linear least squares multiple regression with an estimated intercept term, the r squared equals the square of the Pearson correlation coefficient
    '''
    return calculate_r(obs, exp)**2

def calculate_FAC2(obs, exp):
    '''Calculate fraction of experiment values within a factor of two of observed values (FAC2)'''
    frac = exp/obs
    return (100.0/len(frac)) * len(frac[(frac >= 0.5) & (frac <= 2.0)])
    
def calculate_UPA(obs, exp):
    '''calculate unpaired peak accuracy (UPA)'''
    obs_max = np.max(obs)
    exp_max = np.max(exp)
    return (exp_max - obs_max) - obs_max

#define dictionary storing experiment bias evaluation statistics that can be plotted
experiment_bias_stats_dict = {'MAE':  {'function':calculate_MAE,       'order':0,  'label':'MAE',     'arguments':{},                            'minimum_bias':[0.0],   'vmin':0.0,                'colourbar':sequential_colourmap_warm},
                              'NMAE': {'function':calculate_MAE,       'order':1,  'label':'NMAE',    'arguments':{'normalisation_type':'mean'}, 'minimum_bias':[0.0],   'vmin':0.0,                'colourbar':sequential_colourmap_warm},
                              'MBE':  {'function':calculate_MBE,       'order':2,  'label':'MBE',     'arguments':{},                            'minimum_bias':[0.0]},
                              'NMBE': {'function':calculate_MBE,       'order':3,  'label':'NMBE',    'arguments':{'normalisation_type':'mean'}, 'minimum_bias':[0.0]},
                              'RMSE': {'function':calculate_RMSE,      'order':4,  'label':'RMSE',    'arguments':{},                            'minimum_bias':[0.0]},
                              'NRMSE':{'function':calculate_RMSE,      'order':5,  'label':'NRMSE',   'arguments':{'normalisation_type':'mean'}, 'minimum_bias':[0.0]},
                              'ABPE': {'function':calculate_APBE,      'order':6,  'label':'ABPE',    'arguments':{},                            'minimum_bias':[0.0],   'vmin':0.0,  'vmax':100.0, 'colourbar':sequential_colourmap_warm},
                              'PBE':  {'function':calculate_PBE,       'order':7,  'label':'PBE',     'arguments':{},                            'minimum_bias':[0.0]},
                              'COE':  {'function':calculate_COE,       'order':8,  'label':'COE',     'arguments':{},                            'minimum_bias':[1.0],                'vmax':1.0},
                              'FAC2': {'function':calculate_FAC2,      'order':9,  'label':'FAC2',    'arguments':{},                            'minimum_bias':[100.0], 'vmin':0.0,  'vmax':100.0, 'colourbar':sequential_colourmap_warm},
                              'IOA':  {'function':calculate_IOA,       'order':10, 'label':'IOA',     'arguments':{},                            'minimum_bias':[1.0],   'vmin':-1.0, 'vmax':1.0},
                              'r':    {'function':calculate_r,         'order':11, 'label':'r',       'arguments':{},                            'minimum_bias':[1.0],   'vmin':-1.0, 'vmax':1.0},
                              'r2':   {'function':calculate_r_squared, 'order':12, 'label':'r$^{2}$', 'arguments':{},                            'minimum_bias':[1.0],   'vmin':0.0,  'vmax':1.0,   'colourbar':sequential_colourmap_warm},
                              'UPA':  {'function':calculate_UPA,       'order':13, 'label':'UPA',     'arguments':{},                            'minimum_bias':[0.0]}}

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

#define dictionary for mapping days of week/months as integers to equivalent strings for writing on axes
temporal_axis_mapping_dict = {'dayofweek':{0:'M', 1:'T', 2:'W', 3:'T', 4:'F', 5:'S', 6:'S'},
                              'month':    {1:'J', 2:'F', 3:'M', 4:'A', 5:'M', 6:'J', 7:'J', 8:'A', 9:'S', 10:'O', 11:'N', 12:'D'}}

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

#return defined standardised data flag code associated with a specific standardised data flag name

#**Note: These standardised data flags refer specifically to information/or lack of information provided by the data provider 
#and not to post-processing/quality control of observations (given separately by the 'qa' field).

standard_data_flag_codes = {

#Basic Data Flags
#-----------------------------------------------------
'Valid Data': 0,

'Preliminary Data': 1,

'Missing Data': 2, 

'Invalid Data - Unspecified': 3,

'Un-Flagged Data': 4,
#-----------------------------------------------------


#Estimated Data Flags 
#-----------------------------------------------------
'Estimated Data - Unspecified': 10,

'Estimated Data - Negative Value Detected': 11,

'Estimated Data - No Value Detected': 12,

'Estimated Data - Value Below Detection Limit': 13,

'Estimated Data - Value Above Detection Limit': 14,

'Estimated Data - Value Substituted from Secondary Monitor': 15,

'Estimated Data - Multiple Parameters Aggregated': 16,
#-----------------------------------------------------


#Extreme/Irregular Data Flags
#-----------------------------------------------------
'Extreme/Irregular Data - Unspecified': 20,

'Data Does Not Meet Internal Network Quality Control Criteria': 21,

'High Variability of Data': 22,

'Irregular Data Manually Screened and Accepted': 23,

'Irregular Data Manually Screened and Rejected': 24,

'Negative Value': 25,

'No Value Detected': 26,
    
'Reconstructed/Recalculated Data': 27,

'Value Close to Detection Limit': 28,

'Value Below Acceptable Range': 29,

'Value Above Acceptable Range': 30,

'Value Below Detection Limit': 31,

'Value Above Detection Limit': 32,
#-----------------------------------------------------


#Measurement Issue Data Flags
#-----------------------------------------------------
'Measurement Issue - Unspecified': 40,

'Chemical Issue': 41,

'Erroneous Sampling Operation': 42,

'Extreme Internal Instrument Meteorological Conditions': 43,

'Extreme Ambient Laboratory Meteorological Conditions': 44,

'Extreme External Meteorological Conditions': 45,

'Extreme Sample Transport Conditions': 46,

'Invalid Flow Rate': 47,

'Human Error': 48,

'Low Data Capture': 49,

'Matrix Effect': 50,

'Mechanical Issue/Non-Operational Equipment': 51,

'No Technician': 52,

'Operational Maintenance Check Issue': 53,

'Physical Issue With Filter': 54,

'Power Failure': 55,

'Sample Diluted for Analysis': 56,

'Unmeasured Key Meteorological Parameter': 57,
#-----------------------------------------------------


#Operational Maintenance Data Flags
#-----------------------------------------------------
'Operational Maintenance - Unspecified': 70,

'Calibration': 71,

'Accuracy Check': 72,

'Blank Check': 73,

'Detection Limits Check': 74,

'Precision Check': 75,

'Retention Time Check': 76,

'Span Check': 77,

'Zero Check': 78,

'Instrumental Inspection': 79,

'Instrumental Repair': 80,

'Quality Control Audit': 81,
#-----------------------------------------------------


#Data Formatting/Processing Issue Data Flags
#-----------------------------------------------------
'Data Formatting/Processing Issue': 90,

'Corrected Data Formatting/Processing Issue': 91,
#-----------------------------------------------------


#Local Contamination Data Flags
#-----------------------------------------------------
'Local Contamination - Unspecified': 100,

'Agricultural Contamination': 101,

'Bird-Dropping Contamination': 102,

'Construction Contamination': 103,

'Dust Contamination': 104,

'Fire/Wood Burning Contamination': 105,

'Industrial Contamination': 106,

'Internal Laboratory/Instrument Contamination': 107,

'Insect Contamination': 108,

'Pollen/Leaf Contamination': 109,

'Sea-Salt Contamination': 110,

'Traffic Contamination': 111,
#-----------------------------------------------------


#Exceptional Event Data Flags
#-----------------------------------------------------
'Exceptional Event - Unspecified': 120,

#Natural Events
#---------------------------
'Dust Event': 121,

'Heavy Rain/Snowfall Shower (Squall)': 122,

'High Winds': 123,    

'Seismic Activity': 124,

'Station Inside Cloud': 125,

'Storm': 126,

'Stratospheric Ozone Intrusion': 127,

'Tropical Cyclone (Cyclone/Hurricane/Typhoon)': 128,

'Volcanic Eruptions': 129, 

'Wildfire': 130,

#Anthropogenically Induced Events
#---------------------------
'Chemical Spill/Industrial Accident': 131,

'Cleanup After a Major Disaster': 132,

'Demolition': 133,

'Fireworks': 134,
    
'Infrequent Large Gathering': 135,

'Terrorist Act': 136,
#-----------------------------------------------------


#Aggregation/Representation Flags
#-----------------------------------------------------
'Aggregation/Representation Issue - Unspecified': 150,

'Data Window Completeness < 90%': 151,

'Data Window Completeness < 75%': 152,

'Data Window Completeness < 66%': 153,

'Data Window Completeness < 50%': 154,

'Data Window Completeness < 25%': 155,

'>= 75% of Measurements in Window Below Detection Limit': 156,

'>= 50% of Measurements in Window Below Detection Limit': 157
#-----------------------------------------------------
}

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

#defined standardised qa flag code associated with specific data quality assurance checks
#can also be forced to return N flag codes rather than a specific flag code, by setting get_N_flags to True

standard_qa_flag_codes = {

#Basic QA Flags
#-----------------------------------------------------
#Missing Measurement
#Measurement is missing (i.e. NaN).
'Missing Measurement': 0,

#Infinite Value
#Value is infinite -- happens when data values are outside of the range that the float32 data type can handle (-3.4E+38 to +3.4E+38).
'Infinite Value': 1,

#Negative Measurement
#Measurement is negative in absolute terms.
'Negative Measurement': 2,

#Zero Measurement
#Have measurement equal to zero.
'Zero Measurement': 3,

#Invalid Data Provider Flags - GHOST Decreed
#Measurements are associated with data quality flags given by the data provider which have been decreed by the GHOST project architects to suggest the measurements are associated with substantial uncertainty/bias
'Invalid Data Provider Flags - GHOST Decreed': 4,

#Invalid Data Provider Flags - Network Decreed
#Measurements are associated with data quality flags given by the data provider which have been decreed by the reporting network to suggest the measurements are associated with substantial uncertainty/bias 
'Invalid Data Provider Flags - Network Decreed': 5,
#-----------------------------------------------------


#Duplicate / Overlapping Time Flags
#-----------------------------------------------------
#Duplicate Time - GHOST Decreed Valid Flagged Values Kept
#Multiple measurements reported for the same temporal window - all windows with GHOST decreed flagged invalid data are dropped
'Duplicate Time - GHOST Decreed Valid Flagged Values Kept': 10,

#Duplicate Time - First Value Kept
#Multiple measurements reported for the same temporal window - the first time window is kept preferentially
'Duplicate Time - First Value Kept': 11,

#Overlapping Time - GHOST Decreed Valid Flagged Values Kept
#Multiple measurements with overlapping temporal windows - all windows with GHOST decreed flagged invalid data are dropped
'Overlapping Time - GHOST Decreed Valid Flagged Values Kept': 12,

#Overlapping Time - Finest Temporal Resolutions Kept
#Measurements reported with overlapping temporal windows - only windows with a temporal resolution equal to the finest temporal resolution across the overlapping windows are kept 
'Overlapping Time - Finest Temporal Resolutions Kept': 13,

#Overlapping Time - First Value Kept
#Measurements reported with overlapping temporal windows - the first time window is kept preferentially.
'Overlapping Time - First Value Kept': 14,
#-----------------------------------------------------


#Significantly Shifting Measurement Position Flags
#-----------------------------------------------------
#Significant Latitude Shift
#Latitude shifts by >= 0.0001 degrees from previously reported latitude 
'Significant Latitude Shift': 20, 
    
#Significant Longitude Shift 
#Longitude shifts by >= 0.0001 degrees from previously reported latitude 
'Significant Longitude Shift': 21, 

#Significant Measurement Altitude Shift   
#Measurement Altitude shifts by >= 11m from previously reported measurement altitude
'Significant Measurement Altitude Shift': 22, 
#-----------------------------------------------------


#Duplicate station
#-----------------------------------------------------
#Station has been decreed to be a duplicate (i.e. reporting the same data as from another network, but the data from another network has been preferred)   
'Duplicate Station': 23,
#-----------------------------------------------------


#Metadata Assumption Flags
#-----------------------------------------------------
#Assumed Gas Volume
#Have assumed gas volume when converting between mass density and volume mixing ratio. (Do not have either temperature or pressure, or both).  
'Assumed Gas Volume': 30,

#No Latitude Metadata
#Latitude metadata field is absent, the most recent valid past latitude is assumed to be still valid. 
'No Latitude - Took Most Recent Valid Value': 31,

#Latitude metadata field is absent, the next valid latitude in the time record is assumed to be valid for this period (no available past valid latitudes).  
'No Latitude - Took Next Valid Value': 32,

#Latitude metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Latitude - Used Manually Compiled Metadata': 33,

#No Longitude Metadata
#Longitude metadata field is absent, the most recent valid past longitude is assumed to be still valid.  
'No Longitude - Took Most Recent Valid Value': 34,

#Longitude metadata field is absent, the next valid longitude in the time record is assumed to be valid for this period (no available past valid longitudes).  
'No Longitude - Took Next Valid Value': 35,

#Longitude metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Longitude - Used Manually Compiled Metadata': 36,

#No Altitude Metadata
#Altitude metadata field is absent, the most recent valid past altitude is assumed to be still valid.  
'No Altitude - Took Most Recent Valid Value': 37,

#Altitude Metadata field is absent, the next valid altitude in the time record is assumed to be valid for this period (no available past valid altitudes).  
'No Altitude - Took Next Valid Value': 38,

#Altitude metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Altitude - Used Manually Compiled Metadata': 39,

#Altitude metadata field is absent through entire time record, used ETOPO1 globally gridded altitudes to fill altitude for entire time record.
'No Altitude - Used ETOPO1 to Estimate Altitude': 40,

#No Sampling Height Metadata
#Sampling height metadata field is absent, the most recent valid past sampling height is assumed to be still valid.  
'No Sampling Height - Took Most Recent Valid Value': 41,

#Sampling Height metadata field is absent, the next valid sampling height in the time record is assumed to be valid for this period (no available past valid sampling heights).  
'No Sampling Height - Took Next Valid Value': 42,

#Sampling Height metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Sampling Height - Used Manually Compiled Metadata': 43,

#No Standardised Network Provided Area Classification Metadata
#Standardised Network Provided Area Classification metadata field is absent, the most recent valid past Standardised Network Provided Area Classification is assumed to be still valid.  
'No Area Classification - Took Most Recent Valid Value': 44,

#Standardised Network Provided Area Classification metadata field is absent, the next valid Standardised Network Provided Area Classification in the time record is assumed to be valid for this period (no available past Standardised Network Provided Area Classifications).
'No Area Classification - Took Next Valid Value': 45,

#Standardised Network Provided Area Classification metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Area Classification - Used Manually Compiled Metadata': 46,

#No Standardised Network Provided Station Classification Metadata
#Standardised Network Provided Station Classification metadata field is absent, the most recent valid past Standardised Network Provided Station Classification is assumed to be still valid.  
'No Station Classification - Took Most Recent Valid Value': 47,

#Standardised Network Provided Station Classification metadata field is absent, the next valid Standardised Network Provided Station Classification in the time record is assumed to be valid for this period (no available past Standardised Network Provided Station Classifications).
'No Station Classification - Took Next Valid Value': 48,

#Standardised Network Provided Station Classification metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Station Classification - Used Manually Compiled Metadata': 49,

#No Standardised Network Provided Main Emission Source Metadata
#Standardised Network Provided Main Emission Source metadata field is absent, the most recent valid past Standardised Network Provided Main Emission Source is assumed to be still valid.  
'No Main Emission Source - Took Most Recent Valid Value': 50,

#Standardised Network Provided Main Emission Source metadata field is absent, the next valid Standardised Network Provided Main Emission Source in the time record is assumed to be valid for this period (no available past Standardised Network Provided Main Emission Sources).
'No Main Emission Source - Took Next Valid Value': 51,

#Standardised Network Provided Main Emission Source metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Main Emission Source - Used Manually Compiled Metadata': 52,

#No Standardised Network Provided Land Use Metadata
#Standardised Network Provided Land Use metadata field is absent, the most recent valid past Standardised Network Provided Land Use is assumed to be still valid.  
'No Land Use - Took Most Recent Valid Value': 53,

#Standardised Network Provided Land Use metadata field is absent, the next valid Standardised Network Provided Land Use in the time record is assumed to be valid for this period (no available past Standardised Network Provided Land Uses).
'No Land Use - Took Next Valid Value': 54,

#Standardised Network Provided Land Use metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Land Use - Used Manually Compiled Metadata': 55,

#No Standardised Network Provided Terrain Metadata
#Standardised Network Provided Terrain metadata field is absent, the most recent valid past Standardised Network Provided Terrain is assumed to be still valid.  
'No Terrain - Took Most Recent Valid Value': 56,

#Standardised Network Provided Terrain metadata field is absent, the next valid Standardised Network Provided Terrain in the time record is assumed to be valid for this period (no available past Standardised Network Provided Terrains).
'No Terrain - Took Next Valid Value': 57,

#Standardised Network Provided Terrain metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Terrain - Used Manually Compiled Metadata': 58,

#No Standardised Network Provided Measurement Scale Metadata
#Standardised Network Provided Measurement Scale metadata field is absent, the most recent valid past Standardised Network Provided Measurement Scale is assumed to be still valid.  
'No Measurement Scale - Took Most Recent Valid Value': 59,

#Standardised Network Provided Measurement Scale metadata field is absent, the next valid Standardised Network Provided Measurement Scale in the time record is assumed to be valid for this period (no available past Standardised Network Provided Measurement Scales).
'No Measurement Scale - Took Next Valid Value': 60,

#Standardised Network Provided Measurement Scale metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Measurement Scale - Used Manually Compiled Metadata': 61,

#No Representative Radius Metadata
#Representative Radius metadata field is absent, the most recent valid past Representative Radius is assumed to be still valid.  
'No Representative Radius - Took Most Recent Valid Value': 62,

#Representative Radius metadata field is absent, the next valid Representative Radius in the time record is assumed to be valid for this period (no available past Representative Radii).
'No Representative Radius - Took Next Valid Value': 63,

#Representative Radius metadata field is absent through entire time record, used manually compiled metadata to fill field through entire time record.
'No Representative Radius - Used Manually Compiled Metadata': 64,
#-----------------------------------------------------

    
#Recurring Value
#-----------------------------------------------------
#Do check for persistently recurring values. check is done by using a moving window of 9 measurements. If 7/9 of values in the window are valid.
'Recurring Value': 70,
#-----------------------------------------------------


#Non-Integer Local Timezone (relative to UTC)
#-----------------------------------------------------
#Local timezone has been determined to be non-integer for time of measurement, relative to UTC.
'Non-Integer Local Timezone (relative to UTC)': 71,
#-----------------------------------------------------


#Extreme Data Flags
#-----------------------------------------------------
#Extreme Data - Scientifically Non-Feasible
#Data is greater than a scientifically feasible limit (variable by parameter).
'Extreme Data - Scientifically Non-Feasible': 80,

#Extreme Data - Distributional Outlier
#Data is screened through adjusted boxplot to determine distributional outliers
'Extreme Data - Distributional Outlier': 81,

#Extreme Data - Manually Decreed
#Data has been found and decreed manually to be extreme, a select section of data is flagged
'Extreme Data - Manually Decreed': 82,
#-----------------------------------------------------

    
#Insufficient Measurement Resolution Flags
#-----------------------------------------------------
#Insufficient Measurement Resolution - Documented
#The documented resolution of measurement is coarser than a set limit (variable by measured parameter).
'Insufficient Measurement Resolution - Documented': 83,

#Insufficient Measurement Resolution - Reported
#The reported resolution of measurement is coarser than a set limit (variable by measured parameter).
'Insufficient Measurement Resolution - Reported': 84,

#Insufficient Measurement Resolution - Preferential
#The preferential resolution of measurement (reported, and then documented) is coarser than a set limit (variable by measured parameter).
'Insufficient Measurement Resolution - Preferential': 85,

#No Documented/Reported Measurement Resolution Metadata
#No measurement resolution metadata is available, either reported or documented. 
'No Documented/Reported Measurement Resolution Metadata': 86,

#Insufficient Measurement Resolution - Empirical
#The resolution of measurement is analysed month by month. If the minimum difference between observations is coarser than a set limit (variable by measured parameter), measurements are flagged.
'Insufficient Measurement Resolution - Empirical': 87,
#-----------------------------------------------------


#Limit of Detection Flags
#-----------------------------------------------------
#Below Documented Lower Limit of Detection
#Measurement is below or equal to the instrumental documented lower limit of detection.
'Below Documented Lower Limit of Detection': 88,

#Below Reported Lower Limit of Detection
#Measurement is below or equal to the network reported lower limit of detection. 
'Below Reported Lower Limit of Detection': 89,

#Below Preferential Lower Limit of Detection
#Measurement is below or equal to the preferential lower limit of detection (reported, and then documented). 
'Below Preferential Lower Limit of Detection': 90,

#No Documented/Reported Lower Limit of Detection Metadata
#No lower limit of detection metadata is available, either reported or documented. 
'No Documented/Reported Lower Limit of Detection Metadata': 91,

#Above Documented Upper Limit of Detection
#Measurement is above or equal to the instrumental documented upper limit of detection. 
'Above Documented Upper Limit of Detection': 92,

#Above Reported Upper Limit of Detection
#Measurement is above or equal to the network reported upper limit of detection. 
'Above Reported Upper Limit of Detection': 93,

#Above Preferential Upper Limit of Detection
#Measurement is above or equal to the preferential upper limit of detection (reported, and then documented). 
'Above Preferential Upper Limit of Detection': 94,

#No Documented/Reported Upper Limit of Detection Metadata
#No upper limit of detection metadata is available, either reported or documented. 
'No Documented/Reported Upper Limit of Detection Metadata': 95,
#-----------------------------------------------------


#Measurement Methodology Flags
#-----------------------------------------------------
#Methodology Not Mapped
#The measurement methodology used has not yet been mapped to standardised dictionaries of measurement methodologies.
'Methodology Not Mapped': 100,

#Assumed Primary Sampling
#A level of assumption has been made in determining the primary sampling type.
'Assumed Primary Sampling': 101,

#Assumed Sample Preparation
#A level of assumption has been made in determining the sample preparation.
'Assumed Sample Preparation': 102,

#Assumed Measurement Methodology
#A level of assumption has been made in determining the measurement methodology. 
'Assumed Measurement Methodology': 103,

#Unknown Primary Sampling Instrument
#The specific name of the primary sampling instrument is unknown.
'Unknown Primary Sampling Instrument': 104,

#Unknown Measuring Instrument
#The specific name of measuring instrument is unknown. 
'Unknown Measuring Instrument': 105,

#Erroneous Primary Sampling
#The primary sampling is not appropriate to prepare the specific parameter for subsequent measurement.
'Erroneous Primary Sampling': 106,

#Erroneous Sample Preparation
#The sample preparation is not appropriate to prepare the specific parameter for subsequent measurement.
'Erroneous Sample Preparation': 107,

#Erroneous Measurement Method
#The measurement methodology used is not known to be able to measure the specific parameter. Only do check when known (or have assumed method).
'Erroneous Measurement Methodology': 108,

#Invalid QA Measurement Method 
#The specific measurement methodology has been decreed not to conform to QA standards as the method is not sufficiently proven/ subject to substantial biases/uncertainty. Only do check when known (or have assumed method).
'Invalid QA Measurement Methodology': 109,
#-----------------------------------------------------



#Hourly Temporal Representation Flags
#-----------------------------------------------------
#Two types are checks are done to flag the representativity of hourly data periods (starting and ending at NN:00 UTC).
#First is Data Completeness, i.e. what percentage across the hourly period is represented by valid data (after screening by key QA flags)?
#Second is Maximum Data Gap, i.e. what is the maximum data gap percentage in the provided valid data across the hourly period (after screening by key QA flags)? 
#Multiple separate limits are evaluated by, for each of the checks, to give flexibility in the definitions of temporal representativity.

#Hourly Window Data Completeness < 90%'
'Hourly Window Data Completeness < 90%': 150,

#Hourly Window Data Completeness < 75%'
'Hourly Window Data Completeness < 75%': 151,

#Hourly Window Data Completeness < 66%'
'Hourly Window Data Completeness < 66%': 152,

#Hourly Window Data Completeness < 50%'
'Hourly Window Data Completeness < 50%': 153,

#Hourly Window Data Completeness < 25%'
'Hourly Window Data Completeness < 25%': 154,

#Hourly Window Maximum Data Gap >= 5%'
'Hourly Window Maximum Data Gap >= 5%': 155,

#Hourly Window Maximum Data Gap >= 10%'
'Hourly Window Maximum Data Gap >= 10%': 156,

#Hourly Window Maximum Data Gap >= 15%'
'Hourly Window Maximum Data Gap >= 15%': 157,

#Hourly Window Maximum Data Gap >= 20%'
'Hourly Window Maximum Data Gap >= 20%': 158,

#Hourly Window Maximum Data Gap >= 25%'
'Hourly Window Maximum Data Gap >= 25%': 159,

#-----------------------------------------------------


#Daily Temporal Representation Flags
#-----------------------------------------------------
#Two types are checks are done to flag the representativity of daily data periods (starting and ending at NN 00:00 UTC).
#First is the Data Completeness, i.e. what percentage across the daily period is represented by valid data (after screening by key QA flags)?
#Second is the Maximum Data Gap, i.e. what is the maximum data gap percentage in the provided valid data across the daily period (after screening by key QA flags)? 
#Multiple separate limits are evaluated by, for each of the checks, to give flexibility in the definitions of temporal representativity.

#Daily Window Data Completeness < 90%'
'Daily Window Data Completeness < 90%': 160,

#Daily Window Data Completeness < 75%'
'Daily Window Data Completeness < 75%': 161,

#Daily Window Data Completeness < 66%'
'Daily Window Data Completeness < 66%': 162,

#Daily Window Data Completeness < 50%'
'Daily Window Data Completeness < 50%': 163,

#Daily Window Data Completeness < 25%'
'Daily Window Data Completeness < 25%': 164,

#Daily Window Maximum Data Gap >= 5%'
'Daily Window Maximum Data Gap >= 5%': 165,

#Daily Window Maximum Data Gap >= 10%'
'Daily Window Maximum Data Gap >= 10%': 166,

#Daily Window Maximum Data Gap >= 15%'
'Daily Window Maximum Data Gap >= 15%': 167,

#Daily Window Maximum Data Gap >= 20%'
'Daily Window Maximum Data Gap >= 20%': 168,

#Daily Window Maximum Data Gap >= 25%'
'Daily Window Maximum Data Gap >= 25%': 169,

#-----------------------------------------------------


#Weekly Temporal Representation Flags
#-----------------------------------------------------
#Two types are checks are done to flag the representativity of weekly data periods (starts and ends defined by UTC isocalendar).
#First is Data Completeness, i.e. what percentage across the weekly period is represented by valid data (after screening by key QA flags)?
#Second is Maximum Data Gap, i.e. what is the maximum data gap percentage in the provided valid data across the weekly period (after screening by key QA flags)? 
#Multiple separate limits are evaluated by, for each of the checks, to give flexibility in the definitions of temporal representativity.

#Weekly Window Data Completeness < 90%'
'Weekly Window Data Completeness < 90%': 170,

#Weekly Window Data Completeness < 75%'
'Weekly Window Data Completeness < 75%': 171,

#Weekly Window Data Completeness < 66%'
'Weekly Window Data Completeness < 66%': 172,

#Weekly Window Data Completeness < 50%'
'Weekly Window Data Completeness < 50%': 173,

#Weekly Window Data Completeness < 25%'
'Weekly Window Data Completeness < 25%': 174,

#Weekly Window Maximum Data Gap >= 5%'
'Weekly Window Maximum Data Gap >= 5%': 175,

#Weekly Window Maximum Data Gap >= 10%'
'Weekly Window Maximum Data Gap >= 10%': 176,

#Weekly Window Maximum Data Gap >= 15%'
'Weekly Window Maximum Data Gap >= 15%': 177,

#Weekly Window Maximum Data Gap >= 20%'
'Weekly Window Maximum Data Gap >= 20%': 178,

#Weekly Window Maximum Data Gap >= 25%'
'Weekly Window Maximum Data Gap >= 25%': 179,

#-----------------------------------------------------


#Monthly Temporal Representation Flags
#-----------------------------------------------------
#Two types are checks are done to flag the representativity of monthly data periods (starting and ending at NN:01 00:00 UTC).
#First is Data Completeness, i.e. what percentage across the monthly period is represented by valid data (after screening by key QA flags)?
#Second is Maximum Data Gap, i.e. what is the maximum data gap percentage in the provided valid data across the monthly period (after screening by key QA flags)? 
#Multiple separate limits are evaluated by, for each of the checks, to give flexibility in the definitions of temporal representativity.

#Monthly Window Data Completeness < 90%'
'Monthly Window Data Completeness < 90%': 180,

#Monthly Window Data Completeness < 75%'
'Monthly Window Data Completeness < 75%': 181,

#Monthly Window Data Completeness < 66%'
'Monthly Window Data Completeness < 66%': 182,

#Monthly Window Data Completeness < 50%'
'Monthly Window Data Completeness < 50%': 183,

#Monthly Window Data Completeness < 25%'
'Monthly Window Data Completeness < 25%': 184,

#Monthly Window Maximum Data Gap >= 5%'
'Monthly Window Maximum Data Gap >= 5%': 185,

#Monthly Window Maximum Data Gap >= 10%'
'Monthly Window Maximum Data Gap >= 10%': 186,

#Monthly Window Maximum Data Gap >= 15%'
'Monthly Window Maximum Data Gap >= 15%': 187,

#Monthly Window Maximum Data Gap >= 20%'
'Monthly Window Maximum Data Gap >= 20%': 188,

#Monthly Window Maximum Data Gap >= 25%'
'Monthly Window Maximum Data Gap >= 25%': 189,

#-----------------------------------------------------


#Annual Temporal Representation Flags
#-----------------------------------------------------
#Two types are checks are done to flag the representativity of annual data periods (starting and ending at NNNN:01:01 00:00 UTC).
#First is Data Completeness, i.e. what percentage across the annual period is represented by valid data (after screening by key QA flags)?
#Second is Maximum Data Gap, i.e. what is the maximum data gap percentage in the provided valid data across the annual period (after screening by key QA flags)? 
#Multiple separate limits are evaluated by, for each of the checks, to give flexibility in the definitions of temporal representativity.

#Annual Window Data Completeness < 90%'
'Annual Window Data Completeness < 90%': 190,

#Annual Window Data Completeness < 75%'
'Annual Window Data Completeness < 75%': 191,

#Annual Window Data Completeness < 66%'
'Annual Window Data Completeness < 66%': 192,

#Annual Window Data Completeness < 50%'
'Annual Window Data Completeness < 50%': 193,

#Annual Window Data Completeness < 25%'
'Annual Window Data Completeness < 25%': 194,

#Annual Window Maximum Data Gap >= 5%'
'Annual Window Maximum Data Gap >= 5%': 195,

#Annual Window Maximum Data Gap >= 10%'
'Annual Window Maximum Data Gap >= 10%': 196,

#Annual Window Maximum Data Gap >= 15%'
'Annual Window Maximum Data Gap >= 15%': 197,

#Annual Window Maximum Data Gap >= 20%'
'Annual Window Maximum Data Gap >= 20%': 198,

#Annual Window Maximum Data Gap >= 25%'
'Annual Window Maximum Data Gap >= 25%': 199,

#-----------------------------------------------------


#No Valid Data 
#-----------------------------------------------------
#After screening by key QA flags, no valid data remains.
'No Valid Data': 210

#-----------------------------------------------------    
}

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

#return defined standardised classification flag code associated with specific classification algorithms
#can also be forced to return N flag codes rather than a specific flag code, by setting get_N_flags to True
 
standard_classification_flag_codes = {

#Station Classifications
#-----------------------------------------------------
#High Altitude - Metadata Altitude
#Station determined to be measuring at an altitude >= 1500 metres relative to mean sea level, from altitudes + sampling heights taken from network provided metadata. 
'High Altitude - Metadata Altitude': 0,

#High Altitude - Metadata Derived
#Station determined to be a 'mountain' station, and therefore 'high altitude', derived from standardised 'terrain' network provided metadata.
'High Altitude - Metadata Derived': 1,

#High Altitude - ETOPO1 
#Station determined to be measuring at an altitude >= 1500 metres relative to sea level datum, from altitudes taken from ETOPO1 digital elevation model + sampling heights taken from network provided metadata. 
'High Altitude - ETOPO1': 2,

#High Altitude - Iwahashi Global Landform Classification 
#Station determined to be measuring at a high altitude, derived from the European Soil Data Centre Iwahashi Global Landform Classification. 
'High Altitude - Iwahashi Global Landform Classification': 3,

#High Altitude - Meybeck Global Landform Classification 
#Station determined to be measuring at a high altitude, derived from the European Soil Data Centre Meybeck Global Landform Classification. 
'High Altitude - Meybeck Global Landform Classification': 4,

#Near Coast - Metadata Derived
#Station determined to be located near the coast - derived from standardised 'terrain' network provided metadata. 
'Near Coast - Metadata Derived': 5,

#Near Coast - GSFC
#Station determined to be located < 50km of the coast (either over land or sea) - using GSFC nearest to coastline dataset (0.01 degree grid).
'Near Coast - GSFC': 6,

#Rural Station - Lenient Metadata Derived
#Station determined to be 'rural', using standardised network provided metadata, following lenient classifications. 
'Rural Station - Lenient Metadata Derived': 7,

#Urban Station - Lenient Metadata Derived
#Station determined to be 'urban', using standardised network provided metadata, following lenient classifications. 
'Urban Station - Lenient Metadata Derived': 8,

#Unclassified Station - Lenient Metadata Derived
#Station determined to be 'unclassified', using standardised network provided metadata, following lenient classifications.
'Unclassified Station - Lenient Metadata Derived': 9,

#Rural Station - Strict Metadata Derived
#Station determined to be 'rural', using standardised network provided metadata, following strict classifications. 
'Rural Station - Strict Metadata Derived': 10,

#Urban Station - Strict Metadata Derived
#Station determined to be 'urban', using standardised network provided metadata, following strict classifications. 
'Urban Station - Strict Metadata Derived': 11,

#Unclassified Station - Strict Metadata Derived
#Station determined to be 'unclassified', using standardised network provided metadata, following strict classifications.
'Unclassified Station - Strict Metadata Derived': 12,

#Rural Station - Anthrome (Native Resolution)
#Rural station as defined by using the UMBC Anthrome gridded classification dataset at native resolution (0.0833 degree grid).
'Rural Station - Anthrome (Native Resolution)': 13,

#Urban Station - Anthrome (Native Resolution)
#Urban station as defined by using the modal classification from the UMBC Anthrome gridded classification dataset at native resolution (0.0833 degree grid).
'Urban Station - Anthrome (Native Resolution)': 14,

#Rural Station - Anthrome (Mode in 5km Perimeter)
#Rural station as defined by using the modal classification from the UMBC Anthrome gridded classification dataset in a 5km perimeter around the station location. 
'Rural Station - Anthrome (Mode in 5km Perimeter)': 15,

#Urban Station - Anthrome (Mode in 5km Perimeter)
#Urban station as defined by using the modal classification from the UMBC Anthrome gridded classification dataset in a 5km perimeter around the station location. 
'Urban Station - Anthrome (Mode in 5km Perimeter)': 16,

#Rural Station - Anthrome (Mode in 25km Perimeter)
#Rural station as defined by using the modal classification from the UMBC Anthrome gridded classification dataset in a 25km perimeter around the station location. 
'Rural Station - Anthrome (Mode in 25km Perimeter)': 17,

#Urban Station - Anthrome (Mode in 25km Perimeter)
#Urban station as defined by using the modal classification from the UMBC Anthrome gridded classification dataset in a 25km perimeter around the station location. 
'Urban Station - Anthrome (Mode in 25km Perimeter)': 18,

#Rural Station - TOAR
#Rural station as defined by using a TOAR approach to classification (Tropospheric Ozone Assessment Report).
'Rural Station - TOAR': 19,

#Urban Station - TOAR
#Urban station as defined by using a TOAR approach to classification (Tropospheric Ozone Assessment Report).
'Urban Station - TOAR': 20,

#Unclassified Station - TOAR
#Unclassified station as defined using a TOAR approach to classification (Tropospheric Ozone Assessment Report).
'Unclassified Station - TOAR': 21,

#Rural Station - Joly-Peuch 
#Rural station as defined using a Joly-Peuch approach to classification
'Rural Station - Joly-Peuch': 22,

#Unclassified Station - Joly-Peuch 
#Unclassified station as defined using a Joly-Peuch approach to classification
'Unclassified Station - Joly-Peuch': 23,
#-----------------------------------------------------


#Temporal Period Flags - for classifying periods of time measurements are made within
#-----------------------------------------------------
#Daytime
#Time of measurement is daytime. Done by calculating the solar elevation angle for a latitude/longitude/measurement height at a certain timestamp.
'Daytime': 50,

#Nightime
#Time of measurement is nighttime. Done by calculating the solar elevation angle for a latitude/longitude/measurement height at a certain timestamp.
'Nighttime': 51,

#Weekday
#Time of measurement is weekday (by local time).
'Weekday': 52,

#Weekend
#Time of measurement is weekend (by local time).
'Weekend': 53,

#winter
#Time of measurement is northern hemisphere winter (measurement UTC time in months of December, January or February).
'Winter': 54,    

#spring
#Time of measurement is northern hemisphere spring (measurement UTC time in months of March, April or May).
'Spring': 55,    

#summer
#Time of measurement is northern hemisphere spring (measurement UTC time in months of June, July or August).
'Summer': 56,    

#autumn
#Time of measurement is northern hemisphere spring (measurement UTC time in months of September, October or November).
'Autumn': 57   
#-----------------------------------------------------
}

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

#define dictionary of all parameters with associated key information 

parameter_dictionary = {
'sconco3':         {'long_parameter_name':'ozone',                           'matrix':'GAS',     'standard_units':'nmol mol$^{-1}$', 'chemical_formula':'O3',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcno':         {'long_parameter_name':'nitrogen monoxide',               'matrix':'GAS',     'standard_units':'nmol mol$^{-1}$', 'chemical_formula':'NO',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcno2':        {'long_parameter_name':'nitrogen dioxide',                'matrix':'GAS',     'standard_units':'nmol mol$^{-1}$', 'chemical_formula':'NO2',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcco':         {'long_parameter_name':'carbon monoxide',                 'matrix':'GAS',     'standard_units':'nmol mol$^{-1}$', 'chemical_formula':'CO',      'extreme_lower_limit':0.0, 'extreme_upper_limit':10000.0},
'sconcisop':       {'long_parameter_name':'isoprene',                        'matrix':'GAS',     'standard_units':'nmol mol$^{-1}$', 'chemical_formula':'C5H8',    'extreme_lower_limit':0.0, 'extreme_upper_limit':500.0  },
'sconcso2':        {'long_parameter_name':'sulphur dioxide',                 'matrix':'GAS',     'standard_units':'nmol mol$^{-1}$', 'chemical_formula':'SO2',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcnh3':        {'long_parameter_name':'ammonia',                         'matrix':'GAS',     'standard_units':'nmol mol$^{-1}$', 'chemical_formula':'NH3',     'extreme_lower_limit':0.0, 'extreme_upper_limit':500.0  },
'sconchno3':       {'long_parameter_name':'nitric acid',                     'matrix':'GAS',     'standard_units':'nmol mol$^{-1}$', 'chemical_formula':'HNO3',    'extreme_lower_limit':0.0, 'extreme_upper_limit':500.0  },
'sconcpan':        {'long_parameter_name':'peroxyacetyl nitrate',            'matrix':'GAS',     'standard_units':'nmol mol$^{-1}$', 'chemical_formula':'C2H3NO5', 'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10':            {'long_parameter_name':'PM10 mass',                       'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'',        'extreme_lower_limit':0.0, 'extreme_upper_limit':2000.0 },
'pm2p5':           {'long_parameter_name':'PM2.5 mass',                      'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'',        'extreme_lower_limit':0.0, 'extreme_upper_limit':2000.0 },
'pm1':             {'long_parameter_name':'PM1 mass',                        'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'',        'extreme_lower_limit':0.0, 'extreme_upper_limit':2000.0 },
'sconcal':         {'long_parameter_name':'aluminium',                       'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Al',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcas':         {'long_parameter_name':'arsenic',                         'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'As',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcbc':         {'long_parameter_name':'black carbon',                    'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcc':          {'long_parameter_name':'carbon',                          'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconccorrected':  {'long_parameter_name':'carbon: corrected',               'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcca':         {'long_parameter_name':'calcium',                         'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Ca',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconccd':         {'long_parameter_name':'cadmium',                         'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cd',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconccl':         {'long_parameter_name':'chloride',                        'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cl',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconccobalt':     {'long_parameter_name':'cobalt',                          'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Co',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconccr':         {'long_parameter_name':'chromium',                        'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cr',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconccu':         {'long_parameter_name':'copper',                          'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cu',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcec':         {'long_parameter_name':'elemental carbon',                'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcfe':         {'long_parameter_name':'iron',                            'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Fe',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconchg':         {'long_parameter_name':'mercury',                         'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Hg',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconck':          {'long_parameter_name':'potassium',                       'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'K',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcmg':         {'long_parameter_name':'magnesium',                       'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Mg',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcmn':         {'long_parameter_name':'manganese',                       'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Mn',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 }, 
'sconcna':         {'long_parameter_name':'sodium',                          'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Na',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcnh4':        {'long_parameter_name':'ammonium',                        'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'NH4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcni':         {'long_parameter_name':'nickel',                          'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Ni',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcno3':        {'long_parameter_name':'nitrate',                         'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'NO3',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcoc':         {'long_parameter_name':'organic carbon',                  'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcoccorrected':{'long_parameter_name':'organic carbon: corrected',       'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcpb':         {'long_parameter_name':'lead',                            'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Pb',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcse':         {'long_parameter_name':'selenium',                        'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Se',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcso4':        {'long_parameter_name':'sulphate',                        'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'SO4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcso4nss':     {'long_parameter_name':'sulphate: non-sea salt',          'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'SO4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconcv':          {'long_parameter_name':'vanadium',                        'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'V',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'sconczn':         {'long_parameter_name':'zinc',                            'matrix':'AEROSOL', 'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Zn',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10al':          {'long_parameter_name':'PM10 aluminium',                  'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Al',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10as':          {'long_parameter_name':'PM10 arsenic',                    'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'As',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10bc':          {'long_parameter_name':'PM10 black carbon',               'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10c':           {'long_parameter_name':'PM10 carbon',                     'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10ccorrected':  {'long_parameter_name':'PM10 carbon: corrected',          'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10ca':          {'long_parameter_name':'PM10 calcium',                    'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Ca',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10cd':          {'long_parameter_name':'PM10 cadmium',                    'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cd',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10cl':          {'long_parameter_name':'PM10 chloride',                   'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cl',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10cobalt':      {'long_parameter_name':'PM10 cobalt',                     'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Co',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10cr':          {'long_parameter_name':'PM10 chromium',                   'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cr',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10cu':          {'long_parameter_name':'PM10 copper',                     'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cu',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10ec':          {'long_parameter_name':'PM10 elemental carbon',           'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10fe':          {'long_parameter_name':'PM10 iron',                       'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Fe',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10hg':          {'long_parameter_name':'PM10 mercury',                    'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Hg',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10k':           {'long_parameter_name':'PM10 potassium',                  'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'K',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10mg':          {'long_parameter_name':'PM10 magnesium',                  'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Mg',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10mn':          {'long_parameter_name':'PM10 manganese',                  'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Mn',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10na':          {'long_parameter_name':'PM10 sodium',                     'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Na',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10nh4':         {'long_parameter_name':'PM10 ammonium',                   'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'NH4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10ni':          {'long_parameter_name':'PM10 nickel',                     'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Ni',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10no3':         {'long_parameter_name':'PM10 nitrate',                    'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'NO3',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10oc':          {'long_parameter_name':'PM10 organic carbon',             'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10occorrected': {'long_parameter_name':'PM10 organic carbon: corrected',  'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10pb':          {'long_parameter_name':'PM10 lead',                       'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Pb',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10se':          {'long_parameter_name':'PM10 selenium',                   'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Se',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10so4':         {'long_parameter_name':'PM10 sulphate',                   'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'SO4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10so4nss':      {'long_parameter_name':'PM10 sulphate : non-sea salt',    'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'SO4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10v':           {'long_parameter_name':'PM10 vanadium',                   'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'V',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm10zn':          {'long_parameter_name':'PM10 zinc',                       'matrix':'PM10',    'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Zn',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5al':         {'long_parameter_name':'PM2.5 aluminium',                 'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Al',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5as':         {'long_parameter_name':'PM2.5 arsenic',                   'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'As',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5bc':         {'long_parameter_name':'PM2.5 black carbon',              'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5c':          {'long_parameter_name':'PM2.5 carbon',                    'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5ccorrected': {'long_parameter_name':'PM2.5 carbon: corrected',         'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5ca':         {'long_parameter_name':'PM2.5 calcium',                   'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Ca',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5cd':         {'long_parameter_name':'PM2.5 cadmium',                   'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cd',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },  
'pm2p5cl':         {'long_parameter_name':'PM2.5 chloride',                  'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cl',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5cobalt':     {'long_parameter_name':'PM2.5 cobalt',                    'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Co',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5cr':         {'long_parameter_name':'PM2.5 chromium',                  'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cr',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5cu':         {'long_parameter_name':'PM2.5 copper',                    'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cu',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5ec':         {'long_parameter_name':'PM2.5 elemental carbon',          'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5fe':         {'long_parameter_name':'PM2.5 iron',                      'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Fe',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5hg':         {'long_parameter_name':'PM2.5 mercury',                   'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Hg',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5k':          {'long_parameter_name':'PM2.5 potassium',                 'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'K',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5mg':         {'long_parameter_name':'PM2.5 magnesium',                 'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Mg',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5mn':         {'long_parameter_name':'PM2.5 manganese',                 'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Mn',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5na':         {'long_parameter_name':'PM2.5 sodium',                    'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Na',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5nh3':        {'long_parameter_name':'PM2.5 ammonia',                   'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'NH3',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5nh4':        {'long_parameter_name':'PM2.5 ammonium',                  'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'NH4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5ni':         {'long_parameter_name':'PM2.5 nickel',                    'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Ni',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5no3':        {'long_parameter_name':'PM2.5 nitrate',                   'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'NO3',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5oc':         {'long_parameter_name':'PM2.5 organic carbon',            'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5occorrected':{'long_parameter_name':'PM2.5 organic carbon: corrected', 'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5pb':         {'long_parameter_name':'PM2.5 lead',                      'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Pb',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5se':         {'long_parameter_name':'PM2.5 selenium',                  'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Se',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5so4':        {'long_parameter_name':'PM2.5 sulphate',                  'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'SO4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5so4nss':     {'long_parameter_name':'PM2.5 sulphate : non-sea salt',   'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'SO4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5v':          {'long_parameter_name':'PM2.5 vanadium',                  'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'V',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm2p5zn':         {'long_parameter_name':'PM2.5 zinc',                      'matrix':'PM2.5',   'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Zn',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1al':           {'long_parameter_name':'PM1 aluminium',                   'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Al',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1as':           {'long_parameter_name':'PM1 arsenic',                     'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'As',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1bc':           {'long_parameter_name':'PM1 black carbon',                'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1c':            {'long_parameter_name':'PM1 carbon',                      'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1ccorrected':   {'long_parameter_name':'PM1 carbon: corrected',           'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1ca':           {'long_parameter_name':'PM1 calcium',                     'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Ca',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1cd':           {'long_parameter_name':'PM1 cadmium',                     'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cd',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1cl':           {'long_parameter_name':'PM1 chloride',                    'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cl',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1cobalt':       {'long_parameter_name':'PM1 cobalt',                      'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Co',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1cr':           {'long_parameter_name':'PM1 chromium',                    'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cr',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1cu':           {'long_parameter_name':'PM1 copper',                      'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Cu',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1ec':           {'long_parameter_name':'PM1 elemental carbon',            'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1fe':           {'long_parameter_name':'PM1 iron',                        'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Fe',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1hg':           {'long_parameter_name':'PM1 mercury',                     'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Hg',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1k':            {'long_parameter_name':'PM1 potassium',                   'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'K',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1mg':           {'long_parameter_name':'PM1 magnesium',                   'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Mg',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1mn':           {'long_parameter_name':'PM1 manganese',                   'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Mn',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1na':           {'long_parameter_name':'PM1 sodium',                      'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Na',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1nh4':          {'long_parameter_name':'PM1 ammonium',                    'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'NH4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1ni':           {'long_parameter_name':'PM1 nickel',                      'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Ni',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1no3':          {'long_parameter_name':'PM1 nitrate',                     'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'NO3',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1oc':           {'long_parameter_name':'PM1 organic carbon',              'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1occorrected':  {'long_parameter_name':'PM1 organic carbon: corrected',   'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'C',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1pb':           {'long_parameter_name':'PM1 lead',                        'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Pb',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1se':           {'long_parameter_name':'PM1 selenium',                    'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Se',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1so4':          {'long_parameter_name':'PM1 sulphate',                    'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'SO4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1so4nss':       {'long_parameter_name':'PM1 sulphate : non-sea salt',     'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'SO4',     'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1v':            {'long_parameter_name':'PM1 vanadium',                    'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'V',       'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 },
'pm1zn':           {'long_parameter_name':'PM10 zinc',                       'matrix':'PM1',     'standard_units':'$\mu$g m$^{-3}$', 'chemical_formula':'Zn',      'extreme_lower_limit':0.0, 'extreme_upper_limit':1000.0 }
} 

#------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#

#generate GHOST interactive dashboard

qApp = QtWidgets.QApplication(sys.argv)
qApp.setStyle("Fusion")
GHOST_dash = generate_GHOST_interactive_dashboard('parallel')
sys.exit(qApp.exec_())
