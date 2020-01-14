#WRITTEN BY DENE BOWDALO

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

#providentia.py

#executable of the Providentia tool

###------------------------------------------------------------------------------------###
###IMPORT CONFIGURATION FILE VARIABLES
###------------------------------------------------------------------------------------###

from configuration import *

###------------------------------------------------------------------------------------###
###IMPORT BUILT-IN MODULES
###------------------------------------------------------------------------------------###
import bisect
from collections import OrderedDict
import copy
import datetime
from dateutil.relativedelta import relativedelta
from functools import partial
import gc
import glob
import multiprocessing
import os
import sys
from textwrap import wrap
import time
from weakref import WeakKeyDictionary

###------------------------------------------------------------------------------------###
###IMPORT THIRD-PARTY MODULES
###------------------------------------------------------------------------------------###
from netCDF4 import Dataset, num2date
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
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from PyQt5 import QtCore, QtWidgets, QtGui, Qt
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)
import scipy.stats
import seaborn as sns

###------------------------------------------------------------------------------------###
###IMPORT GHOST STANDARDS 
###------------------------------------------------------------------------------------###
sys.path.insert(1, '{}/GHOST_standards/{}'.format(obs_root,GHOST_version))
from GHOST_standards import standard_parameters, get_standard_metadata, standard_data_flag_name_to_data_flag_code, standard_QA_name_to_QA_code
#modify standard parameter dictionary to have BSC standard parameter names as keys (rather than GHOST)
parameter_dictionary = {}
for param, param_dict in standard_parameters.items():
    parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict  
#get standard metadata dictionary
standard_metadata = get_standard_metadata({'standard_units':''})
#create list of metadata variables to read (make global)
metadata_vars_to_read = [key for key in standard_metadata.keys() if pd.isnull(standard_metadata[key]['metadata_type']) == False]
metadata_dtype = [(key,standard_metadata[key]['data_type']) for key in metadata_vars_to_read]

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

#setup dictionary characterising formats for all GUI window objects (i.e. buttons, titles etc.)

formatting_dict = {'title_menu':               {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20, 'underline':True},
                   'label_menu':               {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'button_menu':              {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'checkbox_menu':            {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'combobox_menu':            {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'lineedit_menu':            {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'title_popup':              {'font':QtGui.QFont("SansSerif", 13),  'colour':'black', 'height':23, 'underline':True},
                   'button_popup':             {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'navigation_button_popup':  {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'checkbox_popup':           {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':18},
                   'rangebox_label_popup':     {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':18},
                   'rangebox_popup':           {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':18, 'max_width':90},
                   'column_header_label_popup':{'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':18, 'italic':True},
                   'tooltip':                  {'font':QtGui.QFont("SansSerif", 8)}            
                   }

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

def set_formatting(PyQt5_obj, formats):

    '''function that takes a PyQt5 object and applies some defined formatting'''

    #first get defined font for object
    defined_font = formats['font']

    #iterate through formats dictionary and apply defined font modifiers/object formatting values
    for format_name, format_val in formats.items():
        if format_name == 'bold':
            defined_font.setBold(format_val) 
       
        if format_name == 'italic':
            defined_font.setItalic(format_val)

        elif format_name == 'underline':
            defined_font.setUnderline(format_val)

        elif format_name == 'height':
            PyQt5_obj.setFixedHeight(format_val)

        elif format_name == 'width':
            PyQt5_obj.setFixedWidth(format_val)

        elif format_name == 'min_height':
            PyQt5_obj.setMinimumHeight(format_val)

        elif format_name == 'min_width':
            PyQt5_obj.setMinimumWidth(format_val)

        elif format_name == 'max_height':
            PyQt5_obj.setMaximumHeight(format_val)

        elif format_name == 'max_width':
            PyQt5_obj.setMaximumWidth(format_val)

        elif format_name == 'colour':
            PyQt5_obj.setStyleSheet("color: {};".format(format_val))


    #now apply font to object
    PyQt5_obj.setFont(defined_font)
        
    #return modified PyQt5 object
    return PyQt5_obj

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

def wrap_tooltip_text(tooltip_text, max_width):

    '''function which takes the text for a tooltip and wraps it by the screen pixel width.
       It does this by estimating the pixel width of the tooltip text (as formatted),
       and then gets the ratio exceedance over the screen pixel width.
       If there is an exceedance (i.e. > 1), the text is then broken into n max_char pieces 
       based on the position of the first exceedance in the text 
       (i.e. the part of the text which first exceeds the screen pixel width)
    ''' 

    tooltip_label = set_formatting(QtWidgets.QLabel(text = tooltip_text), formatting_dict['tooltip'])    
    tooltip_width = tooltip_label.fontMetrics().boundingRect(tooltip_label.text()).width()
    if tooltip_width > max_width:
        ratio = tooltip_width/max_width
        max_char = int(np.floor((len(tooltip_text)/ratio)*1.0))
        tooltip_text = '\n'.join(wrap(tooltip_text, max_char))

    return tooltip_text

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###                             }

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

class ComboBox(QtWidgets.QComboBox):

    '''modify default class of PyQT5 combobox to always dropdown from fixed position box postion, stopping truncation of data'''

    def showPopup(self):
        QtWidgets.QComboBox.showPopup(self)
        self.view().parent().move(self.mapToGlobal(QtCore.QPoint()))

class QVLine(QtWidgets.QFrame):

    '''define class that generates vertical separator line'''

    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QtWidgets.QFrame.VLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

class pop_up_window(QtWidgets.QWidget):

    '''define class that generates generalised pop-up window menu'''

    def __init__(self, menu_root, menu_levels, main_window_geometry):
        super(pop_up_window, self).__init__()
        
        #add input arguments to self
        self.menu_root = menu_root
        self.menu_levels = menu_levels        
        self.main_window_geometry = main_window_geometry
        self.menu_current = menu_root
        for menu_level_ii, menu_level in enumerate(menu_levels):
            self.menu_current = self.menu_current[menu_level]  
        
        #generate GUI window for root page in menu
        self.generate_window()

        #define stylesheet for tooltips
        self.setStyleSheet("QToolTip { font: %spt %s}"%(formatting_dict['tooltip']['font'].pointSizeF(), formatting_dict['tooltip']['font'].family())) 


    def generate_window(self):

        '''generate GUI window for current menu level'''
        
        #get current menu level keys
        menu_current_keys = list(self.menu_current.keys())

        #set window title
        self.setWindowTitle(self.menu_current['window_title'])

        #create parent layout 
        parent_layout = QtWidgets.QVBoxLayout()
        parent_layout.setAlignment(QtCore.Qt.AlignTop) 

        #define spacing/margin variables
        self.layout_spacing = 10
        parent_layout.setSpacing(self.layout_spacing)
        self.page_margin = 5
        parent_layout.setContentsMargins(self.page_margin,self.page_margin,self.page_margin,self.page_margin)

        #set page title 
        title_label = set_formatting(QtWidgets.QLabel(self, text = self.menu_current['page_title']), formatting_dict['title_popup'])
        title_label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignTop)

        #add title to parent frame
        parent_layout.addWidget(title_label)

        #create layout for placing buttons horizontally (aligned left)
        button_row = QtWidgets.QHBoxLayout()
        button_row.setAlignment(QtCore.Qt.AlignLeft)
        #initialise variable to update if have any buttons in button row
        self.have_buttons = False

        #check if need to create home/previous button
        if len(self.menu_levels) >= 1:
            self.have_buttons = True
            #create home button
            root_button = set_formatting(QtWidgets.QPushButton("HOME"), formatting_dict['button_popup']) 
            root_button.setStyleSheet("color: blue;")
            root_button.setToolTip('Return to home menu page')
            root_button.setFixedWidth(80)
            button_row.addWidget(root_button)
            root_button.clicked.connect(self.root_page)
            #create previous button
            previous_button = set_formatting(QtWidgets.QPushButton("PREVIOUS"), formatting_dict['button_popup']) 
            previous_button.setStyleSheet("color: green;")
            previous_button.setToolTip('Return to previous menu page')
            previous_button.setFixedWidth(80)
            button_row.addWidget(previous_button)
            previous_button.clicked.connect(self.previous_page)

        #check if need to create checkbox selection buttons
        if 'select_buttons' in menu_current_keys:
            self.have_buttons = True
            #need to create "Select All" button?
            if 'all' in self.menu_current['select_buttons']:
                select_all_button = set_formatting(QtWidgets.QPushButton("Select All"), formatting_dict['button_popup'])
                select_all_button.setFixedWidth(100)
                button_row.addWidget(select_all_button)
                select_all_button.clicked.connect(self.select_all)
            #need to create "Clear All" button? 
            if 'clear' in self.menu_current['select_buttons']:
                clear_all_button = set_formatting(QtWidgets.QPushButton("Clear All"), formatting_dict['button_popup'])
                clear_all_button.setFixedWidth(100)
                button_row.addWidget(clear_all_button)
                clear_all_button.clicked.connect(self.clear_all)  
            #need to create "Select Default" button? 
            if 'default' in self.menu_current['select_buttons']:
                select_default_button = set_formatting(QtWidgets.QPushButton("Select Default"), formatting_dict['button_popup'])
                select_default_button.setFixedWidth(100)
                button_row.addWidget(select_default_button)
                select_default_button.clicked.connect(self.select_all_default)

        #add button row to parent layout (if have some buttons)
        if self.have_buttons == True:
            parent_layout.addLayout(button_row)

        #create dictionary to store current page status
        self.page_memory = {}

        #gather list of all different menu types that need to be plotted
        menu_types = []
        if 'navigation_buttons' in menu_current_keys:
            menu_types.append('navigation_buttons')
        if 'rangeboxes' in menu_current_keys:
            menu_types.append('rangeboxes')
        if 'checkboxes' in menu_current_keys:
            menu_types.append('checkboxes')
        #have at least1 menu type?
        if len(menu_types) > 0:
            #create grid of menu types (horizontally concatenating different menu type grids)
            grid = self.create_grid(menu_types)     
            #add grid to parent layout
            parent_layout.addLayout(grid)

        #------------------------------------------------------------------------#
        #set finalised layout
        self.setLayout(parent_layout) 

        #set geometry to match that of main window
        self.setGeometry(self.main_window_geometry)

        #show pop-up window
        self.show()

        #------------------------------------------------------------------------#

        #setup event to get selected checkbox indices when closing window
        quit = QtWidgets.QAction("Quit", self)
        quit.triggered.connect(self.closeEvent)

    #------------------------------------------------------------------------#
    #------------------------------------------------------------------------#

    def create_grid(self, menu_types):
        
        '''create grid for each needed checkbox/rangebox/navigation button menu types, that wrap vertically
           and concatenate them horizontally together 
        '''
    
        #create horizontal layout to place all menu types within 
        horizontal_parent = QtWidgets.QHBoxLayout()
        #set spacing between different menu type grids
        horizontal_parent.setSpacing(25)
        #align grids to centre and top
        horizontal_parent.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignTop)
        
        #order appearance of menu types grids in menu (from left to right)
        menu_type_order_dict = {'navigation_buttons':1, 'rangeboxes':2, 'checkboxes':3}
        menu_types = sorted(menu_types, key=menu_type_order_dict.__getitem__) 

        #iterate through passed menu types
        for menu_type in menu_types:

            #create empty grid
            grid = QtWidgets.QGridLayout()
            #align grid to top
            grid.setAlignment(QtCore.Qt.AlignTop)

            #initialise indices used for indexing objects in grid
            start_row_n = 0
            row_n = 0
            column_n = 0

            #get dictionary nested inside current menu for menu type
            menu_current_type = self.menu_current[menu_type]

            #if have no labels for current menu type, continue to next menu type
            if len(menu_current_type['labels']) == 0:
                continue

            #get keys associated with current menu type in current menu level
            current_menu_keys = list(menu_current_type.keys())

            #create dictionary to store variables that save page status per menu type (also set formatting dict for each row in grid, and vertical spacing)
            if menu_type == 'checkboxes':
                row_format_dict = formatting_dict['checkbox_popup']
                grid_vertical_spacing = 0
                if ('keep_selected' in current_menu_keys) & ('remove_selected' in current_menu_keys):
                    self.page_memory['checkboxes'] = {'keep_selected':[],'remove_selected':[], 'n_column_consumed':3, 'ordered_elements':['keep_selected','remove_selected'], 'widget':QtWidgets.QCheckBox}
                elif 'keep_selected' in current_menu_keys:
                    self.page_memory['checkboxes'] = {'keep_selected':[],'n_column_consumed':2, 'ordered_elements':['keep_selected'], 'widget':QtWidgets.QCheckBox}
                elif 'remove_selected' in current_menu_keys:
                    self.page_memory['checkboxes'] = {'remove_selected':[],'n_column_consumed':2, 'ordered_elements':['remove_selected'], 'widget':QtWidgets.QCheckBox}
            elif menu_type == 'rangeboxes':
                row_format_dict = formatting_dict['rangebox_popup']
                grid_vertical_spacing = 3
                if ('current_lower' in current_menu_keys) & ('current_upper' in current_menu_keys):
                    self.page_memory['rangeboxes'] = {'current_lower':[],'current_upper':[], 'n_column_consumed':3, 'ordered_elements':['current_lower','current_upper'], 'widget':QtWidgets.QLineEdit}
                elif 'current_lower' in current_menu_keys:
                    self.page_memory['rangeboxes'] = {'current_lower':[],'n_column_consumed':2, 'ordered_elements':['current_lower'], 'widget':QtWidgets.QLineEdit}
                elif 'current_upper' in current_menu_keys:
                    self.page_memory['rangeboxes'] = {'current_upper':[],'n_column_consumed':2, 'ordered_elements':['current_upper'], 'widget':QtWidgets.QLineEdit}
            elif menu_type == 'navigation_buttons':
                row_format_dict = formatting_dict['navigation_button_popup']
                grid_vertical_spacing = 3
                self.page_memory['navigation_buttons'] = {'buttons':[], 'n_column_consumed':1, 'ordered_elements':['buttons'], 'widget':QtWidgets.QPushButton}
        
            #if have more than 1 column per label, need column headers
            if len(self.page_memory[menu_type]['ordered_elements']) > 1: 
                have_column_headers = True
                start_row_n = 1
            else:
                have_column_headers = False

            #set vertical/horizontal grid spacing
            grid.setHorizontalSpacing(3)
            grid.setVerticalSpacing(grid_vertical_spacing)

            #calculate currently occupied vertical space
            occupied_vertical_space_before_grid = self.page_margin + formatting_dict['title_popup']['height'] + self.layout_spacing + self.page_margin
            if self.have_buttons == True:
                occupied_vertical_space_before_grid += (formatting_dict['button_popup']['height'] + self.layout_spacing)
            if have_column_headers == True:
                occupied_vertical_space_before_grid += (formatting_dict['column_header_label_popup']['height'] + grid_vertical_spacing)

            #initialise variable for tracking available vertical space when appending rows of grid
            currently_occupied_vertical_space = copy.deepcopy(occupied_vertical_space_before_grid)

            #iterate through all grid labels
            for label_ii, label in enumerate(menu_current_type['labels']):

                #evaluate if all available vertical space has been consumed
                row_available_space = self.main_window_geometry.height() - currently_occupied_vertical_space
                #if available space <= than row height, force a new column to be started
                if row_available_space <= (row_format_dict['height']):
                    column_n+=self.page_memory[menu_type]['n_column_consumed']
                    row_n = 0
                    currently_occupied_vertical_space = copy.deepcopy(occupied_vertical_space_before_grid)

                #if menu type == 'rangeboxes', then need to add label to left of rangeboxes, also add a tooltip (if exist)
                if (menu_type == 'checkboxes') or (menu_type == 'rangeboxes'):
                    rangebox_label = set_formatting(QtWidgets.QLabel(self, text = label), formatting_dict['rangebox_label_popup']) 
                    if menu_type == 'rangeboxes':
                        if len(menu_current_type['tooltips']) > 0:
                            rangebox_label.setToolTip(wrap_tooltip_text(menu_current_type['tooltips'][label_ii], self.main_window_geometry.width()))
                    grid.addWidget(rangebox_label, start_row_n+row_n, column_n, QtCore.Qt.AlignLeft) 

                #create all elements in column, per row
                for element_ii, element in enumerate(self.page_memory[menu_type]['ordered_elements']):
                    #if menu type == 'rangeboxes' then add 1 to element ii, because placed a label in first column
                    if (menu_type == 'checkboxes') or (menu_type == 'rangeboxes'):
                        element_label = ''
                        element_ii+=1
                    else:
                        element_label = label

                    #append widget to page memory dictionary 
                    self.page_memory[menu_type][element].append(set_formatting(self.page_memory[menu_type]['widget'](element_label), row_format_dict))  
                    
                    #put text label to left of keep checkbox/rangebox (rather than to right)
                    if menu_type in ['checkboxes','rangeboxes']:
                        #check if checkbox is currently selected, if so select it again
                        #map checkbox value first if necessary
                        if menu_type == 'checkboxes':                   
                            if 'map_vars' in current_menu_keys:
                                var_to_check = menu_current_type['map_vars'][label_ii]
                            else:
                                var_to_check = copy.deepcopy(label)   
                            if var_to_check in menu_current_type[element]: 
                                self.page_memory[menu_type][element][label_ii].setCheckState(QtCore.Qt.Checked)
                        #set rangeboxes to previous set value (if any)
                        elif menu_type == 'rangeboxes':
                            self.page_memory[menu_type][element][label_ii].setText(menu_current_type[element][label_ii])

                    #if menu type == navigation_buttons, add connectivity to buttons, and also add tooltip
                    elif menu_type == 'navigation_buttons':
                        if len(menu_current_type['tooltips']) > 0:
                            self.page_memory[menu_type][element][label_ii].setToolTip(wrap_tooltip_text(menu_current_type['tooltips'][label_ii], self.main_window_geometry.width()))
                        self.page_memory[menu_type][element][label_ii].clicked.connect(self.open_new_page)

                    #add element to grid (aligned left)
                    grid.addWidget(self.page_memory[menu_type][element][label_ii], start_row_n+row_n, column_n+element_ii, QtCore.Qt.AlignLeft)                               

                #iterate row_n
                row_n +=1

                #add row vertical space to total occupied space
                currently_occupied_vertical_space += (row_format_dict['height'] + grid_vertical_spacing)

            #add column headers to menu type grid if needed
            if have_column_headers == True:
                for column_number in np.arange(0, column_n+1, self.page_memory[menu_type]['n_column_consumed']):
                    if menu_type == 'checkboxes':
                        grid.addWidget(set_formatting(QtWidgets.QLabel(self, text = 'K'), formatting_dict['column_header_label_popup']), 0, column_number+1, QtCore.Qt.AlignCenter) 
                        grid.addWidget(set_formatting(QtWidgets.QLabel(self, text = 'R'), formatting_dict['column_header_label_popup']), 0, column_number+2, QtCore.Qt.AlignCenter)
                    elif menu_type == 'rangeboxes':
                        grid.addWidget(set_formatting(QtWidgets.QLabel(self, text = 'Lower Bound'), formatting_dict['column_header_label_popup']), 0, column_number+1, QtCore.Qt.AlignCenter) 
                        grid.addWidget(set_formatting(QtWidgets.QLabel(self, text = 'Upper Bound'), formatting_dict['column_header_label_popup']), 0, column_number+2, QtCore.Qt.AlignCenter)  

            #add menu type grid to horizontal layout
            horizontal_parent.addLayout(grid)

        #return horizontally concatenated menu type grids
        return horizontal_parent

    #------------------------------------------------------------------------#
    #------------------------------------------------------------------------#
    #functions that handle callbacks upon clicking on buttons 

    def open_new_page(self):
        '''function to open new page in pop-up window'''
        
        #get selected navigation button text
        selected_navigation_button = self.sender().text()
        #add selected navigation button text to menu levels list
        self.menu_levels.append(selected_navigation_button)
        #create new pop-up page for selected navigation button
        self.new_window = pop_up_window(self.menu_root, self.menu_levels, self.main_window_geometry)
        #sleep briefly to allow new page to be generated
        time.sleep(0.1)
        #close current pop-up page
        self.close()

    def root_page(self):
        '''function that returns pop-up window to root menu level page'''

        #create new pop-up page for root menu level
        self.new_window = pop_up_window(self.menu_root, [], self.main_window_geometry)
        #sleep briefly to allow new page to be generated
        time.sleep(0.1)
        #close current pop-up page
        self.close()

    def previous_page(self):
        '''function that returns pop-up window to previous menu level page'''

        #create new pop-up page for previous menu level
        self.new_window = pop_up_window(self.menu_root, self.menu_levels[:-1], self.main_window_geometry)
        #sleep briefly to allow new page to be generated
        time.sleep(0.1)
        #close current pop-up page
        self.close()

    def select_all(self):
        '''function to select all checkboxes'''
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for checkbox_ii, checkbox in enumerate(self.page_memory['checkboxes'][element]):
                self.page_memory['checkboxes'][element][checkbox_ii].setCheckState(QtCore.Qt.Checked) 
    
    def clear_all(self):
        '''function to clear all checkboxes'''
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for checkbox_ii, checkbox in enumerate(self.page_memory['checkboxes'][element]):
                self.page_memory['checkboxes'][element][checkbox_ii].setCheckState(QtCore.Qt.Unchecked) 
                
    def select_all_default(self):
        '''function to select all default selected checkboxes'''
        #unselect all checkboxes first 
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for checkbox_ii, checkbox in enumerate(self.page_memory['checkboxes'][element]):
                self.page_memory['checkboxes'][element][checkbox_ii].setCheckState(QtCore.Qt.Unchecked) 
        
        #map default variables to positional indices
        if 'map_vars' in list(self.menu_current['checkboxes'].keys()):
            default_inds = [np.where(self.menu_current['checkboxes']['map_vars'] == default_var)[0][0] for default_var in self.menu_current['checkboxes']['remove_default']]
        else:
            default_inds = [np.where(self.menu_current['checkboxes']['labels'] == default_var)[0][0] for default_var in self.menu_current['checkboxes']['remove_default']]

        #now select only desired default checkboxes
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for default_ind in default_inds:
                self.page_memory['checkboxes'][element][default_ind].setCheckState(QtCore.Qt.Checked)

    #------------------------------------------------------------------------#
    #------------------------------------------------------------------------#

    def closeEvent(self, event):

        '''function to get status of current page upon closing of pop-up window'''

        #take everything from page memory dictionary and put it back into menu level object 
        
        #iterate through menu types
        for menu_type in list(self.page_memory.keys()):
            for element in self.page_memory[menu_type]['ordered_elements']:
                #if menu type == 'checkboxes', get variable names of all checkboxes ticked
                if menu_type == 'checkboxes':
                    selected_vars = []
                    for checkbox_ii, checkbox in enumerate(self.page_memory[menu_type][element]):
                        if checkbox.checkState() == QtCore.Qt.Checked:
                            #map selected position index to variable name
                            if 'map_vars' in list(self.menu_current[menu_type].keys()):
                                selected_vars.append(self.menu_current[menu_type]['map_vars'][checkbox_ii])       
                            else:
                                selected_vars.append(self.menu_current[menu_type]['labels'][checkbox_ii])
                    #update previous selected variable
                    if element == 'keep_selected':
                        self.menu_current[menu_type]['previous_keep_selected'] = copy.deepcopy(self.menu_current[menu_type]['keep_selected'])
                    elif element == 'remove_selected':
                        self.menu_current[menu_type]['previous_keep_selected'] = copy.deepcopy(self.menu_current[menu_type]['remove_selected'])
                    #update selected variable
                    self.menu_current[menu_type][element] = selected_vars 
                #if menu type == 'rangeboxes', get current values of all rangeboxes
                if menu_type == 'rangeboxes':
                    set_vals = []
                    for rangebox in self.page_memory[menu_type][element]:
                        set_vals.append(rangebox.text())
                    #update previous set variable
                    if element == 'current_lower':
                        self.menu_current[menu_type]['previous_lower'] = copy.deepcopy(self.menu_current[menu_type]['current_lower'])
                    elif element == 'current_upper':
                        self.menu_current[menu_type]['previous_upper'] = copy.deepcopy(self.menu_current[menu_type]['current_upper'])
                    #update set value
                    self.menu_current[menu_type][element] = set_vals

#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#
class Providentia_main_window(QtWidgets.QWidget):

    '''define class that generates Providentia main window'''     
 
    #create signals that are fired upon resizing/moving of main Providentia window
    resized = QtCore.pyqtSignal()
    move = QtCore.pyqtSignal()   
 
    def __init__(self, read_type):
        super(Providentia_main_window, self).__init__()
        
        #put read_type into self
        self.read_type = read_type

        #create UI
        self.initUI()
    
        #setup callback events upon resizing/moving of Providentia window 
        self.resized.connect(self.get_geometry)
        self.move.connect(self.get_geometry)

    #------------------------------------------------------------------------#

    def resizeEvent(self, event):
        '''Function to overwrite default PyQt5 resizeEvent function --> for calling get_geometry'''
        self.resized.emit()
        return super(Providentia_main_window, self).resizeEvent(event)

    def moveEvent(self, event):
        '''Function to overwrite default PyQt5 moveEvent function --> for calling get_geometry'''
        self.move.emit()
        return super(Providentia_main_window, self).moveEvent(event)

    def get_geometry(self):
        '''Get current geometry of main Providentia window'''
        self.main_window_geometry = copy.deepcopy(self.geometry())
 
    def initUI(self):

        '''initialise Providentia main window user interface'''

        #set window title
        self.window_title = "Providentia"
        self.setWindowTitle(self.window_title)

        #create parent layout to pull together a configuration bar, a MPL navigation toolbar, and a MPL canvas of plots
        parent_layout = QtWidgets.QVBoxLayout()
        
        #define spacing/margin variables
        parent_layout.setSpacing(0)
        parent_layout.setContentsMargins(0,0,0,0)

        #define stylesheet for tooltips
        self.setStyleSheet("QToolTip { font: %spt %s}"%(formatting_dict['tooltip']['font'].pointSizeF(), formatting_dict['tooltip']['font'].family())) 

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
        self.lb_data_selection = set_formatting(QtWidgets.QLabel(self, text = "Data Selection"), formatting_dict['title_menu'])
        self.lb_data_selection.setToolTip('Setup configuration of data to read into memory')
        self.bu_read = set_formatting(QtWidgets.QPushButton('READ', self), formatting_dict['button_menu'])
        self.bu_read.setFixedWidth(40)
        self.bu_read.setStyleSheet("color: red;")
        self.bu_read.setToolTip('Read selected configuration of data into memory')
        self.ch_colocate = set_formatting(QtWidgets.QCheckBox("Colocate"), formatting_dict['checkbox_menu'])
        self.ch_colocate.setToolTip('Temporally colocate observational/experiment data')
        self.cb_network = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_network.setFixedWidth(95)
        self.cb_network.setToolTip('Select providing observational data network')
        self.cb_resolution = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_resolution.setFixedWidth(95)
        self.cb_resolution.setToolTip('Select temporal resolution of data')
        self.cb_matrix = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_matrix.setFixedWidth(95)
        self.cb_matrix.setToolTip('Select data matrix')
        self.cb_species = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_species.setFixedWidth(95)
        self.cb_species.setToolTip('Select species')
        self.le_start_date = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_start_date.setFixedWidth(70)
        self.le_start_date.setToolTip('Set data start date: YYYYMMDD')
        self.le_end_date = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_end_date.setFixedWidth(70)
        self.le_end_date.setToolTip('Set data end date: YYYYMMDD')
        self.bu_QA = set_formatting(QtWidgets.QPushButton('QA', self), formatting_dict['button_menu'])
        self.bu_QA.setFixedWidth(46)
        self.bu_QA.setToolTip('Select standardised quality assurance flags to filter by')
        self.bu_flags = set_formatting(QtWidgets.QPushButton('FLAGS', self), formatting_dict['button_menu'])
        self.bu_flags.setFixedWidth(46)
        self.bu_flags.setToolTip('Select standardised data reporter provided flags to filter by')
        self.bu_experiments = set_formatting(QtWidgets.QPushButton('EXPS', self), formatting_dict['button_menu'])
        self.bu_experiments.setFixedWidth(46)
        self.bu_experiments.setToolTip('Select experiment/s data to read')
        self.vertical_splitter_1 = QVLine()
        self.vertical_splitter_1.setMaximumWidth(20)
        self.lb_data_filter = set_formatting(QtWidgets.QLabel(self, text = "Data Filter"), formatting_dict['title_menu'])
        self.lb_data_filter.setFixedWidth(65)
        self.lb_data_filter.setToolTip('Select criteria to filter data by')
        self.bu_rep = set_formatting(QtWidgets.QPushButton('% REP', self), formatting_dict['button_menu'])
        self.bu_rep.setFixedWidth(46)
        self.bu_rep.setToolTip('Select % desired representativity in data across whole record and for specific temporal periods')
        self.bu_meta = set_formatting(QtWidgets.QPushButton('META', self), formatting_dict['button_menu'])
        self.bu_meta.setFixedWidth(46)
        self.bu_meta.setToolTip('Select metadata to filter by')
        self.bu_period = set_formatting(QtWidgets.QPushButton('PERIOD', self), formatting_dict['button_menu'])
        self.bu_period.setFixedWidth(50)
        self.bu_period.setToolTip('Select data in specific periods')
        self.bu_screen = set_formatting(QtWidgets.QPushButton('FILTER', self), formatting_dict['button_menu'])
        self.bu_screen.setFixedWidth(50)
        self.bu_screen.setStyleSheet("color: blue;")
        self.bu_screen.setToolTip('Filter data')
        self.lb_data_bounds = set_formatting(QtWidgets.QLabel(self, text = "Bounds"), formatting_dict['label_menu'])
        self.lb_data_bounds.setFixedWidth(47)
        self.lb_data_bounds.setToolTip('Set lower/upper bounds of data')
        self.le_minimum_value = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_minimum_value.setFixedWidth(60)
        self.le_minimum_value.setToolTip('Set lower bound of data')
        self.le_maximum_value = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_maximum_value.setFixedWidth(60)
        self.le_maximum_value.setToolTip('Set upper bound of data')
        self.vertical_splitter_2 = QVLine()
        self.vertical_splitter_2.setMaximumWidth(20)
        self.lb_z = set_formatting(QtWidgets.QLabel(self, text = "Map Z"), formatting_dict['title_menu'])
        self.lb_z.setToolTip('Set map Z statistic')
        self.cb_z_stat = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_z_stat.setFixedWidth(80)
        self.cb_z_stat.setToolTip('Select map Z statistic')
        self.cb_z1 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_z1.setFixedWidth(125)
        self.cb_z1.setToolTip('Select Z1 dataset')
        self.cb_z2 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_z2.setFixedWidth(125)
        self.cb_z2.setToolTip('Select Z2 dataset')
        self.vertical_splitter_3 = QVLine()
        self.vertical_splitter_3.setMaximumWidth(20)
        self.lb_experiment_bias = set_formatting(QtWidgets.QLabel(self, text = "Exp. Bias"), formatting_dict['title_menu'])
        self.lb_experiment_bias.setToolTip('Set experiment bias statistic')
        self.cb_experiment_bias_type = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_experiment_bias_type.setFixedWidth(100)
        self.cb_experiment_bias_type.setToolTip('Select experiment bias type')
        self.cb_experiment_bias_stat = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_experiment_bias_stat.setFixedWidth(100)
        self.cb_experiment_bias_stat.setToolTip('Select experiment bias statistic')
        self.vertical_splitter_4 = QVLine()
        self.vertical_splitter_4.setMaximumWidth(20)
        self.lb_station_selection = set_formatting(QtWidgets.QLabel(self, text = "Site Select"), formatting_dict['title_menu'])
        self.lb_station_selection.setToolTip('Select stations')
        self.ch_select_all = set_formatting(QtWidgets.QCheckBox("All"), formatting_dict['checkbox_menu'])
        self.ch_select_all.setToolTip('Select all stations')
        self.ch_intersect = set_formatting(QtWidgets.QCheckBox("Intersect"), formatting_dict['checkbox_menu'])
        self.ch_intersect.setToolTip('Select stations that intersect with all loaded model domains')

        #position objects on gridded configuration bar
        config_bar.addWidget(self.lb_data_selection, 0, 0, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_colocate, 0, 1, QtCore.Qt.AlignCenter)
        config_bar.addWidget(self.bu_read, 0, 2, QtCore.Qt.AlignCenter)  
        config_bar.addWidget(self.cb_network, 1, 0)
        config_bar.addWidget(self.cb_resolution, 2, 0)
        config_bar.addWidget(self.cb_matrix, 1, 1)
        config_bar.addWidget(self.cb_species, 2, 1)
        config_bar.addWidget(self.le_start_date, 1, 2)
        config_bar.addWidget(self.le_end_date, 2, 2) 
        config_bar.addWidget(self.bu_QA, 0, 3) 
        config_bar.addWidget(self.bu_flags, 1, 3) 
        config_bar.addWidget(self.bu_experiments, 2, 3)
        config_bar.addWidget(self.vertical_splitter_1, 0, 5, 3, 1)
        config_bar.addWidget(self.lb_data_filter, 0, 6, 1, 2, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_rep, 1, 6)
        config_bar.addWidget(self.bu_meta, 2, 6)
        config_bar.addWidget(self.bu_period, 1, 7)
        config_bar.addWidget(self.bu_screen, 2, 7) 
        config_bar.addWidget(self.lb_data_bounds, 0, 8, QtCore.Qt.AlignCenter)
        config_bar.addWidget(self.le_minimum_value, 1, 8)
        config_bar.addWidget(self.le_maximum_value, 2, 8)
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

        #setup pop-up window menu tree for flags
        self.flag_menu = {'window_title':'FLAGS', 'page_title':'Select standardised data reporter provided flags to filter by', 'checkboxes':{}}
        self.flag_menu['checkboxes']['labels'] = np.array(sorted(standard_data_flag_name_to_data_flag_code, key=standard_data_flag_name_to_data_flag_code.get))
        self.flag_menu['checkboxes']['remove_default'] = np.array([1, 2, 3, 10, 11, 12, 13, 14, 15, 16, 20, 21, 24, 25, 26, 29, 30, 31, 32, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 90, 150, 154, 155, 156, 157], dtype=np.uint8)
        self.flag_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
        self.flag_menu['checkboxes']['map_vars'] = np.sort(list(standard_data_flag_name_to_data_flag_code.values()))
        self.flag_menu['select_buttons'] = ['all','clear','default']

        #setup pop-up window menu tree for qa
        self.qa_menu = {'window_title':'QA', 'page_title':'Select standardised data reporter provided flags to filter by', 'checkboxes':{}}
        self.qa_menu['checkboxes']['labels'] = np.array(sorted(standard_QA_name_to_QA_code, key=standard_QA_name_to_QA_code.get))
        self.qa_menu['checkboxes']['remove_default'] = np.array([0, 1, 2, 3, 4, 5, 7, 8, 10, 12, 13, 14, 17, 18, 22, 25, 30, 40, 41, 42], dtype=np.uint8)
        self.qa_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
        self.qa_menu['checkboxes']['map_vars'] = np.sort(list(standard_QA_name_to_QA_code.values()))
        self.qa_menu['select_buttons'] = ['all','clear','default'] 

        #setup pop-up window menu tree for experiments
        self.experiments_menu = {'window_title':'EXPERIMENTS', 'page_title':'Select Experiment/s', 'checkboxes':{}}
        self.experiments_menu['checkboxes']['labels'] = []
        self.experiments_menu['checkboxes']['keep_default'] = []
        self.experiments_menu['checkboxes']['keep_selected'] = []
        self.experiments_menu['checkboxes']['map_vars'] = []
        self.experiments_menu['select_buttons'] = ['all','clear','default']   

        #setup pop-up window menu tree for metadata
        self.metadata_types =  {'STATION POSITION':'Filter stations by measurement position', 
                                'STATION CLASSIFICATIONS':'Filter stations by station provided classifications', 
                                'STATION MISCELLANEOUS':'Filter stations by miscellaneous station provided metadata', 
                                'GLOBALLY GRIDDED CLASSIFICATIONS':'Filter stations by globally gridded classifications', 
                                'MEASUREMENT PROCESS INFORMATION':'Filter stations by measurement process information'}
        self.metadata_menu = {'window_title':'METADATA', 'page_title':'Select metadata type to filter stations by', 'navigation_buttons':{}}
        self.metadata_menu['navigation_buttons']['labels'] = list(self.metadata_types.keys())
        self.metadata_menu['navigation_buttons']['tooltips'] = [self.metadata_types[key] for key in self.metadata_menu['navigation_buttons']['labels']]       
        for metadata_type_ii, metadata_type in enumerate(self.metadata_menu['navigation_buttons']['labels']):
            self.metadata_menu[metadata_type] = {'window_title':metadata_type,'page_title':self.metadata_menu['navigation_buttons']['tooltips'][metadata_type_ii], 'navigation_buttons':{}, 'rangeboxes':{}}
            self.metadata_menu[metadata_type]['navigation_buttons']['labels'] = [metadata_name for metadata_name in standard_metadata.keys() if (standard_metadata[metadata_name]['metadata_type'] == metadata_type) & (standard_metadata[metadata_name]['data_type'] == np.object)]
            self.metadata_menu[metadata_type]['navigation_buttons']['tooltips'] = [standard_metadata[metadata_name]['description'] for metadata_name in self.metadata_menu[metadata_type]['navigation_buttons']['labels']]
            for label in self.metadata_menu[metadata_type]['navigation_buttons']['labels']:
                self.metadata_menu[metadata_type][label] = {'window_title':label,'page_title':'Filter stations by unique {} metadata'.format(label), 'checkboxes':{}}  
                self.metadata_menu[metadata_type][label]['checkboxes']['labels'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['keep_selected'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['remove_selected'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['previous_keep_selected'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['previous_remove_selected'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['keep_default'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['remove_default'] = []
            #self.metadata_menu[metadata_type]['rangeboxes']['previous_labels'] = []
            self.metadata_menu[metadata_type]['rangeboxes']['labels'] = [metadata_name for metadata_name in standard_metadata.keys() if (standard_metadata[metadata_name]['metadata_type'] == metadata_type) & (standard_metadata[metadata_name]['data_type'] != np.object)]
            self.metadata_menu[metadata_type]['rangeboxes']['tooltips'] = [standard_metadata[metadata_name]['description'] for metadata_name in self.metadata_menu[metadata_type]['rangeboxes']['labels']]
            self.metadata_menu[metadata_type]['rangeboxes']['current_lower'] = [''] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels']) 
            self.metadata_menu[metadata_type]['rangeboxes']['current_upper'] = [''] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels']) 
            self.metadata_menu[metadata_type]['rangeboxes']['previous_lower'] = [''] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels']) 
            self.metadata_menu[metadata_type]['rangeboxes']['previous_upper'] = [''] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels']) 
            self.metadata_menu[metadata_type]['rangeboxes']['lower_default'] = [''] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels']) 
            self.metadata_menu[metadata_type]['rangeboxes']['upper_default'] = [''] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels']) 

        #setup pop-up window menu tree for % data representativity 
        self.representativity_menu = {'window_title':'% DATA REPRESENTATIVITY', 'page_title':'Select % Data Representativity Bounds', 'rangeboxes':{}}
        #self.representativity_menu['rangeboxes']['previous_labels'] = [] 
        self.representativity_menu['rangeboxes']['labels'] = [] 
        self.representativity_menu['rangeboxes']['tooltips'] = [] 
        self.representativity_menu['rangeboxes']['current_lower'] = []   
        #self.representativity_menu['rangeboxes']['current_upper'] = [] 

        #setup pop-up window menu tree for data periods
        self.period_menu = {'window_title':'DATA PERIOD', 'page_title':'Select Data Periods', 'checkboxes':{}}
        self.period_menu['checkboxes']['labels'] = []
        self.period_menu['checkboxes']['keep_selected'] = []
        self.period_menu['checkboxes']['remove_selected'] = []
        
        #enable pop up configuration windows
        self.bu_flags.clicked.connect(partial(self.generate_pop_up_window,self.flag_menu))
        self.bu_QA.clicked.connect(partial(self.generate_pop_up_window,self.qa_menu))
        self.bu_experiments.clicked.connect(partial(self.generate_pop_up_window,self.experiments_menu)) 
        self.bu_meta.clicked.connect(partial(self.generate_pop_up_window,self.metadata_menu)) 
        self.bu_rep.clicked.connect(partial(self.generate_pop_up_window,self.representativity_menu)) 
        self.bu_period.clicked.connect(partial(self.generate_pop_up_window,self.period_menu)) 
        
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

    def generate_pop_up_window(self, menu_root):

        '''generate pop up window'''
        
        self.pop_up_window = pop_up_window(menu_root, [], self.main_window_geometry)

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

            #set selected/active values of other fields to be initially None
            self.selected_network = None
            self.active_network = None
            self.selected_resolution = None
            self.active_resolution = None
            self.selected_matrix = None
            self.active_matrix = None
            self.selected_species = None
            self.active_species = None    

            #set active values of variables associated with pop up windows to be empty lists
            self.active_experiment_grids = []  
            self.active_qa = []
            self.active_flags = []  

            #set initial time array to be None
            self.time_array = None      

            self.relevant_yearmonths = None

            #set initial station references to be empty list
            self.station_references = []
        
            #--------------------------------------------------------------------------------------#
            #gather all observational data
            #create nested dictionary storing all observational species data by species matrix, by temporal resolution, by network, associated with list of start YYYYMM yearmonths of data files

            self.all_observation_data = {}

            #set all available networks
            available_networks = ['EBAS','EEA_AQ_eReporting','NCDC_ISD']
        
            #set all available temporal resolutions
            available_resolutions = ['hourly','hourly_instantaneous','daily','monthly']

            #iterate through available networks
            for network in available_networks:

                #check if directory for network exists
                #if not, continue
                if not os.path.exists('%s/%s/%s'%(obs_root,network,GHOST_version)):       
                    continue
    
                #write empty dictionary for network
                self.all_observation_data[network] = {}

                #iterate through available resolutions
                for resolution in available_resolutions:

                    #check if directory for resolution exists
                    #if not, continue
                    if not os.path.exists('%s/%s/%s/%s'%(obs_root,network,GHOST_version,resolution)):       
                        continue

                    #write nested empty dictionary for resolution
                    self.all_observation_data[network][resolution] = {}

                    #get available species for network/resolution
                    available_species = os.listdir('%s/%s/%s/%s'%(obs_root, network, GHOST_version, resolution))

                    #iterate through available files per species
                    for species in available_species:   

                        #get all netCDF monthly files per species
                        species_files = os.listdir('%s/%s/%s/%s/%s'%(obs_root, network, GHOST_version, resolution, species))

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
        else:
            self.selected_network = self.cb_network.currentText()
        
        #update resolution field
        available_resolutions = list(self.available_observation_data[self.cb_network.currentText()].keys())
        #manually force order of available resolutions
        resolution_order_dict = {'hourly':1, 'hourly_instantaneous':2, 'daily':3, 'monthly':4}
        available_resolutions = sorted(available_resolutions, key=resolution_order_dict.__getitem__) 
        self.cb_resolution.addItems(available_resolutions)
        if self.selected_resolution in available_resolutions:
            self.cb_resolution.setCurrentText(self.selected_resolution)
        else:
            self.selected_resolution = self.cb_resolution.currentText()

        #update matrix field
        available_matrices = sorted(self.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()])
        self.cb_matrix.addItems(available_matrices)
        if self.selected_matrix in available_matrices:
            self.cb_matrix.setCurrentText(self.selected_matrix) 
        else:
            self.selected_matrix = self.cb_matrix.currentText()

        #update species field
        available_species = sorted(self.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()][self.cb_matrix.currentText()])
        self.cb_species.addItems(available_species)
        if self.selected_species in available_species:
            self.cb_species.setCurrentText(self.selected_species) 
        else:   
            self.selected_species = self.cb_species.currentText()  

        #update available experiment data dictionary 
        self.get_valid_experiment_files_in_date_range()
        #update selected experiments -- keeping previously selected experiments if available
        #set selected indices as previously selected indices in current available list of experiments
        self.experiments_menu['checkboxes']['keep_selected'] = [previous_selected_experiment for previous_selected_experiment in self.experiments_menu['checkboxes']['keep_selected'] if previous_selected_experiment in self.experiments_menu['checkboxes']['map_vars']]

        #if initialising config bar then check default selection of QA checkboxes
        if self.config_bar_initialisation == True:
            self.qa_menu['checkboxes']['remove_selected'] = copy.deepcopy(self.qa_menu['checkboxes']['remove_default'])

        #unset variable to allow interactive handling from now
        self.block_config_bar_handling_updates = False

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#
    
    def get_valid_obs_files_in_date_range(self):
        
        '''define function that iterates through observational dictionary tree and returns a dictionary of available data in the selected date range'''

        #create dictionary to store available observational data
        self.available_observation_data = {}

        #check if start/end date are valid values, if not, return with no valid obs. files
        selected_start_date = self.le_start_date.text()
        selected_end_date = self.le_end_date.text()
        if (self.valid_date(selected_start_date)) & (self.valid_date(selected_end_date)):
            self.date_range_has_changed = True
            self.selected_start_date = int(selected_start_date)
            self.selected_end_date = int(selected_end_date)
            self.selected_start_date_firstdayofmonth = int(str(self.selected_start_date)[:6]+'01')
        else:
            return

        #check end date is > start date, if not, return with no valid obs. files
        if self.selected_start_date >= self.selected_end_date:
            return

        #check start date and end date are both within if valid date range (19000101 - 20500101), if not, return with no valid obs. files
        if (self.selected_start_date < 19000101) or (self.selected_end_date < 19000101) or (self.selected_start_date >= 20500101) or (self.selected_end_date >= 20500101):
            return 
        
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

    def get_valid_experiment_files_in_date_range(self):

        '''define function which gathers available experiment data for selected network/resolution/species.
           a dictionary is created storing available experiment-grid names associated with valid files in set date range.
        '''
    
        #create dictionary to store available experiment information
        self.available_experiment_data = {}

        #get all different experiment names
        available_experiments = os.listdir('%s/%s'%(exp_root,GHOST_version))          

        #iterate through available experiments
        for experiment in available_experiments:      
 
            #get all available grid types by experiment 
            available_grids = os.listdir('%s/%s/%s'%(exp_root,GHOST_version,experiment))

            #iterate through all available grids
            for grid in available_grids:
            
                #test first if interpolated directory exists before trying to get files from it
                #if it does not exit, continue
                if not os.path.exists('%s/%s/%s/%s/%s/%s/%s'%(exp_root,GHOST_version,experiment,grid,self.selected_resolution,self.selected_species,self.selected_network)):       
                    continue
                else:
                    #get all experiment netCDF files by experiment/grid/selected resolution/selected species/selected network
                    network_files = os.listdir('%s/%s/%s/%s/%s/%s/%s'%(exp_root,GHOST_version,experiment,grid,self.selected_resolution,self.selected_species,self.selected_network))
                    #get start YYYYMM yearmonths of data files
                    network_files_yearmonths = [int(f.split('_')[-1][:6]+'01') for f in network_files] 
                    #limit data files to just those within date range
                    valid_network_files_yearmonths = [ym for ym in network_files_yearmonths if (ym >= self.selected_start_date_firstdayofmonth) & (ym < self.selected_end_date)]
                    #if have some valid data files for experiment-grid, add experiment grid (with associated yearmonths) to dictionary
                    if len(valid_network_files_yearmonths) > 0:
                        self.available_experiment_data['%s-%s'%(experiment,grid)] = valid_network_files_yearmonths

        #get list of available experiment-grid names
        self.experiments_menu['checkboxes']['labels'] = np.array(sorted(list(self.available_experiment_data.keys())))
        self.experiments_menu['checkboxes']['map_vars'] = copy.deepcopy(self.experiments_menu['checkboxes']['labels'])

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

            #check if start date/end date have changed
            if (event_source == self.le_start_date) or (event_source == self.le_end_date):
                self.date_range_has_changed = True
                    
            #update configuration bar fields
            self.update_configuration_bar_fields()  

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#

    def valid_date(self, date_text):
 
        '''define function that determines if a date string is in the correct format'''

        try:
            datetime.datetime.strptime(date_text, '%Y%m%d')
            return True
        except:
            return False

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
        self.previous_active_qa = self.active_qa
        self.previous_active_flags = self.active_flags
        
        #set all currently selected variables as active variables
        self.active_network = self.selected_network
        self.active_resolution = self.selected_resolution
        self.active_matrix = self.selected_matrix
        self.active_species = self.selected_species
        self.active_start_date = self.selected_start_date
        self.active_end_date = self.selected_end_date
        self.active_experiment_grids = copy.deepcopy(self.experiments_menu['checkboxes']['keep_selected'])
        self.active_qa = copy.deepcopy(self.qa_menu['checkboxes']['remove_selected'])
        self.active_flags = copy.deepcopy(self.flag_menu['checkboxes']['remove_selected'])
        
        #--------------------------------------------------------------------#
        #determine what data (if any) needs to be read

        #set variables that inform what data needs to be read (set all initially as False)
        read_all = False
        read_left = False
        read_right = False
        cut_left = False
        cut_right = False
        
        #determine if any of the key variables have changed 
        #(network, resolution, species, qa, flags)
        #if any have changed, observations and any selected experiments have to be re-read entirely
        if (self.active_network != self.previous_active_network) or (self.active_resolution != self.previous_active_resolution) or (self.active_species != self.previous_active_species) or (np.array_equal(self.active_qa,self.previous_active_qa) == False) or (np.array_equal(self.active_flags,self.previous_active_flags) == False):
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

            st = time.time()

            #set new active time array/unique station references/longitudes/latitudes
            #adjust data arrays to account for potential changing number of stations 
            self.read_setup()

            print(1, time.time()-st)

            #need to re-read all observations/experiments?
            if read_all == True:
                #reset data in memory dictionary
                self.data_in_memory = {}
                self.plotting_params = {}
                self.metadata_inds_to_fill = np.arange(len(self.relevant_yearmonths))  
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

                    new_metadata_array = np.full((len(self.station_references), len(self.previous_relevant_yearmonths)), np.NaN, dtype=metadata_dtype) 
                    new_metadata_array[new_station_inds,:] = self.metadata_in_memory[old_station_inds,:]
                    self.metadata_in_memory = new_metadata_array

                    #iterate through all keys in data in memory dictionary
                    for data_label in list(self.data_in_memory.keys()):
                        #create new data array in shape of current station references array
                        if data_label == 'observations':
                            new_data_array = np.full((len(self.station_references),len(self.previous_time_array)), np.NaN, dtype=data_dtype)                        
                        else:
                            new_data_array = np.full((len(self.station_references),len(self.previous_time_array)), np.NaN, dtype=data_dtype[:1])
                        #put the old data into new array in the correct positions
                        new_data_array[new_station_inds,:] = self.data_in_memory[data_label][old_station_inds,:]
                        #overwrite data array with reshaped version
                        self.data_in_memory[data_label] = new_data_array

            #need to cut edges?
            if (cut_left == True) or (cut_right == True):

                #set default edge limits as current edges
                data_left_edge_ind = 0
                data_right_edge_ind = len(self.previous_time_array)

                metadata_left_edge_ind = 0
                metadata_right_edge_ind = len(self.previous_relevant_yearmonths)

                #need to cut on left data edge?                                                                                      
                if cut_left == True:
                    data_left_edge_ind = np.where(self.previous_time_array == self.time_array[0])[0][0]
                    new_relative_delta = relativedelta(datetime.datetime(int(self.relevant_yearmonths[0][:4]), int(self.relevant_yearmonths[0][4:6]), 1, 0, 0), datetime.datetime(int(self.previous_relevant_yearmonths[0][:4]), int(self.previous_relevant_yearmonths[0][4:6]), 1, 0, 0))
                    metadata_left_edge_ind = (new_relative_delta.years*12)+new_relative_delta.months               

                #need to cut on right data edge?
                if cut_right == True:
                    data_right_edge_ind = np.where(self.previous_time_array == self.time_array[-1])[0][0]+1
                    new_relative_delta = relativedelta(datetime.datetime(int(self.relevant_yearmonths[-1][:4]), int(self.relevant_yearmonths[-1][4:6]), 1, 0, 0), datetime.datetime(int(self.previous_relevant_yearmonths[-1][:4]), int(self.previous_relevant_yearmonths[-1][4:6]), 1, 0, 0))
                    metadata_right_edge_ind = (new_relative_delta.years*12)+new_relative_delta.months                 

                #iterate through all keys in data in memory dictionary and cut edges of the associated arrays appropriately 
                for data_label in list(self.data_in_memory.keys()):
                    self.data_in_memory[data_label] = self.data_in_memory[data_label][:,data_left_edge_ind:data_right_edge_ind]
                    self.metadata_in_memory = self.metadata_in_memory[:,metadata_left_edge_ind:metadata_right_edge_ind]

            #need to read on left edge?
            if read_left == True:
                #get n number of new elements on left edge
                n_new_left_data_inds = np.where(self.time_array == self.previous_time_array[0])[0][0]
                new_relative_delta = relativedelta(datetime.datetime(int(self.previous_relevant_yearmonths[0][:4]), int(self.previous_relevant_yearmonths[0][4:6]), 1, 0, 0), datetime.datetime(int(self.relevant_yearmonths[0][:4]), int(self.relevant_yearmonths[0][4:6]), 1, 0, 0))
                n_new_left_metadata_inds = (new_relative_delta.years*12)+new_relative_delta.months
                self.metadata_inds_to_fill = np.arange(0, n_new_left_metadata_inds)   
                self.metadata_in_memory = np.concatenate((np.full((len(self.station_references),n_new_left_metadata_inds),np.NaN, dtype=metadata_dtype), self.metadata_in_memory), axis=1)

                #iterate through all keys in data in memory dictionary and insert read data on left edge of the associated arrays 
                for data_label in list(self.data_in_memory.keys()):
                    #add space on left edge to insert new read data
                    if data_label == 'observations':
                        self.data_in_memory[data_label] = np.concatenate((np.full((len(self.station_references),n_new_left_data_inds),np.NaN, dtype=data_dtype), self.data_in_memory[data_label]), axis=1)
                    else:
                        self.data_in_memory[data_label] = np.concatenate((np.full((len(self.station_references),n_new_left_data_inds),np.NaN, dtype=data_dtype[:1]), self.data_in_memory[data_label]), axis=1)
                    self.read_data(data_label, self.active_start_date, self.previous_active_start_date)

            #need to read on right edge?
            if read_right == True:
                #get n number of new elements on right edge
                n_new_right_data_inds = (len(self.time_array) - 1) - np.where(self.time_array == self.previous_time_array[-1])[0][0]
                new_relative_delta = relativedelta(datetime.datetime(int(self.relevant_yearmonths[-1][:4]), int(self.relevant_yearmonths[-1][4:6]), 1, 0, 0), datetime.datetime(int(self.previous_relevant_yearmonths[-1][:4]), int(self.previous_relevant_yearmonths[-1][4:6]), 1, 0, 0))
                n_new_right_metadata_inds = (new_relative_delta.years*12)+new_relative_delta.months
                self.metadata_inds_to_fill = np.arange(-n_new_right_metadata_inds, 0)                 
                self.metadata_in_memory = np.concatenate((self.metadata_in_memory,np.full((len(self.station_references),n_new_right_metadata_inds),np.NaN, dtype=metadata_dtype)), axis=1)

                #iterate through all keys in data in memory dictionary and insert read data on right edge of the associated arrays
                for data_label in list(self.data_in_memory.keys()):
                    if data_label == 'observations':
                        self.data_in_memory[data_label] = np.concatenate((self.data_in_memory[data_label], np.full((len(self.station_references),n_new_right_data_inds),np.NaN, dtype=data_dtype)), axis=1)
                    else:
                        self.data_in_memory[data_label] = np.concatenate((self.data_in_memory[data_label], np.full((len(self.station_references),n_new_right_data_inds),np.NaN, dtype=data_dtype[:1])), axis=1)
                    self.read_data(data_label, self.previous_active_end_date, self.active_end_date)
        
            print(2, time.time()-st)

            #update menu object fields
            self.update_metadata_fields()
            self.update_representativity_fields()
            self.update_period_fields()

            print(3, time.time()-st)

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
           also, all station metadata is read.
        '''        

        #force garbage collection (to avoid memory issues)
        gc.collect()

        #set current time array, as previous time array
        self.previous_time_array = self.time_array    

        #set current station references, as previous station references 
        self.previous_station_references = self.station_references

        self.previous_relevant_yearmonths = self.relevant_yearmonths

        #get N time chunks between desired start date and end date to set time array
        if (self.active_resolution == 'hourly') or (self.active_resolution == 'hourly_instantaneous'):
            self.active_frequency_code = 'H'
        elif self.active_resolution == 'daily':
            self.active_frequency_code = 'D'
        elif self.active_resolution == 'monthly':
            self.active_frequency_code = 'MS'
        str_active_start_date = str(self.active_start_date)
        str_active_end_date = str(self.active_end_date)
        self.time_array = pd.date_range(start = datetime.datetime(int(str_active_start_date[:4]), int(str_active_start_date[4:6]), int(str_active_start_date[6:8])),end = datetime.datetime(int(str_active_end_date[:4]), int(str_active_end_date[4:6]), int(str_active_end_date[6:8])), freq = self.active_frequency_code)[:-1]

        #get all relevant observational files
        file_root = '%s/%s/%s/%s/%s/%s_'%(obs_root, self.active_network, GHOST_version, self.active_resolution, self.active_species, self.active_species)
        self.relevant_yearmonths = sorted([str(yyyymm) for yyyymm in self.available_observation_data[self.active_network][self.active_resolution][self.active_matrix][self.active_species]])
        relevant_files = sorted([file_root+yyyymm[:6]+'.nc' for yyyymm in self.relevant_yearmonths])       
        #self.monthly_inds = np.array([np.arange(len(self.time_array))[np.all([self.time_array >= datetime.datetime.strptime(start_yyyymm, '%Y%m%d'), self.time_array < datetime.datetime.strptime(self.relevant_yearmonths[month_ii+1], '%Y%m%d')], axis=0)] if month_ii != (len(self.relevant_yearmonths)-1) else np.arange(len(self.time_array))[self.time_array >= datetime.datetime.strptime(start_yyyymm, '%Y%m%d')] for month_ii, start_yyyymm in enumerate(self.relevant_yearmonths)])
        self.N_inds_per_month = np.array([np.count_nonzero(np.all([self.time_array >= datetime.datetime.strptime(start_yyyymm, '%Y%m%d'), self.time_array < datetime.datetime.strptime(self.relevant_yearmonths[month_ii+1], '%Y%m%d')], axis=0)) if month_ii != (len(self.relevant_yearmonths)-1) else np.count_nonzero(self.time_array >= datetime.datetime.strptime(start_yyyymm, '%Y%m%d')) for month_ii, start_yyyymm in enumerate(self.relevant_yearmonths)])

        #redefine some key variables globally (for access by parallel netCDF reading functions)
        global time_array, active_species, selected_qa, selected_flags
        time_array = self.time_array
        active_species = self.active_species   
        selected_qa = self.active_qa
        selected_flags = self.active_flags

        #add all metadata variables to read to station metadata dictionary with associated empty numpy arrays of the appropriate type
        global station_references
        station_references = []
        self.station_longitudes = []
        self.station_latitudes = []
        for relevant_file in relevant_files:
            nCDF_root = Dataset(relevant_file) 
            station_references=np.append(station_references, nCDF_root['station_reference'][:])
            self.station_longitudes=np.append(self.station_longitudes, nCDF_root['longitude'][:])
            self.station_latitudes=np.append(self.station_latitudes, nCDF_root['latitude'][:])
            nCDF_root.close()
        station_references, station_unique_indices = np.unique(station_references, return_index=True)
        self.station_references = station_references
        self.station_longitudes = self.station_longitudes[station_unique_indices]
        self.station_latitudes = self.station_latitudes[station_unique_indices]

        #update measurement units for species (take standard units from parameter dictionary)
        self.measurement_units = parameter_dictionary[self.active_species]['standard_units']

        #set data variables to read (dependent on active data resolution)
        global data_vars_to_read, data_dtype
        if (self.active_resolution == 'hourly') or (self.active_resolution == 'hourly_instantaneous'):
            data_vars_to_read = [self.active_species, 'hourly_native_representativity_percent' ,'daily_native_representativity_percent' ,'monthly_native_representativity_percent', 'annual_native_representativity_percent', 'hourly_native_max_gap_percent', 'daily_native_max_gap_percent', 'monthly_native_max_gap_percent', 'annual_native_max_gap_percent', 'day_night_code', 'weekday_weekend_code', 'season_code']
        elif self.active_resolution == 'daily':
            data_vars_to_read = [self.active_species, 'daily_native_representativity_percent' ,'monthly_native_representativity_percent', 'annual_native_representativity_percent', 'daily_native_max_gap_percent', 'monthly_native_max_gap_percent', 'annual_native_max_gap_percent', 'weekday_weekend_code', 'season_code']
        elif self.active_resolution == 'monthly':
            data_vars_to_read = [self.active_species, 'monthly_native_representativity_percent', 'annual_native_representativity_percent', 'monthly_native_max_gap_percent', 'annual_native_max_gap_percent', 'season_code']
  
        #set data dtype
        #data_dtype = [(key,np.float32) if key == active_species else (key,np.uint8) for key in data_vars_to_read]
        data_dtype = [(key,np.float32) for key in data_vars_to_read]
    
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
            file_root = '%s/%s/%s/%s/%s/%s_'%(obs_root, self.active_network, GHOST_version, self.active_resolution, self.active_species, self.active_species)
            relevant_file_start_dates = sorted(self.available_observation_data[self.active_network][self.active_resolution][self.active_matrix][self.active_species])  
        else:
            process_type = 'experiment'
            experiment_grid_split = data_label.split('-')
            active_experiment = experiment_grid_split[0]
            active_grid = experiment_grid_split[1]  
            file_root = '%s/%s/%s/%s/%s/%s/%s/%s_'%(exp_root, GHOST_version, active_experiment, active_grid, self.active_resolution, self.active_species, self.active_network, self.active_species)
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
            
            if process_type == 'observations':
                self.plotting_params['observations'] = {}
                self.data_in_memory[data_label] = np.full((len(self.station_references),len(self.time_array)), np.NaN, dtype=data_dtype)
                self.metadata_in_memory = np.full((len(self.station_references), len(self.relevant_yearmonths)), np.NaN, dtype=metadata_dtype)

            #if process_type is experiment, get experiment specific grid edges from first relevant file, and save to data in memory dictionary
            if process_type == 'experiment':
                self.data_in_memory[data_label] = np.full((len(self.station_references),len(self.time_array)), np.NaN, dtype=data_dtype[:1])
                self.plotting_params[data_label] = {}
                exp_nc_root = Dataset(relevant_files[0])
                self.plotting_params[data_label]['grid_edge_longitude'] = exp_nc_root['grid_edge_longitude'][:]
                self.plotting_params[data_label]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
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
                self.data_in_memory[data_label][self.read_instance.active_species][full_array_station_indices[np.newaxis,:], time_indices[:,np.newaxis]] = file_data
     
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
            for file_data_ii, file_data in enumerate(all_file_data):
                self.data_in_memory[data_label][file_data[2][:,np.newaxis], file_data[1][np.newaxis,:]] = file_data[0]
                if process_type == 'observations':
                    self.metadata_in_memory[file_data[2][:,np.newaxis], self.metadata_inds_to_fill[file_data_ii]] = file_data[3] 

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_metadata_fields(self):

        '''update the metadata menu object with metadata associated with newly read data
           for non-numeric metadata gets all the unique fields per metadata variable, and sets the available fields as such,
           and for numeric gets the minimum and maximum boundaries of each metadata variable, and sets initial boundaries as such
        '''

        #iterate through metadata variables
        for meta_var in metadata_vars_to_read:

            meta_var_field = self.metadata_in_memory[meta_var]

            #get metadata variable type/data type
            metadata_type = standard_metadata[meta_var]['metadata_type']
            metadata_data_type = standard_metadata[meta_var]['data_type']
            
            #remove NaNs from field
            meta_var_field_nan_removed = meta_var_field[~pd.isnull(meta_var_field)]
            
            #update pop-up metadata menu object with read metadata values
            #for non-numeric metadata gets all the unique fields per metadata variable, and sets the available fields as such
            if metadata_data_type == np.object:
                self.metadata_menu[metadata_type][meta_var]['checkboxes']['labels'] = np.unique(meta_var_field_nan_removed)  
                self.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_default'] = []
                self.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_default'] = []   
            #for numeric fields get the minimum and maximum boundaries of each metadata variable, and sets initial boundaries as such
            else:
                meta_var_index = self.metadata_menu[metadata_type]['rangeboxes']['labels'].index(meta_var)
                if len(meta_var_field_nan_removed) > 0:
                    min_val = str(np.min(meta_var_field_nan_removed))
                    max_val = str(np.max(meta_var_field_nan_removed))
                    self.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = min_val
                    self.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index] = min_val
                    self.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = max_val
                    self.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index] = max_val
                else:
                    self.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = 'nan'
                    self.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index] = 'nan'
                    self.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = 'nan'
                    self.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index] = 'nan'

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 
    
    def update_representativity_fields(self):

        '''update the data representativity menu -> 1D list of rangebox values
           dependent on the temporal resolution, some fields will appear or not
        '''

        #get previous set labels
        previous_labels = copy.deepcopy(self.representativity_menu['rangeboxes']['labels'])

        #get previously set rangebox values
        previous_lower = copy.deepcopy(self.representativity_menu['rangeboxes']['current_lower'])

        #hourly temporal resolution?
        if (self.active_resolution == 'hourly') or (self.active_resolution == 'hourly_instantaneous'):
            self.representativity_menu['rangeboxes']['labels'] = ['hourly_native_representativity_percent', 'hourly_native_max_gap_percent', 'daily_native_representativity_percent', 'daily_representativity_percent', 'daily_native_max_gap_percent', 'daily_max_gap_percent', 'monthly_native_representativity_percent', 'monthly_representativity_percent', 'monthly_native_max_gap_percent', 'monthly_max_gap_percent', 'all_representativity_percent', 'all_max_gap_percent']
        #daily temporal resolution?
        elif self.active_resolution == 'daily':
            self.representativity_menu['rangeboxes']['labels'] = ['daily_native_representativity_percent', 'daily_native_max_gap_percent', 'monthly_native_representativity_percent', 'monthly_representativity_percent', 'monthly_native_max_gap_percent', 'monthly_max_gap_percent', 'all_representativity_percent', 'all_max_gap_percent']
        #monthly temporal resolution?
        elif self.active_resolution == 'monthly':
            self.representativity_menu['rangeboxes']['labels'] = ['monthly_native_representativity_percent', 'monthly_native_max_gap_percent', 'all_representativity_percent', 'all_max_gap_percent']
        
        #initialise rangebox values --> for data representativity fields the default is 0%, for max gap fields % the default is 100%
        self.representativity_menu['rangeboxes']['current_lower'] = []
        for label_ii, label in enumerate(self.representativity_menu['rangeboxes']['labels']):
            if 'max_gap' in label:
                self.representativity_menu['rangeboxes']['current_lower'].append('100')
            else:
                self.representativity_menu['rangeboxes']['current_lower'].append('0') 

        #set new rangebox values (using previous values where possible)
        for label_ii, label in enumerate(self.representativity_menu['rangeboxes']['labels']):
            #label previously existed?
            if label in previous_labels:
                self.representativity_menu['rangeboxes']['current_lower'][label_ii] = previous_lower[previous_labels.index(label)]

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_period_fields(self):

        '''update the data period menu -> list of checkboxes
           dependent on the temporal resolution, some fields will appear or not
        '''

        #hourly temporal resolution?
        if (self.active_resolution == 'hourly') or (self.active_resolution == 'hourly_instantaneous'):
            self.period_menu['checkboxes']['labels'] = ['Daytime', 'Nighttime', 'Weekday', 'Weekend', 'Spring', 'Summer', 'Autumn', 'Winter']
        #daily temporal resolution?
        elif self.active_resolution == 'daily':
            self.period_menu['checkboxes']['labels'] = ['Weekday', 'Weekend', 'Spring', 'Summer', 'Autumn', 'Winter']
        #monthly temporal resolution?
        elif self.active_resolution == 'monthly':
            self.period_menu['checkboxes']['labels'] = ['Spring', 'Summer', 'Autumn', 'Winter']

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------#    

    def update_plotting_parameters(self):

        '''function that updates plotting parameters (colour and zorder) for each selected data array'''

        #assign a colour/zorder to all selected data arrays
        #define observations colour to be 'black'
        self.plotting_params['observations']['colour'] = 'black'
        #define zorder of observations to be 5
        self.plotting_params['observations']['zorder'] = 5
        
        #generate a list of RGB tuples for number of experiments there are
        sns.reset_orig()
        clrs = sns.color_palette('husl', n_colors=len(list(self.data_in_memory.keys()))-1)

        #iterate through sorted experiment names, assigning each experiment a new RGB colour tuple, and zorder
        experiment_ind = 1
        for experiment in sorted(list(self.data_in_memory.keys())):
            if experiment != 'observations':
                #define colour for experiment
                self.plotting_params[experiment]['colour'] = clrs[experiment_ind-1]
                #define zorder for experiment (obs zorder + experiment_ind)
                self.plotting_params[experiment]['zorder'] = self.plotting_params['observations']['zorder'] + experiment_ind
                #update count of experiments
                experiment_ind +=1

#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#

def read_netCDF_data(relevant_file):

    '''function that handles reading of observational/experiment netCDF data
       also handles filtering of observational data based on selected qa/data provider flags
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

    
    #file_data = np.full((len(current_file_station_indices), len(valid_file_time_indices)), np.NaN, data_dtype)
    #for data_var in data_vars_to_read:
    #    read_data = nCDF_root[data_var][current_file_station_indices,valid_file_time_indices]
        #get masked data 
        #data_mask = read_data.mask
        #set masked data as NaN
        #read_data[data_mask] = np.NaN    

    #for observations, set species data based on selected qa flags/standard data provider flags to retain or remove as NaN
    if process_type == 'observations':
        file_data = np.full((len(current_file_station_indices), len(valid_file_time_indices)), np.NaN, dtype=data_dtype)
        for data_var in data_vars_to_read:
            file_data[data_var][:,:] = nCDF_root[data_var][current_file_station_indices,valid_file_time_indices]

        #if some qa flags selected then screen
        if len(selected_qa) > 0:
            #screen out observations which are associated with any of the selected qa flags
            file_data[active_species][np.isin(nCDF_root['qa'][:,valid_file_time_indices,:], selected_qa).any(axis=2)] = np.NaN
        #if some data provider flags selected then screen
        if len(selected_flags) > 0:
            #screen out observations which are associated with any of the selected data provider flags
            file_data[active_species][np.isin(nCDF_root['flag'][:,valid_file_time_indices,:], selected_flags).any(axis=2)] = np.NaN
        
        #get file metadata
        file_metadata = np.full((len(file_station_references), 1), np.NaN, dtype=metadata_dtype)
        for meta_var in metadata_vars_to_read:
            file_metadata[meta_var][current_file_station_indices,0] = nCDF_root[meta_var][:]

    else:
        file_data = np.full((len(current_file_station_indices), len(valid_file_time_indices)), np.NaN, dtype=data_dtype[:1])
        file_data[data_vars_to_read[0]][:,:] = nCDF_root[data_vars_to_read[0]][current_file_station_indices,valid_file_time_indices]    
    
    #close netCDF
    nCDF_root.close()
    
    #return valid species data, time indices relative to active full time array, file station indices relative to all unique station references array 
    if process_type == 'observations':
        return file_data, full_array_time_indices, full_array_station_indices, file_metadata
    else:
        return file_data, full_array_time_indices, full_array_station_indices   

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
        
        #map_ax =              self.gridspec.new_subplotspec((0, 0),   rowspan=45, colspan=45)
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
        #hide all axes  
        #self.map_ax.axis('off')
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

            #add coastlines and gridlines
            land_polygon_resolutions = {'low':'110m', 'medium':'50m', 'high':'10m'}
            feature = cfeature.NaturalEarthFeature('physical', 'land', land_polygon_resolutions[map_coastline_resolution], facecolor='0.85')
            self.map_ax.add_feature(feature)
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
        if (self.read_instance.active_resolution == 'hourly') or (self.read_instance.active_resolution == 'hourly_instantaneous'):
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

        '''function which handles updates data filtering by selected lower/upper limit bounds, selected metadata and selected minimum data availability %'''

        #get set variables names representing percentage data availability (native and non-native)
        active_data_availablity_vars = self.read_instance.representativity_menu['rangeboxes']['labels']
        
        #get set lower/upper data bounds
        selected_lower_bound = self.read_instance.le_minimum_value.text()
        selected_upper_bound = self.read_instance.le_maximum_value.text()

        #check set data availability percent variable bounds and selected lower/upper bounds are numbers
        try:
            data_availability_lower_bounds = []
            data_availability_upper_bounds = []
            for var_ii, var in enumerate(active_data_availablity_vars):
                data_availability_lower_bounds.append(np.float32(self.read_instance.representativity_menu['rangeboxes']['current_lower'][var_ii]))
                #data_availability_upper_bounds.append(np.float32(self.read_instance.representativity_menu['rangeboxes']['current_upper'][var_ii]))
            selected_lower_bound = np.float32(selected_lower_bound)
            selected_upper_bound = np.float32(selected_upper_bound)
        #if any of the fields are not numbers, return from function
        except ValueError:
            return

        #reset data arrays to be un-filtered
        self.read_instance.data_in_memory_filtered = copy.deepcopy(self.read_instance.data_in_memory)

        #filter all observational data out of bounds of lower/upper limits
        inds_out_of_bounds = np.logical_or(self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species]<selected_lower_bound, self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species]>selected_upper_bound)
        self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][inds_out_of_bounds] = np.NaN  

        #filter/limit data for periods selected
        if len(self.read_instance.period_menu['checkboxes']['keep_selected']) > 0:
            if 'Daytime' in self.read_instance.period_menu['checkboxes']['keep_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['day_night_code'] != 0] = np.NaN
            if 'Nighttime' in self.read_instance.period_menu['checkboxes']['keep_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['day_night_code'] != 1] = np.NaN
            if 'Weekday' in self.read_instance.period_menu['checkboxes']['keep_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['weekday_weekend_code'] != 0] = np.NaN
            if 'Weekend' in self.read_instance.period_menu['checkboxes']['keep_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['weekday_weekend_code'] != 1] = np.NaN
            if 'Spring' in self.read_instance.period_menu['checkboxes']['keep_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['season_code'] != 0] = np.NaN
            if 'Summer' in self.read_instance.period_menu['checkboxes']['keep_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['season_code'] != 1] = np.NaN
            if 'Autumn' in self.read_instance.period_menu['checkboxes']['keep_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['season_code'] != 2] = np.NaN
            if 'Winter' in self.read_instance.period_menu['checkboxes']['keep_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['season_code'] != 3] = np.NaN

        if len(self.read_instance.period_menu['checkboxes']['remove_selected']) > 0:
            if 'Daytime' in self.read_instance.period_menu['checkboxes']['remove_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['day_night_code'] == 0] = np.NaN
            if 'Nighttime' in self.read_instance.period_menu['checkboxes']['remove_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['day_night_code'] == 1] = np.NaN
            if 'Weekday' in self.read_instance.period_menu['checkboxes']['remove_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['weekday_weekend_code'] == 0] = np.NaN
            if 'Weekend' in self.read_instance.period_menu['checkboxes']['remove_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['weekday_weekend_code'] == 1] = np.NaN
            if 'Spring' in self.read_instance.period_menu['checkboxes']['remove_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['season_code'] == 0] = np.NaN
            if 'Summer' in self.read_instance.period_menu['checkboxes']['remove_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['season_code'] == 1] = np.NaN
            if 'Autumn' in self.read_instance.period_menu['checkboxes']['remove_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['season_code'] == 2] = np.NaN
            if 'Winter' in self.read_instance.period_menu['checkboxes']['remove_selected']:
                self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations']['season_code'] == 3] = np.NaN

        #filter all obersvational data out of set bounds of native percentage data availability variables
        for var_ii, var in enumerate(active_data_availablity_vars):
            if 'native' in var:
                #max gap variable?
                if 'max_gap' in var:
                    self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations'][var] > data_availability_lower_bounds[var_ii]] = np.NaN
                #data representativity variable?
                else:
                    self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][self.read_instance.data_in_memory_filtered['observations'][var] < data_availability_lower_bounds[var_ii]] = np.NaN

        #filter all obersvational data out of set bounds of non-native percentage data availability variables
        for var_ii, var in enumerate(active_data_availablity_vars):
            if 'native' not in var:
                #get period associate with variable
                period = var.split('_')[0]
                period_inds = np.arange(self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species].shape[1])
                #daily variable?
                if period == 'daily':
                    period_inds_split = np.array_split(period_inds, [24*i for i in range(1, int(np.ceil(len(period_inds)/24)))])
                #monthly variable?
                elif period == 'monthly':
                    period_inds_split = np.array_split(period_inds, np.cumsum(self.read_instance.N_inds_per_month))
                #whole record variable?
                else:
                    period_inds_split = [period_inds]

                #iterate through indices associated with periodic chunks for current period
                for period_inds in period_inds_split:
                    if len(period_inds) > 0:
                        #max gap variable?
                        if 'max_gap' in var:
                            max_gap_percent = max_repeated_NaNs(self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][:, period_inds])
                            self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][max_gap_percent > data_availability_lower_bounds[var_ii]] = np.NaN
                        #data representativity variable?
                        else:
                            data_availability_percent = calculate_data_availability_fraction(self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][:, period_inds])
                            self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][data_availability_percent < data_availability_lower_bounds[var_ii]] = np.NaN

        #fiter all observational data by selected metadata
        #iterate through all metadata
        for meta_var in metadata_vars_to_read:
            metadata_type = standard_metadata[meta_var]['metadata_type']
            metadata_data_type = standard_metadata[meta_var]['data_type']

            #handle non-numeric metadata
            if metadata_data_type == np.object:
                #if any of the keep checkboxes have been selected, and some checkboxes have been changed from previous update, filter out data by fields that have not been selected
                current_keep = self.read_instance.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected']
                previous_keep = self.read_instance.metadata_menu[metadata_type][meta_var]['checkboxes']['previous_keep_selected']
                if (len(current_keep) > 0) & (set(previous_keep) != set(current_keep)):
                    invalid_keep = np.repeat(np.isin(self.read_instance.metadata_in_memory[meta_var][:,:], current_keep, invert=True), self.read_instance.N_inds_per_month, axis=1)
                    self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][invalid_keep] = np.NaN
                #if any of the remove checkboxes have been selected, and some checkboxes have been changed from previous update, filter out data by these selected fields
                current_remove = self.read_instance.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_selected']
                previous_remove = self.read_instance.metadata_menu[metadata_type][meta_var]['checkboxes']['previous_remove_selected']
                if (len(current_remove) > 0) & (set(previous_remove) != set(current_remove)):
                    invalid_remove = np.repeat(np.isin(self.read_instance.metadata_in_memory[meta_var][:,:], current_remove), self.read_instance.N_inds_per_month, axis=1)
                    self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][invalid_remove] = np.NaN 
            #handle numeric metadata
            else:
                meta_var_index = self.read_instance.metadata_menu[metadata_type]['rangeboxes']['labels'].index(meta_var)
                #filter out data with metadata < current lower value (if this is numeric)
                try:
                    current_lower = np.float32(self.read_instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index])
                    previous_lower = np.float32(self.read_instance.metadata_menu[metadata_type]['rangeboxes']['previous_lower'][meta_var_index])
                    #if current lower value is non-NaN, and different from previous value, then filter out data with metadata < current lower value
                    if (pd.isnull(current_lower) == False) & (current_lower != previous_lower):
                        invalid_below = np.repeat(self.read_instance.metadata_in_memory[meta_var][:,:] < current_lower, self.read_instance.N_inds_per_month, axis=1)  
                        self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][invalid_below] = np.NaN
                except:
                    pass
                #filter out data with metadata > current upper value (if this is numeric)
                try:
                    current_upper = np.float32(self.read_instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index])
                    previous_upper = np.float32(self.read_instance.metadata_menu[metadata_type]['rangeboxes']['previous_upper'][meta_var_index])
                    #if current upper value is non-NaN, and different from previous value, then filter out data with metadata > current upper value
                    if (pd.isnull(current_upper) == False) & (current_upper != previous_upper):
                        invalid_above = np.repeat(self.read_instance.metadata_in_memory[meta_var][:,:] > current_upper, self.read_instance.N_inds_per_month, axis=1)  
                        self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species][invalid_above] = np.NaN
                except:
                    pass

        #colocate data (if necessary) 
        self.colocate_data()

        #get indices of stations with > 1 valid measurements ,in all observational data arrays (colocated and non-colocated)
        #iterate through all data arrays
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()): 

            #check if data array is an observational data array
            if data_label.split('_')[0] == 'observations':
                            
                #get absolute data availability number per station in observational data array
                station_data_availability_number = calculate_data_availability_number(self.read_instance.data_in_memory_filtered[data_label][self.read_instance.active_species])                
                
                #get indices of stations with > 1 available measurements    
                #save valid station indices with data array
                self.read_instance.plotting_params[data_label]['valid_station_inds'] = np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]

        #write valid station indices calculated for observations across to associated experimental data arrays
        #iterate through all data arrays
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()): 
        
            #check if data array is not an observational data array
            if data_label.split('_')[0] != 'observations':
        
                #handle colocated experimental arrays
                if '_colocatedto_' in data_label:
                    exp_name = data_label.split('_colocatedto_')[0]
                    self.read_instance.plotting_params[data_label]['valid_station_inds'] = copy.deepcopy(self.read_instance.plotting_params['observations_colocatedto_%s'%(exp_name)]['valid_station_inds'])
                #handle non-located experimental arrays
                else:
                    self.read_instance.plotting_params[data_label]['valid_station_inds'] = copy.deepcopy(self.read_instance.plotting_params['observations']['valid_station_inds'])

        #after subsetting by pre-written associated observational valid stations, get indices of stations with > 1 valid measurements in all experiment data arrays (colocated and non-colocated)
        #iterate through all data arrays
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()): 
    
            #check if data array is not an observational data array
            if data_label.split('_')[0] != 'observations':
                #get indices of associated observational data array valid stations (pre-written to experiment data arrays)
                valid_station_inds = self.read_instance.plotting_params[data_label]['valid_station_inds']
                #get absolute data availability number per station in experiment data array, after subsetting valid observational stations (i.e. number of non-NaN measurements)
                station_data_availability_number = calculate_data_availability_number(self.read_instance.data_in_memory_filtered[data_label][self.read_instance.active_species][valid_station_inds,:])
                #get indices of stations with > 1 available measurements
                valid_station_inds = valid_station_inds[np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]]
                #overwrite previous written valid station indices (now at best a subset of those indices)
                self.read_instance.plotting_params[data_label]['valid_station_inds'] = valid_station_inds       
            
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
            nan_obs = np.isnan(self.read_instance.data_in_memory_filtered['observations'][self.read_instance.active_species])

            #create array for finding instances where have 0 valid values across all experiments
            #initialise as being all True, set as False on the occasion there is a valid value in an experiment 
            exps_all_nan = np.full(nan_obs.shape, True)

            #iterate through experiment data arrays in data in memory dictionary
            for data_label in list(self.read_instance.data_in_memory.keys()):
                if data_label != 'observations':
                    #get all instances experiment are NaN
                    nan_exp = np.isnan(self.read_instance.data_in_memory_filtered[data_label][self.read_instance.active_species])
                    #get all instances where either the observational array or experiment array are NaN at a given time
                    nan_instances = np.any([nan_obs,nan_exp],axis=0)
                    #create new observational array colocated to experiment
                    obs_data = copy.deepcopy(self.read_instance.data_in_memory_filtered['observations'])
                    obs_data[nan_instances] = np.NaN
                    self.read_instance.data_in_memory_filtered['observations_colocatedto_%s'%(data_label)] = obs_data
                    self.read_instance.plotting_params['observations_colocatedto_%s'%(data_label)] = {'colour':self.read_instance.plotting_params['observations']['colour'], 'zorder':self.read_instance.plotting_params['observations']['zorder']} 
                    #create new experiment array colocated to observations
                    exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[data_label])
                    exp_data[nan_instances] = np.NaN
                    self.read_instance.data_in_memory_filtered['%s_colocatedto_observations'%(data_label)] = exp_data 
                    self.read_instance.plotting_params['%s_colocatedto_observations'%(data_label)] = {'colour':self.read_instance.plotting_params[data_label]['colour'], 'zorder':self.read_instance.plotting_params[data_label]['zorder']}
                    #update exps_all_nan array, making False all instances where have valid experiment data
                    exps_all_nan = np.all([exps_all_nan, nan_exp], axis=0) 
                    
            #create observational data array colocated to be non-NaN whenever there is a valid data in at least 1 experiment
            exps_all_nan = np.any([nan_obs,exps_all_nan],axis=0)
            obs_data = copy.deepcopy(self.read_instance.data_in_memory_filtered['observations'])
            obs_data[exps_all_nan] = np.NaN
            self.read_instance.data_in_memory_filtered['observations_colocatedto_experiments'] = obs_data
            self.read_instance.plotting_params['observations_colocatedto_experiments'] = {'colour':self.read_instance.plotting_params['observations']['colour'], 'zorder':self.read_instance.plotting_params['observations']['zorder']}   

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
            self.map_points = self.map_ax.scatter(self.read_instance.station_longitudes[self.active_map_valid_station_inds],self.read_instance.station_latitudes[self.active_map_valid_station_inds], s=map_unselected_station_marker_size, c=self.z_statistic, vmin=self.z_vmin, vmax=self.z_vmax, cmap=self.z_colourmap, picker = 1, zorder=3, transform=self.datacrs, linewidth=0.0, alpha=None)     
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

        #reset alphas of all plotted stations (if have some stations on map)
        if len(self.active_map_valid_station_inds) > 0:
            self.rgba_tuples[:,-1] = 1.0
            marker_sizes = np.full(len(self.z_statistic), map_unselected_station_marker_size)
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
                marker_sizes[self.absolute_selected_station_inds] = map_selected_station_marker_size    
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
            #self.update_selected_station_metadata()

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
                grid_edge_outline_poly = Polygon(np.vstack((self.read_instance.plotting_params[experiment]['grid_edge_longitude'], self.read_instance.plotting_params[experiment]['grid_edge_latitude'])).T, edgecolor=self.read_instance.plotting_params[experiment]['colour'], linewidth=1, linestyle='--', fill=False, zorder=2, transform=self.datacrs)
                #plot grid edge polygon on map
                self.grid_edge_polygons.append(self.map_ax.add_patch(grid_edge_outline_poly))

    #--------------------------------------------------------------------------------#
    #--------------------------------------------------------------------------------# 

    def update_legend(self):

        '''function that updates legend'''

        #create legend elements
        #add observations element
        legend_elements = [Line2D([0], [0], marker='o', color='white', markerfacecolor=self.read_instance.plotting_params['observations']['colour'], markersize=legend_marker_size, label='observations')]
        #add element for each experiment
        for experiment_ind, experiment in enumerate(sorted(list(self.read_instance.data_in_memory.keys()))):
            if experiment != 'observations':
                #add experiment element
                legend_elements.append(Line2D([0], [0], marker='o', color='white', markerfacecolor=self.read_instance.plotting_params[experiment]['colour'], markersize=legend_marker_size, label=experiment))
        
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
                data_array = self.read_instance.data_in_memory_filtered[data_label][self.read_instance.active_species][self.relative_selected_station_inds,:]
            #experiment arrays
            else:          
                #get intersect between selected station indices and valid available indices for experiment data array
                valid_selected_station_indices = np.intersect1d(self.relative_selected_station_inds, self.read_instance.plotting_params[data_label]['valid_station_inds'])
                #get data for valid selected stations
                data_array = self.read_instance.data_in_memory_filtered[data_label][self.read_instance.active_species][valid_selected_station_indices,:]

            #if data array has no valid data for selected stations, do not create a pandas dataframe
            #data array has valid data?
            if data_array.size:
                #add nested dictionary for data array name to selection station data dictionary
                self.selected_station_data[data_label] = {}
                #take cross station median of selected data for data array, and place it in a pandas dataframe -->  add to selected station data dictionary
                self.selected_station_data[data_label]['pandas_df'] = pd.DataFrame(np.nanmedian(data_array, axis=0), index=self.read_instance.time_array, columns=[self.read_instance.active_species])
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
                grouped_data = [g[self.read_instance.active_species].dropna() for n, g in self.selected_station_data[data_label]['pandas_df'].groupby(getattr(self.selected_station_data[data_label]['pandas_df'].index, temporal_aggregation_resolution))]
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
            self.data_array_ts = self.ts_ax.plot(self.selected_station_data[data_label]['pandas_df'].dropna(), color=self.read_instance.plotting_params[data_label]['colour'], marker = 'o', markeredgecolor = None, mew = 0, markersize = time_series_marker_size, linestyle = 'None', zorder=self.read_instance.plotting_params[data_label]['zorder'])
              
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
        if (self.read_instance.active_resolution == 'hourly') or (self.read_instance.active_resolution == 'hourly_instantaneous'):
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
                    patch.set_facecolor(self.read_instance.plotting_params[data_label]['colour'])
                    patch.set_zorder(self.read_instance.plotting_params[data_label]['zorder'])
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
                median_zorder = (self.read_instance.plotting_params['observations']['zorder']+len(list(aggregation_dict[temporal_aggregation_resolution]['plots'].keys())) - 1) + self.read_instance.plotting_params[data_label]['zorder']

                #get xticks (all valid aggregated time indexes) and medians to plot
                xticks = aggregation_dict[temporal_aggregation_resolution]['xticks']
                medians = self.selected_station_data[data_label][temporal_aggregation_resolution]['p50']
                #split arrays if there are any temporal gaps to avoid line drawn being interpolated across missing values
                inds_to_split = np.where(np.diff(xticks) > 1)[0]
                if len(inds_to_split) == 0:
                    aggregation_dict[temporal_aggregation_resolution]['ax'].plot(xticks, medians, marker='o', color=self.read_instance.plotting_params[data_label]['colour'], markersize=temporally_aggregated_marker_size, linewidth=0.5, zorder=median_zorder)
                else:
                    inds_to_split += 1
                    start_ind = 0
                    for end_ind in inds_to_split:
                        aggregation_dict[temporal_aggregation_resolution]['ax'].plot(xticks[start_ind:end_ind], medians[start_ind:end_ind], marker='o', color=self.read_instance.plotting_params[data_label]['colour'], markersize=temporally_aggregated_marker_size, linewidth=0.5, zorder=median_zorder)
                        start_ind = end_ind   
                    aggregation_dict[temporal_aggregation_resolution]['ax'].plot(xticks[start_ind:], medians[start_ind:], marker='o', color=self.read_instance.plotting_params[data_label]['colour'], markersize=temporally_aggregated_marker_size, linewidth=0.5, zorder=median_zorder) 

        #------------------------------------------------------------------------------------------------#
        #plot title (with units)
            
        #if selected data resolution is 'hourly', plot the title on off the hourly aggregation axis 
        if (self.read_instance.active_resolution == 'hourly') or (self.read_instance.active_resolution == 'hourly_instantaneous'):
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
        if (self.read_instance.active_resolution == 'hourly') or (self.read_instance.active_resolution == 'hourly_instantaneous'):
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
        if (self.read_instance.active_resolution == 'hourly') or (self.read_instance.active_resolution == 'hourly_instantaneous'):
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
                    aggregation_dict[temporal_aggregation_resolution]['plots'][data_label] = aggregation_dict[temporal_aggregation_resolution]['ax'].plot(aggregation_dict[temporal_aggregation_resolution]['xticks'], self.selected_station_data[data_label][temporal_aggregation_resolution][selected_experiment_bias_stat], color=self.read_instance.plotting_params[data_label]['colour'], marker = 'o', zorder=self.read_instance.plotting_params[data_label]['zorder'], markersize=temporally_aggregated_experiment_bias_marker_size, linewidth=0.5)            

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
                                    'sample_preparation_types':'Sample Preparation',
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
                                     'primary_sampling_type', 'sample_preparation_types']

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
        
        #get relevant observational array (dependent on colocation)
        if self.colocate_active == False:
            obs_array = self.read_instance.plotting_params['observations']['valid_station_inds']
        else:
            obs_array = self.read_instance.plotting_params['observations_colocatedto_experiments']['valid_station_inds']

        #before doing anything check if have any valid station data for observations, if not update active map valid station indices to be empty list and return
        if len(obs_array) == 0:
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
            self.active_map_valid_station_inds = self.read_instance.plotting_params[z1_array_to_read]['valid_station_inds']
        else:
            #if have z2 array, get intersection of z1 and z2 valid station indices
            self.active_map_valid_station_inds = np.intersect1d(self.read_instance.plotting_params[z1_array_to_read]['valid_station_inds'], self.read_instance.plotting_params[z2_array_to_read]['valid_station_inds'])

        #update absolute selected plotted station indices with respect to new active map valid station indices
        self.absolute_selected_station_inds = np.array([np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind in self.relative_selected_station_inds if selected_ind in self.active_map_valid_station_inds], dtype=np.int)

        #-------------------------------------------------#

        #read z1 data
        z1_array_data = self.read_instance.data_in_memory_filtered[z1_array_to_read][self.read_instance.active_species][self.active_map_valid_station_inds,:]
        #drop NaNs and reshape to object list of station data arrays (if not checking data %)
        if z_statistic_name != 'Data %':
            z1_array_data = drop_NaNs(z1_array_data)
        else:
            z1_array_data.tolist()

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
            z2_array_data = self.read_instance.data_in_memory_filtered[z2_array_to_read][self.read_instance.active_species][self.active_map_valid_station_inds,:]
            #drop NaNs and reshape to object list of station data arrays (if not checking data %)
            if z_statistic_name != 'Data %':
                z2_array_data = drop_NaNs(z2_array_data)
            else:
                z2_array_data = z2_array_data.tolist()          
            #is the difference statistic basic (i.e. mean)?
            if z_statistic_type == 'basic':

                #load default selected z statistic arguments and make separate arguments dictionaries for z1/z2 calculations (as doing 2 separate calculations for z1/z2 and subtracting) 
                function_arguments_z1 = stats_dict['arguments']
                function_arguments_z2 = copy.deepcopy(function_arguments_z1)

                #iterate through stations calculating statistic
                for z_ii in range(len(self.z_statistic)):

                    #set station z1/z2 arrays as arguments in argument dictionaries
                    function_arguments_z1[self.read_instance.active_species] = z1_array_data[z_ii] 
                    function_arguments_z2[self.read_instance.active_species] = z2_array_data[z_ii]

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
                            intersect_lists.append(self.read_instance.plotting_params[data_label]['valid_station_inds'])
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

        '''function that handles single station selection/de-selection upon mouse click'''
        
        #update variable to inform lasso handler that map as already been updated (to not redraw map)
        #the on_click function is only called when a station index has been selected  
        #the variable will be reset by lasso handler (which is always called after on_click)
        self.map_already_updated = True

        #if the mouse click is not recognised as left(1)/right(3) then return from function
        if event.mouseevent.button not in [1,3]:
            return

        #make copy of current full array relative/absolute selected stations indices, before selecting new ones
        self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)
        self.previous_absolute_selected_station_inds = copy.deepcopy(self.absolute_selected_station_inds)

        #get new absolute selected index of station on map 
        absolute_selected_station_inds = np.array([event.ind[0]], dtype=np.int)

        #get new relative selected index with respect to all available stations
        relative_selected_station_inds = self.map_selected_station_inds_to_all_available_inds(absolute_selected_station_inds)

        #if left click (code of 1) --> select station 
        if event.mouseevent.button == 1:
            self.absolute_selected_station_inds = absolute_selected_station_inds
            self.relative_selected_station_inds = relative_selected_station_inds
        #if right click (code of 3) --> unselect station (if currently selected), select station (if currently unselected)  
        elif event.mouseevent.button == 3:
            #if len(self.previous_relative_selected_station_inds) > 0:
            relative_index = np.where(self.previous_relative_selected_station_inds == relative_selected_station_inds)[0]
            if len(relative_index) > 0:
                self.relative_selected_station_inds = np.delete(self.relative_selected_station_inds, relative_index)
            else:
                self.relative_selected_station_inds = np.append(self.relative_selected_station_inds, relative_selected_station_inds)

            absolute_index = np.where(self.previous_absolute_selected_station_inds == absolute_selected_station_inds)[0]
            if len(absolute_index) > 0:
                self.absolute_selected_station_inds = np.delete(self.absolute_selected_station_inds, absolute_index)
            else:
                self.absolute_selected_station_inds = np.append(self.absolute_selected_station_inds, absolute_selected_station_inds)    

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
        #get absolute selected indices of stations on map (the station coordinates contained within lasso)
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
    data = data.tolist()
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
    return (100./data.shape[-1]) * (np.count_nonzero(~np.isnan(data),axis=-1))

def calculate_data_availability_number(data):
    '''calculate data availability absolute number (i.e. number of total data measurements not equal to NaN)'''
    return np.count_nonzero(~np.isnan(data),axis=-1)

def max_repeated_NaNs(data):
    '''get maximum number of consecutive NaNs in array'''
    max_gap_pc = []

    for station_data in data:
        mask = np.concatenate(([False],np.isnan(station_data),[False]))
        if ~mask.any():
            max_gap_pc.append(0)
        else:
            idx = np.nonzero(mask[1:] != mask[:-1])[0]
            max_gap_pc.append((idx[1::2] - idx[::2]).max())
        
    return np.array(max_gap_pc)*(100./data.shape[1])

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
    '''
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

#initialise Providentia main window
qApp = QtWidgets.QApplication(sys.argv)
qApp.setStyle("Fusion")
Providentia_dash = Providentia_main_window('parallel')
sys.exit(qApp.exec_())