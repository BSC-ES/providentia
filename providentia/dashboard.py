""" Module which provides main window """
from .configuration import ProvConfiguration
from .configuration import split_options
from .init_standards import InitStandards
from .canvas import MPLCanvas
from .toolbar import NavigationToolbar
from .toolbar import save_data, conf_dialogs
from .dashboard_aux import ComboBox
from .dashboard_aux import QVLine
from .dashboard_aux import PopUpWindow
from .dashboard_aux import formatting_dict
from .dashboard_aux import set_formatting
from .read import DataReader
from providentia import aux

import copy
import datetime
import os
import os.path
import json
import sys
from functools import partial
from collections import OrderedDict

from PyQt5 import QtCore, QtWidgets, QtGui
import numpy as np
import pandas as pd

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))

class ProvidentiaMainWindow(QtWidgets.QWidget, ProvConfiguration, InitStandards):
    """Define class that generates Providentia dashboard"""

    # create signals that are fired upon resizing/moving of main Providentia window
    resized = QtCore.pyqtSignal()
    move = QtCore.pyqtSignal()

    def __init__(self, **kwargs):
        super(ProvidentiaMainWindow, self).__init__()
        ProvConfiguration.__init__(self, **kwargs)

        # store options to be restored at the end
        dconf_path = (os.path.join(CURRENT_PATH, 'conf/default.conf'))

        # update from config file (if available)
        if ('config' in kwargs) and (os.path.exists(kwargs['config'])):
            if 'section' in kwargs:
                # config and section defined 
                self.from_conf = True
                self.from_section = True
                aux.load_conf(self, fpath=kwargs['config'])
                if kwargs['section'] in self.all_sections:
                    vars(self).update({(k, self.parse_parameter(k, val)) for k, val in self.sub_opts[kwargs['section']].items()})
                else:
                    error = 'Error: The section specified in the command line does not exist.'
                    tip = 'Tip: For subsections, add the name of the parent section followed by a vertical bar (|) before the subsection name (e.g. SECTIONA|Spain).'
                    sys.exit(error + '\n' + tip)

            elif 'section' not in kwargs:
                # config defined, section undefined
                self.from_conf = True
                self.from_section = False
                aux.load_conf(self, fpath=kwargs['config'])    
                all_sections = self.sub_opts.keys()
                selected_section, okpressed = QtWidgets.QInputDialog.getItem(self, 'Sections',
                                                                             'Select section to load',  
                                                                             all_sections, 0, False)
                if okpressed:
                    vars(self).update({(k, self.parse_parameter(k, val)) for k, val in self.sub_opts[selected_section].items()})
        elif ('config' in kwargs) and (not os.path.exists(kwargs['config'])):     
            error = 'Error: The configuration path specified in the command line does not exist.'
            sys.exit(error)
        else:
            if os.path.isfile(dconf_path):
                # config undefined
                self.from_conf = False
                self.from_section = False
                aux.load_conf(self, fpath=dconf_path)
                vars(self).update({(k, self.parse_parameter(k, val)) for k, val in self.sub_opts['default'].items()})

        # update from command line
        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in kwargs.items()})

        # load characteristics per plot type
        # check for self defined plot characteristics file
        if self.plot_characteristics_filename == '':
            self.plot_characteristics_filename = os.path.join(
                CURRENT_PATH, 'conf/plot_characteristics_dashboard.json')
        self.plot_characteristics_templates = json.load(open(self.plot_characteristics_filename))

        # arguments are only local
        self.main_window_geometry = None
        
        # init GHOST standards
        InitStandards.__init__(self, ghost_root=self.ghost_root,
                               ghost_version=self.ghost_version)

        # create dictionary of all available observational GHOST data
        self.all_observation_data = aux.get_ghost_observational_tree(self)

        # load dictionary with non-GHOST esarchive files to read
        nonghost_observation_data_json = json.load(open(os.path.join(CURRENT_PATH, 'conf/nonghost_files.json')))
        # and merge to existing GHOST observational data dict if we have the path
        if self.nonghost_root is not None:
            nonghost_observation_data = aux.get_nonghost_observational_tree(self, nonghost_observation_data_json)
            self.all_observation_data = {**self.all_observation_data, **nonghost_observation_data}

        # initialise DataReader
        self.datareader = DataReader(self)

        #initialise UI
        self.init_ui(**kwargs)

        # setup callback events upon resizing/moving of Providentia window
        self.resized.connect(self.get_geometry)
        self.move.connect(self.get_geometry)

    def resizeEvent(self, event):
        """Function to overwrite default PyQt5 resizeEvent function --> for calling get_geometry"""
        self.resized.emit()
        return super(ProvidentiaMainWindow, self).resizeEvent(event)

    def moveEvent(self, event):
        """Function to overwrite default PyQt5 moveEvent function --> for calling get_geometry"""
        self.move.emit()
        return super(ProvidentiaMainWindow, self).moveEvent(event)

    def get_geometry(self):
        """Get current geometry of main Providentia window and buttons"""

        # get geometry of main window
        self.main_window_geometry = copy.deepcopy(self.geometry())

        # get and update geometry of settings menus
        self.update_buttons_geometry()

    def update_buttons_geometry(self):
        """Update current geometry of buttons"""

        for i, position in enumerate([self.position_1, self.position_2, self.position_3, 
                                      self.position_4, self.position_5]):
            for button, element in zip(self.mpl_canvas.buttons, self.mpl_canvas.elements):

                # get plot type
                plot_type = button.objectName().split('_menu')[0]
                if plot_type == 'periodic_violin':
                    plot_type = 'periodic-violin'
                
                if position == plot_type:
                    # calculate proportional geometry of buttons respect main window
                    x = (self.mpl_canvas.plot_characteristics_templates['general']['settings_menu']['position_'+str(i+1)]['x'] * 
                        self.main_window_geometry.width()) / 1848
                    y = (self.mpl_canvas.plot_characteristics_templates['general']['settings_menu']['position_'+str(i+1)]['y'] * 
                        self.main_window_geometry.height()) / 1016
                    
                    # get geometries
                    old_button_geometry = QtCore.QRect(button.x(), button.y(), 20, 20)
                    new_button_geometry = QtCore.QRect(x, y, 20, 20)
                    
                    # apply new geometry to button
                    button.setGeometry(new_button_geometry)

                    # apply new geometry to container elements
                    for sub_element in element:
                        if isinstance(sub_element, dict):
                            for sub_sub_element in sub_element.values():
                                sub_sub_element.setGeometry(sub_sub_element.x() - old_button_geometry.x() + 
                                                            new_button_geometry.x(), 
                                                            sub_sub_element.y() - old_button_geometry.y() + 
                                                            new_button_geometry.y(),
                                                            sub_sub_element.width(), sub_sub_element.height())
                        else:
                            sub_element.setGeometry(sub_element.x() - old_button_geometry.x() + new_button_geometry.x(), 
                                                    sub_element.y() - old_button_geometry.y() + new_button_geometry.y(),
                                                    sub_element.width(), sub_element.height())
                        
                    continue

    def init_ui(self, **kwargs):
        """Initialise user interface"""

        print("Starting Providentia online...")

        # set window title
        self.window_title = "Providentia"
        self.setWindowTitle(self.window_title)

        # create parent layout to pull together a configuration bar,
        # a MPL navigation toolbar, and a MPL canvas of plots
        parent_layout = QtWidgets.QVBoxLayout()

        # define spacing/margin variables
        parent_layout.setSpacing(0)
        parent_layout.setContentsMargins(0, 0, 0, 0)

        # define stylesheet for tooltips
        self.setStyleSheet("QToolTip { font: %spt %s}" % (formatting_dict['tooltip']['font'].pointSizeF(),
                                                          formatting_dict['tooltip']['font'].family()))

        # setup configuration bar with combo boxes, input boxes and buttons
        # use a gridded layout to place objects
        config_bar = QtWidgets.QGridLayout()

        # define spacing/margin variables
        config_bar.setHorizontalSpacing(3)
        config_bar.setVerticalSpacing(1)
        config_bar.setContentsMargins(5, 0, 0, 0)
        config_bar.setAlignment(QtCore.Qt.AlignLeft)

        # add one more horizontal layout
        hbox = QtWidgets.QHBoxLayout()

        # define all configuration box objects (labels, comboboxes etc.)
        self.lb_data_selection = set_formatting(QtWidgets.QLabel(self, text="Data Selection"),
                                                formatting_dict['title_menu'])
        self.lb_data_selection.setToolTip('Setup configuration of data to read into memory')
        self.bu_read = set_formatting(QtWidgets.QPushButton('READ', self), formatting_dict['button_menu'])
        self.bu_read.setFixedWidth(80)
        self.bu_read.setStyleSheet("color: green;")
        self.bu_read.setToolTip('Read selected configuration of data into memory')
        self.ch_colocate = set_formatting(QtWidgets.QCheckBox("Colocate"), formatting_dict['checkbox_menu'])
        self.ch_colocate.setToolTip('Temporally colocate observational/experiment data')
        self.cb_network = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_network.setFixedWidth(100)
        self.cb_network.setToolTip('Select providing observational data network. '
                                   'Names starting with * indicate non-GHOST datasets')
        self.cb_resolution = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_resolution.setFixedWidth(100)
        self.cb_resolution.setToolTip('Select temporal resolution of data')
        self.cb_matrix = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_matrix.setFixedWidth(100)
        self.cb_matrix.setToolTip('Select data matrix')
        self.cb_species = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_species.setStyleSheet("QComboBox { combobox-popup: 0; }")
        self.cb_species.setFixedWidth(100)
        self.cb_species.setToolTip('Select species')
        self.cb_species.view().setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.le_start_date = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_start_date.setFixedWidth(100)
        self.le_start_date.setToolTip('Set data start date: YYYYMMDD')
        self.le_end_date = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_end_date.setFixedWidth(100)
        self.le_end_date.setToolTip('Set data end date: YYYYMMDD')
        self.bu_QA = set_formatting(QtWidgets.QPushButton('QA', self), formatting_dict['button_menu'])
        self.bu_QA.setFixedWidth(80)
        self.bu_QA.setToolTip('Select standardised quality assurance flags to filter by')
        self.bu_flags = set_formatting(QtWidgets.QPushButton('FLAGS', self), formatting_dict['button_menu'])
        self.bu_flags.setFixedWidth(80)
        self.bu_flags.setToolTip('Select standardised data reporter provided flags to filter by')
        self.bu_experiments = set_formatting(QtWidgets.QPushButton('EXPS', self), formatting_dict['button_menu'])
        self.bu_experiments.setFixedWidth(80)
        self.bu_experiments.setToolTip('Select experiment/s data to read')
        self.vertical_splitter_1 = QVLine()
        self.vertical_splitter_1.setMaximumWidth(20)
        self.lb_data_filter = set_formatting(QtWidgets.QLabel(self, text="Filters"), formatting_dict['title_menu'])
        self.lb_data_filter.setFixedWidth(65)
        self.lb_data_filter.setToolTip('Select criteria to filter data by')
        self.bu_rep = set_formatting(QtWidgets.QPushButton('% REP', self), formatting_dict['button_menu'])
        self.bu_rep.setFixedWidth(80)
        self.bu_rep.setToolTip('Select % desired representativity in data across '
                               'whole record and for specific temporal periods')
        self.bu_meta = set_formatting(QtWidgets.QPushButton('META', self), formatting_dict['button_menu'])
        self.bu_meta.setFixedWidth(80)
        self.bu_meta.setToolTip('Select metadata to filter by')
        self.bu_reset = set_formatting(QtWidgets.QPushButton('RESET', self), formatting_dict['button_menu'])
        self.bu_reset.setFixedWidth(80)
        self.bu_reset.setToolTip('Reset filter fields to initial values')
        self.bu_reset.setStyleSheet("color: red;")
        self.bu_period = set_formatting(QtWidgets.QPushButton('PERIOD', self), formatting_dict['button_menu'])
        self.bu_period.setFixedWidth(80)
        self.bu_period.setToolTip('Select data in specific periods')
        self.bu_screen = set_formatting(QtWidgets.QPushButton('FILTER', self), formatting_dict['button_menu'])
        self.bu_screen.setFixedWidth(80)
        self.bu_screen.setStyleSheet("color: blue;")
        self.bu_screen.setToolTip('Filter data')
        self.lb_data_bounds = set_formatting(QtWidgets.QLabel(self, text="Bounds"), formatting_dict['label_menu'])
        self.lb_data_bounds.setFixedWidth(80)
        self.lb_data_bounds.setToolTip('Set lower/upper bounds of data')
        self.le_minimum_value = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_minimum_value.setFixedWidth(80)
        self.le_minimum_value.setToolTip('Set lower bound of data')
        self.le_maximum_value = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_maximum_value.setFixedWidth(80)
        self.le_maximum_value.setToolTip('Set upper bound of data')
        self.vertical_splitter_2 = QVLine()
        self.vertical_splitter_2.setMaximumWidth(20)
        self.lb_z = set_formatting(QtWidgets.QLabel(self, text="Map Stat"), formatting_dict['title_menu'])
        self.lb_z.setToolTip('Set plotted map statistic')
        self.cb_z_stat = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_z_stat.setFixedWidth(136)
        self.cb_z_stat.setToolTip('Select map statistic')
        self.cb_z1 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_z1.setFixedWidth(136)
        self.cb_z1.setToolTip('Select map dataset 1')
        self.cb_z2 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_z2.setFixedWidth(136)
        self.cb_z2.setToolTip('Select map dataset 2')
        self.vertical_splitter_3 = QVLine()
        self.vertical_splitter_3.setMaximumWidth(20)
        self.lb_experiment_bias = set_formatting(QtWidgets.QLabel(self, text="Exp Bias"),
                                                 formatting_dict['title_menu'])
        self.lb_experiment_bias.setToolTip('Set experiment bias statistic')
        self.cb_experiment_bias_type = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_experiment_bias_type.setFixedWidth(136)
        self.cb_experiment_bias_type.setToolTip('Select experiment bias type')
        self.cb_experiment_bias_stat = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_experiment_bias_stat.setFixedWidth(136)
        self.cb_experiment_bias_stat.setToolTip('Select experiment bias statistic')
        self.vertical_splitter_4 = QVLine()
        self.vertical_splitter_4.setMaximumWidth(20)
        self.lb_station_selection = set_formatting(QtWidgets.QLabel(self, text="Site Select"),
                                                   formatting_dict['title_menu'])
        self.lb_station_selection.setToolTip('Select stations')
        self.ch_select_all = set_formatting(QtWidgets.QCheckBox("All"), formatting_dict['checkbox_menu'])
        self.ch_select_all.setToolTip('Select all stations')
        self.ch_select_all.setFixedWidth(80)
        self.ch_intersect = set_formatting(QtWidgets.QCheckBox("Intersect"), formatting_dict['checkbox_menu'])
        self.ch_intersect.setToolTip('Select stations that intersect with all loaded model domains')
        self.ch_intersect.setFixedWidth(80)
        self.ch_extent = set_formatting(QtWidgets.QCheckBox("Extent"), formatting_dict['checkbox_menu'])
        self.ch_extent.setToolTip('Select stations that are within the map extent')
        self.ch_extent.setFixedWidth(80)
        self.vertical_splitter_5 = QVLine()
        self.vertical_splitter_5.setMaximumWidth(20)
        self.lb_layout_selection = set_formatting(QtWidgets.QLabel(self, text="Layout"),
                                                  formatting_dict['title_menu'])
        self.cb_position_1 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_position_1.setFixedWidth(100)
        self.cb_position_1.setToolTip('Select plot type in top left position')
        self.cb_position_2 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_position_2.setFixedWidth(100)
        self.cb_position_2.setToolTip('Select plot type in top right position')
        self.cb_position_3 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_position_3.setFixedWidth(100)
        self.cb_position_3.setToolTip('Select plot type in bottom left position')
        self.cb_position_4 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_position_4.setFixedWidth(100)
        self.cb_position_4.setToolTip('Select plot type in bottom centre position')
        self.cb_position_5 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_position_5.setFixedWidth(100)
        self.cb_position_5.setToolTip('Select plot type in bottom right position')

        # position objects on gridded configuration bar
        config_bar.addWidget(self.lb_data_selection, 0, 0, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_network, 1, 0)
        config_bar.addWidget(self.cb_resolution, 2, 0)
        config_bar.addWidget(self.cb_matrix, 1, 1)
        config_bar.addWidget(self.cb_species, 1, 2)
        config_bar.addWidget(self.le_start_date, 2, 1)
        config_bar.addWidget(self.le_end_date, 2, 2)
        config_bar.addWidget(self.bu_QA, 1, 3)
        config_bar.addWidget(self.bu_flags, 2, 3)
        config_bar.addWidget(self.bu_experiments, 1, 4)
        config_bar.addWidget(self.ch_colocate, 2, 4)
        config_bar.addWidget(self.bu_read, 3, 4, QtCore.Qt.AlignRight)
        config_bar.addWidget(self.vertical_splitter_1, 0, 5, 4, 1)
        config_bar.addWidget(self.lb_data_filter, 0, 6, 1, 2, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.lb_data_bounds, 1, 6, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.le_minimum_value, 1, 7)
        config_bar.addWidget(self.le_maximum_value, 1, 8)
        config_bar.addWidget(self.bu_rep, 2, 6)
        config_bar.addWidget(self.bu_period, 2, 7)
        config_bar.addWidget(self.bu_meta, 2, 8)
        config_bar.addWidget(self.bu_reset, 3, 7)
        config_bar.addWidget(self.bu_screen, 3, 8)
        config_bar.addWidget(self.vertical_splitter_2, 0, 9, 4, 1)
        config_bar.addWidget(self.lb_z, 0, 10, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_z_stat, 1, 10)
        config_bar.addWidget(self.cb_z1, 2, 10)
        config_bar.addWidget(self.cb_z2, 3, 10)
        config_bar.addWidget(self.vertical_splitter_3, 0, 11, 4, 1)
        config_bar.addWidget(self.lb_experiment_bias, 0, 12, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_experiment_bias_type, 1, 12)
        config_bar.addWidget(self.cb_experiment_bias_stat, 2, 12)
        config_bar.addWidget(self.vertical_splitter_4, 0, 13, 4, 1)
        config_bar.addWidget(self.lb_station_selection, 0, 14, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_select_all, 1, 14)
        config_bar.addWidget(self.ch_intersect, 2, 14)
        config_bar.addWidget(self.ch_extent, 3, 14)
        config_bar.addWidget(self.vertical_splitter_5, 0, 15, 4, 1)
        config_bar.addWidget(self.lb_layout_selection, 0, 16, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_position_1, 1, 16)
        config_bar.addWidget(self.cb_position_2, 1, 17)
        config_bar.addWidget(self.cb_position_3, 2, 16)
        config_bar.addWidget(self.cb_position_4, 2, 17)
        config_bar.addWidget(self.cb_position_5, 2, 18)

        # enable dynamic updating of configuration bar fields which filter data files
        self.cb_network.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_resolution.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_matrix.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_species.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.le_start_date.textChanged.connect(self.handle_config_bar_params_change)
        self.le_end_date.textChanged.connect(self.handle_config_bar_params_change)

        # setup pop-up window menu tree for flags, qa, experiments, 
        # % data representativity, data periods and metadata
        aux.init_flags(self)
        aux.init_qa(self)
        aux.init_experiments(self)
        aux.init_representativity(self)
        aux.init_period(self)
        self.metadata_vars_to_read = []
        aux.init_metadata(self)

        # enable pop up configuration windows
        self.bu_flags.clicked.connect(partial(self.generate_pop_up_window, self.flag_menu))
        self.bu_QA.clicked.connect(partial(self.generate_pop_up_window, self.qa_menu))
        self.bu_experiments.clicked.connect(partial(self.generate_pop_up_window, self.experiments_menu))
        self.bu_meta.clicked.connect(partial(self.generate_pop_up_window, self.metadata_menu))
        self.bu_rep.clicked.connect(partial(self.generate_pop_up_window, self.representativity_menu))
        self.bu_period.clicked.connect(partial(self.generate_pop_up_window, self.period_menu))

        # Setup MPL canvas of plots
        # set variable that blocks updating of MPL canvas until some data has been read
        self.block_MPL_canvas_updates = True
        self.mpl_canvas = MPLCanvas(self)

        # initialise layout options
        self.mpl_canvas.layout_options = ['None', 'distribution', 'metadata', 'periodic', 
                                          'periodic-violin', 'scatter', 'timeseries']
                                          
        # initialise configuration bar fields
        self.config_bar_initialisation = True
        self.update_configuration_bar_fields()
        self.config_bar_initialisation = False

        # Enable interactivity of functions which update MPL canvas
        # enable READ button
        self.bu_read.clicked.connect(self.handle_data_selection_update)
        # enable RESET button
        self.bu_reset.clicked.connect(self.reset_options)
        # enable interactivity of temporal colocation checkbox
        self.ch_colocate.stateChanged.connect(self.mpl_canvas.handle_temporal_colocate_update)
        # enable FILTER button
        self.bu_screen.clicked.connect(self.mpl_canvas.handle_data_filter_update)

        # enable updating of map z statistic
        self.cb_z_stat.currentTextChanged.connect(self.mpl_canvas.handle_map_z_statistic_update)
        self.cb_z1.currentTextChanged.connect(self.mpl_canvas.handle_map_z_statistic_update)
        self.cb_z2.currentTextChanged.connect(self.mpl_canvas.handle_map_z_statistic_update)

        # enable updating of experiment bias statistic
        self.cb_experiment_bias_type.currentTextChanged.connect(self.mpl_canvas.handle_experiment_bias_update)
        self.cb_experiment_bias_stat.currentTextChanged.connect(self.mpl_canvas.handle_experiment_bias_update)

        # enable interactivity of station selection checkboxes
        self.ch_select_all.stateChanged.connect(self.mpl_canvas.select_all_stations)
        self.ch_intersect.stateChanged.connect(self.mpl_canvas.select_intersect_stations)
        self.ch_extent.stateChanged.connect(self.mpl_canvas.select_extent_stations)

        # enable updating of plots layout
        self.cb_position_1.currentTextChanged.connect(self.handle_layout_update)
        self.cb_position_2.currentTextChanged.connect(self.handle_layout_update)
        self.cb_position_3.currentTextChanged.connect(self.handle_layout_update)
        self.cb_position_4.currentTextChanged.connect(self.handle_layout_update)
        self.cb_position_5.currentTextChanged.connect(self.handle_layout_update)

        # Generate MPL navigation toolbar
        self.navi_toolbar = NavigationToolbar(self.mpl_canvas, self)

        # add more buttons on the toolbar, next to the navi_toolbar
        self.savebutton = QtWidgets.QPushButton()
        self.savebutton.setFlat(True)
        self.savebutton.setToolTip("Save current instance of data and metadata")
        self.savebutton.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_data.png")))
        self.savebutton.setIconSize(QtCore.QSize(31, 35))
        self.savebutton.setStyleSheet("QPushButton { border: none;} QPushButton:hover "
                                      "{ border-width: 1px; border-style: solid; border-color: darkgrey; "
                                      "border-radius: 4px; background-color : white; }")
        self.savebutton.clicked.connect(self.savebutton_func)

        # add more buttons on the toolbar, next to the navi_toolbar
        self.conf_load = QtWidgets.QPushButton()
        self.conf_load.setFlat(True)
        self.conf_load.setToolTip("Load toolbar selections from configuration file")
        self.conf_load.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/conf_icon.png")))
        self.conf_load.setIconSize(QtCore.QSize(31, 35))
        self.conf_load.setStyleSheet("QPushButton { border: none;} QPushButton:hover "
                                      "{ border-width: 1px; border-style: solid; border-color: darkgrey; "
                                      "border-radius: 4px; background-color : white; }")
        self.conf_load.clicked.connect(self.conf_load_func)

        # position config bar, navigation toolbar and MPL canvas and elements in parent layout`
        hbox.addWidget(self.savebutton)
        hbox.addWidget(self.conf_load)
        hbox.addWidget(self.navi_toolbar)

        # add config bar and hbox to parent frame
        parent_layout.addLayout(config_bar)
        parent_layout.addLayout(hbox)

        # add MPL canvas of plots to parent frame
        parent_layout.addWidget(self.mpl_canvas)

        # hide plots if needed
        self.mpl_canvas.avoid_hiding_plots = False

        # update variable to inform plotting functions whether to use colocated data/or not
        check_state = self.ch_colocate.checkState()
        if check_state == QtCore.Qt.Checked:
            self.temporal_colocation = True
        else:
            self.temporal_colocation = False

        # starting from a configuration file?
        if self.from_conf: 
            # read
            self.handle_data_selection_update()
            # set fields available for filtering
            aux.representativity_conf(self)
            aux.period_conf(self)
            aux.metadata_conf(self)
            # filter
            self.mpl_canvas.handle_data_filter_update()

        # set finalised layout
        self.setLayout(parent_layout)
        # plot whole dashboard
        self.show()
        # maximise window to fit screen
        self.showMaximized()

    def savebutton_func(self):
        save_data(self.mpl_canvas)

    def conf_load_func(self):
        conf_dialogs(self)

    def generate_pop_up_window(self, menu_root):
        """Generate pop up window"""
        self.pop_up_window = PopUpWindow(menu_root, [], self.main_window_geometry)

    def update_configuration_bar_fields(self):
        """Define function that initialises/updates configuration bar fields"""

        # set variable to block interactive handling while updating config bar parameters
        self.block_config_bar_handling_updates = True

        #turn off some features if using non-GHOST data
        if aux.check_for_ghost(self.network):
            self.enable_ghost_buttons()
        else:
            self.disable_ghost_buttons()

        # set some default configuration values when initialising config bar
        if self.config_bar_initialisation:

            # parse initial config variables (checking for presence of key variables)
            aux.get_parameters(self)

            # set initial selected start-end date as default
            self.le_start_date.setText(str(self.start_date))
            self.le_end_date.setText(str(self.end_date))
            self.date_range_has_changed = False

            # set temporal colocation tickbox
            if self.temporal_colocation:
                self.ch_colocate.setCheckState(QtCore.Qt.Checked)
            else:
                self.ch_colocate.setCheckState(QtCore.Qt.Unchecked)

            # set initial time array, yearmonths and data_labels to be None 
            self.time_array = None
            self.yearmonths = None
            self.data_labels = None
            # set initial station references to be empty dict
            self.station_references = {}

            # set initial selected config variables as set .conf files or defaults
            self.selected_network = copy.deepcopy(self.network[0])
            self.selected_resolution = copy.deepcopy(self.resolution)
            self.selected_matrix = self.parameter_dictionary[self.species[0]]['matrix']
            self.selected_species = copy.deepcopy(self.species[0])

            # set variable stating first read
            self.first_read = True

            # create dictionary of available observational data inside date range
            aux.get_valid_obs_files_in_date_range(self, self.le_start_date.text(), self.le_end_date.text())

            # set qa / flags
            self.flags = aux.which_flags(self)
            self.qa = aux.which_qa(self)
            self.flag_menu['checkboxes']['remove_selected'] = copy.deepcopy(self.flags)
            self.qa_menu['checkboxes']['remove_selected'] = copy.deepcopy(self.qa)

        # if date range has changed then update available observational data dictionary
        if self.date_range_has_changed:
            aux.get_valid_obs_files_in_date_range(self, self.le_start_date.text(), self.le_end_date.text())

        # initialise/update fields - maintain previously selected values wherever possible
        # clear fields
        self.cb_network.clear()
        self.cb_resolution.clear()
        self.cb_matrix.clear()
        self.cb_species.clear()

        # if have no available observational data, return from function, updating variable informing that have no data
        if len(self.available_observation_data) == 0:
            self.no_data_to_read = True
            # unset variable to allow interactive handling from now
            self.block_config_bar_handling_updates = False
            return
        else:
            self.no_data_to_read = False

        # update network field
        available_networks = list(self.available_observation_data.keys())
        self.cb_network.addItems(available_networks)
        if self.selected_network in available_networks:
            self.cb_network.setCurrentText(self.selected_network)
        else:
            self.selected_network = self.cb_network.currentText()

        # update resolution field
        available_resolutions = list(self.available_observation_data[self.selected_network].keys())
        # set order of available resolutions
        available_resolutions = sorted(available_resolutions, key=aux.temporal_resolution_order_dict().__getitem__)
        self.cb_resolution.addItems(available_resolutions)
        if self.selected_resolution in available_resolutions:
            self.cb_resolution.setCurrentText(self.selected_resolution)
        else:
            self.selected_resolution = self.cb_resolution.currentText()

        # update matrix field
        available_matrices = sorted(self.available_observation_data[self.selected_network][self.selected_resolution])
        self.cb_matrix.addItems(available_matrices)
        if self.selected_matrix in available_matrices:
            self.cb_matrix.setCurrentText(self.selected_matrix)
        else:
            self.selected_matrix = self.cb_matrix.currentText()

        # update species field
        available_species = sorted(self.available_observation_data[self.selected_network][self.selected_resolution][self.selected_matrix])
        self.cb_species.addItems(available_species)
        if self.selected_species in available_species:
            self.cb_species.setCurrentText(self.selected_species)
        else:
            self.selected_species = self.cb_species.currentText()

        # update available experiments for selected fields
        aux.get_valid_experiments(self, self.le_start_date.text(), self.le_end_date.text(), self.selected_resolution,
                                  [self.selected_network], [self.selected_species])
        
        # update experiments -- keeping previously selected experiments if available
        if self.config_bar_initialisation:   
            experiments = aux.get_experiments(self)
            self.experiments_menu['checkboxes']['keep_selected'] = [experiment for experiment in experiments
                                                                    if experiment in self.experiments_menu['checkboxes']['map_vars']]
            self.experiments = {experiment:experiment_alias for experiment, experiment_alias in experiments.items()
                                if experiment in self.experiments_menu['checkboxes']['map_vars']}

        self.experiments_menu['checkboxes']['keep_selected'] = [previous_selected_experiment for
                                                                previous_selected_experiment in
                                                                self.experiments_menu['checkboxes']['keep_selected']
                                                                if previous_selected_experiment in
                                                                self.experiments_menu['checkboxes']['map_vars']]
        previous_experiments = self.experiments
        self.experiments = {exp:previous_experiments[exp] if exp in previous_experiments else exp for exp in self.experiments_menu['checkboxes']['keep_selected']}
        
        # update default qa
        default_qa = aux.which_qa(self, return_defaults=True)
        self.qa_menu['checkboxes']['remove_default'] = default_qa

        # initialise/update fields - maintain previously selected values wherever possible
        # clear fields
        self.cb_position_1.clear()
        self.cb_position_2.clear()
        self.cb_position_3.clear()
        self.cb_position_4.clear()
        self.cb_position_5.clear()

        # update position 1 in layout (always map)
        self.cb_position_1.addItems(['map'])
        if self.position_1 in ['map']:
            self.cb_position_1.setCurrentText(self.position_1)
        else:
            self.position_1 = self.cb_position_1.currentText()
        
        # update position 2 in layout
        self.cb_position_2.addItems(self.mpl_canvas.layout_options)
        if self.position_2 in self.mpl_canvas.layout_options:
            self.cb_position_2.setCurrentText(self.position_2)
        else:
            self.position_2 = self.cb_position_2.currentText()
        
        # update position 3 in layout
        self.cb_position_3.addItems(self.mpl_canvas.layout_options)
        if self.position_3 in self.mpl_canvas.layout_options:
            self.cb_position_3.setCurrentText(self.position_3)
        else:
            self.position_3 = self.cb_position_3.currentText()

        # update position 4 in layout
        self.cb_position_4.addItems(self.mpl_canvas.layout_options)
        if self.position_4 in self.mpl_canvas.layout_options:
            self.cb_position_4.setCurrentText(self.position_4)
        else:
            self.position_4 = self.cb_position_4.currentText()

        # update position 5 in layout
        self.cb_position_5.addItems(self.mpl_canvas.layout_options)
        if self.position_5 in self.mpl_canvas.layout_options:
            self.cb_position_5.setCurrentText(self.position_5)
        else:
            self.position_5 = self.cb_position_5.currentText()

        # unset variable to allow interactive handling from now
        self.block_config_bar_handling_updates = False

    def handle_config_bar_params_change(self, changed_param):
        """Define function which handles interactive updates of combo box fields"""

        if (changed_param != '') & (not self.block_config_bar_handling_updates):

            # get event origin source
            event_source = self.sender()
            # if network, resolution, matrix or species have changed then respective
            # current selection for the changed param
            if event_source == self.cb_network:
                self.selected_network = changed_param
            elif event_source == self.cb_resolution:
                self.selected_resolution = changed_param
            elif event_source == self.cb_matrix:
                self.selected_matrix = changed_param
                self.selected_species = sorted(list(
                    self.available_observation_data[self.selected_network][self.selected_resolution][selected_matrix].keys()))[0]
            elif event_source == self.cb_species:
                self.selected_species = changed_param

            # set variable to check if date range changes
            self.date_range_has_changed = False

            # check if start date/end date have changed
            if (event_source == self.le_start_date) or (event_source == self.le_end_date):
                self.date_range_has_changed = True

            # update configuration bar fields
            self.update_configuration_bar_fields()

    def handle_layout_update(self, changed_param):
        
        if (changed_param != '') & (not self.block_config_bar_handling_updates):
            
            # clear all previously plotted artists from selected station plots and hide axes 
            self.mpl_canvas.previous_selected_station_plots = self.mpl_canvas.selected_station_plots
            for plot_type in self.mpl_canvas.previous_selected_station_plots:
                ax = self.mpl_canvas.plot_axes[plot_type]
                if type(ax) == dict:
                    for sub_ax in ax.values():
                        self.mpl_canvas.remove_axis_elements(sub_ax, plot_type)
                else:
                    self.mpl_canvas.remove_axis_elements(ax, plot_type)

            # get event origin source
            event_source = self.sender()

            print('Updating layout')
            if event_source == self.cb_position_2:
                self.position_2 = changed_param
                [self.position_3, self.position_4, self.position_5] = self.avoid_duplicate_plots(self.position_2, 
                                                                                                 [self.position_3, 
                                                                                                  self.position_4, 
                                                                                                  self.position_5])
            elif event_source == self.cb_position_3:
                self.position_3 = changed_param
                [self.position_2, self.position_4, self.position_5] = self.avoid_duplicate_plots(self.position_3, 
                                                                                                 [self.position_2, 
                                                                                                  self.position_4, 
                                                                                                  self.position_5])
            elif event_source == self.cb_position_4:
                self.position_4 = changed_param
                [self.position_2, self.position_3, self.position_5] = self.avoid_duplicate_plots(self.position_4, 
                                                                                                 [self.position_2, 
                                                                                                  self.position_3, 
                                                                                                  self.position_5])
            elif event_source == self.cb_position_5:
                self.position_5 = changed_param
                [self.position_2, self.position_3, self.position_4] = self.avoid_duplicate_plots(self.position_5, 
                                                                                                 [self.position_2, 
                                                                                                  self.position_3, 
                                                                                                  self.position_4])

            # update layout
            self.update_plot_axes()

            # remove empty plots if READ has never been clicked
            if not hasattr(self.mpl_canvas, 'relative_selected_station_inds'):
                for plot_type in self.mpl_canvas.selected_station_plots:
                    ax = self.mpl_canvas.plot_axes[plot_type]
                    if type(ax) == dict:
                        for sub_ax in ax.values():
                            self.mpl_canvas.remove_axis_elements(sub_ax, plot_type)
                    else:
                        self.mpl_canvas.remove_axis_elements(ax, plot_type)
            
            # update characteristics (send warnings if some plots cannot be created)
            self.mpl_canvas.plot.set_plot_characteristics(self.mpl_canvas.selected_station_plots)

            # make plots
            self.mpl_canvas.avoid_hiding_plots = True
            self.mpl_canvas.update_associated_selected_station_plots()
            self.mpl_canvas.avoid_hiding_plots = False

            # update configuration bar fields
            self.update_configuration_bar_fields()

            # update buttons geometry
            self.update_buttons_geometry()
            
            # draw changes
            self.mpl_canvas.figure.canvas.draw()

        return None

    def avoid_duplicate_plots(self, current_plot_type, plots):
        
        if current_plot_type in plots:
            for i, plot_type in enumerate(plots):
                if plot_type == current_plot_type and plot_type != 'None':
                    plots[i] = 'None'

        return plots

    def update_plot_axes(self):

        # update plots
        self.mpl_canvas.selected_station_plots = [self.position_2, 
                                                  self.position_3, 
                                                  self.position_4, 
                                                  self.position_5]
       
        # drop empty plots from list
        self.mpl_canvas.selected_station_plots = list(filter(lambda val: val != 'None', 
                                                             self.mpl_canvas.selected_station_plots))

        # create fixed plot axes on grid
        self.mpl_canvas.plot_axes = {}
        self.mpl_canvas.plot_axes['map'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((2, 0), 
                                                                              rowspan=44, colspan=42), projection=self.mpl_canvas.plotcrs)
        self.mpl_canvas.plot_axes['cb'] = self.mpl_canvas.figure.add_axes([0.0375, 0.53, 0.375, 0.02])
        self.mpl_canvas.plot_axes['legend'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((0, 47), 
                                                                                 rowspan=8,  colspan=53))
        
        # create swappable plot axes on grid
        # position 2 (top right)
        if (self.position_2 == 'periodic') or (self.position_2 == 'periodic-violin'):
            self.mpl_canvas.plot_axes[self.position_2] = {}
            self.mpl_canvas.plot_axes[self.position_2]['hour'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((10, 50), rowspan=17, colspan=50))
            self.mpl_canvas.plot_axes[self.position_2]['dayofweek'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((32, 82), rowspan=17, colspan=18))
            self.mpl_canvas.plot_axes[self.position_2]['month'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((32, 50), rowspan=17, colspan=30))
        elif self.position_2 != 'None':
            self.mpl_canvas.plot_axes[self.position_2] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((12, 50), rowspan=34, colspan=50))
        
        # position 3 (bottom left)
        if (self.position_3 == 'periodic') or (self.position_3 == 'periodic-violin'):
            self.mpl_canvas.plot_axes[self.position_3] = {}
            self.mpl_canvas.plot_axes[self.position_3]['hour'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((56, 0), rowspan=20, colspan=30))
            self.mpl_canvas.plot_axes[self.position_3]['dayofweek'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((80, 20), rowspan=20, colspan=10))
            self.mpl_canvas.plot_axes[self.position_3]['month'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((80, 0), rowspan=20, colspan=18))
        elif self.position_3 != 'None':
            self.mpl_canvas.plot_axes[self.position_3] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((56, 0),  rowspan=44, colspan=30))
        
        # position 4 (bottom centre)
        if (self.position_4 == 'periodic') or (self.position_4 == 'periodic-violin'):
            self.mpl_canvas.plot_axes[self.position_4] = {}
            self.mpl_canvas.plot_axes[self.position_4]['hour'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((56, 35), rowspan=20, colspan=30))
            self.mpl_canvas.plot_axes[self.position_4]['dayofweek'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((80, 55), rowspan=20, colspan=10))
            self.mpl_canvas.plot_axes[self.position_4]['month'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((80, 35), rowspan=20, colspan=18))
        elif self.position_4 != 'None':
            self.mpl_canvas.plot_axes[self.position_4] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((56, 35),  rowspan=44, colspan=30))
        
        # position 5 (bottom right)
        if (self.position_5 == 'periodic') or (self.position_5 == 'periodic-violin'):
            self.mpl_canvas.plot_axes[self.position_5] = {}
            self.mpl_canvas.plot_axes[self.position_5]['hour'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((56, 70), rowspan=20, colspan=30))
            self.mpl_canvas.plot_axes[self.position_5]['dayofweek'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((80, 90), rowspan=20, colspan=10))
            self.mpl_canvas.plot_axes[self.position_5]['month'] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((80, 70), rowspan=20, colspan=18))
        elif self.position_5 != 'None':
            self.mpl_canvas.plot_axes[self.position_5] = self.mpl_canvas.figure.add_subplot(self.mpl_canvas.gridspec.new_subplotspec((56, 70),  rowspan=44, colspan=30))

    def handle_data_selection_update(self):
        """Define function which handles update of data selection
        and MPL canvas upon pressing of READ button
        """

        # if have no data to read, then do not read any data
        if self.no_data_to_read:
            return

        # Update mouse cursor to a waiting cursor
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        # set variable that blocks updating of MPL canvas until all data has been updated
        self.block_MPL_canvas_updates = True

        #set previous active variables
        self.previous_start_date = self.start_date
        self.previous_end_date = self.end_date
        self.previous_network = self.network
        self.previous_resolution = self.resolution
        self.previous_species = self.species
        self.previous_experiments = self.experiments
        self.previous_qa = self.qa
        self.previous_flags = self.flags
        self.previous_data_labels = self.data_labels

        #set new active variables as selected variables from menu
        self.start_date = int(self.le_start_date.text())
        self.end_date = int(self.le_end_date.text())
        self.network = [self.selected_network]
        self.resolution = self.selected_resolution
        self.species = [self.selected_species]  
        self.experiments = {exp:self.previous_experiments[exp] if exp in self.previous_experiments else exp for exp in self.experiments_menu['checkboxes']['keep_selected']}
        self.qa = copy.deepcopy(self.qa_menu['checkboxes']['remove_selected'])
        self.flags = copy.deepcopy(self.flag_menu['checkboxes']['remove_selected'])
        self.data_labels = ['observations'] + list(self.experiments.keys())
        self.networkspeci = '{}-{}'.format(self.network[0],self.species[0])

        #set read operations to be empty list initially
        read_operations = []

        # if first read then need to read all data
        if self.first_read:
            read_operations = ['reset']

        # determine if any of the key variables have changed 
        # (network, resolution, species, qa, flags)
        # if any have changed, observations and any selected experiments have to be re-read entirely
        elif (self.network[0] != self.previous_network[0]) or (
                self.resolution != self.previous_resolution) or (
                self.species[0] != self.previous_species[0]) or (
                np.array_equal(self.qa, self.previous_qa) == False) or (
                np.array_equal(self.flags, self.previous_flags) == False):
            read_operations = ['reset']

        # key variables have not changed, has start/end date?
        else:
            # determine if start date/end date have changed
            if (self.start_date != self.previous_start_date) or (
                    self.end_date != self.previous_end_date):
                # if date range has changed then determine type of overlap with previous date range
                # no overlap (i.e. start date >= than previous end date, or end date <= than previous start date)?
                if (self.start_date >= self.previous_end_date) or (
                        self.end_date <= self.previous_start_date):
                    read_operations = ['reset']
                # data range fully inside previous data range (i.e. start date later and end date earlier)?
                elif (self.start_date > self.previous_start_date) & (
                        self.end_date < self.previous_end_date):
                    read_operations = ['cut_left','cut_right']
                # need to read data on left edge and right edge of previous date range
                # (i.e. start date earlier and end date later)?
                elif (self.start_date < self.previous_start_date) & (
                        self.end_date > self.previous_end_date):
                    read_operations = ['read_left','read_right']
                # need to read data on left edge and cut on right edge of previous date range
                # (i.e. start date earlier and end date earlier)?
                elif (self.start_date < self.previous_start_date) & (
                        self.end_date < self.previous_end_date):
                    read_operations = ['read_left','cut_right']
                # need to cut data on left edge and read data on right edge of previous date range
                # (i.e. start date later and end date later)?
                elif (self.start_date > self.previous_start_date) & (
                        self.end_date > self.previous_end_date):
                    read_operations = ['cut_left','read_right']
                # need to read data on left edge of previous date range (i.e. start date earlier)?
                elif self.start_date < self.previous_start_date:
                    read_operations = ['read_left']
                # need to read data on right edge of previous date range (i.e. end date later)?
                elif self.end_date > self.previous_end_date:
                    read_operations = ['read_right']
                # need to cut data on left edge of previous date range (i.e. start date later)?
                elif self.start_date > self.previous_start_date:
                    read_operations = ['cut_left']
                # need to cut data on right edge of previous date range (i.e. end date earlier)?
                elif self.end_date < self.previous_end_date:
                    read_operations = ['cut_right']

        # determine if any experiments need removing or reading 
        experiments_to_remove = [experiment for experiment in self.previous_experiments if experiment not in self.experiments]
        experiments_to_read = [experiment for experiment in self.experiments if experiment not in self.previous_experiments]
        if 'reset' not in read_operations:
            if len(experiments_to_remove) > 0:
                read_operations.append('remove_exp')
            if len(experiments_to_read) > 0:
                read_operations.append('read_exp')

        # has date range changed?
        if len(read_operations) > 0:

            # set current time array, as previous time array
            self.previous_time_array = self.time_array
            # set current station references, as previous station references
            self.previous_station_references = self.station_references
            # set current relevant yearmonths, as previous relevant yearmonths
            self.previous_yearmonths = self.yearmonths

            # read data
            self.datareader.read_setup(read_operations, experiments_to_remove=experiments_to_remove, experiments_to_read=experiments_to_read)

            #clear canvas entirely if have no valid data
            if self.clear_canvas:
                # clear axes
                for plot_type, ax in self.mpl_canvas.plot_axes.items():
                    if type(ax) == dict:
                        for sub_ax in ax.values():
                            self.mpl_canvas.remove_axis_elements(sub_ax, plot_type)
                    else:
                        self.mpl_canvas.remove_axis_elements(ax, plot_type)  
                # update MPL canvas
                self.mpl_canvas.figure.canvas.draw()  
                # restore mouse cursor to normal
                QtWidgets.QApplication.restoreOverrideCursor()    
                return

            # update fields available for filtering
            aux.update_representativity_fields(self)
            aux.update_period_fields(self)
            aux.update_metadata_fields(self)
            
            #update relevant temporal resolutions 
            self.relevant_temporal_resolutions = aux.get_relevant_temporal_resolutions(self.resolution)

        # if species has changed, or first read, update default species specific lower/upper limits
        if (self.first_read) or (self.species[0] != self.previous_species[0]):
            # update default lower/upper species specific limits and filter data outside limits
            species_lower_limit, species_upper_limit = aux.which_bounds(self, self.species[0])
            # set default limits
            self.le_minimum_value.setText(str(species_lower_limit))
            self.le_maximum_value.setText(str(species_upper_limit))

        # run function to filter data outside lower/upper limits, not using desired
        # measurement methods, and < desired minimum data availability
        self.mpl_canvas.handle_data_filter_update()

        # update map z combobox fields based on data in memory
        # generate lists of basic and basis+bias statistics for using in the z statistic combobox
        self.basic_z_stats = np.array(list(
            OrderedDict(sorted(basic_stats.items(), key=lambda x: x[1]['order'])).keys()))
        self.basic_and_bias_z_stats = np.append(self.basic_z_stats, list(
            OrderedDict(sorted(expbias_stats.items(), key=lambda x: x[1]['order'])).keys()))

        # generate list of sorted z1/z2 data arrays names in memory, putting observations
        # before experiments, and empty string item as first element in z2 array list
        # (for changing from 'difference' statistics to 'absolute')
        if len(self.data_labels) == 1:  
            self.z1_arrays = np.array(['observations'])
        else:
            self.z1_arrays = np.array(self.data_labels)
        self.z2_arrays = np.append([''], self.z1_arrays)

        # initialise map z statistic comboboxes
        self.mpl_canvas.handle_map_z_statistic_update()

        # update experiment bias combobox fields based on data in memory
        # if have no experiment data, all fields are empty
        self.block_config_bar_handling_updates = True
        if len(self.data_labels) == 1: 
            self.cb_experiment_bias_type.clear()
            self.cb_experiment_bias_stat.clear()
            self.block_config_bar_handling_updates = False
        # else, generate combobox lists
        else:
            # set all experiment bias types
            self.cb_experiment_bias_type.clear()
            self.cb_experiment_bias_type.addItems(np.array(['Station Median']))
            self.block_config_bar_handling_updates = False
            self.mpl_canvas.handle_experiment_bias_update()
    
        # reset station select checkboxes to be unchecked
        self.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
        self.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
        self.ch_extent.setCheckState(QtCore.Qt.Unchecked)

        # unset variable to allow updating of MPL canvas
        self.block_MPL_canvas_updates = False

        # update MPL canvas
        self.mpl_canvas.update_MPL_canvas()

        # if first read, then set this now to be False
        if self.first_read:
            self.first_read = False

        # Restore mouse cursor to normal
        QtWidgets.QApplication.restoreOverrideCursor()

    def reset_options(self):
        """Resets all filter fields to initial values"""

        if self.block_MPL_canvas_updates:
            return

        # set mouse cursor to hourglass
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        # reset representativity fields        
        aux.init_representativity(self)
        aux.update_representativity_fields(self)

        # reset period fields 
        aux.init_period(self)
        aux.update_period_fields(self)

        # reset metadata
        aux.init_metadata(self)
        aux.update_metadata_fields(self)

        # reset bounds
        species_lower_limit = np.float32(self.parameter_dictionary[self.species[0]]['extreme_lower_limit'])
        species_upper_limit = np.float32(self.parameter_dictionary[self.species[0]]['extreme_upper_limit'])
        
        # set default limits
        self.le_minimum_value.setText(str(species_lower_limit))
        self.le_maximum_value.setText(str(species_upper_limit))
        
        # unfilter data
        self.mpl_canvas.handle_data_filter_update()
        
        # Restore mouse cursor to normal
        QtWidgets.QApplication.restoreOverrideCursor()

    def disable_ghost_buttons(self):
        """Disable button related only to ghost data"""
        
        # change background-color to indicate that it's nonusable
        self.bu_flags.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        self.bu_QA.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        self.bu_period.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        
        # disable buttons
        self.bu_flags.setEnabled(False)
        self.bu_QA.setEnabled(False)
        self.bu_period.setEnabled(False)
        
    def enable_ghost_buttons(self):
        """Enable button related only to ghost data"""

        # enable buttons        
        self.bu_flags.setEnabled(True)
        self.bu_QA.setEnabled(True)
        self.bu_period.setEnabled(True)

# generate Providentia dashboard
def main(**kwargs):
    """Main function"""
    q_app = QtWidgets.QApplication(sys.argv)
    q_app.setStyle("Fusion")
    ProvidentiaMainWindow(**kwargs)
    sys.exit(q_app.exec_())