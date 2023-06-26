""" Module which provides main window """
from .configuration import ProvConfiguration
from .canvas import MPLCanvas
from .toolbar import NavigationToolbar
from .dashboard_aux import ComboBox, QVLine, PopUpWindow, InputDialog
from .dashboard_aux import set_formatting
from .read import DataReader
from .read_aux import get_default_qa, get_frequency_code
from providentia import aux

import os
import copy
import datetime
import json
import sys
import matplotlib
import mpl_toolkits.axisartist.floating_axes as fa
from matplotlib.projections import PolarAxes
from functools import partial
from collections import OrderedDict
from weakref import WeakKeyDictionary

from PyQt5 import QtCore, QtWidgets, QtGui
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))
formatting_dict = json.load(open(os.path.join(CURRENT_PATH, 'conf/stylesheet.json')))


class ProvidentiaMainWindow(QtWidgets.QWidget):
    """ Class that generates Providentia dashboard. """

    # create signals that are fired upon resizing/moving of main Providentia window
    resized = QtCore.pyqtSignal()
    move = QtCore.pyqtSignal()

    def __init__(self, **kwargs):

        # allow access to methods of parent class QtWidgets.QWidget
        super(ProvidentiaMainWindow, self).__init__()

        # initialise default configuration variables
        # modified by commandline arguments, if given
        provconf = ProvConfiguration(self, **kwargs)

        # update variables from config file (if available)
        self.from_conf = False
        self.current_config = {}
        
        if ('config' in kwargs) and (os.path.exists(kwargs['config'])):
            if 'section' in kwargs:
                # config and section defined 
                aux.load_conf(self, fpath=kwargs['config'])
                if kwargs['section'] in self.all_sections:
                    self.from_conf = True
                    self.current_config = self.sub_opts[kwargs['section']]
                else:
                    error = 'Error: The section specified in the command line ({0}) does not exist.'.format(kwargs['section'])
                    tip = 'Tip: For subsections, add the name of the parent section followed by an interpunct (·) '
                    tip += 'before the subsection name (e.g. SECTIONA·Spain).'
                    sys.exit(error + '\n' + tip)

            elif 'section' not in kwargs:
                # config defined, section undefined
                aux.load_conf(self, fpath=kwargs['config'])    
                all_sections = self.sub_opts.keys()
                
                if len(all_sections) == 1:
                    okpressed = False
                    selected_section = list(all_sections)[0]
                else:
                    title = 'Sections'
                    msg = 'Select section to load'
                    dialog = InputDialog(self, title, msg, all_sections)
                    selected_section, okpressed = dialog.selected_option, dialog.okpressed

                if okpressed or (len(all_sections) == 1):
                    self.from_conf = True
                    self.current_config = self.sub_opts[selected_section]

        elif ('config' in kwargs) and (not os.path.exists(kwargs['config'])):     
            error = 'Error: The path to the configuration file specified in the command line does not exist.'
            sys.exit(error)
        
        # set initial filter species
        self.previous_filter_species = {}
        self.filter_species = {}

        # update variables from defined config file
        if self.current_config:
            for k, val in self.current_config.items():
                setattr(self, k, provconf.parse_parameter(k, val))

        # now all variables have been parsed, check validity of those, throwing errors where necessary
        provconf.check_validity()
        
        # load characteristics per plot type
        # check for self defined plot characteristics file
        if self.plot_characteristics_filename == '':
            self.plot_characteristics_filename = os.path.join(
                CURRENT_PATH, 'conf/plot_characteristics_dashboard.json')
        self.plot_characteristics_templates = json.load(open(self.plot_characteristics_filename))

        # error when using wrong custom plot characteristics path to launch dashboard
        if 'header' in self.plot_characteristics_templates.keys():
            msg = 'It is not possible to use the offline plot characteristics path to launch the dashboard. Consider adding another path to plot_characteristics_filename, as in: '
            msg += 'plot_characteristics_filename = dashboard:/path/plot_characteristics_dashboard.json, offline:/path/plot_characteristics_offline.json.'
            sys.exit(msg)

        # arguments are only local
        self.main_window_geometry = None

        # create dictionary of all available observational GHOST data
        self.all_observation_data = aux.get_ghost_observational_tree(self)

        # load dictionary with non-GHOST esarchive files to read
        nonghost_observation_data_json = json.load(open(os.path.join(CURRENT_PATH, 'conf/nonghost_files.json')))
        # merge to existing GHOST observational data dict if we have the path
        if self.nonghost_root is not None:
            nonghost_observation_data = aux.get_nonghost_observational_tree(self, nonghost_observation_data_json)
            self.all_observation_data = {**self.all_observation_data, **nonghost_observation_data}

        # initialise DataReader
        self.datareader = DataReader(self)

        # initialise UI
        self.init_ui(**kwargs)

        # setup callback events upon resizing/moving of Providentia window
        self.resized.connect(self.get_geometry)
        self.move.connect(self.get_geometry)

    def resizeEvent(self, event):
        """ Function to overwrite default PyQt5 resizeEvent function --> for calling get_geometry. """
        
        self.resized.emit()
        
        return super(ProvidentiaMainWindow, self).resizeEvent(event)

    def moveEvent(self, event):
        """ Function to overwrite default PyQt5 moveEvent function --> for calling get_geometry. """
        
        self.move.emit()
        
        return super(ProvidentiaMainWindow, self).moveEvent(event)

    def get_geometry(self):
        """ Update current geometry of main Providentia window and buttons. """

        # get geometry of main window
        self.main_window_geometry = copy.deepcopy(self.geometry())

        # update geometry of setting menus
        self.update_buttons_geometry()

    def update_buttons_geometry(self):
        """ Update current geometry of buttons. """
        
        for position_ii, position in enumerate([self.position_1, self.position_2, self.position_3, 
                                                self.position_4, self.position_5]):
            for menu_button, save_button, element in zip(self.mpl_canvas.menu_buttons, 
                                                         self.mpl_canvas.save_buttons, 
                                                         self.mpl_canvas.elements):

                # get plot type
                plot_type = menu_button.objectName().split('_menu')[0]
                if position == 'periodic-violin':
                    position = 'periodic_violin'

                if position == plot_type or position == 'None':
                    
                    # calculate proportional geometry of buttons respect main window
                    x = (self.mpl_canvas.plot_characteristics_templates['general']['settings_menu']['position_'
                         + str(position_ii+1)]['x'] * self.main_window_geometry.width()) / 1848
                    y = (self.mpl_canvas.plot_characteristics_templates['general']['settings_menu']['position_' 
                         + str(position_ii+1)]['y'] * self.main_window_geometry.height()) / 1016
                    
                    # get geometries
                    old_button_geometry = QtCore.QRect(menu_button.x(), menu_button.y(), 18, 18)
                    new_button_geometry = QtCore.QRect(x, y, 18, 18)
                    
                    # apply new geometry to menu and save buttons
                    menu_button.setGeometry(new_button_geometry)
                    save_button.setGeometry(menu_button.x()-25, menu_button.y(), 20, 20)

                    # apply new geometry to layout buttons
                    if position_ii > 0:
                        cb_position = getattr(self, 'cb_position_{}'.format(position_ii+1))
                        # layout selector at position 2 needs to be further from menu
                        if (position_ii + 1) == 2:
                            width_diff = 915
                        else:
                            width_diff = 570
                        height_diff = 1
                        width = (width_diff * self.main_window_geometry.width()) / 1848
                        height = (height_diff * self.main_window_geometry.height()) / 1848
                        cb_position.move(menu_button.x()-width, menu_button.y()+height)

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
        """ Initialise user interface. """

        print("Starting Providentia dashboard...")

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
        self.setStyleSheet("QToolTip { font: %spt %s}" % (formatting_dict['tooltip']['font']['size'],
                                                          formatting_dict['tooltip']['font']['style']))

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
        # data selection section
        self.lb_data_selection = set_formatting(QtWidgets.QLabel(self, text="Data Selection"),
                                                formatting_dict['title_menu'])
        self.lb_data_selection.setToolTip('Setup configuration of data to read into memory')
        self.cb_network = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_network.setToolTip('Select providing observational data network. '
                                   'Names starting with * indicate non-GHOST datasets')
        self.cb_resolution = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_resolution.setToolTip('Select temporal resolution of data')
        self.cb_matrix = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_matrix.setToolTip('Select data matrix')
        self.cb_species = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_species.setToolTip('Select species')
        self.le_start_date = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_start_date.setToolTip('Set data start date: YYYYMMDD')
        self.le_end_date = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_end_date.setToolTip('Set data end date: YYYYMMDD')
        self.bu_QA = set_formatting(QtWidgets.QPushButton('QA', self), formatting_dict['button_menu'])
        self.bu_QA.setToolTip('Select standardised quality assurance flags to filter by')
        self.bu_flags = set_formatting(QtWidgets.QPushButton('FLAGS', self), formatting_dict['button_menu'])
        self.bu_flags.setToolTip('Select standardised data reporter provided flags to filter by')
        self.bu_experiments = set_formatting(QtWidgets.QPushButton('EXPS', self), formatting_dict['button_menu'])
        self.bu_experiments.setToolTip('Select experiment/s data to read')
        self.bu_multispecies = set_formatting(QtWidgets.QPushButton('MULTI', self), formatting_dict['button_menu'])
        self.bu_multispecies.setToolTip('Select data to filter by')
        self.bu_read = set_formatting(QtWidgets.QPushButton('READ', self), formatting_dict['button_menu'])
        self.bu_read.setStyleSheet("color: green;")
        self.bu_read.setToolTip('Read selected configuration of data into memory')
        self.vertical_splitter_1 = QVLine()
        self.vertical_splitter_1.setMaximumWidth(20)

        # filters section
        self.lb_data_filter = set_formatting(QtWidgets.QLabel(self, text="Filters"), formatting_dict['title_menu'])
        self.lb_data_filter.setToolTip('Select criteria to filter data by')
        self.bu_rep = set_formatting(QtWidgets.QPushButton('% REP', self), formatting_dict['button_menu'])
        self.bu_rep.setToolTip('Select % desired representativity in data across '
                               'whole record and for specific temporal periods')
        self.bu_meta = set_formatting(QtWidgets.QPushButton('META', self), formatting_dict['button_menu'])
        self.bu_meta.setToolTip('Select metadata to filter by')
        self.bu_reset = set_formatting(QtWidgets.QPushButton('RESET', self), formatting_dict['button_menu'])
        self.bu_reset.setToolTip('Reset filter fields to initial values')
        self.bu_reset.setStyleSheet("color: red;")
        self.bu_period = set_formatting(QtWidgets.QPushButton('PERIOD', self), formatting_dict['button_menu'])
        self.bu_period.setToolTip('Select data in specific periods')
        self.bu_filter = set_formatting(QtWidgets.QPushButton('FILTER', self), formatting_dict['button_menu'])
        self.bu_filter.setStyleSheet("color: blue;")
        self.bu_filter.setToolTip('Filter data')
        self.lb_data_bounds = set_formatting(QtWidgets.QLabel(self, text="Bounds"), formatting_dict['label_menu'])
        self.lb_data_bounds.setToolTip('Set lower/upper bounds of data')
        self.le_minimum_value = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_minimum_value.setToolTip('Set lower bound of data')
        self.le_maximum_value = set_formatting(QtWidgets.QLineEdit(self), formatting_dict['lineedit_menu'])
        self.le_maximum_value.setToolTip('Set upper bound of data')
        self.vertical_splitter_2 = QVLine()
        self.vertical_splitter_2.setMaximumWidth(20)

        # station aggregation section
        self.lb_aggregation = set_formatting(QtWidgets.QLabel(self, text="Aggregation"),
                                             formatting_dict['title_menu'])
        self.lb_aggregation.setToolTip('Select the statistic to aggregate the stations data')
        self.cb_aggregation_statistic = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_aggregation_statistic.setToolTip('Select statistic')
        self.vertical_splitter_3 = QVLine()
        self.vertical_splitter_3.setMaximumWidth(20)
        
        # colocation section
        self.lb_colocate = set_formatting(QtWidgets.QLabel(self, text="Colocation"), formatting_dict['title_menu'])
        self.lb_colocate.setToolTip('Set colocation')
        self.ch_colocate = set_formatting(QtWidgets.QCheckBox("Temporal"), formatting_dict['checkbox_menu'])
        self.ch_colocate.setToolTip('Temporally colocate observational/experiment data')
        self.vertical_splitter_4 = QVLine()
        self.vertical_splitter_4.setMaximumWidth(20)

        # resampling section
        self.lb_resampling = set_formatting(QtWidgets.QLabel(self, text="Resampling"), formatting_dict['title_menu'])
        self.lb_resampling.setToolTip('Set resampling options')
        self.cb_resampling_resolution = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_resampling_resolution.setToolTip('Select temporal resolution to resample the data to')
        self.vertical_splitter_5 = QVLine()
        self.vertical_splitter_5.setMaximumWidth(20)

        # station selection section
        self.lb_station_selection = set_formatting(QtWidgets.QLabel(self, text="Site Selection"),
                                                   formatting_dict['title_menu'])
        self.lb_station_selection.setToolTip('Select stations')
        self.ch_select_all = set_formatting(QtWidgets.QCheckBox("All"), formatting_dict['checkbox_menu'])
        self.ch_select_all.setToolTip('Select all stations')
        self.ch_intersect = set_formatting(QtWidgets.QCheckBox("Intersect"), formatting_dict['checkbox_menu'])
        self.ch_intersect.setToolTip('Select stations that intersect with all loaded model domains')
        self.ch_extent = set_formatting(QtWidgets.QCheckBox("Extent"), formatting_dict['checkbox_menu'])
        self.ch_extent.setToolTip('Select stations that are within the map extent')

        # position objects on gridded configuration bar
        # data selection section
        config_bar.addWidget(self.lb_data_selection, 0, 0, 1, 2, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_network, 1, 0, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_resolution, 2, 0, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_matrix, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_species, 1, 2, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.le_start_date, 2, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.le_end_date, 2, 2, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_QA, 1, 3, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_flags, 2, 3, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_experiments, 1, 4, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_multispecies, 2, 4, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_read, 3, 4, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.vertical_splitter_1, 0, 5, 4, 1, QtCore.Qt.AlignLeft)

        # filters section
        config_bar.addWidget(self.lb_data_filter, 0, 6, 1, 2, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.lb_data_bounds, 1, 6, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.le_minimum_value, 1, 7, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.le_maximum_value, 1, 8, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_rep, 2, 6, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_period, 2, 7, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_meta, 2, 8, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_reset, 3, 7, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.bu_filter, 3, 8, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.vertical_splitter_2, 0, 9, 4, 1, QtCore.Qt.AlignLeft)

        # station aggregation section
        config_bar.addWidget(self.lb_aggregation, 0, 10, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_aggregation_statistic, 1, 10, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.vertical_splitter_3, 0, 11, 4, 1, QtCore.Qt.AlignLeft)

        # colocation section
        config_bar.addWidget(self.lb_colocate, 0, 12, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_colocate, 1, 12, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.vertical_splitter_4, 0, 13, 4, 1, QtCore.Qt.AlignLeft)

        # resampling section
        config_bar.addWidget(self.lb_resampling, 0, 14, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_resampling_resolution, 1, 14, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.vertical_splitter_5, 0, 15, 4, 1, QtCore.Qt.AlignLeft)

        # station selection section
        config_bar.addWidget(self.lb_station_selection, 0, 16, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_select_all, 1, 16, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_intersect, 2, 16, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_extent, 3, 16, QtCore.Qt.AlignLeft)

        # enable dynamic updating of configuration bar fields which filter data files
        self.cb_network.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_resolution.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_matrix.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_species.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.le_start_date.textChanged.connect(self.handle_config_bar_params_change)
        self.le_end_date.textChanged.connect(self.handle_config_bar_params_change)
        self.cb_resampling_resolution.currentTextChanged.connect(self.handle_config_bar_params_change)

        # setup pop-up window menu tree for flags, qa, experiments, 
        # % data representativity, data periods and metadata
        aux.init_flags(self)
        aux.init_qa(self)
        aux.init_experiments(self)
        aux.init_multispecies(self)
        aux.init_representativity(self)
        aux.init_period(self)
        self.metadata_vars_to_read = []
        aux.init_metadata(self)

        # Setup MPL canvas of plots
        # set variable that blocks updating of MPL canvas until some data has been read
        self.block_MPL_canvas_updates = True
        self.mpl_canvas = MPLCanvas(self)
              
        # initialise configuration bar fields
        self.config_bar_initialisation = True
        self.update_configuration_bar_fields()
        self.config_bar_initialisation = False

        # initialise multispecies tab
        self.multispecies_initialisation = True

        # launch with configuriton file?
        if self.from_conf: 
            
            # read
            self.handle_data_selection_update()

            # set filtered multispecies if any
            aux.multispecies_conf(self)

            # set fields available for filtering
            aux.representativity_conf(self)
            aux.period_conf(self)
            aux.metadata_conf(self)
            self.mpl_canvas.handle_data_filter_update()

        # enable pop up configuration windows
        self.bu_flags.clicked.connect(partial(self.generate_pop_up_window, self.flag_menu))
        self.bu_QA.clicked.connect(partial(self.generate_pop_up_window, self.qa_menu))
        self.bu_experiments.clicked.connect(partial(self.generate_pop_up_window, self.experiments_menu))
        self.bu_multispecies.clicked.connect(partial(self.generate_pop_up_window, self.multispecies_menu))
        self.bu_meta.clicked.connect(partial(self.generate_pop_up_window, self.metadata_menu))
        self.bu_rep.clicked.connect(partial(self.generate_pop_up_window, self.representativity_menu))
        self.bu_period.clicked.connect(partial(self.generate_pop_up_window, self.period_menu))

        # Enable interactivity of functions which update MPL canvas
        # enable READ button
        self.bu_read.clicked.connect(self.handle_data_selection_update)

        # enable RESET button
        self.bu_reset.clicked.connect(self.reset_options)

        # enable FILTER button
        self.bu_filter.clicked.connect(self.mpl_canvas.handle_data_filter_update)

        # enable aggregation by changing the statistic
        self.cb_aggregation_statistic.currentTextChanged.connect(self.mpl_canvas.handle_aggregation_update)

        # enable interactivity of temporal colocation checkbox
        self.ch_colocate.stateChanged.connect(self.mpl_canvas.handle_temporal_colocate_update)

        # enable resampling by changing the temporal resolution
        self.cb_resampling_resolution.currentTextChanged.connect(self.mpl_canvas.handle_resampling_update)

        # enable interactivity of station selection checkboxes
        self.ch_select_all.stateChanged.connect(self.mpl_canvas.select_all_stations)
        self.ch_intersect.stateChanged.connect(self.mpl_canvas.select_intersect_stations)
        self.ch_extent.stateChanged.connect(self.mpl_canvas.select_extent_stations)

        # Generate MPL navigation toolbar
        self.navi_toolbar = NavigationToolbar(read_instance=self, canvas_instance=self.mpl_canvas)
        self.navi_toolbar._nav_stack.push(
            WeakKeyDictionary({self.mpl_canvas.plot_axes['map']: (self.mpl_canvas.plot_axes['map']._get_view(), 
                                                                  (self.mpl_canvas.plot_axes['map'].get_position(True).frozen(), 
                                                                   self.mpl_canvas.plot_axes['map'].get_position().frozen()))}))

        # position navigation toolbar in parent layout
        hbox.addWidget(self.navi_toolbar)

        # add config bar and hbox to parent frame
        parent_layout.addLayout(config_bar)
        parent_layout.addLayout(hbox)

        # add MPL canvas of plots to parent frame
        parent_layout.addWidget(self.mpl_canvas)

        # set finalised layout
        self.setLayout(parent_layout)

        # plot whole dashboard
        self.show()

        # maximise window to fit screen
        self.showMaximized()

    def generate_pop_up_window(self, menu_root):
        """ Generate pop up window. """
        
        self.pop_up_window = PopUpWindow(self, menu_root, [], self.main_window_geometry)

    def update_configuration_bar_fields(self):
        """ Function that initialises/updates configuration bar fields. """

        # set variable to block interactive handling while updating config bar parameters
        self.block_config_bar_handling_updates = True

        # set variable to avoid updating the canvas while updating config bar parameters
        self.block_MPL_canvas_updates = True

        # set some default configuration values when initialising config bar
        if self.config_bar_initialisation:

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
            self.selected_resampling_resolution = copy.deepcopy(self.resampling_resolution)
            self.selected_resampling = copy.deepcopy(self.resampling)
            self.selected_aggregation_statistic = copy.deepcopy(self.aggregation_statistic)

            # set initial filter species in widgets as empty dictionaries
            self.selected_widget_network = dict()
            self.selected_widget_matrix = dict()
            self.selected_widget_species = dict()
            self.selected_widget_lower = dict()
            self.selected_widget_upper = dict()
            self.selected_widget_filter_species_fill_value = dict()
            self.selected_widget_apply = dict()
            
            # set variable stating first read
            self.first_read = True

            # create dictionary of available observational data inside date range
            aux.get_valid_obs_files_in_date_range(self, self.le_start_date.text(), self.le_end_date.text())

            # update qa / flags checkboxes 
            self.flag_menu['checkboxes']['remove_selected'] = copy.deepcopy(self.flags)
            self.qa_menu['checkboxes']['remove_selected'] = copy.deepcopy(self.qa_per_species[self.selected_species])
             
        # if date range has changed then update available observational data dictionary
        if self.date_range_has_changed:
            aux.get_valid_obs_files_in_date_range(self, self.le_start_date.text(), self.le_end_date.text())

        # initialise/update fields - maintain previously selected values wherever possible
        # clear fields
        self.cb_network.clear()
        self.cb_resolution.clear()
        self.cb_matrix.clear()
        self.cb_species.clear()
        self.cb_resampling_resolution.clear()
        self.cb_aggregation_statistic.clear()

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

        # turn off some features if using non-GHOST data
        if aux.check_for_ghost(self.selected_network):
            self.enable_ghost_buttons()
        else:
            self.disable_ghost_buttons()

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

        # update aggregation statistic field
        available_aggregation_statistics = ['Mean', 'Median', 'p1', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95', 'p99']
        self.cb_aggregation_statistic.addItems(available_aggregation_statistics)
        if self.selected_aggregation_statistic in available_aggregation_statistics:
            self.cb_aggregation_statistic.setCurrentText(self.selected_aggregation_statistic)
        else:
            self.selected_aggregation_statistic = self.cb_aggregation_statistic.currentText()

        # get available resampling resolutions
        available_resampling_resolutions = copy.deepcopy(available_resolutions)[available_resolutions.index(self.selected_resolution)+1:]
        available_resampling_resolutions.append('yearly')

        # remove resolutions if resampled data would be less than 2 timesteps
        resampling_resolutions = copy.deepcopy(available_resampling_resolutions)
        for resampling_resolution in resampling_resolutions:
            # get active frequency code
            active_frequency_code = get_frequency_code(resampling_resolution)

            # get time array
            start_date = self.le_start_date.text()
            end_date = self.le_end_date.text()
            time_array = pd.date_range(start=datetime.datetime(int(start_date[:4]), int(start_date[4:6]),
                                                               int(start_date[6:8])),
                                       end=datetime.datetime(int(end_date[:4]), int(end_date[4:6]), int(end_date[6:8])),
                                       freq=active_frequency_code)[:-1]

            # show warning when the data consists only of less than 2 timesteps
            if len(time_array) < 2:
                available_resampling_resolutions.remove(resampling_resolution)
            
        # update resampling resolution field
        available_resampling_resolutions = ['None',] + available_resampling_resolutions
        self.cb_resampling_resolution.addItems(available_resampling_resolutions)
        if self.selected_resampling_resolution in available_resampling_resolutions:
            self.cb_resampling_resolution.setCurrentText(self.selected_resampling_resolution)
        else:
            self.selected_resampling_resolution = self.cb_resampling_resolution.currentText()

        # update available experiments for selected fields
        aux.get_valid_experiments(self, self.le_start_date.text(), self.le_end_date.text(), self.selected_resolution,
                                  [self.selected_network], [self.selected_species])
        
        # update experiments -- keeping previously selected experiments if available
        if self.config_bar_initialisation:   
            self.experiments_menu['checkboxes']['keep_selected'] = [experiment for experiment in self.experiments
                                                                    if experiment in 
                                                                    self.experiments_menu['checkboxes']['map_vars']]
            self.experiments = {experiment:experiment_alias for experiment, experiment_alias in self.experiments.items()
                                if experiment in self.experiments_menu['checkboxes']['map_vars']}

        self.experiments_menu['checkboxes']['keep_selected'] = [previous_selected_experiment for
                                                                previous_selected_experiment in
                                                                self.experiments_menu['checkboxes']['keep_selected']
                                                                if previous_selected_experiment in
                                                                self.experiments_menu['checkboxes']['map_vars']]
        previous_experiments = self.experiments
        self.experiments = {exp:previous_experiments[exp] if exp in previous_experiments else exp 
                            for exp in self.experiments_menu['checkboxes']['keep_selected']}

        # update default qa
        default_qa = get_default_qa(self, self.selected_species)
        self.qa_menu['checkboxes']['remove_default'] = default_qa

        # update layout fields
        self.update_layout_fields(self.mpl_canvas)

        # unset variable to allow interactive handling from now
        self.block_config_bar_handling_updates = False

        # unset variable to allow updating the canvas
        self.block_MPL_canvas_updates = False

    def update_layout_fields(self, canvas_instance):
        """ Function which updates layout fields. """

        # set variable to block interactive handling while updating config bar parameters
        self.block_config_bar_handling_updates = True

        # update layout field buttons
        # clear fields
        self.cb_position_2.clear()
        self.cb_position_3.clear()
        self.cb_position_4.clear()
        self.cb_position_5.clear()

        # TODO: For Taylor diagrams, replace this piece of code for the one below when Matplotlib 3.7.2 is available
        # # remove plot types that need active temporal colocation and experiments data
        # for plot_type in ['scatter', 'taylor']:
        #     if ((not self.temporal_colocation) 
        #         or ((self.temporal_colocation) and (len(self.experiments) == 0))): 
        #         if plot_type in canvas_instance.layout_options:
        #             canvas_instance.layout_options.remove(plot_type)
        #     else:
        #         if plot_type not in canvas_instance.layout_options:
        #             canvas_instance.layout_options.append(plot_type)          

        # remove plot types that need active temporal colocation and experiments data
        if 'taylor' in canvas_instance.layout_options:
            canvas_instance.layout_options.remove('taylor')
        for plot_type in ['scatter']:
            if ((not self.temporal_colocation) 
                or ((self.temporal_colocation) and (len(self.experiments) == 0))): 
                if plot_type in canvas_instance.layout_options:
                    canvas_instance.layout_options.remove(plot_type)
            else:
                if plot_type not in canvas_instance.layout_options:
                    canvas_instance.layout_options.append(plot_type)          
 
        # order alphabetically
        layout_options = sorted(canvas_instance.layout_options)

        # update position 2 in layout
        self.cb_position_2.addItems(layout_options)
        if self.position_2 in layout_options:
            self.cb_position_2.setCurrentText(self.position_2)
        else:
            self.block_config_bar_handling_updates = False
            self.cb_position_2.setCurrentText('None')
            self.block_config_bar_handling_updates = True

        # update position 3 in layout
        self.cb_position_3.addItems(layout_options)
        if self.position_3 in layout_options:
            self.cb_position_3.setCurrentText(self.position_3)
        else:
            self.block_config_bar_handling_updates = False
            self.cb_position_3.setCurrentText('None')
            self.block_config_bar_handling_updates = True

        # update position 4 in layout
        self.cb_position_4.addItems(layout_options)
        if self.position_4 in layout_options:
            self.cb_position_4.setCurrentText(self.position_4)
        else:
            self.block_config_bar_handling_updates = False
            self.cb_position_4.setCurrentText('None')
            self.block_config_bar_handling_updates = True

        # update position 5 in layout
        self.cb_position_5.addItems(layout_options)
        if self.position_5 in layout_options:
            self.cb_position_5.setCurrentText(self.position_5)
        else:
            self.block_config_bar_handling_updates = False
            self.cb_position_5.setCurrentText('None')
            self.block_config_bar_handling_updates = True

        # unset variable to allow interactive handling from now
        self.block_config_bar_handling_updates = False

    def handle_config_bar_params_change(self, changed_param):
        """ Function which handles interactive updates of combo box fields. """

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
                    self.available_observation_data[self.selected_network][self.selected_resolution][self.selected_matrix].keys()))[0]
            
            elif event_source == self.cb_species:
                self.selected_species = changed_param

            elif event_source == self.cb_resampling_resolution:
                self.selected_resampling_resolution = changed_param
                self.selected_resampling = False

            # set variable to check if date range changes
            self.date_range_has_changed = False

            # check if start date/end date have changed
            if (event_source == self.le_start_date) or (event_source == self.le_end_date):
                self.date_range_has_changed = True

            # initalise multispecies tab if network, resolution, matrix or species change
            if event_source in [self.cb_network, self.cb_resolution, self.cb_matrix, self.cb_species]:
                aux.init_multispecies(self)

            # update configuration bar fields
            self.update_configuration_bar_fields()

    def handle_layout_update(self, changed_plot_type):
        """ Function which handles update of layout. """
        
        if (changed_plot_type != '') & (not self.block_config_bar_handling_updates):
            
            # get event origin source
            event_source = self.sender()

            # update selected station plots, avoiding duplicates
            if event_source == self.cb_position_2:
                previous_plot_type = copy.deepcopy(self.position_2)
                self.position_2 = copy.deepcopy(changed_plot_type)
                
            elif event_source == self.cb_position_3:
                previous_plot_type = copy.deepcopy(self.position_3)
                self.position_3 = copy.deepcopy(changed_plot_type)

            elif event_source == self.cb_position_4:
                previous_plot_type = copy.deepcopy(self.position_4)
                self.position_4 = copy.deepcopy(changed_plot_type)

            elif event_source == self.cb_position_5:
                previous_plot_type = copy.deepcopy(self.position_5)
                self.position_5 = copy.deepcopy(changed_plot_type)

            # if changed plot type is selected elsewhere, then set that field to None
            if (changed_plot_type == self.position_2) & (event_source != self.cb_position_2):
                self.position_2 = 'None'

            elif (changed_plot_type == self.position_3) & (event_source != self.cb_position_3):
                self.position_3 = 'None'

            if (changed_plot_type == self.position_4) & (event_source != self.cb_position_4):
                self.position_4 = 'None'

            elif (changed_plot_type == self.position_5) & (event_source != self.cb_position_5):
                self.position_5 = 'None'

            # remove axis elements for previous plot type, and from active_dashboard_plots
            if (previous_plot_type in self.active_dashboard_plots) & (previous_plot_type in self.mpl_canvas.plot_axes):
                ax = self.mpl_canvas.plot_axes[previous_plot_type]
                self.mpl_canvas.remove_axis_elements(ax, previous_plot_type)
                if isinstance(ax, dict):
                    for sub_ax in ax.values():
                        sub_ax.remove()
                else:
                    ax.remove()
                self.active_dashboard_plots.remove(previous_plot_type)

            # if changed_plot_type already axis on another axis then remove those axis elements
            if (changed_plot_type in self.active_dashboard_plots) & (changed_plot_type in self.mpl_canvas.plot_axes):
                ax = self.mpl_canvas.plot_axes[changed_plot_type]
                self.mpl_canvas.remove_axis_elements(ax, changed_plot_type)
                if isinstance(ax, dict):
                    for sub_ax in ax.values():
                        sub_ax.remove()
                else:
                    ax.remove()

            # otherwise add plot_type to active_dashboard_plots
            elif changed_plot_type != 'None': 
                self.active_dashboard_plots.append(changed_plot_type)

            # update plot axis for new plot type
            self.update_plot_axis(self.mpl_canvas, event_source, changed_plot_type)

            # update plot if changed_plot_type != None
            if changed_plot_type != 'None':

                # format axis                
                self.mpl_canvas.plot.format_axis(self.mpl_canvas.plot_axes[changed_plot_type], 
                                                 changed_plot_type, 
                                                 self.mpl_canvas.plot_characteristics[changed_plot_type])
                
                # make plot
                self.mpl_canvas.update_associated_active_dashboard_plot(changed_plot_type)

            # update layout fields
            self.update_layout_fields(self.mpl_canvas)

            # update buttons geometry
            self.update_buttons_geometry()
            
            # draw changes
            self.mpl_canvas.figure.canvas.draw()

        return None

    def update_plot_axis(self, canvas_instance, changed_position, changed_plot_type):
        """ Update plot axis position from layout options."""

        # since we have no data, we need to use a random reference to initialize the polar axes
        # the axis will be updated when we create the taylor diagram
        if changed_plot_type == 'taylor':
            reference_stddev = 7.5
            plot_characteristics = canvas_instance.plot_characteristics['taylor']
            ghelper = canvas_instance.plot.get_taylor_diagram_ghelper(reference_stddev, plot_characteristics)

        # position 2 (top right)
        if changed_position == self.cb_position_2 or changed_position == 2:
            if (changed_plot_type == 'periodic') or (changed_plot_type == 'periodic-violin'):
                canvas_instance.plot_axes[changed_plot_type] = {}
                canvas_instance.plot_axes[changed_plot_type]['hour'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((12, 50), rowspan=15, colspan=49))
                canvas_instance.plot_axes[changed_plot_type]['dayofweek'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((31, 81), rowspan=15, colspan=18))
                canvas_instance.plot_axes[changed_plot_type]['month'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((31, 50), rowspan=15, colspan=30))
            elif changed_plot_type == 'statsummary':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((12, 65), rowspan=34, colspan=34))
            elif changed_plot_type == 'taylor':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((9, 65), rowspan=36, colspan=18),
                                                                                                  axes_class=fa.FloatingAxes, grid_helper=ghelper)
            elif changed_plot_type != 'None':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((12, 50), rowspan=34, colspan=49))
            
        # position 3 (bottom left)
        if changed_position == self.cb_position_3 or changed_position == 3:
            if (changed_plot_type == 'periodic') or (changed_plot_type == 'periodic-violin'):
                canvas_instance.plot_axes[changed_plot_type] = {}
                canvas_instance.plot_axes[changed_plot_type]['hour'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 0), rowspan=18, colspan=29))
                canvas_instance.plot_axes[changed_plot_type]['dayofweek'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 19), rowspan=18, colspan=10))
                canvas_instance.plot_axes[changed_plot_type]['month'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 0), rowspan=18, colspan=18))
            elif changed_plot_type == 'statsummary':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 10), rowspan=38, colspan=19))
            elif changed_plot_type == 'taylor':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 5), rowspan=38, colspan=18),
                                                                                                  axes_class=fa.FloatingAxes, grid_helper=ghelper)
            elif changed_plot_type != 'None':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 0), rowspan=38, colspan=29))

        # position 4 (bottom centre)
        if changed_position == self.cb_position_4 or changed_position == 4:
            if (changed_plot_type == 'periodic') or (changed_plot_type == 'periodic-violin'):
                canvas_instance.plot_axes[changed_plot_type] = {}
                canvas_instance.plot_axes[changed_plot_type]['hour'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 35), rowspan=18, colspan=29))
                canvas_instance.plot_axes[changed_plot_type]['dayofweek'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 54), rowspan=18, colspan=10))
                canvas_instance.plot_axes[changed_plot_type]['month'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 35), rowspan=18, colspan=18))
            elif changed_plot_type == 'statsummary':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 45), rowspan=38, colspan=19))
            elif changed_plot_type == 'taylor':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 40), rowspan=38, colspan=18),
                                                                                                  axes_class=fa.FloatingAxes, grid_helper=ghelper)
            elif changed_plot_type != 'None':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 35), rowspan=38, colspan=29))
            
        # position 5 (bottom right)
        if changed_position == self.cb_position_5 or changed_position == 5:
            if (changed_plot_type == 'periodic') or (changed_plot_type == 'periodic-violin'):
                canvas_instance.plot_axes[changed_plot_type] = {}
                canvas_instance.plot_axes[changed_plot_type]['hour'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 70), rowspan=18, colspan=29))
                canvas_instance.plot_axes[changed_plot_type]['dayofweek'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 89), rowspan=18, colspan=10))
                canvas_instance.plot_axes[changed_plot_type]['month'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 70), rowspan=18, colspan=18))
            elif changed_plot_type == 'statsummary':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 80), rowspan=38, colspan=19))
            elif changed_plot_type == 'taylor':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 70), rowspan=38, colspan=18),
                                                                                                  axes_class=fa.FloatingAxes, grid_helper=ghelper)
            elif changed_plot_type != 'None':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 70), rowspan=38, colspan=29))

        # additional changes needed when defining the layout in the configuration changing active_dashboard_plots
        if changed_plot_type == 'timeseries':

            # setup timeseries annotation
            canvas_instance.create_timeseries_annotation()
            canvas_instance.create_timeseries_annotation_vline()
            canvas_instance.lock_timeseries_annotation = False
            canvas_instance.timeseries_annotation_event = canvas_instance.figure.canvas.mpl_connect('motion_notify_event', 
                                                                                                    canvas_instance.hover_timeseries_annotation)

        # additional changes needed when defining the layout in the configuration changing active_dashboard_plots
        elif changed_plot_type == 'scatter':

            # setup scatter annotation
            canvas_instance.create_scatter_annotation()
            canvas_instance.lock_scatter_annotation = False
            canvas_instance.scatter_annotation_event = canvas_instance.figure.canvas.mpl_connect('motion_notify_event', 
                                                                                                 canvas_instance.hover_scatter_annotation)

        # additional changes needed when defining the layout in the configuration changing active_dashboard_plots
        elif changed_plot_type == 'distribution':

            # setup distribution annotation
            canvas_instance.create_distribution_annotation()
            canvas_instance.create_distribution_annotation_vline()
            canvas_instance.lock_distribution_annotation = False
            canvas_instance.distribution_annotation_event = canvas_instance.figure.canvas.mpl_connect('motion_notify_event', 
                                                                                                      canvas_instance.hover_distribution_annotation)

        # additional changes needed when defining the layout in the configuration changing active_dashboard_plots
        elif changed_plot_type == 'periodic':
            
            # setup periodic annotation
            canvas_instance.create_periodic_annotation()
            canvas_instance.create_periodic_annotation_vline()
            canvas_instance.lock_periodic_annotation = dict()
            for resolution in canvas_instance.plot_axes['periodic'].keys():
                canvas_instance.lock_periodic_annotation[resolution] = False
            canvas_instance.periodic_annotation_event = canvas_instance.figure.canvas.mpl_connect('motion_notify_event', 
                                                                                                  canvas_instance.hover_periodic_annotation)

            # update periodic statistic combobox
            self.block_config_bar_handling_updates = False
            
            # update map z combobox fields based on data in memory
            # generate lists of basic and basis+bias statistics for using in the z statistic combobox
            if not hasattr(self, 'basic_z_stats'):
                self.basic_z_stats = np.array(list(
                    OrderedDict(sorted(basic_stats.items(), key=lambda x: x[1]['order'])).keys()))
            if not hasattr(self, 'basic_and_bias_z_stats'):
                self.basic_and_bias_z_stats = np.append(self.basic_z_stats, list(
                    OrderedDict(sorted(expbias_stats.items(), key=lambda x: x[1]['order'])).keys()))

            # update periodic statistic in dashboard
            canvas_instance.handle_periodic_statistic_update()

        # additional changes needed when defining the layout in the configuration changing active_dashboard_plots
        elif changed_plot_type == 'periodic-violin':
            
            # setup periodic violin annotation
            canvas_instance.create_periodic_violin_annotation()
            canvas_instance.create_periodic_violin_annotation_vline()
            canvas_instance.lock_periodic_violin_annotation = dict()
            for resolution in canvas_instance.plot_axes['periodic-violin'].keys():
                canvas_instance.lock_periodic_violin_annotation[resolution] = False
            canvas_instance.periodic_violin_annotation_event = canvas_instance.figure.canvas.mpl_connect('motion_notify_event', 
                                                                                                         canvas_instance.hover_periodic_violin_annotation)

        # additional changes needed when defining the layout in the configuration changing active_dashboard_plots
        elif changed_plot_type == 'taylor':
            
            # initialise polar axis
            canvas_instance.plot.taylor_polar_relevant_axis = canvas_instance.plot_axes[changed_plot_type].get_aux_axes(PolarAxes.PolarTransform())

            # setup taylor annotation
            canvas_instance.lock_taylor_annotation = False
            canvas_instance.taylor_annotation_event = canvas_instance.figure.canvas.mpl_connect('motion_notify_event', 
                                                                                                canvas_instance.hover_taylor_annotation)

    def handle_data_selection_update(self):
        """ Function which handles update of data selection
            and MPL canvas upon pressing of READ button.
        """

        # if have no data to read, then do not read any data
        if self.no_data_to_read:
            return

        # update mouse cursor to a waiting cursor
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        # clear previously selected relative/absolute station indices
        self.mpl_canvas.relative_selected_station_inds = np.array([], dtype=np.int)
        self.mpl_canvas.absolute_selected_station_inds = np.array([], dtype=np.int)
        self.mpl_canvas.absolute_non_selected_station_inds = np.array([], dtype=np.int)
        
        # clear and then hide all axes 
        for plot_type, ax in self.mpl_canvas.plot_axes.items():
            self.mpl_canvas.remove_axis_elements(ax, plot_type)
        
        # update MPL canvas
        self.mpl_canvas.figure.canvas.draw()  
        self.mpl_canvas.figure.canvas.flush_events()

        # set variable that blocks updating of MPL canvas until all data has been updated
        self.block_MPL_canvas_updates = True
        
        # set previous active variables
        self.previous_start_date = self.start_date
        self.previous_end_date = self.end_date
        self.previous_network = self.network
        self.previous_resolution = self.resolution
        self.previous_species = self.species
        self.previous_experiments = self.experiments
        self.previous_qa = self.qa
        self.previous_flags = self.flags
        self.previous_data_labels = self.data_labels
        self.previous_filter_species = self.previous_filter_species
        self.previous_plot_options = {}
        for plot_type in self.mpl_canvas.all_plots:
            self.previous_plot_options[plot_type] = []
        self.previous_statsummary_stats = {}
        self.previous_statsummary_stats['None'] = self.mpl_canvas.plot_characteristics['statsummary']['basic']
        for periodic_cycle in ['Diurnal', 'Weekly', 'Monthly']:
            self.previous_statsummary_stats[periodic_cycle] = []

        # set new active variables as selected variables from menu
        self.start_date = int(self.le_start_date.text())
        self.end_date = int(self.le_end_date.text())
        self.network = [self.selected_network]
        self.resolution = self.selected_resolution
        self.species = [self.selected_species]  
        self.experiments = {exp:self.previous_experiments[exp] if exp in self.previous_experiments else exp 
                            for exp in self.experiments_menu['checkboxes']['keep_selected']}
        self.qa = copy.deepcopy(self.qa_menu['checkboxes']['remove_selected'])
        self.qa_per_species[self.selected_species] = copy.deepcopy(self.qa)
        self.flags = copy.deepcopy(self.flag_menu['checkboxes']['remove_selected'])
        self.networkspecies = ['{}|{}'.format(network,speci) for network, speci in zip(self.network, self.species)]
        self.networkspeci = self.networkspecies[0]
        self.data_labels = ['observations'] + list(self.experiments.keys())
        self.current_plot_options = {}
        for plot_type in self.mpl_canvas.all_plots:
            self.current_plot_options[plot_type] = []
        self.current_statsummary_stats = {}
        for periodic_cycle in ['None', 'Diurnal', 'Weekly', 'Monthly']:
            self.current_statsummary_stats[periodic_cycle] = []

        # if spatial_colocation is not active, force filter_species to be empty dict if it is not already
        # inform user of this
        if (self.filter_species) and (not self.spatial_colocation):
            self.filter_species = {} 
            msg = '"spatial_colocation" must be set to True if wanting to use "filter_species" option.'
            aux.show_message(self.read_instance, msg)

        # set read operations to be empty list initially
        read_operations = []

        # if first read then need to read all data
        if self.first_read:
            read_operations = ['reset']

        # determine if any of the key variables have changed 
        # (network, resolution, species, qa, flags, filter_species)
        # if any have changed, observations and any selected experiments have to be re-read entirely
        elif (self.network[0] != self.previous_network[0]) or (
                self.resolution != self.previous_resolution) or (
                self.species[0] != self.previous_species[0]) or (
                np.array_equal(self.qa, self.previous_qa) == False) or (
                np.array_equal(self.flags, self.previous_flags) == False) or (
                self.filter_species != self.previous_filter_species):
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
            
            # inactivate resampling
            self.resampling = False

            # set current time array, as previous time array
            self.previous_time_array = self.time_array

            # set current station references, as previous station references
            self.previous_station_references = self.station_references
            
            # set current relevant yearmonths, as previous relevant yearmonths
            self.previous_yearmonths = self.yearmonths

            # read data
            self.datareader.read_setup(read_operations, experiments_to_remove=experiments_to_remove, 
                                       experiments_to_read=experiments_to_read)
            
            # restore mouse cursor to normal if have no valid data after read
            if self.invalid_read:
                QtWidgets.QApplication.restoreOverrideCursor()
                return

            # update fields available for filtering
            aux.update_representativity_fields(self)
            aux.update_period_fields(self)
            aux.update_metadata_fields(self)
            
            # update relevant temporal resolutions 
            self.relevant_temporal_resolutions = aux.get_relevant_temporal_resolutions(self.resolution)

        # if species has changed, or first read, update species specific lower/upper limits
        if (self.first_read) or (self.species[0] != self.previous_species[0]):
            # get default GHOST limits
            self.lower_bound[self.species[0]] = np.float32(self.parameter_dictionary[self.species[0]]['extreme_lower_limit']) 
            self.upper_bound[self.species[0]] = np.float32(self.parameter_dictionary[self.species[0]]['extreme_upper_limit']) 
            self.le_minimum_value.setText(str(self.lower_bound[self.species[0]]))
            self.le_maximum_value.setText(str(self.upper_bound[self.species[0]]))

        # run function to update filter
        self.mpl_canvas.handle_data_filter_update()
        
        # update map z combobox fields based on data in memory
        # generate lists of basic and basis+bias statistics for using in the z statistic combobox
        if not hasattr(self, 'basic_z_stats'):
            self.basic_z_stats = np.array(list(
                OrderedDict(sorted(basic_stats.items(), key=lambda x: x[1]['order'])).keys()))
        if not hasattr(self, 'basic_and_bias_z_stats'):
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

        # update map z statistic comboboxes
        self.mpl_canvas.handle_map_z_statistic_update()

        # update periodic statistic combobox
        self.mpl_canvas.handle_periodic_statistic_update()

        # update taylor diagram statistic combobox
        self.mpl_canvas.handle_taylor_correlation_statistic_update()
        
        # unselect all/intersect/extent checkboxes
        self.mpl_canvas.unselect_map_checkboxes()
        
        # reset resampling
        self.cb_resampling_resolution.setCurrentText('None')

        # unset variable to allow updating of MPL canvas
        self.block_MPL_canvas_updates = False

        # update MPL canvas
        self.mpl_canvas.update_MPL_canvas()

        # if first read, then set this now to be False
        # also, if colocate checkbox is ticked, then apply temporal colocation
        if self.first_read:
            self.first_read = False
            if self.ch_colocate.checkState() == QtCore.Qt.Checked:
                self.mpl_canvas.handle_temporal_colocate_update()

        # restore mouse cursor to normal
        QtWidgets.QApplication.restoreOverrideCursor()

    def reset_options(self):
        """ Reset all filter fields to initial values. """

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
        """ Disable button related only to ghost data. """
        
        # change background-color to indicate that it's nonusable
        self.bu_flags.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        self.bu_QA.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        self.bu_period.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        
        # disable buttons
        self.bu_flags.setEnabled(False)
        self.bu_QA.setEnabled(False)
        self.bu_period.setEnabled(False)
        
    def enable_ghost_buttons(self):
        """ Enable button related only to ghost data. """

        # enable buttons        
        self.bu_flags.setEnabled(True)
        self.bu_QA.setEnabled(True)
        self.bu_period.setEnabled(True)

# generate Providentia dashboard
def main(**kwargs):
    """ Main function. """
    
    q_app = QtWidgets.QApplication(sys.argv)
    q_app.setStyle("Fusion")
    ProvidentiaMainWindow(**kwargs)
    sys.exit(q_app.exec_())