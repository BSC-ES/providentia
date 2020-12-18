""" Module which provides main window """
from .configuration import ProvConfiguration
from .configuration import parse_path
from .reading import read_netcdf_data
from .reading import read_netcdf_nonghost
from .reading import get_yearmonths_to_read
from .prov_canvas import MPLCanvas
from .toolbar import NavigationToolbar
from .toolbar import save_data
from .prov_dashboard_aux import ComboBox
from .prov_dashboard_aux import QVLine
from .prov_dashboard_aux import PopUpWindow
from .prov_dashboard_aux import formatting_dict
from .prov_dashboard_aux import set_formatting

import copy
import bisect
import datetime
import gc
import multiprocessing
import os
import os.path
import json
import sys
from functools import partial
from collections import OrderedDict

from PyQt5 import QtCore, QtWidgets, QtGui
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import seaborn as sns
from dateutil.relativedelta import relativedelta

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class ProvidentiaMainWindow(QtWidgets.QWidget, ProvConfiguration):
    """Define class that generates Providentia dashboard"""

    # create signals that are fired upon resizing/moving of main Providentia window
    resized = QtCore.pyqtSignal()
    move = QtCore.pyqtSignal()

    def __init__(self, read_type='parallel', **kwargs):
        super(ProvidentiaMainWindow, self).__init__()
        ProvConfiguration.__init__(self, **kwargs)

        # put read_type into self
        self.read_type = read_type

        # store options to be restored at the end
        # self.localvars = copy.deepcopy(vars(self))
        dconf_path = (os.path.join(CURRENT_PATH, 'conf/default.conf'))
        # update from config file
        if ('config' in kwargs) and ('section' in kwargs):
            self.load_conf(kwargs['section'], kwargs['config'])
            self.from_conf = True
        elif os.path.isfile(dconf_path):
            self.load_conf('default', dconf_path)
            self.from_conf = True
        # update from command line
        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in kwargs.items()})
        # arguments are only local

        self.main_window_geometry = None

        self.parameter_dictionary = {}
        self.standard_metadata = {}
        self.metadata_vars_to_read = {}
        self.metadata_dtype = {}
        self.standard_data_flag_name_to_data_flag_code = {}
        self.standard_QA_name_to_QA_code = {}
        self.init_standards()

        # load necessary dictionaries
        self.basic_stats_dict = json.load(open(os.path.join(CURRENT_PATH,
                                                            'conf/basic_stats_dict.json')))
        self.expbias_dict = json.load(open(os.path.join(CURRENT_PATH,
                                                        'conf/experiment_bias_stats_dict.json')))
        # create UI
        self.init_ui()

        # setup callback events upon resizing/moving of Providentia window
        self.resized.connect(self.get_geometry)
        self.move.connect(self.get_geometry)

    def init_standards(self):
        """ Read from ghost standards """
        sys.path.insert(1, '{}/GHOST_standards/{}'.format(self.obs_root, self.ghost_version))
        from GHOST_standards import standard_parameters, \
            get_standard_metadata, standard_data_flag_name_to_data_flag_code, \
            standard_QA_name_to_QA_code
        # modify standard parameter dictionary to have BSC standard parameter names as
        # keys (rather than GHOST)
        for _, param_dict in standard_parameters.items():
            self.parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict
        # get standard metadata dictionary
        self.standard_metadata = get_standard_metadata({'standard_units':''})
        # create list of metadata variables to read (make global)
        self.metadata_vars_to_read = [key for key in self.standard_metadata.keys()
                                      if pd.isnull(self.standard_metadata[key]['metadata_type'])
                                      == False]
        self.metadata_dtype = [(key, self.standard_metadata[key]['data_type']) for key in
                               self.metadata_vars_to_read]
        self.standard_data_flag_name_to_data_flag_code = \
            standard_data_flag_name_to_data_flag_code
        self.standard_QA_name_to_QA_code = standard_QA_name_to_QA_code
        self.qa_exceptions = ['dir10', 'spd10', 'rho2', 'acprec', 'acsnow', 'si',
                              'cldbot', 'vdist', 'ccovmean', 'cfracmean']
        self.get_qa_codes()

    def get_qa_codes(self):
        """Retrieve QA codes from GHOST_standards using the qa flags' names.

        Specific flags are defined for the following species:
        ['WND_DIR_10M','WND_SPD_10M','RH_2M','PREC_ACCUM','SNOW_ACCUM',
        'SNOW_DEPTH','CEILING_HEIGHT','VIS_DIST','CLOUD_CVG','CLOUD_CVG_FRAC']"""

        # get names from json files
        specific_qa_names = json.load(open("providentia/conf/default_flags.json"))['specific_qa']
        general_qa_names = json.load(open("providentia/conf/default_flags.json"))['general_qa']
        # get codes
        self.specific_qa = [self.standard_QA_name_to_QA_code[qa_name] for qa_name in specific_qa_names]
        self.general_qa = [self.standard_QA_name_to_QA_code[qa_name] for qa_name in general_qa_names]
        # get difference of flags, needed later for updating default selection
        self.qa_diff = list(set(self.general_qa) - set(self.specific_qa))

    def which_qa(self):
        """Checks if the species we currently have selected belongs to the ones
        that have specific qa flags selected as default"""

        if hasattr(self, 'qa'):
            # return subset the user has selected in conf
            return eval(self.qa)
        if self.selected_species in self.qa_exceptions:
            return self.specific_qa
        else:
            return self.general_qa

    def which_flags(self):
        """if there are flags coming from a config file, select those"""

        if hasattr(self, 'flags'):
            return eval(self.flags)
        else:
            return []

    def which_bounds(self):
        """if there are bounds defined in a config file, fill that value,
        if it is withing the feasible bounds of the species"""

        lower = np.float32(self.parameter_dictionary[self.active_species]['extreme_lower_limit'])
        upper = np.float32(self.parameter_dictionary[self.active_species]['extreme_upper_limit'])

        if hasattr(self, 'lower_bound'):
            if self.lower_bound >= lower:
                lower = self.lower_bound

        if hasattr(self, 'upper_bound'):
            if self.upper_bound <= upper:
                upper = self.upper_bound

        return np.float32(lower), np.float32(upper)

    def resizeEvent(self, event):
        '''Function to overwrite default PyQt5 resizeEvent function --> for calling get_geometry'''
        self.resized.emit()
        return super(ProvidentiaMainWindow, self).resizeEvent(event)

    def moveEvent(self, event):
        '''Function to overwrite default PyQt5 moveEvent function --> for calling get_geometry'''
        self.move.emit()
        return super(ProvidentiaMainWindow, self).moveEvent(event)

    def get_geometry(self):
        '''Get current geometry of main Providentia window'''
        self.main_window_geometry = copy.deepcopy(self.geometry())

    def init_ui(self):
        """Initialise user interface"""

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
        self.setStyleSheet("QToolTip { font: %spt %s}"%(formatting_dict['tooltip']['font'].pointSizeF(), formatting_dict['tooltip']['font'].family()))

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
        self.lb_data_selection = set_formatting(QtWidgets.QLabel(self, text="Data Selection"), formatting_dict['title_menu'])
        self.lb_data_selection.setToolTip('Setup configuration of data to read into memory')
        self.bu_read = set_formatting(QtWidgets.QPushButton('READ', self), formatting_dict['button_menu'])
        self.bu_read.setFixedWidth(40)
        self.bu_read.setStyleSheet("color: green;")
        self.bu_read.setToolTip('Read selected configuration of data into memory')
        self.ch_colocate = set_formatting(QtWidgets.QCheckBox("Colocate"), formatting_dict['checkbox_menu'])
        self.ch_colocate.setToolTip('Temporally colocate observational/experiment data')
        self.cb_network = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_network.setFixedWidth(95)
        self.cb_network.setToolTip('Select providing observational data network. Names starting with * indicate non-GHOST datasets')
        self.cb_resolution = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_resolution.setFixedWidth(95)
        self.cb_resolution.setToolTip('Select temporal resolution of data')
        self.cb_matrix = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_matrix.setFixedWidth(136)
        self.cb_matrix.setToolTip('Select data matrix')
        self.cb_species = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_species.setStyleSheet("QComboBox { combobox-popup: 0; }")
        self.cb_species.setFixedWidth(136)
        self.cb_species.setToolTip('Select species')
        self.cb_species.view().setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
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
        self.lb_data_filter = set_formatting(QtWidgets.QLabel(self, text="Filters"), formatting_dict['title_menu'])
        self.lb_data_filter.setFixedWidth(65)
        self.lb_data_filter.setToolTip('Select criteria to filter data by')
        self.bu_rep = set_formatting(QtWidgets.QPushButton('% REP', self), formatting_dict['button_menu'])
        self.bu_rep.setFixedWidth(46)
        self.bu_rep.setToolTip('Select % desired representativity in data across whole record and for specific temporal periods')
        self.bu_meta = set_formatting(QtWidgets.QPushButton('META', self), formatting_dict['button_menu'])
        self.bu_meta.setFixedWidth(46)
        self.bu_meta.setToolTip('Select metadata to filter by')
        self.bu_reset = set_formatting(QtWidgets.QPushButton('RESET', self), formatting_dict['button_menu'])
        self.bu_reset.setFixedWidth(50)
        self.bu_reset.setToolTip('Reset filter fields to initial values')
        self.bu_reset.setStyleSheet("color: red;")
        self.bu_period = set_formatting(QtWidgets.QPushButton('PERIOD', self), formatting_dict['button_menu'])
        self.bu_period.setFixedWidth(50)
        self.bu_period.setToolTip('Select data in specific periods')
        self.bu_screen = set_formatting(QtWidgets.QPushButton('FILTER', self), formatting_dict['button_menu'])
        self.bu_screen.setFixedWidth(50)
        self.bu_screen.setStyleSheet("color: blue;")
        self.bu_screen.setToolTip('Filter data')
        self.lb_data_bounds = set_formatting(QtWidgets.QLabel(self, text="Bounds"), formatting_dict['label_menu'])
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
        self.lb_z = set_formatting(QtWidgets.QLabel(self, text="Map Z"), formatting_dict['title_menu'])
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
        self.lb_experiment_bias = set_formatting(QtWidgets.QLabel(self, text="Exp. Bias"), formatting_dict['title_menu'])
        self.lb_experiment_bias.setToolTip('Set experiment bias statistic')
        self.cb_experiment_bias_type = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_experiment_bias_type.setFixedWidth(100)
        self.cb_experiment_bias_type.setToolTip('Select experiment bias type')
        self.cb_experiment_bias_stat = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_experiment_bias_stat.setFixedWidth(100)
        self.cb_experiment_bias_stat.setToolTip('Select experiment bias statistic')
        self.vertical_splitter_4 = QVLine()
        self.vertical_splitter_4.setMaximumWidth(20)
        self.lb_station_selection = set_formatting(QtWidgets.QLabel(self, text="Site Select"), formatting_dict['title_menu'])
        self.lb_station_selection.setToolTip('Select stations')
        self.ch_select_all = set_formatting(QtWidgets.QCheckBox("All"), formatting_dict['checkbox_menu'])
        self.ch_select_all.setToolTip('Select all stations')
        self.ch_intersect = set_formatting(QtWidgets.QCheckBox("Intersect"), formatting_dict['checkbox_menu'])
        self.ch_intersect.setToolTip('Select stations that intersect with all loaded model domains')

        # position objects on gridded configuration bar
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
        config_bar.addWidget(self.bu_reset, 0, 7)
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

        # enable dynamic updating of configuration bar fields which filter data files
        self.cb_network.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.cb_resolution.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.cb_matrix.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.cb_species.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.le_start_date.textChanged.connect(self.config_bar_params_change_handler)
        self.le_end_date.textChanged.connect(self.config_bar_params_change_handler)

        #setup pop-up window menu tree for flags
        self.flag_menu = {'window_title':'FLAGS', 'page_title':'Select standardised data reporter provided flags to filter by', 'checkboxes':{}}
        self.flag_menu['checkboxes']['labels'] = np.array(sorted(self.standard_data_flag_name_to_data_flag_code, key=self.standard_data_flag_name_to_data_flag_code.get))
        self.flag_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
        self.flag_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
        self.flag_menu['checkboxes']['map_vars'] = np.sort(list(self.standard_data_flag_name_to_data_flag_code.values()))
        self.flag_menu['select_buttons'] = ['all', 'clear', 'default']

        #setup pop-up window menu tree for qa
        self.qa_menu = {'window_title':'QA', 'page_title':'Select standardised quality assurance flags to filter by', 'checkboxes':{}}
        self.qa_menu['checkboxes']['labels'] = np.array(sorted(self.standard_QA_name_to_QA_code, key=self.standard_QA_name_to_QA_code.get))
        self.qa_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
        self.qa_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
        self.qa_menu['checkboxes']['map_vars'] = np.sort(list(self.standard_QA_name_to_QA_code.values()))
        self.qa_menu['select_buttons'] = ['all', 'clear', 'default']

        #setup pop-up window menu tree for experiments
        self.experiments_menu = {'window_title':'EXPERIMENTS', 'page_title':'Select Experiment/s', 'checkboxes':{}}
        self.experiments_menu['checkboxes']['labels'] = []
        self.experiments_menu['checkboxes']['keep_default'] = []
        self.experiments_menu['checkboxes']['keep_selected'] = []
        self.experiments_menu['checkboxes']['map_vars'] = []
        self.experiments_menu['select_buttons'] = ['all', 'clear']

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
            self.metadata_menu[metadata_type] = {'window_title':metadata_type, 'page_title':self.metadata_menu['navigation_buttons']['tooltips'][metadata_type_ii], 'navigation_buttons':{}, 'rangeboxes':{}}
            self.metadata_menu[metadata_type]['navigation_buttons']['labels'] = [metadata_name for metadata_name in self.standard_metadata.keys() if (self.standard_metadata[metadata_name]['metadata_type'] == metadata_type) & (self.standard_metadata[metadata_name]['data_type'] == np.object)]
            self.metadata_menu[metadata_type]['navigation_buttons']['tooltips'] = [self.standard_metadata[metadata_name]['description'] for metadata_name in self.metadata_menu[metadata_type]['navigation_buttons']['labels']]
            for label in self.metadata_menu[metadata_type]['navigation_buttons']['labels']:
                self.metadata_menu[metadata_type][label] = {'window_title':label, 'page_title':'Filter stations by unique {} metadata'.format(label), 'checkboxes':{}}
                self.metadata_menu[metadata_type][label]['checkboxes']['labels'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['keep_selected'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['remove_selected'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['keep_default'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['remove_default'] = []
            self.metadata_menu[metadata_type]['rangeboxes']['labels'] = [metadata_name for metadata_name in self.standard_metadata.keys() if (self.standard_metadata[metadata_name]['metadata_type'] == metadata_type) & (self.standard_metadata[metadata_name]['data_type'] != np.object)]
            self.metadata_menu[metadata_type]['rangeboxes']['tooltips'] = [self.standard_metadata[metadata_name]['description'] for metadata_name in self.metadata_menu[metadata_type]['rangeboxes']['labels']]
            self.metadata_menu[metadata_type]['rangeboxes']['current_lower'] = ['nan'] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels'])
            self.metadata_menu[metadata_type]['rangeboxes']['current_upper'] = ['nan'] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels'])
            self.metadata_menu[metadata_type]['rangeboxes']['lower_default'] = ['nan'] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels'])
            self.metadata_menu[metadata_type]['rangeboxes']['upper_default'] = ['nan'] * len(self.metadata_menu[metadata_type]['rangeboxes']['labels'])

        #setup pop-up window menu tree for % data representativity
        self.representativity_menu = {'window_title':'% DATA REPRESENTATIVITY', 'page_title':'Select % Data Representativity Bounds', 'rangeboxes':{}}
        self.representativity_menu['rangeboxes']['labels'] = []
        self.representativity_menu['rangeboxes']['tooltips'] = []
        self.representativity_menu['rangeboxes']['current_lower'] = []

        # setup pop-up window menu tree for data periods
        self.period_menu = {'window_title':'DATA PERIOD', 'page_title':'Select Data Periods', 'checkboxes':{}}
        self.period_menu['checkboxes']['labels'] = []
        self.period_menu['checkboxes']['keep_selected'] = []
        self.period_menu['checkboxes']['remove_selected'] = []

        # enable pop up configuration windows
        self.bu_flags.clicked.connect(partial(self.generate_pop_up_window, self.flag_menu))
        self.bu_QA.clicked.connect(partial(self.generate_pop_up_window, self.qa_menu))
        self.bu_experiments.clicked.connect(partial(self.generate_pop_up_window, self.experiments_menu))
        self.bu_meta.clicked.connect(partial(self.generate_pop_up_window, self.metadata_menu))
        self.bu_rep.clicked.connect(partial(self.generate_pop_up_window, self.representativity_menu))
        self.bu_period.clicked.connect(partial(self.generate_pop_up_window, self.period_menu))

        # initialise configuration bar fields
        self.config_bar_initialisation = True
        self.update_configuration_bar_fields()
        self.config_bar_initialisation = False

        # Setup MPL canvas of plots
        # set variable that blocks updating of MPL canvas until some data has been read
        self.block_MPL_canvas_updates = True
        self.mpl_canvas = MPLCanvas(self)

        # Enable interactivity of functions which update MPL canvas
        # enable READ button
        self.bu_read.clicked.connect(self.handle_data_selection_update)

        # enable RESET button
        self.bu_reset.clicked.connect(self.reset_options)

        # enable interactivity of colocation checkbox
        self.ch_colocate.stateChanged.connect(self.mpl_canvas.handle_colocate_update)

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

        # Generate MPL navigation toolbar
        self.navi_toolbar = NavigationToolbar(self.mpl_canvas, self)

        # add more buttons on the toolbar, next to the navi_toolbar
        self.savebutton = QtWidgets.QPushButton()
        self.savebutton.setFlat(True)
        self.savebutton.setToolTip("Save current instance of data and metadata")
        self.savebutton.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/conf_icon.png")))
        self.savebutton.setIconSize(QtCore.QSize(31, 35))
        self.savebutton.setStyleSheet("QPushButton { border: none;} QPushButton:hover "
                                      "{ border-width: 1px; border-style: solid; border-color: darkgrey; "
                                      "border-radius: 4px; background-color : white; }")
        self.savebutton.clicked.connect(self.savebutton_func)

        # add more buttons on the toolbar, next to the navi_toolbar
        self.conf_load = QtWidgets.QPushButton()
        self.conf_load.setFlat(True)
        self.conf_load.setToolTip("Load toolbar selections from configuration file")
        self.conf_load.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_data.png")))
        self.conf_load.setIconSize(QtCore.QSize(31, 35))
        self.conf_load.setStyleSheet("QPushButton { border: none;} QPushButton:hover "
                                      "{ border-width: 1px; border-style: solid; border-color: darkgrey; "
                                      "border-radius: 4px; background-color : white; }")
        self.savebutton.clicked.connect(self.savebutton_func)

        # position config bar, navigation toolbar and MPL canvas and elements in parent layout
        hbox.addWidget(self.savebutton)
        hbox.addWidget(self.conf_load)
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

    def load_conf(self):
        """if user has specified a conf file on startup, load
        the session according to that configuration"""

    def savebutton_func(self):
        save_data(self.mpl_canvas)

    def generate_pop_up_window(self, menu_root):
        '''generate pop up window'''

        self.pop_up_window = PopUpWindow(menu_root, [], self.main_window_geometry)

    def update_configuration_bar_fields(self):
        """Define function that initialises/updates configuration bar fields"""

        # set variable to block interactive handling while updating config bar parameters
        self.block_config_bar_handling_updates = True

        self.check_for_ghost()

        # set some default configuration values when initialising config bar
        if self.config_bar_initialisation:
            # set initially selected/active start-end date as default 201601-201701
            self.le_start_date.setText(self.start_date)
            self.le_end_date.setText(self.end_date)
            self.selected_start_date = int(self.le_start_date.text())
            self.selected_end_date = int(self.le_end_date.text())
            self.selected_start_date_firstdayofmonth = int(str(self.selected_start_date)[:6]+'01')
            self.active_start_date = int(self.le_start_date.text())
            self.active_end_date = int(self.le_end_date.text())
            self.date_range_has_changed = False

            # set selected/active values of other fields to be initially None
            # TODO: no need to initialize to None, as we're initializing from conf
            # self.selected_network = None
            self.active_network = None
            self.active_resolution = None
            self.active_matrix = None
            self.active_species = None

            # set selected/active values of variables associated
            # with pop up windows to be empty lists
            self.active_experiment_grids = []
            self.active_qa = []
            self.active_flags = []

            # set initial time array to be None
            self.time_array = None

            self.relevant_yearmonths = None

            # set initial station references to be empty list
            self.station_references = []

            # Gather all observational data
            # create nested dictionary storing all observational species
            # data by species matrix, by temporal resolution, by network,
            # associated with list of start YYYYMM yearmonths of data files

            self.all_observation_data = {}

            # set all available networks
            available_networks = eval(self.available_networks)

            # set all available temporal resolutions
            available_resolutions = ['hourly', '3hourly', '6hourly', 'hourly_instantaneous',
                                     '3hourly_instantaneous', '6hourly_instantaneous',
                                     'daily', 'monthly']

            # iterate through available networks
            for network in available_networks:

                # check if directory for network exists
                # if not, continue
                if not os.path.exists('%s/%s/%s' % (self.obs_root, network, self.ghost_version)):
                    continue

                # write empty dictionary for network
                self.all_observation_data[network] = {}

                # iterate through available resolutions
                for resolution in available_resolutions:

                    # check if directory for resolution exists
                    # if not, continue
                    if not os.path.exists('%s/%s/%s/%s' % (self.obs_root, network, self.ghost_version, resolution)):
                        continue

                    # write nested empty dictionary for resolution
                    self.all_observation_data[network][resolution] = {}

                    # get available species for network/resolution
                    available_species = os.listdir('%s/%s/%s/%s' % (self.obs_root, network, self.ghost_version, resolution))

                    # iterate through available files per species
                    for species in available_species:

                        # get all netCDF monthly files per species
                        species_files = os.listdir(
                            '%s/%s/%s/%s/%s' % (self.obs_root, network, self.ghost_version, resolution, species))

                        # get monthly start date (YYYYMM) of all species files
                        species_files_yearmonths = \
                            [int(f.split('_')[-1][:6]+'01') for f in species_files if f != 'temporary']

                        # get matrix for current species
                        matrix = self.parameter_dictionary[species]['matrix']

                        if matrix not in list(self.all_observation_data[network][resolution].keys()):
                            # write nested empty dictionary for matrix
                            self.all_observation_data[network][resolution][matrix] = {}

                        # write nested dictionary for species, with associated file yearmonths
                        self.all_observation_data[network][resolution][matrix][species] = species_files_yearmonths

            # load dictionary with esarchive files
            esarchive_files = json.load(open(os.path.join(CURRENT_PATH, 'conf/esarchive_files.json')))
            # and merge to existing dict if we have the path
            if self.nonghost_root is not None:
                self.all_observation_data = {**self.all_observation_data, **esarchive_files}
            # create dictionary of observational data inside date range
            self.get_valid_obs_files_in_date_range()

            # check which flags to select, depending if we have conf file or no
            self.flag_menu['checkboxes']['remove_selected'] = self.which_flags()

        # if date range has changed then update available observational data dictionary
        if self.date_range_has_changed:
            self.get_valid_obs_files_in_date_range()

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
        available_resolutions = list(self.available_observation_data[self.cb_network.currentText()].keys())
        # manually force order of available resolutions
        resolution_order_dict = {'hourly':1, '3hourly':2, '6hourly':3, 'hourly_instantaneous':4,
                                 '3hourly_instantaneous':5, '6hourly_instantaneous':6,
                                 'daily':7, 'monthly':8}
        available_resolutions = sorted(available_resolutions, key=resolution_order_dict.__getitem__)
        self.cb_resolution.addItems(available_resolutions)
        if self.selected_resolution in available_resolutions:
            self.cb_resolution.setCurrentText(self.selected_resolution)
        else:
            self.selected_resolution = self.cb_resolution.currentText()

        # update matrix field
        available_matrices = sorted(
            self.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()])
        self.cb_matrix.addItems(available_matrices)
        if self.selected_matrix in available_matrices:
            self.cb_matrix.setCurrentText(self.selected_matrix)
        else:
            self.selected_matrix = self.cb_matrix.currentText()

        # update species field
        available_species = sorted(self.available_observation_data[self.cb_network.currentText()][
                                       self.cb_resolution.currentText()][self.cb_matrix.currentText()])
        self.cb_species.addItems(available_species)
        if self.selected_species in available_species:
            self.cb_species.setCurrentText(self.selected_species)
        else:
            self.selected_species = self.cb_species.currentText()

        # update available experiment data dictionary
        self.get_valid_experiment_files_in_date_range()
        # update selected indices for experiments -- keeping previously selected experiments if available
        # set selected indices as previously selected indices in current available list of experiments
        if self.config_bar_initialisation and hasattr(self, 'experiments'):
            self.experiments_menu['checkboxes']['keep_selected'] = [experiment for experiment in eval(self.experiments)
                                                                    if experiment in
                                                                    self.experiments_menu['checkboxes']['map_vars']]
        self.experiments_menu['checkboxes']['keep_selected'] = [previous_selected_experiment for
                                                                previous_selected_experiment in
                                                                self.experiments_menu['checkboxes']['keep_selected']
                                                                if previous_selected_experiment in
                                                                self.experiments_menu['checkboxes']['map_vars']]

        # since a selection has changed, update also the qa flags
        qa_to_select = self.which_qa()  # first check which flags
        self.qa_menu['checkboxes']['remove_default'] = qa_to_select
        if self.config_bar_initialisation:
            self.qa_menu['checkboxes']['remove_selected'] = qa_to_select
        else:
            # if the selected species has specific qa flags, ensure that none of the
            # inapplicable is selected
            if self.selected_species in self.qa_exceptions:
                self.qa_menu['checkboxes']['remove_selected'] = list(set(
                    self.qa_menu['checkboxes']['remove_selected']) - set(self.qa_diff))

        # if self.config_bar_initialisation:
        #     self.flag_menu['checkboxes']['remove_selected'] = self.which_flags()

        # unset variable to allow interactive handling from now
        self.block_config_bar_handling_updates = False

    def get_valid_obs_files_in_date_range(self):
        """Define function that iterates through observational dictionary tree
        and returns a dictionary of available data in the selected date
        range"""

        # create dictionary to store available observational data
        self.available_observation_data = {}

        # check if start/end date are valid values, if not, return with no valid obs. files
        selected_start_date = self.le_start_date.text()
        selected_end_date = self.le_end_date.text()
        if (self.valid_date(selected_start_date)) & (self.valid_date(selected_end_date)):
            self.date_range_has_changed = True
            self.selected_start_date = int(selected_start_date)
            self.selected_end_date = int(selected_end_date)
            self.selected_start_date_firstdayofmonth = int(str(self.selected_start_date)[:6] + '01')
        else:
            return

        # check end date is > start date, if not, return with no valid obs. files
        if self.selected_start_date >= self.selected_end_date:
            return

        # check start date and end date are both within if valid date range (19000101 - 20500101),
        # if not, return with no valid obs. files
        if (self.selected_start_date < 19000101) or (self.selected_end_date < 19000101) or (
                self.selected_start_date >= 20500101) or (self.selected_end_date >= 20500101):
            return

        # iterate through networks
        for network in list(self.all_observation_data.keys()):
            # iterate through resolutions
            for resolution in list(self.all_observation_data[network].keys()):
                # iterate through matrices
                for matrix in list(self.all_observation_data[network][resolution].keys()):
                    # iterate through species
                    for species in list(self.all_observation_data[network][resolution][matrix].keys()):
                        # get all file yearmonths associated with species
                        species_file_yearmonths = self.all_observation_data[network][resolution][matrix][species]
                        # get file yearmonths within date range
                        valid_species_files_yearmonths = [ym for ym in species_file_yearmonths if
                                                          (ym >= self.selected_start_date_firstdayofmonth) & (
                                                                      ym < self.selected_end_date)]
                        if len(valid_species_files_yearmonths) > 0:
                            # if network not in dictionary yet, add it
                            if network not in list(self.available_observation_data.keys()):
                                self.available_observation_data[network] = {}
                            # if resolution not in dictionary yet, add it
                            if resolution not in list(self.available_observation_data[network].keys()):
                                self.available_observation_data[network][resolution] = {}
                            # if matrix not in dictionary yet, add it
                            if matrix not in list(self.available_observation_data[network][resolution].keys()):
                                self.available_observation_data[network][resolution][matrix] = {}
                            # add species with associated list of file start yearmonths
                            self.available_observation_data[network][resolution][matrix][
                                species] = valid_species_files_yearmonths

    def get_valid_experiment_files_in_date_range(self):
        """Define function which gathers available experiment
        data for selected network/resolution/species.
        A dictionary is created storing available experiment-grid
        names associated with valid files in set date range.
        """

        # create dictionary to store available experiment information
        self.available_experiment_data = {}

        # get all different experiment names
        available_experiments = os.listdir('%s/%s' % (self.exp_root, self.ghost_version))

        # iterate through available experiments
        for experiment in available_experiments:

            # get all available grid types by experiment
            available_grids = os.listdir('%s/%s/%s' % (self.exp_root, self.ghost_version, experiment))

            # iterate through all available grids
            for grid in available_grids:

                # get all available ensemble member numbers
                available_member_numbers = os.listdir('%s/%s/%s/%s' % (self.exp_root, self.ghost_version, experiment, grid))

                # iterate through all available ensemble member numbers
                for member in available_member_numbers:

                    # test first if interpolated directory exists before trying to get files from it
                    # if it does not exit, continue
                    if not os.path.exists(
                            '%s/%s/%s/%s/%s/%s/%s/%s' % (self.exp_root, self.ghost_version, experiment, grid, member,
                                                      self.selected_resolution, self.selected_species,
                                                      self.selected_network)):
                        continue
                    else:
                        # get all experiment netCDF files by experiment/grid/selected
                        # resolution/selected species/selected network
                        network_files = os.listdir(
                            '%s/%s/%s/%s/%s/%s/%s/%s' % (self.exp_root, self.ghost_version,
                                                         experiment, grid, member, self.selected_resolution,
                                                         self.selected_species, self.selected_network))
                        # get start YYYYMM yearmonths of data files
                        network_files_yearmonths = [int(f.split('_')[-1][:6]+'01') for f in network_files]
                        # limit data files to just those within date range
                        valid_network_files_yearmonths = \
                            [ym for ym in network_files_yearmonths if (ym >= self.selected_start_date_firstdayofmonth) &
                             (ym < self.selected_end_date)]

                        # if have some valid data files for experiment-grid-member, add experiment-grid-member key
                        # (with associated yearmonths) to dictionary
                        if len(valid_network_files_yearmonths) > 0:
                            self.available_experiment_data['%s-%s-%s' % (experiment, grid, member)] = valid_network_files_yearmonths

        # get list of available experiment-grid names
        self.experiments_menu['checkboxes']['labels'] = np.array(
            sorted(list(self.available_experiment_data.keys())))
        self.experiments_menu['checkboxes']['map_vars'] = copy.deepcopy(
            self.experiments_menu['checkboxes']['labels'])

    def config_bar_params_change_handler(self, changed_param):
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
                    self.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()][
                        self.cb_matrix.currentText()].keys()))[0]
            elif event_source == self.cb_species:
                self.selected_species = changed_param

            # set variable to check if date range changes
            self.date_range_has_changed = False

            # check if start date/end date have changed
            if (event_source == self.le_start_date) or (event_source == self.le_end_date):
                self.date_range_has_changed = True

            # update configuration bar fields
            self.update_configuration_bar_fields()

            # if we're reading nonghost files, then disable fields again
            self.check_for_ghost()

    def check_for_ghost(self):
        """ It checks whether the selected network comes from GHOST or not.
        In case of non-ghost, it disables ghost-related fields"""

        # if we're reading nonghost files, then disable fields
        if '*' in self.cb_network.currentText():
            self.disable_ghost_buttons()
            self.reading_nonghost = True
        else:
            self.reading_nonghost = False
            self.enable_ghost_buttons()

    def valid_date(self, date_text):
        """define function that determines if a date string is in the correct format"""

        try:
            datetime.datetime.strptime(date_text, '%Y%m%d')
            return True
        except:
            return False

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

        # set all previously active variables as past_active variables
        self.previous_active_network = self.active_network
        self.previous_active_resolution = self.active_resolution
        self.previous_active_matrix = self.active_matrix
        self.previous_active_species = self.active_species
        self.previous_active_start_date = self.active_start_date
        self.previous_active_end_date = self.active_end_date
        self.previous_active_experiment_grids = self.active_experiment_grids
        self.previous_active_qa = self.active_qa
        self.previous_active_flags = self.active_flags

        # set all currently selected variables as active variables
        self.active_network = self.selected_network
        self.active_resolution = self.selected_resolution
        self.active_matrix = self.selected_matrix
        self.active_species = self.selected_species
        self.active_start_date = self.selected_start_date
        self.active_end_date = self.selected_end_date
        self.active_experiment_grids = copy.deepcopy(self.experiments_menu['checkboxes']['keep_selected'])
        self.active_qa = copy.deepcopy(self.qa_menu['checkboxes']['remove_selected'])
        self.active_flags = copy.deepcopy(self.flag_menu['checkboxes']['remove_selected'])

        # --------------------------------------------------------------------#
        # determine what data (if any) needs to be read

        # set variables that inform what data needs to be read (set all initially as False)
        read_all = False
        read_left = False
        read_right = False
        cut_left = False
        cut_right = False

        # determine if any of the key variables have changed
        # (network, resolution, species, qa, flags)
        # if any have changed, observations and any selected experiments have to be re-read entirely
        if (self.active_network != self.previous_active_network) or (
                self.active_resolution != self.previous_active_resolution) or (
                self.active_species != self.previous_active_species) or (
                np.array_equal(self.active_qa, self.previous_active_qa) == False) or (
                np.array_equal(self.active_flags, self.previous_active_flags) == False):
            read_all = True
        # key variables have not changed, has start/end date?
        else:
            # determine if start date/end date have changed
            if (self.active_start_date != self.previous_active_start_date) or (
                    self.active_end_date != self.previous_active_end_date):
                # if date range has changed then determine type of overlap with previous date range
                # no overlap (i.e. start date >= than previous end date, or end date <= than previous start date)?
                if (self.active_start_date >= self.previous_active_end_date) or (
                        self.active_end_date <= self.previous_active_start_date):
                    read_all = True
                # data range fully inside previous data range (i.e. start date later and end date earlier)?
                elif (self.active_start_date > self.previous_active_start_date) & (
                        self.active_end_date < self.previous_active_end_date):
                    cut_left = True
                    cut_right = True
                # need to read data on left edge and right edge of previous date range
                # (i.e. start date earlier and end date later)?
                elif (self.active_start_date < self.previous_active_start_date) & (
                        self.active_end_date > self.previous_active_end_date):
                    read_left = True
                    read_right = True
                # need to read data on left edge and cut on right edge of previous date range
                # (i.e. start date earlier and end date earlier)?
                elif (self.active_start_date < self.previous_active_start_date) & (
                        self.active_end_date < self.previous_active_end_date):
                    read_left = True
                    cut_right = True
                # need to cut data on left edge and read data on right edge of previous date range
                # (i.e. start date later and end date later)?
                elif (self.active_start_date > self.previous_active_start_date) & (
                        self.active_end_date > self.previous_active_end_date):
                    cut_left = True
                    read_right = True
                # need to read data on left edge of previous date range (i.e. start date earlier)?
                elif (self.active_start_date < self.previous_active_start_date):
                    read_left = True
                # need to read data on right edge of previous date range (i.e. end date later)?
                elif (self.active_end_date > self.previous_active_end_date):
                    read_right = True
                # need to cut data on left edge of previous date range (i.e. start date later)?
                elif (self.active_start_date > self.previous_active_start_date):
                    cut_left = True
                # need to cut data on right edge of previous date range (i.e. end date earlier)?
                elif (self.active_end_date < self.previous_active_end_date):
                    cut_right = True

        # determine if any of the active experiments have changed
        # remove experiments that are no longer selected from data_in_memory dictionary
        experiments_to_remove = [experiment for experiment in self.previous_active_experiment_grids if
                                 experiment not in self.active_experiment_grids]
        for experiment in experiments_to_remove:
            del self.data_in_memory[experiment]
        # any new experiments will need completely re-reading
        experiments_to_read = [experiment for experiment in self.active_experiment_grids if
                               experiment not in self.previous_active_experiment_grids]

        # has date range changed?
        if read_all or read_left or read_right or cut_left or cut_right :

            # set new active time array/unique station references/longitudes/latitudes
            # adjust data arrays to account for potential changing number of stations
            self.read_setup()

            # need to re-read all observations/experiments?
            if read_all:
                # reset data in memory dictionary
                self.data_in_memory = {}
                self.plotting_params = {}
                if not self.reading_nonghost:
                    self.metadata_inds_to_fill = np.arange(len(self.relevant_yearmonths))
                # read observations
                self.read_data('observations', self.active_start_date, self.active_end_date)
                # read selected experiments (iterate through)
                for data_label in self.active_experiment_grids:
                    self.read_data(data_label, self.active_start_date, self.active_end_date)
                    # if experiment in experiments_to_read list, remove it (as no longer need to read it)
                    if data_label in experiments_to_read:
                        experiments_to_read.remove(data_label)
            else:
                # if station references array has changed then as cutting/appending to
                # existing data need to rearrange existing data arrays accordingly
                if not np.array_equal(self.previous_station_references, self.station_references):
                    # get indices of stations in previous station references array in current station references array
                    old_station_inds = np.where(np.in1d(self.previous_station_references, self.station_references))[0]
                    # get indices of stations in current station references array
                    # that were in previous station references array
                    new_station_inds = np.where(np.in1d(self.station_references, self.previous_station_references))[0]

                    new_metadata_array = np.full((len(self.station_references), len(self.previous_relevant_yearmonths)),
                                                 np.NaN, dtype=self.metadata_dtype)
                    new_metadata_array[new_station_inds, :] = self.metadata_in_memory[old_station_inds, :]
                    self.metadata_in_memory = new_metadata_array

                    # iterate through all keys in data in memory dictionary
                    for data_label in list(self.data_in_memory.keys()):
                        # create new data array in shape of current station references array
                        if data_label == 'observations':
                            new_data_array = np.full((len(self.station_references), len(self.previous_time_array)),
                                                     np.NaN, dtype=self.data_dtype)
                        else:
                            new_data_array = np.full((len(self.station_references), len(self.previous_time_array)),
                                                     np.NaN, dtype=self.data_dtype[:1])
                        # put the old data into new array in the correct positions
                        new_data_array[new_station_inds, :] = self.data_in_memory[data_label][old_station_inds, :]
                        # overwrite data array with reshaped version
                        self.data_in_memory[data_label] = new_data_array

            # need to cut edges?
            if cut_left or cut_right:

                # set default edge limits as current edges
                data_left_edge_ind = 0
                data_right_edge_ind = len(self.previous_time_array)

                metadata_left_edge_ind = 0
                metadata_right_edge_ind = len(self.previous_relevant_yearmonths)

                # need to cut on left data edge?
                if cut_left:
                    data_left_edge_ind = np.where(self.previous_time_array == self.time_array[0])[0][0]
                    str_first_relevant_yearmonth = str(self.relevant_yearmonths[0])
                    str_previous_first_relevant_yearmonth = str(self.previous_relevant_yearmonths[0])
                    monthly_relative_delta = relativedelta(
                        datetime.datetime(int(str_first_relevant_yearmonth[:4]), int(str_first_relevant_yearmonth[4:6]),
                                          1, 0, 0), datetime.datetime(int(str_previous_first_relevant_yearmonth[:4]),
                                                                      int(str_previous_first_relevant_yearmonth[4:6]), 1,
                                                                      0, 0))
                    metadata_left_edge_ind = (monthly_relative_delta.years * 12) + monthly_relative_delta.months

                # need to cut on right data edge?
                if cut_right:
                    data_right_edge_ind = np.where(self.previous_time_array == self.time_array[-1])[0][0] + 1
                    str_last_relevant_yearmonth = str(self.relevant_yearmonths[-1])
                    str_previous_last_relevant_yearmonth = str(self.previous_relevant_yearmonths[-1])
                    monthly_relative_delta = relativedelta(
                        datetime.datetime(int(str_previous_last_relevant_yearmonth[:4]), int(str_previous_last_relevant_yearmonth[4:6]),
                                          1, 0, 0), datetime.datetime(int(str_last_relevant_yearmonth[:4]),
                                                                      int(str_last_relevant_yearmonth[4:6]),
                                                                      1, 0, 0))
                    metadata_right_edge_ind = metadata_right_edge_ind - ((monthly_relative_delta.years * 12) + monthly_relative_delta.months)

                # do metadata array cut
                if metadata_left_edge_ind == metadata_right_edge_ind:
                    self.metadata_in_memory = self.metadata_in_memory[:, [metadata_left_edge_ind]]
                else:
                    self.metadata_in_memory = self.metadata_in_memory[:, metadata_left_edge_ind:metadata_right_edge_ind]

                # iterate through all keys in data in memory dictionary and
                # cut edges of the associated arrays appropriately
                for data_label in list(self.data_in_memory.keys()):
                    self.data_in_memory[data_label] = self.data_in_memory[data_label][:,
                                                      data_left_edge_ind:data_right_edge_ind]

            # need to read on left edge?
            if read_left:
                # get n number of new elements on left edge
                n_new_left_data_inds = np.where(self.time_array == self.previous_time_array[0])[0][0]

                # get list of yearmonths to read
                yearmonths_to_read = get_yearmonths_to_read(self.relevant_yearmonths, self.active_start_date,
                                                            self.previous_active_start_date)
                # check which yearmonths_to_read in previous matrix
                yearmonths_in_old_matrix = np.isin(yearmonths_to_read,self.previous_relevant_yearmonths)
                # get yearmonths not currently accounted for in matrix
                new_yearmonths = yearmonths_to_read[~yearmonths_in_old_matrix]

                self.metadata_inds_to_fill = np.arange(0, len(yearmonths_to_read))
                self.metadata_in_memory = np.concatenate((np.full(
                    (len(self.station_references), len(new_yearmonths)), np.NaN, dtype=self.metadata_dtype),
                                                          self.metadata_in_memory), axis=1)

                # iterate through all keys in data in memory dictionary and
                # insert read data on left edge of the associated arrays
                for data_label in list(self.data_in_memory.keys()):
                    # add space on left edge to insert new read data
                    if data_label == 'observations':
                        self.data_in_memory[data_label] = np.concatenate((np.full(
                            (len(self.station_references), n_new_left_data_inds), np.NaN, dtype=self.data_dtype),
                                                                          self.data_in_memory[data_label]), axis=1)
                    else:
                        self.data_in_memory[data_label] = np.concatenate((np.full(
                            (len(self.station_references), n_new_left_data_inds), np.NaN, dtype=self.data_dtype[:1]),
                                                                          self.data_in_memory[data_label]), axis=1)
                    self.read_data(data_label, self.active_start_date, self.previous_active_start_date)

            # need to read on right edge?
            if read_right:
                # get n number of new elements on right edge
                n_new_right_data_inds = (len(self.time_array) - 1) - \
                                        np.where(self.time_array == self.previous_time_array[-1])[0][0]

                # get list of yearmonths to read
                yearmonths_to_read = get_yearmonths_to_read(self.relevant_yearmonths, self.previous_active_end_date,
                                                            self.active_end_date)
                # check which yearmonths_to_read in previous matrix
                yearmonths_in_old_matrix = np.isin(yearmonths_to_read,self.previous_relevant_yearmonths)
                # get yearmonths not currently accounted for in matrix
                new_yearmonths = yearmonths_to_read[~yearmonths_in_old_matrix]

                self.metadata_inds_to_fill = np.arange(-len(yearmonths_to_read), 0)
                self.metadata_in_memory = np.concatenate((self.metadata_in_memory, np.full(
                    (len(self.station_references), len(new_yearmonths)), np.NaN, dtype=self.metadata_dtype)), axis=1)

                # iterate through all keys in data in memory dictionary and
                # insert read data on right edge of the associated arrays
                for data_label in list(self.data_in_memory.keys()):
                    if data_label == 'observations':
                        self.data_in_memory[data_label] = np.concatenate((self.data_in_memory[data_label], np.full(
                            (len(self.station_references), n_new_right_data_inds), np.NaN, dtype=self.data_dtype)), axis=1)
                    else:
                        self.data_in_memory[data_label] = np.concatenate((self.data_in_memory[data_label], np.full(
                            (len(self.station_references), n_new_right_data_inds), np.NaN, dtype=self.data_dtype[:1])),
                                                                         axis=1)
                    self.read_data(data_label, self.previous_active_end_date, self.active_end_date)

            # update menu object fields
            self.update_metadata_fields()
            self.update_representativity_fields()
            self.update_period_fields()

        # if have new experiments to read, then read them now
        if len(experiments_to_read) > 0:
            for data_label in experiments_to_read:
                self.read_data(data_label, self.active_start_date, self.active_end_date)

                # --------------------------------------------------------------------#
        # if species has changed, update default species specific lower/upper limits
        if (self.active_species != self.previous_active_species):
            # update default lower/upper species specific limits and filter data outside limits
            # species_lower_limit = np.float32(self.parameter_dictionary[self.active_species]['extreme_lower_limit'])
            # species_upper_limit = np.float32(self.parameter_dictionary[self.active_species]['extreme_upper_limit'])
            species_lower_limit, species_upper_limit = self.which_bounds()
            # set default limits
            self.le_minimum_value.setText(str(species_lower_limit))
            self.le_maximum_value.setText(str(species_upper_limit))

        # --------------------------------------------------------------------#
        # update dictionary of plotting parameters (colour and zorder etc.) for each data array
        self.update_plotting_parameters()

        # --------------------------------------------------------------------#
        # run function to filter data outside lower/upper limits, not using desired
        # measurement methods, and < desired minimum data availability
        self.mpl_canvas.handle_data_filter_update()

        # --------------------------------------------------------------------#
        # update map z combobox fields based on data in memory

        # generate lists of basic and basis+bias statistics for using in the z statistic combobox
        self.basic_z_stats = np.array(list(
            OrderedDict(sorted(self.basic_stats_dict.items(), key=lambda x: x[1]['order'])).keys()))

        self.basic_and_bias_z_stats = np.append(self.basic_z_stats, list(
            OrderedDict(sorted(self.expbias_dict.items(), key=lambda x: x[1]['order'])).keys()))

        # generate list of sorted z1/z2 data arrays names in memory, putting observations
        # before experiments, and empty string item as first element in z2 array list
        # (for changing from 'difference' statistics to 'absolute')
        if len(list(self.data_in_memory.keys())) == 1:
            self.z1_arrays = np.array(['observations'])
        else:
            data_array_labels = np.array(list(self.data_in_memory.keys()))
            self.z1_arrays = np.append(['observations'],
                                       np.delete(data_array_labels, np.where(data_array_labels == 'observations')))
        self.z2_arrays = np.append([''], self.z1_arrays)

        # initialise map z statistic comboboxes
        self.mpl_canvas.handle_map_z_statistic_update()

        # --------------------------------------------------------------------#
        # update experiment bias combobox fields based on data in memory

        # if have no experiment data, all fields are empty
        if len(list(self.data_in_memory.keys())) == 1:
            self.experiment_bias_types = np.array([])
        # else, generate combobox lists
        else:
            # set all experiment bias types
            self.experiment_bias_types = np.array(['Aggregated'])

            # initialise experiment bias comboboxes
            self.mpl_canvas.handle_experiment_bias_update()

        # --------------------------------------------------------------------#
        # reset station select checkboxes to be unchecked
        self.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
        self.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
        # --------------------------------------------------------------------#

        # unset variable to allow updating of MPL canvas
        self.block_MPL_canvas_updates = False

        # update MPL canvas
        self.mpl_canvas.update_MPL_canvas()

        # Restore mouse cursor to normal
        QtWidgets.QApplication.restoreOverrideCursor()

    def reset_options(self):
        """Resets all metadata fields to initial values"""

        if self.block_MPL_canvas_updates:
            return

        #set mouse cursor to hourglass
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        # set rep fields to empty lists and initialize again
        self.representativity_menu['rangeboxes']['labels'] = []
        self.representativity_menu['rangeboxes']['current_lower'] = []
        self.update_representativity_fields()

        # set period fields to empty and initiliaze them
        self.period_menu['checkboxes']['keep_selected'] = []
        self.period_menu['checkboxes']['remove_selected'] = []
        self.update_period_fields()

        # reset metadata
        for metadata_type_ii, metadata_type in enumerate(self.metadata_menu['navigation_buttons']['labels']):
            for label in self.metadata_menu[metadata_type]['navigation_buttons']['labels']:
                self.metadata_menu[metadata_type][label]['checkboxes']['labels'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['keep_selected'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['remove_selected'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['keep_default'] = []
                self.metadata_menu[metadata_type][label]['checkboxes']['remove_default'] = []
            self.metadata_menu[metadata_type]['rangeboxes']['current_lower'] = ['nan'] * len(
                self.metadata_menu[metadata_type]['rangeboxes']['labels'])
            self.metadata_menu[metadata_type]['rangeboxes']['current_upper'] = ['nan'] * len(
                self.metadata_menu[metadata_type]['rangeboxes']['labels'])
            self.metadata_menu[metadata_type]['rangeboxes']['lower_default'] = ['nan'] * len(
                self.metadata_menu[metadata_type]['rangeboxes']['labels'])
            self.metadata_menu[metadata_type]['rangeboxes']['upper_default'] = ['nan'] * len(
                self.metadata_menu[metadata_type]['rangeboxes']['labels'])
        self.update_metadata_fields()

        # reset bounds
        species_lower_limit = np.float32(self.parameter_dictionary[self.active_species]['extreme_lower_limit'])
        species_upper_limit = np.float32(self.parameter_dictionary[self.active_species]['extreme_upper_limit'])
        # set default limits
        self.le_minimum_value.setText(str(species_lower_limit))
        self.le_maximum_value.setText(str(species_upper_limit))

        # unfilter data
        self.mpl_canvas.handle_data_filter_update()

        # Restore mouse cursor to normal
        QtWidgets.QApplication.restoreOverrideCursor()

    def read_setup(self):
        """Function that setups key variables for new read of
        observational/experiment data a time array and arrays
        of unique station references/longitudes/latitudes are created.
        """

        # force garbage collection (to avoid memory issues)
        gc.collect()

        # set current time array, as previous time array
        self.previous_time_array = self.time_array
        # set current station references, as previous station references
        self.previous_station_references = self.station_references
        # set current relevant yearmonths, as previous relevant yearmonths
        self.previous_relevant_yearmonths = self.relevant_yearmonths

        # get N time chunks between desired start date and end date to set time array
        if (self.active_resolution == 'hourly') or (self.active_resolution == 'hourly_instantaneous'):
            self.active_frequency_code = 'H'
        elif (self.active_resolution == '3hourly') or (self.active_resolution == '3hourly_instantaneous'):
            self.active_frequency_code = '3H'
        elif (self.active_resolution == '6hourly') or (self.active_resolution == '6hourly_instantaneous'):
            self.active_frequency_code = '6H'
        elif self.active_resolution == 'daily':
            self.active_frequency_code = 'D'
        elif self.active_resolution == 'monthly':
            self.active_frequency_code = 'MS'
        str_active_start_date = str(self.active_start_date)
        str_active_end_date = str(self.active_end_date)
        self.time_array = pd.date_range(start=datetime.datetime(int(str_active_start_date[:4]),
                                                                int(str_active_start_date[4:6]),
                                                                int(str_active_start_date[6:8])),
                                        end=datetime.datetime(int(str_active_end_date[:4]),
                                                              int(str_active_end_date[4:6]),
                                                              int(str_active_end_date[6:8])),
                                        freq=self.active_frequency_code)[:-1]

        if not self.reading_nonghost:
            # get all relevant observational files
            file_root = '%s/%s/%s/%s/%s/%s_' % (self.obs_root,
                                                self.active_network, self.ghost_version,
                                                self.active_resolution, self.active_species, self.active_species)
        else:
            # get files from nonghost path
            file_root = '%s/%s/%s/%s/%s/%s_' % (self.nonghost_root, self.active_network[1:].lower(),
                                                self.selected_matrix, self.active_resolution, self.active_species,
                                                self.active_species)

        self.relevant_yearmonths = np.sort([yyyymm for yyyymm in
                                           self.available_observation_data[self.active_network][self.active_resolution][
                                               self.active_matrix][self.active_species]])
        relevant_files = sorted([file_root+str(yyyymm)[:6]+'.nc' for yyyymm in self.relevant_yearmonths])
        self.N_inds_per_month = np.array([np.count_nonzero(np.all(
            [self.time_array >= datetime.datetime.strptime(str(start_yyyymm), '%Y%m%d'),
             self.time_array < datetime.datetime.strptime(str(self.relevant_yearmonths[month_ii + 1]), '%Y%m%d')],
            axis=0)) if month_ii != (len(self.relevant_yearmonths) - 1) else np.count_nonzero(
            self.time_array >= datetime.datetime.strptime(str(start_yyyymm), '%Y%m%d')) for month_ii, start_yyyymm in
                                          enumerate(self.relevant_yearmonths)])

        self.station_references = []
        self.station_longitudes = []
        self.station_latitudes = []
        if not self.reading_nonghost:
            for relevant_file in relevant_files:
                ncdf_root = Dataset(relevant_file)
                self.station_references = np.append(self.station_references, ncdf_root['station_reference'][:])
                self.station_longitudes = np.append(self.station_longitudes, ncdf_root['longitude'][:])
                self.station_latitudes = np.append(self.station_latitudes, ncdf_root['latitude'][:])
                ncdf_root.close()
            self.station_references, station_unique_indices = np.unique(self.station_references, return_index=True)
            self.station_longitudes = self.station_longitudes[station_unique_indices]
            self.station_latitudes = self.station_latitudes[station_unique_indices]
        else:
            # first, try to take the data files and handle in case of daily files
            if os.path.exists(relevant_files[0]):
                ncdf_root = Dataset(relevant_files[0])
            else:
                relevant_files = sorted([file_root + str(yyyymm)[:8] + '.nc' for yyyymm in self.relevant_yearmonths])
                ncdf_root = Dataset(relevant_files[0])
            self.station_references = np.array(
                [st_name.tostring().decode('ascii').replace('\x00', '') for st_name in ncdf_root['station_name'][:]],
                dtype=np.str)
            # get staion refs
            if "latitude" in ncdf_root.variables:
                self.station_longitudes = np.append(self.station_longitudes, ncdf_root['longitude'][:])
                self.station_latitudes = np.append(self.station_latitudes, ncdf_root['latitude'][:])
            else:
                self.station_longitudes = np.append(self.station_longitudes, ncdf_root['lon'][:])
                self.station_latitudes = np.append(self.station_latitudes, ncdf_root['lat'][:])
            ncdf_root.close()

        # update measurement units for species (take standard units from parameter dictionary)
        self.measurement_units = self.parameter_dictionary[self.active_species]['standard_units']

        # set data variables to read (dependent on active data resolution)
        if (self.active_resolution == 'hourly') or (self.active_resolution == 'hourly_instantaneous'):
            self.data_vars_to_read = [self.active_species, 'hourly_native_representativity_percent',
                                      'daily_native_representativity_percent',
                                      'monthly_native_representativity_percent',
                                      'annual_native_representativity_percent', 'hourly_native_max_gap_percent',
                                      'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                      'annual_native_max_gap_percent', 'day_night_code', 'weekday_weekend_code',
                                      'season_code', 'time']
        elif (self.active_resolution == 'daily') or (self.active_resolution == '3hourly') or \
                (self.active_resolution == '6hourly') or (self.active_resolution == '3hourly_instantaneous') or \
                (self.active_resolution == '6hourly_instantaneous'):
            self.data_vars_to_read = [self.active_species, 'daily_native_representativity_percent',
                                      'monthly_native_representativity_percent',
                                      'annual_native_representativity_percent',
                                      'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                      'annual_native_max_gap_percent', 'weekday_weekend_code', 'season_code', 'time']
        elif self.active_resolution == 'monthly':
            self.data_vars_to_read = [self.active_species, 'monthly_native_representativity_percent',
                                      'annual_native_representativity_percent', 'monthly_native_max_gap_percent',
                                      'annual_native_max_gap_percent', 'season_code', 'time']

        # set data dtype
        self.data_dtype = [(key, np.float32) for key in self.data_vars_to_read]

    def get_yearmonths_to_read(self, yearmonths, start_date_to_read, end_date_to_read):
        """Function that returns the yearmonths of the files needed to be read.
           This is done by limiting a list of relevant yearmonths by a start/end date
        """

        first_valid_file_ind = bisect.bisect_right(yearmonths, int(start_date_to_read))
        if first_valid_file_ind != 0:
            first_valid_file_ind = first_valid_file_ind - 1
        last_valid_file_ind = bisect.bisect_left(yearmonths, int(end_date_to_read))
        if first_valid_file_ind == last_valid_file_ind:
            return [yearmonths[first_valid_file_ind]]
        else:
            return yearmonths[first_valid_file_ind:last_valid_file_ind]

    def read_data(self, data_label, start_date_to_read, end_date_to_read):
        """Function that handles reading of observational/experiment data"""

        # force garbage collection (to avoid memory issues)
        gc.collect()

        # set process_type variable to self (for access by parallel read function)
        # also get relevant file start dates
        if data_label == 'observations':
            self.process_type = 'observations'
            if not self.reading_nonghost:
                file_root = '%s/%s/%s/%s/%s/%s_' % (self.obs_root,
                                                    self.active_network, self.ghost_version,
                                                    self.active_resolution, self.active_species, self.active_species)
                relevant_file_start_dates = \
                    sorted(self.available_observation_data[self.active_network]
                           [self.active_resolution][self.active_matrix][self.active_species])
            else:
                # get files from nonghost path
                file_root = '%s/%s/%s/%s/%s/%s_' % (self.nonghost_root, self.active_network[1:].lower(),
                                                    self.selected_matrix, self.active_resolution, self.active_species,
                                                    self.active_species)
                relevant_file_start_dates = \
                    sorted(self.available_observation_data[self.active_network]
                           [self.active_resolution][self.active_matrix][self.active_species])

        else:
            self.process_type = 'experiment'
            experiment_grid_split = data_label.split('-')
            active_experiment = experiment_grid_split[0]
            active_grid = experiment_grid_split[1]
            active_member = experiment_grid_split[2]
            file_root = \
                '%s/%s/%s/%s/%s/%s/%s/%s/%s_' % (self.exp_root, self.ghost_version, active_experiment,
                                              active_grid, active_member, self.active_resolution,
                                              self.active_species, self.active_network, self.active_species)
            relevant_file_start_dates = sorted(self.available_experiment_data[data_label])

        # get data files in required date range to read, taking care not to re-read what has already been read
        yearmonths_to_read = get_yearmonths_to_read(relevant_file_start_dates, start_date_to_read, end_date_to_read)
        relevant_files = [file_root+str(yyyymm)[:6]+'.nc' for yyyymm in yearmonths_to_read]

        if not os.path.exists(relevant_files[0]):
            relevant_files = sorted([file_root + str(yyyymm)[:8] + '.nc' for yyyymm in self.relevant_yearmonths])

        # check if data label in data in memory dictionary
        if data_label not in list(self.data_in_memory.keys()):
            # if not create empty array (filled with NaNs) to store species data and place it in the dictionary

            if self.process_type == 'observations':
                self.plotting_params['observations'] = {}
                if not self.reading_nonghost:
                    self.data_in_memory[data_label] = np.full((len(self.station_references), len(self.time_array)),
                                                          np.NaN, dtype=self.data_dtype)
                else:
                    self.data_in_memory[data_label] = np.full((len(self.station_references), len(self.time_array)),
                                                          np.NaN, dtype=self.data_dtype[:1])
                self.metadata_in_memory = np.full((len(self.station_references), len(self.relevant_yearmonths)),
                                                  np.NaN, dtype=self.metadata_dtype)
                if self.reading_nonghost:
                    tmp_ncdf = Dataset(relevant_files[0])
                    # create separate structure of nonghost metadata
                    nonghost_mdata_dtype = [('station_name', np.object), ('latitude', np.float),
                                            ('longitude', np.float), ('altitude', np.float)]
                    if "station_code" in tmp_ncdf.variables:
                        nonghost_mdata_dtype.append(('station_reference', np.object))
                    self.nonghost_metadata = np.full((len(self.station_references)),
                                                     np.NaN, dtype=nonghost_mdata_dtype)

            # if process_type is experiment, get experiment specific grid edges from
            # first relevant file, and save to data in memory dictionary
            if self.process_type == 'experiment':
                self.data_in_memory[data_label] = np.full((len(self.station_references), len(self.time_array)),
                                                          np.NaN, dtype=self.data_dtype[:1])
                self.plotting_params[data_label] = {}
                exp_nc_root = Dataset(relevant_files[0])
                self.plotting_params[data_label]['grid_edge_longitude'] = exp_nc_root['grid_edge_longitude'][:]
                self.plotting_params[data_label]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
                exp_nc_root.close()

        # iterate and read species data in all relevant netCDF files (either in serial/parallel)

        # read serially
        if self.read_type == 'serial':

            # iterate through relevant netCDF files
            for relevant_file in relevant_files:
                # create argument tuple of function
                tuple_arguments = relevant_file, self.time_array, self.station_references, \
                                  self.active_species, self.process_type,\
                                  self.active_qa, self.active_flags, \
                                  self.data_dtype, self.data_vars_to_read, \
                                  self.metadata_dtype, self.metadata_vars_to_read
                # read file
                file_data, time_indices, full_array_station_indices = read_netcdf_data(tuple_arguments)
                # place read data into big array as appropriate
                self.data_in_memory[data_label]['data'][full_array_station_indices[np.newaxis, :],
                                                        time_indices[:, np.newaxis]] = file_data

        # read in parallel
        elif self.read_type == 'parallel':

            # setup pool of N workers on N CPUs
            pool = multiprocessing.Pool(self.n_cpus)

            # read netCDF files in parallel
            if not self.reading_nonghost:
                tuple_arguments = [(file_name, self.time_array, self.station_references, self.active_species,
                                    self.process_type, self.active_qa, self.active_flags,
                                    self.data_dtype, self.data_vars_to_read,
                                    self.metadata_dtype, self.metadata_vars_to_read) for
                                   file_name in relevant_files]
                all_file_data = pool.map(read_netcdf_data, tuple_arguments)
            else:
                tuple_arguments = [
                    (file_name, self.time_array, self.station_references, self.active_species, self.process_type) for
                    file_name in relevant_files]
                all_file_data = pool.map(read_netcdf_nonghost, tuple_arguments)

            # will not submit more files to pool, so close access to it
            pool.close()
            # wait for worker processes to terminate before continuing
            pool.join()

            # iterate through read file data and place data into data array as appropriate
            for file_data_ii, file_data in enumerate(all_file_data):
                try:
                    # some file_data might be none, in case the file did not exist
                    self.data_in_memory[data_label][file_data[2][:, np.newaxis], file_data[1][np.newaxis, :]] = \
                        file_data[0]
                except Exception as e:
                    continue
                if self.process_type == 'observations':
                    if not self.reading_nonghost:
                        self.metadata_in_memory[file_data[2][:, np.newaxis],
                                                self.metadata_inds_to_fill[file_data_ii]] = file_data[3]
                    else:
                        self.nonghost_metadata[file_data[2][:, np.newaxis]] = file_data[3]

    def update_metadata_fields(self):

        """update the metadata menu object with metadata associated with newly read data
           for non-numeric metadata gets all the unique fields per metadata variable,
           and sets the available fields as such, and for numeric gets the minimum and maximum
           boundaries of each metadata variable. 

           If previously metadata settings for a field deviate from the default, then if the same field still 
           exists then the settings (i.e. bounds or checkbox selection) are copied across, rather than setting to the default.  
        """

        # iterate through metadata variables
        for meta_var in self.metadata_vars_to_read:

            meta_var_field = self.metadata_in_memory[meta_var]

            # get metadata variable type/data type
            metadata_type = self.standard_metadata[meta_var]['metadata_type']
            metadata_data_type = self.standard_metadata[meta_var]['data_type']

            # remove NaNs from field
            meta_var_field_nan_removed = meta_var_field[~pd.isnull(meta_var_field)]

            # update pop-up metadata menu object with read metadata values
            # for non-numeric metadata gets all the unique fields per metadata variable
            # and sets the available fields as such
            if metadata_data_type == np.object:
                #get previous fields
                previous_fields = copy.deepcopy(self.metadata_menu[metadata_type][meta_var]['checkboxes']['labels'])
                #update new labels
                self.metadata_menu[metadata_type][meta_var]['checkboxes']['labels'] = np.unique(
                    meta_var_field_nan_removed)
                #if field previously existed, then copy across checkbox settings for field
                #else set initial checkboxes to be all blank
                previous_keep = copy.deepcopy(self.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected'])
                previous_remove = copy.deepcopy(self.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_selected'])
                self.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected'] = []
                self.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_selected'] = []
                for field in self.metadata_menu[metadata_type][meta_var]['checkboxes']['labels']:
                    if field in previous_fields:
                        if field in previous_keep:
                            self.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected'].append(field)
                        if field in previous_remove:
                            self.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_selected'].append(field)
                #set defaults to be empty
                self.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_default'] = []
                self.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_default'] = []
            # for numeric fields get the minimum and maximum boundaries of each metadata variable
            # if previous set values vary from min/max boundaries, copy across the values
            # set as min/max as nan if have no numeric metadata for variable
            else:
                meta_var_index = self.metadata_menu[metadata_type]['rangeboxes']['labels'].index(meta_var)
                #have some numeric values for metadata variable?
                if len(meta_var_field_nan_removed) > 0:
                    min_val = str(np.min(meta_var_field_nan_removed))
                    max_val = str(np.max(meta_var_field_nan_removed))
                    #get previous lower/upper extents and defaults
                    previous_lower_default = copy.deepcopy(self.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index])
                    previous_upper_default = copy.deepcopy(self.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index])
                    previous_lower = copy.deepcopy(self.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index])
                    previous_upper = copy.deepcopy(self.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index])
                    #if previous lower > previous default lower bound then copy across (and also not 'nan') 
                    #initially set as min extent
                    self.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = min_val
                    if (previous_lower != 'nan') & (previous_lower_default != 'nan'):   
                        if previous_lower > previous_lower_default:
                            self.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = copy.deepcopy(previous_lower)
                    #if previous upper < previous default upper bound then copy across (and also not 'nan')
                    #initially set as max extent
                    self.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = max_val
                    if (previous_upper != 'nan') & (previous_upper_default != 'nan'):
                        if previous_upper < previous_upper_default:
                            self.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = copy.deepcopy(previous_upper)
                    #set defaults to min/max extents
                    self.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index] = min_val
                    self.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index] = max_val
                #do not have some numeric values for metadata variable so set as 'nan'
                else:
                    self.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = 'nan'
                    self.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index] = 'nan'
                    self.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = 'nan'
                    self.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index] = 'nan'

    def update_representativity_fields(self):

        '''update the data representativity menu -> 1D list of rangebox values
           dependent on the temporal resolution, some fields will appear or not
        '''

        # get previous set labels
        previous_labels = copy.deepcopy(self.representativity_menu['rangeboxes']['labels'])

        # get previously set rangebox values
        previous_lower = copy.deepcopy(self.representativity_menu['rangeboxes']['current_lower'])

        # hourly temporal resolution?
        if (self.active_resolution == 'hourly') or (self.active_resolution == 'hourly_instantaneous'):
            self.representativity_menu['rangeboxes']['labels'] = ['hourly_native_representativity_percent',
                                                                  'hourly_native_max_gap_percent',
                                                                  'daily_native_representativity_percent',
                                                                  'daily_representativity_percent',
                                                                  'daily_native_max_gap_percent',
                                                                  'daily_max_gap_percent',
                                                                  'monthly_native_representativity_percent',
                                                                  'monthly_representativity_percent',
                                                                  'monthly_native_max_gap_percent',
                                                                  'monthly_max_gap_percent',
                                                                  'all_representativity_percent', 'all_max_gap_percent']
        # daily temporal resolution?
        elif (self.active_resolution == 'daily') or (self.active_resolution == '3hourly') or \
                (self.active_resolution == '6hourly') or (self.active_resolution == '3hourly_instantaneous') or \
                (self.active_resolution == '6hourly_instantaneous'):
            self.representativity_menu['rangeboxes']['labels'] = ['daily_native_representativity_percent',
                                                                  'daily_native_max_gap_percent',
                                                                  'monthly_native_representativity_percent',
                                                                  'monthly_representativity_percent',
                                                                  'monthly_native_max_gap_percent',
                                                                  'monthly_max_gap_percent',
                                                                  'all_representativity_percent', 'all_max_gap_percent']
        # monthly temporal resolution?
        elif self.active_resolution == 'monthly':
            self.representativity_menu['rangeboxes']['labels'] = ['monthly_native_representativity_percent',
                                                                  'monthly_native_max_gap_percent',
                                                                  'all_representativity_percent', 'all_max_gap_percent']

        # initialise rangebox values --> for data representativity fields
        # the default is 0%, for max gap fields % the default is 100%
        self.representativity_menu['rangeboxes']['current_lower'] = []
        for label_ii, label in enumerate(self.representativity_menu['rangeboxes']['labels']):
            if 'max_gap' in label:
                self.representativity_menu['rangeboxes']['current_lower'].append('100')
            else:
                self.representativity_menu['rangeboxes']['current_lower'].append('0')

            # label previously existed?
            if label in previous_labels:
                self.representativity_menu['rangeboxes']['current_lower'][label_ii] = previous_lower[
                    previous_labels.index(label)]

    def update_period_fields(self):

        '''update the data period menu -> list of checkboxes
           dependent on the temporal resolution, some fields will appear or not
        '''

        # hourly temporal resolution?
        if 'hourly' in self.active_resolution:
            self.period_menu['checkboxes']['labels'] = ['Daytime', 'Nighttime', 'Weekday', 'Weekend', 'Spring',
                                                        'Summer', 'Autumn', 'Winter']
        # daily temporal resolution?
        elif self.active_resolution == 'daily':
            self.period_menu['checkboxes']['labels'] = ['Weekday', 'Weekend', 'Spring', 'Summer', 'Autumn', 'Winter']
        # monthly temporal resolution?
        elif self.active_resolution == 'monthly':
            self.period_menu['checkboxes']['labels'] = ['Spring', 'Summer', 'Autumn', 'Winter']

    def update_plotting_parameters(self):
        """Function that updates plotting parameters (colour
        and zorder) for each selected data array
        """

        # assign a colour/zorder to all selected data arrays

        # define observations colour to be 'black'
        self.plotting_params['observations']['colour'] = 'black'
        # define zorder of observations to be 5
        self.plotting_params['observations']['zorder'] = 5

        # generate a list of RGB tuples for number of experiments there are
        sns.reset_orig()
        clrs = sns.color_palette('husl', n_colors=len(list(self.data_in_memory.keys()))-1)

        # iterate through sorted experiment names, assigning each experiment a new RGB colour tuple, and zorder
        experiment_ind = 1
        for experiment in sorted(list(self.data_in_memory.keys())):
            if experiment != 'observations':
                # define colour for experiment
                self.plotting_params[experiment]['colour'] = clrs[experiment_ind-1]
                # define zorder for experiment (obs zorder + experiment_ind)
                self.plotting_params[experiment]['zorder'] = \
                    self.plotting_params['observations']['zorder'] + experiment_ind
                # update count of experiments
                experiment_ind += 1

    def load_conf(self, section=None, fpath=None):
        """ Load existing configurations from file. """

        from .config import read_conf

        if fpath is None:
            fpath = parse_path(self.config_dir, self.config_file)

        if not os.path.isfile(fpath):
            print(("Error %s" % fpath))
            return

        opts = read_conf(section, fpath)
        if section is None:
            return opts

        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in opts.items()})

    def disable_ghost_buttons(self):
        """Disable button related only to ghost data"""
        # TODO: add all ghost related fields to a list
        # and set to False in a list-comprehension way
        # change background-color to indicate that it's nonusable
        self.bu_flags.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        self.bu_rep.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        self.bu_meta.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        self.bu_period.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        self.bu_QA.setStyleSheet("""QPushButton:disabled {background-color:#DCDCDC;}""")
        # and disable
        self.bu_QA.setEnabled(False)
        self.bu_flags.setEnabled(False)
        self.bu_rep.setEnabled(False)
        self.bu_meta.setEnabled(False)
        self.bu_period.setEnabled(False)

    def enable_ghost_buttons(self):
        """Enable button related only to ghost data"""
        self.bu_QA.setEnabled(True)
        self.bu_flags.setEnabled(True)
        self.bu_rep.setEnabled(True)
        self.bu_meta.setEnabled(True)
        self.bu_period.setEnabled(True)


# generate Providentia dashboard
def main(**kwargs):
    """Main function"""
    q_app = QtWidgets.QApplication(sys.argv)
    q_app.setStyle("Fusion")
    ProvidentiaMainWindow(**kwargs)
    sys.exit(q_app.exec_())
