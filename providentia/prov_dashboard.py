from .configuration import ProvConfiguration
from .configuration import parse_path
from .reading import read_netcdf_station
from .reading import read_netcdf_data
from .prov_canvas import MPLCanvas
from .prov_canvas import NavigationToolbar
from .prov_dashboard_aux import ComboBox
from .prov_dashboard_aux import QVLine
# from .prov_dashboard_aux import QHLine
from .prov_dashboard_aux import PopUpWindow

import bisect
import datetime
import gc
import multiprocessing
import os
import json
from collections import OrderedDict
import sys

from netCDF4 import Dataset
import numpy as np
import pandas as pd
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui
import seaborn as sns


class GenerateProvidentiaDashboard(QtWidgets.QWidget, ProvConfiguration):
    """Define class that generates Providentia dashboard"""

    def __init__(self, read_type='parallel', **kwargs):
        super(GenerateProvidentiaDashboard, self).__init__()
        ProvConfiguration.__init__(self, **kwargs)

        # put read_type into self
        self.read_type = read_type

        # store options to be restored at the end
        # self.localvars = copy.deepcopy(vars(self))
        # update from config file
        if ('config' in kwargs) and ('section' in kwargs):
            self.load_conf(kwargs['section'], kwargs['config'])
        # update from command line
        vars(self).update({(k, self.parse_parameter(val)) for k, val in kwargs.items()})
        # arguments are only local

        # create UI
        self.init_ui()

#    def __setattr__(self, key, value):
#        super(GenerateProvidentiaDashboard, self).__setattr__(key, parse_parameter(key, value))

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

        # setup configuration bar with combo boxes, input boxes and buttons
        # use a gridded layout to place objects
        config_bar = QtWidgets.QGridLayout()

        # define spacing/margin variables
        config_bar.setHorizontalSpacing(3)
        config_bar.setVerticalSpacing(1)
        config_bar.setContentsMargins(5, 0, 0, 0)
        config_bar.setAlignment(QtCore.Qt.AlignLeft)

        # define all configuration box objects (labels, comboboxes etc.)
        title_font = QtGui.QFont()
        title_font.setUnderline(True)
        self.lb_data_selection = QtWidgets.QLabel(self, text="Data Selection")
        self.lb_data_selection.setFont(title_font)
        self.lb_data_selection.setToolTip('Setup configuration of data to read into memory')
        self.bu_read = QtWidgets.QPushButton('READ', self)
        self.bu_read.setMinimumWidth(40)
        self.bu_read.setMaximumWidth(40)
        self.bu_read.setStyleSheet("color: red;")
        self.bu_read.setToolTip('Read selected configuration of data into memory')
        self.ch_colocate = QtWidgets.QCheckBox("Colocate")
        self.ch_colocate.setToolTip('')
        # self.cb_network = QtWidgets.QComboBox(self)
        self.cb_network = ComboBox(self)
        self.cb_network.setMinimumWidth(95)
        self.cb_network.setMaximumWidth(95)
        self.cb_network.setToolTip('Select providing observational data network')
        # self.cb_resolution = QtWidgets.QComboBox(self)
        self.cb_resolution = ComboBox(self)
        self.cb_resolution.setMinimumWidth(95)
        self.cb_resolution.setMaximumWidth(95)
        self.cb_resolution.setToolTip('Select temporal resolution of data')
        # self.cb_matrix = QtWidgets.QComboBox(self)
        self.cb_matrix = ComboBox(self)
        self.cb_matrix.setMinimumWidth(95)
        self.cb_matrix.setMaximumWidth(95)
        self.cb_matrix.setToolTip('Select data matrix')
        # self.cb_species = QtWidgets.QComboBox(self)
        self.cb_species = ComboBox(self)
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
        self.lb_data_filter = QtWidgets.QLabel(self, text="Data Filter")
        self.lb_data_filter.setMinimumWidth(65)
        self.lb_data_filter.setMaximumWidth(65)
        self.lb_data_filter.setToolTip('Select criteria to filter data by')
        self.lb_data_filter.setFont(title_font)
        self.bu_screen = QtWidgets.QPushButton('FILTER', self)
        self.bu_screen.setMinimumWidth(46)
        self.bu_screen.setMaximumWidth(46)
        self.bu_screen.setStyleSheet("color: blue;")
        self.bu_screen.setToolTip('')
        self.lb_data_bounds = QtWidgets.QLabel(self, text="Bounds")
        self.lb_data_bounds.setMinimumWidth(47)
        self.lb_data_bounds.setMaximumWidth(47)
        self.lb_data_bounds.setToolTip('')
        self.le_minimum_value = QtWidgets.QLineEdit(self)
        self.le_minimum_value.setMinimumWidth(60)
        self.le_minimum_value.setMaximumWidth(60)
        self.le_minimum_value.setToolTip('')
        self.le_maximum_value = QtWidgets.QLineEdit(self)
        self.le_maximum_value.setMinimumWidth(60)
        self.le_maximum_value.setMaximumWidth(60)
        self.le_maximum_value.setToolTip('')
        self.lb_minimum_data_availability = QtWidgets.QLabel(self, text="% Min.")
        self.lb_minimum_data_availability.setMinimumWidth(47)
        self.lb_minimum_data_availability.setMaximumWidth(47)
        self.lb_minimum_data_availability.setToolTip('')
        self.le_minimum_data_availability = QtWidgets.QLineEdit(self)
        self.le_minimum_data_availability.setMinimumWidth(45)
        self.le_minimum_data_availability.setMaximumWidth(45)
        self.le_minimum_data_availability.setToolTip('')
        self.bu_methods = QtWidgets.QPushButton('METHOD', self)
        self.bu_methods.setMinimumWidth(60)
        self.bu_methods.setMaximumWidth(60)
        self.bu_methods.setToolTip('')
        self.vertical_splitter_2 = QVLine()
        self.vertical_splitter_2.setMaximumWidth(20)
        self.lb_z = QtWidgets.QLabel(self, text="Map Z")
        self.lb_z.setFont(title_font)
        self.lb_z.setToolTip('')
        # self.cb_z_stat = QtWidgets.QComboBox(self)
        self.cb_z_stat = ComboBox(self)
        self.cb_z_stat.setMinimumWidth(80)
        self.cb_z_stat.setMaximumWidth(80)
        self.cb_z_stat.setToolTip('')
        # self.cb_z1 = QtWidgets.QComboBox(self)
        self.cb_z1 = ComboBox(self)
        self.cb_z1.setMinimumWidth(125)
        self.cb_z1.setMaximumWidth(125)
        self.cb_z1.setToolTip('')
        # self.cb_z2 = QtWidgets.QComboBox(self)
        self.cb_z2 = ComboBox(self)
        self.cb_z2.setMinimumWidth(125)
        self.cb_z2.setMaximumWidth(125)
        self.cb_z2.setToolTip('')
        self.vertical_splitter_3 = QVLine()
        self.vertical_splitter_3.setMaximumWidth(20)
        self.lb_experiment_bias = QtWidgets.QLabel(self, text="Exp. Bias")
        self.lb_experiment_bias.setFont(title_font)
        self.lb_experiment_bias.setToolTip('')
        # self.cb_experiment_bias_type = QtWidgets.QComboBox(self)
        self.cb_experiment_bias_type = ComboBox(self)
        self.cb_experiment_bias_type.setMinimumWidth(100)
        self.cb_experiment_bias_type.setMaximumWidth(100)
        self.cb_experiment_bias_type.setToolTip('')
        # self.cb_experiment_bias_stat = QtWidgets.QComboBox(self)
        self.cb_experiment_bias_stat = ComboBox(self)
        self.cb_experiment_bias_stat.setMinimumWidth(100)
        self.cb_experiment_bias_stat.setMaximumWidth(100)
        self.cb_experiment_bias_stat.setToolTip('')
        self.vertical_splitter_4 = QVLine()
        self.vertical_splitter_4.setMaximumWidth(20)
        self.lb_station_selection = QtWidgets.QLabel(self, text="Site Select")
        self.lb_station_selection.setFont(title_font)
        self.lb_station_selection.setToolTip('')
        self.ch_select_all = QtWidgets.QCheckBox("All")
        self.ch_select_all.setToolTip('')
        self.ch_intersect = QtWidgets.QCheckBox("Intersect")
        self.ch_intersect.setToolTip('')

        # position objects on gridded configuration bar
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

        # enable dynamic updating of configuration bar fields which filter data files
        self.cb_network.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.cb_resolution.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.cb_matrix.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.cb_species.currentTextChanged.connect(self.config_bar_params_change_handler)
        self.le_start_date.textChanged.connect(self.config_bar_params_change_handler)
        self.le_end_date.textChanged.connect(self.config_bar_params_change_handler)

        # enable pop up configuration windows
        self.bu_experiments.clicked.connect(self.handle_pop_up_experiments_window)
        self.bu_flags.clicked.connect(self.handle_pop_up_flags_window)
        self.bu_QA.clicked.connect(self.handle_pop_up_qa_window)
        self.bu_classifications.clicked.connect(self.handle_pop_up_classifications_window)
        self.bu_methods.clicked.connect(self.handle_pop_up_methods_window)

        # define all parameters
        self.parameter_dictionary = json.load(open('providentia/conf/parameter.json'))

        # define data provider flags
        self.standard_data_flag_codes = \
            json.load(open('providentia/conf/standard_data_flag_codes.json'))

        self.flag_names = np.array(sorted(self.standard_data_flag_codes,
                                          key=self.standard_data_flag_codes.get))
        self.flag_codes = np.sort(list(self.standard_data_flag_codes.values()))
        self.flag_default_codes = np.array([1, 2, 3, 10, 11, 12, 13, 14, 15, 16, 20,
                                            21, 24, 25, 26, 29, 30, 31, 32, 40, 41, 42,
                                            43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54,
                                            55, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                                            80, 81, 90, 150, 154, 155, 156, 157], dtype=np.uint8)
        self.flag_default_inds = np.array(
            [np.where(self.flag_codes == code)[0][0] for code in self.flag_default_codes], dtype=np.uint8)

        # define qa flags
        self.standard_qa_flag_codes = json.load(open('providentia/conf/standard_qa_flag_codes.json'))

        self.qa_names = np.array(sorted(
            self.standard_qa_flag_codes, key=self.standard_qa_flag_codes.get))
        self.qa_codes = np.sort(list(self.standard_qa_flag_codes.values()))
        self.qa_default_codes = np.array([0, 1, 2, 3, 4, 5, 70, 80, 81, 82, 85, 87, 90, 94, 100, 106, 107, 108],
                                         dtype=np.uint8)
        self.qa_default_inds = np.array([np.where(self.qa_codes == code)[0][0] for code in self.qa_default_codes],
                                        dtype=np.uint8)

        # define classification flags
        self.standard_classification_flag_codes = \
            json.load(open('providentia/conf/standard_classification_flag_codes.json'))

        self.classification_names = np.array(sorted(
            self.standard_classification_flag_codes, key=self.standard_classification_flag_codes.get))

        self.classification_codes = np.sort(list(self.standard_classification_flag_codes.values()))
        self.classification_default_codes_to_retain = np.array([5], dtype=np.uint8)
        self.classification_default_codes_to_remove = np.array([0, 4, 46, 49], dtype=np.uint8)
        self.classification_default_inds_to_retain = np.array([np.where(self.classification_codes == code)[0][0]
                                                               for code in self.classification_default_codes_to_retain],
                                                              dtype=np.uint8)
        self.classification_default_inds_to_remove = np.array([np.where(self.classification_codes == code)[0][0]
                                                               for code in self.classification_default_codes_to_remove],
                                                              dtype=np.uint8)

        # create dictionary to hold indices of selected values in pop-up windows
        self.selected_indices = {
            'EXPERIMENTS': [[]], 'FLAGS': [[]], 'QA': [[]], 'CLASSIFICATIONS': [[], []], 'METHODS': [[]]}

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

        # ------------------------------------------------------------------------# 
        # position config bar, navigation toolbar and MPL canvas and elements in parent layout

        # add config bar to parent frame
        parent_layout.addLayout(config_bar)

        # add MPL navigation toolbar to parent frame
        parent_layout.addWidget(self.navi_toolbar)

        # add MPL canvas of plots to parent frame
        parent_layout.addWidget(self.mpl_canvas)

        # set finalised layout
        self.setLayout(parent_layout)

        # plot whole dashboard
        self.show()

        # maximise window to fit screen
        self.showMaximized()

    def update_configuration_bar_fields(self):
        """Define function that initialises/updates configuration bar fields"""

        # set variable to block interactive handling while updating config bar parameters
        self.block_config_bar_handling_updates = True

        # set some default configuration values when initialising config bar
        if self.config_bar_initialisation is True:
            # set initially selected/active start-end date as default 201601-201701
            self.le_start_date.setText('20160101')
            self.le_end_date.setText('20170101')
            self.selected_start_date = int(self.le_start_date.text())
            self.selected_end_date = int(self.le_end_date.text())
            self.selected_start_date_firstdayofmonth = int(str(self.selected_start_date)[:6]+'01')
            self.active_start_date = int(self.le_start_date.text())
            self.active_end_date = int(self.le_end_date.text())
            self.date_range_has_changed = False

            # set initially selected minimum data availability % to 0.0
            self.le_minimum_data_availability.setText('0.0')

            # set selected/active values of other fields to be initially None
            self.selected_network = None
            self.active_network = None
            self.selected_resolution = None
            self.active_resolution = None
            self.selected_matrix = None
            self.active_matrix = None
            self.selected_species = None
            self.active_species = None

            # set selected/active values of variables associated
            # with pop up windows to be empty lists
            self.active_experiment_grids = []
            self.active_qa_inds = []
            self.active_flag_inds = []
            self.active_classifications_to_retain_inds = []
            self.active_classifications_to_remove_inds = []

            # set initial time array to be None
            self.time_array = None

            # set initial station references to be empty list
            self.station_references = []
            # set initial unique station methods to be empty list
            self.station_unique_methods = []

            # Gather all observational data
            # create nested dictionary storing all observational species
            # data by species matrix, by temporal resolution, by network,
            # associated with list of start YYYYMM yearmonths of data files

            self.all_observation_data = {}

            # set all available networks
            available_networks = ['EBAS', 'EEA_AQ_eReporting']

            # set all available temporal resolutions
            available_resolutions = ['hourly', 'daily', 'monthly']

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

            # create dictionary of observational data inside date range
            self.get_valid_obs_files_in_date_range()

        # if date range has changed then update available observational data dictionary
        if self.date_range_has_changed is True:
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
        available_networks = sorted(list(self.available_observation_data.keys()))
        self.cb_network.addItems(available_networks)
        if self.selected_network in available_networks:
            self.cb_network.setCurrentText(self.selected_network)
        else:
            self.selected_network = self.cb_network.currentText()

        # update resolution field
        available_resolutions = list(self.available_observation_data[self.cb_network.currentText()].keys())
        # manually force order of available resolutions
        resolution_order_dict = {'hourly': 1, 'daily': 2, 'monthly': 3}
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
        available_species = sorted(self.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()][self.cb_matrix.currentText()])
        self.cb_species.addItems(available_species)
        if self.selected_species in available_species:
            self.cb_species.setCurrentText(self.selected_species)
        else:
            self.selected_species = self.cb_species.currentText()

        # update available experiment data dictionary
        self.get_valid_experiment_files_in_date_range()
        # update selected indices for experiments -- keeping previously selected experiments if available
        if len(self.selected_indices['EXPERIMENTS'][0]) > 0:
            previous_selected_experiments = \
                self.previous_available_experiment_grids[self.selected_indices['EXPERIMENTS'][0]]
        else:
            previous_selected_experiments = []
        # set selected indices as previously selected indices in current available list of experiments
        selected_experiments = \
            [previous_selected_experiment for previous_selected_experiment in previous_selected_experiments
             if previous_selected_experiment in self.available_experiment_grids]

        selected_experiment_inds = \
            np.array([np.where(self.available_experiment_grids == selected_experiment)[0][0]
                      for selected_experiment in selected_experiments], dtype=np.uint8)

        self.selected_indices['EXPERIMENTS'] = [selected_experiment_inds]
        # set previous available experiments variable
        self.previous_available_experiment_grids = np.array(self.available_experiment_grids)

        # update selected indices for QA
        # if initialising config bar then check default selection
        if self.config_bar_initialisation is True:
            self.selected_indices['QA'] = [self.qa_default_inds]

        # unset variable to allow interactive handling from now
        self.block_config_bar_handling_updates = False

    def get_valid_obs_files_in_date_range(self):
        """Define function that iterates through observational dictionary tree
        and returns a dictionary of available data in the selected date
        range"""

        # create dictionary to store available observational data
        self.available_observation_data = {}

        # check end date is > start date, if not, return with no valid obs. files
        if self.selected_start_date >= self.selected_end_date:
            return

        # check start date and end date are both within if valid date range (19000101 - 20500101)
        # if not, return with no valid obs. files
        if (self.selected_start_date < 19000101) or \
            (self.selected_end_date < 19000101) or \
            (self.selected_start_date >= 20500101) or \
                (self.selected_end_date >= 20500101):
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
                        valid_species_files_yearmonths = [ym for ym in species_file_yearmonths if (ym >= self.selected_start_date_firstdayofmonth) & (ym < self.selected_end_date)]
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
                            self.available_observation_data[network][resolution][matrix][species] = \
                                valid_species_files_yearmonths

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

                # test first if interpolated directory exists before trying to get files from it
                # if it does not exit, continue
                if not os.path.exists(
                        '%s/%s/%s/%s/%s/%s/%s' % (self.exp_root, self.ghost_version, experiment, grid,
                                                  self.selected_resolution, self.selected_species,
                                                  self.selected_network)):
                    continue
                else:
                    # get all experiment netCDF files by experiment/grid/selected
                    # resolution/selected species/selected network
                    network_files = os.listdir(
                        '%s/%s/%s/%s/%s/%s/%s' % (exp_root, ghost_version, experiment, grid, self.selected_resolution,
                                                  self.selected_species, self.selected_network))
                    # get start YYYYMM yearmonths of data files
                    network_files_yearmonths = [int(f.split('_')[-1][:6]+'01') for f in network_files]
                    # limit data files to just those within date range
                    valid_network_files_yearmonths = \
                        [ym for ym in network_files_yearmonths if (ym >= self.selected_start_date_firstdayofmonth) &
                         (ym < self.selected_end_date)]

                    # if have some valid data files for experiment-grid, add experiment grid
                    # (with associated yearmonths) to dictionary
                    if len(valid_network_files_yearmonths) > 0:
                        self.available_experiment_data['%s-%s' % (experiment, grid)] = valid_network_files_yearmonths

        # get list of available experiment-grid names
        self.available_experiment_grids = np.array(sorted(list(self.available_experiment_data.keys())))


    def config_bar_params_change_handler(self, changed_param):
        """Define function which handles interactive updates of combo box fields"""

        if (changed_param != '') & (self.block_config_bar_handling_updates is False):

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
                self.selected_species = sorted(list(self.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()][self.cb_matrix.currentText()].keys()))[0]
            elif event_source == self.cb_species:
                self.selected_species = changed_param

            # set variable to check if date range changes
            self.date_range_has_changed = False

            # if have start date/end date have changed, make sure both have 8 characters (YYYYMMDD),
            # and are both numbers, before doing update of selection/fields
            if (event_source == self.le_start_date) or (event_source == self.le_end_date):
                valid_date = False
                selected_start_date = self.le_start_date.text()
                selected_end_date = self.le_end_date.text()
                if (len(selected_start_date) == 8) & (len(selected_end_date) == 8):
                    if (selected_start_date.isdigit() is True) & (selected_end_date.isdigit() is True):
                        self.date_range_has_changed = True
                        self.selected_start_date = int(selected_start_date)
                        self.selected_end_date = int(selected_end_date)
                        self.selected_start_date_firstdayofmonth = int(str(self.selected_start_date)[:6]+'01')
                    else:
                        return

            # update configuration bar fields
            self.update_configuration_bar_fields()

    # --------------------------------------------------------------------------------# 
    # --------------------------------------------------------------------------------# 
    # define functions which generate pop up configuration windows for some fields

    def handle_pop_up_experiments_window(self):
        # setup pop up window
        self.experiments_window = \
            PopUpWindow(window_type='EXPERIMENTS', window_titles=['Select Experiment/s'],
                        checkbox_labels=[self.available_experiment_grids], default_checkbox_selection=[[]],
                        selected_indices=self.selected_indices)

    def handle_pop_up_flags_window(self):
        # setup pop up window
        self.qa_window = \
            PopUpWindow(window_type='FLAGS',
                        window_titles=['Select standardised data reporter provided flags to filter by'],
                        checkbox_labels=[self.flag_names], default_checkbox_selection=[self.flag_default_inds],
                        selected_indices=self.selected_indices)

    def handle_pop_up_qa_window(self):
        # setup pop up window
        self.qa_window = \
            PopUpWindow(window_type='QA', window_titles=['Select standardised QA flags to filter by'],
                        checkbox_labels=[self.qa_names], default_checkbox_selection=[self.qa_default_inds],
                        selected_indices=self.selected_indices)

    def handle_pop_up_classifications_window(self):
        # setup pop up window
        self.qa_window = \
            PopUpWindow(window_type='CLASSIFICATIONS',
                        window_titles=['Select standardised classifications to retain',
                                       'Select standardised classifications to remove'],
                        checkbox_labels=[self.classification_names, self.classification_names],
                        default_checkbox_selection=[self.classification_default_inds_to_retain,
                                                    self.classification_default_inds_to_remove],
                        selected_indices=self.selected_indices)

    def handle_pop_up_methods_window(self):
        # only proceed if have some valid stations in memory
        if len(self.station_references) > 0:
            # setup pop up window
            self.qa_window = \
                PopUpWindow(window_type='METHODS',
                            window_titles=['Select standardised measurement methodologies to retain'],
                            checkbox_labels=[self.station_unique_methods],
                            default_checkbox_selection=[np.arange(len(self.station_unique_methods), dtype=np.int)],
                            selected_indices=self.selected_indices)

    # --------------------------------------------------------------------------------# 
    # --------------------------------------------------------------------------------# 
    def handle_data_selection_update(self):
        """Define function which handles update of data selection
        and MPL canvas upon pressing of READ button
        """

        # if have no data to read, then do not read any data
        if self.no_data_to_read is True:
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
        self.previous_active_qa_inds = self.active_qa_inds
        self.previous_active_flag_inds = self.active_flag_inds
        self.previous_active_classifications_to_retain_inds = self.active_classifications_to_retain_inds
        self.previous_active_classifications_to_remove_inds = self.active_classifications_to_remove_inds

        # set all currently selected variables as active variables
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

        # --------------------------------------------------------------------# 
        # determine what data (if any) needs to be read

        # set variables that inform what data needs to be read (set all initially as False)
        read_all = False
        read_left = False
        read_right = False
        cut_left = False
        cut_right = False

        # determine if any of the key variables have changed
        # (network, resolution, species, qa, flags, classifications_to_retain, classifications_to_remove)
        # if any have changed, observations and any selected experiments have to be re-read entirely
        if (self.active_network != self.previous_active_network) or \
                (self.active_resolution != self.previous_active_resolution) or \
                (self.active_species != self.previous_active_species) or \
                (np.array_equal(self.active_qa_inds, self.previous_active_qa_inds) is False) or \
                (np.array_equal(self.active_flag_inds, self.previous_active_flag_inds) is False) or \
                (np.array_equal(self.active_classifications_to_retain_inds,
                                self.previous_active_classifications_to_retain_inds) is False) or \
                (np.array_equal(self.active_classifications_to_remove_inds,
                                self.previous_active_classifications_to_remove_inds) is False):
            read_all = True
        # key variables have not changed, has start/end date?
        else:
            # determine if start date/end date have changed
            if (self.active_start_date != self.previous_active_start_date) or \
                    (self.active_end_date != self.previous_active_end_date):
                # if date range has changed then determine type of overlap with previous date range
                # no overlap (i.e. start date >= than previous end date, or end date <= than previous start date)?
                if (self.active_start_date >= self.previous_active_end_date) or \
                        (self.active_end_date <= self.previous_active_start_date):
                    read_all = True
                # data range fully inside previous data range (i.e. start date later and end date earlier)?
                elif (self.active_start_date > self.previous_active_start_date) & \
                        (self.active_end_date < self.previous_active_end_date):
                    cut_left = True
                    cut_right = True
                # need to read data on left edge and right edge of previous date range
                # (i.e. start date earlier and end date later)?
                elif (self.active_start_date < self.previous_active_start_date) & \
                        (self.active_end_date > self.previous_active_end_date):
                    read_left = True
                    read_right = True
                # need to read data on left edge and cut on right edge of previous date range
                # (i.e. start date earlier and end date earlier)?
                elif (self.active_start_date < self.previous_active_start_date) & \
                        (self.active_end_date < self.previous_active_end_date):
                    read_left = True
                    cut_right = True
                # need to cut data on left edge and read data on right edge of previous date range
                # (i.e. start date later and end date later)?
                elif (self.active_start_date > self.previous_active_start_date) & \
                        (self.active_end_date > self.previous_active_end_date):
                    cut_left = True
                    read_right = True
                # need to read data on left edge of previous date range (i.e. start date earlier)?
                elif self.active_start_date < self.previous_active_start_date:
                    read_left = True
                # need to read data on right edge of previous date range (i.e. end date later)?
                elif self.active_end_date > self.previous_active_end_date:
                    read_right = True
                # need to cut data on left edge of previous date range (i.e. start date later)?
                elif self.active_start_date > self.previous_active_start_date:
                    cut_left = True
                # need to cut data on right edge of previous date range (i.e. end date earlier)?
                elif self.active_end_date < self.previous_active_end_date:
                    cut_right = True

        # ---------------------------------------# 

        # determine if any of the active experiments have changed
        # remove experiments that are no longer selected from data_in_memory dictionary
        experiments_to_remove = [experiment for experiment in self.previous_active_experiment_grids
                                 if experiment not in self.active_experiment_grids]

        for experiment in experiments_to_remove:
            del self.data_in_memory[experiment]

        # any new experiments will need completely re-reading
        experiments_to_read = [experiment for experiment in self.active_experiment_grids
                               if experiment not in self.previous_active_experiment_grids]

        # has date range changed?
        if (read_all is True) or (read_left is True) or (read_right is True) or \
                (cut_left is True) or (cut_right is True):
            # set new active time array/unique station references/longitudes/latitudes
            # adjust data arrays to account for potential changing number of stations
            self.read_setup()

            # need to re-read all observations/experiments?
            if read_all:
                # reset data in memory dictionary
                self.data_in_memory = {}
                # read observations
                self.read_data('observations', self.active_start_date, self.active_end_date)
                # read selected experiments (iterate through)
                for data_label in self.active_experiment_grids:
                    self.read_data(data_label, self.active_start_date, self.active_end_date)
                    # if experiment in experiments_to_read list, remove it (as no longer need to read it)
                    if data_label in experiments_to_read:
                        experiments_to_read.remove(data_label)
            else:
                # if station references array has changed then as cutting/appending to existing data
                # need to rearrange existing data arrays accordingly
                if np.array_equal(self.previous_station_references, self.station_references) is False:
                    # get indices of stations in previous station references array in current station references array
                    old_station_inds = np.where(np.in1d(self.previous_station_references, self.station_references))[0]
                    # get indices of stations in current station references array
                    # that were in previous station references array
                    new_station_inds = np.where(np.in1d(self.station_references, self.previous_station_references))[0]

                    # iterate through all keys in data in memory dictionary
                    for data_label in list(self.data_in_memory.keys()):
                        # create new data array in shape of current station references array,
                        # putting the old data into new array in the correct positions
                        new_data_array = np.full((len(self.station_references),
                                                  len(self.previous_time_array)),
                                                 np.NaN, dtype=np.float32)
                        new_data_array[new_station_inds, :] = \
                            self.data_in_memory[data_label]['data'][old_station_inds, :]
                        # overwrite data array with reshaped version
                        self.data_in_memory[data_label]['data'] = new_data_array

            # need to cut edges?
            if (cut_left is True) or (cut_right is True):

                # set default edge limits as current edges
                left_edge_ind = 0
                right_edge_ind = len(self.previous_time_array)

                # need to cut on left data edge?
                if cut_left is True:
                    left_edge_ind = np.where(self.previous_time_array == self.time_array[0])[0][0]

                # need to cut on right data edge?
                if cut_right is True:
                    right_edge_ind = np.where(self.previous_time_array == self.time_array[-1])[0][0]+1

                # iterate through all keys in data in memory dictionary
                # and cut edges of the associated arrays appropriately
                for data_label in list(self.data_in_memory.keys()):
                    self.data_in_memory[data_label]['data'] = \
                        self.data_in_memory[data_label]['data'][:, left_edge_ind:right_edge_ind]

            # need to read on left edge?
            if read_left is True:
                # get n number of new elements on left edge
                n_new_left_inds = np.where(self.time_array == self.previous_time_array[0])[0][0]
                # iterate through all keys in data in memory dictionary and insert read data on
                # left edge of the associated arrays
                for data_label in list(self.data_in_memory.keys()):
                    # add space on left edge to insert new read data
                    self.data_in_memory[data_label]['data'] = np.insert(
                        self.data_in_memory[data_label]['data'], 0,
                        np.full((len(self.station_references), n_new_left_inds), np.NaN), axis=0)
                    self.read_data(data_label, self.active_start_date, self.previous_active_start_date)

            # need to read on right edge?
            if read_right is True:
                # get n number of new elements on right edge
                n_new_right_inds = \
                    (len(self.time_array) - 1) - np.where(self.time_array == self.previous_time_array[-1])[0][0]
                # iterate through all keys in data in memory dictionary
                # and insert read data on right edge of the associated arrays
                for data_label in list(self.data_in_memory.keys()):
                    self.data_in_memory[data_label]['data'] = \
                        np.append(self.data_in_memory[data_label]['data'],
                                  np.full((len(self.station_references),
                                           n_new_right_inds), np.NaN), axis=0)
                    self.read_data(data_label, self.previous_active_end_date, self.active_end_date)

        # if have new experiments to read, then read them now
        if len(experiments_to_read) > 0:
            for data_label in experiments_to_read:
                self.read_data(data_label, self.active_start_date, self.active_end_date)

        # --------------------------------------------------------------------# 
        # if species has changed, update default species specific lower/upper limits
        if self.active_species != self.previous_active_species:
            # update default lower/upper species specific limits and filter data outside limits
            species_lower_limit = np.float32(self.parameter_dictionary[self.active_species]['extreme_lower_limit'])
            species_upper_limit = np.float32(self.parameter_dictionary[self.active_species]['extreme_upper_limit'])
            # set default limits
            self.le_minimum_value.setText(str(species_lower_limit))
            self.le_maximum_value.setText(str(species_upper_limit))

        # --------------------------------------------------------------------# 
        # update dictionary of plotting parameters (colour and zorder etc.) for each data array
        self.update_plotting_parameters()

        # --------------------------------------------------------------------# 
        # run function to filter data outside lower/upper limits, not using desired measurement methods
        # and < desired minimum data availability
        self.mpl_canvas.handle_data_filter_update()

        # --------------------------------------------------------------------# 
        # update map z combobox fields based on data in memory

        # generate lists of basic and basis+bias statistics for using in the z statistic combobox
        basic_stats_dict = json.load(open('providentia/conf/basic_stats_dict.json'))
        self.basic_z_stats = np.array(list(
            OrderedDict(sorted(basic_stats_dict.items(), key=lambda x: x[1]['order'])).keys()))

        # load experiment bias dictionary from configuration
        expbias_dict = json.load(open('providentia/conf/experiment_bias_stats_dict.json'))
        self.basic_and_bias_z_stats = np.append(self.basic_z_stats, list(
            OrderedDict(sorted(expbias_dict.items(), key=lambda x: x[1]['order'])).keys()))

        # generate list of sorted z1/z2 data arrays names in memory,
        # putting observations before experiments and
        # empty string item as first element in z2 array list
        # (for changing from 'difference' statistics to 'absolute')
        if len(list(self.data_in_memory.keys())) == 1:
            self.z1_arrays = np.array(['observations'])
        else:
            data_array_labels = np.array(list(self.data_in_memory.keys()))
            self.z1_arrays = np.append(['observations'], np.delete(
                data_array_labels, np.where(data_array_labels == 'observations')))
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

        # get N time chunks between desired start date and end date to set time array
        if self.active_resolution == 'hourly':
            self.active_frequency_code = 'H'
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

        # get all relevant observational files
        file_root = '%s/%s/%s/%s/%s/%s_' % (obs_root, self.active_network, ghost_version,
                                            self.active_resolution, self.active_species, self.active_species)
        relevant_files = \
            sorted([file_root+str(yyyymm)[:6]+'.nc' for yyyymm in self.available_observation_data[
                self.active_network][self.active_resolution][self.active_matrix][self.active_species]])

        # redefine some key variables used by parallel netCDF reading functions
        self.selected_qa = self.qa_codes[self.active_qa_inds]
        self.selected_flags = self.flag_codes[self.active_flag_inds]
        self.selected_classifications_to_retain = self.classification_codes[self.active_classifications_to_retain_inds]
        self.selected_classifications_to_remove = self.classification_codes[self.active_classifications_to_remove_inds]

        # Iterate through all relevant observational files
        # and read station references/longitudes/latitudes
        # (either in serial/parallel)

        # define dictionary to store all read metadata (i.e. per file)
        all_read_metadata = {}

        # define dictionary with metadata variables to read, with associated data types
        metadata_dict = {'station_reference': np.object,
                         'station_name': np.object,
                         'latitude': np.float32,
                         'longitude': np.float32,
                         'measurement_altitude': np.float32,
                         'country': np.object,
                         'network': np.object,
                         'standardised_network_provided_area_classification': np.object,
                         'standardised_network_provided_station_classification': np.object,
                         'standardised_network_provided_main_emission_source': np.object,
                         'standardised_network_provided_land_use': np.object,
                         'standardised_network_provided_terrain': np.object,
                         'standardised_network_provided_measurement_scale': np.object,
                         'representative_radius': np.float32,
                         'GSFC_coastline_proximity': np.float32,
                         'primary_sampling_type': np.object,
                         'sample_preparation_types': np.object,
                         'measurement_methodology': np.object,
                         'measuring_instrument_name': np.object,
                         'measuring_instrument_sampling_type': np.object}

        # create list of metadata variables to read (make global)
        metadata_vars_to_read = list(metadata_dict.keys())

        # add all metadata variables to read to station metadata dictionary with associated empty
        # numpy arrays of the appropriate type
        for meta_var in metadata_vars_to_read:
            all_read_metadata[meta_var] = np.array([], dtype=metadata_dict[meta_var])

        # read serially
        if self.read_type == 'serial':
            # iterate through relevant files
            for relevant_file in relevant_files:
                file_metadata = read_netcdf_station((relevant_file, metadata_vars_to_read))
                for meta_var in metadata_vars_to_read:
                    all_read_metadata[meta_var] = np.append(all_read_metadata[meta_var], file_metadata[meta_var])

        # read in parallel
        elif self.read_type == 'parallel':

            # setup pool of N workers on N CPUs
            pool = multiprocessing.Pool(n_CPUs)

            # read netCDF files in parallel
            tuples_list = [(file_name, metadata_vars_to_read) for file_name in relevant_files]
            all_file_metadata = pool.map(read_netcdf_station, tuples_list)
            # will not submit more files to pool, so close access to it
            pool.close()
            # wait for worker processes to terminate before continuing
            pool.join()

            # iterate through relevant file data
            for file_metadata in all_file_metadata:
                for meta_var in metadata_vars_to_read:
                    all_read_metadata[meta_var] = np.append(all_read_metadata[meta_var], file_metadata[meta_var])

        # define dictionary to store read metadata, by station, by metadata variable
        self.station_metadata = {}

        # get unique sorted station references across all files (make global, and also add to self)
        self.station_references, unique_indices = np.unique(all_read_metadata['station_reference'], return_index=True)

        # get standard measurement methodology per station
        self.station_methods = np.array([ref.split('_')[-1].split('--')[0] for ref in self.station_references])

        # for latitude/longitude/measurement altitude/GSFC coastline proximity variables, set static metadata
        # variables (taking the first available value per station in the given time window)
        self.station_longitudes = all_read_metadata['longitude'][unique_indices]
        self.station_latitudes = all_read_metadata['latitude'][unique_indices]
        self.station_measurement_altitudes = all_read_metadata['measurement_altitude'][unique_indices]
        self.station_GSFC_coastline_proximities = all_read_metadata['GSFC_coastline_proximity'][unique_indices]

        # iterate through each unique station
        for station_reference in self.station_references:
            station_inds = np.where(all_read_metadata['station_reference'] == station_reference)[0]
            # create empty dictionary ro store station specific metadata
            self.station_metadata[station_reference] = {}
            # iterate through read metadata variables
            for meta_var in metadata_vars_to_read:
                # get all unique metadata for variable, by station
                unique_metadata = np.unique(all_read_metadata[meta_var][station_inds])
                # append all and unique metadata by station
                self.station_metadata[station_reference][meta_var] = \
                    {'all': all_read_metadata[meta_var][station_inds], 'unique': unique_metadata}

        # update measurement units for species (take standard units from parameter dictionary)
        self.measurement_units = self.parameter_dictionary[self.active_species]['standard_units']


    def read_data(self, data_label, start_date_to_read, end_date_to_read):
        """Function that handles reading of observational/experiment data"""

        # force garbage collection (to avoid memory issues)
        gc.collect()

        # set process_type variable to self (for access by parallel read function)
        # also get relevant file start dates
        if data_label == 'observations':
            self.process_type = 'observations'
            file_root = '%s/%s/%s/%s/%s/%s_' % (obs_root, self.active_network, ghost_version,
                                                self.active_resolution, self.active_species, self.active_species)
            relevant_file_start_dates = \
                sorted(self.available_observation_data[self.active_network]
                       [self.active_resolution][self.active_matrix][self.active_species])

        else:
            self.process_type = 'experiment'
            experiment_grid_split = data_label.split('-')
            active_experiment = experiment_grid_split[0]
            active_grid = experiment_grid_split[1]
            file_root = \
                '%s/%s/%s/%s/%s/%s/%s/%s_' % (exp_root, ghost_version, active_experiment,
                                              active_grid, self.active_resolution, self.active_species,
                                              self.active_network, self.active_species)
            relevant_file_start_dates = sorted(self.available_experiment_data[data_label])

        # create list of relevant files to read
        relevant_files = [file_root+str(yyyymm)[:6]+'.nc' for yyyymm in relevant_file_start_dates]

        # limit data files to required date range to read (i.e. taking care not to re-read what has already been read)
        first_valid_file_ind = bisect.bisect_right(relevant_file_start_dates, int(start_date_to_read))
        if first_valid_file_ind != 0:
            first_valid_file_ind = first_valid_file_ind - 1
        last_valid_file_ind = bisect.bisect_left(relevant_file_start_dates, int(end_date_to_read))
        if first_valid_file_ind == last_valid_file_ind:
            relevant_files = [relevant_files[first_valid_file_ind]]
        else:
            relevant_files = relevant_files[first_valid_file_ind:last_valid_file_ind]

        # check if data label in data in memory dictionary
        if data_label not in list(self.data_in_memory.keys()):
            # if not create empty array (filled with NaNs) to store species data and place it in the dictionary
            self.data_in_memory[data_label] = {'data': np.full((len(self.station_references),
                                                                len(self.time_array)), np.NaN, dtype=np.float32)}

            # if process_type is experiment, get experiment specific grid edges from first relevant file,
            # and save to data in memory dictionary
            if self.process_type == 'experiment':
                exp_nc_root = Dataset(relevant_files[0])
                self.data_in_memory[data_label]['grid_edge_longitude'] = exp_nc_root['grid_edge_longitude'][:]
                self.data_in_memory[data_label]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
                exp_nc_root.close()

        # iterate and read species data in all relevant netCDF files (either in serial/parallel)

        # read serially
        if self.read_type == 'serial':

            # iterate through relevant netCDF files
            for relevant_file in relevant_files:
                # create argument tuple of function
                tuple_arguments = relevant_file, self.time_array, self.station_references, \
                                  self.active_species, self.process_type,\
                                  self.selected_qa, self.selected_flags, \
                                  self.selected_classifications_to_retain,\
                                  self.selected_classifications_to_remove
                # read file
                file_data, time_indices, full_array_station_indices = read_netcdf_data(tuple_arguments)
                # place read data into big array as appropriate
                self.data_in_memory[data_label]['data'][full_array_station_indices[np.newaxis, :],
                                                        time_indices[:, np.newaxis]] = file_data

        # read in parallel
        elif self.read_type == 'parallel':

            # setup pool of N workers on N CPUs
            pool = multiprocessing.Pool(n_CPUs)

            # read netCDF files in parallel
            tuple_arguments = [(file_name, self.time_array, self.station_references, self.active_species,
                                self.process_type, self.selected_qa, self.selected_flags,
                                self.selected_classifications_to_retain,
                                self.selected_classifications_to_remove) for
                               file_name in relevant_files]
            all_file_data = pool.map(read_netcdf_data, tuple_arguments)
            # will not submit more files to pool, so close access to it
            pool.close()
            # wait for worker processes to terminate before continuing
            pool.join()

            # iterate through read file data and place data into data array as appropriate
            for file_data in all_file_data:
                self.data_in_memory[data_label]['data'][file_data[2][:, np.newaxis],
                                                        file_data[1][np.newaxis, :]] = file_data[0]

    def update_plotting_parameters(self):
        """Function that updates plotting parameters (colour
        and zorder) for each selected data array
        """

        # assign a colour/zorder to all selected data arrays

        # define observations colour to be 'black'
        self.data_in_memory['observations']['colour'] = 'black'
        # define zorder of observations to be 5
        self.data_in_memory['observations']['zorder'] = 5

        # generate a list of RGB tuples for number of experiments there are
        sns.reset_orig()
        clrs = sns.color_palette('husl', n_colors=len(list(self.data_in_memory.keys()))-1)

        # iterate through sorted experiment names, assigning each experiment a new RGB colour tuple, and zorder
        experiment_ind = 1
        for experiment in sorted(list(self.data_in_memory.keys())):
            if experiment != 'observations':
                # define colour for experiment
                self.data_in_memory[experiment]['colour'] = clrs[experiment_ind-1]
                # define zorder for experiment (obs zorder + experiment_ind)
                self.data_in_memory[experiment]['zorder'] = \
                    self.data_in_memory['observations']['zorder'] + experiment_ind
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

        vars(self).update({(k, self.parse_parameter(val)) for k, val in opts.items()})


# generate Providentia dashboard
def main(**kwargs):
    """Main function"""
    q_app = QtWidgets.QApplication(sys.argv)
    q_app.setStyle("Fusion")
    GenerateProvidentiaDashboard(**kwargs)
    sys.exit(q_app.exec_())
