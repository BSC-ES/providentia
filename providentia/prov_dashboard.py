""" Module which provides main window """
from .configuration import ProvConfiguration
from .init_standards import InitStandards
from .configuration import parse_path
from .config import split_options
from .reading import get_yearmonths_to_read
from .prov_canvas import MPLCanvas
from .toolbar import NavigationToolbar
from .toolbar import save_data, conf_dialogs
from .prov_dashboard_aux import ComboBox
from .prov_dashboard_aux import QVLine
from .prov_dashboard_aux import PopUpWindow
from .prov_dashboard_aux import formatting_dict
from .prov_dashboard_aux import set_formatting
from .prov_read import DataReader
from .prov_offline import ProvidentiaOffline
from providentia import aux

import copy
import datetime
import os
import os.path
import json
import sys
from glob import glob
from functools import partial
from collections import OrderedDict

from PyQt5 import QtCore, QtWidgets, QtGui
import numpy as np
from dateutil.relativedelta import relativedelta

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class ProvidentiaMainWindow(QtWidgets.QWidget, ProvConfiguration, InitStandards):
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
        dconf_path = (os.path.join(CURRENT_PATH, 'conf/default.conf'))

        #config = configparser.ConfigParser()
        #config.read(conf_to_load)
        #all_sections = config.sections()

        # update from config file (if available)
        #config and section defined 
        self.from_conf = False
        if ('config' in kwargs) and ('section' in kwargs):
            self.load_conf(section=kwargs['section'], fpath=kwargs['config'])
        #just config defined (and only 1 section in file)
        elif 'config' in kwargs: 
            self.load_conf(section=None, fpath=kwargs['config'])    
        #if no config has been loaded 
        if (not self.from_conf) & (os.path.isfile(dconf_path)):
            self.load_conf(section='default', fpath=dconf_path)
            self.from_conf = False
        # update from command line
        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in kwargs.items()})
        # arguments are only local
        self.main_window_geometry = None
        InitStandards.__init__(self, obs_root=self.obs_root,
                               ghost_version=self.ghost_version)

        # load necessary dictionaries
        self.basic_stats_dict = json.load(open(os.path.join(CURRENT_PATH,
                                                            'conf/basic_stats_dict.json')))
        self.expbias_dict = json.load(open(os.path.join(CURRENT_PATH,
                                                        'conf/experiment_bias_stats_dict.json')))
        # initialize DataReader
        self.datareader = DataReader(self)
        if self.offline:
            ProvidentiaOffline(self)
        else:
            self.init_ui()
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
        """Get current geometry of main Providentia window"""
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
        self.bu_read.setFixedWidth(40)
        self.bu_read.setStyleSheet("color: green;")
        self.bu_read.setToolTip('Read selected configuration of data into memory')
        self.ch_colocate = set_formatting(QtWidgets.QCheckBox("Colocate"), formatting_dict['checkbox_menu'])
        self.ch_colocate.setToolTip('Temporally colocate observational/experiment data')
        self.cb_network = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_network.setFixedWidth(95)
        self.cb_network.setToolTip('Select providing observational data network. '
                                   'Names starting with * indicate non-GHOST datasets')
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
        self.bu_rep.setToolTip('Select % desired representativity in data across '
                               'whole record and for specific temporal periods')
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
        self.lb_experiment_bias = set_formatting(QtWidgets.QLabel(self, text="Exp. Bias"),
                                                 formatting_dict['title_menu'])
        self.lb_experiment_bias.setToolTip('Set experiment bias statistic')
        self.cb_experiment_bias_type = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_experiment_bias_type.setFixedWidth(100)
        self.cb_experiment_bias_type.setToolTip('Select experiment bias type')
        self.cb_experiment_bias_stat = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.cb_experiment_bias_stat.setFixedWidth(100)
        self.cb_experiment_bias_stat.setToolTip('Select experiment bias statistic')
        self.vertical_splitter_4 = QVLine()
        self.vertical_splitter_4.setMaximumWidth(20)
        self.lb_station_selection = set_formatting(QtWidgets.QLabel(self, text="Site Select"),
                                                   formatting_dict['title_menu'])
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

        # setup pop-up window menu tree for flags
        self.flag_menu = {'window_title':'FLAGS', 'page_title':'Select standardised data reporter provided flags to filter by', 'checkboxes':{}}
        self.flag_menu['checkboxes']['labels'] = np.array(sorted(self.standard_data_flag_name_to_data_flag_code, key=self.standard_data_flag_name_to_data_flag_code.get))
        self.flag_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
        self.flag_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
        self.flag_menu['checkboxes']['map_vars'] = np.sort(list(self.standard_data_flag_name_to_data_flag_code.values()))
        self.flag_menu['select_buttons'] = ['all', 'clear', 'default']

        # setup pop-up window menu tree for qa
        self.qa_menu = {'window_title':'QA', 'page_title':'Select standardised quality assurance flags to filter by', 'checkboxes':{}}
        self.qa_menu['checkboxes']['labels'] = np.array(sorted(self.standard_QA_name_to_QA_code, key=self.standard_QA_name_to_QA_code.get))
        self.qa_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
        self.qa_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
        self.qa_menu['checkboxes']['map_vars'] = np.sort(list(self.standard_QA_name_to_QA_code.values()))
        self.qa_menu['select_buttons'] = ['all', 'clear', 'default']

        # setup pop-up window menu tree for experiments
        self.experiments_menu = {'window_title': 'EXPERIMENTS', 'page_title': 'Select Experiment/s',
                                 'checkboxes': {'labels': [],
                                                'keep_default': [],
                                                'keep_selected': [],
                                                'map_vars': [],
                                                'select_buttons': ['all', 'clear']}}

        # setup pop-up window menu tree for metadata
        self.metadata_types, self.metadata_menu = aux.init_metadata(self)
        # setup pop-up window menu tree for % data representativity
        self.representativity_menu = {'window_title': '% DATA REPRESENTATIVITY',
                                      'page_title': 'Select % Data Representativity Bounds',
                                      'rangeboxes': {'labels': [],
                                                     'tooltips': [],
                                                     'current_lower': []}}
        # setup pop-up window menu tree for data periods
        self.period_menu = {'window_title': 'DATA PERIOD', 'page_title': 'Select Data Periods',
                            'checkboxes': {'labels': [],
                                           'keep_selected': [],
                                           'remove_selected': []}}

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

        # update variable to inform plotting functions whether to use colocated data/or not
        check_state = self.ch_colocate.checkState()
        if check_state == QtCore.Qt.Checked:
            self.colocate_active = True
        else:
            self.colocate_active = False

        # if we're starting from a configuration file, read first the setup
        if self.from_conf: 
            self.handle_data_selection_update()
            # then see if we have fields that require to be se (meta, rep, period)
            aux.representativity_conf(self)
            if hasattr(self, 'period'):
                self.period_conf()
            # if there are there are metadata reported in configuratoin
            if set([m.lower() for m in self.metadata_vars_to_read]).intersection(vars(self).keys()):
                aux.meta_from_conf(self)
            # call function to apply changes (filter)
            self.mpl_canvas.handle_data_filter_update()

        # set finalised layout
        self.setLayout(parent_layout)
        # plot whole dashboard
        self.show()
        # maximise window to fit screen
        self.showMaximized()

    def period_conf(self):
        keeps, removes = split_options(self.period)
        self.period_menu['checkboxes']['keep_selected'] += keeps
        self.period_menu['checkboxes']['remove_selected'] += removes

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

        self.reading_nonghost = aux.check_for_ghost(self.selected_network)
        if self.reading_nonghost:
            self.disable_ghost_buttons()
        else:
            self.enable_ghost_buttons()

        # set some default configuration values when initialising config bar
        if self.config_bar_initialisation:
            # set initially selected/active start-end date as default 201601-201701
            self.le_start_date.setText(str(self.start_date))
            self.le_end_date.setText(str(self.end_date))
            self.selected_start_date = int(self.le_start_date.text())
            self.selected_end_date = int(self.le_end_date.text())
            self.selected_start_date_firstdayofmonth = int(str(self.selected_start_date)[:6]+'01')
            self.active_start_date = int(self.le_start_date.text())
            self.active_end_date = int(self.le_end_date.text())
            self.date_range_has_changed = False

            # set selected/active values of other fields to be initially None
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
            esarchive_files_empty = json.load(open(os.path.join(CURRENT_PATH, 'conf/esarchive_files.json')))
            # and merge to existing dict if we have the path
            if self.nonghost_root is not None:
                esarchive_files = self.get_esarchive_yearmonth(esarchive_files_empty)
                self.all_observation_data = {**self.all_observation_data, **esarchive_files}
            # create dictionary of observational data inside date range
            self.datareader.get_valid_obs_files_in_date_range(self.le_start_date.text(),
                                                              self.le_end_date.text())

            # check which flags to select, depending if we have conf file or no
            self.flag_menu['checkboxes']['remove_selected'] = aux.which_flags(self)

        # if date range has changed then update available observational data dictionary
        if self.date_range_has_changed:
            self.datareader.get_valid_obs_files_in_date_range(self.le_start_date.text(),
                                                              self.le_end_date.text())

        # initialise/update fields - maintain previously selected values wherever possible
        # clear fields
        self.cb_network.clear()
        self.cb_resolution.clear()
        self.cb_matrix.clear()
        self.cb_species.clear()

        # if have no available observational data, return from function, updating variable informing that have no data
        if len(self.datareader.available_observation_data) == 0:
            self.no_data_to_read = True
            # unset variable to allow interactive handling from now
            self.block_config_bar_handling_updates = False
            return
        else:
            self.no_data_to_read = False

        # update network field
        available_networks = list(self.datareader.available_observation_data.keys())
        self.cb_network.addItems(available_networks)
        if self.selected_network in available_networks:
            self.cb_network.setCurrentText(self.selected_network)
        else:
            self.selected_network = self.cb_network.currentText()

        # update resolution field
        available_resolutions = list(self.datareader.available_observation_data[self.cb_network.currentText()].keys())
        # manually force order of available resolutions
        resolution_order_dict = {'hourly': 1, '3hourly': 2, '6hourly': 3, 'hourly_instantaneous': 4,
                                 '3hourly_instantaneous': 5, '6hourly_instantaneous': 6,
                                 'daily': 7, 'monthly': 8}
        available_resolutions = sorted(available_resolutions, key=resolution_order_dict.__getitem__)
        self.cb_resolution.addItems(available_resolutions)
        if self.selected_resolution in available_resolutions:
            self.cb_resolution.setCurrentText(self.selected_resolution)
        else:
            self.selected_resolution = self.cb_resolution.currentText()

        # update matrix field
        available_matrices = sorted(
            self.datareader.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()])
        self.cb_matrix.addItems(available_matrices)
        if self.selected_matrix in available_matrices:
            self.cb_matrix.setCurrentText(self.selected_matrix)
        else:
            self.selected_matrix = self.cb_matrix.currentText()

        # update species field
        available_species = sorted(self.datareader.available_observation_data[self.cb_network.currentText()][
                                       self.cb_resolution.currentText()][self.cb_matrix.currentText()])
        self.cb_species.addItems(available_species)
        if self.selected_species in available_species:
            self.cb_species.setCurrentText(self.selected_species)
        else:
            self.selected_species = self.cb_species.currentText()

        # update available experiment data dictionary
        self.datareader.get_valid_experiment_files_in_date_range()
        self.datareader.get_valid_obs_files_in_date_range(self.le_start_date.text(),
                                                          self.le_end_date.text())
        # update selected indices for experiments -- keeping previously selected experiments if available
        # set selected indices as previously selected indices in current available list of experiments
        if self.config_bar_initialisation and hasattr(self, 'experiments'):
            conf_experiments = [exp.strip() for exp in self.experiments.split(",")]
            self.experiments_menu['checkboxes']['keep_selected'] = [experiment for experiment in conf_experiments
                                                                    if experiment in
                                                                    self.experiments_menu['checkboxes']['map_vars']]
        self.experiments_menu['checkboxes']['keep_selected'] = [previous_selected_experiment for
                                                                previous_selected_experiment in
                                                                self.experiments_menu['checkboxes']['keep_selected']
                                                                if previous_selected_experiment in
                                                                self.experiments_menu['checkboxes']['map_vars']]

        # since a selection has changed, update also the qa flags
        qa_to_select = aux.which_qa(self)  # first check which flags
        self.qa_menu['checkboxes']['remove_default'] = aux.which_qa(self, return_defaults=True)
        if self.config_bar_initialisation:
            self.qa_menu['checkboxes']['remove_selected'] = qa_to_select
        else:
            # if the selected species has specific qa flags, ensure that none of the
            # inapplicable is selected
            if self.selected_species in self.qa_exceptions:
                self.qa_menu['checkboxes']['remove_selected'] = list(set(
                    self.qa_menu['checkboxes']['remove_selected']) - set(self.qa_diff))

        # unset variable to allow interactive handling from now
        self.block_config_bar_handling_updates = False

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
                    self.datareader.available_observation_data[self.cb_network.currentText()][self.cb_resolution.currentText()][
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
            if aux.check_for_ghost(self.selected_network):
                self.disable_ghost_buttons()
            else:
                self.enable_ghost_buttons()

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
        previous_active_network = self.active_network
        previous_active_resolution = self.active_resolution
        previous_active_species = self.active_species
        previous_active_start_date = self.active_start_date
        previous_active_end_date = self.active_end_date
        previous_active_experiment_grids = self.active_experiment_grids
        previous_active_qa = self.active_qa
        previous_active_flags = self.active_flags

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
        if (self.active_network != previous_active_network) or (
                self.active_resolution != previous_active_resolution) or (
                self.active_species != previous_active_species) or (
                np.array_equal(self.active_qa, previous_active_qa) == False) or (
                np.array_equal(self.active_flags, previous_active_flags) == False):
            read_all = True
        # key variables have not changed, has start/end date?
        else:
            if (self.active_start_date != previous_active_start_date) or (
                    self.active_end_date != previous_active_end_date):
                # if date range has changed then determine type of overlap with previous date range
                # no overlap (i.e. start date >= than previous end date, or end date <= than previous start date)?
                if (self.active_start_date >= previous_active_end_date) or (
                        self.active_end_date <= previous_active_start_date):
                    read_all = True
                # data range fully inside previous data range (i.e. start date later and end date earlier)?
                elif (self.active_start_date > previous_active_start_date) & (
                        self.active_end_date < previous_active_end_date):
                    cut_left = True
                    cut_right = True
                # need to read data on left edge and right edge of previous date range
                # (i.e. start date earlier and end date later)?
                elif (self.active_start_date < previous_active_start_date) & (
                        self.active_end_date > previous_active_end_date):
                    read_left = True
                    read_right = True
                # need to read data on left edge and cut on right edge of previous date range
                # (i.e. start date earlier and end date earlier)?
                elif (self.active_start_date < previous_active_start_date) & (
                        self.active_end_date < previous_active_end_date):
                    read_left = True
                    cut_right = True
                # need to cut data on left edge and read data on right edge of previous date range
                # (i.e. start date later and end date later)?
                elif (self.active_start_date > previous_active_start_date) & (
                        self.active_end_date > previous_active_end_date):
                    cut_left = True
                    read_right = True
                # need to read data on left edge of previous date range (i.e. start date earlier)?
                elif self.active_start_date < previous_active_start_date:
                    read_left = True
                # need to read data on right edge of previous date range (i.e. end date later)?
                elif self.active_end_date > previous_active_end_date:
                    read_right = True
                # need to cut data on left edge of previous date range (i.e. start date later)?
                elif self.active_start_date > previous_active_start_date:
                    cut_left = True
                # need to cut data on right edge of previous date range (i.e. end date earlier)?
                elif self.active_end_date < previous_active_end_date:
                    cut_right = True

        # determine if any of the active experiments have changed
        # remove experiments that are no longer selected from data_in_memory dictionary
        experiments_to_remove = [experiment for experiment in previous_active_experiment_grids if
                                 experiment not in self.active_experiment_grids]
        for experiment in experiments_to_remove:
            del self.datareader.data_in_memory[experiment]
        # any new experiments will need completely re-reading
        experiments_to_read = [experiment for experiment in self.active_experiment_grids if
                               experiment not in previous_active_experiment_grids]

        # has date range changed?
        if read_all or read_left or read_right or cut_left or cut_right:

            # set new active time array/unique station references/longitudes/latitudes
            # adjust data arrays to account for potential changing number of stations
            self.datareader.read_setup(self.active_resolution, self.active_start_date,
                                       self.active_end_date, self.active_network,
                                       self.active_species, self.active_matrix)
                                       
            # clear all elements in dashboard if there are less than 2 timesteps
            if self.datareader.clear_canvas:
                self.mpl_canvas.update_MPL_canvas()
                self.active_network = None
                return

            # need to re-read all observations/experiments?
            if read_all:

                # reset data in memory dictionary
                self.datareader.reset_data_in_memory()

                if not self.reading_nonghost:
                    self.metadata_inds_to_fill = np.arange(len(self.relevant_yearmonths))
                # read observations
                self.datareader.read_data('observations', self.active_start_date,
                                          self.active_end_date, self.active_network,
                                          self.active_resolution, self.active_species,
                                          self.active_matrix)
                # read selected experiments (iterate through)
                for data_label in self.active_experiment_grids:
                    self.datareader.read_data(data_label, self.active_start_date,
                                              self.active_end_date, self.active_network,
                                              self.active_resolution, self.active_species,
                                              self.active_matrix)
                    # if experiment in experiments_to_read list, remove it (as no longer need to read it)
                    if data_label in experiments_to_read:
                        experiments_to_read.remove(data_label)

            else:
                
                # if station references array has changed then as cutting/appending to
                # existing data need to rearrange existing data arrays accordingly
                if not np.array_equal(self.previous_station_references, self.station_references):
                    # get indices of stations in previous station references array in current station references array
                    old_station_inds = np.where(np.in1d(self.previous_station_references,
                                                        self.station_references))[0]
                    # get indices of stations in current station references array
                    # that were in previous station references array
                    new_station_inds = np.where(np.in1d(self.station_references,
                                                        self.previous_station_references))[0]

                    new_metadata_array = np.full((len(self.station_references),
                                                  len(self.previous_relevant_yearmonths)),
                                                 np.NaN, dtype=self.metadata_dtype)

                    # self.datareader.metadata_in_memory = self.datareader.metadata_in_memory[:,:-1]
                    new_metadata_array[new_station_inds, :] = self.datareader.metadata_in_memory[old_station_inds, :]
                    self.datareader.metadata_in_memory = new_metadata_array

                    # iterate through all keys in data in memory dictionary
                    for data_label in list(self.datareader.data_in_memory.keys()):
                        # create new data array in shape of current station references array
                        if data_label == 'observations':
                            new_data_array = np.full((len(self.station_references),
                                                      len(self.previous_time_array)),
                                                     np.NaN, dtype=self.datareader.data_dtype)
                        else:
                            new_data_array = np.full((len(self.station_references),
                                                      len(self.previous_time_array)),
                                                     np.NaN, dtype=self.datareader.data_dtype[:1])
                        # put the old data into new array in the correct positions
                        new_data_array[new_station_inds, :] = self.datareader.data_in_memory[
                                                                  data_label][old_station_inds, :]
                        # overwrite data array with reshaped version
                        self.datareader.data_in_memory[data_label] = new_data_array
            
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
                                                                      int(str_previous_first_relevant_yearmonth[4:6]),
                                                                      1, 0, 0))
                    metadata_left_edge_ind = (monthly_relative_delta.years * 12) + monthly_relative_delta.months

                # need to cut on right data edge?
                if cut_right:
                    data_right_edge_ind = np.where(self.previous_time_array == self.time_array[-1])[0][0] + 1
                    str_last_relevant_yearmonth = str(self.relevant_yearmonths[-1])
                    str_previous_last_relevant_yearmonth = str(self.previous_relevant_yearmonths[-1])
                    monthly_relative_delta = relativedelta(
                        datetime.datetime(int(str_previous_last_relevant_yearmonth[:4]),
                                          int(str_previous_last_relevant_yearmonth[4:6]),
                                          1, 0, 0), datetime.datetime(int(str_last_relevant_yearmonth[:4]),
                                                                      int(str_last_relevant_yearmonth[4:6]), 1, 0, 0))
                    metadata_right_edge_ind = \
                        metadata_right_edge_ind - ((monthly_relative_delta.years * 12) + monthly_relative_delta.months)

                # do metadata array cut
                if metadata_left_edge_ind == metadata_right_edge_ind:
                    self.datareader.metadata_in_memory = self.datareader.metadata_in_memory[:, [metadata_left_edge_ind]]
                else:
                    self.datareader.metadata_in_memory = \
                        self.datareader.metadata_in_memory[:, metadata_left_edge_ind:metadata_right_edge_ind]

                # iterate through all keys in data in memory dictionary and
                # cut edges of the associated arrays appropriately
                for data_label in list(self.datareader.data_in_memory.keys()):
                    self.datareader.data_in_memory[data_label] = \
                        self.datareader.data_in_memory[data_label][:, data_left_edge_ind:data_right_edge_ind]

            # need to read on left edge?
            if read_left:

                # get n number of new elements on left edge
                if self.previous_time_array.size > 0:
                    n_new_left_data_inds = np.where(self.time_array == self.previous_time_array[0])[0][0]
                else:
                    n_new_left_data_inds = len(self.time_array)

                # get list of yearmonths to read
                yearmonths_to_read = get_yearmonths_to_read(self.relevant_yearmonths, self.active_start_date,
                                                            previous_active_start_date, self.active_resolution)
                                                        
                # check which yearmonths_to_read in previous matrix
                yearmonths_in_old_matrix = np.isin(yearmonths_to_read, self.previous_relevant_yearmonths)

                # get yearmonths not currently accounted for in matrix
                if isinstance(yearmonths_to_read, list):
                    yearmonths_to_read = np.asarray(yearmonths_to_read)
                new_yearmonths = yearmonths_to_read[~yearmonths_in_old_matrix]

                if new_yearmonths.size > 0:

                    self.metadata_inds_to_fill = np.arange(0, len(yearmonths_to_read))
                    self.datareader.metadata_in_memory = np.concatenate((np.full(
                        (len(self.station_references), len(new_yearmonths)), np.NaN, dtype=self.metadata_dtype),
                                                            self.datareader.metadata_in_memory), axis=1)

                    # iterate through all keys in data in memory dictionary and
                    # insert read data on left edge of the associated arrays
                    for data_label in list(self.datareader.data_in_memory.keys()):
                        
                        # add space on left edge to insert new read data
                        if data_label == 'observations':
                            self.datareader.data_in_memory[data_label] = np.concatenate((np.full(
                                (len(self.station_references), n_new_left_data_inds), np.NaN,
                                dtype=self.datareader.data_dtype), self.datareader.data_in_memory[data_label]), axis=1)
                        else:
                            self.datareader.data_in_memory[data_label] = np.concatenate((np.full(
                                (len(self.station_references), n_new_left_data_inds), np.NaN,
                                dtype=self.datareader.data_dtype[:1]), self.datareader.data_in_memory[data_label]), axis=1)
                        self.datareader.read_data(data_label, self.active_start_date, previous_active_start_date,
                                                self.active_network, self.active_resolution,
                                                self.active_species, self.active_matrix)

            # need to read on right edge?
            if read_right:
            
                # get n number of new elements on right edge
                if self.previous_time_array.size > 0:
                    n_new_right_data_inds = (len(self.time_array) - 1) - \
                                            np.where(self.time_array == self.previous_time_array[-1])[0][0]
                else:
                    n_new_right_data_inds = (len(self.time_array))
                # get list of yearmonths to read
                yearmonths_to_read = get_yearmonths_to_read(self.relevant_yearmonths, previous_active_end_date,
                                                            self.active_end_date, self.active_resolution)

                # check which yearmonths_to_read in previous matrix
                yearmonths_in_old_matrix = np.isin(yearmonths_to_read, self.previous_relevant_yearmonths)

                # get yearmonths not currently accounted for in matrix
                if isinstance(yearmonths_to_read, list):
                    yearmonths_to_read = np.asarray(yearmonths_to_read)
                new_yearmonths = yearmonths_to_read[~yearmonths_in_old_matrix]

                if new_yearmonths.size > 0:
                    
                    self.metadata_inds_to_fill = np.arange(-len(yearmonths_to_read), 0)
                    self.datareader.metadata_in_memory = np.concatenate((self.datareader.metadata_in_memory, np.full(
                        (len(self.station_references), len(new_yearmonths)), np.NaN, dtype=self.metadata_dtype)), axis=1)
                    
                    # iterate through all keys in data in memory dictionary and
                    # insert read data on right edge of the associated arrays
                    for data_label in list(self.datareader.data_in_memory.keys()):
                        if data_label == 'observations':
                            self.datareader.data_in_memory[data_label] = np.concatenate((self.datareader.data_in_memory[data_label], np.full(
                                (len(self.station_references), n_new_right_data_inds), np.NaN, dtype=self.datareader.data_dtype)), axis=1)
                        else:
                            self.datareader.data_in_memory[data_label] = np.concatenate((self.datareader.data_in_memory[data_label], np.full(
                                (len(self.station_references), n_new_right_data_inds), np.NaN, dtype=self.datareader.data_dtype[:1])),
                                                                            axis=1)
                        self.datareader.read_data(data_label, previous_active_end_date, self.active_end_date,
                                                  self.active_network, self.active_resolution,
                                                  self.active_species, self.active_matrix)

            # update menu object fields
            aux.update_metadata_fields(self)
            self.representativity_menu = aux.representativity_fields(self, self.active_resolution)
            aux.update_period_fields(self.active_resolution, self.period_menu)

        # if have new experiments to read, then read them now
        if len(experiments_to_read) > 0:
            for data_label in experiments_to_read:
                self.datareader.read_data(data_label, self.active_start_date, self.active_end_date,
                                          self.active_network, self.active_resolution,
                                          self.active_species, self.active_matrix)

        # if species has changed, update default species specific lower/upper limits
        if self.active_species != previous_active_species:
            # update default lower/upper species specific limits and filter data outside limits
            species_lower_limit, species_upper_limit = aux.which_bounds(self, self.active_species)
            # set default limits
            self.le_minimum_value.setText(str(species_lower_limit))
            self.le_maximum_value.setText(str(species_upper_limit))

        # --------------------------------------------------------------------#
        # update dictionary of plotting parameters (colour and zorder etc.) for each data array
        self.datareader.update_plotting_parameters()
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
        if len(list(self.datareader.data_in_memory.keys())) == 1:
            self.z1_arrays = np.array(['observations'])
        else:
            data_array_labels = np.array(list(self.datareader.data_in_memory.keys()))
            self.z1_arrays = np.append(['observations'],
                                       np.delete(data_array_labels, np.where(data_array_labels == 'observations')))
        self.z2_arrays = np.append([''], self.z1_arrays)

        # initialise map z statistic comboboxes
        self.mpl_canvas.handle_map_z_statistic_update()

        # update experiment bias combobox fields based on data in memory
        # if have no experiment data, all fields are empty
        if len(list(self.datareader.data_in_memory.keys())) == 1:
            self.experiment_bias_types = np.array([])
        # else, generate combobox lists
        else:
            # set all experiment bias types
            self.experiment_bias_types = np.array(['Aggregated'])

            # initialise experiment bias comboboxes
            self.mpl_canvas.handle_experiment_bias_update()

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

        # set mouse cursor to hourglass
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        # set rep fields to empty lists and initialize again
        self.representativity_menu['rangeboxes']['labels'] = []
        self.representativity_menu['rangeboxes']['current_lower'] = []
        aux.representativity_fields(self, self.selected_resolution)

        # set period fields to empty and initiliaze them
        self.period_menu['checkboxes']['keep_selected'] = []
        self.period_menu['checkboxes']['remove_selected'] = []
        aux.update_period_fields(self.active_resolution, self.period_menu)

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
        aux.update_metadata_fields(self)

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

    def load_conf(self, section=None, fpath=None):
        """ Load existing configurations from file. """

        from .config import read_conf

        if fpath is None:
            fpath = parse_path(self.config_dir, self.config_file)

        if not os.path.isfile(fpath):
            print(("Error %s" % fpath))
            return

        opts = read_conf(section, fpath)
        #if were unable to read file then return
        if opts is None:
            return

        self.opts = opts
        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in opts.items()})
        self.from_conf = True

    def disable_ghost_buttons(self):
        """Disable button related only to ghost data"""
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

    def get_esarchive_yearmonth(self, esfiles):
        """
        Returns the esarchive_files json filles with yearmonth
        from esarchive available files

        :esfiles: contains structure of esarchives to read
        :type esfiles: json
        """
        for n in esfiles:
            network = n[1:].lower()
            for r in esfiles[n]:
                resolution = r
                for d in esfiles[n][r]:
                    detail = d
                    for s in esfiles[n][r][d]:
                        species = s
                        path = "{}/{}/{}/{}/{}".format(self.nonghost_root, network,
                                                       detail, resolution, species)
                        if os.path.exists(path):
                            species_files = glob(path + '*/*_??????.nc')
                            species_files_yearmonths = [int(f.split('_')[-1][:6] + '01') for f in species_files]
                            esfiles[n][r][d][s] = species_files_yearmonths
        return esfiles


# generate Providentia dashboard
def main(**kwargs):
    """Main function"""
    q_app = QtWidgets.QApplication(sys.argv)
    q_app.setStyle("Fusion")
    ProvidentiaMainWindow(**kwargs)
    sys.exit(q_app.exec_())
