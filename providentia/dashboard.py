""" Module which provides main window """

from collections import OrderedDict
import copy
import datetime
from functools import partial
import json
import os
import sys
import time
from weakref import WeakKeyDictionary
import yaml

import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import numpy as np
from packaging.version import Version
import pandas as pd
from PyQt5 import QtCore, QtWidgets, QtGui

from .canvas import MPLCanvas
from .configuration import load_conf
from .configuration import ProvConfiguration
from .dashboard_elements import ComboBox, QVLine, InputDialog
from .dashboard_elements import set_formatting
from .dashboard_interactivity import HoverAnnotation
from .fields_menus import (init_experiments, init_flags, init_qa, update_qa, init_metadata, init_multispecies, init_period, 
                           init_representativity, metadata_conf, multispecies_conf, representativity_conf, period_conf, 
                           update_metadata_fields, update_period_fields, update_representativity_fields)
from .plot_aux import get_taylor_diagram_ghelper
from .plot_formatting import format_axis
from .pop_up_window import PopUpWindow
from .read import DataReader
from .read_aux import (check_for_ghost, get_default_qa, get_frequency_code, generate_file_trees, 
                       get_valid_experiments, get_valid_obs_files_in_date_range,
                       get_nonrelevant_temporal_resolutions, get_relevant_temporal_resolutions,
                       temporal_resolution_order_dict, get_lower_resolutions)
from .toolbar import NavigationToolbar
from .warnings_prv import show_message

from providentia.auxiliar import CURRENT_PATH, join, expand_plot_characteristics


# set font DPI for uniform dashboard appearance across systems
os.environ["QT_FONT_DPI"] = "96"
# enable high DPI pixmaps
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])

class ProvidentiaMainWindow(QtWidgets.QWidget):
    """ Class that generates Providentia dashboard. """

    # create signals that are fired upon resizing/moving of main Providentia window
    resized = QtCore.pyqtSignal()
    move = QtCore.pyqtSignal()

    def __init__(self, **kwargs):

        # allow access to methods of parent class QtWidgets.QWidget
        super(ProvidentiaMainWindow, self).__init__()

        # load statistical yamls
        self.basic_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/basic_stats.yaml')))
        self.expbias_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/experiment_bias_stats.yaml')))

        # load representativity information
        self.representativity_info = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/internal/representativity.yaml')))

        # save warnings that appear next in to show them after the UI is initialised
        self.delay = True
        self.delayed_warnings = []

        # initialise default configuration variables
        # modified by commandline arguments, if given
        self.provconf = ProvConfiguration(self, **kwargs)

        # update variables from config file (if available)
        self.from_conf = False
        self.current_config = {}

        if self.config != '': 
            read_conf = False
            if os.path.exists(self.config):
                read_conf = True
            else:
                if os.path.exists(join(self.config_dir, self.config)):
                    self.config = join(self.config_dir, self.config)
                    read_conf = True

            if read_conf:
                if 'section' in kwargs:
                    # config and section defined 
                    load_conf(self, fpath=self.config)
                    if kwargs['section'] in self.all_sections:
                        self.from_conf = True
                        self.current_config = self.sub_opts[kwargs['section']]
                        self.section = kwargs['section']
                    else:
                        error = 'Error: The section specified in the command line ({0}) does not exist.'.format(kwargs['section'])
                        tip = 'Tip: For subsections, add the name of the parent section followed by an interpunct (·) '
                        tip += 'before the subsection name (e.g. SECTIONA·Spain).'
                        error = error + '\n' + tip
                        self.logger.error(error)
                        sys.exit(1)
                else:
                    # config defined, section undefined
                    load_conf(self, fpath=self.config)    
                    all_sections = self.sub_opts.keys()

                    # if no parent section names are found throw an error
                    if len(all_sections) == 0:
                        error = "Error: No sections were found in the configuration file, make sure to name them using square brackets."
                        self.logger.error(error)
                        sys.exit(1)
                    # if there is only one section, do not ask the user    
                    elif len(all_sections) == 1:
                        okpressed = False
                        self.section = list(all_sections)[0]
                    # ask the user for the section
                    else:
                        title = 'Sections'
                        msg = 'Select section to load'
                        dialog = InputDialog(self, title, msg, all_sections)
                        self.section, okpressed = dialog.selected_option, dialog.okpressed

                    if okpressed or (len(all_sections) == 1):
                        self.from_conf = True
                        self.current_config = self.sub_opts[self.section]
            
            else:
                # have config, but path does not exist
                error = 'Error: The path to the configuration file specified in the command line does not exist.'
                self.logger.error(error)
                sys.exit(1)

        # set initial filter species
        self.previous_filter_species = {}
        self.filter_species = {}

        # set initial calibration factor
        self.previous_calibration_factor = {}
        self.calibration_factor = {}

        # update variables from defined config file
        if self.current_config:
            for k, val in self.current_config.items():
                setattr(self, k, self.provconf.parse_parameter(k, val))

        # now all variables have been parsed, check validity of those, throwing errors where necessary
        self.provconf.check_validity()

        # get operating system specific formatting
        if self.operating_system == 'Mac':
            self.formatting_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/internal/stylesheet_mac.yaml')))
        elif self.operating_system == 'Linux':
            self.formatting_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/internal/stylesheet_linux.yaml')))
        elif self.operating_system == 'Windows':
            self.formatting_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/internal/stylesheet_windows.yaml')))

        # load characteristics per plot type
        # check for self defined plot characteristics file
        if self.plot_characteristics_filename == '':
            self.plot_characteristics_filename = join(
                PROVIDENTIA_ROOT, 'settings/plot_characteristics.yaml')
        plot_characteristics = yaml.safe_load(open(self.plot_characteristics_filename))
        self.plot_characteristics_templates = expand_plot_characteristics(plot_characteristics, 'dashboard')

        # arguments are only local
        self.full_window_geometry = None

        # generate file trees if needed
        generate_file_trees(self)

        # initialise DataReader
        self.datareader = DataReader(self)

        # update map z combobox fields based on data in memory
        # generate lists of basic and basis+bias statistics for using in the z statistic combobox
        if not hasattr(self, 'basic_z_stats'):
            self.basic_z_stats = np.array(list(
                OrderedDict(sorted(self.basic_stats.items(), key=lambda x: x[1]['order'])).keys()))
        if not hasattr(self, 'basic_and_bias_z_stats'):
            self.basic_and_bias_z_stats = np.append(self.basic_z_stats, list(
                OrderedDict(sorted(self.expbias_stats.items(), key=lambda x: x[1]['order'])).keys()))

        # initialise UI
        self.init_ui(**kwargs)

        # setup callback events upon resizing/moving of Providentia window
        self.resized.connect(self.get_geometry)
        self.move.connect(self.get_geometry)
        
        # show delayed warnings
        time.sleep(0.1)
        self.delay = False
        for msg in self.delayed_warnings:
            show_message(self, msg)

    def resizeEvent(self, event):
        """ Function to overwrite default PyQt5 resizeEvent function --> for calling get_geometry. """

        self.resized.emit()

        return super(ProvidentiaMainWindow, self).resizeEvent(event)

    def moveEvent(self, event):
        """ Function to overwrite default PyQt5 moveEvent function --> for calling get_geometry. """
        
        self.move.emit()
        
        return super(ProvidentiaMainWindow, self).moveEvent(event)

    @QtCore.pyqtSlot()
    def get_geometry(self):
        """ Update current geometry of main Providentia window and buttons. """

        # get geometry of full window
        #self.full_window_geometry = copy.deepcopy(self.geometry())
        self.full_window_geometry = copy.deepcopy(self.frameGeometry())

        # update geometry of qt elements
        self.update_qt_elements_geometry(resize=True)

    def update_qt_elements_geometry(self, plot_types='ALL', positions = [1,2,3,4,5], resize=False):
        """ Update current geometry of canvas qt elements"""

        # get dashboard and canvas pixels
        full_window_width = self.full_window_geometry.width()
        full_window_height = self.full_window_geometry.height()
        canvas_width = self.mpl_canvas.frameGeometry().width()
        canvas_height = self.mpl_canvas.frameGeometry().height()
        header_height = full_window_height - canvas_height

        if plot_types == 'ALL':
            plot_types = [self.position_1, self.position_2, self.position_3, self.position_4, self.position_5]
            show_buttons = False
        else:
            show_buttons = True

        #iterate through active dashboard plots
        for position, plot_type in zip(positions, plot_types):

            # gather menu, save buttons and elements for plot type 
            for menu_button, save_button, element in zip(self.mpl_canvas.menu_buttons, 
                                                         self.mpl_canvas.save_buttons, 
                                                         self.mpl_canvas.elements):

                menu_plot_type = menu_button.objectName().split('_menu')[0]
                if plot_type in ['periodic-violin','fairmode-target','fairmode-statsummary']:
                    plot_type = plot_type.replace('-','_')

                # proceed once have objects for plot type
                if plot_type == menu_plot_type:
                    
                    # get position of menu button (set in 1848 x 1016 resolution)
                    x = self.mpl_canvas.plot_characteristics_templates['general']['settings_menu']['position_'
                        + str(position)]['x']
                    y = self.mpl_canvas.plot_characteristics_templates['general']['settings_menu']['position_'
                        + str(position)]['y']

                    # calculate proportional position for different screen resolution
                    x = int((x * canvas_width) / 1848)
                    y = int((y * canvas_height) / 1016)
                    
                    # get geometries (old and new)
                    old_button_geometry = QtCore.QRect(menu_button.x(), menu_button.y(), 18, 18)
                    new_button_geometry = QtCore.QRect(x, y, 18, 18)
                    
                    # apply new geometry to menu and save buttons
                    menu_button.setGeometry(new_button_geometry)
                    save_button.setGeometry(int(menu_button.x() - ((30 * canvas_width) / 1848)), int(menu_button.y()), 20, 20)

                    # show buttons if active
                    if show_buttons:
                        menu_button.show()
                        save_button.show()

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

                    # apply new geometry to layout button and canvas covers (if are resizing)
                    if resize:

                        # apply new geometry to full canvas covers
                        if position == 1:
                            self.mpl_canvas.canvas_cover.setGeometry(0, 0, canvas_width, canvas_height)

                        else:
                            # apply new geometry to layout button
                            cb_position = getattr(self, 'cb_position_{}'.format(position))
                            # layout selector at position 2 needs to be further from menu
                            if position == 2:
                                width_diff = 910
                            else:
                                width_diff = 560
                            height_diff = 1
                            new_x = int(menu_button.x() - ((width_diff * canvas_width) / 1848))
                            new_y = int(menu_button.y() + ((height_diff * canvas_height) / 1016))
                            cb_position.move(new_x, new_y)

                            # apply new geometry to partial canvas covers
                            if position == 2:
                                canvas_x = int(new_x - ((75 * canvas_width) / 1848))
                                self.mpl_canvas.top_right_canvas_cover.setGeometry(canvas_x, new_y,
                                                                                   canvas_width-canvas_x, 
                                                                                   canvas_height-new_y)
                            elif position == 3:
                                self.mpl_canvas.lower_canvas_cover.setGeometry(0, new_y, 
                                                                               canvas_width, 
                                                                               canvas_height-new_y)
              
                    break

    def init_ui(self, **kwargs):
        """ Initialise user interface. """

        self.logger.info("Starting Providentia dashboard...")

        # set window title
        self.window_title = "Providentia"
        self.setWindowTitle(self.window_title)

        # add logo as icon
        self.setWindowIcon(QtGui.QIcon(join(PROVIDENTIA_ROOT, 'assets/logo.png')))

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

        # add one more horizontal layout
        hbox = QtWidgets.QHBoxLayout()

        # define all configuration box objects (labels, comboboxes etc.)
        # data selection section
        self.lb_data_selection = set_formatting(QtWidgets.QLabel(self, text="Data Selection"),
                                                self.formatting_dict['menu_title'])
        self.lb_data_selection.setToolTip('Setup configuration of data to read into memory')
        self.cb_network = set_formatting(ComboBox(self), self.formatting_dict['menu_combobox'])
        self.cb_network.setToolTip('Select providing observational data network. '
                                   'Names starting with * indicate non-GHOST datasets')
        self.cb_resolution = set_formatting(ComboBox(self), self.formatting_dict['menu_combobox'])
        self.cb_resolution.setToolTip('Select temporal resolution of data')
        self.cb_matrix = set_formatting(ComboBox(self), self.formatting_dict['menu_combobox'])
        self.cb_matrix.setToolTip('Select data matrix')
        self.cb_species = set_formatting(ComboBox(self), self.formatting_dict['menu_combobox'])
        self.cb_species.setToolTip('Select species')
        self.le_start_date = set_formatting(QtWidgets.QLineEdit(self), self.formatting_dict['menu_lineedit'])
        self.le_start_date.setToolTip('Set data start date: YYYYMMDD')
        self.le_end_date = set_formatting(QtWidgets.QLineEdit(self), self.formatting_dict['menu_lineedit'])
        self.le_end_date.setToolTip('Set data end date: YYYYMMDD')
        self.bu_QA = set_formatting(QtWidgets.QPushButton('QA', self), self.formatting_dict['menu_button'])
        self.bu_QA.setToolTip('Select standardised quality assurance flags to filter by')
        self.bu_flags = set_formatting(QtWidgets.QPushButton('FLAGS', self), self.formatting_dict['menu_button'])
        self.bu_flags.setToolTip('Select standardised data reporter provided flags to filter by')
        self.bu_experiments = set_formatting(QtWidgets.QPushButton('EXPS', self), self.formatting_dict['menu_button'])
        self.bu_experiments.setToolTip('Select experiment/s data to read')
        self.bu_multispecies = set_formatting(QtWidgets.QPushButton('MULTI', self), self.formatting_dict['menu_button'])
        self.bu_multispecies.setToolTip('Select data to filter by')
        self.bu_read = set_formatting(QtWidgets.QPushButton('READ', self), self.formatting_dict['menu_button'],
                                      extra_arguments={'color': 'green'})
        self.bu_read.setToolTip('Read selected configuration of data into memory')
        self.vertical_splitter_1 = QVLine()
        self.vertical_splitter_1.setMaximumWidth(20)

        # filters section
        self.lb_data_filter = set_formatting(QtWidgets.QLabel(self, text="Filters"), self.formatting_dict['menu_title'])
        self.lb_data_filter.setToolTip('Select criteria to filter data by')
        self.bu_rep = set_formatting(QtWidgets.QPushButton('% REP', self), self.formatting_dict['menu_button'])
        self.bu_rep.setToolTip('Select % desired representativity in data across '
                               'whole record and for specific temporal periods')
        self.bu_meta = set_formatting(QtWidgets.QPushButton('META', self), self.formatting_dict['menu_button'])
        self.bu_meta.setToolTip('Select metadata to filter by')
        self.bu_reset = set_formatting(QtWidgets.QPushButton('RESET', self), self.formatting_dict['menu_button'],
                                       extra_arguments={'color': 'red'})
        self.bu_reset.setToolTip('Reset filter fields to initial values')
        self.bu_period = set_formatting(QtWidgets.QPushButton('PERIOD', self), self.formatting_dict['menu_button'])
        self.bu_period.setToolTip('Select data in specific periods')
        self.bu_filter = set_formatting(QtWidgets.QPushButton('FILTER', self), self.formatting_dict['menu_button'],
                                        extra_arguments={'color': 'blue'})
        self.bu_filter.setToolTip('Filter data')
        self.lb_data_bounds = set_formatting(QtWidgets.QLabel(self, text="Bounds"), self.formatting_dict['menu_label'])
        self.lb_data_bounds.setToolTip('Set lower/upper bounds of data')
        self.le_minimum_value = set_formatting(QtWidgets.QLineEdit(self), self.formatting_dict['menu_lineedit'])
        self.le_minimum_value.setToolTip('Set lower bound of data')
        self.le_maximum_value = set_formatting(QtWidgets.QLineEdit(self), self.formatting_dict['menu_lineedit'])
        self.le_maximum_value.setToolTip('Set upper bound of data')
        self.vertical_splitter_2 = QVLine()
        self.vertical_splitter_2.setMaximumWidth(20)

        # statistical calculation section
        self.lb_statistic = set_formatting(QtWidgets.QLabel(self, text="Statistics"),
                                             self.formatting_dict['menu_title'])
        self.lb_statistic.setToolTip('Select the type of statistical calculation')
        self.lb_statistic_mode = set_formatting(QtWidgets.QLabel(self, text="Mode"), self.formatting_dict['menu_label'])
        self.lb_statistic_mode.setToolTip('Select statistical calculation mode')
        self.cb_statistic_mode = set_formatting(ComboBox(self), self.formatting_dict['menu_combobox'])
        self.cb_statistic_mode.setToolTip('Select statistical calculation mode')
        self.lb_statistic_aggregation = set_formatting(QtWidgets.QLabel(self, text="Aggregation"), self.formatting_dict['menu_label'])
        self.lb_statistic_aggregation.setToolTip('Select statistic for spatial aggregation')
        self.cb_statistic_aggregation = set_formatting(ComboBox(self), self.formatting_dict['menu_combobox'])
        self.cb_statistic_aggregation.setToolTip('Select statistic for spatial aggregation')
        self.vertical_splitter_3 = QVLine()
        self.vertical_splitter_3.setMaximumWidth(20)
        
        # colocation section
        self.lb_colocate = set_formatting(QtWidgets.QLabel(self, text="Colocation"), self.formatting_dict['menu_title'])
        self.lb_colocate.setToolTip('Set colocation')
        self.ch_colocate = set_formatting(QtWidgets.QCheckBox("Temporal"), self.formatting_dict['menu_checkbox'])
        self.ch_colocate.setToolTip('Temporally colocate observational/experiment data')
        self.vertical_splitter_4 = QVLine()
        self.vertical_splitter_4.setMaximumWidth(20)

        # resampling section
        self.lb_resampling = set_formatting(QtWidgets.QLabel(self, text="Resampling"), self.formatting_dict['menu_title'])
        self.lb_resampling.setToolTip('Set resampling options')
        self.cb_resampling_resolution = set_formatting(ComboBox(self), self.formatting_dict['menu_combobox'])
        self.cb_resampling_resolution.setToolTip('Select temporal resolution to resample the data to')
        self.vertical_splitter_5 = QVLine()
        self.vertical_splitter_5.setMaximumWidth(20)

        # station selection section
        self.lb_station_selection = set_formatting(QtWidgets.QLabel(self, text="Site Selection"),
                                                   self.formatting_dict['menu_title'])
        self.lb_station_selection.setToolTip('Select stations')
        self.ch_select_all = set_formatting(QtWidgets.QCheckBox("All"), self.formatting_dict['menu_checkbox'])
        self.ch_select_all.setToolTip('Select all stations')
        self.ch_intersect = set_formatting(QtWidgets.QCheckBox("Intersect"), self.formatting_dict['menu_checkbox'])
        self.ch_intersect.setToolTip('Select stations that intersect with all loaded model domains')
        self.ch_extent = set_formatting(QtWidgets.QCheckBox("Extent"), self.formatting_dict['menu_checkbox'])
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
        config_bar.addWidget(self.lb_statistic, 0, 10, 1, 2, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.lb_statistic_mode, 1, 10, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.lb_statistic_aggregation, 2, 10, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_statistic_mode, 1, 11, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_statistic_aggregation, 2, 11, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.vertical_splitter_3, 0, 12, 4, 1, QtCore.Qt.AlignLeft)

        # colocation section
        config_bar.addWidget(self.lb_colocate, 0, 13, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_colocate, 1, 13, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.vertical_splitter_4, 0, 14, 4, 1, QtCore.Qt.AlignLeft)

        # resampling section
        config_bar.addWidget(self.lb_resampling, 0, 15, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.cb_resampling_resolution, 1, 15, 1, 1, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.vertical_splitter_5, 0, 16, 4, 1, QtCore.Qt.AlignLeft)

        # station selection section
        config_bar.addWidget(self.lb_station_selection, 0, 17, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_select_all, 1, 17, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_intersect, 2, 17, QtCore.Qt.AlignLeft)
        config_bar.addWidget(self.ch_extent, 3, 17, QtCore.Qt.AlignLeft)

        # enable dynamic updating of specific configuration bar fields 
        self.cb_network.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_resolution.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_matrix.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_species.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.le_start_date.textChanged.connect(self.handle_config_bar_params_change)
        self.le_end_date.textChanged.connect(self.handle_config_bar_params_change)
        self.cb_statistic_mode.currentTextChanged.connect(self.handle_config_bar_params_change)
        self.cb_resampling_resolution.currentTextChanged.connect(self.handle_config_bar_params_change)

        # setup pop-up window menu tree for flags, qa, experiments, 
        # % data representativity, data periods and metadata
        init_flags(self)
        init_qa(self)
        init_experiments(self)
        init_multispecies(self)
        init_representativity(self)
        init_period(self)
        self.metadata_vars_to_read = []
        init_metadata(self)

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

        # launch with configuration file?
        if self.from_conf: 
            
            # read
            self.handle_data_selection_update()

            # set filtered multispecies if any
            multispecies_conf(self)

            # set fields available for filtering
            representativity_conf(self)
            period_conf(self)
            metadata_conf(self)
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

        # enable statistical calculation mode changing
        self.cb_statistic_mode.currentTextChanged.connect(self.mpl_canvas.handle_statistic_mode_update)

        # enable statistical aggregation update 
        self.cb_statistic_aggregation.currentTextChanged.connect(self.mpl_canvas.handle_statistic_aggregation_update)

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

        # set minimum height to avoid error 'box_aspect' and 'fig_aspect' must be positive
        self.setMinimumHeight(650)

        # show dashboard. How to do this is different per system 
        if self.operating_system == 'Mac':
            self.showMaximized()
            self.get_geometry()
        elif self.operating_system == 'Linux':
            self.show()
            self.showMaximized()
        elif self.operating_system == 'Windows':
            self.show()
            self.showMaximized()

    @QtCore.pyqtSlot(dict)
    def generate_pop_up_window(self, menu_root):
        """ Generate pop up window. """
        
        self.pop_up_window = PopUpWindow(self, menu_root, [], self.full_window_geometry)

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
            self.data_labels_raw = None

            # set initial station references to be empty dict
            self.station_references = {}

            # set initial selected config variables as set .conf files or defaults
            self.selected_network = copy.deepcopy(self.network[0])
            self.selected_resolution = copy.deepcopy(self.resolution)
            self.selected_matrix = self.parameter_dictionary[self.species[0]]['matrix']
            self.selected_species = copy.deepcopy(self.species[0])
            self.selected_resampling_resolution = copy.deepcopy(self.resampling_resolution)
            self.selected_statistic_mode = copy.deepcopy(self.statistic_mode)
            self.selected_statistic_aggregation = copy.deepcopy(self.statistic_aggregation)
            self.selected_periodic_statistic_aggregation = copy.deepcopy(self.periodic_statistic_aggregation)
            self.selected_periodic_statistic_mode = copy.deepcopy(self.periodic_statistic_mode)
            self.selected_timeseries_statistic_aggregation = copy.deepcopy(self.timeseries_statistic_aggregation)

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
            get_valid_obs_files_in_date_range(self, self.le_start_date.text(), self.le_end_date.text())

            # update qa / flags checkboxes 
            self.flag_menu['checkboxes']['remove_selected'] = copy.deepcopy(self.flags)
            self.qa_menu['checkboxes']['remove_selected'] = copy.deepcopy(self.qa_per_species[self.selected_species])
             
        # if date range has changed then update available observational data dictionary
        if self.date_range_has_changed:
            get_valid_obs_files_in_date_range(self, self.le_start_date.text(), self.le_end_date.text())

        # initialise/update fields - maintain previously selected values wherever possible
        # clear fields
        self.cb_network.clear()
        self.cb_resolution.clear()
        self.cb_matrix.clear()
        self.cb_species.clear()
        self.cb_resampling_resolution.clear()
        self.cb_statistic_mode.clear()
        self.cb_statistic_aggregation.clear()
        self.mpl_canvas.statsummary_periodic_aggregation.clear()
        self.mpl_canvas.statsummary_periodic_mode.clear()
        self.mpl_canvas.timeseries_stat.clear()

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

        # turn off some features if using non-GHOST data (on for ACTRIS)
        if check_for_ghost(self.selected_network) or self.selected_network == 'actris/actris':
            self.enable_ghost_buttons()
        else:
            self.disable_ghost_buttons()
        
        # update resolution field
        available_resolutions = list(self.available_observation_data[self.selected_network].keys())
        # set order of available resolutions
        available_resolutions = sorted(available_resolutions, key=temporal_resolution_order_dict().__getitem__)
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

        # update statistic mode field
        available_statistic_modes = ['Flattened', 'Spatial|Temporal', 'Temporal|Spatial']
        self.cb_statistic_mode.addItems(available_statistic_modes)
        if self.selected_statistic_mode in available_statistic_modes:
            self.cb_statistic_mode.setCurrentText(self.selected_statistic_mode)
        else:
            self.selected_statistic_mode = self.cb_statistic_mode.currentText()

        # update statistic aggregation field
        if self.selected_statistic_mode == 'Flattened':
            available_aggregation_statistics = []
        else:    
            available_aggregation_statistics = ['Mean', 'Median', 'p1', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95', 'p99']
            if self.selected_statistic_aggregation == '':
                self.selected_statistic_aggregation = available_aggregation_statistics[1]
        self.cb_statistic_aggregation.addItems(available_aggregation_statistics)
        if self.selected_statistic_aggregation in available_aggregation_statistics:
            self.cb_statistic_aggregation.setCurrentText(self.selected_statistic_aggregation)
        else:
            self.selected_statistic_aggregation = self.cb_statistic_aggregation.currentText()

        # update statsummary periodic statistic aggregation field
        available_periodic_statistics = ['Mean', 'Median', 'p1', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95', 'p99']
        self.mpl_canvas.statsummary_periodic_aggregation.addItems(available_periodic_statistics)
        if self.selected_periodic_statistic_aggregation in available_periodic_statistics:
            self.mpl_canvas.statsummary_periodic_aggregation.setCurrentText(self.selected_periodic_statistic_aggregation)
        else:
            self.selected_periodic_statistic_aggregation = self.mpl_canvas.statsummary_periodic_aggregation.currentText()
        
        # update statsummary periodic statistic mode field
        available_periodic_modes = ['Independent', 'Cycle']
        self.mpl_canvas.statsummary_periodic_mode.addItems(available_periodic_modes)
        if self.selected_periodic_statistic_mode in available_periodic_modes:
            self.mpl_canvas.statsummary_periodic_mode.setCurrentText(self.selected_periodic_statistic_mode)
        else:
            self.selected_periodic_statistic_mode = self.mpl_canvas.statsummary_periodic_mode.currentText()

        # update timeseries statistic field
        available_timeseries_statistics = ['Mean', 'Median', 'p1', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95', 'p99']
        self.mpl_canvas.timeseries_stat.addItems(available_timeseries_statistics)
        if self.selected_timeseries_statistic_aggregation in available_timeseries_statistics:
            self.mpl_canvas.timeseries_stat.setCurrentText(self.selected_timeseries_statistic_aggregation)
        else:
            self.selected_timeseries_statistic_aggregation = self.mpl_canvas.timeseries_stat.currentText()

        # get available resampling resolutions
        available_resampling_resolutions = get_lower_resolutions(self.selected_resolution)

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
        get_valid_experiments(self, self.le_start_date.text(), self.le_end_date.text(), self.selected_resolution,
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
        previous_default_qa = copy.deepcopy(self.qa_menu['checkboxes']['remove_default']) 
        self.qa_menu['checkboxes']['remove_default'] = default_qa

        # update selected qa if previous selected qa was default (to new default)
        if set(self.qa_menu['checkboxes']['remove_selected']) == set(previous_default_qa):
            self.qa_menu['checkboxes']['remove_selected'] = default_qa

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
    
        # remove plot types that need active temporal colocation and experiments data
        if Version(matplotlib.__version__) < Version("3.8"):
            if 'taylor' in canvas_instance.layout_options:
                canvas_instance.layout_options.remove('taylor')
        else:
            for plot_type in ['scatter', 'taylor', 'fairmode-target', 'fairmode-statsummary']:
                if ((not self.temporal_colocation) 
                    or ((self.temporal_colocation) and (len(self.experiments) == 0))): 
                    if plot_type in canvas_instance.layout_options:
                        canvas_instance.layout_options.remove(plot_type)
                else:
                    if plot_type not in canvas_instance.layout_options:
                        canvas_instance.layout_options.append(plot_type)       

            if 'fairmode-target' in canvas_instance.layout_options and hasattr(self, 'selected_species'):
                if self.selected_species not in ['sconco3', 'sconcno2', 'pm10', 'pm2p5']:
                    canvas_instance.layout_options.remove('fairmode-target')   
                if ((self.selected_species in ['sconco3', 'sconcno2'] and self.selected_resolution != 'hourly') 
                    or (self.selected_species in ['pm10', 'pm2p5'] and (self.selected_resolution not in ['hourly', 'daily']))):
                    canvas_instance.layout_options.remove('fairmode-target')
            
        # order alphabetically
        layout_options = sorted(canvas_instance.layout_options)

        # update position 2 in layout
        self.cb_position_2.clear()
        self.cb_position_2.addItems(layout_options)
        if self.position_2 in layout_options:
            self.cb_position_2.setCurrentText(self.position_2)
        else:
            self.block_config_bar_handling_updates = False
            self.cb_position_2.setCurrentText('None')
            self.block_config_bar_handling_updates = True

        # update position 3 in layout
        self.cb_position_3.clear()
        self.cb_position_3.addItems(layout_options)
        if self.position_3 in layout_options:
            self.cb_position_3.setCurrentText(self.position_3)
        else:
            self.block_config_bar_handling_updates = False
            self.cb_position_3.setCurrentText('None')
            self.block_config_bar_handling_updates = True

        # update position 4 in layout
        self.cb_position_4.clear()
        self.cb_position_4.addItems(layout_options)
        if self.position_4 in layout_options:
            self.cb_position_4.setCurrentText(self.position_4)
        else:
            self.block_config_bar_handling_updates = False
            self.cb_position_4.setCurrentText('None')
            self.block_config_bar_handling_updates = True

        # update position 5 in layout
        self.cb_position_5.clear()
        self.cb_position_5.addItems(layout_options)
        if self.position_5 in layout_options:
            self.cb_position_5.setCurrentText(self.position_5)
        else:
            self.block_config_bar_handling_updates = False
            self.cb_position_5.setCurrentText('None')
            self.block_config_bar_handling_updates = True

        # unset variable to allow interactive handling from now
        self.block_config_bar_handling_updates = False

    @QtCore.pyqtSlot(str)
    def handle_config_bar_params_change(self, changed_param):
        """ Function which handles interactive updates of combo box fields. """

        if (changed_param != '') & (not self.block_config_bar_handling_updates):

            # get event origin source
            event_source = self.sender()

            # if network, resolution, matrix, species, aggregation mode or resampling resolution have changed 
            # then alter respective current selection for the 
            # changed param
            if event_source == self.cb_network:
                self.selected_network = changed_param
                # ensure that QA defaults have been updated if network has changed (i.e. to or from ACTRIS)
                update_qa(self)

            elif event_source == self.cb_resolution:
                self.selected_resolution = changed_param

            elif event_source == self.cb_matrix:
                self.selected_matrix = changed_param
                self.selected_species = sorted(list(
                    self.available_observation_data[self.selected_network][self.selected_resolution][self.selected_matrix].keys()))[0]
            
            elif event_source == self.cb_species:
                self.selected_species = changed_param
                if ((self.selected_species not in ['sconco3', 'sconcno2', 'pm10', 'pm2p5'])
                    or (self.selected_species in ['sconco3', 'sconcno2'] and self.selected_resolution != 'hourly')
                    or (self.selected_species in ['pm10', 'pm2p5'] and self.selected_resolution not in ['hourly', 'daily'])):
                    for plot_type in ['fairmode-target', 'fairmode-statsummary']:
                        if plot_type in self.mpl_canvas.layout_options:
                            self.mpl_canvas.layout_options.remove(plot_type)
                            self.update_layout_fields(self.mpl_canvas)

            elif event_source == self.cb_statistic_mode:
                self.selected_statistic_mode = changed_param

            elif event_source == self.cb_resampling_resolution:
                self.selected_resampling_resolution = changed_param

            # set variable to check if date range changes
            self.date_range_has_changed = False

            # check if start date/end date have changed
            if (event_source == self.le_start_date) or (event_source == self.le_end_date):
                self.date_range_has_changed = True

            # initalise multispecies tab if network, resolution, matrix or species change
            if event_source in [self.cb_network, self.cb_resolution, self.cb_matrix, self.cb_species]:
                init_multispecies(self)

            # update configuration bar fields
            self.update_configuration_bar_fields()

    def handle_layout_update(self, changed_plot_type, sender=None):
        """ Function which handles update of layout. """
        
        if (changed_plot_type != '') & (not self.block_config_bar_handling_updates):
            
            # get event origin source if not given
            if sender is not None:
                if sender == 2:
                    event_source = self.cb_position_2
                elif sender == 3:
                    event_source = self.cb_position_3    
                elif sender == 4:
                    event_source = self.cb_position_4
                elif sender == 5:
                    event_source = self.cb_position_5
            else:
                event_source = self.sender()

            # update selected station plots, avoiding duplicates
            if event_source == self.cb_position_2:
                previous_plot_type = copy.deepcopy(self.position_2)
                self.position_2 = copy.deepcopy(changed_plot_type)
                changed_position = 2
                
            elif event_source == self.cb_position_3:
                previous_plot_type = copy.deepcopy(self.position_3)
                self.position_3 = copy.deepcopy(changed_plot_type)
                changed_position = 3

            elif event_source == self.cb_position_4:
                previous_plot_type = copy.deepcopy(self.position_4)
                self.position_4 = copy.deepcopy(changed_plot_type)
                changed_position = 4

            elif event_source == self.cb_position_5:
                previous_plot_type = copy.deepcopy(self.position_5)
                self.position_5 = copy.deepcopy(changed_plot_type)
                changed_position = 5

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
                elif isinstance(ax, list):
                    for sub_ax in ax:
                        sub_ax.remove()
                else:
                    ax.remove()
                self.active_dashboard_plots.remove(previous_plot_type)

                # hide qt elements for previous plot type
                for menu_button, save_button, element in zip(self.mpl_canvas.menu_buttons, 
                                                             self.mpl_canvas.save_buttons, 
                                                             self.mpl_canvas.elements):

                    menu_plot_type = menu_button.objectName().split('_menu')[0]
                    if menu_plot_type in ['periodic_violin','fairmode_target','fairmode_statsummary']:
                        menu_plot_type = menu_plot_type.replace('_','-')

                    if previous_plot_type == menu_plot_type:
                        menu_button.hide()
                        save_button.hide()
                        if previous_plot_type in ['periodic-violin','fairmode-target','fairmode-statsummary']:
                            previous_plot_type = previous_plot_type.replace('-','_')
                        for element in getattr(self.mpl_canvas, previous_plot_type + '_elements'):
                            if isinstance(element, dict):
                                for sub_element in element.values():
                                    sub_element.hide()
                            else:
                                element.hide()
                        break

            # if changed_plot_type already axis on another axis then remove those axis elements
            if (changed_plot_type in self.active_dashboard_plots) & (changed_plot_type in self.mpl_canvas.plot_axes):
                ax = self.mpl_canvas.plot_axes[changed_plot_type]
                self.mpl_canvas.remove_axis_elements(ax, changed_plot_type)
                if isinstance(ax, dict):
                    for sub_ax in ax.values():
                        sub_ax.remove()
                elif isinstance(ax, list):
                    for sub_ax in ax:
                        sub_ax.remove()
                else:
                    ax.remove()
            # otherwise add plot_type to active_dashboard_plots
            elif changed_plot_type != 'None': 
                self.active_dashboard_plots.append(changed_plot_type)

            # update plot axis for new plot type
            self.update_plot_axis(self.mpl_canvas, event_source, changed_plot_type)

            # update plot if changed_plot_type != 'None'
            if changed_plot_type != 'None':

                # format axis                
                format_axis(self.mpl_canvas, self.mpl_canvas.read_instance, 
                            self.mpl_canvas.plot_axes[changed_plot_type], 
                            changed_plot_type, self.mpl_canvas.plot_characteristics[changed_plot_type],
                            map_extent=self.map_extent)
                
                # make plot
                self.mpl_canvas.update_associated_active_dashboard_plot(changed_plot_type)

                # update qt elements geometry for changed plot type
                self.update_qt_elements_geometry(plot_types=[changed_plot_type], positions=[changed_position])

            # update layout fields
            self.update_layout_fields(self.mpl_canvas)

            # draw changes
            self.mpl_canvas.figure.canvas.draw_idle()

        return None

    def update_plot_axis(self, canvas_instance, changed_position, changed_plot_type):
        """ Update plot axis position from layout options."""

        # since we have no data, we need to use a random reference to initialize the polar axes
        # the axis will be updated when we create the taylor diagram
        if changed_plot_type == 'taylor':
            reference_stddev = 7.5
            plot_characteristics = canvas_instance.plot_characteristics['taylor']
            ghelper = get_taylor_diagram_ghelper(reference_stddev, plot_characteristics)
        elif changed_plot_type == "fairmode-statsummary":
            # get number of rows and columns
            ncols = 4
            nrows = 8 if self.species[0] in ["sconco3", "sconcno2", "pm10"] else 7

        # position 2 (top right)
        if changed_position == self.cb_position_2 or changed_position == 2:
            if changed_plot_type in ['periodic', 'periodic-violin']:
                canvas_instance.plot_axes[changed_plot_type] = {}
                canvas_instance.plot_axes[changed_plot_type]['hour'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((15, 53), rowspan=15, colspan=49))
                canvas_instance.plot_axes[changed_plot_type]['dayofweek'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((34, 83), rowspan=15, colspan=19))
                canvas_instance.plot_axes[changed_plot_type]['month'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((34, 53), rowspan=15, colspan=28))
            elif changed_plot_type == 'taylor':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((9, 65), rowspan=36, colspan=18),
                                                                                                  axes_class=fa.FloatingAxes, grid_helper=ghelper)
            elif changed_plot_type in ['statsummary', 'metadata']:
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((15, 50), rowspan=36, colspan=50))
            elif changed_plot_type == 'fairmode-statsummary':
                inner_gs = canvas_instance.gridspec.new_subplotspec((14, 63), rowspan=36, colspan=36).subgridspec(nrows, ncols,**canvas_instance.plot_characteristics["fairmode-statsummary"]["gridspec_kw"])
            elif changed_plot_type != 'None':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((15, 53), rowspan=34, colspan=51))
            
        # position 3 (bottom left)
        if changed_position == self.cb_position_3 or changed_position == 3:
            if changed_plot_type in ['periodic', 'periodic-violin']:
                canvas_instance.plot_axes[changed_plot_type] = {}
                canvas_instance.plot_axes[changed_plot_type]['hour'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 4), rowspan=18, colspan=28))
                canvas_instance.plot_axes[changed_plot_type]['dayofweek'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 22), rowspan=18, colspan=10))
                canvas_instance.plot_axes[changed_plot_type]['month'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 4), rowspan=18, colspan=17))
            elif changed_plot_type == 'taylor':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 5), rowspan=38, colspan=18),
                                                                                                  axes_class=fa.FloatingAxes, grid_helper=ghelper)
            elif changed_plot_type in ['statsummary', 'metadata']:
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 1), rowspan=38, colspan=28))
            elif changed_plot_type == 'fairmode-statsummary':
                inner_gs = canvas_instance.gridspec.new_subplotspec((61, 8), rowspan=38, colspan=24).subgridspec(nrows, ncols,**canvas_instance.plot_characteristics["fairmode-statsummary"]["gridspec_kw"])
            elif changed_plot_type != 'None':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 4), rowspan=38, colspan=28))
        # position 4 (bottom centre)
        if changed_position == self.cb_position_4 or changed_position == 4:
            if changed_plot_type in ['periodic', 'periodic-violin']:
                canvas_instance.plot_axes[changed_plot_type] = {}
                canvas_instance.plot_axes[changed_plot_type]['hour'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 38), rowspan=18, colspan=28))
                canvas_instance.plot_axes[changed_plot_type]['dayofweek'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 56), rowspan=18, colspan=10))
                canvas_instance.plot_axes[changed_plot_type]['month'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 38), rowspan=18, colspan=17))
            elif changed_plot_type == 'taylor':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 40), rowspan=38, colspan=18),
                                                                                                  axes_class=fa.FloatingAxes, grid_helper=ghelper)
            elif changed_plot_type in ['statsummary', 'metadata']:
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 35), rowspan=38, colspan=28))
            elif changed_plot_type == 'fairmode-statsummary':
                inner_gs = canvas_instance.gridspec.new_subplotspec((61, 41), rowspan=38, colspan=25).subgridspec(nrows, ncols,**canvas_instance.plot_characteristics["fairmode-statsummary"]["gridspec_kw"])
            elif changed_plot_type != 'None':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 38), rowspan=38, colspan=28))
            
        # position 5 (bottom right)
        if changed_position == self.cb_position_5 or changed_position == 5:
            if changed_plot_type in ['periodic', 'periodic-violin']:
                canvas_instance.plot_axes[changed_plot_type] = {}
                canvas_instance.plot_axes[changed_plot_type]['hour'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 72), rowspan=18, colspan=28))
                canvas_instance.plot_axes[changed_plot_type]['dayofweek'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 90), rowspan=18, colspan=10))
                canvas_instance.plot_axes[changed_plot_type]['month'] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((82, 72), rowspan=18, colspan=17))
            elif changed_plot_type == 'taylor':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 72), rowspan=38, colspan=18),
                                                                                                  axes_class=fa.FloatingAxes, grid_helper=ghelper)
            elif changed_plot_type in ['statsummary', 'metadata']:
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 69), rowspan=38, colspan=28))
            elif changed_plot_type == 'fairmode-statsummary':
                inner_gs = canvas_instance.gridspec.new_subplotspec((61, 75), rowspan=38, colspan=25).subgridspec(nrows, ncols,**canvas_instance.plot_characteristics["fairmode-statsummary"]["gridspec_kw"])
            elif changed_plot_type != 'None':
                canvas_instance.plot_axes[changed_plot_type] = canvas_instance.figure.add_subplot(canvas_instance.gridspec.new_subplotspec((60, 72), rowspan=38, colspan=28))

        # initialise polar axis for Taylor plots
        if changed_plot_type == 'taylor':
            canvas_instance.plot.taylor_polar_relevant_axis = canvas_instance.plot_axes[changed_plot_type].get_aux_axes(PolarAxes.PolarTransform())
        
        elif changed_plot_type == "fairmode-statsummary":
            # create gridspec and add it to a list
            canvas_instance.plot_axes[changed_plot_type] = [canvas_instance.figure.add_subplot(inner_gs[i, j]) for i in range(nrows) for j in range(ncols)]

    @QtCore.pyqtSlot()
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
        self.mpl_canvas.relative_selected_station_inds = np.array([], dtype=np.int64)
        self.mpl_canvas.absolute_selected_station_inds = np.array([], dtype=np.int64)
        self.mpl_canvas.absolute_non_selected_station_inds = np.array([], dtype=np.int64)

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
        self.previous_data_labels_raw = self.data_labels_raw
        self.mpl_canvas.previous_plot_options = copy.deepcopy(self.mpl_canvas.current_plot_options) 

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
        self.data_labels = [self.observations_data_label] + list(self.experiments.values())
        self.data_labels_raw = [self.observations_data_label] + list(self.experiments.keys())
        # remove bias plot options if have no experiments loaded
        if len(self.data_labels) == 1:
            for plot_type in self.mpl_canvas.all_plots:
                if 'bias' in self.mpl_canvas.current_plot_options[plot_type]:
                    self.mpl_canvas.current_plot_options[plot_type].remove('bias')
                    self.mpl_canvas.update_plot_options(plot_type)
                    if plot_type == 'statsummary':
                        self.block_config_bar_handling_updates = True
                        self.mpl_canvas.check_statsummary_stats()
                        self.block_config_bar_handling_updates = False

        # if spatial_colocation is not active, force filter_species to be empty dict if it is not already
        # inform user of this
        if (self.filter_species) and (not self.spatial_colocation):
            self.filter_species = {} 
            msg = '"spatial_colocation" must be set to True if wanting to use "filter_species" option.'
            show_message(self.read_instance, msg)

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
                str(dict(sorted(self.filter_species.items()))) != str(dict(sorted(self.previous_filter_species.items())))) or (
                str(dict(sorted(self.calibration_factor.items()))) != str(dict(sorted(self.previous_calibration_factor.items())))):
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
        experiments_to_remove = [experiment for experiment in self.previous_experiments.values() if experiment not in list(self.experiments.values())]
        experiments_to_read = [experiment for experiment in self.experiments.values() if experiment not in list(self.previous_experiments.values())]
        if 'reset' not in read_operations:
            if len(experiments_to_remove) > 0:
                read_operations.append('remove_exp')
            if len(experiments_to_read) > 0:
                read_operations.append('read_exp')

        # has date range changed?
        if len(read_operations) > 0:
            
            # if reading/cutting observations then cover canvas to do updates gracefully
            if ('reset' in read_operations) or ('cut_left' in read_operations) or ('cut_right' in read_operations) or\
               ('read_left' in read_operations) or ('read_right' in read_operations):
                self.mpl_canvas.canvas_cover.show()
            # otherwise, just cover plotting axes as are adding/removing experiments
            else:
                self.mpl_canvas.top_right_canvas_cover.show() 
                self.mpl_canvas.lower_canvas_cover.show()
            self.mpl_canvas.figure.canvas.draw_idle()  
            self.mpl_canvas.figure.canvas.flush_events()

            # clear all axes elements 
            for plot_type, ax in self.mpl_canvas.plot_axes.items():
                self.mpl_canvas.remove_axis_elements(ax, plot_type)

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
            update_representativity_fields(self)
            update_period_fields(self)

            # for non-GHOST delete valid station indices variables because we do not want to 
            # remove the stations with 0 valid measurements before the filter has been updated, 
            # this will happen later
            if hasattr(self, 'valid_station_inds') and (not self.reading_ghost):
                delattr(self, 'valid_station_inds')
                delattr(self, 'valid_station_inds_temporal_colocation')

            update_metadata_fields(self)
            
            # update relevant/nonrelevant temporal resolutions 
            self.relevant_temporal_resolutions = get_relevant_temporal_resolutions(self.resolution)
            self.nonrelevant_temporal_resolutions = get_nonrelevant_temporal_resolutions(self.resolution)

            # if species has changed, or first read, update species specific lower/upper limits
            if (self.first_read) or (self.species[0] != self.previous_species[0]):
                # get default GHOST limits
                self.lower_bound[self.species[0]] = np.float32(self.parameter_dictionary[self.species[0]]['extreme_lower_limit']) 
                self.upper_bound[self.species[0]] = np.float32(self.parameter_dictionary[self.species[0]]['extreme_upper_limit']) 
                self.le_minimum_value.setText(str(self.lower_bound[self.species[0]]))
                self.le_maximum_value.setText(str(self.upper_bound[self.species[0]]))

            # run function to update filter
            self.mpl_canvas.handle_data_filter_update()

            # for non-GHOST, we call update_metadata_fields again to remove the stations that have
            # 0 valid measurements, to do this we need to have valid_station_inds, which is obtained 
            # after filtering
            if not self.reading_ghost:
                update_metadata_fields(self)

            # generate list of sorted z1/z2 data arrays names in memory, putting observations
            # before experiments, and empty string item as first element in z2 array list
            # (for changing from 'difference' statistics to 'absolute')
            if len(self.data_labels) == 1:  
                self.z1_arrays = np.array([self.observations_data_label])
            else:
                self.z1_arrays = np.array(self.data_labels)
            self.z2_arrays = np.append([''], self.z1_arrays)

            # update map z statistic comboboxes
            self.mpl_canvas.handle_map_z_statistic_update()

            # update timeseries statistic combobox
            self.mpl_canvas.handle_timeseries_statistic_update()

            # update resampling resolution
            self.mpl_canvas.handle_resampling_update()
            
            # update timeseries chunk statistic and resolution comboboxes
            self.mpl_canvas.handle_timeseries_chunk_statistic_update()

            # update periodic statistic combobox
            self.mpl_canvas.handle_periodic_statistic_update()

            # update taylor diagram statistic combobox
            self.mpl_canvas.handle_taylor_correlation_statistic_update()

            # update statsummary statistic comboboxes
            self.mpl_canvas.handle_statsummary_statistics_update()
            self.mpl_canvas.handle_statsummary_cycle_update()
            self.mpl_canvas.handle_statsummary_periodic_aggregation_update()
            self.mpl_canvas.handle_statsummary_periodic_mode_update()

            # update fairmode target combobox
            self.mpl_canvas.handle_fairmode_target_classification_update()

            # unselect all/intersect/extent checkboxes
            self.mpl_canvas.unselect_map_checkboxes()

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

    @QtCore.pyqtSlot()
    def reset_options(self):
        """ Reset all filter fields to initial values. """

        if self.block_MPL_canvas_updates:
            return

        # set mouse cursor to hourglass
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        # reset representativity fields        
        init_representativity(self)
        update_representativity_fields(self)

        # reset period fields 
        init_period(self)
        update_period_fields(self)

        # reset metadata
        init_metadata(self)
        update_metadata_fields(self)

        # reset bounds
        species_lower_limit = np.float32(self.parameter_dictionary[self.species[0]]['extreme_lower_limit'])
        species_upper_limit = np.float32(self.parameter_dictionary[self.species[0]]['extreme_upper_limit'])
        
        # set default limits
        self.le_minimum_value.setText(str(species_lower_limit))
        self.le_maximum_value.setText(str(species_upper_limit))
        
        # unfilter data
        self.mpl_canvas.handle_data_filter_update()
        
        # for non-GHOST, we call update_metadata_fields again to remove the stations that have
        # 0 valid measurements, to do this we need to have valid_station_inds, which is obtained 
        # after filtering
        if not self.reading_ghost:
            update_metadata_fields(self)

        # Restore mouse cursor to normal
        QtWidgets.QApplication.restoreOverrideCursor()

    def disable_ghost_buttons(self):
        """ Disable button related only to ghost data. """
        
        # change background-color to indicate that it's nonusable
        self.bu_flags = set_formatting(self.bu_flags, self.formatting_dict['menu_button_disabled'], disabled=True)
        self.bu_QA = set_formatting(self.bu_QA, self.formatting_dict['menu_button_disabled'], disabled=True)
        self.bu_period = set_formatting(self.bu_period, self.formatting_dict['menu_button_disabled'], disabled=True)
        
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
    
    if sys.platform.startswith('darwin'):
        # Set app name, if PyObjC is installed
        # Python 2 has PyObjC preinstalled
        # Python 3: pip3 install pyobjc-framework-Cocoa
        from Foundation import NSBundle
        bundle = NSBundle.mainBundle()
        if bundle:
            app_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
            app_info = bundle.localizedInfoDictionary() or bundle.infoDictionary()
            if app_info:       
                app_info['CFBundleName'] = app_name

    # pause briefly to allow QT modules time to correctly initilise
    time.sleep(0.1)

    # create application
    q_app = QtWidgets.QApplication(sys.argv)

    # set Fusion style
    q_app.setStyle("Fusion")

    # explicitely set colour palette to avoid issues with dark modes (e.g. on Mac)
    p = q_app.palette()
    dcp = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/internal/dashboard_colour_palette.yaml')))
    p.setColor(QtGui.QPalette.Dark, QtGui.QColor(*dcp['Dark']))
    p.setColor(QtGui.QPalette.Light, QtGui.QColor(*dcp['Light']))
    p.setColor(QtGui.QPalette.Window, QtGui.QColor(*dcp['Window']))
    p.setColor(QtGui.QPalette.WindowText, QtGui.QColor(*dcp['WindowText']))
    p.setColor(QtGui.QPalette.Base, QtGui.QColor(*dcp['Base']))
    p.setColor(QtGui.QPalette.AlternateBase, QtGui.QColor(*dcp['AlternateBase']))
    p.setColor(QtGui.QPalette.ToolTipBase, QtGui.QColor(*dcp['ToolTipBase']))
    p.setColor(QtGui.QPalette.ToolTipText, QtGui.QColor(*dcp['ToolTipText']))
    p.setColor(QtGui.QPalette.Button, QtGui.QColor(*dcp['Button']))
    p.setColor(QtGui.QPalette.ButtonText, QtGui.QColor(*dcp['ButtonText']))
    p.setColor(QtGui.QPalette.BrightText, QtGui.QColor(*dcp['BrightText']))
    p.setColor(QtGui.QPalette.Highlight, QtGui.QColor(*dcp['Highlight']))
    p.setColor(QtGui.QPalette.HighlightedText, QtGui.QColor(*dcp['HighlightedText']))
    p.setColor(QtGui.QPalette.Link, QtGui.QColor(*dcp['Link']))
    p.setColor(QtGui.QPalette.LinkVisited, QtGui.QColor(*dcp['LinkVisited']))
    p.setColor(QtGui.QPalette.Mid, QtGui.QColor(*dcp['Mid']))
    p.setColor(QtGui.QPalette.Midlight, QtGui.QColor(*dcp['Midlight']))
    p.setColor(QtGui.QPalette.PlaceholderText, QtGui.QColor(*dcp['PlaceholderText']))
    p.setColor(QtGui.QPalette.Shadow, QtGui.QColor(*dcp['Shadow']))
    p.setColor(QtGui.QPalette.Text, QtGui.QColor(*dcp['Text']))
    q_app.setPalette(p)
    
    q_app.setWindowIcon(QtGui.QIcon(join(PROVIDENTIA_ROOT, 'assets/logo.icns')))
    q_app.setApplicationName("Providentia")
    q_app.setApplicationDisplayName("Providentia")
    q_app.setDesktopFileName("Providentia")

    # open Providentia
    ProvidentiaMainWindow(**kwargs)
    sys.exit(q_app.exec_())