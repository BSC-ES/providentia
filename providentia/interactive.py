""" Class to create interactive session"""

import copy
import datetime
import json
import os
import sys

import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib.projections import PolarAxes
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as fa
import numpy as np
import pandas as pd

from .configuration import load_conf
from .configuration import ProvConfiguration
from .fields_menus import (init_metadata, init_period, init_representativity, metadata_conf,
                           update_metadata_fields, update_period_fields, update_representativity_fields,
                           period_conf, representativity_conf)
from .filter import DataFilter
from .plot import Plot
from .read import DataReader
from .read_aux import (get_ghost_observational_tree, get_nonghost_observational_tree, 
                       get_nonrelevant_temporal_resolutions, get_relevant_temporal_resolutions, 
                       get_valid_experiments, get_valid_obs_files_in_date_range)
from .statistics import calculate_statistic, generate_colourbar, get_selected_station_data, get_z_statistic_info

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, '../settings/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, '../settings/experiment_bias_stats.json')))


class Interactive:
    """ Class to create interactive Providentia session."""

    def __init__(self, **kwargs):

        #set mode as offline
        kwargs['offline'] = True

        #set configuration variables, as well as any other defined variables
        self.set_config(**kwargs)

        # create dictionary of all available observational GHOST data
        self.all_observation_data = get_ghost_observational_tree(self)

        # load dictionary with non-GHOST esarchive files to read
        nonghost_observation_data_json = json.load(open(os.path.join(CURRENT_PATH, '../settings/nonghost_files.json')))
        # merge to existing GHOST observational data dict if we have the path
        if self.nonghost_root is not None:
            nonghost_observation_data = get_nonghost_observational_tree(self, nonghost_observation_data_json)
            self.all_observation_data = {**self.all_observation_data, **nonghost_observation_data}

        # initialise DataReader class
        self.datareader = DataReader(self)

        # initialise Plot class
        self.plot = Plot(read_instance=self, canvas_instance=self)

        # check for self defined plot characteristics file
        if self.plot_characteristics_filename == '':
            self.plot_characteristics_filename = os.path.join(CURRENT_PATH, '../settings/plot_characteristics_interactive.json')
        self.plot_characteristics_templates = json.load(open(self.plot_characteristics_filename))

        # set some key configuration variables
        self.relevant_temporal_resolutions = get_relevant_temporal_resolutions(self.resolution)
        self.nonrelevant_temporal_resolutions = get_nonrelevant_temporal_resolutions(self.resolution)
        self.data_labels = ['observations'] + list(self.experiments.keys())
        self.networkspecies = ['{}|{}'.format(network,speci) for network, speci in zip(self.network, self.species)]

        # get valid observations in date range
        get_valid_obs_files_in_date_range(self, self.start_date, self.end_date)

        # update available experiments for selected fields
        get_valid_experiments(self, self.start_date, self.end_date, self.resolution,
                                self.network, self.species)

        # if have no experiments, force temporal colocation to be False
        if len(self.experiments) == 0:
            self.temporal_colocation = False  

        # read data
        self.read()  

        # update fields available for filtering
        init_representativity(self)
        update_representativity_fields(self)
        representativity_conf(self)
        init_period(self)
        update_period_fields(self)
        period_conf(self)
        init_metadata(self)
        update_metadata_fields(self)
        metadata_conf(self)

        # filter data
        self.filter()

        def read(self):
            """wrapper method to read data"""

            self.datareader.read_setup(['reset'])

            if self.invalid_read:
                print('No valid data to read')
                return

        def filter(self):
            """wrapper method to filter data"""

            DataFilter(self)

        def reset_filter(self):
            """wrapper method to reset filter data"""

        def plot(self, plot=''):
            """wrapper method to make a Providentia plot"""

            self.plot.set_plot_characteristics(plot)

            func = getattr(self.plot, 'make_{}'.format(plot))

            if base_plot_type == 'metadata':
                func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                        plot_options=plot_options, station_inds=station_inds)
            elif base_plot_type == 'periodic':
                func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                        zstat=zstat, plot_options=plot_options)    
            elif base_plot_type == 'distribution':
                func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                    plot_options=plot_options, data_range_min=data_range_min, data_range_max=data_range_max) 
            elif base_plot_type == 'taylor':
                func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                    plot_options=plot_options, stddev_max=stddev_max)
            else:
                func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                    plot_options=plot_options) 

            self.set_plot_title()

        def set_plot_title(self):
            """wrapper method to set plot title"""

        def calculate_statistic(self, statistic=''):
            """wrapper method to calculate statistic/s"""

            return statistic

        def set_config(self, config='', section='', subsection='', **kwargs):
            """wrapper method to set configuration variables"""

            # initialise default configuration variables
            # modified by passed arguments, if given
            provconf = ProvConfiguration(self, **kwargs)

            # update variables from config file
            if ('config' in kwargs) and (os.path.exists(kwargs['config'])):
                load_conf(self, kwargs['config'])
                self.from_conf = True
            elif ('config' in kwargs) and (not os.path.exists(kwargs['config'])):     
                error = 'Error: The path to the configuration file specified in the command line does not exist.'
                sys.exit(error)
            else:
                error = "Error: No configuration file found. The path to the config file must be added as an argument."
                sys.exit(error)

            # parse first available configuration section
            if len(self.parent_section_names) > 1:
                print('Warning: In interactive mode, only the first defined section will be read.')
            self.section = self.parent_section_names[0]

            # update self with section variables
            self.section_opts = self.sub_opts[self.section]
            for k, val in self.section_opts.items():
                setattr(self, k, provconf.parse_parameter(k, val))

            # get subsection names
            self.child_subsection_names = [subsection_name for subsection_name in self.subsection_names 
                                            if self.section == subsection_name.split('·')[0]]
            if len(self.child_subsection_names) > 0:
                if len(self.child_subsection_names) > 1:
                    print('Warning: In interactive mode, only the first defined subsection will be read.')

                    self.subsection = self.child_subsection_names[0]

                # get subsection variables
                self.subsection_opts = self.sub_opts[self.subsection]

                # ensure all fixed section variables defined in subsection have same value as current section variables
                self.subsection_opts = {k: (self.section_opts[k] if k in self.fixed_section_vars else val) 
                                        for (k, val) in self.subsection_opts.items()}

                # update subsection variables
                for k, val in self.subsection_opts.items():
                    setattr(self, k, provconf.parse_parameter(k, val))
            else:
                self.subsection = [self.section]

            # now all variables have been parsed, check validity of those, throwing errors where necessary
            provconf.check_validity()

        def select_station(self, stations=''):
            """wrapper method to select specific station/s"""

        def save(self, fname='', format='nc'):
            """wrapper method to save current data/ metadata in memory"""


        def get_data(self, format='nc'):
            """wrapper method return data / metadata in specific format"""

            return data

        def get_var(self, var=''):
            """wrapper method to return specific data / metadata variable"""

            return var_data





