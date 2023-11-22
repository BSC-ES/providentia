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
from .read import DataReader
from .read_aux import (get_ghost_observational_tree, get_nonghost_observational_tree, 
                       get_nonrelevant_temporal_resolutions, get_relevant_temporal_resolutions, 
                       get_valid_experiments, get_valid_obs_files_in_date_range)
from .statistics import calculate_statistic, generate_colourbar, get_selected_station_data, get_z_statistic_info

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, '../settings/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, '../settings/experiment_bias_stats.json')))


class ProvidentiaInteractive:
    """ Class to create interactive Providentia session."""

    def __init__(self, **kwargs):

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

        # check for self defined plot characteristics file
        if self.plot_characteristics_filename == '':
            self.plot_characteristics_filename = os.path.join(CURRENT_PATH, '../settings/plot_characteristics_offline.json')
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
        print("Reading data\n")
        self.datareader.read_setup(['reset'])

        if self.invalid_read:
            print('No valid data to read')
            return

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

        # set previous QA, flags and filter species as subsection
        #self.previous_qa = copy.deepcopy(self.qa)
        #self.previous_flags = copy.deepcopy(self.flags)
        #self.previous_filter_species = copy.deepcopy(self.filter_species)

        # filter dataset for current subsection
        print('Filtering data\n')
        DataFilter(self)

def read(**kwargs):

    kwargs['offline'] = True

    # initialise interactive class 
    ProvInt = ProvidentiaInteractive(**kwargs)

    return ProvInt

