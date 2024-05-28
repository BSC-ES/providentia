""" Providentia Configuration Module """

import configparser
import copy
import json
import os
import platform
import re
import socket
import subprocess
import sys
import yaml

import matplotlib
import numpy as np
from packaging.version import Version
import pandas as pd

from .read_aux import check_for_ghost, get_default_qa
from .warnings import show_message

MACHINE = os.environ.get('BSC_MACHINE', '')
CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
data_paths = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/data_paths.yaml')))

# set MACHINE to be the hub, workstation or local machine
if MACHINE not in ['power', 'mn4', 'nord3v2', 'mn5']:
    hostname = os.environ.get('HOSTNAME', '')
    ip = socket.gethostbyname(socket.gethostname())
    if "bscearth" in hostname:
        MACHINE = "workstation"
    elif "bscesdust02.bsc.es" in hostname:
        MACHINE = "dust"
    elif ip == "84.88.185.48":
        MACHINE = "hub"
    else:
        MACHINE = "local"


def parse_path(dir, f):
    if os.path.isabs(f):
        return f
    else:
        return os.path.join(dir, f)

class ProvConfiguration:
    """ Class that handles the configuration parameters definitions. """

    def __init__(self, read_instance, **kwargs):
        
        self.read_instance = read_instance 

        # set variable defaults
        self.var_defaults = {
            'ghost_version': '1.5',
            'conf': '',
            'config': '',
            'config_dir': os.path.join(PROVIDENTIA_ROOT, 'configurations/'),
            'operating_system': '',
            'machine': '',
            'cartopy_data_dir': '',
            'available_cpus': '',
            'n_cpus': '',
            'ghost_root': '',
            'nonghost_root': '',
            'exp_root': '',
            'generate_file_tree': False,
            'offline': False,
            'interactive': False,
            'download':False,
            'available_resolutions': ['hourly', '3hourly', '6hourly', 'hourly_instantaneous',
                                      '3hourly_instantaneous', '6hourly_instantaneous',
                                      'daily', 'monthly'],
            'available_networks': ['GHOST','AERONET_v3_lev1.5','AERONET_v3_lev2.0','BJMEMC','CANADA_NAPS','CAPMoN',
                                   'CHILE_SINCA','CNEMC','EANET','EBAS','EBAS-ACTRIS','EBAS-AMAP','EBAS-CAMP', 
                                   'EBAS_COLOSSAL','EBAS-EMEP','EBAS-EUCAARI','EBAS-EUSAAR','EBAS-GUAN','EBAS-HELCOM','EBAS-HTAP',
                                   'EBAS-IMPACTS','EBAS-IMPROVE','EBAS-Independent','EBAS-MOE','EBAS-NILU',
                                   'EBAS-NOAA_ESRL','EBAS-NOAA_GGGRN','EBAS-OECD','EBAS-UK_DECC','EBAS-WMO_WDCA', 'EBAS-WMO_WDCRG',
                                   'EEA_AIRBASE','EEA_AQ_eReporting','JAPAN_NIES','MEXICO_CDMX','MITECO',
                                   'NOAA_ISD','NOAA_ISD_EU','NOAA_ISD_IP','NOAA_ISD_NA','SEARCH','UK_AIR',
                                   'US_EPA_AirNow_DOS','US_EPA_AQS','US_EPA_CASTNET','US_NADP_AMNet','US_NADP_AMoN','WMO_WDCGG'], 
            'network': None,
            'species': None,
            'resolution': None,
            'start_date': None,
            'end_date': None,
            'download_source': 'ghost',
            'observations_data_label': 'observations',
            'experiments': {},
            'qa': None,
            'flags': None,
            'temporal_colocation': False,
            'spatial_colocation': True,
            'map_extent': None, 
            'filter_species': {},
            'calibration_factor': {},
            'lower_bound': None,
            'upper_bound': None,
            'report_type': 'standard',
            'report_summary': True,
            'report_stations': False,
            'report_title': 'Providentia Offline Report',
            'report_filename': 'PROVIDENTIA_Report',
            'active_dashboard_plots': None,
            'resampling_resolution': None,
            'statistic_mode': None,
            'statistic_aggregation': None,
            'periodic_statistic_mode': None,
            'periodic_statistic_aggregation': None,
            'timeseries_statistic_aggregation': None,
            'plot_characteristics_filename': '',
            'harmonise_summary': True,
            'harmonise_stations': True,
            'remove_extreme_stations': None,
            'fixed_section_vars':  ['ghost_version', 'config_dir', 'operating_system', 'cartopy_data_dir', 'available_cpus', 
                                    'n_cpus', 'ghost_root', 'nonghost_root', 'exp_root', 'offline', 'interactive',
                                    'available_resolutions', 'available_networks',
                                    'network', 'species', 'resolution', 'start_date', 'end_date', 
                                    'observations_data_label', 'experiments', 'temporal_colocation', 'spatial_colocation', 
                                    'report_type', 'report_summary', 'report_stations', 'report_title', 
                                    'report_filename', 'plot_characteristics_filename', 
                                    'harmonise_summary', 'harmonise_stations']
        }

        # if variable is given by command line, set that value, otherwise set as default value 
        for k, val in self.var_defaults.items():
            val = kwargs.get(k, val)
            setattr(self.read_instance, k, self.parse_parameter(k, val))

    def parse_parameter(self, key, value):
        """ Parse a parameter. """

        # parse config file name
        if key == 'conf':
            if value != '':
                self.read_instance.config = value
            return ''

        elif key == 'config':
            if hasattr(self.read_instance, 'config'):
                if self.read_instance.config != '':
                    return self.read_instance.config
                else:
                    return value
            else:
                return value
    
        elif key == 'operating_system':
            # get operating system
            operating_system = platform.system()
            if operating_system == 'Darwin':
                operating_system = 'Mac'
            elif operating_system == 'Linux':
                operating_system = 'Linux'
            elif operating_system in ['Windows','MINGW32_NT','MINGW64_NT']:
                operating_system = 'Windows'
            else:
                error = 'Error: The OS cannot be detected.'
                sys.exit(error)
            return operating_system
        
        elif key == 'machine':
            return MACHINE

        elif key == 'available_cpus':
            # get available N CPUs
            if MACHINE in ['power', 'mn4', 'nord3v2', 'mn5']:
                # handle cases where are testing briefly on login nodes (1 cpu capped)
                try:
                    return int(os.getenv('SLURM_CPUS_PER_TASK'))
                except:
                    return 1
            else:
                if self.read_instance.operating_system == 'Linux':
                    return len(os.sched_getaffinity(0))
                else:
                    return os.cpu_count()

        elif key == 'cartopy_data_dir':
            # set cartopy data directory (needed on CTE-POWER/MN4/Nord3v2 as has no external
            # internet connection)
            if MACHINE in ['power', 'mn4', 'nord3v2']:
                return '/gpfs/projects/bsc32/software/rhel/7.5/ppc64le/POWER9/software/Cartopy/0.17.0-foss-2018b-Python-3.7.0/lib/python3.7/site-packages/Cartopy-0.17.0-py3.7-linux-ppc64le.egg/cartopy/data'
            # set directory in MN5 to avoid network issues
            elif MACHINE == 'mn5':
                return '/gpfs/projects/bsc32/software/rhel/7.5/ppc64le/POWER9/software/Cartopy/0.17.0-foss-2018b-Python-3.7.0/lib/python3.7/site-packages/Cartopy-0.17.0-py3.7-linux-ppc64le.egg/cartopy/data'
            # on all other machines pull from internet

        elif key == 'n_cpus':
            # define number of CPUs to process on (leave empty to automatically
            # utilise all available CPUs) NOTE: if this value is set higher than the
            # actual number of CPUs available, then the max number of CPUs is used.

            if (value == '') or (int(value) > self.read_instance.available_cpus):
                return self.read_instance.available_cpus

        elif key == 'ghost_root':
            # define GHOST observational root data directory (if undefined it is
            # automatically taken from the BSC machine the tool is ran on)

            # set default if left undefined
            if value == '':
                ghost_root = data_paths[MACHINE]["ghost_root"]
                return os.path.expanduser(ghost_root[0])+ghost_root[1:]

        elif key == 'nonghost_root':
            # define non-GHOST observational root data directory (if undefined it is
            # automatically taken from the BSC machine the tool is ran on)

            # set default if left undefined
            if value == '':
                nonghost_root = data_paths[MACHINE]["nonghost_root"]
                return os.path.expanduser(nonghost_root[0])+nonghost_root[1:]

        elif key == 'exp_root':
            # define experiment root data directory
            # set experiment root data directory if left undefined
            if value == '':
                exp_root = data_paths[MACHINE]["exp_root"]
                return os.path.expanduser(exp_root[0])+exp_root[1:]
            
        elif key == 'ghost_version':
            # parse GHOST version

            # import GHOST standards 
            sys.path = [path for path in sys.path if 'dependencies/GHOST_standards/' not in path]            
            sys.path.insert(1, os.path.join(CURRENT_PATH, 'dependencies/GHOST_standards/{}'.format(value)))
            if 'GHOST_standards' in sys.modules:
                del sys.modules['GHOST_standards']
            from GHOST_standards import standard_parameters
            from GHOST_standards import get_standard_metadata
            from GHOST_standards import standard_data_flag_name_to_data_flag_code
            from GHOST_standards import standard_QA_name_to_QA_code

            # modify standard parameter dictionary to have BSC standard parameter names as keys (rather than GHOST)
            self.read_instance.parameter_dictionary = dict()
            for _, param_dict in standard_parameters.items():
                self.read_instance.parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict
            
            # get standard metadata dictionary
            self.read_instance.standard_metadata = get_standard_metadata({'standard_units': ''})
            
            # create list of GHOST metadata variables to read
            self.read_instance.ghost_metadata_vars_to_read = [key for key in self.read_instance.standard_metadata.keys() if
                                                              pd.isnull(self.read_instance.standard_metadata[key]['metadata_type']) == False]
            self.read_instance.ghost_metadata_dtype = [(key, self.read_instance.standard_metadata[key]['data_type']) 
                                                       for key in self.read_instance.ghost_metadata_vars_to_read]
            self.read_instance.standard_data_flag_name_to_data_flag_code = standard_data_flag_name_to_data_flag_code
            self.read_instance.standard_QA_name_to_QA_code = standard_QA_name_to_QA_code

            return str(value)

        elif key == 'network':
            # parse network

            if isinstance(value, str):
                # parse multiple networks
                if ',' in value:
                    return [network.strip() for network in value.split(',')]
                else:
                    return [value.strip()]

        elif key == 'species':
            # parse species

            if isinstance(value, str):
                # parse multiple species
                if ',' in value:
                    return [speci.strip() for speci in value.split(',')]
                else:
                    return [value.strip()]

        elif key == 'resolution':
            # parse resolution
            
            if isinstance(value, str):
                # resolution is a list in download mode
                if self.read_instance.download:
                    # parse multiple resolutions
                    if ',' in value:
                        return [speci.strip() for speci in value.split(',')]
                    else:
                        return [value.strip()]
                else:
                    # parse one single resolution
                    return value.strip()

        elif key == 'start_date':
            # parse start_date

            if (isinstance(value, str)) or (isinstance(value, int)):
                # throw error if start_date is empty str
                value = str(value)
                return value.strip()

        elif key == 'end_date':
            # parse end_date

            if (isinstance(value, str)) or (isinstance(value, int)):
                # throw error if start_date is empty str
                value = str(value)
                return value.strip()
        
        elif key == 'qa':
            # parse qa

            # set default qa codes (can differ per GHOST version)
            from GHOST_standards import providentia_defaults
            self.read_instance.default_qa_standard = [self.read_instance.standard_QA_name_to_QA_code[qa_name] 
                                                      for qa_name in providentia_defaults['qa_standard']]
            self.read_instance.default_qa_non_negative = [self.read_instance.standard_QA_name_to_QA_code[qa_name] 
                                                          for qa_name in providentia_defaults['qa_non_negative']]

            # if not None then set QA by that given
            if value is not None:
                # if conf has only 1 QA
                if isinstance(value, int):
                    return [value]
                # empty string
                elif value == "":
                    return []
                # if the QAs are written with their names
                elif isinstance(value, str):
                    return sorted([self.read_instance.standard_QA_name_to_QA_code[q.strip()] for q in value.split(",")])
                # list of integer codes
                else:
                    return sorted(list(value))
            # otherwise, set default QA per species (set later)
            else:
                # set qa to be empty dict (to be later filled)
                return {}

        elif key == 'flags':
            # parse flags

            # if not None then set flags by that given
            if value is not None:
                # if conf has only one flag
                if isinstance(value, int):
                    return [value]
                # empty string
                elif value == "":
                    return []
                # if the flags are written with their names
                elif isinstance(value, str):
                    return [self.read_instance.standard_data_flag_name_to_data_flag_code[f.strip()]
                            for f in value.split(",")]
                # list of integer codes
                else:
                    return sorted(list(value))
            # otherwise, set default (empty list)
            else:
                return []

        elif key == 'experiments':
            # parse experiments

            if isinstance(value, str):
                # empty string
                if value == "":
                    return {}
                # split experiments
                else:
                    # have alternative experiment names for the legend, then parse them?
                    if ('(' in value) & (')' in value):
                        exps = [exp.strip() for exp in value.split('(')[0].strip().split(",")]
                        exps_legend = [exp_legend.strip() for exp_legend in value.split('(')[1].split(')')[0].strip().split(",")]
                    # otherwise set legend names as given experiment names in full
                    else: 
                        exps = [exp.strip() for exp in value.split(",")]
                        exps_legend = copy.deepcopy(exps)
                    return {exp:exp_legend for exp,exp_legend in zip(exps,exps_legend)}

        elif key == 'map_extent':
            # parse map extent

            # map_extent empty?
            if not value:
                # if dashboard, if map extent not defined then fix to global default extent
                if (not self.read_instance.offline) and (not self.read_instance.interactive):
                    return [-180, 180, -90, 90]
            # otherwise parse it
            else:
                if isinstance(value, str):
                    return [float(c.strip()) for c in value.split(',')]

        elif key == 'filter_species':
            # parse filter species

            # per networkspecies to filter by, save in dict as networkspecies:[lower_limit, upper_limit] 
            if isinstance(value, str):

                # strip all whitespace
                value_strip = "".join(value.split())

                # return empty dict if empty str
                if value_strip == '':
                    return {}
                else:
                    # split per networkspecies
                    networkspecies_split = value_strip.split('),')

                    # iterate through networkspecies, saving list of limits per networkspecies
                    filter_networkspecies_dict = {}
                    for networkspeci_split in networkspecies_split:
                        
                        networkspeci_split_2 = networkspeci_split.split('(')

                        # get networkspeci
                        networkspeci = networkspeci_split_2[0].replace(':','|')

                        # get lower and upper limits
                        networkspeci_split_3 = networkspeci_split_2[1].split(',')
                        lower_limit = networkspeci_split_3[0]

                        # if it has fill value
                        if len(networkspeci_split_3) > 2:
                            upper_limit = networkspeci_split_3[1]
                            filter_species_fill_value = float(networkspeci_split_3[2].replace(')',''))
                        # only bounds, fill value will be nan
                        else:
                            upper_limit = networkspeci_split_3[1].replace(')','')
                            filter_species_fill_value = np.nan

                        # save limits per networkspecies
                        if networkspeci in filter_networkspecies_dict:
                             filter_networkspecies_dict[networkspeci].append([lower_limit, upper_limit, 
                                                                              filter_species_fill_value])
                        else:
                            filter_networkspecies_dict[networkspeci] = [[lower_limit, upper_limit, 
                                                                         filter_species_fill_value]]
                    
                    return filter_networkspecies_dict

        elif key == 'lower_bound':
            # parse lower_bound
            
            # if not None then set lower_bound by that given
            # make sure it is a list of values
            if value is not None:
                if value == "":
                    return []
                elif isinstance(value, str):
                    return [np.float32(c.strip()) for c in value.split(',')]
                elif (isinstance(value, int)) or (isinstance(value, float)):
                    return [value]
                #otherwise must be a already a list of values
                else:
                    return value
            # lower_bound empty?
            # then set lower bound using GHOST extreme lower limit for all species in memory (set later)
            else:
                return {}

        elif key == 'upper_bound':
            # parse upper bound

            # if not None then set upper_bound by that given
            # make sure it is a list of values
            if value is not None:
                if value == "":
                    return []
                elif isinstance(value, str):
                    return [np.float32(c.strip()) for c in value.split(',')]
                elif (isinstance(value, int)) or (isinstance(value, float)):
                    return [value]
                # otherwise must be a already a list of values
                else:
                    return value
            # upper_bound empty?
            # then set upper bound using GHOST extreme upper limit for all species in memory (set later)
            else:
                return {}

        elif key == 'active_dashboard_plots':
            # parse active_dashboard_plots

            if isinstance(value, str):
                # parse multiple active_dashboard_plots
                if ',' in value:
                    return [plot.strip() for plot in value.split(',')]
                else:
                    return [value.strip()]

        elif key == 'resampling_resolution':
            # parse resampling resolution

            if isinstance(value, str):
                return value.strip()

        elif key == 'plot_characteristics_filename':
            # parse plot characteristics filename
    
            if isinstance(value, str):
                if value != "":
                    # various paths were provided
                    if "," in value:
                        if ("dashboard:" in value) and ((not self.read_instance.offline) and (not self.read_instance.interactive)):
                            return value.split("dashboard:")[1].split(',')[0]
                        elif ("offline:" in value) and (self.read_instance.offline):
                            return value.split("offline:")[1].split(',')[0]
                        elif ("interactive:" in value) and (self.read_instance.interactive):
                            return value.split("interactive:")[1].split(',')[0]
                        else:
                            msg = 'It is necessary to include the words dashboard, offline or interactive to set different plot characteristics filenames, as in: '
                            msg += 'plot_characteristics_filename = dashboard:/path/plot_characteristics_dashboard.yaml, offline:/path/plot_characteristics_offline.yaml.'
                            sys.exit(msg)
                    # one path was provided
                    else:
                        return value

        elif key == 'calibration_factor':
            # parse calibration factor
            
            if isinstance(value, str):

                # strip all whitespace
                value_strip = "".join(value.split())

                calibration_by_experiment = False
                for experiment in self.read_instance.experiments.keys():
                    if experiment in value_strip:
                        calibration_by_experiment = True
                        break
                
                if calibration_by_experiment:
                    calibration_factor_dict = {}
                    for i, experiment in enumerate(self.read_instance.experiments.keys()):
                        calibration_factor_exp = value_strip.split('(')[i+1].split(')')[0]
                        calibration_factor_dict[experiment] = calibration_factor_exp
                    return calibration_factor_dict
                else:
                    if np.issubdtype(type(value), np.number):
                        return str(value)

        # if no special parsing treatment for variable, simply return value
        return value

    def check_validity(self):
        """ Check validity of set variables after parsing. """
        
        # get non-default fields on config file if launching from a config file
        if hasattr(self.read_instance, "sub_opts"):
            self.read_instance.fields_per_section = {}
            for field_name, fields in self.read_instance.sub_opts.items():
                if field_name in self.read_instance.subsection_names:
                    section_field_name = field_name.split('·')[0]
                    self.read_instance.fields_per_section[field_name] = \
                        fields.keys() - set(self.read_instance.fields_per_section[section_field_name])
                else:
                    self.read_instance.fields_per_section[field_name] = set(fields.keys())
            self.read_instance.non_default_fields_per_section = {
                field_name:fields-set(self.var_defaults) 
                for field_name, fields in self.read_instance.fields_per_section.items()}
       
        # check have network information, 
        # if offline, throw message, stating are using default instead
        if not self.read_instance.network:
            #default = ['GHOST']
            default = ['EBAS']
            msg = "Network (network) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.network = default

        # check have species information, 
        # if offline, throw message, stating are using default instead
        if not self.read_instance.species:
            default = ['sconco3']
            msg = "Species (species) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.species = default

        # if number of networks and species is not the same,
        # and len of one of network or species == 1,
        # then duplicate respestive network/species
        if len(self.read_instance.network) != len(self.read_instance.species):

            # 1 network?
            if len(self.read_instance.network) == 1:
                # duplicate network to match species len
                self.read_instance.network = self.read_instance.network * len(self.read_instance.species)

            # 1 species?
            elif len(self.read_instance.species) == 1:
                # duplicate species to match network len
                self.read_instance.species = self.read_instance.species * len(self.read_instance.network)

            # otherwise throw error
            else:
                error = 'Error: The number of "network" and "species" fields is not the same.'
                sys.exit(error)

        # throw error if one of networks are non all GHOST or non-GHOST
        for network_ii, network in enumerate(self.read_instance.network):
            if network_ii == 0:
                previous_is_ghost = check_for_ghost(network)
            else:
                is_ghost = check_for_ghost(network)
                if is_ghost != previous_is_ghost:
                    error = 'Error: "network" must be all GHOST or non-GHOST'
                    sys.exit(error)
                previous_is_ghost = is_ghost

        # check have resolution information, 
        # if offline, throw message, stating are using default instead
        if not self.read_instance.resolution:
            #default = 'monthly'
            default = 'hourly'
            msg = "Resolution (resolution) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.resolution = default

        # check have start_date information, 
        # if offline, throw message, stating are using default instead
        if not self.read_instance.start_date:
            default = '20180101'
            msg = "Start date (start_date) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.start_date = default

        # check have end_date information, 
        # if offline, throw message, stating are using default instead
        if not self.read_instance.end_date:
            default = '20190101'
            msg = "End date (end_date) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.end_date = default

        # check have statistic_mode information,
        # if offline, throw message, stating are using default instead
        if not self.read_instance.statistic_mode:
            default = 'Temporal|Spatial'
            msg = "Statistic mode (statistic_mode) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.statistic_mode = default

        # check have statistic_aggregation information,
        # if offline, throw message, stating are using default instead
        if not self.read_instance.statistic_aggregation:
            if self.read_instance.statistic_mode == 'Flattened':
                default = ''
            else:    
                default = 'Median'
                msg = "Statistic aggregation (statistic_aggregation) was not defined in the configuration file. Using '{}' as default.".format(default)
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.statistic_aggregation = default
        # if statistic_aggregation is defined ensure that it matches with the statistic_mode
        else:
            if (self.read_instance.statistic_mode == 'Flattened') & (self.read_instance.statistic_aggregation != ''):
                msg = "statistic_mode is set to be 'Flattened', therefore statistic_aggregation must be empty, not '{}'. Setting to be empty.".format(self.read_instance.statistic_aggregation)                
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                self.read_instance.statistic_aggregation = ''
            elif (self.read_instance.statistic_mode != 'Flattened') & (self.read_instance.statistic_aggregation == ''):
                default = 'Median'
                msg = "statistic_mode is set to be '{}', therefore statistic_aggregation must not be empty. Setting to be '{}'.".format(self.read_instance.statistic_mode, default)                
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                self.read_instance.statistic_aggregation = default

        # check have periodic_statistic_mode information,
        # if offline, throw message, stating are using default instead
        if not self.read_instance.periodic_statistic_mode:
            #default = 'Cycle'
            default = 'Independent'
            msg = "Periodic statistic mode (periodic_statistic_mode) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.periodic_statistic_mode = default

        # check have periodic_statistic_aggregation information,
        # if offline, throw message, stating are using default instead
        if not self.read_instance.periodic_statistic_aggregation:
            default = 'Median'
            msg = "Periodic statistic aggregation (periodic_statistic_aggregation) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.periodic_statistic_aggregation = default

        # check have timeseries_statistic_aggregation information,
        # if offline, throw message, stating are using default instead
        if not self.read_instance.timeseries_statistic_aggregation:
            default = 'Median'
            msg = "Timeseries statistic aggregation (timeseries_statistic_aggregation) was not "
            msg += "defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.timeseries_statistic_aggregation = default
        else:
            if ((self.read_instance.statistic_mode == 'Spatial|Temporal')
                and (self.read_instance.timeseries_statistic_aggregation != self.read_instance.statistic_aggregation)):
                msg = "Aggregation statistic and timeseries aggregation statistic are not "
                msg += "the same and Spatial|Temporal mode is active."
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)

        # check have correct active_dashboard_plots information, 
        # should have 4 plots if non-empty, throw error if using dashboard if not
        if not self.read_instance.active_dashboard_plots:
            default = ['timeseries', 'statsummary', 'distribution', 'periodic']
            self.read_instance.active_dashboard_plots = default
        # TODO: For Taylor diagrams, remove this piece of code when we stop using Matplotlib 3.3
        else:
            if Version(matplotlib.__version__) < Version("3.8"):
                if 'taylor' in self.read_instance.active_dashboard_plots:
                    error = 'It is not possible to create Taylor diagrams yet, please remove.'
                    sys.exit(error)

        if (len(self.read_instance.active_dashboard_plots) != 4) and (not self.read_instance.offline) & (not self.read_instance.interactive):
            error = 'Error: there must be 4 "active_dashboard_plots"'
            sys.exit(error)
        
        # if filter_species is active, and spatial_colocation is not active, then cannot filter by species
        # set filter_species to empty dict and advise user of this
        if (self.read_instance.filter_species) and (not self.read_instance.spatial_colocation):
            self.read_instance.filter_species = {}
            msg = 'Spatial colocation (spatial_colocation) must be set to True if wanting to filter by species.'
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)

        # map to multiple species if have * wildcard
        # also duplicate out associated network
        # remove any species for which there exists no data
        new_species = copy.deepcopy(self.read_instance.species)
        for speci_ii, speci in enumerate(self.read_instance.species): 
            if '*' in speci:
                mapped_species = self.multispecies_mapping(speci)
                del new_species[speci_ii]
                new_species[speci_ii:speci_ii] = mapped_species
                network_to_duplicate = self.read_instance.network[speci_ii]
                del self.read_instance.network[speci_ii]
                self.read_instance.network[speci_ii:speci_ii] = [network_to_duplicate]*len(mapped_species)
        self.read_instance.species = copy.deepcopy(new_species)

        # create variable for all unique species (plus filter species)
        filter_species = []
        species_plus_filter_species = copy.deepcopy(self.read_instance.species)
        if self.read_instance.filter_species:
            for networkspeci in self.read_instance.filter_species:
                speci = networkspeci.split('|')[1]
                if speci not in self.read_instance.species:
                    filter_species.append(speci)
                    species_plus_filter_species.append(speci)
                    
        # set lower_bound and upper_bound as dicts with limits per species (including filter species)
        # if type is dict then set bound using GHOST extreme limits per species

        # lower_bound
        # type is dict, then set as default limit per species using GHOST limits
        if isinstance(self.read_instance.lower_bound, dict):
            self.read_instance.lower_bound = {speci:
                                              np.float32(self.read_instance.parameter_dictionary[speci]['extreme_lower_limit']) 
                                              for speci in species_plus_filter_species}
        # otherwise set list values to dict, saving limits per species
        # if have just 1 limit apply for all read species, but if have multiple, set limits per species
        # throw error if have multiple lower bounds, but not equal to number of species to read  
        else:      
            lower_bound_dict = {}
            if len(self.read_instance.lower_bound) == 1:
                for speci in species_plus_filter_species:
                    lower_bound_dict[speci] = self.read_instance.lower_bound[0]
            elif len(self.read_instance.lower_bound) > 1:
                if len(self.read_instance.species) != len(self.read_instance.lower_bound):
                    error = 'Error: "lower_bound" variable must be same length as number of species read.'
                    sys.exit(error)
                else:
                    for speci_ii, speci in enumerate(self.read_instance.species):
                        lower_bound_dict[speci] = self.read_instance.lower_bound[speci_ii] 
                    # add filter_species (using GHOST limits)
                    for speci in filter_species:
                        lower_bound_dict[speci] = np.float32(self.read_instance.parameter_dictionary[speci]['extreme_lower_limit'])
            self.read_instance.lower_bound = lower_bound_dict

        # upper_bound
        # type is dict, then set as default limit per species using GHOST limits
        if isinstance(self.read_instance.upper_bound, dict):
            self.read_instance.upper_bound = {speci:np.float32(self.read_instance.parameter_dictionary[speci]['extreme_upper_limit']) 
                                              for speci in species_plus_filter_species}
        # otherwise set list values to dict, saving limits per species
        # if have just 1 limit apply for all read species, but if have multiple, set limits per species
        # throw error if have multiple upper bounds, but not equal to number of species to read  
        else:      
            upper_bound_dict = {}
            if len(self.read_instance.upper_bound) == 1:
                for speci in species_plus_filter_species:
                    upper_bound_dict[speci] = self.read_instance.upper_bound[0]
            elif len(self.read_instance.upper_bound) > 1:
                if len(self.read_instance.species) != len(self.read_instance.upper_bound):
                    error = 'Error: "upper_bound" variable must be same length as number of species read.'
                    sys.exit(error)
                else:
                    for speci_ii, speci in enumerate(self.read_instance.species):
                        upper_bound_dict[speci] = self.read_instance.upper_bound[speci_ii] 
                    # add filter_species (using GHOST limits)
                    for speci in filter_species:
                        upper_bound_dict[speci] = np.float32(self.read_instance.parameter_dictionary[speci]['extreme_upper_limit'])
            self.read_instance.upper_bound = upper_bound_dict

        # create a variable to set qa per species (including filter species)
        if isinstance(self.read_instance.qa, dict):
            self.read_instance.qa_per_species = {speci:get_default_qa(self.read_instance, speci) 
                                                 for speci in species_plus_filter_species}
            # set qa to be first of qa per species pairs
            self.read_instance.qa = self.read_instance.qa_per_species[list(self.read_instance.qa_per_species.keys())[0]]
        else:
            self.read_instance.qa_per_species = {speci:self.read_instance.qa for speci in species_plus_filter_species}

        # if are using dashboard then just take first network/species pair, as multivar not supported yet
        if ((len(self.read_instance.network) > 1) and (len(self.read_instance.species) > 1) and 
            (not self.read_instance.offline) and (not self.read_instance.interactive)):
             
            msg = 'Multiple networks/species are not supported in the dashboard. First ones will be taken.'
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)

            self.read_instance.network = [self.read_instance.network[0]]
            self.read_instance.species = [self.read_instance.species[0]]
        
        # check bounds inside filter_species
        if self.read_instance.filter_species:
            for networkspeci in self.read_instance.filter_species: 
                for networkspeci_limit_ii, networkspeci_limit in enumerate(self.read_instance.filter_species[networkspeci]):
                    
                    # get bounds
                    lower_limit = networkspeci_limit[0]
                    upper_limit = networkspeci_limit[1]
                    filter_species_fill_value = networkspeci_limit[2]
                    
                    # modify lower bound to be :, or contain > or >=
                    if ('<' in lower_limit):
                        msg = 'Lower bound ({}) for {} cannot contain < or <=. '.format(lower_limit, networkspeci)
                        lower_limit = '>=' + lower_limit.replace('<', '').replace('=', '')
                        msg += 'Setting it to be {}.'.format(lower_limit)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                    elif (':' not in lower_limit) and ('>' not in lower_limit):
                        msg = 'Lower bound ({}) for {} should contain > or >=. '.format(lower_limit, networkspeci)
                        lower_limit = '>=' + lower_limit
                        msg += 'Setting it to be {}.'.format(lower_limit)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)

                    # modify upper bound to be :, or contain < or <=
                    if ('>' in upper_limit):
                        msg = 'Upper bound ({}) for {} cannot contain > or >=. '.format(upper_limit, networkspeci)
                        upper_limit = '<=' + upper_limit.replace('>', '').replace('=', '')
                        msg += 'Setting it to be {}.'.format(upper_limit)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                    elif (':' not in upper_limit) and ('<' not in upper_limit):
                        msg = 'Upper bound ({}) for {} should contain < or <=. '.format(upper_limit, networkspeci)
                        upper_limit = '<=' + upper_limit
                        msg += 'Setting it to be {}.'.format(upper_limit)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                    
                    # update symbols next to values
                    self.read_instance.filter_species[networkspeci][networkspeci_limit_ii] = [lower_limit, upper_limit, 
                                                                                              filter_species_fill_value]

    def multispecies_mapping(self, species):
        """ Map species special case str to multiple species names. """

        multi_species_map = {'vconcaerobin*':['vconcaerobin1','vconcaerobin2','vconcaerobin3','vconcaerobin4',
                            'vconcaerobin5','vconcaerobin6','vconcaerobin7','vconcaerobin8','vconcaerobin9',
                            'vconcaerobin10','vconcaerobin11','vconcaerobin12','vconcaerobin13','vconcaerobin14',
                            'vconcaerobin15','vconcaerobin16','vconcaerobin17','vconcaerobin18','vconcaerobin19',
                            'vconcaerobin20','vconcaerobin21','vconcaerobin22']}

        return multi_species_map[species]


def read_conf(fpath=None):
    """ Read configuration files. """

    res = {}
    res_sub = {}
    config = {}
    all_sections = []
    all_sections_modified = []
    all_sections_commented = []
    repeated_subsections = []
    repeated_subsections_modified = {}
    subsections = []
    subsections_modified = []
    parent_sections = []
    filenames = []

    # get section names (e.g. [SECTIONA], [[Spain]]) and modified names (e.g. SECTIONA, SECTIONA-Spain)
    with open(fpath) as file:
        for line in file:
            if '[' in line and ']' in line and '[[' not in line and ']]' not in line:
                section = line.strip()
                #if first character is comment do not parse section
                if section[0] != '#':
                    section_modified = section.split('[')[1].split(']')[0]
                    if section_modified not in all_sections_modified:
                        parent_sections.append(section_modified)
                        all_sections_modified.append(section_modified)
                    else:
                        error = 'Error: It is not possible to have two sections with the same name.'
                        sys.exit(error)
            elif '[[' in line and ']]' in line:
                subsection = line.strip()
                # if first character is comment do not parse subsection
                if subsection[0] != '#':
                    subsection_modified = section_modified + '·' + line.split('[[')[1].split(']]')[0]
                    subsections.append(subsection)
                    subsections_modified.append(subsection_modified)
                    all_sections_modified.append(subsection_modified)

            if '[' in line and ']' in line:
                # if first character is comment then add section to list to avoid parsing
                if line.strip()[0] == '#':
                    all_sections_commented.append(line.strip())
                else:
                    all_sections.append(line.strip())

    # get repeated elements
    repetition_counts = {section:subsections_modified.count(section) for section in subsections_modified}
    for section, counts in repetition_counts.items():
        if counts > 1:
            repeated_subsections.append(section)
            repeated_subsections_modified[section] = [x for x in all_sections_modified 
                                                        if section.split('[[')[1].split(']]')[0] in x]

    # get attributes for each section and store in dict
    for (i, section), section_modified in zip(enumerate(all_sections), all_sections_modified):
        
        repetition = 0
        copy = False
        config[section_modified] = {}
        
        with open(fpath) as file:
            for line in file:
                line_strip = line.strip()
                # parsing all but last section 
                if section_modified != all_sections_modified[-1]:
                    # start of relevant section 
                    if line_strip == all_sections[i]:
                        if line_strip in repeated_subsections:
                            position = repeated_subsections_modified[section].index(section_modified)
                            if position == repetition:
                                copy = True
                            else:
                                copy = False
                            repetition += 1
                        else:
                            copy = True
                        continue
                    # start of next section, or commented section
                    elif (line_strip == all_sections[i+1]) or (line_strip in all_sections_commented):
                        copy = False
                        continue
                
                # parsing last section
                else:
                    # start of relevant section ?
                    if line_strip == all_sections[-1]:
                        if line_strip in repeated_subsections:
                            position = repeated_subsections_modified[section].index(section_modified)
                            if position == repetition:
                                copy = True
                            else:
                                copy = False
                            repetition += 1
                        else:
                            copy = True
                        continue

                    # start of commented section
                    elif line_strip in all_sections_commented:
                        copy = False
                        continue
                
                # set section attributes
                if copy:
                    if line_strip != '':
                        # initial definition of parameter - value
                        if '=' in line_strip:
                            key = line_strip.split('=', 1)[0].strip()
                            value = line_strip.split('=', 1)[1].strip()
                            config[section_modified][key] = value
                        
                        # lines after adding line breaks for long values
                        # make sure it is not a section name
                        elif ('=' not in line_strip) and (line_strip not in all_sections):
                            # get last key and add current value to values from last key
                            last_key = list(config[section_modified].keys())[-1]
                            value = config[section_modified][last_key] + line_strip.strip()

                            # update values for last key in dict
                            del config[section_modified][last_key]
                            config[section_modified][last_key] = value

    # add section attributes to subsection if do not exist there (e.g. add SECTIONA values to SECTIONA-Spain)
    for section_modified in all_sections_modified:
        
        # determine if subsection or not
        if '·' in section_modified:
            is_subsection = True
            par_section = section_modified.split('·')[0]
            # add attributes from parent section
            for par_k, par_val in config[par_section].items():
                # if first character of key is comment character (#), do not parse this attribute
                if par_k[0] == '#':
                    continue
                try:
                    res_sub[par_k] = eval(par_val)
                except:
                    res_sub[par_k] = par_val
        else:
            is_subsection = False

        # add attributes from current section/subsection
        for k, val in config[section_modified].items():
            if not is_subsection:
                # store filename
                if k == 'report_filename':
                    filenames.append(val)

            # if first character of key is comment character (#), do not parse this attribute
            if k[0] == '#':
                continue
            # overwrite attributes from current subsection
            try:
                res_sub[k] = eval(val)
            except:
                res_sub[k] = val

        # store pairs into res variable
        res[section_modified] = res_sub

        # reset res variable
        res_sub = {}

    return res, all_sections_modified, parent_sections, subsections_modified, filenames
   

def write_conf(section, subsection, fpath, opts):
    """ Write configurations on file. """

    config = configparser.RawConfigParser()

    # update configuration
    for section, section_name in zip(['section', 'subsection'], [section, subsection]):
        if opts[section]:
            config.add_section(section_name)
            for item in opts[section]:
                val = opts[section][item]
                config.set(section_name, item, val)

    # write configuration
    with open(fpath, 'w') as configfile:
        config.write(configfile)


def load_conf(self, fpath=None):
    """ Load existing configurations from file
        for running offline Providentia.
    """

    from .configuration import read_conf

    if fpath is None:
        print("No configuration file found")
        sys.exit(1)

    # if DEFAULT is not present, then return
    if not os.path.isfile(fpath):
        print(("Error %s" % fpath))
        return

    self.sub_opts, self.all_sections, self.parent_section_names, self.subsection_names, self.filenames = read_conf(fpath)


def split_options(read_instance, conf_string, separator="||"):
    """ For the options in the configuration that define the keep and remove
        options. Returns the values in two lists, the keeps and removes.
    """
    
    keeps, removes = [], []

    if separator not in conf_string:
        if ("keep:" in conf_string) and ("remove:" not in conf_string):
            keep_start = conf_string.find("keep:")
            keeps = conf_string[keep_start+5:]
            keeps = keeps.split(",")
            keeps = [k.strip() for k in keeps]
        elif ("keep:" not in conf_string) and ("remove:" in conf_string):
            remove_start = conf_string.find("remove:")
            removes = conf_string[remove_start+7:]
            removes = removes.split(",")
            removes = [r.strip() for r in removes]
        elif ("keep:" in conf_string) and ("remove:" in conf_string):
            msg = 'In order to define the keep and remove options, these must be separated by ||.'
            show_message(msg, from_conf=read_instance.from_conf)
    else:
        if "keep:" in conf_string:
            keep_start, keep_end = conf_string.find("keep:"), conf_string.find(separator)
            keeps = conf_string[keep_start+5:keep_end]
            keeps = keeps.split(",")
            keeps = [k.strip() for k in keeps]
        if "remove:" in conf_string:
            remove_start = conf_string.find("remove:")
            removes = conf_string[remove_start+7:]
            removes = removes.split(",")
            removes = [r.strip() for r in removes]
    
    return keeps, removes
