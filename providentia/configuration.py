""" Providentia Configuration Module """

import configparser
import copy
from datetime import datetime
import os
import platform
import socket
import sys
import yaml

import logging
import numpy as np
from packaging.version import Version
import pandas as pd
import ast

from providentia.auxiliar import CURRENT_PATH, join

MACHINE = os.environ.get('BSC_MACHINE', 'local')

# get current path and providentia root path
PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
data_paths = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'data_paths.yaml')))
default_values = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'prov_defaults.yaml')))
multispecies_map = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'multispecies_shortcurts.yaml')))
interp_experiments = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'interp_experiments.yaml')))

# set MACHINE to be the hub, workstation or local machine
if MACHINE not in ['nord3v2', 'mn5', 'nord4']:
    hostname = os.environ.get('HOSTNAME', '')
    ip = socket.gethostbyname(socket.gethostname())
    if "bscearth" in hostname:
        MACHINE = "workstation"
    elif "transfer" in hostname:
        MACHINE = "storage5"
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
        return join(dir, f)

class ProvConfiguration:
    """ Class that handles the configuration parameters definitions. """

    def __init__(self, read_instance, **kwargs):
        
        self.read_instance = read_instance 
        
        # set variable defaults
        self.var_defaults = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'init_prov_dev.yaml')))
        self.var_defaults['config_dir'] = join(PROVIDENTIA_ROOT, self.var_defaults['config_dir'])
        modifiable_var_defaults = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'init_prov.yaml')))
        self.var_defaults.update(modifiable_var_defaults)

        # if variable is given by command line, set that value, otherwise set as default value 
        for k, val in self.var_defaults.items():
            val = kwargs.get(k, val)
            setattr(self.read_instance, k, self.parse_parameter(k, val))

        # direct output to file/screen
        if hasattr(self.read_instance, 'logger') is False:
            self.switch_logging()

    def parse_parameter(self, key, value, deactivate_warning=False):
        """ Parse a parameter. """
        
        # import show_mesage from warnings
        sys.path.append(join(PROVIDENTIA_ROOT, 'providentia'))
        from warnings_prv import show_message
        
        # make sure we don't pass strings instead of booleans for true and false
        if value == 'true':
            value = True
        elif value == 'false':
            value = False

        # when the value is default then it is as if it was a blank value
        if value == 'default':
            return self.var_defaults[key]

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
                self.read_instance.logger.error(error)
                sys.exit(1)
            return operating_system
        
        elif key == 'machine':
            # set filetree type
            if MACHINE in ['nord3v2', 'mn5', 'dust', 'nord4']:
                self.read_instance.filetree_type = 'remote'
            else:
                self.read_instance.filetree_type = 'local'
            return MACHINE

        elif key == 'available_cpus':
            # get available N CPUs
            if MACHINE in ['nord3v2', 'mn5', 'nord4']:
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
            # set cartopy data directory (needed on Nord3v2/MN5 as has no external
            # internet connection)
            if MACHINE == 'nord3v2':
                return '/gpfs/projects/bsc32/software/rhel/7.5/ppc64le/POWER9/software/Cartopy/0.17.0-foss-2018b-Python-3.7.0/lib/python3.7/site-packages/Cartopy-0.17.0-py3.7-linux-ppc64le.egg/cartopy/data'
            elif MACHINE in ['nord4', 'mn5']:
                return '/gpfs/projects/bsc32/software/rhel/9.2/software/Cartopy/0.23.0-foss-2023b-Python-3.11.5/lib/python3.11/site-packages/cartopy/data'
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

        elif key == 'exp_to_interp_root':
            # define experiment root data directory
            # set experiment root data directory if left undefined
            if value == '':
                exp_to_interp_root = data_paths[MACHINE]["exp_to_interp_root"]
                if exp_to_interp_root != '': 
                    return os.path.expanduser(exp_to_interp_root[0])+exp_to_interp_root[1:]       
        
        elif key == 'ghost_version':
            # parse GHOST version

            # import GHOST standards 
            sys.path = [path for path in sys.path if 'dependencies/GHOST_standards/' not in path]            
            sys.path.insert(1, join(CURRENT_PATH, 'dependencies/GHOST_standards/{}'.format(value)))
            if 'GHOST_standards' in sys.modules:
                del sys.modules['GHOST_standards']
            from GHOST_standards import standard_parameters
            from GHOST_standards import get_standard_metadata
            from GHOST_standards import standard_data_flag_name_to_data_flag_code
            from GHOST_standards import standard_QA_name_to_QA_code
            from GHOST_standards import standard_networks
            from GHOST_standards import standard_temporal_resolutions

            # get ghost_version list
            self.read_instance.possible_ghost_versions = os.listdir(join(CURRENT_PATH,'dependencies', 'GHOST_standards'))
            
            # get GHOST networks
            self.read_instance.ghost_available_networks = list(standard_networks.keys())

            # get GHOST resolutions
            self.read_instance.ghost_available_resolutions = [resolution_dict['temporal_resolution_path'] for resolution_dict in standard_temporal_resolutions.values()]

            # modify standard parameter dictionary to have BSC standard parameter names as keys (rather than GHOST)
            self.read_instance.parameter_dictionary = dict()
            for _, param_dict in standard_parameters.items():
                self.read_instance.parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict

            # get available species
            self.read_instance.available_species = list(self.read_instance.parameter_dictionary.keys())
            
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
                # treat leaving the field blank as default
                if value == '':
                    return self.var_defaults[key]
                # parse multiple networks
                if ',' in value:
                    return [network.strip() for network in value.split(',')]
                else:
                    return [value.strip()]

        elif key == 'species':
            # parse species

            if isinstance(value, str):
                # treat leaving the field blank as default
                if value == '':
                    return self.var_defaults[key]
                # parse multiple species
                if ',' in value:
                    return [speci.strip() for speci in value.split(',')]
                else:
                    return [value.strip()]

        elif key == 'resolution':
            # parse resolution
            
            if isinstance(value, str):
                # treat leaving the field blank as default
                if value == '':
                    return self.var_defaults[key]
                # parse multiple species only in interpolation and download
                if self.read_instance.interpolation or self.read_instance.download:
                    return [res.strip() for res in value.split(',')]
                else:        
                    return value.strip()

        elif key == 'start_date':
            # parse start_date

            if (isinstance(value, str)) or (isinstance(value, int)):
                # treat leaving the field blank as default
                if value == '':
                    return self.var_defaults[key]
                # throw error if start_date is empty str
                value = str(value)
                return value.strip()

        elif key == 'end_date':
            # parse end_date

            if (isinstance(value, str)) or (isinstance(value, int)):
                # treat leaving the field blank as default
                if value == '':
                    return self.var_defaults[key]
                # throw error if start_date is empty str
                value = str(value)
                return value.strip()
        
        elif key == 'qa':
            # parse qa

            from GHOST_standards import providentia_defaults

            # set default qa codes (can differ per GHOST version)
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

            from GHOST_standards import providentia_defaults

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
                    return sorted([self.read_instance.standard_data_flag_name_to_data_flag_code[f.strip()] for f in value.split(",")])
                # list of integer codes
                else:
                    return sorted(list(value))
            # otherwise, set default flags
            else:
                return sorted([self.read_instance.standard_data_flag_name_to_data_flag_code[flag_name] for flag_name in providentia_defaults['flag']])

        elif key in ['add_qa','subtract_qa']:
            # parse add/subtract qa

            # if not None then set QA by that given
            if value is not None:
                # if conf has only one QA
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
            # otherwise, return empty list
            else:
                return []

        elif key in ['add_flags','subtract_flags']:
            # parse add/subtract flags

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
                    # check if all the flags appear in the GHOST_standards of the current version
                    for flag in value.split(","):
                        if flag.strip() not in self.read_instance.standard_data_flag_name_to_data_flag_code:
                            error = (f"Error: Flag '{flag}' not in this GHOST version ({self.read_instance.ghost_version}).")
                            self.read_instance.logger.error(error)
                            sys.exit(1)
                    return sorted([self.read_instance.standard_data_flag_name_to_data_flag_code[f.strip()] for f in value.split(",")])
                # list of integer codes
                else:
                    return sorted(list(value))
            # otherwise, return empty list
            else:
                return []

        elif key == 'domain': # TODO maybe there's no need of domain because it is included on experiments
            # parse domain

            if value is not None:
                # treat leaving the field blank as default
                if value == '':
                    return self.var_defaults[key]
                # split list, if only one domain, then creates list of one element
                domains = [dom.strip() for dom in value.split(",")]      
                return domains
            else:
                return []

        elif key == 'ensemble_options': # TODO maybe there's no need of ensemble num because it is included on experiments
            # parse ensemble_options

            if value is not None:
                # treat leaving the field blank as default
                if value == '':
                    return self.var_defaults[key]

                # split list, if only one ensemble_options, then creates list of one element
                ensemble_opts = []
                for opt in str(value).split(","):
                    # if it is a number, then make it 3 digits, if not it stays as it is
                    if opt.isdigit():
                        opt = opt.strip().zfill(3)
                    # check that it does not start with stat
                    elif opt.startswith('stat'):
                        error = "Error: ensemble option cannot start with 'stat'.\n" \
                        "For ensemble statistics, simply define them based on the stat name provided in the filename, such as:\n" \
                        "   · 'av' for ensemble average\n" \
                        "   · 'av_an' for ensemble analysis average"
                        self.read_instance.logger.error(error)
                        sys.exit(1)
                    
                    ensemble_opts.append(opt)

                ensemble_opts = [opt.strip().zfill(3) if opt.isdigit() else opt for opt in str(value).split(",")]  
                return ensemble_opts
            else:
                return []
            
        elif key == 'experiments':
            # parse experiments

            if isinstance(value, str):
                # empty string
                if value == "":
                    return []
                # split experiments
                else:
                    # have alternative experiment names for the legend, then parse them?
                    if ('(' in value) & (')' in value):
                        exps = [exp.strip() for exp in value.split('(')[0].strip().split(",")]
                        self.read_instance.alias = [exp_legend.strip() for exp_legend in value.split('(')[1].split(')')[0].strip().split(",")]
                    # otherwise set legend names as given experiment names in full
                    else: 
                        exps = [exp.strip() for exp in value.split(",")]
                        self.read_instance.alias = []

                    return exps

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

        elif key == 'calibration_factor':
            # parse calibration factor

            if isinstance(value, (str)):

                # convert to string if not
                if np.issubdtype(type(value), np.number):
                    value = str(value)

                # strip all whitespace
                if ',' in value:
                    return [calibration_factor.strip() for calibration_factor in value.split(',')]
                else:
                    return [value.strip()]
        
        elif key in ['statistic_mode','statistic_aggregation','periodic_statistic_mode','periodic_statistic_aggregation',
                     'timeseries_statistic_aggregation','interp_n_neighbours','interp_reverse_vertical_orientation',
                     'interp_chunk_size','interp_job_array_limit', 'interp_multiprocessing']:
            # treat leaving the field blank as default
            if value == '':
                return self.var_defaults[key]
            
        # if no special parsing treatment for variable, simply return value
        return value

    def decompose_experiments(self):
        """ Get experiment components (experiment-domain-ensemble_options) and fill the class variables with their value."""

        # get possible domains
        possible_domains = default_values["domain"]

        # get separated experiment parts list
        split_experiments = [exp.split("-") for exp in self.read_instance.experiments]

        # check if all the experiments are written in the same way
        if len(set([len(exp) for exp in split_experiments])) > 1:
            error = (f"Error: All the experiments have to follow the same structure: [expID], [expID]-[domain], [expID]-[ensembleNum] or [expID]-[domain]-[ensembleNum].")
            self.read_instance.logger.error(error)
            sys.exit(1)

        # get original domain and ensemble options as passed in the configuration file
        config_domain = copy.deepcopy(self.read_instance.domain) 
        config_ensemble_options = copy.deepcopy(self.read_instance.ensemble_options)

        # initialize experiment id, domain and ensemble options for each of the experiments
        exp_id, exp_dom, exp_ens = None, None, None

        # initialize lists to hold domains/ensemble options inside the experiments
        exp_domains_list = []
        exp_ensemble_options_list = []
        exp_ids_list = []

        # iterate through all the experiments
        for exp_i, split_experiment in enumerate(split_experiments):
            # get experiment name and save into list
            exp_id = split_experiment[0]
            exp_ids_list.append(exp_id)

            # [expID]-[domain] or [expID]-[ensembleNum] 
            if len(split_experiment) == 2:
                end_experiment = split_experiment[-1]
                
                # [expID]-[domain]
                if end_experiment in possible_domains: 
                    # other experiment goes by the format [expID]-[ensembleNum]
                    if exp_ens:
                        error = (f"Error: All the experiments have to follow the same structure: [expID], [expID]-[domain], [expID]-[ensembleNum] or [expID]-[domain]-[ensembleNum].")
                        self.read_instance.logger.error(error)
                        sys.exit(1)
                    exp_dom = end_experiment
                    exp_domains_list.append(exp_dom)
                # [expID]-[ensembleNum]   
                else:
                    # other experiment goes by the format [expID]-[domain]
                    if exp_dom:
                        error = (f"Error: All the experiments have to follow the same structure: [expID], [expID]-[domain], [expID]-[ensembleNum] or [expID]-[domain]-[ensembleNum].")
                        self.read_instance.logger.error(error)
                        sys.exit(1)

                    exp_ens = end_experiment
                    # if it is a number, then make it 3 digits, if not it stays as it is
                    if exp_ens.isdigit():
                        exp_ens = exp_ens.strip().zfill(3)
                    # check that it does not start with stat
                    elif exp_ens.startswith('stat'):
                            error = f"Error: ensemble option {exp_ens} cannot start with 'stat'.\n" \
                            "For ensemble statistics, simply define them based on the stat name provided in the filename, such as:\n" \
                            "   · 'av' for ensemble average\n" \
                            "   · 'av_an' for ensemble analysis average"
                            self.read_instance.logger.error(error)
                            sys.exit(1)
                    exp_ensemble_options_list.append(exp_ens)

            # [expID]-[domain]-[ensembleNum]
            elif len(split_experiment) == 3:               
                exp_dom, exp_ens = split_experiment[1], split_experiment[2]
                exp_domains_list.append(exp_dom)
                # if it is a number, then make it 3 digits, if not it stays as it is
                if exp_ens.isdigit():
                    exp_ens = exp_ens.strip().zfill(3)
                # check that it does not start with stat
                elif exp_ens.startswith('stat'):
                        error = f"Error: ensemble option {exp_ens} cannot start with 'stat'.\n" \
                        "For ensemble statistics, simply define them based on the stat name provided in the filename, such as:\n" \
                        "   · 'av' for ensemble average\n" \
                        "   · 'av_an' for ensemble analysis average"
                        self.read_instance.logger.error(error)
                        sys.exit(1)
                exp_ensemble_options_list.append(exp_ens)
                        
            # if experiment is composed by more than 3 parts, exit
            elif len(split_experiment) > 3:
                error = 'Invalid experiment format, experiments have to consist of three elements maximum.'
                self.read_instance.logger.error(error)
                sys.exit(1)

            # throw error if domain has been defined in configuration file and in experiment name
            if exp_dom and config_domain:
                error = f"Error: Unable to set domain(s) as {', '.join(config_domain)} because the "
                error += f"experiment {self.read_instance.experiments[exp_i]} already contains the domain."
                self.read_instance.logger.error(error)
                sys.exit(1)
            # if there is no domain, fill it with the list from the experiments names
            elif not config_domain:
                self.read_instance.domain = exp_domains_list
            
            # throw error if ensemble options has been defined in configuration file and in experiment name
            if exp_ens and config_ensemble_options:
                error = f"Error: Unable to set ensemble option(s) as {', '.join(config_ensemble_options)} because the "
                error +=  f"experiment {self.read_instance.experiments[exp_i]} already contains the ensemble option."                  
                self.read_instance.logger.error(error)
                sys.exit(1)
            # if there is no ensemble option, fill it with the list from the experiments names
            elif not config_ensemble_options:
                self.read_instance.ensemble_options = exp_ensemble_options_list
            
            # add experiment id to the experiment ids list
            self.read_instance.exp_ids = exp_ids_list

        # when there's no domain/ensemble opt passed in the config file or got from the experiment, then set the default option to true
        if not self.read_instance.domain:
            self.default_domain = True
        if not self.read_instance.ensemble_options:
            self.default_ensemble_options = True

        # set the bool which tells you if domain/ensemble_options have to be combined as it was done in interpolation mode or not
        self.combined_domain, self.combined_ensemble_options = bool(config_domain or self.default_domain), bool(config_ensemble_options or self.default_ensemble_options)

    def check_experiment(self, full_experiment, deactivate_warning):
        # TODO Check if i can only import one time
        from warnings_prv import show_message

        """ Check individual experiment and get list of options. """
        # split full experiment
        experiment, domain, ensemble_option = full_experiment.split('-')
        
        # get all possible experiments
        exp_path = join(self.read_instance.exp_root,self.read_instance.ghost_version)
        self.possible_experiments = [] if not os.path.exists(exp_path) else os.listdir(exp_path)

        # initialise list of possible ghost versions
        available_ghost_versions = []

        # remove possible ghost versions if they are not really in the directories
        possible_ghost_versions = list(set(os.listdir(self.read_instance.exp_root)) & set(self.read_instance.possible_ghost_versions))

        # if ensemble options is allmembers, get all the possible ensemble options
        if ensemble_option == "allmembers":
            exp_found = list(filter(lambda x:x.startswith(experiment+'-'+domain), self.possible_experiments))
           
            # search for other ghost versions
            if not exp_found:
                for ghost_version in possible_ghost_versions:
                    ghost_exp_found = list(filter(lambda x:x.startswith(experiment+'-'+domain), os.listdir(join(self.read_instance.exp_root,ghost_version))))
                    if ghost_exp_found:
                        available_ghost_versions.append(ghost_version)
        # if it is a concrete ensemble option, then just get the experiment from the list
        else:
            exp_found = [full_experiment] if full_experiment in self.possible_experiments else []

            # search for other ghost versions
            if not exp_found:
                available_ghost_versions = list(filter(lambda x:full_experiment in os.listdir(join(self.read_instance.exp_root,x)), possible_ghost_versions))
        
        # if not found because of the ghost version, tell the user
        if available_ghost_versions and ('/' not in self.read_instance.network[0]):
            msg = f"There is no data available for {full_experiment} experiment for the current"
            msg += f" GHOST version ({self.read_instance.ghost_version}). Please check one of the available versions:"
            msg += f" {', '.join(sorted(available_ghost_versions))}"
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)

        return bool(exp_found), exp_found
    
    # TODO use inheritance in the future 
    # TODO add more checking, for now only this is enough
    def check_experiment_interpolation(self, full_experiment, deactivate_warning):
        # TODO Check if i can only import one time
        from warnings_prv import show_message

        """ Checks if experiment, domain and ensemble option combination works 
        for interpolation or the download of non interpolated experiments
        Returns if experiment if valid and the experiment type (if there is one) """
        
        # get the splitted experiment
        experiment, domain, ensemble_option = full_experiment.split('-')

        # accept asterisk to download all experiments
        if experiment == '*':
            return True, experiment
        
        # search if the experiment id is in the interp_experiments file
        # initialize experiment search variables
        experiment_exists = False
        msg = ""

        # for HPC machines, search in interp_experiments
        # if it's local interpolation, don't enter
        if not (self.read_instance.machine == "local" and self.read_instance.interpolation is True):
            for experiment_type, experiment_dict in interp_experiments.items():
                if experiment in experiment_dict["experiments"]:
                    experiment_exists = True
                    break
            
            msg += f"Cannot find the experiment ID '{experiment}' in '{join('settings', 'interp_experiments.yaml')}'. Please add it to the file. "

        # get directory from data_paths if it doesn't exists in the interp_experiments file 
        # if executed from the hpc machines and want to do a download, don't enter
        if experiment_exists is False and not (self.read_instance.machine != "local" and self.read_instance.download is True):
            # get the path to the non interpolated experiments
            # in the current machine if it is an intepolation
            if self.read_instance.interpolation is True:
                exp_to_interp_path = join(self.read_instance.exp_to_interp_root, experiment)
                if os.path.exists(exp_to_interp_path):
                    experiment_exists = True
            # in the remote machine if it is a local download
            else:
                # connect to the remote machine
                self.read_instance.connect()        
                # get all possible experiments
                exp_to_interp_path = join(self.read_instance.exp_to_interp_remote_path,experiment,domain)
                try:
                    self.read_instance.sftp.stat(exp_to_interp_path)
                    experiment_exists = True
                except FileNotFoundError:
                    pass     
            
            msg += f"Cannot find the {experiment} experiment with the {domain} domain in '{self.read_instance.exp_to_interp_root}'."

        # if experiment does not exist, exit
        # supressed warning deactivation
        if experiment_exists is False:
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)

        return experiment_exists, [full_experiment]
    
    # TODO maybe remove this one and keep the download check since its much cleaner
    def check_experiment_download(self, full_experiment, deactivate_warning):
        """ Check individual experiment and get list of options."""

        # TODO Check if i can only import one time
        from warnings_prv import show_message

        # split full experiment
        experiment, domain, ensemble_option = full_experiment.split('-')

        # accept asterisk to download all experiments
        if experiment == '*':
            return True, experiment
        
        # all experiments pass this check because the real one is in the remote machine
        exp_found = [full_experiment]
        
        # connect to the remote machine
        self.read_instance.connect()        
        
        # get all possible experiments
        exp_path = join(self.read_instance.exp_remote_path,self.read_instance.ghost_version)
        self.possible_experiments = self.read_instance.sftp.listdir(exp_path)

        # TODO repeated code, put this into a method in the future?
        # if ensemble options is allmembers, get all the possible ensemble options
        if ensemble_option == "allmembers":
            exp_found = list(sorted(filter(lambda x:x.startswith(experiment+'-'+domain), self.possible_experiments)))
           
            if not exp_found:
                # initialise list of possible ghost versions
                available_ghost_versions = []
                
                # search for other ghost versions
                for ghost_version in self.read_instance.possible_ghost_versions:
                    ghost_exp_found = list(filter(lambda x:x.startswith(experiment+'-'+domain), self.read_instance.sftp.listdir(join(self.read_instance.exp_remote_path,ghost_version))))
                    if ghost_exp_found:
                        available_ghost_versions.append(ghost_version)

                msg = f"There is no experiment {experiment}-{domain} data for the current ghost version ({self.read_instance.ghost_version})." 
                if available_ghost_versions:
                    msg += f" Please, check one of the available versions: {', '.join(sorted(available_ghost_versions))}"
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)

        return bool(exp_found), exp_found        

    def check_validity(self, deactivate_warning=False):
        """ Check validity of set variables after parsing. """
       
        # import matplotlib for Taylor Diagrams
        import matplotlib

        # import check_for_ghost and get_default_qa
        from read_aux import check_for_ghost, get_default_qa

        # import show_mesage from warnings
        from warnings_prv import show_message

        # check if species is valid
        if self.read_instance.species:
            for speci in self.read_instance.species:
                if '*' not in speci and speci not in self.read_instance.parameter_dictionary:
                    error = f'Error: species "{speci}" not valid.'
                    self.read_instance.logger.error(error)
                    sys.exit(1)

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
        # in download mode of non interpolated experiments is allowed to not have network, so continue
        if not self.read_instance.network and not (self.read_instance.download and not self.read_instance.interpolated):
            # default = ['GHOST']
            if self.read_instance.interpolation:
                default = self.read_instance.ghost_available_networks
            else:
                default = default_values['network']
            msg = "Network (network) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
            self.read_instance.network = default

        # check have species information, TODO REFACTOR INTERPOLATION STUFF
        # if offline, throw message, stating are using default instead
        # in download mode is allowed to not pass any species, so continue
        if (not self.read_instance.species and not self.read_instance.download and not self.read_instance.interpolation):
            if self.read_instance.interpolation:
                default = self.read_instance.available_species
            else:
                default = default_values['species']
            msg = "Species (species) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
            self.read_instance.species = default

        # if number of networks and species is not the same,
        # and len of one of network or species == 1,
        # then duplicate respestive network/species
        # in download mode is allowed to not have a different number, so continue
        # TODO in download mode and interpolation separate this somehow.
        if len(self.read_instance.network) != len(self.read_instance.species) and not (self.read_instance.download or self.read_instance.interpolation):

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
                self.read_instance.logger.error(error)
                sys.exit(1)

        # throw error if one of networks are non all GHOST or non-GHOST
        # in download mode it is allowed to have mixed networks #TODO Change this
        if not self.read_instance.download and not self.read_instance.interpolation:
            for network_ii, network in enumerate(self.read_instance.network):
                if network_ii == 0:
                    previous_is_ghost = check_for_ghost(network)
                else:
                    is_ghost = check_for_ghost(network)
                    if is_ghost != previous_is_ghost:
                        error = 'Error: "network" must be all GHOST or non-GHOST'
                        self.read_instance.logger.error(error)
                        sys.exit(1)
                    previous_is_ghost = is_ghost

        # check have resolution information, TODO when refactoring init change this way of checking defaults
        # if offline, throw message, stating are using default instead
        if (not self.read_instance.resolution and not self.read_instance.download):
            if self.read_instance.interpolation:
                default = ["hourly", "hourly_instantaneous", "daily", "monthly"]
            else:
                #default = 'monthly'
                default = default_values['resolution']
            msg = "Resolution (resolution) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
            self.read_instance.resolution = default

        # check start_date format, TODO START DATE IS DIFFERENT IN INTERPOLATION (check in the refactoring)
        # if offline, throw message, stating are using default instead
        if not self.read_instance.start_date:
            if self.read_instance.interpolation:
                default = '201801'
            else:
                default = default_values['start_date']
            msg = "Start date (start_date) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
            self.read_instance.start_date = default
        else:
            len_start_date = len(self.read_instance.start_date)
            if self.read_instance.interpolation:
                if len_start_date != 6:
                    if len_start_date == 8:
                        self.read_instance.start_date = self.read_instance.start_date[:-2]
                        msg = "Start date (start_date) was defined as YYYYMMDD, changing it to YYYYMM. Using '{}'.".format(self.read_instance.start_date)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
                    else:
                        error = "Error: Format of Start date (start_date) not correct, please change it to YYYYMM."
                        self.read_instance.logger.error(error)
                        sys.exit(1)
            else:
                if len_start_date != 8:
                    if len_start_date == 6:
                        self.read_instance.start_date = self.read_instance.start_date + "01"
                        msg = "Start date (start_date) was defined as YYYYMM, changing it to YYYYMMDD. Using '{}'.".format(self.read_instance.start_date)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
                    else:
                        error = "Error: The format of Start date (start_date) is not correct, please change it to YYYYMMDD."
                        self.read_instance.logger.error(error)
                        sys.exit(1)

        # check end_date  format, TODO START DATE IS DIFFERENT IN INTERPOLATION (check in the refactoring)
        # if offline, throw message, stating are using default instead
        if not self.read_instance.end_date:
            if self.read_instance.interpolation:
                default = '201901'
            else:
                default = default_values['end_date']
            msg = "End date (end_date) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
            self.read_instance.end_date = default
        else:
            len_end_date = len(self.read_instance.end_date)
            if self.read_instance.interpolation:
                if len_end_date != 6:
                    if len_end_date == 8:
                        self.read_instance.end_date = self.read_instance.end_date[:-2]
                        msg = "End Date (end_date) was defined as YYYYMMDD, changing it to YYYYMM. Using '{}'.".format(self.read_instance.end_date)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
                    else:
                        error = "Error: Format of End Date (end_date) not correct, please change it to YYYYMM."
                        self.read_instance.logger.error(error)
                        sys.exit(1)
            else:
                if len_end_date != 8:
                    if len_end_date == 6:
                        self.read_instance.end_date = self.read_instance.end_date + "01"
                        msg = "End Date (end_date) was defined as YYYYMM, changing it to YYYYMMDD. Using '{}'.".format(self.read_instance.end_date)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
                    else:
                        error = "Error: The format of End Date (end_date) is not correct, please change it to YYYYMMDD."
                        self.read_instance.logger.error(error)
                        sys.exit(1)

        # check have interp_n_neighbours information, TODO ONLY FOR INTERPOLATION
        # if offline, throw message, stating are using default instead
        # TODO CHANGE THE MESSAGE WHEN ITS DEFAULT AND WHEN ITS BECAUSE I DIDNT PUT THE NAME
        if self.read_instance.interpolation and not self.read_instance.interp_n_neighbours:
            default = default_values['interp_n_neighbours']
            msg = "Number of neighbours (interp_n_neighbours) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
            self.read_instance.interp_n_neighbours = default
    
        # TODO MAYBE CHANGE THIS initialization to somewhere else or take it from another place
        # TODO and should it be provconf or no????
        # initialise possible domains
        self.default_domain = False
        self.default_ensemble_options = False

        # make sure there are experiments in interpolation
        if self.read_instance.interpolation and (len(self.read_instance.experiments) == 0):
            error = 'Error: No experiments were provided in the configuration file.'
            self.read_instance.logger.error(error)
            sys.exit(1)

        # get domain, ensemble options, experiment ids and flag to get the default values of these variables
        self.decompose_experiments()

        # check have domain information, TODO ONLY FOR INTERPOLATION
        # TODO think if we need one separated variable for this one because it is already included on experiments
        # if offline, throw message, stating are using default instead
        if self.read_instance.experiments and self.default_domain:
            default = default_values['domain']
            msg = "Domain (domain) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.domain = default

        # check have ensemble_options information, TODO ONLY FOR INTERPOLATION
        # TODO think if we need one separated variable for this one because it is already included on experiments
        # if offline, throw message, stating are using default instead
        # TODO maybe think this a bit better, if i dont pass it it should check better if i already passed it in experiments and so
        if self.read_instance.experiments and self.default_ensemble_options:
            if self.read_instance.interpolation:
                default = ["000"]
            else:
                default = default_values['ensemble_options']
            msg = "Ensemble options (ensemble_options) was not defined in the configuration file. Using '{}' as default.".format(default)
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
            self.read_instance.ensemble_options = default

        # check if alias can be set (in case there is an alias)
        if self.read_instance.alias:
            # to set an alias is mandatory to have the same number of experiments and legends
            if len(self.read_instance.experiments)==len(self.read_instance.alias): 
                # if all experiments are full length ([expID]-[domain]-[ensembleNum]) or there's only one experiment with only one possible combination,
                # then they can be set as alias     
                if all([len(exp.split("-"))==3 for exp in self.read_instance.experiments]) or \
                    (len(self.read_instance.experiments) == 1 and len(self.read_instance.domain) == 1 and len(self.read_instance.ensemble_options) == 1):
                    self.read_instance.alias_flag = True
            
            # show warning if alias not possible to be set
            if not deactivate_warning and not self.read_instance.alias_flag:
                msg = "Experiment alias could not be set."
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)

        # before checking the experiment check that the remote download has the interpolated tag as False, if not exit
        if self.read_instance.download and MACHINE in ["storage5", "nord3v2", "nord4"] and self.read_instance.interpolated is True:
            error = F"Error: Nothing from the {self.read_instance.section} section was copied to gpfs, change the interpolated field to 'False'."
            self.read_instance.logger.error(error)
            sys.exit(1)

        # get correct check experiment function
        # TODO do it using heritage
        # if the current mode is interpolation or the experiment i want to download is not interpolated
        if self.read_instance.interpolation or (self.read_instance.download and self.read_instance.interpolated is False):
            check_experiment_fun = self.check_experiment_interpolation
        elif self.read_instance.download:
            check_experiment_fun = self.check_experiment_download
        else:
            check_experiment_fun = self.check_experiment

        # temp dictionary to store experiments
        # TODO change names
        final_experiments = []
        correct_experiments = {}

        # TODO keep the check only in download or configuration
        # in case of zenodo download don't even check
        if not (self.read_instance.download and self.read_instance.bsc_download_choice == "n"): 
            # join experiments
            for exp_i, experiment in enumerate(self.read_instance.exp_ids):
                # experiment, domain, ensemble_options
                if self.combined_domain and self.combined_ensemble_options:
                    final_experiments += [f'{experiment}-{domain}-{ens_opt}' for domain in self.read_instance.domain for ens_opt in self.read_instance.ensemble_options]
                else:
                    if self.combined_domain or self.combined_ensemble_options:
                        # experiment-ensemble_options, domain
                        if self.combined_domain:
                            final_experiments += [f'{experiment}-{domain}-{self.read_instance.ensemble_options[exp_i]}' for domain in self.read_instance.domain]
                        # experiment-domain, ensemble_options 
                        else:
                            final_experiments += [f'{experiment}-{self.read_instance.domain[exp_i]}-{ens_opt}' for ens_opt in self.read_instance.ensemble_options]
                    # experiment-domain-ensemble_options
                    else:
                        final_experiments.append(f'{experiment}-{self.read_instance.domain[exp_i]}-{self.read_instance.ensemble_options[exp_i]}')

            for exp_i, experiment in enumerate(final_experiments):
                # TODO change boolean name
                exp_is_valid, valid_exp_list = check_experiment_fun(experiment, deactivate_warning)
                if exp_is_valid:
                    for valid_exp in valid_exp_list:
                        if self.read_instance.alias_flag:
                            correct_experiments[valid_exp] = self.read_instance.alias[exp_i]
                        else:
                            correct_experiments[valid_exp] = valid_exp
        
        # if experiments were passed and there's no valid experiment, show warning
        if self.read_instance.experiments != [] and correct_experiments == {}:
            msg = 'No experiment data available.'
            if self.read_instance.interpolation:
                error = "Error: " + msg
                self.read_instance.logger.error(error)
                sys.exit(1)
            else:
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)

        # replace experiments by new ones found
        self.read_instance.experiments = correct_experiments

        # check calibration factor
        if self.read_instance.calibration_factor:

            # detect if calibration factor is passed by experiment
            calibration_by_experiment = not self.read_instance.calibration_factor[0][0] in ['+', '-', '*', '/']

            # control that calibration factor not by experiment can only be one element
            if not calibration_by_experiment and len(self.read_instance.calibration_factor) > 1:
                error = "Error: When calibration factor is not provided by the experiment, only one value can be passed."
                self.read_instance.logger.error(error)
                sys.exit(1)

            # create dictionary per experiment
            calibration_factor_dict = {}

            # if calibration is by experiment
            if calibration_by_experiment:
                for i, experiment in enumerate(self.read_instance.experiments):
                    for calibration_factor in self.read_instance.calibration_factor:
                        if experiment in calibration_factor:
                            calibration_factor_exp = calibration_factor.split("(")[1][:-1]
                            calibration_factor_dict[experiment] = calibration_factor_exp
            # if the same calibration is applied to all experiments
            else:
                calibration_factor_dict = {experiment:self.read_instance.calibration_factor[0] for experiment in self.read_instance.experiments}                 

            # replace calibration factors by new dictionary
            self.read_instance.calibration_factor = calibration_factor_dict

        # check have statistic_mode information,
        # if offline, throw message, stating are using default instead
        # TODO not needed in interpolation 
        if not self.read_instance.statistic_mode and not self.read_instance.interpolation:
            default = default_values['statistic_mode']
            self.read_instance.statistic_mode = default

        # check have statistic_aggregation information,
        # if offline, throw message, stating are using default instead
        # TODO not needed in interpolation 
        if not self.read_instance.interpolation:
            default = default_values['statistic_aggregation'][self.read_instance.statistic_mode]
            if not self.read_instance.statistic_aggregation:  
                self.read_instance.statistic_aggregation = default
            # if statistic_aggregation is defined ensure that it matches with the statistic_mode
            else:
                if self.read_instance.statistic_mode == 'Flattened':
                    self.read_instance.statistic_aggregation = default

        # check have periodic_statistic_mode information,
        # if offline, throw message, stating are using default instead
        # TODO not needed in interpolation 
        if not self.read_instance.periodic_statistic_mode and not self.read_instance.interpolation:
            #default = 'Cycle'
            default = default_values['periodic_statistic_mode']
            self.read_instance.periodic_statistic_mode = default

        # check have periodic_statistic_aggregation information,
        # if offline, throw message, stating are using default instead
        # TODO not needed in interpolation 
        if not self.read_instance.periodic_statistic_aggregation and not self.read_instance.interpolation:
            default = default_values['periodic_statistic_aggregation']
            self.read_instance.periodic_statistic_aggregation = default

        # check have timeseries_statistic_aggregation information,
        # if offline, throw message, stating are using default instead
        # TODO not needed in interpolation 
        if not self.read_instance.timeseries_statistic_aggregation and not self.read_instance.interpolation:
            default = default_values['timeseries_statistic_aggregation']
            self.read_instance.timeseries_statistic_aggregation = default

        # check have correct active_dashboard_plots information, 
        # should have 4 plots if non-empty, throw error if using dashboard if not
        if not self.read_instance.active_dashboard_plots:
            default = default_values['active_dashboard_plots']
            self.read_instance.active_dashboard_plots = default
        # TODO: For Taylor diagrams, remove this piece of code when we stop using Matplotlib 3.3
        else:
            if Version(matplotlib.__version__) < Version("3.8"):
                if 'taylor' in self.read_instance.active_dashboard_plots:
                    error = 'It is not possible to create Taylor diagrams yet, please remove from settings/report_plots.yaml.'
                    self.read_instance.logger.error(error)
                    sys.exit(1)

        if (len(self.read_instance.active_dashboard_plots) != 4) and (not self.read_instance.offline) & (not self.read_instance.interactive):
            error = 'Error: there must be 4 "active_dashboard_plots"'
            self.read_instance.logger.error(error)
            sys.exit(1)
        
        # if filter_species is active, and spatial_colocation is not active, then cannot filter by species
        # set filter_species to empty dict and advise user of this
        if (self.read_instance.filter_species) and (not self.read_instance.spatial_colocation):
            self.read_instance.filter_species = {}
            msg = 'Spatial colocation (spatial_colocation) must be set to True if wanting to filter by species.'
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)

        # map to multiple species if have * wildcard
        # also duplicate out associated network
        # remove any species for which there exists no data
        new_species = copy.deepcopy(self.read_instance.species)
        for speci_ii, speci in enumerate(self.read_instance.species): 
            if '*' in speci:
                # throw mapping error if species not ablet to be mapped
                if speci not in multispecies_map:
                    error = f'Error: not able to map species "{speci}".'
                    self.read_instance.logger.error(error)
                    sys.exit(1)
                
                mapped_species = multispecies_map[speci]
                del new_species[speci_ii]
                new_species[speci_ii:speci_ii] = mapped_species
                # in download mode is not necessary to duplicate the networks
                if not self.read_instance.download:
                    network_to_duplicate = self.read_instance.network[speci_ii]
                    del self.read_instance.network[speci_ii]
                    self.read_instance.network[speci_ii:speci_ii] = [network_to_duplicate]*len(mapped_species)
        self.read_instance.species = copy.deepcopy(new_species)

        # get species and filter species which are not on the current ghost version
        invalid_species = set(self.read_instance.species) - set(self.read_instance.available_species)
        invalid_filter_species = set(map(lambda x:x.split('|')[1], self.read_instance.filter_species)) - set(self.read_instance.available_species)                          
        
        # check species, remove the ones that are not on the ghost version       
        if invalid_species:                                                            
            msg = f'Removing invalid species {", ".join(invalid_species)} for the current GHOST version ({self.read_instance.ghost_version})'
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
            for inv_species in invalid_species:
                self.read_instance.species.remove(inv_species)
            # exit if there are no valid species left
            if not self.read_instance.species:
                error = f"Error: No valid species for the current GHOST version ({self.read_instance.ghost_version})"
                self.read_instance.logger.error(error)
                sys.exit(1)

        # check filter species, remove the ones that are not on the ghost version     
        if invalid_filter_species:
            msg = f'Removing invalid filter species {", ".join(invalid_filter_species)} for the current GHOST version ({self.read_instance.ghost_version})'
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
            # remove them from the filter_species atribute
            for filter_species in self.read_instance.filter_species.keys():
                if filter_species.split('|')[1] in invalid_filter_species:
                    del self.read_instance.filter_species[filter_species]

        # TODO change this in the refactoring, not do this in download mode
        if not self.read_instance.download: 
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
                        self.read_instance.logger.error(error)
                        sys.exit(1)
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
                        self.read_instance.logger.error(error)
                        sys.exit(1)
                    else:
                        for speci_ii, speci in enumerate(self.read_instance.species):
                            upper_bound_dict[speci] = self.read_instance.upper_bound[speci_ii] 
                        # add filter_species (using GHOST limits)
                        for speci in filter_species:
                            upper_bound_dict[speci] = np.float32(self.read_instance.parameter_dictionary[speci]['extreme_upper_limit'])
                self.read_instance.upper_bound = upper_bound_dict

            # create a variable to set qa per species (including filter species), setting defaults in the process
            if isinstance(self.read_instance.qa, dict):
                self.read_instance.qa_per_species = {speci:get_default_qa(self.read_instance, speci) 
                                                    for speci in species_plus_filter_species}
                # set qa to be first of qa per species pairs
                self.read_instance.qa = self.read_instance.qa_per_species[list(self.read_instance.qa_per_species.keys())[0]]
            else:
                self.read_instance.qa_per_species = {speci:self.read_instance.qa for speci in species_plus_filter_species}

        # add to qa
        if self.read_instance.add_qa:
            for qa_flag_to_add in self.read_instance.add_qa:
                if qa_flag_to_add not in self.read_instance.qa:
                    self.read_instance.qa.append(qa_flag_to_add)
                for speci in self.read_instance.qa_per_species:
                    if qa_flag_to_add not in self.read_instance.qa_per_species[speci]:
                        self.read_instance.qa_per_species[speci].append(qa_flag_to_add)

            self.read_instance.qa = sorted(self.read_instance.qa) 
            for speci in self.read_instance.qa_per_species:
                self.read_instance.qa_per_species[speci] = sorted(self.read_instance.qa_per_species[speci])

        # subtract from qa
        if self.read_instance.subtract_qa:
            for qa_flag_to_remove in self.read_instance.subtract_qa:
                if qa_flag_to_remove in self.read_instance.qa:
                    self.read_instance.qa.remove(qa_flag_to_remove)
                for speci in self.read_instance.qa_per_species:
                    if qa_flag_to_remove in self.read_instance.qa_per_species[speci]:
                        self.read_instance.qa_per_species[speci].remove(qa_flag_to_remove)

            self.read_instance.qa = sorted(self.read_instance.qa) 
            for speci in self.read_instance.qa_per_species:
                self.read_instance.qa_per_species[speci] = sorted(self.read_instance.qa_per_species[speci])

        # add to flags
        if self.read_instance.add_flags:
            for flag_to_add in self.read_instance.add_flags:
                if flag_to_add not in self.read_instance.flags:
                    self.read_instance.flags.append(flag_to_add)

            self.read_instance.flags = sorted(self.read_instance.flags) 

        # subtract from flags
        if self.read_instance.subtract_flags:
            for flag_to_remove in self.read_instance.subtract_flags:
                if flag_to_remove in self.read_instance.flags:
                    self.read_instance.flags.remove(flag_to_remove)

            self.read_instance.flags = sorted(self.read_instance.flags) 

        # if are using dashboard then just take first network/species pair, as multivar not supported yet
        if ((len(self.read_instance.network) > 1) and (len(self.read_instance.species) > 1) and 
            (not self.read_instance.offline) and (not self.read_instance.interactive) and (not self.read_instance.download) and (not self.read_instance.interpolation)):
             
            msg = 'Multiple networks/species are not supported in the dashboard. First ones will be taken.'
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)

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
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
                    elif (':' not in lower_limit) and ('>' not in lower_limit):
                        msg = 'Lower bound ({}) for {} should contain > or >=. '.format(lower_limit, networkspeci)
                        lower_limit = '>=' + lower_limit
                        msg += 'Setting it to be {}.'.format(lower_limit)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)

                    # modify upper bound to be :, or contain < or <=
                    if ('>' in upper_limit):
                        msg = 'Upper bound ({}) for {} cannot contain > or >=. '.format(upper_limit, networkspeci)
                        upper_limit = '<=' + upper_limit.replace('>', '').replace('=', '')
                        msg += 'Setting it to be {}.'.format(upper_limit)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
                    elif (':' not in upper_limit) and ('<' not in upper_limit):
                        msg = 'Upper bound ({}) for {} should contain < or <=. '.format(upper_limit, networkspeci)
                        upper_limit = '<=' + upper_limit
                        msg += 'Setting it to be {}.'.format(upper_limit)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
                    
                    # update symbols next to values
                    self.read_instance.filter_species[networkspeci][networkspeci_limit_ii] = [lower_limit, upper_limit, 
                                                                                              filter_species_fill_value]

        if (MACHINE == 'local') and (not self.read_instance.interp_multiprocessing) and (self.read_instance.interpolation):
            msg = 'During interpolation multiprocessing must be turned on for local runs, activating.'
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf, deactivate=deactivate_warning)
            self.read_instance.interp_multiprocessing = True
    
    def switch_logging(self):
        # create logger
        self.read_instance.logger = logging.getLogger("")
        self.read_instance.logger.setLevel(logging.INFO) 

        # remove previous handlers
        while self.read_instance.logger.handlers:
            self.read_instance.logger.removeHandler(self.read_instance.logger.handlers[0])

        # interpolation does not use this feature
        if self.read_instance.logfile != False and self.read_instance.interpolation is False:
            # default path, default name
            if self.read_instance.logfile == True:
                # get log filename and filepath
                timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
                filename = f"{timestamp}.log"
            # default path, custom name
            else:
                filename = str(self.read_instance.logfile)

            # custom path, custom name (no need of creating the file path)
            if type(self.read_instance.logfile) == str and os.sep in self.read_instance.logfile:
                file_path = self.read_instance.logfile
            # default path (the file name depends on the mode)
            else:
                # get the mode being used right now
                mode_list = ["offline", "interactive", "download", "interpolation"] 
                mode = "dashboard"
                for temp_mode in mode_list:
                    if getattr(self.read_instance, temp_mode) is True:
                        mode = temp_mode if temp_mode != "interactive" else "notebook"
                file_path = join(PROVIDENTIA_ROOT, 'logs', mode, filename)

            # redirect output to a file
            handler = logging.FileHandler(file_path)
            print(f"Output redirected to {file_path}")
        else:
            # redirect output to terminal
            handler = logging.StreamHandler(sys.stdout)

        formatter = logging.Formatter("%(message)s")
        handler.setFormatter(formatter)
        self.read_instance.logger.addHandler(handler)

        # suppress paramiko logs in download
        if self.read_instance.download is True:
            logging.getLogger("paramiko").setLevel(logging.WARNING)

def read_conf(self, fpath=None):
    """ Read configuration files. """

    res = {}
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
                        self.read_instance.logger.error(error)
                        sys.exit(1)
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
                # allow # after first character to partially comment lines
                line_strip = line.strip().split('#')[0]
                
                # get current section                        
                if '[' in line and ']' in line and '[[' not in line and ']]' not in line:
                    current_section = line_strip.replace('[', '').replace(']', '')
                
                # get current subsection
                if ('[[' in line_strip) and (']]' in line_strip):
                    current_modified_subsection = current_section + "·" + line_strip.replace('[[', '').replace(']]', '')
                
                # parsing all but last section 
                if section_modified != all_sections_modified[-1]:
                    # start of relevant section 
                    if line_strip == all_sections[i]:
                        # if subsection, make sure its section corresponds to current section to avoid
                        # problems with repeated subsection names (SECTIONA·SUBSECTION, SECTIONB·SUBSECTION)
                        if ('[[' in line_strip) and (']]' in line_strip):
                            if (current_modified_subsection != section_modified):
                                copy = False
                            else:
                                copy = True
                        elif line_strip in repeated_subsections:
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
                        # if subsection, make sure its section corresponds to current section to avoid
                        # problems with repeated subsection names (SECTIONA·SUBSECTION, SECTIONB·SUBSECTION)
                        if ('[[' in line_strip) and (']]' in line_strip):
                            if (current_modified_subsection != section_modified):
                                copy = False
                            else:
                                copy = True
                        elif line_strip in repeated_subsections:
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
                    # if lines are not empty and not commented
                    if (line_strip != '') and ('#' not in line_strip):
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
        
        # reset res variable
        res_sub = {}

        # determine if subsection or not
        if '·' in section_modified:
            is_subsection = True
            par_section = section_modified.split('·')[0]
            # add attributes from parent section
            for par_k, par_val in config[par_section].items():
                # if first character of key is comment character (#), do not parse this attribute
                if par_k[0] == '#':
                    continue
                # transform str booleans, ints, floats etc. to their real type if possible
                try:
                    par_val = ast.literal_eval(par_val)
                except:
                    pass
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
            # transform str booleans, ints, floats etc. to their real type if possible
            try:
                val = ast.literal_eval(val)
            except:
                pass
            res_sub[k] = val

        # store pairs into res variable
        res[section_modified] = res_sub

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
    sys.path.append(join(PROVIDENTIA_ROOT, 'providentia'))
    from configuration import read_conf

    if fpath is None:
        self.read_instance.logger.info("No configuration file found")
        sys.exit(1)

    # if DEFAULT is not present, then return
    if not os.path.isfile(fpath):
        self.read_instance.logger.error(f"Error {fpath}")
        return

    self.sub_opts, self.all_sections, self.parent_section_names, self.subsection_names, self.filenames = read_conf(self, fpath)

def split_options(read_instance, conf_string, separator="||"):
    """ For the options in the configuration that define the keep and remove
        options. Returns the values in two lists, the keeps and removes.
    """
    # import show_mesage from warnings
    from .warnings_prv import show_message
    
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
