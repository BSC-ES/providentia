import copy
import glob
import inspect
import os
import json
import random
import socket
import subprocess
import sys
import time
from pydoc import locate
import numpy as np

MACHINE = os.environ.get('BSC_MACHINE', '')

# get current path and providentia root path
CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = os.path.dirname(os.path.dirname(CURRENT_PATH))

sys.path.append(os.path.join(PROVIDENTIA_ROOT, 'providentia', 'interpolation'))
sys.path.append(os.path.join(PROVIDENTIA_ROOT, 'providentia'))

from aux import get_aeronet_bin_radius_from_bin_variable, check_for_ghost
from mapping_species import mapping_species
from configuration import ProvConfiguration, load_conf

# load the data_paths for the different machines and the default values jsons
data_paths = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings', 'data_paths.json')))
config_format = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings', 'prov_interp_defaults.json')))
bin_vars = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings', 'multispecies_shortcurts.json')))

# load the defined experiments paths and agrupations jsons
experiment_paths = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings', 'experiment_paths.json')))
experiment_names = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings', 'experiment_names.json')))

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

class SubmitInterpolation(object):
    """ Class that handles the interpolation submission. """

    def __init__(self):

        # start timer
        self.start = time.time()

        # get current BSC machine
        self.machine = MACHINE

        # define current working directory and
        # arguments/submit/interpolation log subdirectories
        self.working_directory = os.getcwd()
        if self.machine == "nord3":
            self.working_directory = self.working_directory[17:]
        self.arguments_dir = '{}/arguments'.format(self.working_directory)
        self.submit_dir = '{}/submit'.format(self.working_directory)
        self.interpolation_log_dir = '{}/interpolation_logs'.format(self.working_directory)

        # set SLURM job ID as unique ID for tracking tasks
        # defined to process in the configuration file
        self.unique_ID = sys.argv[1]

        json_kwargs = "".join(sys.argv[2:]).replace(":", '":"').replace(",", '","').replace("{", '{"').replace("}", '"}')
        kwargs = json.loads(json_kwargs)

        # initialize commandline arguments, if given
        provconf = ProvConfiguration(self, **kwargs)

        # update variables from config file
        if self.config != '':  
            read_conf = False
            if os.path.exists(self.config):
                read_conf = True
            else: 
                if os.path.exists(os.path.join(self.config_dir, self.config)):
                    self.config = os.path.join(self.config_dir, self.config)
                    read_conf = True
            if read_conf:
                load_conf(self, self.config)
                self.from_conf = True
            else:
                error = 'Error: The path to the configuration file specified in the command line does not exist.'
                sys.exit(error)
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            sys.exit(error)

        # get main section args
        self.config_dict = self.sub_opts[self.parent_section_names[0]]
        self.config_args = self.config_dict.keys()

        # dictionary that stores utilized interpolation variables
        self.interpolation_variables = {}

        # put configuration variables into self, assigning defaults where neccessary 
        self.set_configuration_defaults(vars_to_set=['ghost_version'])
        self.set_configuration_defaults(vars_not_to_set=['ghost_version'])

        # print variables used        
        print("\nVariables used for the interpolation:")
        for arg,value in self.interpolation_variables.items():
            print(f"{arg}: {value}")

        # define the QOS (Quality of Service) used to manage jobs on the SLURM system
        if self.machine == 'mn5':
            self.qos = 'gp_bsces'
        else:
            self.qos = 'bsc_es'

        # import unit converter
        sys.path.append(os.path.join(PROVIDENTIA_ROOT, 'providentia', 'dependencies','unit-converter'))
        import unit_converter

    def set_configuration_defaults(self, vars_to_set=None, vars_not_to_set=None):
        ''' Fill relevant fields appropriately with defaults, where 'default' is set in configuration file.
            Set some key variables for job submission and check formatting of configuration file variables are correct, 
            if not throw errors.
            
            :param vars_to_set: list with configuration variables to be set
            :type binvar: list
            :param vars_not_to_set: list with configuration variables not to be set
            :type vars_not_to_set: list
        '''
  
        # if have defined variables to set, do just that
        # otherwise set all variables in config file
        if not vars_to_set:
            vars_to_set = list(config_format.keys())
        
        # remove variables that do not want to set
        if vars_not_to_set:
            for var_not_to_set in vars_not_to_set:
                vars_to_set.remove(var_not_to_set)
            
        # check each variable to set

        for var_to_set in vars_to_set:

            # get format information for variable
            var_config_format = list(config_format[var_to_set].keys())

            # if variable not in configuration file, then set default value for field if available
            # otherwise throw error
            if not var_to_set in self.config_args: 
                if 'default' not in var_config_format:
                    sys.exit('CONFIGURATION FILE DOES NOT CONTAIN REQUIRED ARGUMENT: {}'.format(var_to_set))
                else:
                    if var_to_set == 'species_to_process':
                        self.config_dict[var_to_set] =  [self.standard_parameters[param]['bsc_parameter_name'] 
                                                          for param in self.standard_parameters.keys()] 
                    else:
                        self.config_dict[var_to_set] = config_format[var_to_set]['default']
                        
            # otherwise if variable in configuration file, make sure formatting is correct
            else:
                # set default for config argument (if required)
                if 'default' in var_config_format:
                    if self.config_dict[var_to_set] == 'default':
                        if var_to_set == 'species_to_process':
                            self.config_dict[var_to_set] = [self.standard_parameters[param]['bsc_parameter_name'] 
                                                              for param in self.standard_parameters.keys()]
                        else:
                            self.config_dict[var_to_set] =  config_format[var_to_set]['default']
                    # transform ensemble_options to string when ensemble id is passed
                    elif var_to_set == 'ensemble_options' and type(self.config_dict[var_to_set]) == int:
                        ensembleid = str(self.config_dict[var_to_set])
                        self.config_dict[var_to_set] = '0' * (3-len(ensembleid)) + ensembleid
                
                # transform non-default atributes which are still in str format to list
                if config_format[var_to_set]['type'] == 'list' and locate(config_format[var_to_set]['type']) != type(self.config_dict[var_to_set]):
                    if type(config_format[var_to_set]['type']) == str:
                        self.config_dict[var_to_set] = self.config_dict[var_to_set].split("(")[0].replace(" ", "").split(",")
                    if config_format[var_to_set]['subtype'] == 'int': # it is an int
                        self.config_dict[var_to_set] = [int(var) for var in self.config_dict[var_to_set]]

                # check primary typing is correct
                if locate(config_format[var_to_set]['type']) != type(self.config_dict[var_to_set]):
                    sys.exit('CONFIGURATION FILE ARGUMENT: {} NEEDS TO BE A {} TYPE'.format(
                        var_to_set, config_format[var_to_set]['type']))

                # check subtyping is correct
                if 'subtype' in var_config_format:
                    for var in self.config_dict[var_to_set]:
                        if locate(config_format[var_to_set]['subtype']) != type(var):
                            sys.exit('CONFIGURATION FILE ARGUMENT: {} NEEDS TO BE A LIST CONTAINING {} TYPES'.format(
                                var_to_set, config_format[var_to_set]['subtype']))

            # set some extra variables for ghost_version
            if var_to_set == 'ghost_version':
                # import GHOST standards
                sys.path.insert(1, data_paths[self.machine]["ghost_root"] + '/GHOST_standards/{}'.format(
                    self.config_dict['ghost_version']))
                from GHOST_standards import standard_parameters
                self.standard_parameters = standard_parameters

            # fill in bin wildcard parameters
            if var_to_set == 'species_to_process':
                if self.config_dict[var_to_set] != 'default':
                    for arg_var in self.config_dict[var_to_set]:
                        if arg_var in list(bin_vars.keys()):
                            self.config_dict[var_to_set] = bin_vars[arg_var]

            # add config argument to self
            setattr(self, var_to_set, self.config_dict[var_to_set])

            # add config argument to used interpolation variable list
            self.interpolation_variables[var_to_set] = self.config_dict[var_to_set]
 
    def gather_arguments(self):
        ''' Gather list of arguments for all unique tasks to process, as defined in the configuration file. '''

        # create arguments list
        self.arguments = []

        # create list to hold root directory locations of .out log files generated for each task processed
        self.output_log_roots = []

        # iterate through desired experiment IDs
        for experiment_to_process in self.experiments:
            
            print('\nEXPERIMENT: {0}'.format(experiment_to_process))

            # create arguments list
            exp_arguments = []

            # initialise files lists
            obs_files = []
            exp_files = []
            obs_files_dates = []
            exp_files_dates = []

            # get experiment type and specific directory 
            experiment_exists = False
            for self.experiment_type in experiment_names:
                if experiment_to_process in experiment_names[self.experiment_type]:
                    exp_dict = experiment_paths[self.experiment_type]
                    experiment_exists = True
                    break
            
            # if experiment id is not defined, exit
            if not experiment_exists:
                msg = f"THE EXPERIMENT ID '{experiment_to_process}' DOESN'T EXIST"
                sys.exit(msg)

            # take gpfs directory preferentially over esarchive directory
            if 'gpfs' in list(exp_dict.keys()):
                exp_dir = exp_dict['gpfs']
            else:
                exp_dir = exp_dict['esarchive']
            if ('gpfs' in exp_dir) and (self.machine == 'mn5'):
                exp_dir = exp_dir.replace('/gpfs/', '/gpfs/tapes/MN4/')

            # add file to directory path
            exp_dir += f"{experiment_to_process}/" 

            # get all grid type subdirectories for current experiment
            available_grid_types = [name for name in os.listdir(exp_dir) if os.path.isdir("{}/{}".format(
                exp_dir,name))]
            
            # get intersection between desired grid types to process and grid types available in directory
            available_grid_types = [x for x in available_grid_types if x in self.domain]

            # iterate through grid types to process
            for grid_type in available_grid_types:

                # get all temporal resolution subdirectories for current experiment/grid_type
                available_model_temporal_resolutions = [name for name in os.listdir("{}/{}/".format(exp_dir, grid_type)) 
                                                        if os.path.isdir("{}/{}/{}".format(exp_dir, grid_type, name))]

                # iterate through species to process
                for speci_ii, speci_to_process in enumerate(self.species):

                    original_speci_to_process = copy.deepcopy(speci_to_process)

                    # iterate through temporal_resolutions to output
                    for temporal_resolution_to_output in self.resolution:

                        experiment_species_ensemblestat = []

                        # only keep available temporal resolutions which are equal or finer in resolution to that 
                        # wanted to output
                        # in case temporal_resolution_to_output is instantaneous, but this not available in model, 
                        # attempt to take non-instantaneous resolution 
                        # in case temporal_resolution_to_output is non-instantaneous, but this not available in model, 
                        # attempt to take instantaneous resolution

                        if temporal_resolution_to_output == 'hourly':
                            resolutions_to_keep = ['hourly', 'hourly_instantaneous']
                        elif temporal_resolution_to_output == 'hourly_instantaneous':
                            resolutions_to_keep = ['hourly_instantaneous', 'hourly']
                        elif temporal_resolution_to_output == '3hourly':
                            resolutions_to_keep = ['3hourly','hourly', '3hourly_instantaneous', 'hourly_instantaneous']
                        elif temporal_resolution_to_output == '3hourly_instantaneous':
                            resolutions_to_keep = ['3hourly_instantaneous', 'hourly_instantaneous', '3hourly', 'hourly']
                        elif temporal_resolution_to_output == '6hourly':  
                            resolutions_to_keep = ['6hourly', '3hourly', 'hourly', '6hourly_instantaneous',
                                                '3hourly_instantaneous', 'hourly_instantaneous']
                        elif temporal_resolution_to_output == '6hourly_instantaneous':
                            resolutions_to_keep = ['6hourly_instantaneous', '3hourly_instantaneous', 'hourly_instantaneous',
                                                '6hourly', '3hourly', 'hourly']
                        elif temporal_resolution_to_output == 'daily':
                            resolutions_to_keep = ['daily', '6hourly', '3hourly', 'hourly', '6hourly_instantaneous',
                                                '3hourly_instantaneous', 'hourly_instantaneous']
                        elif temporal_resolution_to_output == 'monthly':
                            resolutions_to_keep = ['monthly', 'daily', '6hourly', '3hourly', 'hourly', 
                                                '6hourly_instantaneous', '3hourly_instantaneous', 'hourly_instantaneous']
                        
                        # iterate through resolutions_to_keep until find one for which have speci_to_process (or mapped speci)
                        have_valid_resolution = False
                        for model_temporal_resolution in resolutions_to_keep:
                            # test if have directory for current speci_to_process
                            if os.path.isdir("{}/{}/{}/{}".format(exp_dir, grid_type, model_temporal_resolution, speci_to_process)):
                                have_valid_resolution = True
                                break

                            # test if have speci directory in ensemble-stats
                            elif os.path.isdir("{}/{}/{}/ensemble-stats".format(exp_dir, grid_type, model_temporal_resolution)):

                                # get all ensemble-stats species
                                experiment_species_ensemblestat = list(np.unique([name.split('_')[0] 
                                                                    for name in os.listdir("{}/{}/{}/ensemble-stats".format(
                                                                      exp_dir,grid_type,model_temporal_resolution)) 
                                                                        if os.path.isdir("{}/{}/{}/ensemble-stats/{}".format(
                                                                          exp_dir,grid_type,model_temporal_resolution,name))]))
                                
                                # test if have speci_to_process in experiment_species_ensemblestat
                                if speci_to_process in experiment_species_ensemblestat:
                                    have_valid_resolution = True
                                    break
                                # test if have mapped speci_to_process in experiment_species_ensemblestat
                                elif speci_to_process in mapping_species:
                                    for speci_to_map in mapping_species[speci_to_process]:
                                        if speci_to_map in experiment_species_ensemblestat:
                                            # if have a binned size distribution variable to map, first check if bin radius is within model's bin extents
                                            # if not, do not process species
                                            if 'vconcaerobin' in speci_to_map:
                                                # model is MONARCH?
                                                if self.experiment_type == 'monarch':
                                                    # check if bin radius is within MONARCH's bin extents (0.2-20.0 um)
                                                    bin_radius = get_aeronet_bin_radius_from_bin_variable(speci_to_map)
                                                    if (bin_radius >= 0.2) & (bin_radius <= 20.0):
                                                        speci_to_process = copy.deepcopy(speci_to_map)
                                                        have_valid_resolution = True
                                                        break
                                            else:
                                                speci_to_process = copy.deepcopy(speci_to_map)
                                                have_valid_resolution = True
                                                break
                                    
                                    if have_valid_resolution:
                                        break

                            # for some variables it is possible to extract the variable information by mapping to a different variable name, with a higher dimensionality
                            # this currently is implemented for 2 cases:
                            # -- 4D binned size distribution
                            # -- 4D gas variables
                            # first check if speci_to_process can be mapped
                            # then check if the variable to map to exists for the experiment/grid_type/resolution (these can be multiple, list order sets the priority)
                            else:
                                # get species that can do mapping for
                                available_species_to_map_from = list(mapping_species.keys())
                                # check if speci_to_process can be mapped
                                if speci_to_process in available_species_to_map_from:

                                    # if it can be then check then if the variable to map to exists for the experiment/grid_type/resolution 
                                    # (these can be multiple, list order sets the priority)
                                    for speci_to_map in mapping_species[speci_to_process]:
                                        if os.path.isdir("{}/{}/{}/{}".format(exp_dir, grid_type, model_temporal_resolution, speci_to_map)):

                                            # if have a binned size distribution variable to map, first check if bin radius is within model's bin extents
                                            # if not, do not process species
                                            if 'vconcaerobin' in speci_to_map:
                                                # model is MONARCH?
                                                if self.experiment_type == 'monarch':
                                                    # check if bin radius is within MONARCH's bin extents (0.2-20.0 um)
                                                    bin_radius = get_aeronet_bin_radius_from_bin_variable(speci_to_map)
                                                    if (bin_radius >= 0.2) & (bin_radius <= 20.0):
                                                        speci_to_process = copy.deepcopy(speci_to_map)
                                                        have_valid_resolution = True
                                                        break
                                            else:
                                                speci_to_process = copy.deepcopy(speci_to_map)
                                                have_valid_resolution = True
                                                break
                                
                                    if have_valid_resolution:
                                        break

                        # only proceed if have valid resolution
                        if have_valid_resolution:

                            # iterate through observational networks to interpolate against
                            for network_to_interpolate_against in self.network:

                                # define if network is in GHOST format
                                self.reading_ghost = check_for_ghost(network_to_interpolate_against)

                                # get all relevant observational files
                                # GHOST
                                if self.reading_ghost:
                                    obs_files = np.sort(glob.glob(
                                        data_paths[self.machine]['ghost_root'] + '/{}/{}/{}/{}/{}*.nc'.format(
                                            network_to_interpolate_against, self.ghost_version, 
                                            temporal_resolution_to_output, original_speci_to_process, 
                                            original_speci_to_process)))
                                # non-GHOST
                                else:
                                    obs_files = np.sort(glob.glob(
                                        data_paths[self.machine]['nonghost_root'] + '/{}/{}/{}/{}*.nc'.format(
                                            network_to_interpolate_against, temporal_resolution_to_output,
                                            original_speci_to_process, original_speci_to_process)))

                                # if have no observational files then continue
                                if len(obs_files) == 0:
                                    continue

                                # iterate through ensemble options
                                for ensemble_option in self.ensemble_options:
                                
                                    # check if ensemble option is ensemble stat and get all relevant experiment files
                                    if 'stat_' in ensemble_option:
                                        ensemble_stat = True
                                        exp_files_all = np.sort(glob.glob('{}/{}/{}/ensemble-stats/{}_{}/{}*{}.nc'.format(
                                            exp_dir, grid_type, model_temporal_resolution, speci_to_process, 
                                            ensemble_option[5:], speci_to_process, ensemble_option[5:])))
                                    else:
                                        ensemble_stat = False
                                        exp_files_all = np.sort(glob.glob('{}/{}/{}/{}/{}*.nc'.format(
                                            exp_dir, grid_type, model_temporal_resolution, speci_to_process, 
                                            speci_to_process)))    
                                        
                                        # drop all analysis files ending with '_an.nc' which are not in ensemble-stats
                                        exp_files_all = np.array([f for f in exp_files_all if '_an.nc' not in f])

                                    # if have no relevant experiment files then continue
                                    if len(exp_files_all) == 0:
                                        continue

                                    # ensemble stat?
                                    if ensemble_stat:
                                        available_ensemble_options = [ensemble_option[5:]]    

                                    # not ensemble stat?
                                    else:                                
                                        
                                        # determine if simulation generates files with ensemble member numbers or not 
                                        # (test first file)
                                        if exp_files_all[0].split('/')[-1].rsplit('_', 1)[0] != speci_to_process:
                                            have_ensemble_members = True
                                            # if have ensemble members in filename, get all unique numbers
                                            unique_ensemble_members = np.unique([f.split('/{}-'.format(speci_to_process))[-1][:3] 
                                                                                for f in exp_files_all])
                                            # get intersection between desired ensemble members to process and those 
                                            # available in directory
                                            # if no members defined explicitly to process, process them all 
                                            if ensemble_option == 'allmembers':
                                                available_ensemble_options = unique_ensemble_members
                                            else:
                                                if ensemble_option in unique_ensemble_members:
                                                    available_ensemble_options = [ensemble_option]
                                                else:
                                                    continue
                                        else:
                                            have_ensemble_members = False
                                            # if have defined ensemble members to process, then continue as no files in this 
                                            # directory have ensemble member number
                                            if ensemble_option != 'allmembers':
                                                continue
                                            # otherwise, proceed (tag files as ensemble member '000' for sake of operation)
                                            else:
                                                available_ensemble_options = ['000']
                                    
                                    # iterate through available ensemble options to process
                                    for available_ensemble_option in available_ensemble_options:
                                
                                        # limit experiment files to be just those for specific ensemble member 
                                        # (where neccessary) 
                                        exp_file_speci = copy.deepcopy(speci_to_process)
                                        exp_files = copy.deepcopy(exp_files_all)
                                        if ensemble_stat == False:
                                            if have_ensemble_members == True:
                                                exp_file_speci = '{}-{}'.format(speci_to_process, available_ensemble_option)
                                                exp_files = np.sort([f for f in exp_files_all if '{}_'.format(exp_file_speci) 
                                                                    in f])
                                        
                                        # get all observational file start dates (year and month)
                                        obs_files_dates = []
                                        for obs_file in obs_files:
                                            obs_file_date = obs_file.split('{}_'.format(original_speci_to_process))[-1].split('_')[0].split('.nc')[0]
                                            obs_files_dates=np.append(obs_files_dates,obs_file_date[:6])

                                        # get all experiment file start dates (year and month)
                                        exp_files_dates = []   
                                        for exp_file in exp_files:
                                            exp_file_date = exp_file.split('{}_'.format(exp_file_speci))[-1].split('_')[0].split('.nc')[0]
                                            exp_files_dates=np.append(exp_files_dates,exp_file_date[:6])

                                        # remove observational files outside date ranges 
                                        obs_files_ii = np.array([obs_files_ii for obs_files_ii, obs_files_date in enumerate(obs_files_dates) 
                                                                if ((int(obs_files_date) >= int(self.start_date)) 
                                                                    and (int(obs_files_date) < int(self.end_date)))])       
                                        if len(obs_files_ii) == 0:
                                            obs_files = []
                                            obs_files_dates = []
                                            continue
                                        obs_files = obs_files[obs_files_ii]
                                        obs_files_dates = obs_files_dates[obs_files_ii]

                                        # remove experiment files outside date ranges 
                                        exp_files_ii = np.array([exp_files_ii for exp_files_ii, exp_files_date in enumerate(exp_files_dates) 
                                                                if ((int(exp_files_date) >= int(self.start_date)) 
                                                                    and (int(exp_files_date) < int(self.end_date)))])      
                                        if len(exp_files_ii) == 0:
                                            exp_files = []
                                            exp_files_dates = []
                                            continue
                                        exp_files = exp_files[exp_files_ii]
                                        exp_files_dates = exp_files_dates[exp_files_ii]
                                        
                                        # get intersection of file yearmonths between observations and experiment
                                        intersect_yearmonths = np.intersect1d(obs_files_dates, exp_files_dates)

                                        # if have no intersecting months, continue
                                        if len(intersect_yearmonths) == 0:
                                            continue

                                        # create Providentia experiment code (expid-region-ensembleoption)
                                        prov_exp_code = '{}-{}-{}'.format(experiment_to_process, grid_type, 
                                                                        available_ensemble_option)
        
                                        # create directories to store slurm output/error logs for interpolation task of 
                                        # specific combination of iterated variables (if does not already exist)
                                        if not os.path.exists('{}/{}/{}/{}/{}'.format(self.interpolation_log_dir, 
                                                                                      prov_exp_code, 
                                                                                      original_speci_to_process, 
                                                                                      network_to_interpolate_against, 
                                                                                      temporal_resolution_to_output)):
                                            os.makedirs('{}/{}/{}/{}/{}'.format(self.interpolation_log_dir, prov_exp_code, 
                                                                                original_speci_to_process, 
                                                                                network_to_interpolate_against, 
                                                                                temporal_resolution_to_output))

                                        # iterate through intersecting yearmonths and write all current variable arguments 
                                        # to arguments file
                                        for yearmonth in intersect_yearmonths:

                                            # append current iterative arguments to arguments list               
                                            self.arguments.append("{} {} {} {} {} {} {}".format(prov_exp_code, 
                                                                                                model_temporal_resolution, 
                                                                                                speci_to_process, 
                                                                                                network_to_interpolate_against, 
                                                                                                temporal_resolution_to_output, 
                                                                                                yearmonth, 
                                                                                                original_speci_to_process))

                                            # append root name of .out file that will be output for each processed task
                                            self.output_log_roots.append('{}/{}/{}/{}/{}/{}'.format(self.interpolation_log_dir, 
                                                                                                    prov_exp_code, 
                                                                                                    original_speci_to_process, 
                                                                                                    network_to_interpolate_against, 
                                                                                                    temporal_resolution_to_output, 
                                                                                                    yearmonth))

                                            # remove previous output logs
                                            previous_logs = glob.glob('{}/{}/{}/{}/{}/{}*'.format(self.interpolation_log_dir, 
                                                                                                  prov_exp_code, 
                                                                                                  original_speci_to_process, 
                                                                                                  network_to_interpolate_against, 
                                                                                                  temporal_resolution_to_output, 
                                                                                                  yearmonth))
                                            for previous_log in previous_logs:
                                                os.remove(previous_log) 


            # if have no arguments after iteration, return message stating that
            exp_arguments = copy.deepcopy(self.arguments)
            if len(exp_arguments) == 0:
                msg = 'NO INTERSECTING OBSERVATIONS/EXPERIMENT DATA FOR INTERPOLATION. \n' 
                if self.start_date == self.end_date:
                    msg += 'If you want to interpolate data for one month, '
                    msg += 'you need to set the end date to be the next one. \n'
                    msg += 'e.g. For November 2018, this is 201811 to 201812.'
                else:
                    msg += 'TRYING TO MATCH DATA FROM:\n'
                    if len(obs_files) == 0:
                        msg += f'Observations data between {self.start_date} and {self.end_date} cannot be found.'
                    else:
                        msg += f'Observation files: {obs_files}\n'
                        msg += f'Observation dates: {obs_files_dates}\n'
                    if len(exp_files) == 0:
                        msg += f'Experiment data between {self.start_date} and {self.end_date} cannot be found.'
                    else:
                        msg += f'Experiment files: {exp_files}\n'
                        msg += f'Experiment dates: {exp_files_dates}'
                print(msg)
                continue
        
        # if have no arguments for all experiments, return message stating that
        if len(self.arguments) == 0:
            if len(self.experiments) > 1:
                msg = 'INTERPOLATION CANNOT BE DONE FOR ANY EXPERIMENT'
                sys.exit(msg)
            else:
                sys.exit()

        # randomise the order of the arguments list
        random.shuffle(self.arguments)     


    def create_greasy_arguments_file(self):
        ''' Create greasy arguments text file storing all different tasks to run by greasy. '''

        # define list to store chunked argument files (to be submitted using greasy)
        argument_files = []

        # get the CPU chunk size -- set initially as miniumum number of arguments per file 
        N_arguments_per_file_minimum = copy.deepcopy(self.chunk_size)
        N_arguments_per_file = copy.deepcopy(self.chunk_size)
        
        # divide the number of arguments by the CPU chunk size, to determine how many argument files will be needed to 
        # submit all jobs
        N_submit_files = int(np.ceil(len(self.arguments)/self.chunk_size))
        
        # set argument remainder as 0 initially
        argument_remainder = 0

        # if the number of argument files is greater than the job array limit (i.e the limit on the number of argument 
        # files that can be processed simultaneously)
        # then add adjust minimum N arguments per file appropriately (i.e. split extra arguments across the maximum 
        # number of argument files evenly)        
        if N_submit_files > self.job_array_limit:
            # update the minimum number of arguments per file
            N_arguments_per_file_minimum = int(np.floor(len(self.arguments)/self.job_array_limit))
            N_arguments_per_file = copy.deepcopy(N_arguments_per_file_minimum)

            # if the number of extra arguments does not divide equally between all files get the remainder
            argument_remainder = int(len(self.arguments)%self.job_array_limit)
            
            # if have argument remainder then update N_arguments_per_file variable to be 1 greater than minimum for 
            # first file written (and for all files thereafter until  remainder is accounted for)
            if argument_remainder > 0:
                N_arguments_per_file = N_arguments_per_file_minimum + 1
                # subtract 1 from the argument remainder
                argument_remainder -= 1
            
            # set N submit files as N of job array limit
            N_submit_files = copy.deepcopy(self.job_array_limit)

        # create file which will store a list of all chunked argument filenames
        greasy_file = open('{}/{}.grz'.format(self.arguments_dir, self.unique_ID), 'w')
        
        # create all chunked argument filenames
        for ii in range(N_submit_files):
            argument_fname = '{}/{}_{}.txt'.format(self.arguments_dir,self.unique_ID,ii)
            argument_files.append(argument_fname)
            greasy_file.write('{}\n'.format(argument_fname))
        greasy_file.close()

        # initialise variables to count N argument files written and N lines written
        file_ii = 0
        n_lines_written = 0
        
        # define special characters that will need to be escaped in string written        
        special_characters = ['(',')']

        # open first arguments file for writing
        arguments_file = open(argument_files[file_ii], 'w')

        # iterate through different arguments, writing line by line to arguments file
        for arguments_ii, str_to_write in enumerate(self.arguments):

            # escape certain special characters in str_to_write
            for ch in special_characters:
                str_to_write = str_to_write.replace(ch,'\{}'.format(ch))
        
            # write arguments str to current file
            arguments_file.write('python -u {}/experiment_interpolation.py {}\n'.format(self.working_directory, 
                                                                                        str_to_write))

            # iterate n lines written    
            n_lines_written += 1
    
            # if have written sufficient arguments to file (and not currently processing last argument)
            # then close current file and open another
            if (n_lines_written == N_arguments_per_file) & (arguments_ii < (len(self.arguments)-1)):
                # reset n lines written
                n_lines_written = 0
                
                # close current arguments file
                arguments_file.close()
                
                # iterate n arguments files written
                file_ii += 1
                
                # if have argument remainder, write an extra argument to the next file
                if argument_remainder > 0:
                    N_arguments_per_file = N_arguments_per_file_minimum + 1
                    # subtract 1 from the argument remainder
                    argument_remainder -= 1
                # otherwise N arguments to write to next file should be minimum
                else:
                    N_arguments_per_file = copy.deepcopy(N_arguments_per_file_minimum)
                
                # open new arguments file
                arguments_file = open(argument_files[file_ii], 'w')
                
        # close current arguments file
        arguments_file.close()


    def create_slurm_submission_script(self):
        ''' Write a slurm submission shell script that submits a greasy job. '''

        # create job_fname (unique_ID + 'sh.')
        self.job_fname = self.unique_ID+'.sh'
    
        # get all argument files
        argument_files = sorted(glob.glob('{}/{}_*.txt'.format(self.arguments_dir,self.unique_ID)))

        # read how many lines are in first arguments file
        with open(argument_files[0]) as f: 
            for ii, line in enumerate(f):
                pass
            N_arguments = ii + 1

        # cap the number of simultaneously running tasks to be the defined CPU chunk size  
        max_tasks = copy.deepcopy(self.chunk_size)

        # if the number of arguments is > capped max tasks, then set N simultaneous tasks to be the max tasks permitted
        if N_arguments > max_tasks:
            n_simultaneous_tasks = copy.deepcopy(max_tasks) 
        # else, if the number of arguments is <= capped max tasks, then set N simultaneous tasks to be N arguments
        else:
            n_simultaneous_tasks = copy.deepcopy(N_arguments)

        # create slurm submission script
        submit_file = open(self.submit_dir+'/'+self.job_fname,"w")

        submit_file.write("#!/bin/bash\n")
        submit_file.write("\n")
        submit_file.write("#SBATCH --job-name=PRVI_{}\n".format(self.unique_ID))
        submit_file.write("#SBATCH --ntasks={}\n".format(n_simultaneous_tasks))
        # fix number of nodes to be 1 (for faster execution)
        submit_file.write("#SBATCH --nodes=1\n")
        submit_file.write("#SBATCH --time=02:00:00\n")
        submit_file.write("#SBATCH --array=1-{}\n".format(len(argument_files)))
        # if machine is power9, then there are 4 threads per physical CPU, which show as 4 separate CPUs. 
        # if multithreading is enabled, ensure 4 virtual cores are actually seen as 1 CPU.
        if (self.machine == 'power') & (self.multithreading == True):
            submit_file.write("#SBATCH --cpus-per-task=4\n")
        submit_file.write("#SBATCH --qos={}\n".format(self.qos))
        submit_file.write("##SBATCH --output=/dev/null\n")
        submit_file.write("##SBATCH --error=/dev/null\n")
        submit_file.write("\n")
        if self.machine == 'mn5':
            submit_file.write("#SBATCH --account=bsc32\n")  
            submit_file.write("#SBATCH --ntasks-per-node=1\n")
            submit_file.write("#SBATCH --cpus-per-task={}\n".format(n_simultaneous_tasks+1))
        else:
            submit_file.write("source {}/load_modules.sh\n".format(self.working_directory))
        submit_file.write("export GREASY_NWORKERS=$SLURM_NPROCS\n") 
        submit_file.write("export GREASY_LOGFILE={}/{}_$SLURM_ARRAY_TASK_ID.log\n".format(self.submit_dir, 
                                                                                          self.unique_ID))
        submit_file.write("export SLURM_CPU_BIND=none\n")
        submit_file.write("arguments_store={}/{}.grz\n".format(self.arguments_dir, self.unique_ID))
        submit_file.write("argument_file=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')\n")
        submit_file.write("\n")
        submit_file.write("greasy $argument_file")

        # close submit file
        submit_file.close()


    def create_lsf_submission_script(self):
        """ Write a lsf submission shell script that submits a greasy job. """

        # create job_fname (unique_ID + 'sh.')
        self.job_fname = self.unique_ID + '.sh'

        # get all argument files
        argument_files = sorted(glob.glob('{}/{}_*.txt'.format(self.arguments_dir, self.unique_ID)))

        # read how many lines are in first arguments file
        with open(argument_files[0]) as f:
            for ii, line in enumerate(f):
                pass
            N_arguments = ii + 1

        # cap the number of simultaneously running tasks to be the defined CPU chunk size
        max_tasks = copy.deepcopy(self.chunk_size)

        # if the number of arguments is > capped max tasks,
        # then set N simultaneous tasks to be the max tasks permitted
        if N_arguments > max_tasks:
            n_simultaneous_tasks = copy.deepcopy(max_tasks)
            # else, if the number of arguments is <= capped max tasks,
            # then set N simultaneous tasks to be N arguments
        else:
            n_simultaneous_tasks = copy.deepcopy(N_arguments)

        # create slurm submission script
        submit_file = open(self.submit_dir + '/' + self.job_fname, "w")

        submit_file.write("#!/bin/bash\n")
        submit_file.write("\n")
        submit_file.write("#BSUB -n {}\n".format(n_simultaneous_tasks))
        submit_file.write("#BSUB -W 01:00\n")
        submit_file.write("#BSUB -J PRVI_{}[1-{}]\n".format(self.unique_ID, len(argument_files)))
        submit_file.write("#BSUB -q {}\n".format(self.qos))
        submit_file.write("#BSUB -oo /dev/null\n")
        submit_file.write("#BSUB -eo /dev/null\n")
        submit_file.write("\n")

        submit_file.write("source {}/load_modules.sh\n".format(self.working_directory))
        submit_file.write("export GREASY_NWORKERS=$LSB_DJOB_NUMPROC\n")
        submit_file.write("export GREASY_LOGFILE={}/{}_$LSB_JOBINDEX.log\n".format(self.submit_dir, self.unique_ID))
        submit_file.write("arguments_store={}/{}.grz\n".format(self.arguments_dir, self.unique_ID))
        submit_file.write("argument_file=$(cat $arguments_store | awk -v var=$LSB_JOBINDEX 'NR==var {print $1}')\n")
        submit_file.write("\n")
        submit_file.write("greasy $argument_file")

        # close submit file
        submit_file.close()

    
    def submit_job(self):

        # time start of interpolation jobs
        interpolation_start = time.time()

        # submit slurm script
        submit_complete = False
        while submit_complete == False:

            if self.machine == "nord3":
                submit_process = subprocess.Popen(['bsub'], stdout=subprocess.PIPE,
                                              stdin=open('{}/{}'.format(self.submit_dir, self.job_fname), 'r'))
            else:
                submit_process = subprocess.Popen(['sbatch', self.job_fname],
                                              cwd=self.submit_dir, stdout=subprocess.PIPE)
            submit_status = submit_process.communicate()[0]
            submit_return_code = submit_process.returncode

            # if sbatch fails, time out for 60 seconds and then try again
            if submit_return_code != 0:
                time.sleep(60)
                continue
            else:
                submit_complete = True

            # take a 1 second pause between submittals (to help slurm scheduler)
            time.sleep(1)

        # monitor number of jobs in queue (every 10 seconds) until there are 0 left in the squeue
        all_tasks_finished = False
        # flag for monitoring jobs that were submitted in n3
        job_entered = False 

        while all_tasks_finished == False:
            if self.machine == "nord3":
                # cmd = ['bjobs', '-noheader', '-J', 'PRVI_{}[1]'.format(self.unique_ID)]
                cmd = ['bjobs', '-noheader']
            else:
                cmd = ['squeue', '-h', '-n', 'PRVI_{}'.format(self.unique_ID)]
            squeue_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding='utf8')
            squeue_status = squeue_process.communicate()[0]
            n_jobs_in_queue = len(squeue_status.split('\n')[:-1])
            # if number of jobs in queue > 0, then sleep for 10
            # seconds and then check again how many jobs there are in queue
            if (self.machine in ('mn4', 'power','nord3v2', 'amd', 'mn5')) and (n_jobs_in_queue > 0):
                time.sleep(10)
                continue
            elif self.machine == 'nord3':
                # has submitted jobs entered the queue?
                if self.unique_ID[1:] in squeue_status.split('\n')[1:][0]:
                    job_entered = True
                    time.sleep(10)
                    continue
                # if not, then wait
                if not job_entered:
                    time.sleep(10)
                    continue
                else:
                    # if submitted job entered but still running, wait
                    if n_jobs_in_queue > 1:
                        time.sleep(10)
                        continue

            # if no more jobs in the squeue, then now check the outcome of all the jobs
            # if any jobs have failed/not finished, write them out to file
            failed_tasks=[]
            not_finished_tasks=[]
            process_times=[]
            for output_log_root in self.output_log_roots:
                output_log_file = glob.glob('{}_*'.format(output_log_root))
                # have an output log file? (i.e. job has finished)
                if len(output_log_file) > 0:
                    output_log_file = output_log_file[0]
                    process_code = int(output_log_file.split('_')[-1].split('.out')[0])
                    # if process code == 0, job completed successfully
                    if process_code == 0:
                        # append interpolation job process time
                        process_times.append(float(subprocess.check_output(
                            ['tail', '-1', output_log_file], encoding='utf8').strip()))
                    # otherwise, job failed --> append failed log file
                    else:
                        failed_tasks.append(output_log_file)
                # no output log file, therefore append to not finished list
                else:
                    not_finished_tasks.append(output_log_root)

            # break out of while loop
            all_tasks_finished = True

        # get interpolation time
        interpolation_time = (time.time() - interpolation_start)/60.

        # stop timer
        total_time = (time.time()-self.start)/60.

        # have 0 failed/non-finished tasks?
        if (len(failed_tasks) == 0) & (len(not_finished_tasks) == 0):
            # get queue time
            process_time = np.max(process_times)
            queue_time = interpolation_time - process_time
            overhead_time = total_time - interpolation_time
            print('ALL {} INTERPOLATION TASKS COMPLETED SUCCESSFULLY IN {:.2f} MINUTES\n'
                  '({:.2f} MINUTES PROCESING, {:.2f} MINUTES QUEUING, {:.2f} MINUTES ON OVERHEADS)'.
                  format(len(self.output_log_roots), total_time, process_time, queue_time, overhead_time))
        else:
            print('{}/{} INTERPOLATION TASKS FINISHED SUCCESSFULLY IN {:.2f} MINUTES'.format(
                  len(self.output_log_roots)-(len(not_finished_tasks)+len(failed_tasks)), len(self.output_log_roots), 
                  total_time))
            if len(failed_tasks) > 0:
                print('THE FOLLOWING INTERPOLATION TASKS FAILED: {}'.format(failed_tasks))
            if len(not_finished_tasks) > 0:
                print('THE FOLLOWING INTERPOLATION TASKS DID NOT FINISH: {}'.format(not_finished_tasks))


if __name__ == "__main__":

    # initialise SubmitInterpolation object
    SI = SubmitInterpolation()

    # get all unique arguments to process interpolation tasks
    SI.gather_arguments()
    
    # create greasy arguments file
    SI.create_greasy_arguments_file()

    # create submission script according to machine
    if SI.machine == "nord3":
        SI.create_lsf_submission_script()
    else:
        SI.create_slurm_submission_script()
    
    # submit interpolation jobs
    SI.submit_job()
