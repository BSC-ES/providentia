import copy
import glob
import os
import yaml
import random
import subprocess
import sys
import time
from pydoc import locate
import numpy as np
import multiprocessing

MACHINE = os.environ.get('BSC_MACHINE', 'local')

# get current path and providentia root path
CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = os.path.dirname(os.path.dirname(CURRENT_PATH))

sys.path.append(os.path.join(PROVIDENTIA_ROOT, 'providentia', 'interpolation'))
sys.path.append(os.path.join(PROVIDENTIA_ROOT, 'providentia'))

from aux import get_aeronet_bin_radius_from_bin_variable, get_model_bin_radii, check_for_ghost
from mapping_species import mapping_species
from configuration import ProvConfiguration, load_conf

# load the defined experiments paths and agrupations jsons
interp_experiments = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings', 'interp_experiments.yaml')))

class SubmitInterpolation(object):
    """ Class that handles the interpolation submission. """

    def __init__(self,**kwargs):

        # start timer
        self.start = time.time()

        # define current working directory and
        # arguments/greasy/interpolation log subdirectories
        self.working_directory = CURRENT_PATH     
        self.arguments_dir = os.path.join(PROVIDENTIA_ROOT, 'logs/interpolation/arguments')
        self.submit_dir = os.path.join(PROVIDENTIA_ROOT, 'logs/interpolation/greasy_logs')
        self.interpolation_log_dir = os.path.join(PROVIDENTIA_ROOT, 'logs/interpolation/interpolation_logs')

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

        # get section args
        if "section" in kwargs:
            section = kwargs["section"]
            if section in self.parent_section_names:
                self.current_config = self.sub_opts[section]
            else:
                error = f"Error: Defined section '{section}' does not exist in configuration file."
                sys.exit(error)
        else:
            # if no parent section names are found throw an error
            if len(self.parent_section_names) == 0:
                error = "Error: No sections were found in the configuration file, make sure to name them using square brackets."
                sys.exit(error)
            self.current_config = self.sub_opts[self.parent_section_names[0]]
            print(f"Taking first defined section ({self.parent_section_names[0]}) to be read.")

        # dictionary that stores utilized interpolation variables
        self.interpolation_variables = {}

        # update variables from defined config file
        if self.current_config:
            for k, val in sorted(self.current_config.items()):
                setattr(self, k, provconf.parse_parameter(k, val))

        # now all variables have been parsed, check validity of those, throwing errors where necessary
        provconf.check_validity()

        # TODO Hardcoded
        interp_print_variables = ['ghost_version', 'start_date', 'end_date', 'experiments', 
                                  'species', 'network', 'resolution', 'forecast', 'forecast_day', 
                                  'interp_n_neighbours', 'interp_reverse_vertical_orientation', 
                                  'interp_chunk_size', 'interp_job_array_limit', 'exp_root', 
                                  'ghost_root', 'nonghost_root', 'exp_to_interp_root', 'interp_multiprocessing']

        # print variables used, if all species are used print "All Species"        
        print("\nVariables used for the interpolation:\n")
        for arg in interp_print_variables:
            if arg != "experiments":
                print(f"{arg}: {getattr(self, arg)}")
            else:
                print(f"{arg}:")
                for exp, alias in getattr(self, arg).items():
                    if self.alias_flag:
                        print(f" - {exp} ({alias})")
                    else:
                        print(f" - {exp}")

        # define the QOS (Quality of Service) used to manage jobs on the SLURM system
        if MACHINE == 'mn5':
            self.qos = 'gp_bsces'
        else:
            self.qos = 'bsc_es'

        # import unit converter
        sys.path.append(os.path.join(PROVIDENTIA_ROOT, 'providentia', 'dependencies','unit-converter'))
        import unit_converter

    def gather_arguments(self):
        ''' Gather list of arguments for all unique tasks to process, as defined in the configuration file. '''

        # create arguments list
        self.arguments = []

        # create list to hold root directory locations of .out log files generated for each task processed
        self.output_log_roots = []  
        
        # create list of arguments to compare between the iterations
        last_arguments = None

        # iterate through desired experiment IDs and its types
        for (exp_dom_ens, alias), experiment_type in zip(self.experiments.items(), self.experiment_types):
            
            experiment_to_process, grid_type, ensemble_option = exp_dom_ens.split("-") 

            print('\nEXPERIMENT: {0}'.format(alias))

            # create arguments list
            exp_arguments = []

            # initialise files lists
            obs_files = []
            exp_files = []
            obs_files_dates = []
            exp_files_dates = []

            exp_dir = None
            # if local machine, get directory from data_paths
            if MACHINE == 'local':  
                exp_to_interp_path = os.path.join(self.exp_to_interp_root, experiment_to_process)
                if os.path.exists(exp_to_interp_path):
                    exp_dir = exp_to_interp_path
                
                # if none of the paths are in this current machine, break
                if exp_dir is None:
                    error = f"Error: None of the experiment paths in {self.exp_to_interp_root} are available in this machine ({MACHINE})."
                    sys.exit(error)
            # for HPC machines
            else:
                # get experiment type and specific directories
                exp_dir_list = interp_experiments[experiment_type]["paths"]
                for temp_exp_dir in exp_dir_list:
                    if os.path.exists(os.path.join(temp_exp_dir,experiment_to_process)):
                        exp_dir = temp_exp_dir
                        break

                # take first functional directory 
                if exp_dir is None:
                    error = f"Error: None of the experiment paths in {os.path.join('settings', 'interp_experiments.yaml')} are available in this machine ({MACHINE})."
                    sys.exit(error)
                
                # add file to directory path
                exp_dir += f"{experiment_to_process}/"

            # get model bin edges
            r_edges, rho_bins = get_model_bin_radii(experiment_type)

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
                                        if ('vconcaerobin' in speci_to_process) or ('vconcaerobin' in speci_to_map):
                                            # check if bin radius is within model's bin extents
                                            if 'vconcaerobin' in speci_to_process:
                                                speci_to_check = copy.deepcopy(speci_to_process)
                                            elif 'vconcaerobin' in speci_to_map:
                                                speci_to_check = copy.deepcopy(speci_to_map)
                                            bin_radius = get_aeronet_bin_radius_from_bin_variable(speci_to_check)
                                            if (bin_radius >= r_edges[0]) & (bin_radius <= r_edges[-1]):
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
                                        if ('vconcaerobin' in speci_to_process) or ('vconcaerobin' in speci_to_map):
                                            # check if bin radius is within model's bin extents
                                            if 'vconcaerobin' in speci_to_process:
                                                speci_to_check = copy.deepcopy(speci_to_process)
                                            elif 'vconcaerobin' in speci_to_map:
                                                speci_to_check = copy.deepcopy(speci_to_map)
                                            bin_radius = get_aeronet_bin_radius_from_bin_variable(speci_to_check)
                                            if (bin_radius >= r_edges[0]) & (bin_radius <= r_edges[-1]):
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
                                    self.ghost_root + '/{}/{}/{}/{}/{}*.nc'.format(
                                        network_to_interpolate_against, self.ghost_version, 
                                        temporal_resolution_to_output, original_speci_to_process, 
                                        original_speci_to_process)))
                            # non-GHOST
                            else:
                                obs_files = np.sort(glob.glob(
                                    self.nonghost_root + '/{}/{}/{}/{}*.nc'.format(
                                        network_to_interpolate_against, temporal_resolution_to_output,
                                        original_speci_to_process, original_speci_to_process)))

                            # if have no observational files then continue
                            if len(obs_files) == 0:
                                continue
                            
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
                                # if the ensemble number is in the nc file name next to species
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
                                # if there's no ensemble number in the file name
                                else:
                                    have_ensemble_members = False
                                    # if have defined ensemble members to process, then continue as no files in this 
                                    # directory have ensemble member number
                                    if ensemble_option not in ['allmembers', '000']:
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
                                    self.arguments.append("{} {} {} {} {} {} {} {}".format(prov_exp_code, 
                                                                                        model_temporal_resolution, 
                                                                                        speci_to_process, 
                                                                                        network_to_interpolate_against, 
                                                                                        temporal_resolution_to_output, 
                                                                                        yearmonth, 
                                                                                        original_speci_to_process,
                                                                                        self.slurm_job_id))

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

            # boolean that says if there was any new valid data in the iteration
            new_arguments = last_arguments != self.arguments

            # if list is empty or have no arguments after iteration, return message stating that
            if len(self.arguments) == 0 or not new_arguments:
                msg = 'NO INTERSECTING OBSERVATIONAL AND EXPERIMENT DATA FOR INTERPOLATION. \n' 
                if self.start_date == self.end_date:
                    msg += 'If you want to interpolate data for one month, '
                    msg += 'you need to set the end date to be the next one. \n'
                    msg += 'e.g. For November 2018, this is 201811 to 201812.'
                else:
                    if len(exp_files) == 0:
                        msg += f'Experiment data between {self.start_date} and {self.end_date} cannot be found.'
                    else:
                        msg += f'Experiment files: {exp_files}\n'
                        msg += f'Experiment dates: {exp_files_dates}'
                        if len(obs_files) == 0:
                            msg += f'Observational data between {self.start_date} and {self.end_date} cannot be found.'
                        else:
                            msg += f'Observational files: {obs_files}\n'
                            msg += f'Observational dates: {obs_files_dates}\n'
            else:
                msg = '***INTERSECTING OBSERVATIONAL AND EXPERIMENTAL DATA IS AVAILABLE FOR INTERPOLATION.***' 

            print(msg)

            # get the arguments from the last iteration
            last_arguments = copy.deepcopy(self.arguments)
        
        # if have no arguments for all experiments, return message stating that
        if len(self.arguments) == 0:
            if len(self.experiments) > 1:
                error = 'INTERPOLATION CANNOT BE DONE FOR ANY EXPERIMENT'
                sys.exit(error)
            else:
                sys.exit()

        # randomise the order of the arguments list
        random.shuffle(self.arguments)     

    def create_greasy_arguments_file(self):
        ''' Create greasy arguments text file storing all different tasks to run by greasy. '''

        # define list to store chunked argument files (to be submitted using greasy)
        argument_files = []

        # get the CPU chunk size -- set initially as miniumum number of arguments per file 
        N_arguments_per_file_minimum = copy.deepcopy(int(self.interp_chunk_size))
        N_arguments_per_file = copy.deepcopy(int(self.interp_chunk_size))
        
        # divide the number of arguments by the CPU chunk size, to determine how many argument files will be needed to 
        # submit all jobs
        N_submit_files = int(np.ceil(len(self.arguments)/int(self.interp_chunk_size)))
        
        # set argument remainder as 0 initially
        argument_remainder = 0

        # if the number of argument files is greater than the job array limit (i.e the limit on the number of argument 
        # files that can be processed simultaneously)
        # then add adjust minimum N arguments per file appropriately (i.e. split extra arguments across the maximum 
        # number of argument files evenly)        
        if N_submit_files > int(self.interp_job_array_limit):
            # update the minimum number of arguments per file
            N_arguments_per_file_minimum = int(np.floor(len(self.arguments)/int(self.interp_job_array_limit)))
            N_arguments_per_file = copy.deepcopy(N_arguments_per_file_minimum)

            # if the number of extra arguments does not divide equally between all files get the remainder
            argument_remainder = int(len(self.arguments)%int(self.interp_job_array_limit))
            
            # if have argument remainder then update N_arguments_per_file variable to be 1 greater than minimum for 
            # first file written (and for all files thereafter until  remainder is accounted for)
            if argument_remainder > 0:
                N_arguments_per_file = N_arguments_per_file_minimum + 1
                # subtract 1 from the argument remainder
                argument_remainder -= 1
            
            # set N submit files as N of job array limit
            N_submit_files = copy.deepcopy(int(self.interp_job_array_limit))

        # create file which will store a list of all chunked argument filenames
        greasy_file = open('{}/{}.grz'.format(self.arguments_dir, self.slurm_job_id), 'w')
        
        # create all chunked argument filenames
        for ii in range(N_submit_files):
            argument_fname = '{}/{}_{}.txt'.format(self.arguments_dir,self.slurm_job_id,ii)
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
        self.job_fname = self.slurm_job_id+'.sh'
    
        # get all argument files
        argument_files = sorted(glob.glob('{}/{}_*.txt'.format(self.arguments_dir,self.slurm_job_id)))

        # read how many lines are in first arguments file
        with open(argument_files[0]) as f: 
            for ii, line in enumerate(f):
                pass
            N_arguments = ii + 1

        # cap the number of simultaneously running tasks to be the defined CPU chunk size  
        max_tasks = copy.deepcopy(int(self.interp_chunk_size))

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
        submit_file.write("#SBATCH --job-name=PRVI_{}\n".format(self.slurm_job_id))
        submit_file.write("#SBATCH --ntasks={}\n".format(n_simultaneous_tasks))
        # fix number of nodes to be 1 (for faster execution)
        submit_file.write("#SBATCH --nodes=1\n")
        submit_file.write("#SBATCH --time=48:00:00\n")
        submit_file.write("#SBATCH --array=1-{}\n".format(len(argument_files)))
        submit_file.write("#SBATCH --qos={}\n".format(self.qos))
        # submit_file.write("#SBATCH --output=/dev/null\n") # decomment when debugging
        # submit_file.write("#SBATCH --error=/dev/null\n")
        if MACHINE == 'mn5': # TODO when checking if debug works check this
            submit_file.write("#SBATCH --account=bsc32\n")  
            submit_file.write("#SBATCH --ntasks-per-node={}\n".format(n_simultaneous_tasks))
            submit_file.write("#SBATCH --cpus-per-task=1\n")
            submit_file.write("\n")
        else:
            submit_file.write("\n")
            submit_file.write("source {}/bin/load_modules.sh\n".format(PROVIDENTIA_ROOT))
        submit_file.write("export GREASY_NWORKERS=$SLURM_NPROCS\n") 
        submit_file.write("export GREASY_LOGFILE={}/{}_$SLURM_ARRAY_TASK_ID.log\n".format(self.submit_dir, 
                                                                                          self.slurm_job_id))
        submit_file.write("export SLURM_CPU_BIND=none\n")
        submit_file.write("arguments_store={}/{}.grz\n".format(self.arguments_dir, self.slurm_job_id))
        submit_file.write("argument_file=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')\n")
        submit_file.write("\n")
        submit_file.write("greasy $argument_file")

        # close submit file
        submit_file.close()

    def create_lsf_submission_script(self):
        """ Write a lsf submission shell script that submits a greasy job. """

        # create job_fname (unique_ID + 'sh.')
        self.job_fname = self.slurm_job_id + '.sh'

        # get all argument files
        argument_files = sorted(glob.glob('{}/{}_*.txt'.format(self.arguments_dir, self.slurm_job_id)))

        # read how many lines are in first arguments file
        with open(argument_files[0]) as f:
            for ii, line in enumerate(f):
                pass
            N_arguments = ii + 1

        # cap the number of simultaneously running tasks to be the defined CPU chunk size
        max_tasks = copy.deepcopy(int(self.interp_chunk_size))

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
        submit_file.write("#BSUB -J PRVI_{}[1-{}]\n".format(self.slurm_job_id, len(argument_files)))
        submit_file.write("#BSUB -q {}\n".format(self.qos))
        submit_file.write("#BSUB -oo /dev/null\n")
        submit_file.write("#BSUB -eo /dev/null\n")
        submit_file.write("\n")

        submit_file.write("source {}/bin/load_modules.sh\n".format(PROVIDENTIA_ROOT))
        submit_file.write("export GREASY_NWORKERS=$LSB_DJOB_NUMPROC\n")
        submit_file.write("export GREASY_LOGFILE={}/{}_$LSB_JOBINDEX.log\n".format(self.submit_dir, self.slurm_job_id))
        submit_file.write("arguments_store={}/{}.grz\n".format(self.arguments_dir, self.slurm_job_id))
        submit_file.write("argument_file=$(cat $arguments_store | awk -v var=$LSB_JOBINDEX 'NR==var {print $1}')\n")
        submit_file.write("\n")
        submit_file.write("greasy $argument_file")

        # close submit file
        submit_file.close()

    def submit_job_greasy(self):

        # time start of interpolation jobs
        interpolation_start = time.time()

        # submit slurm script
        submit_complete = False
        while submit_complete == False:

            if MACHINE == "nord3":
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
            if MACHINE == "nord3":
                # cmd = ['bjobs', '-noheader', '-J', 'PRVI_{}[1]'.format(self.slurm_job_id)]
                cmd = ['bjobs', '-noheader']
            else:
                cmd = ['squeue', '-h', '-n', 'PRVI_{}'.format(self.slurm_job_id)]
            squeue_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding='utf8')
            squeue_status = squeue_process.communicate()[0]
            n_jobs_in_queue = len(squeue_status.split('\n')[:-1])
            # if number of jobs in queue > 0, then sleep for 10
            # seconds and then check again how many jobs there are in queue
            if (MACHINE in ('nord3v2', 'amd', 'mn5')) and (n_jobs_in_queue > 0):
                time.sleep(10)
                continue
            elif MACHINE == 'nord3':
                # has submitted jobs entered the queue?
                if self.slurm_job_id[1:] in squeue_status.split('\n')[1:][0]:
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
            print('\nALL {} INTERPOLATION TASKS COMPLETED SUCCESSFULLY IN {:.2f} MINUTES\n'
                  '({:.2f} MINUTES PROCESSING, {:.2f} MINUTES QUEUING, {:.2f} MINUTES ON OVERHEADS)'.
                  format(len(self.output_log_roots), total_time, process_time, queue_time, overhead_time))
        else:
            print('\n{}/{} INTERPOLATION TASKS FINISHED SUCCESSFULLY IN {:.2f} MINUTES'.format(
                  len(self.output_log_roots)-(len(not_finished_tasks)+len(failed_tasks)), len(self.output_log_roots), 
                  total_time))
            if len(failed_tasks) > 0:
                print('THE FOLLOWING INTERPOLATION TASKS FAILED: {}'.format(failed_tasks))
            if len(not_finished_tasks) > 0:
                print('THE FOLLOWING INTERPOLATION TASKS DID NOT FINISH: {}'.format(not_finished_tasks))
    
    def run_command(self, commands):
        arguments_list = commands.strip().split()
        subprocess.run(arguments_list, capture_output=True, text=True)

    def submit_job_multiprocessing(self):
        
        # launch interpolation
        commands = [f'python -u {CURRENT_PATH}/experiment_interpolation.py '
                    + argument for argument in self.arguments]
        with multiprocessing.Pool(processes=self.n_cpus) as pool:
            pool.map(self.run_command, commands)

        # stop timer
        total_time = (time.time()-self.start)/60.

        # if no more jobs in the squeue, then now check the outcome of all the jobs
        # if any jobs have failed/not finished, write them out to file
        failed_tasks=[]
        not_finished_tasks=[]
        for output_log_root in self.output_log_roots:
            output_log_file = glob.glob('{}_*'.format(output_log_root))
            # have an output log file? (i.e. job has finished)
            if len(output_log_file) > 0:
                output_log_file = output_log_file[0]
                process_code = int(output_log_file.split('_')[-1].split('.out')[0])
                # if process code is not 0, job failed
                if process_code != 0:
                    failed_tasks.append(output_log_file)
            # no output log file, therefore append to not finished list
            else:
                not_finished_tasks.append(output_log_root)

        # have 0 failed/non-finished tasks?
        if (len(failed_tasks) == 0) & (len(not_finished_tasks) == 0):
            print('\nALL {} INTERPOLATION TASKS COMPLETED SUCCESSFULLY IN {:.2f} MINUTES\n'.
                  format(len(self.output_log_roots), total_time))
        else:
            print('\n{}/{} INTERPOLATION TASKS FINISHED SUCCESSFULLY IN {:.2f} MINUTES'.format(
                  len(self.output_log_roots)-(len(not_finished_tasks)+len(failed_tasks)), len(self.output_log_roots), 
                  total_time))
            if len(failed_tasks) > 0:
                print('THE FOLLOWING INTERPOLATION TASKS FAILED: {}'.format(failed_tasks))
            if len(not_finished_tasks) > 0:
                print('THE FOLLOWING INTERPOLATION TASKS DID NOT FINISH: {}'.format(not_finished_tasks))
    

def main(**kwargs):

    # initialise SubmitInterpolation object
    SI = SubmitInterpolation(**kwargs)

    # get all unique arguments to process interpolation tasks
    SI.gather_arguments()

    # create greasy arguments file
    SI.create_greasy_arguments_file()

    # submit interpolation jobs
    if SI.interp_multiprocessing:
        SI.submit_job_multiprocessing()
    else:
        if SI.machine == 'local':
            error = 'Error: It is not possible to interpolate locally without using multiprocessing.'
            sys.exit(error)
        else:
            # create submission script according to machine
            if SI.machine == "nord3":
                SI.create_lsf_submission_script()
            else:
                SI.create_slurm_submission_script()
            SI.submit_job_greasy()