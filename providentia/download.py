import sys
import os
import shutil

from base64 import decodebytes
import cdsapi
import copy
from dotenv import dotenv_values
import json 
from netCDF4 import Dataset
import numpy as np
import paramiko 
import re 
import requests
import signal
import subprocess
import tarfile
import time
from tqdm import tqdm
import xarray as xr
import yaml
import zipfile   

from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from getpass import getpass

from .actris import (get_files_path, temporally_average_data, get_data,
                     get_files_per_var, is_wavelength_var, get_files_to_download,
                     get_files_info, parameters_dict)
from providentia.auxiliar import CURRENT_PATH, join
from .configuration import ProvConfiguration, load_conf
from .read_aux import check_for_ghost
from .warnings_prv import show_message


PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)
REMOTE_MACHINE = "storage5"

# load the defined experiments paths, agrupations yaml and mapping species
data_paths = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/data_paths.yaml')))
interp_experiments = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'interp_experiments.yaml')))
mapping_species =  yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'mapping_species.yaml')))
cams_options = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'cams_dataset.yaml')))
cams_variables_level = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'cams_variables_level.yaml')))

class Download(object):
    def __init__(self, **kwargs):

        # get providentia start time
        self.prov_start_time = time.time()

        # get ssh user and password 
        env = dotenv_values(join(PROVIDENTIA_ROOT, ".env"))

        # get ssh user and password 
        self.prv_user = env.get("PRV_USER")
        self.prv_password = env.get("PRV_PWD")

        # get user preference over GHOST download, overwrite and origin update (ACTRIS)
        self.overwrite_choice = env.get("OVERWRITE")
        self.origin_update_choice = env.get("ORIGIN_UPDATE")

        # set timeout limit
        self.timeoutLimit = 3 * 60

        # initialise default configuration variables
        # modified by commandline arguments, if given
        self.provconf = ProvConfiguration(self, **kwargs)

        self.logger.info("Starting Providentia download...")

        # update variables from config file
        if self.config != '':  
            read_conf = False
            if os.path.exists(self.config):
                read_conf = True
            else: 
                if os.path.exists(join(self.config_dir, self.config)):
                    self.config = join(self.config_dir, self.config)
                    read_conf = True
            if read_conf:
                load_conf(self, self.config)
                self.from_conf = True
                # get the section in case it was passed on the command line                
                if 'section' in kwargs:
                    # config and section defined 
                    if kwargs['section'] in self.all_sections:
                        self.sections = [kwargs['section']]
                    else:
                        error = 'Error: The section specified in the command line ({0}) does not exist.'.format(kwargs['section'])
                        tip = 'Tip: For subsections, add the name of the parent section followed by an interpunct (·) '
                        tip += 'before the subsection name (e.g. SECTIONA·Spain).'
                        error = error + '\n' + tip
                        self.logger.error(error)
                        sys.exit(1)
                # if no section passed, then get all the parent sections
                else:
                    # if no parent section names are found throw an error
                    if len(self.parent_section_names) == 0:
                        error = "Error: No sections were found in the configuration file, make sure to name them using square brackets."
                        self.logger.error(error)
                        sys.exit(1)
                    self.sections = self.parent_section_names
            else:
                error = 'Error: The path to the configuration file specified in the command line does not exist.'
                self.logger.error(error)
                sys.exit(1)
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            self.logger.error(error)
            sys.exit(1)

        # variable that saves whether some experiments/observations were downloaded before
        self.overwritten_files_flag = False

        # initialize the necessary things in local
        if self.machine == "local":

            # TODO move some of these variables to configuration.py
            # initialise zenodo url
            self.ghost_url = 'https://zenodo.org/records/10637450'

            # initialise remote hostname
            self.remote_hostname = "transfer1.bsc.es"

            # create empty directories for all the paths if they don't exist
            for path in [self.nonghost_root,self.ghost_root,self.exp_root,self.exp_to_interp_root]:
                if not os.path.exists(path):
                    try:
                        os.makedirs(path)
                    except PermissionError as error:
                        os.system(f"sudo mkdir -p {path}")
                        os.system(f"sudo chmod o+w {path}")

            # initialise ssh 
            self.ssh = None

            # initialise boolean thath indicates whether remote machine changed 
            self.switched_remote = False

    def run(self):
        for section_ind, section in enumerate(self.sections):
            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # update self with section variables
            for k, val in self.section_opts.items():
                setattr(self, k, self.provconf.parse_parameter(k, val))

            # if there's experiments, ask the user whether they want interpolated or non-interpolated
            if self.experiments:
                # self.interpolated = input("Experiments were detected in the configuration file. Do you want to download the interpolated versions? (Otherwise, the non-interpolated experiments will be downloaded) ([y]/n): ")
                # while self.interpolated.lower() not in ['','y','n']:
                #     self.interpolated = input("Experiments were detected in the configuration file. Do you want to download the interpolated versions? (Otherwise, the non-interpolated experiments will be downloaded) ([y]/n): ")

                # # set the interpolated parameter  
                self.interpolated = False

            # now all variables have been parsed, check validity of those, throwing errors where necessary
            self.provconf.check_validity(deactivate_warning=True)

            # from here generate control if user stopped execution
            signal.signal(signal.SIGINT, self.sighandler)
            
            # only the local download iterates through the networks
            if self.machine in "local":
                # if networks is none and is not the non interpolated mode, raise error
                if not self.network and self.interpolated is True:
                    error = "Error: No networks were passed."
                    self.logger.error(error)
                    sys.exit(1)
                
                # when one of those symbols is passed, get all networks
                if self.network == ["*"]:
                    self.get_all_networks()

                # networks
                if self.network:
                    # combine all networks and species combinations to download (for network and filter species)
                    combined_networks = [(network, None) for network in self.network] + \
                                        [(network_specie.split('|')[0], network_specie.split('|')[1]) for network_specie in self.filter_species]
                    
                    # save main species
                    main_species = copy.deepcopy(self.species)
                    
                    # download network observations with species and filter_species
                    for network, filter_species in combined_networks:
                        # change species when turn of filter species
                        if filter_species is not None:
                            self.species = [filter_species]

                        # get the files to be downloaded, check if they were already downloaded 
                        # GHOST
                        if check_for_ghost(network):
                            # ask whether the user wants to download from the zenodo or bsc machine
                            bsc_download = input("GHOST network detected. Download from the BSC remote machine? (Otherwise, it will be retrieved from Zenodo) ([y]/n): ")
                            while bsc_download.lower() not in ['','y','n']:
                                bsc_download = input("GHOST network detected. Download from the BSC remote machine? (Otherwise, it will be retrieved from Zenodo) ([y]/n): ")
                            
                            # get the zenodo or bsc machine depending on the user answer
                            download_ghost = self.download_ghost_network_sftp if bsc_download.lower() in ['','y'] else self.download_ghost_network_zenodo
                            
                            # download GHOST network
                            initial_check_nc_files = download_ghost(network, initial_check=True)
                            files_to_download = self.select_files_to_download(initial_check_nc_files)
                            if not initial_check_nc_files or files_to_download:
                                download_ghost(network, initial_check=False, files_to_download=files_to_download)
                        # ACTRIS
                        elif network == 'actris/actris':
                            self.download_actris_network()
                        # non-GHOST
                        else:
                            initial_check_nc_files = self.download_nonghost_network(network, initial_check=True)
                            files_to_download = self.select_files_to_download(initial_check_nc_files)
                            if not initial_check_nc_files or files_to_download:
                                self.download_nonghost_network(network, initial_check=False, files_to_download=files_to_download)

                    # get orignal species back
                    self.species = main_species

            # when one of those symbols is passed, get all experiments
            if self.experiments == {'*' : '*'}:
                self.get_all_experiments()

            if self.experiments:
                # remote machine experiment download
                if self.machine in ["storage5", "nord3v2", "nord4"]:
                    # get function to download experiment depending on the configuration file field
                    download_experiment_fun = self.copy_non_interpolated_experiment
                # local experiment download
                else:
                    for experiment in self.experiments:
                        # CAMS experiment
                        if experiment.startswith(tuple(cams_options.keys())):
                            self.download_cams_experiment(experiment)
                        # BSC machines
                        else:
                            download_experiment_fun = self.download_experiment if self.interpolated else self.download_non_interpolated_experiment
                    
                            # iterate the experiments download
                            for experiment in self.experiments.keys():
                                initial_check_nc_files = download_experiment_fun(experiment, initial_check=True)
                                files_to_download = self.select_files_to_download(initial_check_nc_files)
                                if not initial_check_nc_files or files_to_download:
                                    download_experiment_fun(experiment, initial_check=False, files_to_download=files_to_download)

            # remove section variables from memory
            for k in self.section_opts:
                try:
                    vars(self).pop(k)
                except:
                    pass

            # reset domain and ensemble options for new section
            self.domain = []
            self.ensemble_options = []

        # show message in case experiments or observations were ignored
        if self.overwritten_files_flag == True:
            self.logger.info("\nSome experiments/observations were found but were not downloaded because the OVERWRITE option is set to 'n'.")

        if self.machine == "local":
            # close connection, if it exists
            if self.ssh is not None:
                self.ssh.close() 
                self.sftp.close()

    def connect(self):
        # declare that we are using the remote machine
        global REMOTE_MACHINE
        
        # initialise the paths
        self.ghost_remote_obs_path = data_paths[REMOTE_MACHINE]["ghost_root"]
        self.nonghost_remote_obs_path = data_paths[REMOTE_MACHINE]["nonghost_root"]
        self.exp_remote_path = data_paths[REMOTE_MACHINE]["exp_root"]
        self.exp_to_interp_remote_path = data_paths[REMOTE_MACHINE]["exp_to_interp_root"]

        # get public remote machine public key and add it to ssh object
        _, output = subprocess.getstatusoutput(f"ssh-keyscan -t ed25519 {self.remote_hostname}")

        ed25519_key = output.split()[-1].encode()
        key = paramiko.Ed25519Key(data=decodebytes(ed25519_key))
        
        # encode the output public key if possible
        try:
            ed25519_key = output.split()[-1].encode()
            key = paramiko.Ed25519Key(data=decodebytes(ed25519_key))
        
        # in case transfer broke
        except:
            msg = f"Remote machine {REMOTE_MACHINE} not working right now."

            # if the remote machine has not been changed
            if not self.switched_remote:
                # change remote machine and hostname
                self.remote_hostname, REMOTE_MACHINE = "glogin4.bsc.es", "mn5" 
                msg += f" Changing it to {REMOTE_MACHINE}..."

            show_message(self, msg)
            
            # if the remote machine has not been changed
            if not self.switched_remote:
                # connect but with the new machine
                self.switched_remote = True
                return self.connect()
            # if it has been changed already, exit
            else:
                error = "Error: None of the machines are working right now. Try later."
                self.logger.error(error)
                sys.exit(1)

        self.ssh = paramiko.SSHClient()
        hostkeys = self.ssh.get_host_keys().add(self.remote_hostname, 'ed25519', key)

        # initialise temporal variables
        prv_user, prv_password = None, None

        # if couldn't get user, ask for it
        if self.prv_user is None:
            prv_user = ''
            while prv_user == '':
                prv_user = input(f"\nInsert BSC {REMOTE_MACHINE} ssh user: ")
            self.prv_user = prv_user
        
        # if couldn't get user, check if you have to ask for it
        if self.prv_password is None:
            # check if user needs a password
            try:
                self.ssh.connect(self.remote_hostname, username=self.prv_user, password='placeholder')
            # if authentication error, that means that the user or and the password are wrong
            except paramiko.ssh_exception.AuthenticationException:
                # if name was not changed, then user in .env is not valid
                if prv_user is None:
                    error = f"Authentication failed. Please, check if PRV_USER on {join(PROVIDENTIA_ROOT, '.env')} aligns with your BSC {REMOTE_MACHINE} ssh user."
                    error += "\nIf it does not, change the user to the correct one. If it does, delete the whole PRV_USER row and execute again."
                    self.logger.error(error)
                    sys.exit(1)
                else:
                    prv_password = getpass("Insert password: ")
                    self.prv_password = prv_password

        # catch identification method
        try:
            # connect through ssh and create Secure File Transfer Protocol object
            self.ssh.connect(self.remote_hostname, username=self.prv_user, password=self.prv_password)
            self.sftp = self.ssh.open_sftp()
            
        # if credentials are invalid, throw an error
        except paramiko.ssh_exception.AuthenticationException:
            error = "Authentication failed."
            # if user or password were taken from .env (did not change), tell the user to check .env
            if prv_user is None:
                error += f" Please, check your credentials on {join(PROVIDENTIA_ROOT, '.env')}"
            self.logger.error(error)
            sys.exit(1)

        # if pwd or user changed, ask if user wants to remember credentials
        if (prv_user is not None) or (prv_password is not None):
            # ask user if they want their credentials saved
            remind_txt = input("\nRemember credentials (y/n)? ")
            while remind_txt.lower() not in ['y','n']:
                remind_txt = input("\nRemember credentials (y/n)? ")
            
            # create .env with the input user and/or password
            if remind_txt.lower() == 'y':
                with open(join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                    if prv_user is not None:
                        f.write(f"PRV_USER={self.prv_user}\n")
                    if prv_password is not None:
                        f.write(f"PRV_PWD={self.prv_password}\n")

                self.logger.info(f"\nRemote machine credentials saved on {join(PROVIDENTIA_ROOT, '.env')}\n")
                    
    def select_files_to_download(self, nc_files_to_download):
        """ Returns the files that are not already downloaded. """
        # initialise list of non-downloaded files
        not_downloaded_files = []
        if nc_files_to_download:
            # get the downloaded and not downloaded files
            not_downloaded_files = list(filter(lambda x:not os.path.exists(x), nc_files_to_download))
            downloaded_files = list(filter(lambda x:os.path.exists(x), nc_files_to_download))
            
            # get the files that were downloaded before the execution
            downloaded_before_execution_files = list(filter(lambda x:self.prov_start_time > os.path.getctime(x), downloaded_files))

            # if there was any file downloaded before the execution    
            if downloaded_before_execution_files:
                # make the user choose between overwriting or not overwriting
                if self.overwrite_choice not in ['y','n']:
                    # ask if user wants to overwrite
                    while self.overwrite_choice not in ['y','n']:
                        self.overwrite_choice = input("\nThere are some files that were already downloaded in a previous download, do you want to overwrite them (y/n)? ").lower() 
                    # ask if user wants to remember the decision
                    remind_txt = None
                    while remind_txt not in ['y','n']:
                        remind_txt = input("\nDo you want to remember your decision for future downloads (y/n)? ").lower() 
                    # save the decision
                    if remind_txt == 'y':
                        with open(join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                            f.write(f"OVERWRITE={self.overwrite_choice}\n")
                # if user wants to overwrite then add the the files that were 
                # downloaded before the execution as if they were never download
                if self.overwrite_choice == 'y':
                    not_downloaded_files += downloaded_before_execution_files
                # change overwritten files boolean to True to indicate that some files were ignored
                else:
                    self.overwritten_files_flag = True

        return not_downloaded_files

    def download_nonghost_network(self, network, initial_check, files_to_download=None):
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect() 
        
        if not initial_check:
            # print current_network
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nDownloading non-GHOST {network} network data from {REMOTE_MACHINE}...")

        # if not valid network, check if user put the network on init_prov 
        # TODO Move to configuration.py
        if network not in self.nonghost_available_networks:
            msg = f"The {network} network could not be found on {join(PROVIDENTIA_ROOT,'settings','init_prov.yaml')} nonghost_available_networks list."
            msg += "\nPlease, add the network to the list and execute again."
            show_message(self, msg, deactivate=initial_check)
            return
        
        # check if nonghost network exists in directory
        # TODO: Change this to somewhere in configuration, the one up too
        try:
            self.sftp.stat(join(self.nonghost_remote_obs_path,network))
        except FileNotFoundError:
            msg = f"There is no data available in {REMOTE_MACHINE} for {network} network."
            show_message(self, msg, deactivate=initial_check)
            return

        # check if all resolutions are in init_prov, if not warning and delete the not correct ones
        # TODO move to configuration.py
        not_available_resolutions = set(self.resolution) - set(self.nonghost_available_resolutions)
        if not_available_resolutions:
            available_resolutions = set(self.resolution) - not_available_resolutions
            msg = f"The resolution/s {', '.join(available_resolutions)} could not be found on {join(PROVIDENTIA_ROOT,'settings','init_prov.yaml')} nonghost_available_resolutions list."
            msg += "\nPlease, add the necessary resolutions to the list and execute again."
            show_message(self, msg, deactivate=initial_check)
            return

        # get resolution and species combinations
        res_spec_dir = []

        sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(join(self.nonghost_remote_obs_path,network))).intersection(self.nonghost_available_resolutions)
        for resolution in sftp_resolutions:
            try:
                sftp_species = self.species if self.species else set(self.sftp.listdir(join(self.nonghost_remote_obs_path,network,resolution))).intersection(self.available_species)
            except FileNotFoundError:
                msg = f"There is no data available in {REMOTE_MACHINE} for {network} network at {resolution} resolution"
                show_message(self, msg, deactivate=initial_check)
                continue
            for species in sftp_species: 
                res_spec_dir.append(join(self.nonghost_remote_obs_path,network,resolution,species))
        
        if res_spec_dir:
            
            # initialise list with all the nc files to be downloaded
            initial_check_nc_files = []

            if not initial_check:
                self.logger.info(f"\n{network} observations to download ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range within the specie and resolution combination
            for remote_dir in res_spec_dir:

                local_dir = join(self.nonghost_root,remote_dir.split('/',6)[-1])
                species = remote_dir.split('/')[-1]
                resolution = remote_dir.split('/')[-2]

                #  print the species, resolution and network combinations that are going to be downloaded
                if not initial_check:
                    self.logger.info(f"\n  - {local_dir}, source: {remote_dir} ({REMOTE_MACHINE})")

                try:
                    nc_files = self.sftp.listdir(remote_dir)
                except FileNotFoundError:
                    msg = f"There is no data available in {REMOTE_MACHINE} for {network} network for {species} species at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue
                
                # get the nc files in the date range
                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if network + species + resolution + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available in {REMOTE_MACHINE} from {self.start_date} to {self.end_date} for {network} network {species} species at {resolution} resolution."
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # download the valid resolution specie date combinations
                else:
        
                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # sort nc_files
                    valid_nc_files.sort() 

                    if not initial_check:
                        # get the ones that are not already downloaded
                        valid_nc_files = list(filter(lambda x:join(local_dir,x) in files_to_download, valid_nc_files))
                        if not valid_nc_files:
                            msg = "Files were already downloaded."
                            show_message(self, msg, deactivate=initial_check)     
                            continue  

                    if not initial_check and not self.logfile:
                        # print the tqdm bar if output goes to screen       
                        valid_nc_files_iter = tqdm(valid_nc_files,bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                    else:
                        # do not print the bar if it is the initial check
                        valid_nc_files_iter = valid_nc_files

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files_iter:
                        local_path = join(local_dir,nc_file)
                        if initial_check:
                            initial_check_nc_files.append(local_path)
                        else:
                            # get last downloaded file in case there was a keyboard interrupt
                            self.latest_nc_file_path = local_path

                            # initialize the timeout and get the file
                            self.ncfile_dl_start_time = time.time()
                            remote_path = join(remote_dir,nc_file)
                            self.sftp.get(remote_path, local_path, callback=self.check_time)

            return initial_check_nc_files
        
    def download_ghost_network_sftp(self, network, initial_check, files_to_download=None):
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect() 

        if not initial_check:
            # print current_network
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nDownloading GHOST {network} network data from {REMOTE_MACHINE}...")

        # if not valid network, next
        if network not in self.sftp.listdir(self.ghost_remote_obs_path):
            msg = f"There is no data available in {REMOTE_MACHINE} for {network} network."
            show_message(self, msg, deactivate=initial_check)
            return 
        
        # if not valid combination of GHOST version and network, next 
        elif self.ghost_version not in self.sftp.listdir(join(self.ghost_remote_obs_path,network)):
            msg = f"There is no data available in {REMOTE_MACHINE} for {network} network for the current ghost version ({self.ghost_version})."
            
            available_ghost_versions = set(self.sftp.listdir(join(self.ghost_remote_obs_path,network))).intersection(self.possible_ghost_versions)

            # list that saves the GHOST versions with valid nc files
            valid_available_ghost_versions = []
            
            # check for combinations of species, resolutions, network, and day in the available versions
            if available_ghost_versions:
                # iterate the different GHOST versions
                for possible_ghost_version in available_ghost_versions:
                    remote_dir_ghost_version = join(self.ghost_remote_obs_path, network, possible_ghost_version)
                    
                    # iterate the different resolutions
                    sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir_ghost_version)).intersection(self.ghost_available_resolutions)
                    for resolution in sftp_resolutions:
                        try:
                            species_list = self.species if self.species else self.sftp.listdir(join(remote_dir_ghost_version, resolution))
                        except FileNotFoundError:
                            continue

                        # iterate the different species
                        for species in species_list:
                            # look for valid nc files in the date range
                            try:
                                nc_files = self.sftp.listdir(join(remote_dir_ghost_version, resolution, species))
                                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)
                                if valid_nc_files:
                                    valid_available_ghost_versions.append(possible_ghost_version)
                                    break
                            except FileNotFoundError:
                                continue

            if valid_available_ghost_versions:
                msg += f" Please check one of the available versions: {', '.join(sorted(valid_available_ghost_versions))}"
            elif available_ghost_versions:
                msg += " There are no other versions available at the moment with this configuration."
            else:
                msg += " There are no other versions available at the moment."

            show_message(self, msg, deactivate=initial_check)
            return

        # get resolution and species combinations
        res_spec_dir = []

        remote_dir = join(self.ghost_remote_obs_path,network,self.ghost_version)

        sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir)).intersection(self.ghost_available_resolutions)
        for resolution in sftp_resolutions:
            try:
                sftp_species = self.species if self.species else set(self.sftp.listdir(join(remote_dir,resolution))).intersection(self.available_species)
            except FileNotFoundError:
                msg = f"There is no data available in {REMOTE_MACHINE} for {network} network at {resolution} resolution"
                show_message(self, msg, deactivate=initial_check)
                continue
            for species in sftp_species: 
                res_spec_dir.append(join(remote_dir,resolution,species))
        
        # print the species, resolution and network combinations that are going to be downloaded
        if res_spec_dir:            
            
            # initialise list with all the nc files to be downloaded
            initial_check_nc_files = []
            
            if not initial_check:
                self.logger.info(f"\n{network} observations to download ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range within the specie and resolution combination
            for remote_dir in res_spec_dir:

                local_dir = join(self.ghost_root,remote_dir.split('/',7)[-1])
                species = remote_dir.split('/')[-1]
                resolution = remote_dir.split('/')[-2]

                #  print the species, resolution and network combinations that are going to be downloaded
                if not initial_check:
                    self.logger.info(f"\n  - {local_dir}, source: {remote_dir} ({REMOTE_MACHINE})")

                try:
                    nc_files = self.sftp.listdir(remote_dir)
                except FileNotFoundError:
                    msg = f"There is no data available in {REMOTE_MACHINE} for {network} network for {species} species at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue
                
                # get the nc files in the date range
                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if network + species + resolution + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available in {REMOTE_MACHINE} from {self.start_date} to {self.end_date} for {network} network {species} species at {resolution} resolution."
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # download the valid resolution specie date combinations
                else:

                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # sort nc_files
                    valid_nc_files.sort() 

                    if not initial_check:
                        # get the ones that are not already downloaded
                        valid_nc_files = list(filter(lambda x:join(local_dir,x) in files_to_download, valid_nc_files))
                        if not valid_nc_files:
                            msg = "Files were already downloaded."
                            show_message(self, msg, deactivate=initial_check)     
                            continue  

                    if not initial_check and not self.logfile:
                        # print the tqdm bar if output goes to screen          
                        valid_nc_files_iter = tqdm(valid_nc_files,bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                    else:
                        # do not print the bar if it is the initial check
                        valid_nc_files_iter = valid_nc_files

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files_iter:
                        local_path = join(local_dir,nc_file)
                        if initial_check:
                            initial_check_nc_files.append(local_path)
                        else:
                            # get last downloaded file in case there was a keyboard interrupt
                            self.latest_nc_file_path = local_path

                            # initialize the timeout and get the file
                            self.ncfile_dl_start_time = time.time()
                            remote_path = join(remote_dir,nc_file)
                            self.sftp.get(remote_path, local_path, callback=self.check_time)
                       
            return initial_check_nc_files

    def download_ghost_network_zenodo(self, network, initial_check, files_to_download=None):
        # import remotezip
        from remotezip import RemoteZip
        
        if not initial_check:
            # print current_network
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nDownloading GHOST {network} network data from Zenodo...")

        # if first time reading a GHOST network, get current zips urls in zenodo page
        if not hasattr(self,"zenodo_ghost_available_networks"): 
            self.fetch_zenodo_networks()

        # if not valid network, next
        if network not in self.zenodo_ghost_available_networks:
            msg = f"There is no data available in Zenodo for {network} network."
            show_message(self, msg, deactivate=initial_check)
            return

        # get url to download the zip file for the current network
        zip_file_urls = self.zenodo_ghost_available_networks[network]

        # get resolution and/or species combinations
        # network
        if not self.resolution and not self.species:
            res_spec_combinations = [""]
        # network + resolution 
        elif not self.species:
            res_spec_combinations = [f"/{resolution}/" for resolution in self.resolution]
        # network + species 
        elif not self.resolution:
            res_spec_combinations = [f"/{species}.tar.xz"  for species in self.species]   
        # network + resolution + species 
        else:
            res_spec_combinations = [f"/{resolution}/{species}.tar.xz" for species in self.species for resolution in self.resolution]

        # get all the species tar files which fulfill the combination condition
        res_spec_dir_tail = []
        with RemoteZip(zip_file_urls) as zip:
            for combi in res_spec_combinations:
                res_spec_dir_tail += list(filter(lambda x: combi in x, zip.namelist()))
            res_spec_dir_tail = list(filter(lambda x: x[-7:] == '.tar.xz', res_spec_dir_tail))

            # warning if network + species + resolution combination is gets no matching results
            if not res_spec_dir_tail:
                print_spec = f'{",".join(self.species)} species' if self.species else ""
                print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
                msg = f"There is no data available in Zenodo for {network} network {print_spec} {print_res}"
                show_message(self, msg, deactivate=initial_check)

            # check if there's any possible combination with user's network, resolution and species
            else:
                # initialise list with all the nc files to be downloaded
                initial_check_nc_files = []

                # print the species, resolution and network combinations that are going to be downloaded
                if not initial_check:
                    self.logger.info(f"\n{network} observations to download:")
                
                # initialise only possible GHOST version
                # TODO in the future there will be various versions in zenodo
                last_ghost_version = '1.5' 
                
                for remote_dir_tail in res_spec_dir_tail:
                    resolution, specie = remote_dir_tail.split("/")[1:]
                    specie = specie[:-7]
                    local_dir = join(self.ghost_root,network,last_ghost_version,resolution,specie)

                    #  print the species, resolution and network combinations that are going to be downloaded
                    if not initial_check:
                        self.logger.info(f"\n  - {local_dir}")

                    # create temporal dir to store the middle tar file with its directories
                    temp_dir = join(self.ghost_root,'.temp')
                    if not os.path.exists(temp_dir):
                        os.mkdir(temp_dir)

                    zip.extract(remote_dir_tail,temp_dir)
                    
                    # get path and the name of the directory of the tar file
                    tar_path = join(self.ghost_root, network, str(last_ghost_version), *remote_dir_tail.split("/")[1:])
                    temp_path = join(temp_dir,remote_dir_tail)

                    # extract nc file from tar file
                    with tarfile.open(temp_path) as tar_file:
                        # get the nc files that are between the start and end date
                        tar_names = tar_file.getnames()
                        valid_nc_file_names = self.get_valid_nc_files_in_date_range(tar_names)
                        
                        # warning if network + species + resolution + date range combination gets no matching results                        
                        if not valid_nc_file_names:
                            print_spec = f'{",".join(self.species)} species' if self.species else ""
                            print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
                            msg = f"There is no data available in Zenodo from {self.start_date} to {self.end_date} for {network} network {print_spec} {print_res}"
                            show_message(self, msg, deactivate=initial_check)
                            continue
                        # download the valid resolution specie date combinations
                        else:                 
                            # create directories if they don't exist
                            tar_dir = os.path.dirname(tar_path) 
                            if not os.path.exists(tar_dir):
                                os.makedirs(tar_dir)
                            
                            if not initial_check:
                                # get the ones that are not already downloaded
                                valid_nc_file_names = list(filter(lambda x:join(local_dir,x.split('/')[-1]) in files_to_download, valid_nc_file_names))
                                if not valid_nc_file_names:
                                    msg = "Files were already downloaded."
                                    show_message(self, msg, deactivate=initial_check)     
                                    continue  
                                
                            valid_nc_files = list(filter(lambda x:x.name in valid_nc_file_names, tar_file.getmembers()))
                          
                            # sort nc_files
                            valid_nc_files.sort(key = lambda x:x.name)   

                            if not initial_check and not self.logfile:
                                # print the tqdm bar if output goes to screen   
                                valid_nc_files_iter = tqdm(valid_nc_files,bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                            else:
                                # do not print the bar if it is the initial check
                                valid_nc_files_iter = valid_nc_files

                            for nc_file in valid_nc_files_iter:
                                local_path = join(tar_dir,nc_file.name)
                                if initial_check:
                                    initial_check_nc_files.append(local_path)
                                else:
                                    # get last downloaded file in case there was a keyboard interrupt
                                    self.latest_nc_file_path = local_path

                                    # extract the file
                                    tar_file.extract(member = nc_file, path = tar_dir)
                
                # remove the temp directory
                shutil.rmtree(os.path.dirname(os.path.dirname(temp_path)))
                
                return initial_check_nc_files             
                                
    def download_experiment(self, experiment, initial_check, files_to_download=None):
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect()  
        
        if not initial_check:
            # print current experiment
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nDownloading {experiment} experiment data from {REMOTE_MACHINE}...")
            
        # get resolution and species combinations
        res_spec_dir = []
        
        # domain and ensemble option are part of the experiment name, all united by dash (-)
        experiment_new = experiment

        # domain and ensemble option are directories
        experiment_old = experiment.replace("-","/")
        
        # get remote directory format depending on the GHOST version
        experiment = experiment_old if self.ghost_version in ["1.2", "1.3", "1.3.1"] else experiment_new

        # get remote directory
        remote_dir = join(self.exp_remote_path,self.ghost_version,experiment)

        # check if experiment exists
        try:
            self.sftp.stat(remote_dir)
        except FileNotFoundError:
            msg = f"There is no data available in {REMOTE_MACHINE} for {experiment_new} experiment for the current GHOST version ({self.ghost_version})."

            # get possible GHOST versions from the combination of GHOST_standards and the real avaibles in the experiment remote machine path
            possible_ghost_versions = set(self.sftp.listdir(self.exp_remote_path)).intersection(set(self.possible_ghost_versions))
            
            # get available experiments in other GHOST versions (considering different formats)
            available_ghost_versions = []

            for possible_ghost_version in possible_ghost_versions:
                try:
                    # get experiment path depending on the GHOST version
                    remote_dir_ghost_version = join(self.exp_remote_path, possible_ghost_version, experiment_old if possible_ghost_version in ["1.2", "1.3", "1.3.1"] else experiment_new)
                    
                    # check if directory exists
                    self.sftp.stat(remote_dir_ghost_version)

                    # if it doesn't break, the experiment exists in this version
                    available_ghost_versions.append(possible_ghost_version)

                except FileNotFoundError:
                    continue
                            
            # list that saves the GHOST versions with valid nc files
            valid_available_ghost_versions = []
            
            # check for combinations of species, resolutions, network, and day in the available versions
            if available_ghost_versions:

                # iterate the different GHOST versions
                for possible_ghost_version in available_ghost_versions:
                    remote_dir_ghost_version = join(self.exp_remote_path, possible_ghost_version, experiment_old if possible_ghost_version in ["1.2", "1.3", "1.3.1"] else experiment_new)
                    
                    # iterate the different resolutions
                    sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir)).intersection(self.nonghost_available_resolutions)
                    for resolution in sftp_resolutions:                        
                        try:
                            species_list = self.species if self.species else self.sftp.listdir(join(remote_dir_ghost_version, resolution))
                        except FileNotFoundError:
                            continue
                        
                        # iterate the different species
                        for species in species_list:
                            try:
                                network_list = self.network if self.network else self.sftp.listdir(join(remote_dir_ghost_version, resolution, species))
                            except FileNotFoundError:
                                continue
                            
                            # iterate the different networks
                            for network in network_list:
                                # look for valid nc files in the date range
                                try:
                                    nc_files = self.sftp.listdir(join(remote_dir_ghost_version, resolution, species, network))
                                    valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)
                                    if valid_nc_files:
                                        valid_available_ghost_versions.append(possible_ghost_version)
                                        break
                                except FileNotFoundError:
                                    continue
                            
            if valid_available_ghost_versions and check_for_ghost(network):
                msg += f" Please check one of the available versions: {', '.join(sorted(valid_available_ghost_versions))}"
            elif available_ghost_versions:
                msg += " There are no other versions available at the moment with this configuration."
            else:
                msg += " There are no other versions available at the moment."

            show_message(self, msg, deactivate=initial_check)
            return

        sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir)).intersection(self.nonghost_available_resolutions)
        for resolution in sftp_resolutions:
            try:
                sftp_species = self.species if self.species else set(self.sftp.listdir(join(remote_dir,resolution))).intersection(self.available_species)
            except FileNotFoundError:
                msg = f"There is no data available in {REMOTE_MACHINE} for {experiment_new} experiment at {resolution} resolution"
                show_message(self, msg, deactivate=initial_check)
                continue
            for species in sftp_species: 
                try:
                    sftp_network = self.network if self.network else self.sftp.listdir(join(remote_dir,resolution,species))
                except FileNotFoundError:
                    msg = f"There is no data available in {REMOTE_MACHINE} for {experiment_new} experiment for {species} species at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue
                for network in sftp_network:
                    # if network is nonghost, change the slashes to dashes
                    if not check_for_ghost(network):
                        network = network.replace("/", "-")
                    res_spec_dir.append(join(remote_dir,resolution,species,network))
        
        # print the species, resolution and experiment combinations that are going to be downloaded
        if res_spec_dir:

            # initialise list with all the nc files to be downloaded
            initial_check_nc_files = []

            if not initial_check:
                self.logger.info(f"\n{experiment_new} experiment data to download ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range
            for remote_dir in res_spec_dir:
                if not initial_check:
                    local_path = remote_dir.split('/',7)[-1]
                    if self.ghost_version in ["1.2", "1.3", "1.3.1"]:
                        self.logger.info(f"\n  - {join(self.exp_root,self.ghost_version,'-'.join(local_path.split('/')[1:4]),*local_path.split('/')[4:])}, source: {remote_dir} ({REMOTE_MACHINE})")
                    else:
                        self.logger.info(f"\n  - {join(self.exp_root,local_path)}, source: {remote_dir} ({REMOTE_MACHINE})")
            
                network = remote_dir.split('/')[-1]
                species = remote_dir.split('/')[-2]
                resolution = remote_dir.split('/')[-3]
                
                # get nc files if directory is found
                try:
                    nc_files = self.sftp.listdir(remote_dir)
                except FileNotFoundError:
                    msg = f"There is no data available in {REMOTE_MACHINE} for {experiment_new} experiment for {species} species {network} network at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # get the nc files in the date range       
                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if experiment + species + resolution + network + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available in {REMOTE_MACHINE} from {self.start_date} to {self.end_date} for {experiment_new} experiment {species} species {network} network at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # download the valid resolution specie date combinations
                else:
                    # create local directory (always with experiments on the new format)
                    local_dir = join(self.exp_root,self.ghost_version,experiment_new,resolution,species,network)
                    
                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # sort nc_files
                    valid_nc_files.sort() 

                    if not initial_check:
                        # get the ones that are not already downloaded
                        valid_nc_files = list(filter(lambda x:join(local_dir,x) in files_to_download, valid_nc_files))
                        if not valid_nc_files:
                            msg = "Files were already downloaded."
                            show_message(self, msg, deactivate=initial_check)     
                            continue

                    if not initial_check and not self.logfile:
                        # print the tqdm bar if output goes to screen            
                        valid_nc_files_iter = tqdm(valid_nc_files, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                    else:
                        # do not print the bar if it is the initial check
                        valid_nc_files_iter = valid_nc_files

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files_iter:
                        local_path = join(local_dir,nc_file)
                        if initial_check:
                            initial_check_nc_files.append(local_path)
                        else:
                            # get last downloaded file in case there was a keyboard interrupt
                            self.latest_nc_file_path = local_path

                            # initialize the timeout and get the file
                            self.ncfile_dl_start_time = time.time()
                            remote_path = join(remote_dir,nc_file)
                            self.sftp.get(remote_path, local_path, callback=self.check_time) 
            
            return initial_check_nc_files

    def download_non_interpolated_experiment(self, experiment, initial_check, files_to_download=None):
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect()  
        
        if not initial_check:
            # print current experiment
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nDownloading {experiment} non-interpolated experiment data from {REMOTE_MACHINE}...")
            
        # get resolution and species combinations
        res_spec_dir = []

        # get experiment id and the domain
        exp_id, domain, ensemble_options = experiment.split("-")

        # initialise warning message and experiment exists boolean
        msg = ""
        experiment_exists = False

        # see if the experiment is any of the interp_experiment.yaml lists
        for experiment_type, experiment_dict in interp_experiments.items():
            if exp_id in experiment_dict["experiments"]:
                experiment_exists = True
                break
        
        # if it is in the list, check if the paths work
        if experiment_exists is True:
            # get boolean to False again until the paths works
            experiment_exists = False

            # get all paths that work
            # if there is none, show a warning
            exp_dir_functional_list = []    
            for exp_dir in experiment_dict["paths"]:
                # esarchive in transfer5 is located inside gpfs
                if "/esarchive/" == exp_dir[:11]:
                    exp_dir = join("/gpfs/archive/bsc32/",exp_dir[1:])
                # check if directory exists in the remote machine
                try:
                    self.sftp.stat(exp_dir)
                    exp_dir_functional_list.append(exp_dir)      
                except FileNotFoundError:
                    pass

            # if none of the paths are in this current machine, break
            if not exp_dir_functional_list:
                msg += f"None of the paths specified in {join('settings', 'interp_experiments.yaml')} are available on the remote machine ({REMOTE_MACHINE}). "
            # if any path works, get the first one that has the experiment
            else:
                # get first functional directory  
                for exp_dir in exp_dir_functional_list:
                    remote_dir = join(exp_dir,exp_id,domain)
                    # check if remote experiment and domain directories exist in the remote machine
                    try:
                        self.sftp.stat(remote_dir)
                        experiment_exists = True
                        break
                    except FileNotFoundError:
                        pass

                # if the experiment-domain combination is not possible, show the warning
                if experiment_exists is False:
                    msg += f"There is no data available for the {exp_id} experiment with the {domain} domain in none of the paths specified in {join('settings', 'interp_experiments.yaml')} in the remote machine ({REMOTE_MACHINE}). "

        # if experiment was not in the list, or any of the paths were available
        # or there was no valid path experiment combination then search in the gpfs directory
        if experiment_exists is False:
            # get all possible experiments
            exp_to_interp_path = join(self.exp_to_interp_remote_path,exp_id,domain)
            try:
                self.sftp.stat(exp_to_interp_path)
                remote_dir = exp_to_interp_path
                experiment_exists = True
            except FileNotFoundError:
                pass 
            
            # add to the message if experiment was not found in the gpfs remote directory
            msg += f"Cannot find the {exp_id} experiment with the {domain} domain in '{self.exp_to_interp_remote_path}'."    
        
        # if the experiment-domain combination is not possible, break
        if experiment_exists is False:
            show_message(self, msg, deactivate=initial_check)
            return

        # get all the resolutions available in the remote directory
        sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir)).intersection(self.nonghost_available_resolutions)

        # iterate through the resolutions
        for resolution in sftp_resolutions:
            try:
                # get available species ("normal" and mapped)
                available_species = self.available_species+[spec[0] for spec in mapping_species.values()]
                sftp_species = self.species if self.species else set(self.sftp.listdir(join(remote_dir,resolution))).intersection(available_species)
            except FileNotFoundError:
                msg = f"There is no data available in {REMOTE_MACHINE} for the {exp_id} experiment with the {domain} domain at {resolution} resolution"
                show_message(self, msg, deactivate=initial_check)
                continue

            # iterate through the species
            for speci_to_process in sftp_species: 
                # initialize boolean that saves whether species was found
                species_exists = False
                species = speci_to_process
                # first try with the original species
                try:
                    # if it is an ensemble member
                    if ensemble_options.isdigit() or ensemble_options == 'allmembers':
                        res_spec = join(remote_dir,resolution,species)
                    # if it is an ensemble statistic
                    else:
                        res_spec = join(remote_dir,resolution,"ensemble-stats",species+"_"+ensemble_options)
  
                    self.sftp.stat(res_spec)
                    species_exists = True
                # if there are none, try with the mapped species
                except FileNotFoundError:
                    # change species name to the species to map
                    if speci_to_process in mapping_species:
                        for mapping_speci in mapping_species[speci_to_process]:
                            try:
                                # if it is an ensemble member
                                if ensemble_options.isdigit() or ensemble_options == 'allmembers':
                                    res_spec = join(remote_dir,resolution, mapping_speci)
                                # if it is an ensemble statistic
                                else:
                                    res_spec = join(remote_dir,resolution, "ensemble-stats", species + "_" + ensemble_options)
  
                                self.sftp.stat(res_spec)  
                                species_exists = True
                                break
                            except FileNotFoundError:
                                pass
                
                # if no species were found, then show the message
                if species_exists is False:
                    msg = f"There is no data available in {REMOTE_MACHINE} for the {exp_id} experiment with the {domain} domain for {species} species at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # add the path with the resolution and species combination to the list
                res_spec_dir.append(res_spec)
                        
        # print the species, resolution and experiment combinations that are going to be downloaded
        if res_spec_dir:

            # initialise list with all the nc files to be downloaded
            initial_check_nc_files = []

            if not initial_check:
                self.logger.info(f"\n{experiment} experiment data to download ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range
            for remote_dir in res_spec_dir:
                if not initial_check:
                    local_path = remote_dir.split('/', 6)[-1]
                    self.logger.info(f"\n  - {join(self.exp_to_interp_root,local_path)}, source: {remote_dir} ({REMOTE_MACHINE})")
                         
                # get nc files
                nc_files = self.sftp.listdir(remote_dir)

                if nc_files:
                    # if it is an ensemble member
                    if ensemble_options.isdigit() or ensemble_options == 'allmembers':
                        # get the domain, resolution and species from the path
                        domain, resolution, species = remote_dir.split('/')[-3:]

                        # identify format of the directory
                        # the format is a tuple of how many - and how many _ are there
                        # the directory format is choosen by popularity
                        formats_list = [(file.count("-"), file.count("_")) for file in nc_files]
                        number_of_formats_dict = {format: formats_list.count(format) for format in set(formats_list)}
                        format = max(number_of_formats_dict, key=number_of_formats_dict.get)
                        
                        # filter and get only the files that follow the format (number of dashes and hyphens and end of file)
                        nc_files = list(filter(lambda x:(x.count("-"),x.count("_")) == format and x.endswith(".nc"),nc_files))
                        
                        # example: od550du_2019040212.nc (0,1)
                        if format == (0,1):
                            # when there is no ensemble option in the name only allmembers and 000 are valid
                            if ensemble_options == '000' or ensemble_options == 'allmembers':
                                nc_files = list(filter(lambda x:x.split("_")[0] == species, nc_files))

                        # example: od550du-000_2021020812.nc (1,1)
                        elif format == (1,1):
                            # filter by ensemble option in case that ensemble option is not allmembers
                            if ensemble_options != 'allmembers':
                                nc_files = list(filter(lambda x:x.split("_")[0] == species+'-'+ensemble_options,nc_files))
                           
                        else:
                            # TODO delete this in the future
                            error = "It is not possible to download this nc file type yet. Please, contact the developers.", nc_files
                            self.logger.error(error)
                            sys.exit(1)
                    
                    # if it is an ensemble statistic
                    else:
                        # get the domain, resolution and species from the path
                        domain, resolution, _, species = remote_dir.split('/')[-4:]
                        species = species.split("_",1)[0]

                        # filter the nc files to only get the ones that have the correct species and stats
                        nc_files = list(filter(lambda x:x.split("_")[0] == species and "_".join(x[:-3].split("_")[2:]) == ensemble_options, nc_files))
                
                # if there is no options with the ensemble option, tell the user
                if nc_files == []:
                    msg = f"There is no data available in {REMOTE_MACHINE} for the {exp_id} experiment with the {domain} domain with the {ensemble_options} ensemble option."
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # get the nc files in the date range        
                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if experiment + species + resolution + network + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available in {REMOTE_MACHINE} from {self.start_date} to {self.end_date} for {experiment} experiment {species} species at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # download the valid resolution specie date combinations
                else:
                    # create local directory 
                    # if it is an ensemble member
                    if ensemble_options.isdigit() or ensemble_options == 'allmembers':
                        local_dir = join(self.exp_to_interp_root,exp_id,domain,resolution,species)
                    else:
                        local_dir = join(self.exp_to_interp_root,exp_id,domain,resolution,"ensemble-stats",species+"_"+ensemble_options)
                    
                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # sort nc_files
                    valid_nc_files.sort() 

                    if not initial_check:
                        # get the ones that are not already downloaded
                        valid_nc_files = list(filter(lambda x:join(local_dir,x) in files_to_download, valid_nc_files))
                        if not valid_nc_files:
                            msg = "Files were already downloaded."
                            show_message(self, msg, deactivate=initial_check)     
                            continue    

                    if not initial_check and not self.logfile:
                        # print the tqdm bar if output goes to screen        
                        valid_nc_files_iter = tqdm(valid_nc_files, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                    else:
                        # do not print the bar if it is the initial check
                        valid_nc_files_iter = valid_nc_files

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files_iter:
                        local_path = join(local_dir,nc_file)
                        if initial_check:
                            initial_check_nc_files.append(local_path)
                        else:
                            # get last downloaded file in case there was a keyboard interrupt
                            self.latest_nc_file_path = local_path

                            # initialize the timeout and get the file
                            self.ncfile_dl_start_time = time.time()
                            remote_path = join(remote_dir, nc_file)
                            self.sftp.get(remote_path, local_path, callback=self.check_time) 
            
            return initial_check_nc_files

        # tell the user if not valid resolution specie date combinations
        else:
            msg = "There are no available observations to be downloaded."
            show_message(self, msg, deactivate=initial_check)
            
    def copy_non_interpolated_experiment(self, experiment, initial_check, files_to_download=None):
        if not initial_check:
            # print current experiment
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nCopying {experiment} non-interpolated experiment data from esarchive to gpfs in {self.machine}...")
            
        # get resolution and species combinations
        res_spec_dir = []

        # get experiment id and the domain
        exp_id, domain, ensemble_options = experiment.split("-")

        # get experiment type
        for experiment_type, experiment_dict in interp_experiments.items():
            if exp_id in experiment_dict["experiments"]:
                break
        
        # get experiment specific directories list
        exp_dir_list = experiment_dict["paths"]

        # take all functional directories
        exp_dir_functional_list = []    
        for exp_dir in exp_dir_list:
            # make sure that it comes from esarchive
            if "/esarchive/" in exp_dir:
                # esarchive in transfer5 is located inside gpfs
                if "/esarchive/" == exp_dir[:11] and self.machine == "storage5":
                    exp_dir = join("/gpfs/archive/bsc32/",exp_dir[1:])
                # check if directory exists in esarchive
                if os.path.exists(exp_dir):
                    exp_dir_functional_list.append(exp_dir)     
            
        # if none of the paths are in this current machine, break
        if not exp_dir_functional_list:
            msg = f"None of the paths specified in {join('settings', 'interp_experiments.yaml')} are available on esarchive."
            show_message(self, msg, deactivate=initial_check)
            return
        
        # take first functional directory  
        esarchive_dir = None
        for exp_dir in exp_dir_functional_list:
            temp_esarchive_dir = join(exp_dir,exp_id,domain)
            # check if experiment and domain directories exist in esarchive machine
            if os.path.exists(temp_esarchive_dir): 
                esarchive_dir = temp_esarchive_dir
                break
        
        # if the experiment-domain combination is not possible, break
        if esarchive_dir is None:
            msg = f"There is no data available for the {exp_id} experiment with the {domain} domain in none of the paths specified in {join('settings', 'interp_experiments.yaml')} in esarchive."
            show_message(self, msg, deactivate=initial_check)
            return

        # get all the resolutions available in the esarchive directory
        sftp_resolutions = self.resolution if self.resolution else set(os.listdir(esarchive_dir)).intersection(self.nonghost_available_resolutions)

        # iterate through the resolutions
        for resolution in sftp_resolutions:
            try:
                # get available species ("normal" and mapped)
                available_species = self.available_species+[spec[0] for spec in mapping_species.values()]
                sftp_species = self.species if self.species else set(os.listdir(join(esarchive_dir,resolution))).intersection(available_species)
            except FileNotFoundError:
                msg = f"There is no data available in esarchive for the {exp_id} experiment with the {domain} domain at {resolution} resolution"
                show_message(self, msg, deactivate=initial_check)
                continue

            # iterate through the species
            for speci_to_process in sftp_species: 
                # initialize boolean that saves whether species was found
                species_exists = False
                species = speci_to_process

                # if it is an ensemble member
                if ensemble_options.isdigit() or ensemble_options == 'allmembers':
                    res_spec = join(esarchive_dir,resolution,species)
                # if it is an ensemble statistic
                else:
                    res_spec = join(esarchive_dir,resolution,"ensemble-stats",species+"_"+ensemble_options)
                species_exists = os.path.exists(res_spec)
                # if there are none, try with the mapped species
                if species_exists is False:
                    # change species name to the species to map
                    if speci_to_process in mapping_species:
                        for species in mapping_species[speci_to_process]:
                            # if it is an ensemble member
                            if ensemble_options.isdigit() or ensemble_options == 'allmembers':
                                res_spec = join(esarchive_dir,resolution,species)
                            # if it is an ensemble statistic
                            else:
                                res_spec = join(esarchive_dir,resolution,"ensemble-stats",species+"_"+ensemble_options)
  
                            species_exists = os.path.exists(res_spec)
                
                # if no species were found, then show the message
                if species_exists is False:
                    msg = f"There is no data available in esarchive for the {exp_id} experiment with the {domain} domain for {species} species at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # add the path with the resolution and species combination to the list
                res_spec_dir.append(res_spec)
                        
        # print the species, resolution and experiment combinations that are going to be copied
        if res_spec_dir:

            # initialise list with all the nc files to be copied
            initial_check_nc_files = []

            if not initial_check:
                self.logger.info(f"\n{experiment} experiment data to copy ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range
            for esarchive_dir in res_spec_dir:
                if not initial_check:
                    self.logger.info(f"\n  - {join(self.exp_to_interp_root,'/'.join(esarchive_dir.split('/')[-4:]))}, source: {esarchive_dir} ({self.machine})")
                         
                # get nc files
                nc_files = os.listdir(esarchive_dir)

                if nc_files:
                    # if it is an ensemble member
                    if ensemble_options.isdigit() or ensemble_options == 'allmembers':
                        # get the domain, resolution and species from the path
                        domain, resolution, species = esarchive_dir.split('/')[-3:]

                        # identify format of the directory
                        # the format is a tuple of how many - and how many _ are there
                        # the directory format is choosen by popularity
                        formats_list = [(file.count("-"), file.count("_")) for file in nc_files]
                        number_of_formats_dict = {format: formats_list.count(format) for format in set(formats_list)}
                        format = max(number_of_formats_dict, key=number_of_formats_dict.get)
                        
                        # filter and get only the files that follow the format (number of dashes and hyphens and end of file)
                        nc_files = list(filter(lambda x:(x.count("-"),x.count("_")) == format and x.endswith(".nc"), nc_files))
                        
                        # example: od550du_2019040212.nc (0,1)
                        if format == (0,1):
                            # when there is no ensemble option in the name only allmembers and 000 are valid
                            if ensemble_options == '000' or ensemble_options == 'allmembers':
                                nc_files = list(filter(lambda x:x.split("_")[0] == species, nc_files))
                        
                        # example: od550du-000_2021020812.nc (1,1)
                        elif format == (1,1):
                            # filter by ensemble option in case that ensemble option is not allmembers
                            if ensemble_options != 'allmembers':
                                nc_files = list(filter(lambda x:x.split("_")[0] == species+'-'+ensemble_options,nc_files))
                           
                        else:
                            error = "It is not possible to copy this nc file type yet. Please, contact the developers.", nc_files
                            self.logger.error(error)
                            sys.exit(1)
                    
                    # if it is an ensemble statistic
                    else:
                        # get the domain, resolution and species from the path
                        domain, resolution, _, species = esarchive_dir.split('/')[-4:]
                        species = species.split("_",1)[0]

                        # filter the nc files to only get the ones that have the correct species and stats
                        nc_files = list(filter(lambda x:x.split("_")[0] == species and "_".join(x[:-3].split("_")[2:]) == ensemble_options, nc_files))
                        
                # if there is no options with the ensemble option, tell the user
                if nc_files == []:
                    msg = f"There is no data available in esarchive for the {exp_id} experiment with the {domain} domain with the {ensemble_options} ensemble option."
                    show_message(self, msg, deactivate=initial_check)
                    continue
                
                # get the nc files in the date range        
                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if experiment + species + resolution + network + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available in esarchive from {self.start_date} to {self.end_date} for {experiment} experiment {species} species at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # copy the valid resolution specie date combinations
                else:
                    # if it is an ensemble member
                    if ensemble_options.isdigit() or ensemble_options == 'allmembers':
                        gpfs_dir = join(self.exp_to_interp_root,exp_id,domain,resolution,species)
                    else:
                        gpfs_dir = join(self.exp_to_interp_root,exp_id,domain,resolution,"ensemble-stats",species+"_"+ensemble_options)
                    
                    # create directories if they don't exist
                    if not os.path.exists(gpfs_dir):
                        os.makedirs(gpfs_dir) 
                        # give to each directory 770 permissions and make group owner bsc32 
                        temp_gpfs_dir = gpfs_dir
                        for i in range(4):
                            os.system(f"chmod 770 {temp_gpfs_dir}; chgrp bsc32 {temp_gpfs_dir}")
                            temp_gpfs_dir = os.path.dirname(temp_gpfs_dir)

                    # sort nc_files
                    valid_nc_files.sort() 

                    if not initial_check:
                        # get the ones that are not already copied
                        valid_nc_files = list(filter(lambda x:join(gpfs_dir,x) in files_to_download, valid_nc_files))
                        if not valid_nc_files:
                            msg = "Files were already copied."
                            show_message(self, msg, deactivate=initial_check)     
                            continue   

                    if not initial_check and not self.logfile:
                        # print the tqdm bar if output goes to screen         
                        valid_nc_files_iter = tqdm(valid_nc_files, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Copying files from esarchive to gpfs ({len(valid_nc_files)})")
                    else:
                        # do not print the bar if it is the initial check
                        valid_nc_files_iter = valid_nc_files

                    # copy each individual nc file using sftp protocol
                    for nc_file in valid_nc_files_iter:
                        gpfs_path = join(gpfs_dir,nc_file)
                        
                        if initial_check:
                            initial_check_nc_files.append(gpfs_path)
                        
                        else:
                            # get last downloaded file in case there was a keyboard interrupt
                            self.latest_nc_file_path = gpfs_path
                            
                            # TODO: remove if the permissions are not needed
                            # check if the file already exists
                            # new_file = not os.path.isfile(gpfs_path)
                           
                            # get rsync command depending on which machine it was ran
                            rsync_command = "dtrsync" if self.machine == "storage5" else "rsync"
                            
                            # copy file
                            esarchive_path = join(esarchive_dir, nc_file)
                            try:
                                with open(os.devnull, 'wb') as devnull:
                                    subprocess.check_call([rsync_command, esarchive_path, gpfs_path], stdout=devnull, stderr=subprocess.STDOUT)
                            except subprocess.CalledProcessError:
                                error = f'Failed to copy the files. Try later.'
                                self.logger.error(error)
                                sys.exit(1)

                            # TODO: fails when creating a new file, wait until users use the copy option to see if it is really needed
                            # give to each new file 770 permissions to directory and make group owner bsc32 
                            # if new_file:                        
                            #     os.system(f"chmod 770 {gpfs_path}; chgrp bsc32 {gpfs_path}")

                    # dtrsync generates output files that are generated in 3 seconds 
                    if not initial_check and self.machine == "storage5":
                        time.sleep(3)
            
                        # remove the output files frojm dtrsync
                        for file in os.listdir(PROVIDENTIA_ROOT):
                            if file.startswith("dtrsync_"):
                                os.remove(join(PROVIDENTIA_ROOT,file)) 
            
            return initial_check_nc_files

        # tell the user if not valid resolution specie date combinations
        else:
            msg = "There are no available observations to be copied."
            show_message(self, msg, deactivate=initial_check)
            
    def fetch_zenodo_networks(self):
        # Get urls from zenodo to get GHOST zip files url

        # initialize dictionary to store possible networks
        self.zenodo_ghost_available_networks = {} 

        response = requests.get(self.ghost_url)

        # Check if the request was successful (status code 200)
        if response.status_code != 200:
            error = f'Failed to retrieve the webpage. Status code: {response.status_code}'
            self.logger.error(error)
            sys.exit(1)
        
        # fill network dictionary with its corresponding zip url
        for line in response.text.split(">"):
            if '<link rel="alternate" type="application/zip" href=' in line:
                zip_file_url = line.split('href="')[-1][:-1]
                zip_network = line.split("/")[-1][:-5]
                self.zenodo_ghost_available_networks[zip_network] = zip_file_url

    def fetch_cams_dates(self, url):
        
        # send HTTP GET request and get
        response = requests.get(url)

        # get the info for all the fields
        fields_info = re.findall(r'(\{"name":".*?".*?)(?=\,{"name":)', response.text, re.DOTALL)

        # get the date field minimum and maximum
        for field in fields_info:
            if '"name":"date"' in field:
                # transform the str to a dict
                field_dict = json.loads(field)

                # convert to datetime format
                minstart = datetime.strptime(field_dict["details"]["minStart"], '%Y-%m-%d')
                maxend = datetime.strptime(field_dict["details"]["maxEnd"], '%Y-%m-%d')

                # get the mininimum start date and maximum end date
                return minstart, maxend

    def get_all_networks(self):
        # get user input to know which kind of network wants
        download_source = None
        while download_source is None:
            download_source = input("\nDo you want to download GHOST, non-GHOST or all networks? (g/n/a) ").lower()
            download_source = download_source if download_source in ["g","n","a"] else None

        if download_source in ["g","a"]:
            if self.bsc_download_choice == 'y':
                self.network = self.ghost_available_networks
            elif self.bsc_download_choice == 'n':
                if not hasattr(self,"zenodo_ghost_available_networks"): 
                    self.fetch_zenodo_networks()
                    self.network = list(self.zenodo_ghost_available_networks.keys())
            else:
                error = "Download option not valid, check your .env file."
                self.logger.error(error)
                sys.exit(1)

        if download_source in ["n","a"]:
            self.network = self.nonghost_available_networks
    
    def get_all_experiments(self):
        # download all interpolated experiments
        if self.interpolated is True:
            # check if ssh exists and check if still active, connect if not
            if (self.ssh is None) or (self.ssh.get_transport().is_active()):
                self.connect()  

            # get directory content and format it as the experiments       
            experiment_list = self.sftp.listdir(join(self.exp_remote_path,self.ghost_version))
        # download all non interpolated experiments
        else:
            # get all the experiments id
            experiments = []
            for experiment_dict in interp_experiments.values():
                experiments += experiment_dict["experiments"]
            # get all the domain and ensemble options combinations 
            experiment_list = []
            # TODO hardcoded
            for domain in ["ip", "d03", "d01", "regional", "eu", "reg", "ex", "bcn", "cat", "d02", "global","regional_i01", "regional_i02", "regional_i03"]:
                for exp in experiments:
                    experiment_list.append(exp+"-"+domain+"-allmembers")

        self.experiments = dict(zip(experiment_list,experiment_list))

    def get_valid_nc_files_in_date_range(self, nc_files):
        valid_nc_files = []
        for nc_file in sorted(nc_files):
            if ".nc" in nc_file:
                ym = nc_file[:-3].split("_")[1]
                # from yyyymm to yyyymmdd
                if len(ym) == 6:
                    ym = '{}01'.format(ym)
                # from yyyymmddhh to yyyymmdd
                elif len(ym) == 10:
                    ym = ym[:-2]
                # get the date range
                if int(ym) >= int(self.start_date) and int(ym) < int(self.end_date):
                    valid_nc_files.append(nc_file)
                    
        return valid_nc_files        
    
    def download_actris_network(self):
        
        resolution = self.resolution[0]
        target_start_date = datetime(int(self.start_date[:4]), int(self.start_date[4:6]), int(self.start_date[6:8]), 0)
        target_end_date = datetime(int(self.end_date[:4]), int(self.end_date[4:6]), int(self.end_date[6:8]), 23, 59, 59) - timedelta(days=1)

        for var in self.species:
            
            # check if variable name is available
            if var not in parameters_dict.keys():
                self.logger.info(f'Data for {var} cannot be downloaded')
                continue
            else:
                actris_parameter = parameters_dict[var]

            # get files that were already downloaded
            initial_check_nc_files = get_files_to_download(self.nonghost_root, target_start_date, target_end_date, resolution, var)
            files_to_download = self.select_files_to_download(initial_check_nc_files)
            if not files_to_download:
                msg = f"\nFiles were already downloaded for {var} at {resolution} "
                msg += f"resolution between {target_start_date} and {target_end_date}."
                show_message(self, msg, deactivate=False)     
                continue 
            
            # get files info path
            path = get_files_path(var)

            # if file does not exist
            if not os.path.isfile(path):
                # get files information
                self.logger.info(f'\nFile containing information of the files available in Thredds for {var} ({path}) does not exist, creating.')
                combined_data = get_files_per_var(var)
                all_files = combined_data[var]['files']
                files_info = get_files_info(all_files, var, path)
                    
            # if file exists
            else:
                # ask if user wants to update file information from NILU Thredds
                if self.origin_update_choice not in ['y','n']:
                    while self.origin_update_choice not in ['y','n']:
                        self.origin_update_choice = input(f"\nFile containing information of the files available in Thredds for {var} ({path}) already exists. Do you want to update it (y/n)? ").lower() 
                    # ask if user wants to remember the decision
                    remind_txt = None
                    while remind_txt not in ['y','n']:
                        remind_txt = input("\nDo you want to remember your decision for future downloads (y/n)? ").lower() 
                    # save the decision
                    if remind_txt == 'y':
                        with open(join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                            f.write(f"ORIGIN_UPDATE={self.origin_update_choice}\n")
                if self.origin_update_choice == 'n':
                    # get files information
                    files_info = yaml.safe_load(open(join(CURRENT_PATH, path)))
                    files_info = {k: v for k, v in files_info.items() if k.strip() and v}
                else:
                    # get files information
                    combined_data = get_files_per_var(var)
                    all_files = combined_data[var]['files']
                    files_info = get_files_info(all_files, var, path)
            
            # go to next variable if no data is found
            if len(files_info) == 0:
                continue

            # filter files by resolution and dates
            self.logger.info('Filtering files by resolution and dates...')
            files = {}
            for file, attributes in files_info.items():
                if attributes["resolution"] == resolution:
                    start_date = datetime.strptime(attributes["start_date"], "%Y-%m-%dT%H:%M:%S UTC")
                    end_date = datetime.strptime(attributes["end_date"], "%Y-%m-%dT%H:%M:%S UTC")
                    for file_to_download in files_to_download:
                        file_to_download_yearmonth = file_to_download.split(f'{var}_')[1].split('.nc')[0]
                        file_to_download_start_date = datetime.strptime(file_to_download_yearmonth, "%Y%m")
                        file_to_download_end_date = datetime(file_to_download_start_date.year, file_to_download_start_date.month, 1) + relativedelta(months=1, seconds=-1)
                        if file_to_download_start_date <= end_date and file_to_download_end_date >= start_date:
                            # from filtered files, save those that are provided multiple times
                            station = attributes["station_reference"]
                            if station not in files:
                                files[station] = []
                            if file not in files[station]:
                                files[station].append(file)

            if len(files) != 0:

                # get data and metadata for each file within period
                start = time.time()
                combined_ds_list, metadata, wavelength = get_data(files, var, actris_parameter, resolution, 
                                                                  target_start_date, target_end_date)
                end = time.time()
                elapsed_minutes = (end - start) / 60
                print(f"Time to read data: {elapsed_minutes:.2f} minutes")

                # check if there is data after reading available files
                if len(combined_ds_list) == 0:
                    self.logger.info('No data were found')
                    continue

                # get flag dimension per station
                N_flag_codes_dims = []
                for ds in combined_ds_list:
                    N_flag_codes_dims.append(ds.dims['N_flag_codes'])
                
                # get maximum number of flags across all stations
                N_flag_codes_max = max(N_flag_codes_dims)
                
                # recreate flag variable so that all stations have the same dimension and can be concatenated, leave nan for unknown values
                combined_ds_list_corrected_flag = []
                for ds in combined_ds_list:
                    flag_data = ds['flag']
                    da_flag = xr.DataArray(
                            np.full((flag_data.sizes['station'], flag_data.sizes['time'], N_flag_codes_max), np.nan),
                            dims=["station", "time", "N_flag_codes"],
                            coords={
                                "time": flag_data.coords["time"],
                            },
                            name="flag"
                        )
                    da_flag[:, :, :flag_data.values.shape[-1]] = flag_data.values
                    ds = ds.drop_vars('flag')
                    ds['flag'] = da_flag
                    combined_ds_list_corrected_flag.append(ds)
            
                # combine and create new dataset
                self.logger.info('Combining files...')
                start = time.time()
                combined_ds = temporally_average_data(combined_ds_list_corrected_flag, resolution, var, self.ghost_version, target_start_date, target_end_date)
                end = time.time()
                elapsed_minutes = (end - start) / 60
                print(f"Time to temporally average: {elapsed_minutes:.2f} minutes")

                # add metadata
                for key, value in metadata[resolution].items():
                    if key in ['latitude', 'longitude']:
                        value = [float(val) for val in value]
                    elif key in ['altitude', 'sampling_height']:
                        value = [float(val.replace('m', '').strip()) if isinstance(val, str) else val for val in value]
                    combined_ds[key] = xr.Variable(data=value, dims=('station'))

                # calculate measurement_altitude if altitude and sampling_height exist
                if ('altitude' in combined_ds.keys()) and ('sampling_height' in combined_ds.keys()):
                    value = combined_ds['altitude'].values + combined_ds['sampling_height'].values
                    combined_ds['measurement_altitude'] = xr.Variable(data=value, dims=('station'))

                # add units for lat and lon
                # TODO: Check attrs geospatial_lat_units and geospatial_lon_units
                combined_ds.latitude.attrs['units'] = 'degrees_north'
                combined_ds.longitude.attrs['units'] = 'degrees_east'

                # add general attrs
                combined_ds.attrs['data_license'] = 'BSD-3-Clause. Copyright 2025 Alba Vilanova Cortezón'
                combined_ds.attrs['source'] = 'Observations'
                combined_ds.attrs['institution'] = 'Barcelona Supercomputing Center'
                combined_ds.attrs['creator_name'] = 'Alba Vilanova Cortezón'
                combined_ds.attrs['creator_email'] = 'alba.vilanova@bsc.es'
                combined_ds.attrs['application_area'] = 'Monitoring atmospheric composition'
                combined_ds.attrs['domain'] = 'Atmosphere'
                combined_ds.attrs['observed_layer'] = 'Land surface'

                # save data per year and month
                path = join(self.nonghost_root, f'actris/actris/{resolution}/{var}')
                if not os.path.isdir(path):
                    os.makedirs(path, exist_ok=True)
                saved_files = 0
                for year, ds_year in combined_ds.groupby('time.year'):
                    for month, ds_month in ds_year.groupby('time.month'):
                        filename = f"{path}/{var}_{year}{month:02d}.nc"
                        if filename in files_to_download:
                            combined_ds_yearmonth = combined_ds.sel(time=f"{year}-{month:02d}")

                            # add title to attrs
                            extra_info = ''
                            wavelength_var = is_wavelength_var(actris_parameter)
                            if wavelength_var and wavelength is not None:
                                extra_info = f' at {wavelength}nm'
                            combined_ds_yearmonth.attrs['title'] = f'Surface {parameters_dict[var]}{extra_info} in the ACTRIS network in {year}-{month:02d}.'

                            # order attrs
                            custom_order = ['title', 'institution', 'creator_name', 'creator_email',
                                            'source', 'application_area', 'domain', 'observed_layer',
                                            'data_license']
                            ordered_attrs = {key: combined_ds_yearmonth.attrs[key] 
                                            for key in custom_order 
                                            if key in combined_ds_yearmonth.attrs}
                            combined_ds_yearmonth.attrs = ordered_attrs

                            # remove stations if all variable data is nan
                            # previous_n_stations = len(combined_ds_yearmonth.station)
                            combined_ds_yearmonth = combined_ds_yearmonth.dropna(dim="station", subset=[var], how="all")
                            combined_ds_yearmonth = combined_ds_yearmonth.assign_coords(station=range(len(combined_ds_yearmonth.station)))
                            # current_n_stations = len(combined_ds_yearmonth.station)
                            # n_stations_diff = previous_n_stations - current_n_stations
                            # if n_stations_diff > 0:
                            #     self.logger.info(f'Data for {n_stations_diff} stations was removed because all data was NaN during {month}-{year}.')
                            
                            # remove file if it exists
                            if os.path.isfile(filename):
                                os.system("rm {}".format(filename))
                                
                            # save file
                            combined_ds_yearmonth.to_netcdf(filename)

                            # change permissions
                            os.system("chmod 777 {}".format(filename))
                            self.logger.info(f"Saved: {filename}")
                            saved_files += 1
                            
                self.logger.info(f'Total number of saved files: {saved_files}')

            else:
                self.logger.info('No files were found')

    def download_cams_experiment(self, experiment): 
        # print current_experiment
        self.logger.info('\n'+'-'*40)
        self.logger.info(f"\nDownloading {experiment} experiment data from the Atmosphere Data Store...")

        # get experiment id and the domain
        config_expid, domain, ensemble_options = experiment.split("-")
        
        # get the CAMS dataset and the experiment name (if there's one)
        exp_id = None
        if len(config_expid.split('_')) == 2:
            prefix = config_expid
        else:
            prefix, exp_id = config_expid.rsplit('_', 1)

        # check if the domain is the correct one for the dataset
        if domain not in cams_options[prefix]:
            possible_domains = "', '".join(cams_options[prefix])
            msg = (
            f"The current domain '{domain}' is not valid for the CAMS dataset."
            f"It must be '{possible_domains}'.")            
            show_message(self, msg)
            return

        # get prefix - domain dictionary
        cams_dict = cams_options[prefix][domain]

        # get CAMS dataset
        dataset = cams_dict["dataset"]

        # make sure the experiment is available in the dataset
        if 'experiments' in cams_dict and exp_id not in cams_dict["experiments"]:
            msg = f"Cannot find the {exp_id} experiment in the {dataset} dataset."    
            show_message(self, msg)
            return

        # download an experiment
        if 'experiments' not in cams_dict and exp_id is not None:
            msg = f"The {dataset} does not admit experiments, change the experiment in the configuration file to '{prefix}'."    
            show_message(self, msg)
            return

        # only ensemble options allmembers and 000 are valid
        if ensemble_options not in ['000', 'allmembers']:
            msg = (
            f"The current ensemble option '{ensemble_options}' is not valid for the CAMS '{dataset}' dataset."
            f"It must be '000' or 'allmembers'.")            
            show_message(self, msg)
            return
        
        # get minimum and maximum possible dates
        min_start_date, max_end_date = self.fetch_cams_dates(cams_dict['url'])

        # convert the selected dates to datetetime
        cams_start_date = datetime.strptime(self.start_date, "%Y%m%d")
        cams_end_date = datetime.strptime(self.end_date, "%Y%m%d")

        # if the minimum date is over the end date
        if min_start_date > cams_end_date or max_end_date < cams_start_date:
            msg = f"The selected dates are unavailable. Please choose dates between {min_start_date.strftime('%Y-%m-%d')} and {max_end_date.strftime('%Y-%m-%d')}."
            show_message(self, msg)
            return

        # check if the start date is within limits
        if min_start_date > cams_start_date:
            cams_start_date = min_start_date
            
        # check if the end date is within limits
        if max_end_date < cams_end_date:
            cams_end_date = max_end_date

        # get the ghost to cams vocabulary mapping
        parameters_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'ghost_cams_variables.yaml')))
        
        # create cdsapirc file in case it was not created
        cdsapirc_path = join(os.getenv("HOME"),'.cdsapirc')
        if not os.path.isfile(cdsapirc_path):
            
            # ask the user whether they want to crete the file in the home directory
            create_file = input(f"'.cdsapirc' file not found. Creating it at {cdsapirc_path}. Do you agree? ([y]/n)").lower()
            while create_file not in ['','y','n']:
                create_file = input(f"'.cdsapirc' file not found. Creating it at {cdsapirc_path}. Do you agree? ([y]/n)").lower()

            if create_file in ['', 'y']: 
                with open(cdsapirc_path, "w") as f:
                    f.write("url: https://ads.atmosphere.copernicus.eu/api\n")
                    f.write("key: 6101c82f-6d0b-4278-ac89-82743a81502c\n") # TODO get the user key
            else:
                self.logger.error("Error: Cannot proceed without '.cdsapirc'. CAMS experiment data download requires this file.")
                sys.exit(1)

        # iterate throught the resolutions
        for resolution in self.resolution:

            # check if the resolution is the correct one for the dataset
            if resolution != cams_dict["resolution"]:
                msg = (
                f"The current resolution '{resolution}' is not valid for the CAMS '{dataset}' dataset. "
                f"It must be '{cams_dict['resolution']}'.")            
                show_message(self, msg)
                continue
            
            # iterate through the species
            for species in self.species: 
                # check if species is in the ghost_cams_variables file
                if species not in parameters_dict:
                    msg = f"The species '{species}' is not available in CAMS."
                    show_message(self, msg)
                    continue
            
                # get the species in the cams vocabulary
                cams_species = parameters_dict[species]

                # check if the mapped species are available in the dataset
                if cams_species not in cams_dict['variable']:
                    msg = f"Mapped species '{cams_species}' for input species '{species}' is not available in the CAMS '{dataset}' dataset."          
                    show_message(self, msg)
                    continue

                # create client
                self.logger.info('')
                client = cdsapi.Client(retry_max=1, quiet=True)

                # get directory structure
                dir_tail = join(config_expid, domain, resolution, species)

                # get temporal and final dir
                temp_dir = join(self.exp_to_interp_root,'.temp', dir_tail)
                final_dir = join(self.exp_to_interp_root, dir_tail)

                # create temporal and final dirs to store the middle zip file with its directories
                os.makedirs(temp_dir, exist_ok=True)
                os.makedirs(final_dir, exist_ok=True)

                # iterate through the dates
                current_cams_date = cams_start_date
                while current_cams_date <= cams_end_date:

                    # add one day or one month depending if it is forecast or analysis
                    if cams_dict['forecast'] is True:
                        next_cams_date = current_cams_date.replace(day=1) + relativedelta(months=1)
                    else:
                        next_cams_date =  current_cams_date + timedelta(days=1) 

                    # get the date in cams str format
                    current_cams_date_str = current_cams_date.strftime('%Y-%m-%d')

                    # create the request
                    request = {
                    "variable": [cams_species],
                    "time": ["00:00"],
                    "data_format": "netcdf_zip"
                    }

                    # add leadtime hour to the request if the dataset has it
                    if 'leadtime_hour' in cams_dict:
                        request["leadtime_hour"] = list(range(0,cams_dict['leadtime_hour']+1))

                    # add type to the request if the dataset has it
                    if 'type' in cams_dict:
                        request["type"] = cams_dict['type']

                    # if it's forecast one file per day, analysis one file per month
                    if cams_dict['forecast'] is True:  # monthly
                        request["date"] = [f"{current_cams_date_str}/{(next_cams_date - timedelta(days=1)).strftime('%Y-%m-%d')}"]
                    else: # daily
                        request["date"] = [f"{current_cams_date_str}/{current_cams_date_str}"]

                    # add the experiment if models are available in the dataset
                    if 'experiments' in cams_dict:
                        request["model"] = [exp_id]

                    # if species is multi level, get the 0
                    if cams_species in cams_variables_level['multi']:
                        request["level"] = ["0"]

                    # get file name and final path
                    file_name = f"{species}-000_{current_cams_date.strftime('%Y%m%d')}.nc"
                    final_path = join(final_dir, file_name)

                    # get temporal path
                    temp_path = join(temp_dir, 'zip_file')

                    # print the request
                    self.logger.info(f"Dataset -> {cams_dict['dataset']}")
                    self.logger.info('Request -> {')
                    for k,v in request.items():
                        self.logger.info(f"{k} : {v}")
                    self.logger.info('}\n')

                    # make the request
                    try:
                        self.logger.info(f"Downloading {final_path}") # TODO change message
                        client.retrieve(dataset, request, target=temp_path)
                    except requests.exceptions.HTTPError as err:
                        # bad request
                        if err.response.status_code == 400: 
                            self.logger.info("\nBad request (400): The server could not understand the request.")
                            self.logger.info(f"Details: {err}")
                        # connection error
                        elif err.response.status_code == 500: 
                            self.logger.info("\nServer error (500): The server encountered an error while processing the request.")
                            self.logger.info(f"Details: {err}")
                            self.logger.info("Please try again later.")
                            return
                        else:
                            self.logger.info(f"\nUnexpected error ({err.response.status_code}):")
                            self.logger.info(f"Details: {err}")
                        # do the next download
                        continue

                    # extract file 
                    with zipfile.ZipFile(temp_path, 'r') as zip_ref:
                        zip_file_name = zip_ref.namelist()[0]
                        zip_ref.extractall(temp_dir)

                    # format the cams files and move them to the corresponding folder
                    self.logger.info(f"Formatting {final_path}\n") 
                    self.format_cams(join(temp_dir,zip_file_name), final_path, cams_species, species)

                    # add one day to the date
                    current_cams_date = next_cams_date    

                # remove the temp directory tail
                shutil.rmtree(join(self.exp_to_interp_root,'.temp'))

    def format_cams(self, input_filepath, output_filepath, cams_species, species):  
        # open original netcdf file      
        og_grp = Dataset(input_filepath, 'r', format="NETCDF4") 

        # extract date 
        date_str = og_grp['time'].long_name.split()[-1]

        # create new netcdf file
        root_grp = Dataset(output_filepath, 'w', format="NETCDF4") 
        root_grp.set_auto_mask(True)	
        
        # give permission to others 

        # copy global attributes
        root_grp.setncatts({k: og_grp.getncattr(k) for k in og_grp.ncattrs()})

        # copy dimensions 
        for name, dim in og_grp.dimensions.items():
            # skip level
            if name == 'level': 
                continue
            root_grp.createDimension(name, len(dim))

        # copy variables
        for name, og_var in og_grp.variables.items():
            # skip level
            if name == "level":
                continue

            # change species name
            if name not in ["longitude","latitude","time"]:
                name = species

            # get dimensions without level
            dims = tuple(dim for dim in og_var.dimensions if dim != "level")
            
            # create the variable
            var = root_grp.createVariable(name, og_var.datatype, dims)

            # add atribures
            if name == "time":
                # add specific atributes to time
                var.units = f'hours since {date_str[:4]}-{date_str[4:6]}-{date_str[-2:]} 00:00:00'                
                var.calendar = 'standard'
            else:
                # copy atributes from the original file
                var.setncatts({k: og_var.getncattr(k) for k in og_var.ncattrs() if k != '_FillValue'})

            # get the data from the original file
            data = og_var[:]

            # remove level dimension from species
            if name == species:
                data = np.squeeze(data, axis=1)
            
            # add the data to the variable
            var[:] = data

        # add grid_mapping
        root_grp[species].setncattr('grid_mapping', 'crs')
        
        # add coordinates
        root_grp[species].setncattr('coordinates', 'latitude longitude')

        # add crs
        crs_var = root_grp.createVariable('crs', 'u1')  
        crs_var.setncatts({
            'grid_mapping_name': 'latitude_longitude',
            'semi_major_axis': 6371000.0,
            'inverse_flattening': 0.0
        })
        
        # close the original and new netcdf files
        root_grp.close()
        og_grp.close()           

    def check_time(self, size, file_size):
        if (time.time() - self.ncfile_dl_start_time) > self.timeoutLimit:
            error = 'Download timeout, try later.'
            self.logger.error(error)
            sys.exit(1)
            
    def sighandler(self, *unused):
        self.logger.info('\nKeyboard Interrupt. Stopping execution.')
        
        # close connection, if it exists
        if hasattr(self, 'ssh'): 
            if self.ssh is not None:
                self.logger.info("\nClosing ssh connection...")
                self.ssh.close()
                if hasattr(self, 'sftp'):
                    self.sftp.close()

        # delete the las downloaded nc file to avoid corrupted files
        if hasattr(self, 'latest_nc_file_path'):
            self.logger.info(f"\nDeleting last file to avoid corruption: {self.latest_nc_file_path}...")
            if os.path.isfile(self.latest_nc_file_path):
                os.remove(self.latest_nc_file_path)

        # remove the output files from dtrsync in case it was a download from storage5
        if self.machine == "storage5":
            time.sleep(3)
                
            for file in os.listdir(PROVIDENTIA_ROOT):
                if file.startswith("dtrsync_"):
                    os.remove(join(PROVIDENTIA_ROOT,file)) 

        # delete temp dir if necessary
        temp_dir = join(self.ghost_root,'.temp')
        if os.path.exists(temp_dir):
            self.logger.info(f"\nDeleting {temp_dir}")
            shutil.rmtree(temp_dir)
        
        self.logger.info("\nExiting...")
        sys.exit()


def main(**kwargs):
    """ Main function when running download function. """

    download = Download(**kwargs)
    download.run()