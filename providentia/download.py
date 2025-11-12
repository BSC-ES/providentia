import sys
import os
import shutil

from base64 import decodebytes
import copy
from dotenv import dotenv_values
from getpass import getpass
import paramiko 
from remotezip import RemoteZip
import requests
from remotezip import RemoteZip
import signal
import subprocess
import tarfile
import time
from tqdm import tqdm
import xarray as xr
import yaml 

from .actris import Actris
from .cams import (Cams, cams_options)
from providentia.auxiliar import CURRENT_PATH, join
from .configuration import ProvConfiguration, load_conf
from .read_aux import check_for_ghost
from .warnings_prv import show_message
from .zenodo import Zenodo

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)
REMOTE_MACHINE = "storage5"

# load the defined experiments paths, agrupations yaml and mapping species
data_paths = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/data_paths.yaml')))
interp_experiments = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'interp_experiments.yaml')))
mapping_species =  yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'mapping_species.yaml')))

class Download(object):
    def __init__(self, **kwargs):

        # get providentia start time
        self.prov_start_time = time.time()

        # get ssh user and password 
        env = dotenv_values(join(PROVIDENTIA_ROOT, ".env"))

        # get ssh user and password 
        self.prv_user = env.get("PRV_USER")
        self.prv_password = env.get("PRV_PWD")

        # get origin update (ACTRIS)
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
                        msg = 'Error: The section specified in the command line ({0}) does not exist.'.format(kwargs['section'])
                        msg += '\nTip: For subsections, add the name of the parent section followed by an interpunct (·) '
                        msg += 'before the subsection name (e.g. SECTIONA·Spain). Available: {0}'.format(self.all_sections)
                        self.logger.error(msg)
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

            # initialise remote hostname
            self.remote_hostname = "transfer1.bsc.es"

            # initialise ssh 
            self.ssh = None

            # initialise boolean thath indicates whether remote machine changed 
            self.switched_remote = False

    def run(self, **kwargs):
        for section_ind, section in enumerate(self.sections):
            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # update self with section variables
            for k, val in self.section_opts.items():
                if k not in kwargs:
                    setattr(self, k, self.provconf.parse_parameter(k, val))
            
            # now all variables have been parsed, check validity of those, throwing errors where necessary
            self.provconf.check_validity()

            # from here generate control if user stopped execution
            signal.signal(signal.SIGINT, self.sighandler)
            
            # only the local download iterates through the networks
            if self.machine in "local":

                # networks
                if self.network:

                    if self.network == ["*"]:
                        # get user input to know which kind of network wants
                        download_source = input("\nDo you want to download all the GHOST networks? (Otherwise all the non-GHOST networks will be downloaded) ([y]/n): ")
                        while download_source.lower() not in ['','y','n']:
                            download_source = input("\nDo you want to download all the GHOST networks? (Otherwise all the non-GHOST networks will be downloaded) ([y]/n): ")
                        self.reading_ghost = download_source.lower() in ['','y']

                    # if there are GHOST networks, ask the user whether they want to download it from zenodo or HPC machines
                    if self.reading_ghost:
                        # ask whether the user wants to download from the zenodo or bsc machine
                        self.bsc_download = input("\nDo you want to download from the BSC remote machine? (Otherwise, GHOST data will be retrieved from Zenodo) ([y]/n): ")
                        while self.bsc_download.lower() not in ['','y','n']:
                            self.bsc_download = input("\nDo you want to download from the BSC remote machine? (Otherwise, GHOST data will be retrieved from Zenodo) ([y]/n): ")

                        # initialise the Zenodo object if user chose a Zenodo download
                        if self.bsc_download.lower() not in ['', 'y']:
                            self.zenodo = Zenodo(self)        

                    # get all networks if wildcard is passed
                    if self.network == ["*"]:
                        self.get_all_networks()

                    # combine all networks and species combinations to download (for network and filter species)
                    combined_networks = [(network, None) for network in self.network] + \
                                        [(network_specie.split('|')[0], network_specie.split('|')[1]) for network_specie in self.filter_species]
                    
                    # save main species
                    main_species = copy.deepcopy(self.species)

                    # get the download function 
                    download_fun = (
                    self.download_ghost_network_sftp if self.reading_ghost and self.bsc_download.lower() in ['', 'y']
                    else self.zenodo.download_ghost_network_zenodo if self.reading_ghost
                    else self.download_nonghost_network)
                    
                    # download network observations with species and filter_species
                    for network, filter_species in combined_networks:
                        # change species when turn of filter species
                        if filter_species is not None:
                            self.species = [filter_species]

                        # ACTRIS
                        if network == 'actris/actris':   
                            for resolution in self.resolution:
                                actris = Actris(self, resolution)
                                actris.download_actris_network()
                        # GHOST and non-GHOST 
                        else:
                            # download GHOST network
                            initial_check_nc_files = download_fun(network, initial_check=True)
                            files_to_download = self.select_files_to_download(initial_check_nc_files)
                            if not initial_check_nc_files or files_to_download:
                                download_fun(network, initial_check=False, files_to_download=files_to_download)

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
                            self.cams = Cams(self)
                            self.cams.download_cams_experiment(experiment)
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

            # reset domain and ensemble for new section
            self.domain = []
            self.ensemble = []

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
            remind_txt = input("\nRemember credentials ([y]/n)? ")
            while remind_txt.lower() not in ['y','n']:
                remind_txt = input("\nRemember credentials ([y]/n)? ")
            
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
                if not isinstance(self.dl_overwrite, bool):
                    # ask if user wants to overwrite
                    overwrite_choice = None
                    while overwrite_choice not in ['y','n','']:
                        overwrite_choice = input("\nThere are some files that were already downloaded in a previous download, do you want to overwrite them ([y]/n)? ").lower() 
                    self.dl_overwrite = overwrite_choice != 'n'

                # if user wants to overwrite then add the files downloaded before the execution as if they were never downloaded
                if self.dl_overwrite:
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
            msg = f"The {network} network could not be found on {join(PROVIDENTIA_ROOT,'settings','available_inputs.yaml')} nonghost_available_networks list."
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
            msg = f"The resolution/s {', '.join(available_resolutions)} could not be found on {join(PROVIDENTIA_ROOT,'settings','available_inputs.yaml')} nonghost_available_resolutions list."
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
            msg = f"There is no data available in {REMOTE_MACHINE} for {network} network for the current GHOST version ({self.ghost_version})."
            
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
                msg += f" Please check one of the available versions for the current species, resolution and dates: {', '.join(sorted(valid_available_ghost_versions))}"
            elif available_ghost_versions:
                msg += " There are no other versions available at the moment for the current species, resolution and dates with this configuration."
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
        
        # domain and ensemble are part of the experiment name, all united by dash (-)
        experiment_new = experiment

        # domain and ensemble are directories
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
                msg += f" Please check one of the available versions for the current species, resolution and dates: {', '.join(sorted(valid_available_ghost_versions))}"
            elif available_ghost_versions:
                msg += " There are no other versions available at the moment for the current species, resolution and dates with this configuration."
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
                
                # get network, species and resolution
                network = remote_dir.split('/')[-1]
                species = remote_dir.split('/')[-2]
                resolution = remote_dir.split('/')[-3]

                # get local directory 
                local_dir = join(self.exp_root,self.ghost_version,experiment_new,resolution,species,network)

                # print source and destination  
                if not initial_check:
                    self.logger.info(f"\n  - {local_dir}, source: {remote_dir} ({REMOTE_MACHINE})")
            
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
        exp_id, domain, ensemble = experiment.split("-")

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
                    if ensemble.isdigit() or ensemble in ['allmembers', 'av_an']:
                        res_spec = join(remote_dir,resolution,species)
                    # if it is an ensemble statistic
                    else:
                        res_spec = join(remote_dir,resolution,"ensemble-stats",species+"_"+ensemble)
  
                    self.sftp.stat(res_spec)
                    species_exists = True
                # if there are none, try with the mapped species
                except FileNotFoundError:
                    # change species name to the species to map
                    if speci_to_process in mapping_species:
                        for mapping_speci in mapping_species[speci_to_process]:
                            try:
                                # if it is an ensemble member
                                if ensemble.isdigit() or ensemble in ['allmembers', 'av_an']:
                                    res_spec = join(remote_dir,resolution, mapping_speci)
                                # if it is an ensemble statistic
                                else:
                                    res_spec = join(remote_dir,resolution, "ensemble-stats", species + "_" + ensemble)
  
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
                         
                # get nc files
                nc_files = self.sftp.listdir(remote_dir)

                if nc_files:
                    # if it is an ensemble member
                    if ensemble.isdigit() or ensemble in ['allmembers', 'av_an']:
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
                        if format == (0, 1):
                            # when there is no ensemble in the name only allmembers and 000 are valid
                            if ensemble == '000' or ensemble == 'allmembers':
                                nc_files = list(filter(lambda x:x.split("_")[0] == species, nc_files))

                        # example: od550du-000_2021020812.nc (1,1)
                        elif format == (1, 1):
                            # filter by ensemble in case that ensemble is not allmembers
                            if ensemble != 'allmembers':
                                nc_files = list(filter(lambda x:x.split("_")[0] == species+'-'+ensemble,nc_files))
                        
                        # example: od550du_2018011700_av_an.nc
                        elif format == (0, 3):
                            nc_files = list(filter(lambda x:x.split("_")[0] == species, nc_files))

                        # unknown format
                        else:
                            error = f"It is not possible to download this nc file type yet. Please, contact the developers. Files to download: {nc_files}"
                            self.logger.error(error)
                            sys.exit(1)
                    
                    # if it is an ensemble statistic
                    else:
                        # get the domain, resolution and species from the path
                        domain, resolution, _, species = remote_dir.split('/')[-4:]
                        species = species.split("_",1)[0]

                        # filter the nc files to only get the ones that have the correct species and stats
                        nc_files = list(filter(lambda x:x.split("_")[0] == species and "_".join(x[:-3].split("_")[2:]) == ensemble, nc_files))
                
                # add ensemble-stats directory if it is an ensemble member
                if ensemble.isdigit() or ensemble in ['allmembers', 'av_an']:
                    local_dir = join(self.exp_to_interp_root,exp_id,domain,resolution,species)
                else:
                    local_dir = join(self.exp_to_interp_root,exp_id,domain,resolution,"ensemble-stats",species+"_"+ensemble)

                # print source and destination
                if not initial_check:
                    self.logger.info(f"\n  - {local_dir}, source: {remote_dir} ({REMOTE_MACHINE})")
                    
                # if there is no options with the ensemble, tell the user
                if nc_files == []:
                    msg = f"There is no data available in {REMOTE_MACHINE} for the {exp_id} experiment with the {domain} domain with the {ensemble} ensemble."
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
        exp_id, domain, ensemble = experiment.split("-")

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
                if ensemble.isdigit() or ensemble in ['allmembers', 'av_an']:
                    res_spec = join(esarchive_dir,resolution,species)
                # if it is an ensemble statistic
                else:
                    res_spec = join(esarchive_dir,resolution,"ensemble-stats",species+"_"+ensemble)
                species_exists = os.path.exists(res_spec)
                # if there are none, try with the mapped species
                if species_exists is False:
                    # change species name to the species to map
                    if speci_to_process in mapping_species:
                        for species in mapping_species[speci_to_process]:
                            # if it is an ensemble member
                            if ensemble.isdigit() or ensemble in ['allmembers', 'av_an']:
                                res_spec = join(esarchive_dir,resolution,species)
                            # if it is an ensemble statistic
                            else:
                                res_spec = join(esarchive_dir,resolution,"ensemble-stats",species+"_"+ensemble)
  
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
                    if ensemble.isdigit() or ensemble in ['allmembers', 'av_an']:
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
                            # when there is no ensemble in the name only allmembers and 000 are valid
                            if ensemble == '000' or ensemble == 'allmembers':
                                nc_files = list(filter(lambda x:x.split("_")[0] == species, nc_files))
                        
                        # example: od550du-000_2021020812.nc (1,1)
                        elif format == (1,1):
                            # filter by ensemble in case that ensemble is not allmembers
                            if ensemble != 'allmembers':
                                nc_files = list(filter(lambda x:x.split("_")[0] == species+'-'+ensemble,nc_files))
                           
                        # example: od550du_2018011700_av_an.nc
                        elif format == (0, 3):
                            nc_files = list(filter(lambda x:x.split("_")[0] == species, nc_files))

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
                        nc_files = list(filter(lambda x:x.split("_")[0] == species and "_".join(x[:-3].split("_")[2:]) == ensemble, nc_files))
                        
                # if there is no options with the ensemble, tell the user
                if nc_files == []:
                    msg = f"There is no data available in esarchive for the {exp_id} experiment with the {domain} domain with the {ensemble} ensemble."
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
                    if ensemble.isdigit() or ensemble in ['allmembers', 'av_an']:
                        gpfs_dir = join(self.exp_to_interp_root,exp_id,domain,resolution,species)
                    else:
                        gpfs_dir = join(self.exp_to_interp_root,exp_id,domain,resolution,"ensemble-stats",species+"_"+ensemble)
                    
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

                            # get last downloaded file in case there was a keyboard interrupt
                            self.latest_nc_file_path = gpfs_path

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

    def get_all_networks(self): 
        if self.reading_ghost:
            if self.bsc_download.lower() in ['', 'y']:
                self.network = self.ghost_available_networks
            elif not hasattr(self,"zenodo_ghost_available_networks"): 
                if not self.zenodo.fetch_zenodo_networks():
                    return
                self.network = list(self.zenodo_ghost_available_networks.keys())
        else:
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
            # get all the domain and ensemble combinations 
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
            if os.path.isfile(self.latest_nc_file_path):
                self.logger.info(f"\nDeleting last file to avoid corruption: {self.latest_nc_file_path}...")
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
    download.run(**kwargs)