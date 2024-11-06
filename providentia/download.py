import sys
import os
import shutil

import requests
from io import BytesIO
import subprocess
import yaml
import json
from dotenv import dotenv_values
import paramiko 
from base64 import decodebytes
import signal
import copy
import time

# urlparse
from tqdm import tqdm
from remotezip import RemoteZip
import tarfile
from datetime import datetime, timedelta
from getpass import getpass

from .configuration import ProvConfiguration, load_conf
from .read_aux import check_for_ghost
from .warnings_prv import show_message

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)
REMOTE_MACHINE = "storage5"

# load the defined experiments paths, agrupations yaml and mapping species
data_paths = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/data_paths.yaml')))
interp_experiments = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings', 'interp_experiments.yaml')))
mapping_species =  yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings', 'internal', 'mapping_species.yaml')))

def check_time(size, file_size):
    if (time.time() - download.ncfile_dl_start_time) > download.timeoutLimit:
        error = 'Download timeout, try later.'
        sys.exit(error)
        
def sighandler(*unused):
    print('\nKeyboard Interrupt. Stopping execution.')
    
    # close connection, if it exists
    if download.ssh is not None:
        print("\nClosing ssh connection...")
        download.ssh.close()
        download.sftp.close()

    # delete the las downloaded nc file to avoid corrupted files
    if hasattr(download,'latest_nc_file_path'):
        print(f"\nDeleting last file to avoid corruption: {download.latest_nc_file_path}...")
        os.remove(download.latest_nc_file_path)

    # delete temp dir if necessary
    temp_dir = os.path.join(download.ghost_root,'.temp')
    if os.path.exists(temp_dir):
        print(f"\nDeleting {temp_dir}")
        shutil.rmtree(temp_dir)
    
    print("\nExiting...")
    sys.exit()

class ProvidentiaDownload(object):
    def __init__(self, **kwargs):
        # TODO move some of these variables to configuration.py
        # initialise zenodo url
        self.ghost_url = 'https://zenodo.org/records/10637450'

        # get providentia start time
        self.prov_start_time = time.time()

        # initialise remote hostname
        self.remote_hostname = "transfer1.bsc.es"

        # in case transfer broke
        # global REMOTE_MACHINE
        # self.remote_hostname, REMOTE_MACHINE = "glogin4.bsc.es", "mn5" 

        # get ssh user and password 
        env = dotenv_values(os.path.join(PROVIDENTIA_ROOT, ".env"))

        # get ssh user and password 
        self.prv_user = env.get("PRV_USER")
        self.prv_password = env.get("PRV_PWD")

        # get user preference over GHOST download and overwrite
        self.bsc_download_choice = env.get("BSC_DL_CHOICE")
        self.overwrite_choice = env.get("OVERWRITE")

        # set timeout limit
        self.timeoutLimit = 3 * 60

        # initialise default configuration variables
        # modified by commandline arguments, if given
        self.provconf = ProvConfiguration(self, **kwargs)

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
                # get the section in case it was passed on the command line                
                if 'section' in kwargs:
                    # config and section defined 
                    if kwargs['section'] in self.all_sections:
                        self.sections = [kwargs['section']]
                    else:
                        error = 'Error: The section specified in the command line ({0}) does not exist.'.format(kwargs['section'])
                        tip = 'Tip: For subsections, add the name of the parent section followed by an interpunct (·) '
                        tip += 'before the subsection name (e.g. SECTIONA·Spain).'
                        sys.exit(error + '\n' + tip)
                # if no section passed, then get all the parent sections
                else:
                    # if no parent section names are found throw an error
                    if len(self.parent_section_names) == 0:
                        error = "Error: No sections were found in the configuration file, make sure to name them using square brackets."
                        sys.exit(error)
                    self.sections = self.parent_section_names
            else:
                error = 'Error: The path to the configuration file specified in the command line does not exist.'
                sys.exit(error)
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            sys.exit(error)

        # create empty directories for all the paths if they don't exist
        for path in [self.nonghost_root,self.ghost_root,self.exp_root,self.exp_to_interp_root]:
            if not os.path.exists(path):
                try:
                    os.makedirs(path)
                except PermissionError as error:
                    os.system(f"sudo mkdir -p {path}")
                    os.system(f"sudo chmod o+w {path}")

        # initialise type of download
        if not self.bsc_download_choice:
            self.confirm_bsc_download()

        # initialise ssh 
        self.ssh = None

        # initialise boolean thath indicates whether remote machine changed 
        self.switched_remote = False

        # variable that saves whether some experiments/observations were downloaded before
        self.overwritten_files_flag = False

    def run(self):
        for section_ind, section in enumerate(self.sections):
            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # update self with section variables
            for k, val in self.section_opts.items():
                setattr(self, k, self.provconf.parse_parameter(k, val))

            # now all variables have been parsed, check validity of those, throwing errors where necessary
            self.provconf.check_validity(deactivate_warning=True)

            # if networks is none and is not the non interpolated mode, raise error
            if not self.network and self.interpolated is True:
                error = "Error: No networks were passed."
                sys.exit(error)
            
            # when one of those symbols is passed, get all networks
            if self.network == ["*"]:
                self.get_all_networks()

            # from here generate control if user stopped execution
            signal.signal(signal.SIGINT, sighandler)

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

                    # get the files to be downlaoded, check if they files were already downlaoded and download if not
                    # download from the remote machine
                    if self.bsc_download_choice == 'y':
                        # GHOST
                        if check_for_ghost(network):
                            initial_check_nc_files = self.download_ghost_network_sftp(network, initial_check=True)
                            files_to_download = self.select_files_to_download(initial_check_nc_files)
                            if not initial_check_nc_files or files_to_download:
                                self.download_ghost_network_sftp(network, initial_check=False, files_to_download=files_to_download)
                        # non-GHOST
                        else:
                            initial_check_nc_files = self.download_nonghost_network(network, initial_check=True)
                            files_to_download = self.select_files_to_download(initial_check_nc_files)
                            if not initial_check_nc_files or files_to_download:
                                self.download_nonghost_network(network, initial_check=False, files_to_download=files_to_download)

                    # download from the zenodo webpage
                    elif self.bsc_download_choice == 'n':
                        # GHOST
                        if check_for_ghost(network):
                            initial_check_nc_files = self.download_ghost_network_zenodo(network, initial_check=True)
                            files_to_download = self.select_files_to_download(initial_check_nc_files)
                            if not initial_check_nc_files or files_to_download:
                                self.download_ghost_network_zenodo(network, initial_check=False, files_to_download=files_to_download)
                        # non-GHOST
                        else:
                            error = f"Error: It is not possible to download files from the non-GHOST network {network} from the zenodo webpage."
                            sys.exit(error)
                    
                    # download option invalid
                    else:
                        error = "Error: Download option not valid, check your .env file."
                        sys.exit(error)

                # get orignal species back
                self.species = main_species

            # when one of those symbols is passed, get all experiments
            if self.experiments == {'*': '*'}:
                self.get_all_experiments()

            # experiment
            # download from the remote machine
            if self.experiments:
                if self.bsc_download_choice == 'y':
                    # get function to download experiment depending on the configuration file field
                    self.download_experiment_fun = self.download_experiment if self.interpolated is True else self.download_non_interpolated_experiment
                    # iterate the experiments download
                    for experiment in self.experiments.keys():
                        initial_check_nc_files = self.download_experiment_fun(experiment, initial_check=True)
                        files_to_download = self.select_files_to_download(initial_check_nc_files)
                        if not initial_check_nc_files or files_to_download:
                            self.download_experiment_fun(experiment, initial_check=False, files_to_download=files_to_download)
                # download from the zenodo webpage
                else:
                    error = f"Error: It is not possible to download experiments from the zenodo webpage."
                    sys.exit(error)

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
            print("\nSome experiments/observations were found but were not downloaded because the OVERWRITE option is set to 'n'.")

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
                sys.exit(error)

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
                    error = f"Authentication failed. Please, check if PRV_USER on {os.path.join(PROVIDENTIA_ROOT, '.env')} aligns with your BSC {REMOTE_MACHINE} ssh user."
                    error += "\nIf it does not, change the user to the correct one. If it does, delete the whole PRV_USER row and execute again."
                    sys.exit(error)
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
                error += f" Please, check your credentials on {os.path.join(PROVIDENTIA_ROOT, '.env')}"
            sys.exit(error)

        # if pwd or user changed, ask if user wants to remember credentials
        if (prv_user is not None) or (prv_password is not None):
            # ask user if they want their credentials saved
            remind_txt = input("\nRemember credentials (y/n)? ")
            while remind_txt.lower() not in ['y','n']:
                remind_txt = input("\nRemember credentials (y/n)? ")
            
            # create .env with the input user and/or password
            if remind_txt.lower() == 'y':
                with open(os.path.join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                    if prv_user is not None:
                        f.write(f"PRV_USER={self.prv_user}\n")
                    if prv_password is not None:
                        f.write(f"PRV_PWD={self.prv_password}\n")

                print(f"\nRemote machine credentials saved on {os.path.join(PROVIDENTIA_ROOT, '.env')}")

    def confirm_bsc_download(self):
        # get user choice regarding bsc downloads
        self.bsc_download_choice = None
        while self.bsc_download_choice not in ['y','n']:
            self.bsc_download_choice = input("\nDo you want to download from BSC remote machine (y/n)? ").lower()

        remind_txt = None
        while remind_txt not in ['y','n']:
            remind_txt = input("\nDo you want to remember your decision (y/n)? ").lower() 
        
        if remind_txt == 'y':
            with open(os.path.join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                f.write(f"BSC_DL_CHOICE={self.bsc_download_choice}\n")
            
            print(f"\nBSC download choice saved on {os.path.join(PROVIDENTIA_ROOT, '.env')}")
                    
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
                        with open(os.path.join(PROVIDENTIA_ROOT, ".env"),"a") as f:
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
            print('\n'+'-'*40)
            print(f"\nDownloading non-GHOST {network} network data from {REMOTE_MACHINE}...")

        # if not valid network, check if user put the network on init_prov 
        # TODO Move to configuration.py
        if network not in self.nonghost_available_networks:
            msg = f"The {network} network could not be found on {os.path.join(PROVIDENTIA_ROOT,'settings','init_prov.yaml')} nonghost_available_networks list."
            msg += "\nPlease, add the network to the list and execute again."
            show_message(self, msg, deactivate=initial_check)
            return
        
        # check if nonghost network exists in directory
        # TODO: Change this to somewhere in configuration, the one up too
        try:
            self.sftp.stat(os.path.join(self.nonghost_remote_obs_path,network))
        except FileNotFoundError:
            msg = f"There is no data available in {REMOTE_MACHINE} for {network} network."
            show_message(self, msg, deactivate=initial_check)
            return

        # check if all resolutions are in init_prov, if not warning and delete the not correct ones
        # TODO move to configuration.py
        not_available_resolutions = set(self.resolution) - set(self.nonghost_available_resolutions)
        if not_available_resolutions:
            available_resolutions = set(self.resolution) - not_available_resolutions
            msg = f"The resolution/s {', '.join(not_available_resolutions)} could not be found on {os.path.join(PROVIDENTIA_ROOT,'settings','init_prov.yaml')} nonghost_available_resolutions list."
            msg += "\nPlease, add the necessary resolutions to the list and execute again."
            show_message(self, msg, deactivate=initial_check)
            return

        # get resolution and species combinations
        res_spec_dir = []

        sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(os.path.join(self.nonghost_remote_obs_path,network))).intersection(self.nonghost_available_resolutions)
        for resolution in sftp_resolutions:
            try:
                sftp_species = self.species if self.species else set(self.sftp.listdir(os.path.join(self.nonghost_remote_obs_path,network,resolution))).intersection(self.available_species)
            except FileNotFoundError:
                msg = f"There is no data available in {REMOTE_MACHINE} for {network} network at {resolution} resolution"
                show_message(self, msg, deactivate=initial_check)
                continue
            for species in sftp_species: 
                res_spec_dir.append(os.path.join(self.nonghost_remote_obs_path,network,resolution,species))
        
        if res_spec_dir:
            
            # initialise list with all the nc files to be downloaded
            initial_check_nc_files = []

            if not initial_check:
                print(f"\n{network} observations to download ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range within the specie and resolution combination
            for remote_dir in res_spec_dir:

                local_dir = os.path.join(self.nonghost_root,remote_dir.split('/',6)[-1])
                species = remote_dir.split('/')[-1]
                resolution = remote_dir.split('/')[-2]

                #  print the species, resolution and network combinations that are going to be downloaded
                if not initial_check:
                    print(f"\n  - {local_dir}")

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
                        valid_nc_files = list(filter(lambda x:os.path.join(local_dir,x) in files_to_download, valid_nc_files))
                        if not valid_nc_files:
                            msg = "Files were already downloaded."
                            show_message(self, msg, deactivate=initial_check)     
                            continue         
                        valid_nc_files_iter = tqdm(valid_nc_files,bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                    else:
                        # do not print the bar if it is the initial check
                        valid_nc_files_iter = valid_nc_files

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files_iter:
                        local_path = os.path.join(local_dir,nc_file)
                        if initial_check:
                            initial_check_nc_files.append(local_path)
                        else:
                            self.latest_nc_file_path = local_path
                            remote_path = os.path.join(remote_dir,nc_file)
                            self.ncfile_dl_start_time = time.time()
                            self.sftp.get(remote_path,local_path, callback=check_time)

            return initial_check_nc_files
        
    def download_ghost_network_sftp(self, network, initial_check, files_to_download=None):
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect() 

        if not initial_check:
            # print current_network
            print('\n'+'-'*40)
            print(f"\nDownloading GHOST {network} network data from {REMOTE_MACHINE}...")

        # if not valid network, next
        if network not in self.sftp.listdir(self.ghost_remote_obs_path):
            msg = f"There is no data available in {REMOTE_MACHINE} for {network} network."
            show_message(self, msg, deactivate=initial_check)
            return 
        
        # if not valid combination of GHOST version and network, next 
        elif self.ghost_version not in self.sftp.listdir(os.path.join(self.ghost_remote_obs_path,network)):
            msg = f"There is no data available in {REMOTE_MACHINE} for {network} network for the current ghost version ({self.ghost_version})."
            
            available_ghost_versions = set(self.sftp.listdir(os.path.join(self.ghost_remote_obs_path,network))).intersection(self.possible_ghost_versions)

            # list that saves the GHOST versions with valid nc files
            valid_available_ghost_versions = []
            
            # check for combinations of species, resolutions, network, and day in the available versions
            if available_ghost_versions:
                # iterate the different GHOST versions
                for possible_ghost_version in available_ghost_versions:
                    remote_dir_ghost_version = os.path.join(self.ghost_remote_obs_path, network, possible_ghost_version)
                    
                    # iterate the different resolutions
                    sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir_ghost_version)).intersection(self.ghost_available_resolutions)
                    for resolution in sftp_resolutions:
                        try:
                            species_list = self.species if self.species else self.sftp.listdir(os.path.join(remote_dir_ghost_version, resolution))
                        except FileNotFoundError:
                            continue

                        # iterate the different species
                        for species in species_list:
                            # look for valid nc files in the date range
                            try:
                                nc_files = self.sftp.listdir(os.path.join(remote_dir_ghost_version, resolution, species))
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

        remote_dir = os.path.join(self.ghost_remote_obs_path,network,self.ghost_version)

        sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir)).intersection(self.ghost_available_resolutions)
        for resolution in sftp_resolutions:
            try:
                sftp_species = self.species if self.species else set(self.sftp.listdir(os.path.join(remote_dir,resolution))).intersection(self.available_species)
            except FileNotFoundError:
                msg = f"There is no data available in {REMOTE_MACHINE} for {network} network at {resolution} resolution"
                show_message(self, msg, deactivate=initial_check)
                continue
            for species in sftp_species: 
                res_spec_dir.append(os.path.join(remote_dir,resolution,species))
        
        # print the species, resolution and network combinations that are going to be downloaded
        if res_spec_dir:            
            
            # initialise list with all the nc files to be downloaded
            initial_check_nc_files = []
            
            if not initial_check:
                print(f"\n{network} observations to download ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range within the specie and resolution combination
            for remote_dir in res_spec_dir:

                local_dir = os.path.join(self.ghost_root,remote_dir.split('/',7)[-1])
                species = remote_dir.split('/')[-1]
                resolution = remote_dir.split('/')[-2]

                #  print the species, resolution and network combinations that are going to be downloaded
                if not initial_check:
                    print(f"\n  - {local_dir}")

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
                        valid_nc_files = list(filter(lambda x:os.path.join(local_dir,x) in files_to_download, valid_nc_files))
                        if not valid_nc_files:
                            msg = "Files were already downloaded."
                            show_message(self, msg, deactivate=initial_check)     
                            continue         
                        valid_nc_files_iter = tqdm(valid_nc_files,bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                    else:
                        # do not print the bar if it is the initial check
                        valid_nc_files_iter = valid_nc_files

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files_iter:
                        local_path = os.path.join(local_dir,nc_file)
                        if initial_check:
                            initial_check_nc_files.append(local_path)
                        else:
                            self.latest_nc_file_path = local_path
                            remote_path = os.path.join(remote_dir,nc_file)
                            self.ncfile_dl_start_time = time.time()
                            self.sftp.get(remote_path,local_path, callback=check_time)
                       
            return initial_check_nc_files

    def download_ghost_network_zenodo(self, network, initial_check, files_to_download=None):
        if not initial_check:
            # print current_network
            print('\n'+'-'*40)
            print(f"\nDownloading GHOST {network} network data from Zenodo...")

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
                    print(f"\n{network} observations to download:")
                
                # initialise only possible GHOST version
                # TODO in the future there will be various versions in zenodo
                last_ghost_version = '1.5' 
                
                for remote_dir_tail in res_spec_dir_tail:
                    resolution, specie = remote_dir_tail.split("/")[1:]
                    specie = specie[:-7]
                    local_dir = os.path.join(self.ghost_root,network,last_ghost_version,resolution,specie)

                    #  print the species, resolution and network combinations that are going to be downloaded
                    if not initial_check:
                        print(f"\n  - {local_dir}")

                    # create temporal dir to store the middle tar file with its directories
                    temp_dir = os.path.join(self.ghost_root,'.temp')
                    if not os.path.exists(temp_dir):
                        os.mkdir(temp_dir)

                    zip.extract(remote_dir_tail,temp_dir)
                    
                    # get path and the name of the directory of the tar file
                    tar_path = os.path.join(self.ghost_root, network, str(last_ghost_version), *remote_dir_tail.split("/")[1:])
                    temp_path = os.path.join(temp_dir,remote_dir_tail)

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
                                valid_nc_file_names = list(filter(lambda x:os.path.join(local_dir,x.split('/')[-1]) in files_to_download, valid_nc_file_names))
                                if not valid_nc_file_names:
                                    msg = "Files were already downloaded."
                                    show_message(self, msg, deactivate=initial_check)     
                                    continue  
                                
                            valid_nc_files = list(filter(lambda x:x.name in valid_nc_file_names, tar_file.getmembers()))
                          
                            # sort nc_files
                            valid_nc_files.sort(key = lambda x:x.name)   

                            if not initial_check:    
                                valid_nc_files_iter = tqdm(valid_nc_files,bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                            else:
                                # do not print the bar if it is the initial check
                                valid_nc_files_iter = valid_nc_files

                            for nc_file in valid_nc_files_iter:
                                local_path = os.path.join(tar_dir,nc_file.name)
                                if initial_check:
                                    initial_check_nc_files.append(local_path)
                                else:
                                    self.latest_nc_file_path = local_path
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
            print('\n'+'-'*40)
            print(f"\nDownloading {experiment} experiment data from {REMOTE_MACHINE}...")
            
        # get resolution and species combinations
        res_spec_dir = []
        
        # domain and ensemble option are part of the experiment name, all united by dash (-)
        experiment_new = experiment

        # domain and ensemble option are directories
        experiment_old = experiment.replace("-","/")
        
        # get remote directory format depending on the GHOST version
        experiment = experiment_old if self.ghost_version in ["1.2", "1.3", "1.3.1"] else experiment_new

        # get remote directory
        remote_dir = os.path.join(self.exp_remote_path,self.ghost_version,experiment)

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
                    remote_dir_ghost_version = os.path.join(self.exp_remote_path, possible_ghost_version, experiment_old if possible_ghost_version in ["1.2", "1.3", "1.3.1"] else experiment_new)
                    
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
                    remote_dir_ghost_version = os.path.join(self.exp_remote_path, possible_ghost_version, experiment_old if possible_ghost_version in ["1.2", "1.3", "1.3.1"] else experiment_new)
                    
                    # iterate the different resolutions
                    sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir)).intersection(self.nonghost_available_resolutions)
                    for resolution in sftp_resolutions:                        
                        try:
                            species_list = self.species if self.species else self.sftp.listdir(os.path.join(remote_dir_ghost_version, resolution))
                        except FileNotFoundError:
                            continue
                        
                        # iterate the different species
                        for species in species_list:
                            try:
                                network_list = self.network if self.network else self.sftp.listdir(os.path.join(remote_dir_ghost_version, resolution, species))
                            except FileNotFoundError:
                                continue
                            
                            # iterate the different networks
                            for network in network_list:
                                # look for valid nc files in the date range
                                try:
                                    nc_files = self.sftp.listdir(os.path.join(remote_dir_ghost_version, resolution, species, network))
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
                sftp_species = self.species if self.species else set(self.sftp.listdir(os.path.join(remote_dir,resolution))).intersection(self.available_species)
            except FileNotFoundError:
                msg = f"There is no data available in {REMOTE_MACHINE} for {experiment_new} experiment at {resolution} resolution"
                show_message(self, msg, deactivate=initial_check)
                continue
            for species in sftp_species: 
                try:
                    sftp_network = self.network if self.network else self.sftp.listdir(os.path.join(remote_dir,resolution,species))
                except FileNotFoundError:
                    msg = f"There is no data available in {REMOTE_MACHINE} for {experiment_new} experiment for {species} species at {resolution} resolution"
                    show_message(self, msg, deactivate=initial_check)
                    continue
                for network in sftp_network:
                    # if network is nonghost, change the slashes to dashes
                    if not check_for_ghost(network):
                        network = network.replace("/", "-")
                    res_spec_dir.append(os.path.join(remote_dir,resolution,species,network))
        
        # print the species, resolution and experiment combinations that are going to be downloaded
        if res_spec_dir:

            # initialise list with all the nc files to be downloaded
            initial_check_nc_files = []

            if not initial_check:
                print(f"\n{experiment_new} experiment data to download ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range
            for remote_dir in res_spec_dir:
                if not initial_check:
                    local_path = remote_dir.split('/',7)[-1]
                    if self.ghost_version in ["1.2", "1.3", "1.3.1"]:
                        print(f"\n  - {os.path.join(self.exp_root,self.ghost_version,'-'.join(local_path.split('/')[1:4]),*local_path.split('/')[4:])}")
                    else:
                        print(f"\n  - {os.path.join(self.exp_root,local_path)}")
            
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
                    local_dir = os.path.join(self.exp_root,self.ghost_version,experiment_new,resolution,species,network)
                    
                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # sort nc_files
                    valid_nc_files.sort() 

                    if not initial_check:
                        # get the ones that are not already downloaded
                        valid_nc_files = list(filter(lambda x:os.path.join(local_dir,x) in files_to_download, valid_nc_files))
                        if not valid_nc_files:
                            msg = "Files were already downloaded."
                            show_message(self, msg, deactivate=initial_check)     
                            continue         
                        valid_nc_files_iter = tqdm(valid_nc_files, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                    else:
                        # do not print the bar if it is the initial check
                        valid_nc_files_iter = valid_nc_files

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files_iter:
                        local_path = os.path.join(local_dir,nc_file)
                        if initial_check:
                            initial_check_nc_files.append(local_path)
                        else:
                            self.latest_nc_file_path = local_path
                            remote_path = os.path.join(remote_dir,nc_file)
                            self.ncfile_dl_start_time = time.time()
                            self.sftp.get(remote_path,local_path, callback=check_time) 
            
            return initial_check_nc_files

    def download_non_interpolated_experiment(self, experiment, initial_check, files_to_download=None):
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect()  
        
        if not initial_check:
            # print current experiment
            print('\n'+'-'*40)
            print(f"\nDownloading {experiment} non-interpolated experiment data from {REMOTE_MACHINE}...")
            
        # get resolution and species combinations
        res_spec_dir = []

        # get experiment id and the domain
        exp_id, domain, ensemble_options = experiment.split("-")

        # get stat if it is an ensemble statistic
        if ensemble_options.startswith("stat_"): 
            stat = ensemble_options.split("_",1)[-1]

        # get experiment type
        for experiment_type, experiment_dict in interp_experiments.items():
            if exp_id in experiment_dict["experiments"]:
                break
        
        # get experiment specific directories list
        exp_dir_list = experiment_dict["paths"]

        # take all functional directories
        exp_dir_functional_list = []    
        for exp_dir in exp_dir_list:
            # esarchive in transfer5 is located inside gpfs
            if "/esarchive/" == exp_dir[:11]:
                exp_dir = os.path.join("/gpfs/archive/bsc32/",exp_dir[1:])
            # check if directory exists in the remote machine
            try:
                self.sftp.stat(exp_dir)
                exp_dir_functional_list.append(exp_dir)      
                break
            except FileNotFoundError:
                pass

        # if none of the paths are in this current machine, break
        if not exp_dir_functional_list:
            msg = f"None of the paths specified in {os.path.join('settings', 'interp_experiments.yaml')} are available on the remote machine ({REMOTE_MACHINE})."
            show_message(self, msg, deactivate=initial_check)
            return
        
        # take first functional directory  
        remote_dir = None
        for exp_dir in exp_dir_functional_list:
            temp_remote_dir = os.path.join(exp_dir,exp_id,domain)
            # check if remote experiment and domain directories exist in the remote machine
            try:
                self.sftp.stat(temp_remote_dir)
                remote_dir = temp_remote_dir
                break
            except FileNotFoundError:
                pass
        
        # if the experiment-domain combination is possible, break
        if remote_dir is None:
            msg = f"There is no data available for the {exp_id} experiment with the {domain} domain in none of the paths specified in {os.path.join('settings', 'interp_experiments.yaml')} in the remote machine ({REMOTE_MACHINE})."
            show_message(self, msg, deactivate=initial_check)
            return

        # get all the resolutions available in the remote directory
        sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir)).intersection(self.nonghost_available_resolutions)

        # iterate through the resolutions
        for resolution in sftp_resolutions:
            try:
                # get available species ("normal" and mapped)
                available_species = self.available_species+[spec[0] for spec in mapping_species.values()]
                sftp_species = self.species if self.species else set(self.sftp.listdir(os.path.join(remote_dir,resolution))).intersection(available_species)
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
                    if not ensemble_options.startswith("stat_"):
                        res_spec = os.path.join(remote_dir,resolution,species)
                    # if it is an ensemble statistic
                    else:
                        res_spec = os.path.join(remote_dir,resolution,"ensemble-stats",species+"_"+stat)
  
                    self.sftp.stat(res_spec)
                    species_exists = True
                # if there are none, try with the mapped species
                except FileNotFoundError:
                    # change species name to the species to map
                    if speci_to_process in mapping_species:
                        for species in mapping_species[speci_to_process]:
                            try:
                                # if it is an ensemble member
                                if not ensemble_options.startswith("stat_"):
                                    res_spec = os.path.join(remote_dir,resolution,species)
                                # if it is an ensemble statistic
                                else:
                                    res_spec = os.path.join(remote_dir,resolution,"ensemble-stats",species+"_"+stat)
  
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
                print(f"\n{experiment} experiment data to download ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range
            for remote_dir in res_spec_dir:
                if not initial_check:
                    local_path = remote_dir.split('/',7)[-1]
                    print(f"\n  - {os.path.join(self.exp_to_interp_root,local_path)}")
                         
                # get nc files
                nc_files = self.sftp.listdir(remote_dir)

                if nc_files:
                        # if it is an ensemble member
                        if not ensemble_options.startswith("stat_"):
                            # get the domain, resolution and species from the path
                            domain, resolution, species = remote_dir.split('/')[-3:]

                            # identify format of the directory
                            # the format is a tuple of how many - and how many _ are there
                            # the directory format is choosen by popularity
                            formats_list = [(file.count("-"), file.count("_")) for file in nc_files]
                            number_of_formats_dict = {format: formats_list.count(format) for format in set(formats_list)}
                            format = max(number_of_formats_dict, key=number_of_formats_dict.get)
                            
                            # filter and get only the files that follow the format
                            nc_files = list(filter(lambda x:(x.count("-"),x.count("_")) == format,nc_files))
                            
                            # example: od550du_2019040212.nc (0,1)
                            if format == (0,1):
                                # when there is no ensemble option in the name only allmembers and 000 are valid
                                if ensemble_options == '000' or ensemble_options == 'allmembers':
                                    nc_files = list(filter(lambda x:x.split("_")[0] == species, nc_files))
                                else:
                                    msg = f"There is no data available in {REMOTE_MACHINE} for the {exp_id} experiment with the {domain} domain with the {ensemble_options} ensemble option."
                                    show_message(self, msg, deactivate=initial_check)
                                    continue
                            # example: od550du-000_2021020812.nc (1,1)
                            elif format == (1,1):
                                # filter by ensemble option in case that ensemble option is not allmembers
                                if ensemble_options != 'allmembers':
                                    nc_files = list(filter(lambda x:x.split("_")[0] == species+'-'+ensemble_options,nc_files))
                                # if there is no options with the ensemble option, tell the user
                                if nc_files is []:
                                    msg = f"There is no data available in {REMOTE_MACHINE} for the {exp_id} experiment with the {domain} domain with the {ensemble_options} ensemble option."
                                    show_message(self, msg, deactivate=initial_check)
                                    continue
                            else:
                                # TODO delete this in the future
                                error = "It is not possible to download this nc file type yet. Please, contact the developers.", nc_files
                                sys.exit(error)
                        
                        # if it is an ensemble statistic
                        else:
                            # get the domain, resolution and species from the path
                            domain, resolution, _, species = remote_dir.split('/')[-4:]
                            species = species.split("_",1)[0]

                            # filter the nc files to only get the ones that have the correct species and stats
                            nc_files = list(filter(lambda x:x.split("_")[0] == species and "_".join(x[:-3].split("_")[2:]) == stat, nc_files))
                
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
                    if not ensemble_options.startswith("stat_"):
                        local_dir = os.path.join(self.exp_to_interp_root,exp_id,domain,resolution,species)
                    else:
                        local_dir = os.path.join(self.exp_to_interp_root,exp_id,domain,resolution,"ensemble-stats",species+"_"+stat)
                    
                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # sort nc_files
                    valid_nc_files.sort() 

                    if not initial_check:
                        # get the ones that are not already downloaded
                        valid_nc_files = list(filter(lambda x:os.path.join(local_dir,x) in files_to_download, valid_nc_files))
                        if not valid_nc_files:
                            msg = "Files were already downloaded."
                            show_message(self, msg, deactivate=initial_check)     
                            continue         
                        valid_nc_files_iter = tqdm(valid_nc_files, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                    else:
                        # do not print the bar if it is the initial check
                        valid_nc_files_iter = valid_nc_files

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files_iter:
                        local_path = os.path.join(local_dir,nc_file)
                        if initial_check:
                            initial_check_nc_files.append(local_path)
                        else:
                            self.latest_nc_file_path = local_path
                            remote_path = os.path.join(remote_dir,nc_file)
                            self.ncfile_dl_start_time = time.time()
                            self.sftp.get(remote_path,local_path,callback=check_time) 
            
            return initial_check_nc_files

        # tell the user if not valid resolution specie date combinations
        else:
            msg = "There are no available observations to be downloaded."
            show_message(self, msg, deactivate=initial_check)
            
    def fetch_zenodo_networks(self):
        # Get urls from zenodo to get GHOST zip files url

        # initialize dictionary to store possible networks
        self.zenodo_ghost_available_networks = {} 

        response = requests.get(self.ghost_url)

        # Check if the request was successful (status code 200)
        if response.status_code != 200:
            error = f'Failed to retrieve the webpage. Status code: {response.status_code}'
            sys.exit(error)
        
        # fill network dictionary with its corresponding zip url
        for line in response.text.split(">"):
            if '<link rel="alternate" type="application/zip" href=' in line:
                zip_file_url = line.split('href="')[-1][:-1]
                zip_network = line.split("/")[-1][:-5]
                self.zenodo_ghost_available_networks[zip_network] = zip_file_url

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
                sys.exit(error)

        if download_source in ["n","a"]:
            self.network = self.nonghost_available_networks
    
    def get_all_experiments(self):
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect()  

        # download all interpolated experiments
        if self.interpolated is True:
            # get directory content and format it as the experiments       
            experiment_list = self.sftp.listdir(os.path.join(self.exp_remote_path,self.ghost_version))
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
        
def main(**kwargs):
    """ Main function when running download function. """
    # initialise break blocker
    global download  
    download = ProvidentiaDownload(**kwargs)
    download.run()