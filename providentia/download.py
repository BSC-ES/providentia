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

# urlparse
from tqdm import tqdm
from remotezip import RemoteZip
import tarfile
from datetime import datetime, timedelta
from getpass import getpass

from .configuration import ProvConfiguration, load_conf
from .read_aux import check_for_ghost
from .warnings import show_message
# TODO delete when sure
# from .read_aux import get_ghost_observational_tree, get_nonghost_observational_tree

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)
REMOTE_MACHINE = "storage5"

data_paths = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/data_paths.yaml')))

def sighandler(*unused):
    print('\nKeyboard Interrupt. Stopping execution.')
    
    # close connection, if it exists
    if download.ssh != None:
        print("\nClosing ssh connection...")
        download.ssh.close()
        download.sftp.close()

    # delete temp dir if necessary
    temp_dir = os.path.join(download.ghost_root,'.temp')
    if os.path.exists(temp_dir):
        dot_temp_contents = os.listdir(temp_dir)
        if dot_temp_contents:
            print(f"\nDeleting {temp_dir} contents...")
            for path in dot_temp_contents:
                shutil.rmtree(os.path.join(temp_dir,path))
    
    # TODO delete when sure
    # download.update_filetrees()
    
    print("\nExiting...")
    sys.exit()

class ProvidentiaDownload(object):
    def __init__(self,**kwargs):
        # TODO move some of these variables to configuration.py
        # initialise zenodo url
        self.ghost_url = 'https://zenodo.org/records/10637450'

        # initialize dictionaries to store possible networks
        self.ghost_zip_files = {} # TODO: Maybe juntar este con el del zenodo

        # initialise remote hostname
        self.remote_hostname = "transfer1.bsc.es"

        # get ghost_version list
        self.possible_ghost_versions = os.listdir(os.path.join(CURRENT_PATH,'dependencies', 'GHOST_standards'))

        # in case transfer broke
        # global REMOTE_MACHINE
        # self.remote_hostname, REMOTE_MACHINE = "glogin4.bsc.es", "mn5" 

        self.ghost_remote_obs_path = data_paths[REMOTE_MACHINE]["ghost_root"]
        self.nonghost_remote_obs_path = data_paths[REMOTE_MACHINE]["nonghost_root"]
        self.exp_remote_path = data_paths[REMOTE_MACHINE]["exp_root"]

        # get ssh user and password 
        env = dotenv_values(os.path.join(PROVIDENTIA_ROOT, ".env"))

        # get ssh user and password 
        self.prv_user = env.get("PRV_USER")
        self.prv_password = env.get("PRV_PWD")

        # get ssh user and password 
        self.bsc_download_choice = env.get("BSC_DL_CHOICE")

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
            else:
                error = 'Error: The path to the configuration file specified in the command line does not exist.'
                sys.exit(error)
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            sys.exit(error)

        # create empty directories for all the paths if they don't exist
        for path in [self.nonghost_root,self.ghost_root,self.exp_root]:
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

    def run(self):
        for section_ind, section in enumerate(self.parent_section_names):
            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # update self with section variables
            for k, val in self.section_opts.items():
                setattr(self, k, self.provconf.parse_parameter(k, val))

            # now all variables have been parsed, check validity of those, throwing errors where necessary
            self.provconf.check_validity(deactivate_warning=True)

            # if networks is none, raise error
            if not self.network:
                error = "Error: No networks were passed."
                sys.exit(error)
            
            # when one of those symbols is passed, get all networks
            if self.network == ["*"] or self.network == ["default"]:
                self.get_all_networks()

            # from here generate control if user stopped execution
            signal.signal(signal.SIGINT, sighandler)

            # combine all networks and species combinations to download (for network and filter species)
            combined_networks = [(network, None) for network in self.network] + \
                                [(network_specie.split('|')[0], network_specie.split('|')[1]) for network_specie in self.filter_species]
            
            for network, filter_species in combined_networks:
                # change species when turn of filter species
                if filter_species != None:
                    self.species = [filter_species]

                # ghost
                if check_for_ghost(network):
                    if self.bsc_download_choice == 'y':
                        self.download_ghost_network_sftp(network)
                    elif self.bsc_download_choice == 'n':
                        self.download_ghost_network(network)
                    else:
                        error = "Download option not valid, check your .env file."
                        sys.exit(error)
                # non-ghost
                else:
                    self.download_nonghost_network(network)
            
            # when one of those symbols is passed, get all experiments
            if self.experiments == {'*': '*'} or self.experiments == {'default': 'default'}:
                self.get_all_experiments()

            # experiment
            for experiment in self.experiments.keys():
                self.download_experiment(experiment)
            
            # TODO delete when sure
            # update filetrees
            # self.update_filetrees()

        # close connection, if it exists
        if self.ssh != None:
            self.ssh.close() 
            self.sftp.close()

    def connect(self):
        # get public remote machine public key and add it to ssh object
        _, output = subprocess.getstatusoutput(f"ssh-keyscan -t ed25519 {self.remote_hostname}")
        ed25519_key = output.split()[-1].encode()
        key = paramiko.Ed25519Key(data=decodebytes(ed25519_key))

        self.ssh = paramiko.SSHClient()
        hostkeys = self.ssh.get_host_keys().add(self.remote_hostname, 'ed25519', key)

        # initialise temporal variables
        remind = False
        prv_user, prv_password = None, None

        # if couldn't get user, ask for it
        if self.prv_user == None:
            prv_user = ''
            while prv_user == '':
                prv_user = input(f"\nInsert BSC {REMOTE_MACHINE} ssh user: ")
            self.prv_user = prv_user
        
        # if couldn't get user, check if you have to ask for it
        if self.prv_password == None:
            # check if user needs a password
            try:
                self.ssh.connect(self.remote_hostname, username=self.prv_user, password='placeholder')
            # if authentication error, that means that the user or and the password are wrong
            except paramiko.ssh_exception.AuthenticationException:
                # if name was not changed, then user in .env is not valid
                if prv_user == None:
                    error = f"Authentication failed. Please, check if PRV_USER on {os.path.join(PROVIDENTIA_ROOT, '.env')} aligns with your BSC {REMOTE_MACHINE} ssh user."
                    error += "\nIf it does not, change the user to the correct one. If it does, delete the whole PRV_USER row and execute again."
                    sys.exit(error)
                else:
                    prv_password = getpass("Insert password: ")
                    self.prv_password = prv_password
        
        # if pwd or user changed, ask for credentials
        if prv_user != None or prv_password != None:
            # ask user if they want their credentials saved
            remind_txt = input("\nRemember credentials (y/n)? ")
            while remind_txt.lower() not in ['y','n']:
                remind_txt = input("\nRemember credentials (y/n)? ")

            remind = remind_txt.lower() == 'y'

        # catch identification method
        try:
            # connect through ssh and create Secure File Transfer Protocol object
            self.ssh.connect(self.remote_hostname, username=self.prv_user, password=self.prv_password)
            self.sftp = self.ssh.open_sftp()
            
            # create .env with the input user and/or password
            if remind:
                with open(os.path.join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                    if prv_user != None:
                        f.write(f"PRV_USER={self.prv_user}\n")
                    if prv_password != None:
                        f.write(f"PRV_PWD={self.prv_password}\n")

                print(f"\nRemote machine credentials saved on {os.path.join(PROVIDENTIA_ROOT, '.env')}")
        # if credentials are invalid, throw an error
        except paramiko.ssh_exception.AuthenticationException:
            error = "Authentication failed."
            # if user or password were taken from .env (did not change), tell the user to check .env
            if prv_user == None:
                error += f" Please, check your credentials on {os.path.join(PROVIDENTIA_ROOT, '.env')}"
            sys.exit(error)

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
                    
    def download_nonghost_network(self,network):
        # check if ssh exists and check if still active, connect if not
        if self.ssh == None or self.ssh.get_transport().is_active():
            self.connect() 

        # print current_network
        print('\n'+'-'*40)
        print(f"\nDownloading non-GHOST network: {network}...")

        # if not valid network, check if user put the network on init_prov 
        if network not in self.nonghost_available_networks:
            msg = f"The network {network} could not be found on {os.path.join(PROVIDENTIA_ROOT,'settings','init_prov.yaml')} nonghost_available_networks list."
            msg += "\nPlease, add the network to the list and execute again."
            show_message(self, msg)
            return 
        
        # check if nonghost network exists in directory
        # TODO: Change this to somewhere in configuration, the one up too
        try:
            self.sftp.stat(os.path.join(self.nonghost_remote_obs_path,network))
        except FileNotFoundError:
            msg = f"There is no data available for {network} network."
            show_message(self, msg)
            return 

        # check if all resolutions are in init_prov, if not warning and delete the not correct ones
        not_available_resolutions = set(self.resolution) - set(self.nonghost_available_resolutions)
        if not_available_resolutions:
            available_resolutions = set(self.resolution) - not_available_resolutions
            msg = f"The resolution/s {', '.join(not_available_resolutions)} could not be found on {os.path.join(PROVIDENTIA_ROOT,'settings','init_prov.yaml')} nonghost_available_resolutions list."
            msg += "\nPlease, add the necessary resolutions to the list and execute again."
            show_message(self, msg)
            return 

        # get resolution and species combinations
        res_spec_dir = []

        sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(os.path.join(self.nonghost_remote_obs_path,network))).intersection(self.ghost_available_resolutions)
        for resolution in sftp_resolutions:
            try:
                sftp_species  = self.species if self.species else self.sftp.listdir(os.path.join(self.nonghost_remote_obs_path,network,resolution))
            except FileNotFoundError:
                msg = f"There is no data available for {network} network at {resolution} resolution"
                show_message(self, msg)
                continue
            for species in sftp_species: 
                res_spec_dir.append(os.path.join(self.nonghost_remote_obs_path,network,resolution,species))
        
        # print the species, resolution and network combinations that are going to be downloaded
        if res_spec_dir:
            print(f"\n{network} observations to download: ({len(res_spec_dir)})")
            for combi_print in res_spec_dir:
                print(f"  - {os.path.join(self.nonghost_root,combi_print.split('/',6)[-1])}")

            valid_res_spec_dir_nc_files = []
            # get all the nc files in the date range within the specie and resolution combination
            for remote_dir in res_spec_dir:
                local_dir = os.path.join(self.nonghost_root,remote_dir.split('/',6)[-1])
                species = remote_dir.split('/')[-1]
                resolution = remote_dir.split('/')[-2]

                try:
                    nc_files = self.sftp.listdir(remote_dir)
                except FileNotFoundError:
                    msg = f"There is no data available for {network} network {species} species at {resolution} resolution"
                    show_message(self, msg)
                    continue

                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if network + species + resolution + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available from {self.start_date} to {self.end_date} for {network} network {species} species at {resolution} resolution."
                    show_message(self, msg)
                    continue
                
                unique_valid_nc_files = copy.deepcopy(valid_nc_files)
                valid_res_spec_dir_nc_files.append((remote_dir,local_dir,unique_valid_nc_files))
            
            print()

            # download the valid resolution specie date combinations
            if valid_res_spec_dir_nc_files:
                for remote_dir,local_dir,valid_nc_files in tqdm(valid_res_spec_dir_nc_files,ascii=True, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"Downloading valid observations ({len(valid_res_spec_dir_nc_files)})"):
                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files:
                        remote_path = os.path.join(remote_dir,nc_file)
                        local_path = os.path.join(local_dir,nc_file)
                        self.sftp.get(remote_path,local_path)

                print(f"\n{network} observations downloaded: ({len(valid_res_spec_dir_nc_files)})")
                for _,local_dir,_ in valid_res_spec_dir_nc_files:
                    print(f"  - {os.path.join(local_dir)}")
            
            # tell the user if not valid resolution specie date combinations
            else:
                print("There are no valid observations to be downloaded.")
        
    def download_ghost_network_sftp(self,network):
        # check if ssh exists and check if still active, connect if not
        if self.ssh == None or self.ssh.get_transport().is_active():
            self.connect() 

        # print current_network
        print('\n'+'-'*40)
        print(f"\nDownloading GHOST network: {network} from {REMOTE_MACHINE}...")

        # If not valid network, next
        if network not in self.sftp.listdir(self.ghost_remote_obs_path):
            msg = f"There is no data available for {network} network."
            show_message(self, msg)
            return 
        
        # If not valid combination of ghost version and network, next 
        elif self.ghost_version not in self.sftp.listdir(os.path.join(self.ghost_remote_obs_path,network)):
            msg = f"There is no data available for {network} network for the current ghost version ({self.ghost_version})."
            
            available_ghost_versions = set(self.sftp.listdir(os.path.join(self.ghost_remote_obs_path,network))).intersection(self.possible_ghost_versions)

            if available_ghost_versions:
                msg += f" Please check one of the available versions: {', '.join(sorted(available_ghost_versions))}"
            else:
                msg += f" There are no other versions available at the moment."
            
            show_message(self, msg)
            return

        # get resolution and species combinations
        res_spec_dir = []

        remote_dir = os.path.join(self.ghost_remote_obs_path,network,self.ghost_version)

        sftp_resolutions = self.resolution if self.resolution else set(self.sftp.listdir(remote_dir)).intersection(self.ghost_available_resolutions)
        for resolution in sftp_resolutions:
            try:
                sftp_species  = self.species if self.species else self.sftp.listdir(os.path.join(remote_dir,resolution))
            except FileNotFoundError:
                msg = f"There is no data available for {network} network at {resolution} resolution"
                show_message(self, msg)
                continue
            for species in sftp_species: 
                res_spec_dir.append(os.path.join(remote_dir,resolution,species))
        
        # print the species, resolution and network combinations that are going to be downloaded
        if res_spec_dir:
            print(f"\n{network} observations to download: ({len(res_spec_dir)})")
            for combi_print in res_spec_dir:
                print(f"  - {os.path.join(self.ghost_root,combi_print.split('/',7)[-1])}")

            valid_res_spec_dir_nc_files = []
            # get all the nc files in the date range within the specie and resolution combination
            for remote_dir in res_spec_dir:
                local_dir = os.path.join(self.ghost_root,remote_dir.split('/',7)[-1])
                species = remote_dir.split('/')[-1]
                resolution = remote_dir.split('/')[-2]

                try:
                    nc_files = self.sftp.listdir(remote_dir)
                except FileNotFoundError:
                    msg = f"There is no data available for {network} network {species} species at {resolution} resolution"
                    show_message(self, msg)
                    continue

                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if network + species + resolution + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available from {self.start_date} to {self.end_date} for {network} network {species} species at {resolution} resolution."
                    show_message(self, msg)
                    continue
                
                unique_valid_nc_files = copy.deepcopy(valid_nc_files)
                valid_res_spec_dir_nc_files.append((remote_dir,local_dir,unique_valid_nc_files))
            
            print()

            # download the valid resolution specie date combinations
            if valid_res_spec_dir_nc_files:
                for remote_dir,local_dir,valid_nc_files in tqdm(valid_res_spec_dir_nc_files,ascii=True, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"Downloading valid observations ({len(valid_res_spec_dir_nc_files)})"):
                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files:
                        remote_path = os.path.join(remote_dir,nc_file)
                        local_path = os.path.join(local_dir,nc_file)
                        self.sftp.get(remote_path,local_path)

                print(f"\n{network} observations downloaded: ({len(valid_res_spec_dir_nc_files)})")
                for _,local_dir,_ in valid_res_spec_dir_nc_files:
                    print(f"  - {os.path.join(local_dir)}")

            # tell the user if not valid resolution specie date combinations
            else:
                print("There are no valid observations to be downloaded.")

    def download_ghost_network(self,network):
        # print current_network
        print('\n'+'-'*40)
        print(f"\nDownloading GHOST network: {network} from Zenodo...")

        # if first time reading a ghost network, get current zips urls in zenodo page
        if not self.ghost_zip_files: 
            self.get_ghost_zip_files()

        if network not in self.ghost_zip_files:
            msg = f"There is no data available for {network} network."
            show_message(self, msg)
            return

        # get url to download the zip file for the current network
        zip_url = self.ghost_zip_files[network]

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
        res_spec_final_combinations = []
        with RemoteZip(zip_url) as zip:
            for combi in res_spec_combinations:
                res_spec_final_combinations += list(filter(lambda x: combi in x, zip.namelist()))
            res_spec_final_combinations = list(filter(lambda x: x[-7:] == '.tar.xz', res_spec_final_combinations))
            
            # warning if network + species + resolution combination is gets no matching results
            # TODO change warning as it is on sftp mode
            if not res_spec_final_combinations:
                print_spec = f'{",".join(self.species)} species' if self.species else ""
                print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
                msg = f"There is no data available for {network} network {print_spec} {print_res}"
                show_message(self, msg)

            # check if there's any possible combination with user's network, resolution and species
            else:
                # Print the species, resolution and network combinations that are going to be downloaded
                print(f"\n{network} observations to download:")
                last_ghost_version = 1.5 #TODO In the future there will be various versions in zenodo
                for species in res_spec_final_combinations:
                    hourly_specie_print = "/".join(species.split("/")[1:])
                    print(f"  - {self.ghost_root}/{network}/{last_ghost_version}/{hourly_specie_print[:-7]}")
                
                # create temporal dir to store the middle tar file with its directories
                temp_dir = os.path.join(self.ghost_root,'.temp')
                if not os.path.exists(temp_dir):
                    os.mkdir(temp_dir)

                # extract species from zip files
                for specie_to_get in tqdm(res_spec_final_combinations,desc="Downloading Observations"):
                    zip.extract(specie_to_get,temp_dir)
                    
                    # get path and the name of the directory of the tar file
                    tar_path = os.path.join(self.ghost_root, network, str(last_ghost_version), "/".join(specie_to_get.split("/")[1:]))
                    temp_path = os.path.join(temp_dir,specie_to_get)
                    
                    tar_dir = os.path.dirname(tar_path) 

                    if not os.path.exists(tar_dir):
                        os.makedirs(tar_dir) #TODO maybe add the chmod

                    # extract nc file from tar file
                    with tarfile.open(temp_path) as tar_file:
                        tar_names = tar_file.getnames()

                        # get the nc files that are between the start and end date
                        valid_nc_file_names = self.get_valid_nc_files_in_date_range(tar_names)
                        valid_nc_files = [tar_member for tar_member in tar_file.getmembers() if tar_member.name in valid_nc_file_names]

                        tar_file.extractall(path = tar_dir, members = valid_nc_files)
                    
                    # remove the temp directory
                    shutil.rmtree(os.path.dirname(os.path.dirname(temp_path)))

                    # warning if network + species + resolution + date range combination gets no matching results                        
                    if not valid_nc_files:
                        print_spec = f'{",".join(self.species)} species' if self.species else ""
                        print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
                        msg = f"There is no data available from {self.start_date} to {self.end_date} for {network} network {print_spec} {print_res}"
                        show_message(self, msg)
                            
    def download_experiment(self, experiment):
        # check if ssh exists and check if still active, connect if not
        if self.ssh == None or self.ssh.get_transport().is_active():
            self.connect()  

        # print current experiment
        print('\n'+'-'*40)
        print(f"\nDownloading Experiment: {experiment} from {REMOTE_MACHINE}...")
        
        # get resolution and species combinations
        res_spec_dir = []

        remote_dir = os.path.join(self.exp_remote_path,self.ghost_version,experiment)

        # TODO, should check if experiment is checked
        if experiment not in self.sftp.listdir(os.path.join(self.exp_remote_path,self.ghost_version)):
            msg = f"There is no data available for {experiment} experiment for the current ghost version ({self.ghost_version})."
            
            possible_ghost_versions = set(self.sftp.listdir(self.exp_remote_path)).intersection(set(self.possible_ghost_versions))
            available_ghost_versions = list(filter(lambda x:experiment in self.sftp.listdir(os.path.join(self.exp_remote_path,x)),possible_ghost_versions))
            
            if available_ghost_versions:
                msg += f" Please check one of the available versions: {', '.join(sorted(available_ghost_versions))}"
            else:
                msg += f" There are no other versions available at the moment."
            show_message(self, msg)
            return

        sftp_resolutions = self.resolution if self.resolution else self.sftp.listdir(remote_dir)
        for resolution in sftp_resolutions:
            try:
                sftp_species  = self.species if self.species else self.sftp.listdir(os.path.join(remote_dir,resolution))
            except FileNotFoundError:
                msg = f"There is no data available for {experiment} experiment at {resolution} resolution"
                show_message(self, msg)
                continue
            for species in sftp_species: 
                try:
                    sftp_network = self.network if self.network else self.sftp.listdir(os.path.join(remote_dir,resolution,species))
                except FileNotFoundError:
                    msg = f"There is no data available for {experiment} experiment {species} species at {resolution} resolution"
                    show_message(self, msg)
                    continue
                for network in sftp_network:
                    res_spec_dir.append(os.path.join(remote_dir,resolution,species,network))
        
        # print the species, resolution and experiment combinations that are going to be downloaded
        if res_spec_dir:
            print(f"\n{experiment} experiments to download: ({len(res_spec_dir)})")
            for combi_print in res_spec_dir:
                print(f"  - {os.path.join(self.exp_root,combi_print.split('/',7)[-1])}")

            valid_res_spec_dir_nc_files = []
            # get all the nc files in the date range
            for remote_dir in res_spec_dir:
                local_dir = os.path.join(self.exp_root,remote_dir.split('/',7)[-1])

                network = remote_dir.split('/',11)[-1]
                species = remote_dir.split('/',11)[-2]
                resolution = remote_dir.split('/',11)[-3]

                try:
                    nc_files = self.sftp.listdir(remote_dir)
                except FileNotFoundError:
                    msg = f"There is no data available for {experiment} experiment {species} species {network} network at {resolution} resolution"
                    show_message(self, msg)
                    continue

                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if experiment + species + resolution + network + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available from {self.start_date} to {self.end_date} for {experiment} experiment {species} species {network} network at {resolution} resolution"
                    show_message(self, msg)
                    continue

                unique_valid_nc_files = copy.deepcopy(valid_nc_files)
                valid_res_spec_dir_nc_files.append((remote_dir,local_dir,unique_valid_nc_files))
            
            print()

            # download the valid resolution specie date combinations
            if valid_res_spec_dir_nc_files:
                for remote_dir,local_dir,valid_nc_files in tqdm(valid_res_spec_dir_nc_files,ascii=True, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"Downloading valid experiments ({len(valid_res_spec_dir_nc_files)})"):
                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files:
                        remote_path = os.path.join(remote_dir,nc_file)
                        local_path = os.path.join(local_dir,nc_file)
                        self.sftp.get(remote_path,local_path)

                print(f"\n{experiment} experiments downloaded: ({len(valid_res_spec_dir_nc_files)})")
                for _,local_dir,_ in valid_res_spec_dir_nc_files:
                    print(f"  - {os.path.join(local_dir)}")

            # tell the user if not valid resolution specie date combinations
            else:
                print("There are no valid experiments to be downloaded.")

        # tell the user if not valid resolution specie date combinations
        else:
            print("There are no valid experiments to be downloaded.")

    def get_ghost_zip_files(self):
        # Get urls from zenodo to get ghost zip files url
        # Send an HTTP GET request to the URL
        # TODO CHANGE METHOD NAME
        response = requests.get(self.ghost_url)

        # Check if the request was successful (status code 200)
        if response.status_code != 200:
            msg = f'Failed to retrieve the webpage. Status code: {response.status_code}'
            show_message(self,msg)
        
        # fill network dictionary with its corresponding zip url
        for line in response.text.split(">"):
            if '<link rel="alternate" type="application/zip" href=' in line:
                zip_url = line.split('href="')[-1][:-1]
                zip_network = line.split("/")[-1][:-5]
                self.ghost_zip_files[zip_network] = zip_url

    def get_all_networks(self):
        # get user input to know which kind of network wants
        download_source = None
        while not download_source:
            download_source = input("\nDo you want to download GHOST, non-GHOST or all networks? (g/n/a) ").lower()
            download_source = download_source if download_source in ["g","n","a"] else None

        if download_source in ["g","a"]:
            if self.bsc_download_choice == 'y':
                self.network = self.ghost_available_networks
            elif self.bsc_download_choice == 'n':
                if not self.ghost_zip_files: 
                    self.get_ghost_zip_files()
                    self.network = list(self.ghost_zip_files.keys())
            else:
                error = "Download option not valid, check your .env file."
                sys.exit(error)

        if download_source in ["n","a"]:
            self.network = self.nonghost_available_networks
    
    def get_all_experiments(self):
        # check if ssh exists and check if still active, connect if not
        if self.ssh == None or self.ssh.get_transport().is_active():
            self.connect()  

        # get directory content and format it as the experiments       
        experiment_list = self.sftp.listdir(os.path.join(self.exp_remote_path,self.ghost_version))
        self.experiments = dict(zip(experiment_list,experiment_list))

    def get_valid_nc_files_in_date_range(self, nc_files):
        valid_nc_files = []
        for nc_file in nc_files:
            if ".nc" in nc_file:
                ym = nc_file.split("_")[-1].split(".nc")[0]
                if int('{}01'.format(ym)) >= int(self.start_date) and int('{}01'.format(ym)) < int(self.end_date):
                    valid_nc_files.append(nc_file)
                    
        return valid_nc_files        

    # TODO: delete when sure
    # def update_filetrees(self):
    #     print('\nUpdating filetrees...')
    #     get_ghost_observational_tree(download)
    #     get_nonghost_observational_tree(download)
        
def main(**kwargs):
    """ Main function when running download function. """
    # initialise break blocker
    global download  
    download = ProvidentiaDownload(**kwargs)
    download.run()