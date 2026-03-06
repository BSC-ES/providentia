""" Class for Providentia download mode """

from base64 import decodebytes
import copy
from getpass import getpass
import os
import shutil
import signal
import subprocess
import sys
import time

from dotenv import dotenv_values
import paramiko 
from tqdm import tqdm
import yaml 

from .actris import Actris
from .cams import (Cams, cams_options)
from providentia.auxiliar import CURRENT_PATH, join
from .configuration import ProvConfiguration, load_conf
from .read_aux import check_for_ghost
from .warnings_prv import show_message
from .zenodo import Zenodo

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

# load the defined models paths, agrupations yaml and mapping species
data_paths = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/data_paths.yaml')))
interp_models = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'interp_models.yaml')))
mapping_species =  yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'mapping_species.yaml')))
dl_hpc = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'dl_hpc.yaml')))
temporal_resolution_map = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'temporal_resolution_map.yaml')))

class Download(object):
    """Main class responsible for handling Providentia data downloads."""

    def __init__(self, **kwargs):
        """
        Initialize Providentia download mode.

        Parameters
        ----------
        **kwargs : dict
            Optional command-line arguments that override default configuration values.
        """

        # update self with command line arguments
        self.commandline_arguments = copy.deepcopy(kwargs)

        # get providentia start time
        self.prov_start_time = time.time()

        # get ssh user and password 
        env = dotenv_values(join(PROVIDENTIA_ROOT, ".env"))

        # get ssh user and password 
        self.prv_user = env.get("PRV_USER")
        self.prv_password = env.get("PRV_PWD")

        # get origin update (ACTRIS)
        self.origin_update_choice = env.get("ORIGIN_UPDATE")

        # initialise default configuration variables
        # modified by commandline arguments, if given
        self.provconf = ProvConfiguration(self, **self.commandline_arguments)

        self.logger.info("Starting Providentia download...")

        # update variables from config file
        if self.config != '':  
            read_conf = False
            if os.path.exists(self.config):
                read_conf = True
            elif os.path.exists(join(self.config_dir, self.config)):
                    self.config = join(self.config_dir, self.config)
                    read_conf = True

            if read_conf:
                load_conf(self, self.config)
                self.from_conf = True
                # get the section in case it was passed on the command line                
                if 'section' in self.commandline_arguments:
                    # config and section defined 
                    if self.commandline_arguments['section'] in self.all_sections:
                        self.sections = [self.commandline_arguments['section']]
                    else:
                        msg = 'Error: The section specified in the command line ({0}) does not exist.'.format(self.commandline_arguments['section'])
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

        # variable that saves whether some models/observations were downloaded before
        self.overwritten_files_flag = False

        # initialize the necessary things in local
        if self.machine == "local":

            # index to know which remote machine has been tried already
            self.switched_remote = 0

            # initialise remote hostname and machine
            self.remote_hostname, self.remote_machine = dl_hpc[self.switched_remote]

            # initialise ssh 
            self.ssh = None

    def run(self):
        """Execute the Providentia download process for all configured sections."""
        
        for section_ind, section in enumerate(self.sections):
            self.logger.info('Starting to download for {} section'.format(section))

            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # update self with section variables (if not passed via command line)
            for k, val in self.section_opts.items():
                if k not in self.commandline_arguments:
                    setattr(self, k, self.provconf.parse_parameter(k, val))
            
            # now all variables have been parsed, check validity of those, throwing errors where necessary
            self.provconf.check_validity()

            # TODO: make it work directly with the asterisk
            # transform asterisk fields to empty since it was originally coded this way
            for field in ['species', 'resolution']:
                if getattr(self, field) == ['*']:
                    setattr(self, field, '')

            # TODO: remove the filters instead
            # transform asterisk fields to a low and high date
            for field, date_num in {'start_date' : '0', 'end_date' : '9'}.items():
                if getattr(self, field) == '*':
                    setattr(self, field, date_num * 8)
            
            # from here generate control if user stopped execution
            signal.signal(signal.SIGINT, self.sighandler)
            
            # only the local download iterates through the networks
            if self.machine in "local":                        
                # networks
                if self.network and self.dl_mode != 'mod':

                    if self.network == ["*"]:
                        if self.network_type:
                            # exit if the value is neither BSC or zenodo
                            network_type = str(self.network_type).lower()
                            if network_type not in ['ghost', 'non-ghost']:
                                error = f"Error: Invalid 'network_type': '{self.network_type}'. Expected 'ghost' or 'non-ghost'."
                                self.logger.error(error)
                                sys.exit(1)

                            self.reading_ghost = network_type == 'ghost'  
                        else:
                            # get user input to know which kind of network wants
                            while True:
                                download_source = input("\nDo you want to download all the GHOST networks? (Otherwise all the non-GHOST networks will be downloaded) ([y]/n): ").lower()
                                if download_source in ['', 'y', 'n']:
                                    break
                            
                            # get the boolean value from the answer of the user
                            self.reading_ghost = download_source in ['', 'y']

                    # if there are GHOST networks, ask the user whether they want to download it from zenodo or HPC machines
                    if self.reading_ghost:
                        if not self.dl_ghost_source:
                            # ask whether the user wants to download from the zenodo or bsc machine
                            while True: 
                                dl_ghost_source = input("\nDo you want to download observational data from the BSC remote machine? (Otherwise, GHOST observational data will be retrieved from Zenodo) ([y]/n): ").lower()
                                if dl_ghost_source in ['', 'y', 'n']:
                                    break

                            self.dl_ghost_source = 'bsc' if dl_ghost_source in ['', 'y'] else 'zenodo'
                        else:
                            # exit if the value is neither bsc or zenodo
                            if str(self.dl_ghost_source).lower() not in ['bsc', 'zenodo']:
                                error = f"Error: Invalid 'dl_ghost_source': '{self.dl_ghost_source}'. Expected 'bsc' or 'zenodo'."
                                self.logger.error(error)
                                sys.exit(1)

                            self.dl_ghost_source = self.dl_ghost_source.lower()

                        # initialise the Zenodo object if user chose a Zenodo download
                        if self.dl_ghost_source == 'zenodo':
                            # warn the user about Zenodo download speed
                            msg = "Downloading from Zenodo can take a little time, please be patient."
                            show_message(self, msg)
                            
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
                    self.download_ghost_network_sftp if self.reading_ghost and self.dl_ghost_source == 'bsc'
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
                                actris.download_actris_data()
                        # GHOST and non-GHOST 
                        else:
                            # download GHOST network
                            initial_check_nc_files = download_fun(network, initial_check=True)
                            files_to_download = self.select_files_to_download(initial_check_nc_files)
                            if not initial_check_nc_files or files_to_download:
                                download_fun(network, initial_check=False, files_to_download=files_to_download)

                    # get orignal species back
                    self.species = main_species

            # when one of those symbols is passed, get all models
            if self.experiments == {'*' : '*'}:
                self.get_all_models()

            if self.experiments and self.dl_mode != 'obs':
                # remote machine model download
                if self.machine in ["storage5", "nord4"]:
                    # get function to download model depending on the configuration file field
                    download_model_fun = self.copy_non_interpolated_model
                # local model download
                else:
                    for model in self.experiments:
                        # CAMS model
                        if model.startswith(tuple(cams_options.keys())):
                            self.cams = Cams(self)
                            initial_check_nc_files = self.cams.download_cams_model(model, initial_check=True)
                            files_to_download = self.select_files_to_download(initial_check_nc_files)
                            if not initial_check_nc_files or files_to_download:
                                self.cams.download_cams_model(model, initial_check=False, files_to_download=files_to_download)
                        # BSC machines
                        else:
                            download_model_fun = self.download_model if self.dl_interpolated else self.download_non_interpolated_model
                    
                            # iterate the models download
                            for model in self.experiments.keys():
                                initial_check_nc_files = download_model_fun(model, initial_check=True)
                                files_to_download = self.select_files_to_download(initial_check_nc_files)
                                if not initial_check_nc_files or files_to_download:
                                    download_model_fun(model, initial_check=False, files_to_download=files_to_download)

            # remove section variables from memory
            for k in self.section_opts:
                try:
                    vars(self).pop(k)
                except:
                    pass

            # reset domain and ensemble for new section
            self.domain = []
            self.ensemble = []

            # reinitialise default configuration variables
            # modified by commandline arguments, if given
            self.provconf = ProvConfiguration(self, **self.commandline_arguments)

        # show message in case models or observations were ignored
        if self.overwritten_files_flag == True:
            self.logger.info("\nSome models/observations were found but were not downloaded because the user chose not to overwrite or because 'dl_overwrite' is set to False.")

        if self.machine == "local":
            # close connection, if it exists
            if self.ssh is not None:
                self.ssh.close() 
                self.sftp.close()
        
    def connect(self): 
        """Establish an SSH connection to a remote BSC machine for downloading data."""

        # initialise the paths
        self.ghost_remote_obs_path = data_paths[self.remote_machine]["ghost_root"]
        self.nonghost_remote_obs_path = data_paths[self.remote_machine]["nonghost_root"]
        self.mod_remote_path = data_paths[self.remote_machine]["mod_root"]
        self.mod_to_interp_remote_path = data_paths[self.remote_machine]["mod_to_interp_root"]

        # get public remote machine public key and add it to ssh object
        _, output = subprocess.getstatusoutput(f"ssh-keyscan -t ed25519 {self.remote_hostname}")
        
        # encode the output public key if possible
        try:
            ed25519_key = output.split()[-1].encode()
            key = paramiko.Ed25519Key(data=decodebytes(ed25519_key))
        
        # in case transfer broke
        except:
            # if the remote machine has not been changed
            if self.switched_remote < len(dl_hpc) - 1:

                # change remote machine and hostname
                self.switched_remote += 1
                previous_remote_hostname = self.remote_hostname
                self.remote_hostname, self.remote_machine = dl_hpc[self.switched_remote]
                
                msg = f"Remote machine {previous_remote_hostname} not working right now. Changing it to {self.remote_hostname}..."
                show_message(self, msg)
 
                # connect with the new machine
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
        counter = 0

        # iterate through the different user password combinations
        while True:

            # ask for user if not in .env
            if self.prv_user is None:
                prv_user = ''
                while prv_user == '':
                    prv_user = input(f"\nInsert BSC {self.remote_machine} ssh user: ")
                self.prv_user = prv_user
            
            # ask for password if not in .env
            if self.prv_password is None:
                try:
                    self.ssh.connect(self.remote_hostname, username=self.prv_user, password='placeholder')
                # if authentication error, that means that the user or and the password are wrong
                except paramiko.ssh_exception.AuthenticationException:
                    # if name was not changed, then user in .env is not valid
                    if prv_user is None:
                        error = f"Authentication failed. Please, check if PRV_USER on {join(PROVIDENTIA_ROOT, '.env')} aligns with your BSC {self.remote_machine} ssh user."
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

                # exit the loop when the connection was succesful
                break
                
            # if credentials are invalid, notify the user and retry
            except paramiko.ssh_exception.AuthenticationException:
                msg = "Authentication failed. Please re-enter your credentials."
                show_message(self, msg)

                self.prv_user, self.prv_password = None, None
                counter += 1

                if counter == 2:
                    error = "Error: Maximum number of authentication attempts reached. Please close and reopen the terminal, then try again"
                    self.logger.error(error)
                    sys.exit(1)

        # if pwd or user changed, ask if user wants to remember credentials
        if (prv_user is not None) or (prv_password is not None):
            # ask user if they want their credentials saved
            while True:
                remind_txt = input("\nRemember credentials ([y]/n)? ").lower()
                if remind_txt in ['y', 'n', '']:
                    break
            
            # create .env with the input user and/or password
            if remind_txt in ['y', '']:
                with open(join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                    if prv_user is not None:
                        f.write(f"PRV_USER={self.prv_user}\n")
                    if prv_password is not None:
                        f.write(f"PRV_PWD={self.prv_password}\n")

                self.logger.info(f"\nRemote machine credentials saved on {join(PROVIDENTIA_ROOT, '.env')}")
                    
    def select_files_to_download(self, nc_files_to_download):
        """
        Checks a list of files against the local filesystem and returns
        the subset of files that are not already downloaded. 

        Parameters
        ----------
        nc_files_to_download : list of str
            A list of file paths intended for download.

        Returns
        -------
        not_downloaded_files : list of str
            A list of file paths that are either not present locally or should be 
            re-downloaded depending on the `dl_overwrite` attribute.
        """
        
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
                    while True:
                        dl_overwrite = input("\nThere are some files that were already downloaded in a previous download, do you want to overwrite them ([y]/n)? ").lower() 
                        if dl_overwrite in ['y', 'n', '']:
                            break
                    
                    # get the boolean value
                    self.dl_overwrite = dl_overwrite != 'n'

                # if user wants to overwrite then add the files downloaded before the execution as if they were never downloaded
                if self.dl_overwrite:
                    not_downloaded_files += downloaded_before_execution_files
                # change overwritten files boolean to True to indicate that some files were ignored
                else:
                    self.overwritten_files_flag = True

        return not_downloaded_files

    def download_nonghost_network(self, network, initial_check, files_to_download=None):
        """
        Download non-GHOST network data from a remote machine via SFTP.

        Parameters
        ----------
        network : str
            Name of the non-GHOST network to download.
        initial_check : bool
            If True, only performs a check of available files without downloading.
            If False, downloads files and displays progress.
        files_to_download : list of str, optional
            Specific file paths to download. Only files in this list are considered.

        Returns
        -------
        initial_check_nc_files : list of str
            A list of file paths intended for download.
        """
                
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect() 
        
        if not initial_check:
            # print current_network
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nDownloading non-GHOST {network} network data from {self.remote_machine}...")

        # if not valid network, check if user put the network on init.yaml 
        if network not in self.nonghost_available_networks:
            msg = f"The {network} network could not be found on {join(PROVIDENTIA_ROOT,'settings','available_inputs.yaml')} nonghost_available_networks list."
            msg += "\nPlease, add the network to the list and execute again."
            show_message(self, msg, deactivate=initial_check)
            return
        
        # check if nonghost network exists in directory
        try:
            self.sftp.stat(join(self.nonghost_remote_obs_path, network))
        except FileNotFoundError:
            msg = f"There is no data available in {self.remote_machine} for {network} network."
            show_message(self, msg, deactivate=initial_check)
            return

        # check if all resolutions are in init.yaml, if not warning and delete the not correct ones
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
                msg = f"There is no data available in {self.remote_machine} for {network} network at {resolution} resolution."
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
                    self.logger.info(f"\n  - {local_dir}, source: {remote_dir} ({self.remote_machine})")

                try:
                    nc_files = self.sftp.listdir(remote_dir)
                except FileNotFoundError:
                    msg = f"There is no data available in {self.remote_machine} for {network} network for {species} species at {resolution} resolution."
                    show_message(self, msg, deactivate=initial_check)
                    continue
                
                # get the nc files in the date range
                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if network + species + resolution + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available in {self.remote_machine} from {self.start_date} to {self.end_date} for {network} network {species} species at {resolution} resolution."
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
        """
        Download GHOST network data from a remote machine via SFTP.

        Parameters
        ----------
        network : str
            Name of the GHOST network to download.
        initial_check : bool
            If True, only performs a check of available files without downloading.
            If False, downloads files and displays progress.
        files_to_download : list of str, optional
            Specific file paths to download. Only files in this list are considered.

        Returns
        -------
        initial_check_nc_files : list of str
            A list of file paths intended for download.
        """
        
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect() 

        if not initial_check:
            # print current_network
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nDownloading GHOST {network} network data from {self.remote_machine}...")

        # if not valid network, next
        if network not in self.sftp.listdir(self.ghost_remote_obs_path):
            msg = f"There is no data available in {self.remote_machine} for {network} network."
            show_message(self, msg, deactivate=initial_check)
            return 
        
        # if not valid combination of GHOST version and network, next 
        elif self.ghost_version not in self.sftp.listdir(join(self.ghost_remote_obs_path,network)):
            msg = f"There is no data available in {self.remote_machine} for {network} network for the current GHOST version ({self.ghost_version})."
            
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
                msg = f"There is no data available in {self.remote_machine} for {network} network at {resolution} resolution."
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
                    self.logger.info(f"\n  - {local_dir}, source: {remote_dir} ({self.remote_machine})")

                try:
                    nc_files = self.sftp.listdir(remote_dir)
                except FileNotFoundError:
                    msg = f"There is no data available in {self.remote_machine} for {network} network for {species} species at {resolution} resolution."
                    show_message(self, msg, deactivate=initial_check)
                    continue
                
                # get the nc files in the date range
                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if network + species + resolution + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available in {self.remote_machine} from {self.start_date} to {self.end_date} for {network} network {species} species at {resolution} resolution."
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

    def download_model(self, model, initial_check, files_to_download=None):
        """
        Download interpolated model data from a remote machine via SFTP.

        Parameters
        ----------
        model : str
            Name of the model to download.
        initial_check : bool
            If True, only performs a check of available files without downloading.
            If False, downloads files and displays progress.
        files_to_download : list of str, optional
            Specific file paths to download. Only files in this list are considered.

        Returns
        -------
        initial_check_nc_files : list of str
            A list of file paths intended for download.
        """
         
        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect()  
        
        if not initial_check:
            # print current model
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nDownloading {model} model data from {self.remote_machine}...")
            
        # get resolution and species combinations
        res_spec_dir = []
        
        # domain and ensemble are part of the model name, all united by dash (-)
        model_new = model

        # domain and ensemble are directories
        model_old = model.replace("-","/")
        
        # get remote directory format depending on the GHOST version
        model = model_old if self.ghost_version in ["1.2", "1.3", "1.3.1"] else model_new

        # get remote directory
        remote_dir = join(self.mod_remote_path,self.ghost_version,model)

        # check if model exists
        try:
            self.sftp.stat(remote_dir)
        except FileNotFoundError:
            msg = f"There is no data available in {self.remote_machine} for {model_new} model for the current GHOST version ({self.ghost_version})."

            # get possible GHOST versions from the combination of GHOST_standards and the real avaibles in the model remote machine path
            possible_ghost_versions = set(self.sftp.listdir(self.mod_remote_path)).intersection(set(self.possible_ghost_versions))
            
            # get available models in other GHOST versions (considering different formats)
            available_ghost_versions = []

            for possible_ghost_version in possible_ghost_versions:
                try:
                    # get model path depending on the GHOST version
                    remote_dir_ghost_version = join(self.mod_remote_path, possible_ghost_version, model_old if possible_ghost_version in ["1.2", "1.3", "1.3.1"] else model_new)
                    
                    # check if directory exists
                    self.sftp.stat(remote_dir_ghost_version)

                    # if it doesn't break, the model exists in this version
                    available_ghost_versions.append(possible_ghost_version)

                except FileNotFoundError:
                    continue
                            
            # list that saves the GHOST versions with valid nc files
            valid_available_ghost_versions = []
            
            # check for combinations of species, resolutions, network, and day in the available versions
            if available_ghost_versions:

                # iterate the different GHOST versions
                for possible_ghost_version in available_ghost_versions:
                    remote_dir_ghost_version = join(self.mod_remote_path, possible_ghost_version, model_old if possible_ghost_version in ["1.2", "1.3", "1.3.1"] else model_new)
                    
                    # iterate the different resolutions
                    sftp_resolutions = self.model_resolution if self.model_resolution else set(self.sftp.listdir(remote_dir)).intersection(self.nonghost_available_resolutions)
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

        sftp_resolutions = self.model_resolution if self.model_resolution else set(self.sftp.listdir(remote_dir)).intersection(self.nonghost_available_resolutions)
        for resolution in sftp_resolutions:
            try:
                sftp_species = self.species if self.species else set(self.sftp.listdir(join(remote_dir,resolution))).intersection(self.available_species)
            except FileNotFoundError:
                msg = f"There is no data available in {self.remote_machine} for {model_new} model at {resolution} resolution."
                show_message(self, msg, deactivate=initial_check)
                continue
            for species in sftp_species: 
                try:
                    sftp_network = self.network if self.network else self.sftp.listdir(join(remote_dir,resolution,species))
                except FileNotFoundError:
                    msg = f"There is no data available in {self.remote_machine} for {model_new} model for {species} species at {resolution} resolution."
                    show_message(self, msg, deactivate=initial_check)
                    continue
                for network in sftp_network:
                    # if network is nonghost, change the slashes to dashes
                    if not check_for_ghost(network):
                        network = network.replace("/", "-")
                    res_spec_dir.append(join(remote_dir,resolution,species,network))
        
        # print the species, resolution and model combinations that are going to be downloaded
        if res_spec_dir:

            # initialise list with all the nc files to be downloaded
            initial_check_nc_files = []

            if not initial_check:
                self.logger.info(f"\n{model_new} model data to download ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range
            for remote_dir in res_spec_dir:
                
                # get network, species and resolution
                network = remote_dir.split('/')[-1]
                species = remote_dir.split('/')[-2]
                resolution = remote_dir.split('/')[-3]

                # get local directory 
                local_dir = join(self.mod_root, self.ghost_version, model_new, resolution, species, network)

                # print source and destination  
                if not initial_check:
                    self.logger.info(f"\n  - {local_dir}, source: {remote_dir} ({self.remote_machine})")
            
                # get nc files if directory is found
                try:
                    nc_files = self.sftp.listdir(remote_dir)
                except FileNotFoundError:
                    msg = f"There is no data available in {self.remote_machine} for {model_new} model for {species} species {network} network at {resolution} resolution."
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # get the nc files in the date range       
                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if model + species + resolution + network + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available in {self.remote_machine} from {self.start_date} to {self.end_date} for {model_new} model {species} species {network} network at {resolution} resolution."
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
        
    def find_model(self, mod_id, domain, initial_check):
        # initialise warning message and model exists boolean
        msg = ""
        model_exists = False
        model_paths = []
        remote_dir = None

        # check model is in any of the interp_models.yaml lists
        for model_type, model_dict in interp_models.items():
            if mod_id in model_dict["models"]:
                model_paths = model_dict["paths"]
                break
        
        # if in the list, check if the paths work
        if model_paths:
            
            # stores all working paths
            mod_dir_functional_list = []    

            for mod_dir in model_paths:
                # esarchive in transfer5 is located inside gpfs
                if "/esarchive/" == mod_dir[:11] and self.remote_hostname.startswith('transfer'):
                    mod_dir = join("/gpfs/archive/bsc32/", mod_dir[1:])
                
                # add directory if it exists in the remote machine
                try:
                    self.sftp.stat(mod_dir)
                    mod_dir_functional_list.append(mod_dir)      
                except FileNotFoundError:
                    pass

            # if none of the paths are in this current machine, break
            if not mod_dir_functional_list:
                msg += f"None of the paths specified in {join('settings', 'interp_models.yaml')} are available on the remote machine ({self.remote_machine}). "
            
            else:
                # get first functional directory that has the model 
                for mod_dir in mod_dir_functional_list:                   
                    try:
                        remote_dir = join(mod_dir, mod_id, domain)
                        self.sftp.stat(remote_dir)
                        model_exists = True
                        break
                    except FileNotFoundError:
                        pass

                # if the model-domain combination is not possible, show the warning
                if model_exists is False:
                    msg += f"There is no data available for the {mod_id} model with the {domain} domain in none of the paths specified in {join('settings', 'interp_models.yaml')} in the remote machine ({self.remote_machine}). "

        # if no valid path model, search in the data_paths.yaml directory
        if model_exists is False:
            # get all possible models
            try:
                remote_dir = join(self.mod_to_interp_remote_path, mod_id, domain)
                self.sftp.stat(remote_dir)
                model_exists = True
            except FileNotFoundError:
                pass 
            
        # if the model-domain combination is not possible, break
        if model_exists is False:
            # add to the message if model was not found in the gpfs remote directory
            msg += f"Cannot find the {mod_id} model with the {domain} domain in '{self.mod_to_interp_remote_path}'."    
            show_message(self, msg, deactivate=initial_check)
        
        return model_exists, remote_dir
    
    def find_available_resolutions(self, remote_dir, mod_id, domain, initial_check):
        resolution_not_found = False
        valid_resolutions = []

        # iterate through the resolutions
        for resolution in self.model_resolution:
            # check existence of resolution directory
            try:
                self.sftp.stat(join(remote_dir, resolution))
            except FileNotFoundError:
                # show resolution not found
                msg = (f"There is no data available in {self.remote_machine} for the {mod_id} model with the "
                f"{domain} domain at {resolution} resolution.")
                show_message(self, msg, deactivate=initial_check)

                if not resolution_not_found:
                    # obtain the possible resolutions
                    possible_resolutions = self.sftp.listdir(remote_dir)

                    # get the order of the possible resolutions from temporal_resolution_map
                    resolution_order = temporal_resolution_map[resolution]

                    # sort the possible resolutions in the temporal_resolution_map order
                    sorted_resolutions = []
                    for res in resolution_order:
                        if res in possible_resolutions:
                            sorted_resolutions.append(res)

                    # get resolution input from user
                    if sorted_resolutions:
                        while True:
                            user_resolution = input(f"\nNo non-interpolated model data found for {resolution} resolution. "
                            f"You can specify any of the temporal resolutions found for the model: {', '.join(sorted_resolutions)}. "
                            "If left empty, the first available option will be selected automatically: ").lower()

                            if user_resolution in possible_resolutions or user_resolution == '':
                                if user_resolution != '':
                                    possible_resolutions = [user_resolution]
                                break
                        
                        # add resolution to the return list
                        valid_resolutions.append(possible_resolutions[0])
                        
                        # block input from the user happening again
                        resolution_not_found = True
                    else: 
                        msg = (f"There is no data available in {self.remote_machine} for the {mod_id} model " 
                                f"with the {domain} domain at {resolution} resolution.")
                        show_message(self, msg, deactivate=initial_check)
                        continue
            else:
                # add resolution to the return list
                valid_resolutions.append(resolution)
        
        return valid_resolutions
    
    def find_available_species(self, valid_resolutions, remote_dir, ensemble, mod_id, domain, initial_check):
        valid_resolutions_species_dir = []

        for resolution in valid_resolutions: 
            species_exists = False
                
            # iterate through the species
            for original_species in self.species: 
                species_to_process = [original_species]

                # iterate through the original species and if not found through the mapped species
                if original_species in mapping_species:
                    species_to_process = species_to_process + mapping_species[original_species]

                for species in species_to_process:
                    # check existance of species directory
                    try:
                        # ensemble member
                        if ensemble.isdigit() or ensemble == 'allmembers':
                            resolution_species_dir = join(remote_dir, resolution, species)
                        # ensemble statistic
                        else:
                            resolution_species_dir = join(remote_dir, resolution, "ensemble-stats", species + "_" + ensemble)

                        self.sftp.stat(resolution_species_dir)

                        # store the resolution species pair
                        valid_resolutions_species_dir.append([resolution_species_dir, mod_id, domain, resolution, species, ensemble])
                        
                        # do not show warning
                        species_exists = True

                        # get the first species found
                        break
                    
                    # check existance of mapped species directory
                    except FileNotFoundError:
                        pass
                    
            # if no species were found, then show the message
            if species_exists is False:
                msg = (f"There is no data available in {self.remote_machine} for the {mod_id} model with the "
                       f"{domain} domain for {species} species at {resolution} resolution.")
                show_message(self, msg, deactivate=initial_check)
                continue
        
        return valid_resolutions_species_dir
    
    def find_available_nc_files(self, valid_resolutions_species_dir, initial_check):
        # initialise list with all the nc files to be downloaded
        path_files_dict = {}

        # get all the nc files in the date range
        for resolution_species_dir, mod_id, domain, resolution, species, ensemble in valid_resolutions_species_dir:
                        
            # get nc files
            nc_files = self.sftp.listdir(resolution_species_dir)

            if nc_files:
                # ensemble member
                if ensemble.isdigit() or ensemble == 'allmembers':

                    # identify format of the directory
                    # the format is a tuple of how many '-' and how many '_' are there, e.g.: (0,1)
                    # the directory format is choosen by popularity
                    formats_list = [(file.count("-"), file.count("_")) for file in nc_files]
                    number_of_formats_dict = {format: formats_list.count(format) for format in set(formats_list)}
                    format = max(number_of_formats_dict, key=number_of_formats_dict.get)
                    
                    # filter and get only the files that follow the format choosen
                    nc_files = list(filter(lambda x:(x.count("-"),x.count("_")) == format and x.endswith(".nc"), nc_files))
                    
                    # example: od550du_2019040212.nc (0,1)
                    if format == (0, 1):
                        # if no ensemble in the name only allmembers and 000 are valid
                        if ensemble == '000' or ensemble == 'allmembers':
                            nc_files = list(filter(lambda x:x.split("_")[0] == species, nc_files))

                    # example: od550du-000_2021020812.nc (1,1)
                    elif format == (1, 1):
                        # filter by ensemble in case that ensemble is not allmembers
                        if ensemble != 'allmembers':
                            nc_files = list(filter(lambda x:x.split("_")[0] == species + '-' + ensemble, nc_files))

                    # unknown format
                    else:
                        msg = f"It is not possible to download this nc file type yet. Please, contact the developers. Files to download: {nc_files}"
                        show_message(self, msg, deactivate=initial_check)
                        continue

                    local_dir = join(self.mod_to_interp_root, mod_id, domain, resolution, species)

                # ensemble statistic
                else:
                    # example: sconco3_2025120300_av_an.nc
                    nc_files = list(filter(lambda x:x.split("_")[0] == species and "_".join(x[:-3].split("_")[2:]) == ensemble, nc_files))

                    local_dir = join(self.mod_to_interp_root, mod_id, domain, resolution, "ensemble-stats", species + "_" + ensemble)

                # store the available files for each path
                if nc_files:
                    path_files_dict[local_dir] = {"nc_files" : nc_files,
                                                  "mod_id" : mod_id, 
                                                  "species" : species, 
                                                  "resolution" : resolution}
                
                else:
                    msg = f"There is no data available in {self.remote_machine} for the {mod_id} model with the {domain} domain with the {ensemble} ensemble."
                    show_message(self, msg, deactivate=initial_check)

        return path_files_dict
    
    def build_nc_file_paths_in_range(self, path_files_dict, initial_check):
        initial_check_nc_files = []

        for path, model_dict in path_files_dict.items():
            nc_files = model_dict["nc_files"]

            valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files) 
            valid_nc_files.sort()   

            # warning if model + species + resolution + network + date range combination gets no matching results       
            if not valid_nc_files:                 
                msg = (f"There is no data available in {self.remote_machine} from {self.start_date} to {self.end_date} "
                       f"for {model_dict["mod_id"]} model {model_dict["species"]} species at {model_dict["resolution"]} resolution.")
                show_message(self, msg, deactivate=initial_check)
                continue 

            for nc_file in valid_nc_files:
                initial_check_nc_files.append(join(path, nc_file))

        return initial_check_nc_files

    def download_non_interpolated_model(self, model, initial_check, files_to_download=None):
        """
        Download non-interpolated model data from a remote machine via SFTP.

        Parameters
        ----------
        model : str
            Name of the model to download.
        initial_check : bool
            If True, only performs a check of available files without downloading.
            If False, downloads files and displays progress.
        files_to_download : list of str, optional
            Specific file paths to download. Only files in this list are considered.

        Returns
        -------
        initial_check_nc_files : list of str
            A list of file paths intended for download.
        """

        # check if ssh exists and check if still active, connect if not
        if (self.ssh is None) or (self.ssh.get_transport().is_active()):
            self.connect()  

        # get model id and the domain
        mod_id, domain, ensemble = model.split("-")
        
        if initial_check:
            # print current model
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nDownloading {model} non-interpolated model data from {self.remote_machine}...")

            # look for the model in the remote machine
            model_exists, remote_dir = self.find_model(mod_id, domain, initial_check)

            if not model_exists:
                return
            
            valid_resolutions = self.find_available_resolutions(remote_dir, mod_id, domain, initial_check)
            
            valid_resolutions_species_dir = self.find_available_species(valid_resolutions, remote_dir, ensemble, mod_id, domain, initial_check)

            path_files_dict = self.find_available_nc_files(valid_resolutions_species_dir, initial_check)

            initial_check_nc_files = self.build_nc_file_paths_in_range(path_files_dict, initial_check)

            return initial_check_nc_files
                        
        else:
            self.logger.info(f"\n{model} model data to download ({len(res_spec_dir)}):")

            # print source and destination
            self.logger.info(f"\n  - {local_dir}, source: {remote_dir} ({self.remote_machine})")
 
            # create directories if they don't exist
            if not os.path.exists(local_dir):
                os.makedirs(local_dir) 

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

            # tell the user if not valid resolution specie date combinations
            msg = "There is no available model output to be downloaded."
            show_message(self, msg, deactivate=initial_check)

    # TODO         
    def copy_non_interpolated_model(self, model, initial_check, files_to_download=None):
        """
        Copy from esarchive to gpfs interpolated model data from a remote machine via SFTP.

        Parameters
        ----------
        model : str
            Name of the model to download.
        initial_check : bool
            If True, only performs a check of available files without downloading.
            If False, downloads files and displays progress.
        files_to_download : list of str, optional
            Specific file paths to download. Only files in this list are considered.

        Returns
        -------
        initial_check_nc_files : list of str
            A list of file paths intended for download.
        """

        if not initial_check:
            # print current model
            self.logger.info('\n'+'-'*40)
            self.logger.info(f"\nCopying {model} non-interpolated model data from esarchive to gpfs in {self.machine}...")
            
        # get resolution and species combinations
        res_spec_dir = []

        # get model id and the domain
        mod_id, domain, ensemble = model.split("-")

        # get model type
        for model_type, model_dict in interp_models.items():
            if mod_id in model_dict["models"]:
                break
        
        # get model specific directories list
        mod_dir_list = model_dict["paths"]

        # take all functional directories
        mod_dir_functional_list = []    
        for mod_dir in mod_dir_list:
            # make sure that it comes from esarchive
            if "/esarchive/" in mod_dir:
                # esarchive in transfer5 is located inside gpfs
                if "/esarchive/" == mod_dir[:11] and self.machine == "storage5":
                    mod_dir = join("/gpfs/archive/bsc32/",mod_dir[1:])
                # check if directory exists in esarchive
                if os.path.exists(mod_dir):
                    mod_dir_functional_list.append(mod_dir)     
            
        # if none of the paths are in this current machine, break
        if not mod_dir_functional_list:
            msg = f"None of the paths specified in {join('settings', 'interp_models.yaml')} are available on esarchive."
            show_message(self, msg, deactivate=initial_check)
            return
        
        # take first functional directory  
        esarchive_dir = None
        for mod_dir in mod_dir_functional_list:
            temp_esarchive_dir = join(mod_dir,mod_id,domain)
            # check if model and domain directories exist in esarchive machine
            if os.path.exists(temp_esarchive_dir): 
                esarchive_dir = temp_esarchive_dir
                break
        
        # if the model-domain combination is not possible, break
        if esarchive_dir is None:
            msg = f"There is no data available for the {mod_id} model with the {domain} domain in none of the paths specified in {join('settings', 'interp_models.yaml')} in esarchive."
            show_message(self, msg, deactivate=initial_check)
            return

        # get all the resolutions available in the esarchive directory
        sftp_resolutions = self.model_resolution if self.model_resolution else set(os.listdir(esarchive_dir)).intersection(self.nonghost_available_resolutions)

        # iterate through the resolutions
        for resolution in sftp_resolutions:
            try:
                # get available species ("normal" and mapped)
                available_species = self.available_species+[spec[0] for spec in mapping_species.values()]
                sftp_species = self.species if self.species else set(os.listdir(join(esarchive_dir,resolution))).intersection(available_species)
            except FileNotFoundError:
                msg = f"There is no data available in esarchive for the {mod_id} model with the {domain} domain at {resolution} resolution."
                show_message(self, msg, deactivate=initial_check)
                continue

            # iterate through the species
            for speci_to_process in sftp_species: 
                # initialize boolean that saves whether species was found
                species_exists = False
                species = speci_to_process

                # if it is an ensemble member
                if ensemble.isdigit() or ensemble == 'allmembers':
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
                            if ensemble.isdigit() or ensemble == 'allmembers':
                                res_spec = join(esarchive_dir,resolution,species)
                            # if it is an ensemble statistic
                            else:
                                res_spec = join(esarchive_dir,resolution,"ensemble-stats",species+"_"+ensemble)
  
                            species_exists = os.path.exists(res_spec)
                
                # if no species were found, then show the message
                if species_exists is False:
                    msg = f"There is no data available in esarchive for the {mod_id} model with the {domain} domain for {species} species at {resolution} resolution."
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # add the path with the resolution and species combination to the list
                res_spec_dir.append(res_spec)
                        
        # print the species, resolution and model combinations that are going to be copied
        if res_spec_dir:

            # initialise list with all the nc files to be copied
            initial_check_nc_files = []

            if not initial_check:
                self.logger.info(f"\n{model} model data to copy ({len(res_spec_dir)}):")
            
            # get all the nc files in the date range
            for esarchive_dir in res_spec_dir:
                if not initial_check:
                    self.logger.info(f"\n  - {join(self.mod_to_interp_root,'/'.join(esarchive_dir.split('/')[-4:]))}, source: {esarchive_dir} ({self.machine})")
                         
                # get nc files
                nc_files = os.listdir(esarchive_dir)

                if nc_files:
                    # if it is an ensemble member
                    if ensemble.isdigit() or ensemble == 'allmembers':
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
                    msg = f"There is no data available in esarchive for the {mod_id} model with the {domain} domain with the {ensemble} ensemble."
                    show_message(self, msg, deactivate=initial_check)
                    continue
                
                # get the nc files in the date range        
                valid_nc_files = self.get_valid_nc_files_in_date_range(nc_files)

                # warning if model + species + resolution + network + date range combination gets no matching results       
                if not valid_nc_files:                 
                    msg = f"There is no data available in esarchive from {self.start_date} to {self.end_date} for {model} model {species} species at {resolution} resolution."
                    show_message(self, msg, deactivate=initial_check)
                    continue

                # copy the valid resolution specie date combinations
                else:
                    # if it is an ensemble member
                    if ensemble.isdigit() or ensemble == 'allmembers':
                        gpfs_dir = join(self.mod_to_interp_root,mod_id,domain,resolution,species)
                    else:
                        gpfs_dir = join(self.mod_to_interp_root,mod_id,domain,resolution,"ensemble-stats",species+"_"+ensemble)
                    
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
            msg = "There is no available model output to be copied."
            show_message(self, msg, deactivate=initial_check)

    def get_all_networks(self): 
        """Populate the `self.network` attribute with all the available networks."""

        if self.reading_ghost:
            if self.dl_ghost_source == 'bsc':
                self.network = self.ghost_available_networks
            elif not hasattr(self,"zenodo_ghost_available_networks"): 
                if not self.zenodo.fetch_zenodo_networks():
                    return
                self.network = list(self.zenodo_ghost_available_networks.keys())
        else:
            self.network = self.nonghost_available_networks
    
    def get_all_models(self):
        """Populate the `self.experiments` attribute with all the available models."""
        
        # download all interpolated models
        if self.dl_interpolated is True:
            # check if ssh exists and check if still active, connect if not
            if (self.ssh is None) or (self.ssh.get_transport().is_active()):
                self.connect()  

            # get directory content and format it as the models       
            model_list = self.sftp.listdir(join(self.mod_remote_path,self.ghost_version))
        # download all non interpolated models
        else:
            # get all the models id
            models = []
            for model_dict in interp_models.values():
                models += model_dict["models"]
            # get all the domain and ensemble combinations 
            model_list = []
            # TODO hardcoded
            for domain in ["ip", "d03", "d01", "regional", "eu", "reg", "ex", "bcn", "cat", "d02", "global","regional_i01", "regional_i02", "regional_i03"]:
                for mod in models:
                    model_list.append(mod+"-"+domain+"-allmembers")

        self.experiments = dict(zip(model_list,model_list))

    def get_valid_nc_files_in_date_range(self, nc_files):
        """
        Filter NetCDF files by the object's date range.

        Parameters
        ----------
        nc_files : list of str
            List of NetCDF filenames to filter.

        Returns
        -------
        valid_nc_files : list of str
            Subset of `nc_files` whose dates fall within the date range.
        """

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
        """
        Monitor download time and abort if a timeout is exceeded.

        Parameters
        ----------
        size : int
            Number of bytes transferred so far.
        file_size : int
            Total size of the file being downloaded.
        """
         
        if (time.time() - self.ncfile_dl_start_time) > self.dl_timeout:
            error = 'Download timeout, try later.'
            self.logger.error(error)
            sys.exit(1)
            
    def sighandler(self, *unused):
        """
        Handle keyboard interrupts and perform cleanup before exiting.

        Parameters
        ----------
        *unused : tuple
            Positional arguments required by the signal handler interface.
            These arguments are not used.
        """

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

        # delete Zenodo and CAMS temp dirs if necessary
        for root in ['mod_to_interp_root', 'ghost_root']:
            temp_dir = join(getattr(self, root),'.temp')
            if os.path.exists(temp_dir):
                self.logger.info(f"\nDeleting {temp_dir}")
                shutil.rmtree(temp_dir)
        
        self.logger.info("\nExiting...")
        sys.exit()

def main(**kwargs):
    """
    Initialize the download environment and launch the download workflow.
    
    Parameters
    ----------
    **kwargs
        Optional command-line arguments that override default configuration values.
    """

    download = Download(**kwargs)
    download.run()