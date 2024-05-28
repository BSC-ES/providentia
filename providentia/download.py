import sys
import os

import requests
from io import BytesIO
from pexpect import pxssh
import subprocess
import yaml
from dotenv import load_dotenv

# urlparse
from tqdm import tqdm
from remotezip import RemoteZip
import tarfile
from datetime import datetime, timedelta

from .configuration import ProvConfiguration, load_conf
from .read_aux import check_for_ghost

# MACHINE = os.environ.get('BSC_MACHINE', '')
MACHINE = 'local' #it will always be local

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

data_paths = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/data_paths.yaml')))

class ProvidentiaDownload(object):
    def __init__(self,**kwargs):
        # initialise zenodo url
        self.ghost_url = 'https://zenodo.org/records/10637450'

        # initialize dictionary to store zenodo's ghost zips
        self.ghost_zip_files = {}

        # TODO CHANGE TO DEFAULTS
        # This is an atribute which is used in case there's no networks, if no networks then all the networks from this
        self.download_source = "ghost"

        # initialise remote hostname
        self.remote_hostname = "transfer1.bsc.es"

        # get home user
        self.home_user = os.getlogin()
        
        # flag to indicate if user wants their user and password saved 
        remind = False

        # TODO: Put if exp or nonghost network only
        # get user and password if they are in .env
        load_dotenv()
        user = os.getenv("PRV_USER")
        password = os.getenv("PRV_PWD")

        # if couldn't get variables from .env, ask user for them
        if not user or not password:
            user = input("Introduce BSC user: ")
            password = input("Introduce password: ")
            remind_txt = input("Remember password (y/n)? ")
            remind = remind_txt.lower() == 'y'

        # initialise default configuration variables
        # modified by commandline arguments, if given
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

        # initialise data_paths for local
        # TODO CHANGE THIS IN CONFIGURATION

        self.nonghost_root = os.path.join('/home',self.home_user,data_paths[MACHINE]['nonghost_root'])
        self.ghost_root = os.path.join('/home',self.home_user,data_paths[MACHINE]['ghost_root']) 
        self.exp_root = os.path.join('/home',self.home_user,data_paths[MACHINE]['exp_root']) 

        # create empty directories for all the paths if they don't exist
        for path in [self.nonghost_root,self.ghost_root,self.exp_root]:
            if not os.path.exists(path):
                os.makedirs(path)

        for section_ind, section in enumerate(self.parent_section_names):
            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # update self with section variables
            for k, val in self.section_opts.items():
                setattr(self, k, provconf.parse_parameter(k, val))

            if not self.network:
                self.get_all_networks()

            # network
            for network in self.network:
                # ghost
                if check_for_ghost(network):
                    self.download_ghost_network(network)
                # non-ghost
                else:
                    self.download_nonghost_network(network)

            # print(self.nonghost_root,self.network,self.resolution,self.species)

            # experiment
            for experiment in self.experiments:
                pass

            # print(self.nonghost_root,self.network,self.resolution,self.species)
        
        # print(self.exp_root)

    def download_nonghost_network(self,network):
        # network
        if not self.resolution and not self.species:
            wanted_directories = ['']
        # network + resolution 
        elif not self.species:
            wanted_directories = [self.resolution]
        # network + species 
        elif not self.resolution:
            wanted_directories = self.species
        # network + resolution + species 
        else:
            wanted_directories = [os.path.join(self.resolution,specie) for specie in self.species]
        
        # print(wanted_directories)

        # list_subfolders_with_paths = [f.path for f in os.scandir(self.nonghost_root) if f.is_dir()]
        # print(list_subfolders_with_paths)
        for wanted_directory in wanted_directories:
            # print("\npaula",os.path.join(self.nonghost_root,network,wanted_directory))
            directory = os.path.join(self.nonghost_root,network)
            # print(wanted_directories,"\n")
            for (root, directories, files) in os.walk(directory):
                if wanted_directory in root:
                    print("root",root)#,"\n")

       # dirs_to_get = [directory for directory in  if directory[-7:] == '.tar.xz' and wanted_directories in directory]

    def download_ghost_network(self,network):
        # if first time reading a ghost network, get current zips urls in zenodo page
        if not self.ghost_zip_files: 
            self.get_ghost_zip_files()
        
        # get resolution and/or species combinations
        # network
        if not self.resolution and not self.species:
            res_spec_combinations = [""]
        # network + resolution 
        elif not self.species:
            res_spec_combinations = [f"/{resolution}/" for resolution in self.resolution]
        # network + species 
        elif not self.resolution:
            res_spec_combinations = self.species    
        # network + resolution + species 
        else:
            res_spec_combinations = [f"/{resolution}/{specie}" for specie in self.species for resolution in self.resolution]

        # get all the possible networks from network e.g. EBAS gets EBAS-ACTRIS_oyktp, EBAS-AMAP_hcxwm, EBAS-CAMP_daczg, etc
        current_networks = list(filter(lambda x: network in x, self.ghost_zip_files))

        for curr_net in current_networks:
            # get url to download the zip file for the current network
            zip_url = self.ghost_zip_files[curr_net]

            #TODO: Add comprovation in case is already downloaded

            # get all the species tar files which fulfill the combination condition
            species_to_download = []
            with RemoteZip(zip_url) as zip:
                for combi in res_spec_combinations:
                    species_to_download += list(filter(lambda x: combi in x, zip.namelist()))
                species_to_download = list(filter(lambda x: x[-7:] == '.tar.xz', species_to_download))

                # Print the specie, resolution and network combinations that are going to be downloaded
                print("Observations to download:")
                for specie in species_to_download:
                    print(f"  - {specie}")

                # extract species from zip files
                for specie_to_get in tqdm(species_to_download,unit="specie",desc="Downloading species"):
                    zip.extract(specie_to_get,self.ghost_root)
                    
                    # get path and the name of the directory of the tar file
                    tar_path = os.path.join(self.ghost_root, specie_to_get)
                    tar_dir = os.path.dirname(tar_path)

                    # extract nc file from tar file
                    with tarfile.open(tar_path) as tar_file:
                        tar_names = tar_file.getnames()

                        # TODO Maybe valid dates is gonna change to the month generator, see what is said on providentia meeting
                        # get the nc files that are between the start and end date
                        valid_nc_file_names = self.get_valid_dates(tar_names)
                        valid_nc_files = [tar_member for tar_member in tar_file.getmembers() if tar_member.name in valid_nc_file_names]

                        tar_file.extractall(path = tar_dir, members = valid_nc_files)
                    
                    # remove the tar file
                    os.remove(tar_path)

    def download_exp(self):
        pass    

    def get_ghost_zip_files(self):
        # Get urls from zenodo to get ghost zip files url
        # Send an HTTP GET request to the URL
        # TODO CHANGE NAME
        response = requests.get(self.ghost_url)

        # Check if the request was successful (status code 200)
        if response.status_code != 200:
            msg = f'Failed to retrieve the webpage. Status code: {response.status_code}'
        
        # fill network dictionary with its corresponding zip url
        for line in response.text.split(">"):
            if '<link rel="alternate" type="application/zip" href=' in line:
                zip_url = line.split('href="')[-1][:-1]
                zip_network = line.split("/")[-1][:-5]
                self.ghost_zip_files[zip_network] = zip_url

    def get_all_networks(self):
        # if not network passed, then get all networks

        if self.download_source.lower() == "ghost":
            if not self.ghost_zip_files: 
                self.get_ghost_zip_files()
            
            self.networks = list(self.ghost_zip_files.keys())

        elif self.download_source.lower() == "nonghost":
            pass

        else:
            print(f"Error download source type {self.download_source} not valid")

    def generate_months(start_date, end_date):
        # TODO en proceso maybe quitar
        start = datetime.strptime(start_date, "%Y%m%d")
        end = datetime.strptime(end_date, "%Y%m%d")
        
        month_list = []
        
        current_month = start
        while current_month < end:
            formatted_month = current_month.strftime("%Y%m")
            current_month += timedelta(days=31)
            month_list.append(current_month)
        
        formatted_end_month = end.strftime("%Y%m") # si no quiere que incluya el ultimo borrar esta linea
        month_list.append(formatted_end_month)
        return month_list
    
    def get_valid_dates(self, tar_names):
        nc_files = []

        for nc_file in tar_names:
            if ".nc" in nc_file:
                ym = nc_file.split("_")[-1].split(".nc")[0]
                if int('{}01'.format(ym)) >= int(self.start_date) and int('{}01'.format(ym)) < int(self.end_date):
                    nc_files.append(nc_file)
                    
        return nc_files        

def main(**kwargs):
    """ Main function when running download function. """
   
    ProvidentiaDownload(**kwargs)