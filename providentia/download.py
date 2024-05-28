import sys
import os

import requests
from io import BytesIO
import subprocess
import yaml
from dotenv import dotenv_values
import paramiko 
from base64 import decodebytes

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
        # TODO move some of these variables to default
        # initialise zenodo url
        self.ghost_url = 'https://zenodo.org/records/10637450'

        # initialize dictionaries to store possible networks
        self.ghost_zip_files = {}
        self.nonghost_observation_data = {}

        # initialise remote hostname and esarchive path in storage5
        # TODO MAYBE MOVE TO DATAPATHS
        self.remote_hostname = "transfer1.bsc.es"
        self.nonghost_remote_obs_path = "/gpfs/archive/bsc32/esarchive/obs"

        # get home user
        self.home_user = os.getlogin()

        # get ssh user and password if they are in .env
        env = dotenv_values(os.path.join(PROVIDENTIA_ROOT, ".env"))
  
        self.prv_user = env.get("PRV_USER")
        self.prv_password = env.get("PRV_PWD")

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

        # create empty directories for all the paths if they don't exist
        for path in [self.nonghost_root,self.ghost_root,self.exp_root]:
            if not os.path.exists(path):
                try:
                    os.makedirs(path)
                except PermissionError as error:
                    os.system(f"sudo mkdir -p {path}")
                    os.system(f"chmod o+w {path}")

        for section_ind, section in enumerate(self.parent_section_names):
            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # update self with section variables
            for k, val in self.section_opts.items():
                setattr(self, k, provconf.parse_parameter(k, val))

            # if networks is none, then get all possible networks 
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

            # experiment
            for experiment in self.experiments:
                pass

    def download_nonghost_network(self,network):
        # if first time reading a nonghost network, get current networks from yaml
        if not self.nonghost_observation_data: 
            self.nonghost_observation_data = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/nonghost_networks.yaml')))

        # If not valid network, next
        if network not in self.nonghost_observation_data:
            return

        # flag to indicate if user wants their user and password saved 
        remind = False

        # if couldn't get credentials, ask user for them
        if not self.prv_user  or not self.prv_password:
            self.prv_user  = input("Introduce BSC user: ")
            self.prv_password  = input("Introduce password: ")
            remind_txt = input("Remember password (y/n)? ")
            remind = remind_txt.lower() == 'y'

        # get public remote machine public key and add it to ssh object
        _, output = subprocess.getstatusoutput("ssh-keyscan -t ed25519 transfer1.bsc.es")
        ed25519_key = output.split()[-1].encode()
        key = paramiko.Ed25519Key(data=decodebytes(ed25519_key))

        ssh = paramiko.SSHClient()
        hostkeys = ssh.get_host_keys().add(self.remote_hostname, 'ed25519', key)
        
        # connect to machine
        ssh.connect(self.remote_hostname, username=self.prv_user, password=self.prv_password)

        # create .env with the input user and password
        if remind:
            with open(os.path.join(PROVIDENTIA_ROOT, ".env"),"w") as f:
                f.write(f"PRV_USER={self.prv_user}\n")
                f.write(f"PRV_PWD={self.prv_password}")

        # get resolution and/or species combinations
        res_spec_combinations = []
        # network
        if not self.resolution and not self.species:
            for res, specie_list in self.nonghost_observation_data[network].items():
                for spec in specie_list:
                    res_spec_combinations.append(os.path.join(network,res,spec))
            
        # network + resolution 
        elif not self.species:
            for res in self.resolution:
                for spec in self.nonghost_observation_data[network].get(res,[]):
                    res_spec_combinations.append(os.path.join(network,res,spec))
        # network + species 
        elif not self.resolution:
            for spec in self.species:
                for res, specie_list in self.nonghost_observation_data[network].items():
                    if spec in specie_list:
                        res_spec_combinations.append(os.path.join(network,res,spec))
        # network + resolution + species 
        else:
            for spec in self.species:
                for res in self.resolution:
                    if spec in self.nonghost_observation_data[network].get(res,[]):
                        res_spec_combinations.append(os.path.join(network,res,spec))
        
        # create Secure File Transfer Protocol object
        sftp = ssh.open_sftp()

        # print the specie, resolution and network combinations that are going to be downloaded
        if res_spec_combinations:
            print(f"{network} observations to download:")
            for combi_print in res_spec_combinations:
                print(f"  - {self.nonghost_root}/{combi_print}")

            # get all the nc files in the date range
            for combi in tqdm(res_spec_combinations,desc="Downloading Observations"):
                remote_dir = os.path.join(self.nonghost_remote_obs_path,combi)
                local_dir = os.path.join(self.nonghost_root,combi)

                valid_nc_files = self.get_valid_dates(sftp.listdir(remote_dir))

                if valid_nc_files:

                    # create directories if they don't exist
                    if not os.path.exists(local_dir):
                        os.makedirs(local_dir) 

                    # download each individual nc file using sftp protocol
                    for nc_file in valid_nc_files:
                        remote_path = os.path.join(remote_dir,nc_file)
                        local_path = os.path.join(local_dir,nc_file)
                        sftp.get(remote_path,local_path)

        # close conections
        sftp.close()
        ssh.close()

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

            # get all the species tar files which fulfill the combination condition
            res_spec_final_combinations = []
            with RemoteZip(zip_url) as zip:
                for combi in res_spec_combinations:
                    res_spec_final_combinations += list(filter(lambda x: combi in x, zip.namelist()))
                res_spec_final_combinations = list(filter(lambda x: x[-7:] == '.tar.xz', res_spec_final_combinations))

                # TODO Deberia printear cuanto va a ocupar pero no se sabe porque para cada uno no se va a descargar todas las fechas
                # SE PUEDE PONER DESPUES

                # Print the specie, resolution and network combinations that are going to be downloaded
                print(f"{network} observations to download:")
                for specie in res_spec_final_combinations:
                    print(f"  - {self.nonghost_root}{specie[:-7]}")

                # extract species from zip files
                for specie_to_get in tqdm(res_spec_final_combinations,desc="Downloading Observations:"):
                    zip.extract(specie_to_get,self.ghost_root)
                    
                    # get path and the name of the directory of the tar file
                    tar_path = os.path.join(self.ghost_root, specie_to_get)
                    tar_dir = os.path.dirname(tar_path)

                    # extract nc file from tar file
                    with tarfile.open(tar_path) as tar_file:
                        tar_names = tar_file.getnames()

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
        # TODO CHANGE METHOD'S NAME
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
        # TODO maybe change this prints to warnings
        # if not network passed, then get all networks
        print(f"Not networks found, downloading all {self.download_source.lower()} networks.")

        if self.download_source.lower() == "ghost":
            if not self.ghost_zip_files: 
                self.get_ghost_zip_files()
            
            self.network = list(self.ghost_zip_files.keys())

        elif self.download_source.lower() == "nonghost":
            self.nonghost_observation_data = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/nonghost_networks.yaml')))
            self.network = list(self.nonghost_observation_data.keys())

        else:
            #TODO PONER ERROR DE VERDAD
            print(f"Error download source type {self.download_source} not valid")

    def generate_months(start_date, end_date):
        # TODO en proceso maybe quitar
        # DEPRECATED DEPRECATED DEPRECATED
        start = datetime.strptime(start_date, "%Y%m%d")
        end = datetime.strptime(end_date, "%Y%m%d")
        
        month_list = []
        
        current_month = start
        while current_month < end:
            formatted_month = current_month.strftime("%Y%m")
            current_month += timedelta(days=31)
            month_list.append(current_month)
        
        formatted_end_month = end.strftime("%Y%m") 
        month_list.append(formatted_end_month)
        return month_list
    
    def get_valid_dates(self, tar_names):
        # get valid nc files inside the date range
        # TODO CAMBIAR NOMBRE Y COMENTAR 
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