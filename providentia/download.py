import sys
import os

import requests
from io import BytesIO
import subprocess
import yaml
import json
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

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

data_paths = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/data_paths.yaml')))

class ProvidentiaDownload(object):
    def __init__(self,**kwargs):
        # TODO move some of these variables to default
        # initialise zenodo url
        self.ghost_url = 'https://zenodo.org/records/10637450'

        # initialize dictionaries to store possible networks
        self.ghost_zip_files = {} # TODO: Maybe juntar este con el del zenodo
        self.nonghost_observation_data = {}
        self.ghost_observation_data = {}

        # initialise remote hostname and path in storage5
        # TODO MAYBE MOVE TO DATAPATHS
        self.remote_hostname = "transfer1.bsc.es"
        REMOTE_MACHINE = "storage5"

        self.ghost_remote_obs_path = data_paths[REMOTE_MACHINE]["ghost_root"]
        self.nonghost_remote_obs_path = data_paths[REMOTE_MACHINE]["nonghost_root"]
        self.exp_remote_path = data_paths[REMOTE_MACHINE]["exp_root"]

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

            # now all variables have been parsed, check validity of those, throwing errors where necessary
            provconf.check_validity()

            # if networks is none, raise error
            if not self.network:
                error = "Error: No networks were passed."
                sys.exit(error)
            
            # when one of those symbols is passed, get all networks
            if self.network == ["*"] or self.network == ["default"]:
                self.get_all_networks()
            
            for network in self.network:
                # ghost
                if check_for_ghost(network):
                    self.download_ghost_network(network)
                # non-ghost
                else:
                    self.download_nonghost_network(network)

            # filter species networks
            for i,network_specie in enumerate(self.filter_species):
                # get network and species from filter_species
                network, species = network_specie.split('|')
                self.species = [species]

                # ghost
                if check_for_ghost(network):
                    self.download_ghost_network(network)
                # non-ghost
                else:
                    self.download_nonghost_network(network)

                # TODO experiment
                for experiment in self.experiments:
                    pass

    def connect(self):
        # flag to indicate if user wants their user and password saved 
        remind = False

        # if couldn't get credentials, ask user for them
        if not self.prv_user or not self.prv_password:
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
        
        return ssh

    def download_nonghost_network(self,network):
        # print current_network
        print(f"\nDownloading non-GHOST network: {network}...")

        # if first time reading a nonghost network, get current networks from yaml
        if not self.nonghost_observation_data: 
            self.nonghost_observation_data = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/nonghost_networks.yaml')))

        # If not valid network, next
        if network not in self.nonghost_observation_data:
            print(f"There is no data available for {network} network.")
            return

        ssh = self.connect()

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
        
        # warning if network + species + resolution combination is gets no matching results
        if not res_spec_combinations:
            print_spec = f'{",".join(self.species)} species' if self.species else ""
            print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
            print(f"There is no data available for {network} network {print_spec} {print_res}")
            return

        # print the species, resolution and network combinations that are going to be downloaded
        else:
            print(f"{network} observations to download:")
            for combi_print in res_spec_combinations:
                print(f"  - {self.nonghost_root}/{combi_print}")
            
            # get all the nc files in the date range
            for combi in tqdm(res_spec_combinations,desc="Downloading Observations"):
                remote_dir = os.path.join(self.nonghost_remote_obs_path,combi)
                local_dir = os.path.join(self.nonghost_root,combi)

                valid_nc_files = self.get_valid_dates(sftp.listdir(remote_dir))

                # warning if network + species + resolution + date range combination gets no matching results       
                if not valid_nc_files:                 
                    print_spec = f'{",".join(self.species)} species' if self.species else ""
                    print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
                    print(f"There is no data available from {self.start_date} to {self.end_date} for {network} network {print_spec} {print_res}")

                else:
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
        # print current_network
        print(f"\nDownloading GHOST network: {network}...")

        # if first time reading a ghost network, get current networks from yaml
        if not self.ghost_observation_data: 
            temp_ghost_dict = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/ghost_filetree.json'))) 

            self.ghost_observation_data = {}
            for flt_network, resolution_dict in temp_ghost_dict.items():
                self.ghost_observation_data[flt_network] = {}
                for flt_resolution, FAKE_species_dict in resolution_dict.items():
                    for flt_FAKE_specie, species_dict in FAKE_species_dict.items():
                        self.ghost_observation_data[flt_network][flt_resolution] = list(species_dict.keys())
            
        # # TODO PREGUNTAR QUE SON LOS VALORES INTERMEDIOS
        # for flt_network, resolution_dict in self.ghost_observation_data.items():
        #     print(resolution_dict.keys())
        #     for flt_resolution, species_dict in resolution_dict.items():
        #         print(species_dict.keys())
        #         for nose,species2_dict in species_dict.items():
        #             print(species2_dict.keys())
        #             break
        #         break
        #     break
   
        # If not valid network, next
        if network not in self.ghost_observation_data:
            print(f"There is no data available for {network} network.")
            return

        ssh = self.connect()

        # get resolution and/or species combinations
        res_spec_combinations = []
        
        # network
        if not self.resolution and not self.species:
            for res, specie_list in self.ghost_observation_data[network].items():
                for spec in specie_list:
                    res_spec_combinations.append(os.path.join(network,self.ghost_version,res,spec))    
        # network + resolution 
        elif not self.species:
            for res in self.resolution:
                for spec in self.ghost_observation_data[network].get(res,[]):
                    res_spec_combinations.append(os.path.join(network,self.ghost_version,res,spec))
        # network + species 
        elif not self.resolution:
            for spec in self.species:
                for res, specie_list in self.ghost_observation_data[network].items():
                    if spec in specie_list:
                        res_spec_combinations.append(os.path.join(network,self.ghost_version,res,spec))
        # network + resolution + species 
        else:
            for spec in self.species:
                for res in self.resolution:
                    if spec in self.ghost_observation_data[network].get(res,[]):
                        res_spec_combinations.append(os.path.join(network,self.ghost_version,res,spec))
        
        # create Secure File Transfer Protocol object
        sftp = ssh.open_sftp()
        
        # warning if network + species + resolution combination is gets no matching results
        if not res_spec_combinations:
            print_spec = f'{",".join(self.species)} species' if self.species else ""
            print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
            print(f"There is no data available for {network} network {print_spec} {print_res}")
            return

        # print the species, resolution and network combinations that are going to be downloaded
        else:
            print(f"{network} observations to download:")
            for combi_print in res_spec_combinations:
                print(f"  - {self.ghost_root}/{combi_print}")

            # get all the nc files in the date range
            for combi in tqdm(res_spec_combinations,desc="Downloading Observations"):
                remote_dir = os.path.join(self.ghost_remote_obs_path,combi)
                local_dir = os.path.join(self.ghost_root,combi)

                valid_nc_files = self.get_valid_dates(sftp.listdir(remote_dir))

                # warning if network + species + resolution + date range combination gets no matching results       
                if not valid_nc_files:                 
                    print_spec = f'{",".join(self.species)} species' if self.species else ""
                    print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
                    print(f"There is no data available from {self.start_date} to {self.end_date} for {network} network {print_spec} {print_res}")

                else:
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

    def download_ghost_network_zenodo(self,network):
        # print current_network
        print(f"\nDownloading {network}...")

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
            res_spec_combinations = [f"/{species}.tar.xz"  for species in self.species]   
        # network + resolution + species 
        else:
            res_spec_combinations = [f"/{resolution}/{species}.tar.xz" for species in self.species for resolution in self.resolution]

       
        # get all the possible networks from network e.g. EBAS gets EBAS-ACTRIS_oyktp, EBAS-AMAP_hcxwm, EBAS-CAMP_daczg, etc
        if network in ["GHOST","EBAS"]:
            current_networks = list(filter(lambda x: network in x, self.ghost_zip_files))
        else:
            current_networks = [network] if network in self.ghost_zip_files else []
            # if long EBAS or similar was passed (due to getting all networks) then path used would be EBAS
            for short_network in ["GHOST","EBAS"]:
                if short_network in network:
                    network = short_network
                    break

        # TODO MAYBE CHANGE TO WARNING
        if not current_networks:
            print(f"There is no data available for {network} network.")
            return

        for curr_net in current_networks:
            # get url to download the zip file for the current network
            zip_url = self.ghost_zip_files[curr_net]

            # get all the species tar files which fulfill the combination condition
            res_spec_final_combinations = []
            with RemoteZip(zip_url) as zip:
                for combi in res_spec_combinations:
                    res_spec_final_combinations += list(filter(lambda x: combi in x, zip.namelist()))
                res_spec_final_combinations = list(filter(lambda x: x[-7:] == '.tar.xz', res_spec_final_combinations))
                
                # warning if network + species + resolution combination is gets no matching results
                if not res_spec_final_combinations:
                    print_spec = f'{",".join(self.species)} species' if self.species else ""
                    print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
                    print(f"There is no data available for {network} network {print_spec} {print_res}")
                    return
                
                # TODO Deberia printear cuanto va a ocupar pero no se sabe porque para cada uno no se va a descargar todas las fechas
                # SE PUEDE PONER DESPUES

                # check if there's any possible combination with user's network, resolution and species
                if res_spec_final_combinations:
                    # Print the species, resolution and network combinations that are going to be downloaded
                    print(f"{network} observations to download:")
                    last_ghost_version = 1.5
                    for species in res_spec_final_combinations:
                        hourly_specie_print = "/".join(species.split("/")[1:])
                        print(f"  - {self.ghost_root}/{network}/{last_ghost_version}/{hourly_specie_print[:-7]}")

                    # extract species from zip files
                    for specie_to_get in tqdm(res_spec_final_combinations,desc="Downloading Observations"):
                        zip.extract(specie_to_get,self.ghost_root)
                        
                        # get path and the name of the directory of the tar file
                        tar_path = os.path.join(self.ghost_root, network, str(last_ghost_version), "/".join(specie_to_get.split("/")[1:]))
                        # tar_path = os.path.join(self.ghost_root, specie_to_get)
                        
                        # TODO ESTOY AKI POR ALGUNA RAZON NO EXISTE PERO DESPUES SI, VER XQ
                        # NO SE PUEDE PEDIR A DENE QUE PONGA EL 1.5 EN LA RUTA
                        # A VER SE PUEDE PERO TENDIRAMOS QUE CREAR Y BORARR LA CARPETA
                        tar_dir = os.path.dirname(tar_path) 
                        print(os.path.exists(tar_path))
                        print(list(tar_dir))

                        if not os.path.exists(tar_dir):
                            os.makedirs(tar_dir)

                        print(tar_path)

                        # extract nc file from tar file
                        with tarfile.open(tar_path) as tar_file:
                            tar_names = tar_file.getnames()

                            # get the nc files that are between the start and end date
                            valid_nc_file_names = self.get_valid_dates(tar_names)
                            valid_nc_files = [tar_member for tar_member in tar_file.getmembers() if tar_member.name in valid_nc_file_names]

                            tar_file.extractall(path = tar_dir, members = valid_nc_files)
                        
                        # remove the tar file
                        os.remove(tar_path)

                        # warning if network + species + resolution + date range combination gets no matching results                        
                        if not valid_nc_files:
                            print_spec = f'{",".join(self.species)} species' if self.species else ""
                            print_res = f'at {",".join(self.resolution)} resolutions' if self.resolution else ""
                            print(f"There is no data available from {self.start_date} to {self.end_date} for {network} network {print_spec} {print_res}")
                            
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
        # get user input to know which kind of network wants
        download_source = None
        while not download_source:
            download_source = input("Do you want to download GHOST, non-GHOST or all networks? G/N/A ")
            download_source = download_source.upper() if download_source.upper() in ["G","N","A"] else None

        # TODO change this to new way of getting ghost networks
        if download_source in ["G","A"]:
            if not self.ghost_zip_files: 
                self.get_ghost_zip_files()
            
            self.network = list(self.ghost_zip_files.keys())

        if download_source in ["N","A"]:
            self.nonghost_observation_data = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/nonghost_networks.yaml')))
            self.network = list(self.nonghost_observation_data.keys())

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