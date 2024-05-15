import sys
import os

import requests, zipfile
from io import BytesIO

# urlparse
#from tqdm import tqdm
# from remotezip import RemoteZip


from .configuration import ProvConfiguration, load_conf
from .read_aux import check_for_ghost

class ProvidentiaDownload(object):
    def __init__(self,**kwargs):
        # initialise zenodo url
        self.ghost_url = 'https://zenodo.org/records/10637450'

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

        for section_ind, section in enumerate(self.parent_section_names):
            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # update self with section variables
            for k, val in self.section_opts.items():
                setattr(self, k, provconf.parse_parameter(k, val))

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
            wanted_directories = [""]
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
        # network
        if not self.resolution and not self.species:
            wanted_directories = [""]
        # network + resolution 
        elif not self.species:
            wanted_directories = [f"{network}/{self.resolution}/"]
        # network + species 
        elif not self.resolution:
            wanted_directories = self.species
        # network + resolution + species 
        else:
            wanted_directories = [f"{network}/{self.resolution}/{specie}" for specie in self.species]

        ### TODO: Copiar todo lo de en sucio cuando tenga las librerias que faltan.

    def download_exp(self):
        pass    

def main(**kwargs):
    """ Main function when running download function. """
   
    ProvidentiaDownload(**kwargs)