""" Class for downloading and formatting data from Zenodo """

import os
import shutil
import sys
import time

import requests
import tarfile
from tqdm import tqdm
import yaml 
from zipfile import ZipFile

from providentia.auxiliar import CURRENT_PATH, join
from .warnings_prv import show_message

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)
zenodo_dois = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'zenodo', 'zenodo_dois.yaml')))

class Zenodo:
    """Class responsible for handling GHOST network data downloads from Zenodo."""

    def __init__(self, download_instance):
        """
        Initialize the Zenodo download helper with a download instance.

        Parameters
        ----------
        download_instance : Download
            Stores the instance of the 'download' class.
        """
        
        self.download_instance = download_instance

        # get url for the zenodo GHOST repository 
        if self.download_instance.ghost_version not in zenodo_dois:
            error = (f"Error: Current GHOST version ({self.download_instance.ghost_version}) is not available on Zenodo. "
                    f"Please choose one of the available versions: {tuple(zenodo_dois.keys())}.")            
            self.download_instance.logger.error(error)
            sys.exit(1)

        # load zenodo artifact mapping
        self.artifact_mapping = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'zenodo', f'zenodo_{self.download_instance.ghost_version}.yaml')))

        # create session
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Providentia/3.0 (contact: paula.serrano@bsc.es)",
            "Accept": "application/json",
        })

    def get_zip(self):
        """
        Perform a HTTP GET request against the Zenodo REST API.

        The response content is read in chunks using 
        `Response.iter_content`, which is useful for large files.

        Returns
        -------
        requests.Response
            HTTP response object for the request.
        """

        # initialize retry and timing controls
        retries = 0
        min_interval = 0.2
        last_request = 0

        # keep trying until request succeeds or fails definitively
        while True:
            # respect the minimum interval between requests
            delta = time.time() - last_request
            if delta < min_interval:
                time.sleep(min_interval - delta)

            # send GET request with streaming enabled
            resp = self.session.get(self.url, stream=True, timeout=30)
            last_request = time.time()

            # if request succeeds, return the response
            if resp.status_code == 200:
                return resp
            
            # if server asks us to wait, retry after sleeping
            if resp.status_code in (429, 502, 503) and retries < 3:
                retries += 1
                time.sleep(2 ** retries)
                continue
            
            # for other errors, print response content and raise exception
            print(f"Error {resp.status_code}: {resp.text}")
            resp.raise_for_status() 

    def download_zip(self, network, artifact_network):
        """
        Download the ZIP file from Zenodo in chunks and a 
        tqdm progress bar shows the download progress.

        Parameters
        ----------
        network : str
            Name of the network to download.
        artifact_network : str
            Name of the network to download with the Zenodo artifact.

        Returns
        -------
        zip_path : str
            Absolute path to the zip file.
        """

        # get the record id for the current GHOST version
        record_id = zenodo_dois[self.download_instance.ghost_version]  

        # get url to download the zip file for the current network
        self.url = f"https://zenodo.org/api/records/{record_id}/files/{artifact_network}.zip/content"

        # get path were zip file is going to get downloaded
        zip_path = f"{self.temp_dir}/{network}.zip"

        # open streaming GET request to the file
        with self.get_zip() as r:
            # extract total file size from headers for progress bar
            total = int(r.headers.get("Content-Length", 0))

            # create the folder if it does not exist
            os.makedirs(os.path.dirname(zip_path), exist_ok=True)

            # open the file and create the tqdm object
            with open(zip_path, "wb") as f, tqdm(
                total=total, unit="B", desc=f"    Downloading {network}.zip", unit_scale=True) as pbar:
                # write the file content in chunks (64 KB each)
                for chunk in r.iter_content(chunk_size=65536):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))

        return zip_path

    def extract_zip(self, network, zip_path):
        """
        Extract the contents of a ZIP archive into a temporary directory.

        Parameters
        ----------
        network : str
            Name of the network to download.
        zip_path : str
            Absolute path to the zip file.
        """
        
        # create directory where the contents of the zip file will go to
        extracted_zip_path = f"{self.temp_dir}/{network}"

        # open the ZIP file
        with ZipFile(zip_path, "r") as zipf:            
            zipf.extractall(extracted_zip_path)

    def extract_tar(self):
        pass

    def download_ghost_network_zenodo(self, network, initial_check, files_to_download=None):
        """
        Download GHOST network data from Zenodo.

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

        if not initial_check:
            # print current_network
            self.download_instance.logger.info('\n'+'-'*40)
            self.download_instance.logger.info(f"\nDownloading GHOST {network} network data from Zenodo...")
        
        # exit if network is not uploaded
        if network not in self.artifact_mapping.keys():
            msg = f"Network '{network}' is not available for GHOST version {self.download_instance.ghost_version} on Zenodo."
            show_message(self.download_instance, msg, deactivate=initial_check)
            return

        # get the GHOST artifact value for the corresponding network
        artifact_network = self.artifact_mapping[network]    

        # create temporal dir to store the zip file and its tar components
        self.temp_dir = os.path.join(self.download_instance.ghost_root, ".temp")
        os.makedirs(self.temp_dir, exist_ok=True)

        # download zip on the temporal directory
        zip_path = self.download_zip(network, artifact_network)

        # extract zip on the temporal directory
        self.extract_zip(network, zip_path)

        # extract tar on the temporal directory
        self.extract_tar()

        # get resolution and/or species combinations
        # network
        if not self.download_instance.resolution and not self.download_instance.species:
            res_spec_combinations = [""]
        # network + resolution 
        elif not self.download_instance.species:
            res_spec_combinations = [f"/{resolution}/" for resolution in self.download_instance.resolution]
        # network + species 
        elif not self.download_instance.resolution:
            res_spec_combinations = [f"/{species}.tar.xz"  for species in self.download_instance.species]   
        # network + resolution + species 
        else:
            res_spec_combinations = [f"/{resolution}/{species}.tar.xz" for species in self.download_instance.species for resolution in self.download_instance.resolution]

        # get all the species tar files which fulfill the combination condition
        res_spec_dir_tail = []
        with RemoteZip(zip_file_url) as zip:
            for combi in res_spec_combinations:
                res_spec_dir_tail += list(filter(lambda x: combi in x, zip.namelist()))
            res_spec_dir_tail = list(filter(lambda x: x[-7:] == '.tar.xz', res_spec_dir_tail))

            # warning if network + species + resolution combination is gets no matching results
            if not res_spec_dir_tail:
                print_spec = f'{",".join(self.download_instance.species)} species' if self.download_instance.species else ""
                print_res = f'at {",".join(self.download_instance.resolution)} resolution(s)' if self.download_instance.resolution else ""
                msg = f"There is no data available in Zenodo for {network} network for {print_spec} {print_res}"
                show_message(self.download_instance, msg, deactivate=initial_check)

            # check if there's any possible combination with user's network, resolution and species
            else:
                # initialise list with all the nc files to be downloaded
                initial_check_nc_files = []

                # print the species, resolution and network combinations that are going to be downloaded
                if not initial_check:
                    self.download_instance.logger.info(f"\n{network} observations to download:")
                
                for remote_dir_tail in res_spec_dir_tail:
                    resolution, species = remote_dir_tail.split("/")[1:]
                    species = species[:-7]
                    local_dir = join(self.download_instance.ghost_root, network, self.download_instance.ghost_version, resolution, species)

                    #  print the species, resolution and network combinations that are going to be downloaded
                    if not initial_check:
                        self.download_instance.logger.info(f"\n  - {local_dir}")

                    # create temporal dir to store the middle tar file with its directories
                    self.temp_dir = join(self.download_instance.ghost_root,'.temp')
                    if not os.path.exists(self.temp_dir):
                        os.mkdir(self.temp_dir)

                    zip.extract(remote_dir_tail,self.temp_dir)
                    
                    # get path and the name of the directory of the tar file
                    tar_path = join(self.download_instance.ghost_root, network, str(self.download_instance.ghost_version), *remote_dir_tail.split("/")[1:])
                    temp_path = join(self.temp_dir,remote_dir_tail)

                    # extract nc file from tar file
                    with tarfile.open(temp_path) as tar_file:
                        # get the nc files that are between the start and end date
                        tar_names = tar_file.getnames()
                        valid_nc_file_names = self.download_instance.get_valid_nc_files_in_date_range(tar_names)
                        
                        # warning if network + species + resolution + date range combination gets no matching results                        
                        if not valid_nc_file_names:
                            print_spec = f'{",".join(self.download_instance.species)} species' if self.download_instance.species else ""
                            print_res = f'at {",".join(self.download_instance.resolution)} resolutions' if self.download_instance.resolution else ""
                            msg = f"There is no data available in Zenodo from {self.download_instance.start_date} to {self.download_instance.end_date} for {network} network for {print_spec} species at {print_res} resolution"
                            show_message(self.download_instance, msg, deactivate=initial_check)
                            continue
                        # download the valid resolution species date combinations
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
                                    show_message(self.download_instance, msg, deactivate=initial_check)     
                                    continue  
                                
                            valid_nc_files = list(filter(lambda x:x.name in valid_nc_file_names, tar_file.getmembers()))
                            
                            # sort nc_files
                            valid_nc_files.sort(key = lambda x:x.name)   

                            if not initial_check and not self.download_instance.logfile:
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
                                    self.download_instance.latest_nc_file_path = local_path

                                    # extract the file
                                    tar_file.extract(member = nc_file, path = tar_dir)
                
                # remove the temp directory
                shutil.rmtree(os.path.dirname(os.path.dirname(temp_path)))
                
                return initial_check_nc_files