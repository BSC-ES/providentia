""" Class for downloading and formatting data from Zenodo """

import os
import shutil
import sys
import time

from datetime import datetime
import json
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
            error = f"Error {resp.status_code}: {resp.text}"
            self.download_instance.logger.info(error)
            sys.exit(1)

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
        zip_path = join(self.temp_dir, f"{network}.zip")

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

    def extract_zip(self, file_paths, zip_path, initial_check):
        """
        Extract certain tar.xz of a ZIP archive into a temporary directory.

        Parameters
        ----------
        file_paths : list of str
            List of file paths used to determine which archives must be extracted.
        zip_path : str
            Absolute path to the ZIP file.
        initial_check : bool
            Flag indicating whether missing archives should deactivate the download
            process during the initial validation step.

        Returns
        -------
        valid_files_info : dict
            Dictionary where each key is the final directory, and
            each value is a dictionary containing metadata for each successfully
            extracted nc file. The metadata includes network, resolution, species,
            filenames in the format of `species/species_YYYYMMM.nc` and tar_member
            with relative path of the `.tar.xz` archive inside the ZIP file.
        """

        files_info = {}

        # fill the dictionary with the metadata
        for path in file_paths:
            dir = join('/', *path.split('/')[:-3], *path.split('/')[-3:-1])

            if dir not in files_info:
                files_info[dir] = {
                    "network": path.split('/')[-4],
                    "resolution": path.split('/')[-3],
                    "species": path.split('/')[-2],
                    "tar_member": f"{path.split('/')[-5]}/{path.split('/')[-3]}/{path.split('/')[-2]}.tar.xz",
                    "filenames": [join(*path.split('/')[-2:])],
                }
            else:
                files_info[dir]["filenames"].append(join(*path.split('/')[-2:]))

        valid_files_info = {}

        # open zip file
        with ZipFile(zip_path) as zipf:
            zip_members = set(zipf.namelist())

            # iterate through the different nc files inside the zip and validate each one of them
            for dir, file_info_dict in files_info.items():
                for filename in file_info_dict['filenames']:
                    
                    tar_member = file_info_dict["tar_member"]

                    if tar_member in zip_members:
                        
                        # include the tar files that were on the ZIP file on the final directory
                        if dir not in valid_files_info:
                            valid_files_info[dir] = {
                                "network": file_info_dict["network"],
                                "resolution": file_info_dict["resolution"],
                                "species": file_info_dict["species"],
                                "tar_path": join(self.temp_dir, tar_member),
                                "filenames": [filename],
                            }
                            
                            # extract valid tar member 
                            zipf.extract(tar_member, self.temp_dir)
                        
                        else:
                            valid_files_info[dir]["filenames"].append(filename)
                    
                    else:
                        # throw a warning if nc file is not found in the zip file
                        msg = f"Missing archive in ZIP: {tar_member}"
                        show_message(self.download_instance, msg, deactivate=initial_check)

        return valid_files_info

    def extract_tar(self, files_info):
        """
        Extract selected files from previously extracted `.tar.xz` archives.

        Parameters
        ----------
        files_info : dict
            Dictionary where each key is the final directory, and
            each value is a dictionary containing metadata for each successfully
            extracted nc file. The metadata includes network, resolution, species,
            filenames in the format of `species/species_YYYYMMM.nc` and tar_path
            with the absolute path of the `.tar.xz` archive inside the ZIP file.
        """
        self.download_instance.logger.info("\n    Extracting TAR files to:")

        # iterate through each file
        for dir, file_info_dict in files_info.items():
            self.download_instance.logger.info(f"\n    - {dir}")

            for filename in tqdm(file_info_dict["filenames"], bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}', desc=f"      TAR extraction in progress... ({len(file_info_dict['filenames'])})"):
                # open tar file
                with tarfile.open(file_info_dict["tar_path"]) as tar:
                    try:
                        member = tar.getmember(filename)
                    except KeyError:
                        # throw a warning if nc file is not found in the tar file
                        msg = f"Missing file in TAR: {filename}"
                        show_message(self.download_instance, msg)
                        continue

                    # create directory where the tar file will be extracted
                    os.makedirs(dir, exist_ok=True)
                    
                    # extract valid nc file
                    with tar.extractfile(member) as src:
                        with open(join(os.path.dirname(dir), filename), "wb") as dst:
                            shutil.copyfileobj(src, dst)

    def filter_dates(self, dates):
        """
        Filter a list of dates based on the start and end date of the download instance.

        Parameters
        ----------
        dates : list of str
            List of dates in 'YYYYMMDD' string format to be filtered.

        Returns
        -------
        valid_dates : list of str
            List of dates that are within the range specified by 
            `self.download_instance.start_date` and `self.download_instance.end_date`.
        """

        # transform YYYYMMDD into datetime format
        start = datetime.strptime(self.download_instance.start_date, "%Y%m%d")
        end = datetime.strptime(self.download_instance.end_date, "%Y%m%d")

        valid_dates = []
        for date in dates:
            # transform YYYMM into datetime and append '01'
            d = datetime.strptime(date + "01", "%Y%m%d")
            if start <= d < end:
                valid_dates.append(date)

        return valid_dates
        
    def check_filetrees(self, network):
        """
        Check the filetree for valid resolutions, species and dates for a given network.

        Parameters
        ----------
        network : str
            The network identifier used to select the subset of the filetree 
            JSON corresponding to that network.

        Returns
        -------
        valid_filetree : dict of dict
            Nested dictionary where keys are valid resolutions and values are 
            dictionaries mapping species names to lists of valid dates. 
            The structure is `{network: resolution: {species: [date]}}}`
        """

        # load filetree JSON for the current GHOST version
        filetree_path = join(PROVIDENTIA_ROOT, 'settings', 'internal', 'zenodo', f'zenodo_ghost_filetree_{self.download_instance.ghost_version}.json')
        with open(filetree_path) as f:
            zenodo_filetrees = json.load(f)

        valid_filetree = {network : {}}

        # determine which resolution to process
        valid_resolutions = self.download_instance.resolution if self.download_instance.resolution != '' else zenodo_filetrees[network].keys()
        for resolution in valid_resolutions:
            if resolution in zenodo_filetrees[network]:
                # determine which species to process
                valid_species = self.download_instance.species if self.download_instance.species != '' else zenodo_filetrees[network][resolution].keys()
                for species in valid_species:
                    if species in zenodo_filetrees[network][resolution]:
                        # obtain dates in the date range
                        dates = zenodo_filetrees[network][resolution][species]
                        valid_dates = self.filter_dates(dates)
                        
                        # add resolution, species and dates to the output if needed
                        if valid_dates:
                            if resolution not in valid_filetree[network]:
                                valid_filetree[network][resolution] = {}
                        else:
                            msg =  f"No valid dates for species '{species}', resolution '{resolution}', "
                            f"network '{network}' in the requested date range."
                            show_message(self.download_instance, msg)

                        valid_filetree[network][resolution][species] = valid_dates
                    else:
                        msg = f"Species '{species}' not found for network '{network}', resolution '{resolution}'. Skipping."
                        show_message(self.download_instance, msg)
            else:
                msg = f"Resolution '{resolution}' not available for network '{network}'. Skipping."
                show_message(self.download_instance, msg)
        
        return valid_filetree

    def filetree_to_paths(self, filetree):
        """
        Convert the filetree dictionary into a list of file paths.

        Parameters
        ----------
        filetree : dict
            Nested dictionary where keys are valid resolutions and values are 
            dictionaries mapping species names to lists of valid dates. 
            The structure is `{network: resolution: {species: [date]}}}`

        Returns
        -------
        paths : list of str
            List of absolute file paths.
        """

        paths = []

        for network, resolutions in filetree.items():
            for resolution, species_dict in resolutions.items():
                for species, dates in species_dict.items():
                    for date in dates:
                        filename = f"{species}_{date}.nc"
                        path = join(self.download_instance.ghost_root, network, self.download_instance.ghost_version , resolution, species, filename)
                        paths.append(path)

        return paths

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
   
        # exit if network is not uploaded
        if network not in self.artifact_mapping.keys():
            msg = f"Network '{network}' is not available for GHOST version {self.download_instance.ghost_version} on Zenodo."
            show_message(self.download_instance, msg, deactivate=initial_check)
            return
 
        if initial_check:
            # print current network
            self.download_instance.logger.info('\n'+'-'*40)
            self.download_instance.logger.info(f"\nDownloading GHOST {network} network data from Zenodo...")
   
            # obtain the filetree that matches with the configuration file
            valid_filetree = self.check_filetrees(network)

            # convert the filetree to absolute paths
            initial_check_nc_files = self.filetree_to_paths(valid_filetree)

            return initial_check_nc_files
        
        elif files_to_download:
                # get the GHOST artifact value for the corresponding network
                artifact_network = self.artifact_mapping[network]    

                # create temporal dir to store the zip file and its tar components
                self.temp_dir = os.path.join(self.download_instance.ghost_root, ".temp")
                os.makedirs(self.temp_dir, exist_ok=True)

                # download zip on the temporal directory
                zip_path = self.download_zip(network, artifact_network)

                # extract zip on the temporal directory
                valid_files_info = self.extract_zip(files_to_download, zip_path, initial_check)

                # extract tar on the temporal directory
                self.extract_tar(valid_files_info)                    
    
                # remove the temporal directory
                shutil.rmtree(self.temp_dir)