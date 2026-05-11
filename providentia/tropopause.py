from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import requests
import re
import os
from tqdm import tqdm
import yaml
from netCDF4 import Dataset
import numpy as np
import shutil

import requests

from urllib.parse import urljoin

from .warnings_prv import show_message

from providentia.auxiliar import CURRENT_PATH, join
from .warnings_prv import show_message

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

cams_species_units = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'cams_species_units.yaml')))
ghost_tropopause_variables = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'tropopause.yaml')))

class Tropopause(object):
    """
    Class that manages the interaction with Juelich
    to retrieve ERA5 datasets and convert downloaded 
    NetCDF files into the Providentia-compatible format.
    """

    def __init__(self, download_instance):
        """
        Parameters
        ----------
        download_instance : providentia.Download
            Download controller instance. 
        """

        self.download_instance = download_instance
    
    def format_data(self, input_filepath, output_filepath, species): 
        """
        Reformat a raw CAMS NetCDF file into a standardized
        Providentia-compatible NetCDF.

        Parameters
        ----------
        input_filepath : str
            Path to the input CAMS NetCDF file.
        output_filepath : str
            Path where the formatted NetCDF file will be written.
        species : str
            Providentia species name.
        prefix : str
            CAMS dataset prefix (e.g. 'cams_forecast', 'cams_reanalysis').
        domain : str
            Spatial domain, 'global' or 'regional'.
        resolution : str
            Temporal resolution of the data.
        cams_species : str
            CAMS variable name corresponding to the Providentia species.
        url : str
            CAMS dataset URL
        """

        # open original netcdf file      
        input_file = Dataset(input_filepath, 'r', format="NETCDF4") 

        # create new netcdf file
        output_file = Dataset(output_filepath, 'w', format="NETCDF4") 
        output_file.set_auto_mask(True)	

        # change the last downloaded file
        self.download_instance.latest_nc_file_path = output_filepath

        # create the dimensions
        for input_dim_name, output_dim_name in input_file.dimensions.items():
            output_file.createDimension(input_dim_name, output_dim_name.size)

        # copy variables
        for input_var_name in input_file.variables:
            # get the output var name    
            if input_var_name in ['lat', 'lon', 'time']:
                if input_var_name == 'lat':
                    output_var_name = 'latitude'
                elif input_var_name == 'lon':
                    output_var_name = 'longitude'
                else:
                    output_var_name = input_var_name 

            elif input_var_name == ghost_tropopause_variables[species]:
                output_var_name = species

            else:
                continue
            
            # get the variable
            input_var = input_file[input_var_name]
            
            # create the variable
            output_var = output_file.createVariable(output_var_name, input_var.datatype, input_var.dimensions)

            # add calendar and units attributes to the time variable
            if output_var_name == "time":
                output_var.setncattr('calendar', 'standard')
                output_var.setncattr('units', input_var.units)

            # add coordinates, grid_mapping and units to the species variable
            elif output_var_name == species:               
                output_var.setncattr('coordinates', 'latitude longitude')
                output_var.setncattr('grid_mapping', 'crs')
                output_var.setncattr('units', cams_species_units[input_var.units])
            
            if output_var_name == "time":  
                data = np.arange(len(input_var[:]))
            else:
                data = input_var[:]

            # add the data to the variable
            output_var[:] = data
        
        # add grid_mapping
        output_file[species].setncattr('grid_mapping', 'crs')
        
        # add coordinates
        output_file[species].setncattr('coordinates', 'lat lon')

        # add crs
        crs_var = output_file.createVariable('crs', 'u1')  
        crs_var.setncatts({
            'grid_mapping_name': 'latitude_longitude',
            'semi_major_axis': 6371000.0,
            'inverse_flattening': 0.0
        })
        
        # close the original and new netcdf files
        output_file.close()
        input_file.close()     

    def download_tropopause_model(self, model, initial_check, files_to_download=None): 
        """
        Download and process tropopause model data. Validates the
        requested model configuration, checks date availability
        downloads CAMS NetCDF files and reformats them into
        Providentia-compatible NetCDF outputs.

        Parameters
        ----------
        model : str
            CAMS model specification string in the form
            '<dataset>_<model>_<stream>-<domain>-<ensemble>'.
        initial_check : bool
            If True, perform validation and file discovery without downloading
            or formatting data.
        files_to_download : list of str, optional
            Subset of files to download, used to restrict downloads to specific
            dates or files.

        Returns
        -------
        initial_check_nc_files : list of str or None
            List of expected output file paths when 'initial_check' is True.
            Returns None otherwise.
        """

        if not initial_check:
            # print current model
            self.download_instance.logger.info('\n'+'-'*40)
            self.download_instance.logger.info(f"\nDownloading {model} model data from SDL Climate Science...")

        # get model id and the domain
        config_modid, domain, ensemble_options = model.split("-")

        # make the necessary checks to the dates
        # TODO HARD CODED
        min_start_date = datetime(2000, 1, 1)
        max_end_date = datetime(2020, 10, 31)

        cams_start_date = datetime.strptime(self.download_instance.start_date, "%Y%m%d")
        cams_end_date = datetime.strptime(self.download_instance.end_date, "%Y%m%d") - timedelta(days=1)

        # if the minimum date is over the end date
        if min_start_date > cams_end_date or max_end_date < cams_start_date:
            msg = f"The selected dates are unavailable. Please choose dates between {min_start_date.strftime('%Y-%m-%d')} and {max_end_date.strftime('%Y-%m-%d')}."
            show_message(self.download_instance, msg, deactivate=initial_check)
            return 
        
        # initialise list with all the nc files to be downloaded
        initial_check_nc_files = []

        # get model resolution
        resolution_list = self.download_instance.model_resolution if self.download_instance.model_resolution else self.download_instance.resolution

        # iterate through the resolutions
        for resolution in resolution_list:

            # get the resolution for the cams dataset
            correct_resolution = "hourly" # TODO HARDCODED
            
            # check if the resolution is the correct one for the dataset TODO change for the finer resolutions
            if resolution != correct_resolution:
                msg = f"The current resolution '{resolution}' is not valid. It must be '{correct_resolution}'."            
                show_message(self.download_instance, msg, deactivate=initial_check)
                continue

            valid_species = []

            # iterate through the species
            for ghost_species in self.download_instance.species: 

                # check if species is in the ghost_cams_variables file
                if ghost_species not in ghost_tropopause_variables.keys():
                    msg = f"The species '{ghost_species}' is not available in the ERA5 tropopause dataset."
                    show_message(self.download_instance, msg, deactivate=initial_check)
                    continue

                valid_species.append(ghost_species)

            # initialize iterators controlers
            next_cams_date = cams_start_date.replace(day=1) + relativedelta(months=1) - timedelta(days=1)

            all_dates = [cams_start_date + timedelta(days=i) for i in range((next_cams_date - cams_start_date).days + 1)]

            # get directory structure
            dir_tail = join(config_modid, domain, resolution)

            # temporal directory for the zip file
            temp_root_dir = join(self.download_instance.mod_to_interp_root,'.temp')
            temp_dir = join(temp_root_dir, dir_tail)

            if valid_species:
                
                # print progress depending on initial check 
                if initial_check:
                    dates_iterator = all_dates
                else:
                    dates_iterator = tqdm(all_dates, bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"\nDownloading and formating {', '.join(valid_species)} ({len(all_dates)})")

                # iterate through all dates to format each of the day files
                for date in dates_iterator:
                    # create temporal dir to store the middle zip file with its directories
                    os.makedirs(temp_dir, exist_ok=True)

                    temp_file = f"era5_{date.year}_{str(date.month).zfill(2)}_{str(date.day).zfill(2)}.nc"
                    
                    # get temporal path
                    temp_path = join(temp_dir, temp_file)

                    if not initial_check:
                        # create url
                        base_url = "https://datapub.fz-juelich.de/slcs/tropopause/data/v2/era5/"
                        file_url = f"{base_url}{date.year}/{temp_file}"

                        # submit the request
                        r = requests.get(file_url)
                        with open(temp_path, "wb") as out:
                            out.write(r.content)

                    for species in valid_species:

                        # final directory for the nc files
                        final_dir = join(self.download_instance.mod_to_interp_root, dir_tail, species)

                        # get final path
                        file_name = f"{species}_{date.strftime('%Y%m%d')}.nc"
                        final_path = join(final_dir, file_name)
                    
                        # add the ncfiles to the files to download list
                        if initial_check:
                            initial_check_nc_files.append(final_path)

                        else:
                            # create dir
                            os.makedirs(final_dir, exist_ok=True)

                            self.format_data(temp_path, final_path, species)
                        
                        # change the last downloaded file
                        self.download_instance.latest_nc_file_path = "/path/to/file"

        # remove the temp directory tail
        if os.path.exists(temp_root_dir):
            shutil.rmtree(temp_root_dir)
                    
        return initial_check_nc_files