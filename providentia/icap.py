from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import requests
import os
from tqdm import tqdm
import yaml
from netCDF4 import Dataset
import numpy as np
import shutil
import pandas as pd

from .warnings_prv import show_message

from providentia.auxiliar import CURRENT_PATH, join

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

cams_species_units = yaml.safe_load(
    open(
        join(
            PROVIDENTIA_ROOT, "settings", "internal", "cams", "cams_species_units.yaml"
        )
    )
)

mapping_species = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "mapping_species.yaml"))
)

class ICAP(object):
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

    def control_domain(self, domain, initial_check):
        # TODO HARDCODED
        correct_resolution = "6hourly"
        
        # check if the domain is the correct one for the dataset
        if domain not in cams_options[prefix]:
            possible_domains = "', '".join(cams_options[prefix])
            msg = (
                f"The current domain '{domain}' is not valid for the CAMS dataset. "
                f"It must be '{possible_domains}'."
            )
            show_message(self.download_instance, msg, deactivate=initial_check)
            return

    def control_resolutions(self, initial_check):
        resolution_list = (
            self.download_instance.model_resolution
            if self.download_instance.model_resolution
            else self.download_instance.resolution
        )
        
        # TODO HARDCODED
        correct_resolution = "6hourly"

        final_resolution_list = []

        for resolution in resolution_list:
            if correct_resolution == resolution:
                final_resolution_list.append(resolution)
            else:
                msg = f"The current resolution '{resolution}' is not valid. It must be '{correct_resolution}'."
                show_message(
                    self.download_instance, msg, deactivate=initial_check
                )

        return final_resolution_list

    def control_species(self, initial_check):

        # TODO HARDCODED
        correct_species = "od550aero"

        final_species_list = []

        for species in self.download_instance.species:
            if correct_species == species or correct_species in mapping_species[correct_species]:
                final_species_list.append(correct_species)
            else:
                msg = f"The species '{species}' is not available in CAMS."
                show_message(
                    self.download_instance, msg, deactivate=initial_check
                )

        return final_species_list
        
    def control_dates(self):
        pass

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
        input_file = Dataset(input_filepath, "r", format="NETCDF4")

        # create new netcdf file
        output_file = Dataset(output_filepath, "w", format="NETCDF4")
        output_file.set_auto_mask(True)

        # change the last downloaded file
        self.download_instance.latest_nc_file_path = output_filepath

        # create the dimensions
        for input_dim_name, output_dim_name in input_file.dimensions.items():
            output_file.createDimension(input_dim_name, output_dim_name.size)

        # copy variables
        for input_var_name in input_file.variables:
            # get the output var name
            if input_var_name in ["lat", "lon", "time"]:
                if input_var_name == "lat":
                    output_var_name = "latitude"
                elif input_var_name == "lon":
                    output_var_name = "longitude"
                else:
                    output_var_name = input_var_name

            elif input_var_name == ghost_tropopause_variables[species]:
                output_var_name = species

            else:
                continue

            # get the variable
            input_var = input_file[input_var_name]

            # create the variable
            output_var = output_file.createVariable(
                output_var_name, input_var.datatype, input_var.dimensions
            )

            # add calendar and units attributes to the time variable
            if output_var_name == "time":
                output_var.setncattr("calendar", "standard")
                output_var.setncattr("units", input_var.units)

            # add coordinates, grid_mapping and units to the species variable
            elif output_var_name == species:
                output_var.setncattr("coordinates", "latitude longitude")
                output_var.setncattr("grid_mapping", "crs")
                output_var.setncattr("units", cams_species_units[input_var.units])

            if output_var_name == "time":
                data = np.arange(len(input_var[:]))
            else:
                data = input_var[:]

            # add the data to the variable
            output_var[:] = data

        # add grid_mapping
        output_file[species].setncattr("grid_mapping", "crs")

        # add coordinates
        output_file[species].setncattr("coordinates", "lat lon")

        # add crs
        crs_var = output_file.createVariable("crs", "u1")
        crs_var.setncatts(
            {
                "grid_mapping_name": "latitude_longitude",
                "semi_major_axis": 6371000.0,
                "inverse_flattening": 0.0,
            }
        )

        # close the original and new netcdf files
        output_file.close()
        input_file.close()

    def build_nc_file_paths_in_range(self):
        # look for model before download

        model = "cams_icap"
        domain = "global"

        path = join(self.download_instance.mod_to_interp_root, model, domain, resolution, species)

        months = pd.date_range(start=self.download_instance.start_date, end=self.download_instance.end_date, freq="MS")

        initial_check_nc_files = {
            path: {
                "nc_files": [f"{species}_{date}.nc" for date in months.strftime("%Y%m")]
            }
        }

        return initial_check_nc_files
    
    def find_available_dates(self):
        # TODO Hardcoded
        first_start_date = datetime(2014, 11, 1)
        last_end_date = datetime.today()
        
    # TODO change name
    def download_ICAP(self, files_to_download_info):
        # TODO hardcoded
        ORIG = "https://nrlgodae1.nrlmry.navy.mil/ftp/outgoing/nrl/ICAP-MME/"

        # get directory structure
        dir_tail = join(config_modid, domain, resolution, ghost_species)

        # temporal directory for the zip file
        temp_root_dir = join(self.download_instance.mod_to_interp_root, ".temp")
        temp_dir = join(temp_root_dir, dir_tail)

        ODIR = "./recon/aemet/icap/original_files"
        MEANDIR = "./recon/aemet/icap/6hourly/mean"

        for path in [ODIR, MEANDIR]:
            os.makedirs(path, exist_ok=True)

        current = self.download_instance.start_date
        
        while current <= end:

            date = current.strftime("%Y%m%d")
            date1 = current.strftime("%Y-%m-%d")

            Y = current.strftime("%Y")
            M = current.strftime("%m")

            nam0 = f"icap_{date}00_C4_dustaod550.nc"
            nam1 = f"{date}_ICAP-MME_MEAN.nc"

            # Download source file
            source_file = os.path.join(ODIR, nam0)

            if not os.path.exists(source_file):

                url = f"{ORIG}{Y}/{Y}{M}/{nam0}"
                print(f"Downloading {url}")

                response = requests.get(url, stream=True)
                response.raise_for_status()

                with open(source_file, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)

            else:
                print("File already downloaded")

            print("Running formatting scripts")

            # Mean file
            mean_file = os.path.join(MEANDIR, nam1)

            if not os.path.exists(mean_file):

                fmt(date, str(ODIR), str(MEANDIR), nam1, date1)

            else:
                print("Mean file already formatted")

            current += timedelta(days=1)

    def download_ICAP_model(
        self, model, initial_check, files_to_download=None
    ):
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

        # get model id and the domain
        mod_id, domain, ensemble = model.split("-")

        if initial_check:
            # print current model
            self.download_instance.logger.info("\n" + "-" * 40)
            self.download_instance.logger.info(
                f"\nDownloading {model} model data from the U.S. Naval Research Laboratory..."
            )

            self.control_domain(initial_check)

            valid_resolutions = self.control_resolutions(initial_check)

            if not valid_resolutions:
                return

            valid_resolutions_species_dir = self.control_species()

            cams_start_date, cams_end_date = self.control_dates()

            initial_check_nc_files = self.build_nc_file_paths_in_range()

            return initial_check_nc_files

        elif files_to_download:
            self.download_instance.logger.info(
                f"\n{model} model data to download ({len(files_to_download)}):"
            )

            self.download_ICAP(files_to_download)

        else:
            # tell the user if not valid resolution specie date combinations
            msg = "There is no available model output to be downloaded."
            show_message(self, msg, deactivate=initial_check)
