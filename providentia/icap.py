from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import requests
import os
from tqdm import tqdm
import yaml
from netCDF4 import Dataset
import numpy as np
import shutil

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
ghost_tropopause_variables = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "internal", "cams", "tropopause.yaml"))
)


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

    def find_mode(self):
        # look for model before download
        pass



    def download_non_interpolated_model(
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
            self.logger.info("\n" + "-" * 40)
            self.logger.info(
                f"\nDownloading {model} model data from the U.S. Naval Research Laboratory..."
            )

            # check if ssh exists and check if still active, connect if not
            if (self.ssh is None) or (self.ssh.get_transport().is_active()):
                self.connect()

            # look for the model in the remote machine
            model_exists, remote_dir = find_model(self, mod_id, domain, initial_check)

            if not model_exists:
                return

            # valid_resolutions = find_available_resolutions(
            #     self, remote_dir, mod_id, domain, initial_check
            # )

            # valid_resolutions_species_dir = find_available_species(
            #     self,
            #     valid_resolutions,
            #     remote_dir,
            #     ensemble,
            #     mod_id,
            #     domain,
            #     initial_check,
            # )

            path_files_dict = find_available_nc_files(
                self, valid_resolutions_species_dir, initial_check
            )

            initial_check_nc_files = build_nc_file_paths_in_range(
                self, path_files_dict, initial_check
            )

            return initial_check_nc_files

        elif files_to_download:
            self.logger.info(
                f"\n{model} model data to download ({len(files_to_download)}):"
            )

            download_non_interpolated_sftp(self, files_to_download)

        else:
            # tell the user if not valid resolution specie date combinations
            msg = "There is no available model output to be downloaded."
            show_message(self, msg, deactivate=initial_check)

        return initial_check_nc_files
