from datetime import date, datetime, timedelta, time
import requests
import os
from tqdm import tqdm
import yaml
from netCDF4 import Dataset
import numpy as np
import pandas as pd

from .warnings_prv import show_message

from providentia.auxiliar import CURRENT_PATH, join

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

mapping_species = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "mapping_species.yaml"))
)

class BDRC(object):
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

        self.bdrc_fixed_values = yaml.safe_load(    
                                    open(join(PROVIDENTIA_ROOT, "settings", "internal",  "bdrc", "bdrc.yaml"))
                                    )

        self.bdrc_password = yaml.safe_load(    
                                    open(join(PROVIDENTIA_ROOT, "settings", "internal",  "bdrc", "bdrc_password.yaml"))
                                    )

        self.bdrc_model = yaml.safe_load(
                                    open(join(PROVIDENTIA_ROOT, "settings", "internal",  "bdrc", "bdrc_model.yaml"))
                                )

        # transform dates into datetime
        for k, v in self.bdrc_fixed_values.items():
            if "date" in k:
                self.bdrc_fixed_values[k] = self.convert_date_to_datetime(v)
    
    def convert_date_to_datetime(self, dt):
        """
        Convert dates to datetimes.

        Supports individual dates, lists of dates, dictionaries containing
        dates and the string ``today``.

        Parameters
        ----------
        dt : datetime.date, str, list, dict
            Value to convert. A string equal to ``today`` is converted to
            the current datetime. Other strings are converted to ``None``.

        Returns
        -------
        datetime.datetime, list, dict or None
            The converted value.
        """

        if isinstance(dt, date) and not isinstance(dt, datetime):
            return datetime.combine(dt, time(0, 0))

        if isinstance(dt, str):
            return datetime.today() if dt == "today" else None

        if isinstance(dt, list):
            return [
                datetime.combine(d, time(0, 0))
                for d in dt
            ]

        if isinstance(dt, dict):
            return {
                model: datetime.combine(d, time(0, 0))
                for model, d in dt.items()
            }

    def control_mod_id(self, mod_id):
        """
        Validate the requested model id for BDRC.

        Unsupported model id is reported to the user.

        Parameters
        ----------
        mod_id : str
            Model ID requested by the user.

        Returns
        -------
        bool
            ``True`` if the requested model ID is supported by BDRC,
            otherwise ``False``.
        """

        # get model IDs from the yaml file
        available_models = self.bdrc_model.keys()
        
        # check if the model ID is the correct one for the dataset
        if mod_id not in available_models:
            msg = (
                f"The current model ID '{mod_id}' is not valid for the BDRC dataset. "
                f"It must be '{', '.join(available_models)}'."
            )
            show_message(self.download_instance, msg)

            return False
        
        return True

    def control_domain(self, domain):
        """
        Validate the requested domain for BDRC.

        Unsupported domain is reported to the user.

        Parameters
        ----------
        domain : str
            Domain requested by the user.

        Returns
        -------
        bool
            ``True`` if the requested domain is supported by BDRC,
            otherwise ``False``.
        """

        # get the domain from the yaml file
        correct_domain = self.bdrc_fixed_values["domain"]
        
        # check if the domain is the correct one for the dataset
        if domain != correct_domain:
            msg = (
                f"The current domain '{domain}' is not valid for the Barcelona Dust Regional Center."
                f"It must be '{correct_domain}'."
            )
            show_message(self.download_instance, msg)

            return False
        
        return True

    def control_resolutions(self):
        """
        Validate the requested resolutions for BDRC.

        If a model-specific resolution is provided, it takes priority over the
        general resolution requested by the user. 

        Unsupported resolutions are reported to the user and are not included
        in the returned list.

        Returns
        -------
        list of str
            List of valid resolutions supported by the BDRC dataset.
        """

        # priorize model resolution in case user passed it
        resolution_list = (
            self.download_instance.model_resolution
            if self.download_instance.model_resolution
            else self.download_instance.resolution
        )
        
        # get the domain from the yaml file
        correct_resolution = self.bdrc_fixed_values["resolution"]

        final_resolution_list = []

        # check if the resolution is the correct one for the dataset
        for resolution in resolution_list:
            if correct_resolution == resolution:
                final_resolution_list.append(resolution)
            else:
                msg = f"The current resolution '{resolution}' is not valid for the Barcelona Dust Regional Center. It must be '{correct_resolution}'."
                show_message(
                    self.download_instance, msg
                )

        return final_resolution_list

    def control_species(self):
        """
        Validate and normalize the requested species for BDRC.

        Some model species may be mapped using ``mapping_species`` 
        and are converted to yaml file BDRC species name.

        Unsupported species are reported to the user and are not included
        in the returned list.

        Returns
        -------
        list of str
            List of unique BDRC species name for the requested
            species.
        """

        # create dictionary with possible species and its mappings
        correct_species = {}
        for species in self.bdrc_fixed_values["species"]:
            correct_species[species] = mapping_species.get(species, [])

            correct_species[species].extend(
                k for k, v in mapping_species.items()
                if species in v
            )

        # make sure that species is either the exact species or mapped
        final_species_list = []
        for species in self.download_instance.species:
            for correct_s, species_list in correct_species.items():
                if species == correct_s or species in species_list:
                    final_species_list.append(correct_s)
                    break
            else:
                show_message(
                    self.download_instance,
                    f"The species '{species}' is not available in the Barcelona Dust Regional Center."
                )

        # remove duplicates
        final_species_list = list(set(final_species_list))

        return final_species_list
        
    def control_dates(self, mod_id):
        """
        Adjust the requested date range to the BDRC limits.

        If the requested period does not overlap with the available BDRC
        period, an error message is displayed and ``(None, None)`` is returned.

        If the requested period partially falls outside the BDRC availability
        period, the dates are clipped to the supported range.

        Parameters
        ----------
        mod_id : str
            Model ID requested by the user.

        Returns
        -------
        start_date : datetime or None
            Validated start date, adjusted to the minimum BDRC date if
            necessary.
        end_date : datetime or None
            Validated end date, adjusted to the maximum BDRC date if
            necessary.
        """

        # get the start and end dates from the yaml file
        min_start_date = self.bdrc_fixed_values["start_date"][mod_id]
        max_end_date = self.bdrc_fixed_values["end_date"]

        end_date = datetime.strptime(self.download_instance.end_date, "%Y%m%d") - timedelta(days=1)
        start_date = datetime.strptime(self.download_instance.start_date, "%Y%m%d")
        
        # if the minimum date is over the end date
        if min_start_date > end_date or max_end_date < start_date:
            msg = f"The selected dates are unavailable. Please choose dates between {min_start_date.strftime('%Y-%m-%d')} and {max_end_date.strftime('%Y-%m-%d')}."
            show_message(self.download_instance, msg)
            return None, None

        # check if the start date is within limits
        if min_start_date > start_date:
            start_date = min_start_date

        # check if the end date is within limits
        if max_end_date < end_date:
            end_date = max_end_date
        
        return start_date, end_date 

    def build_nc_file_paths_in_range(self, start_date, end_date, mod_id, domain, resolutions_list, species_list):
        """
        Build the expected NetCDF file paths for a given date range.

        For each combination of resolution and species, this method creates
        the corresponding local output directory and temporary directory.
        It also generates the expected NetCDF filenames for every day in
        the requested date range.

        Parameters
        ----------
        start_date : datetime-like
            Start date of the requested period.
        end_date : datetime-like
            End date of the requested period.
        mod_id : str
            Model identifier.
        domain : str
            Geographic domain for the model data.
        resolutions_list : list of str
            List of spatial resolutions to process.
        species_list : list of str
            List of species to process.

        Returns
        -------
        dict
            Dictionary keyed by the local directory path. Each entry contains:
            - nc_files : list of str
                NetCDF filenames, one for each day in the date range.
            - temp_dir : str
                Directory where temporary NetCDF files are stored.
        """
            
        # create list with all the DAYS in the date range
        date_list = pd.date_range(start=start_date, end=end_date, freq="D")

        initial_check_nc_files = {}

        # create dictionary with the paths and files
        for resolution in resolutions_list:
            for species in species_list:

                dir_tail = join(
                                    mod_id,
                                    domain,
                                    resolution,
                                    species,
                                )

                local_dir = join(
                                    self.download_instance.mod_to_interp_root,
                                    dir_tail
                                )

                temp_dir = join(
                                    self.download_instance.mod_to_interp_root, 
                                    ".temp",
                                    dir_tail
                                )
                
                initial_check_nc_files[local_dir] = {
                        "nc_files" : [f"{species}_{date}.nc" for date in date_list.strftime("%Y%m%d")],
                        "temp_dir" : temp_dir
                    }

        return initial_check_nc_files

    def extract_info_from_ncfile(self, nc_file):
        """
        Extract the species name and date from an BDRC NetCDF filename.

        Parameters
        ----------
        nc_file : str
            BDRC NetCDF filename containing the species and date information.

        Returns
        -------
        species : str
            Species name extracted from the filename.
        date : datetime.datetime
            Date extracted from the filename and converted to a datetime object.
        """

        # the species name is the first component of the filename.
        species = nc_file.split("_")[0]

        # extract the date from the filename
        date = datetime.strptime(nc_file.split("_")[-1][:-3], "%Y%m%d")

        return species, date

    def get_session(self):
        """
        Create and configure an authenticated HTTP session
        using the BDRC username and password.
        """

        USERNAME = self.bdrc_password["USERNAME"]
        PASSWORD = self.bdrc_password["PASSWORD"]

        self.session = requests.Session()
        self.session.auth = (USERNAME, PASSWORD)
        
    def download(self, temp_dir, date, mod_id):
        """
        Download the BDRC NetCDF file for the given date to a temporary directory.

        Parameters
        ----------
        temp_dir : str
            Directory where the downloaded BDRC NetCDF file will be stored.
        date : datetime.datetime
            Date and time used to determine the BDRC file to download.
        mod_id : str
            Model ID requested by the user.

        Returns
        -------
        temp_path : str
            Path to the downloaded BDRC NetCDF file.
        """

        # build the BDRC filename using the date and product name
        filename = f"{date.strftime('%Y%m%d')}{self.bdrc_model[mod_id]['ending']}_{self.bdrc_model[mod_id]['format']}"

        # get model ID with the capitalization of the url
        url_mod_id = self.bdrc_model[mod_id]["url"]

        # download source file
        temp_path = os.path.join(temp_dir, filename)

        # BDRC files are organised on the server by year and year/month
        url = f"{self.bdrc_fixed_values['url']}{url_mod_id}/{date.year}/{date.month:02d}/{filename}"

        # control not available dates
        try:
            response = self.session.get(
                        url,
                        verify=False,
                        stream=True,
                        timeout=300,
                    )
            response.raise_for_status()
        except requests.exceptions.HTTPError as exc:
            if response.status_code == 404:
                self.download_instance.logger.warning("No data available for day %s", date)
                return None

        # write the downloaded content to the temporary file in chunks
        with open(temp_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        return temp_path     

    def format_data(self, input_filepath, output_dir, nc_file, species):
        """
        Reformat a raw BDRC NetCDF file into a standardized
        Providentia-compatible NetCDF.

        Parameters
        ----------
        input_filepath : str
            Path where the input BDRC NetCDF file is located.
        output_dir : str
            Directory where the formatted NetCDF file will be written.
        nc_file : str
            Final name of the BDRC NetCDF file.
        species : str
            Providentia species name.
        """

        output_filepath = join(output_dir, nc_file)

        # get last downloaded file in case there was a keyboard interrupt
        self.download_instance.latest_nc_file_path = output_filepath  

        # open original netcdf file
        input_file = Dataset(input_filepath, "r", format="NETCDF4")

        # create new netcdf file
        output_file = Dataset(output_filepath, "w", format="NETCDF4")
        output_file.set_auto_mask(True)

        # create the dimensions
        for input_dim_name, output_dim_name in input_file.dimensions.items():
            output_file.createDimension(input_dim_name, output_dim_name.size)

        name_map = {
        "latitude": "latitude",
        "longitude": "longitude",
        "time": "time",
        }

        species_map = {
                        "od550du" : {"OD550_DUST": "od550du"},
                        "sconcdu" : {"SCONC_DUST": "sconcdu"}
                                    }

        unit_map = {
                    "od550du" : "unitless",
                    "sconcdu" : "ug m-3"
                                }

        name_map.update(species_map[species])

        # copy variables
        for input_var_name in input_file.variables:

            if input_var_name not in name_map:
                continue

            # get the output var name
            output_var_name = name_map[input_var_name]

            # get the variable
            input_var = input_file[input_var_name]

            # create the variable
            output_var = output_file.createVariable(
                output_var_name, input_var.datatype, input_var.dimensions
            )

            # add atributes for the species
            if output_var_name == "time":
                output_var.setncattr("calendar", "proleptic_gregorian")
                output_var.setncattr("units", input_var.units)

            elif output_var_name == species:
                output_file[species].setncattr("grid_mapping", "crs")
                output_file[species].setncattr("coordinates", "latitude longitude")
                output_file[species].setncattr("units", unit_map[species])

            # get the data from
            if output_var_name == "time":
                data = np.arange(len(input_var[:]))
                data *= 6
            else:
                data = input_var[:]

            # add the data to the variable
            output_var[:] = data

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

        # change the last downloaded file
        self.download_instance.latest_nc_file_path = "/path/to/file"
    
    def download_BDRC_model(
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
                f"\nDownloading {model} model data from the Barcelona Dust Regional Center..."
            )

            correct_mod_id = self.control_mod_id(mod_id)

            if not correct_mod_id:
                return

            correct_domain = self.control_domain(domain)

            if not correct_domain:
                return

            resolutions_list = self.control_resolutions()

            if not resolutions_list:
                return

            species_list = self.control_species()

            if not species_list:
                return
            
            start_date, end_date = self.control_dates(mod_id)

            if start_date is None and end_date is None:
                return

            initial_check_nc_files = self.build_nc_file_paths_in_range(start_date, end_date, mod_id, domain, resolutions_list, species_list)

            return initial_check_nc_files

        elif files_to_download:
            self.download_instance.logger.info(
                f"\n{model} model data to download ({len(files_to_download)}):"
            )

            for local_dir, files_to_download_dict in files_to_download.items():

                # create directories 
                os.makedirs(files_to_download_dict["temp_dir"], exist_ok=True)
                os.makedirs(local_dir, exist_ok=True)

                # obtain request session
                self.get_session()

                self.download_instance.logger.info(
                                f"\n  - {local_dir}"
                            )
        
                # iterate through each individual nc file
                for nc_file in tqdm(
                            files_to_download_dict['nc_files'],
                            bar_format="{l_bar}{bar}|{n_fmt}/{total_fmt}",
                            desc=f"    Downloading files ({len(files_to_download_dict['nc_files'])})",
                            ):

                    species, date = self.extract_info_from_ncfile(nc_file)

                    # peform de download and formatting only if the date is available      
                    temp_path = self.download(files_to_download_dict["temp_dir"], date, mod_id)

                    # continue with loop if day is not available
                    if not temp_path:
                        continue
            
                    self.format_data(temp_path, local_dir, nc_file, species)

        else:
            # tell the user if not valid resolution specie date combinations
            msg = "There is no available model output to be downloaded."
            show_message(self.download_instance, msg, deactivate=initial_check)
