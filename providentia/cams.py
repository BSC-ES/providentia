""" Class for downloading and formatting CAMS data """

from datetime import datetime, timedelta, timezone
import os
import re
import shutil
import sys
import zipfile

import cdsapi
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
import numpy as np
import requests
from tqdm import tqdm
import yaml

from providentia.auxiliar import CURRENT_PATH, join
from .warnings_prv import show_message

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

cams_options = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "internal", "cams", "cams_dataset.yaml"))
)
cams_variables_level = yaml.safe_load(
    open(
        join(
            PROVIDENTIA_ROOT,
            "settings",
            "internal",
            "cams",
            "cams_variables_level.yaml",
        )
    )
)
cams_formatting = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "internal", "cams", "cams_formatting.yaml"))
)
ghost_cams_variables = yaml.safe_load(
    open(
        join(
            PROVIDENTIA_ROOT,
            "settings",
            "internal",
            "cams",
            "ghost_cams_variables.yaml",
        )
    )
)
cams_stream = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "internal", "cams", "cams_stream.yaml"))
)
cams_species_units = yaml.safe_load(
    open(
        join(
            PROVIDENTIA_ROOT, "settings", "internal", "cams", "cams_species_units.yaml"
        )
    )
)
cdsapirc_urls = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "internal", "cams", "cdsapirc_urls.yaml"))
)


class Cams(object):
    """
    Class that manages the interaction with the Copernicus Atmosphere Data Store
    to retrieve CAMS and ERA5 datasets, validate requested models and dates, and
    convert downloaded NetCDF files into the Providentia-compatible format.
    """

    def __init__(self, download_instance):
        """
        Parameters
        ----------
        download_instance : providentia.Download
            Download controller instance.
        """

        self.download_instance = download_instance

    def get_model_info(self, config_modid):
        # count underscores to determine format
        u_count = config_modid.count("_")

        # get the prefix
        prefix = config_modid.rsplit("_", u_count - 1)[0]

        return u_count, prefix

    def control_domain(self, domain, prefix):
        # check if the domain is the correct one for the dataset
        if domain not in cams_options[prefix]:
            possible_domains = "', '".join(cams_options[prefix])
            msg = (
                f"The current domain '{domain}' is not valid for the CAMS dataset. "
                f"It must be '{possible_domains}'."
            )
            show_message(self.download_instance, msg)

            return False

        return True

    def control_species(self, cams_dict, dataset): 
        final_species_list = []
             
        # iterate through the species
        for species in self.download_instance.species:
            # if are interpolating between species, then only get model part.
            if '@' in species:
                species = species.split("@")[0]
            # check if species is in the ghost_cams_variables file
            if species not in ghost_cams_variables:
                msg = f"The species '{species}' is not available in CAMS."
                show_message(self.download_instance, msg)
                continue

            # get the CAMS species mapped to GHOST
            mapped_cams_species = ghost_cams_variables[species]

            # get the CAMS and GHOST list in [[cams_species, ghost_species]] format
            if type(mapped_cams_species) is list:
                cams_ghost_species_list = [
                    [ghost_cams_variables[ghost_species], ghost_species]
                    for ghost_species in mapped_cams_species
                ]
                cams_ghost_species_list.append([None, species])
            else:
                cams_ghost_species_list = [[mapped_cams_species, species]]

            for cams_species, ghost_species in cams_ghost_species_list:
                if ghost_species not in ["dir10", "spd10", "cldaf", "clddf", "photi"]:
                    # check if the mapped species are available in the dataset
                    if cams_species not in cams_dict["variable"]:
                        msg = f"Mapped species '{cams_species}' for input species '{ghost_species}' is not available in the CAMS '{dataset}' dataset."
                        show_message(
                            self.download_instance, msg
                        )
                        continue

                final_species_list.append([cams_species, ghost_species])

        return final_species_list

    def get_level(self, species_list, url):
        level_list = []

        for cams_species, ghost_species in species_list:
            # get the species' level
            level_list.append(
                "multi"
                if "multi" in cams_variables_level[url]
                and cams_species in cams_variables_level[url]["multi"]
                else "single"
            )

        return level_list

    def control_resolution(self, cams_dict, species_list, level_list):
        # get model resolution
        resolution_list = (
            self.download_instance.model_resolution
            if self.download_instance.model_resolution
            else self.download_instance.resolution
        )

        final_specie_resolution_dict = {}

        for (cams_species, ghost_species), level in zip(species_list, level_list):

            resolution_not_found = False

            # iterate through the resolutions
            for resolution in resolution_list:
                if not resolution_not_found:

                    # get the resolution for the cams dataset
                    correct_resolution = (
                        cams_dict["resolution"][level]
                        if type(cams_dict["resolution"]) is dict
                        else cams_dict["resolution"]
                    )

                    final_specie_resolution_dict[ghost_species] = []

                    # check if the resolution is the correct one for the dataset
                    if resolution == correct_resolution:
                        final_specie_resolution_dict[ghost_species].append(resolution)
                    else:
                        while True:
                            user_resolution = input(
                            f"\nNo non-interpolated model data is available for {resolution} resolution and {ghost_species} species. "
                            f"Available temporal resolution for this model: {correct_resolution}. "
                            f"Please specify the desired resolution or press Enter to automatically select "
                            f"'{correct_resolution}'. "
                            f"Otherwise enter 'n' if you do not want to download this model at another resolution. "
                            ).lower()

                            if (user_resolution == correct_resolution or user_resolution in ["", "n"]):
                                if user_resolution == "":
                                    final_specie_resolution_dict[ghost_species].append(correct_resolution)
                                elif user_resolution != "n":
                                    final_specie_resolution_dict[ghost_species].append(user_resolution)

                                # block input from the user happening again
                                resolution_not_found = True

                                break     

                    # remove in case there is no resolution for the species
                    if not final_specie_resolution_dict[ghost_species]:
                        del final_specie_resolution_dict[ghost_species]

        return final_specie_resolution_dict

    # DEPRECATED not sure if its necessary anymore, check if a bug appears
    def control_file_dates(self):
        # check if any of the files to download is in the current date range
        
        date_in_range = False

        for file in files_to_download:
            # get date from the file path
            date_str = (
                os.path.basename(file)
                .split("_")[1]
                .split(".")[0]
            )
            date = datetime.strptime(date_str, date_format)

            # normalize date to UTC
            date = (
                date.astimezone(timezone.utc)
                if date.tzinfo
                else date.replace(tzinfo=timezone.utc)
            )

            # break if any file is in the date range
            if cams_dict["forecast"]:
                in_range = (
                    current_cams_date <= date <= next_cams_date
                )
            else:
                in_range = (
                    date.year == current_cams_date.year
                    and date.month == current_cams_date.month
                )

            if in_range:
                date_in_range = True
                break

        # continue if there is no file in the date range
        if not date_in_range:
            current_cams_date = next_cams_date + timedelta(
                days=1
            )
            pass
            # continue
            
    def build_nc_file_paths_in_range(self, config_modid, domain, mod_id, stream, resolutions_list, species_list, level_list, cams_dict, cams_start_date, cams_end_date):
        # temporal directory for the zip file
        temp_root_dir = join(
            self.download_instance.mod_to_interp_root, ".temp"
        )

        initial_check_nc_files = {}

        for resolution in resolutions_list:
            for (cams_species, ghost_species), level in zip(species_list, level_list):

                # get directory structure
                dir_tail = join(config_modid, domain, resolution, ghost_species)

                # create specific final and temporal directory
                temp_dir = join(temp_root_dir, dir_tail)
                final_dir = join(
                    self.download_instance.mod_to_interp_root, dir_tail
                )

                # initialize iterators controlers
                current_cams_date = cams_start_date
                next_cams_date = (
                    cams_start_date.replace(day=1)
                    + relativedelta(months=1)
                    - timedelta(days=1)
                )

                #  if the lookahead days made the month change
                if cams_dict["lookahead_days"] > 0 and self.month_change:
                    next_cams_date = next_cams_date + relativedelta(months=1)

                # add the ncfiles to the files to download list
                initial_check_nc_files[final_dir] = {
                                "nc_files" : [],
                                "temp_dir" : temp_dir,
                                "requests" : {},
                                "ghost_species" : ghost_species,
                                "cams_species" : cams_species,
                                "resolution" : resolution,
                                "level" : level, 
                                "url" : cams_dict["url"],
                                "split" : cams_dict["split"],
                                                    }

                if cams_dict["forecast"] is True:
                    initial_check_nc_files[final_dir]["dated_file_format"] = cams_dict["dated_file_format"]
                
                if "dated_file_format" in cams_dict and len(all_dates) != 1:
                    initial_check_nc_files[final_dir]["zip_files"] = []

                while current_cams_date <= cams_end_date:
                    # add one month
                    next_cams_date = (
                        cams_end_date
                        if next_cams_date > cams_end_date
                        else next_cams_date
                    )

                    # get filename depending whether it is a download for the whole month or just a day
                    date_format = (
                        "%Y%m" if cams_dict["forecast"] is False else "%Y%m%d"
                    )

                    # if ghost_species not in [
                    #     "dir10",
                    #     "spd10",
                    #     "cldaf",
                    #     "clddf",
                    #     "photi",
                    # ]:
                    # create dictionary to do the request
                    request, is_valid = self.create_request(
                        cams_species,
                        cams_dict,
                        current_cams_date,
                        next_cams_date,
                        level,
                        stream,
                        mod_id,
                    )

                    # get request string to print later
                    request_str = self.get_request_str(cams_dict["dataset"], request)

                    # link the request to the files
                    initial_check_nc_files[final_dir]["requests"][request_str] = {"dict_request" : request,
                                                                                  "nc_files" : [],
                                                                                  "dataset" : cams_dict["dataset"]
                                                                                }

                    # jump to the next date
                    if not is_valid:
                        current_cams_date = next_cams_date + timedelta(days=1)
                        continue

                    # self.control_file_dates() ## ?? TODO

                    # iterate through the different days of the month if forecast
                    if "dated_file_format" in cams_dict:
                        all_dates = [
                            current_cams_date + timedelta(days=i)
                            for i in range(
                                (next_cams_date - current_cams_date).days + 1
                            )
                        ]
                    else:
                        all_dates = [current_cams_date]

                    initial_check_nc_files[final_dir]["requests"][request_str]["all_dates"] = all_dates

                    # iterate through all dates to format each of the day files
                    for date in all_dates:
                        if (
                            ghost_species
                            not in ["dir10", "spd10", "cldaf", "clddf", "photi"] ## ?
                        ):
                            # get the file format
                            if  "zip_files" in initial_check_nc_files[final_dir]:
                                zip_file_name = (
                                    cams_dict["dated_file_format"]
                                    .replace("yyyy", f"{date.year:04d}")
                                    .replace("mm", f"{date.month:02d}")
                                    .replace("dd", f"{date.day:02d}")
                                )

                                initial_check_nc_files[final_dir]["zip_files"].append(zip_file_name)

                        # format the file name
                        file_name = f"{ghost_species}_{date.strftime(date_format)}.nc"

                        initial_check_nc_files[final_dir]["nc_files"].append(file_name)
                        initial_check_nc_files[final_dir]["requests"][request_str]["nc_files"].append(file_name)

                        # add one day to the date
                        current_cams_date = next_cams_date + timedelta(days=1)
    
                        # prepare next cams date for the next iteration
                        next_cams_date = (
                            current_cams_date.replace(day=1)
                            + relativedelta(months=1)
                            - timedelta(days=1)
                        )

        return initial_check_nc_files

    def retrieve_request(self, temp_path, dataset, request, cdsapirc_path):
        client = cdsapi.Client(retry_max=1, quiet=True)

        try:
            client.retrieve(dataset, request, target=temp_path)
        except requests.exceptions.HTTPError as err:
            # invalid credential on .cdsapirc
            if err.response.status_code == 401:
                self.download_instance.logger.info(
                    "\nBad request (401): Client Error. Invalid credentials in the .cdsapirc file. "
                    "Removing authentication file...\n"
                    "Please run the program again so Providentia can recreate the file automatically, "
                    "or manually create a new .cdsapirc file by following the instructions at: "
                    "https://cds.climate.copernicus.eu/how-to-api"
                )
                os.remove(cdsapirc_path)
                return
            # bad request
            if err.response.status_code == 400:
                self.download_instance.logger.info(
                    "\nBad request (400): The server could not understand the request."
                )
                self.download_instance.logger.info(
                    f"Details: {err}"
                )
            # connection error
            elif err.response.status_code == 500:
                self.download_instance.logger.info(
                    "\nServer error (500): The server encountered an error while processing the request."
                )
                self.download_instance.logger.info(
                    f"Details: {err}"
                )
                self.download_instance.logger.info(
                    "Please try again later."
                )
                return
            else:
                self.download_instance.logger.info(
                    f"\nUnexpected error ({err.response.status_code}):"
                )
                self.download_instance.logger.info(
                    f"Details: {err}"
                )

    def extract_zip(self, temp_dir, temp_path):
        with zipfile.ZipFile(temp_path, "r") as zip_ref:
            zip_file_name = zip_ref.namelist()[0]
            zip_ref.extractall(temp_dir)

    def get_request_str(self, dataset, request):
        """
        Get formatted string representation of a CAMS request dictionary.

        Parameters
        ----------
        dataset : str
            Name of the dataset for which the request is made.
        request : dict
            Dictionary containing request parameters.
        Returns
        -------
        request_str : str
            String that represents a CAMS request dictionary.
        """

        request_str = ""

        request_str += f"dataset = '{dataset}'\n"
        request_str += "request = {\n"
        for k, v in request.items():
            if type(v) is str:
                v = f"'{v}'"
            request_str += f"'{k}' : {v},\n"
        request_str += "}\n"

        return request_str

    def fetch_cams_dates(self, url, cams_dict):
        """
        Extract the minimum and maximum available dates for a CAMS dataset by webscrapping.

        Parameters
        ----------
        url : str
            URL of the CAMS dataset.
        cams_dict : dict
            CAMS dataset configuration dictionary.

        Returns
        -------
        minstart : datetime.datetime
            Minimum available date of the dataset.
        maxend : datetime.datetime
            Maximum available date of the dataset.
        """

        # send HTTP GET request and get the text
        try:
            response = requests.get(url, timeout=20)
            response.raise_for_status()
            r = response.text
        except requests.exceptions.ReadTimeout:
            minstart, maxend = datetime.strptime("19990101", "%Y%m%d"), datetime.now()
            msg = f"Request timed out when accessing {url}. Using default period: {minstart.strftime('%Y-%m-%d')} and {maxend.strftime('%Y-%m-%d')} as minimum and maximum dates."
            show_message(self.download_instance, msg, print=True)
            return minstart, maxend

        # do the webscrapping depending if there is whole dates or only month
        if cams_dict["months_list"] is False:
            # get the minstart and maxend dictionary
            minstart_dict = re.findall(r'"minStart":".*?"', r, re.DOTALL)
            maxend_dict = re.findall(r'"maxEnd":".*?"', r, re.DOTALL)

            # get the value from the dictionary
            minstart = datetime.strptime(minstart_dict[0].split('"')[-2], "%Y-%m-%d")
            maxend = datetime.strptime(maxend_dict[0].split('"')[-2], "%Y-%m-%d")
        else:
            # get the interval dictionary
            match = re.search(r'"interval":\[\["(.*?)","(.*?)"\]\]', r)

            # get the date value
            minstart = datetime.strptime(match.group(1), "%Y-%m-%dT%H:%M:%S%z")
            maxend = datetime.strptime(match.group(2), "%Y-%m-%dT%H:%M:%S%z")

        return minstart, maxend

    def control_dates(self, url, cams_dict, initial_check=False):
        """
        Check that the configured start and end dates are valid,
        determines the minimum and maximum available dates from the CAMS service,
        and adjusts the requested period accordingly.

        Parameters
        ----------
        url : str
            CAMS dataset metadata URL used to query date availability.
        cams_dict : dict
            Dictionary describing the CAMS dataset configuration.
        initial_check : bool, optional
            If True, do not show warnings.

        Returns
        -------
        cams_start_date : datetime.datetime or None
            Start date of the valid CAMS download period.
        cams_end_date : datetime.datetime or None
            End date of the valid CAMS download period.
        """

        # get minimum and maximum possible dates
        min_start_date, max_end_date = self.fetch_cams_dates(url, cams_dict)

        # convert the selected dates to datetetime
        cams_start_date = datetime.strptime(self.download_instance.start_date, "%Y%m%d")
        cams_end_date = datetime.strptime(
            self.download_instance.end_date, "%Y%m%d"
        ) - timedelta(days=1)

        # warn the user that download is going to be for N days before
        if cams_dict["lookahead_days"] > 0:
            # download N days ahead for forecast
            lookahead_cams_start_date = cams_start_date - timedelta(days=cams_dict["lookahead_days"])

            # bool to control if the lookahead days change make the start date on the previous month
            self.month_change = lookahead_cams_start_date.month != cams_start_date.month

            cams_start_date = lookahead_cams_start_date

            msg = f"Model data will be downloaded {cams_dict['lookahead_days']} day(s) in advance relative to the configured date."
            show_message(self.download_instance, msg, deactivate=initial_check)

        # normalize all to UTC
        min_start_date = (
            min_start_date.astimezone(timezone.utc)
            if min_start_date.tzinfo
            else min_start_date.replace(tzinfo=timezone.utc)
        )
        max_end_date = (
            max_end_date.astimezone(timezone.utc)
            if max_end_date.tzinfo
            else max_end_date.replace(tzinfo=timezone.utc)
        )
        cams_start_date = (
            cams_start_date.astimezone(timezone.utc)
            if cams_start_date.tzinfo
            else cams_start_date.replace(tzinfo=timezone.utc)
        )
        cams_end_date = (
            cams_end_date.astimezone(timezone.utc)
            if cams_end_date.tzinfo
            else cams_end_date.replace(tzinfo=timezone.utc)
        )

        # if the minimum date is over the end date
        if min_start_date > cams_end_date or max_end_date < cams_start_date:
            msg = f"The selected dates are unavailable. Please choose dates between {min_start_date.strftime('%Y-%m-%d')} and {max_end_date.strftime('%Y-%m-%d')}."
            show_message(self.download_instance, msg, deactivate=initial_check)
            return None, None

        # check if the start date is within limits
        if min_start_date > cams_start_date:
            cams_start_date = min_start_date

        # check if the end date is within limits
        if max_end_date < cams_end_date:
            cams_end_date = max_end_date

        return cams_start_date, cams_end_date

    def manage_cdsapirc(self, model):
        # create cdsapirc file in case it was not created
        cdsapirc_path = join(os.getenv("HOME"), ".cdsapirc")

        # get url for .cdsapirc
        for model_type in cdsapirc_urls:
            if model_type in model:
                break

        self.cdsapirc_url = cdsapirc_urls[model_type]

        # create csapirc file necessary for the download
        if not os.path.isfile(cdsapirc_path):
            self.create_cdsapirc(cdsapirc_path)
        else:
            self.change_cdsapirc(cdsapirc_path)

        return cdsapirc_path

    def create_request(
        self,
        cams_species,
        cams_dict,
        current_cams_date,
        next_cams_date,
        level,
        stream,
        mod_id,
    ):
        """
        Build the request required by the Copernicus Atmosphere Data Store API.

        Parameters
        ----------
        cams_species : str
            CAMS variable name to request.
        cams_dict : dict
            CAMS dataset configuration dictionary.
        current_cams_date : datetime.datetime
            Start date of the request period.
        next_cams_date : datetime.datetime
            End date of the request period.
        level : str
            Variable level type, 'single' or 'multi'.
        stream : str or None
            CAMS data stream, 'validated' or 'interim'.
        mod_id : str or None
            Model identifier for multi-model datasets.

        Returns
        -------
        dict : request
            Request dictionary.
        is_valid : bool
            Indicator of a valid request.
        """

        # initialise request and boolean
        request = {"variable": cams_species}
        is_valid = True

        # add area the request if the dataset has it and user asked for it
        if cams_dict["area"] and hasattr(self.download_instance, "longitude") and hasattr(self.download_instance, "latitude"):
            request["area"] = [self.download_instance.latitude[1], self.download_instance.longitude[0], 
             self.download_instance.latitude[0], self.download_instance.longitude[1]]
            
        # add leadtime hour to the request if the dataset has it
        if "leadtime_hour" in cams_dict and level in cams_dict["leadtime_hour"]:
            request["leadtime_hour"] = cams_dict["leadtime_hour"][level]

        # pass numerical date if the request allows it, pass month as a list
        if cams_dict["months_list"] is False:
            request[
                "date"
            ] = f"{current_cams_date.strftime('%Y-%m-%d')}/{next_cams_date.strftime('%Y-%m-%d')}"
        else:
            request["year"] = str(current_cams_date.year)
            request["month"] = str(current_cams_date.strftime("%m"))

        # pass days as a list if the request needs it
        if cams_dict["days_list"] is True:
            request["day"] = [
                f"{i:02d}"
                for i in range(1, (next_cams_date - current_cams_date).days + 2)
            ]

        # add interim and or validated stream to the request
        if cams_dict["stream"] is True:
            # check whether the stream is available for the year
            if stream not in cams_stream[request["year"]]:
                msg = f"The current stream '{stream}' is not available for the current date: {request['month']}-{request['year']}. Continuing..."
                show_message(self.download_instance, msg)
                is_valid = False
                return request, is_valid

            request["type"] = stream

        # add the model if models are available in the dataset
        if "model" in cams_dict:
            request["model"] = mod_id

        # get the level and apply it if the species is multi level
        level_variable = "level" if "level" in cams_dict else "model_level"
        if level == "multi":
            # choose the maximum level if there was a level increase at some point
            if "level_boundary" in cams_dict:
                boundary_date = datetime.combine(
                    cams_dict["level_boundary"], datetime.min.time()
                )
                boundary_date = (
                    boundary_date.astimezone(timezone.utc)
                    if boundary_date.tzinfo
                    else boundary_date.replace(tzinfo=timezone.utc)
                )

                if current_cams_date <= boundary_date:
                    request[level_variable] = cams_dict[level_variable][
                        "before_increase"
                    ]
                else:
                    request[level_variable] = cams_dict[level_variable][
                        "after_increase"
                    ]

            # without a level increase, just get the level
            else:
                request[level_variable] = cams_dict[level_variable]

        # copy shared fields from config into the request
        shared_variables = [
            "time",
            "type",
            "data_format",
            "product_type",
            "download_format",
        ]

        for shared_variable in shared_variables:
            if shared_variable in cams_dict:
                request[shared_variable] = cams_dict[shared_variable]

        return request, is_valid

    def create_cdsapirc(self, cdsapirc_path):
        """
        Create the `.cdsapirc` configuration file required for the CAMS API.

        Parameters
        ----------
        cdsapirc_path : str
            Absolute path where the `.cdsapirc` file will be created.
        """

        # ask the user whether they want to create the file in the home directory
        while True:
            create_file = input(
                f"\n'.cdsapirc' file not found. Creating it at {cdsapirc_path}. Do you agree? ([y]/n) "
            ).lower()
            if create_file in ["", "y", "n"]:
                break
        # create file if user agreed with it
        if create_file in ["", "y"]:
            # ask the user for the personal access token
            personal_access_token = input(
                "\nEnter your personal access token, which you can find at https://ads.atmosphere.copernicus.eu/profile after login: "
            )
            # create the .cdsapirc file with the user's acces token
            with open(cdsapirc_path, "w") as f:
                f.write(f"url: {self.cdsapirc_url}\n")
                f.write(f"key: {personal_access_token}\n")
        else:
            self.download_instance.logger.error(
                "Error: Cannot proceed without '.cdsapirc'. CAMS model data download requires this file."
            )
            sys.exit(1)

    def change_cdsapirc(self, cdsapirc_path):
        """
        Change the `.cdsapirc` configuration file required for the CAMS API.

        Parameters
        ----------
        cdsapirc_path : str
            Absolute path where the `.cdsapirc` file will be created.
        """

        # get .cdsapirc file contents
        with open(cdsapirc_path, "r") as f:
            data = f.readlines()

        # change url if necessary
        for i, line in enumerate(data):
            if line.startswith("url"):
                if self.cdsapirc_url not in line:
                    data[i] = f"url: {self.cdsapirc_url}\n"
                    with open(cdsapirc_path, "w") as f:
                        f.truncate(0)
                        f.writelines(data)
                break

    def get_model(
        self,
        cams_dict,
        u_count,
        config_modid,
        dataset,
        ensemble_option,
        initial_check=False,
    ):
        """
        Determine model ID and stream from CAMS configuration.

        Parameters
        ----------
        cams_dict : dict
            CAMS dataset configuration dictionary.
        u_count : int
            Number of underscores in the model identifier defined in the configuration file.
        config_modid : str
            Exact model identifier as specified in the configuration file.
        dataset : str
            Name of the dataset to which the model belongs.
        ensemble_option : str
            Ensemble member specified in the configuration file.
        initial_check : bool, optional
            If True, suppress warnings related to unavailable streams.

        Returns
        -------
        mod_id : str or None
            Model identifier extracted from the configuration, None if no model applies.
        stream : str or None
            Stream extracted from the configuration, None if no stream applies.
        error : bool
            True if the configuration is invalid, False if it is valid.
        """

        # determine model ID or stream
        mod_id, stream, error = None, None, True

        if u_count == 1:
            # e.g. cams_forecast
            if "model" in cams_dict:
                msg = f"The model '{config_modid}' is missing the model. Please add one (e.g. '{config_modid}_ensemble')."
                show_message(self.download_instance, msg, deactivate=initial_check)
                return mod_id, stream, error

        elif u_count == 2:
            # e.g. cams_analysis_ensemble
            last_element = config_modid.rsplit("_", 1)[1]

            if cams_dict.get("stream") and "model" in cams_dict:
                msg = f"The '{dataset}' dataset needs a model and a stream. Please add one (e.g. 'cams_reanalysis_ensemble_interim')."
                show_message(self.download_instance, msg, deactivate=initial_check)
                return mod_id, stream, error

            if "model" in cams_dict:
                if last_element not in cams_dict["model"]:
                    msg = f"Cannot find the {last_element} model in the '{dataset}' dataset."
                    show_message(self.download_instance, msg, deactivate=initial_check)
                    return mod_id, stream, error
                mod_id = last_element
            else:
                msg = f"The '{dataset}' dataset does not admit models or streams."
                show_message(self.download_instance, msg, deactivate=initial_check)
                return mod_id, stream, error

        elif u_count == 3:
            # e.g. cams_reanalysis_ensemble_interim
            if not (cams_dict.get("stream") and "model" in cams_dict):
                msg = f"The '{dataset}' dataset does not admit models and streams, change the model in the configuration file."
                show_message(self.download_instance, msg, deactivate=initial_check)
                return mod_id, stream, error

            # extract model and stream
            _, temp_mod_id, temp_stream = config_modid.rsplit("_", 2)

            if temp_mod_id not in cams_dict["model"]:
                msg = f"Cannot find the {temp_mod_id} model in the '{dataset}' dataset."
                show_message(self.download_instance, msg, deactivate=initial_check)
                return mod_id, stream, error

            if temp_stream not in ["validated", "interim"]:
                msg = f"'{temp_stream}' is not a valid stream. Availabe streams: validated, interim."
                show_message(self.download_instance, msg, deactivate=initial_check)
                return mod_id, stream, error

            mod_id, stream = temp_mod_id, temp_stream

            # add reanalysis suffix
            stream += "_reanalysis"

        else:
            msg = f"The '{config_modid}' format is not valid."
            show_message(self.download_instance, msg, deactivate=initial_check)
            return mod_id, stream, error

        # only ensemble options allmembers and 000 are valid
        if ensemble_option not in ["000", "allmembers"]:
            msg = (
                f"The current ensemble option '{ensemble_option}' is not valid for the CAMS '{dataset}' dataset. "
                f"It must be '000' or 'allmembers'."
            )
            show_message(self.download_instance, msg, deactivate=initial_check)
            return mod_id, stream, error

        # successful parsing
        error = False
        return mod_id, stream, error

    def extract_date(self, input_file, prefix, domain):
        """
        Extract the date from a CAMS NetCDF file depending on model type and domain.

        Parameters
        ----------
        input_file : netCDF4.Dataset or similar
            NetCDF dataset object representing the original CAMS file.
        prefix : str
            Model type prefix, e.g. 'cams_analysis', 'cams_forecast' or 'cams_reanalysis'.
        domain : str
            Domain of the model, either 'regional' or 'global'.

        Returns
        -------
        time_str : str
            Extracted date as a string in 'YYYY-MM-DD' format.
        """

        if prefix in ["cams_analysis", "cams_forecast"] and domain == "regional":
            time = input_file["time"].long_name.split()[-1]
            time = datetime.strptime(time, "%Y%m%d")
        elif prefix == "cams_reanalysis" and domain == "regional":
            time = input_file["time"].units.split()[2]
            time = datetime.strptime(time, "%Y-%m-%d")
        elif prefix == "cams_forecast" and domain == "global":
            time = input_file["forecast_reference_time"][0]
            time = datetime.fromtimestamp(int(time))
        elif prefix in ["cams_reanalysis", "era5_reanalysis"] and domain == "global":
            time = input_file["valid_time"][0]
            time = datetime.fromtimestamp(int(time))

        time_str = time.strftime("%Y-%m-%d")

        return time_str

    def format_data(
        self,
        input_filepath,
        output_filepath,
        species,
        prefix,
        domain,
        resolution,
        cams_species,
        url,
    ):
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

        self.download_instance.logger.info(f"Formatting {output_filepath}")

        # get file formatting
        cams_providentia_map = cams_formatting[prefix][domain]

        # open original netcdf file
        input_file = Dataset(input_filepath, "r", format="NETCDF4")

        # extract date
        date_str = self.extract_date(input_file, prefix, domain)

        # create new netcdf file
        output_file = Dataset(output_filepath, "w", format="NETCDF4")
        output_file.set_auto_mask(True)

        # change the last downloaded file
        self.download_instance.latest_nc_file_path = output_filepath

        for input_dim_name, output_dim_name in cams_providentia_map.items():
            # skip single species
            if output_dim_name == "level" and (
                "single" in cams_variables_level[url]
                and cams_species in cams_variables_level[url]["single"]
            ):
                continue
            # skip dimensions not used
            if not output_dim_name:
                continue
            # get dimension
            dim = input_file.dimensions[input_dim_name]
            # create the dimension with the new name
            output_file.createDimension(output_dim_name, len(dim))

        # copy variables
        for input_var_name in input_file.variables:
            # get the output var name
            if input_var_name in cams_providentia_map:
                output_var_name = cams_providentia_map[input_var_name]
                # skip variables not used
                if not output_var_name:
                    continue
            else:
                output_var_name = species

            # get the variable
            input_var = input_file[input_var_name]

            # change the name of the dimensions into the providentia name
            output_var_dims = [
                cams_providentia_map.get(name, name)
                for name in input_var.dimensions
                if name in cams_providentia_map and cams_providentia_map[name]
            ]
            # create the variable
            output_var = output_file.createVariable(
                output_var_name, input_var.datatype, output_var_dims
            )

            # add calendar and units attributes to the time variable
            if output_var_name == "time":
                output_var.setncattr(
                    "calendar", cams_options[prefix][domain]["calendar"]
                )
                output_var.setncattr("units", f"hours since {date_str}")

            # add to level the priority of the level
            elif output_var_name == "level":
                output_var.positive = "up"

            # add coordinates, grid_mapping and units to the species variable
            elif output_var_name == species:
                output_var.setncattr("coordinates", "latitude longitude")
                output_var.setncattr("grid_mapping", "crs")
                output_var.setncattr("units", cams_species_units[input_var.units])

            # get the data from
            if output_var_name == "time":
                data = np.arange(len(input_var[:]))
                if resolution == "3hourly":
                    data *= 3
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

    def calculate_dir10(self, u10_data, v10_data):
        """
        Calculate the wind direction from the u and v
        components of wind.

        Parameters
        ----------
        u10_data : np.array
            Array that stores the u component of wind.
        v10_data : np.array
            Array that stores the v component of wind.
        """
        # create empty array to store the operation result
        dir10_data = np.empty_like(u10_data, dtype="float32")

        # calculate direction in chunks to not run out of memory
        for i in range(u10_data.shape[0]):
            u = u10_data[i]
            v = v10_data[i]
            dir10_data[i] = np.sqrt(u * u + v * v)

        return dir10_data

    def calculate_spd10(self, u10_data, v10_data):
        """
        Calculate the wind speed from the u and v
        components of wind.

        Parameters
        ----------
        u10_data : np.array
            Array that stores the u component of wind.
        v10_data : np.array
            Array that stores the v component of wind.
        """
        # create empty array to store the operation result
        dir10_data = np.empty_like(u10_data, dtype="float32")

        # calculate speed in chunks to not run out of memory
        for i in range(u10_data.shape[0]):
            u = u10_data[i]
            v = v10_data[i]
            degrees = np.degrees(np.arctan2(-u, -v))

            # transform negative degrees to positive
            degrees[degrees < 0] += 360
            dir10_data[i] = degrees

        return dir10_data

    def calculate_radiation(self, radiation1_data, radiation2_data):
        """
        Calculate the radiation from two different radiations.

        Parameters
        ----------
        radiation1_data : np.array
            Array that stores the dividend of the radiation.
        radiation2_data : np.array
            Array that stores the divisor of the radiation.
        """
        # create empty array to store the operation result
        radiation_result_data = np.empty_like(radiation1_data, dtype="float32")

        # calculate speed in chunks to not run out of memory
        for i in range(radiation1_data.shape[0]):
            radiation1 = radiation1_data[i]
            radiation2 = radiation2_data[i]
            radiation_result_data[i] = radiation1 / radiation2

        return radiation_result_data

    def format_product(self, output_filepath, species, ghost_cams_variables):
        """
        Combine u and v components of wind to create
        a standardized Providentia-compatible NetCDF
        for wind direction and speed.

        Parameters
        ----------
        output_filepath : str
            Path where the formatted NetCDF file will be written.
        species : str
            Providentia species name.
        """

        self.download_instance.logger.info(f"Formatting {output_filepath}")

        var1, var2 = ghost_cams_variables[species]

        # open var1 netcdf file
        input_filepath_var1 = output_filepath.replace(species, var1)
        input_file_var1 = Dataset(input_filepath_var1, "r", format="NETCDF4")

        # open var2 netcdf file
        input_filepath_var2 = output_filepath.replace(species, var2)
        input_file_var2 = Dataset(input_filepath_var2, "r", format="NETCDF4")

        # get the speed/direction/radiation calculation function
        if species == "dir10":
            calculate_product = self.calculate_dir10
        elif species == "spd10":
            calculate_product = self.calculate_spd10
        else:
            calculate_product = self.calculate_radiation

        # create new netcdf file
        output_file = Dataset(output_filepath, "w", format="NETCDF4")
        output_file.set_auto_mask(True)

        # change the last downloaded file
        self.download_instance.latest_nc_file_path = output_filepath

        # create the dimensions
        for input_dim_name, output_dim_name in input_file_var1.dimensions.items():
            output_file.createDimension(input_dim_name, output_dim_name.size)

        # copy variables
        for input_var_name in input_file_var1.variables:
            # get the variable
            input_var_var1 = input_file_var1[input_var_name]

            # get the speed/direction from the u and v components of wind
            if input_var_name == var1:
                input_var_var2 = input_file_var2[var2]
                data = calculate_product(input_var_var1[:], input_var_var2[:])
                output_var_name = species
            # copy the data of u for the rest of the variables
            else:
                data = input_var_var1[:]
                output_var_name = input_var_name

            # create the variable
            output_var = output_file.createVariable(
                output_var_name, input_var_var1.datatype, input_var_var1.dimensions
            )

            # add attributes to the variable
            output_var.setncatts(input_var_var1.__dict__)

            # change the direction units to degrees
            if input_var_name == var1:
                if species == "dir10":
                    output_var.setncattr("units", "degrees")
                elif species in ["cldaf", "clddf", "photi"]:
                    output_var.setncattr("units", "unitless")

            # add the data to the variable
            output_var[:] = data

        # close the original and new netcdf files
        output_file.close()
        input_file_var1.close()
        input_file_var2.close()

    def split_nc_file(
        self, input_file_name, all_dates, dated_file_format, temp_dir, prefix, domain, level
    ):
        """
        Split a multi-day CAMS NetCDF forecast file into daily files.

        Parameters
        ----------
        input_file_name : str
            Name of the input NetCDF file to split.
        all_dates : list of datetime.datetime
            List of dates corresponding to forecast slices.
        dated_file_format : str
            String with the NCfile date format.
        temp_dir : str
            Temporary directory containing the input file and output slices.
        prefix : str
            CAMS dataset prefix.
        domain : str
            Spatial domain,'global' or 'regional'.
        level : str
            Variable level type, 'single' or 'multi'.
        """

        # get file formatting
        cams_providentia_map = cams_formatting[prefix][domain]

        # read the input netcdf file
        input_filepath = join(temp_dir, input_file_name)
        input_file = Dataset(input_filepath, "r")

        # set available dimensions
        available_dimensions = ["forecast_period", "latitude", "longitude"]
        if level == "multi":
            available_dimensions.append("model_level")

        # create tqdm iterator
        all_dates_iter = tqdm(
            all_dates,
            bar_format="{l_bar}{bar}|{n_fmt}/{total_fmt}",
            desc=f"Splitting {input_file_name} file ({len(all_dates)})",
        )

        # loop through the possible dates
        for i, date in enumerate(all_dates_iter):
            # create a new file for each slice
            output_file_name = (
                dated_file_format
                .replace("yyyy", f"{date.year:04d}")
                .replace("mm", f"{date.month:02d}")
                .replace("dd", f"{date.day:02d}")
            )
            output_filepath = join(temp_dir, output_file_name)
            output_file = Dataset(output_filepath, "w", format="NETCDF4")

            # copy all the dimensions to the new file, leave forecas_reference_time as one
            for dim in available_dimensions:
                output_file.createDimension(dim, input_file.dimensions[dim].size)
            output_file.createDimension("forecast_reference_time", 1)

            # create the variable in the new dataset and assign it the slice
            output_var = output_file.createVariable(
                "forecast_reference_time",
                input_file["forecast_reference_time"].datatype,
                input_file["forecast_reference_time"].dimensions,
            )
            output_var[:] = input_file["forecast_reference_time"][i]

            # copy all other variables from the original dataset
            for input_var_name in input_file.variables:
                # skip forecast_reference_time
                if input_var_name == "forecast_reference_time":
                    continue
                input_var = input_file[input_var_name]
                # create the same variable in the new file
                output_var = output_file.createVariable(
                    input_var_name, input_var.datatype, input_var.dimensions
                )

                # copy attributes
                output_var.setncatts(
                    {attr: input_var.getncattr(attr) for attr in input_var.ncattrs()}
                )

                # copy data, get the specific dimensions depending on the variable
                if input_var_name == "valid_time":
                    output_var[:] = input_var[i, :]
                elif input_var_name not in cams_providentia_map:
                    if level == "multi":
                        output_var[:] = input_var[:, i, :, :, :]
                    else:
                        output_var[:] = input_var[:, i, :, :]
                else:
                    output_var[:] = input_var[:]

            # close new dataset
            output_file.close()

        self.download_instance.logger.info("")

        # close original dataset
        input_file.close()

    def download_cams_model(self, model, initial_check, files_to_download=None):
        """
        Download and process CAMS model data. Validates the
        requested model configuration, checks date availability, builds ADS
        requests, downloads CAMS NetCDF files and reformats them into
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

        if initial_check:
            # print current model
            self.download_instance.logger.info("\n" + "-" * 40)
            self.download_instance.logger.info(
                f"\nDownloading {model} model data from the Atmosphere/Climate Data Store..."
            )

            # get model id and the domain
            config_modid, domain, ensemble_options = model.split("-")

            u_count, prefix = self.get_model_info(config_modid)

            correct_domain = self.control_domain(domain, prefix)

            if not correct_domain:
                return

            # get the dictionary with the dataset characteristics
            cams_dict = cams_options[prefix][domain]

            # make the necessary checks to the model
            mod_id, stream, invalid_model = self.get_model(
                cams_dict, u_count, config_modid, cams_dict["dataset"], 
                ensemble_options, initial_check
            )

            if invalid_model:
                return

            # make the necessary checks to the dates
            cams_start_date, cams_end_date = self.control_dates(
                cams_dict["url"], cams_dict, initial_check
            )

            if cams_start_date is None and cams_end_date is None:
                return

            species_list = self.control_species(cams_dict, cams_dict["dataset"])

            if not species_list:
                return

            # get level depending on the species
            level_list = self.get_level(species_list, cams_dict["url"])

            resolutions_list = self.control_resolution(cams_dict, species_list, level_list)

            if not resolutions_list:
                return

            initial_check_nc_files = self.build_nc_file_paths_in_range(config_modid, domain, mod_id, stream, resolutions_list, species_list, level_list, cams_dict, cams_start_date, cams_end_date)

            return initial_check_nc_files

        elif files_to_download:

                self.download_instance.logger.info(
                                f"\n{model} model data to download ({len(files_to_download)}):"
                            )

                # create / change the cdsapirc
                cdsapirc_path = self.manage_cdsapirc(model)

                for local_dir, files_to_download_dict in files_to_download.items():
    
                    # create directories
                    os.makedirs(files_to_download_dict["temp_dir"], exist_ok=True)
                    os.makedirs(local_dir, exist_ok=True)

                    self.download_instance.logger.info("")

                    for request_str, request_dict in files_to_download_dict["requests"].items():

                        # get files the files that passed the overwrite filter
                        nc_files = list(sorted(set(files_to_download_dict["nc_files"]) & set(request_dict["nc_files"])))

                        if nc_files:

                            # print the request
                            self.download_instance.logger.info(request_str)

                            # get zip file path
                            temp_path = join(files_to_download_dict["temp_dir"], "zip_file")

                            # get zip file from ADS / CDS
                            # self.retrieve_request(temp_path, request_dict["dataset"], request_dict["dict_request"], cdsapirc_path)
                            
                            # extract file
                            zip_file_name = self.extract_zip(files_to_download_dict["temp_dir"], temp_path)

                            if files_to_download_dict["ghost_species"] not in [
                                                    "dir10",
                                                    "spd10",
                                                    "cldaf",
                                                    "clddf",
                                                    "photi",
                                                ]:
                                # split the forecast file
                                if request_dict["split"] is True:
                                    self.split_nc_file(
                                        zip_file_name,
                                        request_dict["all_dates"],
                                        request_dict["dated_file_format"],
                                        files_to_download_dict["temp_dir"],
                                        prefix,
                                        domain,
                                        files_to_download_dict["level"],
                                    )

                            for nc_file in nc_files:

                                # obtain middle and final path
                                final_path = join(local_dir, nc_file)
                                input_filepath = join(files_to_download_dict["temp_dir"], nc_file)

                                # format the ncfiles
                                if files_to_download_dict["ghost_species"] in [
                                    "dir10",
                                    "spd10",
                                    "cldaf",
                                    "clddf",
                                    "photi",
                                ]:
                                    self.format_product(
                                        final_path, files_to_download_dict["ghost_species"], ghost_cams_variables
                                    )
                                else:
                                    self.format_data(
                                        input_filepath,
                                        final_path,
                                        files_to_download_dict["ghost_species"],
                                        prefix,
                                        domain,
                                        files_to_download_dict["resolution"],
                                        files_to_download_dict["cams_species"],
                                        files_to_download_dict["url"],
                                    )

                                # change the last downloaded file
                                self.download_instance.latest_nc_file_path = "/path/to/file"

