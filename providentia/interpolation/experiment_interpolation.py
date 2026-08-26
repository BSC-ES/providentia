""" Class for the spatial and temporal interpolation of model data to observational data"""

import ast
from calendar import monthrange
import copy
import datetime
import glob
import os
import pathlib
import shutil
import subprocess
import sys
import time
import traceback

import cartopy.crs as ccrs
import cftime
import dateutil.relativedelta as relativedelta
from netCDF4 import Dataset, num2date, chartostring
import numpy as np
from packaging.version import Version
import pandas as pd
import pyproj
import scipy
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, Point
import xarray as xr
import yaml

from aux_interp import (
    check_for_ghost,
    findMiddle,
    check_directory_existence,
    set_file_permissions_ownership,
    get_aeronet_bin_radius_from_bin_variable,
    get_aeronet_model_bin,
    get_model_to_aeronet_bin_transform_factor,
)

sys.path.append(str(pathlib.Path(__file__).resolve().parents[2]))
from providentia.auxiliar import (
    CURRENT_PATH,
    join,
    get_machine,
    get_conversion_factor,
    get_standard_parameters_by_speci,
)

# set current MACHINE
MACHINE = get_machine()

# get current path and providentia root path
PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

# load the defined models paths and agrupations jsons
interp_models = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "interp_models.yaml"))
)


class ModelInterpolation(object):
    """Class which handles interpolation of model data to surface data, both spatially and temporally."""

    def __init__(self, submit_args):
        """
        Sets up the interpolation environment by reading submission arguments,
        retrieving model and observational files, handling ensemble or forecast
        data, applying GHOST standards, and computing forecast information.

        Parameters
        ----------
        submit_args : dict
            Dictionary of job submission parameters, including:
            - model_temporal_resolution : str
            - mod_speci_to_process : str
            - network_to_interpolate_against : str
            - temporal_resolution_to_output : str
            - yearmonth : str
            - speci_to_process : str
            - job_id : str
            - prov_mod_code : str
        """

        self.log_file_str = "STARTING INTERPOLATION\n"

        self.log_file_str += (
            "\n".join(f"{k}: {v}" for k, v in submit_args.items()) + "\n"
        )

        # set variables from input keywords
        self.model_temporal_resolution = submit_args["model_temporal_resolution"]
        self.mod_speci_to_process = submit_args["mod_speci_to_process"]
        self.network_to_interpolate_against = submit_args[
            "network_to_interpolate_against"
        ]
        self.temporal_resolution_to_output = submit_args[
            "temporal_resolution_to_output"
        ]
        self.yearmonth = submit_args["yearmonth"]
        self.speci_to_process = submit_args["speci_to_process"]
        self.unique_id = submit_args["job_id"]
        self.prov_mod_code = submit_args["prov_mod_code"]
        (
            self.experiment_to_process,
            self.grid_type,
            self.ensemble,
        ) = self.prov_mod_code.split("-")

        # determine if are interpolating to another species (e.g. sconco3@t2) or not
        if '@' in self.speci_to_process:
            speci_to_process_split = self.speci_to_process.split('@')
            self.original_mod_speci_to_process = speci_to_process_split[0]
            self.obs_speci_to_process = speci_to_process_split[1]
        else:
            self.original_mod_speci_to_process = copy.deepcopy(self.speci_to_process)
            self.obs_speci_to_process = copy.deepcopy(self.speci_to_process)

        # get year/month string
        self.year = self.yearmonth[:4]
        self.month = self.yearmonth[4:]

        # determine if ensemble is member or emsemble stat
        self.ensemble_member = self.ensemble.isdigit()

        # dictionary to save utilized interpolation variables
        self.interpolation_variables = {}

        # from file in management_logs, get and set the arguments which were not passed as a paremeter
        submission_file = join(
            PROVIDENTIA_ROOT,
            "logs/interpolation/management_logs",
            f"{self.unique_id}.out",
        )
        with open(submission_file, "r") as f:
            submission_file_txt = f.read().split()

        # join unclosed lists together in 1 element
        submission_file_txt_joined = []
        cur = ""
        n = 0
        for t in submission_file_txt:
            cur += " " + t if cur else t
            n += t.count("[") - t.count("]")
            if n == 0:
                submission_file_txt_joined.append(cur)
                cur = ""

        # get configuration variables from the management_logs
        for variable_key in [
            "ghost_version",
            "forecast",
            "interp_spinup_timesteps",
            "interp_model_downsampling",
            "interp_model_upsampling",
            "interp_n_neighbours",
            "interp_reverse_vertical_orientation",
            "interp_cleanup",
            "mod_root",
            "ghost_root",
            "mod_to_interp_root",
            "nonghost_root",
        ]:
            variable_val_idx = submission_file_txt_joined.index(variable_key + ":") + 1
            variable_val = submission_file_txt_joined[variable_val_idx]
            # make sure that boolean, None and lists are set correctly
            literals = {"true": True, "false": False, "none": None, "[]": []}
            if variable_val.lower() in literals:
                setattr(self, variable_key, literals[variable_val.lower()])
            else:
                setattr(self, variable_key, variable_val)

        # import GHOST standards
        self.standard_parameter_speci = get_standard_parameters_by_speci(
            self.original_mod_speci_to_process, self.ghost_version
        )

        # get GHOST lower/upper limits for variable
        self.GHOST_speci_lower_limit = self.standard_parameter_speci[
            "extreme_lower_limit"
        ]
        self.GHOST_speci_upper_limit = self.standard_parameter_speci[
            "extreme_upper_limit"
        ]

        mod_dir = None
        # for HPC machines, search in interp_models
        if MACHINE != "local":
            # get model type and specific directories
            for model_type, model_dict in interp_models.items():
                if self.experiment_to_process in model_dict["models"]:
                    # take first functional directory
                    for temp_mod_dir in model_dict["paths"]:
                        temp_mod_dir = join(temp_mod_dir, self.experiment_to_process)
                        if os.path.exists(temp_mod_dir):
                            mod_dir = temp_mod_dir
                            break
                    break

        # if local machine or if not mod_dir, get directory from data_paths
        if mod_dir is None:
            mod_to_interp_path = join(
                self.mod_to_interp_root, self.experiment_to_process
            )
            if os.path.exists(mod_to_interp_path):
                mod_dir = mod_to_interp_path

        # define if network is in GHOST format
        self.reading_ghost = check_for_ghost(self.network_to_interpolate_against)

        # get relevant observational file
        # GHOST
        if self.reading_ghost:
            self.obs_file = glob.glob(
                self.ghost_root
                + "/{}/{}/{}/{}/{}_{}*.nc".format(
                    self.network_to_interpolate_against,
                    self.ghost_version,
                    self.temporal_resolution_to_output,
                    self.obs_speci_to_process,
                    self.obs_speci_to_process,
                    self.yearmonth,
                )
            )[0]
        # non-GHOST
        else:
            self.obs_file = glob.glob(
                self.nonghost_root
                + "/{}/{}/{}/{}_{}*.nc".format(
                    self.network_to_interpolate_against,
                    self.temporal_resolution_to_output,
                    self.obs_speci_to_process,
                    self.obs_speci_to_process,
                    self.yearmonth,
                )
            )[0]

        # set obs file for getting units and boolean setting if can convert units or not
        self.obs_file_units = copy.deepcopy(self.obs_file)
           
        # get relevant observational file for unit conversion if are interpolating between species (e.g. sconco3@t2)
        if '@' in self.speci_to_process:
            # GHOST
            if self.reading_ghost:
                obs_files = glob.glob(self.ghost_root + '/{}/{}/{}/{}/{}_{}*.nc'\
                              .format(self.network_to_interpolate_against, self.ghost_version,
                                      self.temporal_resolution_to_output, self.original_mod_speci_to_process,
                                      self.original_mod_speci_to_process, self.yearmonth))
            # non-GHOST
            else:
                obs_files = glob.glob(self.nonghost_root + '/{}/{}/{}/{}_{}*.nc'\
                             .format(self.network_to_interpolate_against, self.temporal_resolution_to_output,
                                     self.original_mod_speci_to_process, self.original_mod_speci_to_process, self.yearmonth))
                
            # if have valid obs file then take it
            if len(obs_files) > 0:
                self.obs_file_units = obs_files[0]
            else:
                self.obs_file_units = None

        # get relevant model files
        if self.ensemble_member:
            all_model_files = np.sort(
                glob.glob(
                    "{}/{}/{}/{}/{}*{}*.nc".format(
                        mod_dir,
                        self.grid_type,
                        self.model_temporal_resolution,
                        self.mod_speci_to_process,
                        self.mod_speci_to_process,
                        self.yearmonth,
                    )
                )
            )

            # drop all analysis files ending with '_an.nc' which are not in ensemble-stats
            all_model_files = [f for f in all_model_files if "_an.nc" not in f]

            # isolate model files to be only those associated with relevant ensemble
            self.model_files = np.sort(
                [
                    f
                    for f in all_model_files
                    if "{}-{}_".format(self.mod_speci_to_process, self.ensemble) in f
                ]
            )

            # if the number of remaining model files is 0, this is because the files
            # do not contain an ensemble member number, therefore take all model files
            if len(self.model_files) == 0:
                self.model_files = all_model_files

        else:
            self.model_files = np.sort(
                glob.glob(
                    "{}/{}/{}/ensemble-stats/{}_{}/{}*{}*{}.nc".format(
                        mod_dir,
                        self.grid_type,
                        self.model_temporal_resolution,
                        self.mod_speci_to_process,
                        self.ensemble,
                        self.mod_speci_to_process,
                        self.yearmonth,
                        self.ensemble,
                    )
                )
            )

        # if have forecast data, then get the max number of complete forecast days across model files
        if self.forecast:
            # ititalise max forecast days as 0
            max_forecast_days = 0

            # iterate and read chunked model files
            for model_ii, model_file in enumerate(self.model_files):
                # read in chunked netcdf file
                mod_nc_root = Dataset(model_file)

                # check if have time dimension in daily file, if do not, do not process file
                if "time" not in list(mod_nc_root.dimensions.keys()):
                    mod_nc_root.close()
                    continue

                # get date from filename
                if not self.ensemble_member:
                    file_date = model_file.replace("_" + self.ensemble, "").split("_")[
                        -1
                    ][:-3]
                else:
                    file_date = model_file.split("_")[-1][:-3]

                # get file chunk type
                if len(file_date) == 6:
                    chunk_type = "monthly"
                elif len(file_date) in [8, 10]:
                    chunk_type = "daily"
                else:
                    mod_nc_root.close()
                    continue

                # cannot have chunked monthly file which has forecast data, so throw an error
                if chunk_type == "monthly":
                    self.log_file_str += "File {} is monthly, for which forecast data cannot be processed. Terminating process.".format(
                        model_file
                    )
                    create_output_logfile(1, self.log_file_str)
                # cannot have monthly resolution file which has forecast data, so throw an error
                elif self.model_temporal_resolution == "monthly":
                    self.log_file_str += "File {} resolution is monthly, for which forecast data cannot be processed. Terminating process.".format(
                        model_file
                    )
                    create_output_logfile(1, self.log_file_str)

                # get number of timesteps in file
                n_timesteps = mod_nc_root.dimensions["time"].size

                # get number of complete forecast days
                if self.model_temporal_resolution in ["hourly", "hourly_instantaneous"]:
                    forecast_days = int(n_timesteps / 24.0)
                elif self.model_temporal_resolution in [
                    "3hourly",
                    "3hourly_instantaneous",
                ]:
                    forecast_days = int(n_timesteps / (24.0 / 3.0))
                elif self.model_temporal_resolution in [
                    "6hourly",
                    "6hourly_instantaneous",
                ]:
                    forecast_days = int(n_timesteps / (24.0 / 6.0))
                elif self.model_temporal_resolution == "daily":
                    forecast_days = int(n_timesteps)

                # determine if current file has more complete forecast days than the previous maximum
                if forecast_days > max_forecast_days:
                    max_forecast_days = forecast_days

                # close model netCDF root
                mod_nc_root.close()

            # if max forecast days is 1 or less then forecast is not active, as have no forecast data, so deactivate it
            if max_forecast_days <= 1:
                self.forecast = None
                self.forecast_days = 1
            # otherwise set forecast days  as max forecast days
            else:
                self.forecast_days = max_forecast_days

        # else, set forecast days as 1 (as have no forecast)
        else:
            self.forecast_days = 1

        # if have hour in model fname, and it is not 0 then insert previous file also to get times from previous yearmonth that are offset
        # if forecast is set also get times from previous yearmonth to get prior forecasts for first days of the month
        if not self.ensemble_member:
            first_file_date = (
                self.model_files[0].replace("_" + self.ensemble, "").split("_")[-1][:-3]
            )
        else:
            first_file_date = self.model_files[0].split("_")[-1][:-3]

        if (len(first_file_date) == 10) or (self.forecast):
            if len(first_file_date) == 10:
                model_start_hour = int(first_file_date[-2:])
            else:
                model_start_hour = 0

            # model start hour is not 0, or forecast is set
            if (model_start_hour != 0) or (self.forecast):
                # get previous month
                prev_yearmonth = (
                    datetime.datetime.strptime(self.yearmonth, "%Y%m")
                    - relativedelta.relativedelta(months=1)
                ).strftime("%Y%m")
                # get previous month files
                if self.ensemble_member:
                    prev_month_files = np.sort(
                        glob.glob(
                            "{}/{}/{}/{}/{}*{}*.nc".format(
                                mod_dir,
                                self.grid_type,
                                self.model_temporal_resolution,
                                self.mod_speci_to_process,
                                self.mod_speci_to_process,
                                prev_yearmonth,
                            )
                        )
                    )
                else:
                    prev_month_files_final = np.sort(
                        glob.glob(
                            "{}/{}/{}/ensemble-stats/{}_{}/{}*{}*{}.nc".format(
                                mod_dir,
                                self.grid_type,
                                self.model_temporal_resolution,
                                self.mod_speci_to_process,
                                self.ensemble,
                                self.mod_speci_to_process,
                                prev_yearmonth,
                                self.ensemble,
                            )
                        )
                    )

                if self.ensemble_member:
                    # drop all analysis files ending with '_an.nc' which are not in ensemble-stats
                    prev_month_files = [
                        f for f in prev_month_files if "_an.nc" not in f
                    ]

                    # isolate model files to be only those associated with relevant ensemble
                    prev_month_files_final = np.sort(
                        [
                            f
                            for f in prev_month_files
                            if "{}-{}_".format(self.mod_speci_to_process, self.ensemble)
                            in f
                        ]
                    )

                    # if the number of remaining model files is 0, this is because the files
                    # do not contain an ensemble member number, therefore take all model files
                    if len(prev_month_files_final) == 0:
                        prev_month_files_final = prev_month_files

                # if have previous month files, then if forecast insert the previous N forecast day files
                # otherwise, insert the last one as first file in all_model_files
                if len(prev_month_files_final) > 0:
                    valid_prev_month_files = []
                    if self.forecast:
                        first_file_YYYYMMDD = datetime.datetime.strptime(
                            first_file_date[:8], "%Y%m%d"
                        )
                        for prev_month_file in prev_month_files_final:
                            if not self.ensemble_member:
                                file_YYYYMMDD = datetime.datetime.strptime(
                                    prev_month_file.replace(
                                        "_" + self.ensemble, ""
                                    ).split("_")[-1][:-3][:8],
                                    "%Y%m%d",
                                )
                            else:
                                file_YYYYMMDD = datetime.datetime.strptime(
                                    self.model_files[0].split("_")[-1][:-3][:8],
                                    "%Y%m%d",
                                )

                            if model_start_hour != 0:
                                N_prev_days = self.forecast_days + 1
                            else:
                                N_prev_days = self.forecast_days

                            if (
                                (first_file_YYYYMMDD - file_YYYYMMDD).days
                            ) <= N_prev_days:
                                valid_prev_month_files.append(prev_month_file)
                    else:
                        valid_prev_month_files.append(prev_month_files_final[-1])

                    # add valid previous month files to start of model files array
                    self.model_files = np.insert(
                        self.model_files, 0, valid_prev_month_files
                    )

        # if have no valid model files, exit
        if len(self.model_files) == 0:
            self.log_file_str += "No valid model files in {}. Skipping month.".format(
                self.yearmonth
            )
            create_output_logfile(1, self.log_file_str)

    def get_model_information(self):
        """
        Take first valid model file in month and get grid dimension/coordinate information.
        Put initial object read in a  try/except to handle reading of corrupted files.
        Iterate through files until have read a valid file.
        If do not read a valid file, skip month.
        """

        for model_file_ii, model_file in enumerate(self.model_files):
            try:
                # load instance of model file netCDF
                mod_nc_root = Dataset(model_file)

                # get instance of species variable
                mod_speci_obj = mod_nc_root[self.mod_speci_to_process]

                # get species units
                if hasattr(mod_speci_obj, "units"):
                    self.mod_speci_units = mod_speci_obj.units
                else:
                    self.log_file_str += f"Missing 'units' attribute for variable '{self.mod_speci_to_process}' in file {model_file}\n"
                    create_output_logfile(1, self.log_file_str)

                # get model grid type
                if hasattr(mod_speci_obj, "grid_mapping"):
                    self.mod_grid_type = mod_speci_obj.grid_mapping
                else:
                    self.log_file_str += f"Missing 'grid_mapping' attribute for variable '{self.mod_speci_to_process}' in file {model_file}\n"
                    create_output_logfile(1, self.log_file_str)

                # get individual dimension variable names
                # standard (no vertical dimension)
                self.have_vertical_dimension = False
                self.have_bin_dimension = False
                if len(mod_speci_obj.shape) == 3:
                    self.x_varname = mod_speci_obj.dimensions[2]
                    self.y_varname = mod_speci_obj.dimensions[1]
                # mapped size distribution variable, with bin dimension
                elif (len(mod_speci_obj.shape) == 4) and (
                    "vconcaerobin" in self.original_mod_speci_to_process
                ):
                    self.have_bin_dimension = True
                    self.x_varname = mod_speci_obj.dimensions[3]
                    self.y_varname = mod_speci_obj.dimensions[2]
                # with vertical dimension
                elif len(mod_speci_obj.shape) == 4:
                    self.have_vertical_dimension = True
                    self.x_varname = mod_speci_obj.dimensions[3]
                    self.y_varname = mod_speci_obj.dimensions[2]
                    self.z_varname = mod_speci_obj.dimensions[1]

                    # check if vertical dimension goes up or down to get correct index for surface
                    mod_vert_obj = mod_nc_root[self.z_varname]
                    direction = mod_vert_obj.positive

                    # if direction == 'up', surface index is location of mininum value in z var
                    if direction == "up":
                        if self.interp_reverse_vertical_orientation:
                            self.z_index = np.argmax(mod_vert_obj[:])
                        else:
                            self.z_index = np.argmin(mod_vert_obj[:])
                    # if direction == 'down', surface index is location of maximum value in z var
                    elif direction == "down":
                        if self.interp_reverse_vertical_orientation:
                            self.z_index = np.argmin(mod_vert_obj[:])
                        else:
                            self.z_index = np.argmax(mod_vert_obj[:])
                    # if cannot determine a surface index, terminate process
                    else:
                        self.log_file_str += "Cannot determine surface index in vertical dimension. Terminating process."
                        create_output_logfile(1, self.log_file_str)

                # check if species grid dimensions are named correctly, and in correct BSC standard order
                # if not terminate process
                # this is done by checking the variable names of the x, y (and z if required) dimensions
                # X dimension is valid if 'lon' is contained within name, or is == 'x'
                if ("lon" not in self.x_varname) and (self.x_varname != "x"):
                    self.log_file_str += (
                        "X dimension incorrectly named. Terminating process."
                    )
                    create_output_logfile(1, self.log_file_str)
                # Y dimension is valid if 'lat' is contained within name, or is == 'y'
                if ("lat" not in self.y_varname) and (self.y_varname != "y"):
                    self.log_file_str += (
                        "Y dimension incorrectly named. Terminating process."
                    )
                    create_output_logfile(1, self.log_file_str)
                # Z dimension is valid if == 'z' or 'lev' or 'alt' or 'height'
                if self.have_vertical_dimension:
                    if (
                        (self.z_varname != "lev")
                        and (self.z_varname != "z")
                        and (self.z_varname != "alt")
                        and (self.z_varname != "height")
                        and (self.z_varname != "level")
                    ):
                        self.log_file_str += (
                            "Z dimension incorrectly named. Terminating process."
                        )
                        create_output_logfile(1, self.log_file_str)

                # get instances of x/y grid dimension variables
                mod_lon_obj = mod_nc_root[self.x_varname]
                mod_lat_obj = mod_nc_root[self.y_varname]

                # get size of x/y grid dimensions
                self.x_N = mod_lon_obj.size
                self.y_N = mod_lat_obj.size

                # get 1D x/y values
                self.x = mod_lon_obj[:]
                self.y = mod_lat_obj[:]

                # get name of longitude/latitude grid centre coordinate variables
                if hasattr(mod_speci_obj, "coordinates"):
                    grid_centre_coordinates = mod_speci_obj.coordinates.split(" ")
                else:
                    self.log_file_str += f"Missing 'coordinates' attribute for variable '{self.mod_speci_to_process}' in file {model_file}\n"
                    create_output_logfile(1, self.log_file_str)

                lon_centre_varname = grid_centre_coordinates[1]
                lat_centre_varname = grid_centre_coordinates[0]

                # check if species grid centre coordinate are named correctly, and in correct BSC standard order
                # if not terminate process
                # longitude coordinate is valid if 'lon' is contained within name
                if "lon" not in lon_centre_varname:
                    self.log_file_str += "Longitude grid centre coordinate incorrectly named. Set it by defining variable coordinates. Terminating process."
                    create_output_logfile(1, self.log_file_str)
                # latitude coordinate is valid if 'lat' is contained within name
                if "lat" not in lat_centre_varname:
                    self.log_file_str += "Latitude grid centre coordinate incorrectly named. Set it by defining variable coordinates. Terminating process."
                    create_output_logfile(1, self.log_file_str)

                # get longitude and latitude grid centre values
                self.mod_lons_centre = np.float32(mod_nc_root[lon_centre_varname][:])
                self.mod_lats_centre = np.float32(mod_nc_root[lat_centre_varname][:])

                # check for unusual ranges
                lon_min = float(self.mod_lons_centre.min())
                lon_max = float(self.mod_lons_centre.max())
                lat_min = float(self.mod_lats_centre.min())
                lat_max = float(self.mod_lats_centre.max())
                if (lon_min < -180) or (lon_max > 360):
                    self.log_file_str += "Longitude range is unusual ({}, {}) - needs manual inspection. Terminating process.".format(
                        lon_min, lon_max
                    )
                    create_output_logfile(1, self.log_file_str)
                if (lat_min < -90) or (lat_max > 180):
                    self.log_file_str += "Latitude range is unusual ({}, {}) - needs manual inspection. Terminating process.".format(
                        lat_min, lat_max
                    )
                    create_output_logfile(1, self.log_file_str)

                # check if need to order or shift coordinates

                # need to shift coordinates
                # longitudes are 0-360? (0 centred over Greenwich)
                self.shift_lon = lon_max > 180
                # latitudes are 0-180? (0 = South Pole, 90 = Equator, 180 = North Pole)
                self.shift_lat = (lat_min >= 0) and (lat_max > 90)

                # shift coordinates
                if self.shift_lon:
                    self.log_file_str += (
                        "Shifting longitudes from 0-360 to -180-180. \n"
                    )
                    self.mod_lons_centre = ((self.mod_lons_centre + 180) % 360) - 180
                    if self.mod_grid_type == "crs":
                        self.x = ((self.x + 180) % 360) - 180
                if self.shift_lat:
                    self.log_file_str += "Shifting latitudes from 0-180 to -90-90. \n"
                    self.mod_lats_centre = self.mod_lats_centre - 90
                    if self.mod_grid_type == "crs":
                        self.y = self.y - 90

                # initialise tolerance    
                tol = 1e-10

                # need to order coordinates?
                # geographic grid?
                if self.mod_grid_type == "crs":
                    x_diff = np.diff(self.x)
                    y_diff = np.diff(self.y)
                    increasing_x = np.all(x_diff > tol)
                    increasing_y = np.all(y_diff > tol)
                    self.order_lon = not increasing_x
                    self.order_lat = not increasing_y
                # projected grid?
                else:
                    lon_test = self.mod_lons_centre[0, :]
                    lat_test = self.mod_lats_centre[:, 0]
                    lon_diff = np.diff(lon_test)
                    lat_diff = np.diff(lat_test)
                    increasing_lon = np.all(lon_diff > tol)
                    increasing_lat = np.all(lat_diff > tol)
                    self.order_lon = not increasing_lon
                    self.order_lat = not increasing_lat

                # order longitude
                if self.order_lon:
                    self.log_file_str += "Ordering longitudes. \n"
                    # geographic grid?
                    if self.mod_grid_type == "crs":
                        self.x_idx = np.argsort(self.x)
                        self.x = self.x[self.x_idx]
                        self.mod_lons_centre = self.mod_lons_centre[self.x_idx]
                    # projected grid?
                    else:
                        self.x_idx = np.argsort(self.mod_lons_centre[0, :])
                        self.mod_lons_centre = self.mod_lons_centre[:, self.x_idx]
                        self.mod_lats_centre = self.mod_lats_centre[:, self.x_idx]
                        self.x = self.x[self.x_idx]

                # order latitude
                if self.order_lat:
                    self.log_file_str += "Ordering latitudes. \n"
                    # geographic grid?
                    if self.mod_grid_type == "crs":
                        self.y_idx = np.argsort(self.y)
                        self.y = self.y[self.y_idx]
                        self.mod_lats_centre = self.mod_lats_centre[self.y_idx]
                    # projected grid?
                    else:
                        self.y_idx = np.argsort(self.mod_lats_centre[:, 0])
                        self.mod_lons_centre = self.mod_lons_centre[self.y_idx, :]
                        self.mod_lats_centre = self.mod_lats_centre[self.y_idx, :]
                        self.y = self.y[self.y_idx]

                # Check for monotonicity before ordering coordinates if needed (with tolerance for small numerical issues)
                # TODO: Review if there was any reason why prior to this commit the monotonicity check was done before shift of coordinates
                # There are cases (as seen in interpolation of cams_forecast_ensemble) where we need to shift the coordinates
                # before the monotonicity check, as the original coordinates are in 0-360 format and therefore not monotonic,
                # but once shifted to -180-180 they are monotonic.
                
                x_diff = np.diff(self.x)
                y_diff = np.diff(self.y)
                increasing_x = np.all(x_diff > tol)
                decreasing_x = np.all(x_diff < -tol)
                increasing_y = np.all(y_diff > tol)
                decreasing_y = np.all(y_diff < -tol)
                
                if not (increasing_x or decreasing_x):
                    # remove scientific notation
                    np.set_printoptions(suppress=True)
                    self.log_file_str += (
                        "Longitude is not monotonic. Terminating process."
                    )
                    self.log_file_str += (
                        "\nLongitudes: \n{0}.\nDifferences: \n{1}".format(
                            self.x, x_diff
                        )
                    )
                    create_output_logfile(1, self.log_file_str)
                if not (increasing_y or decreasing_y):
                    # remove scientific notation
                    np.set_printoptions(suppress=True)
                    self.log_file_str += (
                        "Latitude is not monotonic. Terminating process."
                    )
                    self.log_file_str += (
                        "\nLatitudes: \n{0}.\nDifferences: \n{1}".format(self.y, y_diff)
                    )
                    create_output_logfile(1, self.log_file_str)

                # close model netCDF root
                mod_nc_root.close()

                # get the last valid model file
                self.valid_model_file = copy.deepcopy(model_file)

                # break out of for loop, now that have read a valid model file in the month
                break

            except Exception as e:
                # in case of exception close model netCDF root
                mod_nc_root.close()

                # if have got to last file of month and that is corrupted, return from function
                if model_file_ii == (len(self.model_files) - 1):
                    self.log_file_str += (
                        "All model files corrupted in {}. Skipping month. \n".format(
                            self.yearmonth
                        )
                    )
                    self.log_file_str += "Error: {}".format(e)
                    create_output_logfile(1, self.log_file_str)
                # else, continue to next file in month
                else:
                    self.log_file_str += f"File {model_file} failed with error: {e}. Trying to read next file.\n"
                    continue

        # correct units in case of having bin dimension
        if self.have_bin_dimension:
            self.mod_speci_units = self.standard_parameter_speci["standard_units"]

    def create_grid_domain_edge_polygon(self):
        """
        Create grid domain edge polygon from model netCDF file.
        This is handled differently for regular grids (i.e. following lines of longitude/latitude),
        and non-regular grids (e.g. lambert-conformal).
        """
        # load instance of model file netCDF
        mod_nc_root = Dataset(self.valid_model_file)

        # if grid type == 'crs, then is regular grid
        if self.mod_grid_type == "crs":
            # set general grid type flag
            self.general_grid_type = "regular"

            # set longitude/latitude centre coordinates
            x_centre = self.mod_lons_centre
            y_centre = self.mod_lats_centre

        # here handle non-regular grids
        # currently only rotated pole and lambert conformal are supported
        elif self.mod_grid_type in [
            "rotated_pole",
            "Lambert_conformal",
            "Lambert_Conformal",
        ]:
            # set general grid type flag
            self.general_grid_type = "non-regular"

            # get instance of variable which defines key variables associated with the non-regular grid
            non_regular_grid_type_obj = mod_nc_root[self.mod_grid_type]

            # get non-regular grid-specific variables used for defining coordinate reference system
            if self.mod_grid_type == "rotated_pole":
                pole_longitude = np.float32(
                    non_regular_grid_type_obj.grid_north_pole_longitude
                )
                pole_latitude = np.float32(
                    non_regular_grid_type_obj.grid_north_pole_latitude
                )

                # add workaround to fix interpolation in the southern hemisphere
                # if the centre gridbox value (average of 2, if even number) of the geographic latitude variable is <0:
                # then add 'central_rotated_longitude = 180.0'
                lat_centre_i = findMiddle(self.mod_lats_centre.shape[0])
                lon_centre_i = findMiddle(self.mod_lats_centre.shape[1])
                centre_lat = np.average(
                    self.mod_lats_centre[lat_centre_i, lon_centre_i]
                )
                if centre_lat < 0.0:
                    central_rotated_longitude = 180.0
                else:
                    central_rotated_longitude = 0.0

            elif self.mod_grid_type in ["Lambert_conformal", "Lambert_Conformal"]:
                standard_parallels = np.float32(
                    non_regular_grid_type_obj.standard_parallel.split(",")
                )
                central_longitude = np.float32(
                    non_regular_grid_type_obj.longitude_of_central_meridian
                )
                central_latitude = np.float32(
                    non_regular_grid_type_obj.latitude_of_projection_origin
                )

            # read in non-regular grid x/y grid centre coordinates
            x_centre = np.float32(mod_nc_root[self.x_varname][:])
            y_centre = np.float32(mod_nc_root[self.y_varname][:])

        # the grid type cannot be handled, therefore terminate process
        else:
            self.log_file_str += (
                "Cannot handle grid of type: {}. Terminating process".format(
                    self.mod_grid_type
                )
            )
            create_output_logfile(1, self.log_file_str)

        # get x/y grid resolution (taken from average of increment between x/y grid centres)
        x_res = np.mean(np.diff(x_centre))
        y_res = np.mean(np.diff(y_centre))

        # if either of the x/y latitude increment resolutions are negative, this is because they run unconventionally
        # from east to west/south to north
        # set indices for indexing accordingly
        if x_res < 0:
            x_left_ind = -1
            x_right_ind = 0
        else:
            x_left_ind = 0
            x_right_ind = -1
        if y_res < 0:
            y_bottom_ind = -1
            y_top_ind = 0
        else:
            y_bottom_ind = 0
            y_top_ind = -1

        # force x/y resolutions increments to be positive
        x_res = np.abs(x_res)
        y_res = np.abs(y_res)

        # get x/y coordinates of grid edges
        x_edge = x_centre - (x_res / 2.0)
        if x_right_ind == -1:
            x_edge = np.append(x_edge, x_edge[x_right_ind] + x_res)
        else:
            x_edge = np.insert(x_edge, 0, x_edge[x_right_ind] + x_res)
        y_edge = y_centre - (y_res / 2.0)
        if y_top_ind == -1:
            y_edge = np.append(y_edge, y_edge[y_top_ind] + y_res)
        else:
            y_edge = np.insert(y_edge, 0, y_edge[y_top_ind] + y_res)

        # get x/y coordinates around grid edges
        left_edge_x = np.full(len(y_edge), x_edge[x_left_ind])
        right_edge_x = np.full(len(y_edge), x_edge[x_right_ind])
        if x_left_ind == 0:
            top_edge_x = x_edge[1:-1]
            bottom_edge_x = x_edge[-2::-1]
        else:
            top_edge_x = x_edge[-2::-1]
            bottom_edge_x = x_edge[1:-1]
        x_grid_edge = np.concatenate(
            (left_edge_x, top_edge_x, right_edge_x, bottom_edge_x)
        )

        bottom_edge_y = np.full(len(x_edge) - 1, y_edge[y_bottom_ind])
        top_edge_y = np.full(len(x_edge) - 2, y_edge[y_top_ind])
        if y_bottom_ind == 0:
            left_edge_y = y_edge[:]
            right_edge_y = y_edge[::-1]
        else:
            left_edge_y = y_edge[::-1]
            right_edge_y = y_edge[:]
        y_grid_edge = np.concatenate(
            (left_edge_y, top_edge_y, right_edge_y, bottom_edge_y)
        )

        # regular grid type?
        if self.general_grid_type == "regular":
            # stack longitude/latitude bounding edge coordinate pairs
            self.model_grid_outline = np.vstack((x_grid_edge, y_grid_edge)).T

            # get 2D mesh of model longitude/latitude gridcell centres
            self.mod_lons_centre, self.mod_lats_centre = np.meshgrid(
                self.mod_lons_centre, self.mod_lats_centre
            )

        # non-regular grid type?
        elif self.general_grid_type == "non-regular":
            # create cartopy coordinate reference system for the specific type non-standard grid, on WGS84 ellipsoid
            if self.mod_grid_type == "rotated_pole":
                non_regular_grid_crs = ccrs.RotatedPole(
                    pole_longitude=pole_longitude,
                    pole_latitude=pole_latitude,
                    central_rotated_longitude=central_rotated_longitude,
                )
            elif self.mod_grid_type in ["Lambert_conformal", "Lambert_Conformal"]:
                non_regular_grid_crs = ccrs.LambertConformal(
                    central_longitude=central_longitude,
                    central_latitude=central_latitude,
                    standard_parallels=standard_parallels,
                )

            # define a regular gridded coordinate reference system (Plate Carree), on WGS84 ellipsoid
            plate_carree = ccrs.PlateCarree()

            # convert bounding coordinates of box defined in non-regular grid coordinates to regular grid coordinates
            self.model_grid_outline = plate_carree.transform_points(
                non_regular_grid_crs, x_grid_edge, y_grid_edge
            )[:, :2]

        # convert longitude/latitude bounding coordinate pairs to a shapely polygon
        self.model_grid_outline_poly = Polygon(self.model_grid_outline)

        # close model netCDF root - now have all neccessary grid information
        mod_nc_root.close()

    def get_observational_objects(self):
        """Get necessary observational objects."""

        # get observational file netCDF root
        obs_nc_root = Dataset(self.obs_file)

        # get measured observational variable units
        # if are interpolating between species, and have no observational file to get units from, 
        # then take from GHOST standards
        if self.obs_file_units is None:
            self.obs_units = self.standard_parameter_speci['standard_units']
        # if are interpolating between species, and have an observational file to get units from, 
        # then take it from that
        elif self.obs_file != self.obs_file_units:
            obs_nc_root_units = Dataset(self.obs_file_units)
            obs_measured_var_obj = obs_nc_root_units[self.original_mod_speci_to_process]
            self.obs_units = obs_measured_var_obj.units
            obs_nc_root_units.close()
        # otherwise take from standard observational file
        else:
            # get measured observational variable object 
            obs_measured_var_obj = obs_nc_root[self.obs_speci_to_process]
            self.obs_units = obs_measured_var_obj.units

        # station object
        # for GHOST always is "station_reference"
        # for non-GHOST, try for "station_reference", then "station_code", then "station_name".
        if self.reading_ghost:
            obs_station_reference_obj = obs_nc_root["station_reference"]
        else:
            if "station_reference" in obs_nc_root.variables:
                obs_station_reference_obj = obs_nc_root["station_reference"]
            elif "station_code" in obs_nc_root.variables:
                obs_station_reference_obj = obs_nc_root["station_code"]
            elif "station_name" in obs_nc_root.variables:
                obs_station_reference_obj = obs_nc_root["station_name"]

        # get station references atributes
        if self.reading_ghost:
            self.obs_station_reference_units = obs_station_reference_obj.units
            self.obs_station_reference_standard_name = (
                obs_station_reference_obj.standard_name
            )
            self.obs_station_reference_description = (
                obs_station_reference_obj.description
            )
            self.obs_station_reference_long_name = obs_station_reference_obj.long_name

        # lon/lat objects
        if "latitude" in obs_nc_root.variables:
            obs_lon_obj = obs_nc_root["longitude"]
            obs_lat_obj = obs_nc_root["latitude"]
        else:
            obs_lon_obj = obs_nc_root["lon"]
            obs_lat_obj = obs_nc_root["lat"]

        # get lon lat atributes
        self.obs_lat_obj_units, self.obs_lon_obj_units = (
            obs_lat_obj.units,
            obs_lon_obj.units,
        )
        if self.reading_ghost:
            self.obs_lat_obj_standard_name, self.obs_lon_obj_standard_name = (
                obs_lat_obj.standard_name,
                obs_lon_obj.standard_name,
            )
            self.obs_lat_obj_long_name, self.obs_lon_obj_long_name = (
                obs_lat_obj.long_name,
                obs_lon_obj.long_name,
            )
            self.obs_lat_obj_description, self.obs_lon_obj_description = (
                obs_lat_obj.description,
                obs_lon_obj.description,
            )

        # get station data
        # GHOST
        if self.reading_ghost:
            self.station_references = obs_station_reference_obj[:]
        # non-GHOST
        else:
            meta_shape = obs_station_reference_obj.shape
            self.station_references = obs_station_reference_obj[:]
            meta_val_dtype = np.array([self.station_references[0]]).dtype

            if len(meta_shape) == 2:
                if meta_val_dtype == np.dtype(object):
                    self.station_references = np.array(
                        ["".join(val) for val in self.station_references]
                    )
                else:
                    self.station_references = chartostring(self.station_references)

        # get indices of all non-NaN stations (can be NaN for some non-GHOST files)
        non_nan_station_indices = np.array(
            [
                ref_ii
                for ref_ii, ref in enumerate(self.station_references)
                if ref.lower() != "nan"
            ]
        )
        self.station_references = self.station_references[non_nan_station_indices]

        # get valid longitude and latitudes
        self.obs_lons = obs_lon_obj[non_nan_station_indices]
        self.obs_lats = obs_lat_obj[non_nan_station_indices]

        # close the nc file
        obs_nc_root.close()

    def get_resampling_direction(self):
        """
        Determine the resampling direction between model input and output temporal resolutions.

        Returns
        -------
        str or None
            'downsampling' if input is finer than output,
            'upsampling' if input is coarser than output,
            None if input and output resolutions are equal.
        """

        freq_map = {
            "hourly": 1,
            "hourly_instantaneous": 1,
            "3hourly": 3,
            "3hourly_instantaneous": 3,
            "6hourly": 6,
            "6hourly_instantaneous": 6,
            "daily": 24,
            "monthly": 24 * 30,
        }

        in_freq = freq_map[self.model_temporal_resolution]
        out_freq = freq_map[self.temporal_resolution_to_output]

        # input resolution finer than output resolution (downsampling)
        if in_freq < out_freq:
            return "downsampling"
        # input resolution coarser than output resolution (upsampling)
        elif in_freq > out_freq:
            return "upsampling"
        # no resampling
        else:
            return

    def get_monthly_model_data(self):
        """Read all relevant model data in yearmonth into memory."""

        # determine resampling direction (upsampling, downsampling or None)
        resampling_direction = self.get_resampling_direction()

        # get number of days in month processing
        self.days_in_month = monthrange(int(self.year), int(self.month))[1]

        # create monthly datetime array at resolution to output
        start_month_dt = datetime.datetime(
            year=int(self.year), month=int(self.month), day=1, hour=0, minute=0
        )
        end_month_dt = start_month_dt + relativedelta.relativedelta(months=1)
        if self.temporal_resolution_to_output in ["hourly", "hourly_instantaneous"]:
            self.yearmonth_time = np.arange(0, self.days_in_month * 24.0)
            if Version(pd.__version__) >= Version("2.2"):
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="h", inclusive="left"
                )
            elif Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="H", inclusive="left"
                )
            else:
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="H", closed="left"
                )
            self.descriptive_temporal_resolution = "hours"
            self.temporal_resolution_to_output_code = "H"
        elif self.temporal_resolution_to_output in ["3hourly", "3hourly_instantaneous"]:
            self.yearmonth_time = np.arange(0, self.days_in_month * 24.0, 3.0)
            if Version(pd.__version__) >= Version("2.2"):
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="3h", inclusive="left"
                )
            elif Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="3H", inclusive="left"
                )
            else:
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="3H", closed="left"
                )
            self.descriptive_temporal_resolution = "hours"
            self.temporal_resolution_to_output_code = "3H"
        elif self.temporal_resolution_to_output in ["6hourly", "6hourly_instantaneous"]:
            self.yearmonth_time = np.arange(0, self.days_in_month * 24.0, 6.0)
            if Version(pd.__version__) >= Version("2.2"):
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="6h", inclusive="left"
                )
            elif Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="6H", inclusive="left"
                )
            else:
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="6H", closed="left"
                )
            self.descriptive_temporal_resolution = "hours"
            self.temporal_resolution_to_output_code = "6H"
        elif self.temporal_resolution_to_output == "daily":
            self.yearmonth_time = np.arange(0, self.days_in_month)
            if Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="D", inclusive="left"
                )
            else:
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="D", closed="left"
                )
            self.descriptive_temporal_resolution = "days"
            self.temporal_resolution_to_output_code = "D"
        elif self.temporal_resolution_to_output == "monthly":
            self.yearmonth_time = np.arange(0, 1)
            if Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="MS", inclusive="left"
                )
            else:
                self.yearmonth_dt = pd.date_range(
                    start_month_dt, end_month_dt, freq="MS", closed="left"
                )
            self.descriptive_temporal_resolution = "months"
            self.temporal_resolution_to_output_code = "MS"

        # create array for storing daily monthly model data
        # if have forecast data, then have an extra dimension for the number of forecast days
        if self.forecast:
            self.monthly_model_data = np.full(
                (len(self.yearmonth_time), self.forecast_days, self.y_N, self.x_N),
                np.nan,
                dtype=np.float32,
            )
        else:
            self.monthly_model_data = np.full(
                (len(self.yearmonth_time), self.y_N, self.x_N), np.nan, dtype=np.float32
            )

        # if have mapped size distribution variable:
        # -- get aeronet bin radius for original variable
        # -- get index of model bin to use for transformation of bin types (and rmin/rmax)
        # -- get transform factor to go between model bin and aeronet bin
        if self.have_bin_dimension:
            aeronet_bin_radius = get_aeronet_bin_radius_from_bin_variable(
                self.original_mod_speci_to_process
            )
            bin_index, rmin, rmax, rho_bin = get_aeronet_model_bin(
                self.model_name, aeronet_bin_radius
            )
            bin_transform_factor = get_model_to_aeronet_bin_transform_factor(
                self.model_name, rmin, rmax
            )

        # iterate and read chunked model files
        failed_files = 0
        for model_ii, model_file in enumerate(self.model_files):
            # put model file read in try/except to catch corrupt model files
            try:
                # read in chunked netcdf file
                mod_nc_root = Dataset(model_file)

                # check if have time dimension in daily file, if do not, do not process file
                if "time" not in list(mod_nc_root.dimensions.keys()):
                    mod_nc_root.close()
                    self.log_file_str += "File {} is corrupt. Skipping.\n".format(
                        model_file
                    )
                    failed_files += 1
                    continue

                # get date from filename
                if not self.ensemble_member:
                    file_date = model_file.replace("_" + self.ensemble, "").split("_")[
                        -1
                    ][:-3]
                else:
                    file_date = model_file.split("_")[-1][:-3]

                # get file time (handle monthly resolution data differently to hourly/daily
                # as num2date does not support 'months since' units)
                file_time = mod_nc_root["time"][:]
                time_units = mod_nc_root["time"].units
                time_calendar = mod_nc_root["time"].calendar

                # if have >= spinup timesteps than the number of timesteps in the file, then throw error
                self.interp_spinup_timesteps = int(self.interp_spinup_timesteps)
                if self.interp_spinup_timesteps >= len(file_time):
                    self.log_file_str += "File {} has <= number of timesteps of 'interp_spinup_timesteps':{}. Terminating process.".format(
                        model_file, self.interp_spinup_timesteps
                    )
                    create_output_logfile(1, self.log_file_str)

                # convert file time to datetime
                if "months" in time_units:
                    monthly_start_date = time_units.split(" ")[2]
                    file_time_dt = pd.date_range(
                        start=monthly_start_date, periods=1, freq="MS"
                    )
                else:
                    file_time_dt = num2date(
                        file_time, units=time_units, calendar=time_calendar
                    )

                    # convert to pandas datetime
                    if Version(cftime.__version__) <= Version("1.0.3.4"):
                        # remove microseconds
                        file_time_dt = pd.to_datetime(
                            [t.replace(microsecond=0) for t in file_time_dt]
                        )
                    else:
                        # bug fix for newer versions of cftime
                        file_time_dt = file_time_dt.astype("datetime64[s]")
                        file_time_dt = pd.to_datetime([t for t in file_time_dt])

                # get file start and end datetime
                # monthly file
                if len(file_date) == 6:
                    start_file_dt = datetime.datetime(
                        year=int(file_date[:4]),
                        month=int(file_date[4:6]),
                        day=1,
                        hour=0,
                        minute=0,
                    )
                # daily file
                elif len(file_date) == 8:
                    start_file_dt = datetime.datetime(
                        year=int(file_date[:4]),
                        month=int(file_date[4:6]),
                        day=int(file_date[6:8]),
                        hour=0,
                        minute=0,
                    )
                # daily file
                elif len(file_date) == 10:
                    start_file_dt = datetime.datetime(
                        year=int(file_date[:4]),
                        month=int(file_date[4:6]),
                        day=int(file_date[6:8]),
                        hour=int(file_date[8:10]),
                        minute=0,
                    )
                else:
                    mod_nc_root.close()
                    self.log_file_str += 'Resolution could not be detected in {}, check the date in the filename as now it shows as "{}".\n'.format(
                        model_file, file_date
                    )
                    failed_files += 1
                    continue

                # for forecast data, get valid indices of file time per forecast day
                if self.forecast:
                    # get indices of file time within each forecast day (excluding spinup timsesteps)
                    valid_file_time_inds = []
                    for forecast_day in range(self.forecast_days):
                        forecast_start_dt = start_file_dt + datetime.timedelta(
                            days=forecast_day
                        )
                        forecast_end_dt = start_file_dt + datetime.timedelta(
                            days=forecast_day + 1
                        )
                        valid_inds = np.where(
                            (file_time_dt >= forecast_start_dt)
                            & (file_time_dt < forecast_end_dt)
                        )[0]
                        valid_file_time_inds.append(
                            np.array(
                                [
                                    valid_ind
                                    for valid_ind in valid_inds
                                    if valid_ind >= self.interp_spinup_timesteps
                                ],
                                dtype=np.int32,
                            )
                        )

                # for non-forecast runs, get all timesteps (excluding spinup timesteps)
                else:
                    # get indices of file time inside yearmonth
                    valid_inds = np.where(
                        (file_time_dt >= start_month_dt) & (file_time_dt < end_month_dt)
                    )[0]
                    valid_file_time_inds = [
                        np.array(
                            [
                                valid_ind
                                for valid_ind in valid_inds
                                if valid_ind >= self.interp_spinup_timesteps
                            ],
                            dtype=np.int32,
                        )
                    ]

                # set flag to indicate if have bad time format in file
                bad_time_format = False

                # iterate through forecast days
                for forecast_day_ii in range(self.forecast_days):
                    # cut file time dt for only valid times
                    # for forecast data ensure that this is the refers to the first day of forecast to get the correct time indices
                    valid_file_time_dt = file_time_dt[valid_file_time_inds[0]]

                    # read valid data from file for valid indices
                    # have bin dimension?
                    if self.have_bin_dimension:
                        read_data = mod_nc_root[self.mod_speci_to_process][
                            valid_file_time_inds[forecast_day_ii], bin_index, :, :
                        ]

                        # convert model units from kg m-2 to um3/um2
                        read_data = 1e6 * read_data / rho_bin

                        # transform model bin data to aeronet bin
                        read_data = read_data * bin_transform_factor
                    # has vertical dimension
                    elif self.have_vertical_dimension:
                        read_data = mod_nc_root[self.mod_speci_to_process][
                            valid_file_time_inds[forecast_day_ii], self.z_index, :, :
                        ]
                    # has no vertical dimension?
                    else:
                        read_data = mod_nc_root[self.mod_speci_to_process][
                            valid_file_time_inds[forecast_day_ii], :, :
                        ]

                    # convert model units to observational units
                    read_data = read_data * self.conversion_factor

                    # set any model values outside GHOST extreme limits for variable to be NaN
                    read_data[
                        (read_data < self.GHOST_speci_lower_limit)
                        | (read_data > self.GHOST_speci_upper_limit)
                    ] = np.nan

                    # reorder data coordinates if needed
                    if self.order_lon:
                        read_data = np.take(read_data, self.x_idx, axis=-1)
                    if self.order_lat:
                        read_data = np.take(read_data, self.y_idx, axis=-2)

                    # create xarray for resampling
                    xr_data = xr.DataArray(
                        dims=("time", "latitude", "longitude"),
                        data=read_data,
                        coords=dict(
                            time=valid_file_time_dt, latitude=self.y, longitude=self.x
                        ),
                    )

                    # do resampling
                    if resampling_direction:
                        # downsampling (finer to coarser)
                        if resampling_direction == "downsampling":
                            # mean
                            if self.interp_model_downsampling == "mean":
                                xr_data = xr_data.resample(
                                    time=self.temporal_resolution_to_output_code
                                ).mean()

                            # median
                            elif self.interp_model_downsampling == "median":
                                xr_data = xr_data.resample(
                                    time=self.temporal_resolution_to_output_code
                                ).median()

                        # upsampling (coarser to finer)
                        elif resampling_direction == "upsampling":
                            # convert original data frequency code to pandas offset
                            offset = pd.tseries.frequencies.to_offset(
                                pd.infer_freq(xr_data.time.to_index())
                            )

                            # extend the last timestamp to cover the full final period
                            start_time = pd.Timestamp(xr_data.time.values[0])
                            last_time = pd.Timestamp(xr_data.time.values[-1])
                            end_extended = last_time + offset - pd.Timedelta(seconds=1)

                            # create new continuous index at frequency to output
                            new_index = pd.date_range(
                                start=start_time,
                                end=end_extended,
                                freq=self.temporal_resolution_to_output_code,
                            )

                            # fill gaps (repeating values between measurements)
                            if self.interp_model_upsampling == "fill":
                                xr_data = xr_data.reindex(
                                    time=new_index, method="ffill"
                                )

                            # leave gaps (setting nans between measurements)
                            elif self.interp_model_upsampling == "gaps":
                                xr_data = xr_data.reindex(time=new_index)

                    # get indices in yearmonth to fill with read data. Adjust for forecast day if neccessary.
                    if (self.forecast) and (forecast_day_ii != 0):
                        file_time = xr_data.time.values + pd.Timedelta(
                            days=forecast_day_ii
                        )
                    else:
                        file_time = xr_data.time.values
                    inds_to_fill = np.isin(self.yearmonth_dt, file_time)

                    # check if any file fimes are not in monthly yearmonth time array
                    # this should not be the case unless time format is bad, or forecast data is being read
                    if (not any(inds_to_fill)) and (self.forecast):
                        continue
                    elif (not any(inds_to_fill)) and (not self.forecast):
                        if not bad_time_format:
                            self.log_file_str += "Time in model dataset {} is not in standard format: \n {}".format(
                                model_file, file_time
                            )
                            failed_files += 1
                            bad_time_format = True
                        continue

                    # fill in data array
                    if self.forecast:
                        self.monthly_model_data[
                            inds_to_fill, forecast_day_ii, :, :
                        ] = xr_data.values
                    else:
                        self.monthly_model_data[inds_to_fill, :, :] = xr_data.values

                # close model netCDF root
                mod_nc_root.close()

            except Exception:
                mod_nc_root.close()
                self.log_file_str += "File {} is corrupt. Skipping.\n{}".format(
                    model_file, traceback.format_exc()
                )
                failed_files += 1

        if failed_files == len(self.model_files):
            self.log_file_str += "No model dataset could be interpolated."
            create_output_logfile(1, self.log_file_str)

    def n_nearest_neighbour_inverse_distance_weights(self):
        """
        Calculate N nearest neighbour inverse distance weights (and indices) of model gridcells centres
        from an array of geographic observational station coordinates. Both observational and model geographic
        longitude/latitude coordinates are first converted to cartesian ECEF (Earth Centred, Earth Fixed)
        coordinates, before calculating distances. Weights returned for obervational stations not contained
        within model grid extents are all zero.
        """

        # for each pair of observational station geographic coordinates, test if station is inside model grid
        obs_inside_model_grid = np.array(
            [
                self.model_grid_outline_poly.contains(Point(obs_lon, obs_lat))
                for obs_lon, obs_lat in zip(self.obs_lons, self.obs_lats)
            ]
        )

        # flatten model centre lon/lat arrays
        # generate equal length altitude arrays
        obs_alts = np.zeros(len(self.obs_lons))
        flat_mod_lons_centre = self.mod_lons_centre.ravel()
        flat_mod_lats_centre = self.mod_lats_centre.ravel()
        flat_mod_alts = np.zeros(len(flat_mod_lats_centre))

        # convert observational/model geographic longitude/latitude coordinates to cartesian ECEF (Earth Centred,
        # Earth Fixed) coordinates, assuming WGS84 datum and ellipsoid, and that all heights = 0.
        # ECEF coordiantes represent positions (in meters) as X, Y, Z coordinates, approximating the earth surface
        # as an ellipsoid of revolution.
        # This conversion is for the subsequent calculation of euclidean distances of the model gridcell centres
        # from each observational station.
        # Defining the distance between two points on the earth's surface as simply the euclidean distance
        # between the two lat/lon pairs could lead to inaccurate results depending on the distance
        # between two points (i.e. 1 deg. of longitude varies with latitude).
        if Version(pyproj.__version__) >= Version("2.0.0"):
            lla = {"proj": "latlong", "ellps": "WGS84", "datum": "WGS84"}
            ecef = {"proj": "geocent", "ellps": "WGS84", "datum": "WGS84"}
            transformer = pyproj.Transformer.from_crs(lla, ecef)
            obs_x, obs_y, obs_z = transformer.transform(
                self.obs_lons, self.obs_lats, obs_alts, radians=False
            )
            mod_x, mod_y, mod_z = transformer.transform(
                flat_mod_lons_centre, flat_mod_lats_centre, flat_mod_alts, radians=False
            )
        else:
            lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
            ecef = pyproj.Proj(proj="geocent", ellps="WGS84", datum="WGS84")
            obs_x, obs_y, obs_z = pyproj.transform(
                lla, ecef, self.obs_lons, self.obs_lats, obs_alts, radians=False
            )
            mod_x, mod_y, mod_z = pyproj.transform(
                lla,
                ecef,
                flat_mod_lons_centre,
                flat_mod_lats_centre,
                flat_mod_alts,
                radians=False,
            )

        # stack converted cartesian coordinates for preparation of calculation of nearest neighbour distances
        # using Scipy cKDTree
        obs_lonlat = np.column_stack((obs_x, obs_y, obs_z))
        mod_lonlat = np.column_stack((mod_x, mod_y, mod_z))

        # generate cKDtree
        tree = cKDTree(mod_lonlat)

        # get n-neighbour nearest distances/indices (ravel form) of model gridcell centres from each observational station
        if Version(scipy.__version__) < Version("1.6.0"):
            dists, idx = tree.query(obs_lonlat, k=int(self.interp_n_neighbours))
        else:
            dists, idx = tree.query(
                obs_lonlat, k=int(self.interp_n_neighbours), workers=-1
            )

        # for n neighbours == 1, do reshaping of array so doesn't break
        if len(dists.shape) == 1:
            dists = np.reshape(dists, (dists.shape[0], 1))
            idx = np.reshape(idx, (idx.shape[0], 1))

        # if any dists are < 1m, then set as 1m (i.e. very small distance), to avoid issue calculating weight
        dists[dists < 1.0] = 1.0

        # get nearest neighbour indices
        self.nearest_neighbour_inds = np.column_stack(
            np.unravel_index(idx, self.mod_lons_centre.shape)
        )

        # take the reciprocals of the nearest neighbours distances from the observational points
        self.inverse_dists = np.reciprocal(dists)

        # set reciprocal distances for all observational points outside model grid extent to be 0
        self.inverse_dists[~obs_inside_model_grid, :] = 0.0

    def write_netCDF(self):
        """
        Write netCDF, returning interpolated model data to observational surface stations
        for yearmonth.
        """

        # set output directory where observational interpolated monthly model netcdf will be saved
        # GHOST
        if self.reading_ghost:
            network_name = self.network_to_interpolate_against
        # non-GHOST
        else:
            # as it appears in PRV (e.g. nasa-aeronet/oneill_v3-lev15 -> nasa-aeronet-oneill_v3-lev15)
            network_name = self.network_to_interpolate_against.replace("/", "-")
        output_dir = join(
            self.mod_root,
            self.ghost_version,
            self.prov_mod_code,
            self.temporal_resolution_to_output,
            self.speci_to_process,
            network_name,
        )

        # check if need to create any directories in path
        check_directory_existence(output_dir, self.mod_root)

        # create netCDF dataset
        netCDF_fname = "{}/{}_{}.nc".format(
            output_dir, self.speci_to_process, self.yearmonth
        )
        root_grp = Dataset(netCDF_fname, "w", format="NETCDF4")

        # auto mask arrays
        root_grp.set_auto_mask(True)

        # file contents
        msg = "Inverse distance weighting ({} neighbours) interpolated ".format(
            self.interp_n_neighbours
        )
        msg += "{} model data for the component {} ".format(
            self.experiment_to_process, self.original_mod_speci_to_process
        )
        if '@' in self.speci_to_process:
            msg += 'with reference to {} measurement stations in the '.format(self.obs_speci_to_process)
        else:
            msg += 'with reference to the measurement stations in the '
        msg += "{} network ".format(self.network_to_interpolate_against)
        msg += "in {}-{}.".format(self.year, self.month)
        root_grp.title = msg
        root_grp.institution = "Barcelona Supercomputing Center"
        root_grp.source = "Model {}".format(self.experiment_to_process)
        root_grp.creator_name = "Dene R. Bowdalo"
        root_grp.creator_email = "dene.bowdalo@bsc.es"

        # iterate through each station to process
        try:
            for ii, station_reference in enumerate(self.station_references):
                if ii == 0:
                    # netcdf dimensions
                    root_grp.createDimension("station", len(self.obs_lons))
                    root_grp.createDimension("time", len(self.yearmonth_time))
                    root_grp.createDimension("model_longitude", self.x_N)
                    root_grp.createDimension("model_latitude", self.y_N)
                    root_grp.createDimension(
                        "grid_edge", self.model_grid_outline.shape[0]
                    )
                    # if have forecast data then create forecast data dimension
                    if self.forecast:
                        root_grp.createDimension("forecast_day", self.forecast_days)

                    # create time variable
                    time_var = root_grp.createVariable("time", "u4", ("time"))
                    time_var.standard_name = "time"
                    time_var.long_name = "time"
                    time_var.units = "{} since {}-{}-01 00:00:00".format(
                        self.descriptive_temporal_resolution, self.year, self.month
                    )
                    msg = "Time in {} since {}-{}-01 00:00 UTC. ".format(
                        self.descriptive_temporal_resolution, self.year, self.month
                    )
                    msg += "Time given refers to the start of the time window the measurement "
                    msg += "is representative of (temporal resolution)"
                    time_var.description = msg
                    time_var.axis = "T"
                    time_var.calendar = "standard"
                    time_var.tz = "UTC"

                    # create observational equivalent station reference variable
                    station_reference_var = root_grp.createVariable(
                        "station_reference", str, ("station")
                    )
                    if self.reading_ghost:
                        station_reference_var.standard_name = (
                            self.obs_station_reference_standard_name
                        )
                        station_reference_var.long_name = (
                            self.obs_station_reference_long_name
                        )
                        station_reference_var.units = self.obs_station_reference_units
                        station_reference_var.description = (
                            self.obs_station_reference_description
                        )

                    # create observational equivalent longitude/latitude variables
                    longitude_var = root_grp.createVariable(
                        "longitude", "f8", ("station")
                    )
                    if self.reading_ghost:
                        longitude_var.standard_name = self.obs_lon_obj_standard_name
                        longitude_var.units = self.obs_lon_obj_units
                        longitude_var.long_name = self.obs_lon_obj_long_name
                        longitude_var.description = self.obs_lon_obj_description
                    longitude_var.axis = "X"

                    latitude_var = root_grp.createVariable(
                        "latitude", "f8", ("station")
                    )
                    if self.reading_ghost:
                        latitude_var.standard_name = self.obs_lat_obj_standard_name
                        latitude_var.units = self.obs_lat_obj_units
                        latitude_var.long_name = self.obs_lat_obj_long_name
                        latitude_var.description = self.obs_lat_obj_description
                    latitude_var.axis = "Y"

                    # create 2D meshed longitude/latitude gridcell centre variables
                    model_centre_longitude_var = root_grp.createVariable(
                        "model_centre_longitude",
                        "f8",
                        ("model_latitude", "model_longitude"),
                    )
                    model_centre_longitude_var.standard_name = "model centre longitude"
                    model_centre_longitude_var.long_name = "model centre longitude"
                    model_centre_longitude_var.units = self.obs_lon_obj_units
                    msg = "2D meshed grid centre longitudes with "
                    msg += "{} longitudes in {} bands of latitude".format(
                        self.x_N, self.y_N
                    )
                    model_centre_longitude_var.description = msg
                    model_centre_longitude_var.axis = "X"

                    model_centre_latitude_var = root_grp.createVariable(
                        "model_centre_latitude",
                        "f8",
                        ("model_latitude", "model_longitude"),
                    )
                    model_centre_latitude_var.standard_name = "model centre latitude"
                    model_centre_latitude_var.long_name = "model centre latitude"
                    model_centre_latitude_var.units = self.obs_lat_obj_units
                    msg = "2D meshed grid centre longitudes with "
                    msg += "{} longitudes in {} bands of latitude".format(
                        self.y_N, self.x_N
                    )
                    model_centre_latitude_var.description = msg
                    model_centre_latitude_var.axis = "Y"

                    # create grid domain edge longitude/latitude variables
                    grid_edge_longitude_var = root_grp.createVariable(
                        "grid_edge_longitude", "f8", ("grid_edge")
                    )
                    grid_edge_longitude_var.standard_name = "grid edge longitude"
                    grid_edge_longitude_var.long_name = "grid edge longitude"
                    grid_edge_longitude_var.units = self.obs_lon_obj_units
                    msg = "Longitude coordinate along edge of grid domain "
                    msg += "(going clockwise around grid boundary from bottom-left corner)."
                    grid_edge_longitude_var.description = msg
                    grid_edge_longitude_var.axis = "X"

                    grid_edge_latitude_var = root_grp.createVariable(
                        "grid_edge_latitude", "f8", ("grid_edge")
                    )
                    grid_edge_latitude_var.standard_name = "grid edge latitude"
                    grid_edge_latitude_var.long_name = "grid edge latitude"
                    grid_edge_latitude_var.units = self.obs_lat_obj_units
                    msg = "Latitude coordinate along edge of grid domain "
                    msg += "(going clockwise around grid boundary from bottom-left corner)."
                    grid_edge_latitude_var.description = msg
                    grid_edge_latitude_var.axis = "Y"

                    # create measured variable
                    if self.forecast:
                        measured_var = root_grp.createVariable(
                            self.speci_to_process,
                            "f4",
                            ("station", "time", "forecast_day"),
                        )
                    else:
                        measured_var = root_grp.createVariable(
                            self.speci_to_process, "f4", ("station", "time")
                        )
                    # add attributes
                    measured_var.standard_name = self.speci_to_process
                    measured_var.units = self.obs_units
                    if '@' in self.speci_to_process:
                        measured_var.description = 'Interpolated value of {} from the model {} with reference to {} measurement stations in the {} network'.format(
                                                  self.original_mod_speci_to_process, self.experiment_to_process, self.obs_speci_to_process, self.network_to_interpolate_against)
                    else:
                        measured_var.description = 'Interpolated value of {} from the model {} with reference to the measurement stations in the {} network'.format(
                                                  self.original_mod_speci_to_process, self.experiment_to_process, self.network_to_interpolate_against)

                    # write to variables
                    time_var[:] = self.yearmonth_time
                    station_reference_var[:] = self.station_references
                    longitude_var[:] = self.obs_lons
                    latitude_var[:] = self.obs_lats
                    model_centre_longitude_var[:] = self.mod_lons_centre
                    model_centre_latitude_var[:] = self.mod_lats_centre
                    grid_edge_longitude_var[:] = self.model_grid_outline[:, 0]
                    grid_edge_latitude_var[:] = self.model_grid_outline[:, 1]

                # iterate through observational stations
                # use calculated interpolated weights per station to model grid to produce model reciprocal output
                # if all station weights are 0, station is outside model grid domain
                # set all values to be NaN
                station_weights = self.inverse_dists[ii, :]
                if np.all(station_weights == 0):
                    if self.forecast:
                        interp_vals = np.full(
                            (len(self.yearmonth_time), self.forecast_days),
                            np.nan,
                            dtype=np.float32,
                        )
                    else:
                        interp_vals = np.full(
                            len(self.yearmonth_time), np.nan, dtype=np.float32
                        )
                else:
                    # get reciprocal model data at N nearest neighbours to observational station
                    if self.forecast:
                        cut_model_data = self.monthly_model_data[
                            :,
                            :,
                            self.nearest_neighbour_inds[
                                ii, : int(self.interp_n_neighbours)
                            ],
                            self.nearest_neighbour_inds[
                                ii, int(self.interp_n_neighbours) :
                            ],
                        ]
                    else:
                        cut_model_data = self.monthly_model_data[
                            :,
                            self.nearest_neighbour_inds[
                                ii, : int(self.interp_n_neighbours)
                            ],
                            self.nearest_neighbour_inds[
                                ii, int(self.interp_n_neighbours) :
                            ],
                        ]
                    # create mask where data == NaN or infinite
                    invalid_mask = ~np.isfinite(cut_model_data)

                    # create masked array
                    cut_model_data = np.ma.MaskedArray(
                        cut_model_data, mask=invalid_mask
                    )
                    # interpolate masked array across time dimension using interpolated weights per station
                    interp_vals = np.ma.average(
                        cut_model_data, weights=station_weights, axis=-1
                    )

                # write measured variable
                measured_var[ii, :] = interp_vals

            # close writing to netCDF
            root_grp.close()

        except Exception as e:
            self.log_file_str += "File {} could not be written. Error: {}.\n".format(
                netCDF_fname, e
            )
            create_output_logfile(1, self.log_file_str)

        # compress netCDF file
        try:
            subprocess.Popen(
                ["ncks", "-O", "--dfl_lvl", "1", netCDF_fname, netCDF_fname],
                stdout=subprocess.PIPE,
            )
        except:
            self.log_file_str += (
                "NCO could not be found, please install it in your system "
            )
            self.log_file_str += (
                "with conda install -c conda-forge nco --override-channels.\n"
            )

        # give 770 permissions for file and make owner bsc32 if machine isn't local
        if MACHINE != "local":
            set_file_permissions_ownership(netCDF_fname)

        # copy file to esarchive (if have access and is being written to exp_interp on gpfs)
        if MACHINE == "nord4":
            # only do copy if file has being written to exp_interp default gpfs directory
            if "exp_interp" in netCDF_fname:
                # set esarchive output dir
                esarchive_output_dir = "/esarchive/recon/prov_interp/{}".format(
                    "/".join(netCDF_fname.split("/exp_interp/")[1].split("/")[:-1])
                )

                # check if need to create any directories in path
                check_directory_existence(
                    esarchive_output_dir, "/esarchive/recon/prov_interp"
                )

                # set esarchive fname
                esarchive_netCDF_fname = "{}/{}".format(
                    esarchive_output_dir, netCDF_fname.split("/")[-1]
                )

                try:
                    # copy file (without permissions)
                    shutil.copyfile(netCDF_fname, esarchive_netCDF_fname)

                    # give 770 permissions for file and make owner bsc32
                    set_file_permissions_ownership(esarchive_netCDF_fname)
                except (PermissionError, FileNotFoundError) as e:
                    self.log_file_str += "Interpolated file/s could not be copied to esarchive. Error: {}".format(
                        e
                    )
                    create_output_logfile(1, self.log_file_str)

    def cleanup(self):
        """
        Remove non-inteprolated model files used for interpolation.
        """

        self.log_file_str += "\nRemoving non-interpolated model files:\n"

        for model_file in self.model_files:
            try:
                os.remove(model_file)
                self.log_file_str += "Model file {} removed.\n".format(model_file)
                
            except Exception as e:
                self.log_file_str += "Model file {} could not be removed. Error: {}.\n".format(
                    model_file, e
                )
                create_output_logfile(1, self.log_file_str)

def create_output_logfile(process_code, log_file_str):
    """
    Create a logfile for stating outcome of interpolation job.

    Parameters
    ----------
    process_code : int
        Code indicating the outcome of the interpolation process:
        - 0 : Process completed without issue
        - 1 : Caught error in process
        - 2 : Uncaught error in process

    log_file_str : str
        Content to write in the logfile.
    """

    output_logfile_dir = (
        f"{join(PROVIDENTIA_ROOT, 'logs/interpolation/interpolation_logs/')}"
        f"{submit_args['prov_mod_code']}/"
        f"{submit_args['speci_to_process']}/"
        f"{submit_args['network_to_interpolate_against']}/"
        f"{submit_args['temporal_resolution_to_output']}/"
        f"{submit_args['yearmonth']}"
    )

    f = open(f"{output_logfile_dir}_{process_code}.out", "w")
    f.write(str(log_file_str))
    f.close()

    # exit from current process after writing logfile
    sys.exit()


if __name__ == "__main__":
    try:
        # time start of yearmonth interpolation
        interpolation_start = time.time()

        # get arguments passed from submittal script --> put into dict
        submit_args = {
            "prov_mod_code": sys.argv[1],
            "model_temporal_resolution": sys.argv[2],
            "mod_speci_to_process": sys.argv[3],
            "network_to_interpolate_against": sys.argv[4],
            "temporal_resolution_to_output": sys.argv[5],
            "yearmonth": sys.argv[6],
            "speci_to_process": sys.argv[7],
            "job_id": sys.argv[8],
        }

        # initialise ModelInterpolation object
        EI = ModelInterpolation(submit_args)

        # read model domain information
        EI.get_model_information()

        # create polygon along edge of model domain
        EI.create_grid_domain_edge_polygon()

        # get necessary data objects from the observational file
        EI.get_observational_objects()

        # get unit conversion factor between observations and model data
        EI.conversion_factor = get_conversion_factor(
            EI.mod_speci_units, EI.obs_units, EI.standard_parameter_speci
        )
        EI.log_file_str += (
            EI.mod_speci_units
            + " to "
            + EI.obs_units
            + " conversion factor is: "
            + str(EI.conversion_factor)
            + "\n"
        )
        if isinstance(EI.conversion_factor, str):
            EI.log_file_str += EI.conversion_factor
            create_output_logfile(1, EI.log_file_str)

        # read relevant monthly model data into memory
        EI.get_monthly_model_data()

        # get interpolation weights of model grid to observational stations
        # (using inverse distance weighting interpolation)
        EI.n_nearest_neighbour_inverse_distance_weights()

        # write out netCDF for yearmonth, interpolating model data to surface observational stations
        EI.write_netCDF()

        # remove non-inteprolated model files used for interpolation
        if EI.interp_cleanup:
            EI.cleanup()

        # get total time of interpolation
        interpolation_time = time.time() - interpolation_start

        # return valid process logfile (0)
        EI.log_file_str += str((time.time() - interpolation_start) / 60.0)
        create_output_logfile(0, EI.log_file_str)

    # write error log file if have uncaught internal error
    except Exception:
        log_file_str = "STARTING INTERPOLATION\n"
        log_file_str += str(traceback.format_exc())
        create_output_logfile(2, log_file_str)
