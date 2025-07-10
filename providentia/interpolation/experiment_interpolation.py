import ast
from calendar import monthrange
import cartopy.crs as ccrs
import cftime
import copy
import glob
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
import yaml
import pathlib
from pydoc import locate

import datetime
import dateutil.relativedelta as relativedelta
from netCDF4 import Dataset, num2date, date2num, chartostring
import numpy as np
from packaging.version import Version
import pandas as pd
import pyproj
import scipy
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, Point
import xarray as xr
from aux_interp import (check_for_ghost, findMiddle, check_directory_existence, set_file_permissions_ownership,
                 get_aeronet_bin_radius_from_bin_variable, get_aeronet_model_bin, 
                 get_model_to_aeronet_bin_transform_factor)

sys.path.append(str(pathlib.Path(__file__).resolve().parents[2]))
from providentia.auxiliar import CURRENT_PATH, join

MACHINE = os.environ.get('BSC_MACHINE', 'local')

# get current path and providentia root path
PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

# load the defined experiments paths and agrupations jsons
interp_experiments = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'interp_experiments.yaml')))

# add unit-converter submodule to python load path
sys.path.append(join(PROVIDENTIA_ROOT,'providentia','dependencies','unit-converter'))
import unit_converter

class ExperimentInterpolation(object):
    """ Class which handles interpolation of experiment data to surface observations. """

    def __init__(self,submit_args):
        
        self.log_file_str = 'STARTING INTERPOLATION\n'

        self.log_file_str += '\n'.join(f"{k}: {v}" for k, v in submit_args.items()) + '\n'

        # set variables from input keywords
        self.model_temporal_resolution            = submit_args['model_temporal_resolution']
        self.speci_to_process                     = submit_args['speci_to_process']
        self.network_to_interpolate_against       = submit_args['network_to_interpolate_against']
        self.temporal_resolution_to_output        = submit_args['temporal_resolution_to_output']
        self.yearmonth                            = submit_args['yearmonth']
        self.original_speci_to_process            = submit_args['original_speci_to_process']
        self.unique_id                            = submit_args['job_id']
        self.prov_exp_code                        = submit_args['prov_exp_code']
        self.experiment_to_process, self.grid_type, self.ensemble_option = self.prov_exp_code.split('-')
        
        # get year/month string
        self.year = self.yearmonth[:4]
        self.month = self.yearmonth[4:]

        # determine if ensemble option is member or emsemble stat
        self.ensemble_member = self.ensemble_option.isdigit()

        # dictionary to save utilized interpolation variables
        self.interpolation_variables = {}
        
        # from file in management_logs, get and set the arguments which were not passed as a paremeter
        submission_file = join(PROVIDENTIA_ROOT, 'logs/interpolation/management_logs',f'{self.unique_id}.out')
        with open(submission_file, 'r') as f:
            submission_file_txt = f.read().split()

        # get configuration variables from the management_logs
        for variable_key in ["ghost_version", "interp_n_neighbours", 
                             "interp_reverse_vertical_orientation", 
                             "exp_root", "ghost_root", "exp_to_interp_root", 
                             "nonghost_root", "forecast"]:
            variable_val_idx = submission_file_txt.index(variable_key+":")+1
            variable_val = submission_file_txt[variable_val_idx]
            # make sure vertical orientiation is a boolean!
            if variable_key in ["interp_reverse_vertical_orientation"]:
                setattr(self, variable_key, ast.literal_eval(variable_val))
            else:
                setattr(self, variable_key, variable_val)

        # import GHOST standards
        sys.path.insert(1, join(PROVIDENTIA_ROOT, 'providentia/dependencies/GHOST_standards/{}'.format(self.ghost_version)))
        from GHOST_standards import standard_parameters
        self.standard_parameters = standard_parameters
        
        # get cut of standard parameters for original speci to process
        for standard_parameter in self.standard_parameters.keys():
            if self.standard_parameters[standard_parameter]['bsc_parameter_name'] == self.original_speci_to_process:
                self.standard_parameter_speci = self.standard_parameters[standard_parameter]
                break
        
        # get GHOST lower/upper limits for variable
        self.GHOST_speci_lower_limit = self.standard_parameter_speci['extreme_lower_limit']
        self.GHOST_speci_upper_limit = self.standard_parameter_speci['extreme_upper_limit']

        exp_dir = None
        # for HPC machines, search in interp_experiments
        if MACHINE != "local":
            # get experiment type and specific directories
            for experiment_type, experiment_dict in interp_experiments.items():
                if self.experiment_to_process in experiment_dict["experiments"]:
                    # take first functional directory 
                    for temp_exp_dir in experiment_dict["paths"]:
                        temp_exp_dir = join(temp_exp_dir, self.experiment_to_process)
                        if os.path.exists(temp_exp_dir):
                            exp_dir = temp_exp_dir
                            break
                    break

        # if local machine or if not exp_dir, get directory from data_paths
        if exp_dir is None: 
            exp_to_interp_path = join(self.exp_to_interp_root, self.experiment_to_process)
            if os.path.exists(exp_to_interp_path):
                exp_dir = exp_to_interp_path

        # define if network is in GHOST format
        self.reading_ghost = check_for_ghost(self.network_to_interpolate_against)
        
        # get relevant observational file
        # GHOST
        if self.reading_ghost:
            self.obs_file = glob.glob(self.ghost_root + '/{}/{}/{}/{}/{}_{}*.nc'\
                                      .format(self.network_to_interpolate_against, self.ghost_version,
                                              self.temporal_resolution_to_output, self.original_speci_to_process,
                                              self.original_speci_to_process, self.yearmonth))[0]
        # non-GHOST
        else:
            self.obs_file = glob.glob(self.nonghost_root + '/{}/{}/{}/{}_{}*.nc'\
                .format(self.network_to_interpolate_against, self.temporal_resolution_to_output,
                        self.original_speci_to_process, self.original_speci_to_process, self.yearmonth))[0]
           
        # get relevant model files
        if self.ensemble_member:
            all_model_files = np.sort(glob.glob('{}/{}/{}/{}/{}*{}*.nc'\
                                                .format(exp_dir, self.grid_type,
                                                        self.model_temporal_resolution,
                                                        self.speci_to_process, self.speci_to_process, self.yearmonth)))
             
            # drop all analysis files ending with '_an.nc' which are not in ensemble-stats
            all_model_files = [f for f in all_model_files if '_an.nc' not in f] 

            # isolate model files to be only those associated with relevant ensemble option
            self.model_files = np.sort([f for f in all_model_files if '{}-{}_'.format(
                self.speci_to_process,self.ensemble_option) in f])
            
            # if the number of remaining model files is 0, this is because the files
            # do not contain an ensemble member number, therefore take all model files
            if len(self.model_files) == 0:
                self.model_files = all_model_files
     
        else:            
            self.model_files = np.sort(glob.glob('{}/{}/{}/ensemble-stats/{}_{}/{}*{}*{}.nc'\
                                                .format(exp_dir, self.grid_type,
                                                        self.model_temporal_resolution,
                                                        self.speci_to_process, self.ensemble_option, 
                                                        self.speci_to_process, self.yearmonth, self.ensemble_option)))
            
        # if have hour in model fname, and it is not 0 then insert previous file also to get times from previous yearmonth that are offset
        if not self.ensemble_member:
            first_file_date = self.model_files[0].replace('_' + self.ensemble_option, '').split('_')[-1][:-3]
        else:
            first_file_date = self.model_files[0].split('_')[-1][:-3]  
        if len(first_file_date) == 10:
            model_start_hour = int(first_file_date[-2:])
            if model_start_hour != 0:
                # get previous month
                prev_yearmonth = (datetime.datetime.strptime(self.yearmonth, '%Y%m') - relativedelta.relativedelta(months=1)).strftime('%Y%m')
                # get previous month files
                if self.ensemble_member:
                    prev_month_files = np.sort(glob.glob('{}/{}/{}/{}/{}*{}*.nc'\
                                                        .format(exp_dir, self.grid_type,
                                                                self.model_temporal_resolution,
                                                                self.speci_to_process, self.speci_to_process, prev_yearmonth)))
                else:
                    prev_month_files_final = np.sort(glob.glob('{}/{}/{}/ensemble-stats/{}_{}/{}*{}*{}.nc'\
                                                        .format(exp_dir, self.grid_type,
                                                                self.model_temporal_resolution,
                                                                self.speci_to_process, self.ensemble_option, 
                                                                self.speci_to_process, prev_yearmonth, self.ensemble_option)))

                if self.ensemble_member:
                    # drop all analysis files ending with '_an.nc' which are not in ensemble-stats
                    prev_month_files = [f for f in prev_month_files if '_an.nc' not in f] 

                    # isolate model files to be only those associated with relevant ensemble option
                    prev_month_files_final = np.sort([f for f in prev_month_files if '{}-{}_'.format(
                        self.speci_to_process,self.ensemble_option) in f])
                    
                    # if the number of remaining model files is 0, this is because the files
                    # do not contain an ensemble member number, therefore take all model files
                    if len(prev_month_files_final) == 0:
                        prev_month_files_final = prev_month_files

                # if have previous month files, insert last one as first file in all_model_files
                if len(prev_month_files_final) > 0:
                    self.model_files = np.insert(self.model_files, 0, prev_month_files_final[-1])

        # if have no valid model files, exit
        if len(self.model_files) == 0:
            self.log_file_str += 'No valid model files in {}. Skipping month.'.format(self.yearmonth)
            create_output_logfile(1, self.log_file_str)

    def get_model_information(self):
        """ Take first valid model file in month and get grid dimension/coordinate information.
            Put initial object read in a  try/except to handle reading of corrupted files.
            Iterate through files until have read a valid file.
            If do not read a valid file, skip month.
        """

        for model_file_ii, model_file in enumerate(self.model_files):
            
            try:
                # load instance of model file netCDF
                mod_nc_root = Dataset(model_file)

                # get all variable names in file
                mod_nc_varnames = list(mod_nc_root.variables.keys())

                # get instance of species variable
                mod_speci_obj = mod_nc_root[self.speci_to_process]
    
                # get species units
                self.mod_speci_units = mod_speci_obj.units
                 
                # get model grid type
                self.mod_grid_type = mod_speci_obj.grid_mapping
        
                # get x/y grid dimension variable names
                grid_dimensions = mod_speci_obj.dimensions

            except Exception as e:
                # in case of exception close model netCDF root
                mod_nc_root.close()
                
                # if have got to last file of month and that is corrupted, return from function
                if model_file_ii == (len(self.model_files)-1):
                    self.log_file_str += 'All model files corrupted in {}. Skipping month. \n'.format(self.yearmonth)
                    self.log_file_str += 'Error: {}'.format(e)
                    create_output_logfile(1, self.log_file_str)
                # else, continue to next file in month
                else:
                    continue 

            # get indivudual dimension variable names  
            # standard (no vertical dimension)
            self.have_vertical_dimension = False
            self.have_bin_dimension = False
            if len(mod_speci_obj.shape) == 3:
                self.x_varname = mod_speci_obj.dimensions[2]
                self.y_varname = mod_speci_obj.dimensions[1]
            # mapped size distribution variable, with bin dimension
            elif (len(mod_speci_obj.shape) == 4) & ('vconcaerobin' in self.original_speci_to_process):
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
                if direction == 'up':
                    if self.interp_reverse_vertical_orientation:
                        self.z_index = np.argmax(mod_vert_obj[:])
                    else:
                        self.z_index = np.argmin(mod_vert_obj[:])
                # if direction == 'down', surface index is location of maximum value in z var
                elif direction == 'down':
                    if self.interp_reverse_vertical_orientation:
                        self.z_index = np.argmin(mod_vert_obj[:])
                    else:
                        self.z_index = np.argmax(mod_vert_obj[:])
                # if cannot determine a surface index, terminate process
                else: 
                    self.log_file_str += 'Cannot determine surface index in vertical dimension. Terminating process.'
                    create_output_logfile(1, self.log_file_str)

            # check if species grid dimensions are named correctly, and in correct BSC standard order
            # if not terminate process
            # this is done by checking the variable names of the x, y (and z if required) dimensions
            # X dimension is valid if 'lon' is contained within name, or is == 'x'
            if ('lon' not in self.x_varname) & (self.x_varname != 'x'):
                self.log_file_str += 'X dimension incorrectly named. Terminating process.'
                create_output_logfile(1, self.log_file_str)
            # Y dimension is valid if 'lat' is contained within name, or is == 'y'
            if ('lat' not in self.y_varname) & (self.y_varname != 'y'):
                self.log_file_str += 'Y dimension incorrectly named. Terminating process.'
                create_output_logfile(1, self.log_file_str)
            # Z dimension is valid if == 'z' or 'lev' or 'alt' or 'height'
            if self.have_vertical_dimension:
                if ((self.z_varname != 'lev') & (self.z_varname != 'z') & (self.z_varname != 'alt')
                    & (self.z_varname != 'height')):
                    self.log_file_str += 'Z dimension incorrectly named. Terminating process.'
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
            grid_centre_coordinates = mod_speci_obj.coordinates.split(' ')
            lon_centre_varname = grid_centre_coordinates[1] 
            lat_centre_varname = grid_centre_coordinates[0]

            # check if species grid centre coordinate are named correctly, and in correct BSC standard order
            # if not terminate process
            # longitude coordinate is valid if 'lon' is contained within name
            if ('lon' not in lon_centre_varname):
                self.log_file_str += 'Longitude grid centre coordinate incorrectly named. Terminating process.'
                create_output_logfile(1, self.log_file_str)
            # latitude coordinate is valid if 'lat' is contained within name
            if ('lat' not in lat_centre_varname):
                self.log_file_str += 'Latitude grid centre coordinate incorrectly named. Terminating process.'
                create_output_logfile(1, self.log_file_str)

            # get longitude and latitude grid centre values
            self.mod_lons_centre = np.float32(mod_nc_root[lon_centre_varname][:])
            self.mod_lats_centre = np.float32(mod_nc_root[lat_centre_varname][:])
            
            # check if there are coordinates to remap
            lons_to_remap = self.mod_lons_centre[np.where(self.mod_lons_centre > 180.0)]
            lats_to_remap = self.mod_lats_centre[np.where(self.mod_lats_centre > 90.0)]
            self.coords_remapping = False
            if len(lons_to_remap) > 0 or len(lats_to_remap) > 0:
                if self.mod_grid_type == 'crs':
                    # correct model grid longitude/latitude centre variables to be between -180:180 and -90:90 respectively
                    self.mod_lons_centre[np.where(self.mod_lons_centre > 180.0)] = self.mod_lons_centre[np.where(self.mod_lons_centre > 180.0)] - 360.0
                    self.mod_lats_centre[np.where(self.mod_lats_centre > 90.0)] = self.mod_lats_centre[np.where(self.mod_lats_centre > 90.0)] - 180.0
                    
                    # sort coordinates
                    self.mod_lons_centre = sorted(self.mod_lons_centre)
                    self.mod_lats_centre = sorted(self.mod_lats_centre) 

                    # set variable for values reassignation later
                    self.coords_remapping = True    
                else: 
                    self.log_file_str += 'Cannot handle grid of type: {} with these coordinates. Please remap. Terminating process'.format(self.mod_grid_type)
                    create_output_logfile(1, self.log_file_str)
            
            # close model netCDF root
            mod_nc_root.close()

            # get the last valid model file
            self.valid_model_file = copy.deepcopy(model_file)

            # break out of for loop, now that have read a valid model file in the month
            break

        # correct units in case of having bin dimension
        if self.have_bin_dimension:
            self.mod_speci_units = self.standard_parameter_speci['standard_units']

    def create_grid_domain_edge_polygon(self):
        """ Create grid domain edge polygon from model netCDF file.
            This is handled differently for regular grids (i.e. following lines of longitude/latitude), 
            and non-regular grids (e.g. lambert-conformal).
        """
        # load instance of model file netCDF
        mod_nc_root = Dataset(self.valid_model_file)

        # if grid type == 'crs, then is regular grid
        if self.mod_grid_type == 'crs':
            
            # set general grid type flag
            self.general_grid_type = 'regular'

            # set longitude/latitude centre coordinates
            x_centre = self.mod_lons_centre
            y_centre = self.mod_lats_centre

        # here handle non-regular grids
        # currently only rotated pole and lambert conformal are supported
        elif self.mod_grid_type in ['rotated_pole', 'Lambert_conformal', 'Lambert_Conformal']:

            # set general grid type flag
            self.general_grid_type = 'non-regular'

            # get instance of variable which defines key variables associated with the non-regular grid
            non_regular_grid_type_obj = mod_nc_root[self.mod_grid_type]

            # get non-regular grid-specific variables used for defining coordinate reference system
            if self.mod_grid_type == 'rotated_pole':
                pole_longitude = np.float32(non_regular_grid_type_obj.grid_north_pole_longitude)
                pole_latitude = np.float32(non_regular_grid_type_obj.grid_north_pole_latitude)
                
                # add workaround to fix interpolation in the southern hemisphere
                # if the centre gridbox value (average of 2, if even number) of the geographic latitude variable is <0:
                # then add 'central_rotated_longitude = 180.0'
                lat_centre_i = findMiddle(self.mod_lats_centre.shape[0])
                lon_centre_i = findMiddle(self.mod_lats_centre.shape[1])
                centre_lat = np.average(self.mod_lats_centre[lat_centre_i,lon_centre_i])
                if centre_lat < 0.0:
                    central_rotated_longitude = 180.0
                else:
                    central_rotated_longitude = 0.0

            elif self.mod_grid_type in ['Lambert_conformal', 'Lambert_Conformal']:
                standard_parallels = np.float32(non_regular_grid_type_obj.standard_parallel.split(','))
                central_longitude = np.float32(non_regular_grid_type_obj.longitude_of_central_meridian)
                central_latitude  = np.float32(non_regular_grid_type_obj.latitude_of_projection_origin)

            # read in non-regular grid x/y grid centre coordinates 
            x_centre = np.float32(mod_nc_root[self.x_varname][:])
            y_centre = np.float32(mod_nc_root[self.y_varname][:])

        # the grid type cannot be handled, therefore terminate process
        else:
            self.log_file_str += 'Cannot handle grid of type: {}. Terminating process'.format(self.mod_grid_type)
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
        x_edge = x_centre - (x_res/2.)
        if x_right_ind == -1:
            x_edge = np.append(x_edge, x_edge[x_right_ind]+x_res)
        else:
            x_edge = np.insert(x_edge, 0, x_edge[x_right_ind]+x_res)
        y_edge = y_centre - (y_res/2.)
        if y_top_ind == -1:
            y_edge = np.append(y_edge, y_edge[y_top_ind]+y_res)
        else:
            y_edge = np.insert(y_edge, 0, y_edge[y_top_ind]+y_res) 

        # get x/y coordinates around grid edges
        left_edge_x = np.full(len(y_edge), x_edge[x_left_ind])
        right_edge_x = np.full(len(y_edge), x_edge[x_right_ind])
        if x_left_ind == 0:
            top_edge_x = x_edge[1:-1]
            bottom_edge_x = x_edge[-2::-1]
        else:
            top_edge_x = x_edge[-2::-1]
            bottom_edge_x = x_edge[1:-1]
        x_grid_edge = np.concatenate((left_edge_x,top_edge_x,right_edge_x,bottom_edge_x))

        bottom_edge_y = np.full(len(x_edge)-1, y_edge[y_bottom_ind])
        top_edge_y = np.full(len(x_edge)-2, y_edge[y_top_ind])
        if y_bottom_ind == 0:
            left_edge_y = y_edge[:]
            right_edge_y = y_edge[::-1]
        else:
            left_edge_y = y_edge[::-1]
            right_edge_y = y_edge[:]
        y_grid_edge = np.concatenate((left_edge_y,top_edge_y,right_edge_y,bottom_edge_y))

        # regular grid type?
        if self.general_grid_type == 'regular':
            # stack longitude/latitude bounding edge coordinate pairs 
            self.model_grid_outline = np.vstack((x_grid_edge,y_grid_edge)).T
        
            # get 2D mesh of model longitude/latitude gridcell centres
            self.mod_lons_centre, self.mod_lats_centre = np.meshgrid(self.mod_lons_centre, self.mod_lats_centre)
            
        # non-regular grid type?
        elif self.general_grid_type == 'non-regular':
            # create cartopy coordinate reference system for the specific type non-standard grid, on WGS84 ellipsoid
            if self.mod_grid_type == 'rotated_pole':
                non_regular_grid_crs = ccrs.RotatedPole(pole_longitude=pole_longitude, pole_latitude=pole_latitude, 
                                                        central_rotated_longitude=central_rotated_longitude)
            elif self.mod_grid_type in ['Lambert_conformal', 'Lambert_Conformal']:
                non_regular_grid_crs = ccrs.LambertConformal(central_longitude=central_longitude, 
                                                             central_latitude=central_latitude, 
                                                             standard_parallels=standard_parallels)
            
            # define a regular gridded coordinate reference system (Plate Carree), on WGS84 ellipsoid
            plate_carree = ccrs.PlateCarree()

            # convert bounding coordinates of box defined in non-regular grid coordinates to regular grid coordinates
            self.model_grid_outline = plate_carree.transform_points(non_regular_grid_crs, x_grid_edge, y_grid_edge)[:,:2]

        # convert longitude/latitude bounding coordinate pairs to a shapely polygon 
        self.model_grid_outline_poly = Polygon(self.model_grid_outline)

        # close model netCDF root - now have all neccessary grid information
        mod_nc_root.close()

    def get_observational_objects(self):
        """ Get necessary observational objects. """

        # get observational file netCDF root
        obs_nc_root = Dataset(self.obs_file)

        # get measured observational variable object 
        obs_measured_var_obj = obs_nc_root[self.original_speci_to_process]

        # get the metadata
        self.obs_units = obs_measured_var_obj.units
        if self.reading_ghost:
            self.obs_long_name = obs_measured_var_obj.long_name
            self.obs_standard_name = obs_measured_var_obj.standard_name

        # station object
        # for GHOST always is "station_reference" 
        # for non-GHOST, try for "station_reference", then "station_code", then "station_name".
        if self.reading_ghost:
            obs_station_reference_obj = obs_nc_root['station_reference']
        else:
            if 'station_reference' in obs_nc_root.variables:
                obs_station_reference_obj = obs_nc_root['station_reference']
            elif 'station_code' in obs_nc_root.variables:
                obs_station_reference_obj = obs_nc_root['station_code']
            elif 'station_name' in obs_nc_root.variables:
                obs_station_reference_obj = obs_nc_root['station_name']

        # get station references atributes
        if self.reading_ghost:
            self.obs_station_reference_units = obs_station_reference_obj.units
            self.obs_station_reference_standard_name = obs_station_reference_obj.standard_name
            self.obs_station_reference_description = obs_station_reference_obj.description
            self.obs_station_reference_long_name = obs_station_reference_obj.long_name

        # lon/lat objects
        if "latitude" in obs_nc_root.variables:
            obs_lon_obj = obs_nc_root['longitude']
            obs_lat_obj = obs_nc_root['latitude']
        else:
            obs_lon_obj = obs_nc_root['lon']
            obs_lat_obj = obs_nc_root['lat']

        # get lon lat atributes
        self.obs_lat_obj_units, self.obs_lon_obj_units = obs_lat_obj.units, obs_lon_obj.units
        if self.reading_ghost:
            self.obs_lat_obj_standard_name, self.obs_lon_obj_standard_name = obs_lat_obj.standard_name, obs_lon_obj.standard_name
            self.obs_lat_obj_long_name, self.obs_lon_obj_long_name = obs_lat_obj.long_name, obs_lon_obj.long_name
            self.obs_lat_obj_description, self.obs_lon_obj_description = obs_lat_obj.description, obs_lon_obj.description

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
                    self.station_references = np.array([''.join(val) for val in self.station_references])
                else:
                    self.station_references = chartostring(self.station_references)
    
        # get indices of all non-NaN stations (can be NaN for some non-GHOST files)
        non_nan_station_indices = np.array([ref_ii for ref_ii, ref in enumerate(self.station_references) if ref.lower() != 'nan'])
        self.station_references = self.station_references[non_nan_station_indices]

        # get valid longitude and latitudes
        self.obs_lons = obs_lon_obj[non_nan_station_indices]
        self.obs_lats = obs_lat_obj[non_nan_station_indices]

        # close the nc file
        obs_nc_root.close()

    def get_conversion_factor(self):
        """ Get conversion factor between observations and experiment data. """

        # get units conversion factor between model and observations (go from model to observational units)    
        obs_speci_units = self.obs_units

        # if units are unitless, then no need for conversion (i.e. conversion factor = 1.0)   
        if (obs_speci_units == 'unitless') or (obs_speci_units == '-') or (obs_speci_units == '1'):
            self.conversion_factor = 1.0
            return

        # unit converter module does not produce conversion factor for temperature, but both observational and model 
        # units should be Kelvin (i.e. conversion factor = 1.0) 
        # if model not in K then return error
        elif obs_speci_units == 'K': 
            if self.mod_speci_units == 'K':
                self.conversion_factor = 1.0
                return
            else:
                self.log_file_str += "Experiment units should be 'K', but are set as '{}'".format(self.mod_speci_units)
                create_output_logfile(1, self.log_file_str)

        # unit converter module does not produce conversion factor for angular degrees, but both observational and model 
        # units should be in angular degrees (i.e. conversion factor = 1.0) 
        # if model not in K then return error
        elif obs_speci_units == 'angular degrees': 
            if ((self.mod_speci_units.lower() == 'angular degrees') or (self.mod_speci_units.lower() == 'degrees') 
                or (self.mod_speci_units.lower() == '°')):
                self.conversion_factor = 1.0
                return
            else:
                self.log_file_str += "Experiment units should be 'angular degrees', but are set as '{}'".format(
                    self.mod_speci_units)
                create_output_logfile(1, self.log_file_str)
        
        # otherwise check if the unit quantities are equal
        conv_obj = unit_converter.convert_units(obs_speci_units,obs_speci_units,1)
        obs_quantity = conv_obj.output_represented_quantity
        conv_obj = unit_converter.convert_units(self.mod_speci_units,self.mod_speci_units,1)
        model_quantity = conv_obj.output_represented_quantity

        # observations and model quantities not equal (convert to observational units, standard_temperature=293.15, 
        # standard_pressure=1013.25)
        if obs_quantity != model_quantity:
            # determine chemical formula of species 
            speci_chemical_formula = self.standard_parameter_speci['chemical_formula']
            # if cannot determine chemical formula of species, then terminate process
            if speci_chemical_formula == '':
                self.log_file_str += 'Cannot determine speci chemical formula needed for unit conversion. Terminating process.'
                create_output_logfile(1, self.log_file_str)            
            input_units ={'temperature':'K', 'pressure':'hPa', 'molar_mass':'kg mol-1', model_quantity:self.mod_speci_units}
            input_values = {'temperature':293.15, 'pressure':1013.25, 'molar_mass':unit_converter.get_molecular_mass(speci_chemical_formula), 
                                                                                                                     model_quantity:1.0}
            conv_obj = unit_converter.convert_units(input_units, obs_speci_units, input_values, 
                                                    conversion_input_quantity=model_quantity)
            self.conversion_factor = conv_obj.conversion_factor
        else:
            conv_obj = unit_converter.convert_units(self.mod_speci_units, obs_speci_units, 1.0) 
            self.conversion_factor = conv_obj.conversion_factor

    def get_monthly_model_data(self):
        """ Read all relevant model data in yearmonth into memory. """

        # temporal resolution to output is coarser than model resolution, therefore will need to modify temporal 
        # resolution and perform resampling?
        resampling = False
        if self.temporal_resolution_to_output != self.model_temporal_resolution: 
            resampling = True

        # get number of days in month processing
        self.days_in_month = monthrange(int(self.year),int(self.month))[1]
        
        # create monthly time/dy variables
        start_month_dt = datetime.datetime(year=int(self.year), month=int(self.month), day=1, hour=0, minute=0)
        end_month_dt = start_month_dt + relativedelta.relativedelta(months=1)
        if self.temporal_resolution_to_output in ['hourly', 'hourly_instantaneous']:
            self.yearmonth_time = np.arange(0,self.days_in_month*24.0)
            if Version(pd.__version__) >= Version("2.2"):
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='h', inclusive='left')
            elif Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='H', inclusive='left')
            else:
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='H', closed='left')
            self.descriptive_temporal_resolution = 'hours'
            self.temporal_resolution_to_output_code = 'H'
        elif self.temporal_resolution_to_output in ['3hourly', '3hourly_instantaneous']:
            self.yearmonth_time = np.arange(0,self.days_in_month*24.0,3.0)
            if Version(pd.__version__) >= Version("2.2"):
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='3h', inclusive='left')
            elif Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='3H', inclusive='left')
            else:
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='3H', closed='left')
            self.descriptive_temporal_resolution = 'hours'
            self.temporal_resolution_to_output_code = '3H'
        elif self.temporal_resolution_to_output in ['6hourly', '6hourly_instantaneous']:
            self.yearmonth_time = np.arange(0,self.days_in_month*24.0,6.0)
            if Version(pd.__version__) >= Version("2.2"):
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='6h', inclusive='left')
            elif Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='6H', inclusive='left')
            else:
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='6H', closed='left')
            self.descriptive_temporal_resolution = 'hours'
            self.temporal_resolution_to_output_code = '6H'
        elif self.temporal_resolution_to_output == 'daily':
            self.yearmonth_time = np.arange(0,self.days_in_month)
            if Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='D', inclusive='left')
            else:
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='D', closed='left')
            self.descriptive_temporal_resolution = 'days'
            self.temporal_resolution_to_output_code = 'D'
        elif self.temporal_resolution_to_output == 'monthly':
            self.yearmonth_time = np.arange(0,1)
            if Version(pd.__version__) >= Version("1.4"):
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='MS', inclusive='left')
            else:
                self.yearmonth_dt = pd.date_range(start_month_dt, end_month_dt, freq='MS', closed='left')
            self.descriptive_temporal_resolution = 'months'
            self.temporal_resolution_to_output_code = 'MS'

        # if have forecast data, then get the max number of complete forecast days across model files
        if self.forecast:

            # ititalise max forecast days as 0
            max_forecast_days = 0

            # iterate and read chunked model files
            for model_ii, model_file in enumerate(self.model_files):
                # read in chunked netcdf file
                mod_nc_root = Dataset(model_file)

                # check if have time dimension in daily file, if do not, do not process file
                if 'time' not in list(mod_nc_root.dimensions.keys()):
                    mod_nc_root.close()
                    continue

                # get date from filename
                if not self.ensemble_member:
                    file_date = model_file.replace('_' + self.ensemble_option, '').split('_')[-1][:-3]
                else:
                    file_date = model_file.split('_')[-1][:-3]    

                # get file chunk type
                if len(file_date) == 6:
                    chunk_type = 'monthly'
                elif len(file_date) in [8,10]:
                    chunk_type = 'daily'
                else:
                    mod_nc_root.close()
                    continue

                # cannot have chunked monthly file which has forecast data, so throw an error
                if chunk_type == 'monthly':
                    log_file_str += 'File {} is monthly, for which forecast data cannot be processed. Terminating process.'.format(model_file,self.forecast_day)
                    create_output_logfile(1)
                # cannot have monthly resolution file which has forecast data, so throw an error
                elif self.model_temporal_resolution == 'monthly':
                    log_file_str += 'File {} resolution is monthly, for which forecast data cannot be processed. Terminating process.'.format(model_file,self.forecast_day)
                    create_output_logfile(1)

                # get number of timesteps in file
                n_timesteps = mod_nc_root.dimensions['time'].size

                # get number of complete forecast days
                if self.model_temporal_resolution in ['hourly', 'hourly_instantaneous']:
                    forecast_days = int(n_timesteps / 24.0)
                elif self.model_temporal_resolution in ['3hourly', '3hourly_instantaneous']:
                    forecast_days = int(n_timesteps / (24.0/3.0))
                elif self.model_temporal_resolution in ['6hourly', '6hourly_instantaneous']:
                    forecast_days = int(n_timesteps / (24.0/6.0))
                elif self.model_temporal_resolution == 'daily':
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

        # create array for storing daily monthly model data
        # if have forecast data, then have an extra dimension for the number of forecast days
        if self.forecast:
            self.monthly_model_data = np.full((len(self.yearmonth_time), self.forecast_days, self.y_N, self.x_N), np.NaN, dtype=np.float32)
        else:
            self.monthly_model_data = np.full((len(self.yearmonth_time), self.y_N, self.x_N), np.NaN, dtype=np.float32)

        # if have mapped size distribution variable: 
        # -- get aeronet bin radius for original variable 
        # -- get index of model bin to use for transformation of bin types (and rmin/rmax)
        # -- get transform factor to go between model bin and aeronet bin
        if self.have_bin_dimension:
            aeronet_bin_radius = get_aeronet_bin_radius_from_bin_variable(self.original_speci_to_process)
            bin_index, rmin, rmax, rho_bin = get_aeronet_model_bin(self.model_name, aeronet_bin_radius)
            bin_transform_factor = get_model_to_aeronet_bin_transform_factor(self.model_name, rmin, rmax)
              
        # iterate and read chunked model files
        failed_files = 0
        for model_ii, model_file in enumerate(self.model_files):

            # put model file read in try/except to catch corrupt model files
            try:
                # read in chunked netcdf file
                mod_nc_root = Dataset(model_file)
            
                # check if have time dimension in daily file, if do not, do not process file
                if 'time' not in list(mod_nc_root.dimensions.keys()):
                    mod_nc_root.close()
                    self.log_file_str += 'File {} is corrupt. Skipping.\n'.format(model_file)
                    failed_files += 1
                    continue 
                
                # get date from filename
                if not self.ensemble_member:
                    file_date = model_file.replace('_' + self.ensemble_option, '').split('_')[-1][:-3]
                else:
                    file_date = model_file.split('_')[-1][:-3]                
                
                # get file time (handle monthly resolution data differently to hourly/daily
                # as num2date does not support 'months since' units)
                file_time = mod_nc_root['time'][:] 
                time_units = mod_nc_root['time'].units
                time_calendar = mod_nc_root['time'].calendar
                
                # if have >= spinup timesteps than the number of timesteps in the file, then throw error
                self.interp_spinup_timesteps = int(self.interp_spinup_timesteps)
                if self.interp_spinup_timesteps >= len(file_time):
                    self.log_file_str += "File {} has <= number of timesteps of 'interp_spinup_timesteps':{}. Terminating process.".format(model_file, self.interp_spinup_timesteps)
                    create_output_logfile(1, self.log_file_str)

                # convert file time to datetime
                if 'months' in time_units:
                    monthly_start_date = time_units.split(' ')[2]
                    file_time_dt = pd.date_range(start=monthly_start_date, periods=1, freq='MS')
                else:
                    file_time_dt = num2date(file_time, units=time_units, calendar=time_calendar)

                    # convert to pandas datetime
                    if Version(cftime.__version__) <= Version("1.0.3.4"):
                        # remove microseconds
                        file_time_dt = pd.to_datetime([t.replace(microsecond=0) for t in file_time_dt])
                    else:
                        # bug fix for newer versions of cftime
                        file_time_dt = file_time_dt.astype('datetime64[s]')
                        file_time_dt = pd.to_datetime([t for t in file_time_dt])

                # get file start and end datetime
                if len(file_date) == 6:
                    chunk_type = 'monthly'
                    start_file_dt = datetime.datetime(year=int(file_date[:4]), month=int(file_date[4:6]), day=1, hour=0,
                                                      minute=0)
                    end_file_dt = start_file_dt + relativedelta.relativedelta(months=1)
                elif len(file_date) == 8:
                    chunk_type = 'daily'
                    start_file_dt = datetime.datetime(year=int(file_date[:4]), month=int(file_date[4:6]), 
                                                      day=int(file_date[6:8]), hour=0, minute=0)
                    end_file_dt = start_file_dt + datetime.timedelta(days=1)
                elif len(file_date) == 10:
                    chunk_type = 'daily'
                    start_file_dt = datetime.datetime(year=int(file_date[:4]), month=int(file_date[4:6]), 
                                                      day=int(file_date[6:8]), hour=int(file_date[8:10]), minute=0)
                    end_file_dt = start_file_dt + datetime.timedelta(days=1)
                else:
                    mod_nc_root.close()
                    self.log_file_str += 'Resolution could not be detected in {}, check the date in the filename as now it shows as "{}".\n'.format(model_file, file_date)
                    failed_files += 1
                    continue

                # for forceast data, get valid indices of file time per forecast day
                if self.forecast:

                    # get indices of file time within each forecast day (excluding spinup timsesteps)
                    valid_file_time_inds = []
                    for forecast_day in range(self.forecast_days):
                        forecast_start_dt = start_file_dt + datetime.timedelta(days=forecast_day)
                        forecast_end_dt = start_file_dt + datetime.timedelta(days=forecast_day+1)
                        if forecast_day == 0:
                            valid_inds = np.where((file_time_dt >= forecast_start_dt) & (file_time_dt < forecast_end_dt) & (file_time_dt >= start_month_dt) & (file_time_dt < end_month_dt))[0]
                        else:
                            valid_inds = np.where((file_time_dt >= forecast_start_dt) & (file_time_dt < forecast_end_dt))[0]
                        valid_file_time_inds.append(np.array([valid_ind for valid_ind in valid_inds if valid_ind >= self.interp_spinup_timesteps]), dtype=np.int32)

                # for non-forecast runs, get all timesteps (excluding spinup timesteps)
                else:
                    # get indices of file time inside yearmonth
                    valid_inds = np.where((file_time_dt >= start_month_dt) & (file_time_dt < end_month_dt))[0]
                    valid_file_time_inds = np.array([valid_ind for valid_ind in valid_inds if valid_ind >= self.interp_spinup_timesteps], dtype=np.int32)

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
                        read_data = mod_nc_root[self.speci_to_process][valid_file_time_inds[forecast_day_ii],bin_index,:,:]
                        
                        # convert model units from kg m-2 to um3/um2
                        read_data = 1e6 * read_data / rho_bin
                        
                        # transform model bin data to aeronet bin
                        read_data = read_data * bin_transform_factor
                    # has vertical dimension 
                    elif self.have_vertical_dimension:
                        read_data = mod_nc_root[self.speci_to_process][valid_file_time_inds[forecast_day_ii],self.z_index,:,:] 
                    # has no vertical dimension?
                    else:
                        read_data = mod_nc_root[self.speci_to_process][valid_file_time_inds[forecast_day_ii],:,:]

                    # convert model units to observational units
                    read_data = read_data * self.conversion_factor
                    
                    # set any model values outside GHOST extreme limits for variable to be NaN
                    read_data[(read_data < self.GHOST_speci_lower_limit) | (read_data > self.GHOST_speci_upper_limit)] = np.NaN 

                    # create xarray for resampling
                    xr_data = xr.DataArray(dims=("time","latitude","longitude"),
                                           data=read_data,
                                           coords=dict(time=valid_file_time_dt, latitude=self.y, longitude=self.x))
                    
                    # reassign values to correct coordinates if model centre coordinates have been remapped 
                    # (only in the case of regular grids)
                    if self.coords_remapping:
                        # correct 1D arrays
                        self.x = [x if x <= 180 else x - 360 for x in self.x]
                        self.y = [y if y <= 90 else y - 180 for y in self.y]

                        # assign and sort coordinates
                        xr_data = xr_data.assign_coords(longitude=self.x).sortby('longitude')
                        xr_data = xr_data.assign_coords(latitude=self.y).sortby('latitude')

                    # do resampling (taking mean) to coarser temporal resolution if neccessary
                    if resampling:
                        xr_data = xr_data.resample(time=self.temporal_resolution_to_output_code).mean()     
                    
                    # get indices in yearmonth to fill with read data
                    inds_to_fill = np.isin(self.yearmonth_dt, xr_data.time.values)
                    if not any(inds_to_fill):
                        if not bad_time_format:
                            self.log_file_str += 'Time in model dataset {} is not in standard format: \n {}'.format(model_file, xr_data.time.values)
                            failed_files += 1
                            bad_time_format = True
                        continue
                    
                    # fill in data array
                    if self.forecast:
                        self.monthly_model_data[inds_to_fill,forecast_day_ii,:,:] = xr_data.values
                    else:
                        self.monthly_model_data[inds_to_fill,:,:] = xr_data.values

                # close model netCDF root
                mod_nc_root.close()

            except Exception as e:
                mod_nc_root.close()
                self.log_file_str += 'File {} is corrupt. Skipping.\n{}'.format(model_file, traceback.format_exc())
                failed_files += 1
        
        if failed_files == len(self.model_files):
            self.log_file_str += 'No model dataset could be interpolated.'
            create_output_logfile(1, self.log_file_str)
        
    def n_nearest_neighbour_inverse_distance_weights(self):
        """ Calculate N nearest neighbour inverse distance weights (and indices) of model gridcells centres 
            from an array of geographic observational station coordinates. Both observational and model geographic 
            longitude/latitude coordinates are first converted to cartesian ECEF (Earth Centred, Earth Fixed) 
            coordinates, before calculating distances. Weights returned for obervational stations not contained 
            within model grid extents are all zero.
        """

        # for each pair of observational station geographic coordinates, test if station is inside model grid
        obs_inside_model_grid = np.array([self.model_grid_outline_poly.contains(Point(obs_lon, obs_lat)) 
                                          for obs_lon, obs_lat in zip(self.obs_lons,self.obs_lats)]) 

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
            obs_x, obs_y, obs_z = transformer.transform(self.obs_lons, self.obs_lats, obs_alts, radians=False)
            mod_x, mod_y, mod_z = transformer.transform(flat_mod_lons_centre, flat_mod_lats_centre, flat_mod_alts, 
                                                        radians=False)
        else:
            lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
            ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
            obs_x, obs_y, obs_z = pyproj.transform(lla, ecef, self.obs_lons, self.obs_lats, obs_alts, radians=False)
            mod_x, mod_y, mod_z = pyproj.transform(lla, ecef, flat_mod_lons_centre, flat_mod_lats_centre, flat_mod_alts, 
                                                radians=False)

        # stack converted cartesian coordinates for preparation of calculation of nearest neighbour distances 
        # using Scipy cKDTree
        obs_lonlat = np.column_stack((obs_x,obs_y,obs_z))
        mod_lonlat = np.column_stack((mod_x,mod_y,mod_z))
    
        # generate cKDtree
        tree = cKDTree(mod_lonlat)
    
        # get n-neighbour nearest distances/indices (ravel form) of model gridcell centres from each observational station  
        if Version(scipy.__version__) < Version("1.6.0"):
            dists,idx = tree.query(obs_lonlat,k=int(self.interp_n_neighbours))
        else:
            dists,idx = tree.query(obs_lonlat,k=int(self.interp_n_neighbours),workers=-1)
        
        # for n neighbours == 1, do reshaping of array so doesn't break
        if len(dists.shape) == 1:
            dists = np.reshape(dists,(dists.shape[0],1))
            idx = np.reshape(idx,(idx.shape[0],1))
            
        # if any dists are < 1m, then set as 1m (i.e. very small distance), to avoid issue calculating weight
        dists[dists < 1.0] = 1.0

        # get nearest neighbour indices
        self.nearest_neighbour_inds = np.column_stack(np.unravel_index(idx, self.mod_lons_centre.shape))

        # take the reciprocals of the nearest neighbours distances from the observational points
        self.inverse_dists = np.reciprocal(dists)

        # set reciprocal distances for all observational points outside model grid extent to be 0
        self.inverse_dists[~obs_inside_model_grid,:] = 0.0

    def write_netCDF(self):
        """ Write netCDF, returning interpolated model data to observational surface stations 
            for yearmonth.
        """

        # set output directory where observational interpolated monthly model netcdf will be saved
        # GHOST
        if self.reading_ghost:
            network_name = self.network_to_interpolate_against
        # non-GHOST
        else:
            # as it appears in PRV (e.g. nasa-aeronet/oneill_v3-lev15 -> nasa-aeronet-oneill_v3-lev15)
            network_name = self.network_to_interpolate_against.replace('/', '-')
        output_dir = join(self.exp_root, self.ghost_version, self.prov_exp_code, 
                                  self.temporal_resolution_to_output, self.original_speci_to_process, network_name)

        # check if need to create any directories in path 
        check_directory_existence(output_dir, self.exp_root)

        # create netCDF dataset
        netCDF_fname = '{}/{}_{}.nc'.format(output_dir, self.original_speci_to_process, self.yearmonth)
        root_grp = Dataset(netCDF_fname, 'w', format="NETCDF4") 

        # auto mask arrays
        root_grp.set_auto_mask(True)	

        # file contents
        msg = 'Inverse distance weighting ({} neighbours) interpolated '.format(self.interp_n_neighbours)
        msg += '{} experiment data for the component {} '.format(self.experiment_to_process, 
                                                                self.original_speci_to_process)
        msg += 'with reference to the measurement stations in the '
        msg += '{} network '.format(self.network_to_interpolate_against)
        msg += 'in {}-{}.'.format(self.year, self.month)
        root_grp.title = msg
        root_grp.institution = 'Barcelona Supercomputing Center'
        root_grp.source = 'Experiment {}'.format(self.experiment_to_process)
        root_grp.creator_name = 'Dene R. Bowdalo'
        root_grp.creator_email = 'dene.bowdalo@bsc.es'

        # iterate through each station to process
        try:
            for ii, station_reference in enumerate(self.station_references):
                if ii == 0:
                    # netcdf dimensions
                    root_grp.createDimension('station', len(self.obs_lons))
                    root_grp.createDimension('time', len(self.yearmonth_time))
                    root_grp.createDimension('model_longitude', self.x_N)
                    root_grp.createDimension('model_latitude', self.y_N)
                    root_grp.createDimension('grid_edge', self.model_grid_outline.shape[0])
                    # if have forecast data then create forecast data dimension
                    if self.forecast:
                        root_grp.createDimension('forecast_day', self.forecast_days)

                    # create time variable
                    time_var = root_grp.createVariable('time', 'u4', ('time'))
                    time_var.standard_name = 'time'
                    time_var.long_name = 'time'
                    time_var.units = '{} since {}-{}-01 00:00:00'.format(self.descriptive_temporal_resolution, self.year, 
                                                                        self.month)
                    msg = 'Time in {} since {}-{}-01 00:00 UTC. '.format(self.descriptive_temporal_resolution, self.year, 
                                                                        self.month)
                    msg += 'Time given refers to the start of the time window the measurement '
                    msg += 'is representative of (temporal resolution)'
                    time_var.description = msg
                    time_var.axis = 'T'
                    time_var.calendar = 'standard'
                    time_var.tz = 'UTC'

                    # create observational equivalent station reference variable
                    station_reference_var = root_grp.createVariable('station_reference', str, ('station'))
                    if self.reading_ghost:
                        station_reference_var.standard_name = self.obs_station_reference_standard_name
                        station_reference_var.long_name = self.obs_station_reference_long_name
                        station_reference_var.units = self.obs_station_reference_units
                        station_reference_var.description = self.obs_station_reference_description

                    # create observational equivalent longitude/latitude variables
                    longitude_var = root_grp.createVariable('longitude', 'f8', ('station'))
                    if self.reading_ghost:
                        longitude_var.standard_name = self.obs_lon_obj_standard_name
                        longitude_var.units = self.obs_lon_obj_units
                        longitude_var.long_name = self.obs_lon_obj_long_name
                        longitude_var.description = self.obs_lon_obj_description
                    longitude_var.axis = 'X'

                    latitude_var = root_grp.createVariable('latitude', 'f8', ('station'))
                    if self.reading_ghost:
                        latitude_var.standard_name = self.obs_lat_obj_standard_name
                        latitude_var.units = self.obs_lat_obj_units
                        latitude_var.long_name = self.obs_lat_obj_long_name
                        latitude_var.description = self.obs_lat_obj_description
                    latitude_var.axis = 'Y'

                    # create 2D meshed longitude/latitude gridcell centre variables
                    model_centre_longitude_var = root_grp.createVariable('model_centre_longitude', 'f8', 
                                                                        ('model_latitude','model_longitude'))
                    model_centre_longitude_var.standard_name = 'model centre longitude'
                    model_centre_longitude_var.long_name = 'model centre longitude'
                    model_centre_longitude_var.units = self.obs_lon_obj_units
                    msg = '2D meshed grid centre longitudes with '
                    msg += '{} longitudes in {} bands of latitude'.format(self.x_N, self.y_N)
                    model_centre_longitude_var.description = msg
                    model_centre_longitude_var.axis = 'X'            
        
                    model_centre_latitude_var = root_grp.createVariable('model_centre_latitude', 'f8', 
                                                                        ('model_latitude','model_longitude'))
                    model_centre_latitude_var.standard_name = 'model centre latitude'
                    model_centre_latitude_var.long_name = 'model centre latitude'
                    model_centre_latitude_var.units = self.obs_lat_obj_units
                    msg = '2D meshed grid centre longitudes with '
                    msg += '{} longitudes in {} bands of latitude'.format(self.y_N, self.x_N)
                    model_centre_latitude_var.description = msg
                    model_centre_latitude_var.axis = 'Y'

                    # create grid domain edge longitude/latitude variables
                    grid_edge_longitude_var = root_grp.createVariable('grid_edge_longitude', 'f8', ('grid_edge'))   
                    grid_edge_longitude_var.standard_name = 'grid edge longitude'
                    grid_edge_longitude_var.long_name = 'grid edge longitude'
                    grid_edge_longitude_var.units = self.obs_lon_obj_units
                    msg = 'Longitude coordinate along edge of grid domain '
                    msg += '(going clockwise around grid boundary from bottom-left corner).'
                    grid_edge_longitude_var.description = msg
                    grid_edge_longitude_var.axis = 'X'

                    grid_edge_latitude_var = root_grp.createVariable('grid_edge_latitude', 'f8', ('grid_edge'))
                    grid_edge_latitude_var.standard_name = 'grid edge latitude'
                    grid_edge_latitude_var.long_name = 'grid edge latitude'
                    grid_edge_latitude_var.units = self.obs_lat_obj_units
                    msg = 'Latitude coordinate along edge of grid domain '
                    msg += '(going clockwise around grid boundary from bottom-left corner).'
                    grid_edge_latitude_var.description = msg
                    grid_edge_latitude_var.axis = 'Y'

                    # create measured variable
                    if self.forecast:
                        measured_var = root_grp.createVariable(self.original_speci_to_process, 'f4', ('station', 'forecast_day', 'time'))
                    else:
                        measured_var = root_grp.createVariable(self.original_speci_to_process, 'f4', ('station', 'time'))
                    # GHOST
                    if self.reading_ghost:
                        measured_var.long_name = self.obs_long_name
                        measured_var.units = self.obs_units
                        measured_var.standard_name = self.obs_standard_name
                        measured_var.description = 'Interpolated value of {} from the experiment {} ' \
                                            'with reference to the measurement stations in the {} network'.format(
                            self.obs_standard_name, self.experiment_to_process, 
                            self.network_to_interpolate_against)
                    # non-GHOST
                    else:
                        measured_var.standard_name = self.original_speci_to_process
                        measured_var.description = 'Interpolated value of {} from the experiment {} with reference to the measurement stations in the {} network'.format(
                            self.original_speci_to_process, self.experiment_to_process,
                            self.network_to_interpolate_against)

                    # write to variables
                    time_var[:] = self.yearmonth_time
                    station_reference_var[:] = self.station_references
                    longitude_var[:] = self.obs_lons
                    latitude_var[:] = self.obs_lats                                                           
                    model_centre_longitude_var[:] = self.mod_lons_centre
                    model_centre_latitude_var[:] = self.mod_lats_centre
                    grid_edge_longitude_var[:] = self.model_grid_outline[:,0]
                    grid_edge_latitude_var[:] = self.model_grid_outline[:,1]

                # iterate through observational stations
                # use calculated interpolated weights per station to model grid to produce model reciprocal output
                # if all station weights are 0, station is outside model grid domain
                # set all values to be NaN
                station_weights = self.inverse_dists[ii,:]
                if np.all(station_weights == 0):
                    interp_vals = np.full(len(self.yearmonth_time), np.NaN, dtype=np.float32)
                else:
                    # get reciprocal model data at N nearest neighbours to observational station 
                    if self.forecast:
                        cut_model_data = self.monthly_model_data[:, :, self.nearest_neighbour_inds[ii,:int(self.interp_n_neighbours)],
                                                                       self.nearest_neighbour_inds[ii,int(self.interp_n_neighbours):]]
                    else:
                        cut_model_data = self.monthly_model_data[:, self.nearest_neighbour_inds[ii,:int(self.interp_n_neighbours)],
                                                                    self.nearest_neighbour_inds[ii,int(self.interp_n_neighbours):]]
                    
                    # create mask where data == NaN or infinite
                    invalid_mask = ~np.isfinite(cut_model_data)
                    
                    # create masked array
                    cut_model_data = np.ma.MaskedArray(cut_model_data, mask=invalid_mask)
                    
                    # interpolate masked array across time dimension using interpolated weights per station
                    interp_vals = np.ma.average(cut_model_data, weights=station_weights, axis=-1)

                # write measured variable 
                measured_var[ii,:] = interp_vals

            # close writing to netCDF
            root_grp.close() 
        
        except Exception as e:
            self.log_file_str += str(interp_vals.shape)
            self.log_file_str += 'File {} could not be written. Error: {}.\n'.format(netCDF_fname, e)
            create_output_logfile(1, self.log_file_str)

        # compress netCDF file
        try:
            compress_process = subprocess.Popen(['ncks', '-O', '--dfl_lvl', '1', netCDF_fname, netCDF_fname], 
                                                stdout=subprocess.PIPE)
            compress_status = compress_process.communicate()[0]
            compress_return_code = compress_process.returncode
        except:
            self.log_file_str += 'NCO could not be found, please install it in your system ' 
            self.log_file_str += 'with sudo apt install nco (Debian/Ubuntu) or brew install nco (macOS).'
    
        # give 770 permissions for file and make owner bsc32 if machine isn't local
        if MACHINE != 'local':
            set_file_permissions_ownership(netCDF_fname)

        # copy file to esarchive (if have access)
        if MACHINE in ['nord3v2', 'nord4']:
            
            # set esarchive output dir
            esarchive_output_dir = '/esarchive/recon/prov_interp/{}'.format('/'.join(netCDF_fname.split('/exp_interp/')[1].split('/')[:-1]))

            # check if need to create any directories in path
            check_directory_existence(esarchive_output_dir, '/esarchive/recon/prov_interp')

            # set esarchive fname
            esarchive_netCDF_fname = '{}/{}'.format(esarchive_output_dir, netCDF_fname.split('/')[-1])
            
            try:
                # copy file (without permissions)
                shutil.copyfile(netCDF_fname, esarchive_netCDF_fname)

                # give 770 permissions for file and make owner bsc32
                set_file_permissions_ownership(esarchive_netCDF_fname)
            except (PermissionError, FileNotFoundError) as e:
                self.log_file_str += 'Interpolated file/s could not be copied to esarchive. Error: {}'.format(e)
                create_output_logfile(1, self.log_file_str)

def create_output_logfile(process_code, log_file_str):
    """ Create a logfile for stating outcome of interpolation job'
        the filename is prefixed with a code referencing the job outcome.
        The process codes are:
        0: Process completed without issue
        1: Caught error in process
        2: Uncaught error in process

        :param process_code: interpolation outcome code
        :type process_code: int
    """
    output_logfile_dir = (f"{join(PROVIDENTIA_ROOT, 'logs/interpolation/interpolation_logs/')}"
    f"{submit_args['prov_exp_code']}/"
    f"{submit_args['original_speci_to_process']}/"
    f"{submit_args['network_to_interpolate_against']}/"
    f"{submit_args['temporal_resolution_to_output']}/"
    f"{submit_args['yearmonth']}")

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
        submit_args = {'prov_exp_code': sys.argv[1], 
                       'model_temporal_resolution': sys.argv[2], 
                       'speci_to_process': sys.argv[3], 
                       'network_to_interpolate_against': sys.argv[4], 
                       'temporal_resolution_to_output': sys.argv[5], 
                       'yearmonth': sys.argv[6], 
                       'original_speci_to_process': sys.argv[7],
                       'job_id': sys.argv[8]
                       }   

        # initialise ExperimentInterpolation object
        EI = ExperimentInterpolation(submit_args)
 
        # read model domain information
        EI.get_model_information()

        # create polygon along edge of model domain
        EI.create_grid_domain_edge_polygon()

        # get necessary data objects from the observational file
        EI.get_observational_objects()

        # get unit conversion factor between observations and experiment data
        EI.get_conversion_factor()

        # read relevant monthly model data into memory
        EI.get_monthly_model_data()

        # get interpolation weights of model grid to observational stations 
        # (using inverse distance weighting interpolation) 
        EI.n_nearest_neighbour_inverse_distance_weights()

        # write out netCDF for yearmonth, interpolating model data to surface observational stations
        EI.write_netCDF()

        # get total time of interpolation
        interpolation_time = time.time() - interpolation_start

        # return valid process logfile (0)
        EI.log_file_str += str((time.time() - interpolation_start)/60.)
        create_output_logfile(0, EI.log_file_str)

    # write error log file if have uncaught internal error
    except Exception as e:
        log_file_str = 'STARTING INTERPOLATION\n'
        log_file_str += str(traceback.format_exc())
        create_output_logfile(2, log_file_str)
