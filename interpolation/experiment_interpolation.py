#WRITTEN BY DENE BOWDALO

###--------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------###

#experiment_interpolation.py

#module which interpolates experiment output to surface observations 

###--------------------------------------------------------------------------------------------------###
###IMPORT MODULES
###--------------------------------------------------------------------------------------------------###

import experiment_interpolation_functions
import glob
from netCDF4 import Dataset
import inspect
import numpy as np
import os
import sys
import time
import traceback

#change current working directory from submit directory to import global configuration file
working_directory = os.getcwd().split('/submit')[0]
os.chdir(working_directory)

#Read configuration file
from configuration import *

###--------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------###

#initialise log file string
log_file_str = 'STARTING INTERPOLATION\n'

try:
    #time start of yearmonth interpolation
    interpolation_start = time.time()

    #get arguments passed from submittal script
    experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]

    #create output logfile directory
    output_logfile_dir = '{}/interpolation_logs/{}/{}/{}/{}/{}/{}/{}'.format(working_directory, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth)

    #read defined experiments dictionary
    from defined_experiments import defined_experiments_dictionary

    #get experiment specific directory (take gpfs experiment directory preferentially over esarchive directory)
    exp_dict = defined_experiments_dictionary[experiment_to_process]
    if 'gpfs' in list(exp_dict.keys()):
        exp_dir = exp_dict['gpfs']
    else:
        exp_dir = exp_dict['esarchive']

    #get relevant observational file
    obs_file = glob.glob('/gpfs/projects/bsc32/AC_cache/obs/ghost/{}/{}/{}/{}/{}_{}*.nc'.format(GHOST_network_to_interpolate_against, GHOST_version, temporal_resolution_to_output, speci_to_process, speci_to_process, yearmonth))[0]

    #get relevant model files
    model_files = np.sort(glob.glob('{}/{}/{}/{}/{}_{}*.nc'.format(exp_dir, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, speci_to_process, yearmonth)))

    #read model domain information
    mod_nc_root, mod_grid_type, mod_speci_units, mod_lons_centre, mod_lats_centre, x_N, y_N, x_varname, y_varname, z_index = experiment_interpolation_functions.get_model_information(model_files, speci_to_process, yearmonth, log_file_str, output_logfile_dir)

    #create polygon along edge of model domain
    model_grid_outline, model_grid_outline_poly, mod_lons_centre, mod_lats_centre = experiment_interpolation_functions.create_grid_domain_edge_polygon(mod_nc_root, mod_grid_type, mod_lons_centre, mod_lats_centre, x_varname, y_varname, log_file_str, output_logfile_dir)

    #read relevant monthly model data into memory
    yearmonth_time, monthly_model_data, log_file_str = experiment_interpolation_functions.get_monthly_model_data(speci_to_process, yearmonth, temporal_resolution_to_output, model_files, x_N, y_N, z_index, log_file_str)

    #get observational file netCDF root
    obs_nc_root = Dataset(obs_file)

    #get interpolation weights of model grid to observational stations (using inverse distance weighting interpolation) 
    nearest_neighbour_inds, inverse_dists = experiment_interpolation_functions.n_nearest_neighbour_inverse_distance_weights(obs_nc_root['longitude'][:], obs_nc_root['latitude'][:], mod_lons_centre, mod_lats_centre, model_grid_outline_poly, n_neighbours=int(n_neighbours_to_find))

    #write out yearmonth netCDF, interpolating model data to surface observational stations
    experiment_interpolation_functions.write_yearmonth_netCDF(obs_nc_root, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth, int(n_neighbours_to_find), mod_speci_units, mod_lons_centre, mod_lats_centre, x_N, y_N, model_grid_outline, yearmonth_time, monthly_model_data, nearest_neighbour_inds, inverse_dists, GHOST_version)

    #get total time of interpolation
    interpolation_time = time.time() - interpolation_start

    #return valid process logfile (0)
    log_file_str += str((time.time() - interpolation_start)/60.)
    experiment_interpolation_functions.create_output_logfile(0,log_file_str,output_logfile_dir)

#write error log file if have uncaught internal error
except Exception as e:
    try:
        log_file_str = inspect.trace()[-1][0].f_locals['log_file_str']
    except:
        pass
    log_file_str += str(traceback.format_exc())
    experiment_interpolation_functions.create_output_logfile(2,log_file_str,output_logfile_dir)