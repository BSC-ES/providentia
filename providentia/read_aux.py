""" Auxiliary reading functions """

import datetime
from glob import glob
import json
import os
import sys

import bisect
import cftime
from netCDF4 import Dataset, num2date, chartostring
import numpy as np
from packaging.version import Version
import pandas as pd

from providentia.auxiliar import CURRENT_PATH, join
from providentia.warnings_prv import show_message

# initialise dictionary for storing pointers to shared memory variables in read step 
shared_memory_vars = {}

PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])


def drop_nans(data):
    """
    Returns a 1D numpy array with all NaN values removed.

    Parameters
    ----------
    data : numpy.ndarray or pandas.Series
        The input data containing potential NaN values.

    Returns
    -------
    numpy.ndarray
        A filtered 1D array containing only valid numerical values.
    """
    return data[~pd.isnull(data)]


def init_shared_vars_read_netcdf_data(data_in_memory, data_in_memory_shape, ghost_data_in_memory, 
                                      ghost_data_in_memory_shape, timestamp_array, qa, flags):
    """
    Initialises worker processes by granting access to shared memory variables before parallel NetCDF reading.

    Parameters
    ----------
    data_in_memory : multiprocessing.RawArray
        Shared memory buffer for storing the primary scientific data.
    data_in_memory_shape : tuple
        The dimensions of the primary data array.
    ghost_data_in_memory : multiprocessing.RawArray
        Shared memory buffer for storing metadata or GHOST-specific processed data.
    ghost_data_in_memory_shape : tuple
        The dimensions of the GHOST data array.
    timestamp_array : multiprocessing.RawArray
        Shared memory buffer containing the temporal coordinates for the data.
    qa : multiprocessing.RawArray
        Shared memory buffer for Quality Assurance scores.
    flags : multiprocessing.RawArray
        Shared memory buffer for data quality flags.
    """

    shared_memory_vars['data_in_memory'] = data_in_memory
    shared_memory_vars['data_in_memory_shape'] = data_in_memory_shape
    shared_memory_vars['ghost_data_in_memory'] = ghost_data_in_memory
    shared_memory_vars['ghost_data_in_memory_shape'] = ghost_data_in_memory_shape
    shared_memory_vars['timestamp_array'] = timestamp_array
    shared_memory_vars['qa'] = qa
    shared_memory_vars['flag'] = flags


def read_netcdf_data(tuple_arguments):
    """
    Handles the reading and filtering of observational or model netCDF data.

    Parameters
    ----------
    tuple_arguments : tuple
        A collection of arguments including file paths, station references, species names, 
        shared memory handles, and quality control settings.

    Returns
    -------
    numpy.ndarray or list or None
        Returns a structured metadata array if reading observations; an empty list if no 
        stations are found; or None for model data or missing files.
    """

    # assign arguments from tuple to variables
    relevant_file, station_references, station_names, speci,\
    observations_data_label, data_label, data_labels, reading_ghost, ghost_data_vars_to_read,\
    metadata_dtype, metadata_vars_to_read, logger, default_qa, filter_read, network, forecast_indices = tuple_arguments

    # wrap shared arrays as numpy arrays to more easily manipulate the data
    data_in_memory = np.frombuffer(shared_memory_vars['data_in_memory'], dtype=np.float32).reshape(shared_memory_vars['data_in_memory_shape'][:])
    if (reading_ghost or network == 'actris/actris') & (data_label == observations_data_label): 
        qa = np.frombuffer(shared_memory_vars['qa'], dtype=np.uint8)
        flags = np.frombuffer(shared_memory_vars['flag'], dtype=np.uint8)
        if reading_ghost:
            if not filter_read:
                ghost_data_in_memory = np.frombuffer(shared_memory_vars['ghost_data_in_memory'], 
                                                    dtype=np.float32).reshape(shared_memory_vars['ghost_data_in_memory_shape'][:])
    timestamp_array = np.frombuffer(shared_memory_vars['timestamp_array'], dtype=np.int64)

    # read netCDF frame
    ncdf_root = Dataset(relevant_file)

    # get file time (handle monthly resolution data differently to hourly/daily
    # as num2date does not support 'months since' units)
    file_time = ncdf_root['time'][:] 
    time_units = ncdf_root['time'].units
    
    if 'months' in time_units:
        monthly_start_date = time_units.split(' ')[2]
        file_time_dt = pd.date_range(start=monthly_start_date, periods=1, freq='MS')
    else:
        file_time_dt = num2date(file_time, units=time_units)
        
        # convert to pandas datetime
        if Version(cftime.__version__) <= Version("1.0.3.4"):
            # remove microseconds
            file_time_dt = pd.to_datetime([t.replace(microsecond=0) for t in file_time_dt])
        else:
            # bug fix for newer versions of cftime
            file_time_dt = file_time_dt.astype('datetime64[s]')
            file_time_dt = pd.to_datetime([t for t in file_time_dt])

    # get file time as integer timestamp
    file_timestamp = file_time_dt.asi8

    # get valid file time indices (i.e. those times in active full time array)
    valid_file_time_indices = np.where(np.logical_and(file_timestamp>=timestamp_array[0], 
                                                      file_timestamp<=timestamp_array[-1]))[0]

    # get indices relative to active full timestamp array
    full_array_time_indices = np.searchsorted(timestamp_array, file_timestamp[valid_file_time_indices])

    # get all station references in file (do little extra work to get non-GHOST observational station references)
    if (not reading_ghost) & (data_label == observations_data_label):
        if 'station_reference' in ncdf_root.variables:
            station_reference_var = 'station_reference'
        elif 'station_code' in ncdf_root.variables:
            station_reference_var = 'station_code'
        elif 'station_name' in ncdf_root.variables:
            station_reference_var = 'station_name'
        else: 
            error = 'Error: {} cannot be read because it has no station_name.'.format(relevant_file)
            logger.error(error)
            sys.exit(1)

        meta_shape = ncdf_root[station_reference_var].shape
        file_station_references = ncdf_root[station_reference_var][:]
        meta_val_dtype = np.array([file_station_references[0]]).dtype

        if len(meta_shape) == 2:
            if meta_val_dtype == np.dtype(object):
                file_station_references = np.array([''.join(val) for val in file_station_references])
            else:
                file_station_references = chartostring(file_station_references)

    # GHOST and interpolated model data
    else:
        file_station_references = ncdf_root['station_reference'][:]

    # get indices of all non-NaN stations (can be NaN for some non-GHOST files)
    non_nan_station_indices = np.array([ref_ii for ref_ii, ref in enumerate(file_station_references) if ref.lower() != 'nan'])
    file_station_references = file_station_references[non_nan_station_indices]

    # get indices of file station station references that are contained in all unique station references array
    current_file_station_indices = np.where(np.in1d(file_station_references, station_references))[0]

    # for all unique station references that are contained within file station references array
    # get the index of the station reference in the unique station references array 
    index = np.argsort(station_references)
    sorted_station_references = station_references[index]
    sorted_index = np.searchsorted(sorted_station_references, file_station_references[current_file_station_indices])
    full_array_station_indices = np.take(index, sorted_index, mode="clip")

    # if have zero current_file_station_indices in all unique station references, 
    # then check if it is because of old-style of Providentia-interpolation output, 
    # where all station_references were for 'station_name'  
    if (data_label != observations_data_label) & (len(current_file_station_indices) == 0):

        # get indices of file station station references that are contained in all unique station references array
        current_file_station_indices = np.where(np.in1d(file_station_references, station_names))[0]

        # for all unique station references that are contained within file station references array
        # get the index of the station reference in the unique station references array 
        index = np.argsort(station_names)
        sorted_station_names = station_names[index]
        sorted_index = np.searchsorted(sorted_station_names, file_station_references[current_file_station_indices])
        full_array_station_indices = np.take(index, sorted_index, mode="clip")

    # if still have zero current_file_station_indices in all unique station references (can happen due to station colocation)
    # then return from function without reading
    if len(current_file_station_indices) == 0:
        # return empty metadata list if reading observations
        if (data_label == observations_data_label) & (not filter_read):
            return []
        else:
            return 

    # read observations
    if data_label == observations_data_label:

        # read species variable
        # GHOST
        if reading_ghost:
            # if need to filter by qa load non-filtered array, otherwise load prefiltered array (if available)
            if (default_qa) & ('{}_prefiltered_defaultqa'.format(speci) in list(ncdf_root.variables.keys())):
                species_data = ncdf_root['{}_prefiltered_defaultqa'.format(speci)][current_file_station_indices, 
                                                                                valid_file_time_indices]
                # set qa to None as not filtering by them
                qa = None
            else:
                species_data = ncdf_root[speci][current_file_station_indices, valid_file_time_indices]
        # non-GHOST
        else:
            # transpose array to swap station and time dimensions
            if ncdf_root[speci].dimensions == ('time', 'station'):
                species_data = ncdf_root[speci][valid_file_time_indices, current_file_station_indices].T
            # do not transpose
            else:
                species_data = ncdf_root[speci][current_file_station_indices, valid_file_time_indices]
        
        # reading GHOST data?
        if reading_ghost:

            # read GHOST data variables
            if not filter_read:
                for ghost_data_var_ii, ghost_data_var in enumerate(ghost_data_vars_to_read):
                    ghost_data_in_memory[ghost_data_var_ii, full_array_station_indices[:, np.newaxis], 
                                        full_array_time_indices[np.newaxis, :]] = \
                        ncdf_root[ghost_data_var][current_file_station_indices, valid_file_time_indices]
        
        if (reading_ghost) or (network == 'actris/actris'):
            # if some qa flags selected then screen observations
            if qa is not None:
                if len(qa) > 0:
                    # screen out observations which are associated with any of the selected qa flags
                    species_data[np.isin(ncdf_root['qa'][current_file_station_indices, valid_file_time_indices, :], 
                                        qa).any(axis=2)] = np.nan
                
            # if some data provider flags selected then screen observations
            if len(flags) > 0:
                # screen out observations which are associated with any of the selected data provider flags
                species_data[np.isin(ncdf_root['flag'][current_file_station_indices, valid_file_time_indices, :], 
                                    flags).any(axis=2)] = np.nan

        # write filtered species data to shared file data
        data_in_memory[data_labels.index(observations_data_label), full_array_station_indices[:, np.newaxis], 
                    full_array_time_indices[np.newaxis, :]] = species_data

        # get file metadata
        if not filter_read:
            file_metadata = np.full((len(station_references), 1), np.nan, dtype=metadata_dtype)
            for meta_var in metadata_vars_to_read:
                # do extra work for non-GHOST data 
                if not reading_ghost:
                    # get correct variable name for .nc
                    if meta_var == 'longitude':
                        if "longitude" in ncdf_root.variables:
                            meta_var_nc = 'longitude'
                        else:
                            meta_var_nc = 'lon'
                    elif meta_var == 'latitude':
                        if "latitude" in ncdf_root.variables:
                            meta_var_nc = 'latitude'
                        else:
                            meta_var_nc = 'lat'
                    elif meta_var == 'altitude':
                        if "altitude" in ncdf_root.variables:
                            meta_var_nc = 'altitude'
                        else:
                            meta_var_nc = 'alt'
                    elif meta_var == 'station_reference':
                        if 'station_reference' in ncdf_root.variables:
                            meta_var_nc = 'station_reference'
                        elif 'station_code' in ncdf_root.variables:
                            meta_var_nc = 'station_code'
                        elif 'station_name' in ncdf_root.variables:
                            meta_var_nc = 'station_name'
                    elif meta_var == 'area_classification':
                        if 'area_classification' in ncdf_root.variables:
                            meta_var_nc = 'area_classification'
                        else:
                            meta_var_nc = 'station_area'
                    elif meta_var == 'station_classification':
                        if 'station_classification' in ncdf_root.variables:
                            meta_var_nc = 'station_classification'
                        else:
                            meta_var_nc = 'station_type'
                    else:
                        meta_var_nc = meta_var
        
                    # check meta variable is in netCDF
                    if meta_var_nc not in ncdf_root.variables:
                        continue

                    meta_shape = ncdf_root[meta_var_nc].shape
                    meta_val = ncdf_root[meta_var_nc][current_file_station_indices]
                    meta_val_dtype = np.array([meta_val[0]]).dtype

                    # do str formatting where neccessary
                    if meta_val_dtype not in [np.int8, np.int16, np.int32, np.int64, 
                                              np.uint8, np.uint16, np.uint32, np.uint64,
                                              np.float16, np.float32, np.float64]:

                        if len(meta_shape) == 2:
                            if meta_val_dtype == np.dtype(object):
                                meta_val = np.array([''.join(val) for val in meta_val])
                            else:
                                meta_val = chartostring(meta_val)
                    
                    # do str formatting (capitalization) to the metadata
                    if isinstance(meta_val,str):
                        meta_val = np.char.capitalize(meta_val)

                # GHOST metadata
                else:
                    meta_var_nc = meta_var

                    # check meta variable is in netCDF
                    if meta_var_nc not in ncdf_root.variables:
                        continue

                    meta_val = ncdf_root[meta_var_nc][current_file_station_indices]

                # put metadata in array
                file_metadata[meta_var][full_array_station_indices, 0] = meta_val

    # model data
    else:

        # determine if data is structured as forecast data or not
        if 'forecast_day' in ncdf_root[speci].dimensions:
            have_forecast = True
            #if so, check what type of forecast data
            # check if have daily forecast set
            daily_forecast = np.any([True for data_label in data_labels if '-daily' in data_label])    
            # check if have combined forecast set
            combined_forecast = np.any([True for data_label in data_labels if '-combined' in data_label]) 
            # check if have day forecast set
            day_forecast = np.any([True for data_label in data_labels if '-day' in data_label])
        # data is not structured as forecast
        else:
            have_forecast = False
        
        # if have no passed forecast indices, then take first day preferentially (if data is structured as forecast data)
        if len(forecast_indices) == 0:
            forecast_indices = np.array([0], dtype=np.int32)

        # iterate through forecast indices
        for forecast_index in forecast_indices:

            # if want a specific forecast day, and it is available in the netCDF, then take it 
            if have_forecast:
                # daily forecast
                if daily_forecast:
                    data_label_forecast = '{}-daily{}'.format(data_label, forecast_index+1)
                # combined forecast
                elif combined_forecast:
                    data_label_forecast = '{}-combined{}'.format(data_label, forecast_index+1)
                # N day forecast
                elif day_forecast:
                    data_label_forecast = '{}-day{}'.format(data_label, forecast_index+1)
                else:
                    data_label_forecast = data_label
                relevant_data = ncdf_root[speci][current_file_station_indices, valid_file_time_indices, forecast_index]
            # else if forecast day not available in the netCDF, then just take the data as it is
            else:
                relevant_data = ncdf_root[speci][current_file_station_indices, valid_file_time_indices]
                data_label_forecast = data_label

            # mask out fill values for parameter field
            relevant_data[relevant_data.mask] = np.nan

            # put data in array
            data_in_memory[data_labels.index(data_label_forecast), full_array_station_indices[:, np.newaxis], 
                           full_array_time_indices[np.newaxis, :]] = relevant_data

    # close netCDF
    ncdf_root.close()

    # return metadata if reading observations
    if (data_label == observations_data_label) & (not filter_read):
        return file_metadata


def read_netcdf_metadata(tuple_arguments):
    """
    Extracts basic metadata from an observational NetCDF file.

    Parameters
    ----------
    tuple_arguments : tuple
        A collection containing the file path to be read, a boolean indicating if it is a GHOST 
        format file, and a logger instance for error reporting.

    Returns
    -------
    list
        A list of NumPy arrays or lists containing the extracted metadata for 
        'station_reference', 'longitude', 'latitude', 'station_name', and 'measurement_altitude'.
    """

    # assign arguments from tuple to variables
    relevant_file, reading_ghost, logger = tuple_arguments

    # read netCDF frame
    ncdf_root = Dataset(relevant_file)

    # set metadata variables to read
    metadata_vars_to_read = ['station_reference', 'longitude', 'latitude', 'station_name', 'measurement_altitude']
    metadata_read = []

    # iterate though metadata variables to read
    for meta_var in metadata_vars_to_read:

        # do extra work for non-GHOST data 
        if not reading_ghost:
        
            # station reference
            if meta_var == 'station_reference':
                if 'station_reference' in ncdf_root.variables:
                    station_reference_var = 'station_reference'
                elif 'station_code' in ncdf_root.variables:
                    station_reference_var = 'station_code'
                elif 'station_name' in ncdf_root.variables:
                    station_reference_var = 'station_name'
                else: 
                    error = 'Error: {} cannot be read because it has no station_name.'.format(relevant_file)
                    logger.error(error)
                    sys.exit(1)

                meta_shape = ncdf_root[station_reference_var].shape
                meta_val = ncdf_root[station_reference_var][:]
                meta_val_dtype = np.array([meta_val[0]]).dtype
                if len(meta_shape) == 2:
                    if meta_val_dtype == np.dtype(object):
                        meta_val = np.array([''.join(val) for val in meta_val])
                    else:
                        meta_val = chartostring(meta_val)

                # get indices of all non-NaN stations (can be NaN for some non-GHOST files)
                non_nan_station_indices = np.array([ref_ii for ref_ii, ref in enumerate(meta_val) if ref.lower() != 'nan'])
                meta_val = meta_val[non_nan_station_indices]

            # longitude
            elif meta_var == 'longitude':
                if "longitude" in ncdf_root.variables:
                    meta_val = ncdf_root['longitude'][non_nan_station_indices]
                else:
                    meta_val = ncdf_root['lon'][non_nan_station_indices]
            
            # latitude
            elif meta_var == 'latitude':
                if "latitude" in ncdf_root.variables:
                    meta_val = ncdf_root['latitude'][non_nan_station_indices]
                else:
                    meta_val = ncdf_root['lat'][non_nan_station_indices]

            # station name
            elif meta_var == 'station_name':
                if "station_name" in ncdf_root.variables:
                    meta_shape = ncdf_root['station_name'].shape
                    meta_val = ncdf_root['station_name'][non_nan_station_indices]
                    meta_val_dtype = np.array([meta_val[0]]).dtype
                    if len(meta_shape) == 2:
                        if meta_val_dtype == np.dtype(object):
                            meta_val = np.array([''.join(val) for val in meta_val])
                        else:
                            meta_val = chartostring(meta_val)
                else:
                    meta_val = []

            # measurement altitude
            elif meta_var == 'measurement_altitude':
                if meta_var in ncdf_root.variables:
                    meta_val = ncdf_root[meta_var][non_nan_station_indices]
                else:
                    meta_val = []

        # GHOST metadata
        else:
            meta_val = ncdf_root[meta_var][:]

        # append read metadata 
        metadata_read.append(meta_val)

    # close netCDF
    ncdf_root.close()

    # return read metadata
    return metadata_read


def check_forecast_dimension(tuple_arguments):
    """
    Determines if a NetCDF file contains a forecast dimension and returns the number of lead days.

    Parameters
    ----------
    tuple_arguments : str
        The file path to the NetCDF file being inspected.

    Returns
    -------
    int
        The size of the 'forecast_day' dimension, or 0 if the dimension is not present.
    """

    # assign arguments from tuple to variables
    relevant_file = tuple_arguments

    # read netCDF frame
    ncdf_root = Dataset(relevant_file)

    # check if forecast day dimension exists
    if 'forecast_day' in ncdf_root.dimensions:
        n_forecast_days = ncdf_root.dimensions['forecast_day'].size
    else:
        n_forecast_days = 0

    # close netCDF
    ncdf_root.close()

    return n_forecast_days


def get_yearmonths_to_read(available_yearmonths, start_date_to_read, end_date_to_read, resolution):
    """
    Determines which monthly data files fall within a specified date range based on temporal resolution.

    Parameters
    ----------
    available_yearmonths : list of str
        List of available file strings in 'YYYYMM' format.
    start_date_to_read : str or int
        The beginning of the requested period in 'YYYYMMDD' format.
    end_date_to_read : str or int
        The end of the requested period in 'YYYYMMDD' format.
    resolution : str
        The temporal resolution of the data (e.g., 'monthly', 'daily', or 'hourly').

    Returns
    -------
    list of str
        A subset of available_yearmonths representing the files required to cover the date range.
    """
    
    available_yearmonthdays = [int(yearmonth+'01') for yearmonth in available_yearmonths]

    first_valid_file_ind = bisect.bisect_right(available_yearmonthdays, int(start_date_to_read))
    if first_valid_file_ind != 0:
        first_valid_file_ind = first_valid_file_ind - 1
    last_valid_file_ind = bisect.bisect_left(available_yearmonthdays, int(end_date_to_read))

    # read only complete months
    if resolution == 'monthly':
        if str(end_date_to_read)[6:8] != '01':
            if str(end_date_to_read)[0:6] == str(available_yearmonths[-1]):
                last_valid_file_ind -= 1
        if str(start_date_to_read)[6:8] != '01':
            if str(start_date_to_read)[0:6] == str(available_yearmonths[0]):
                first_valid_file_ind += 1

    if first_valid_file_ind == last_valid_file_ind:
        return [available_yearmonths[first_valid_file_ind]]
    else:
        return available_yearmonths[first_valid_file_ind:last_valid_file_ind]


def get_default_qa(instance, speci):
    """
    Returns the standard GHOST quality assurance flags based on a species' physical limits.

    Parameters
    ----------
    instance : object
        An instance of the application class containing parameter metadata.
    speci : str
        The name of the chemical species or parameter to evaluate.

    Returns
    -------
    list
        A sorted list of integer codes representing the default QA flags.
    """

    if instance.parameter_dictionary[speci]['extreme_lower_limit'] < 0.0:
        return sorted(instance.default_qa_non_negative)
    else:
        return sorted(instance.default_qa_standard)
    

def get_frequency_code(resolution):
    """
    Maps human-readable temporal resolutions to Pandas offset aliases.

    Parameters
    ----------
    resolution : str
        The descriptive name of the temporal resolution (e.g., 'hourly', 'daily').

    Returns
    -------
    active_frequency_code : str
        The corresponding Pandas frequency string used for resampling and date ranges.
    """
    
    if resolution in ['hourly', 'hourly_instantaneous']:
        active_frequency_code = 'h'
    elif resolution in ['3hourly', '3hourly_instantaneous']:
        active_frequency_code = '3h'
    elif resolution in ['6hourly', '6hourly_instantaneous']:
        active_frequency_code = '6h'
    elif resolution == 'daily':
        active_frequency_code = 'D'
    elif resolution == 'monthly':
        active_frequency_code = 'MS'
    elif resolution == 'annual':
        active_frequency_code = 'YS'

    return active_frequency_code

   
def get_chunk_size(active_resolution, chunk_resolution):
    """
    Calculates the number of data steps required to fill a specific temporal chunk.

    Parameters
    ----------
    active_resolution : str
        The base temporal resolution of the data (e.g., 'hourly', '3hourly').
    chunk_resolution : str
        The target resolution for data partitioning (e.g., 'daily').

    Returns
    -------
    int
        The ratio of the target duration to the base duration.
    """

    # convert resolutions to hours
    hours_dict = {
        "hourly": 1,
        "3hourly": 3,
        "6hourly": 6,
        "daily": 24,
    }

    base_hours = hours_dict[active_resolution]
    target_hours = hours_dict[chunk_resolution]

    return target_hours // base_hours


def check_for_ghost(network_name):
    """
    Determines if a given monitoring network comes from GHOST or not.

    Parameters
    ----------
    network_name : str
        The name of the atmospheric monitoring network to verify.

    Returns
    -------
    bool
        True if the network is a GHOST network, False otherwise.
    """

    if '/' in network_name:
        return False
    else:
        return True


def get_ghost_observational_tree(instance):
    """
    Creates a nested dictionary representing the GHOST observational data tree and exports it to JSON.

    Parameters
    ----------
    instance : object
        An instance of the application class containing GHOST configuration and versioning.

    Returns
    -------
    dict
        A nested dictionary structure: {network: {resolution: {matrix: {species: [YYYYMM, ...]}}}}.
    """

    # create dictionary for storing filetree
    ghost_observation_data = {}

    # iterate through available networks
    for network in instance.ghost_available_networks:

        # check if directory for network exists
        # if not, continue
        if not os.path.exists('%s/%s/%s' % (instance.ghost_root, network, instance.ghost_version)):
            continue

        # write empty dictionary for network
        ghost_observation_data[network] = {}

        # iterate through available resolutions
        for resolution in instance.ghost_available_resolutions:

            # check if directory for resolution exists
            # if not, continue
            if not os.path.exists('%s/%s/%s/%s' % (instance.ghost_root, network, instance.ghost_version, resolution)):
                continue

            # write nested empty dictionary for resolution
            ghost_observation_data[network][resolution] = {}

            # get available species for network/resolution
            available_species = os.listdir('%s/%s/%s/%s' % (instance.ghost_root, network, instance.ghost_version, resolution))

            # iterate through available files per species
            for speci in available_species:

                # get all available netCDF files
                available_files = os.listdir(
                    '%s/%s/%s/%s/%s' % (instance.ghost_root, network, instance.ghost_version, resolution, speci))
                
                # continue if have no files
                if len(available_files) == 0:
                    continue

                # get monthly start date (YYYYMM) of all files
                file_yearmonths = sorted([f.split('_')[-1][:6] for f in available_files if f != 'temporary'])
                
                # get matrix for current species
                if speci in instance.parameter_dictionary:
                    matrix = instance.parameter_dictionary[speci]['matrix']
                    if matrix not in ghost_observation_data[network][resolution]:
                        # write nested empty dictionary for matrix
                        ghost_observation_data[network][resolution][matrix] = {}

                    # write nested dictionary for species, with associated file yearmonths
                    ghost_observation_data[network][resolution][matrix][speci] = file_yearmonths

    # save file tree out to yaml
    with open(join(PROVIDENTIA_ROOT, 'settings/internal/ghost_filetree_{}.json'.format(instance.ghost_version)), 'w') as json_file:
        json.dump(ghost_observation_data, json_file, indent=4)

    return ghost_observation_data


def get_nonghost_observational_tree(instance):
    """
    Creates a nested dictionary representing the non-GHOST observational data tree and exports it to JSON.

    Parameters
    ----------
    instance : object
        An instance of the application class containing non-GHOST configuration and directory roots.

    Returns
    -------
    dict
        A nested dictionary structure: {network: {resolution: {matrix: {species: [YYYYMM, ...]}}}}.
    """

    # create dictionary for storing filetree
    nonghost_observation_data = {}

    # iterate through networks
    for network in instance.nonghost_available_networks:

        # check if directory for network exists
        # if not, continue
        if not os.path.exists('%s/%s' % (instance.nonghost_root, network)):
            continue

        # write empty dictionary for network
        nonghost_observation_data[network] = {}

        # iterate through available resolutions
        for resolution in instance.nonghost_available_resolutions:

            # check if directory for resolution exists
            # if not, continue
            if not os.path.exists('%s/%s/%s' % (instance.nonghost_root, network, resolution)):
                continue

            # write nested empty dictionary for resolution
            nonghost_observation_data[network][resolution] = {}

            # get available species for network/resolution
            available_species = os.listdir('%s/%s/%s' % (instance.nonghost_root, network, resolution))

            # iterate through available files per species
            for speci in available_species:

                # get all available netCDF files 
                available_files = glob('%s/%s/%s/%s/%s_??????.nc' % (instance.nonghost_root, network, resolution, speci, speci))

                # continue if have no files
                if len(available_files) == 0:
                    continue

                # get monthly start date (YYYYMM) of all files
                file_yearmonths = sorted([f.split('_')[-1][:6] for f in available_files])

                # get matrix for current species
                if speci in instance.parameter_dictionary:
                    matrix = instance.parameter_dictionary[speci]['matrix']
                    if matrix not in nonghost_observation_data[network][resolution]:
                        # write nested empty dictionary for matrix
                        nonghost_observation_data[network][resolution][matrix] = {}

                    # write nested dictionary for species, with associated file yearmonths
                    nonghost_observation_data[network][resolution][matrix][speci] = file_yearmonths
        
    # save file tree out to yaml
    with open(join(PROVIDENTIA_ROOT, 'settings/internal/nonghost_filetree.json'), 'w') as json_file:
        json.dump(nonghost_observation_data, json_file, indent=4)

    return nonghost_observation_data


def get_valid_obs_files_in_date_range(instance, start_date, end_date):
    """
    Filters the full observational data tree to identify files available within a specific date range.

    Parameters
    ----------
    instance : object
        An instance of the application class containing the full observational data tree.
    start_date : str
        The start date in 'YYYYMMDD' format.
    end_date : str
        The end date in 'YYYYMMDD' format.
    """

    # create dictionary to store available observational data
    instance.available_observation_data = {}

    # check if start/end date are valid values, if not, return with no valid obs. files
    if (not valid_date(start_date)) or (not valid_date(end_date)):
        msg = f'One of start date ({start_date}) or end date ({end_date}) are not valid.'
        show_message(instance, msg, print=True)
        return

    # check end date is > start date, if not, return with no valid obs. files
    if int(start_date) >= int(end_date):
        msg = f'Start date ({start_date}) exceeds end date ({end_date}).'
        show_message(instance, msg, print=True)
        return

    # check start date and end date are both within if valid date range (19000101 - 20500101),
    # if not, return with no valid obs. files
    if (int(start_date) < 19000101) or (int(end_date) < 19000101) or (int(start_date) >= 20500101) or (int(end_date) >= 20500101):
        msg = f'One of start date ({start_date}) or end date ({end_date}) are not valid.'
        show_message(instance, msg, print=True)
        return 

    # get start date on first of month
    start_date_firstdayofmonth = int(str(start_date)[:6] + '01')

    # iterate through observational dictionary
    for network in instance.all_observation_data:
        for resolution in instance.all_observation_data[network]:
            for matrix in instance.all_observation_data[network][resolution]:
                for speci in instance.all_observation_data[network][resolution][matrix]:
                    
                    # get available files
                    file_yearmonths = instance.all_observation_data[network][resolution][matrix][speci]

                    # get file yearmonths within date range
                    valid_file_yearmonths = sorted([ym for ym in file_yearmonths if
                                                    (int('{}01'.format(ym)) >= start_date_firstdayofmonth) & (int('{}01'.format(ym)) < int(end_date))])
                    
                    # add yearmonths to available observation data dict 
                    if len(valid_file_yearmonths) > 0:
                        if network not in instance.available_observation_data:
                            instance.available_observation_data[network] = {}
                        if resolution not in instance.available_observation_data[network]:
                            instance.available_observation_data[network][resolution] = {}
                        if matrix not in instance.available_observation_data[network][resolution]:
                            instance.available_observation_data[network][resolution][matrix] = {}
                        instance.available_observation_data[network][resolution][matrix][speci] = valid_file_yearmonths


def get_valid_models(instance, start_date, end_date, resolution, networks, species):
    """
    Identifies models within a given date range and chosen networks and species. 

    Parameters
    ----------
    instance : object
        An instance of the application class containing directory roots and menu configurations.
    start_date : str
        The start date in 'YYYYMMDD' format.
    end_date : str
        The end date in 'YYYYMMDD' format.
    resolution : str
        The temporal resolution (e.g. 'hourly', 'daily').
    networks : list of str
        The monitoring networks to match against model data.
    species : list of str
        The chemical species or parameters to match against model data.
    """

    # get all different model names (from providentia-interpolation output dir)
    available_models = []
    if os.path.exists(join(instance.mod_root,instance.ghost_version)):
        available_models = os.listdir('%s/%s' % (instance.mod_root, instance.ghost_version))

    # create dictionary to store available model data
    instance.available_model_data = {}

    #list for saving models to add to models pop-up 
    models_to_add = []

    # get start date on first of month
    start_date_firstdayofmonth = int(str(start_date)[:6] + '01')

    # iterate through networks and species
    for network, speci in zip(networks, species):

        # iterate through available models
        for model in available_models:

            # get folder where interpolated models are saved
            if '/' not in network:           
                files_directory = '%s/%s/%s/%s/%s/%s' % (instance.mod_root, instance.ghost_version, 
                                                         model, resolution, speci, network)
            else:
                files_directory = '%s/%s/%s/%s/%s/%s' % (instance.mod_root, instance.ghost_version, 
                                                          model, resolution, speci,
                                                          network.replace('/', '-'))
                
            # test if interpolated directory exists for model
            # if it does not exit, continue
            if not os.path.exists(files_directory):
                continue
            else:
                # get all available netCDF files (handling potential permissions issues)
                try:
                    available_files = os.listdir(files_directory)
                except PermissionError as error:
                    continue

            # get monthly start date (YYYYMM) of all files
            file_yearmonths = sorted([f.split('_')[-1][:6] for f in available_files])

            # write nested dictionary for model, with associated file yearmonths
            if len(file_yearmonths) > 0:

                # get file yearmonths within date range
                valid_file_yearmonths = sorted([ym for ym in file_yearmonths if 
                                                (int('{}01'.format(ym)) >= start_date_firstdayofmonth) & (int('{}01'.format(ym)) < int(end_date))])

                #if have valid files, then add model to pop-up menu, 
                # and add yearmonths to available model data
                if len(valid_file_yearmonths) > 0:
                    models_to_add.append(model)

                    if network not in instance.available_model_data:
                        instance.available_model_data[network] = {}
                    if resolution not in instance.available_model_data[network]:
                        instance.available_model_data[network][resolution] = {}
                    if speci not in instance.available_model_data[network][resolution]:
                        instance.available_model_data[network][resolution][speci] = {}
                    if model not in instance.available_model_data[network][resolution][speci]:
                        instance.available_model_data[network][resolution][speci][model] = valid_file_yearmonths

    # set list of model names to add on models pop-up
    if instance.mode not in ['report', 'library']:
        models_to_add = np.array(sorted(models_to_add))
        instance.models_menu['models']['labels'] = models_to_add
        instance.models_menu['models']['map_vars'] = models_to_add

def get_possible_temporal_resolutions():
    """
    Returns a comprehensive list of all supported temporal resolutions for data processing.

    Returns
    -------
    list
        A list of strings representing temporal resolutions from hourly to annual.
    """

    # define possible temporal resolutions
    possible_resolutions = ['hourly', 'hourly_instantaneous', '3hourly', '3hourly_instantaneous', 
                            '6hourly', '6hourly_instantaneous', 'daily', 'monthly', 'annual']
    
    return possible_resolutions


def get_temporal_resolution_order():
    """
    Provides a dictionary mapping temporal resolution names to their preferred display order in menus.

    Returns
    -------
    dict
        A dictionary where keys are resolution strings and values are integer ranks.
    """

    resolution_order_dict = {'hourly': 1, '3hourly': 2, '6hourly': 3, 'hourly_instantaneous': 4,
                             '3hourly_instantaneous': 5, '6hourly_instantaneous': 6,
                             'daily': 7, 'monthly': 8}

    return resolution_order_dict


def get_possible_resampling_resolutions(resolution, daily_forecast=False):
    """
    Identifies permissible lower temporal resolutions for resampling based on the base resolution.

    Parameters
    ----------
    resolution : str
        The current temporal resolution of the data (e.g. 'hourly', 'daily').
    daily_forecast : bool, optional
        Flag indicating if the data belongs to a daily forecast cycle (default is False).

    Returns
    -------
    list of str
        A list of target resolutions that are coarser than or equal to the input resolution.
    """
    
    if daily_forecast:
        if resolution in ['hourly', 'hourly_instantaneous']:
            resolutions = ['hourly', '3hourly', '6hourly', 'daily']
        elif resolution in ['3hourly', '3hourly_instantaneous']:
            resolutions = ['3hourly', '6hourly', 'daily']
        elif resolution in ['6hourly', '6hourly_instantaneous']:
            resolutions = ['6hourly', 'daily']
        elif resolution == 'daily':
            resolutions = ['daily']
        elif resolution == 'monthly':
            resolutions = []
        elif resolution == 'annual':
            resolutions = []
    else:
        if resolution in ['hourly', 'hourly_instantaneous']:
            resolutions = ['hourly', '3hourly', '6hourly', 'daily', 'monthly', 'annual']
        elif resolution in ['3hourly', '3hourly_instantaneous']:
            resolutions = ['3hourly', '6hourly', 'daily', 'monthly', 'annual']
        elif resolution in ['6hourly', '6hourly_instantaneous']:
            resolutions = ['6hourly', 'daily', 'monthly', 'annual']
        elif resolution == 'daily':
            resolutions = ['daily', 'monthly', 'annual']
        elif resolution == 'monthly':
            resolutions = ['monthly', 'annual']
        elif resolution == 'annual':
            resolutions = ['annual']

    return resolutions


def get_periodic_relevant_temporal_resolutions(resolution):        
    """
    Identifies the appropriate time cycles for periodic analysis based on a base temporal resolution.

    Parameters
    ----------
    resolution : str
        The name of the selected temporal resolution (e.g. 'hourly', 'daily').

    Returns
    -------
    relevant_temporal_resolutions : list
        A list of periodic groupings (e.g. 'hour', 'dayofweek') valid for the input resolution.
    """

    if 'hourly' in resolution:
        relevant_temporal_resolutions = ['hour', 'dayofweek', 'month']
    elif resolution == 'daily':
        relevant_temporal_resolutions = ['dayofweek', 'month']
    elif resolution == 'monthly':
        relevant_temporal_resolutions = ['month']
    else:
        relevant_temporal_resolutions = []
        
    return relevant_temporal_resolutions


def get_periodic_nonrelevant_temporal_resolutions(resolution):        
    """
    Identifies time cycles that cannot be analysed based on a base temporal resolution.

    Parameters
    ----------
    resolution : str
        The name of the selected temporal resolution (e.g. 'hourly', 'daily').

    Returns
    -------
    nonrelevant_temporal_resolutions : list
        A list of periodic groupings that are scientifically invalid for the input resolution.
    """

    if 'hourly' in resolution:
        nonrelevant_temporal_resolutions = []
    elif resolution == 'daily':
        nonrelevant_temporal_resolutions = ['hour']
    elif resolution == 'monthly':
        nonrelevant_temporal_resolutions = ['hour', 'dayofweek']
    else:
        nonrelevant_temporal_resolutions = ['hour', 'dayofweek', 'month']
        
    return nonrelevant_temporal_resolutions


def valid_date(date_text):
    """
    Validates whether a given date string or integer conforms to the 'YYYYMMDD' format.

    Parameters
    ----------
    date_text : str or int
        The date value to be validated.

    Returns
    -------
    bool
        True if the date is valid and correctly formatted, False otherwise.
    """

    try:
        datetime.datetime.strptime(str(date_text), '%Y%m%d')
        return True
    except Exception as e:
        return False


def generate_file_trees(instance):
    """
    Handles the dynamic generation or loading of observational data catalogues.

    Parameters
    ----------
    instance : object
        An instance of the application class containing configuration flags and versioning.
    """
    
    # get dictionaries of observational GHOST and non-GHOST filetrees, either created dynamically or loaded
    # if have filetree flags, then these overwrite any defaults
    gft = False
    if instance.generate_file_tree:
        gft = True
    elif instance.disable_file_tree:
        gft = False
    # by default generate filetree on MN5
    elif instance.machine in ['mn5', 'nord4']:
        gft = True
    # by default generate filetree locally
    elif instance.filetree_type == 'local':
        gft = True

    # generate file trees
    ghost_filetree_path = join(PROVIDENTIA_ROOT, 'settings/internal/ghost_filetree_{}.json'.format(instance.ghost_version)) 
    nonghost_filetree_path = join(PROVIDENTIA_ROOT, 'settings/internal/nonghost_filetree.json') 

    # generate file trees for ghost
    if gft or (not os.path.exists(ghost_filetree_path)):
        if not os.path.exists(ghost_filetree_path):
            instance.logger.info(f'Generating file tree {ghost_filetree_path}...')
        else:
            instance.logger.info(f'Updating file tree {ghost_filetree_path}...')
        instance.all_observation_data = get_ghost_observational_tree(instance) 
    # load file trees
    else:
        instance.logger.info(f'Loading file tree {ghost_filetree_path}...')
        try:
            instance.all_observation_data = json.load(open(join(PROVIDENTIA_ROOT, 'settings/internal/ghost_filetree_{}.json'.format(instance.ghost_version)))) 
        except FileNotFoundError as file_error:
            error = "Error: Trying to load 'settings/internal/ghost_filetree_{}.json' but file does not exist. Run with the flag '--gft' to generate this file.".format(instance.ghost_version)
            instance.logger.error(error)
            sys.exit(1)

    if instance.nonghost_root is not None:
        # generate file trees for nonghost
        if gft or (not os.path.exists(nonghost_filetree_path)):
            if not os.path.exists(nonghost_filetree_path):
                instance.logger.info(f'Generating file tree {nonghost_filetree_path}...')
            else:
                instance.logger.info(f'Updating file tree {nonghost_filetree_path}...')
            nonghost_observation_data = get_nonghost_observational_tree(instance)
        # load file trees
        else:
            instance.logger.info(f'Loading file tree {nonghost_filetree_path}...')
            try:
                nonghost_observation_data = json.load(open(join(PROVIDENTIA_ROOT, 'settings/internal/nonghost_filetree.json')))
            except FileNotFoundError as file_error:
                error = "Error: Trying to load 'settings/internal/nonghost_filetree.json' but file does not exist. Run with the flag '--gft' to generate this file."
                instance.logger.error(error)
                sys.exit(1)

        # merge GHOST and non-GHOST filetrees
        instance.all_observation_data = {**instance.all_observation_data, **nonghost_observation_data}


def get_valid_metadata(read_instance, variable, valid_station_idxs, networkspeci):
    """
    Extracts the first non-NaN metadata values for selected stations across all time steps.

    Parameters
    ----------
    read_instance : object
        An instance of the application class responsible for data reading operations.
    variable : str
        The specific metadata field to retrieve (e.g., 'station_name', 'latitude').
    valid_station_idxs : list of int
        Indices of the stations to be processed.
    networkspeci : str
        A combined string representing the network and species identifier.

    Returns
    -------
    valid_metadata : list
        A list of non-null metadata values, one per station.
    """

    stations_metadata = read_instance.canvas_instance.selected_station_metadata[networkspeci][
                            variable][valid_station_idxs, :]
    valid_metadata = []
    for station_metadata in  stations_metadata:
        first_valid_station_metadata = next((
            value for value in station_metadata 
            if value == value), 'nan')
        # if metadata returns an array, get first
        if isinstance(first_valid_station_metadata, (np.ndarray)):
            first_valid_station_metadata = first_valid_station_metadata.item()
        valid_metadata.append(first_valid_station_metadata)

    return valid_metadata