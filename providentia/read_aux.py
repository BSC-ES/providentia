""" Module storing static reading functions """
from netCDF4 import num2date
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import bisect
import time

#initialise dictionary for storing pointers to shared memory variables in read step 
shared_memory_vars = {}

def drop_nans(data):
    """Function that returns numpy object of lists
    of station data with NaNs removed
    """

    # reshape numpy array to have lists of data per station
    data = data.tolist()
    # iterate through each list of station data and remove NaNs
    for station_ii, station_data in enumerate(data):
        data[station_ii] = np.array(station_data)[~np.isnan(station_data)]
    # return numpy object of lists of station data with NaNs removed
    return np.array(data)

def init_shared_vars_read_netcdf_data(file_data, file_data_shape, timestamp_array, qa, flags):
    """Function which called before netCDF read function,
       to initialise each worker process.
       Purpose of this function is to access shared memory variables  
    """

    shared_memory_vars['file_data'] = file_data
    shared_memory_vars['file_data_shape'] = file_data_shape
    shared_memory_vars['timestamp_array'] = timestamp_array
    shared_memory_vars['qa'] = qa
    shared_memory_vars['flag'] = flags


def read_netcdf_data(tuple_arguments):
    """Function that handles reading of observational/experiment
    netCDF data also handles filtering of observational data based
    on selected qa/flag/classification flags. If file does not exist,
    returns None.
    """

    # assign arguments from tuple to variables
    relevant_file, station_references, active_species, process_type, \
    data_vars_to_read, metadata_dtype, metadata_vars_to_read = tuple_arguments

    #wrap shared arrays as numpy arrays to more easily manipulate the data
    file_data_shared = np.frombuffer(shared_memory_vars['file_data'], dtype=np.float32).reshape(shared_memory_vars['file_data_shape'])

    # read netCDF frame, if files doesn't exist, return with None
    try:
        ncdf_root = Dataset(relevant_file)
    except Exception as e:
        return

    # get time units
    time_units = ncdf_root['time'].units

    # get file time (handle monthly resolution data differently to hourly/daily
    # as num2date does not support 'months since' units)
    if 'months' in time_units:
        monthly_start_date = time_units.split(' ')[2]
        file_time = pd.date_range(start=monthly_start_date, periods=1, freq='MS')
    else:
        file_time = num2date(ncdf_root['time'][:], time_units)
        # remove microseconds and convert to integer
        file_time = pd.to_datetime([t.replace(microsecond=0) for t in file_time])
    #get file time as integer timestamp
    file_timestamp = file_time.asi8

    # get valid file time indices (i.e. those times in active full time array)
    valid_file_time_indices = np.where(np.logical_and(file_timestamp>=shared_memory_vars['timestamp_array'][0], file_timestamp<=shared_memory_vars['timestamp_array'][-1]))[0]

    # get indices relative to active full timestamp array
    full_array_time_indices = np.searchsorted(shared_memory_vars['timestamp_array'], file_timestamp[valid_file_time_indices])

    # get all station references in file
    file_station_references = ncdf_root['station_reference'][:]

    # get indices of all unique station references that are contained
    # within file station references array
    full_array_station_indices = \
        np.where(np.in1d(station_references, file_station_references))[0]

    # get indices of file station station references that are
    # contained in all unique station references array
    current_file_station_indices = \
        np.where(np.in1d(file_station_references, station_references))[0]

    # for observations, set species data based on selected qa flags/standard data provider
    # flags/classifications to retain or remove as NaN
    s = time.time()
    if process_type == 'observations':
        for data_var_ii, data_var in enumerate(data_vars_to_read):
            #if species variable, then if need to filter by qa load non-filtered array, otherwise load prefiltered array (if available)
            if data_var == active_species:
                if (not shared_memory_vars['qa']) & ('{}_prefiltered_defaultqa'.format(data_var) in list(ncdf_root.variables.keys())):
                    species_data = ncdf_root['{}_prefiltered_defaultqa'.format(data_var)][current_file_station_indices, valid_file_time_indices]
                else:
                    species_data = ncdf_root[data_var][current_file_station_indices, valid_file_time_indices]
            else:
                file_data_shared[data_var_ii,full_array_station_indices[:, np.newaxis],full_array_time_indices[np.newaxis, :]] =\
                    ncdf_root[data_var][current_file_station_indices, valid_file_time_indices]

        # if some qa flags selected then screen
        if shared_memory_vars['qa']:
            # screen out observations which are associated with any of the selected qa flags
            species_data[np.isin(ncdf_root['qa'][current_file_station_indices, valid_file_time_indices, :], shared_memory_vars['qa']).any(axis=2)] = np.NaN

        # if some data provider flags selected then screen
        if shared_memory_vars['flag']:
            # screen out observations which are associated with any of the selected data provider flags
            species_data[np.isin(ncdf_root['flag'][current_file_station_indices, valid_file_time_indices, :], shared_memory_vars['flag']).any(axis=2)] = np.NaN

        #write filtered species data to shared file data
        file_data_shared[data_vars_to_read.index(active_species),full_array_station_indices[:, np.newaxis],full_array_time_indices[np.newaxis, :]] =\
            species_data

        # get file metadata
        file_metadata = np.full((len(file_station_references), 1), np.NaN, dtype=metadata_dtype)
        for meta_var_ii, meta_var in enumerate(metadata_vars_to_read):
            file_metadata[meta_var][current_file_station_indices,0] = ncdf_root[meta_var][:]

    #experiment data
    else:
        relevant_data = ncdf_root[active_species][current_file_station_indices, valid_file_time_indices]
        # mask out fill values for parameter field
        relevant_data[relevant_data.mask] = np.NaN
        file_data_shared[data_vars_to_read.index(active_species),full_array_station_indices[:, np.newaxis],full_array_time_indices[np.newaxis, :]] = relevant_data

    # close netCDF
    ncdf_root.close()

    # return valid species data, time indices relative to active full time array,
    # file station indices relative to all unique station references array
    if process_type == 'observations':
        return full_array_station_indices, file_metadata
    else:
        return full_array_station_indices


def read_netcdf_nonghost(tuple_arguments):
    """Function to handle reading of non-ghost files. If file to be read
    does not exist, returns None"""

    # assign arguments from tuple to variables
    relevant_file, station_references, active_species, process_type = tuple_arguments
    # nonghost have separate metadata to read
    # read netCDF frame, if files doesn't exist, return with None
    try:
        ncdf_root = Dataset(relevant_file)
    except Exception as e:
        return

    latitude = "latitude"
    longitude = "longitude"
    if "latitude" not in ncdf_root.variables:
        latitude = "lat"
        longitude = "lon"
    metadata_vars_to_read = ['station_name', latitude, longitude, 'altitude']
    metadata_dtype = [('station_name', np.object), (latitude, np.float),
                      (longitude, np.float), ('altitude', np.float)]
    if "station_code" in ncdf_root.variables:
        # rename it to station_reference for consistency
        metadata_vars_to_read.append('station_reference')
        metadata_dtype.append(('station_reference', np.object))
    if "station_type" in ncdf_root.variables:
        metadata_vars_to_read.append('station_type')
        metadata_dtype.append(('station_type', np.object))
    if "station_area" in ncdf_root.variables:
        metadata_vars_to_read.append('station_area')
        metadata_dtype.append(('station_area', np.object))
    # get time units
    time_units = ncdf_root['time'].units

    # get file time (handle monthly resolution data differently to hourly/daily
    # as num2date does not support 'months since' units)
    if 'months' in time_units:
        monthly_start_date = time_units.split(' ')[2]
        file_time = pd.date_range(start=monthly_start_date, periods=1, freq='MS')
    else:
        file_time = num2date(ncdf_root['time'][:], time_units)
        # remove microseconds
        file_time = pd.to_datetime([t.replace(microsecond=0) for t in file_time])
    

    # get valid file time indices (i.e. those times in active full time array)
    valid_file_time_indices = \
        np.array([i for i, val in enumerate(file_time)
                  if (val >= shared_memory_vars['time_array'][0]) & (val <= shared_memory_vars['time_array'][-1])],
                 dtype=np.int)
    # cut file time for valid indices
    file_time = file_time[valid_file_time_indices]

    # get indices relative to active full time array
    full_array_time_indices = np.searchsorted(shared_memory_vars['time_array'], file_time)

    # get all station references in file
    # file_station_references = station_references  #ncdf_root['station_name'][:]
    if process_type == 'observations':
        file_station_references = np.array([st_name.tostring().decode('ascii').replace('\x00', '')
                                        for st_name in ncdf_root['station_name'][:]], dtype=np.str)
    # if we're reading exp, the station references have been handled and exist in
    # in variable station_reference
    else:
        file_station_references = ncdf_root['station_reference'][:]
    # get indices of all unique station references that are contained
    # within file station references array
    full_array_station_indices = \
        np.where(np.in1d(station_references, file_station_references))[0]
    # get indices of file station station references that are
    # contained in all unique station references array
    current_file_station_indices = \
        np.where(np.in1d(file_station_references, station_references))[0]
    # for observations, set species data based on selected qa flags/standard data provider
    # flags/classifications to retain or remove as NaN
    if process_type == 'observations':
        file_data = np.full((len(current_file_station_indices),
                             len(valid_file_time_indices)), np.NaN)
        file_metadata = np.full((len(file_station_references), 1), np.NaN, dtype=metadata_dtype)
        # for data_var in data_vars_to_read:
        file_data[:] = ncdf_root[active_species][valid_file_time_indices, current_file_station_indices].T
        for meta_var in metadata_vars_to_read:
            if meta_var == "station_name":
                file_metadata[meta_var][current_file_station_indices, 0] = station_references
            elif meta_var == "station_reference":
                codes = np.array([st_code.tostring().decode('ascii').replace('\x00', '')
                                  for st_code in ncdf_root['station_code'][:]], dtype=np.str)
                file_metadata[meta_var][current_file_station_indices, 0] = codes
            elif meta_var == "station_type":
                file_station_types = np.array([st_type.tostring().decode('ascii').replace('\x00', '')
                                               for st_type in ncdf_root['station_type'][:]], dtype=np.str)
                file_metadata[meta_var][current_file_station_indices, 0] = file_station_types
            elif meta_var == "station_area":
                file_station_areas = np.array([st_area.tostring().decode('ascii').replace('\x00', '')
                                               for st_area in ncdf_root['station_area'][:]], dtype=np.str)
                file_metadata[meta_var][current_file_station_indices, 0] = file_station_areas
            else:
                file_metadata[meta_var][current_file_station_indices, 0] = ncdf_root[meta_var][:]

    else:
        file_data = np.full((len(current_file_station_indices),
                             len(valid_file_time_indices)), np.NaN)

        relevant_data = ncdf_root[active_species][current_file_station_indices, valid_file_time_indices]
        # mask out fill values for parameter field
        relevant_data[relevant_data.mask] = np.NaN
        file_data[:] = relevant_data

    # close netCDF
    ncdf_root.close()

    # return valid species data, time indices relative to active full time array,
    # file station indices relative to all unique station references array
    if process_type == 'observations':
        return file_data, full_array_time_indices, full_array_station_indices, file_metadata
    else:
        return file_data, full_array_time_indices, full_array_station_indices


def get_yearmonths_to_read(yearmonths, start_date_to_read, end_date_to_read, resolution):
    """Function that returns the yearmonths of the files needed to be read.
       This is done by limiting a list of relevant yearmonths by a start/end date
    """

    first_valid_file_ind = bisect.bisect_right(yearmonths, int(start_date_to_read))
    if first_valid_file_ind != 0:
        first_valid_file_ind = first_valid_file_ind - 1
    last_valid_file_ind = bisect.bisect_left(yearmonths, int(end_date_to_read))

    # read only complete months
    if (resolution == 'monthly') and (str(end_date_to_read)[6:8] != '01'):
        last_valid_file_ind -= 1

    if first_valid_file_ind == last_valid_file_ind:
        return [yearmonths[first_valid_file_ind]]
    else:
        return yearmonths[first_valid_file_ind:last_valid_file_ind]
