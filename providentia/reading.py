""" Module storing reading functions """
from netCDF4 import num2date
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import bisect


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


def read_netcdf_data(tuple_arguments):
    """Function that handles reading of observational/experiment
    netCDF data also handles filtering of observational data based
    on selected qa/flag/classification flags. If file does not exist,
    returns None.
    """

    # assign arguments from tuple to variables
    relevant_file, time_array, station_references, active_species, process_type, \
    selected_qa, selected_flags, data_dtype, data_vars_to_read, \
    metadata_dtype, metadata_vars_to_read = tuple_arguments

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
        # remove microseconds
        file_time = pd.to_datetime([t.replace(microsecond=0) for t in file_time])

    # get valid file time indices (i.e. those times in active full time array)
    valid_file_time_indices = \
        np.array([i for i, val in enumerate(file_time)
                  if (val >= time_array[0]) & (val <= time_array[-1])],
                 dtype=np.int)
    # cut file time for valid indices
    file_time = file_time[valid_file_time_indices]
    # get indices relative to active full time array
    full_array_time_indices = np.searchsorted(time_array, file_time)

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
    if process_type == 'observations':

        file_data = np.full((len(current_file_station_indices),
                             len(valid_file_time_indices)), np.NaN,
                            dtype=data_dtype)
        for data_var in data_vars_to_read:
            if data_var == 'time':
                # len(len(file_data)) = number of stations
                # save time as unix time to avoid dtype issues
                unix_time = np.array([t.value // 10**9 for t in file_time])
                file_data['time'] = unix_time.repeat(len(file_data))\
                    .reshape(file_time.shape[0], len(file_data)).T
            else:
                file_data[data_var][:, :] = ncdf_root[data_var][current_file_station_indices,
                                                                valid_file_time_indices]

        # if some qa flags selected then screen
        if len(selected_qa) > 0:
            # screen out observations which are associated with any of the selected qa flags
            file_data[active_species]\
                    [np.isin(ncdf_root['qa'][:, valid_file_time_indices, :],
                             selected_qa).any(axis=2)] = np.NaN

        # if some data provider flags selected then screen
        if len(selected_flags) > 0:
            # screen out observations which are associated with any of the
            # selected data provider flags
            file_data[active_species]\
                    [np.isin(ncdf_root['flag'][:, valid_file_time_indices, :],
                             selected_flags).any(axis=2)] = np.NaN

        # get file metadata
        file_metadata = np.full((len(file_station_references), 1), np.NaN, dtype=metadata_dtype)
        for meta_var in metadata_vars_to_read:
            file_metadata[meta_var][current_file_station_indices, 0] = ncdf_root[meta_var][:]

    else:
        file_data = np.full((len(current_file_station_indices),
                             len(valid_file_time_indices)), np.NaN,
                            dtype=data_dtype[:1])
        
        relevant_data = ncdf_root[data_vars_to_read[0]][current_file_station_indices, valid_file_time_indices]
        # mask out fill values for parameter field
        relevant_data[relevant_data.mask] = np.NaN
        file_data[data_vars_to_read[0]][:, :] = relevant_data

    # close netCDF
    ncdf_root.close()

    # return valid species data, time indices relative to active full time array,
    # file station indices relative to all unique station references array
    if process_type == 'observations':
        return file_data, full_array_time_indices, full_array_station_indices, file_metadata
    else:
        return file_data, full_array_time_indices, full_array_station_indices


def read_netcdf_nonghost(tuple_arguments):
    """Function to handle reading of non-ghost files. If file to be read
    does not exist, returns None"""

    # assign arguments from tuple to variables
    relevant_file, time_array, station_references, active_species, process_type = tuple_arguments
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
        # remove microseconds
        file_time = pd.to_datetime([t.replace(microsecond=0) for t in file_time])

    # get valid file time indices (i.e. those times in active full time array)
    valid_file_time_indices = \
        np.array([i for i, val in enumerate(file_time)
                  if (val >= time_array[0]) & (val <= time_array[-1])],
                 dtype=np.int)
    # cut file time for valid indices
    file_time = file_time[valid_file_time_indices]

    # get indices relative to active full time array
    full_array_time_indices = np.searchsorted(time_array, file_time)

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
        # for data_var in data_vars_to_read:
        file_data[:] = ncdf_root[active_species][valid_file_time_indices, current_file_station_indices].T

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
    # if process_type == 'observations':
    return file_data, full_array_time_indices, full_array_station_indices


def get_yearmonths_to_read(yearmonths, start_date_to_read, end_date_to_read):
    """Function that returns the yearmonths of the files needed to be read.
       This is done by limiting a list of relevant yearmonths by a start/end date
    """

    first_valid_file_ind = bisect.bisect_right(yearmonths, int(start_date_to_read))
    if first_valid_file_ind != 0:
        first_valid_file_ind = first_valid_file_ind - 1
    last_valid_file_ind = bisect.bisect_left(yearmonths, int(end_date_to_read))
    if first_valid_file_ind == last_valid_file_ind:
        return [yearmonths[first_valid_file_ind]]
    else:
        return yearmonths[first_valid_file_ind:last_valid_file_ind]
