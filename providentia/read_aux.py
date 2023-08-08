""" Module storing static reading functions """
from netCDF4 import Dataset, num2date, chartostring
import numpy as np
import pandas as pd
import bisect
import time
import sys

# initialise dictionary for storing pointers to shared memory variables in read step 
shared_memory_vars = {}

def drop_nans(data):
    """ Function that return 1D numpy array with NaNs removed. """
    return data[~pd.isnull(data)]

def init_shared_vars_read_netcdf_data(data_in_memory, data_in_memory_shape, ghost_data_in_memory, 
                                      ghost_data_in_memory_shape, timestamp_array, qa, flags):
    """ Function which called before netCDF read function,
        to initialise each worker process.
        The purpose of this function is to access shared memory variables.
    """
    shared_memory_vars['data_in_memory'] = data_in_memory
    shared_memory_vars['data_in_memory_shape'] = data_in_memory_shape
    shared_memory_vars['ghost_data_in_memory'] = ghost_data_in_memory
    shared_memory_vars['ghost_data_in_memory_shape'] = ghost_data_in_memory_shape
    shared_memory_vars['timestamp_array'] = timestamp_array
    shared_memory_vars['qa'] = qa
    shared_memory_vars['flag'] = flags

def read_netcdf_data(tuple_arguments):
    """ Function that handles reading of observational/experiment
        netCDF data also handles filtering of observational data based
        on selected qa/flag/classification flags. If file does not exist,
        returns None.
    """

    # assign arguments from tuple to variables
    relevant_file, station_references, station_names, speci,\
    data_label, data_labels, reading_ghost, ghost_data_vars_to_read,\
    metadata_dtype, metadata_vars_to_read, default_qa, filter_read = tuple_arguments

    # wrap shared arrays as numpy arrays to more easily manipulate the data
    data_in_memory = np.frombuffer(shared_memory_vars['data_in_memory'], dtype=np.float32).reshape(shared_memory_vars['data_in_memory_shape'][:])
    if (reading_ghost) & (data_label == 'observations'): 
        qa = np.frombuffer(shared_memory_vars['qa'], dtype=np.uint8)
        flags = np.frombuffer(shared_memory_vars['flag'], dtype=np.uint8)
        if not filter_read:
            ghost_data_in_memory = np.frombuffer(shared_memory_vars['ghost_data_in_memory'], 
                                                 dtype=np.float32).reshape(shared_memory_vars['ghost_data_in_memory_shape'][:])
    timestamp_array = np.frombuffer(shared_memory_vars['timestamp_array'], dtype=np.int64)

    # read netCDF frame
    ncdf_root = Dataset(relevant_file)

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

    # get file time as integer timestamp
    file_timestamp = file_time.asi8

    # get valid file time indices (i.e. those times in active full time array)
    valid_file_time_indices = np.where(np.logical_and(file_timestamp>=timestamp_array[0], 
                                                      file_timestamp<=timestamp_array[-1]))[0]

    # get indices relative to active full timestamp array
    full_array_time_indices = np.searchsorted(timestamp_array, file_timestamp[valid_file_time_indices])

    # get all station references in file (do little extra work to get non-GHOST observational station references)
    if (not reading_ghost) & (data_label == 'observations'):
        if 'station_reference' in ncdf_root.variables:
            station_reference_var = 'station_reference'
        elif 'station_code' in ncdf_root.variables:
            station_reference_var = 'station_code'
        elif 'station_name' in ncdf_root.variables:
            station_reference_var = 'station_name'
        else: 
            print('Error: {} cannot be read because it has no station_name.'.format(relevant_file))
            sys.exit()

        meta_shape = ncdf_root[station_reference_var].shape
        file_station_references = ncdf_root[station_reference_var][:]
        meta_val_dtype = np.array([file_station_references[0]]).dtype

        if len(meta_shape) == 2:
            if meta_val_dtype == np.dtype(object):
                file_station_references = np.array([''.join(val) for val in file_station_references])
            else:
                file_station_references = chartostring(file_station_references)

    # GHOST and interpolated experiment data
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
    if (data_label != 'observations') & (len(current_file_station_indices) == 0):

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
        if (data_label == 'observations') & (not filter_read):
            return []
        else:
            return 

    # read observations
    if data_label == 'observations':

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
        # non-GHOST (transpose array to swap station and time dimensions)
        else:
            species_data = ncdf_root[speci][valid_file_time_indices, current_file_station_indices].T

        # reading GHOST data?
        if reading_ghost:

            # read GHOST data variables
            if not filter_read:
                for ghost_data_var_ii, ghost_data_var in enumerate(ghost_data_vars_to_read):
                    ghost_data_in_memory[ghost_data_var_ii, full_array_station_indices[:, np.newaxis], 
                                        full_array_time_indices[np.newaxis, :]] = \
                        ncdf_root[ghost_data_var][current_file_station_indices, valid_file_time_indices]

            # if some qa flags selected then screen observations
            if qa is not None:
                if len(qa) > 0:
                    # screen out observations which are associated with any of the selected qa flags
                    species_data[np.isin(ncdf_root['qa'][current_file_station_indices, valid_file_time_indices, :], 
                                        qa).any(axis=2)] = np.NaN
                
            # if some data provider flags selected then screen observations
            if len(flags) > 0:
                # screen out observations which are associated with any of the selected data provider flags
                species_data[np.isin(ncdf_root['flag'][current_file_station_indices, valid_file_time_indices, :], 
                                    flags).any(axis=2)] = np.NaN

        # write filtered species data to shared file data
        data_in_memory[data_labels.index('observations'), full_array_station_indices[:, np.newaxis], 
                       full_array_time_indices[np.newaxis, :]] = species_data

        # get file metadata
        if not filter_read:
            file_metadata = np.full((len(station_references), 1), np.NaN, dtype=metadata_dtype)
            for meta_var in metadata_vars_to_read:
                # do extra work for non-GHOST data 
                if not reading_ghost:
                    # get correct variable name for .nc
                    if meta_var == 'longitude':
                        if "longitude" not in ncdf_root.variables:
                            meta_var_nc = 'lon'
                        else:
                            meta_var_nc = 'longitude'
                    elif meta_var == 'latitude':
                        if "latitude" not in ncdf_root.variables:
                            meta_var_nc = 'lat'
                        else:
                            meta_var_nc = 'latitude'
                    elif meta_var == 'altitude':
                        if "altitude" not in ncdf_root.variables:
                            meta_var_nc = 'alt'
                        else:
                            meta_var_nc = 'altitude'
                    elif meta_var == 'station_reference':
                        if "station_reference" not in ncdf_root.variables:
                            meta_var_nc = 'station_code'
                        else:
                            meta_var_nc = 'station_reference'
                    else:
                        meta_var_nc = meta_var
        
                    # check meta variable is in netCDF
                    if meta_var_nc not in ncdf_root.variables:
                        continue

                    meta_shape = ncdf_root[meta_var_nc].shape
                    meta_val = ncdf_root[meta_var_nc][current_file_station_indices]
                    meta_val_dtype = np.array([meta_val[0]]).dtype

                    # some extra str formatting
                    if meta_var in ['station_reference', 'station_name', 'station_classification', 
                                    'area_classification']:

                        if len(meta_shape) == 2:
                            if meta_val_dtype == np.dtype(object):
                                meta_val = np.array([''.join(val) for val in meta_val])
                            else:
                                meta_val = chartostring(meta_val)

                # GHOST metadata
                else:
                    meta_var_nc = meta_var
                    meta_val = ncdf_root[meta_var_nc][current_file_station_indices]

                # put metadata in array
                file_metadata[meta_var][full_array_station_indices, 0] = meta_val

    # experiment data
    else:
        relevant_data = ncdf_root[speci][current_file_station_indices, valid_file_time_indices]

        # mask out fill values for parameter field
        relevant_data[relevant_data.mask] = np.NaN
        
        # put data in array
        data_in_memory[data_labels.index(data_label), full_array_station_indices[:, np.newaxis], 
                       full_array_time_indices[np.newaxis, :]] = relevant_data

    # close netCDF
    ncdf_root.close()

    # return metadata if reading observations
    if (data_label == 'observations') & (not filter_read):
        return file_metadata


def read_netcdf_metadata(tuple_arguments):

    """ Function that handles reading of observational basic metadata from a netCDF"""

    # assign arguments from tuple to variables
    relevant_file, reading_ghost = tuple_arguments

    # read netCDF frame
    ncdf_root = Dataset(relevant_file)

    # set metadata variables to read
    metadata_vars_to_read = ['station_reference', 'station_name', 'longitude', 'latitude', 'measurement_altitude']
    metadata_read = []

    # iterate though metadata variables to read
    for meta_var in metadata_vars_to_read:
        # do extra work for non-GHOST data 
        if not reading_ghost:
            # get correct variable name for .nc
            if meta_var == 'longitude':
                if "longitude" not in ncdf_root.variables:
                    meta_var_nc = 'lon'
                else:
                    meta_var_nc = 'longitude'
            elif meta_var == 'latitude':
                if "latitude" not in ncdf_root.variables:
                    meta_var_nc = 'lat'
                else:
                    meta_var_nc = 'latitude'
            elif meta_var == 'altitude':
                if "altitude" not in ncdf_root.variables:
                    meta_var_nc = 'alt'
                else:
                    meta_var_nc = 'altitude'
            elif meta_var == 'station_reference':
                if "station_reference" not in ncdf_root.variables:
                    meta_var_nc = 'station_code'
                else:
                    meta_var_nc = 'station_reference'
            else:
                meta_var_nc = meta_var

            # check meta variable is in netCDF, otherwise append empty list
            if meta_var_nc not in ncdf_root.variables:
                metadata_read.append([])
                continue

            meta_shape = ncdf_root[meta_var_nc].shape
            meta_val = ncdf_root[meta_var_nc][:]
            meta_val_dtype = np.array([meta_val[0]]).dtype

            # some extra str formatting
            if meta_var in ['station_reference', 'station_name', 'station_classification', 
                            'area_classification']:

                if len(meta_shape) == 2:
                    if meta_val_dtype == np.dtype(object):
                        meta_val = np.array([''.join(val) for val in meta_val])
                    else:
                        meta_val = chartostring(meta_val)

        # GHOST metadata
        else:
            meta_val = ncdf_root[meta_var][:]

        # append read metadata 
        metadata_read.append(meta_val)

    # close netCDF
    ncdf_root.close()

    # return read metadata
    return metadata_read

def get_yearmonths_to_read(available_yearmonths, start_date_to_read, end_date_to_read, resolution):
    """ Function that returns the yearmonths of the files to be read.
        Filters out yearmonths outside given date range.
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
    """ Return the default values according to GHOST standards. 

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :return: QA flags' codes in list
        :rtype: list
    """

    if instance.parameter_dictionary[speci]['extreme_lower_limit'] < 0.0:
        return sorted(instance.default_qa_non_negative)
    else:
        return sorted(instance.default_qa_standard)
    
def get_frequency_code(resolution):
    """ Get frequency code. """
    
    if resolution in ['hourly', 'hourly_instantaneous']:
        active_frequency_code = 'H'
    elif resolution in ['3hourly', '3hourly_instantaneous']:
        active_frequency_code = '3H'
    elif resolution in ['6hourly', '6hourly_instantaneous']:
        active_frequency_code = '6H'
    elif resolution == 'daily':
        active_frequency_code = 'D'
    elif resolution == 'monthly':
        active_frequency_code = 'MS'
    elif resolution == 'yearly':
        active_frequency_code = 'AS'

    return active_frequency_code