""" Module storing static reading functions """

import datetime
from glob import glob
import json
import os
import sys
import time

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
    observations_data_label, data_label, data_labels, reading_ghost, ghost_data_vars_to_read,\
    metadata_dtype, metadata_vars_to_read, logger, default_qa, filter_read, network = tuple_arguments

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
    
    # time_calendar = ncdf_root['time'].calendar
    if 'months' in time_units:
        monthly_start_date = time_units.split(' ')[2]
        file_time_dt = pd.date_range(start=monthly_start_date, periods=1, freq='MS')
    else:
        # file_time_dt = num2date(file_time, units=time_units, calendar=time_calendar)
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
                                        qa).any(axis=2)] = np.NaN
                
            # if some data provider flags selected then screen observations
            if len(flags) > 0:
                # screen out observations which are associated with any of the selected data provider flags
                species_data[np.isin(ncdf_root['flag'][current_file_station_indices, valid_file_time_indices, :], 
                                    flags).any(axis=2)] = np.NaN

        # write filtered species data to shared file data
        data_in_memory[data_labels.index(observations_data_label), full_array_station_indices[:, np.newaxis], 
                       full_array_time_indices[np.newaxis, :]] = species_data

        # get file metadata
        if not filter_read:
            file_metadata = np.full((len(station_references), 1), np.NaN, dtype=metadata_dtype)
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
    if (data_label == observations_data_label) & (not filter_read):
        return file_metadata


def read_netcdf_metadata(tuple_arguments):

    """ Function that handles reading of observational basic metadata from a netCDF"""

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
                if "measurement_altitude" in ncdf_root.variables:
                    meta_val = ncdf_root['measurement_altitude'][non_nan_station_indices]
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
    """ Return the default qa flags according to GHOST standards. 

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


def check_for_ghost(network_name):
    """ Check whether the selected network comes from GHOST or not.
        All non-GHOST networks start with an asterisk at their name.
    """

    if '/' in network_name:
        return False
    else:
        return True


def get_ghost_observational_tree(instance):
    """ Create GHOST observational data tree as a nested dictionary,
        storing a list of start YYYYMM yearmonths per:
        network / resolution / matrix / speci

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :return: GHOST observational tree dictionary
        :rtype: dict
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
    """ Fill non-GHOST observational data tree,
        storing a list of start YYYYMM yearmonths per:
        network / resolution / matrix / speci

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :return: non-GHOST observational tree dictionary
        :rtype: dict
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
    """ Iterate through observational dictionary tree and return 
        a dictionary of available data in the selected daterange

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :param start_date: start date (e.g. "20201101")
        :type start_date: str
        :param end_date: end date (e.g. "20201101")
        :type end_date: str
        :return: available observational tree dictionary
        :rtype: dict
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


def get_valid_experiments(instance, start_date, end_date, resolution, networks, species):
    """ Get valid experiments for daterange, and selected parameters.
        Update experiment pop-up with valid experiments.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :param start_date: start date (e.g. "20201101")
        :type start_date: str
        :param end_date: end date (e.g. "20201231")
        :type end_date: str
        :param resolution: resolution (e.g. "hourly")
        :type resolution: str
        :param networks: list of networks
        :type networks: list
        :param species: list of species
        :type species: list 
    """

    # get all different experiment names (from providentia-interpolation output dir)
    available_experiments = []
    if os.path.exists(join(instance.exp_root,instance.ghost_version)):
        available_experiments = os.listdir('%s/%s' % (instance.exp_root, instance.ghost_version))

    # create dictionary to store available experiment data
    instance.available_experiment_data = {}

    #list for saving experiments to add to experiments pop-up 
    experiments_to_add = []

    # get start date on first of month
    start_date_firstdayofmonth = int(str(start_date)[:6] + '01')

    # iterate through networks and species
    for network, speci in zip(networks, species):

        # iterate through available experiments
        for experiment in available_experiments:

            # get folder where interpolated experiments are saved
            if '/' not in network:           
                files_directory = '%s/%s/%s/%s/%s/%s' % (instance.exp_root, instance.ghost_version, 
                                                         experiment, resolution, speci, network)
            else:
                files_directory = '%s/%s/%s/%s/%s/%s' % (instance.exp_root, instance.ghost_version, 
                                                          experiment, resolution, speci,
                                                          network.replace('/', '-'))
                
            # test if interpolated directory exists for experiment
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

            # write nested dictionary for experiment, with associated file yearmonths
            if len(file_yearmonths) > 0:

                # get file yearmonths within date range
                valid_file_yearmonths = sorted([ym for ym in file_yearmonths if 
                                                (int('{}01'.format(ym)) >= start_date_firstdayofmonth) & (int('{}01'.format(ym)) < int(end_date))])

                #if have valid files, then add experiment to pop-up menu, 
                # and add yearmonths to available experiment data
                if len(valid_file_yearmonths) > 0:
                    experiments_to_add.append(experiment)

                    if network not in instance.available_experiment_data:
                        instance.available_experiment_data[network] = {}
                    if resolution not in instance.available_experiment_data[network]:
                        instance.available_experiment_data[network][resolution] = {}
                    if speci not in instance.available_experiment_data[network][resolution]:
                        instance.available_experiment_data[network][resolution][speci] = {}
                    if experiment not in instance.available_experiment_data[network][resolution][speci]:
                        instance.available_experiment_data[network][resolution][speci][experiment] = valid_file_yearmonths

    # set list of experiment names to add on experiments pop-up
    if (not instance.offline) and (not instance.interactive):
        experiments_to_add = np.array(sorted(experiments_to_add))
        instance.experiments_menu['checkboxes']['labels'] = experiments_to_add
        instance.experiments_menu['checkboxes']['map_vars'] = experiments_to_add


def temporal_resolution_order_dict():
    """ Return temporal resolution order for menus as a dictionary.

        :return: temporal resolution order
        :rtype: dict
    """

    resolution_order_dict = {'hourly': 1, '3hourly': 2, '6hourly': 3, 'hourly_instantaneous': 4,
                             '3hourly_instantaneous': 5, '6hourly_instantaneous': 6,
                             'daily': 7, 'monthly': 8}

    return resolution_order_dict


def get_relevant_temporal_resolutions(resolution):        
    """ Get relevant temporal reolsutions for periodic plots, by selected temporal resolution.

        :param resolution: name of selected temporal resolution
        :type resolution: str
        :return: relevant temporal resolutions
        :rtype: list
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

def get_nonrelevant_temporal_resolutions(resolution):        
    """ Get non-relevant temporal reolsutions for periodic plots, by selected temporal resolution.

        :param resolution: name of selected temporal resolution
        :type resolution: str
        :return: non-relevant temporal resolutions
        :rtype: list
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
    """ Determine if a date string is in the correct format. """

    try:
        datetime.datetime.strptime(str(date_text), '%Y%m%d')
        return True
    except Exception as e:
        return False


def get_lower_resolutions(resolution):
    """ Get available lower resolutions. """

    if resolution in ['hourly', 'hourly_instantaneous', '3hourly', '3hourly_instantaneous', 
                      '6hourly', '6hourly_instantaneous']:
        resolutions = ['daily', 'monthly', 'annual']
    elif resolution == 'daily':
        resolutions = ['monthly', 'annual']
    elif resolution == 'monthly':
        resolutions = ['annual']
    elif resolution == 'annual':
        resolutions = []

    return resolutions


def generate_file_trees(instance, force=False):
    """ Generate file trees. Force if we want to remove depedency on the machine.
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
    if gft or force:
        instance.logger.info('Generating file trees...')
        instance.all_observation_data = get_ghost_observational_tree(instance)
        if instance.nonghost_root is not None:
            nonghost_observation_data = get_nonghost_observational_tree(instance)
    # load file trees
    else:
        try:
            instance.all_observation_data = json.load(open(join(PROVIDENTIA_ROOT, 'settings/internal/ghost_filetree_{}.json'.format(instance.ghost_version)))) 
        except FileNotFoundError as file_error:
            error = "Error: Trying to load 'settings/internal/ghost_filetree_{}.json' but file does not exist. Run with the flag '--gft' to generate this file.".format(instance.ghost_version)
            instance.logger.error(error)
            sys.exit(1)
        if instance.nonghost_root is not None:
            try:
                nonghost_observation_data = json.load(open(join(PROVIDENTIA_ROOT, 'settings/internal/nonghost_filetree.json')))
            except FileNotFoundError as file_error:
                error = "Error: Trying to load 'settings/internal/nonghost_filetree.json' but file does not exist. Run with the flag '--gft' to generate this file."
                instance.logger.error(error)
                sys.exit(1)

    # merge GHOST and non-GHOST filetrees
    if instance.nonghost_root is not None:
        instance.all_observation_data = {**instance.all_observation_data, **nonghost_observation_data}


def get_valid_metadata(read_instance, variable, valid_station_idxs, networkspeci):
    """ Get metadata without nans from the data for each timestep
        If for all timesteps the metadata is nan, set nan as valid metadata
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