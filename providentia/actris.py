import bisect
import copy
import csv
import ctypes
import datetime
import itertools
import os
import requests
import sys
import yaml
import re
import time
from tqdm import tqdm
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
import multiprocessing
from dateutil.relativedelta import relativedelta

from providentia.auxiliar import CURRENT_PATH, join, pad_array
from .warnings_prv import show_message

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

# load ACTRIS mapping files
ghost_actris_variables = yaml.safe_load(open(join(
    PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'ghost_actris_variables.yaml')))
metadata_dict = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'metadata.yaml')))
coverages_dict = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'coverages.yaml')))
variable_mapping = yaml.safe_load(open(join(
    PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'variable_mapping.yaml')))
variable_mapping = {k: v for k,
                    v in variable_mapping.items() if k.strip() and v}
flags_dict = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'flags.yaml')))

# initialise dictionary for storing pointers to shared memory variables in read step 
shared_memory_vars = {}

def create_variable_mapping_file():

    result = {
        value['preferred_term'].replace('"', ''): {'var': key[2], 'units': key[0]}
        for key, value in variable_mapping.items()
    }
    with open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'variable_mapping.yaml'),
              mode='w') as file:
        yaml.dump(result, file, default_flow_style=False)


def create_actris_variables_file():

    with open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'actris_variables.csv'),
              mode="w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        for key in variable_mapping.keys():
            writer.writerow([key, variable_mapping[key]['var']])


def create_ghost_variables_file(ghost_version):

    sys.path.insert(1, join(
        PROVIDENTIA_ROOT, 'providentia/dependencies/GHOST_standards/{}'.format(ghost_version)))
    from GHOST_standards import standard_parameters
    with open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'ghost_variables.csv'),
              mode="w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        for key in standard_parameters.keys():
            writer.writerow([standard_parameters[key]['long_parameter_name'],
                             standard_parameters[key]['bsc_parameter_name'],
                             ', '.join(standard_parameters[key]['ebas_parameter_name'])])


def get_files_per_var(download_instance, var):
    """Get all files available in ACTRIS server per variable

    Parameters
    ----------
    var : str
        Variable

    Returns
    -------
    dict
        Dictionary with files per variable
    """

    files_per_var = {}
    base_url = "https://prod-actris-md.nilu.no/metadata/content"

    if var not in files_per_var:
        files_per_var[var] = {}

    variable_files = []
    page = 0
    while True:
        # set up URL with pagination
        url = f"{base_url}/{ghost_actris_variables[var]}/page/{page}"
        response = requests.get(url)

        # check if the response is valid and contains data
        if response.status_code != 200:
            download_instance.logger.error(
                f"Error fetching page {page}. Status code: {response.status_code}")
            break

        data = response.json()

        # check if there's content in the data
        if not data:
            break

        # loop through each entry in the data and get OPeNDAP URL
        for item in data:
            doi = item.get("md_identification", {}).get(
                "identifier", {}).get("pid")
            opendap_urls = [protocol_dict['dataset_url'] for protocol_dict in item.get(
                'md_distribution_information', []) if protocol_dict.get('protocol') == 'OPeNDAP']

            # print DOI and OPeNDAP URL if both are present
            if doi and opendap_urls:
                variable_files.append(opendap_urls)

        # go to the next page
        page += 1

    files_per_var[var]['files'] = list(
        itertools.chain.from_iterable(variable_files))

    return files_per_var


def get_files_path(var):
    """Get path of file where files per variable are saved

    Parameters
    ----------
    var : str
        Variable

    Returns
    -------
    str
        Path where files per variable are saved
    """

    alpha_var = ''.join(x for x in var if x.isalpha())
    if alpha_var in ['lsco', 'absco', 'lbsco', 'odaero']:
        path = join(PROVIDENTIA_ROOT, 'settings', 'internal',
                    'actris', f'files/{alpha_var}.yaml')
    else:
        path = join(PROVIDENTIA_ROOT, 'settings', 'internal',
                    'actris', f'files/{var}.yaml')

    return path


def get_standard_flags_and_qa(flags, ghost_version):
    """Convert flags from EBAS standards to GHOST standards and get QA

    Parameters
    ----------
    flags : numpy.array
        Array with flags in EBAS standards
    ghost_version : float
        GHOST version
    """

    sys.path.insert(1, join(
        PROVIDENTIA_ROOT, 'providentia/dependencies/GHOST_standards/{}'.format(ghost_version)))
    from GHOST_standards import standard_data_flag_name_to_data_flag_code

    qa = []
    standard_flags = []
    network_decreed_validities = []
    GHOST_decreed_validities = []
    for flag in flags:
        # get standard flag name from GHOST standards
        flag_info = flags_dict[str(int(flag))] if flag != 0 else flags_dict['000']
        standard_flag_name = flag_info["standard_data_flag_name"][0]
        standard_flag = standard_data_flag_name_to_data_flag_code[standard_flag_name]
        standard_flags.append(standard_flag)

        # save validities
        network_decreed_validity = flag_info["network_decreed_validity"][0]
        GHOST_decreed_validity = flag_info["GHOST_decreed_validity"][0]
        network_decreed_validities.append(network_decreed_validity)
        GHOST_decreed_validities.append(GHOST_decreed_validity)

    # get qa if there is any flag that is invalid
    if 'I' in np.unique(GHOST_decreed_validities):
        qa.append(6)  # 'Invalid Data Provider Flags - GHOST Decreed'
    if 'I' in np.unique(network_decreed_validities):
        qa.append(7)  # 'Invalid Data Provider Flags - Network Decreed'

    return np.array(standard_flags, dtype=np.float32), np.array(qa, dtype=np.float32)


def create_time_pairs(time):

    time_pairs = list(zip(time[:-1], time[1:]))
    last_time_pair = (time[-1], time[-1] + (time[-1] - time[-2]))
    time_pairs.append(last_time_pair)

    return time_pairs


def datetime_to_fractional_minutes_from_reference(date):
    """Get datetime in fractional minutes

    Parameters
    ----------
    date : datetime.datetime
        Datetime

    Returns
    -------
    int
        Datetime in fractional minutes
    """

    past_date = datetime.datetime(1,1,1,0,0)
    difference = date - past_date
    
    return difference.total_seconds() / datetime.timedelta(minutes=1).total_seconds()
                    

def get_window_indices(standard_start_date, standard_end_date, valid_start_times, valid_end_times, last_relevant_index):

    window_start_minute = datetime_to_fractional_minutes_from_reference(standard_start_date)
    window_end_minute = datetime_to_fractional_minutes_from_reference(standard_end_date)

    # get index where valid start times >= window start minute
    valid_start_times_begin = bisect.bisect_left(valid_start_times, window_start_minute, 
                                                    lo=last_relevant_index, hi=len(valid_start_times))
    
    # get index where valid start times < window end minute
    valid_start_times_end = bisect.bisect_left(valid_start_times, window_end_minute, 
                                                lo=last_relevant_index, hi=len(valid_start_times))

    # get index of first start time < window start minute
    start_time_index_before_window = valid_start_times_begin - 1

    # get index where valid end times > window start minute
    valid_end_times_begin = bisect.bisect_right(valid_end_times, window_start_minute, 
                                                lo=last_relevant_index, hi=len(valid_end_times))
    
    # get index where valid end times <= window end minute
    valid_end_times_end = bisect.bisect_right(valid_end_times, window_end_minute, 
                                                lo=last_relevant_index, hi=len(valid_end_times))

    # get index of first end time > window end minute
    end_time_index_after_window = copy.deepcopy(valid_end_times_end)

    # get indices of all measurement periods which entirely contain window
    if start_time_index_before_window == end_time_index_after_window:
        window_in_measurement_indices = np.array([start_time_index_before_window])
    else:
        window_in_measurement_indices = np.array([], dtype=np.uint32)

    # get indices of measurements entirely contained within window
    start_time_indices_in_window = np.arange(valid_start_times_begin, valid_start_times_end, dtype=np.uint32)
    end_time_indices_in_window = np.arange(valid_end_times_begin, valid_end_times_end, dtype=np.uint32)
    measurement_in_window_indices = np.intersect1d(start_time_indices_in_window,end_time_indices_in_window, assume_unique=True)
    
    # get indices of measurements which overlap on left/right edges of window
    left_overlap_indices = np.setdiff1d(end_time_indices_in_window, measurement_in_window_indices, assume_unique=True)
    left_overlap_indices = np.setdiff1d(left_overlap_indices, window_in_measurement_indices, assume_unique=True)
    right_overlap_indices = np.setdiff1d(start_time_indices_in_window, measurement_in_window_indices, assume_unique=True)
    right_overlap_indices = np.setdiff1d(right_overlap_indices, window_in_measurement_indices, assume_unique=True)

    # deal with cases where left/right borders align but measurement period entirely contains window
    # add indices of these cases to window_in_measurement_indices
    # remove these indices also from the overlap indices
    left_window_in_measurement_indices = right_overlap_indices[np.where(valid_start_times[right_overlap_indices] == window_start_minute)]
    right_window_in_measurement_indices = left_overlap_indices[np.where(valid_end_times[left_overlap_indices] == window_end_minute)]
    if len(left_window_in_measurement_indices) > 0:
        window_in_measurement_indices = np.concatenate((window_in_measurement_indices, left_window_in_measurement_indices))
        right_overlap_indices = np.setdiff1d(right_overlap_indices, left_window_in_measurement_indices, assume_unique=True)
    if len(right_window_in_measurement_indices) > 0:
        window_in_measurement_indices = np.concatenate((window_in_measurement_indices, right_window_in_measurement_indices))
        left_overlap_indices = np.setdiff1d(left_overlap_indices, right_window_in_measurement_indices, assume_unique=True)
    
    # if there is a left border overlap, get the number of minutes the measurement period overlaps the measurement window
    if len(left_overlap_indices) > 0:
        left_overlap = valid_end_times[left_overlap_indices[0]] - window_start_minute
    # otherwise left overlap == 0
    else:
        left_overlap = 0
    # if there is a right border overlap, get the number of minutes the measurement period overlaps the measurement window
    if len(right_overlap_indices) > 0:
        right_overlap = window_end_minute - valid_start_times[right_overlap_indices[0]]
    # otherwise right overlap == 0
    else:
        right_overlap = 0
    # concatenate all relevant measurements indices in current window
    window_indices = np.sort(np.concatenate((measurement_in_window_indices,window_in_measurement_indices,left_overlap_indices,right_overlap_indices)))
            
    return window_indices, right_overlap, left_overlap


def temporally_average_data(station_ds, var, ghost_version, standard_time_pairs, vfunc):
    """Temporally average data and get unique flags in the valid times (temporally averaged)

    Parameters
    ----------
    ds : xarray.Dataset
        Data per station
    resolution : str
        Temporal resolution
    var : str
        Variable
    ghost_version : float
        GHOST version

    Returns
    -------
    ds : xarray.Dataset
        Data per station with valid times per month
    """

    # initialise averaged data
    station_averaged_data = []
    station_flag_data = []
    station_qa_data = []

    # remove station selection
    station_var = station_ds[var].values
    station_flag = station_ds['flag'].values
    station_time_bnds = station_ds['time_bnds'].values

    # get measurement start and end times
    start_times = station_time_bnds[:, 0]
    end_times = station_time_bnds[:, 1]

    # get timedelta between start and end times
    valid_timedeltas =  np.array([(end_time - start_time).astype('timedelta64[m]').astype(np.float32) 
                                  for (end_time, start_time) in zip(end_times, start_times)])
    
    # get measurement start and end times as integers
    valid_start_times = np.array([datetime_to_fractional_minutes_from_reference(t) 
                                  for t in pd.to_datetime(start_times).to_pydatetime()])
    valid_end_times = np.array([datetime_to_fractional_minutes_from_reference(t) 
                                for t in pd.to_datetime(end_times).to_pydatetime()])

    # initialise variable for finding relevant measurement indices
    last_relevant_index = 0

    for i, (standard_start_date, standard_end_date) in enumerate(standard_time_pairs):
        
        # get window indices
        window_indices, right_overlap, left_overlap = get_window_indices(standard_start_date, 
                                                                         standard_end_date, 
                                                                         valid_start_times, 
                                                                         valid_end_times, 
                                                                         last_relevant_index)

        # only one overlap value
        if len(window_indices) == 1:
            window_var_data = station_var[window_indices][0]
            flag_data = station_flag[window_indices][0]

            # get unique flag values and convert to standard flag names (instead of EBAS) and qa
            flags = np.unique(flag_data)
            valid_flags = flags[~np.isnan(flags)]
            window_flag_data, window_qa_data = get_standard_flags_and_qa(
                valid_flags, ghost_version)
            
            # record last index of current window indices for next iteration
            last_relevant_index = window_indices[-1]

        # multiple overlap values
        elif len(window_indices) > 1:
            flag_data = station_flag[window_indices]
            var_data = station_var[window_indices]

            GHOST_decreed_validities = vfunc(flag_data)
            GHOST_invalid = np.any(
                GHOST_decreed_validities == 'I', axis=1)

            # if there are invalid values in period but some of them are valid, convert invalid ones to nan
            if np.any(GHOST_invalid) and not np.all(GHOST_invalid):
                var_data[GHOST_invalid] = np.nan
                flag_data[GHOST_invalid, :] = np.nan
            
            # weight data by timedeltas for averaging different temporal resolution data
            window_weights = valid_timedeltas[window_indices]

            # modify weights on left and right edges to adjust for actual time sampled within window, 
            # if measurement overlaps edge
            if left_overlap > 0.0:
                window_weights[0] = left_overlap
            if right_overlap > 0.0:
                window_weights[-1] = right_overlap

            window_var_data = np.average(var_data, weights=window_weights)

            # get unique flag values and convert to standard flag names (instead of EBAS) and qa
            flags = np.unique(flag_data)
            valid_flags = flags[~np.isnan(flags)]
            window_flag_data, window_qa_data = get_standard_flags_and_qa(
                valid_flags, ghost_version)
            
            # record last index of current window indices for next iteration
            last_relevant_index = window_indices[-1]

        # no overlap values
        else:
            window_var_data = np.nan
            window_flag_data = [np.nan]
            window_qa_data = [np.nan]

        # save flags, qa and var
        # TODO: Check why flag data is not sorted
        station_flag_data.append(
            pad_array(window_flag_data, length=186))
        station_qa_data.append(pad_array(window_qa_data, length=2))
        station_averaged_data.append(window_var_data)
    
    return station_averaged_data, station_flag_data, station_qa_data


def is_wavelength_var(actris_parameter):
    """Check if ACTRIS parameter depends on wavelength

    Parameters
    ----------
    actris_parameter : str
        ACTRIS parameter

    Returns
    -------
    bool
        Indicates if ACTRIS parameter depends on wavelength
    """

    wavelength_var = False
    if actris_parameter in ['aerosol particle light absorption coefficient',
                            'aerosol particle optical depth',
                            'aerosol particle light hemispheric backscatter coefficient',
                            'aerosol particle light scattering coefficient',
                            'aerosol particle equivalent black carbon mass concentration']:
        wavelength_var = True

    return wavelength_var


def get_files_info(download_instance, files, var, path):
    """Read variables, resolution, start date and end date from all files in ACTRIS server per variable.

    Parameters
    ----------
    files : list
        All files per variable
    var : str
        Variable
    path : str
        Path where files per variable are saved

    Returns
    -------
    dict
        Dictionary with details per file
    """

    files_info = {}
    tqdm_iter = tqdm(files, bar_format='{l_bar}{bar}|{n_fmt}/{total_fmt}',
                     desc=f"Creating information file")
    for file in tqdm_iter:
        # open file
        try:
            ds = Dataset(file)
        except:
            continue

        # get resolution
        coverage = ds.time_coverage_resolution
        try:
            file_resolution = coverages_dict[coverage]
        except:
            file_resolution = f'Unrecognised ({coverage})'

        # save in dict
        files_info[file] = {}
        files_info[file]['resolution'] = file_resolution
        files_info[file]['variables'] = list(ds.variables.keys())
        for var in ['time_coverage_start', 'time_coverage_end', 'ebas_statistics', 
                    'ebas_station_code', 'ebas_station_latitude', 'ebas_station_longitude',
                    'ebas_data_level', 'ebas_station_altitude', 'ebas_measurement_height',
                    'ebas_instrument_name', 'ebas_method_ref', 'ebas_revision_date',
                    'ebas_station_code']:
            if var in ds.ncattrs():
                files_info[file][var] = ds.getncattr(var)

        if "Wavelength" in files_info[file]['variables']:
            wavelengths = ds.variables["Wavelength"][:]
            wavelengths = [float(w) for w in wavelengths]
            files_info[file]['wavelengths'] = wavelengths

    # create file
    datasets = {
        url: data
        for url, data in files_info.items()
    }
    if len(datasets) != 0:
        path_dir = os.path.dirname(path)
        if not os.path.exists(path_dir):
            os.makedirs(path_dir)
        with open(path, 'w') as file:
            yaml.dump(datasets, file, default_flow_style=False)
    else:
        download_instance.logger.error(f'Error: No data could be found for {var}')
        return
    
    return files_info


def get_var_in_file(ds, var, actris_parameter, ebas_component):

    unformatted_units = variable_mapping[actris_parameter]['units']
    units = unformatted_units.replace('/', '_per_').replace(' ', '_')
    units_var = f'{ebas_component}_{units}'
    possible_vars = [ebas_component,
                     f'{ebas_component}_amean',
                     units_var,
                     f'{units_var}_amean']
    if var in ['sconcso4', 'precso4']:
        possible_vars.append(f'sulphate_corrected_{units}')
    da_var_exists = False
    for possible_var in possible_vars:
        if possible_var in ds:
            da_var_exists = True
            break

    # continue to next file if variable cannot be read
    if not da_var_exists:
        return None

    return possible_vars, possible_var

def select_station_file(urls, files_info):

    attrs_dict = {}
    urls = np.array(urls)
    # create dictionary with information about statistics, data level and revision date of available files
    for url in urls:
        attrs = files_info[url]
        for attr in ['ebas_statistics', 'ebas_data_level', 'ebas_revision_date']:
            if attr not in attrs_dict:
                attrs_dict[attr] = np.array([])
            if attr in attrs:
                attrs_dict[attr] = np.append(attrs_dict[attr], attrs[attr])
            else:
                attrs_dict[attr] = np.append(attrs_dict[attr], '')

    urls_statistics = np.unique(attrs_dict['ebas_statistics'])
    if len(urls_statistics) > 1:
        for attr_val in ['arithmetic mean', '']:
            if attr_val in attrs_dict['ebas_statistics']:
                is_attr_val = attrs_dict['ebas_statistics'] == attr_val

                # remove urls if the statistics are not attr_val (arithmetic mean or undefined)
                for attr in ['ebas_statistics', 'ebas_data_level', 'ebas_revision_date']:
                    attrs_dict[attr] = attrs_dict[attr][is_attr_val]
                urls = urls[is_attr_val]
                break
    
    # if all the urls contain statistics that are not arithmetic mean or undefined, 
    # then continue to next station
    elif len(urls_statistics) == 1:
        if urls_statistics[0] not in ['', 'arithmetic mean']:
            return

    # check if we still have files
    if len(urls) > 1:
        urls_data_levels = np.unique(attrs_dict['ebas_data_level'])
        # if we have different data levels
        if len(urls_data_levels) > 1:
            max_level_ind = np.nanargmax(np.float32([np.nan if x == '' else x for x in attrs_dict['ebas_data_level']]))
            is_max = attrs_dict['ebas_data_level'] == attrs_dict['ebas_data_level'][max_level_ind]
            
            # remove urls if the data level is not the maximum of all files
            for attr in ['ebas_statistics', 'ebas_data_level', 'ebas_revision_date']:
                attrs_dict[attr] = attrs_dict[attr][is_max]
            urls = urls[is_max]

    # check if we still have files
    if len(urls) > 1:
        urls_revision_dates = np.unique(attrs_dict['ebas_revision_date'])
        # if we have different revision dates
        if len(urls_revision_dates) > 1:
            is_most_recent = attrs_dict['ebas_revision_date'] == attrs_dict['ebas_revision_date'][np.argmax(np.float32(attrs_dict['ebas_revision_date']))]
            
            # remove urls if they aren't the most recent
            for attr in ['ebas_statistics', 'ebas_data_level', 'ebas_revision_date']:
                attrs_dict[attr] = attrs_dict[attr][is_most_recent]
            urls = urls[is_most_recent]               

    # get first file after checks
    # in case we have multiple files, that would mean all our remaining files have the same revision date, 
    # data level and statistics
    url = urls[0]

    return url


def read_data(args):

    (i, station, urls, data_shape, flag_shape, qa_shape,
     metadata_shape, files_info, var, actris_parameter, ebas_component,
     target_start_date, target_end_date,
     variable_mapping, metadata_dict, standard_time_pairs,
     vfunc, ghost_version) = args

    local_errors = ""
    local_warnings = ""

    # open files
    try:
        if len(urls) > 1:
            url = select_station_file(urls, files_info)
            if url is None:
                return station, local_errors, local_warnings
        else:
            url = urls[0]

        nc = Dataset(url, mode='r')
        ds = xr.open_dataset(xr.backends.NetCDF4DataStore(nc))

        possible_vars, possible_var = get_var_in_file(ds, var, actris_parameter, ebas_component)
        if possible_var is None:
            local_errors = f'No variable name matches for {possible_vars}. Existing keys: {list(ds.data_vars)}.'
            return station, local_errors, local_warnings

    except Exception as error:
        local_errors = f'Opening file: {error}.'
        return station, local_errors, local_warnings
  
    # remove time duplicates if any (keep first)
    ds = ds.sel(time=~ds['time'].to_index().duplicated())
    
    # select data in period range
    ds = ds.sel(time=slice(target_start_date, target_end_date))
    if ds.time.size == 0:
        local_warnings += f'No data available after filtering by time.'
        return station, local_errors, local_warnings

    # rename qc dimension
    ds = ds.rename_dims({f'{possible_var}_qc_flags': 'N_flag_codes'})

    # get lowest level if tower height is in coordinates
    if 'Tower_inlet_height' in list(ds.coords):
        local_warnings += f'Taking data from first height (Tower_inlet_height={min(ds.Tower_inlet_height.values)}). '
        ds = ds.sel(Tower_inlet_height=min(
            ds.Tower_inlet_height.values), drop=True)

    # get data at desired wavelength if wavelength is in coordinates
    if 'Wavelength' in list(ds[possible_var].coords) or 'Wavelengthx' in list(ds[possible_var].coords):
        # Select most common wavelength for black carbon (name does not provide it)
        if var == 'sconcbc':
            wavelength = 880
            local_warnings += f'Wavelength appears in dimensions. Selected wavelength: {wavelength}. '
        # Get wavelength from variable name for other variables
        else:
            wavelength = float(re.findall(r'\d+', var)[0])

        # Select data for wavelength
        found_wavelength = False
        if 'Wavelengthx' in list(ds[possible_var].coords):
            if wavelength in ds.Wavelengthx.values:
                ds = ds.sel(Wavelengthx=wavelength, drop=True)
                found_wavelength = True
            else:
                existing_wavelengths = ds.Wavelengthx.values
        elif 'Wavelength' in list(ds[possible_var].coords):
            if wavelength in ds.Wavelength.values:
                ds = ds.sel(Wavelength=wavelength, drop=True)
                found_wavelength = True
            else:
                existing_wavelengths = ds.Wavelength.values

        if not found_wavelength:
            local_warnings += f'Data at {wavelength}nm could not be found. Existing wavelengths: {existing_wavelengths}. '
            return station, local_errors, local_warnings

    # remove artifact and fraction (sconcoc)
    # TODO: Discuss this
    if 'Artifact' in list(ds.coords):
        local_warnings += f'Taking data from first artifact dimension (Artifact={ds.Artifact.values[0]}). '
        ds = ds.isel(Artifact=0, drop=True)
    if 'Fraction' in list(ds.coords):
        local_warnings += f'Taking data from first fraction dimension (Fraction={ds.Fraction.values[0]}). '
        ds = ds.isel(Fraction=0, drop=True)

    # read variable
    da_var = ds[possible_var]

    # avoid datasets that do not have defined units
    if 'ebas_unit' not in da_var.attrs:
        local_errors = f'No units were defined.'
        return station, local_errors, local_warnings

    # avoid datasets that do not have the same units as in variable mapping
    if da_var.attrs['ebas_unit'] != variable_mapping[actris_parameter]['units']:
        local_errors = f"Units {da_var.attrs['ebas_unit']} do not match those in variable mapping "
        local_errors += f"dictionary ({variable_mapping[actris_parameter]['units']})."
        return station, local_errors, local_warnings

    # remove all attributes except units
    da_var_attrs = copy.deepcopy(da_var.attrs)
    da_var.attrs = {key: value for key,
                    value in da_var.attrs.items() if key in ['ebas_unit', 'ebas_station_code']}
    
    # read quality control data
    try:
        flag_data = ds[f'{possible_var}_qc'].transpose(
            "time", "N_flag_codes")
    except:
        local_errors = "Flag data could not be transposed."
        return station, local_errors, local_warnings

    # rename variable to BSC standards
    station_ds = da_var.to_dataset(name=var)

    # add time bounds
    station_ds['time_bnds'] = ds['time_bnds']

    # add quality control data
    station_ds['flag'] = flag_data

    # temporally average data from original times to standard times
    station_averaged_data, station_flag_data, station_qa_data = temporally_average_data(station_ds, var, ghost_version, standard_time_pairs, vfunc)
    if np.isnan(station_averaged_data).all():
        local_warnings += 'No data after temporal averaging.'
        return station, local_errors, local_warnings

    data_np = np.frombuffer(shared_memory_vars['data'], dtype=np.float32).reshape(data_shape)
    flag_np = np.frombuffer(shared_memory_vars['flag'], dtype=np.float32).reshape(flag_shape)
    qa_np = np.frombuffer(shared_memory_vars['qa'], dtype=np.float32).reshape(qa_shape)

    data_np[i, :] = station_averaged_data
    flag_np[i, :, :] = station_flag_data
    qa_np[i, :, :] = station_qa_data

    # save metadata
    metadata_np = {}
    for ghost_key, ebas_key in metadata_dict.items():
        metadata_np[ghost_key] = np.frombuffer(shared_memory_vars['metadata'][ghost_key], dtype='S75').reshape(metadata_shape)
        if ebas_key in da_var_attrs.keys():
            val =  da_var_attrs[ebas_key]
        elif ebas_key in ds.attrs.keys():
            val = ds.attrs[ebas_key]
        else:
            val = str(np.nan)
        metadata_np[ghost_key][i] = val.encode('utf-8')

    return station, local_errors, local_warnings


def init_shared_vars_read_data(shared_data, shared_flag_data, shared_qa_data, shared_metadata):
    shared_memory_vars['data'] = shared_data
    shared_memory_vars['flag'] = shared_flag_data
    shared_memory_vars['qa'] = shared_qa_data
    shared_memory_vars['metadata'] = {}
    for ghost_key in metadata_dict.keys():
        shared_memory_vars['metadata'][ghost_key] = shared_metadata[ghost_key]


def ghost_validity_flag_mapper(flag):
    if np.isnan(flag):
        return np.nan
    elif flag == 0:
        return flags_dict['000']["GHOST_decreed_validity"][0]
    else:
        return flags_dict[str(int(flag))]["GHOST_decreed_validity"][0]
    

def get_data(download_instance, files, var, actris_parameter, resolution, target_start_date, target_end_date, 
             files_info, ghost_version, n_cpus):
    """Read variable and metadata data standarising dimensions

    Parameters
    ----------
    files : list
        Files per variable
    var : str
        Variable
    actris_parameter : str
        ACTRIS parameter
    resolution : str
        Resolution
    target_start_date : datetime.datetime
        Target start date (defined from configuration file)
    target_end_date : datetime.datetime
        Target end date (defined from configuration file)

    Returns
    -------
    list
        Standarised variable data
    list
        Standarised metadata
    wavelength : None, float, int
        Wavelength
    """

    # get EBAS component
    ebas_component = variable_mapping[actris_parameter]['var']

    # get valid dates frequency
    if resolution == 'hourly':
        frequency = 'h'
    elif resolution == '3hourly':
        frequency = '3h'
    elif resolution == 'daily':
        frequency = 'D'
    # TODO: Review this
    elif resolution == 'monthly':
        frequency = 'MS'

    standard_time = pd.date_range(start=target_start_date, end=target_end_date,
                                  freq=frequency).to_pydatetime()
    
    # get dimension of new arrays
    n_stations = len(files)
    n_times = len(standard_time)
    data_shape = (n_stations, n_times)
    flag_shape = (n_stations, n_times, 186)
    qa_shape = (n_stations, n_times, 2)
    metadata_shape = (n_stations)

    # initialise averaged data
    averaged_data = np.full(data_shape, fill_value=np.nan, dtype=np.float32)
    averaged_flag_data = np.full(flag_shape, fill_value=np.nan, dtype=np.float32)
    averaged_qa_data = np.full(qa_shape, fill_value=np.nan, dtype=np.float32)

    # define vectorize function to get GHOST_decreed_validity by overlap flags later
    vfunc = np.vectorize(ghost_validity_flag_mapper, otypes=['O'])

    # create pairs of valid dates
    standard_time_pairs = create_time_pairs(standard_time)

    # create data specific arrays to share across processes (for parallel multiprocessing use)
    shared_data = multiprocessing.RawArray(ctypes.c_float, int(np.prod(data_shape)))
    shared_flag_data = multiprocessing.RawArray(ctypes.c_float, int(np.prod(flag_shape)))
    shared_qa_data = multiprocessing.RawArray(ctypes.c_float, int(np.prod(qa_shape)))
    shared_metadata = {}
    for ghost_key in metadata_dict.keys():
        shared_metadata[ghost_key] = multiprocessing.RawArray(ctypes.c_char, int(np.prod(metadata_shape)*75))

    # read data and metadata in parallel
    errors = []
    warnings = []
    args_list = [
        (
            i, station, urls, data_shape, flag_shape, qa_shape, metadata_shape,
            files_info, var, actris_parameter, ebas_component,
            target_start_date, target_end_date,
            variable_mapping, metadata_dict, standard_time_pairs,
            vfunc, ghost_version
        )
        for i, (station, urls) in enumerate(files.items())
    ]
    pool = multiprocessing.Pool(
        processes=n_cpus,
        initializer=init_shared_vars_read_data,
        initargs=(shared_data, shared_flag_data, shared_qa_data, shared_metadata)
    )
    for station, error, warning in tqdm(
            pool.imap(read_data, args_list),
            total=len(args_list),
            desc="Reading data",
            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}'
        ):
        if error:
            errors.append(f'{station} - Error: {error}')
        if warning:
            warnings.append(f'{station} - Warning: {warning}')

    # print errors and warnings if any
    if errors:
        download_instance.logger.info("=== ERRORS ===")
        for e in errors:
            download_instance.logger.info(e)

    if warnings:
        download_instance.logger.info("=== WARNINGS ===")
        for w in warnings:
            download_instance.logger.info(w)
            
    pool.close()
    
    # wait for worker processes to terminate before continuing
    pool.join()

    if (len(errors) + len(warnings)) == len(args_list):
        download_instance.logger.info('All datasets have thrown an error or warning, aborting.')
        return None, None

    # get combined data and metadata after read
    averaged_data = np.frombuffer(shared_data, dtype=np.float32).reshape(data_shape)
    averaged_flag_data = np.frombuffer(shared_flag_data, dtype=np.float32).reshape(flag_shape)
    averaged_qa_data = np.frombuffer(shared_qa_data, dtype=np.float32).reshape(qa_shape)
    metadata = {}
    for ghost_key in metadata_dict.keys():
        metadata[ghost_key] = np.frombuffer(shared_metadata[ghost_key], dtype='S75').reshape(metadata_shape)

    # create dataset with averaged data
    units = variable_mapping[actris_parameter]['units']
    
    # case for temperature (t2), unit converter assumes multiple units when there is a space
    if units == 'deg C':
        units = 'degC'
    # case for lsco, lbsco and absco
    elif units == '1/Mm':
        units = 'Mm-1'
    # case for od
    elif units == 'no unit':
        units = 'unitless'

    combined_ds = xr.Dataset(
        data_vars={
            var: (['station', 'time'], averaged_data,
                  {'units': units})
        },
        coords={
            'station': np.arange(averaged_data.shape[0]),
            'time': standard_time
        }
    )

    # add flags variable
    da_flag = xr.DataArray(
        averaged_flag_data,
        dims=["station", "time", "N_flag_codes"],
        coords={
            "time": standard_time,
        },
        name="flag"
    )
    combined_ds['flag'] = da_flag

    # add qa variable
    da_qa = xr.DataArray(
        averaged_qa_data,
        dims=["station", "time", "N_qa_codes"],
        coords={
            "time": standard_time,
        },
        name="qa"
    )
    combined_ds['qa'] = da_qa

    # add metadata
    for key, value in metadata.items():
        if key in ['latitude', 'longitude']:
            value = [float(val) if val != b'' else np.nan for val in value]
        elif key in ['altitude', 'sampling_height']:
            value = [float(val.decode('utf-8').replace('m', '').strip()) 
                     if val != b'' else np.nan for val in value]
        else:
            value = [val.decode('utf-8') for val in value]
        combined_ds[key] = xr.Variable(data=value, dims=('station'))

    # calculate measurement_altitude if altitude and sampling_height exist
    if ('altitude' in combined_ds.keys()) and ('sampling_height' in combined_ds.keys()):
        value = combined_ds['altitude'].values + combined_ds['sampling_height'].values
        combined_ds['measurement_altitude'] = xr.Variable(data=value, dims=('station'))

    # add units for lat and lon
    # TODO: Check attrs geospatial_lat_units and geospatial_lon_units
    combined_ds.latitude.attrs['units'] = 'degrees_north'
    combined_ds.longitude.attrs['units'] = 'degrees_east'

    # add general attrs
    combined_ds.attrs['data_license'] = f'BSD-3-Clause. Copyright {datetime.datetime.now().year} Alba Vilanova Cortezón'
    combined_ds.attrs['source'] = 'Observations'
    combined_ds.attrs['institution'] = 'Barcelona Supercomputing Center'
    combined_ds.attrs['creator_name'] = 'Alba Vilanova Cortezón'
    combined_ds.attrs['creator_email'] = 'alba.vilanova@bsc.es'
    combined_ds.attrs['application_area'] = 'Monitoring atmospheric composition'
    combined_ds.attrs['domain'] = 'Atmosphere'
    combined_ds.attrs['observed_layer'] = 'Land surface'

    return combined_ds


def get_files_to_download(nonghost_root, target_start_date, target_end_date, resolution, var):
    """Get filenames that should be downloaded

    Parameters
    ----------
    nonghost_root : str
        Directory where non-GHOST data is saved
    target_start_date : datetime.datetime
        Target start date (defined from configuration file)
    target_end_date : datetime.datetime
        Target end date (defined from configuration file)
    resolution : str
        Resolution
    var : str
        Variable

    Returns
    -------
    list
        Filenames that should be downloaded
    """

    base_dir = join(nonghost_root, 'actris/actris', resolution, var)
    paths = []
    current_date = copy.deepcopy(target_start_date)
    while current_date <= target_end_date:

        # save path
        path = f"{base_dir}/{var}_{current_date.strftime('%Y%m')}.nc"
        paths.append(path)

        # get following month
        next_month = current_date.month % 12 + 1
        next_year = current_date.year + (current_date.month // 12)
        current_date = current_date.replace(year=next_year, month=next_month)

    return paths


def download_actris_network(instance, resolution):

    target_start_date = datetime.datetime(int(instance.start_date[:4]), 
                                          int(instance.start_date[4:6]), 
                                          int(instance.start_date[6:8]), 0)
    target_end_date = datetime.datetime(int(instance.end_date[:4]), 
                                        int(instance.end_date[4:6]), 
                                        int(instance.end_date[6:8]), 23, 59, 59) - datetime.timedelta(days=1)

    for var in instance.species:
        
        # check if variable name is available
        if var not in ghost_actris_variables.keys():
            instance.logger.info(f'Data for {var} cannot be downloaded')
            continue
        else:
            actris_parameter = ghost_actris_variables[var]

        # get files that were already downloaded
        initial_check_nc_files = get_files_to_download(instance.nonghost_root, target_start_date, target_end_date, 
                                                       resolution, var)
        files_to_download = instance.select_files_to_download(initial_check_nc_files)
        if not files_to_download:
            msg = f"\nFiles were already downloaded for {var} at {resolution} "
            msg += f"resolution between {target_start_date} and {target_end_date}."
            show_message(instance, msg, deactivate=False)     
            continue 
        
        # get files info path
        path = get_files_path(var)

        # if file does not exist
        if not os.path.isfile(path):
            # get files information
            instance.logger.info(f'\nFile containing information of the files available in Thredds for {var} ({path}) does not exist, creating.')
            combined_data = get_files_per_var(instance, var)
            all_files = combined_data[var]['files']
            files_info = get_files_info(instance, all_files, var, path)
                
        # if file exists
        else:
            # ask if user wants to update file information from NILU Thredds
            if instance.origin_update_choice not in ['y','n']:
                while instance.origin_update_choice not in ['y','n']:
                    instance.origin_update_choice = input(f"\nFile containing information of the files available in Thredds for {var} ({path}) already exists. Do you want to update it (y/n)? ").lower() 
                # ask if user wants to remember the decision
                remind_txt = None
                while remind_txt not in ['y','n']:
                    remind_txt = input("\nDo you want to remember your decision for future downloads (y/n)? ").lower() 
                # save the decision
                if remind_txt == 'y':
                    with open(join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                        f.write(f"ORIGIN_UPDATE={instance.origin_update_choice}\n")
            if instance.origin_update_choice == 'n':
                # get files information
                files_info = yaml.safe_load(open(join(CURRENT_PATH, path)))
                files_info = {k: v for k, v in files_info.items() if k.strip() and v}
            else:
                # get files information
                combined_data = get_files_per_var(instance, var)
                all_files = combined_data[var]['files']
                files_info = get_files_info(instance, all_files, var, path)
        
        # go to next variable if no data is found
        if files_info is not None:
            if len(files_info) == 0:
                continue
        else:
            continue

        # get wavelength
        wavelength_var = is_wavelength_var(actris_parameter)
        if wavelength_var:
            # select most common wavelength for black carbon (name does not provide it)
            if var == 'sconcbc':
                wavelength = 880
                instance.logger.info(f'Wavelength appears in dimensions. Selected wavelength: {wavelength}.')
            # get wavelength from variable name for other variables
            else:
                wavelength = float(re.findall(r'\d+', var)[0])
        else:
            wavelength = None

        # filter files by resolution and dates
        instance.logger.info('Filtering files by resolution and dates...')
        files = {}
        for file, attributes in files_info.items():
            if attributes["resolution"] == resolution:
                start_date = datetime.datetime.strptime(attributes["time_coverage_start"], "%Y-%m-%dT%H:%M:%S UTC")
                end_date = datetime.datetime.strptime(attributes["time_coverage_end"], "%Y-%m-%dT%H:%M:%S UTC")
                for file_to_download in files_to_download:
                    file_to_download_yearmonth = file_to_download.split(f'{var}_')[1].split('.nc')[0]
                    file_to_download_start_date = datetime.datetime.strptime(file_to_download_yearmonth, "%Y%m")
                    file_to_download_end_date = datetime.datetime(file_to_download_start_date.year, 
                                                                  file_to_download_start_date.month, 1) + relativedelta(months=1, seconds=-1)
                    if file_to_download_start_date <= end_date and file_to_download_end_date >= start_date:
                        if 'wavelengths' in attributes:
                            if wavelength is None:
                                instance.logger.error(f'Dataset has wavelength in its dimensions but wavelength is None. Revise if ACTRIS parameter ({actris_parameter}) is included in is_wavelength_var function.')
                                break
                            if wavelength not in attributes['wavelengths']:
                                continue
                        # from filtered files, save those that are provided multiple times
                        station = attributes["ebas_station_code"]
                        if station not in files:
                            files[station] = []
                        if file not in files[station]:
                            files[station].append(file)

        if len(files) != 0:

            # get data for each file within period and temporally average to standard times
            # start = time.time()
            combined_ds = get_data(instance, files, var, actris_parameter, resolution, 
                                    target_start_date, target_end_date, files_info,
                                    instance.ghost_version, instance.n_cpus)
            if combined_ds is None:
                continue
            # end = time.time()
            # elapsed_minutes = (end - start) / 60
            # instance.logger.info(f"Time to read data: {elapsed_minutes:.2f} minutes")

            # save data per year and month
            path = join(instance.nonghost_root, f'actris/actris/{resolution}/{var}')
            if not os.path.isdir(path):
                os.makedirs(path, exist_ok=True)
            saved_files = 0
            for year, ds_year in combined_ds.groupby('time.year'):
                for month, ds_month in ds_year.groupby('time.month'):
                    filename = f"{path}/{var}_{year}{month:02d}.nc"
                    if filename in files_to_download:
                        combined_ds_yearmonth = combined_ds.sel(time=f"{year}-{month:02d}")

                        # add title to attrs
                        extra_info = ''
                        if wavelength_var and wavelength is not None:
                            extra_info = f' at {wavelength}nm'
                        combined_ds_yearmonth.attrs['title'] = f'Surface {ghost_actris_variables[var]}{extra_info} in the ACTRIS network in {year}-{month:02d}.'

                        # order attrs
                        custom_order = ['title', 'institution', 'creator_name', 'creator_email',
                                        'source', 'application_area', 'domain', 'observed_layer',
                                        'data_license']
                        ordered_attrs = {key: combined_ds_yearmonth.attrs[key] 
                                        for key in custom_order 
                                        if key in combined_ds_yearmonth.attrs}
                        combined_ds_yearmonth.attrs = ordered_attrs

                        # remove stations if all variable data is nan
                        # previous_n_stations = len(combined_ds_yearmonth.station)
                        combined_ds_yearmonth = combined_ds_yearmonth.dropna(dim="station", subset=[var], how="all")
                        combined_ds_yearmonth = combined_ds_yearmonth.assign_coords(station=range(len(combined_ds_yearmonth.station)))
                        # current_n_stations = len(combined_ds_yearmonth.station)
                        # n_stations_diff = previous_n_stations - current_n_stations
                        # if n_stations_diff > 0:
                        #     instance.logger.info(f'Data for {n_stations_diff} stations was removed because all data was NaN during {month}-{year}.')

                        # remove file if it exists
                        if os.path.isfile(filename):
                            os.system("rm {}".format(filename))

                        # do not save if empty
                        if len(combined_ds_yearmonth[var].values) == 0:
                            continue
                            
                        # save file
                        combined_ds_yearmonth.to_netcdf(filename)

                        # change permissions
                        os.system("chmod 777 {}".format(filename))
                        instance.logger.info(f"Saved: {filename}")
                        saved_files += 1
                        
            instance.logger.info(f'Total number of saved files: {saved_files}')

        else:
            instance.logger.info(f'No files were found at {resolution} resolution for {var}')
