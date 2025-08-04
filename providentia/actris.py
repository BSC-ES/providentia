import bisect
import copy
import csv
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

from providentia.auxiliar import CURRENT_PATH, join, pad_array

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

# load ACTRIS mapping files
parameters_dict = yaml.safe_load(open(join(
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
        url = f"{base_url}/{parameters_dict[var]}/page/{page}"
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
        sys.exit(1)

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
            is_max = attrs_dict['ebas_data_level'] == attrs_dict['ebas_data_level'][np.argmax(np.float32(attrs_dict['ebas_data_level']))]
            
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


def get_data(download_instance, files, var, actris_parameter, resolution, target_start_date, target_end_date, 
             files_info, ghost_version):
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

    # initialise metadata
    metadata = {}

    # get EBAS component
    ebas_component = variable_mapping[actris_parameter]['var']

    # initialise wavelength
    wavelength = None
    
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

    # initialise averaged data
    averaged_data = np.full(
        (len(files.items()), len(standard_time)), fill_value=np.nan, dtype=np.float32)
    averaged_flag_data = np.full(
        (len(files.items()), len(standard_time), 186), fill_value=np.nan, dtype=np.float32)
    averaged_qa_data = np.full(
        (len(files.items()), len(standard_time), 2), fill_value=np.nan, dtype=np.float32)

    # define vectorize function to get GHOST_decreed_validity by overlap flags later
    vfunc = np.vectorize(
        lambda flag: np.nan if np.isnan(flag) 
                    else flags_dict['000']["GHOST_decreed_validity"][0] if flag == 0 
                    else flags_dict[str(int(flag))]["GHOST_decreed_validity"][0],
        otypes=['O']
    )

    # create pairs of valid dates
    standard_time_pairs = create_time_pairs(standard_time)

    errors = {}
    warnings = {}
    tqdm_iter = tqdm(
        files.items(), bar_format='{l_bar}{bar}|{n_fmt}/{total_fmt}', desc=f"Reading data")
    for i, (station, urls) in enumerate(tqdm_iter):

        # open files
        try:
            if len(urls) > 1:
                url = select_station_file(urls, files_info)
                if url is None:
                    continue
            else:
                url = urls[0]

            nc = Dataset(url, mode='r')
            ds = xr.open_dataset(xr.backends.NetCDF4DataStore(nc))

            possible_vars, possible_var = get_var_in_file(ds, var, actris_parameter, ebas_component)
            if possible_var is None:
                errors[
                    station] = f'No variable name matches for {possible_vars}. Existing keys: {list(ds.data_vars)}.'
                continue

        except Exception as error:
            errors[station] = f'Error opening file: {error}.'
            continue

        warnings[station] = ""

        # remove time duplicates if any (keep first)
        ds = ds.sel(time=~ds['time'].to_index().duplicated())
        
        # select data in period range
        ds = ds.sel(time=slice(target_start_date, target_end_date))
        if ds.time.size == 0:
            errors[station] = f'No data available after filtering by time.'
            continue

        # rename qc dimension
        ds = ds.rename_dims({f'{possible_var}_qc_flags': 'N_flag_codes'})

        # get lowest level if tower height is in coordinates
        if 'Tower_inlet_height' in list(ds.coords):
            warnings[station] += f'Taking data from first height (Tower_inlet_height={min(ds.Tower_inlet_height.values)}). '
            ds = ds.sel(Tower_inlet_height=min(
                ds.Tower_inlet_height.values), drop=True)

        # get data at desired wavelength if wavelength is in coordinates
        if 'Wavelength' in list(ds.coords) or 'Wavelengthx' in list(ds.coords):
            # Select most common wavelength for black carbon (name does not provide it)
            if var == 'sconcbc':
                wavelength = 880
                warnings[
                    station] += f'Wavelength appears in dimensions. Selected wavelength: {wavelength}. '
            # Get wavelength from variable name for other variables
            else:
                wavelength = float(re.findall(r'\d+', var)[0])

            # Select data for wavelength
            found_wavelength = False
            if 'Wavelengthx' in list(ds.coords):
                if wavelength in ds.Wavelengthx.values:
                    ds = ds.sel(Wavelengthx=wavelength, drop=True)
                    found_wavelength = True
                else:
                    existing_wavelengths = ds.Wavelengthx.values
            elif 'Wavelength' in list(ds.coords):
                if wavelength in ds.Wavelength.values:
                    ds = ds.sel(Wavelength=wavelength, drop=True)
                    found_wavelength = True
                else:
                    existing_wavelengths = ds.Wavelength.values

            if not found_wavelength:
                warnings[station] += f'Data at {wavelength}nm could not be found. Existing wavelengths: {existing_wavelengths}. '
                continue

        # remove artifact and fraction (sconcoc)
        # TODO: Discuss this
        if 'Artifact' in list(ds.coords):
            warnings[station] += f'Taking data from first artifact dimension (Artifact={ds.Artifact.values[0]}). '
            ds = ds.isel(Artifact=0, drop=True)
        if 'Fraction' in list(ds.coords):
            warnings[station] += f'Taking data from first fraction dimension (Fraction={ds.Fraction.values[0]}). '
            ds = ds.isel(Fraction=0, drop=True)

        # read variable
        da_var = ds[possible_var]

        # avoid datasets that do not have defined units
        if 'ebas_unit' not in da_var.attrs:
            errors[station] = f'No units were defined.'
            continue

        # avoid datasets that do not have the same units as in variable mapping
        if da_var.attrs['ebas_unit'] != variable_mapping[actris_parameter]['units']:
            errors[station] = f"Units {da_var.attrs['ebas_unit']} do not match those in variable mapping "
            errors[station] += f"dictionary ({variable_mapping[actris_parameter]['units']})."
            continue

        # remove all attributes except units
        da_var_attrs = copy.deepcopy(da_var.attrs)
        da_var.attrs = {key: value for key,
                        value in da_var.attrs.items() if key in ['ebas_unit', 'ebas_station_code']}

        # read quality control data
        flag_data = ds[f'{possible_var}_qc'].transpose(
            "time", "N_flag_codes")

        # rename variable to BSC standards
        station_ds = da_var.to_dataset(name=var)

        # add time bounds
        station_ds['time_bnds'] = ds['time_bnds']

        # add quality control data
        station_ds['flag'] = flag_data

        # temporally average data from original times to standard times
        station_averaged_data, station_flag_data, station_qa_data = temporally_average_data(station_ds, var, ghost_version, standard_time_pairs, vfunc)
        if np.isnan(station_averaged_data).all():
            warnings[station] += 'No data after temporal averaging. '
            continue

        # save metadata
        for ghost_key, ebas_key in metadata_dict.items():
            # create key if it does not exist
            if ghost_key not in metadata.keys():
                metadata[ghost_key] = []

            # search value in var attrs
            if ebas_key in da_var_attrs.keys():
                metadata[ghost_key].append(da_var_attrs[ebas_key])
            # search value in ds attrs
            elif ebas_key in ds.attrs.keys():
                metadata[ghost_key].append(ds.attrs[ebas_key])
            # not found -> nan
            else:
                metadata[ghost_key].append(np.nan)

        # save
        averaged_data[i, :] = station_averaged_data
        averaged_flag_data[i, :, :] = station_flag_data
        averaged_qa_data[i, :, :] = station_qa_data

    # drop stations that have nan for all times
    mask = ~np.isnan(averaged_data).all(axis=1)
    averaged_data = averaged_data[mask]
    averaged_flag_data = averaged_flag_data[mask, :]
    averaged_qa_data = averaged_qa_data[mask, :]

    # create dataset with averaged data
    units = variable_mapping[actris_parameter]['units']
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

    # show errors
    if len(errors) > 0:
        download_instance.logger.info(f'Collected errors ({len(errors)}):')
        for file, error in errors.items():
            download_instance.logger.info(f'{file} - Error: {error}')

    # show warnings
    if len(warnings) > 0 and all(len(warning) > 0 for warning in warnings.values()):
        download_instance.logger.info(f'Collected warnings ({len(warnings)}):')
        for file, warning in warnings.items():
            if len(warning) > 0:
                download_instance.logger.info(f'{file} - Warning: {warning}')

    return combined_ds, metadata, wavelength


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
