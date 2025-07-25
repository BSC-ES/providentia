import copy
import csv
import datetime
import itertools
import os
import requests
import sys
import yaml
import re

from tqdm import tqdm
import numpy as np
import pandas as pd
import xarray as xr

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


def get_files_per_var(var):
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
            print(
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
        flag_info = flags_dict[str(int(flag))]
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


def check_overlap(station_ds, var, standard_start_date, standard_end_date, measurement_time_pairs):

    overlap_durations = []
    overlap_var_data = []
    overlap_flag_data = []

    for i, (measurement_start_date, measurement_end_date) in enumerate(measurement_time_pairs):
        overlap_start = max(measurement_start_date, standard_start_date)
        overlap_end = min(measurement_end_date, standard_end_date)

        # Detect and save overlap
        if overlap_start < overlap_end:
            overlap_duration = np.array(pd.Timedelta(
                overlap_end - overlap_start).total_seconds() / 60)
            data = station_ds.sel(time=slice(
                measurement_start_date, measurement_end_date))

            # get values at 0 to remove time dimension (always 1)
            var_data = data[var].values[0]
            flag_data = data['flag'].values[0]

            overlap_durations.append(overlap_duration)
            overlap_var_data.append(var_data)
            overlap_flag_data.append(flag_data)

        # If there is no overlap and there have already been overlaps, go to next period
        elif len(overlap_durations) > 0:
            break

    overlap_durations = np.array(overlap_durations)
    overlap_var_data = np.array(overlap_var_data)
    overlap_flag_data = np.array(overlap_flag_data)

    return overlap_durations, overlap_var_data, overlap_flag_data


def temporally_average_data(combined_ds_list, resolution, var, ghost_version, target_start_date, target_end_date):
    """Temporally average data and get unique flags in the valid times (temporally averaged)

    Parameters
    ----------
    combined_ds_list : xarray.Dataset
        Concatenated data for all stations with original times per month
    resolution : str
        Temporal resolution
    var : str
        Variable
    ghost_version : float
        GHOST version

    Returns
    -------
    combined_ds : xarray.Dataset
        Concatenated data for all stations with valid times per month
    """

    # get valid dates frequency
    if resolution == 'hourly':
        frequency = 'h'
        timedelta = np.timedelta64(30, 'm')
    elif resolution == '3hourly':
        frequency = '3h'
        timedelta = np.timedelta64(90, 'm')
    elif resolution == 'daily':
        frequency = 'D'
        timedelta = np.timedelta64(12, 'h')
    # TODO: Review this
    elif resolution == 'monthly':
        frequency = 'MS'
        timedelta = np.timedelta64(15, 'D')

    standard_time = pd.date_range(start=target_start_date, end=target_end_date,
                                  freq=frequency).to_numpy(dtype='datetime64[ns]')

    # initialise averaged data
    averaged_data = np.empty(
        (len(combined_ds_list), len(standard_time)), dtype=np.float32)
    averaged_flag_data = np.empty(
        (len(combined_ds_list), len(standard_time), 186), dtype=np.float32)
    averaged_qa_data = np.empty(
        (len(combined_ds_list), len(standard_time), 2), dtype=np.float32)

    # define vectorize function to get GHOST_decreed_validity by overlap flags later
    vfunc = np.vectorize(
        lambda flag: np.nan if np.isnan(flag) else flags_dict[str(
            int(flag))]["GHOST_decreed_validity"][0],
        otypes=['O']
    )

    # create pairs of valid dates
    standard_time_pairs = create_time_pairs(standard_time)

    tqdm_iter = tqdm(combined_ds_list, bar_format='{l_bar}{bar}|{n_fmt}/{total_fmt}',
                     desc=f"Temporally averaging data")
    for station_i, station_ds in enumerate(tqdm_iter):

        # initialise averaged data
        station_averaged_data = []
        station_flag_data = []
        station_qa_data = []

        # remove station selection
        station_ds = station_ds.isel(station=0)

        measurement_time_pairs = [(t - timedelta, t + timedelta)
                                  for t in station_ds.time.values]
        for i, (standard_start_date, standard_end_date) in enumerate(standard_time_pairs):
            overlap_durations, overlap_var_data, overlap_flag_data = check_overlap(station_ds, var,
                                                                                   standard_start_date,
                                                                                   standard_end_date,
                                                                                   measurement_time_pairs)

            if len(overlap_var_data) > 0:
                if len(overlap_var_data) > 1:
                    GHOST_decreed_validities = vfunc(overlap_flag_data)
                    GHOST_invalid = np.any(
                        GHOST_decreed_validities == 'I', axis=1)
                    # if there are invalid values in period but some of them are valid, convert invalid ones to nan
                    if np.any(GHOST_invalid) and not np.all(GHOST_invalid):
                        for i in range(len(overlap_var_data)):
                            if GHOST_invalid[i]:
                                overlap_var_data[i] = np.array(np.nan)
                                overlap_flag_data[i, :] = np.array(
                                    [np.nan]*overlap_flag_data.shape[1])

                # remove nan
                valid_mask = ~np.isnan(overlap_var_data)
                valid_durations = overlap_durations[valid_mask]
                valid_values = overlap_var_data[valid_mask]

                # get weighted average if there is more than one value
                if len(valid_values) > 1:
                    mean = np.average(valid_values, weights=valid_durations)
                # get value if there is only one value
                elif len(valid_values) == 1:
                    mean = valid_values[0]
                # otherwise nan
                else:
                    mean = np.nan
                station_averaged_data.append(mean)

                # get unique flag values and convert to standard flag names (instead of EBAS) and qa
                flags = np.unique(overlap_flag_data)
                valid_flags = flags[~np.isnan(flags)]
                standard_flag_codes, qa = get_standard_flags_and_qa(
                    valid_flags, ghost_version)

                # save flags and qa
                station_flag_data.append(
                    pad_array(standard_flag_codes, length=186))
                station_qa_data.append(pad_array(qa, length=2))

            # if file has no data for dates, set it to be nan
            else:
                station_averaged_data.append(np.nan)
                station_flag_data.append([np.nan]*186)
                station_qa_data.append([np.nan]*2)

        averaged_data[station_i, :] = station_averaged_data
        averaged_flag_data[station_i, :, :] = station_flag_data
        averaged_qa_data[station_i, :, :] = station_qa_data

    # create dataset with averaged data
    units = combined_ds_list[0][var].attrs['ebas_unit']
    combined_ds = xr.Dataset(
        data_vars={
            var: (['station', 'time'], averaged_data,
                  {'units': units})
        },
        coords={
            'station': np.arange(len(combined_ds_list)),
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

    return combined_ds


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


def get_files_info(files, var, path):
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
            ds = xr.open_dataset(file)
        except:
            continue

        # get statistics
        if 'ebas_statistics' in ds.attrs:
            file_statistics = ds.attrs['ebas_statistics']
        else:
            file_statistics = 'Unknown'

        # get resolution
        coverage = ds.time_coverage_resolution
        try:
            file_resolution = coverages_dict[coverage]
        except:
            file_resolution = f'Unrecognised ({coverage})'

        file_start_date = ds.time_coverage_start
        file_end_date = ds.time_coverage_end
        file_variables = list(ds.data_vars.keys())
        file_reference = ds.attrs['ebas_station_code']

        # save in dict
        files_info[file] = {}
        files_info[file]['resolution'] = file_resolution
        files_info[file]['start_date'] = file_start_date
        files_info[file]['end_date'] = file_end_date
        files_info[file]['variables'] = file_variables
        files_info[file]['statistics'] = file_statistics
        files_info[file]['station_reference'] = file_reference

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
        print(f'Error: No data could be found for {var}')

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


def get_data(files, var, actris_parameter, resolution, target_start_date, target_end_date):
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

    # combine datasets that have the same variable and resolution
    combined_ds_list = []
    metadata = {}
    metadata[resolution] = {}

    # get EBAS component
    ebas_component = variable_mapping[actris_parameter]['var']

    # initialise wavelength
    wavelength = None
    
    errors = {}
    warnings = {}
    tqdm_iter = tqdm(
        files.items(), bar_format='{l_bar}{bar}|{n_fmt}/{total_fmt}', desc=f"Reading data")
    for i, (station, urls) in enumerate(tqdm_iter):

        # open file
        try:
            if len(urls) == 1:
                ds = xr.open_dataset(urls[0])
                possible_vars, possible_var = get_var_in_file(ds, var, actris_parameter, ebas_component)
                if possible_var is None:
                    errors[
                        station] = f'No variable name matches for {possible_vars}. Existing keys: {list(ds.data_vars)}.'
                    continue
            else:
                url_vars = []
                for url in urls:
                    ds = xr.open_dataset(url)
                    possible_vars, possible_var = get_var_in_file(ds, var, actris_parameter, ebas_component)
                    if possible_var is None:
                        errors[
                            station] = f'No variable name matches for {possible_vars}. Existing keys: {list(ds.data_vars)}.'
                        continue
                    url_vars.append(possible_var)

                # throw error if the variable names across datasets are different for main species
                if len(list(set(url_vars))) > 1:
                    msg = 'There is more than one dataset for the same station in the Thredds and the variable names are not the same. '
                    msg += f'Unique variable names for all files: {list(set(url_vars))}'
                    errors[station] = msg
                    continue
                ds = xr.open_mfdataset(urls, combine='nested', concat_dim='time')
                ds = ds.sortby('time')
        except Exception as error:
            errors[station] = f'{i} - Error opening file: {error}.'
            continue
        
        warnings[station] = ""

        # remove time duplicates if any (keep first)
        ds = ds.sel(time=~ds['time'].to_index().duplicated())

        # select data in period range
        ds = ds.sel(time=slice(target_start_date, target_end_date))

        # assign station code as dimension
        ds = ds.expand_dims(dim={'station': [i]})

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

        # save metadata
        for ghost_key, ebas_key in metadata_dict.items():
            # create key if it does not exist
            if ghost_key not in metadata[resolution].keys():
                metadata[resolution][ghost_key] = []

            # search value in var attrs
            if ebas_key in da_var.attrs.keys():
                metadata[resolution][ghost_key].append(da_var.attrs[ebas_key])
            # search value in ds attrs
            elif ebas_key in ds.attrs.keys():
                metadata[resolution][ghost_key].append(ds.attrs[ebas_key])
            # not found -> nan
            else:
                metadata[resolution][ghost_key].append(np.nan)

        # remove all attributes except units
        da_var.attrs = {key: value for key,
                        value in da_var.attrs.items() if key in ['ebas_unit', 'ebas_station_code']}

        # read quality control data
        flag_data = ds[f'{possible_var}_qc'].transpose(
            "station", "time", "N_flag_codes")

        # rename variable to BSC standards
        ds_station = da_var.to_dataset(name=var)

        # add quality control data
        ds_station['flag'] = flag_data

        # append modified dataset to list
        combined_ds_list.append(ds_station)

    # show errors
    if len(errors) > 0:
        print(f'Collected errors ({len(errors)}):')
        for file, error in errors.items():
            print(f'{file} - Error: {error}')

    # show warnings
    if len(warnings) > 0 and all(len(warning) > 0 for warning in warnings.values()):
        print(f'Collected warnings ({len(warnings)}):')
        for file, warning in warnings.items():
            if len(warning) > 0:
                print(f'{file} - Warning: {warning}')

    return combined_ds_list, metadata, wavelength


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
