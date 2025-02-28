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

from providentia.auxiliar import CURRENT_PATH, join

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

# load ACTRIS mapping files
parameters_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'ghost_actris_variables.yaml')))
metadata_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'metadata.yaml')))
coverages_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'coverages.yaml')))
variable_mapping = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'variable_mapping.yaml')))
variable_mapping = {k: v for k, v in variable_mapping.items() if k.strip() and v}


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
    
    sys.path.insert(1, join(PROVIDENTIA_ROOT, 'providentia/dependencies/GHOST_standards/{}'.format(ghost_version)))
    from GHOST_standards import standard_parameters
    with open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'ghost_variables.csv'), 
              mode="w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        for key in standard_parameters.keys():
            writer.writerow([standard_parameters[key]['long_parameter_name'], 
                             standard_parameters[key]['bsc_parameter_name'], 
                             ', '.join( standard_parameters[key]['ebas_parameter_name'])])


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
        path = join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', f'files/{alpha_var}/files.yaml')
    else:
        path = join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', f'files/{var}/files.yaml')

    return path


def pad_array(arr):
    """Pad array with nan if its length is less than 186."""

    pad_size = max(0, 186 - len(arr))

    return np.pad(arr, (0, pad_size), constant_values=np.nan)


def temporally_average_data(combined_ds, resolution, year, month, var):
    """Temporally average data and get unique flags in the valid times (temporally averaged)

    Parameters
    ----------
    combined_ds : xarray.Dataset
        Concatenated data for all stations with original times per month
    resolution : str
        Temporal resolution
    year : int
        Year
    month : int
        Month
    var : str
        Variable

    Returns
    -------
    combined_ds : xarray.Dataset
        Concatenated data for all stations with valid times per month
    """

    # get valid dates frequency
    if resolution == 'hourly':
        frequency = 'h'
    elif resolution == 'daily':
        frequency = 'D'
    elif resolution == 'monthly':
        frequency = 'MS'

    # get start and end of period to construct valid dates
    time = combined_ds.time.values
    start_date = datetime.datetime(year, month, 1)
    first_day_next_month = datetime.datetime(year, month % 12 + 1, 1) if month != 12 else datetime.datetime(year + 1, 1, 1)
    end_date = (first_day_next_month - datetime.timedelta(days=1)).replace(hour=23)
    valid_dates = pd.date_range(start=start_date, end=end_date, freq=frequency).to_numpy(dtype='datetime64[ns]')
    
    # initialise averaged data
    averaged_data = np.empty((len(combined_ds.station.values), len(valid_dates)))
    flag_data = np.empty((len(combined_ds.station.values), len(valid_dates), 186))

    for station_i, station in enumerate(combined_ds.station.values):

        # initialise averaged data
        station_averaged_data = []

        # read data per station
        data = combined_ds[var].isel(station=station_i).values

        # get indices where values are nan
        valid_idxs = ~np.isnan(data)

        # go to next station if all values are nan
        if np.sum(valid_idxs) == 0:
            continue

        # filter nan values
        valid_time = time[valid_idxs]
        valid_data = data[valid_idxs]
        
        # calculate weighted averages
        if len(valid_data) != 0:
            for date in valid_dates:
            
                # get differences between valid time and actual times in nanoseconds
                time_diffs = (valid_time - date).astype('timedelta64[ns]').astype(float)
            
                # get positive differences and negative differences to differentiate 
                # between the actual times that are earlier than the valid date (negative), and those that are later (positive)
                positive_diffs = time_diffs[time_diffs > 0]
                negative_diffs = time_diffs[time_diffs < 0]
                
                # find the closest actual time after the valid time
                closest_positive = None
                if len(positive_diffs) > 0:
                    closest_positive_idx = np.abs(positive_diffs).argmin()
                    closest_positive = positive_diffs[closest_positive_idx]
                    closest_positive_time = valid_time[time_diffs == positive_diffs[closest_positive_idx]][0]
                    closest_positive_value = valid_data[time_diffs == positive_diffs[closest_positive_idx]][0]
            
                # find the closest actual time before the valid time
                closest_negative = None
                if len(negative_diffs) > 0:
                    closest_negative_idx = np.abs(negative_diffs).argmin()
                    closest_negative = negative_diffs[closest_negative_idx]
                    closest_negative_time = valid_time[time_diffs == negative_diffs[closest_negative_idx]][-1]
                    closest_negative_value = valid_data[time_diffs == negative_diffs[closest_negative_idx]][-1]
            
                # when the valid time only has a value in one direction, get closest value without calculating weights
                if closest_positive is None:
                    value = closest_negative_value
                elif closest_negative is None:
                    value = closest_positive_value
                # in the rest of cases, calculate weights of 2 closest values and make average
                else:
                    # get 2 closest times and make positive to be able to compare differences
                    closest_diffs = np.abs([closest_negative, closest_positive])
            
                    # we do the reverse, since we want the differences to have a heavier weight if these are smaller (nearer the actual time)
                    weights = 1 / closest_diffs
                    
                    # finally we normalize them to have values between 0 and 1
                    weights_normalized = weights / np.sum(weights)
            
                    # get average
                    value = np.average([closest_negative_value, closest_positive_value], weights=weights_normalized)
        
                # save averaged data
                station_averaged_data.append(value)
        
            averaged_data[station_i, :] = station_averaged_data
        else:
            averaged_data[station_i, :] = [np.nan]*len(valid_dates)

        # create pairs of valid dates
        time_pairs = list(zip(valid_dates[:-1], valid_dates[1:]))
        last_time_pair = (valid_dates[-1], valid_dates[-1] + (valid_dates[-1] - valid_dates[-2]))
        time_pairs.append(last_time_pair)
        
        # get data between each valid start date and end date, then get unique values for the available timesteps, and pad
        # to have arrays of the same length (maximum length is 186, when there is less these are going to be nan)
        station_flag_data = []
        for start_date, end_date in time_pairs:
            unique_flag_values_per_pair = np.unique(combined_ds.flag.sel(time=slice(start_date, end_date), station=station).values)
            unique_flag_values_per_pair_nonan = unique_flag_values_per_pair[~np.isnan(unique_flag_values_per_pair)]
            station_flag_data.append(pad_array(unique_flag_values_per_pair_nonan))

        # get station flag data for the valid dates
        flag_data[station_i, :, :] = station_flag_data

    # create new variable with averaged data
    combined_averaged_da = xr.DataArray(
        data=averaged_data,
        coords={'station': combined_ds.station.values, 'time': valid_dates}, 
        dims=['station', 'time'],
        attrs={'units': combined_ds[var].attrs['ebas_unit']})
    
    # drop old variable and associated time
    combined_ds = combined_ds.drop_vars(var)
    combined_ds = combined_ds.drop_dims('time')
    
    # add new variable
    combined_ds[var] = combined_averaged_da

    # add qa variable
    da_flag = xr.DataArray(
            flag_data,
            dims=["station", "time", "N_flag_codes"],
            coords={
                "time": combined_ds.coords['time'],
            },
            name="flag"
        )
    combined_ds['flag'] = da_flag

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
    tqdm_iter = tqdm(files,bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Creating information file ({len(files)})")
    for file in tqdm_iter:
        # open file
        try:
            ds = xr.open_dataset(file)
        except:
            continue

        # get resolution
        coverage = ds.time_coverage_resolution
        try:             
            file_resolution = coverages_dict[coverage]
        except:
            file_resolution = f'Unrecognised ({coverage})'
            
        file_start_date = ds.time_coverage_start
        file_end_date = ds.time_coverage_end
        file_variables = list(ds.data_vars.keys())
        files_info[file] = {}
        files_info[file]['resolution'] = file_resolution
        files_info[file]['start_date'] = file_start_date
        files_info[file]['end_date'] = file_end_date
        files_info[file]['variables'] = file_variables

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
        print(f'    Error: No data could be found for {var}')
    
    return files_info


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
    tqdm_iter = tqdm(files,bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Reading data ({len(files)})")
    for i, file in enumerate(tqdm_iter):
        
        # open file
        try:
            ds = xr.open_dataset(file)
        except:
            errors[file] = 'Error opening file'
            continue

        # remove time duplicates if any (keep first)
        ds = ds.sel(time=~ds['time'].to_index().duplicated())
        
        # select data in period range
        ds = ds.sel(time=slice(target_start_date, target_end_date))
        
        # assign station code as dimension
        ds = ds.expand_dims(dim={'station': [i]})

        # select data for that variable only
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
            errors[file] = f'No variable name matches for {possible_vars}. Existing keys: {list(ds.data_vars)}'
            continue
            
        # rename qc dimension
        ds = ds.rename_dims({f'{possible_var}_qc_flags': 'N_flag_codes'})

        # get lowest level if tower height is in coordinates
        if 'Tower_inlet_height' in list(ds.coords):
            warnings[file] = f'Taking data from first height (Tower_inlet_height={min(ds.Tower_inlet_height.values)})'
            ds = ds.sel(Tower_inlet_height=min(ds.Tower_inlet_height.values), drop=True)

        # get data at desired wavelength if wavelength is in coordinates
        if 'Wavelength' in list(ds.coords) or 'Wavelengthx' in list(ds.coords):
            # Select most common wavelength for black carbon (name does not provide it)
            if var == 'sconcbc':
                wavelength = 880
                warnings[file] = f'Wavelength appears in dimensions. Selected wavelength: {wavelength}'
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
                warnings[file] = f'Data at {wavelength}nm could not be found. Existing wavelengths: {existing_wavelengths}'
                continue             
                
        # remove artifact and fraction (sconcoc)
        # TODO: Discuss this
        if 'Artifact' in list(ds.coords):
            warnings[file] = f'Taking data from first artifact dimension (Artifact={ds.Artifact.values[0]})'
            ds = ds.isel(Artifact=0, drop=True)
        if 'Fraction' in list(ds.coords):
            warnings[file] = f'Taking data from first fraction dimension (Fraction={ds.Fraction.values[0]})'
            ds = ds.isel(Fraction=0, drop=True)

        # read variable
        da_var = ds[possible_var]
        
        # avoid datasets that do not have defined units
        if 'ebas_unit' not in da_var.attrs:
            errors[file] = f'No units were defined'
            continue

        # avoid datasets that do not have the same units as in variable mapping
        if da_var.attrs['ebas_unit'] != variable_mapping[actris_parameter]['units']:
            errors[file] = f"Units {da_var.attrs['ebas_unit']} do not match those in variable mapping "
            errors[file] += f"dictionary ({variable_mapping[actris_parameter]['units']})"
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
        da_var.attrs = {key: value for key, value in da_var.attrs.items() if key == 'ebas_unit'}

        # read quality control data
        flag_data = ds[f'{possible_var}_qc'].transpose("station", "time", "N_flag_codes")

        # rename variable to BSC standards
        ds_station = da_var.to_dataset(name=var)

        # add quality control data
        ds_station['flag'] = flag_data
              
        # append modified dataset to list
        combined_ds_list.append(ds_station)
    
    # show errors
    if len(errors) > 0:
        print(f'\nCollected errors ({len(errors)}):')
        for file, error in errors.items():
            print(f'{file} - Error: {error}')
            
    # show warnings
    if len(warnings) > 0:
        print(f'\nCollected warnings ({len(warnings)}):')
        for file, warning in warnings.items():
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