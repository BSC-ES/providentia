import requests
import itertools
import copy
import os
import yaml
import xarray as xr
import numpy as np
from variable_mapping import variable_mapping
import datetime
import csv
import sys
import time
import pandas as pd
import re
from dotenv import dotenv_values
from dateutil.relativedelta import relativedelta
from tqdm import tqdm

from providentia.auxiliar import CURRENT_PATH, join

PROVIDENTIA_ROOT = '/home/avilanov/software/Providentia/'
nonghost_root = '/home/avilanov/data/providentia/obs/nonghost/'
CURRENT_PATH = os.getcwd()
sys.path = [path for path in sys.path if '/home/avilanov/software/Providentia/providentia/dependencies/GHOST_standards/' not in path]            
sys.path.insert(1, os.path.join(CURRENT_PATH, '/home/avilanov/software/Providentia/providentia/dependencies/GHOST_standards/1.5'))
from GHOST_standards import standard_parameters, get_standard_metadata

coverages_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'coverages.yaml')))
parameters_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'ghost_actris_variables.yaml')))
variable_mapping = yaml.safe_load(open(join(CURRENT_PATH, 'settings', 'internal', 'actris', 'variable_mapping.yaml')))
variable_mapping = {k: v for k, v in variable_mapping.items() if k.strip() and v}
metadata_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'metadata.yaml')))

variables = ['sconco3']
resolution = 'hourly'
target_start_date = datetime.datetime(2018, 1, 1, 0)
target_end_date = datetime.datetime(2019, 1, 1 ,0)

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


def select_files_to_download(prov_start_time, nc_files_to_download):
    """ Returns the files that are not already downloaded. """
    # initialise list of non-downloaded files
    not_downloaded_files = []
    
    # get ssh user and password 
    env = dotenv_values(join(PROVIDENTIA_ROOT, ".env"))
    overwrite_choice = env.get("OVERWRITE")

    if nc_files_to_download:
        # get the downloaded and not downloaded files
        not_downloaded_files = list(filter(lambda x:not os.path.exists(x), nc_files_to_download))
        downloaded_files = list(filter(lambda x:os.path.exists(x), nc_files_to_download))
        
        # get the files that were downloaded before the execution
        downloaded_before_execution_files = list(filter(lambda x:prov_start_time > os.path.getctime(x), downloaded_files))

        # if there was any file downloaded before the execution    
        if downloaded_before_execution_files:
            # make the user choose between overwriting or not overwriting
            if overwrite_choice not in ['y','n']:
                # ask if user wants to overwrite
                while overwrite_choice not in ['y','n']:
                    overwrite_choice = input("\nThere are some files that were already downloaded in a previous download, do you want to overwrite them (y/n)? ").lower() 
                # ask if user wants to remember the decision
                remind_txt = None
                while remind_txt not in ['y','n']:
                    remind_txt = input("\nDo you want to remember your decision for future downloads (y/n)? ").lower() 
                # save the decision
                if remind_txt == 'y':
                    with open(join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                        f.write(f"OVERWRITE={overwrite_choice}\n")
            # if user wants to overwrite then add the the files that were 
            # downloaded before the execution as if they were never download
            if overwrite_choice == 'y':
                not_downloaded_files += downloaded_before_execution_files
            # change overwritten files boolean to True to indicate that some files were ignored
            else:
                overwritten_files_flag = True

    return not_downloaded_files


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


def pad_array(arr):
    """Pad array with 255 if its length is less than 186."""
    pad_size = max(0, 186 - len(arr))
    return np.pad(arr, (0, pad_size), constant_values=255)


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
        # to have arrays of the same length (maximum length is 186, when there is less these are going to be 255)
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


# get providentia start time
prov_start_time = time.time()

env = dotenv_values(join(PROVIDENTIA_ROOT, ".env"))
origin_update_choice = env.get("ORIGIN_UPDATE")

for var in variables:

    # check if variable name is available
    if var not in parameters_dict.keys():
        print(f'Data for {var} cannot be downloaded')
        continue
    else:
        actris_parameter = parameters_dict[var]
    
    # get files that were already downloaded
    initial_check_nc_files = get_files_to_download(nonghost_root, target_start_date, target_end_date, resolution, var)
    files_to_download = select_files_to_download(prov_start_time, initial_check_nc_files)
    if not files_to_download:
        msg = f"\nFiles were already downloaded for {var} at {resolution} "
        msg += f"resolution between {target_start_date} and {target_end_date}."
        print(msg)  
        continue 

    # get files info path
    path = get_files_path(var)
    
    # if file does not exist
    if not os.path.isfile(path):
        # get files information
        print(f'\nFile containing information of the files available in Thredds for {var} ({path}) does not exist, creating.')
        combined_data = get_files_per_var(var)
        all_files = combined_data[var]['files']
        files_info = get_files_info(all_files, var, path)
            
    # if file exists
    else:
        # ask if user wants to update file information from NILU Thredds
        if origin_update_choice not in ['y','n']:
            while origin_update_choice not in ['y','n']:
                origin_update_choice = input(f"\nFile containing information of the files available in Thredds for {var} ({path}) already exists. Do you want to update it (y/n)? ").lower() 
            # ask if user wants to remember the decision
            remind_txt = None
            while remind_txt not in ['y','n']:
                remind_txt = input("\nDo you want to remember your decision for future downloads (y/n)? ").lower() 
            # save the decision
            if remind_txt == 'y':
                with open(join(PROVIDENTIA_ROOT, ".env"),"a") as f:
                    f.write(f"ORIGIN_UPDATE={origin_update_choice}\n")
        if origin_update_choice == 'n':
            # get files information
            files_info = yaml.safe_load(open(join(CURRENT_PATH, path)))
            files_info = {k: v for k, v in files_info.items() if k.strip() and v}
        else:
            # get files information
            combined_data = get_files_per_var(var)
            all_files = combined_data[var]['files']
            files_info = get_files_info(all_files, var, path)
    
    # go to next variable if no data is found
    if len(files_info) == 0:
        continue
        
    # filter files by resolution and dates
    print('    Filtering files by resolution and dates...')
    files = []
    for file, attributes in files_info.items():
        if attributes["resolution"] == resolution:
            start_date = datetime.datetime.strptime(attributes["start_date"], "%Y-%m-%dT%H:%M:%S UTC")
            end_date = datetime.datetime.strptime(attributes["end_date"], "%Y-%m-%dT%H:%M:%S UTC")
            for file_to_download in files_to_download:
                file_to_download_yearmonth = file_to_download.split(f'{var}_')[1].split('.nc')[0]
                file_to_download_start_date = datetime.datetime.strptime(file_to_download_yearmonth, "%Y%m")
                file_to_download_end_date = datetime.datetime(file_to_download_start_date.year, file_to_download_start_date.month, 1) + relativedelta(months=1, seconds=-1)
                if file_to_download_start_date <= end_date and file_to_download_end_date >= start_date:
                    if file not in files:
                        files.append(file)
    
    if len(files) != 0:
            
        # get data and metadata for each file within period
        combined_ds_list, metadata, wavelength = get_data(files, var, actris_parameter, resolution, target_start_date, target_end_date)

        # get flag dimension per station
        N_flag_codes_dims = []
        for ds in combined_ds_list:
            N_flag_codes_dims.append(ds.dims['N_flag_codes'])
        
        # get maximum number of flags across all stations
        N_flag_codes_max = max(N_flag_codes_dims)
        
        # recreate flag variable so that all stations have the same dimension and can be concatenated, leave nan for unknown values
        combined_ds_list_corrected_flag = []
        for ds in combined_ds_list:
            flag_data = ds['flag']
            da_flag = xr.DataArray(
                    np.full((flag_data.sizes['station'], flag_data.sizes['time'], N_flag_codes_max), np.nan),
                    dims=["station", "time", "N_flag_codes"],
                    coords={
                        "time": flag_data.coords["time"],
                    },
                    name="flag"
                )
            da_flag[:, :, :flag_data.values.shape[-1]] = flag_data.values
            ds = ds.drop_vars('flag')
            ds['flag'] = da_flag
            combined_ds_list_corrected_flag.append(ds)

        # combine and create new dataset
        print('    Combining files...')
        try:
            combined_ds = xr.concat(combined_ds_list_corrected_flag, 
                                    dim='station', 
                                    combine_attrs='drop_conflicts')
        except Exception as error:
            print(f'Error: Datasets could not be combined - {error}')
            continue
        
        # add metadata
        for key, value in metadata[resolution].items():
            if key in ['latitude', 'longitude']:
                value = [float(val) for val in value]
            elif key in ['altitude', 'sampling_height']:
                value = [float(val.replace('m', '').strip()) if isinstance(val, str) else val for val in value]
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
        combined_ds.attrs['data_license'] = 'BSD-3-Clause. Copyright 2025 Alba Vilanova Cortezón'
        combined_ds.attrs['source'] = 'Observations'
        combined_ds.attrs['institution'] = 'Barcelona Supercomputing Center'
        combined_ds.attrs['creator_name'] = 'Alba Vilanova Cortezón'
        combined_ds.attrs['creator_email'] = 'alba.vilanova@bsc.es'
        combined_ds.attrs['application_area'] = 'Monitoring atmospheric composition'
        combined_ds.attrs['domain'] = 'Atmosphere'
        combined_ds.attrs['observed_layer'] = 'Land surface'

        # save data per year and month
        path = join(nonghost_root, f'actris/actris/{resolution}/{var}')
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
        saved_files = 0
        for year, ds_year in combined_ds.groupby('time.year'):
            for month, ds_month in ds_year.groupby('time.month'):
                filename = f"{path}/{var}_{year}{month:02d}.nc"
                if filename in files_to_download:
                    combined_ds_yearmonth_unaveraged = combined_ds.sel(time=f"{year}-{month:02d}")
                    combined_ds_yearmonth = temporally_average_data(combined_ds_yearmonth_unaveraged, resolution, year, month, var)
    
                    # add title to attrs
                    extra_info = ''
                    wavelength_var = is_wavelength_var(actris_parameter)
                    if wavelength_var and wavelength is not None:
                        extra_info = f' at {wavelength}nm'
                    combined_ds_yearmonth.attrs['title'] = f'Surface {parameters_dict[var]}{extra_info} in the ACTRIS network in {year}-{month:02d}.'
    
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
                    #     print(f'    Data for {n_stations_diff} stations was removed because all data was NaN during {month}-{year}.')
                    
                    # save file
                    combined_ds_yearmonth.to_netcdf(filename)

                    # change permissions
                    os.system("chmod 777 {}".format(filename))
                    print(f"    Saved: {filename}")
                    saved_files += 1
                    
        print(f'    Total number of saved files: {saved_files}')

    else:
        print('    No files were found')