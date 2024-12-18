
import datetime
import itertools
import os
import requests
import yaml

import numpy as np
import pandas as pd
import xarray as xr

from providentia.auxiliar import CURRENT_PATH, join

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

# load ACTRIS mapping files
parameters_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'ghost_actris_variables.yaml')))
metadata_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'metadata.yaml')))
coverages_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'coverages.yaml')))
units_dict = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'units.yaml')))
variable_mapping = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'variable_mapping.yaml')))
variable_mapping = {k: v for k, v in variable_mapping.items() if k.strip() and v}


def get_files_per_var(var):
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


def get_files_info(var, files):

    files_info = {}
    files_info[var] = {}
    for i, file in enumerate(files):
        try:
            ds = xr.open_dataset(file)
        except:
            print(i, '-', file, '- Error: Could not open dataset')
            continue
        coverage = ds.time_coverage_resolution
        try:
            file_resolution = coverages_dict[coverage]
        except:
            print(i, '-', file, '- Error: Unknown coverage:', coverage)
            continue
        start_date = ds.time_coverage_start
        end_date = ds.time_coverage_end
        variables = list(ds.data_vars.keys())
        files_info[var][file] = {}
        files_info[var][file]['resolution'] = file_resolution
        files_info[var][file]['start_date'] = start_date
        files_info[var][file]['end_date'] = end_date
        files_info[var][file]['variables'] = variables
        print(i, '-', file, '- OK')

    return files_info


def get_files_path(var):

    alpha_var = ''.join(x for x in var if x.isalpha())
    if alpha_var in ['lsco', 'absco', 'lbsco', 'odaero']:
        path = join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', f'files/{alpha_var}/files.yaml')
    else:
        path = join(PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', f'files/{var}/files.yaml')

    return path


def create_files_info_file(var):

    path = get_files_path(var)
    print('Searching for data in Thredds...')
    combined_data = get_files_per_var(var)

    # if file does not exist
    if not os.path.isfile(path):

        print(f'File {path} does not exist, creating.')

        # get files information
        files = combined_data[var]['files']
        print(f'Reading basic information from {len(files)} files')
        files_info = get_files_info(var, files)

        # create file
        datasets = {
            url: data
            for url, data in files_info[var].items()
        }
        if len(datasets) != 0:
            path_dir = os.path.dirname(path)
            if not os.path.exists(path_dir):
                os.makedirs(path_dir)
            with open(path, 'w') as file:
                yaml.dump(datasets, file, default_flow_style=False)

    # if file exists
    else:
        print(
            f'File {path} already exists, checking if it needs an update.')

        # get currently available files
        current_files = combined_data[var]['files']

        # get previously saved files
        file_info_to_update = yaml.safe_load(
            open(join(CURRENT_PATH, path)))
        previous_files = list(file_info_to_update.keys())

        # keep old files only if they are currently available
        files_to_remove = []
        for file in previous_files:
            if file not in current_files:
                files_to_remove.append(file)

        # add new files
        files_to_add = []
        for file in current_files:
            if file not in previous_files:
                files_to_add.append(file)

        # get information from new files
        datasets = []
        if len(files_to_add) > 0:
            print(
                '- Files to add ({0}): {1}'.format(len(files_to_add), files_to_add))
            print('Reading basic information...')
            files_info = get_files_info(var, files_to_add)
            datasets = {
                url: data
                for url, data in files_info[var].items()
            }

            # add new data to dictionary
            if len(datasets) != 0:
                file_info_to_update.update(datasets)

        # remove unavailable data
        if len(files_to_remove) > 0:
            print(
                '- Files to remove ({0}): {1}'.format(len(files_to_remove), files_to_remove))
            for file in files_to_remove:
                print(f'Removing file {file}')
                file_info_to_update.pop(file, None)

        # recreate file info
        if len(datasets) != 0 or len(files_to_remove) > 0:
            print('Updating file...')
            os.remove(path)
            with open(path, 'w') as file:
                yaml.dump(file_info_to_update, file,
                            default_flow_style=False)
        else:
            print('No relevant changes were found.')


def filter_files(var, resolution, target_start_date, target_end_date):
    
    create_files_info_file(var)

    files = []
    path = get_files_path(var)
    files_info = yaml.safe_load(open(join(CURRENT_PATH, path)))
    files_info = {k: v for k, v in files_info.items() if k.strip() and v}
    for file, attributes in files_info.items():
        if attributes["resolution"] == resolution:
            start_date = datetime.datetime.strptime(
                attributes["start_date"], "%Y-%m-%dT%H:%M:%S UTC")
            end_date = datetime.datetime.strptime(
                attributes["end_date"], "%Y-%m-%dT%H:%M:%S UTC")
            if start_date <= target_end_date and end_date >= target_start_date:
                files.append(file)
    return files


def temporally_average_data(combined_ds, resolution, year, month, var):

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
    first_day_next_month = datetime.datetime(
        year, month % 12 + 1, 1) if month != 12 else datetime.datetime(year + 1, 1, 1)
    end_date = first_day_next_month - datetime.timedelta(days=1)
    valid_dates = pd.date_range(
        start=start_date, end=end_date, freq=frequency).to_numpy(dtype='datetime64[ns]')

    # initialise averaged data
    averaged_data = np.empty(
        (len(combined_ds.station.values), len(valid_dates)))

    for station_i, station in enumerate(combined_ds.station.values):
        # initialise averaged data
        station_averaged_data = []

        # read data per station
        data = combined_ds[var].isel(station=station_i).values

        # ignore data (times and values) if the values are nan
        valid_idxs = ~np.isnan(data)
        valid_time = time[valid_idxs]
        valid_data = data[valid_idxs]

        # calculate weighted averages
        if len(valid_data) != 0:
            for date in valid_dates:

                # get differences between valid time and actual times in minutes
                time_diffs = (
                    valid_time - date).astype('timedelta64[ns]').astype(float)

                # get positive differences and negative differences to differentiate
                # between the actual times that are earlier than the valid date (negative), and those that are later (positive)
                positive_diffs = time_diffs[time_diffs > 0]
                negative_diffs = time_diffs[time_diffs < 0]

                # find the closest actual time after the valid time
                closest_positive = None
                if len(positive_diffs) > 0:
                    closest_positive_idx = np.abs(positive_diffs).argmin()
                    closest_positive = positive_diffs[closest_positive_idx]
                    closest_positive_time = valid_time[time_diffs ==
                                                       positive_diffs[closest_positive_idx]][0]
                    closest_positive_value = valid_data[time_diffs ==
                                                        positive_diffs[closest_positive_idx]][0]

                # find the closest actual time before the valid time
                closest_negative = None
                if len(negative_diffs) > 0:
                    closest_negative_idx = np.abs(negative_diffs).argmin()
                    closest_negative = negative_diffs[closest_negative_idx]
                    closest_negative_time = valid_time[time_diffs ==
                                                       negative_diffs[closest_negative_idx]][-1]
                    closest_negative_value = valid_data[time_diffs ==
                                                        negative_diffs[closest_negative_idx]][-1]

                # when the valid time only has a value in one direction, get closest value without calculating weights
                if closest_positive is None:
                    value = closest_negative_value
                elif closest_negative is None:
                    value = closest_positive_value
                # in the rest of cases, calculate weights of 2 closest values and make average
                else:
                    # get 2 closest times and make positive to be able to compare differences
                    closest_diffs = np.abs(
                        [closest_negative, closest_positive])

                    # we do the reverse, since we want the differences in minutes to have a heavier weight if these are smaller (nearer the actual time)
                    weights = 1 / closest_diffs

                    # finally we normalize them to have values between 0 and 1
                    weights_normalized = weights / np.sum(weights)

                    # get average
                    value = np.average(
                        [closest_negative_value, closest_positive_value], weights=weights_normalized)

                # save averaged data
                station_averaged_data.append(value)

            averaged_data[station_i, :] = station_averaged_data
        else:
            averaged_data[station_i, :] = [np.nan]*len(valid_dates)

    # create new variable with averaged data
    combined_averaged_ds = xr.DataArray(
        data=averaged_data,
        coords={'station': combined_ds.station.values, 'time': valid_dates},
        dims=['station', 'time'],
        attrs={'units': combined_ds[var].units})

    # drop old variable and associated time
    combined_ds = combined_ds.drop_vars(var)
    combined_ds = combined_ds.drop_dims('time')

    # add new variable
    combined_ds[var] = combined_averaged_ds

    return combined_ds
