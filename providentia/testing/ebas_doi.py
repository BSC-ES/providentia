import os
import yaml
import xarray as xr
import numpy as np
import datetime
import sys
import time
from dotenv import dotenv_values
from dateutil.relativedelta import relativedelta

from providentia.auxiliar import join
from providentia.actris import (get_files_per_var, get_files_path, get_files_info,
                                get_files_to_download, is_wavelength_var, get_data,
                                temporally_average_data)

PROVIDENTIA_ROOT = '/home/avilanov/software/Providentia/'
nonghost_root = '/home/avilanov/data/providentia/obs/nonghost/'
CURRENT_PATH = os.getcwd()

parameters_dict = yaml.safe_load(open(join(
    PROVIDENTIA_ROOT, 'settings', 'internal', 'actris', 'ghost_actris_variables.yaml')))

variables = ['sconco3']
resolution = 'hourly'
target_start_date = datetime.datetime(2018, 1, 1, 0)
target_end_date = datetime.datetime(2019, 1, 1, 0)


def select_files_to_download(prov_start_time, nc_files_to_download):
    """ Returns the files that are not already downloaded. """
    # initialise list of non-downloaded files
    not_downloaded_files = []

    # get ssh user and password
    env = dotenv_values(join(PROVIDENTIA_ROOT, ".env"))
    overwrite_choice = env.get("OVERWRITE")

    if nc_files_to_download:
        # get the downloaded and not downloaded files
        not_downloaded_files = list(
            filter(lambda x: not os.path.exists(x), nc_files_to_download))
        downloaded_files = list(
            filter(lambda x: os.path.exists(x), nc_files_to_download))

        # get the files that were downloaded before the execution
        downloaded_before_execution_files = list(
            filter(lambda x: prov_start_time > os.path.getctime(x), downloaded_files))

        # if there was any file downloaded before the execution
        if downloaded_before_execution_files:
            # make the user choose between overwriting or not overwriting
            if overwrite_choice not in ['y', 'n']:
                # ask if user wants to overwrite
                while overwrite_choice not in ['y', 'n']:
                    overwrite_choice = input(
                        "\nThere are some files that were already downloaded in a previous download, do you want to overwrite them (y/n)? ").lower()
                # ask if user wants to remember the decision
                remind_txt = None
                while remind_txt not in ['y', 'n']:
                    remind_txt = input(
                        "\nDo you want to remember your decision for future downloads (y/n)? ").lower()
                # save the decision
                if remind_txt == 'y':
                    with open(join(PROVIDENTIA_ROOT, ".env"), "a") as f:
                        f.write(f"OVERWRITE={overwrite_choice}\n")
            # if user wants to overwrite then add the the files that were
            # downloaded before the execution as if they were never download
            if overwrite_choice == 'y':
                not_downloaded_files += downloaded_before_execution_files
            # change overwritten files boolean to True to indicate that some files were ignored
            else:
                overwritten_files_flag = True

    return not_downloaded_files


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
    initial_check_nc_files = get_files_to_download(
        nonghost_root, target_start_date, target_end_date, resolution, var)
    files_to_download = select_files_to_download(
        prov_start_time, initial_check_nc_files)
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
        print(
            f'\nFile containing information of the files available in Thredds for {var} ({path}) does not exist, creating.')
        combined_data = get_files_per_var(var)
        all_files = combined_data[var]['files']
        files_info = get_files_info(all_files, var, path)

    # if file exists
    else:
        # ask if user wants to update file information from NILU Thredds
        if origin_update_choice not in ['y', 'n']:
            while origin_update_choice not in ['y', 'n']:
                origin_update_choice = input(
                    f"\nFile containing information of the files available in Thredds for {var} ({path}) already exists. Do you want to update it (y/n)? ").lower()
            # ask if user wants to remember the decision
            remind_txt = None
            while remind_txt not in ['y', 'n']:
                remind_txt = input(
                    "\nDo you want to remember your decision for future downloads (y/n)? ").lower()
            # save the decision
            if remind_txt == 'y':
                with open(join(PROVIDENTIA_ROOT, ".env"), "a") as f:
                    f.write(f"ORIGIN_UPDATE={origin_update_choice}\n")
        if origin_update_choice == 'n':
            # get files information
            files_info = yaml.safe_load(open(join(CURRENT_PATH, path)))
            files_info = {k: v for k, v in files_info.items()
                          if k.strip() and v}
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
            start_date = datetime.datetime.strptime(
                attributes["start_date"], "%Y-%m-%dT%H:%M:%S UTC")
            end_date = datetime.datetime.strptime(
                attributes["end_date"], "%Y-%m-%dT%H:%M:%S UTC")
            for file_to_download in files_to_download:
                file_to_download_yearmonth = file_to_download.split(f'{var}_')[
                    1].split('.nc')[0]
                file_to_download_start_date = datetime.datetime.strptime(
                    file_to_download_yearmonth, "%Y%m")
                file_to_download_end_date = datetime.datetime(
                    file_to_download_start_date.year, file_to_download_start_date.month, 1) + relativedelta(months=1, seconds=-1)
                if file_to_download_start_date <= end_date and file_to_download_end_date >= start_date:
                    if file not in files:
                        files.append(file)

    if len(files) != 0:

        # get data and metadata for each file within period
        combined_ds_list, metadata, wavelength = get_data(
            files, var, actris_parameter, resolution, target_start_date, target_end_date)

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
                np.full(
                    (flag_data.sizes['station'], flag_data.sizes['time'], N_flag_codes_max), np.nan),
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
                value = [float(val.replace('m', '').strip()) if isinstance(
                    val, str) else val for val in value]
            combined_ds[key] = xr.Variable(data=value, dims=('station'))

        # calculate measurement_altitude if altitude and sampling_height exist
        if ('altitude' in combined_ds.keys()) and ('sampling_height' in combined_ds.keys()):
            value = combined_ds['altitude'].values + \
                combined_ds['sampling_height'].values
            combined_ds['measurement_altitude'] = xr.Variable(
                data=value, dims=('station'))

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
                    combined_ds_yearmonth_unaveraged = combined_ds.sel(
                        time=f"{year}-{month:02d}")
                    combined_ds_yearmonth = temporally_average_data(
                        combined_ds_yearmonth_unaveraged, resolution, year, month, var)

                    # add title to attrs
                    extra_info = ''
                    wavelength_var = is_wavelength_var(actris_parameter)
                    if wavelength_var and wavelength is not None:
                        extra_info = f' at {wavelength}nm'
                    combined_ds_yearmonth.attrs[
                        'title'] = f'Surface {parameters_dict[var]}{extra_info} in the ACTRIS network in {year}-{month:02d}.'

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
                    combined_ds_yearmonth = combined_ds_yearmonth.dropna(
                        dim="station", subset=[var], how="all")
                    combined_ds_yearmonth = combined_ds_yearmonth.assign_coords(
                        station=range(len(combined_ds_yearmonth.station)))
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
