""" Class for downloading and formatting Actris data """

import bisect
import copy
import ctypes
import datetime
import itertools
import multiprocessing
import os
import re
import sys

from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import pycountry
import requests
from tqdm import tqdm
import xarray as xr
import yaml

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

class Actris:

    def __init__(self, download_instance, resolution):
        """
        Initialise the object

        Parameters
        ----------
        download_instance : object
            Instance used to download data
        resolution : str
            Desired output resolution
        """

        self.download_instance = download_instance
        self.resolution = resolution

    def get_files_per_var(self, base_url, var):
        """
        Get all files available in ACTRIS server per variable

        Parameters
        ----------
        base_url : str
            URL to NILU Thredds where data is stored
        var : str
            Variable

        Returns
        -------
        dict
            Dictionary with files per variable
        """

        files_per_var = {}

        if var not in files_per_var:
            files_per_var[var] = {}

        variable_files = []
        page = 0
        while True:
            # set up URL with pagination
            url = f"{base_url}/{ghost_actris_variables[var]}/page/{page}"
            try:
                response = requests.get(url)
            except:
                self.download_instance.logger.error(
                    f"Error connecting to server at {url}")
                sys.exit()

            # check if the response is valid and contains data
            if response.status_code != 200:
                self.download_instance.logger.error(
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

    def get_files_path(self, var):
        """
        Get path of file where files per variable are saved

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

    def get_standard_flags_and_qa(self, flags, ghost_version):
        """
        Convert flags from EBAS standards to GHOST standards and get QA

        Parameters
        ----------
        flags : numpy.array
            Array with flags in EBAS standards
        ghost_version : float
            GHOST version

        Returns
        -------
        numpy.array
            Array with flags in GHOST standards
        numpy.array
            Array with QA in GHOST standards
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

    def create_time_pairs(self, time):
        """
        Build pairs of consecutive time values

        Parameters
        ----------
        time : numpy.array
            Standard datetimes

        Returns
        -------
        list
            Pairs of standard datetimes
        """

        time_pairs = list(zip(time[:-1], time[1:]))
        last_time_pair = (time[-1], time[-1] + (time[-1] - time[-2]))
        time_pairs.append(last_time_pair)

        return time_pairs

    def datetime_to_fractional_minutes_from_reference(self, date):
        """
        Get datetime in fractional minutes (e.g. 10 minutes 45 seconds beacomes 10.45)

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
                        
    def get_window_indices(self, standard_start_date, standard_end_date, valid_start_times, 
                           valid_end_times, last_relevant_index):
        """
        Determines which measurement intervals overlap with a given time window

        Parameters
        ----------
        standard_start_date : datetime.datetime
            Window standard start date
        standard_end_date : datetime.datetime
            Window standard end date
        valid_start_times : numpy.array
            Actual start times
        valid_end_times : numpy.array
            Actual end times
        last_relevant_index : int
            Index to start searching from to ignore earlier timesteps

        Returns
        -------
        numpy.array
            All measurement periods that overlap with the window (fully inside or overlapping edges)
        int
            Overlap time from the standard start date of the window
        int
            Overlap time from the standard end date of the window
        """

        window_start_minute = self.datetime_to_fractional_minutes_from_reference(standard_start_date)
        window_end_minute = self.datetime_to_fractional_minutes_from_reference(standard_end_date)

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
        start_time_indices_in_window = np.arange(valid_start_times_begin, valid_start_times_end, 
                                                 dtype=np.uint32)
        end_time_indices_in_window = np.arange(valid_end_times_begin, valid_end_times_end, dtype=np.uint32)
        measurement_in_window_indices = np.intersect1d(start_time_indices_in_window,
                                                       end_time_indices_in_window, 
                                                       assume_unique=True)
        
        # get indices of measurements which overlap on left/right edges of window
        left_overlap_indices = np.setdiff1d(end_time_indices_in_window, measurement_in_window_indices, 
                                            assume_unique=True)
        left_overlap_indices = np.setdiff1d(left_overlap_indices, window_in_measurement_indices, 
                                            assume_unique=True)
        right_overlap_indices = np.setdiff1d(start_time_indices_in_window, measurement_in_window_indices, 
                                             assume_unique=True)
        right_overlap_indices = np.setdiff1d(right_overlap_indices, window_in_measurement_indices, 
                                             assume_unique=True)

        # deal with cases where left/right borders align but measurement period entirely contains window
        # add indices of these cases to window_in_measurement_indices
        # remove these indices also from the overlap indices
        left_window_in_measurement_indices = right_overlap_indices[np.where(
            valid_start_times[right_overlap_indices] == window_start_minute)]
        right_window_in_measurement_indices = left_overlap_indices[np.where(
            valid_end_times[left_overlap_indices] == window_end_minute)]
        if len(left_window_in_measurement_indices) > 0:
            window_in_measurement_indices = np.concatenate((window_in_measurement_indices, 
                                                            left_window_in_measurement_indices))
            right_overlap_indices = np.setdiff1d(right_overlap_indices, 
                                                 left_window_in_measurement_indices, 
                                                 assume_unique=True)
        if len(right_window_in_measurement_indices) > 0:
            window_in_measurement_indices = np.concatenate((window_in_measurement_indices, 
                                                            right_window_in_measurement_indices))
            left_overlap_indices = np.setdiff1d(left_overlap_indices, 
                                                right_window_in_measurement_indices, 
                                                assume_unique=True)
        
        # if there is a left border overlap, get the number of minutes the measurement period 
        # overlaps the measurement window
        if len(left_overlap_indices) > 0:
            left_overlap = valid_end_times[left_overlap_indices[0]] - window_start_minute
        # otherwise left overlap == 0
        else:
            left_overlap = 0
        
        # if there is a right border overlap, get the number of minutes the measurement period 
        # overlaps the measurement window
        if len(right_overlap_indices) > 0:
            right_overlap = window_end_minute - valid_start_times[right_overlap_indices[0]]
        # otherwise right overlap == 0
        else:
            right_overlap = 0
        
        # concatenate all relevant measurements indices in current window
        window_indices = np.sort(np.concatenate((measurement_in_window_indices,
                                                 window_in_measurement_indices,
                                                 left_overlap_indices,
                                                 right_overlap_indices)))
                
        return window_indices, right_overlap, left_overlap

    def temporally_average_data(self, station_ds, var, ghost_version, standard_time_pairs, vfunc):
        """
        Temporally average data and get unique flags in the valid times (temporally averaged)

        Parameters
        ----------
        station_ds : xarray.Dataset
            Data per station
        var : str
            Variable
        ghost_version : float
            GHOST version
        standard_time_pairs : list
            Pairs of standard datetimes
        vfunc : callable
            Vectorised function to map ghost validity flags
        
        Returns
        -------
        list
            Station temporally averaged data
        list
            Station temporally averaged flag data
        list
            Station temporally averaged QA data
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
        valid_timedeltas = np.array([(end_time - start_time).astype('timedelta64[m]').astype(np.float32) 
                                    for (end_time, start_time) in zip(end_times, start_times)])
        
        # get measurement start and end times as integers
        valid_start_times = np.array([self.datetime_to_fractional_minutes_from_reference(t) 
                                      for t in pd.to_datetime(start_times).to_pydatetime()])
        valid_end_times = np.array([self.datetime_to_fractional_minutes_from_reference(t) 
                                    for t in pd.to_datetime(end_times).to_pydatetime()])

        # initialise variable for finding relevant measurement indices
        last_relevant_index = 0

        for i, (standard_start_date, standard_end_date) in enumerate(standard_time_pairs):
            
            # get window indices
            window_indices, right_overlap, left_overlap = self.get_window_indices(standard_start_date, 
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
                window_flag_data, window_qa_data = self.get_standard_flags_and_qa(
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
                window_flag_data, window_qa_data = self.get_standard_flags_and_qa(
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

    def is_wavelength_var(self, actris_parameter):
        """
        Check if ACTRIS parameter depends on wavelength

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

    def get_files_info(self, files, var, path):
        """
        Read variables, resolution, start date and end date from all files in ACTRIS server 
        per variable

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

        if len(files) == 0:
            return

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
            self.download_instance.logger.error(f'Error: No data could be found for {var}')
            return
        
        return files_info


    def get_var_in_file(self, ds, var, actris_parameter, ebas_component):
        """
        Get variable name from dataset

        Parameters
        ----------
        ds : xarray.Dataset
            Dataset
        var : str
            GHOST variable name
        actris_parameter : str
            ACTRIS vocabulary name
        ebas_component : str
            EBAS vocabulary name

        Returns
        -------
        list
            List of possible variable names
        str
            Actual variable name found in dataset
        """

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

    def select_station_file(self, urls, files_info):
        """
        Select one station file from available ones, depending on 
        statistics, data level and revision date

        Parameters
        ----------
        urls : list
            List of available URLs for the same station
        files_info : dict
            Dictionary with information of all files in Thredds for one variable

        Returns
        -------
        str
            Selected URL
        """

        attrs_dict = {}
        urls = np.array(urls)
        # create dictionary with information about statistics, data level and revision date of 
        # available files
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
                max_level_ind = np.nanargmax(np.float32([np.nan if x == '' else x for x 
                                                         in attrs_dict['ebas_data_level']]))
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

    def str_to_nan(self, val):
        """
        Detect strings that only contain commas

        Parameters
        ----------
        val : str
            Metadata value

        Returns
        -------
        bool
            Whether string contains only commas or not
        """

        stripped = val.replace(" ", "")
        if stripped and all(c == "," for c in stripped):
            return True
        return False

    def read_data(self, args):
        """
        Read data for one station and extract valuable information

        Parameters
        ----------
        args : tuple
            Tuple containing all necessary information to process data for one station

        Returns
        -------
        str
            Selected URL
        str
            Station code
        str
            Station errors
        str
            Station warnings
        """

        (i, station, urls, data_shape, flag_shape, qa_shape,
        metadata_shape, files_info, var, actris_parameter, ebas_component,
        target_start_date, target_end_date,
        variable_mapping, metadata_dict, standard_time_pairs,
        vfunc, ghost_version) = args

        local_errors = ""
        local_warnings = ""

        # open files
        if len(urls) > 1:
            try:
                url = self.select_station_file(urls, files_info)
                if url is None:
                    local_warnings = f"No station file is valid."
                    return urls, station, local_errors, local_warnings
            except Exception as error:
                local_errors = f'Selecting station file: {error}.'
                return urls, station, local_errors, local_warnings
        else:
            url = urls[0]

        try:
            nc = Dataset(url, mode='r')
            ds = xr.open_dataset(xr.backends.NetCDF4DataStore(nc))
        except Exception as error:
            local_errors = f'Opening file: {error}.'
            return url, station, local_errors, local_warnings
        
        possible_vars, possible_var = self.get_var_in_file(ds, var, actris_parameter, ebas_component)
        if possible_var is None:
            local_errors = f'No variable name matches for {possible_vars}. Existing keys: {list(ds.data_vars)}.'
            return url, station, local_errors, local_warnings
    
        # remove time duplicates if any (keep first)
        ds = ds.sel(time=~ds['time'].to_index().duplicated())
        
        # select data in period range
        ds = ds.sel(time=slice(target_start_date, target_end_date))
        if ds.time.size == 0:
            local_warnings += f'No data available after filtering by time.'
            return url, station, local_errors, local_warnings

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
                return url, station, local_errors, local_warnings

        # read variable
        da_var = ds[possible_var]

        # avoid datasets that do not have defined units
        if 'ebas_unit' not in da_var.attrs:
            local_errors = f'No units were defined.'
            return url, station, local_errors, local_warnings

        # avoid datasets that do not have the same units as in variable mapping
        if da_var.attrs['ebas_unit'] != variable_mapping[actris_parameter]['units']:
            local_errors = f"Units {da_var.attrs['ebas_unit']} do not match those in variable mapping "
            local_errors += f"dictionary ({variable_mapping[actris_parameter]['units']})."
            return url, station, local_errors, local_warnings

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
            return url, station, local_errors, local_warnings

        # rename variable to BSC standards
        station_ds = da_var.to_dataset(name=var)

        # add time bounds
        station_ds['time_bnds'] = ds['time_bnds']

        # add quality control data
        station_ds['flag'] = flag_data

        # temporally average data from original times to standard times
        station_averaged_data, station_flag_data, station_qa_data = self.temporally_average_data(
            station_ds, var, ghost_version, standard_time_pairs, vfunc)
        if np.isnan(station_averaged_data).all():
            local_warnings += 'No data after temporal averaging.'
            return url, station, local_errors, local_warnings

        data_np = np.frombuffer(shared_memory_vars['data'], dtype=np.float32).reshape(data_shape)
        flag_np = np.frombuffer(shared_memory_vars['flag'], dtype=np.float32).reshape(flag_shape)
        qa_np = np.frombuffer(shared_memory_vars['qa'], dtype=np.float32).reshape(qa_shape)

        data_np[i, :] = station_averaged_data
        flag_np[i, :, :] = station_flag_data
        qa_np[i, :, :] = station_qa_data

        # get metadata from DOI in API if possible
        doi = ds.attrs['doi']
        metadata = self.get_metadata(doi)
        if len(metadata) != 0:
            facility_metadata = metadata[0]['md_data_identification']['facility']
            specific_metadata = metadata[0]['md_actris_specific']
            contact_metadata = metadata[0]['md_metadata']['contact'][0]
            contact_identification = metadata[0]['md_identification']['contact'][0]
        # else:
        #     local_warnings += f"Metadata cannot be read from DOI {doi} in API, reading from attributes.."

        # save metadata
        metadata_np = {}
        for ghost_key, actris_key in metadata_dict.items():
            
            metadata_np[ghost_key] = np.frombuffer(shared_memory_vars['metadata'][ghost_key], dtype='S75').reshape(metadata_shape)

            actris_api_keys = actris_key['API']
            actris_attr_keys = actris_key['attributes']
            
            # search in API first
            vals = []
            if len(metadata) != 0:
                for actris_api_key in actris_api_keys:
                    val = None
                    if actris_api_key in facility_metadata.keys() and len(actris_api_key) > 0:
                        val = facility_metadata[actris_api_key]
                    elif actris_api_key in specific_metadata.keys() and len(actris_api_key) > 0:
                        val = specific_metadata[actris_api_key]
                    elif ghost_key in ['contact_email_address', 'contact_institution', 'contact_name']:
                        if ghost_key == 'contact_name':
                            val = f"{contact_metadata['first_name']} {contact_metadata['last_name']}" 
                        elif len(actris_api_key) > 0:
                            val = contact_metadata[actris_api_key]
                    elif ghost_key in ['principal_investigator_email_address', 
                                    'principal_investigator_institution',
                                    'principal_investigator_name']:
                        if ghost_key == 'principal_investigator_name':
                            val = f"{contact_identification['first_name']} {contact_identification['last_name']}"
                        elif len(actris_api_key) > 0:
                            val = contact_identification[actris_api_key]
                    if val is not None:
                        vals.append(val)

            # search in attributes if it does not exist in API
            if len(vals) == 0:
                for actris_attr_key in actris_attr_keys:
                    val = None
                    if len(actris_attr_key) > 0:
                        if actris_attr_key in da_var_attrs.keys():
                            val = da_var_attrs[actris_attr_key]
                        elif actris_attr_key in ds.attrs.keys():
                            val = ds.attrs[actris_attr_key]
                    if val is not None:
                        vals.append(val)

            # if not found, make nan
            if len(vals) == 0:
                metadata_val = str(np.nan)
            # if found and only one value
            elif len(vals) == 1:
                vals = vals[0]
                # if it is a string made out of commas, make nan
                if isinstance(vals, str) and self.str_to_nan(vals):
                    vals = str(np.nan)
                # keep everything as string
                elif not isinstance(vals, str):
                    vals = str(vals)
                metadata_val = vals
            # convert lists into strings
            elif len(vals) > 1:
                metadata_val = ", ".join(vals).lstrip(", ")
            
            # remove all leading character ,
            metadata_val = metadata_val.lstrip(", ")

            metadata_np[ghost_key][i] = metadata_val.encode('utf-8')

        return url, station, local_errors, local_warnings

    def init_shared_vars_read_data(self, shared_data, shared_flag_data, shared_qa_data, shared_metadata, 
                                   data_shape):
        """
        Initialise shared arrays across processes

        Parameters
        ----------
        shared_data : multiprocessing.RawArray
            Shared data array
        shared_flag_data : multiprocessing.RawArray
            Shared flag data array
        shared_qa_data : multiprocessing.RawArray
            Shared QA data array
        shared_metadata : multiprocessing.RawArray
            Shared metadata array
        data_shape : tuple
            Data shape
        """

        # multiprocessing.RawArray does not support NaN initialisation and initialised to zeros
        # replace 0 by nan before reading data so that if there are any errors we can later 
        # drop the stations
        shared_data = np.frombuffer(shared_data, dtype=np.float32).reshape(data_shape)
        shared_data[:] = np.nan
        shared_memory_vars['data'] = shared_data
        shared_memory_vars['flag'] = shared_flag_data
        shared_memory_vars['qa'] = shared_qa_data
        shared_memory_vars['metadata'] = {}
        for ghost_key in metadata_dict.keys():
            shared_memory_vars['metadata'][ghost_key] = shared_metadata[ghost_key]

    def ghost_validity_flag_mapper(self, flag):
        """Get GHOST decreed validity from flag

        Parameters
        ----------
        flag : float
            GHOST flag

        Returns
        -------
        str
            GHOST decreed validity
        """

        if np.isnan(flag):
            return np.nan
        elif flag == 0:
            return flags_dict['000']["GHOST_decreed_validity"][0]
        else:
            return flags_dict[str(int(flag))]["GHOST_decreed_validity"][0]

    def get_metadata(self, doi):
        """
        Get metadata dictionary from DOI

        Parameters
        ----------
        doi : str
            NILU's DOI

        Returns
        -------
        dict
            Metadata dictionary
        """

        url = "https://prod-actris-md2.nilu.no/metadata/pid"

        headers = {
        "accept": "application/json",
        "Content-Type": "application/json"
        }
        
        data = {"pid": doi}

        response = requests.post(url, headers=headers, json=data)
        metadata = response.json()
        
        return metadata

    def get_data(self, files, var, actris_parameter, target_start_date, target_end_date, files_info):
        """
        Read variable and metadata data standarising dimensions

        Parameters
        ----------
        files : list
            Files per variable
        var : str
            Variable
        actris_parameter : str
            ACTRIS parameter
        target_start_date : datetime.datetime
            Target start date (defined from configuration file)
        target_end_date : datetime.datetime
            Target end date (defined from configuration file)
        files_info : dict
            Dictionary with information of all files in Thredds for one variable
        
        Returns
        -------
        xarray.Dataset
            Temporally averaged data with flags, QA and metadata
        """

        # get EBAS component
        ebas_component = variable_mapping[actris_parameter]['var']

        # get valid dates frequency
        if self.resolution == 'hourly':
            frequency = 'h'
        elif self.resolution == '3hourly':
            frequency = '3h'
        elif self.resolution == 'daily':
            frequency = 'D'
        # TODO: Review this
        elif self.resolution == 'monthly':
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
        vfunc = np.vectorize(self.ghost_validity_flag_mapper, otypes=['O'])

        # create pairs of valid dates
        standard_time_pairs = self.create_time_pairs(standard_time)

        # create data specific arrays to share across processes (for parallel multiprocessing use)
        shared_data = multiprocessing.RawArray(ctypes.c_float, int(np.prod(data_shape)))
        shared_flag_data = multiprocessing.RawArray(ctypes.c_float, int(np.prod(flag_shape)))
        shared_qa_data = multiprocessing.RawArray(ctypes.c_float, int(np.prod(qa_shape)))
        shared_metadata = {}
        for ghost_key in metadata_dict.keys():
            shared_metadata[ghost_key] = multiprocessing.RawArray(ctypes.c_char, int(np.prod(
                metadata_shape)*75))

        # read data and metadata in parallel
        errors = []
        warnings = []
        args_list = [
            (
                i, station, urls, data_shape, flag_shape, qa_shape, metadata_shape,
                files_info, var, actris_parameter, ebas_component,
                target_start_date, target_end_date,
                variable_mapping, metadata_dict, standard_time_pairs,
                vfunc, self.download_instance.ghost_version
            )
            for i, (station, urls) in enumerate(files.items())
        ]
        pool = multiprocessing.Pool(
            processes=self.download_instance.n_cpus,
            initializer=self.init_shared_vars_read_data,
            initargs=(shared_data, shared_flag_data, shared_qa_data, shared_metadata, data_shape)
        )
        for url, station, error, warning in tqdm(
                pool.imap(self.read_data, args_list),
                total=len(args_list),
                desc="    Processing station data",
                bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}'
            ):
            if error:
                errors.append(f'{url} ({station}): {error}')
            if warning:
                warnings.append(f'{url} ({station}): {warning}')

        # print errors and warnings if any
        if errors:
            self.download_instance.logger.info(f"    === ERRORS ({len(errors)}) ===")
            for e in errors:
                self.download_instance.logger.info(f"    {e}")

        if warnings:
            self.download_instance.logger.info(f"    === WARNINGS ({len(warnings)}) ===")
            for w in warnings:
                self.download_instance.logger.info(f"    {w}")
                
        pool.close()
        
        # wait for worker processes to terminate before continuing
        pool.join()

        if len(errors) == len(args_list):
            self.download_instance.logger.info('All datasets have thrown an error, aborting.')
            return

        # get combined data and metadata after read
        averaged_data = np.frombuffer(shared_data, dtype=np.float32).reshape(data_shape)
        averaged_flag_data = np.frombuffer(shared_flag_data, dtype=np.float32).reshape(flag_shape)
        averaged_qa_data = np.frombuffer(shared_qa_data, dtype=np.float32).reshape(qa_shape)
        metadata = {}
        for ghost_key in metadata_dict.keys():
            metadata[ghost_key] = np.frombuffer(shared_metadata[ghost_key], dtype='S75').reshape(
                metadata_shape)

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
                        if val != b'' and val != b'nan' else 0 for val in value]
            else:
                value = [val.decode('utf-8', errors="replace") for val in value]
            combined_ds[key] = xr.Variable(data=value, dims=('station'))
        
        # add country
        value = [pycountry.countries.get(alpha_2=val[0:2]).name if val else np.nan 
                 for val in combined_ds['country'].values]
        combined_ds['country'] = xr.Variable(data=value, dims=('station'))

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

    def get_files_to_download(self, target_start_date, target_end_date, var):
        """
        Get filenames that should be downloaded

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

        base_dir = join(self.download_instance.nonghost_root, 'actris/actris', self.resolution, var)
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

    def download_actris_data(self):
        """
        Download ACTRIS data
        """

        target_start_date = datetime.datetime(int(self.download_instance.start_date[:4]), 
                                              int(self.download_instance.start_date[4:6]), 
                                              int(self.download_instance.start_date[6:8]), 0)
        target_end_date = datetime.datetime(int(self.download_instance.end_date[:4]), 
                                            int(self.download_instance.end_date[4:6]), 
                                            int(self.download_instance.end_date[6:8]), 23, 59, 59) - datetime.timedelta(days=1)

        for var in self.download_instance.species:
            
            # check if variable name is available
            if var not in ghost_actris_variables.keys():
                self.download_instance.logger.info(f"Data for {var} cannot be downloaded because it was not mapped in 'settings/internal/actris/ghost_actris_variables.yaml'.")
                continue
            else:
                actris_parameter = ghost_actris_variables[var]

            # get files that were already downloaded
            initial_check_nc_files = self.get_files_to_download(target_start_date, target_end_date, var)
            files_to_download = self.download_instance.select_files_to_download(initial_check_nc_files)
            if not files_to_download:
                msg = f"Files were already downloaded for {var} at {self.resolution} "
                msg += f"resolution between {target_start_date} and {target_end_date}."
                show_message(self.download_instance, msg, deactivate=False)     
                continue 
            
            # get files info path
            info_path = self.get_files_path(var)

            # define NILU path
            base_url = "https://prod-actris-md.nilu.no/metadata/content"
            
            # if file does not exist
            if not os.path.isfile(info_path):
                # get files information
                self.download_instance.logger.info(f'\nFile containing information of the files available in Thredds for {var} ({info_path}) does not exist, creating.')
                combined_data = self.get_files_per_var(base_url, var)
                all_files = combined_data[var]['files']
                files_info = self.get_files_info(all_files, var, info_path)
                    
            # if file exists
            else:
                # make the user wants to update file information from NILU Thredds
                if not isinstance(self.download_instance.dl_thredds_update, bool):
                    # ask if user wants to update
                    while True:
                        dl_thredds_update = input(f"\nFile containing information of the files available in Thredds for {actris_parameter} ({info_path}) already exists. Do you want to update it (y/[n])? ").lower() 
                        if dl_thredds_update in ['y','n','']:
                            break
                    
                    # get the boolean value
                    self.download_instance.dl_thredds_update = dl_thredds_update not in ['n', '']

                if self.download_instance.dl_thredds_update:
                    # get files information
                    combined_data = self.get_files_per_var(base_url, var)
                    all_files = combined_data[var]['files']
                    files_info = self.get_files_info(all_files, var, info_path)
                else:
                    # get files information
                    files_info = yaml.safe_load(open(join(CURRENT_PATH, info_path)))
                    files_info = {k: v for k, v in files_info.items() if k.strip() and v}
            
            # go to next variable if no data is found
            if files_info is not None:
                if len(files_info) == 0:
                    print(f'Warning: No files found for {var}.')
                    continue
            else:
                print(f'Warning: No files found for {var}.')
                continue

            # get wavelength
            wavelength_var = self.is_wavelength_var(actris_parameter)
            if wavelength_var:
                # select most common wavelength for black carbon (name does not provide it)
                if var == 'sconcbc':
                    wavelength = 880
                    self.download_instance.logger.info(f'Wavelength appears in dimensions. Selected wavelength: {wavelength}.')
                # get wavelength from variable name for other variables
                else:
                    wavelength = float(re.findall(r'\d+', var)[0])
            else:
                wavelength = None

            # filter files by resolution and dates
            self.download_instance.logger.info('\n'+'-'*40)
            self.download_instance.logger.info(f"\nDownloading ACTRIS framework data from EBAS DOI...")
            path = join(self.download_instance.nonghost_root, f'actris/actris/{self.resolution}/{var}')
            self.download_instance.logger.info(f"\n  - {path}, source: {base_url}/{ghost_actris_variables[var]}")
            files = {}
            for file, attributes in files_info.items():
                if attributes["resolution"] == self.resolution:
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
                                    self.download_instance.logger.error(f'Dataset has wavelength in its dimensions but wavelength is None. Revise if ACTRIS parameter ({actris_parameter}) is included in is_wavelength_var function.')
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
                combined_ds = self.get_data(files, var, actris_parameter, 
                                            target_start_date, target_end_date, files_info)
                if combined_ds is None:
                    continue

                # save data per year and month
                if not os.path.isdir(path):
                    os.makedirs(path, exist_ok=True)

                valid_nc_files = []
                for year, ds_year in combined_ds.groupby('time.year'):
                    for month, ds_month in ds_year.groupby('time.month'):
                        filename = f"{var}_{year}{month:02d}.nc"
                        filepath = f"{path}/{filename}"
                        if filepath in files_to_download:
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
                            #     self.download_instance.logger.info(f'Data for {n_stations_diff} stations was removed because all data was NaN during {month}-{year}.')
                            
                            # add acknowledgements
                            dois = ', '.join(x for x in combined_ds_yearmonth.doi.values if x not in ('', '[', ']'))
                            combined_ds_yearmonth.attrs['acknowledgements'] = f'This data is compiled by measurements from these DOI: {dois}'

                            # remove file if it exists
                            if os.path.isfile(filepath):
                                os.system("rm {}".format(filepath))

                            # do not save if empty
                            if len(combined_ds_yearmonth[var].values) == 0:
                                continue

                            # get last downloaded file in case there was a keyboard interrupt
                            self.download_instance.latest_nc_file_path = filepath

                            # save file
                            combined_ds_yearmonth.to_netcdf(filepath)

                            # change permissions
                            os.system("chmod 777 {}".format(filepath))
                            valid_nc_files.append(filename)
                
                # print download of valid files
                valid_nc_files_iter = tqdm(valid_nc_files,bar_format= '{l_bar}{bar}|{n_fmt}/{total_fmt}',desc=f"    Downloading files ({len(valid_nc_files)})")
                for filename in valid_nc_files_iter:
                    pass

            else:
                self.download_instance.logger.info(f'No files were found at {self.resolution} resolution for {var}. You can check what is available at {info_path}.')
