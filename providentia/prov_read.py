from .reading import get_yearmonths_to_read, read_netcdf_data, read_netcdf_nonghost

import os
import gc
import copy
import datetime
import multiprocessing

import numpy as np
import pandas as pd
import seaborn as sns
from netCDF4 import Dataset


class DataReader:
    """Class that reads observational/experiment data in memory."""

    def __init__(self, read_instance, read_type='parallel'):
        self.read_instance = read_instance
        self.read_type = read_type

    def read_all(self):

        # check if reading GHOST or non-GHOST files
        self.read_instance.reading_nonghost = self.check_for_ghost()

        # get valid observational files in range
        self.get_valid_obs_files_in_date_range(self.read_instance.le_start_date.text(),
                                               self.read_instance.le_end_date.text())

        # update available experiment data dictionary
        self.get_valid_experiment_files_in_date_range()

        # setup read
        self.read_setup(self.read_instance.active_resolution, self.read_instance.active_start_date,
                        self.read_instance.active_end_date, self.read_instance.active_network,
                        self.read_instance.active_species, self.read_instance.active_matrix)

        # update dictionary of plotting parameters (colour and zorder etc.) for each data array
        self.update_plotting_parameters()

    def reset_data_in_memory(self):
        self.data_in_memory = {}
        self.plotting_params = {}

    def check_for_ghost(self):
        """ It checks whether the selected network comes from GHOST or not.
        In case of non-ghost, it disables ghost-related fields"""

        # if we're reading nonghost files, then disable fields
        if '*' in self.read_instance.cb_network.currentText():
            # self.disable_ghost_buttons()
            return True
        else:
            return False
            # self.enable_ghost_buttons()

    def valid_date(self, date_text):
        """define function that determines if a date string is in the correct format"""

        try:
            datetime.datetime.strptime(str(date_text), '%Y%m%d')
            return True
        except Exception as e:
            return False

    def read_setup(self, resolution, start_date, end_date, network, species, matrix):
        """Setup key variables for new read of observational/experiment
        data a time array and create arrays of unique station
        references/longitudes/latitudes.
        """

        # force garbage collection (to avoid memory issues)
        gc.collect()
        self.read_instance.reading_nonghost = False  # TODO: remoe line, was added for testing
        # series of actions not applicable when --offline
        if not self.read_instance.offline:
            self.read_instance.reading_nonghost = self.check_for_ghost()

            # set current time array, as previous time array
            self.read_instance.previous_time_array = self.read_instance.time_array
            # set current station references, as previous station references
            self.read_instance.previous_station_references = self.read_instance.station_references
            # set current relevant yearmonths, as previous relevant yearmonths
            self.read_instance.previous_relevant_yearmonths = self.read_instance.relevant_yearmonths

        # get N time chunks between desired start date and end date to set time array
        if (resolution == 'hourly') or (resolution == 'hourly_instantaneous'):
            self.active_frequency_code = 'H'
        elif (resolution == '3hourly') or (resolution == '3hourly_instantaneous'):
            self.active_frequency_code = '3H'
        elif (resolution == '6hourly') or (resolution == '6hourly_instantaneous'):
            self.active_frequency_code = '6H'
        elif resolution == 'daily':
            self.active_frequency_code = 'D'
        elif resolution == 'monthly':
            self.active_frequency_code = 'MS'
        str_active_start_date = str(start_date)
        str_active_end_date = str(end_date)
        self.read_instance.time_array = pd.date_range(start=datetime.datetime(int(str_active_start_date[:4]),
                                                                              int(str_active_start_date[4:6]),
                                                                              int(str_active_start_date[6:8])),
                                                      end=datetime.datetime(int(str_active_end_date[:4]),
                                                                            int(str_active_end_date[4:6]),
                                                                            int(str_active_end_date[6:8])),
                                                      freq=self.active_frequency_code)[:-1]

        if not self.read_instance.reading_nonghost:
            # get all relevant observational files
            file_root = '%s/%s/%s/%s/%s/%s_' % (self.read_instance.obs_root, network,
                                                self.read_instance.ghost_version, resolution,
                                                species, species)
        else:
            # get files from nonghost path
            file_root = '%s/%s/%s/%s/%s/%s_' % (self.read_instance.nonghost_root, network[1:].lower(),
                                                self.read_instance.selected_matrix,
                                                resolution, species, species)

        self.read_instance.relevant_yearmonths = np.sort([yyyymm for yyyymm in self.available_observation_data[
            network][resolution][matrix][species]])

        relevant_files = sorted([file_root+str(yyyymm)[:6]+'.nc'
                                 for yyyymm in self.read_instance.relevant_yearmonths])
        self.N_inds_per_month = np.array([np.count_nonzero(np.all(
            [self.read_instance.time_array >= datetime.datetime.strptime(str(start_yyyymm), '%Y%m%d'),
             self.read_instance.time_array < datetime.datetime.strptime(str(self.read_instance.relevant_yearmonths[month_ii + 1]), '%Y%m%d')],
            axis=0)) if month_ii != (len(self.read_instance.relevant_yearmonths) - 1) else np.count_nonzero(
            self.read_instance.time_array >= datetime.datetime.strptime(str(start_yyyymm), '%Y%m%d')) for month_ii, start_yyyymm in
                                          enumerate(self.read_instance.relevant_yearmonths)])

        self.read_instance.station_references = []
        self.station_longitudes = []
        self.station_latitudes = []
        if not self.read_instance.reading_nonghost:
            for relevant_file in relevant_files:
                ncdf_root = Dataset(relevant_file)
                self.read_instance.station_references = np.append(self.read_instance.station_references, ncdf_root['station_reference'][:])
                self.station_longitudes = np.append(self.station_longitudes, ncdf_root['longitude'][:])
                self.station_latitudes = np.append(self.station_latitudes, ncdf_root['latitude'][:])
                ncdf_root.close()
            self.read_instance.station_references, station_unique_indices = np.unique(self.read_instance.station_references, return_index=True)
            self.station_longitudes = self.station_longitudes[station_unique_indices]
            self.station_latitudes = self.station_latitudes[station_unique_indices]
        else:
            # first, try to take the data files and handle in case of daily files
            if os.path.exists(relevant_files[0]):
                ncdf_root = Dataset(relevant_files[0])
            else:
                relevant_files = sorted([file_root + str(yyyymm)[:8] + '.nc'
                                         for yyyymm in self.read_instance.relevant_yearmonths])
                ncdf_root = Dataset(relevant_files[0])
            self.read_instance.station_references = np.array(
                [st_name.tostring().decode('ascii').replace('\x00', '')
                 for st_name in ncdf_root['station_name'][:]], dtype=np.str)
            # get staion refs
            if "latitude" in ncdf_root.variables:
                self.station_longitudes = np.append(self.station_longitudes, ncdf_root['longitude'][:])
                self.station_latitudes = np.append(self.station_latitudes, ncdf_root['latitude'][:])
            else:
                self.station_longitudes = np.append(self.station_longitudes, ncdf_root['lon'][:])
                self.station_latitudes = np.append(self.station_latitudes, ncdf_root['lat'][:])
            ncdf_root.close()

        # update measurement units for species (take standard units from parameter dictionary)
        self.measurement_units = self.read_instance.parameter_dictionary[species]['standard_units']

        # set data variables to read (dependent on active data resolution)
        if not self.read_instance.reading_nonghost:
            if (resolution == 'hourly') or (resolution == 'hourly_instantaneous'):
                self.data_vars_to_read = [species, 'hourly_native_representativity_percent',
                                          'daily_native_representativity_percent',
                                          'monthly_native_representativity_percent',
                                          'annual_native_representativity_percent', 'hourly_native_max_gap_percent',
                                          'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                          'annual_native_max_gap_percent', 'day_night_code', 'weekday_weekend_code',
                                          'season_code', 'time']
            elif (resolution == '3hourly') or \
                    (resolution == '6hourly') or (resolution == '3hourly_instantaneous') or \
                    (resolution == '6hourly_instantaneous'):
                 self.data_vars_to_read = [species, 'daily_native_representativity_percent',
                                           'monthly_native_representativity_percent',
                                           'annual_native_representativity_percent',
                                           'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                           'annual_native_max_gap_percent', 'day_night_code', 'weekday_weekend_code',
                                           'season_code', 'time']
            elif resolution == 'daily':
                self.data_vars_to_read = [species, 'daily_native_representativity_percent',
                                          'monthly_native_representativity_percent',
                                          'annual_native_representativity_percent',
                                          'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                          'annual_native_max_gap_percent', 'weekday_weekend_code', 'season_code', 'time']
            elif resolution == 'monthly':
                self.data_vars_to_read = [species, 'monthly_native_representativity_percent',
                                          'annual_native_representativity_percent', 'monthly_native_max_gap_percent',
                                          'annual_native_max_gap_percent', 'season_code', 'time']
        else:
            self.data_vars_to_read = [species]

        # set data dtype
        self.data_dtype = [(key, np.float32) for key in self.data_vars_to_read]

    def read_data(self, data_label, start_date, end_date,
                  network, resolution, species, matrix):
        """Function that handles reading of observational/experiment data"""

        # force garbage collection (to avoid memory issues)
        gc.collect()

        # get relevant file start dates
        if data_label == 'observations':
            process_type = 'observations'
            if not self.read_instance.reading_nonghost:
                file_root = '%s/%s/%s/%s/%s/%s_' % (self.read_instance.obs_root, network,
                                                    self.read_instance.ghost_version,
                                                    resolution, species, species)
                relevant_file_start_dates = \
                    sorted(self.available_observation_data[network][resolution][matrix][species])
            else:
                # get files from nonghost path
                file_root = '%s/%s/%s/%s/%s/%s_' % (self.read_instance.nonghost_root,
                                                    network[1:].lower(), matrix,
                                                    resolution, species, species)
                relevant_file_start_dates = \
                    sorted(self.available_observation_data[network][resolution][matrix][species])

        else:
            process_type = 'experiment'
            file_root = \
                '%s/%s/%s/%s/%s/%s/%s_' % (self.read_instance.exp_root, self.read_instance.ghost_version, data_label,
                                           resolution, species, network, species)

            relevant_file_start_dates = sorted(self.available_experiment_data[data_label])

        # get data files in required date range to read, taking care not to re-read what has already been read
        yearmonths_to_read = get_yearmonths_to_read(relevant_file_start_dates, start_date, end_date)
        relevant_files = [file_root+str(yyyymm)[:6]+'.nc' for yyyymm in yearmonths_to_read]

        if not os.path.exists(relevant_files[0]):
            relevant_files = sorted([file_root + str(yyyymm)[:8] + '.nc' for yyyymm
                                     in self.read_instance.relevant_yearmonths])

        # check if data label in data in memory dictionary
        if data_label not in list(self.data_in_memory.keys()):
            # if not create empty array (filled with NaNs) to store species data and place it in the dictionary

            if process_type == 'observations':
                self.plotting_params['observations'] = {}
                if not self.read_instance.reading_nonghost:
                    self.data_in_memory[data_label] = np.full((len(self.read_instance.station_references),
                                                               len(self.read_instance.time_array)),
                                                              np.NaN, dtype=self.data_dtype)
                else:
                    self.data_in_memory[data_label] = np.full((len(self.read_instance.station_references),
                                                               len(self.read_instance.time_array)),
                                                              np.NaN, dtype=self.data_dtype[:1])
                self.metadata_in_memory = np.full((len(self.read_instance.station_references),
                                                   len(self.read_instance.relevant_yearmonths)),
                                                  np.NaN, dtype=self.read_instance.metadata_dtype)
                if self.read_instance.reading_nonghost:
                    tmp_ncdf = Dataset(relevant_files[0])
                    # create separate structure of nonghost metadata
                    nonghost_mdata_dtype = [('station_name', np.object), ('latitude', np.float),
                                            ('longitude', np.float), ('altitude', np.float)]
                    if "station_code" in tmp_ncdf.variables:
                        nonghost_mdata_dtype.append(('station_reference', np.object))
                    self.nonghost_metadata = np.full((len(self.read_instance.station_references)),
                                                     np.NaN, dtype=nonghost_mdata_dtype)

            # if process_type is experiment, get experiment specific grid edges from
            # first relevant file, and save to data in memory dictionary
            if process_type == 'experiment':
                self.data_in_memory[data_label] = np.full((len(self.read_instance.station_references),
                                                           len(self.read_instance.time_array)),
                                                          np.NaN, dtype=self.data_dtype[:1])
                self.plotting_params[data_label] = {}
                exp_nc_root = Dataset(relevant_files[0])
                self.plotting_params[data_label]['grid_edge_longitude'] = \
                    exp_nc_root['grid_edge_longitude'][:]
                self.plotting_params[data_label]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
                exp_nc_root.close()

        # iterate and read species data in all relevant netCDF files (either in serial/parallel)

        # read serially
        if self.read_type == 'serial':

            # iterate through relevant netCDF files
            for relevant_file in relevant_files:
                # create argument tuple of function
                tuple_arguments = relevant_file, self.read_instance.time_array, self.read_instance.station_references, \
                                  species, process_type,\
                                  self.read_instance.active_qa, self.read_instance.active_flags, \
                                  self.data_dtype, self.data_vars_to_read, \
                                  self.read_instance.metadata_dtype, self.read_instance.metadata_vars_to_read
                # read file
                file_data, time_indices, full_array_station_indices = read_netcdf_data(tuple_arguments)
                # place read data into big array as appropriate
                self.data_in_memory[data_label]['data'][full_array_station_indices[np.newaxis, :],
                                                        time_indices[:, np.newaxis]] = file_data

        # read in parallel
        elif self.read_type == 'parallel':
            # setup pool of N workers on N CPUs
            pool = multiprocessing.Pool(self.read_instance.n_cpus)
            # read netCDF files in parallel
            if not self.read_instance.reading_nonghost:
                tuple_arguments = [(file_name, self.read_instance.time_array, self.read_instance.station_references,
                                    species, process_type, self.read_instance.active_qa,
                                    self.read_instance.active_flags, self.data_dtype,
                                    self.data_vars_to_read, self.read_instance.metadata_dtype,
                                    self.read_instance.metadata_vars_to_read) for file_name in relevant_files]
                all_file_data = pool.map(read_netcdf_data, tuple_arguments)
            else:
                tuple_arguments = [
                    (file_name, self.read_instance.time_array, self.read_instance.station_references,
                     species, process_type) for
                    file_name in relevant_files]
                all_file_data = pool.map(read_netcdf_nonghost, tuple_arguments)

            pool.close()
            # wait for worker processes to terminate before continuing
            pool.join()

            # iterate through read file data and place data into data array as appropriate
            for file_data_ii, file_data in enumerate(all_file_data):
                try:
                    # some file_data might be none, in case the file did not exist
                    self.data_in_memory[data_label][file_data[2][:, np.newaxis], file_data[1][np.newaxis, :]] = \
                        file_data[0]
                except Exception as e:
                    continue
                if process_type == 'observations':
                    if not self.read_instance.reading_nonghost:
                        self.metadata_in_memory[file_data[2][:, np.newaxis],
                                                self.read_instance.metadata_inds_to_fill[file_data_ii]] = file_data[3]
                    else:
                        self.nonghost_metadata[file_data[2][:, np.newaxis]] = file_data[3]

    def get_valid_obs_files_in_date_range(self, selected_start_date, selected_end_date):
        """Define function that iterates through observational dictionary tree
        and returns a dictionary of available data in the selected date
        range"""

        # create dictionary to store available observational data
        self.available_observation_data = {}

        # check if start/end date are valid values, if not, return with no valid obs. files
        if (self.valid_date(selected_start_date)) & (self.valid_date(selected_end_date)):
            self.read_instance.date_range_has_changed = True
            self.read_instance.selected_start_date = int(selected_start_date)
            self.read_instance.selected_end_date = int(selected_end_date)
            self.read_instance.selected_start_date_firstdayofmonth = \
                int(str(self.read_instance.selected_start_date)[:6] + '01')
        else:
            return

        # check end date is > start date, if not, return with no valid obs. files
        if self.read_instance.selected_start_date >= self.read_instance.selected_end_date:
            return
        # check start date and end date are both within if valid date range (19000101 - 20500101),
        # if not, return with no valid obs. files
        if (self.read_instance.selected_start_date < 19000101) or (self.read_instance.selected_end_date < 19000101) or (
                self.read_instance.selected_start_date >= 20500101) or (self.read_instance.selected_end_date >= 20500101):
            return

        # iterate through networks
        for network in list(self.read_instance.all_observation_data.keys()):
            for resolution in list(self.read_instance.all_observation_data[network].keys()):
                for matrix in list(self.read_instance.all_observation_data[network][resolution].keys()):
                    for species in list(self.read_instance.all_observation_data[network][resolution][matrix].keys()):
                        # get all file yearmonths associated with species
                        species_file_yearmonths = self.read_instance.all_observation_data[network][resolution][matrix][species]
                        # get file yearmonths within date range
                        valid_species_files_yearmonths = [ym for ym in species_file_yearmonths if
                                                          (ym >= self.read_instance.selected_start_date_firstdayofmonth) & (
                                                                      ym < self.read_instance.selected_end_date)]
                        if len(valid_species_files_yearmonths) > 0:
                            # if network/res/matrix/species not in dictionary yet, add it
                            if network not in list(self.available_observation_data.keys()):
                                self.available_observation_data[network] = {}
                            if resolution not in list(self.available_observation_data[network].keys()):
                                self.available_observation_data[network][resolution] = {}
                            if matrix not in list(self.available_observation_data[network][resolution].keys()):
                                self.available_observation_data[network][resolution][matrix] = {}
                            self.available_observation_data[network][resolution][matrix][
                                species] = valid_species_files_yearmonths

    def get_valid_experiment_files_in_date_range(self):
        """Define function which gathers available experiment
        data for selected network/resolution/species.
        A dictionary is created storing available experiment-grid
        names associated with valid files in set date range.
        """

        # create dictionary to store available experiment information
        self.available_experiment_data = {}

        # get all different experiment names
        available_experiments = os.listdir('%s/%s' % (self.read_instance.exp_root, self.read_instance.ghost_version))

        # iterate through available experiments
        for experiment in available_experiments:

            # test first if interpolated directory exists before trying to get files from it
            # if it does not exit, continue
            if not os.path.exists(
                    '%s/%s/%s/%s/%s/%s' % (self.read_instance.exp_root, self.read_instance.ghost_version, experiment,
                                           self.read_instance.selected_resolution, self.read_instance.selected_species,
                                           self.read_instance.selected_network)):
                continue
            else:
                # get all experiment netCDF files by experiment/grid/selected
                # resolution/selected species/selected network
                network_files = os.listdir(
                    '%s/%s/%s/%s/%s/%s' % (self.read_instance.exp_root, self.read_instance.ghost_version,
                                           experiment, self.read_instance.selected_resolution,
                                           self.read_instance.selected_species, self.read_instance.selected_network))
                # get start YYYYMM yearmonths of data files
                network_files_yearmonths = [int(f.split('_')[-1][:6] + '01') for f in network_files]
                # limit data files to just those within date range
                valid_network_files_yearmonths = \
                    [ym for ym in network_files_yearmonths if (ym >= self.read_instance.selected_start_date_firstdayofmonth) &
                     (ym < self.read_instance.selected_end_date)]

                # if have some valid data files for experiment, add experiment key
                # (with associated yearmonths) to dictionary
                if len(valid_network_files_yearmonths) > 0:
                    self.available_experiment_data['%s' % (experiment)] = valid_network_files_yearmonths

        # get list of available experiment-grid names
        if not self.read_instance.offline:
            self.read_instance.experiments_menu['checkboxes']['labels'] = np.array(
                sorted(list(self.available_experiment_data.keys())))
            self.read_instance.experiments_menu['checkboxes']['map_vars'] = copy.deepcopy(
                self.read_instance.experiments_menu['checkboxes']['labels'])

    def update_plotting_parameters(self):
        """Function that updates plotting parameters (colour
        and zorder) for each selected data array
        """

        # assign a colour/zorder to all selected data arrays
        # define observations colour to be 'black'
        self.plotting_params['observations']['colour'] = 'black'
        # define zorder of observations to be 5
        self.plotting_params['observations']['zorder'] = 5

        # generate a list of RGB tuples for number of experiments there are
        sns.reset_orig()
        clrs = sns.color_palette('husl', n_colors=len(list(self.data_in_memory.keys()))-1)

        # iterate through sorted experiment names, assigning each experiment a new RGB colour tuple, and zorder
        experiment_ind = 1
        for experiment in sorted(list(self.data_in_memory.keys())):
            if experiment != 'observations':
                # define colour for experiment
                self.plotting_params[experiment]['colour'] = clrs[experiment_ind-1]
                # define zorder for experiment (obs zorder + experiment_ind)
                self.plotting_params[experiment]['zorder'] = \
                    self.plotting_params['observations']['zorder'] + experiment_ind
                # update count of experiments
                experiment_ind += 1
