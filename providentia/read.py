from .read_aux import get_yearmonths_to_read, init_shared_vars_read_netcdf_data, read_netcdf_data, read_netcdf_nonghost
from providentia import aux

import sys
import os
import copy
import ctypes
import datetime
import multiprocessing
import time

import numpy as np
import pandas as pd
import seaborn as sns
from netCDF4 import Dataset


class DataReader:
    """Class that reads observational/experiment data into memory."""

    def __init__(self, read_instance, read_type='parallel'):
        self.read_instance = read_instance
        self.read_type = read_type

    def reset_data_in_memory(self):
        self.data_in_memory = {}
        self.plotting_params = {}
        
    def read_setup(self, resolution, start_date, end_date, network, species, matrix, reset=False):
        """Setup key variables for new read of observational/experiment
        data a time array and create arrays of unique station
        references/longitudes/latitudes.

        :param resolution: resolution (e.g. "hourly")
        :type resolution: str
        :param start_date: start date (e.g. "20201101")
        :type start_date: str
        :param end_date: end date (e.g. "20201231")
        :type end_date: str
        :param network: network (e.g. "EBAS")
        :type network: str
        :param species: species (e.g. "sconco3")
        :type species: str
        :param matrix: matrix (e.g. "gas")
        :type matrix: str
        """

        # series of actions not applicable when --offline
        if not self.read_instance.offline:
            self.read_instance.reading_nonghost = aux.check_for_ghost(network)

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
        # get time array as integer timestamps
        self.read_instance.timestamp_array = self.read_instance.time_array.asi8

        #get relevant observational files
        if not self.read_instance.reading_nonghost:
            # get all relevant observational files
            file_root = '%s/%s/%s/%s/%s/%s_' % (self.read_instance.obs_root, network,
                                                self.read_instance.ghost_version, resolution,
                                                species, species)
        else:
            # get files from nonghost path
            file_root = '%s/%s/%s/%s/%s/%s_' % (self.read_instance.nonghost_root, network[1:].lower(),
                                                matrix, resolution, species, species)

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

        #get station references, longitudes and latitudes
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
                                          'season_code']
            elif (resolution == '3hourly') or \
                    (resolution == '6hourly') or (resolution == '3hourly_instantaneous') or \
                    (resolution == '6hourly_instantaneous'):
                 self.data_vars_to_read = [species, 'daily_native_representativity_percent',
                                           'monthly_native_representativity_percent',
                                           'annual_native_representativity_percent',
                                           'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                           'annual_native_max_gap_percent', 'day_night_code', 'weekday_weekend_code',
                                           'season_code']
            elif resolution == 'daily':
                self.data_vars_to_read = [species, 'daily_native_representativity_percent',
                                          'monthly_native_representativity_percent',
                                          'annual_native_representativity_percent',
                                          'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                          'annual_native_max_gap_percent', 'weekday_weekend_code', 'season_code']
            elif resolution == 'monthly':
                self.data_vars_to_read = [species, 'monthly_native_representativity_percent',
                                          'annual_native_representativity_percent', 'monthly_native_max_gap_percent',
                                          'annual_native_max_gap_percent', 'season_code']
        else:
            self.data_vars_to_read = [species]

        #need to reset metadata and data arrays?
        if reset:
                    
            #create new data AND metadata arrays 

            #data
            self.data_in_memory[data_label] = np.full((len(self.read_instance.data_labels),
                                                       len(self.data_vars_to_read),
                                                       len(self.read_instance.station_references),
                                                       len(self.read_instance.time_array)),
                                                       np.NaN, dtype=np.float32)

            #metadata
            if self.read_instance.reading_nonghost:
                tmp_ncdf = Dataset(relevant_files[0])
                meta_dtype = [('station_name', np.object), ('latitude', np.float32),
                            ('longitude', np.float32), ('altitude', np.float32)]
                if "station_code" in tmp_ncdf.variables:
                    meta_dtype.append(('station_reference', np.object))
                if "station_type" in tmp_ncdf.variables:
                    meta_dtype.append(('station_type', np.object))
                if "station_area" in tmp_ncdf.variables:
                    meta_dtype.append(('station_area', np.object))
                tmp_ncdf.close()
            else:
                meta_dtype = self.read_instance.metadata_dtype
            self.metadata_in_memory = np.full((len(self.read_instance.station_references),
                                            len(self.read_instance.relevant_yearmonths)),
                                            np.NaN, dtype=meta_dtype)
                        
            #create new plotting param dictionaries per data label
            # get experiment specific grid edges for exp, from first relevant file
            for data_label in self.read_instance.data_labels:
                self.plotting_params[data_label] = {}
                if data_label != 'observations':
                    exp_nc_root = Dataset(relevant_files[0])
                    self.plotting_params[data_label]['grid_edge_longitude'] = \
                        exp_nc_root['grid_edge_longitude'][:]
                    self.plotting_params[data_label]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
                    exp_nc_root.close()


    def read_data(self, data_labels, start_date, end_date,
                  network, resolution, species, matrix):
        """Function that handles reading of observational/experiment data.

        :param data_labels: label of data array to read (e.g "observations" or expid)
        :type data_labels: list 
        :param start_date: start date (e.g. "20201101")
        :type start_date: str
        :param end_date: end date (e.g. "20201231")
        :type end_date: str
        :param resolution: resolution (e.g. "hourly")
        :type resolution: str
        :param network: network (e.g. "EBAS")
        :type network: str
        :param species: species (e.g. "sconco3")
        :type species: str
        :param matrix: matrix (e.g. "gas")
        :type matrix: str
        """

        print('READ DATA START')

        #iterate through data labels, setting up variables for read
        relevant_files = []
        process_types = []
        for data_label in data_labels:

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
            yearmonths_to_read = get_yearmonths_to_read(relevant_file_start_dates, start_date, end_date, resolution)
            
            files_to_read = [file_root+str(yyyymm)[:6]+'.nc' for yyyymm in yearmonths_to_read]
            if not os.path.exists(files_to_read[0]):
                relevant_files = relevant_files + sorted([file_root + str(yyyymm)[:8] + '.nc' for yyyymm
                                                          in self.read_instance.relevant_yearmonths])

        #if active selected qa == default qa, no need to screen by qa, so set selected qa to None
        if self.read_instance.active_qa == self.read_instance.default_qa:
            qa_to_filter = []
        else:
            qa_to_filter = self.read_instance.active_qa  

        #create arrays to share across processes (for parallel multiprocessing use)
        #this only works for numerical dtypes, i.e. not strings
        file_data_shared_shape = (len(data_labels), len(self.data_vars_to_read), len(self.read_instance.station_references), len(self.read_instance.time_array))
        file_data_shared = multiprocessing.RawArray(ctypes.c_float, file_data_shared_shape[0] * file_data_shared_shape[1] * file_data_shared_shape[2] * file_data_shared_shape[3])  
        timestamp_array_shared = multiprocessing.RawArray(ctypes.c_int64, len(self.read_instance.timestamp_array))
        qa_shared = multiprocessing.RawArray(ctypes.c_uint8, len(qa_to_filter))
        flags_shared = multiprocessing.RawArray(ctypes.c_uint8, len(self.read_instance.active_flags))
        # Wrap file_data_shared as an numpy array so we can easily manipulates its data.
        file_data_shared_np = np.frombuffer(file_data_shared, dtype=np.float32).reshape(file_data_shared_shape)
        #fill arrays
        data_label_indices = [self.read_instance.data_labels.index(data_label) for data_label in data_labels]
        np.copyto(file_data_shared_np, self.data_in_memory[data_label_indices, :, :, :])
        timestamp_array_shared[:] = self.read_instance.timestamp_array
        qa_shared[:] = qa_to_filter
        flags_shared[:] = self.read_instance.active_flags

        # iterate and read species data in all relevant netCDF files (either in serial/parallel)
        s = time.time()
        # read serially
        if self.read_type == 'serial':
            
            # iterate through relevant netCDF files
            for relevant_file_ii, relevant_file in enumerate(relevant_files):
                # create argument tuple of function
                tuple_arguments = relevant_file, self.read_instance.station_references, species, process_type, \
                                  self.data_vars_to_read, self.read_instance.metadata_dtype, self.read_instance.metadata_vars_to_read
                # read file
                file_data, time_indices, full_array_station_indices, file_metadata = read_netcdf_data(tuple_arguments)
                # place read data into big array as appropriate
                self.data_in_memory[data_label_indices[:, np.newaxis, np.newaxis], :
                                    full_array_station_indices[:, np.newaxis],
                                    time_indices[np.newaxis, :]] = file_data
                #place metadata
                if process_type == 'observations':
                    self.metadata_in_memory[full_array_station_indices[:, np.newaxis],
                                           self.read_instance.metadata_inds_to_fill[relevant_file_ii]] = file_metadata

        elif self.read_type == 'parallel':

            # setup pool of N workers on N CPUs
            pool = multiprocessing.Pool(self.read_instance.n_cpus, initializer=init_shared_vars_read_netcdf_data, initargs=(file_data_shared, file_data_shared_shape, timestamp_array_shared, qa_shared, flags_shared))
            # read netCDF files in parallel
            if not self.read_instance.reading_nonghost:
                tuple_arguments = [(file_name, self.read_instance.station_references, species, process_type, 
                                    self.data_vars_to_read, self.read_instance.metadata_dtype, self.read_instance.metadata_vars_to_read) for file_name in relevant_files]
                print('POOLING', time.time() - s)
                returned_data = pool.map(read_netcdf_data, tuple_arguments)
            else:
                tuple_arguments = [
                    (file_name, self.read_instance.station_references,
                     species, process_type) for
                    file_name in relevant_files]
                returned_data = pool.map(read_netcdf_nonghost, tuple_arguments)

            pool.close()
            # wait for worker processes to terminate before continuing
            pool.join()
            
            print('READY TO JOIN', time.time() - s)
            # iterate through read file data and place metadata into full array as appropriate
            for returned_data_ii, returned_data_per_month in enumerate(returned_data):
                if process_type == 'observations':
                    self.metadata_in_memory[returned_data_per_month[0][:, np.newaxis],
                                            self.read_instance.metadata_inds_to_fill[returned_data_ii]] = returned_data_per_month[1]

            print('METADATA PLACED', time.time() - s)

        #overwrite data in memory
        self.data_in_memory[data_label] = file_data_shared_np

        # check if datasets consist of arrays full of -9999.0 or nan values or if they are empty
        if (self.data_in_memory[data_label].size == 0 or
            np.isin(self.data_in_memory[data_label].flatten(), [-9999.0, np.nan]).all()):

            if self.data_in_memory[data_label].size == 0:
                print('Error: The observation or experiment datasets are empty.')
            
            elif np.isin(self.data_in_memory[data_label].flatten(), [-9999.0, np.nan]).all():
                print('Error: The observation or experiment datasets are void.')

            print('Check if the data from the observations was downloaded correctly and')
            print('if the experiments were interpolated at the stations of the network of interest.')
            sys.exit()

        print('READ DATA END', time.time() - s)

    def get_valid_obs_files_in_date_range(self, start_date, end_date):
        """Define function that iterates through observational dictionary tree
        and returns a dictionary of available data in the selected date
        range

        :param start_date: start date (e.g. "20201101")
        :type start_date: str
        :param end_date: end date (e.g. "20201101")
        :type end_date: str
        """

        # create dictionary to store available observational data
        self.available_observation_data = {}

        # check if start/end date are valid values, if not, return with no valid obs. files
        if (not aux.valid_date(start_date)) or (not aux.valid_date(end_date)):
            return False

        # check end date is > start date, if not, return with no valid obs. files
        if start_date >= end_date:
            return False

        # check start date and end date are both within if valid date range (19000101 - 20500101),
        # if not, return with no valid obs. files
        if (int(start_date) < 19000101) or (int(end_date) < 19000101) or (int(start_date) >= 20500101) or (int(end_date) >= 20500101):
            return False

        #get start date at first of month
        start_date_firstdayofmonth = int(start_date[:6] + '01')

        # iterate through networks
        for network in list(self.read_instance.all_observation_data.keys()):
            for resolution in list(self.read_instance.all_observation_data[network].keys()):
                for matrix in list(self.read_instance.all_observation_data[network][resolution].keys()):
                    for species in list(self.read_instance.all_observation_data[network][resolution][matrix].keys()):
                        # get all file yearmonths associated with species
                        species_file_yearmonths = self.read_instance.all_observation_data[network][resolution][matrix][species]
                        # get file yearmonths within date range
                        valid_species_files_yearmonths = [ym for ym in species_file_yearmonths if
                                                          (ym >= start_date_firstdayofmonth) & (ym < int(end_date))]
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

        return True

    def get_valid_experiment_files_in_date_range(self, start_date, end_date, resolution, network, species):
        """Define function which gathers available experiment
        data for selected network/resolution/species.
        A dictionary is created storing available experiment
        names associated with valid files in set date range.

        :param start_date: start date (e.g. "20201101")
        :type start_date: str
        :param end_date: end date (e.g. "20201231")
        :type end_date: str
        :param resolution: resolution (e.g. "hourly")
        :type resolution: str
        :param network: network (e.g. "EBAS")
        :type network: str
        :param species: species (e.g. "sconco3")
        :type species: str
        """

        # create dictionary to store available experiment information
        self.available_experiment_data = {}

        # get all different experiment names
        available_experiments = os.listdir('%s/%s' % (self.read_instance.exp_root, self.read_instance.ghost_version))

        #get start date at first of month
        start_date_firstdayofmonth = int(start_date[:6] + '01')

        # iterate through available experiments
        for experiment in available_experiments:

            # test first if interpolated directory exists before trying to get files from it
            # if it does not exit, continue
            if not os.path.exists(
                    '%s/%s/%s/%s/%s/%s' % (self.read_instance.exp_root, self.read_instance.ghost_version, 
                                           experiment, resolution, species, network)):
                continue
            else:
                # get all experiment netCDF files by experiment/grid/selected
                # resolution/selected species/selected network
                network_files = os.listdir(
                    '%s/%s/%s/%s/%s/%s' % (self.read_instance.exp_root, self.read_instance.ghost_version,
                                           experiment, resolution, species, network))

                # get start YYYYMM yearmonths of data files
                network_files_yearmonths = [int(f.split('_')[-1][:6] + '01') for f in network_files]
                # limit data files to just those within date range
                valid_network_files_yearmonths = \
                    [ym for ym in network_files_yearmonths if (ym >= start_date_firstdayofmonth) & (ym < int(end_date))]

                # if have some valid data files for experiment, add experiment key
                # (with associated yearmonths) to dictionary
                if len(valid_network_files_yearmonths) > 0:
                    self.available_experiment_data['%s' % (experiment)] = valid_network_files_yearmonths

        # get list of available experiment names
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
        # define observations colour
        self.plotting_params['observations']['colour'] = self.read_instance.plot_characteristics_templates['general']['obs_markerfacecolor']
        # define zorder 
        self.plotting_params['observations']['zorder'] = self.read_instance.plot_characteristics_templates['general']['obs_zorder']

        # generate a list of RGB tuples for number of experiments there are
        sns.reset_orig()
        clrs = sns.color_palette(self.read_instance.plot_characteristics_templates['general']['legend_color_palette'], n_colors=len(list(self.data_in_memory.keys()))-1)

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
