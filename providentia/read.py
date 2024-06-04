""" Class that reads observational/experiment data into memory """

import copy
import ctypes
import datetime
import multiprocessing
import os
import sys
import time

from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset, chartostring
import numpy as np
import pandas as pd

from .plot_aux import update_plotting_parameters
from .read_aux import (check_for_ghost, get_default_qa, get_frequency_code, get_yearmonths_to_read, 
                       init_shared_vars_read_netcdf_data, read_netcdf_data, read_netcdf_metadata)
from .spatial_colocation import (resolve_duplicate_spatial_colocation_matches, spatial_colocation_ghost, 
                                 spatial_colocation_nonghost)
from .warnings import show_message

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class DataReader:

    def __init__(self, read_instance):
        self.read_instance = read_instance
        
    def read_setup(self, operations, experiments_to_remove=None, experiments_to_read=None):
        """ Setup structures for read of new observational/experiment data and then perform read.

            :param operations: list of instructions on how to adjust data structures
            :type operation: list
            :param experiments_to_remove: list of experiments to remove from arrays
            :type experiments_to_remove: list
            :param experiments_to_read: list of experiments to add to arrays
            :type experiments_to_read: list
        """

        # changing time dimension ?
        if ('reset' in operations) or ('read_left' in operations) or ('read_right' in operations) or ('cut_left' in operations) or ('cut_right' in operations):

            # turn off read from configuration
            if 'reset' not in operations:
                self.read_instance.from_conf = False

            # determine if reading GHOST or non-GHOST
            self.read_instance.reading_ghost = check_for_ghost(self.read_instance.network[0])

            # get active frequency code
            self.read_instance.active_frequency_code = get_frequency_code(self.read_instance.resolution)

            # get time array
            self.read_instance.time_array = pd.date_range(start=datetime.datetime(int(str(self.read_instance.start_date)[:4]),
                                                                                  int(str(self.read_instance.start_date)[4:6]),
                                                                                  int(str(self.read_instance.start_date)[6:8])),
                                                          end=datetime.datetime(int(str(self.read_instance.end_date)[:4]),
                                                                                int(str(self.read_instance.end_date)[4:6]),
                                                                                int(str(self.read_instance.end_date)[6:8])),
                                                          freq=self.read_instance.active_frequency_code)[:-1]

            # show warning when the data consists only of less than 2 timesteps
            if len(self.read_instance.time_array) < 2:
                self.read_instance.invalid_read = True
                msg = 'Extend the time range or enhance the resolution (e.g. from monthly to daily) to create plots. '
                msg += 'Plots will only be created when period is longer than 2 timesteps.'
                show_message(self.read_instance, msg)
                if (self.read_instance.from_conf) and (not self.read_instance.offline) and (not self.read_instance.interactive):
                    sys.exit('Error: Providentia will not be launched.')
                elif (self.read_instance.offline):
                    sys.exit('Error: Offline report will not be created.')
                elif (self.read_instance.interactive):
                    sys.exit('Error: Data cannot be read')
                else:
                    self.read_instance.first_read = True
                    return
            else:
                # get list of extra networkspecies to read, used for filtering data
                # read only networkspecies not present in current networkspecies to read
                self.read_instance.filter_networkspecies = []
                if self.read_instance.filter_species:
                    for networkspeci in self.read_instance.filter_species:
                        if networkspeci not in self.read_instance.networkspecies:
                            self.read_instance.filter_networkspecies.append(networkspeci)
                    
                    # do not update data bounds for current networkspecies
                    filter_species = copy.deepcopy(self.read_instance.filter_species)
                    for networkspeci in filter_species:
                        if networkspeci in self.read_instance.networkspecies:
                            msg = 'The current network-species ({}) cannot be selected as a filter species. '.format(networkspeci)
                            msg += 'If you want to change its data bounds, use the lower and upper bounds parameters.'
                            show_message(self.read_instance, msg)
                            del self.read_instance.filter_species[networkspeci]

                # get yearmonths in data range (incomplete months are removed for monthly resolution)
                self.read_instance.yearmonths = list(np.unique(['{}0{}'.format(dti.year,dti.month) if len(str(dti.month)) == 1 else '{}{}'.format(dti.year,dti.month) \
                                                                for dti in self.read_instance.time_array]))
                
                # get time array as integer timestamps
                self.read_instance.timestamp_array = self.read_instance.time_array.asi8

                # get N indices per yearmonth
                self.read_instance.N_inds_per_yearmonth = np.array([np.count_nonzero(np.all(
                    [self.read_instance.time_array >= datetime.datetime.strptime(start_yyyymm+'01','%Y%m%d'),
                    self.read_instance.time_array < datetime.datetime.strptime(self.read_instance.yearmonths[month_ii + 1]+'01','%Y%m%d')], axis=0)) 
                    if month_ii != (len(self.read_instance.yearmonths) - 1) else np.count_nonzero(
                    self.read_instance.time_array >= datetime.datetime.strptime(start_yyyymm+'01','%Y%m%d')) 
                    for month_ii, start_yyyymm in enumerate(self.read_instance.yearmonths)])

                # get unique basic metadata across networkspecies
                # for this step include filter networkspecies
                self.read_basic_metadata()

                # iterate through station_references per networkspecies
                # if have 0 valid stations then drop  
                for networkspeci, stn_refs in self.read_instance.station_references.items():
                    if len(stn_refs) == 0:
                        if networkspeci in self.read_instance.networkspecies:
                            self.read_instance.networkspecies.remove(networkspeci)
                            print('Warning: There is no available observational data for the network|species: {}. Dropping.'.format(networkspeci))
                        elif networkspeci in self.read_instance.filter_networkspecies:
                            self.read_instance.filter_networkspecies.remove(networkspeci)
                            print('Warning: There is no available observational data for the filter network|species: {}. Dropping.'.format(networkspeci))

                # if have zero networkspecies left, then return with invalid_read
                if len(self.read_instance.networkspecies) == 0:
                    self.read_instance.invalid_read = True
                    return

                # set invalid_read to be False if have data to read
                self.read_instance.invalid_read = False

        # need to reset all data structures 
        if 'reset' in operations:  

            # uninitialise filter object
            if (not self.read_instance.offline) and (not self.read_instance.interactive):
                self.read_instance.mpl_canvas.filter_data = None

            # data
            self.read_instance.data_in_memory = {networkspeci: 
                                                 np.full((len(self.read_instance.data_labels),
                                                          len(self.read_instance.station_references[networkspeci]),
                                                          len(self.read_instance.time_array)),
                                                          np.NaN, dtype=np.float32) for networkspeci in self.read_instance.networkspecies}

            # filter data (if active)
            if self.read_instance.filter_species:
                self.read_instance.filter_data_in_memory = {networkspeci: 
                                                            np.full((len(self.read_instance.station_references[networkspeci]),
                                                                     len(self.read_instance.time_array)),
                                                                     np.NaN, dtype=np.float32) for networkspeci in self.read_instance.filter_networkspecies}
            else:
                self.read_instance.filter_data_in_memory = {}

            # GHOST data --> data variables which change per measurement (for filtering)
            if self.read_instance.reading_ghost:

                # set ghost data variables to read (dependent on data resolution)
                if (self.read_instance.resolution == 'hourly') or (self.read_instance.resolution == 'hourly_instantaneous'):
                    self.read_instance.ghost_data_vars_to_read = ['hourly_native_representativity_percent',
                                                                  'daily_native_representativity_percent',
                                                                  'monthly_native_representativity_percent',
                                                                  'annual_native_representativity_percent', 
                                                                  'day_night_code', 'weekday_weekend_code',
                                                                  'season_code']
                elif (self.read_instance.resolution == '3hourly') or \
                        (self.read_instance.resolution == '6hourly') or (self.read_instance.resolution == '3hourly_instantaneous') or \
                        (self.read_instance.resolution == '6hourly_instantaneous'):
                    self.read_instance.ghost_data_vars_to_read = ['daily_native_representativity_percent',
                                                                  'monthly_native_representativity_percent',
                                                                  'annual_native_representativity_percent',
                                                                  'day_night_code', 'weekday_weekend_code',
                                                                  'season_code']
                elif self.read_instance.resolution == 'daily':
                    self.read_instance.ghost_data_vars_to_read = ['daily_native_representativity_percent',
                                                                  'monthly_native_representativity_percent',
                                                                  'annual_native_representativity_percent',
                                                                  'weekday_weekend_code', 'season_code']
                elif self.read_instance.resolution == 'monthly':
                    self.read_instance.ghost_data_vars_to_read = ['monthly_native_representativity_percent',
                                                                  'annual_native_representativity_percent', 
                                                                  'season_code']

                self.read_instance.ghost_data_in_memory = {networkspeci:
                                        np.full((len(self.read_instance.ghost_data_vars_to_read),
                                                 len(self.read_instance.station_references[networkspeci]),
                                                 len(self.read_instance.time_array)),
                                                 np.NaN, dtype=np.float32) for networkspeci in self.read_instance.networkspecies} 
            else:
                self.read_instance.ghost_data_in_memory = {}
                self.read_instance.ghost_data_vars_to_read = []

            # metadata 
            # non-GHOST
            if not self.read_instance.reading_ghost:
                self.read_instance.metadata_dtype = [(nonghost_var, self.read_instance.standard_metadata[nonghost_var]['data_type'])
                                                      for nonghost_var in self.read_instance.nonghost_metadata_vars_to_read]
                self.read_instance.metadata_vars_to_read = self.read_instance.nonghost_metadata_vars_to_read
            # GHOST
            else:
                self.read_instance.metadata_dtype = self.read_instance.ghost_metadata_dtype
                self.read_instance.metadata_vars_to_read = self.read_instance.ghost_metadata_vars_to_read

            # show warning when there is a non-defined field if launching from a config file
            if hasattr(self.read_instance, "non_default_fields_per_section"):
                invalid_args = {field_name:fields-set(self.read_instance.metadata_vars_to_read) 
                                for field_name, fields in self.read_instance.non_default_fields_per_section.items()}
                invalid_var = [f"""{i} ('{"', '".join(j)}')""" for i,j in invalid_args.items() if j]
                if invalid_var:
                    msg = f"Invalid field(s) in configuration file {self.read_instance.config.split('/')[-1]}. "
                    msg += f"Section(s) and Field(s): {', '.join(invalid_var)}."
                    show_message(self.read_instance, msg)

            self.read_instance.metadata_in_memory = {networkspeci: 
                                                     np.full((len(self.read_instance.station_references['{}'.format(networkspeci)]),
                                                              len(self.read_instance.yearmonths)),
                                                              np.NaN, dtype=self.read_instance.metadata_dtype) 
                                                              for networkspeci in self.read_instance.networkspecies}

            # get list of yearmonths to read
            yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.start_date,
                                                        self.read_instance.end_date, self.read_instance.resolution)

            # read data 
            self.read_data(yearmonths_to_read, self.read_instance.data_labels)

            # update measurement units for all species (take standard units for each speci from parameter dictionary)
            # non-GHOST
            if not self.read_instance.reading_ghost:
                # import unit converter
                sys.path.insert(1, os.path.join(CURRENT_PATH, 'dependencies/unit-converter'))
                import unit_converter

                # convert non-GHOST units to standard format
                nonghost_standard_units = {}
                for speci in self.read_instance.nonghost_units.keys():
                    input_units = self.read_instance.nonghost_units[speci]
                    if input_units != '-':
                        output_units = copy.deepcopy(input_units)
                        conv_obj = unit_converter.convert_units(input_units, output_units, 1)
                        nonghost_standard_units[speci] = conv_obj.output_standard_units
                    else:
                        nonghost_standard_units[speci] = 'unitless'
                self.read_instance.measurement_units = {speci:nonghost_standard_units[speci] 
                                                        for speci in self.read_instance.species}
            # GHOST
            else:
                self.read_instance.measurement_units = {speci:self.read_instance.parameter_dictionary[speci]['standard_units'] 
                                                        for speci in self.read_instance.species}

            # reset plotting params
            self.read_instance.plotting_params = {}

            # iterate through data labels
            for data_label in self.read_instance.data_labels:
                self.read_instance.plotting_params[data_label] = {}

                # get experiment specific grid edges for exp, from first relevant file
                if data_label != self.read_instance.observations_data_label:
                    # iterate through networkspecies until find one which has valid files to read
                    for valid_networkspeci in self.read_instance.networkspecies:
                        if data_label in self.files_to_read[valid_networkspeci]:
                            exp_nc_root = Dataset(self.files_to_read[valid_networkspeci][data_label][0])
                            self.read_instance.plotting_params[data_label]['grid_edge_longitude'] = exp_nc_root['grid_edge_longitude'][:]
                            self.read_instance.plotting_params[data_label]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
                            exp_nc_root.close()
                            break
            
            # update plotting parameters colours and zorder
            update_plotting_parameters(self.read_instance) 

        # need to read on left / read on right / cut on left / cut on right (for dashboard)
        if ('read_left' in operations) or ('read_right' in operations) or ('cut_left' in operations) or ('cut_right' in operations):  

            # if station references array has changed then as cutting / appending to
            # need to rearrange existing metadata/data arrays accordingly
            if not np.array_equal(self.read_instance.previous_station_references[self.read_instance.networkspecies[0]], \
                                    self.read_instance.station_references[self.read_instance.networkspecies[0]]):

                # get indices of stations in previous station references array in current station references array
                old_station_inds = np.where(np.in1d(self.read_instance.previous_station_references[self.read_instance.networkspecies[0]],
                                                    self.read_instance.station_references[self.read_instance.networkspecies[0]]))[0]
                                                    
                # get indices of stations in current station references array
                # that were in previous station references array
                new_station_inds = np.where(np.in1d(self.read_instance.station_references[self.read_instance.networkspecies[0]],
                                                    self.read_instance.previous_station_references[self.read_instance.networkspecies[0]]))[0]

                # rearrange metadata station dimension
                new_metadata_in_memory = np.full((len(self.read_instance.station_references[self.read_instance.networkspecies[0]]),
                                                  len(self.read_instance.previous_yearmonths)),
                                                  np.NaN, dtype=self.read_instance.metadata_dtype)
                new_metadata_in_memory[new_station_inds, :] = self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]][old_station_inds, :]
                self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]] = new_metadata_in_memory

                # rearrage data array station dimension
                new_data_in_memory = np.full((len(self.read_instance.previous_data_labels),
                                              len(self.read_instance.station_references[self.read_instance.networkspecies[0]]),
                                              len(self.read_instance.previous_time_array)),
                                              np.NaN, dtype=np.float32)
                # put the old data into new array in the correct positions
                new_data_in_memory[:, new_station_inds, :] = self.read_instance.data_in_memory[self.read_instance.networkspecies[0]][:, old_station_inds, :]
                # overwrite data array with reshaped version
                self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = new_data_in_memory

                # rearrage filter data array station dimension
                # iterate through all filter networkspecies
                if self.read_instance.filter_species:    
                    for filter_networkspeci in self.read_instance.filter_networkspecies:            
                        new_filter_data_in_memory = np.full((len(self.read_instance.station_references[filter_networkspeci]),
                                                             len(self.read_instance.previous_time_array)),
                                                             np.NaN, dtype=np.float32)
                        # put the old data into new array in the correct positions
                        new_filter_data_in_memory[new_station_inds, :] = self.read_instance.filter_data_in_memory[filter_networkspeci][old_station_inds, :]
                        # overwrite data array with reshaped version
                        self.read_instance.filter_data_in_memory[filter_networkspeci] = new_filter_data_in_memory

                # rearrage ghost data array station dimension
                if self.read_instance.reading_ghost:
                    new_ghost_data_in_memory = np.full((len(self.read_instance.ghost_data_vars_to_read),
                                                        len(self.read_instance.station_references[self.read_instance.networkspecies[0]]),
                                                        len(self.read_instance.previous_time_array)),
                                                        np.NaN, dtype=np.float32)
                    # put the old ghost data into new array in the correct positions
                    new_ghost_data_in_memory[:, new_station_inds, :] = self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]][:, old_station_inds, :]
                    # overwrite ghost data array with reshaped version
                    self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]] = new_ghost_data_in_memory

            # need to cut on left / cut on right
            if ('cut_left' in operations) or ('cut_right' in operations):

                # set default edge limits as current edges
                data_left_edge_ind = 0
                data_right_edge_ind = len(self.read_instance.previous_time_array)

                metadata_left_edge_ind = 0
                metadata_right_edge_ind = len(self.read_instance.previous_yearmonths)

                # need to cut on left data edge?
                if 'cut_left' in operations:
                    data_left_edge_ind = np.where(self.read_instance.previous_time_array == self.read_instance.time_array[0])[0][0]
                    str_first_yearmonth = str(self.read_instance.yearmonths[0])
                    str_previous_first_yearmonth = str(self.read_instance.previous_yearmonths[0])
                    monthly_relative_delta = relativedelta(
                        datetime.datetime(int(str_first_yearmonth[:4]), int(str_first_yearmonth[4:6]),
                                        1, 0, 0), datetime.datetime(int(str_previous_first_yearmonth[:4]),
                                                                    int(str_previous_first_yearmonth[4:6]),
                                                                    1, 0, 0))
                    metadata_left_edge_ind = (monthly_relative_delta.years * 12) + monthly_relative_delta.months

                # need to cut on right data edge?
                if 'cut_right' in operations:
                    data_right_edge_ind = np.where(self.read_instance.previous_time_array == self.read_instance.time_array[-1])[0][0] + 1
                    str_last_yearmonth = str(self.read_instance.yearmonths[-1])
                    str_previous_last_yearmonth = str(self.read_instance.previous_yearmonths[-1])
                    monthly_relative_delta = relativedelta(
                        datetime.datetime(int(str_previous_last_yearmonth[:4]),
                                        int(str_previous_last_yearmonth[4:6]),
                                        1, 0, 0), datetime.datetime(int(str_last_yearmonth[:4]),
                                                                    int(str_last_yearmonth[4:6]), 1, 0, 0))
                    metadata_right_edge_ind = \
                        metadata_right_edge_ind - ((monthly_relative_delta.years * 12) + monthly_relative_delta.months)

                # cut edges of metadata array
                if metadata_left_edge_ind == metadata_right_edge_ind:
                    self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]] = \
                        self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]][:, [metadata_left_edge_ind]]
                else:
                    self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]] = \
                        self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]][:, metadata_left_edge_ind:metadata_right_edge_ind]

                # cut edges of data array
                self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = \
                    self.read_instance.data_in_memory[self.read_instance.networkspecies[0]][:, :, data_left_edge_ind:data_right_edge_ind]

                # cut edges of filter data array
                if self.read_instance.filter_species:  
                    self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]] = \
                        self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]][:, data_left_edge_ind:data_right_edge_ind]

                # cut edges of ghost data array
                if self.read_instance.reading_ghost:
                    self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]] = \
                        self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]][:, :, data_left_edge_ind:data_right_edge_ind]

            # need to read on left / read on right
            if ('read_left' in operations) or ('read_right' in operations):

                # save list of all yearmonths to read on both edges
                all_yearmonths_to_read = []

                # need to read on left 
                if 'read_left' in operations:

                    # get n number of new elements on left edge
                    if self.read_instance.previous_time_array.size > 0:
                        n_new_left_data_inds = np.where(self.read_instance.time_array == self.read_instance.previous_time_array[0])[0][0]
                    else:
                        n_new_left_data_inds = len(self.read_instance.time_array)

                    # get list of yearmonths to read
                    yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.start_date,
                                                                self.read_instance.previous_start_date, self.read_instance.resolution)
                                                            
                    # check which yearmonths_to_read previously read
                    yearmonths_previously_read = np.isin(yearmonths_to_read, self.read_instance.previous_yearmonths)

                    # get yearmonths not currently accounted for
                    if isinstance(yearmonths_to_read, list):
                        yearmonths_to_read = np.asarray(yearmonths_to_read)
                    yearmonths_to_read = list(yearmonths_to_read[~yearmonths_previously_read])

                    # keep current month to update data array if number of months has not changed but number of days has
                    if len(yearmonths_to_read) == 0:
                        yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.start_date,
                                                                    self.read_instance.previous_start_date, self.read_instance.resolution)

                    # add space for new data on left edge of the metadata array
                    self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]] = \
                        np.concatenate((np.full((len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), len(yearmonths_to_read)), 
                            np.NaN, dtype=self.read_instance.metadata_dtype), self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]]), axis=1)

                    # insert space for new data on left edge of the data array
                    self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = \
                        np.concatenate((np.full((len(self.read_instance.previous_data_labels), len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_left_data_inds), 
                            np.NaN, dtype=np.float32), self.read_instance.data_in_memory[self.read_instance.networkspecies[0]]), axis=2)

                    # add space for new data on left edge of the filter data array
                    if self.read_instance.filter_species:  
                        self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((np.full((len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_left_data_inds), 
                                np.NaN, dtype=np.float32), self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]]), axis=1)

                    # insert space for new ghost data on left edge of the ghost data array
                    if self.read_instance.reading_ghost:
                        self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((np.full((len(self.read_instance.ghost_data_vars_to_read), len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_left_data_inds), 
                                np.NaN, dtype=np.float32), self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]]), axis=2)
                    
                    # add yearmonths_to_read to list for both edges
                    all_yearmonths_to_read.extend(yearmonths_to_read)

                # need to read on right
                if 'read_right' in operations:

                    # get n number of new elements on right edge
                    if self.read_instance.previous_time_array.size > 0:
                        n_new_right_data_inds = (len(self.read_instance.time_array) - 1) - \
                                                np.where(self.read_instance.time_array == self.read_instance.previous_time_array[-1])[0][0]
                    else:
                        n_new_right_data_inds = (len(self.read_instance.time_array))

                    # get list of yearmonths to read
                    yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.previous_end_date,
                                                                self.read_instance.end_date, self.read_instance.resolution)

                    # check which yearmonths_to_read previously read
                    yearmonths_previously_read = np.isin(yearmonths_to_read, self.read_instance.previous_yearmonths)

                    # get yearmonths not currently accounted for
                    if isinstance(yearmonths_to_read, list):
                        yearmonths_to_read = np.asarray(yearmonths_to_read)
                    yearmonths_to_read = list(yearmonths_to_read[~yearmonths_previously_read])

                    # keep current month to update data array if number of months has not changed but number of days has
                    if len(yearmonths_to_read) == 0:
                        yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.previous_end_date,
                                                                    self.read_instance.end_date, self.read_instance.resolution)

                    # add space for new data on right edge of the metadata array
                    self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]] = \
                        np.concatenate((self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]], 
                            np.full((len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), len(yearmonths_to_read)), 
                                np.NaN, dtype=self.read_instance.metadata_dtype)), axis=1)

                    # insert space for new data on right edge of the data array
                    self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = \
                        np.concatenate((self.read_instance.data_in_memory[self.read_instance.networkspecies[0]], 
                            np.full((len(self.read_instance.previous_data_labels), len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_right_data_inds), 
                                np.NaN, dtype=np.float32)), axis=2)

                    # insert space for new data on right edge of the filter data array
                    if self.read_instance.filter_species: 
                        self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]], 
                                np.full((len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_right_data_inds), 
                                    np.NaN, dtype=np.float32)), axis=1)

                    # insert space for new ghost data on right edge of the ghost data array
                    if self.read_instance.reading_ghost:
                        self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]], 
                                np.full((len(self.read_instance.ghost_data_vars_to_read), len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_right_data_inds), 
                                    np.NaN, dtype=np.float32)), axis=2)    

                    # add yearmonths_to_read to list for both edges
                    all_yearmonths_to_read.extend(yearmonths_to_read)

                # read data 
                self.read_data(all_yearmonths_to_read, self.read_instance.data_labels) 

        # need to remove experiment/s ?
        if 'remove_exp' in operations: 

            # get indices of experiments to remove
            experiments_to_remove_inds = [self.read_instance.previous_data_labels.index(experiment) for experiment in experiments_to_remove]
            
            # remove experiment data
            self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = \
                np.delete(self.read_instance.data_in_memory[self.read_instance.networkspecies[0]], experiments_to_remove_inds, axis=0)

            # remove plotting paramaters for experiments removed
            for experiment in experiments_to_remove:
                del self.read_instance.plotting_params[experiment]

        # need to read experiment/s ? 
        if 'read_exp' in operations: 

            # insert space for new experiments in data array
            for experiment in experiments_to_read:
                experiments_to_read_ind = self.read_instance.data_labels.index(experiment) 
                
                self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = \
                    np.insert(self.read_instance.data_in_memory[self.read_instance.networkspecies[0]], 
                                experiments_to_read_ind,
                                np.full((1, len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), len(self.read_instance.time_array)), 
                                np.NaN, dtype=np.float32),                     
                                axis=0)

            # get list of yearmonths to read
            yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.start_date,
                                                        self.read_instance.end_date, self.read_instance.resolution)

            # read data
            self.read_data(yearmonths_to_read, experiments_to_read)       

            # add experiment specific grid edges for exp to plotting params
            for experiment_to_read in experiments_to_read:
                self.read_instance.plotting_params[experiment_to_read] = {}
                exp_nc_root = Dataset(self.files_to_read[self.read_instance.networkspecies[0]][experiment_to_read][0])
                self.read_instance.plotting_params[experiment_to_read]['grid_edge_longitude'] = \
                    exp_nc_root['grid_edge_longitude'][:]
                self.read_instance.plotting_params[experiment_to_read]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
                exp_nc_root.close()

        # need to read or/and remove experiment/s ? 
        if ('read_exp' in operations) or ('remove_exp' in operations): 
            
            # update plotting parameters colours and zorder
            update_plotting_parameters(self.read_instance) 

        # print basic information species
        print('SELECTED SPECIES')
        print('- Main network-species', self.read_instance.networkspecies)
        if self.read_instance.filter_species:
            print('- Filter network-species', self.read_instance.filter_species)

    def read_basic_metadata(self):     
        """ Get basic unique metadata across networkspecies wanting to read
            The basic fields are: station_reference, station_name, longitude, latitude, measurement_altitude, 
            station_classification and area_classification

            If have multiple species, then spatially colocate across species 
            to get matching stations across stations.
        """

        # define dictionaries for storing basic metadata across all species to read
        self.read_instance.station_references = {}
        self.read_instance.station_names = {}
        self.read_instance.station_longitudes = {}
        self.read_instance.station_latitudes = {}
        self.read_instance.station_measurement_altitudes = {}
        self.read_instance.nonghost_units = {}

        # iterate through network, speci pairs
        for networkspeci in (self.read_instance.networkspecies + self.read_instance.filter_networkspecies):
        
            # get indivudual network and species strings
            network = networkspeci.split('|')[0]
            speci = networkspeci.split('|')[1]

            # get species matrix
            matrix = self.read_instance.parameter_dictionary[speci]['matrix']

            # get file root
            # GHOST
            if self.read_instance.reading_ghost:
                file_root = '%s/%s/%s/%s/%s/%s_' % (self.read_instance.ghost_root, network,
                                                    self.read_instance.ghost_version, self.read_instance.resolution,
                                                    speci, speci)
            # non-GHOST
            else:
                file_root = '%s/%s/%s/%s/%s_' % (self.read_instance.nonghost_root, network,
                                                 self.read_instance.resolution, speci, speci)

            # get relevant files
            relevant_files_before_filter = sorted([file_root+str(yyyymm)+'.nc' for yyyymm in self.read_instance.yearmonths])
            relevant_files = copy.deepcopy(relevant_files_before_filter)

            # drop files if they don't exist
            for file in relevant_files_before_filter:
                if not os.path.exists(file):
                    relevant_files.remove(file)

            # if have 0 files to read for networkspeci, then drop networkspeci
            if len(relevant_files) == 0:
                if networkspeci in self.read_instance.networkspecies:
                    self.read_instance.networkspecies.remove(networkspeci)
                    print('Warning: There is no available observational data for the network|species: {}. Dropping.'.format(networkspeci))
                elif networkspeci in self.read_instance.filter_networkspecies:
                    self.read_instance.filter_networkspecies.remove(networkspeci)
                    del self.read_instance.filter_species[networkspeci]
                    print('Warning: There is no available observational data for the filter network|species: {}. Dropping.'.format(networkspeci))
                continue
                
            # get basic metadata for speci
            # GHOST - read in parallel
            if self.read_instance.reading_ghost:
                
                # define arrays for storing speci metadata
                speci_station_references = []
                speci_station_names = []
                speci_station_longitudes = []
                speci_station_latitudes = []
                speci_station_measurement_altitudes = []

                # read metadata in parallel
                tuple_arguments = []
                for fname in relevant_files:
                    tuple_arguments.append((fname, self.read_instance.reading_ghost))
                pool = multiprocessing.Pool(self.read_instance.n_cpus)
                returned_data = pool.map(read_netcdf_metadata, tuple_arguments)
                pool.close()
                pool.join()

                # unzip returned data
                for returned_data_per_month in returned_data:
                    speci_station_references = np.append(speci_station_references, returned_data_per_month[0])
                    speci_station_names = np.append(speci_station_names, returned_data_per_month[1])
                    speci_station_longitudes = np.append(speci_station_longitudes, returned_data_per_month[2])
                    speci_station_latitudes = np.append(speci_station_latitudes, returned_data_per_month[3])
                    speci_station_measurement_altitudes = np.append(speci_station_measurement_altitudes, returned_data_per_month[4])

                speci_station_references, station_unique_indices = np.unique(speci_station_references, return_index=True)
                self.read_instance.station_references[networkspeci] = speci_station_references
                self.read_instance.station_names[networkspeci] = speci_station_names[station_unique_indices]
                self.read_instance.station_longitudes[networkspeci] = speci_station_longitudes[station_unique_indices]
                self.read_instance.station_latitudes[networkspeci] = speci_station_latitudes[station_unique_indices]
                self.read_instance.station_measurement_altitudes[networkspeci] = speci_station_measurement_altitudes[station_unique_indices]

            # non-GHOST - take from first file (metadata will be same across time)
            else:
                self.read_instance.nonghost_metadata_vars_to_read = []
            
                ncdf_root = Dataset(relevant_files[0])
                if 'station_reference' in ncdf_root.variables:
                    station_reference_var = 'station_reference'
                elif 'station_code' in ncdf_root.variables:
                    station_reference_var = 'station_code'
                elif 'station_name' in ncdf_root.variables:
                    station_reference_var = 'station_name'
                else: 
                    error = 'Error: {} cannot be read because it has no station_name.'.format(relevant_file)
                    sys.exit(error)
                self.read_instance.nonghost_metadata_vars_to_read.append('station_reference') 

                meta_shape = ncdf_root[station_reference_var].shape
                self.read_instance.station_references[networkspeci] = ncdf_root[station_reference_var][:]
                meta_val_dtype = np.array([self.read_instance.station_references[networkspeci][0]]).dtype
                if len(meta_shape) == 2:
                    if meta_val_dtype == np.dtype(object):
                        self.read_instance.station_references[networkspeci] = np.array([''.join(val) 
                                                                    for val in self.read_instance.station_references[networkspeci]])
                    else:
                        self.read_instance.station_references[networkspeci] = chartostring(self.read_instance.station_references[networkspeci])

                # get indices of all non-NaN stations (can be NaN for some non-GHOST files)
                non_nan_station_indices = np.array([ref_ii for ref_ii, ref in enumerate(self.read_instance.station_references[networkspeci]) 
                                                    if ref.lower() != 'nan'])
                self.read_instance.station_references[networkspeci] = self.read_instance.station_references[networkspeci][non_nan_station_indices]
                
                if "latitude" in ncdf_root.variables:
                    self.read_instance.station_longitudes[networkspeci] = ncdf_root['longitude'][non_nan_station_indices]
                    self.read_instance.station_latitudes[networkspeci] = ncdf_root['latitude'][non_nan_station_indices]
                else:
                    self.read_instance.station_longitudes[networkspeci] = ncdf_root['lon'][non_nan_station_indices]
                    self.read_instance.station_latitudes[networkspeci] = ncdf_root['lat'][non_nan_station_indices]
                self.read_instance.nonghost_metadata_vars_to_read.append('longitude') 
                self.read_instance.nonghost_metadata_vars_to_read.append('latitude') 

                if "station_name" in ncdf_root.variables:
                    meta_shape = ncdf_root['station_name'].shape
                    self.read_instance.station_names[networkspeci] = ncdf_root['station_name'][non_nan_station_indices]
                    meta_val_dtype = np.array([self.read_instance.station_names[networkspeci][0]]).dtype
                    if len(meta_shape) == 2:
                        if meta_val_dtype == np.dtype(object):
                            self.read_instance.station_names[networkspeci] = np.array([''.join(val) for val in self.read_instance.station_names[networkspeci]])
                        else:
                            self.read_instance.station_names[networkspeci] = chartostring(self.read_instance.station_names[networkspeci])
                    self.read_instance.nonghost_metadata_vars_to_read.append('station_name') 

                if "measurement_altitude" in ncdf_root.variables:
                    self.read_instance.station_measurement_altitudes[networkspeci] = ncdf_root['measurement_altitude'][non_nan_station_indices]
                    self.read_instance.nonghost_metadata_vars_to_read.append('measurement_altitude') 

                # get non-GHOST measurement units
                self.read_instance.nonghost_units[speci] = ncdf_root[speci].units

                # get list of nonghost variables to read (subset of GHOST variables)
                for ghost_metadata_var in self.read_instance.ghost_metadata_vars_to_read:
                    if (ghost_metadata_var in ncdf_root.variables) & (ghost_metadata_var not in self.read_instance.nonghost_metadata_vars_to_read):
                        self.read_instance.nonghost_metadata_vars_to_read.append(ghost_metadata_var) 

                ncdf_root.close()

        # if have more than 1 networkspecies (including filter networkspecies), and spatial_colocation is active,
        # then spatially colocate stations across species
        if (len((self.read_instance.networkspecies + self.read_instance.filter_networkspecies)) > 1) & (self.read_instance.spatial_colocation):
            # get intersecting station indices across species (handle both GHOST and non-GHOST cases)
            if self.read_instance.reading_ghost:
                intersecting_indices = spatial_colocation_ghost(self.read_instance.station_longitudes, 
                                                                self.read_instance.station_latitudes, 
                                                                self.read_instance.station_measurement_altitudes)
            else:
                intersecting_indices = spatial_colocation_nonghost(self.read_instance.station_references, 
                                                                   self.read_instance.station_longitudes, 
                                                                   self.read_instance.station_latitudes)
            
            # iterate through networkspecies specific intersecting indices, setting 
            for ns, ns_intersects in intersecting_indices.items():
                self.read_instance.station_references[ns] = self.read_instance.station_references[ns][ns_intersects]
                self.read_instance.station_longitudes[ns] = self.read_instance.station_longitudes[ns][ns_intersects]
                self.read_instance.station_latitudes[ns] = self.read_instance.station_latitudes[ns][ns_intersects]
                if ns in self.read_instance.station_measurement_altitudes:
                    self.read_instance.station_measurement_altitudes[ns] = self.read_instance.station_measurement_altitudes[ns][ns_intersects]
                if ns in self.read_instance.station_names:
                    self.read_instance.station_names[ns] = self.read_instance.station_names[ns][ns_intersects]

    def read_data(self, yearmonths_to_read, data_labels):
        """ Function that handles reading of observational/experiment data.

            :param yearmonths_to_read: list of yearmonths to read
            :type yearmonths_to_read: list
            :param data_labels: Data arrays to plot
            :type data_labels: list
        """

        # create arrays to share across processes (for parallel multiprocessing use)
        # this only works for numerical dtypes, i.e. not strings
        timestamp_array_shared = multiprocessing.RawArray(ctypes.c_int64, len(self.read_instance.timestamp_array))
        if (self.read_instance.reading_ghost) & (self.read_instance.observations_data_label in data_labels):
            flags_shared = multiprocessing.RawArray(ctypes.c_uint8, len(self.read_instance.flags))
        else:
            flags_shared = None
        # fill arrays
        timestamp_array_shared[:] = self.read_instance.timestamp_array
        if (self.read_instance.reading_ghost) & (self.read_instance.observations_data_label in data_labels):
            flags_shared[:] = self.read_instance.flags

        # create dictionary for saving files to read
        self.files_to_read = {}

        # iterate through networkspecies + filter_networkspecies
        for networkspeci in (self.read_instance.networkspecies + self.read_instance.filter_networkspecies):

            # determine if filter networkspecies or not
            if networkspeci in self.read_instance.filter_networkspecies:
                filter_read = True
            else:
                filter_read = False

            # get indivudual network and species strings
            network = networkspeci.split('|')[0]
            speci = networkspeci.split('|')[1]

            # add dictionary of files to read per network-speci 
            self.files_to_read[networkspeci] = {}

            # iterate through data labels
            for data_label in data_labels:

                # get raw data label (non-alias)
                data_label_raw = self.read_instance.data_labels_raw[self.read_instance.data_labels.index(data_label)]

                # get species matrix
                matrix = self.read_instance.parameter_dictionary[speci]['matrix']

                # get relevant file start dates
                # observations
                if data_label == self.read_instance.observations_data_label:
                    
                    # GHOST
                    if self.read_instance.reading_ghost:
                        file_root = '%s/%s/%s/%s/%s/%s_' % (self.read_instance.ghost_root, network,
                                                            self.read_instance.ghost_version,
                                                            self.read_instance.resolution, speci, speci)
                        try:
                            available_yearmonths = self.read_instance.available_observation_data[network][self.read_instance.resolution][matrix][speci]
                        except KeyError:
                            continue

                    # non-GHOST
                    else:
                        file_root = '%s/%s/%s/%s/%s_' % (self.read_instance.nonghost_root, network, 
                                                        self.read_instance.resolution, speci, speci)
                        try:
                            available_yearmonths = self.read_instance.available_observation_data[network][self.read_instance.resolution][matrix][speci]
                        except KeyError:
                            continue

                # experiments 
                else:
                    # if are reading filter species continue to next data_label
                    if filter_read:
                        continue 

                    elif '/' in network:
                        file_root = \
                            '%s/%s/%s/%s/%s/%s/%s_' % (self.read_instance.exp_root, self.read_instance.ghost_version, 
                                                        data_label_raw, self.read_instance.resolution, speci, 
                                                        network.replace('/', '-'), speci)
                    else:
                        file_root = \
                            '%s/%s/%s/%s/%s/%s/%s_' % (self.read_instance.exp_root, self.read_instance.ghost_version, 
                                                       data_label_raw, self.read_instance.resolution, speci, network, speci)
                    try:
                        available_yearmonths = self.read_instance.available_experiment_data[network][self.read_instance.resolution][speci][data_label_raw]
                    except KeyError:
                        continue

                # get intersection of yearmonths_to_read and available_yearmonths
                yearmonths_to_read_intersect = list(set(yearmonths_to_read) & set(available_yearmonths))
                self.files_to_read[networkspeci][data_label] = sorted([file_root+str(yyyymm)+'.nc' for yyyymm in yearmonths_to_read_intersect])
                
            # if active qa == default qa, no need to screen by QA, so inform reading function of this
            default_qa = get_default_qa(self.read_instance, speci)
            if self.read_instance.qa_per_species[speci] == default_qa:
                default_qa_active = True
            else:
                default_qa_active = False

            # create network/speci specific arrays to share across processes (for parallel multiprocessing use)
            # this only works for numerical dtypes, i.e. not strings
            if not filter_read:
                data_in_memory_shared_shape = (len(data_labels), len(self.read_instance.station_references[networkspeci]), len(self.read_instance.time_array))
            else:
                data_in_memory_shared_shape = (1, len(self.read_instance.station_references[networkspeci]), len(self.read_instance.time_array))
            data_in_memory_shared = multiprocessing.RawArray(ctypes.c_float, data_in_memory_shared_shape[0] * data_in_memory_shared_shape[1] * data_in_memory_shared_shape[2])  
            if (self.read_instance.reading_ghost) & (self.read_instance.observations_data_label in data_labels):
                qa_shared = multiprocessing.RawArray(ctypes.c_uint8, len(self.read_instance.qa_per_species[speci]))
                if not filter_read:
                    ghost_data_in_memory_shared_shape = (len(self.read_instance.ghost_data_vars_to_read), len(self.read_instance.station_references[networkspeci]), len(self.read_instance.time_array))
                    ghost_data_in_memory_shared = multiprocessing.RawArray(ctypes.c_float, ghost_data_in_memory_shared_shape[0] * ghost_data_in_memory_shared_shape[1] * ghost_data_in_memory_shared_shape[2])  
                else:
                    ghost_data_in_memory_shared_shape = None
                    ghost_data_in_memory_shared = None
            else:
                qa_shared = None
                ghost_data_in_memory_shared_shape = None
                ghost_data_in_memory_shared = None

            # wrap data_in_memory_shared and ghost_data_in_memory_shared as numpy arrays so we can easily manipulate the data.
            data_in_memory_shared_np = np.frombuffer(data_in_memory_shared, dtype=np.float32).reshape(data_in_memory_shared_shape)
            if (self.read_instance.reading_ghost) & (self.read_instance.observations_data_label in data_labels) & (not filter_read):
                ghost_data_in_memory_shared_np = np.frombuffer(ghost_data_in_memory_shared, dtype=np.float32).reshape(ghost_data_in_memory_shared_shape)

            # fill arrays
            if not filter_read:
                data_label_indices = [self.read_instance.data_labels.index(data_label) for data_label in data_labels]
                np.copyto(data_in_memory_shared_np, self.read_instance.data_in_memory[networkspeci][data_label_indices, :, :])
            else:
                np.copyto(data_in_memory_shared_np, self.read_instance.filter_data_in_memory[networkspeci][:, :])
            if (self.read_instance.reading_ghost) & (self.read_instance.observations_data_label in data_labels):      
                qa_shared[:] = self.read_instance.qa_per_species[speci]
                if not filter_read:
                    np.copyto(ghost_data_in_memory_shared_np, self.read_instance.ghost_data_in_memory[networkspeci])  

            # iterate and read species data in all relevant netCDF files (either in serial/parallel)

            # read data in parallel
            # setup pool of N workers on N CPUs
            pool = multiprocessing.Pool(self.read_instance.n_cpus, initializer=init_shared_vars_read_netcdf_data, 
                                        initargs=(data_in_memory_shared, data_in_memory_shared_shape, 
                                                  ghost_data_in_memory_shared, ghost_data_in_memory_shared_shape, 
                                                  timestamp_array_shared, qa_shared, flags_shared))
            # read netCDF files in parallel
            tuple_argument_fields = ['filename', 'station_references', 'station_names', 'speci', 
                                     'observations_data_label', 'data_label', 'data_labels', 
                                     'reading_ghost', 'ghost_data_vars_to_read', 
                                     'metadata_dtype', 'metadata_vars_to_read', 'default_qa_active', 'filter_read']
            tuple_arguments = []

            for data_label in self.files_to_read[networkspeci]:
                for fname in self.files_to_read[networkspeci][data_label]:
                    tuple_arguments.append((fname, self.read_instance.station_references[networkspeci], 
                                            self.read_instance.station_names[networkspeci], speci, 
                                            self.read_instance.observations_data_label, data_label, data_labels, 
                                            self.read_instance.reading_ghost, 
                                            self.read_instance.ghost_data_vars_to_read, 
                                            self.read_instance.metadata_dtype, 
                                            self.read_instance.metadata_vars_to_read,
                                            default_qa_active, filter_read))
  
            returned_data = pool.map(read_netcdf_data, tuple_arguments)

            pool.close()
            # wait for worker processes to terminate before continuing
            pool.join()
            
            # do not read data if there are not enough datasets (less than 2 timesteps)
            if not returned_data:
                continue

            # finalise assignment of non-filter species
            if not filter_read:
                # iterate through read file data and place metadata into full array as appropriate
                for returned_data_ii, returned_data_per_month in enumerate(returned_data):
                    returned_filename = tuple_arguments[returned_data_ii][tuple_argument_fields.index('filename')]
                    returned_data_label = tuple_arguments[returned_data_ii][tuple_argument_fields.index('data_label')]
                    returned_yearmonth = returned_filename.split('_')[-1][:6]
                    if returned_data_label == self.read_instance.observations_data_label:
                        # if returned_data_per_month is empty list, do not add
                        if len(returned_data_per_month) > 0:
                            self.read_instance.metadata_in_memory[networkspeci][:, self.read_instance.yearmonths.index(returned_yearmonth)] = returned_data_per_month[:, 0]

                # save to data in memory
                self.read_instance.data_in_memory[networkspeci][data_label_indices, :, :] = data_in_memory_shared_np
                if (self.read_instance.reading_ghost) & (self.read_instance.observations_data_label in data_labels):
                    self.read_instance.ghost_data_in_memory[networkspeci] = ghost_data_in_memory_shared_np

                # set data array for final validation checks
                data_array = self.read_instance.data_in_memory[networkspeci]

            # finalise assignment of filter species
            else:
                # save to filter data in memory
                self.read_instance.filter_data_in_memory[networkspeci][:, :] = data_in_memory_shared_np

                # set data array for final validation checks
                data_array = self.read_instance.filter_data_in_memory[networkspeci]

            # check if read data array consist of arrays full of -9999.0 or nan values or if they are empty
            if (data_array.size == 0) or \
                (np.isin(data_array.flatten(), [-9999.0, np.nan]).all()):

                if data_array.size == 0:
                    error = 'Error: The observation and experiment arrays for {} are empty.'.format(networkspeci)
            
                elif np.isin(data_array.flatten(), [-9999.0, np.nan]).all():
                    error = 'Error: All observation and experiment arrays for {} are void.'.format(networkspeci)

                sys.exit(error)
