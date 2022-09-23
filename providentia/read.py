from .read_aux import get_yearmonths_to_read, init_shared_vars_read_netcdf_data, read_netcdf_data
from .aux import check_for_ghost, get_basic_metadata, update_plotting_parameters

import sys
import os
import copy
import ctypes
import datetime
import multiprocessing
import time

from dateutil.relativedelta import relativedelta
import numpy as np
import pandas as pd
from netCDF4 import Dataset

class DataReader:
    """Class that reads observational/experiment data into memory."""

    def __init__(self, read_instance):
        self.read_instance = read_instance
        
    def read_setup(self, operations, experiments_to_remove=None, experiments_to_read=None):
        """
        Setup structures for read of new observational/experiment data and then perform read.

        :param operations: list of instructions on how to adjust data structures
        :type operation: list
        :param experiments_to_remove: list of experiments to remove from arrays
        :type experiments_to_remove: list
        :param experiments_to_read: list of experiments to add to arrays
        :type experiments_to_read: list
        """

        # changing time dimension ?
        if ('reset' in operations) or ('read_left' in operations) or ('read_right' in operations) or ('cut_left' in operations) or ('cut_right' in operations):

            # determine if reading GHOST or non-GHOST
            self.read_instance.reading_ghost = check_for_ghost(self.read_instance.network[0])

            # set active frequency code
            if (self.read_instance.resolution == 'hourly') or (self.read_instance.resolution == 'hourly_instantaneous'):
                self.active_frequency_code = 'H'
            elif (self.read_instance.resolution == '3hourly') or (self.read_instance.resolution == '3hourly_instantaneous'):
                self.active_frequency_code = '3H'
            elif (self.read_instance.resolution == '6hourly') or (self.read_instance.resolution == '6hourly_instantaneous'):
                self.active_frequency_code = '6H'
            elif self.read_instance.resolution == 'daily':
                self.active_frequency_code = 'D'
            elif self.read_instance.resolution == 'monthly':
                self.active_frequency_code = 'MS'

            # get time array
            self.read_instance.time_array = pd.date_range(start=datetime.datetime(int(str(self.read_instance.start_date)[:4]),
                                                                                int(str(self.read_instance.start_date)[4:6]),
                                                                                int(str(self.read_instance.start_date)[6:8])),
                                                          end=datetime.datetime(int(str(self.read_instance.end_date)[:4]),
                                                                                int(str(self.read_instance.end_date)[4:6]),
                                                                                int(str(self.read_instance.end_date)[6:8])),
                                                          freq=self.active_frequency_code)[:-1]

            # show warning when the data consists only of less than 2 timesteps
            if len(self.read_instance.time_array) < 2:
                self.read_instance.clear_canvas = True
                print('Warning: Extend the time range or enhance the resolution (e.g. from monthly to daily) to create plots.')
                return
            else:
                self.read_instance.clear_canvas = False

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

                # get unique basic metadata across networks / species
                self.read_instance.station_references, self.read_instance.station_longitudes, self.read_instance.station_latitudes =\
                    get_basic_metadata(self.read_instance, self.read_instance.network, self.read_instance.species, self.read_instance.resolution) 

        # need to reset all data structures 
        if 'reset' in operations:  

            # data
            self.read_instance.data_in_memory = {'{}|{}'.format(network,speci): 
                                                 np.full((len(self.read_instance.data_labels),
                                                          len(self.read_instance.station_references['{}|{}'.format(network,speci)]),
                                                          len(self.read_instance.time_array)),
                                                          np.NaN, dtype=np.float32) for network, speci in zip(self.read_instance.network, self.read_instance.species)}

            # GHOST data --> data variables which change per measurement (for filtering)
            if self.read_instance.reading_ghost:

                # set ghost data variables to read (dependent on data resolution)
                if (self.read_instance.resolution == 'hourly') or (self.read_instance.resolution == 'hourly_instantaneous'):
                    self.read_instance.ghost_data_vars_to_read = ['hourly_native_representativity_percent',
                                                                'daily_native_representativity_percent',
                                                                'monthly_native_representativity_percent',
                                                                'annual_native_representativity_percent', 'hourly_native_max_gap_percent',
                                                                'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                                                'annual_native_max_gap_percent', 'day_night_code', 'weekday_weekend_code',
                                                                'season_code']
                elif (self.read_instance.resolution == '3hourly') or \
                        (self.read_instance.resolution == '6hourly') or (self.read_instance.resolution == '3hourly_instantaneous') or \
                        (self.read_instance.resolution == '6hourly_instantaneous'):
                    self.read_instance.ghost_data_vars_to_read = ['daily_native_representativity_percent',
                                                                'monthly_native_representativity_percent',
                                                                'annual_native_representativity_percent',
                                                                'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                                                'annual_native_max_gap_percent', 'day_night_code', 'weekday_weekend_code',
                                                                'season_code']
                elif self.read_instance.resolution == 'daily':
                    self.read_instance.ghost_data_vars_to_read = ['daily_native_representativity_percent',
                                                                'monthly_native_representativity_percent',
                                                                'annual_native_representativity_percent',
                                                                'daily_native_max_gap_percent', 'monthly_native_max_gap_percent',
                                                                'annual_native_max_gap_percent', 'weekday_weekend_code', 'season_code']
                elif self.read_instance.resolution == 'monthly':
                    self.read_instance.ghost_data_vars_to_read = ['monthly_native_representativity_percent',
                                                                'annual_native_representativity_percent', 'monthly_native_max_gap_percent',
                                                                'annual_native_max_gap_percent', 'season_code']

                self.read_instance.ghost_data_in_memory = {'{}|{}'.format(network,speci):
                                        np.full((len(self.read_instance.ghost_data_vars_to_read),
                                        len(self.read_instance.station_references['{}|{}'.format(network,speci)]),
                                        len(self.read_instance.time_array)),
                                        np.NaN, dtype=np.float32) for network, speci in zip(self.read_instance.network, self.read_instance.species)} 

            else:
                self.read_instance.ghost_data_vars_to_read = []

            # metadata 
            #non-GHOST
            if not self.read_instance.reading_ghost:
                self.read_instance.metadata_dtype = [('station_name', np.object), ('latitude', np.float32),
                                                     ('longitude', np.float32), ('altitude', np.float32),
                                                     ('station_reference', np.object), ('station_classification', np.object),
                                                     ('area_classification', np.object)]
                self.read_instance.metadata_vars_to_read = [meta_dtype[0] for meta_dtype in self.read_instance.metadata_dtype]
            #GHOST
            else:
                self.read_instance.metadata_dtype = self.read_instance.ghost_metadata_dtype
                self.read_instance.metadata_vars_to_read = self.read_instance.ghost_metadata_vars_to_read

            self.read_instance.metadata_in_memory = {'{}|{}'.format(network,speci): 
                                                     np.full((len(self.read_instance.station_references['{}|{}'.format(network,speci)]),
                                                              len(self.read_instance.yearmonths)),
                                                              np.NaN, dtype=self.read_instance.metadata_dtype) for network, speci in zip(self.read_instance.network, self.read_instance.species)}

            # get list of yearmonths to read
            yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.start_date,
                                                        self.read_instance.end_date, self.read_instance.resolution)

            # read data 
            self.read_data(self.read_instance.network, self.read_instance.species, self.read_instance.resolution, yearmonths_to_read, self.read_instance.data_labels)

            # update measurement units for all species (take standard units for each speci from parameter dictionary)
            self.read_instance.measurement_units = {speci:self.read_instance.parameter_dictionary[speci]['standard_units'] for speci in self.read_instance.species}

            # reset plotting params
            self.read_instance.plotting_params = {}
            #iterate through data labels
            for data_label in self.read_instance.data_labels:
                self.read_instance.plotting_params[data_label] = {}
                # get experiment specific grid edges for exp, from first relevant file
                if data_label != 'observations':
                    exp_nc_root = Dataset(self.files_to_read['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])][data_label][0])
                    self.read_instance.plotting_params[data_label]['grid_edge_longitude'] = exp_nc_root['grid_edge_longitude'][:]
                    self.read_instance.plotting_params[data_label]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
                    exp_nc_root.close()
            
            # update plotting parameters colours and zorder
            update_plotting_parameters(self.read_instance) 

        # need to read on left / read on right / cut on left / cut on right (for dashboard)
        if ('read_left' in operations) or ('read_right' in operations) or ('cut_left' in operations) or ('cut_right' in operations):  

            # if station references array has changed then as cutting / appending to
            # need to rearrange existing metadata/data arrays accordingly
            if not np.array_equal(self.read_instance.previous_station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])], \
                                    self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]):

                # get indices of stations in previous station references array in current station references array
                old_station_inds = np.where(np.in1d(self.read_instance.previous_station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])],
                                                    self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]))[0]
                                                    
                # get indices of stations in current station references array
                # that were in previous station references array
                new_station_inds = np.where(np.in1d(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])],
                                                    self.read_instance.previous_station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]))[0]

                #rearrange metadata station dimension
                new_metadata_in_memory = np.full((len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]),
                                                  len(self.read_instance.previous_yearmonths)),
                                                  np.NaN, dtype=self.read_instance.metadata_dtype)
                new_metadata_in_memory[new_station_inds, :] = self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])][old_station_inds, :]
                self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = new_metadata_in_memory

                # rearrage data array station dimension
                new_data_in_memory = np.full((len(self.read_instance.previous_data_labels),
                                              len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]),
                                              len(self.read_instance.previous_time_array)),
                                              np.NaN, dtype=np.float32)
                # put the old data into new array in the correct positions
                new_data_in_memory[:, new_station_inds, :] = self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])][:, old_station_inds, :]
                # overwrite data array with reshaped version
                self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = new_data_in_memory

                # rearrage ghost data array station dimension
                if self.read_instance.reading_ghost:
                    new_ghost_data_in_memory = np.full((len(self.read_instance.ghost_data_vars_to_read),
                                                        len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]),
                                                        len(self.read_instance.previous_time_array)),
                                                        np.NaN, dtype=np.float32)
                    # put the old ghost data into new array in the correct positions
                    new_ghost_data_in_memory[:, new_station_inds, :] = self.read_instance.ghost_data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])][:, old_station_inds, :]
                    # overwrite ghost data array with reshaped version
                    self.read_instance.ghost_data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = new_ghost_data_in_memory

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
                    self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                        self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])][:, [metadata_left_edge_ind]]
                else:
                    self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                        self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])][:, metadata_left_edge_ind:metadata_right_edge_ind]

                # cut edges of data array
                self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                    self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])][:, :, data_left_edge_ind:data_right_edge_ind]

                #cut edges of ghost data array
                if self.read_instance.reading_ghost:
                    self.read_instance.ghost_data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                        self.read_instance.ghost_data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])][:, :, data_left_edge_ind:data_right_edge_ind]

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

                    #need to read new yearmonths?
                    if len(yearmonths_to_read) > 0:

                        #add space for new data on left edge of the metadata array
                        self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                            np.concatenate((np.full((len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), len(yearmonths_to_read)), 
                                np.NaN, dtype=self.read_instance.metadata_dtype), self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), axis=1)

                        # insert space for new data on left edge of the data array
                        self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                            np.concatenate((np.full((len(self.read_instance.previous_data_labels), len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), n_new_left_data_inds), 
                                np.NaN, dtype=np.float32), self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), axis=2)

                        # insert space for new ghost data on left edge of the ghost data array
                        if self.read_instance.reading_ghost:
                            self.read_instance.ghost_data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                                np.concatenate((np.full((len(self.read_instance.ghost_data_vars_to_read), len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), n_new_left_data_inds), 
                                    np.NaN, dtype=np.float32), self.read_instance.ghost_data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), axis=2)
                        
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

                    #need to read new yearmonths?
                    if len(yearmonths_to_read) > 0:

                        #add space for new data on right edge of the metadata array
                        self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                            np.concatenate((self.read_instance.metadata_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])], 
                                np.full((len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), len(yearmonths_to_read)), 
                                    np.NaN, dtype=self.read_instance.metadata_dtype)), axis=1)

                        # insert space for new data on right edge of the data array
                        self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                            np.concatenate((self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])], 
                                np.full((len(self.read_instance.previous_data_labels), len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), n_new_right_data_inds), 
                                    np.NaN, dtype=np.float32)), axis=2)

                        # insert space for new ghost data on right edge of the ghost data array
                        if self.read_instance.reading_ghost:
                            self.read_instance.ghost_data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                                np.concatenate((self.read_instance.ghost_data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])], 
                                    np.full((len(self.read_instance.ghost_data_vars_to_read), len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), n_new_right_data_inds), 
                                        np.NaN, dtype=np.float32)), axis=2)    

                        # add yearmonths_to_read to list for both edges
                        all_yearmonths_to_read.extend(yearmonths_to_read)

                # read data 
                self.read_data(self.read_instance.network, self.read_instance.species, self.read_instance.resolution, all_yearmonths_to_read, self.read_instance.previous_data_labels) 

        # need to remove experiment/s ?
        if 'remove_exp' in operations: 

            experiments_to_remove_inds = [self.read_instance.previous_data_labels.index(experiment) for experiment in experiments_to_remove]
            self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                np.delete(self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])], experiments_to_remove_inds, axis=0)

        # need to read experiment/s ? 
        if 'read_exp' in operations: 

            # insert space for new experiments in data array
            self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])] = \
                np.concatenate((self.read_instance.data_in_memory['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])], \
                np.full((len(experiments_to_read), len(self.read_instance.station_references['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])]), len(self.read_instance.time_array)), 
                    np.NaN, dtype=np.float32)), axis=0)

            # get list of yearmonths to read
            yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.start_date,
                                                        self.read_instance.end_date, self.read_instance.resolution)

            # read data 
            self.read_data(self.read_instance.network, self.read_instance.species, self.read_instance.resolution, yearmonths_to_read, experiments_to_read)       

            # add experiment specific grid edges for exp to plotting params
            for experiment_to_read in experiments_to_read:
                self.read_instance.plotting_params[experiment_to_read] = {}
                exp_nc_root = Dataset(self.files_to_read['{}|{}'.format(self.read_instance.network[0],self.read_instance.species[0])][experiment_to_read][0])
                self.read_instance.plotting_params[experiment_to_read]['grid_edge_longitude'] = \
                    exp_nc_root['grid_edge_longitude'][:]
                self.read_instance.plotting_params[experiment_to_read]['grid_edge_latitude'] = exp_nc_root['grid_edge_latitude'][:]
                exp_nc_root.close()

            # update plotting parameters colours and zorder
            update_plotting_parameters(self.read_instance) 

    def read_data(self, networks, species, resolution, yearmonths_to_read, data_labels, filter=False):
        """
        Function that handles reading of observational/experiment data.

        :param networks: list of networks to read
        :type networks: list
        :param species: list of species to read
        :type species: list
        :param resolution: temporal resolution to read
        :type resolution: str
        :param yearmonths_to_read: list of yearmonths to read
        :type yearmonths_to_read: list
        :param data_labels: data labels to read
        :type data_labels: list
        :param filter: boolean informing whether are reading a variable to filter by or not
        :type filter: boolean
        """

        # create arrays to share across processes (for parallel multiprocessing use)
        # this only works for numerical dtypes, i.e. not strings
        timestamp_array_shared = multiprocessing.RawArray(ctypes.c_int64, len(self.read_instance.timestamp_array))
        if (self.read_instance.reading_ghost) & ('observations' in data_labels):
            flags_shared = multiprocessing.RawArray(ctypes.c_uint8, len(self.read_instance.flags))
        else:
            flags_shared = None
        # fill arrays
        timestamp_array_shared[:] = self.read_instance.timestamp_array
        if (self.read_instance.reading_ghost) & ('observations' in data_labels):
            flags_shared[:] = self.read_instance.flags

        # create dictionary for saving files to read
        self.files_to_read = {}

        # iterate through networks and species
        for network, speci in zip(networks, species):

            # add dictionary of files to read per network-speci 
            self.files_to_read['{}|{}'.format(network,speci)] = {}

            # iterate through data labels
            for data_label in data_labels:

                # get species matrix
                matrix = self.read_instance.parameter_dictionary[speci]['matrix']

                # get relevant file start dates
                # observations
                if data_label == 'observations':
                    
                    # GHOST
                    if self.read_instance.reading_ghost:
                        file_root = '%s/%s/%s/%s/%s/%s_' % (self.read_instance.ghost_root, network,
                                                            self.read_instance.ghost_version,
                                                            resolution, speci, speci)
                        try:
                            available_yearmonths = self.read_instance.available_observation_data[network][resolution][matrix][speci]
                        except KeyError:
                            error = 'Error: The folder {0} does not exist.'.format(file_root[:-1])
                            tip = 'Tip: Consider interpolating the network data for the given configuration using Providentia Interpolation.'
                            print(error + '\n' + tip)
                            sys.exit()

                    # non-GHOST
                    else:
                        file_root = '%s/%s/%s/%s/%s_' % (self.read_instance.nonghost_root, network, 
                                                         resolution, speci, speci)
                        try:
                            available_yearmonths = self.read_instance.available_observation_data[network][resolution][matrix][speci]
                        except KeyError:
                            error = 'Error: The folder {0} does not exist.'.format(file_root[:-1])
                            tip = 'Tip: Consider interpolating the network data for the given configuration using Providentia Interpolation.'
                            print(error + '\n' + tip)
                            sys.exit()

                # experiments
                else:
                    file_root = \
                        '%s/%s/%s/%s/%s/%s/%s_' % (self.read_instance.exp_root, self.read_instance.ghost_version, 
                                                   data_label, resolution, speci, network, speci)
                    try:
                        available_yearmonths = self.read_instance.available_experiment_data[network][resolution][speci][data_label]
                    except KeyError:
                        error = 'Error: The folder {0} does not exist.'.format(file_root[:-1])
                        tip = 'Tip: Consider interpolating the network data for the given configuration using Providentia Interpolation.'
                        print(error + '\n' + tip)
                        sys.exit()

                # get intersection of yearmonths_to_read and available_yearmonths
                yearmonths_to_read_intersect = list(set(yearmonths_to_read) & set(available_yearmonths))
                self.files_to_read['{}|{}'.format(network,speci)][data_label] = sorted([file_root+str(yyyymm)+'.nc' for yyyymm in yearmonths_to_read_intersect])

            # if active qa == default qa, no need to screen by qa, so set selected qa to None
            if speci in self.read_instance.met_parameters:
                default_qa = self.read_instance.default_qa_met
            else:
                default_qa = self.read_instance.default_qa_standard
            if self.read_instance.qa == default_qa:
                qa_to_filter = []
            else:
                qa_to_filter = self.read_instance.qa 

            # create network/ speci specific arrays to share across processes (for parallel multiprocessing use)
            # this only works for numerical dtypes, i.e. not strings
            data_in_memory_shared_shape = (len(data_labels), len(self.read_instance.station_references['{}|{}'.format(network,speci)]), len(self.read_instance.time_array))
            data_in_memory_shared = multiprocessing.RawArray(ctypes.c_float, data_in_memory_shared_shape[0] * data_in_memory_shared_shape[1] * data_in_memory_shared_shape[2])  
            if (self.read_instance.reading_ghost) & ('observations' in data_labels):
                ghost_data_in_memory_shared_shape = (len(self.read_instance.ghost_data_vars_to_read), len(self.read_instance.station_references['{}|{}'.format(network,speci)]), len(self.read_instance.time_array))
                ghost_data_in_memory_shared = multiprocessing.RawArray(ctypes.c_float, ghost_data_in_memory_shared_shape[0] * ghost_data_in_memory_shared_shape[1] * ghost_data_in_memory_shared_shape[2])  
                qa_shared = multiprocessing.RawArray(ctypes.c_uint8, len(qa_to_filter))
            else:
                ghost_data_in_memory_shared_shape = None
                ghost_data_in_memory_shared = None
                qa_shared = None
            # wrap data_in_memory_shared and ghost_data_in_memory_shared as numpy arrays so we can easily manipulate the data.
            data_in_memory_shared_np = np.frombuffer(data_in_memory_shared, dtype=np.float32).reshape(data_in_memory_shared_shape)
            if (self.read_instance.reading_ghost) & ('observations' in data_labels):
                ghost_data_in_memory_shared_np = np.frombuffer(ghost_data_in_memory_shared, dtype=np.float32).reshape(ghost_data_in_memory_shared_shape)
            # fill arrays
            data_label_indices = [self.read_instance.data_labels.index(data_label) for data_label in data_labels]
            np.copyto(data_in_memory_shared_np, self.read_instance.data_in_memory['{}|{}'.format(network,speci)][data_label_indices, :, :])
            if (self.read_instance.reading_ghost) & ('observations' in data_labels):
                np.copyto(ghost_data_in_memory_shared_np, self.read_instance.ghost_data_in_memory['{}|{}'.format(network,speci)])        
                qa_shared[:] = qa_to_filter

            # iterate and read species data in all relevant netCDF files (either in serial/parallel)
            s = time.time()

            # read data in parallel
            # setup pool of N workers on N CPUs
            pool = multiprocessing.Pool(self.read_instance.n_cpus, initializer=init_shared_vars_read_netcdf_data, initargs=(data_in_memory_shared, data_in_memory_shared_shape, ghost_data_in_memory_shared, ghost_data_in_memory_shared_shape, timestamp_array_shared, qa_shared, flags_shared))
            # read netCDF files in parallel
            tuple_argument_fields = ['filename','station_references','speci','data_label','data_labels','reading_ghost','ghost_data_vars_to_read','metadata_dtype','metadata_vars_to_read']
            tuple_arguments = []
            for data_label in self.files_to_read['{}|{}'.format(network,speci)]:
                for fname in self.files_to_read['{}|{}'.format(network,speci)][data_label]:
                    tuple_arguments.append((fname, self.read_instance.station_references['{}|{}'.format(network,speci)], speci, data_label, data_labels,
                                            self.read_instance.reading_ghost, self.read_instance.ghost_data_vars_to_read, self.read_instance.metadata_dtype, self.read_instance.metadata_vars_to_read))

            returned_data = pool.map(read_netcdf_data, tuple_arguments)

            pool.close()
            # wait for worker processes to terminate before continuing
            pool.join()
            
            # iterate through read file data and place metadata into full array as appropriate
            for returned_data_ii, returned_data_per_month in enumerate(returned_data):
                returned_filename = tuple_arguments[returned_data_ii][tuple_argument_fields.index('filename')]
                returned_data_label = tuple_arguments[returned_data_ii][tuple_argument_fields.index('data_label')]
                returned_yearmonth = returned_filename.split('_')[-1][:6]
                if returned_data_label == 'observations':
                    self.read_instance.metadata_in_memory['{}|{}'.format(network,speci)][:, self.read_instance.yearmonths.index(returned_yearmonth)] = returned_data_per_month[:, 0]

            # save to data in memory
            self.read_instance.data_in_memory['{}|{}'.format(network,speci)][data_label_indices, :, :] = data_in_memory_shared_np
            if (self.read_instance.reading_ghost) & ('observations' in data_labels):
                self.read_instance.ghost_data_in_memory['{}|{}'.format(network,speci)] = ghost_data_in_memory_shared_np

            # check if datasets consist of arrays full of -9999.0 or nan values or if they are empty
            if (self.read_instance.data_in_memory['{}|{}'.format(network,speci)].size == 0) or \
                (np.isin(self.read_instance.data_in_memory['{}|{}'.format(network,speci)].flatten(), [-9999.0, np.nan]).all()):

                if self.read_instance.data_in_memory['{}|{}'.format(network,speci)].size == 0:
                    error = 'Error: The observation or experiment datasets are empty.'
            
                elif np.isin(self.read_instance.data_in_memory['{}|{}'.format(network,speci)].flatten(), [-9999.0, np.nan]).all():
                    error = 'Error: The observation or experiment datasets are void.'

                tip = 'Tip: Check if the data from the observations was downloaded correctly and if the experiments were interpolated at the stations of the network of interest.'
                sys.exit(error + '\n' + tip)
