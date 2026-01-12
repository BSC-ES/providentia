""" Class for reading data into memory """

import copy
import ctypes
import datetime
import multiprocessing
import os
import sys

from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
import numpy as np
import pandas as pd

from providentia.auxiliar import get_standard_units, get_standard_parameters_by_speci
from .fields_menus import (init_representativity, init_period, init_metadata,
                           update_representativity_fields, update_period_fields, update_metadata_fields,
                           representativity_conf, period_conf, metadata_conf)
from .plot_aux import update_plotting_parameters
from .read_aux import (check_for_ghost, get_default_qa, get_frequency_code, get_yearmonths_to_read, 
                       init_shared_vars_read_netcdf_data, read_netcdf_data, read_netcdf_metadata, 
                       check_forecast_dimension)
from .spatial_colocation import SpatialColocation
from .warnings_prv import show_message


class DataReader:
    """ Class for reading data into memory """

    def __init__(self, read_instance):
        """
        Initialise the DataReader instance.

        Parameters
        ----------
        read_instance : object
            An object instance to be updated with data reading configurations and state.
        """

        self.read_instance = read_instance
        
    def read_setup(self, operations, models_to_remove=None, models_to_read=None):
        """
        Setup structures for read of new observational/model data and then perform read.

        Parameters
        ----------
        operations : list
            A list of instructions on how to adjust data structures (e.g., 'reset', 'read_left', 'cut_right').
        models_to_remove : list, optional
            A list of models to remove from arrays.
        models_to_read : list, optional
            A list of models to add to arrays.
        """

        # changing time dimension ?
        if ('reset' in operations) or ('read_left' in operations) or ('read_right' in operations) or ('cut_left' in operations) or ('cut_right' in operations):

            # if reset is not in operations then cannot be reading from conf file
            if 'reset' not in operations:
                self.read_instance.from_conf = False

            # determine if reading GHOST or non-GHOST
            self.read_instance.reading_ghost = check_for_ghost(self.read_instance.network[0])

            # re-initialise and update fields for active resolution
            init_representativity(self.read_instance)
            init_period(self.read_instance)
            update_representativity_fields(self.read_instance)
            update_period_fields(self.read_instance)
            # if loading from a conf file then load new fields
            if self.read_instance.from_conf:
                representativity_conf(self.read_instance)
                period_conf(self.read_instance)

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
                if (self.read_instance.from_conf) and (self.read_instance.mode not in ['report', 'library']):
                    error = 'Error: Providentia will not be launched.'
                    self.read_instance.logger.error(error)
                    sys.exit(1) 
                elif (self.read_instance.mode == 'report'):
                    error = 'Error: Report will not be created.'
                    self.read_instance.logger.error(error)
                    sys.exit(1) 
                elif (self.read_instance.mode == 'library'):
                    error = 'Error: Data cannot be read.'
                    self.read_instance.logger.error(error)
                    sys.exit(1) 
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

                # get N indices per year
                self.read_instance.N_inds_per_year = np.array([np.count_nonzero(np.all(
                    [self.read_instance.time_array >= datetime.datetime.strptime(str(year)+'0101','%Y%m%d'),
                    self.read_instance.time_array < datetime.datetime.strptime(str(year+1)+'0101','%Y%m%d')], axis=0)) 
                    for year in np.unique([int(str(yyyymm)[:4]) for yyyymm in self.read_instance.yearmonths])])

                # get unique basic metadata across networkspecies
                # for this step include filter networkspecies
                self.read_basic_metadata()

                # iterate through station_references per networkspecies
                # if have 0 valid stations then drop  
                for networkspeci, stn_refs in self.read_instance.station_references.items():
                    if len(stn_refs) == 0:
                        if networkspeci in self.read_instance.networkspecies:
                            self.read_instance.networkspecies.remove(networkspeci)
                            self.read_instance.species.remove(networkspeci.split('|')[1])
                            msg = 'There is no available observational data for the network|species: {}. Dropping.'.format(networkspeci)
                            show_message(self.read_instance, msg)
                        elif networkspeci in self.read_instance.filter_networkspecies:
                            self.read_instance.filter_networkspecies.remove(networkspeci)
                            self.read_instance.filter_species.remove(networkspeci.split('|')[1])
                            msg = 'There is no available observational data for the filter network|species: {}. Dropping.'.format(networkspeci)
                            show_message(self.read_instance, msg)

                # if have zero networkspecies left, then return with invalid_read
                if len(self.read_instance.networkspecies) == 0:
                    self.read_instance.invalid_read = True
                    return

                # set invalid_read to be False if have data to read
                self.read_instance.invalid_read = False

        # need to reset all data structures 
        if 'reset' in operations:  

            # uninitialise filter object
            if self.read_instance.mode not in ['report', 'library']:
                self.read_instance.mpl_canvas.filter_data = None

            # get list of yearmonths to read
            yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.start_date,
                                                        self.read_instance.end_date, self.read_instance.resolution)

            # check if any of the model data has a forecast dimension to handle 
            # only for report / library modes as dashboard handled previously
            if self.read_instance.mode in ['report', 'library']:

                self.read_instance.original_data_labels = copy.deepcopy(self.read_instance.data_labels)
                self.read_instance.original_data_labels_raw = copy.deepcopy(self.read_instance.data_labels_raw)
                self.read_instance.original_models = copy.deepcopy(self.read_instance.experiments)

                self.check_forecast(yearmonths_to_read=yearmonths_to_read)
                self.read_instance.data_labels, self.read_instance.data_labels_raw, self.read_instance.experiments = self.update_forecast_indices(init=True)

            # create data in memory array
            self.read_instance.data_in_memory = {networkspeci: 
                                                np.full((len(self.read_instance.data_labels),
                                                        len(self.read_instance.station_references[networkspeci]),
                                                        len(self.read_instance.time_array)),
                                                        np.nan, dtype=np.float32) for networkspeci in self.read_instance.networkspecies} 

            # filter data (if active)
            if self.read_instance.filter_species:
                self.read_instance.filter_data_in_memory = {networkspeci: 
                                                            np.full((len(self.read_instance.station_references[networkspeci]),
                                                                     len(self.read_instance.time_array)),
                                                                     np.nan, dtype=np.float32) for networkspeci in self.read_instance.filter_networkspecies}
            else:
                self.read_instance.filter_data_in_memory = {}

            # initialise variables and get representativity resolutions
            self.read_instance.ghost_data_in_memory = {}
            self.read_instance.ghost_data_vars_to_read = []
            if self.read_instance.resolution in ['hourly', 'hourly_instantaneous']:
                resolution = 'hourly'
            elif self.read_instance.resolution in ['daily', '3hourly', '6hourly', '3hourly_instantaneous', '6hourly_instantaneous']:
                resolution = 'daily'
            elif self.read_instance.resolution == 'monthly':
                resolution = 'monthly'

            # get data variables which change per measurement (for filtering)
            if self.read_instance.reading_ghost:
                # get representativity fields (only native because non-native and all are calculated on-the-fly)
                self.read_instance.ghost_data_vars_to_read = [var for var 
                                                              in self.read_instance.representativity_info['ghost'][resolution]['map_vars'] 
                                                              if 'native' in var]
                # add annual native representativity percent    
                self.read_instance.ghost_data_vars_to_read.append('annual_native_representativity_percent') 

                # add periodic code variables
                self.read_instance.ghost_data_vars_to_read.append('season_code')
                if self.read_instance.resolution != 'monthly':
                    self.read_instance.ghost_data_vars_to_read.append('weekday_weekend_code')
                if self.read_instance.resolution not in ['monthly', 'daily']:
                    self.read_instance.ghost_data_vars_to_read.append('day_night_code')

                # initialise data in memory for GHOST with NaN for these variables
                self.read_instance.ghost_data_in_memory = {networkspeci:
                                        np.full((len(self.read_instance.ghost_data_vars_to_read),
                                                 len(self.read_instance.station_references[networkspeci]),
                                                 len(self.read_instance.time_array)),
                                                 np.nan, dtype=np.float32) for networkspeci in self.read_instance.networkspecies} 
            
            if self.read_instance.reading_ghost:
                representativity_valid_vars = self.read_instance.representativity_info['ghost'][resolution]['map_vars']
            else:
                representativity_valid_vars = self.read_instance.representativity_info['nonghost'][resolution]['map_vars']

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
                # get all the valid args in one list
                period_set = ['period'] if self.read_instance.reading_ghost else []
                valid_fields = self.read_instance.metadata_vars_to_read + self.read_instance.ghost_data_vars_to_read + \
                    period_set + representativity_valid_vars

                # remove all the valid fields from the invalid field list
                self.read_instance.invalid_fields = {field_name: fields-set(valid_fields)
                    for field_name, fields in self.read_instance.non_default_fields_per_section.items() 
                    if field_name==self.read_instance.section or field_name.startswith(self.read_instance.section+"·")}

                # show warning if there's an invalid field
                invalid_var = [f"""{i} ('{"', '".join(j)}')""" for i,j in self.read_instance.invalid_fields.items() if j]
                if invalid_var:
                    msg = f"Invalid field(s) in configuration file {self.read_instance.config.split('/')[-1]}. "
                    msg += f"Section(s) and Field(s): {', '.join(invalid_var)}."
                    show_message(self.read_instance, msg)

                    # delete from instance all invalid fields from the configuration file
                    for section_invalid_fields in self.read_instance.invalid_fields.values():
                        for k in section_invalid_fields:
                            # control if the atribute exists because in report mode the subsection ones are not set yet
                            if hasattr(self.read_instance, k):
                                delattr(self.read_instance, k)                 

            self.read_instance.metadata_in_memory = {networkspeci: 
                                                     np.full((len(self.read_instance.station_references['{}'.format(networkspeci)]),
                                                              len(self.read_instance.yearmonths)),
                                                              np.nan, dtype=self.read_instance.metadata_dtype) 
                                                              for networkspeci in self.read_instance.networkspecies}

            # read data 
            self.read_data(yearmonths_to_read, self.read_instance.data_labels)

            # update measurement units for all species (take standard units for each speci from parameter dictionary)
            # non-GHOST
            if not self.read_instance.reading_ghost:

                # convert non-GHOST units to standard format
                nonghost_standard_units = {}
                for speci in self.read_instance.nonghost_units.keys():
                    input_units = self.read_instance.nonghost_units[speci]
                    standard_parameter_speci = get_standard_parameters_by_speci(speci, self.read_instance.ghost_version)
                    standard_input_units = get_standard_units(input_units, standard_parameter_speci)
                    if 'Error:' in standard_input_units:
                        self.read_instance.logger.error(standard_input_units)
                        sys.exit(1)
                    nonghost_standard_units[speci] = standard_input_units
                self.read_instance.measurement_units = {speci.split('|')[1]:nonghost_standard_units[speci.split('|')[1]] 
                                                        for speci in self.read_instance.networkspecies}
            # GHOST
            else:
                self.read_instance.measurement_units = {speci.split('|')[1]:self.read_instance.parameter_dictionary[speci.split('|')[1]]['standard_units'] 
                                                        for speci in self.read_instance.networkspecies}

            # update plotting parameters colours, zorder and model grid edges
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
                                                  np.nan, dtype=self.read_instance.metadata_dtype)
                new_metadata_in_memory[new_station_inds, :] = self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]][old_station_inds, :]
                self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]] = new_metadata_in_memory

                # rearrage data array station dimension
                new_data_in_memory = np.full((len(self.read_instance.previous_data_labels),
                                              len(self.read_instance.station_references[self.read_instance.networkspecies[0]]),
                                              len(self.read_instance.previous_time_array)),
                                              np.nan, dtype=np.float32)
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
                                                             np.nan, dtype=np.float32)
                        # put the old data into new array in the correct positions
                        new_filter_data_in_memory[new_station_inds, :] = self.read_instance.filter_data_in_memory[filter_networkspeci][old_station_inds, :]
                        # overwrite data array with reshaped version
                        self.read_instance.filter_data_in_memory[filter_networkspeci] = new_filter_data_in_memory

                # rearrage ghost data array station dimension
                if self.read_instance.reading_ghost:
                    new_ghost_data_in_memory = np.full((len(self.read_instance.ghost_data_vars_to_read),
                                                        len(self.read_instance.station_references[self.read_instance.networkspecies[0]]),
                                                        len(self.read_instance.previous_time_array)),
                                                        np.nan, dtype=np.float32)
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

                    # count number of yermonths not in previous yearmonths
                    n_new_yearmonths = np.count_nonzero(~np.isin(yearmonths_to_read, self.read_instance.previous_yearmonths))

                    # add space for new data on left edge of the metadata array (if needed)
                    if n_new_yearmonths > 0:
                        self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((np.full((len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_yearmonths), 
                                np.nan, dtype=self.read_instance.metadata_dtype), self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]]), axis=1)

                    # insert space for new data on left edge of the data array
                    self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = \
                        np.concatenate((np.full((len(self.read_instance.previous_data_labels), len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_left_data_inds), 
                            np.nan, dtype=np.float32), self.read_instance.data_in_memory[self.read_instance.networkspecies[0]]), axis=2)

                    # add space for new data on left edge of the filter data array
                    if self.read_instance.filter_species:  
                        self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((np.full((len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_left_data_inds), 
                                np.nan, dtype=np.float32), self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]]), axis=1)

                    # insert space for new ghost data on left edge of the ghost data array
                    if self.read_instance.reading_ghost:
                        self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((np.full((len(self.read_instance.ghost_data_vars_to_read), len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_left_data_inds), 
                                np.nan, dtype=np.float32), self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]]), axis=2)
                    
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

                    # count number of yermonths not in previous yearmonths
                    n_new_yearmonths = np.count_nonzero(~np.isin(yearmonths_to_read, self.read_instance.previous_yearmonths))

                    # add space for new data on right edge of the metadata array (if needed)
                    if n_new_yearmonths > 0:
                        self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((self.read_instance.metadata_in_memory[self.read_instance.networkspecies[0]], 
                                np.full((len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_yearmonths), 
                                    np.nan, dtype=self.read_instance.metadata_dtype)), axis=1)
                    
                    # insert space for new data on right edge of the data array
                    self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = \
                        np.concatenate((self.read_instance.data_in_memory[self.read_instance.networkspecies[0]], 
                            np.full((len(self.read_instance.previous_data_labels), len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_right_data_inds), 
                                np.nan, dtype=np.float32)), axis=2)

                    # insert space for new data on right edge of the filter data array
                    if self.read_instance.filter_species: 
                        self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((self.read_instance.filter_data_in_memory[self.read_instance.networkspecies[0]], 
                                np.full((len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_right_data_inds), 
                                    np.nan, dtype=np.float32)), axis=1)

                    # insert space for new ghost data on right edge of the ghost data array
                    if self.read_instance.reading_ghost:
                        self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]] = \
                            np.concatenate((self.read_instance.ghost_data_in_memory[self.read_instance.networkspecies[0]], 
                                np.full((len(self.read_instance.ghost_data_vars_to_read), len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), n_new_right_data_inds), 
                                    np.nan, dtype=np.float32)), axis=2)    

                    # add yearmonths_to_read to list for both edges
                    all_yearmonths_to_read.extend(yearmonths_to_read)

                # read data 
                self.read_data(all_yearmonths_to_read, self.read_instance.data_labels) 

        # need to remove model/s ?
        if 'remove_mod' in operations: 

            # get indices of models to remove
            models_to_remove_inds = [self.read_instance.previous_data_labels.index(model) for model in models_to_remove]
            
            # remove model data
            self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = \
                np.delete(self.read_instance.data_in_memory[self.read_instance.networkspecies[0]], models_to_remove_inds, axis=0)

        # need to read model/s ? 
        if 'read_mod' in operations: 

            # insert space for new models in data array
            for model_to_read in models_to_read:
                models_to_read_ind = self.read_instance.data_labels.index(model_to_read) 
                
                self.read_instance.data_in_memory[self.read_instance.networkspecies[0]] = \
                    np.insert(self.read_instance.data_in_memory[self.read_instance.networkspecies[0]], 
                                models_to_read_ind,
                                np.full((1, len(self.read_instance.station_references[self.read_instance.networkspecies[0]]), len(self.read_instance.time_array)), 
                                np.nan, dtype=np.float32),                     
                                axis=0)

            # get list of yearmonths to read
            yearmonths_to_read = get_yearmonths_to_read(self.read_instance.yearmonths, self.read_instance.start_date,
                                                        self.read_instance.end_date, self.read_instance.resolution)

            # read data
            self.read_data(yearmonths_to_read, models_to_read)       

        # if removing or adding models, update plotting parameters colours, zorder and model grid edges
        if ('remove_mod' in operations) or ('read_mod' in operations): 
            data_labels_to_remove = [mod for mod in models_to_remove if mod in list(self.read_instance.plotting_params.keys())]
            data_labels_to_add = [mod for mod in models_to_read if mod not in list(self.read_instance.plotting_params.keys())]
            update_plotting_parameters(self.read_instance, data_labels_to_remove=data_labels_to_remove, 
                                       data_labels_to_add=data_labels_to_add)

        # for non-GHOST delete valid station indices variables because we do not want to 
        # remove the stations with 0 valid measurements before the filter has been updated, 
        # this will happen later
        if hasattr(self.read_instance, 'valid_station_inds') and (not self.read_instance.reading_ghost):
            delattr(self.read_instance, 'valid_station_inds')
            delattr(self.read_instance, 'valid_station_inds_temporal_colocation')

        # re-initialise and update metadata fields for read metadata
        init_metadata(self.read_instance)
        update_metadata_fields(self.read_instance)
        # if loading from a conf file then load metadata fields
        if self.read_instance.from_conf:
            metadata_conf(self.read_instance)

        # determine if have daily forecast active or not
        self.read_instance.daily_forecast = np.any([True for data_label in self.read_instance.data_labels if '-daily' in data_label])

        # determine if have combined forecast active or not
        self.read_instance.combined_forecast = np.any([True for data_label in self.read_instance.data_labels if '-combined' in data_label])

        # if have reading daily or combined forecast data, make original copy of data labels, models and plotting params, 
        # as will be modified later and may need restoring, for dashboard mode only
        if ((self.read_instance.daily_forecast) or (self.read_instance.combined_forecast)) & (self.read_instance.mode not in ['report', 'library']):
            self.read_instance.original_data_labels = copy.deepcopy(self.read_instance.data_labels)
            self.read_instance.original_data_labels_raw = copy.deepcopy(self.read_instance.data_labels_raw)
            self.read_instance.original_models = copy.deepcopy(self.read_instance.experiments)
            self.read_instance.original_plotting_params = copy.deepcopy(self.read_instance.plotting_params)
        
        # print basic information species
        self.read_instance.logger.info('\nOBSERVATIONS')
        for network in self.read_instance.networkspecies:
            self.read_instance.logger.info(f" - {network}")
        if self.read_instance.filter_species:
            self.read_instance.logger.info('OBSERVATIONS TO FILTER BY')
            for network_model, parameters in self.read_instance.filter_species.items():
                self.read_instance.logger.info(f" - {network_model} {parameters}")                       
        # print models after observations
        if self.read_instance.experiments:
            self.read_instance.logger.info("MODELS")
            mods_printed = []
            for model, alias in self.read_instance.experiments.items():
                if '-daily' in model:
                    model = '{}-daily'.format(model.split('-daily')[0])
                    alias = '{}-daily'.format(alias.split('-daily')[0])
                if '-combined' in model:
                    model = '{}-combined'.format(model.split('-combined')[0])
                    alias = '{}-combined'.format(alias.split('-combined')[0])
                if model not in mods_printed:
                    str_model = f" - {model}"
                    if self.read_instance.alias_flag: 
                        str_model += f" ({alias})"
                    self.read_instance.logger.info(str_model)
                    mods_printed.append(model)
        self.read_instance.logger.info("")

    def read_basic_metadata(self):     
        """
        Extracts unique basic station metadata and handles spatial colocation across multiple networkspecies (if set) in parallel.
        The basic metadata fields are: station_reference, station_name, longitude, latitude, measurement_altitude, 
        station_classification and area_classification.
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
                    self.read_instance.species.remove(networkspeci.split('|')[1])
                    msg = 'There is no available observational data for the network|species: {}. Dropping.'.format(networkspeci)
                    show_message(self.read_instance, msg)
                elif networkspeci in self.read_instance.filter_networkspecies:
                    self.read_instance.filter_networkspecies.remove(networkspeci)
                    if networkspeci.split('|')[1] in self.read_instance.filter_species:
                        self.read_instance.filter_species.remove(networkspeci.split('|')[1])
                    msg = 'There is no available observational data for the filter network|species: {}. Dropping.'.format(networkspeci)
                    show_message(self.read_instance, msg)
                continue
                
            # get basic metadata for speci - read in parallel
                
            # define arrays for storing speci metadata
            speci_station_references = []
            speci_station_names = []
            speci_station_longitudes = []
            speci_station_latitudes = []
            speci_station_measurement_altitudes = []

            # read metadata in parallel
            tuple_arguments = []
            for fname in relevant_files:
                tuple_arguments.append((fname, self.read_instance.reading_ghost, self.read_instance.logger))
            pool = multiprocessing.Pool(self.read_instance.n_cpus)
            returned_data = pool.map(read_netcdf_metadata, tuple_arguments)
            pool.close()
            pool.join()

            # unzip returned data
            for returned_data_per_month in returned_data:
                speci_station_references = np.append(speci_station_references, returned_data_per_month[0])
                speci_station_longitudes = np.append(speci_station_longitudes, returned_data_per_month[1])
                speci_station_latitudes = np.append(speci_station_latitudes, returned_data_per_month[2])
                speci_station_names = np.append(speci_station_names, returned_data_per_month[3])
                speci_station_measurement_altitudes = np.append(speci_station_measurement_altitudes, returned_data_per_month[4])

            # for non-GHOST data, check if station name and measurement altitude are available in files, and
            # get list of additional variables to read that maybe available in files (must be subset of GHOST variables)
            # also get non-GHOST measurement units
            if not self.read_instance.reading_ghost:

                # define base list of non-GHOST metadata variables to read
                self.read_instance.nonghost_metadata_vars_to_read = ['station_reference','longitude','latitude']
                basic_metadata_read = [speci_station_references, speci_station_longitudes, speci_station_latitudes]
                
                # check if station names are available in files
                if len(speci_station_names) > 0:
                    self.read_instance.nonghost_metadata_vars_to_read.append('station_name') 
                    basic_metadata_read.append(speci_station_names)
                    have_station_name = True
                else:
                    have_station_name = False

                # check if measurement altitudes are available in files
                if len(speci_station_measurement_altitudes) > 0:
                    self.read_instance.nonghost_metadata_vars_to_read.append('measurement_altitude') 
                    basic_metadata_read.append(speci_station_measurement_altitudes)
                    have_measurement_altitude = True
                else:
                    have_measurement_altitude = False

                # get non-GHOST measurement units (take from first relevant file)
                ncdf_root = Dataset(relevant_files[0])
                self.read_instance.nonghost_units[speci] = ncdf_root[speci].units

                # get list of nonghost variables to read (subset of GHOST variables)
                for ghost_metadata_var in self.read_instance.ghost_metadata_vars_to_read:
                    if (ghost_metadata_var in ncdf_root.variables) & (ghost_metadata_var not in self.read_instance.nonghost_metadata_vars_to_read):
                        self.read_instance.nonghost_metadata_vars_to_read.append(ghost_metadata_var) 

                # check if area classification is available in files (if not already added)
                if 'area_classification' not in self.read_instance.nonghost_metadata_vars_to_read:
                    if 'station_area' in ncdf_root.variables:
                        self.read_instance.nonghost_metadata_vars_to_read.append('area_classification') 

                # check if station classification is available in files (if not already added)
                if 'station_classification' not in self.read_instance.nonghost_metadata_vars_to_read:
                    if 'station_type' in ncdf_root.variables:
                        self.read_instance.nonghost_metadata_vars_to_read.append('station_classification') 

                # check if DOI (ACTRIS) is available in files (if not already added)
                if 'doi' not in self.read_instance.nonghost_metadata_vars_to_read:
                    if 'doi' in ncdf_root.variables:
                        self.read_instance.nonghost_metadata_vars_to_read.append('doi') 

                # check if national facility (ACTRIS) is available in files (if not already added)
                if 'actris_national_facility' not in self.read_instance.nonghost_metadata_vars_to_read:
                    if 'actris_national_facility' in ncdf_root.variables:
                        self.read_instance.nonghost_metadata_vars_to_read.append('actris_national_facility') 

                # close first relevant file
                ncdf_root.close()
            
            # for GHOST data, just set up some variables for checks
            else:
                basic_metadata_read = [speci_station_references, speci_station_longitudes, speci_station_latitudes, speci_station_names, speci_station_measurement_altitudes]
                have_station_name = True
                have_measurement_altitude = True

            # double check that lengths of read variables are the same
            if len(set(map(len, basic_metadata_read))) != 1:
                error = 'Error: Some metadata variables do not appear in all netCDF files. This should not be the case!'
                self.read_instance.logger.error(error)
                sys.exit(1) 

            # get unique station references and apply unique indices to the other variables        
            speci_station_references, station_unique_indices = np.unique(speci_station_references, return_index=True)
            self.read_instance.station_references[networkspeci] = speci_station_references
            self.read_instance.station_longitudes[networkspeci] = speci_station_longitudes[station_unique_indices]
            self.read_instance.station_latitudes[networkspeci] = speci_station_latitudes[station_unique_indices]
            if have_station_name:
                self.read_instance.station_names[networkspeci] = speci_station_names[station_unique_indices]
            if have_measurement_altitude:
                self.read_instance.station_measurement_altitudes[networkspeci] = speci_station_measurement_altitudes[station_unique_indices]

        # if have more than 1 networkspecies (including filter networkspecies), and spatial_colocation is active,
        # then spatially colocate stations across species
        if (len((self.read_instance.networkspecies + self.read_instance.filter_networkspecies)) > 1) & (self.read_instance.spatial_colocation):
            
            # get intersecting station indices across species (handle both GHOST and non-GHOST cases)
            self.sc = SpatialColocation(self.read_instance)

            # if have zero intersecting indices after spatial colocation, then deactivate spatial colocation
            if len(self.sc.intersecting_indices[self.sc.firstnetworkspeci]) == 0:
                msg = "spatial_colocation is set to False, as have 0 intersecting stations across species."
                show_message(self.read_instance, msg)
                self.read_instance.spatial_colocation = False
            # otherwise, iterate through networkspecies specific intersecting indices, reducing associated variables 
            # for specific intersecting indices per networkspecies
            else:        
                for ns, ns_intersects in self.sc.intersecting_indices.items():
                    self.read_instance.station_references[ns] = self.read_instance.station_references[ns][ns_intersects]
                    self.read_instance.station_longitudes[ns] = self.read_instance.station_longitudes[ns][ns_intersects]
                    self.read_instance.station_latitudes[ns] = self.read_instance.station_latitudes[ns][ns_intersects]
                    if ns in self.read_instance.station_measurement_altitudes:
                        self.read_instance.station_measurement_altitudes[ns] = self.read_instance.station_measurement_altitudes[ns][ns_intersects]
                    if ns in self.read_instance.station_names:
                        self.read_instance.station_names[ns] = self.read_instance.station_names[ns][ns_intersects]


    def check_forecast(self, yearmonths_to_read=None, data_labels=None, data_labels_raw=None, networkspecies=None):
        """
        Identifies the available forecast days for model data by inspecting NetCDF dimensions.

        Parameters
        ----------
        yearmonths_to_read : list, optional
            Target months to check; if None, it is calculated from the current instance date range.
        data_labels : list, optional
            Display labels for the data; defaults to instance data labels.
        data_labels_raw : list, optional
            Raw identifiers for the data; defaults to instance raw data labels.
        networkspecies : list, optional
            List of networkspecies pairs; defaults to instance networkspecies.
        """

        # set data labels if not defined
        if data_labels is None:
            data_labels = copy.deepcopy(self.read_instance.data_labels)

        # set data labels raw if not defined
        if data_labels_raw is None:
            data_labels_raw = copy.deepcopy(self.read_instance.data_labels_raw)

        # set networkspecies if not defined
        if networkspecies is None:
            networkspecies = copy.deepcopy(self.read_instance.networkspecies)

        # set forecast_days_per_data_label to be empty dictionary
        self.read_instance.forecast_days_per_data_label = {}

        # set yearmonths_to_read if not defined
        if yearmonths_to_read is None:
            # gather currently selected variables
            start_date = self.read_instance.le_start_date.text()
            end_date = self.read_instance.le_end_date.text()
            resolution = self.read_instance.selected_resolution
            
            # get active frequency code
            active_frequency_code = get_frequency_code(resolution)

            # get time array
            time_array = pd.date_range(start=datetime.datetime(int(start_date[:4]), int(start_date[4:6]), int(start_date[6:8])),
                                        end=datetime.datetime(int(end_date[:4]), int(end_date[4:6]), int(end_date[6:8])),
                                        freq=active_frequency_code)[:-1]

            # get yearmonths in data range (incomplete months are removed for monthly resolution)
            yearmonths = list(np.unique(['{}0{}'.format(dti.year,dti.month) if len(str(dti.month)) == 1 else '{}{}'.format(dti.year,dti.month) \
                            for dti in time_array]))

            # get yearmonths to read
            yearmonths_to_read = get_yearmonths_to_read(yearmonths, start_date, end_date, resolution)

        # iterate through network, speci pairs
        for networkspeci in networkspecies:

            # get indivudual network and species strings
            network = networkspeci.split('|')[0]
            speci = networkspeci.split('|')[1]

            # add dictionary of files to read
            files_to_read = {}

            # add dictionary for available and selected forecast days per data label for this networkspeci
            if networkspeci not in self.read_instance.forecast_days_per_data_label:
                self.read_instance.forecast_days_per_data_label[networkspeci] = {}

            # iterate through data labels
            for data_label, data_label_raw in zip(data_labels, data_labels_raw):

                # only check model data labels
                if data_label != self.read_instance.observations_data_label:

                    if '/' in network:
                        file_root = \
                            '%s/%s/%s/%s/%s/%s/%s_' % (self.read_instance.mod_root, self.read_instance.ghost_version, 
                                                        data_label_raw, self.read_instance.resolution, speci, 
                                                        network.replace('/', '-'), speci)
                    else:
                        file_root = \
                            '%s/%s/%s/%s/%s/%s/%s_' % (self.read_instance.mod_root, self.read_instance.ghost_version, 
                                                        data_label_raw, self.read_instance.resolution, speci, network, speci)

                    try:
                        available_yearmonths = self.read_instance.available_model_data[network][self.read_instance.resolution][speci][data_label_raw]
                    except KeyError:
                        continue

                    # get intersection of yearmonths_to_read and available_yearmonths
                    yearmonths_to_read_intersect = list(set(yearmonths_to_read) & set(available_yearmonths))
                    files_to_read[data_label] = sorted([file_root+str(yyyymm)+'.nc' for yyyymm in yearmonths_to_read_intersect])[0]

            # read forecast dimension for each model
            # if only one data label, then read forecast dimension directly
            if len(files_to_read) == 0:
                continue
            elif len(files_to_read) == 1:
                returned_data = check_forecast_dimension(files_to_read[list(files_to_read.keys())[0]])
                returned_data = [returned_data]
            # if more than one data label, then read forecast dimension in parallel
            elif len(files_to_read) > 1:
                tuple_arguments = []
                for data_label in files_to_read:
                    tuple_arguments.append((files_to_read[data_label]))
                pool = multiprocessing.Pool(self.read_instance.n_cpus)
                returned_data = pool.map(check_forecast_dimension, tuple_arguments)
                pool.close()
                pool.join()

            # unzip returned data
            n_forecast_days = np.array(returned_data, dtype=np.int32)

            # set n forecast days for each data label
            for data_label, n_forecast_day in zip(list(files_to_read.keys()), n_forecast_days):

                # set n forecast days per data label
                self.read_instance.forecast_days_per_data_label[networkspeci][data_label] = n_forecast_day

 
    def update_forecast_indices(self, data_labels=None, data_labels_raw=None, 
                                selected_data_labels=None, selected_data_labels_raw=None, 
                                networkspecies=None, init=False):
        """
        Updates forecast-related indices, labels, and model configurations based on current settings and selections.

        Parameters
        ----------
        data_labels : list, optional
            Display labels for the data; defaults to instance data labels.
        data_labels_raw : list, optional
            Raw identifiers for the data; defaults to instance raw data labels.
        selected_data_labels : list, optional
            Subset of display labels currently selected; defaults to all data labels.
        selected_data_labels_raw : list, optional
            Subset of raw identifiers currently selected; defaults to all raw labels.
        networkspecies : list, optional
            List of networkspecies pairs; defaults to instance networkspecies.
        init : bool, optional
            Flag indicating if this is an initialisation call to reset menu structures.

        Returns
        -------
        new_data_labels : list
            Updated list of data labels including generated forecast suffixes.
        new_data_labels_raw : list
            Updated list of raw data labels including generated forecast suffixes.
        new_models : dict
            Mapping of updated raw data labels to their display counterparts.
        """

        # set data labels if not defined
        if data_labels is None:
            data_labels = copy.deepcopy(self.read_instance.data_labels)

        # set selected data labels if not defined
        if selected_data_labels is None:
            selected_data_labels = copy.deepcopy(data_labels)

        # set data labels raw if not defined
        if data_labels_raw is None:
            data_labels_raw = copy.deepcopy(self.read_instance.data_labels_raw)

        # set selected data labels raw if not defined
        if selected_data_labels_raw is None:
            selected_data_labels_raw = copy.deepcopy(data_labels_raw)

        # set networkspecies if not defined
        if networkspecies is None:
            networkspecies = copy.deepcopy(self.read_instance.networkspecies)

        # list to store the forecast days that will be processed
        wanted_forecast_days = []
        if len(self.read_instance.forecast) > 0:
            for fct_ii, fct in enumerate(self.read_instance.forecast):
                # on first pass get forecast_type (e.g., 'day', 'daily', 'combined')
                if fct_ii == 0:
                    if 'day' in fct:
                        forecast_type = 'day'
                    elif 'combined' in fct:
                        forecast_type = 'combined'
                    elif 'daily' in fct:
                        forecast_type = 'daily' 

                # remove forecast_type from forecast variable (e.g., "day3" -> "3")
                fct = fct.split(forecast_type)[-1]
                # if no forecast day then set forecast type (e.g., just "day")
                if fct == '':
                    wanted_forecast_days.append(forecast_type)
                # if colon present, get a range of forecast days (e.g., "1:5")
                elif ':' in fct:
                    colon_split = fct.split(':')
                    start_fct_day = int(colon_split[0])
                    end_fct_day = int(colon_split[1])+1
                    # do not allow forecast day to start before 1
                    if start_fct_day < 1:
                        start_fct_day = 1
                    # create array of days (e.g., [1,2,3,4,5])
                    fct_days = np.arange(start_fct_day,end_fct_day,dtype=np.int32)
                    for fct_day in fct_days:
                        wanted_forecast_days.append(fct_day)
                # else it’s a single forecast day (e.g., "3")
                else:
                    wanted_forecast_days.append(int(fct))

        # initialise dictionary to store forecast indices per data label
        self.read_instance.forecast_indices_per_data_label = {}

        # prepare unique lists to track which labels are added/removed across network species
        unique_data_labels_to_add = []
        unique_data_labels_raw_to_add = []
        unique_data_labels_to_remove = []
        unique_data_labels_raw_to_remove = []

        # iterate over each network species
        for networkspeci in networkspecies:

            # temporary lists for each species
            data_labels_to_remove = []
            data_labels_raw_to_remove = []
            data_labels_to_add = []
            data_labels_raw_to_add = []

            # reset models pop-up menu options if in dashboard mode, and are initialising
            if (self.read_instance.mode not in ['report', 'library']) & (init):
                self.read_instance.models_menu['models']['forecast'] = {}
                self.read_instance.models_menu['models']['forecast_days'] = {}
            # otherwise reset selected and disabled forecast variable and day options (to ensure data labels that are no longer selected are cleaned)
            elif (self.read_instance.mode not in ['report', 'library']) & (not init): 
                for data_label in self.read_instance.models_menu['models']['forecast']:
                    self.read_instance.models_menu['models']['forecast'][data_label][1] = []
                    self.read_instance.models_menu['models']['forecast'][data_label][2] = []
                    self.read_instance.models_menu['models']['forecast_days'][data_label][1] = []

            # initialise dictionary for this network species
            self.read_instance.forecast_indices_per_data_label[networkspeci] = {}
            for data_label, data_label_raw in zip(data_labels, data_labels_raw):
                # skip the observation data label
                if data_label == self.read_instance.observations_data_label:
                     continue
                self.read_instance.forecast_indices_per_data_label[networkspeci][data_label] = {}
                # if do not have forecast days for data label, data does not exist for data label, so continue
                if data_label not in self.read_instance.forecast_days_per_data_label[networkspeci]:
                    continue
                # get number of forecast days for this label
                n_forecast_days = self.read_instance.forecast_days_per_data_label[networkspeci][data_label]

                # reset models pop-up menu options if in dashboard mode for specific data label (as now know will reset it)
                if self.read_instance.mode not in ['report', 'library']:
                    if data_label in self.read_instance.models_menu['models']['forecast']:
                        del self.read_instance.models_menu['models']['forecast'][data_label]
                        del self.read_instance.models_menu['models']['forecast_days'][data_label]

                # if no forecast days available, just update the menu with empty entry
                if n_forecast_days == 0:
                    self.update_forecast_menu(networkspeci, data_label, '', n_forecast_days, '')
                else:
                    # if no specific forecast requested, still update empty menu entry
                    if len(wanted_forecast_days) == 0:
                        self.update_forecast_menu(networkspeci, data_label, '', n_forecast_days, '')
                    else:
                        # mark the data label for removal if it was selected (will be replaced by forecast versions)
                        if data_label in selected_data_labels:
                            data_labels_to_remove.append(data_label)
                            data_labels_raw_to_remove.append(data_label_raw)

                        # if forecast type covers all days (combined/daily/day)
                        if wanted_forecast_days[0] in ['combined', 'daily', 'day']:
                            for wanted_forecast_index in range(n_forecast_days):
                                # create new label names with forecast day index
                                new_data_label = '{}-{}{}'.format(data_label, forecast_type, wanted_forecast_index+1)
                                new_data_label_raw = '{}-{}{}'.format(data_label_raw, forecast_type, wanted_forecast_index+1)
                                # map new labels to forecast indices
                                self.read_instance.forecast_indices_per_data_label[networkspeci][data_label][new_data_label] = wanted_forecast_index
                                # update the forecast menu with the new labels
                                self.update_forecast_menu(networkspeci, data_label, new_data_label, n_forecast_days, wanted_forecast_days[0])
                                # if label was selected, add the new one to be displayed
                                if data_label in selected_data_labels:
                                    data_labels_to_add.append(new_data_label)
                                    data_labels_raw_to_add.append(new_data_label_raw)
                                
                        else:
                            # otherwise, handle explicit list of wanted forecast days
                            for wanted_forecast_day in wanted_forecast_days:
                                if wanted_forecast_day <= n_forecast_days:
                                    # create new forecast label for that specific day
                                    new_data_label = '{}-{}{}'.format(data_label, forecast_type, wanted_forecast_day)
                                    new_data_label_raw = '{}-{}{}'.format(data_label_raw, forecast_type, wanted_forecast_day)
                                    # store forecast index (zero-based)
                                    self.read_instance.forecast_indices_per_data_label[networkspeci][data_label][new_data_label] = wanted_forecast_day-1
                                    # update forecast menu
                                    self.update_forecast_menu(networkspeci, data_label, new_data_label, n_forecast_days, wanted_forecast_day)
                                    # if label was selected, add new label to lists
                                    if data_label in selected_data_labels:
                                        data_labels_to_add.append(new_data_label)
                                        data_labels_raw_to_add.append(new_data_label_raw)
                                    
            # ensure unique tracking of labels to add/remove across all species
            for data_label_to_remove in data_labels_to_remove:
                if data_label_to_remove not in unique_data_labels_to_remove:
                    unique_data_labels_to_remove.append(data_label_to_remove)
            for data_label_raw_to_remove in data_labels_raw_to_remove:
                if data_label_raw_to_remove not in unique_data_labels_raw_to_remove:
                    unique_data_labels_raw_to_remove.append(data_label_raw_to_remove)
            for data_label_to_add in data_labels_to_add:
                if data_label_to_add not in unique_data_labels_to_add:
                    unique_data_labels_to_add.append(data_label_to_add)
            for data_label_raw_to_add in data_labels_raw_to_add:
                if data_label_raw_to_add not in unique_data_labels_raw_to_add:
                    unique_data_labels_raw_to_add.append(data_label_raw_to_add)

        # update data labels and models (replace removed labels with new forecasted ones)
        if (len(unique_data_labels_to_remove) > 0) & (len(unique_data_labels_to_add) > 0):

            new_data_labels = []
            new_data_labels_raw = []
            new_models = {}

            # iterate through selected data labels and build updated lists
            for data_label, data_label_raw in zip(selected_data_labels, selected_data_labels_raw):
                # keep original labels if not removed
                if data_label not in unique_data_labels_to_remove:
                    new_data_labels.append(data_label)
                    new_data_labels_raw.append(data_label_raw)
                    if data_label != self.read_instance.observations_data_label:
                        new_models[data_label_raw] = data_label
                else:
                    # replace removed labels with corresponding forecast labels
                    for data_label_to_add, data_label_raw_to_add in zip(unique_data_labels_to_add, unique_data_labels_raw_to_add):
                        # find base label before forecast suffix
                        base_data_label_to_add = data_label_to_add.split('-day')[0].split('-daily')[0].split('-combined')[0]
                        if base_data_label_to_add == data_label:
                            new_data_labels.append(data_label_to_add)
                            new_data_labels_raw.append(data_label_raw_to_add)
                            new_models[data_label_raw_to_add] = data_label_to_add
        else:
            # if no changes were made, keep existing selections
            new_data_labels = copy.deepcopy(selected_data_labels)
            new_data_labels_raw = copy.deepcopy(selected_data_labels_raw)
            new_models = {} 
            for data_label, data_label_raw in zip(new_data_labels, new_data_labels_raw):
                if data_label != self.read_instance.observations_data_label:
                    new_models[data_label_raw] = data_label

        # return updated data label lists and model mapping
        return new_data_labels, new_data_labels_raw, new_models

    def update_forecast_menu(self, networkspeci, data_label, new_data_label, n_forecast_days, wanted_forecast_day):
        """
        Update the forecast menu options and selection states for the dashboard interface.

        Parameters
        ----------
        networkspeci : str
            The specific networkspecies string being processed.
        data_label : str
            The original display label of the model dataset.
        new_data_label : str
            The newly generated label containing the forecast suffix.
        n_forecast_days : int
            Total number of forecast days available in the dataset.
        wanted_forecast_day : str or int
            The specific forecast day or type requested (e.g., 'combined', 'daily', or an integer).
        """

        # Only proceed if are in dashboard mode
        if self.read_instance.mode not in ['report', 'library']:
            
            # Case 1: No forecast days are available
            if n_forecast_days == 0:
                
                # Initialise all forecast-related lists as empty
                available_forecast_vars = []
                available_forecast_day_vars = []
                selected_forecast_vars = []
                selected_forecast_day_var = '' 
                disabled_forecast_vars = []

            else:
                # Case 2: Forecast days exist → build available options

                # Types of forecast data that can be selected
                available_forecast_vars = ['combined', 'daily', 'day']

                # List all available forecast day labels (e.g. 'day 1', 'day 2', ...)
                available_forecast_day_vars = ['day {}'.format(forecast_day+1) for forecast_day in range(n_forecast_days)]

                # If no specific forecast day is requested
                if wanted_forecast_day == '':
                    selected_forecast_vars = []
                    selected_forecast_day_var = '' 
                    disabled_forecast_vars = []

                else:
                    # Extract the forecast type from the new data label (e.g., '-combined', '-daily', '-day...')
                    forecast_extension = new_data_label.split('-')[-1]

                    # Determine which forecast variable is selected and which others should be disabled
                    if 'combined' in forecast_extension:
                        selected_forecast_vars = ['combined']
                        disabled_forecast_vars = ['daily', 'day']
                    elif 'daily' in forecast_extension:
                        selected_forecast_vars = ['daily']
                        disabled_forecast_vars = ['combined', 'day']
                    elif 'day' in forecast_extension:
                        selected_forecast_vars = ['day']
                        disabled_forecast_vars = ['combined', 'daily']

                    # Get the index of the selected forecast day from a mapping of forecast indices
                    selected_forecast_index = self.read_instance.forecast_indices_per_data_label[networkspeci][data_label][new_data_label]

                    # Convert index to a display label (e.g. index 0 → 'day 1')
                    selected_forecast_day_var = 'day {}'.format(selected_forecast_index + 1)

            # If this data_label has no existing entry in the forecast menu, initialize it
            if data_label not in self.read_instance.models_menu['models']['forecast']:
                # Store available, selected, and disabled forecast variable options
                self.read_instance.models_menu['models']['forecast'][data_label] = [
                    available_forecast_vars,
                    selected_forecast_vars,
                    disabled_forecast_vars
                ]

                # Initialise available forecast day options and empty selected list
                self.read_instance.models_menu['models']['forecast_days'][data_label] = [
                    available_forecast_day_vars,
                    []
                ]

            # If a specific forecast day is selected, add it to the selected forecast days list
            if selected_forecast_day_var != '':
                self.read_instance.models_menu['models']['forecast_days'][data_label][1].append(selected_forecast_day_var)


    def read_data(self, yearmonths_to_read, data_labels):
        """
        Manages the parallel reading of observational and model NetCDF data into shared memory arrays.

        Parameters
        ----------
        yearmonths_to_read : list
            List of year-month strings (YYYYMM) identifying the files to be read.
        data_labels : list
            List of data labels representing the datasets to be loaded into memory.
        """

        # create arrays to share across processes (for parallel multiprocessing use)
        # this only works for numerical dtypes, i.e. not strings
        timestamp_array_shared = multiprocessing.RawArray(ctypes.c_int64, len(self.read_instance.timestamp_array))
        if (self.read_instance.reading_ghost or self.read_instance.network[0] == 'actris/actris') & (self.read_instance.observations_data_label in data_labels):
            flags_shared = multiprocessing.RawArray(ctypes.c_uint8, len(self.read_instance.flags))
        else:
            flags_shared = None
        # fill arrays
        timestamp_array_shared[:] = self.read_instance.timestamp_array
        if (self.read_instance.reading_ghost or self.read_instance.network[0] == 'actris/actris') & (self.read_instance.observations_data_label in data_labels):
            flags_shared[:] = self.read_instance.flags

        # create dictionary for saving files to read
        self.read_instance.files_to_read = {}

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
            self.read_instance.files_to_read[networkspeci] = {}

            # iterate through data labels
            for data_label in data_labels:

                # get raw data label (non-alias)
                data_label_raw = self.read_instance.data_labels_raw[self.read_instance.data_labels.index(data_label)]

                # get base data label and data label raw (i.e. without -dayN or -daily or -combined suffix for forecast data)
                base_data_label = data_label.split('-day')[0].split('-daily')[0].split('-combined')[0]
                base_data_label_raw = data_label_raw.split('-day')[0].split('-daily')[0].split('-combined')[0]

                # if already have base data label in files_to_read, then continue
                if base_data_label in self.read_instance.files_to_read[networkspeci]:
                    continue

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

                # models 
                else:
                    # if are reading filter species continue to next data_label
                    if filter_read:
                        continue 

                    elif '/' in network:
                        file_root = \
                            '%s/%s/%s/%s/%s/%s/%s_' % (self.read_instance.mod_root, self.read_instance.ghost_version, 
                                                        base_data_label_raw, self.read_instance.resolution, speci, 
                                                        network.replace('/', '-'), speci)
                    else:
                        file_root = \
                            '%s/%s/%s/%s/%s/%s/%s_' % (self.read_instance.mod_root, self.read_instance.ghost_version, 
                                                       base_data_label_raw, self.read_instance.resolution, speci, network, speci)
                    try:
                        available_yearmonths = self.read_instance.available_model_data[network][self.read_instance.resolution][speci][base_data_label_raw]
                    except KeyError:
                        continue

                # get intersection of yearmonths_to_read and available_yearmonths
                yearmonths_to_read_intersect = list(set(yearmonths_to_read) & set(available_yearmonths))
                self.read_instance.files_to_read[networkspeci][base_data_label] = sorted([file_root+str(yyyymm)+'.nc' for yyyymm in yearmonths_to_read_intersect])
                
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
            if (self.read_instance.reading_ghost or self.read_instance.network[0] == 'actris/actris') & (self.read_instance.observations_data_label in data_labels):
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
            
            if (self.read_instance.reading_ghost or self.read_instance.network[0] == 'actris/actris') & (self.read_instance.observations_data_label in data_labels):      
                qa_shared[:] = self.read_instance.qa_per_species[speci]
                if (self.read_instance.reading_ghost):      
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
                                     'metadata_dtype', 'metadata_vars_to_read', 'default_qa_active', 'filter_read', 
                                     'network', 'forecast_indices']
            tuple_arguments = []

            # iterate through base data labels
            for base_data_label in self.read_instance.files_to_read[networkspeci]:
                # set forecast indices
                # no forecast indices if loading observations
                if base_data_label == self.read_instance.observations_data_label:
                    forecast_indices = np.array([], dtype=np.int32)
                # no forecast indices if loading non-forecast model
                elif len(self.read_instance.forecast_indices_per_data_label[networkspeci][base_data_label]) == 0:
                    forecast_indices = np.array([], dtype=np.int32)
                # get forecast indices if loading forecast model
                else:
                    forecast_indices = np.array([self.read_instance.forecast_indices_per_data_label[networkspeci][base_data_label][data_label] 
                                                 for data_label in data_labels if data_label in self.read_instance.forecast_indices_per_data_label[networkspeci][base_data_label]], dtype=np.int32)
                    
                for fname in self.read_instance.files_to_read[networkspeci][base_data_label]:
                    tuple_arguments.append((fname, self.read_instance.station_references[networkspeci], 
                                            self.read_instance.station_names[networkspeci], speci, 
                                            self.read_instance.observations_data_label, 
                                            base_data_label, data_labels, 
                                            self.read_instance.reading_ghost, 
                                            self.read_instance.ghost_data_vars_to_read, 
                                            self.read_instance.metadata_dtype, 
                                            self.read_instance.metadata_vars_to_read,
                                            self.read_instance.logger,
                                            default_qa_active, filter_read, 
                                            network, forecast_indices))

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
                    error = 'Error: The observation and model arrays for {} are empty.'.format(networkspeci)
            
                elif np.isin(data_array.flatten(), [-9999.0, np.nan]).all():
                    error = 'Error: All observation and model arrays for {} are void.'.format(networkspeci)

                self.read_instance.logger.error(error)
                sys.exit(1) 