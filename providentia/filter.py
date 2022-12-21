from .calculate import Stats
from .configuration import split_options

import copy
import numpy as np
import pandas as pd

from PyQt5 import QtWidgets
from providentia import aux

class DataFilter:
    """ Class that filters observational/experiment data into memory as required. """

    def __init__(self, read_instance):
        self.read_instance = read_instance

        # get indices of some data variables
        self.obs_index = self.read_instance.data_labels.index('observations')
        if self.read_instance.reading_ghost:
            if self.read_instance.resolution != 'daily' and self.read_instance.resolution != 'monthly':
                self.day_night_index = self.read_instance.ghost_data_vars_to_read.index('day_night_code')
            if self.read_instance.resolution != 'monthly':
                self.weekday_weekend_index = self.read_instance.ghost_data_vars_to_read.index('weekday_weekend_code')
            self.season_index = self.read_instance.ghost_data_vars_to_read.index('season_code')

        # apply filtering
        self.filter_all()

    def filter_all(self):
        """ Call methods to start filtering. """

        self.reset_data_filter()
        self.filter_by_species()
        self.filter_data_limits()
        self.filter_by_period()
        self.filter_by_data_availability()
        self.filter_by_metadata()
        self.temporally_colocate_data()
        self.get_valid_stations_after_filtering()

    def reset_data_filter(self):
        """ Resets data arrays to be un-filtered"""

        self.read_instance.data_in_memory_filtered = copy.deepcopy(self.read_instance.data_in_memory)
        self.read_instance.temporal_colocation_nans = {}
        self.read_instance.valid_station_inds = {}
        self.read_instance.valid_station_inds_temporal_colocation = {}

    def filter_by_species(self):
        """ Function which filters read species by other species.
            For N other species a lower and upper limit are set. 
            Where values for each species are outside of these ranges,
            then impose NaNs upon all read species in memory.
            Only filter if spatial colocation is True.
        """

        # filter all read species by set species ranges
        if (self.read_instance.filter_species) and (self.read_instance.spatial_colocation):

            # initialise array to set where temporally to filter species
            # initialse being all False, set as True where data is outside given bounds for species
            inds_to_filter = np.full(self.read_instance.data_in_memory_filtered[self.read_instance.networkspecies[0]][self.obs_index,:,:].shape, False)    

            # iterate through all species to filter by
            for filter_networkspeci, speci_limits in self.read_instance.filter_species.items():
                
                # get lower and upper limits for species
                lower_limit = speci_limits[0]
                upper_limit = speci_limits[1]
                filter_species_fill_value = speci_limits[2]

                # get where data is outside bounds
                if filter_networkspeci in self.read_instance.networkspecies:
                    invalid_inds_per_species = np.logical_or(self.read_instance.data_in_memory_filtered[filter_networkspeci][self.obs_index, :,:] < lower_limit,
                                                             self.read_instance.data_in_memory_filtered[filter_networkspeci][self.obs_index, :,:] > upper_limit)
                else:
                    invalid_inds_per_species = np.logical_or(self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] < lower_limit,
                                                             self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] > upper_limit)

                # update inds_to_filter array, making True all instances where have data outside bounds
                inds_to_filter = np.any([inds_to_filter, invalid_inds_per_species], axis=0)

                # set all inds to filter as NaN for all networkspecies in memory
                for networkspeci in self.read_instance.networkspecies:
                    self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index, inds_to_filter] = filter_species_fill_value      

    def filter_data_limits(self):
        """ Filter out (set to NaN) data which exceed the lower/upper limits. """

        # iterate through networkspecies  
        for networkspeci in self.read_instance.networkspecies:

            # get speci str
            speci = networkspeci.split('|')[1]

            # get lower/upper data bounds
            if self.read_instance.offline:
                lower_bound = self.read_instance.lower_bound[speci]
                upper_bound = self.read_instance.upper_bound[speci]
            else:
                lower_bound = self.read_instance.le_minimum_value.text()
                upper_bound = self.read_instance.le_maximum_value.text()

            # check selected lower/upper bounds are numbers
            try:
                lower_bound = np.float32(lower_bound)
                upper_bound = np.float32(upper_bound)
            # if any of the fields are not numbers, return from function
            except ValueError:
                print("Warning: Data limit fields must be numeric")
                return

            # filter all observational/experiment data out of bounds of lower/upper limits
            inds_out_of_bounds = np.logical_or(self.read_instance.data_in_memory_filtered[networkspeci][:,:,:] < lower_bound,
                                               self.read_instance.data_in_memory_filtered[networkspeci][:,:,:] > upper_bound)
            self.read_instance.data_in_memory_filtered[networkspeci][inds_out_of_bounds] = np.NaN

    def filter_by_period(self):
        """ Filter data for selected periods (keeping or removing data, as defined). """

        keeps, removes = [], []
        if self.read_instance.offline:
            if hasattr(self.read_instance, 'period'):
                keeps, removes = split_options(self.read_instance.period)
        else:
            keeps = self.read_instance.period_menu['checkboxes']['keep_selected']
            removes = self.read_instance.period_menu['checkboxes']['remove_selected']

        # filter/limit data for periods selected
        if len(keeps) > 0:
            day_night_codes_to_keep = []
            if 'Daytime' in keeps:
                day_night_codes_to_keep.append(0)
            if 'Nighttime' in keeps:
                day_night_codes_to_keep.append(1)
            if len(day_night_codes_to_keep) == 1:
                if (self.read_instance.resolution != 'daily') & (self.read_instance.resolution != 'monthly'):
                    # iterate through network / species  
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.day_night_index,:,:], day_night_codes_to_keep, invert=True)
                        self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN

            weekday_weekend_codes_to_keep = []
            if 'Weekday' in keeps:
                weekday_weekend_codes_to_keep.append(0)
            if 'Weekend' in keeps:
                weekday_weekend_codes_to_keep.append(1)
            if len(weekday_weekend_codes_to_keep) == 1:
                if self.read_instance.resolution != 'monthly':
                    # iterate through network / species  
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.weekday_weekend_index,:,:], weekday_weekend_codes_to_keep, invert=True)
                        self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN

            season_codes_to_keep = []
            if 'Spring' in keeps:
                season_codes_to_keep.append(0)
            if 'Summer' in keeps:
                season_codes_to_keep.append(1)
            if 'Autumn' in keeps:
                season_codes_to_keep.append(2)
            if 'Winter' in keeps:
                season_codes_to_keep.append(3)
            if (len(season_codes_to_keep) > 0) & (len(season_codes_to_keep) < 4):
                # iterate through network / species  
                for networkspeci in self.read_instance.networkspecies:
                    inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.season_index,:,:], season_codes_to_keep, invert=True)
                    self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN

        if len(removes) > 0:
            day_night_codes_to_remove = []
            if 'Daytime' in removes:
                day_night_codes_to_remove.append(0)
            if 'Nighttime' in removes:
                day_night_codes_to_remove.append(1)
            if len(day_night_codes_to_remove) > 0:
                if (self.read_instance.resolution != 'daily') & (self.read_instance.resolution != 'monthly'):
                    # iterate through network / species  
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.day_night_index,:,:], day_night_codes_to_remove)
                        self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN

            weekday_weekend_codes_to_remove = []
            if 'Weekday' in removes:
                weekday_weekend_codes_to_remove.append(0)
            if 'Weekend' in removes:
                weekday_weekend_codes_to_remove.append(1)
            if len(weekday_weekend_codes_to_remove) > 0:
                if self.read_instance.resolution != 'monthly':
                    # iterate through network / species  
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.weekday_weekend_index,:,:], weekday_weekend_codes_to_remove)
                        self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN

            season_codes_to_remove = []
            if 'Spring' in removes:
                season_codes_to_remove.append(0)
            if 'Summer' in removes:
                season_codes_to_remove.append(1)
            if 'Autumn' in removes:
                season_codes_to_remove.append(2)
            if 'Winter' in removes:
                season_codes_to_remove.append(3)
            if len(season_codes_to_remove) > 0:
                # iterate through network / species  
                for networkspeci in self.read_instance.networkspecies:
                    inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.season_index,:,:], season_codes_to_remove)
                    self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN

    def filter_by_data_availability(self):
        """ Function which filters data by selected data availability variables. """

        # get set variables names representing percentage data availability (native and non-native)
        active_data_availablity_vars = self.read_instance.representativity_menu['rangeboxes']['map_vars']

        try:
            data_availability_lower_bounds = []
            for var_ii, var in enumerate(active_data_availablity_vars):
                data_availability_lower_bounds.append(
                    np.float32(self.read_instance.representativity_menu['rangeboxes']['current_lower'][var_ii]))
        # if any of the fields are not numbers, return from function
        except ValueError:
            print("Warning: Data availability fields must be numeric")
            return

        # filter observations by native percentage data availability variables (only GHOST data)
        if self.read_instance.reading_ghost:
            for var_ii, var in enumerate(active_data_availablity_vars):
                if 'native' in var:
                    var_index = self.read_instance.ghost_data_vars_to_read.index(var)
                    
                    # iterate through network / species  
                    for networkspeci in self.read_instance.networkspecies:

                        # max gap variable?
                        if 'max_gap' in var:
                            # bound is < 100?:
                            if data_availability_lower_bounds[var_ii] < 100:
                                inds_to_screen = self.read_instance.ghost_data_in_memory[networkspeci][var_index,:,:] > data_availability_lower_bounds[var_ii]
                                self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index, inds_to_screen] = np.NaN
                        # data representativity variable?
                        else:
                            # bound is > 0?
                            if data_availability_lower_bounds[var_ii] > 0:
                                inds_to_screen = self.read_instance.ghost_data_in_memory[networkspeci][var_index,:,:] < data_availability_lower_bounds[var_ii]
                                self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index, inds_to_screen] = np.NaN

        # filter observations and experiment data by non-native percentage data availability variables 
        # (calculated on the fly)
        for var_ii, var in enumerate(active_data_availablity_vars):
            if 'native' not in var:
                # max gap variable?
                if 'max_gap' in var:
                    # bound is == 100?
                    if data_availability_lower_bounds[var_ii] == 100:
                        continue
                # data representativity variable?
                else:
                    # bound == 0?
                    if data_availability_lower_bounds[var_ii] == 0:
                        continue

                # get period associate with variable
                period = var.split('_')[0]
                period_inds = np.arange(len(self.read_instance.time_array))
                # daily variable?
                if period == 'daily':
                    period_inds_split = np.array_split(period_inds,
                                                       [24 * i for i in range(1, int(np.ceil(len(period_inds) / 24)))])
                # monthly variable?
                elif period == 'monthly':
                    period_inds_split = np.array_split(period_inds, np.cumsum(self.read_instance.N_inds_per_yearmonth))
                # whole record variable?
                else:
                    period_inds_split = [period_inds]

                # iterate through indices associated with periodic chunks for current period
                for period_inds in period_inds_split:
                    if len(period_inds) > 0:

                        # iterate through networkspecies  
                        for networkspeci in self.read_instance.networkspecies:

                            # max gap variable?
                            if 'max_gap' in var:
                                max_gap_percent = Stats.max_repeated_nans_fraction(
                                    self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,period_inds])
                                inds_to_screen = np.where(max_gap_percent > data_availability_lower_bounds[var_ii])[0]
                                self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,inds_to_screen[:,np.newaxis],period_inds[np.newaxis,:]] = np.NaN

                            # data representativity variable?
                            else:
                                data_availability_percent = Stats.calculate_data_avail_fraction(
                                    self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,period_inds])
                                inds_to_screen = np.where(data_availability_percent < data_availability_lower_bounds[var_ii])[0]
                                self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,inds_to_screen[:,np.newaxis],period_inds[np.newaxis,:]] = np.NaN

    def filter_by_metadata(self):
        """ Filter data by selected metadata. """

        # validate fields before filtering
        if not self.validate_values():
            return

        # iterate through metadata in memory
        for meta_var in self.read_instance.metadata_vars_to_read:
            
            if meta_var == 'lat':
                meta_var = 'latitude'
            elif meta_var == 'lon':
                meta_var = 'longitude'

            metadata_type = self.read_instance.standard_metadata[meta_var]['metadata_type']
            metadata_data_type = self.read_instance.standard_metadata[meta_var]['data_type']

            # handle non-numeric metadata
            if metadata_data_type == np.object:

                # iterate through networkspecies  
                for networkspeci in self.read_instance.networkspecies:

                    # if any of the keep checkboxes are selected, filter out data by fields that have not been selected
                    current_keep = self.read_instance.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected']
                    if len(current_keep) > 0:
                        invalid_keep = np.repeat(
                            np.isin(self.read_instance.metadata_in_memory[networkspeci][meta_var][:, :], current_keep, invert=True), 
                                    self.read_instance.N_inds_per_yearmonth, axis=1)
                        self.read_instance.data_in_memory_filtered[networkspeci][:,invalid_keep] = np.NaN
                
                    # if any of the remove checkboxes have been selected, filter out data by these selected fields
                    current_remove = self.read_instance.metadata_menu[metadata_type][meta_var]['checkboxes'][
                        'remove_selected']
                    if len(current_remove) > 0:
                        invalid_remove = np.repeat(
                            np.isin(self.read_instance.metadata_in_memory[networkspeci][meta_var][:, :], current_remove),
                                    self.read_instance.N_inds_per_yearmonth, axis=1)
                        self.read_instance.data_in_memory_filtered[networkspeci][:,invalid_remove] = np.NaN
            
            # handle numeric metadata
            else:
                meta_var_index = self.read_instance.metadata_menu[metadata_type]['rangeboxes']['labels'].index(meta_var)
                current_lower = np.float32(
                    self.read_instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index])
                current_upper = np.float32(
                    self.read_instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index])
                
                # get array with selected filters (those with the check Apply on)
                current_apply = self.read_instance.metadata_menu[metadata_type]['rangeboxes']['apply_selected']

                if len(current_apply) > 0:
                    # apply bounds and remove nans if variable has been selected
                    if meta_var in current_apply:

                        # iterate through networkspecies  
                        for networkspeci in self.read_instance.networkspecies:

                            # if current lower value is non-NaN, then filter out data with metadata < current lower value
                            if not pd.isnull(current_lower):
                                lower_default = np.float32(
                                    self.read_instance.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index])
                                if current_lower >= lower_default: 
                                    invalid_below = np.repeat(self.read_instance.metadata_in_memory[networkspeci][meta_var][:, :] < current_lower, 
                                                              self.read_instance.N_inds_per_yearmonth, axis=1)
                                    self.read_instance.data_in_memory_filtered[networkspeci][:,invalid_below] = np.NaN

                            # if current upper < than the maximum extent, then filter out
                            # data with metadata > current upper value (if this is numeric)
                            # if current upper value is non-NaN, then filter out data with metadata > current upper value
                            if not pd.isnull(current_upper):
                                upper_default = np.float32(
                                    self.read_instance.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index])
                                if current_upper <= upper_default: 
                                    invalid_above = np.repeat(self.read_instance.metadata_in_memory[networkspeci][meta_var][:, :] > current_upper, 
                                                              self.read_instance.N_inds_per_yearmonth, axis=1)
                                    self.read_instance.data_in_memory_filtered[networkspeci][:,invalid_above] = np.NaN

                            # remove nans
                            invalid_nan = np.repeat(pd.isnull(self.read_instance.metadata_in_memory[networkspeci][meta_var][:, :]), 
                                                    self.read_instance.N_inds_per_yearmonth, axis=1)
                            self.read_instance.data_in_memory_filtered[networkspeci][:,invalid_nan] = np.NaN

    def validate_values(self):
        """ Validate that field inserted by user is float. """

        # iterate through metadata in memory
        for meta_var in self.read_instance.metadata_vars_to_read:

            metadata_type = self.read_instance.standard_metadata[meta_var]['metadata_type']
            metadata_data_type = self.read_instance.standard_metadata[meta_var]['data_type']

            if metadata_data_type != np.object:
                meta_var_index = self.read_instance.metadata_menu[metadata_type][
                    'rangeboxes']['labels'].index(meta_var)
                try:
                    np.float32(self.read_instance.metadata_menu[metadata_type][
                                   'rangeboxes']['current_lower'][meta_var_index])
                    np.float32(self.read_instance.metadata_menu[metadata_type][
                                   'rangeboxes']['current_upper'][meta_var_index])
                    return True
                except ValueError as e:
                    if self.read_instance.offline:
                        print("Warning: Error in metadata fields. The field of '{}' "
                              "should be numeric, \n{}".format(meta_var, str(e)))
                    else:
                        QtWidgets.QMessageBox.critical(self.read_instance, "Error in metadata fields",
                                                       "The field of '{}' should be numeric, \n{}"
                                                       .format(meta_var, str(e)),
                                                       QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.NoButton)
                    return False

    def temporally_colocate_data(self):
        """ Define function which temporally colocates observational and experiment data.
            This is done across all networks / species if spatial colocation is active,
            otherwise it is done independently per network / species
            This in reality means storing the indices for the temporal colocation.
        """

        # if do not have any experiment data loaded, no colocation is possible.
        if len(self.read_instance.data_labels) == 1:
            return
        else:
            # colocate observational data array to every different experiment array in memory, 
            # and all experiment arrays to observations.
            # wherever there is a NaN at one time in one of the observations/experiment arrays 
            # the other array value is also made NaN.
            # this is done across all networks / species if spatial colocation is active,
            # otherwise it is done inderpendently per network / species
            
            # iterate through network / species  
            for ii, networkspeci in enumerate(self.read_instance.networkspecies):

                # initialise arrays to determine where have NaNs
                if (ii == 0) or (not self.read_instance.spatial_colocation):
                    # create array for finding instances where have 0 valid values across all observations
                    # initialise as being all False (i.e. non-NaN), set as True on the occasion there is a NaN in the observations
                    obs_all_nan = np.full(self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,:].shape, False)

                    # create array for finding instances where have 0 valid values across all experiments
                    # initialise as being all False (i.e. non-NaN), set as True on the occasion there is a NaN in an experiment
                    exps_all_nan = np.full(obs_all_nan.shape, False)

                # get all instances observations is NaN
                nan_obs = np.isnan(self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,:])             

                # update obs_all_nan array, making True all instances where have NaNs
                # if all observations are nan then do not update
                if not np.all(nan_obs):
                    obs_all_nan = np.any([obs_all_nan, nan_obs], axis=0)

                # iterate through experiment data arrays in data in memory dictionary
                # save indices for colocation with observations
                for experiment in self.read_instance.experiments:
                    
                    #get expid data label index
                    exp_data_index = self.read_instance.data_labels.index(experiment)
                    
                    # get all instances experiment is NaN
                    nan_exp = np.isnan(self.read_instance.data_in_memory_filtered[networkspeci][exp_data_index,:,:])

                    # update exps_all_nan array, making True all instances where have NaNs
                    # if all experiment values are nan then do not update for that experiment
                    if not np.all(nan_exp):
                        exps_all_nan = np.any([exps_all_nan, nan_exp], axis=0)

                # if spatial colocation is not active,
                # get indices where one of observations and experiments per network /species is NaN
                if not self.read_instance.spatial_colocation:
                    self.read_instance.temporal_colocation_nans[networkspeci] = np.any([obs_all_nan, exps_all_nan], axis=0)

            # if spatial colocation is active, 
            # get indices where one of observations and experiments across networks / species is NaN
            if self.read_instance.spatial_colocation:
                for networkspeci in self.read_instance.networkspecies:
                    self.read_instance.temporal_colocation_nans[networkspeci] = np.any([obs_all_nan, exps_all_nan], axis=0)

    def get_valid_stations_after_filtering(self):
        """ Get valid station indices after all filtering has been performed.
            These are saved in a dictionary per network/species, per data label. 
            There is an mirror dictionary saved for the temporally colocated version of the data. 
        """

        # iterate through networkspecies  
        for networkspeci in self.read_instance.networkspecies:

            self.read_instance.valid_station_inds[networkspeci] = {}
            self.read_instance.valid_station_inds_temporal_colocation[networkspeci] = {}

            # get observational station indices with > 1 valid measurements 
            for data_label in self.read_instance.data_labels:

                # check if data array is observational data array
                if data_label == 'observations':

                    # get obs data array
                    obs_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[networkspeci][self.read_instance.data_labels.index(data_label),:,:])

                    # get absolute data availability number per station in observational data array
                    station_data_availability_number = Stats.calculate_data_avail_number(obs_data)

                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds[networkspeci][data_label] = \
                        np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]

                    if len(self.read_instance.data_labels) > 1:
                        # get colocated obs data array (if have experiments)
                        obs_data[self.read_instance.temporal_colocation_nans[networkspeci]] = np.NaN

                        # get absolute data availability number per station in observational data array
                        station_data_availability_number = Stats.calculate_data_avail_number(obs_data)

                        # get indices of stations with > 1 available measurements
                        self.read_instance.valid_station_inds_temporal_colocation[networkspeci][data_label] = \
                            np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]

            # get equivalent valid station indices for experimental arrays
            # subset observational valid station indices, with experimental array stations with > 1 valid measurements
            # therefore number of observational valid indices will always be >= experimental valid indices
            for data_label in self.read_instance.data_labels:

                # check if data array is not an observational data array
                if data_label != 'observations':

                    # get indices of valid observational data array stations
                    valid_station_inds = copy.deepcopy(self.read_instance.valid_station_inds[networkspeci]['observations'])

                    # get experimental data array (first subset by valid observational stations)
                    exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[networkspeci][self.read_instance.data_labels.index(data_label),valid_station_inds,:])

                    # get absolute data availability number per station in experiment data array
                    station_data_availability_number = Stats.calculate_data_avail_number(exp_data)

                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds[networkspeci][data_label] = \
                        valid_station_inds[np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]]

                    # get colocated experimental data array (first subset by valid observational stations)
                    exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[networkspeci][self.read_instance.data_labels.index(data_label),:,:])
                    exp_data[self.read_instance.temporal_colocation_nans[networkspeci]] = np.NaN
                    exp_data = exp_data[valid_station_inds,:]

                    # get absolute data availability number per station in experiment data array
                    station_data_availability_number = Stats.calculate_data_avail_number(exp_data)
                    
                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds_temporal_colocation[networkspeci][data_label] = \
                        valid_station_inds[np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]]