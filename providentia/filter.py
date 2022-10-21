from .calculate import Stats
<<<<<<< HEAD
from .config import split_options
=======
from .configuration import split_options
>>>>>>> master

import copy
import numpy as np
import pandas as pd

from PyQt5 import QtWidgets
<<<<<<< HEAD

=======
from providentia import aux
>>>>>>> master

class DataFilter:
    """
    "Class that filters observational/experiment data into memory as required."
    """

    def __init__(self, read_instance):
        self.read_instance = read_instance
<<<<<<< HEAD
=======

        # get indices of some data variables
        self.obs_index = self.read_instance.data_labels.index('observations')
        if self.read_instance.reading_ghost:
            if self.read_instance.resolution != 'daily' and self.read_instance.resolution != 'monthly':
                self.day_night_index = self.read_instance.ghost_data_vars_to_read.index('day_night_code')
            if self.read_instance.resolution != 'monthly':
                self.weekday_weekend_index = self.read_instance.ghost_data_vars_to_read.index('weekday_weekend_code')
            self.season_index = self.read_instance.ghost_data_vars_to_read.index('season_code')

        #apply filtering
>>>>>>> master
        self.filter_all()

    def filter_all(self):
        # call functions to start filtering
        self.reset_data_filter()
        self.filter_data_limits()
        self.filter_by_period()
        self.filter_by_data_availability()
        self.filter_by_metadata()
<<<<<<< HEAD
        self.colocate_data()
=======
        self.temporally_colocate_data()
        self.filter_by_species()
>>>>>>> master
        self.get_valid_stations_after_filtering()

    def reset_data_filter(self):
        """Resets data arrays to be un-filtered"""
<<<<<<< HEAD
        self.read_instance.data_in_memory_filtered = copy.deepcopy(self.read_instance.datareader.data_in_memory)
=======

        self.read_instance.data_in_memory_filtered = copy.deepcopy(self.read_instance.data_in_memory)
        self.read_instance.temporal_colocation_nans = {}
        self.read_instance.valid_station_inds = {}
        self.read_instance.valid_station_inds_temporal_colocation = {}
>>>>>>> master

    def filter_data_limits(self):
        """Filter out (set to NaN) data which exceed the lower/upper limits"""

<<<<<<< HEAD
        # get set lower/upper data bounds
        if self.read_instance.offline:
            species = self.read_instance.selected_species
            selected_lower_bound = self.read_instance.minimum_value
            selected_upper_bound = self.read_instance.maximum_value
        else:
            species = self.read_instance.active_species
            selected_lower_bound = self.read_instance.le_minimum_value.text()
            selected_upper_bound = self.read_instance.le_maximum_value.text()

        # check selected lower/upper bounds are numbers
        try:
            selected_lower_bound = np.float32(selected_lower_bound)
            selected_upper_bound = np.float32(selected_upper_bound)
        # if any of the fields are not numbers, return from function
        except ValueError:
            return

        # filter all observational data out of bounds of lower/upper limits
        inds_out_of_bounds = np.logical_or(self.read_instance.data_in_memory_filtered[
                                               'observations'][species] < selected_lower_bound,
                                           self.read_instance.data_in_memory_filtered[
                                               'observations'][species] > selected_upper_bound)
        self.read_instance.data_in_memory_filtered['observations'][species][inds_out_of_bounds] = np.NaN
=======
        # iterate through network / species  
        for network, speci in zip(self.read_instance.network, self.read_instance.species):

            # get network / species str
            networkspeci = '{}|{}'.format(network,speci)

            # get set lower/upper data bounds
            if self.read_instance.offline:
                lower_bound, upper_bound = aux.which_bounds(self.read_instance, 
                                                            self.read_instance.species[0])
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
>>>>>>> master

    def filter_by_period(self):
        """Filters data for selected periods (keeping or removing data, as defined)"""

        keeps, removes = [], []
        if self.read_instance.offline:
<<<<<<< HEAD
            species = self.read_instance.selected_species
            if hasattr(self.read_instance, 'period'):
                keeps, removes = split_options(self.read_instance.period)
                print(keeps, removes)
        else:
            species = self.read_instance.active_species
=======
            if hasattr(self.read_instance, 'period'):
                keeps, removes = split_options(self.read_instance.period)
        else:
>>>>>>> master
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
<<<<<<< HEAD
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['day_night_code'],
                            day_night_codes_to_keep, invert=True)] = np.NaN
=======
                if (self.read_instance.resolution != 'daily') & (self.read_instance.resolution != 'monthly'):
                    # iterate through network / species  
                    for network, speci in zip(self.read_instance.network, self.read_instance.species):
                        networkspeci = '{}|{}'.format(network,speci)
                        inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.day_night_index,:,:], day_night_codes_to_keep, invert=True)
                        self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN
>>>>>>> master

            weekday_weekend_codes_to_keep = []
            if 'Weekday' in keeps:
                weekday_weekend_codes_to_keep.append(0)
            if 'Weekend' in keeps:
                weekday_weekend_codes_to_keep.append(1)
            if len(weekday_weekend_codes_to_keep) == 1:
<<<<<<< HEAD
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['weekday_weekend_code'],
                            weekday_weekend_codes_to_keep, invert=True)] = np.NaN
=======
                if self.read_instance.resolution != 'monthly':
                    # iterate through network / species  
                    for network, speci in zip(self.read_instance.network, self.read_instance.species):
                        networkspeci = '{}|{}'.format(network,speci)
                        inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.weekday_weekend_index,:,:], weekday_weekend_codes_to_keep, invert=True)
                        self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN
>>>>>>> master

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
<<<<<<< HEAD
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['season_code'],
                            season_codes_to_keep, invert=True)] = np.NaN
=======
                # iterate through network / species  
                for network, speci in zip(self.read_instance.network, self.read_instance.species):
                    networkspeci = '{}|{}'.format(network,speci)
                    inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.season_index,:,:], season_codes_to_keep, invert=True)
                    self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN
>>>>>>> master

        if len(removes) > 0:
            day_night_codes_to_remove = []
            if 'Daytime' in removes:
                day_night_codes_to_remove.append(0)
            if 'Nighttime' in removes:
                day_night_codes_to_remove.append(1)
            if len(day_night_codes_to_remove) > 0:
<<<<<<< HEAD
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['day_night_code'],
                            day_night_codes_to_remove)] = np.NaN
=======
                if (self.read_instance.resolution != 'daily') & (self.read_instance.resolution != 'monthly'):
                    # iterate through network / species  
                    for network, speci in zip(self.read_instance.network, self.read_instance.species):
                        networkspeci = '{}|{}'.format(network,speci)
                        inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.day_night_index,:,:], day_night_codes_to_remove)
                        self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN
>>>>>>> master

            weekday_weekend_codes_to_remove = []
            if 'Weekday' in removes:
                weekday_weekend_codes_to_remove.append(0)
            if 'Weekend' in removes:
                weekday_weekend_codes_to_remove.append(1)
            if len(weekday_weekend_codes_to_remove) > 0:
<<<<<<< HEAD
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['weekday_weekend_code'],
                            weekday_weekend_codes_to_remove)] = np.NaN
=======
                if self.read_instance.resolution != 'monthly':
                    # iterate through network / species  
                    for network, speci in zip(self.read_instance.network, self.read_instance.species):
                        networkspeci = '{}|{}'.format(network,speci)
                        inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.weekday_weekend_index,:,:], weekday_weekend_codes_to_remove)
                        self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN
>>>>>>> master

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
<<<<<<< HEAD
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['season_code'],
                            season_codes_to_remove)] = np.NaN
=======
                # iterate through network / species  
                for network, speci in zip(self.read_instance.network, self.read_instance.species):
                    networkspeci = '{}|{}'.format(network,speci)
                    inds_to_screen = np.isin(self.read_instance.ghost_data_in_memory[networkspeci][self.season_index,:,:], season_codes_to_remove)
                    self.read_instance.data_in_memory_filtered[networkspeci][:, inds_to_screen] = np.NaN
>>>>>>> master

    def filter_by_data_availability(self):
        """Function which filters data by selected data availability variables"""

<<<<<<< HEAD
        if self.read_instance.offline:
            species = self.read_instance.selected_species
        else:
            species = self.read_instance.active_species

        # get set variables names representing percentage data availability (native and non-native)
        active_data_availablity_vars = self.read_instance.representativity_menu['rangeboxes']['labels']
=======
        # get set variables names representing percentage data availability (native and non-native)
        active_data_availablity_vars = self.read_instance.representativity_menu['rangeboxes']['map_vars']
>>>>>>> master

        try:
            data_availability_lower_bounds = []
            for var_ii, var in enumerate(active_data_availablity_vars):
                data_availability_lower_bounds.append(
                    np.float32(self.read_instance.representativity_menu['rangeboxes']['current_lower'][var_ii]))
        # if any of the fields are not numbers, return from function
        except ValueError:
<<<<<<< HEAD
            return

        if not self.read_instance.reading_nonghost:
            for var_ii, var in enumerate(active_data_availablity_vars):
                if 'native' in var:
                    # max gap variable?
                    if 'max_gap' in var:
                        # bound is < 100?:
                        if data_availability_lower_bounds[var_ii] < 100:
                            self.read_instance.data_in_memory_filtered['observations'][species][
                                self.read_instance.data_in_memory_filtered['observations'][var] >
                                data_availability_lower_bounds[var_ii]] = np.NaN
                    # data representativity variable?
                    else:
                        # bound is > 0?
                        if data_availability_lower_bounds[var_ii] > 0:
                            self.read_instance.data_in_memory_filtered['observations'][species][
                                self.read_instance.data_in_memory_filtered['observations'][var] <
                                data_availability_lower_bounds[var_ii]] = np.NaN

        # filter all observational data out of set bounds of non-native percentage data availability variables
=======
            print("Warning: Data availability fields must be numeric")
            return

        # filter observations by native percentage data availability variables (only GHOST data)
        if self.read_instance.reading_ghost:
            for var_ii, var in enumerate(active_data_availablity_vars):
                if 'native' in var:
                    var_index = self.read_instance.ghost_data_vars_to_read.index(var)
                    
                    # iterate through network / species  
                    for network, speci in zip(self.read_instance.network, self.read_instance.species):

                        # get network / species str
                        networkspeci = '{}|{}'.format(network,speci)

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
>>>>>>> master
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
<<<<<<< HEAD
                period_inds = np.arange(
                    self.read_instance.data_in_memory_filtered['observations'][species].shape[
                        1])
=======
                period_inds = np.arange(len(self.read_instance.time_array))
>>>>>>> master
                # daily variable?
                if period == 'daily':
                    period_inds_split = np.array_split(period_inds,
                                                       [24 * i for i in range(1, int(np.ceil(len(period_inds) / 24)))])
                # monthly variable?
                elif period == 'monthly':
<<<<<<< HEAD
                    period_inds_split = np.array_split(period_inds, np.cumsum(self.read_instance.datareader.N_inds_per_month))
=======
                    period_inds_split = np.array_split(period_inds, np.cumsum(self.read_instance.N_inds_per_yearmonth))
>>>>>>> master
                # whole record variable?
                else:
                    period_inds_split = [period_inds]

                # iterate through indices associated with periodic chunks for current period
                for period_inds in period_inds_split:
                    if len(period_inds) > 0:
<<<<<<< HEAD
                        # max gap variable?
                        if 'max_gap' in var:
                            max_gap_percent = Stats.max_repeated_nans_fraction(
                                self.read_instance.data_in_memory_filtered['observations'][species][:, period_inds])
                            self.read_instance.data_in_memory_filtered['observations'][species][
                                max_gap_percent > data_availability_lower_bounds[var_ii]] = np.NaN
                        # data representativity variable?
                        else:
                            data_availability_percent = Stats.calculate_data_avail_fraction(
                                self.read_instance.data_in_memory_filtered['observations'][species][:, period_inds])
                            self.read_instance.data_in_memory_filtered['observations'][species][
                                data_availability_percent < data_availability_lower_bounds[var_ii]] = np.NaN

    def filter_by_metadata(self):
        """Filters data by selected metadata"""

        if self.read_instance.offline:
            species = self.read_instance.selected_species
        else:
            species = self.read_instance.active_species

=======

                        # iterate through network / species  
                        for network, speci in zip(self.read_instance.network, self.read_instance.species):

                            # get network / species str
                            networkspeci = '{}|{}'.format(network,speci)

                            # max gap variable?
                            if 'max_gap' in var:
                                max_gap_percent = Stats.max_repeated_nans_fraction(
                                    self.read_instance.data_in_memory_filtered[networkspeci][:,:,period_inds])
                                inds_to_screen = max_gap_percent > data_availability_lower_bounds[var_ii]
                                self.read_instance.data_in_memory_filtered[networkspeci][inds_to_screen] = np.NaN
                            # data representativity variable?
                            else:
                                data_availability_percent = Stats.calculate_data_avail_fraction(
                                    self.read_instance.data_in_memory_filtered[networkspeci][:,:,period_inds])
                                inds_to_screen = data_availability_percent < data_availability_lower_bounds[var_ii]
                                self.read_instance.data_in_memory_filtered[networkspeci][inds_to_screen] = np.NaN
    
    def filter_by_metadata(self):
        """Filters data by selected metadata"""

>>>>>>> master
        # validate fields before filtering
        if not self.validate_values():
            return

<<<<<<< HEAD
        # iterate through all metadata
        for meta_var in self.read_instance.metadata_vars_to_read:

=======
        # iterate through metadata in memory
        for meta_var in self.read_instance.metadata_vars_to_read:
            
>>>>>>> master
            if meta_var == 'lat':
                meta_var = 'latitude'
            elif meta_var == 'lon':
                meta_var = 'longitude'

            metadata_type = self.read_instance.standard_metadata[meta_var]['metadata_type']
            metadata_data_type = self.read_instance.standard_metadata[meta_var]['data_type']

            # handle non-numeric metadata
            if metadata_data_type == np.object:
<<<<<<< HEAD
                # if any of the keep checkboxes are selected, filter out data by fields that have not been selected
                current_keep = self.read_instance.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected']
                if len(current_keep) > 0:
                    invalid_keep = np.repeat(
                        np.isin(self.read_instance.datareader.metadata_in_memory[meta_var][:, :],
                                current_keep, invert=True), self.read_instance.datareader.N_inds_per_month, axis=1)
                    self.read_instance.data_in_memory_filtered['observations'][species][
                        invalid_keep] = np.NaN
                # if any of the remove checkboxes have been selected, filter out data by these selected fields
                current_remove = self.read_instance.metadata_menu[metadata_type][meta_var]['checkboxes'][
                    'remove_selected']
                if len(current_remove) > 0:
                    invalid_remove = np.repeat(
                        np.isin(self.read_instance.datareader.metadata_in_memory[meta_var][:, :], current_remove),
                        self.read_instance.datareader.N_inds_per_month, axis=1)
                    self.read_instance.data_in_memory_filtered['observations'][species][
                        invalid_remove] = np.NaN
=======

                # iterate through network / species  
                for network, speci in zip(self.read_instance.network, self.read_instance.species):

                    # get network / species str
                    networkspeci = '{}|{}'.format(network,speci)

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
            
>>>>>>> master
            # handle numeric metadata
            else:
                meta_var_index = self.read_instance.metadata_menu[metadata_type]['rangeboxes']['labels'].index(meta_var)
                current_lower = np.float32(
                    self.read_instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index])
                current_upper = np.float32(
                    self.read_instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index])
<<<<<<< HEAD

                # if current lower value is non-NaN, then filter out data with metadata < current lower value
                if not pd.isnull(current_lower):
                    lower_default = np.float32(
                        self.read_instance.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index])
                    if current_lower > lower_default:
                        if not self.read_instance.reading_nonghost:
                            invalid_below = np.repeat(self.read_instance.datareader.metadata_in_memory[meta_var][:, :] <
                                                      current_lower, self.read_instance.datareader.N_inds_per_month, axis=1)
                        else:
                            invalid_below = np.repeat(self.read_instance.datareader.nonghost_metadata[meta_var][:, :] <
                                                      current_lower, self.read_instance.datareader.N_inds_per_month, axis=1)
                        self.read_instance.data_in_memory_filtered['observations'][species][
                            invalid_below] = np.NaN
                # if current upper < than the maximum extent, then filter out
                # data with metadata > current upper value (if this is numeric)
                # if current upper value is non-NaN, then filter out data with metadata > current upper value
                if not pd.isnull(current_upper):
                    upper_default = np.float32(
                        self.read_instance.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index])
                    if current_upper < upper_default:
                        if not self.read_instance.reading_nonghost:
                            invalid_above = np.repeat(self.read_instance.datareader.metadata_in_memory[meta_var][:, :] >
                                                      current_upper, self.read_instance.datareader.N_inds_per_month, axis=1)
                        else:
                            invalid_above = np.repeat(self.read_instance.datareader.nonghost_metadata[meta_var][:, :] >
                                                      current_upper, self.read_instance.datareader.N_inds_per_month, axis=1)
                        self.read_instance.data_in_memory_filtered['observations'][species][
                            invalid_above] = np.NaN
=======
                
                # get array with selected filters (those with the check Apply on)
                current_apply = self.read_instance.metadata_menu[metadata_type]['rangeboxes']['apply_selected']

                if len(current_apply) > 0:
                    # apply bounds and remove nans if variable has been selected
                    if meta_var in current_apply:

                        # iterate through network / species  
                        for network, speci in zip(self.read_instance.network, self.read_instance.species):

                            # get network / species str
                            networkspeci = '{}|{}'.format(network,speci)

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
>>>>>>> master

    def validate_values(self):
        """Validates that field inserted by user is float"""

<<<<<<< HEAD
        # iterate through all metadata
        for meta_var in self.read_instance.metadata_vars_to_read:
=======
        # iterate through metadata in memory
        for meta_var in self.read_instance.metadata_vars_to_read:

>>>>>>> master
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
<<<<<<< HEAD
                        print("Error in metadata fields. The field of '{}' "
=======
                        print("Warning: Error in metadata fields. The field of '{}' "
>>>>>>> master
                              "should be numeric, \n{}".format(meta_var, str(e)))
                    else:
                        QtWidgets.QMessageBox.critical(self.read_instance, "Error in metadata fields",
                                                       "The field of '{}' should be numeric, \n{}"
                                                       .format(meta_var, str(e)),
                                                       QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.NoButton)
                    return False

<<<<<<< HEAD
    def colocate_data(self):
        """Define function which colocates observational and experiment data"""

        if self.read_instance.offline:
            species = self.read_instance.selected_species
        else:
            species = self.read_instance.active_species

        # if do not have any experiment data loaded, no colocation is possible,
        # therefore return from function
        if len(list(self.read_instance.datareader.data_in_memory.keys())) == 1:
            return
        else:
            # otherwise, colocate observational and experiment data, creating new data arrays

            # colocate observational data array to every different experiment array in memory, and vice versa
            # wherever there is a NaN at one time in one of the observations/experiment arrays,
            # the other array value is also made NaN

            # get all instances observations are NaN
            nan_obs = np.isnan(
                self.read_instance.data_in_memory_filtered['observations'][species])

            # create array for finding instances where have 0 valid values across all experiments
            # initialise as being all True, set as False on the occasion there is a valid value in an experiment
            exps_all_nan = np.full(nan_obs.shape, True)

            # get name of all experiment labels in memory
            exp_labels = sorted(list(self.read_instance.datareader.data_in_memory.keys()))
            exp_labels.remove('observations')

            # iterate through experiment data arrays in data in memory dictionary
            exp_nan_dict = {}
            for exp_label in exp_labels:
                # get all instances experiment are NaN
                exp_nan_dict[exp_label] = np.isnan(self.read_instance.data_in_memory_filtered[exp_label][species])
                # get all instances where either the observational array or experiment array are NaN at a given time
                nan_instances = np.any([nan_obs, exp_nan_dict[exp_label]], axis=0)
                # create new observational array colocated to experiment
                obs_data = copy.deepcopy(self.read_instance.data_in_memory_filtered['observations'])
                obs_data[nan_instances] = np.NaN
                self.read_instance.data_in_memory_filtered['observations_colocatedto_{}'.format(exp_label)] = obs_data
                self.read_instance.datareader.plotting_params['observations_colocatedto_{}'.format(exp_label)] = {
                    'colour': self.read_instance.datareader.plotting_params['observations']['colour'],
                    'zorder': self.read_instance.datareader.plotting_params['observations']['zorder']}
                # create new experiment array colocated to observations
                exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[exp_label])
                exp_data[nan_instances] = np.NaN
                self.read_instance.data_in_memory_filtered['{}_colocatedto_observations'.format(exp_label)] = exp_data
                self.read_instance.datareader.plotting_params['{}_colocatedto_observations'.format(exp_label)] = {
                    'colour': self.read_instance.datareader.plotting_params[exp_label]['colour'],
                    'zorder': self.read_instance.datareader.plotting_params[exp_label]['zorder']}
                # update exps_all_nan array, making False all instances where have valid experiment data
                exps_all_nan = np.all([exps_all_nan, exp_nan_dict[exp_label]], axis=0)

            # colocate experiments with all other experiments
            for exp_label_ii, exp_label in enumerate(exp_labels):
                for exp_label_2 in exp_labels[exp_label_ii + 1:]:
                    # get all instances where either of the experiment arrays are NaN at a given time
                    nan_instances = np.any([exp_nan_dict[exp_label], exp_nan_dict[exp_label_2]], axis=0)
                    # create new experiment array for experiment1 colocated to experiment2
                    exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[exp_label])
                    exp_data[nan_instances] = np.NaN
                    self.read_instance.data_in_memory_filtered[
                        '{}_colocatedto_{}'.format(exp_label, exp_label_2)] = exp_data
                    self.read_instance.datareader.plotting_params['{}_colocatedto_{}'.format(exp_label, exp_label_2)] = {
                        'colour': self.read_instance.datareader.plotting_params[exp_label]['colour'],
                        'zorder': self.read_instance.datareader.plotting_params[exp_label]['zorder']}
                    # create new experiment array for experiment2 colocated to experiment1
                    exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[exp_label_2])
                    exp_data[nan_instances] = np.NaN
                    self.read_instance.data_in_memory_filtered[
                        '{}_colocatedto_{}'.format(exp_label_2, exp_label)] = exp_data
                    self.read_instance.datareader.plotting_params['{}_colocatedto_{}'.format(exp_label_2, exp_label)] = {
                        'colour': self.read_instance.datareader.plotting_params[exp_label_2]['colour'],
                        'zorder': self.read_instance.datareader.plotting_params[exp_label_2]['zorder']}

            # create observational data array colocated to be non-NaN whenever
            # there is a valid data in at least 1 experiment
            exps_all_nan = np.any([nan_obs, exps_all_nan], axis=0)
            obs_data = copy.deepcopy(self.read_instance.data_in_memory_filtered['observations'])
            obs_data[exps_all_nan] = np.NaN
            self.read_instance.data_in_memory_filtered['observations_colocatedto_experiments'] = obs_data
            self.read_instance.datareader.plotting_params['observations_colocatedto_experiments'] = {
                'colour': self.read_instance.datareader.plotting_params['observations']['colour'],
                'zorder': self.read_instance.datareader.plotting_params['observations']['zorder']}

    def get_valid_stations_after_filtering(self):

        if self.read_instance.offline:
            species = self.read_instance.selected_species
        else:
            species = self.read_instance.active_species

        # get intersect of indices of stations with >= % minimum data availability percent,
        # and with > 1 valid measurements ,in all observational data arrays (colocated and non-colocated)
        # then subset these indices with standard methods == selected methods,
        # iterate through all data arrays
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()):

            # check if data array is an observational data array
            if data_label.split('_')[0] == 'observations':
                # calculate data availability fraction per station in observational data array
                # get absolute data availability number per station in observational data array
                station_data_availability_number = Stats.calculate_data_avail_number(
                    self.read_instance.data_in_memory_filtered[data_label][species])

                # get indices of stations with > 1 available measurements
                # save valid station indices with data array
                self.read_instance.datareader.plotting_params[data_label]['valid_station_inds'] = \
                    np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]

        # write valid station indices calculated for observations across to associated experimental data arrays
        # iterate through all data arrays
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()):

            # check if data array is not an observational data array
            if data_label.split('_')[0] != 'observations':

                # handle colocated experimental arrays
                if '_colocatedto_' in data_label:
                    exp_name = data_label.split('_colocatedto_')[0]
                    self.read_instance.datareader.plotting_params[data_label]['valid_station_inds'] = \
                        copy.deepcopy(
                            self.read_instance.datareader.plotting_params['observations_colocatedto_{}'.format(exp_name)][
                                'valid_station_inds'])
                # handle non-colocated experimental arrays
                else:
                    self.read_instance.datareader.plotting_params[data_label]['valid_station_inds'] = \
                        copy.deepcopy(self.read_instance.datareader.plotting_params['observations']['valid_station_inds'])

        # after subsetting by pre-written associated observational valid stations, get
        # indices of stations with > 1 valid measurements in all experiment data arrays (colocated and non-colocated)
        # iterate through all data arrays
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()):

            # check if data array is not an observational data array
            if data_label.split('_')[0] != 'observations':
                # get indices of associated observational data array valid stations
                # (pre-written to experiment data arrays)
                valid_station_inds = self.read_instance.datareader.plotting_params[data_label]['valid_station_inds']
                # get absolute data availability number per station in experiment data array
                # after subsetting valid observational stations (i.e. number of non-NaN measurements)
                # update stats object data and call data availability function
                station_data_availability_number = \
                    Stats.calculate_data_avail_number(
                        self.read_instance.data_in_memory_filtered[data_label][species][valid_station_inds, :])
                # get indices of stations with > 1 available measurements
                valid_station_inds = valid_station_inds[np.arange(len(station_data_availability_number), dtype=np.int)[
                    station_data_availability_number > 1]]
                # overwrite previous written valid station indices (now at best a subset of those indices)
                self.read_instance.datareader.plotting_params[data_label]['valid_station_inds'] = valid_station_inds
=======
    def temporally_colocate_data(self):
        """Define function which temporally colocates observational and experiment data.
           This is done across all networks / species if spatial colocation is active,
           otherwise it is done inderpendently per network / species
           This in reality means storing the indices for the temporal colocation.
        """

        # if do not have any experiment data loaded, no colocation is possible.
        if len(self.read_instance.data_labels) == 1:
            return
        else:
            # colocate observational data array to every different experiment array in memory, 
            # and all experiment arrays to observations.
            # wherever there is a NaN at one time in one of the observations/experiment arrays across all networks / species
            # the other array value is also made NaN

            # iterate through network / species  
            for ii, (network, speci) in enumerate(zip(self.read_instance.network, self.read_instance.species)):

                # get network / species str
                networkspeci = '{}|{}'.format(network,speci)

                if (self.read_instance.spatial_colocation) or (ii == 0):
                    # create array for finding instances where have 0 valid values across all observations in all networks / species
                    # initialise as being all False (i.e. non-NaN), set as True on the occasion there is a NaN in the observations
                    obs_all_nan = np.full(self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,:].shape, False)

                    # create array for finding instances where have 0 valid values across all experiments, in all networks / species
                    # initialise as being all False (i.e. non-NaN), set as True on the occasion there is a NaN in an experiment
                    exps_all_nan = np.full(obs_all_nan.shape, False)

                # get all instances observations is NaN
                nan_obs = np.isnan(self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,:])
                
                # update obs_all_nan array, making True all instances where have NaNs
                obs_all_nan = np.any([obs_all_nan, nan_obs], axis=0)

                # iterate through experiment data arrays in data in memory dictionary
                # save indices for colocation with observations
                for experiment in self.read_instance.experiments:
                    
                    #get expid data label index
                    exp_data_index = self.read_instance.data_labels.index(experiment)
                    
                    # get all instances experiment is NaN
                    nan_exp = np.isnan(self.read_instance.data_in_memory_filtered[networkspeci][exp_data_index,:,:])
                
                    # update exps_all_nan array, making True all instances where have NaNs
                    exps_all_nan = np.any([exps_all_nan, nan_exp], axis=0)

                # if spatial colocation is not active,
                # get indices where one of observations and experiments per network /species is NaN
                if not self.read_instance.spatial_colocation:
                    self.read_instance.temporal_colocation_nans[networkspeci] = np.any([obs_all_nan, exps_all_nan], axis=0)

            # get indices where one of observations and experiments across networks / species is NaN
            if self.read_instance.spatial_colocation:
                for network, speci in zip(self.read_instance.network, self.read_instance.species):
                    self.read_instance.temporal_colocation_nans[networkspeci] = np.any([obs_all_nan, exps_all_nan], axis=0)

    def filter_by_species(self):
        """Define function which filters read species by other species."""

        pass

    def get_valid_stations_after_filtering(self):
        """Get valid station indices after all filtering has been performed.
           These are saved in a dictionary per network/species, per data label. 
           There is an mirror dictionary saved for the temporally colocated version of the data. 
        """

        # iterate through networks / species  
        for network, speci in zip(self.read_instance.network, self.read_instance.species):

            # get network / species str
            networkspeci = '{}|{}'.format(network,speci)

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

                    # get colocated experimental data array (if have experiments, first subset by valid observational stations)
                    exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[networkspeci][self.read_instance.data_labels.index(data_label),:,:])
                    exp_data[self.read_instance.temporal_colocation_nans[networkspeci]] = np.NaN
                    exp_data = exp_data[valid_station_inds,:]

                    # get absolute data availability number per station in experiment data array
                    station_data_availability_number = Stats.calculate_data_avail_number(exp_data)
                    
                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds_temporal_colocation[networkspeci][data_label] = \
                        valid_station_inds[np.arange(len(station_data_availability_number), dtype=np.int)[station_data_availability_number > 1]]
>>>>>>> master
