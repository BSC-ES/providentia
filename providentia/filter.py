""" Class that filters observational/experiment data into memory as required """

import copy
import json
import os
import yaml

import numpy as np
import pandas as pd
import ast

from .calculate import Stats, ExpBias
from .configuration import split_options
from .statistics import get_z_statistic_info, exceedance_lim
from .warnings_prv import show_message

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])

class DataFilter:

    def __init__(self, read_instance):
        self.read_instance = read_instance

        # get indices of some data variables
        self.obs_index = self.read_instance.data_labels.index(self.read_instance.observations_data_label)
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
        self.filter_extreme_stations()
        self.temporally_colocate_data()
        self.apply_calibration_factor()
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
            Where values for each species are outside of these ranges, or NaN,
            then impose NaNs upon all read species in memory.
            Only filter if spatial colocation is True.
        """

        # filter all read species by set species ranges
        if (self.read_instance.filter_species) and (self.read_instance.spatial_colocation):

            # iterate through all species to filter by
            for filter_networkspeci, speci_all_limits in self.read_instance.filter_species.items():
                
                # get where data is inside bounds or NaN
                for speci_limit in speci_all_limits:

                    # initialise array to set where temporally to filter species
                    # initialse being all False, set as True where data is inside given bounds for species
                    inds_to_filter = np.full(self.read_instance.data_in_memory_filtered[self.read_instance.networkspecies[0]][self.obs_index,:,:].shape, False)    

                    # get lower and upper limits for species
                    lower_limit = speci_limit[0]
                    upper_limit = speci_limit[1]
                    filter_species_fill_value = speci_limit[2]

                    # remove symbols from limits and transform into float
                    if lower_limit != ':':
                        lower_limit_val = float(lower_limit.replace('>', '').replace('=', ''))
                    if upper_limit != ':':
                        upper_limit_val = float(upper_limit.replace('<', '').replace('=', ''))

                    # get filter conditions
                    if ':' in upper_limit and ':' in lower_limit:
                        print('Upper and lower bounds are :, no data filter will be applied for {0}'.format(filter_networkspeci))
                        return
                    if ':' in upper_limit:
                        if '=' in lower_limit:
                            valid_inds_per_species = (self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] >= lower_limit_val)
                        else:
                            valid_inds_per_species = (self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] > lower_limit_val)
                    elif ':' in lower_limit:
                        if '=' in upper_limit:
                            valid_inds_per_species = (self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] <= upper_limit_val)
                        else:
                            valid_inds_per_species = (self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] < upper_limit_val)
                    else:
                        if '=' in upper_limit and '=' in lower_limit:
                            valid_inds_per_species = np.logical_and.reduce((self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] >= lower_limit_val,
                                                                            self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] <= upper_limit_val))
                        elif '=' in upper_limit and '=' not in lower_limit:
                            valid_inds_per_species = np.logical_and.reduce((self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] > lower_limit_val,
                                                                            self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] <= upper_limit_val))
                        elif '=' not in upper_limit and '=' in lower_limit:
                            valid_inds_per_species = np.logical_and.reduce((self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] >= lower_limit_val,
                                                                            self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] < upper_limit_val))
                        else:
                            valid_inds_per_species = np.logical_and.reduce((self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] > lower_limit_val,
                                                                            self.read_instance.filter_data_in_memory[filter_networkspeci][:,:] < upper_limit_val))

                    valid_inds_per_species = np.logical_or.reduce((valid_inds_per_species, 
                                                                   np.isnan(self.read_instance.filter_data_in_memory[filter_networkspeci][:,:])))
                    
                    # update inds_to_filter array, making True all instances where we have data inside of bounds
                    inds_to_filter = np.any([inds_to_filter, valid_inds_per_species], axis=0)
                    
                    # set all inds to filter as fill value for all networkspecies in memory
                    for networkspeci in self.read_instance.networkspecies:
                        self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index, inds_to_filter] = filter_species_fill_value      

    def filter_data_limits(self):
        """ Filter out (set to NaN) data which exceed the lower/upper limits. """

        # iterate through networkspecies  
        for networkspeci in self.read_instance.networkspecies:

            # get speci str
            speci = networkspeci.split('|')[1]

            # get lower/upper data bounds
            if (self.read_instance.offline) or (self.read_instance.interactive):
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
                msg = 'Data limit fields must be numeric.'
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                return

            # filter all observational/experiment data out of bounds of lower/upper limits
            inds_out_of_bounds = np.logical_or(self.read_instance.data_in_memory_filtered[networkspeci][:,:,:] < lower_bound,
                                               self.read_instance.data_in_memory_filtered[networkspeci][:,:,:] > upper_bound)
            self.read_instance.data_in_memory_filtered[networkspeci][inds_out_of_bounds] = np.NaN

    def filter_by_period(self):
        """ Filter data for selected periods (keeping or removing data, as defined). """

        keeps, removes = [], []
        if (self.read_instance.offline) or (self.read_instance.interactive):
            if hasattr(self.read_instance, 'period'):
                keeps, removes = split_options(self.read_instance, self.read_instance.period)
            else:
                # Get period if apply_filter is used in interactive mode
                if self.read_instance.period_menu['checkboxes']['keep_selected']:
                    keeps = self.read_instance.period_menu['checkboxes']['keep_selected']
                if self.read_instance.period_menu['checkboxes']['remove_selected']:
                    removes = self.read_instance.period_menu['checkboxes']['remove_selected']         
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
            msg = 'Data availability fields must be numeric.'
            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
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

                            # data representativity variable?
                            data_availability_percent = Stats.calculate_data_avail_fraction(
                                self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,period_inds[0]:period_inds[-1]+1])
                            inds_to_screen = np.where(data_availability_percent < data_availability_lower_bounds[var_ii])[0]
                            self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,inds_to_screen[:,np.newaxis],period_inds[np.newaxis,:]] = np.NaN

    def filter_by_metadata(self):
        """ Filter data by selected metadata. """

        # iterate through metadata in memory
        for meta_var in self.read_instance.metadata_vars_to_read:
            
            # validate field before filtering
            if not self.validate_values(meta_var):
                # go to next variable if filter cannot be applied
                continue

            if meta_var == 'lat':
                meta_var = 'latitude'
            elif meta_var == 'lon':
                meta_var = 'longitude'

            metadata_type = self.read_instance.standard_metadata[meta_var]['metadata_type']
            metadata_data_type = self.read_instance.standard_metadata[meta_var]['data_type']

            # handle non-numeric metadata
            if metadata_data_type == object:

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

    def validate_values(self, meta_var):
        """ Validate that field inserted by user is float. """
        
        metadata_type = self.read_instance.standard_metadata[meta_var]['metadata_type']
        metadata_data_type = self.read_instance.standard_metadata[meta_var]['data_type']
        
        if metadata_data_type != object:
            meta_var_index = self.read_instance.metadata_menu[metadata_type][
                'rangeboxes']['labels'].index(meta_var)
            try:
                np.float32(self.read_instance.metadata_menu[metadata_type][
                                'rangeboxes']['current_lower'][meta_var_index])
                np.float32(self.read_instance.metadata_menu[metadata_type][
                                'rangeboxes']['current_upper'][meta_var_index])
                return True
            except ValueError as e:
                msg = "Error in metadata fields. The field of '{}' should be numeric.".format(meta_var)
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                if (not self.read_instance.offline) and (not self.read_instance.interactive):
                    self.read_instance.metadata_menu[metadata_type]['rangeboxes']['apply_selected'].remove(meta_var)
                return False
        else:
            return True

    def filter_extreme_stations(self):
        """ Define function which filters out extreme stations based on set statistical limits.
            There can be multiple limit arguments for a statistic e.g. 'MB': ['<10','>20']
            There can also be limits per species e.g. 'RMSE': {'sconco3': ['<50.0', '>70.0'], 'sconco':[>100.0]}
            An absolute statistic can be set to be a bias statistic by adding '_bias' e.g. 
            'p95_bias': ['<10','>20']

            If statistic is an absolute statistic, then only remove stations based on observations.
            If statistic is a bias statistic, then remove collection of all stations outside limits across all obs-exp
            comparsions.
        """

        # option to remove extreme stations set?
        if self.read_instance.remove_extreme_stations:

            # load yaml of defined stattistics limits
            remove_extreme_stations_fname = os.path.join(PROVIDENTIA_ROOT, 'settings/remove_extreme_stations.yaml')
            stat_defs = yaml.safe_load(open(remove_extreme_stations_fname))

            # get specific set of limits (if available)
            # throw wraning if not
            if self.read_instance.remove_extreme_stations not in stat_defs:
                msg = "'{}' not defined in '{}'. Not removing extreme stations.".format(self.read_instance.remove_extreme_stations, remove_extreme_stations_fname)
                show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                return
            else:
                stat_limits = stat_defs[self.read_instance.remove_extreme_stations]

            # loop through and calculate each statistic per station and remove stations outside statistical limits
            for zstat in stat_limits:
                # get list of statistical limits for specific stat 
                stat_arguments = stat_limits[zstat]
                
                # if have a dict, the limits are specific per species, so limit for species
                speci_specific_limits = False
                if type(stat_arguments) == dict:
                    speci_specific_limits = True                        

                # determine if station is absolute or bias statistic
                zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat) 

                # handle some possible errors

                # if have only observations data then cannot calculate bias statistic, so continue to next stat
                if z_statistic_sign == 'bias':
                    if len(self.read_instance.data_labels) == 1:
                        msg = "Cannot remove extreme stations via calculation of '{}' as no experiment data has been read.".format(zstat)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                        continue
                        
                # if temporal_colocation is not active then cannot calculate ExpBias statistic, so continue to next stat 
                if z_statistic_type == 'expbias':
                    if not self.read_instance.temporal_colocation:
                        msg = "Cannot remove extreme stations via calculation of '{}' as 'temporal_colocation' is not active.".format(zstat)
                        show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                        continue

                # get dictionary containing necessary information for calculation of selected statistic
                if z_statistic_type == 'basic':
                    stats_dict = self.read_instance.basic_stats[base_zstat]
                else:
                    stats_dict = self.read_instance.expbias_stats[base_zstat]

                # load default selected z statistic arguments for passing to statistical function
                function_arguments = stats_dict['arguments']

                # iterate through network / species  
                for ii, networkspeci in enumerate(self.read_instance.networkspecies):

                    # if stat is exceedances then add threshold value (if available)  
                    if base_zstat == 'Exceedances':
                        function_arguments['threshold'] = exceedance_lim(networkspeci)

                    # get list of statistic limits specific for speci (if wanted)
                    if speci_specific_limits:
                        # if have not defined limits for speci, then throw warning and continue to next speci
                        if networkspeci not in stat_arguments:
                            msg = "No statistical limits defined for '{}' in '{}' section of '{}'. Not removing extreme stations for '{}'.".format(
                                  networkspeci, self.read_instance.remove_extreme_stations, remove_extreme_stations_fname, speci)
                            show_message(self.read_instance, msg, from_conf=self.read_instance.from_conf)
                            continue
                        else:
                            specific_stat_arguments = stat_arguments[networkspeci]
                    else:
                        specific_stat_arguments = copy.deepcopy(stat_arguments)

                    # calculate statistic per station and then compare statistic against limits
                    # set stations exceeding limits to NaN
                    # calculate basic stats
                    if (z_statistic_type == 'basic') and (z_statistic_sign != 'bias'):            
                        data_array_a = self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,:]
                        calc_stat = np.array(getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments))
                        non_finite_stat = ~np.isfinite(calc_stat)
                        self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,non_finite_stat,:] = np.NaN
                        for specific_stat_argument in specific_stat_arguments:
                            invalid_stations = ast.literal_eval('calc_stat{}'.format(specific_stat_argument))
                            self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,invalid_stations,:] = np.NaN
                    # calculate basic bias stats and expbias stats
                    else:
                        for data_label in self.read_instance.data_labels:
                            if data_label != self.read_instance.observations_data_label:
                                #get expid data label index
                                exp_data_index = self.read_instance.data_labels.index(data_label)
                                data_array_a = self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,:]
                                data_array_b = self.read_instance.data_in_memory_filtered[networkspeci][exp_data_index,:,:]
                                # calculate basic bias stats
                                if (z_statistic_type == 'basic') and (z_statistic_sign == 'bias'): 
                                    statistic_a = np.array(getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments))
                                    statistic_b = np.array(getattr(Stats, stats_dict['function'])(data_array_b, **function_arguments))
                                    calc_stat = statistic_b - statistic_a
                                # calculate expbias stats
                                elif z_statistic_type == 'expbias':
                                    calc_stat = np.array(getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':data_array_a,'exp':data_array_b}}))
                                non_finite_stat = ~np.isfinite(calc_stat)
                                self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,non_finite_stat,:] = np.NaN
                                for specific_stat_argument in specific_stat_arguments:
                                    invalid_stations = ast.literal_eval('calc_stat{}'.format(specific_stat_argument))
                                    self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,invalid_stations,:] = np.NaN
                                    

    def temporally_colocate_data(self):
        """ Define function which temporally colocates observational and experiment data.
            If spatial colocation is active, then data is also temporally colocated across all network / species,
            otherwise it is done independently per network / species.
            This in reality means storing the indices for the temporal colocation.
        """

        # if do not have any experiment data loaded, no colocation is possible.
        #if not self.read_instance.temporal_colocation:
        #    # iterate through network / species setting all data as being available
        #    for ii, networkspeci in enumerate(self.read_instance.networkspecies):
        #        self.read_instance.temporal_colocation_nans[networkspeci] = np.full(self.read_instance.data_in_memory_filtered[networkspeci][self.obs_index,:,:].shape, False)
        #    return
        #else:
            # colocate observational data array to every different experiment array in memory, 
            # and all experiment arrays to observations.
            # wherever there is a NaN at one time in one of the observations/experiment arrays 
            # the other array value is also made NaN.
            # this is done across all networks / species if spatial colocation is active,
            # otherwise it is done independently per network / species
            
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
            for experiment in self.read_instance.experiments.values():
                
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
                if data_label == self.read_instance.observations_data_label:

                    # get obs data array
                    obs_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[networkspeci][self.read_instance.data_labels.index(data_label),:,:])

                    # get absolute data availability number per station in observational data array
                    if obs_data.size == 0:
                        station_data_availability_number = np.array([])
                    else:
                        station_data_availability_number = Stats.calculate_data_avail_number(obs_data)
                        
                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds[networkspeci][data_label] = \
                        np.arange(len(station_data_availability_number), dtype=np.int64)[station_data_availability_number > 1]

                    # get colocated obs data array if have temporal colocation is active
                    obs_data[self.read_instance.temporal_colocation_nans[networkspeci]] = np.NaN

                    # get absolute data availability number per station in observational data array
                    if obs_data.size == 0:
                        station_data_availability_number = np.array([])
                    else:
                        station_data_availability_number = Stats.calculate_data_avail_number(obs_data)
                    
                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds_temporal_colocation[networkspeci][data_label] = \
                        np.arange(len(station_data_availability_number), dtype=np.int64)[station_data_availability_number > 1]

            # get equivalent valid station indices for experimental arrays
            # subset observational valid station indices, with experimental array stations with > 1 valid measurements
            # therefore number of observational valid indices will always be > experimental valid indices
            for data_label in self.read_instance.data_labels:

                # check if data array is not an observational data array
                if data_label != self.read_instance.observations_data_label:

                    # get indices of valid observational data array stations
                    valid_station_inds = copy.deepcopy(self.read_instance.valid_station_inds[networkspeci][self.read_instance.observations_data_label])

                    # get experimental data array (first subset by valid observational stations)
                    exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[networkspeci][self.read_instance.data_labels.index(data_label),valid_station_inds,:])

                    # get absolute data availability number per station in experiment data array
                    if exp_data.size == 0:
                        station_data_availability_number = np.array([])
                    else:
                        station_data_availability_number = Stats.calculate_data_avail_number(exp_data)
                    
                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds[networkspeci][data_label] = \
                        valid_station_inds[np.arange(len(station_data_availability_number), dtype=np.int64)[station_data_availability_number > 1]]
                    
                    # get colocated experimental data array (first subset by valid observational stations)
                    exp_data = copy.deepcopy(self.read_instance.data_in_memory_filtered[networkspeci][self.read_instance.data_labels.index(data_label),:,:])
                    exp_data[self.read_instance.temporal_colocation_nans[networkspeci]] = np.NaN
                    exp_data = exp_data[valid_station_inds,:]

                    # get absolute data availability number per station in experiment data array
                    if exp_data.size == 0:
                        station_data_availability_number = np.array([])
                    else:
                        station_data_availability_number = Stats.calculate_data_avail_number(exp_data)
                    
                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds_temporal_colocation[networkspeci][data_label] = \
                        valid_station_inds[np.arange(len(station_data_availability_number), dtype=np.int64)[station_data_availability_number > 1]]

    def apply_calibration_factor(self):
        """ Apply calibration factor to add or subtract a number to the experiments, 
            multiply or divide the experiment data by a certain value.
        """

        if self.read_instance.calibration_factor:

            # iterate through networkspecies  
            for networkspeci_ii, networkspeci in enumerate(self.read_instance.networkspecies):      
                
                # remove observations from data labels
                relevant_data_labels = copy.deepcopy(self.read_instance.data_labels)
                relevant_data_labels.remove(self.read_instance.observations_data_label)

                # get calibration factor per experiment
                for data_label_ii, data_label in enumerate(relevant_data_labels):

                    # get calibration factor per experiment
                    if isinstance(self.read_instance.calibration_factor, dict):
                        exp_label = list(self.read_instance.experiments.keys())[
                            list(self.read_instance.experiments.values()).index(data_label)]
                        calibration_factor = self.read_instance.calibration_factor[exp_label]
                    else:
                        calibration_factor = self.read_instance.calibration_factor

                    # get calibration factor per networkspeci
                    if (len(self.read_instance.networkspecies) > 1) and (',' in calibration_factor):
                        calibration_factor = calibration_factor.split(',')[networkspeci_ii]
                    
                    msg = 'Applying calibration factor: '
                    msg += '{0} in {1} to {2}'.format(calibration_factor, data_label, networkspeci)
                    print(msg)
                    
                    # apply calibration factor
                    if '*' in calibration_factor:
                        self.read_instance.data_in_memory_filtered[networkspeci][data_label_ii+1,:,:] *= \
                            float(calibration_factor.replace('*', ''))
                    elif '/' in calibration_factor:
                        self.read_instance.data_in_memory_filtered[networkspeci][data_label_ii+1,:,:] /= \
                            float(calibration_factor.replace('/', ''))
                    elif '-' in calibration_factor:
                        self.read_instance.data_in_memory_filtered[networkspeci][data_label_ii+1,:,:] -= \
                            float(calibration_factor.replace('-', ''))
                    else:
                        self.read_instance.data_in_memory_filtered[networkspeci][data_label_ii+1,:,:] += \
                            float(calibration_factor)
