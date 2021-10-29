from .calculate import Stats
from .config import split_options

import copy
import numpy as np
import pandas as pd

from PyQt5 import QtWidgets


class DataFilter:
    """
    "Class that filters observational/experiment data into memory as required."
    """

    def __init__(self, read_instance):
        self.read_instance = read_instance
        self.filter_all()

    def filter_all(self):
        # call functions to start filtering
        self.reset_data_filter()
        self.filter_data_limits()
        self.filter_by_period()
        self.filter_by_data_availability()
        self.filter_by_metadata()
        self.colocate_data()
        self.get_valid_stations_after_filtering()

    def reset_data_filter(self):
        """Resets data arrays to be un-filtered"""
        self.read_instance.data_in_memory_filtered = copy.deepcopy(self.read_instance.datareader.data_in_memory)

    def filter_data_limits(self):
        """Filter out (set to NaN) data which exceed the lower/upper limits"""

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

    def filter_by_period(self):
        """Filters data for selected periods (keeping or removing data, as defined)"""

        keeps, removes = [], []
        if self.read_instance.offline:
            species = self.read_instance.selected_species
            if hasattr(self.read_instance, 'period'):
                keeps, removes = split_options(self.read_instance.period)
                print(keeps, removes)
        else:
            species = self.read_instance.active_species
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
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['day_night_code'],
                            day_night_codes_to_keep, invert=True)] = np.NaN

            weekday_weekend_codes_to_keep = []
            if 'Weekday' in keeps:
                weekday_weekend_codes_to_keep.append(0)
            if 'Weekend' in keeps:
                weekday_weekend_codes_to_keep.append(1)
            if len(weekday_weekend_codes_to_keep) == 1:
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['weekday_weekend_code'],
                            weekday_weekend_codes_to_keep, invert=True)] = np.NaN

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
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['season_code'],
                            season_codes_to_keep, invert=True)] = np.NaN

        if len(removes) > 0:
            day_night_codes_to_remove = []
            if 'Daytime' in removes:
                day_night_codes_to_remove.append(0)
            if 'Nighttime' in removes:
                day_night_codes_to_remove.append(1)
            if len(day_night_codes_to_remove) > 0:
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['day_night_code'],
                            day_night_codes_to_remove)] = np.NaN

            weekday_weekend_codes_to_remove = []
            if 'Weekday' in removes:
                weekday_weekend_codes_to_remove.append(0)
            if 'Weekend' in removes:
                weekday_weekend_codes_to_remove.append(1)
            if len(weekday_weekend_codes_to_remove) > 0:
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['weekday_weekend_code'],
                            weekday_weekend_codes_to_remove)] = np.NaN

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
                self.read_instance.data_in_memory_filtered['observations'][species][
                    np.isin(self.read_instance.data_in_memory_filtered['observations']['season_code'],
                            season_codes_to_remove)] = np.NaN

    def filter_by_data_availability(self):
        """Function which filters data by selected data availability variables"""

        if self.read_instance.offline:
            species = self.read_instance.selected_species
        else:
            species = self.read_instance.active_species

        # get set variables names representing percentage data availability (native and non-native)
        active_data_availablity_vars = self.read_instance.representativity_menu['rangeboxes']['labels']

        try:
            data_availability_lower_bounds = []
            for var_ii, var in enumerate(active_data_availablity_vars):
                data_availability_lower_bounds.append(
                    np.float32(self.read_instance.representativity_menu['rangeboxes']['current_lower'][var_ii]))
        # if any of the fields are not numbers, return from function
        except ValueError:
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
                period_inds = np.arange(
                    self.read_instance.data_in_memory_filtered['observations'][species].shape[
                        1])
                # daily variable?
                if period == 'daily':
                    period_inds_split = np.array_split(period_inds,
                                                       [24 * i for i in range(1, int(np.ceil(len(period_inds) / 24)))])
                # monthly variable?
                elif period == 'monthly':
                    period_inds_split = np.array_split(period_inds, np.cumsum(self.read_instance.datareader.N_inds_per_month))
                # whole record variable?
                else:
                    period_inds_split = [period_inds]

                # iterate through indices associated with periodic chunks for current period
                for period_inds in period_inds_split:
                    if len(period_inds) > 0:
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

        # validate fields before filtering
        if not self.validate_values():
            return

        # iterate through all metadata
        for meta_var in self.read_instance.metadata_vars_to_read:

            if meta_var == 'lat':
                meta_var = 'latitude'
            elif meta_var == 'lon':
                meta_var = 'longitude'

            metadata_type = self.read_instance.standard_metadata[meta_var]['metadata_type']
            metadata_data_type = self.read_instance.standard_metadata[meta_var]['data_type']

            # handle non-numeric metadata
            if metadata_data_type == np.object:
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
            # handle numeric metadata
            else:
                meta_var_index = self.read_instance.metadata_menu[metadata_type]['rangeboxes']['labels'].index(meta_var)
                current_lower = np.float32(
                    self.read_instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index])
                current_upper = np.float32(
                    self.read_instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index])

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

    def validate_values(self):
        """Validates that field inserted by user is float"""

        # iterate through all metadata
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
                        print("Error in metadata fields. The field of '{}' "
                              "should be numeric, \n{}".format(meta_var, str(e)))
                    else:
                        QtWidgets.QMessageBox.critical(self.read_instance, "Error in metadata fields",
                                                       "The field of '{}' should be numeric, \n{}"
                                                       .format(meta_var, str(e)),
                                                       QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.NoButton)
                    return False

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
