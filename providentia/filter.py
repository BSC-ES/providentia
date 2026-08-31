""" Class for filtering of observational and model data """

import ast
import copy
import sys

import numpy as np
import pandas as pd
import yaml

from providentia.auxiliar import (
    CURRENT_PATH,
    join,
    get_conversion_factor,
    get_standard_parameters_by_speci,
)
from .calculate import Stats, ModBias
from .configuration import split_options
from .fields_menus import update_metadata_fields
from .plot_aux import update_plotting_parameters
from .statistics import merge_forecast_days, get_z_statistic_info, exceedance_lim
from .warnings_prv import show_message

PROVIDENTIA_ROOT = "/".join(CURRENT_PATH.split("/")[:-1])


class DataFilter:
    """Class for filtering of observational and model data"""

    def __init__(self, read_instance):
        """
        Initialise the DataFilter instance.

        Parameters
        ----------
        read_instance : object
            An instance of the data-reading class containing the data arrays and labels to be filtered.
        """

        self.read_instance = read_instance

        # get index of observational data array
        self.obs_index = self.read_instance.data_labels.index(
            self.read_instance.observations_data_label
        )

        # apply filtering
        self.filter_all()

    def filter_all(self):
        """
        Sequentially executes all data filtering procedures.
        """

        self.reset_data_filter()
        self.order_stations()
        self.filter_by_species()
        self.filter_data_limits()
        self.filter_by_period()
        self.filter_by_data_availability()
        self.filter_by_metadata()
        self.filter_extreme_stations()
        self.apply_calibration_factor()
        self.forecast_daily_switch()
        self.temporally_colocate_data()
        self.get_valid_stations()
        self.cap_stations()

        # save time index after filter in case it is overwritten later and can be reset
        self.read_instance.time_index_after_filter = copy.deepcopy(
            self.read_instance.time_index
        )

    def reset_data_filter(self):
        """
        Resets data arrays and temporal indices to their original, unfiltered state.
        """

        # reset data in memory
        # (plain .copy() is safe and much faster than deepcopy here: data_in_memory
        # values are pure float32 ndarrays, with no nested mutable objects to deep-copy)
        self.read_instance.data_in_memory_filtered = {
            networkspeci: data_array.copy()
            for networkspeci, data_array in self.read_instance.data_in_memory.items()
        }

        # reset metadata in memory
        # (plain .copy() is safe here too: the only object-dtype metadata fields hold
        # immutable strings, so a shallow array copy is exactly as independent as deepcopy)
        self.read_instance.metadata_in_memory_filtered = {
            networkspeci: metadata_array.copy()
            for networkspeci, metadata_array in self.read_instance.metadata_in_memory.items()
        }

        # reset basic metadata order
        if self.read_instance.station_order_inds_invert is not None:
            for networkspeci in self.read_instance.networkspecies:
                self.read_instance.station_references[networkspeci] = (
                    self.read_instance.station_references[networkspeci][self.read_instance.station_order_inds_invert]
                )
                self.read_instance.station_names[networkspeci] = (
                    self.read_instance.station_names[networkspeci][self.read_instance.station_order_inds_invert]
                )
                self.read_instance.station_longitudes[networkspeci] = (
                    self.read_instance.station_longitudes[networkspeci][self.read_instance.station_order_inds_invert]
                )
                self.read_instance.station_latitudes[networkspeci] = (
                    self.read_instance.station_latitudes[networkspeci][self.read_instance.station_order_inds_invert]
                )
                self.read_instance.station_measurement_altitudes[networkspeci] = (
                    self.read_instance.station_measurement_altitudes[networkspeci][self.read_instance.station_order_inds_invert]
                )
            self.read_instance.station_order_inds_invert = None
        
        self.read_instance.temporal_colocation_nans = {}
        self.read_instance.valid_station_inds = {}
        self.read_instance.valid_station_inds_temporal_colocation = {}

        # reset time index to original time array
        self.read_instance.time_index = self.read_instance.time_array

        # if are reading daily or combined forecast data restore data labels, models and plotting params
        # to how they were upon read for dashboard mode
        if (
            (self.read_instance.daily_forecast)
            or (self.read_instance.combined_forecast)
        ) & (self.read_instance.mode not in ["report", "library"]):
            self.read_instance.data_labels = copy.deepcopy(
                self.read_instance.original_data_labels
            )
            self.read_instance.data_labels_raw = copy.deepcopy(
                self.read_instance.original_data_labels_raw
            )
            self.read_instance.experiments = copy.deepcopy(
                self.read_instance.original_models
            )
            self.read_instance.plotting_params = copy.deepcopy(
                self.read_instance.original_plotting_params
            )

    def order_stations(self):
        """
        Orders stations based on a specified metadata variable and direction (ascending or descending).
        The station order is defined in the configuration file as 'station_order' in the format 'metadata_variable||direction'.
        """

        if self.read_instance.station_order:

            # check station_order is valid metadata variable or random
            if ((self.read_instance.station_order not in self.read_instance.metadata_vars_to_read)
                and (self.read_instance.station_order != 'random')):
                error = ( 
                    "Error: station_order metadata variable '{}' does not exist.".format(
                    self.read_instance.station_order
                    )
                )
                self.read_instance.logger.error(error)
                sys.exit(1)

            for networkspeci in self.read_instance.networkspecies:

                # get indices of stations ordered by station_order variable
                # this handles ascending and descending order for both numerical and string metadata variables,
                # as well as random sort
                if self.read_instance.station_order == 'random':
                    station_order_inds = np.random.permutation(
                        self.read_instance.metadata_in_memory_filtered[networkspeci].shape[0]
                    )
                else:
                    # metadata array is [station, month]; a station can be NaN for a given
                    # month if it wasn't active then, so take the first valid (non-NaN)
                    # value per station across months, rather than just column 0
                    station_metadata = self.read_instance.metadata_in_memory_filtered[networkspeci][
                        self.read_instance.station_order
                    ]
                    n_stations = station_metadata.shape[0]

                    # mask of valid (non-missing) entries per station/month
                    valid_mask = ~(pd.isna(station_metadata) | np.vectorize(lambda x: isinstance(x, str) and x.strip().lower() == 'nan')(station_metadata))
                    has_any_valid = valid_mask.any(axis=1)

                    # index of first valid month per station (argmax finds first True;
                    # returns 0 for all-NaN rows, but we overwrite those with NaN below)
                    first_valid_month = np.argmax(valid_mask, axis=1)

                    first_valid_values = station_metadata[np.arange(n_stations), first_valid_month]
                    # cast to object so we can safely assign np.nan even if metadata is string dtype
                    first_valid_values = first_valid_values.astype(object)
                    first_valid_values[~has_any_valid] = np.nan

                    # sort, always keeping NaNs at the end regardless of direction
                    ascending = self.read_instance.station_order_direction != "descending"
                    sorted_series = pd.Series(first_valid_values).sort_values(ascending=ascending, na_position="last")
                    station_order_inds = sorted_series.index.to_numpy()

                # get indices to reorder stations to how they originally were
                self.read_instance.station_order_inds_invert = np.argsort(station_order_inds)

                # reorder data_in_memory_filtered and metadata_in_memory_filtered according to station_order_inds
                self.read_instance.data_in_memory_filtered[networkspeci] = (
                    self.read_instance.data_in_memory_filtered[networkspeci][
                        :, station_order_inds, :
                    ]
                )
                for meta_var in self.read_instance.metadata_vars_to_read:
                    self.read_instance.metadata_in_memory_filtered[networkspeci][meta_var] = (
                        self.read_instance.metadata_in_memory_filtered[networkspeci][meta_var][
                            station_order_inds, :
                        ]
                    )

                # reorder basic metadata according to station_order_inds
                self.read_instance.station_references[networkspeci] = (
                    self.read_instance.station_references[networkspeci][station_order_inds]
                )
                self.read_instance.station_names[networkspeci] = (
                    self.read_instance.station_names[networkspeci][station_order_inds]
                )
                self.read_instance.station_longitudes[networkspeci] = (
                    self.read_instance.station_longitudes[networkspeci][station_order_inds]
                )
                self.read_instance.station_latitudes[networkspeci] = (
                    self.read_instance.station_latitudes[networkspeci][station_order_inds]
                )
                self.read_instance.station_measurement_altitudes[networkspeci] = (
                    self.read_instance.station_measurement_altitudes[networkspeci][station_order_inds]
                )

    def filter_by_species(self):
        """
        Filters read species using one or more "filter species".

        For each filter species, one or more lower/upper limit rules are defined.
        Wherever the filter-species data satisfy a given rule OR are NaN, the corresponding
        elements in *all* read species are overwritten with that rule's fill value.

        Filtering is only performed if:
        - filter_species is set
        - spatial_colocation is active
        """

        # ------------------------------------------------------------------
        # Only run if filter species are configured AND spatial colocation is active
        # ------------------------------------------------------------------
        if not (
            self.read_instance.filter_species and self.read_instance.spatial_colocation
        ):
            return

        # Local bindings to reduce repeated deep attribute lookups
        ri = self.read_instance
        networkspecies = ri.networkspecies
        obs_index = self.obs_index
        data_in_memory_filtered = ri.data_in_memory_filtered
        filter_data_in_memory = ri.filter_data_in_memory

        # ------------------------------------------------------------------
        # Helper to parse a bound string such as:
        #   '>10', '>=10', '<5', '<=5', ':'
        #
        # Returns:
        #   is_open      : bool   True if bound is ':' (i.e. no bound)
        #   inclusive    : bool   True if '=' is present
        #   value        : float or None
        # ------------------------------------------------------------------
        def _parse_bound(bound_str):
            is_open = ":" in bound_str
            if is_open:
                return True, False, None

            inclusive = "=" in bound_str
            value = float(bound_str.replace(">", "").replace("<", "").replace("=", ""))
            return False, inclusive, value

        # ------------------------------------------------------------------
        # Iterate over each filter species
        # ------------------------------------------------------------------
        for filter_networkspeci, speci_all_limits in ri.filter_species.items():
            # Data for the current filter species
            # Shape is expected to match the "station x time" layout used below
            data = filter_data_in_memory[filter_networkspeci]

            # Cache NaN mask once per filter species (used for every limit)
            nan_mask = np.isnan(data)

            # ------------------------------------------------------------------
            # Apply each limit independently
            # ------------------------------------------------------------------
            for speci_limit in speci_all_limits:
                lower_limit = speci_limit[0]
                upper_limit = speci_limit[1]
                filter_species_fill_value = speci_limit[2]

                # Parse bounds
                lower_open, lower_inclusive, lower_value = _parse_bound(lower_limit)
                upper_open, upper_inclusive, upper_value = _parse_bound(upper_limit)

                # --------------------------------------------------------------
                # If BOTH bounds are open (':' and ':'), return immediately
                # --------------------------------------------------------------
                if lower_open and upper_open:
                    return

                # --------------------------------------------------------------
                # Build validity mask for current filter species and current bound rule
                # --------------------------------------------------------------
                if upper_open:
                    # Only lower bound is active
                    if lower_inclusive:
                        valid_inds_per_species = data >= lower_value
                    else:
                        valid_inds_per_species = data > lower_value

                elif lower_open:
                    # Only upper bound is active
                    if upper_inclusive:
                        valid_inds_per_species = data <= upper_value
                    else:
                        valid_inds_per_species = data < upper_value

                else:
                    # Both lower and upper bounds are active
                    if lower_inclusive:
                        lower_mask = data >= lower_value
                    else:
                        lower_mask = data > lower_value

                    if upper_inclusive:
                        upper_mask = data <= upper_value
                    else:
                        upper_mask = data < upper_value

                    # Equivalent to np.logical_and.reduce((lower_mask, upper_mask))
                    # but faster and simpler for just two conditions
                    valid_inds_per_species = lower_mask & upper_mask

                # --------------------------------------------------------------
                # NaNs in the filter species are always treated as valid-to-filter
                # --------------------------------------------------------------
                valid_inds_per_species |= nan_mask

                # --------------------------------------------------------------
                # Apply THIS RULE'S fill value immediately to all read species
                # --------------------------------------------------------------
                for networkspeci in networkspecies:
                    data_in_memory_filtered[networkspeci][
                        obs_index, valid_inds_per_species
                    ] = filter_species_fill_value

    def filter_data_limits(self):
        """
        Filter out data which exceed the lower/upper limits.
        """

        # iterate through networkspecies
        for networkspeci in self.read_instance.networkspecies:
            # get speci str
            speci = networkspeci.split("|")[1]

            # get lower/upper data bounds
            if self.read_instance.mode in ["report", "library"]:
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
                msg = "Data limit fields must be numeric."
                show_message(
                    self.read_instance, msg, from_conf=self.read_instance.from_conf
                )
                return

            # convert units of bounds from standard to actual before filtering
            standard_parameter_speci = get_standard_parameters_by_speci(
                speci, self.read_instance.ghost_version
            )
            initial_units = standard_parameter_speci["standard_units"]
            final_units = self.read_instance.measurement_units[speci]
            conversion_factor = get_conversion_factor(
                initial_units, final_units, standard_parameter_speci
            )
            if isinstance(conversion_factor, str):
                self.read_instance.logger.error(conversion_factor)
                sys.exit(1)
            lower_bound *= conversion_factor
            upper_bound *= conversion_factor

            # filter all observational/model data out of bounds of lower/upper limits
            arr = self.read_instance.data_in_memory_filtered[networkspeci]
            inds_out_of_bounds = (arr < lower_bound) | (arr > upper_bound)
            self.read_instance.data_in_memory_filtered[networkspeci][
                inds_out_of_bounds
            ] = np.nan

    def filter_by_period(self):
        """
        Filter data for selected periods (keeping or removing data, as defined).
        """

        # set appropriate data and variable name arrays
        data_array = self.read_instance.ghost_data_in_memory
        varname_array = self.read_instance.ghost_data_vars_to_read

        keeps, removes = [], []
        if self.read_instance.mode in ["report", "library"]:
            if hasattr(self.read_instance, "period"):
                keeps, removes = split_options(
                    self.read_instance, self.read_instance.period
                )
            else:
                # Get period if apply_filter is used in library mode
                if self.read_instance.period_menu["checkboxes"]["keep_selected"]:
                    keeps = self.read_instance.period_menu["checkboxes"][
                        "keep_selected"
                    ]
                if self.read_instance.period_menu["checkboxes"]["remove_selected"]:
                    removes = self.read_instance.period_menu["checkboxes"][
                        "remove_selected"
                    ]
        else:
            keeps = self.read_instance.period_menu["checkboxes"]["keep_selected"]
            removes = self.read_instance.period_menu["checkboxes"]["remove_selected"]

        available_values = [
            "Daytime",
            "Nighttime",
            "Weekday",
            "Weekend",
            "Spring",
            "Summer",
            "Autumn",
            "Winter",
        ]

        # filter/limit data for periods selected
        if len(keeps) > 0:
            for keep in keeps:
                if keep not in available_values:
                    msg = f"\nCannot filter data by period: {keep}, ignoring filter value. "
                    msg += f"Available filters: {available_values}"
                    self.read_instance.logger.info(msg)

            day_night_codes_to_keep = []
            if "Daytime" in keeps:
                day_night_codes_to_keep.append(0)
            if "Nighttime" in keeps:
                day_night_codes_to_keep.append(1)
            if len(day_night_codes_to_keep) == 1:
                if "hourly" in self.read_instance.resolution:
                    day_night_index = varname_array.index("day_night_code")
                    # iterate through network / species
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(
                            data_array[networkspeci][day_night_index, :, :],
                            day_night_codes_to_keep,
                            invert=True,
                        )
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            :, inds_to_screen
                        ] = np.nan

            weekday_weekend_codes_to_keep = []
            if "Weekday" in keeps:
                weekday_weekend_codes_to_keep.append(0)
            if "Weekend" in keeps:
                weekday_weekend_codes_to_keep.append(1)
            if len(weekday_weekend_codes_to_keep) == 1:
                if self.read_instance.resolution not in ["monthly", "annual"]:
                    weekday_weekend_index = varname_array.index("weekday_weekend_code")
                    # iterate through network / species
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(
                            data_array[networkspeci][weekday_weekend_index, :, :],
                            weekday_weekend_codes_to_keep,
                            invert=True,
                        )
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            :, inds_to_screen
                        ] = np.nan

            season_codes_to_keep = []
            if "Spring" in keeps:
                season_codes_to_keep.append(0)
            if "Summer" in keeps:
                season_codes_to_keep.append(1)
            if "Autumn" in keeps:
                season_codes_to_keep.append(2)
            if "Winter" in keeps:
                season_codes_to_keep.append(3)
            if (len(season_codes_to_keep) > 0) & (len(season_codes_to_keep) < 4):
                if self.read_instance.resolution != "annual":
                    season_index = varname_array.index("season_code")
                    # iterate through network / species
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(
                            data_array[networkspeci][season_index, :, :],
                            season_codes_to_keep,
                            invert=True,
                        )
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            :, inds_to_screen
                        ] = np.nan

        if len(removes) > 0:
            for remove in removes:
                if remove not in available_values:
                    msg = f"\nCannot filter data by period: {remove}, ignoring filter value. "
                    msg += f"Available filters: {available_values}"
                    self.read_instance.logger.info(msg)

            day_night_codes_to_remove = []
            if "Daytime" in removes:
                day_night_codes_to_remove.append(0)
            if "Nighttime" in removes:
                day_night_codes_to_remove.append(1)
            if len(day_night_codes_to_remove) > 0:
                if "hourly" in self.read_instance.resolution:
                    day_night_index = varname_array.index("day_night_code")
                    # iterate through network / species
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(
                            data_array[networkspeci][day_night_index, :, :],
                            day_night_codes_to_remove,
                        )
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            :, inds_to_screen
                        ] = np.nan

            weekday_weekend_codes_to_remove = []
            if "Weekday" in removes:
                weekday_weekend_codes_to_remove.append(0)
            if "Weekend" in removes:
                weekday_weekend_codes_to_remove.append(1)
            if len(weekday_weekend_codes_to_remove) > 0:
                if self.read_instance.resolution not in ["monthly", "annual"]:
                    weekday_weekend_index = varname_array.index("weekday_weekend_code")
                    # iterate through network / species
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(
                            data_array[networkspeci][weekday_weekend_index, :, :],
                            weekday_weekend_codes_to_remove,
                        )
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            :, inds_to_screen
                        ] = np.nan

            season_codes_to_remove = []
            if "Spring" in removes:
                season_codes_to_remove.append(0)
            if "Summer" in removes:
                season_codes_to_remove.append(1)
            if "Autumn" in removes:
                season_codes_to_remove.append(2)
            if "Winter" in removes:
                season_codes_to_remove.append(3)
            if len(season_codes_to_remove) > 0:
                if self.read_instance.resolution != "annual":
                    season_index = varname_array.index("season_code")
                    # iterate through network / species
                    for networkspeci in self.read_instance.networkspecies:
                        inds_to_screen = np.isin(
                            data_array[networkspeci][season_index, :, :],
                            season_codes_to_remove,
                        )
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            :, inds_to_screen
                        ] = np.nan

    def filter_by_data_availability(self):
        """
        Filters data based on native and dynamically calculated data completeness thresholds.
        """

        # get set variables names representing percentage data availability (native and non-native)
        if float(".".join(self.read_instance.ghost_version.split(".")[:2])) < 1.6:
            active_data_availablity_vars = self.read_instance.coverage_menu[
                "rangeboxes"
            ]["map_vars_old"]
        else:
            active_data_availablity_vars = self.read_instance.coverage_menu[
                "rangeboxes"
            ]["map_vars"]

        try:
            data_availability_lower_bounds = []
            for var_ii, var in enumerate(active_data_availablity_vars):
                data_availability_lower_bounds.append(
                    np.float32(
                        self.read_instance.coverage_menu["rangeboxes"]["current_lower"][
                            var_ii
                        ]
                    )
                )
        # if any of the fields are not numbers, return from function
        except ValueError:
            msg = "Data availability fields must be numeric."
            show_message(
                self.read_instance, msg, from_conf=self.read_instance.from_conf
            )
            return

        # filter observations by native percentage data availability variables (only GHOST data)
        if self.read_instance.reading_ghost:
            # get appropriate data and variable name arrays
            data_array = self.read_instance.ghost_data_in_memory
            varname_array = self.read_instance.ghost_data_vars_to_read

            # iterate through data availability variables
            for var_ii, var in enumerate(active_data_availablity_vars):
                # variable is GHOST native?
                if "native" in var:
                    var_index = varname_array.index(var)

                    # iterate through network / species
                    for networkspeci in self.read_instance.networkspecies:
                        # max gap variable?
                        if "max_gap" in var:
                            # bound is < 100?:
                            if data_availability_lower_bounds[var_ii] < 100:
                                inds_to_screen = (
                                    data_array[networkspeci][var_index, :, :]
                                    > data_availability_lower_bounds[var_ii]
                                )
                                self.read_instance.data_in_memory_filtered[
                                    networkspeci
                                ][self.obs_index, inds_to_screen] = np.nan
                        # data coverage variable?
                        else:
                            # bound is > 0?
                            if data_availability_lower_bounds[var_ii] > 0:
                                inds_to_screen = (
                                    data_array[networkspeci][var_index, :, :]
                                    < data_availability_lower_bounds[var_ii]
                                )
                                self.read_instance.data_in_memory_filtered[
                                    networkspeci
                                ][self.obs_index, inds_to_screen] = np.nan

        # filter observations and model data by non-native percentage data availability variables
        # (calculated on the fly)
        for var_ii, var in enumerate(active_data_availablity_vars):
            if "native" not in var:
                # max gap variable?
                if "max_gap" in var:
                    # bound is == 100?
                    if data_availability_lower_bounds[var_ii] == 100:
                        continue
                # data coverage variable?
                else:
                    # bound == 0?
                    if data_availability_lower_bounds[var_ii] == 0:
                        continue

                # get period associate with variable
                period = var.split("_")[0]
                period_inds = np.arange(len(self.read_instance.time_array))

                # daily variable?
                if period == "daily":
                    period_inds_split = np.array_split(
                        period_inds,
                        [24 * i for i in range(1, int(np.ceil(len(period_inds) / 24)))],
                    )
                # monthly variable?
                elif period == "monthly":
                    period_inds_split = np.array_split(
                        period_inds, np.cumsum(self.read_instance.N_inds_per_yearmonth)
                    )
                # annual variable?
                elif period == "annual":
                    # if annual variable, then split indices into whole record variable
                    period_inds_split = np.array_split(
                        period_inds, np.cumsum(self.read_instance.N_inds_per_year)
                    )
                # whole record variable?
                else:
                    period_inds_split = [period_inds]

                # iterate through indices associated with periodic chunks for current period
                for period_inds in period_inds_split:
                    if len(period_inds) > 0:
                        # iterate through networkspecies
                        for networkspeci in self.read_instance.networkspecies:
                            # data coverage variable?
                            data_availability_percent = (
                                Stats.calculate_data_avail_fraction(
                                    self.read_instance.data_in_memory_filtered[
                                        networkspeci
                                    ][
                                        self.obs_index,
                                        :,
                                        period_inds[0] : period_inds[-1] + 1,
                                    ]
                                )
                            )
                            inds_to_screen = np.where(
                                data_availability_percent
                                < data_availability_lower_bounds[var_ii]
                            )[0]
                            self.read_instance.data_in_memory_filtered[networkspeci][
                                self.obs_index,
                                inds_to_screen[:, np.newaxis],
                                period_inds[np.newaxis, :],
                            ] = np.nan

    def filter_by_metadata(self):
        """
        Filters data based on categorical and numerical station metadata.
        """

        # iterate through metadata in memory
        for meta_var in self.read_instance.metadata_vars_to_read:
            # validate field before filtering
            if not self.validate_values(meta_var):
                # go to next variable if filter cannot be applied
                continue

            if meta_var == "lat":
                meta_var = "latitude"
            elif meta_var == "lon":
                meta_var = "longitude"

            metadata_type = self.read_instance.standard_metadata[meta_var][
                "metadata_type"
            ]
            metadata_data_type = self.read_instance.standard_metadata[meta_var][
                "data_type"
            ]

            # handle non-numeric metadata
            if metadata_data_type is object:
                # iterate through networkspecies
                for networkspeci in self.read_instance.networkspecies:
                    # if any of the keep checkboxes are selected, filter out data by fields that have not been selected
                    current_keep = self.read_instance.metadata_menu[metadata_type][
                        meta_var
                    ]["checkboxes"]["keep_selected"]
                    if len(current_keep) > 0:
                        invalid_keep = np.repeat(
                            np.isin(
                                self.read_instance.metadata_in_memory_filtered[networkspeci][
                                    meta_var
                                ][:, :],
                                current_keep,
                                invert=True,
                            ),
                            self.read_instance.N_inds_per_yearmonth,
                            axis=1,
                        )
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            :, invalid_keep
                        ] = np.nan
                        for current_keep_value in current_keep:
                            invalid_keep_value = np.repeat(
                                np.isin(
                                    self.read_instance.metadata_in_memory_filtered[networkspeci][
                                        meta_var
                                    ][:, :],
                                    current_keep_value,
                                    invert=True,
                                ),
                                self.read_instance.N_inds_per_yearmonth,
                                axis=1,
                            )
                            if np.all(invalid_keep_value):
                                available_values = np.unique(
                                    self.read_instance.metadata_in_memory_filtered[networkspeci][
                                        meta_var
                                    ].astype(str)
                                )
                                msg = f"Cannot filter by {meta_var}: {current_keep_value} in {networkspeci} data, ignoring filter value. "
                                msg += (
                                    f"Available filters: {available_values.tolist()}\n"
                                )
                                self.read_instance.logger.info(msg)

                    # if any of the remove checkboxes have been selected, filter out data by these selected fields
                    current_remove = self.read_instance.metadata_menu[metadata_type][
                        meta_var
                    ]["checkboxes"]["remove_selected"]
                    if len(current_remove) > 0:
                        invalid_remove = np.repeat(
                            np.isin(
                                self.read_instance.metadata_in_memory_filtered[networkspeci][
                                    meta_var
                                ][:, :],
                                current_remove,
                            ),
                            self.read_instance.N_inds_per_yearmonth,
                            axis=1,
                        )
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            :, invalid_remove
                        ] = np.nan
                        for current_remove_value in current_remove:
                            invalid_remove_value = np.repeat(
                                np.isin(
                                    self.read_instance.metadata_in_memory_filtered[networkspeci][
                                        meta_var
                                    ][:, :],
                                    current_remove_value,
                                ),
                                self.read_instance.N_inds_per_yearmonth,
                                axis=1,
                            )
                            if not np.any(invalid_remove_value):
                                available_values = np.unique(
                                    self.read_instance.metadata_in_memory_filtered[networkspeci][
                                        meta_var
                                    ].astype(str)
                                )
                                msg = f"Cannot filter by {meta_var}: {current_remove_value} in {networkspeci} data, ignoring filter value. "
                                msg += (
                                    f"Available filters: {available_values.tolist()}\n"
                                )
                                self.read_instance.logger.info(msg)

            # handle numeric metadata
            else:
                meta_var_index = self.read_instance.metadata_menu[metadata_type][
                    "rangeboxes"
                ]["labels"].index(meta_var)
                current_lower = np.float32(
                    self.read_instance.metadata_menu[metadata_type]["rangeboxes"][
                        "current_lower"
                    ][meta_var_index]
                )
                current_upper = np.float32(
                    self.read_instance.metadata_menu[metadata_type]["rangeboxes"][
                        "current_upper"
                    ][meta_var_index]
                )

                # get array with selected filters (those with the check Apply on)
                current_apply = self.read_instance.metadata_menu[metadata_type][
                    "rangeboxes"
                ]["apply_selected"]

                if len(current_apply) > 0:
                    # apply bounds and remove nans if variable has been selected
                    if meta_var in current_apply:
                        # iterate through networkspecies
                        for networkspeci in self.read_instance.networkspecies:
                            # if current lower value is non-NaN, then filter out data with metadata < current lower value
                            if not pd.isnull(current_lower):
                                lower_default = np.float32(
                                    self.read_instance.metadata_menu[metadata_type][
                                        "rangeboxes"
                                    ]["lower_default"][meta_var_index]
                                )
                                if current_lower >= lower_default:
                                    invalid_below = np.repeat(
                                        self.read_instance.metadata_in_memory_filtered[
                                            networkspeci
                                        ][meta_var][:, :]
                                        < current_lower,
                                        self.read_instance.N_inds_per_yearmonth,
                                        axis=1,
                                    )
                                    self.read_instance.data_in_memory_filtered[
                                        networkspeci
                                    ][:, invalid_below] = np.nan

                            # if current upper < than the maximum extent, then filter out
                            # data with metadata > current upper value (if this is numeric)
                            # if current upper value is non-NaN, then filter out data with metadata > current upper value
                            if not pd.isnull(current_upper):
                                upper_default = np.float32(
                                    self.read_instance.metadata_menu[metadata_type][
                                        "rangeboxes"
                                    ]["upper_default"][meta_var_index]
                                )
                                if current_upper <= upper_default:
                                    invalid_above = np.repeat(
                                        self.read_instance.metadata_in_memory_filtered[
                                            networkspeci
                                        ][meta_var][:, :]
                                        > current_upper,
                                        self.read_instance.N_inds_per_yearmonth,
                                        axis=1,
                                    )
                                    self.read_instance.data_in_memory_filtered[
                                        networkspeci
                                    ][:, invalid_above] = np.nan

                            # remove nans
                            invalid_nan = np.repeat(
                                pd.isnull(
                                    self.read_instance.metadata_in_memory_filtered[networkspeci][
                                        meta_var
                                    ][:, :]
                                ),
                                self.read_instance.N_inds_per_yearmonth,
                                axis=1,
                            )
                            self.read_instance.data_in_memory_filtered[networkspeci][
                                :, invalid_nan
                            ] = np.nan

    def validate_values(self, meta_var):
        """
        Ensures user-provided metadata range boundaries are valid floating-point numbers.

        Parameters
        ----------
        meta_var : str
            The name of the metadata variable to validate.

        Returns
        -------
        is_valid : bool
            Returns True if the metadata values are numeric or non-numeric objects, False otherwise.
        """

        metadata_type = self.read_instance.standard_metadata[meta_var]["metadata_type"]
        metadata_data_type = self.read_instance.standard_metadata[meta_var]["data_type"]

        if metadata_data_type is not object:
            meta_var_index = self.read_instance.metadata_menu[metadata_type][
                "rangeboxes"
            ]["labels"].index(meta_var)
            try:
                np.float32(
                    self.read_instance.metadata_menu[metadata_type]["rangeboxes"][
                        "current_lower"
                    ][meta_var_index]
                )
                np.float32(
                    self.read_instance.metadata_menu[metadata_type]["rangeboxes"][
                        "current_upper"
                    ][meta_var_index]
                )
                return True
            except ValueError:
                msg = "Error in metadata fields. The field of '{}' should be numeric.".format(
                    meta_var
                )
                show_message(
                    self.read_instance, msg, from_conf=self.read_instance.from_conf
                )
                if self.read_instance.mode not in ["report", "library"]:
                    self.read_instance.metadata_menu[metadata_type]["rangeboxes"][
                        "apply_selected"
                    ].remove(meta_var)
                return False
        else:
            return True

    def filter_extreme_stations(self):
        """
        Filters out extreme stations based on set statistical thresholds.
        There can be multiple limit arguments for a statistic e.g. 'MB': ['<10','>20']
        There can also be limits per species e.g. 'RMSE': {'sconco3': ['<50.0', '>70.0'], 'sconco':[>100.0]}
        An absolute statistic can be set to be a bias statistic by adding '_bias' e.g.
        'p95_bias': ['<10','>20']

        If statistic is an absolute statistic, then only remove stations based on observations.
        If statistic is a bias statistic, then remove collection of all stations outside limits across all obs-mod
        comparsions.
        """

        # option to remove extreme stations set?
        if self.read_instance.remove_extreme_stations:
            # load yaml of defined stattistics limits
            remove_extreme_stations_fname = join(
                PROVIDENTIA_ROOT, "settings/remove_extreme_stations.yaml"
            )
            stat_defs = yaml.safe_load(open(remove_extreme_stations_fname))

            # get specific set of limits (if available)
            # throw wraning if not
            if self.read_instance.remove_extreme_stations not in stat_defs:
                msg = "'{}' not defined in '{}'. Not removing extreme stations.".format(
                    self.read_instance.remove_extreme_stations,
                    remove_extreme_stations_fname,
                )
                show_message(
                    self.read_instance, msg, from_conf=self.read_instance.from_conf
                )
                return
            else:
                stat_limits = stat_defs[self.read_instance.remove_extreme_stations]

            # loop through and calculate each statistic per station and remove stations outside statistical limits
            for zstat in stat_limits:
                # get list of statistical limits for specific stat
                stat_arguments = stat_limits[zstat]

                # if have a dict, the limits are specific per species, so limit for species
                speci_specific_limits = False
                if type(stat_arguments) is dict:
                    speci_specific_limits = True

                # determine if station is absolute or bias statistic
                (
                    zstat,
                    base_zstat,
                    z_statistic_type,
                    z_statistic_sign,
                    z_statistic_period,
                ) = get_z_statistic_info(zstat=zstat)

                # handle some possible errors

                # if have only observations data then cannot calculate bias statistic, so continue to next stat
                if z_statistic_sign == "bias":
                    if len(self.read_instance.data_labels) == 1:
                        msg = "Cannot remove extreme stations via calculation of '{}' as no model data has been read.".format(
                            zstat
                        )
                        show_message(
                            self.read_instance,
                            msg,
                            from_conf=self.read_instance.from_conf,
                        )
                        continue

                # if temporal_colocation is not active then cannot calculate ModBias statistic, so continue to next stat
                if z_statistic_type == "modbias":
                    if not self.read_instance.temporal_colocation:
                        msg = "Cannot remove extreme stations via calculation of '{}' as 'temporal_colocation' is not active.".format(
                            zstat
                        )
                        show_message(
                            self.read_instance,
                            msg,
                            from_conf=self.read_instance.from_conf,
                        )
                        continue

                # get dictionary containing necessary information for calculation of selected statistic
                if z_statistic_type == "basic":
                    stats_dict = self.read_instance.basic_stats[base_zstat]
                else:
                    stats_dict = self.read_instance.modbias_stats[base_zstat]

                # load default selected z statistic arguments for passing to statistical function
                function_arguments = stats_dict["arguments"]

                # iterate through network / species
                for ii, networkspeci in enumerate(self.read_instance.networkspecies):
                    # if stat is exceedances then add threshold value (if available)
                    if base_zstat == "Exceedances":
                        function_arguments["threshold"] = exceedance_lim(
                            self.read_instance, networkspeci
                        )

                    # get list of statistic limits specific for speci (if wanted)
                    if speci_specific_limits:
                        # if have not defined limits for speci, then throw warning and continue to next speci
                        if networkspeci not in stat_arguments:
                            msg = "No statistical limits defined for '{}' in '{}' section of '{}'. Not removing extreme stations for '{}'.".format(
                                networkspeci,
                                self.read_instance.remove_extreme_stations,
                                remove_extreme_stations_fname,
                                networkspeci,
                            )
                            show_message(
                                self.read_instance,
                                msg,
                                from_conf=self.read_instance.from_conf,
                            )
                            continue
                        else:
                            specific_stat_arguments = stat_arguments[networkspeci]
                    else:
                        specific_stat_arguments = copy.deepcopy(stat_arguments)

                    # calculate statistic per station and then compare statistic against limits
                    # set stations exceeding limits to NaN
                    # calculate basic stats
                    if (z_statistic_type == "basic") and (z_statistic_sign != "bias"):
                        data_array_a = self.read_instance.data_in_memory_filtered[
                            networkspeci
                        ][self.obs_index, :, :]
                        calc_stat = np.array(
                            getattr(Stats, stats_dict["function"])(
                                data_array_a, **function_arguments
                            )
                        )
                        non_finite_stat = ~np.isfinite(calc_stat)
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            self.obs_index, non_finite_stat, :
                        ] = np.nan
                        for specific_stat_argument in specific_stat_arguments:
                            invalid_stations = ast.literal_eval(
                                "calc_stat{}".format(specific_stat_argument)
                            )
                            self.read_instance.data_in_memory_filtered[networkspeci][
                                self.obs_index, invalid_stations, :
                            ] = np.nan
                    # calculate basic bias stats and modbias stats
                    else:
                        for data_label in self.read_instance.data_labels:
                            if data_label != self.read_instance.observations_data_label:
                                # get modid data label index
                                mod_data_index = self.read_instance.data_labels.index(
                                    data_label
                                )
                                data_array_a = (
                                    self.read_instance.data_in_memory_filtered[
                                        networkspeci
                                    ][self.obs_index, :, :]
                                )
                                data_array_b = (
                                    self.read_instance.data_in_memory_filtered[
                                        networkspeci
                                    ][mod_data_index, :, :]
                                )
                                # calculate basic bias stats
                                if (z_statistic_type == "basic") and (
                                    z_statistic_sign == "bias"
                                ):
                                    statistic_a = np.array(
                                        getattr(Stats, stats_dict["function"])(
                                            data_array_a, **function_arguments
                                        )
                                    )
                                    statistic_b = np.array(
                                        getattr(Stats, stats_dict["function"])(
                                            data_array_b, **function_arguments
                                        )
                                    )
                                    calc_stat = statistic_b - statistic_a
                                # calculate modbias stats
                                elif z_statistic_type == "modbias":
                                    calc_stat = np.array(
                                        getattr(ModBias, stats_dict["function"])(
                                            **{
                                                **function_arguments,
                                                **{
                                                    "obs": data_array_a,
                                                    "mod": data_array_b,
                                                },
                                            }
                                        )
                                    )
                                non_finite_stat = ~np.isfinite(calc_stat)
                                self.read_instance.data_in_memory_filtered[
                                    networkspeci
                                ][self.obs_index, non_finite_stat, :] = np.nan
                                for specific_stat_argument in specific_stat_arguments:
                                    invalid_stations = eval(
                                        "calc_stat{}".format(specific_stat_argument)
                                    )
                                    self.read_instance.data_in_memory_filtered[
                                        networkspeci
                                    ][self.obs_index, invalid_stations, :] = np.nan

    def apply_calibration_factor(self):
        """
        Adjusts model data arrays through basic arithmetic operations using user-defined calibration factors.
        """

        if self.read_instance.calibration_factor:
            # iterate through networkspecies
            for networkspeci_ii, networkspeci in enumerate(
                self.read_instance.networkspecies
            ):
                # remove observations from data labels
                relevant_data_labels = copy.deepcopy(self.read_instance.data_labels)
                relevant_data_labels.remove(self.read_instance.observations_data_label)

                # get calibration factor per model
                for data_label_ii, data_label in enumerate(relevant_data_labels):
                    # get calibration factor per model
                    if isinstance(self.read_instance.calibration_factor, dict):
                        mod_label = list(self.read_instance.experiments.keys())[
                            list(self.read_instance.experiments.values()).index(
                                data_label
                            )
                        ]
                        if mod_label not in self.read_instance.calibration_factor:
                            msg = f"No calibration factor applied to model {mod_label}."
                            self.read_instance.logger.info(msg)
                            continue
                        calibration_factor = self.read_instance.calibration_factor[
                            mod_label
                        ]
                    else:
                        calibration_factor = self.read_instance.calibration_factor

                    # get calibration factor per networkspeci
                    if (len(self.read_instance.networkspecies) > 1) and (
                        "," in calibration_factor
                    ):
                        calibration_factor = calibration_factor.split(",")[
                            networkspeci_ii
                        ]

                    msg = "Applying calibration factor: "
                    msg += "{0} in {1} model".format(calibration_factor, data_label)
                    self.read_instance.logger.info(msg)

                    # apply calibration factor
                    if (
                        calibration_factor.count("*") == 1
                        and calibration_factor[0] == "*"
                    ):
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            data_label_ii + 1, :, :
                        ] *= float(calibration_factor.replace("*", ""))
                    elif (
                        calibration_factor.count("/") == 1
                        and calibration_factor[0] == "/"
                    ):
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            data_label_ii + 1, :, :
                        ] /= float(calibration_factor.replace("/", ""))
                    elif (
                        calibration_factor.count("-") == 1
                        and calibration_factor[0] == "-"
                    ):
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            data_label_ii + 1, :, :
                        ] -= float(calibration_factor.replace("-", ""))
                    elif (
                        calibration_factor.count("+") == 1
                        and calibration_factor[0] == "+"
                    ):
                        self.read_instance.data_in_memory_filtered[networkspeci][
                            data_label_ii + 1, :, :
                        ] += float(calibration_factor)
                    else:
                        error = f"Error: Invalid format '{calibration_factor}' in calibration factor. Accepted formats are: '+num', '-num', '*num', or '/num'."
                        self.read_instance.logger.error(error)
                        sys.exit(1)

    def forecast_daily_switch(self):
        """
        Restructures forecast data by merging separate forecast days into a single dimension and updating metadata.
        """

        # Exit early if neither daily_forecast nor combined_forecast are active
        if (not self.read_instance.daily_forecast) and (
            not self.read_instance.combined_forecast
        ):
            return

        # Extract base labels by removing '-daily' and '-combined' suffixes from data_labels
        base_data_labels = [
            label.split("-daily")[0].split("-combined")[0]
            for label in self.read_instance.data_labels
        ]
        # Keep unique base labels while preserving order
        unique_base_data_labels = list(dict.fromkeys(base_data_labels))

        # Pass 1: Determine the maximum number of forecast days across all labels
        # and collect which forecast days are actually active
        self.read_instance.max_forecast_days = 0
        self.read_instance.active_forecast_days = []
        for networkspeci in self.read_instance.networkspecies:
            for base_data_label in base_data_labels:
                if base_data_label != self.read_instance.observations_data_label:
                    current_count = 0
                    for data_label in self.read_instance.data_labels:
                        if data_label.startswith(base_data_label):
                            current_count += 1
                            # Check if this label has forecast indices for the current network
                            if (
                                data_label
                                in self.read_instance.forecast_indices_per_data_label[
                                    networkspeci
                                ][base_data_label]
                            ):
                                forecast_day = (
                                    self.read_instance.forecast_indices_per_data_label[
                                        networkspeci
                                    ][base_data_label][data_label]
                                    + 1
                                )
                                # Collect active forecast days, avoiding duplicates
                                if (
                                    forecast_day
                                    not in self.read_instance.active_forecast_days
                                ):
                                    self.read_instance.active_forecast_days.append(
                                        forecast_day
                                    )

                    # Update max_forecast_days if current label has more forecast days
                    if current_count > self.read_instance.max_forecast_days:
                        self.read_instance.max_forecast_days = current_count

        # Sort the active forecast days for consistent ordering
        self.read_instance.active_forecast_days = np.sort(
            self.read_instance.active_forecast_days
        )

        # Pass 2: Rebuild the in-memory data array to merge forecast days separated as different models to same dimension (tiled)
        for networkspeci in self.read_instance.networkspecies:
            # merge forecast days as different models to same dimension (tiled)
            new_data_in_memory = merge_forecast_days(
                self.read_instance,
                networkspeci,
                self.read_instance.data_labels,
                unique_base_data_labels,
                self.read_instance.data_in_memory_filtered[networkspeci],
            )

            # Replace old filtered data with the newly built array
            self.read_instance.data_in_memory_filtered[
                networkspeci
            ] = new_data_in_memory

        # Rebuild the global time index by tiling the original time array for each forecast day
        self.read_instance.time_index = pd.DatetimeIndex(
            np.tile(self.read_instance.time_array, self.read_instance.max_forecast_days)
        )

        # Update model labels and data_labels to reflect daily or combined forecasts
        data_labels_to_remove = []
        data_labels_to_add = set()
        new_models = {}

        for mod_raw, mod in self.read_instance.experiments.items():
            if ("-daily" in mod) or ("-combined" in mod):
                if self.read_instance.daily_forecast:
                    new_mod = f"{mod.split('-daily')[0]}-daily"
                    new_mod_raw = f"{mod_raw.split('-daily')[0]}-daily"
                elif self.read_instance.combined_forecast:
                    new_mod = f"{mod.split('-combined')[0]}-combined"
                    new_mod_raw = f"{mod_raw.split('-combined')[0]}-combined"
                if new_mod_raw not in new_models:
                    new_models[new_mod_raw] = new_mod
                data_labels_to_remove.append(mod)
                data_labels_to_add.add(new_mod)
            else:
                # Keep models without '-daily' or '-combined' unchanged
                new_models[mod_raw] = mod

        # Update the models dictionary
        self.read_instance.experiments = dict(new_models)
        # Rebuild the data_labels and raw data_labels including observations first
        self.read_instance.data_labels = [
            self.read_instance.observations_data_label
        ] + list(self.read_instance.experiments.values())
        self.read_instance.data_labels_raw = [
            self.read_instance.observations_data_label
        ] + list(self.read_instance.experiments.keys())

        # Update plotting parameters to reflect changes in labels and daily forecasts
        data_labels_to_remove = [
            mod
            for mod in data_labels_to_remove
            if mod in list(self.read_instance.plotting_params.keys())
        ]
        data_labels_to_add = [
            mod
            for mod in list(data_labels_to_add)
            if mod not in list(self.read_instance.plotting_params.keys())
        ]
        update_plotting_parameters(
            self.read_instance,
            data_labels_to_remove=data_labels_to_remove,
            data_labels_to_add=data_labels_to_add,
            daily_forecast=True,
        )

    def temporally_colocate_data(self):
        """
        Temporally colocate observational and model data.
        If spatial colocation is active, then data is also temporally colocated across all network / species,
        otherwise it is done independently per network / species.
        This in reality means storing the indices for the temporal colocation.
        """

        # iterate through network / species
        for ii, networkspeci in enumerate(self.read_instance.networkspecies):
            # initialise arrays to determine where have NaNs
            if (ii == 0) or (not self.read_instance.spatial_colocation):
                # create array for finding instances where have 0 valid values across all observations
                # initialise as being all False (i.e. non-NaN), set as True on the occasion there is a NaN in the observations
                obs_all_nan = np.full(
                    self.read_instance.data_in_memory_filtered[networkspeci][
                        self.obs_index, :, :
                    ].shape,
                    False,
                )

                # create array for finding instances where have 0 valid values across all models
                # initialise as being all False (i.e. non-NaN), set as True on the occasion there is a NaN in an model
                mods_all_nan = np.full(obs_all_nan.shape, False)

            # get all instances observations is NaN
            nan_obs = np.isnan(
                self.read_instance.data_in_memory_filtered[networkspeci][
                    self.obs_index, :, :
                ]
            )

            # update obs_all_nan array, making True all instances where have NaNs
            # if all observations are nan then do not update
            if not np.all(nan_obs):
                obs_all_nan = np.any([obs_all_nan, nan_obs], axis=0)

            # iterate through model data arrays in data in memory dictionary
            # save indices for colocation with observations
            for model in self.read_instance.experiments.values():
                # get modid data label index
                mod_data_index = self.read_instance.data_labels.index(model)

                # get all instances model is NaN
                nan_mod = np.isnan(
                    self.read_instance.data_in_memory_filtered[networkspeci][
                        mod_data_index, :, :
                    ]
                )

                # update mods_all_nan array, making True all instances where have NaNs
                # if all model values are nan then do not update for that model
                if not np.all(nan_mod):
                    mods_all_nan = np.any([mods_all_nan, nan_mod], axis=0)

            # if spatial colocation is not active,
            # get indices where one of observations and models per network /species is NaN
            if not self.read_instance.spatial_colocation:
                self.read_instance.temporal_colocation_nans[networkspeci] = np.any(
                    [obs_all_nan, mods_all_nan], axis=0
                )

        # if spatial colocation is active,
        # get indices where one of observations and models across networks / species is NaN
        if self.read_instance.spatial_colocation:
            for networkspeci in self.read_instance.networkspecies:
                self.read_instance.temporal_colocation_nans[networkspeci] = np.any(
                    [obs_all_nan, mods_all_nan], axis=0
                )

    def get_valid_stations(self):
        """
        Get valid station indices before and after all filtering has been performed.
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
                    # (.copy() suffices: pure float32 data, and it is mutated in place
                    # below via the temporal colocation NaN mask, so a copy is required,
                    # but does not need to be a deepcopy)
                    obs_data = self.read_instance.data_in_memory_filtered[networkspeci][
                        self.read_instance.data_labels.index(data_label), :, :
                    ].copy()

                    # get absolute data availability number per station in observational data array
                    if obs_data.size == 0:
                        station_data_availability_number = np.array([])
                    else:
                        station_data_availability_number = (
                            Stats.calculate_data_avail_number(obs_data)
                        )

                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds[networkspeci][
                        data_label
                    ] = np.arange(
                        len(station_data_availability_number), dtype=np.int32
                    )[
                        station_data_availability_number > 1
                    ]

                    # get colocated obs data array if have temporal colocation is active
                    obs_data[
                        self.read_instance.temporal_colocation_nans[networkspeci]
                    ] = np.nan

                    # get absolute data availability number per station in observational data array
                    if obs_data.size == 0:
                        station_data_availability_number = np.array([])
                    else:
                        station_data_availability_number = (
                            Stats.calculate_data_avail_number(obs_data)
                        )

                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds_temporal_colocation[
                        networkspeci
                    ][data_label] = np.arange(
                        len(station_data_availability_number), dtype=np.int32
                    )[
                        station_data_availability_number > 1
                    ]

            # get equivalent valid station indices for model arrays
            # subset observational valid station indices, with model array stations with > 1 valid measurements
            # therefore number of observational valid indices will always be > model valid indices
            for data_label in self.read_instance.data_labels:
                # check if data array is not an observational data array
                if data_label != self.read_instance.observations_data_label:
                    # get indices of valid observational data array stations
                    valid_station_inds = self.read_instance.valid_station_inds[
                        networkspeci
                    ][self.read_instance.observations_data_label].copy()

                    # get model data array (first subset by valid observational stations)
                    mod_data = self.read_instance.data_in_memory_filtered[networkspeci][
                        self.read_instance.data_labels.index(data_label),
                        valid_station_inds,
                        :,
                    ].copy()

                    # get absolute data availability number per station in model data array
                    if mod_data.size == 0:
                        station_data_availability_number = np.array([])
                    else:
                        station_data_availability_number = (
                            Stats.calculate_data_avail_number(mod_data)
                        )

                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds[networkspeci][
                        data_label
                    ] = valid_station_inds[
                        np.arange(
                            len(station_data_availability_number), dtype=np.int32
                        )[station_data_availability_number > 1]
                    ]

                    # get colocated model data array
                    mod_data = self.read_instance.data_in_memory_filtered[networkspeci][
                        self.read_instance.data_labels.index(data_label), :, :
                    ].copy()
                    mod_data[
                        self.read_instance.temporal_colocation_nans[networkspeci]
                    ] = np.nan
                    mod_data = mod_data[valid_station_inds, :]

                    # get absolute data availability number per station in model data array
                    if mod_data.size == 0:
                        station_data_availability_number = np.array([])
                    else:
                        station_data_availability_number = (
                            Stats.calculate_data_avail_number(mod_data)
                        )

                    # get indices of stations with > 1 available measurements
                    self.read_instance.valid_station_inds_temporal_colocation[
                        networkspeci
                    ][data_label] = valid_station_inds[
                        np.arange(
                            len(station_data_availability_number), dtype=np.int32
                        )[station_data_availability_number > 1]
                    ]

    def cap_stations(self):
        """Cap stations to a maximum number of stations"""

        if self.read_instance.station_cap:

            for networkspeci in self.read_instance.networkspecies:

                # cap metadata_in_memory
                #self.read_instance.metadata_in_memory_filtered[networkspeci] = self.read_instance.metadata_in_memory_filtered[networkspeci][:self.read_instance.station_cap, :]

                for data_label in self.read_instance.data_labels:

                    # cap valid station inds (colocated and not)
                    self.read_instance.valid_station_inds[networkspeci][data_label] = (
                        self.read_instance.valid_station_inds[networkspeci][
                            data_label
                        ][:self.read_instance.station_cap]
                    )
                    self.read_instance.valid_station_inds_temporal_colocation[networkspeci][data_label] = (
                        self.read_instance.valid_station_inds_temporal_colocation[networkspeci][
                            data_label
                        ][:self.read_instance.station_cap]
                    )

            # update available metadata fields for capped metadata
            update_metadata_fields(self.read_instance, cap=True)
