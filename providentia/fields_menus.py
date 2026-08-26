""" Functions to initialise and update fields for filtering (flags, period, metadata, etc.) """

import copy

import numpy as np
import pandas as pd


def init_flags(instance):
    """
    Initialise the internal dictionary structure and default values for network quality assurance flags.

    Parameters
    ----------
    instance : object
        An object instance to be updated with flag menu configurations.
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, "flag_menu"):
        instance.flag_menu = {
            "window_title": "Network QA Flags",
            "page_title": "Select standardised data reporter provided flags to filter by",
            "checkboxes": {},
        }
        instance.flag_menu["select_buttons"] = ["all", "clear", "default"]

    # reset fields
    instance.flag_menu["checkboxes"]["labels"] = np.array(
        sorted(
            instance.standard_data_flag_name_to_data_flag_code,
            key=instance.standard_data_flag_name_to_data_flag_code.get,
        )
    )
    instance.flag_menu["checkboxes"]["remove_default"] = np.array([], dtype=np.uint8)
    instance.flag_menu["checkboxes"]["remove_selected"] = np.array([], dtype=np.uint8)
    instance.flag_menu["checkboxes"]["map_vars"] = np.sort(
        list(instance.standard_data_flag_name_to_data_flag_code.values())
    )


def init_qa(instance):
    """
    Initialise the internal dictionary structure and default values for GHOST quality assurance flags.

    Parameters
    ----------
    instance : object
        An object instance to be updated with quality assurance menu configurations.
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, "qa_menu"):
        instance.qa_menu = {
            "window_title": "GHOST QA Flags",
            "page_title": "Select standardised quality assurance flags to filter by",
            "checkboxes": {},
        }
        instance.qa_menu["select_buttons"] = ["all", "clear", "default"]

    # reset fields
    if instance.network == ["actris/actris"]:
        qa_labels = {
            "Invalid Data Provider Flags - GHOST Decreed": 6,
            "Invalid Data Provider Flags - Network Decreed": 7,
        }
    else:
        qa_labels = instance.standard_QA_name_to_QA_code

    qa_key = qa_labels.get
    qa_map_vars = qa_labels.values()
    instance.qa_menu["checkboxes"]["labels"] = np.array(sorted(qa_labels, key=qa_key))
    instance.qa_menu["checkboxes"]["remove_default"] = np.array([], dtype=np.uint8)
    instance.qa_menu["checkboxes"]["remove_selected"] = np.array([], dtype=np.uint8)
    instance.qa_menu["checkboxes"]["map_vars"] = np.sort(list(qa_map_vars))


def init_models(instance):
    """
    Initialise the internal dictionary structure and default values for model selection and forecast options.

    Parameters
    ----------
    instance : object
        An object instance to be updated with model menu configurations.
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, "models_menu"):
        instance.models_menu = {
            "window_title": "Models",
            "page_title": "Select Models",
            "models": {},
        }
        instance.models_menu["select_buttons"] = ["all", "clear"]

    # reset fields
    instance.models_menu["models"]["labels"] = []
    instance.models_menu["models"]["keep_selected"] = {
        "interpolated": [],
        "noninterpolated": [],
    }
    instance.models_menu["models"]["enabled"] = {
        "interpolated": {},
        "noninterpolated": {},
    }
    instance.models_menu["models"]["forecast"] = {}
    instance.models_menu["models"]["forecast_days"] = {}
    instance.models_menu["models"]["map_vars"] = []


def init_multispecies(instance):
    """
    Initialise the internal dictionary structure for multispecies filtering fields and range boundaries.

    Parameters
    ----------
    instance : object
        An object instance to be updated with multispecies menu configurations.
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, "multispecies_menu"):
        instance.multispecies_menu = {
            "window_title": "Species Filtering",
            "page_title": "Select species to filter by",
            "multispecies": {},
        }

    # reset rangeboxes
    instance.multispecies_menu["multispecies"]["labels"] = []
    instance.multispecies_menu["multispecies"]["current_lower"] = {}
    instance.multispecies_menu["multispecies"]["current_upper"] = {}
    instance.multispecies_menu["multispecies"]["current_filter_species_fill_value"] = {}
    instance.multispecies_menu["multispecies"]["apply_selected"] = {}
    instance.multispecies_menu["multispecies"]["previous_lower"] = {}
    instance.multispecies_menu["multispecies"]["previous_upper"] = {}
    instance.multispecies_menu["multispecies"]["previous_apply"] = {}
    instance.multispecies_menu["multispecies"][
        "previous_filter_species_fill_value"
    ] = {}


def init_coverage(instance):
    """
    Initialise the internal dictionary structure and default values for data coverage percentage fields.

    Parameters
    ----------
    instance : object
        An object instance to be updated with coverage menu configurations.
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, "coverage_menu"):
        instance.coverage_menu = {
            "window_title": "% Data Coverage",
            "page_title": "Select Minimum Required % Data Coverage",
            "rangeboxes": {},
        }

    # reset fields
    instance.coverage_menu["rangeboxes"]["tooltips"] = []
    instance.coverage_menu["rangeboxes"]["labels"] = []
    instance.coverage_menu["rangeboxes"]["current_lower"] = []
    instance.coverage_menu["rangeboxes"]["map_vars_old"] = []
    instance.coverage_menu["rangeboxes"]["map_vars"] = []
    instance.coverage_menu["rangeboxes"]["subtitles"] = []
    instance.coverage_menu["rangeboxes"]["subtitle_inds"] = []


def init_period(instance):
    """
    Initialise the internal dictionary structure and default values for data period selection fields.

    Parameters
    ----------
    instance : object
        An object instance to be updated with period menu configurations.
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, "period_menu"):
        instance.period_menu = {
            "window_title": "Data Period",
            "page_title": "Select Data Periods",
            "checkboxes": {},
        }

    # reset fields
    instance.period_menu["checkboxes"]["labels"] = []
    instance.period_menu["checkboxes"]["keep_selected"] = []
    instance.period_menu["checkboxes"]["remove_selected"] = []


def init_metadata(instance):
    """
    Initialise the internal dictionary structure and nested menus for station metadata.

    Parameters
    ----------
    instance : object
        An object instance to be updated with metadata menu configurations.
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, "metadata_menu"):
        instance.metadata_types = {
            "STATION POSITION": "Filter stations by measurement position",
            "STATION CLASSIFICATIONS": "Filter stations by station provided classifications",
            "STATION MISCELLANEOUS": "Filter stations by miscellaneous station provided metadata",
            "GLOBALLY GRIDDED CLASSIFICATIONS": "Filter stations by globally gridded classifications",
            "MEASUREMENT PROCESS INFORMATION": "Filter stations by measurement process information",
        }

        instance.metadata_menu = {
            "window_title": "Metadata",
            "page_title": "Select metadata type to filter stations by",
            "navigation_buttons": {},
        }

        instance.metadata_menu["navigation_buttons"]["labels"] = list(
            instance.metadata_types.keys()
        )
        instance.metadata_menu["navigation_buttons"]["tooltips"] = [
            instance.metadata_types[key]
            for key in instance.metadata_menu["navigation_buttons"]["labels"]
        ]

        for metadata_type_ii, metadata_type in enumerate(
            instance.metadata_menu["navigation_buttons"]["labels"]
        ):
            # setup nested menu
            instance.metadata_menu[metadata_type] = {
                "window_title": metadata_type,
                "page_title": instance.metadata_menu["navigation_buttons"]["tooltips"][
                    metadata_type_ii
                ],
                "navigation_buttons": {},
                "rangeboxes": {},
            }

    # reset fields
    for metadata_type_ii, metadata_type in enumerate(
        instance.metadata_menu["navigation_buttons"]["labels"]
    ):
        # reset rangebox labels
        instance.metadata_menu[metadata_type]["rangeboxes"]["labels"] = [
            metadata_var
            for metadata_var in instance.metadata_vars_to_read
            if (
                instance.standard_metadata[metadata_var]["metadata_type"]
                == metadata_type
            )
            & (instance.standard_metadata[metadata_var]["data_type"] is not object)
        ]

        # reset rangebox tooltips
        instance.metadata_menu[metadata_type]["rangeboxes"]["tooltips"] = [
            instance.standard_metadata[metadata_var]["description"]
            for metadata_var in instance.metadata_menu[metadata_type]["rangeboxes"][
                "labels"
            ]
        ]

        # reset rangeboxes
        instance.metadata_menu[metadata_type]["rangeboxes"]["current_lower"] = [
            "nan"
        ] * len(instance.metadata_menu[metadata_type]["rangeboxes"]["labels"])
        instance.metadata_menu[metadata_type]["rangeboxes"]["current_upper"] = [
            "nan"
        ] * len(instance.metadata_menu[metadata_type]["rangeboxes"]["labels"])
        instance.metadata_menu[metadata_type]["rangeboxes"]["lower_default"] = [
            "nan"
        ] * len(instance.metadata_menu[metadata_type]["rangeboxes"]["labels"])
        instance.metadata_menu[metadata_type]["rangeboxes"]["upper_default"] = [
            "nan"
        ] * len(instance.metadata_menu[metadata_type]["rangeboxes"]["labels"])
        instance.metadata_menu[metadata_type]["rangeboxes"]["apply_selected"] = []

        # reset checkbox labels
        instance.metadata_menu[metadata_type]["navigation_buttons"]["labels"] = [
            metadata_var
            for metadata_var in instance.metadata_vars_to_read
            if (
                instance.standard_metadata[metadata_var]["metadata_type"]
                == metadata_type
            )
            & (instance.standard_metadata[metadata_var]["data_type"] is object)
        ]

        # reset checkbox tooltips
        instance.metadata_menu[metadata_type]["navigation_buttons"]["tooltips"] = [
            instance.standard_metadata[metadata_var]["description"]
            for metadata_var in instance.metadata_menu[metadata_type][
                "navigation_buttons"
            ]["labels"]
        ]

        # reset checkboxes
        for metadata_var in instance.metadata_menu[metadata_type]["navigation_buttons"][
            "labels"
        ]:
            # metadata variable already in dict?
            # then just reset lists
            if metadata_var in instance.metadata_menu[metadata_type]:
                instance.metadata_menu[metadata_type][metadata_var]["checkboxes"][
                    "labels"
                ] = []
                instance.metadata_menu[metadata_type][metadata_var]["checkboxes"][
                    "keep_selected"
                ] = []
                instance.metadata_menu[metadata_type][metadata_var]["checkboxes"][
                    "keep_default"
                ] = []
                instance.metadata_menu[metadata_type][metadata_var]["checkboxes"][
                    "remove_selected"
                ] = []
                instance.metadata_menu[metadata_type][metadata_var]["checkboxes"][
                    "remove_default"
                ] = []
            # otherwise, create infrastructure to store metadata var information
            else:
                instance.metadata_menu[metadata_type][metadata_var] = {
                    "window_title": metadata_var,
                    "page_title": "Filter stations by unique {} metadata".format(
                        metadata_var
                    ),
                    "checkboxes": {},
                }
                instance.metadata_menu[metadata_type][metadata_var]["checkboxes"] = {
                    "labels": [],
                    "keep_selected": [],
                    "keep_default": [],
                    "remove_selected": [],
                    "remove_default": [],
                }

        # remove metadata type checkbox var dicts not in metadata_vars_to_read
        metadata_type_checkbox_vars = [
            metadata_type_var
            for metadata_type_var in instance.metadata_menu[metadata_type].keys()
            if metadata_type_var
            not in ["window_title", "page_title", "navigation_buttons", "rangeboxes"]
        ]
        metadata_type_checkbox_vars_to_remove = [
            metadata_type_checkbox_var
            for metadata_type_checkbox_var in metadata_type_checkbox_vars
            if metadata_type_checkbox_var not in instance.metadata_vars_to_read
        ]
        for (
            metadata_type_checkbox_var_to_remove
        ) in metadata_type_checkbox_vars_to_remove:
            del instance.metadata_menu[metadata_type][
                metadata_type_checkbox_var_to_remove
            ]


def update_qa(instance):
    """
    Update the quality assurance menu labels and mappings.

    Parameters
    ----------
    instance : object
        An object instance to be updated with quality assurance menu configurations.
    """

    # reset fields
    if instance.selected_network == "actris/actris":
        qa_labels = {
            "Invalid Data Provider Flags - GHOST Decreed": 6,
            "Invalid Data Provider Flags - Network Decreed": 7,
        }
    else:
        qa_labels = instance.standard_QA_name_to_QA_code

    qa_key = qa_labels.get
    qa_map_vars = qa_labels.values()
    instance.qa_menu["checkboxes"]["labels"] = np.array(sorted(qa_labels, key=qa_key))
    instance.qa_menu["checkboxes"]["remove_default"] = np.array([], dtype=np.uint8)
    instance.qa_menu["checkboxes"]["remove_selected"] = np.array([], dtype=np.uint8)
    instance.qa_menu["checkboxes"]["map_vars"] = np.sort(list(qa_map_vars))


def update_coverage_fields(instance):
    """
    Update the data coverage menu labels and values based on the active data resolution and network type.

    Parameters
    ----------
    instance : object
        An object instance to be updated with coverage menu configurations.
    """

    # get previously set rangebox labels / values
    previous_mapped_labels = copy.deepcopy(
        instance.coverage_menu["rangeboxes"]["map_vars"]
    )
    previous_mapped_labels_old = copy.deepcopy(
        instance.coverage_menu["rangeboxes"]["map_vars_old"]
    )
    previous_lower = copy.deepcopy(
        instance.coverage_menu["rangeboxes"]["current_lower"]
    )

    # build coverage menu
    if (instance.reading_ghost) & (instance.ghost_features == "max"):
        network_type = "ghost"
    else:
        network_type = "nonghost"

    # set resolution to get coverage information based on active resolution
    if instance.resolution in ["hourly", "hourly_instantaneous"]:
        resolution = "hourly"
    elif instance.resolution in ["3hourly", "3hourly_instantaneous"]:
        resolution = "3hourly"
    elif instance.resolution in ["6hourly", "6hourly_instantaneous"]:
        resolution = "6hourly"
    elif instance.resolution == "daily":
        resolution = "daily"
    elif instance.resolution == "monthly":
        resolution = "monthly"
    elif instance.resolution == "annual":
        resolution = "annual"

    instance.coverage_menu["rangeboxes"]["map_vars"] = instance.coverage_info[
        network_type
    ][resolution]["map_vars"]
    instance.coverage_menu["rangeboxes"]["map_vars_old"] = instance.coverage_info[
        network_type
    ][resolution]["map_vars_old"]
    instance.coverage_menu["rangeboxes"]["labels"] = instance.coverage_info[
        network_type
    ][resolution]["labels"]
    instance.coverage_menu["rangeboxes"]["subtitles"] = instance.coverage_info[
        network_type
    ][resolution]["subtitles"]
    instance.coverage_menu["rangeboxes"]["subtitle_inds"] = instance.coverage_info[
        network_type
    ][resolution]["subtitle_inds"]

    # initialise rangebox values --> for data coverage fields
    # the default is 0%, for max gap fields % the default is 100%
    instance.coverage_menu["rangeboxes"]["current_lower"] = []
    for label_ii, (label_mapped, label_mapped_old) in enumerate(
        zip(
            instance.coverage_menu["rangeboxes"]["map_vars"],
            instance.coverage_menu["rangeboxes"]["map_vars_old"],
        )
    ):
        if "max_gap" in label_mapped:
            instance.coverage_menu["rangeboxes"]["current_lower"].append("100")
        else:
            instance.coverage_menu["rangeboxes"]["current_lower"].append("0")

        # label previously existed?
        if (label_mapped in previous_mapped_labels) or (
            label_mapped_old in previous_mapped_labels_old
        ):
            instance.coverage_menu["rangeboxes"]["current_lower"][
                label_ii
            ] = previous_lower[previous_mapped_labels.index(label_mapped)]


def update_period_fields(instance):
    """
    Update the data period menu labels based on the active temporal resolution.

    Parameters
    ----------
    instance : object
        An object instance to be updated with period menu configurations.
    """

    # hourly/3hourly/6hourly temporal resolution?
    if "hourly" in instance.resolution:
        instance.period_menu["checkboxes"]["labels"] = [
            "Daytime",
            "Nighttime",
            "Weekday",
            "Weekend",
            "Spring",
            "Summer",
            "Autumn",
            "Winter",
        ]
    # daily temporal resolution?
    elif instance.resolution == "daily":
        instance.period_menu["checkboxes"]["labels"] = [
            "Weekday",
            "Weekend",
            "Spring",
            "Summer",
            "Autumn",
            "Winter",
        ]
        # drop selected fields from higher temporal resolutions
        labels_to_remove = ["Daytime", "Nighttime"]
        for label in labels_to_remove:
            if label in instance.period_menu["checkboxes"]["keep_selected"]:
                instance.period_menu["checkboxes"]["keep_selected"].remove(label)
            if label in instance.period_menu["checkboxes"]["remove_selected"]:
                instance.period_menu["checkboxes"]["remove_selected"].remove(label)

    # monthly temporal resolution?
    elif instance.resolution == "monthly":
        instance.period_menu["checkboxes"]["labels"] = [
            "Spring",
            "Summer",
            "Autumn",
            "Winter",
        ]
        # drop selected fields from higher temporal resolutions
        labels_to_remove = ["Daytime", "Nighttime", "Weekday", "Weekend"]
        for label in labels_to_remove:
            if label in instance.period_menu["checkboxes"]["keep_selected"]:
                instance.period_menu["checkboxes"]["keep_selected"].remove(label)
            if label in instance.period_menu["checkboxes"]["remove_selected"]:
                instance.period_menu["checkboxes"]["remove_selected"].remove(label)

    # annual temporal resolution?
    elif instance.resolution == "annual":
        instance.period_menu["checkboxes"]["labels"] = []
        # drop selected fields from higher temporal resolutions
        labels_to_remove = [
            "Daytime",
            "Nighttime",
            "Weekday",
            "Weekend",
            "Spring",
            "Summer",
            "Autumn",
            "Winter",
        ]
        for label in labels_to_remove:
            if label in instance.period_menu["checkboxes"]["keep_selected"]:
                instance.period_menu["checkboxes"]["keep_selected"].remove(label)
            if label in instance.period_menu["checkboxes"]["remove_selected"]:
                instance.period_menu["checkboxes"]["remove_selected"].remove(label)


def update_metadata_fields(instance, cap=False):
    """
    Update metadata menu fields using metadata values currently loaded in memory.

    Non-numeric metadata variables are updated with the sorted unique categorical
    fields found across all active network/species combinations. Numeric metadata
    variables are updated with the minimum and maximum valid bounds found across
    all active network/species combinations.

    Existing user selections are preserved where possible. For categorical
    metadata, previously selected keep/remove checkbox selections are copied
    across if the corresponding fields still exist after the metadata update.
    For numeric metadata, previous lower/upper bounds are retained if they are
    stricter than the previous defaults, matching the existing behaviour.

    Parameters
    ----------
    instance : object
        Object instance containing metadata arrays, metadata menu definitions,
        active network/species configuration, station-validity information, and
        standard metadata definitions. The function updates
        ``instance.metadata_menu`` in place.
    cap: boolean, optional
        optional boolean indicating that have capped stations. 

    Returns
    -------
    None
        The function modifies ``instance.metadata_menu`` in place and does not
        return a value.
    """

    metadata_menu = instance.metadata_menu
    metadata_vars = instance.metadata_vars_to_read
    standard_metadata = instance.standard_metadata
    networkspecies = instance.networkspecies
    if cap:
        metadata_in_memory = instance.metadata_in_memory_filtered
    else:
        metadata_in_memory = instance.metadata_in_memory

    # ------------------------------------------------------------------
    # Local helper: clean object metadata values.
    # ------------------------------------------------------------------
    def _clean_object_unique_values(values):
        """
        Remove invalid categorical metadata values from an already reduced array.

        Values equal to the string ``'nan'``, empty strings, and pandas-style
        missing values are removed. This preserves the behaviour of the current
        implementation for object metadata, but is intended to operate on a
        reduced unique-value array rather than the full station/month metadata
        array.

        Parameters
        ----------
        values : numpy.ndarray
            One-dimensional array of object/categorical metadata values. This
            array is usually already reduced by ``pd.unique``.

        Returns
        -------
        clean_values : numpy.ndarray
            One-dimensional array containing only valid categorical metadata
            values.
        """
        values = np.asarray(values, dtype=object).ravel()

        if values.size == 0:
            return values

        return values[(values != "nan") & (values != "") & (~pd.isna(values))]

    # ------------------------------------------------------------------
    # Local helper: get values for one metadata variable/networkspecies.
    # ------------------------------------------------------------------
    def _get_metadata_values(arr, inds=None):
        """
        Return flattened metadata values, optionally restricted to valid stations.

        Parameters
        ----------
        arr : numpy.ndarray
            Metadata array for one network/species and one metadata variable.
        inds : numpy.ndarray or None, optional
            Optional valid-station indices. If provided, only those stations are
            selected before flattening.

        Returns
        -------
        values : numpy.ndarray
            One-dimensional metadata value array.
        """
        if inds is None:
            return arr.ravel()

        return arr[inds].ravel()

    # ------------------------------------------------------------------
    # 1) Check whether metadata menu needs reset.
    # ------------------------------------------------------------------
    reset_meta = False
    nav_labels = metadata_menu["navigation_buttons"]["labels"]

    for metadata_type in nav_labels:
        required_count = 0

        for meta_var in metadata_vars:
            meta_info = standard_metadata[meta_var]

            if (
                meta_info["metadata_type"] is metadata_type
                and meta_info["data_type"] is not object
            ):
                required_count += 1

        if required_count != len(metadata_menu[metadata_type]["rangeboxes"]["labels"]):
            reset_meta = True
            break

    if reset_meta:
        init_metadata(instance)
        metadata_menu = instance.metadata_menu

    # ------------------------------------------------------------------
    # 2) Normalize network/species to lists once.
    # ------------------------------------------------------------------
    if isinstance(instance.network, str):
        instance.network = [instance.network]

    if isinstance(instance.species, str):
        instance.species = [instance.species]

    # ------------------------------------------------------------------
    # 3) Precompute valid station indices once.
    # ------------------------------------------------------------------
    use_valid_station_inds = (
        hasattr(instance, "valid_station_inds") and not instance.reading_ghost
    )

    valid_indices_by_netspec = {}

    if use_valid_station_inds:
        obs_label = instance.observations_data_label

        if instance.temporal_colocation:
            valid_source = instance.valid_station_inds_temporal_colocation
        else:
            valid_source = instance.valid_station_inds

        for netspec in networkspecies:
            valid_indices_by_netspec[netspec] = valid_source[netspec][obs_label]

    # ------------------------------------------------------------------
    # 4) Precompute numeric rangebox label indices.
    # ------------------------------------------------------------------
    rangebox_index_by_type = {}

    for metadata_type in nav_labels:
        if "rangeboxes" in metadata_menu[metadata_type]:
            labels = metadata_menu[metadata_type]["rangeboxes"]["labels"]
            rangebox_index_by_type[metadata_type] = {
                label: ii for ii, label in enumerate(labels)
            }

    # ------------------------------------------------------------------
    # 5) Main metadata update loop.
    # ------------------------------------------------------------------
    for meta_var in metadata_vars:
        meta_info = standard_metadata[meta_var]
        metadata_type = meta_info["metadata_type"]
        metadata_data_type = meta_info["data_type"]

        # ==============================================================
        # OBJECT / CATEGORICAL METADATA
        # ==============================================================
        if metadata_data_type is object:
            # ----------------------------------------------------------
            # Fast categorical path:
            #
            # Use pd.unique on each flattened chunk first. For repeated
            # station/month metadata this can reduce thousands of values
            # down to a handful before missing-value checks and final sort.
            # ----------------------------------------------------------
            unique_chunks = []

            for netspec in networkspecies:
                arr = metadata_in_memory[netspec][meta_var]

                if use_valid_station_inds:
                    vals = _get_metadata_values(arr, valid_indices_by_netspec[netspec])
                else:
                    vals = _get_metadata_values(arr)

                if vals.size == 0:
                    continue

                # Hash-based unique is usually faster than np.unique on
                # large repeated object arrays.
                vals_unique = pd.unique(vals)
                vals_unique = _clean_object_unique_values(vals_unique)

                if vals_unique.size > 0:
                    unique_chunks.append(vals_unique)

            if len(unique_chunks) == 0:
                new_fields = np.array([], dtype=object)
            elif len(unique_chunks) == 1:
                # Final np.unique preserves the sorted output behaviour of
                # the current implementation.
                new_fields = np.unique(unique_chunks[0])
            else:
                new_fields = np.unique(np.concatenate(unique_chunks))

            checkbox_menu = metadata_menu[metadata_type][meta_var]["checkboxes"]

            previous_fields = checkbox_menu["labels"]
            previous_keep = checkbox_menu["keep_selected"]
            previous_remove = checkbox_menu["remove_selected"]

            previous_field_set = set(previous_fields)
            previous_keep_set = set(previous_keep)
            previous_remove_set = set(previous_remove)

            checkbox_menu["labels"] = new_fields

            keep_selected = []
            remove_selected = []

            for field in new_fields:
                if field in previous_field_set:
                    if field in previous_keep_set:
                        keep_selected.append(field)
                    if field in previous_remove_set:
                        remove_selected.append(field)

            checkbox_menu["keep_selected"] = keep_selected
            checkbox_menu["remove_selected"] = remove_selected
            checkbox_menu["keep_default"] = []
            checkbox_menu["remove_default"] = []

            continue

        # ==============================================================
        # NUMERIC METADATA
        # ==============================================================

        current_min = None
        current_max = None

        for netspec in networkspecies:
            arr = metadata_in_memory[netspec][meta_var]

            if use_valid_station_inds:
                vals = _get_metadata_values(arr, valid_indices_by_netspec[netspec])
            else:
                vals = _get_metadata_values(arr)

            if vals.size == 0:
                continue

            # Same validity rule as current function for numeric metadata:
            # ignore NaNs only.
            #
            # Avoid np.nanmin / np.nanmax warnings for all-NaN chunks.
            nan_vals = np.isnan(vals)

            if nan_vals.all():
                continue

            chunk_min = np.nanmin(vals)
            chunk_max = np.nanmax(vals)

            if current_min is None:
                current_min = chunk_min
                current_max = chunk_max
            else:
                if chunk_min < current_min:
                    current_min = chunk_min
                if chunk_max > current_max:
                    current_max = chunk_max

        rangeboxes = metadata_menu[metadata_type]["rangeboxes"]
        meta_var_index = rangebox_index_by_type[metadata_type][meta_var]

        previous_lower_default = rangeboxes["lower_default"][meta_var_index]
        previous_upper_default = rangeboxes["upper_default"][meta_var_index]
        previous_lower = rangeboxes["current_lower"][meta_var_index]
        previous_upper = rangeboxes["current_upper"][meta_var_index]

        if current_min is not None:
            min_val = str(current_min)
            max_val = str(current_max)

            # Preserve current behaviour:
            # previous bounds are compared as strings, not numerically.
            new_lower = min_val

            if previous_lower not in ("nan", None) and previous_lower_default not in (
                "nan",
                None,
            ):
                if previous_lower > previous_lower_default:
                    new_lower = previous_lower

            new_upper = max_val

            if previous_upper not in ("nan", None) and previous_upper_default not in (
                "nan",
                None,
            ):
                if previous_upper < previous_upper_default:
                    new_upper = previous_upper

            rangeboxes["current_lower"][meta_var_index] = new_lower
            rangeboxes["current_upper"][meta_var_index] = new_upper
            rangeboxes["lower_default"][meta_var_index] = min_val
            rangeboxes["upper_default"][meta_var_index] = max_val

        else:
            rangeboxes["current_lower"][meta_var_index] = "nan"
            rangeboxes["lower_default"][meta_var_index] = "nan"
            rangeboxes["current_upper"][meta_var_index] = "nan"
            rangeboxes["upper_default"][meta_var_index] = "nan"

            if meta_var in rangeboxes["apply_selected"]:
                rangeboxes["apply_selected"].remove(meta_var)

    return None


def multispecies_conf(instance):
    """
    Configure multispecies filtering variables and widget states using settings loaded from a configuration file.

    Parameters
    ----------
    instance : object
        An object instance to be updated with multispecies menu configurations.
    """

    if hasattr(instance, "filter_species"):
        filter_species = copy.deepcopy(instance.filter_species)
        for (networkspeci_ii, networkspeci), networkspeci_bounds in zip(
            enumerate(filter_species.keys()), filter_species.values()
        ):
            for bounds in networkspeci_bounds:
                # update menu_current
                if (
                    "networkspeci_" + str(networkspeci_ii)
                ) not in instance.multispecies_menu["multispecies"]["labels"]:
                    instance.multispecies_menu["multispecies"]["labels"].append(
                        "networkspeci_" + str(networkspeci_ii)
                    )

                # add values
                instance.multispecies_menu["multispecies"]["current_lower"][
                    networkspeci_ii
                ] = bounds[0]
                instance.multispecies_menu["multispecies"]["current_upper"][
                    networkspeci_ii
                ] = bounds[1]
                instance.multispecies_menu["multispecies"][
                    "current_filter_species_fill_value"
                ][networkspeci_ii] = bounds[2]
                instance.multispecies_menu["multispecies"]["apply_selected"][
                    networkspeci_ii
                ] = True

                # set initial selected config variables as set .conf files or defaults
                instance.selected_widget_network.update(
                    {networkspeci_ii: networkspeci.split("|")[0]}
                )
                instance.selected_widget_matrix.update(
                    {
                        networkspeci_ii: instance.parameter_dictionary[
                            networkspeci.split("|")[1]
                        ]["matrix"]
                    }
                )
                instance.selected_widget_species.update(
                    {networkspeci_ii: networkspeci.split("|")[1]}
                )
                instance.selected_widget_lower.update({networkspeci_ii: bounds[0]})
                instance.selected_widget_upper.update({networkspeci_ii: bounds[1]})
                instance.selected_widget_filter_species_fill_value.update(
                    {networkspeci_ii: bounds[2]}
                )
                instance.selected_widget_apply.update({networkspeci_ii: True})

                networkspeci_ii += 1

            # filtering tab is initialized from conf
            instance.multispecies_initialisation = False


def coverage_conf(instance):
    """
    Configure coverage filter variables using settings loaded from a configuration file.

    Parameters
    ----------
    instance : object
        An object instance to be updated with coverage menu configurations.
    """

    for i, (label, label_old) in enumerate(
        zip(
            instance.coverage_menu["rangeboxes"]["map_vars"],
            instance.coverage_menu["rangeboxes"]["map_vars_old"],
        )
    ):
        if hasattr(instance, label):
            instance.coverage_menu["rangeboxes"]["current_lower"][i] = str(
                getattr(instance, label)
            )
        elif hasattr(instance, label_old):
            instance.coverage_menu["rangeboxes"]["current_lower"][i] = str(
                getattr(instance, label_old)
            )


def period_conf(instance):
    """
    Configure period filter variables and checkbox selections using settings loaded from a configuration file.

    Parameters
    ----------
    instance : object
        An object instance to be updated with period menu configurations.
    """

    from .configuration import split_options

    if hasattr(instance, "period"):
        keeps, removes = split_options(instance, instance.period)
        instance.period_menu["checkboxes"]["keep_selected"] = keeps
        instance.period_menu["checkboxes"]["remove_selected"] = removes


def metadata_conf(instance):
    """
    Configure metadata filter ranges and categorical selections using settings loaded from a configuration file.

    Parameters
    ----------
    instance : object
        An object instance to be updated with metadata menu configurations.
    """

    from .configuration import split_options

    for menu_type in instance.metadata_types:
        # treat first ranges
        current_lower = []
        current_upper = []
        apply_selected = []
        for i, label in enumerate(
            instance.metadata_menu[menu_type]["rangeboxes"]["labels"]
        ):
            if hasattr(instance, label):
                current_lower.append(str(getattr(instance, label)[0]))
                current_upper.append(str(getattr(instance, label)[1]))
                apply_selected.append(label)
            else:
                current_lower.append(
                    instance.metadata_menu[menu_type]["rangeboxes"]["current_lower"][i]
                )
                current_upper.append(
                    instance.metadata_menu[menu_type]["rangeboxes"]["current_upper"][i]
                )
        instance.metadata_menu[menu_type]["rangeboxes"]["current_lower"] = current_lower
        instance.metadata_menu[menu_type]["rangeboxes"]["current_upper"] = current_upper
        instance.metadata_menu[menu_type]["rangeboxes"][
            "apply_selected"
        ] = apply_selected

        # and then treat the keep/remove
        for label in instance.metadata_menu[menu_type]["navigation_buttons"]["labels"]:
            if hasattr(instance, label):
                keeps, removes = split_options(instance, getattr(instance, label))
                instance.metadata_menu[menu_type][label]["checkboxes"][
                    "keep_selected"
                ] = keeps
                instance.metadata_menu[menu_type][label]["checkboxes"][
                    "remove_selected"
                ] = removes
