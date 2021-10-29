"""
Contains multiple static functions that are shared by
the dashboard and the offline version. Currently, it
includes function related to the initialization and update
of metadata, checks fields coming from conf files etc.
"""

import copy
import json
import datetime
import numpy as np
import pandas as pd

from .config import split_options


def which_bounds(instance, species):
    """Returns lower/upper bounds of species selected. If
    there are bounds defined in a config file, fill that value,
    if they are within the feasible bounds of the species.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    :param species: The species current selected for evaluation (e.g. sconco3)
    :type species: str
    :return: The lower and upper bound
    :rtype: np.float32
    """

    lower = np.float32(instance.parameter_dictionary[species]['extreme_lower_limit'])
    upper = np.float32(instance.parameter_dictionary[species]['extreme_upper_limit'])

    if hasattr(instance, 'lower_bound'):
        if instance.lower_bound >= lower:
            lower = instance.lower_bound

    if hasattr(instance, 'upper_bound'):
        if instance.upper_bound <= upper:
            upper = instance.upper_bound

    return np.float32(lower), np.float32(upper)


def which_qa(instance, return_defaults=False):
    """Returns QA flags for the species selected. If return_defaults
    is set to true, it will just return the default values according
    to GHOST standards. If there is a config file which has QA defined,
    it will return those.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    :param return_defaults: flag for just returning the default QA
    :type return_defaults: bool
    :return: QA flags' codes in list
    :rtype: list
    """

    if return_defaults or (not hasattr(instance, 'qa')):
        if instance.selected_species in instance.qa_exceptions:
            return instance.specific_qa
        else:
            return instance.general_qa

    if hasattr(instance, 'qa'):
        # if conf has only 1 QA
        if isinstance(instance.qa, int):
            return [instance.qa]
        # if the QAs are written with their names
        elif instance.qa == "":
            return []
        elif isinstance(instance.qa, str):
            return [instance.standard_QA_name_to_QA_code[q.strip()] for q in instance.qa.split(",")]
        # return subset the user has selected in conf
        else:
            return list(instance.qa)


def which_flags(instance):
    """If there are flags coming from a config file,
    select those. Otherwise, return empty list.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    :return: list of flags' codes
    :rtype: list
    """

    if hasattr(instance, 'flags'):
        # if conf has only one flag
        if isinstance(instance.flags, int):
            return [instance.flags]
        # if flags are writtern as strings
        elif instance.flags == "":
            return []
        elif isinstance(instance.flags, str):
            return [instance.standard_data_flag_name_to_data_flag_code[f.strip()]
                    for f in instance.flags.split(",")]
        else:
            return list(instance.flags)
    else:
        return []


def get_qa_codes(instance):
    """Retrieve QA codes from GHOST_standards using the QA flags' names.

    Specific flags are defined for the following species:
    ['WND_DIR_10M','WND_SPD_10M','RH_2M','PREC_ACCUM','SNOW_ACCUM',
    'SNOW_DEPTH','CEILING_HEIGHT','VIS_DIST','CLOUD_CVG','CLOUD_CVG_FRAC']

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    :return: three lists which contain flags' codes
    :rtype: list
    """

    # get names from json files
    specific_qa_names = json.load(open(
        "providentia/conf/default_flags.json"))['specific_qa']
    general_qa_names = json.load(open(
        "providentia/conf/default_flags.json"))['general_qa']
    # get codes
    specific_qa = [instance.standard_QA_name_to_QA_code[qa_name]
                   for qa_name in specific_qa_names]
    general_qa = [instance.standard_QA_name_to_QA_code[qa_name]
                  for qa_name in general_qa_names]
    # get difference of flags, needed later for updating default selection
    qa_diff = list(set(general_qa) - set(specific_qa))

    return specific_qa, general_qa, qa_diff


def exceedance_lim(species):
    """Returns the exceedance limit depending on the species input. If
    species doesn't have a reported limit, returns np.NaN.

    :param species: name of species currently selected (e.g. sconco3)
    :type species: str
    :return: value of exceedance limit
    :rtype: int
    """
    exceedance_limits = {'sconco3': 90.21, 'sconcno2': 106.38}
    if species in exceedance_limits.keys():
        return exceedance_limits[species]
    else:
        return np.NaN


def temp_axis_dict():
    """Returns a temporal mapping as a dictionary used for the plots.

    :return: numbering of months/days
    :rtype: dict
    """
    map_dict = {'dayofweek': {0: 'M', 1: 'T', 2: 'W', 3: 'T', 4: 'F', 5: 'S', 6: 'S'},
                'month': {1: 'J', 2: 'F', 3: 'M', 4: 'A', 5: 'M', 6: 'J',
                          7: 'J', 8: 'A', 9: 'S', 10: 'O', 11: 'N', 12: 'D'}
                }
    return map_dict


def representativity_conf(instance):
    """Comes here if there is a configuration loaded. Checks if there is a
    representative field loaded in the object from the conf and if there is
    assigns the value in the representativity menu.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """
    for i, label in enumerate(instance.representativity_menu['rangeboxes']['labels']):
        if hasattr(instance, label):
            instance.representativity_menu['rangeboxes']['current_lower'][i] = str(getattr(instance, label))


def init_metadata(instance):
    """Initialise internal structure to store metadata.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    # setup pop-up window menu tree for metadata
    metadata_types = {'STATION POSITION': 'Filter stations by measurement position',
                      'STATION CLASSIFICATIONS': 'Filter stations by station provided classifications',
                      'STATION MISCELLANEOUS': 'Filter stations by miscellaneous station provided metadata',
                      'GLOBALLY GRIDDED CLASSIFICATIONS': 'Filter stations by globally gridded classifications',
                      'MEASUREMENT PROCESS INFORMATION': 'Filter stations by measurement process information'}
    metadata_menu = {'window_title': 'METADATA', 'page_title': 'Select metadata type to filter stations by',
                     'navigation_buttons': {}}

    metadata_menu['navigation_buttons']['labels'] = list(metadata_types.keys())
    metadata_menu['navigation_buttons']['tooltips'] = [metadata_types[key] for key in
                                                       metadata_menu['navigation_buttons']['labels']]

    for metadata_type_ii, metadata_type in enumerate(metadata_menu['navigation_buttons']['labels']):
        metadata_menu[metadata_type] = {'window_title': metadata_type,
                                        'page_title': metadata_menu['navigation_buttons']['tooltips'][
                                            metadata_type_ii], 'navigation_buttons': {}, 'rangeboxes': {}}
        metadata_menu[metadata_type]['navigation_buttons']['labels'] = \
            [metadata_name for metadata_name in instance.standard_metadata.keys() if
             (instance.standard_metadata[metadata_name]['metadata_type'] == metadata_type) &
             (instance.standard_metadata[metadata_name]['data_type'] == np.object)]
        metadata_menu[metadata_type]['navigation_buttons']['tooltips'] = \
            [instance.standard_metadata[metadata_name]['description'] for metadata_name in
             metadata_menu[metadata_type]['navigation_buttons']['labels']]

        for label in metadata_menu[metadata_type]['navigation_buttons']['labels']:
            metadata_menu[metadata_type][label] = {'window_title': label,
                                                   'page_title': 'Filter stations by unique {} metadata'.format(
                                                       label), 'checkboxes': {}}
            metadata_menu[metadata_type][label]['checkboxes'] = {'labels': [], 'keep_selected': [],
                                                                 'remove_selected': [], 'keep_default': [],
                                                                 'remove_default': []}

        metadata_menu[metadata_type]['rangeboxes']['labels'] = \
            [metadata_name for metadata_name in instance.standard_metadata.keys()
             if (instance.standard_metadata[metadata_name]['metadata_type'] == metadata_type)
             & (instance.standard_metadata[metadata_name]['data_type'] != np.object)]
        metadata_menu[metadata_type]['rangeboxes']['tooltips'] = \
            [instance.standard_metadata[metadata_name]['description']
             for metadata_name in metadata_menu[metadata_type]['rangeboxes']['labels']]
        metadata_menu[metadata_type]['rangeboxes']['current_lower'] = \
            ['nan'] * len(metadata_menu[metadata_type]['rangeboxes']['labels'])
        metadata_menu[metadata_type]['rangeboxes']['current_upper'] = \
            ['nan'] * len(metadata_menu[metadata_type]['rangeboxes']['labels'])
        metadata_menu[metadata_type]['rangeboxes']['lower_default'] = \
            ['nan'] * len(metadata_menu[metadata_type]['rangeboxes']['labels'])
        metadata_menu[metadata_type]['rangeboxes']['upper_default'] = \
            ['nan'] * len(metadata_menu[metadata_type]['rangeboxes']['labels'])

    return metadata_types, metadata_menu


def meta_from_conf(instance):
    """Comes here if there in a loaded configuration
    there are also metadata fields.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    for menu_type in instance.metadata_types:
        # treat first ranges
        for i, label_cap in enumerate(instance.metadata_menu[menu_type]['rangeboxes']['labels']):
            label = label_cap.lower()
            if hasattr(instance, label):
                instance.metadata_menu[menu_type]['rangeboxes']['current_lower'][i] = str(getattr(instance, label)[0])
                instance.metadata_menu[menu_type]['rangeboxes']['current_upper'][i] = str(getattr(instance, label)[1])
        # and then treat the keep/remove
        for label_cap in instance.metadata_menu[menu_type]['navigation_buttons']['labels']:
            label = label_cap.lower()
            if hasattr(instance, label):
                keeps, removes = split_options(getattr(instance, label))
                instance.metadata_menu[menu_type][label_cap]['checkboxes']['keep_selected'] = keeps
                instance.metadata_menu[menu_type][label_cap]['checkboxes']['remove_selected'] = removes


def update_metadata_fields(instance):
    """Update the metadata menu object with metadata associated with newly read data
    for non-numeric metadata gets all the unique fields per metadata variable,
    and sets the available fields as such, and for numeric gets the minimum and maximum
    boundaries of each metadata variable.
    If previously metadata settings for a field deviate from the default,
    then if the same field still exists then the settings (i.e. bounds or
    checkbox selection) are copied across, rather than setting to the default.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    # iterate through metadata variables
    for meta_var in instance.metadata_vars_to_read:

        meta_var_field = instance.datareader.metadata_in_memory[meta_var]
        if instance.reading_nonghost and meta_var == 'latitude':
            meta_var_field = instance.datareader.nonghost_metadata[meta_var]
        if instance.reading_nonghost and meta_var == 'longitude':
            meta_var_field = instance.datareader.nonghost_metadata[meta_var]

        # get metadata variable type/data type
        metadata_type = instance.standard_metadata[meta_var]['metadata_type']
        metadata_data_type = instance.standard_metadata[meta_var]['data_type']

        # remove NaNs from field
        meta_var_field_nan_removed = meta_var_field[~pd.isnull(meta_var_field)]

        # update pop-up metadata menu object with read metadata values
        # for non-numeric metadata gets all the unique fields per metadata variable
        # and sets the available fields as such
        if metadata_data_type == np.object:
            # get previous fields
            previous_fields = copy.deepcopy(instance.metadata_menu[metadata_type][meta_var]['checkboxes']['labels'])
            # update new labels
            instance.metadata_menu[metadata_type][meta_var]['checkboxes']['labels'] = np.unique(
                meta_var_field_nan_removed)
            # if field previously existed, then copy across checkbox settings for field
            # else set initial checkboxes to be all blank
            previous_keep = copy.deepcopy(instance.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected'])
            previous_remove = copy.deepcopy(instance.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_selected'])
            instance.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected'] = []
            instance.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_selected'] = []
            for field in instance.metadata_menu[metadata_type][meta_var]['checkboxes']['labels']:
                if field in previous_fields:
                    if field in previous_keep:
                        instance.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected'].append(field)
                    if field in previous_remove:
                        instance.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_selected'].append(field)
            # set defaults to be empty
            instance.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_default'] = []
            instance.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_default'] = []
        # for numeric fields get the minimum and maximum boundaries of each metadata variable
        # if previous set values vary from min/max boundaries, copy across the values
        # set as min/max as nan if have no numeric metadata for variable
        else:
            meta_var_index = instance.metadata_menu[metadata_type]['rangeboxes']['labels'].index(meta_var)
            # have some numeric values for metadata variable?
            if len(meta_var_field_nan_removed) > 0:
                min_val = str(np.min(meta_var_field_nan_removed))
                max_val = str(np.max(meta_var_field_nan_removed))
                # get previous lower/upper extents and defaults
                previous_lower_default = copy.deepcopy(instance.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index])
                previous_upper_default = copy.deepcopy(instance.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index])
                previous_lower = copy.deepcopy(instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index])
                previous_upper = copy.deepcopy(instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index])
                # if previous lower > previous default lower bound then copy across (and also not 'nan')
                # initially set as min extent
                instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = min_val
                if (previous_lower != 'nan') & (previous_lower_default != 'nan'):
                    if previous_lower > previous_lower_default:
                        instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = copy.deepcopy(previous_lower)
                # if previous upper < previous default upper bound then copy across (and also not 'nan')
                # initially set as max extent
                instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = max_val
                if (previous_upper != 'nan') & (previous_upper_default != 'nan'):
                    if previous_upper < previous_upper_default:
                        instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = copy.deepcopy(previous_upper)
                # set defaults to min/max extents
                instance.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index] = min_val
                instance.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index] = max_val
            # do not have some numeric values for metadata variable so set as 'nan'
            else:
                instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = 'nan'
                instance.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index] = 'nan'
                instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = 'nan'
                instance.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index] = 'nan'


def representativity_fields(instance, resolution):
    """Update the data representativity menu -> 1D list of rangebox values
    dependent on the temporal resolution, some fields will appear or not.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    :param resolution: selected resolution (e.g. "hourly", "3hourly")
    :type resolution: str
    :return: menu of representativity
    :rtype: dict
    """

    # get previous set labels of dashboard
    if instance.offline:
        repr_menu = {'rangeboxes': {'labels': [], 'current_lower': []}}
        previous_labels = []
    else:
        previous_labels = copy.deepcopy(instance.representativity_menu['rangeboxes']['labels'])
        # get previously set rangebox values
        previous_lower = copy.deepcopy(instance.representativity_menu['rangeboxes']['current_lower'])
        repr_menu = instance.representativity_menu

    # hourly temporal resolution?
    if (resolution == 'hourly') or (resolution == 'hourly_instantaneous'):
        repr_menu['rangeboxes']['labels'] = ['hourly_native_representativity_percent',
                                             'hourly_native_max_gap_percent',
                                             'daily_native_representativity_percent',
                                             'daily_representativity_percent',
                                             'daily_native_max_gap_percent',
                                             'daily_max_gap_percent',
                                             'monthly_native_representativity_percent',
                                             'monthly_representativity_percent',
                                             'monthly_native_max_gap_percent',
                                             'monthly_max_gap_percent',
                                             'all_representativity_percent', 'all_max_gap_percent']
    # daily temporal resolution?
    elif (resolution == 'daily') or (resolution == '3hourly') or \
            (resolution == '6hourly') or (resolution == '3hourly_instantaneous') or \
            (resolution == '6hourly_instantaneous'):
        repr_menu['rangeboxes']['labels'] = ['daily_native_representativity_percent',
                                             'daily_native_max_gap_percent',
                                             'monthly_native_representativity_percent',
                                             'monthly_representativity_percent',
                                             'monthly_native_max_gap_percent',
                                             'monthly_max_gap_percent',
                                             'all_representativity_percent', 'all_max_gap_percent']
    # monthly temporal resolution?
    elif resolution == 'monthly':
        repr_menu['rangeboxes']['labels'] = ['monthly_native_representativity_percent',
                                             'monthly_native_max_gap_percent',
                                             'all_representativity_percent', 'all_max_gap_percent']

    # initialise rangebox values --> for data representativity fields
    # the default is 0%, for max gap fields % the default is 100%
    repr_menu['rangeboxes']['current_lower'] = []
    for label_ii, label in enumerate(repr_menu['rangeboxes']['labels']):
        if 'max_gap' in label:
            repr_menu['rangeboxes']['current_lower'].append('100')
        else:
            repr_menu['rangeboxes']['current_lower'].append('0')

        # label previously existed?
        if label in previous_labels:
            repr_menu['rangeboxes']['current_lower'][label_ii] = \
                previous_lower[previous_labels.index(label)]

    return repr_menu


def update_period_fields(resolution, period_menu):
    """Update the data period menu -> list of checkboxes
    dependent on the temporal resolution, some fields will appear or not.

    :param resolution: selected resolution (e.g. "hourly", "3hourly")
    :type resolution: str
    :param period_menu: contains options of period
    :type period_menu: dict
    :return: updated menu of period
    :rtype: dict
    """
    # hourly temporal resolution?
    if 'hourly' in resolution:
        period_menu['checkboxes']['labels'] = ['Daytime', 'Nighttime', 'Weekday', 'Weekend',
                                               'Spring', 'Summer', 'Autumn', 'Winter']
    # daily temporal resolution?
    elif resolution == 'daily':
        period_menu['checkboxes']['labels'] = ['Weekday', 'Weekend', 'Spring',
                                               'Summer', 'Autumn', 'Winter']

        # drop selected fields from higher temporal resolutions
        labels_to_remove = ['Daytime', 'Nighttime']
        for label in labels_to_remove:
            if label in period_menu['checkboxes']['keep_selected']:
                period_menu['checkboxes']['keep_selected'].remove(label)
            if label in period_menu['checkboxes']['remove_selected']:
                period_menu['checkboxes']['remove_selected'].remove(label)

    # monthly temporal resolution?
    elif resolution == 'monthly':
        period_menu['checkboxes']['labels'] = ['Spring', 'Summer', 'Autumn', 'Winter']
        # drop selected fields from higher temporal resolutions
        labels_to_remove = ['Daytime', 'Nighttime', 'Weekday', 'Weekend']
        for label in labels_to_remove:
            if label in period_menu['checkboxes']['keep_selected']:
                period_menu['checkboxes']['keep_selected'].remove(label)
            if label in period_menu['checkboxes']['remove_selected']:
                period_menu['checkboxes']['remove_selected'].remove(label)


def to_pandas_dataframe(instance, species):
    """Function that takes selected station data within
    data arrays and puts it into a pandas dataframe.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    :param species: The species current selected for evaluation (e.g. sconco3)
    :type species: str
    :return: station data by data array
    :rtype: dict
    """

    # create new dictionary to store selection station data by data array
    selected_station_data = {}

    # iterate through data arrays in data in memory filtered dictionary
    for data_label in list(instance.data_in_memory_filtered.keys()):

        # if colocation is not active, do not convert colocated data arrays to pandas data frames
        if not instance.colocate_active:
            if 'colocated' in data_label:
                continue
        # else, if colocation is active, do not convert non-colocated data arrays to pandas data frames
        elif instance.colocate_active:
            if 'colocated' not in data_label:
                continue

        # observational arrays
        if data_label.split('_')[0] == 'observations':
            # get data for selected stations
            data_array = instance.data_in_memory_filtered[data_label][
                             species][instance.relative_selected_station_inds, :]
        # experiment arrays
        else:
            # get intersect between selected station indices and valid available indices for experiment data array
            valid_selected_station_indices = np.intersect1d(instance.relative_selected_station_inds,
                                                            instance.datareader.plotting_params[
                                                                data_label]['valid_station_inds'])
            # get data for valid selected stations
            data_array = instance.data_in_memory_filtered[data_label][species][valid_selected_station_indices, :]

        # if data array has no valid data for selected stations, do not create a pandas dataframe
        # data array has valid data?
        if data_array.size:
            # add nested dictionary for data array name to selection station data dictionary
            selected_station_data[data_label] = {}
            # take cross station median of selected data for data array, and place it in a pandas
            # dataframe -->  add to selected station data dictionary
            selected_station_data[data_label]['pandas_df'] = pd.DataFrame(np.nanmedian(data_array, axis=0),
                                                                          index=instance.time_array,
                                                                          columns=['data'])

    return selected_station_data


def valid_date(date_text):
    """Determines if a date string is in the correct format."""

    try:
        datetime.datetime.strptime(str(date_text), '%Y%m%d')
        return True
    except Exception as e:
        return False


def check_for_ghost(network_name):
    """It checks whether the selected network comes from GHOST or not.
    All non-GHOST networks start with an asterisk at their name."""

    if '*' in network_name:
        return True
    else:
        return False
