import copy
import json
import numpy as np

from .config import split_options


def which_bounds(instance, species):
    """If there are bounds defined in a config file, fill that value,
    if it is withing the feasible bounds of the species"""

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
    """Checks if the species we currently have selected belongs to the ones
    that have specific qa flags selected as default"""

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
    """if there are flags coming from a config file, select those"""

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
    """Retrieve QA codes from GHOST_standards using the qa flags' names.

    Specific flags are defined for the following species:
    ['WND_DIR_10M','WND_SPD_10M','RH_2M','PREC_ACCUM','SNOW_ACCUM',
    'SNOW_DEPTH','CEILING_HEIGHT','VIS_DIST','CLOUD_CVG','CLOUD_CVG_FRAC']"""

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


def representativity_conf(instance):
    """Comes here if there is a configuration loaded. Checks if there is a
    representative field loaded in the object from the conf and if there is
    assigns the value in the representativity menu"""
    for i, label in enumerate(instance.representativity_menu['rangeboxes']['labels']):
        if hasattr(instance, label):
            instance.representativity_menu['rangeboxes']['current_lower'][i] = str(getattr(instance, label))


def init_metadata(instance):
    """Initialise internal structure to store metadata."""

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
    """Comes here if there in a loaded configuration there are also metadata fields."""

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


def representativity_fields(instance, resolution):
    """Update the data representativity menu -> 1D list of rangebox values
       dependent on the temporal resolution, some fields will appear or not
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
