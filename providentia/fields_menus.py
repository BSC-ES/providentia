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
    if not hasattr(instance, 'flag_menu'):
        instance.flag_menu = {'window_title':'Network QA Flags', 
                              'page_title':'Select standardised data reporter provided flags to filter by', 
                              'checkboxes':{}}
        instance.flag_menu['select_buttons'] = ['all', 'clear', 'default']
    
    # reset fields
    instance.flag_menu['checkboxes']['labels'] = np.array(sorted(instance.standard_data_flag_name_to_data_flag_code, 
                                                                 key=instance.standard_data_flag_name_to_data_flag_code.get))
    instance.flag_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
    instance.flag_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
    instance.flag_menu['checkboxes']['map_vars'] = np.sort(list(instance.standard_data_flag_name_to_data_flag_code.values()))
    

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
    if not hasattr(instance, 'qa_menu'):
        instance.qa_menu = {'window_title':'GHOST QA Flags', 
                            'page_title':'Select standardised quality assurance flags to filter by', 
                            'checkboxes':{}}
        instance.qa_menu['select_buttons'] = ['all', 'clear', 'default']
    
    # reset fields
    if instance.network == ['actris/actris']:
        qa_labels = {'Invalid Data Provider Flags - GHOST Decreed': 6, 
                     'Invalid Data Provider Flags - Network Decreed': 7}
    else:
        qa_labels = instance.standard_QA_name_to_QA_code

    qa_key = qa_labels.get
    qa_map_vars = qa_labels.values()
    instance.qa_menu['checkboxes']['labels'] = np.array(sorted(qa_labels, 
                                                               key=qa_key))
    instance.qa_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
    instance.qa_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
    instance.qa_menu['checkboxes']['map_vars'] = np.sort(list(qa_map_vars))
    

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
    if not hasattr(instance, 'models_menu'):
        instance.models_menu = {'window_title': 'Models', 
                                     'page_title': 'Select Models', 
                                     'models': {}}
        instance.models_menu['select_buttons'] = ['all', 'clear']
    
    # reset fields
    instance.models_menu['models']['labels'] = [] 
    instance.models_menu['models']['keep_selected'] = [] 
    instance.models_menu['models']['forecast'] = {} 
    instance.models_menu['models']['forecast_days'] = {}
    instance.models_menu['models']['map_vars'] = [] 


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
    if not hasattr(instance, 'multispecies_menu'):
        instance.multispecies_menu = {'window_title': 'Species Filtering', 
                                      'page_title': 'Select species to filter by',
                                      'multispecies': {},
                                     }

    # reset rangeboxes
    instance.multispecies_menu['multispecies']['labels'] = []
    instance.multispecies_menu['multispecies']['current_lower'] = {}
    instance.multispecies_menu['multispecies']['current_upper'] = {}
    instance.multispecies_menu['multispecies']['current_filter_species_fill_value'] = {}
    instance.multispecies_menu['multispecies']['apply_selected'] = {}
    instance.multispecies_menu['multispecies']['previous_lower'] = {}
    instance.multispecies_menu['multispecies']['previous_upper'] = {}
    instance.multispecies_menu['multispecies']['previous_apply'] = {}
    instance.multispecies_menu['multispecies']['previous_filter_species_fill_value'] = {}


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
    if not hasattr(instance, 'coverage_menu'):
        instance.coverage_menu = {'window_title': '% Data Coverage', 
                                          'page_title': 'Select Minimum Required % Data Coverage', 
                                          'rangeboxes':{}}
    
    # reset fields
    instance.coverage_menu['rangeboxes']['tooltips'] = []
    instance.coverage_menu['rangeboxes']['labels'] = [] 
    instance.coverage_menu['rangeboxes']['current_lower'] = []    
    instance.coverage_menu['rangeboxes']['map_vars_old'] = []                                              
    instance.coverage_menu['rangeboxes']['map_vars'] = []
    instance.coverage_menu['rangeboxes']['subtitles'] = []
    instance.coverage_menu['rangeboxes']['subtitle_inds'] = []


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
    if not hasattr(instance, 'period_menu'):
        instance.period_menu = {'window_title': 'Data Period', 
                                'page_title': 'Select Data Periods', 
                                'checkboxes':{}}
    
    # reset fields
    instance.period_menu['checkboxes']['labels'] = []
    instance.period_menu['checkboxes']['keep_selected'] = []
    instance.period_menu['checkboxes']['remove_selected'] = []


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
    if not hasattr(instance, 'metadata_menu'):
        instance.metadata_types = {'STATION POSITION': 'Filter stations by measurement position',
                                   'STATION CLASSIFICATIONS': 'Filter stations by station provided classifications',
                                   'STATION MISCELLANEOUS': 'Filter stations by miscellaneous station provided metadata',
                                   'GLOBALLY GRIDDED CLASSIFICATIONS': 'Filter stations by globally gridded classifications',
                                   'MEASUREMENT PROCESS INFORMATION': 'Filter stations by measurement process information'}
            
        instance.metadata_menu = {'window_title': 'Metadata', 
                                  'page_title': 'Select metadata type to filter stations by',
                                  'navigation_buttons': {}}

        instance.metadata_menu['navigation_buttons']['labels'] = list(instance.metadata_types.keys())
        instance.metadata_menu['navigation_buttons']['tooltips'] = [instance.metadata_types[key] for key in
                                                                    instance.metadata_menu['navigation_buttons']['labels']]

        for metadata_type_ii, metadata_type in enumerate(instance.metadata_menu['navigation_buttons']['labels']):
            
            # setup nested menu
            instance.metadata_menu[metadata_type] = {'window_title': metadata_type,
                                                    'page_title': instance.metadata_menu['navigation_buttons']['tooltips'][
                                                     metadata_type_ii], 'navigation_buttons': {}, 'rangeboxes': {}}
        
    # reset fields
    for metadata_type_ii, metadata_type in enumerate(instance.metadata_menu['navigation_buttons']['labels']):

        # reset rangebox labels    
        instance.metadata_menu[metadata_type]['rangeboxes']['labels'] = \
            [metadata_var for metadata_var in instance.metadata_vars_to_read
            if (instance.standard_metadata[metadata_var]['metadata_type'] == metadata_type)
            & (instance.standard_metadata[metadata_var]['data_type'] != object)]

        # reset rangebox tooltips
        instance.metadata_menu[metadata_type]['rangeboxes']['tooltips'] = \
            [instance.standard_metadata[metadata_var]['description']
            for metadata_var in instance.metadata_menu[metadata_type]['rangeboxes']['labels']]

        # reset rangeboxes
        instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'] = \
            ['nan'] * len(instance.metadata_menu[metadata_type]['rangeboxes']['labels'])
        instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'] = \
            ['nan'] * len(instance.metadata_menu[metadata_type]['rangeboxes']['labels'])
        instance.metadata_menu[metadata_type]['rangeboxes']['lower_default'] = \
            ['nan'] * len(instance.metadata_menu[metadata_type]['rangeboxes']['labels'])
        instance.metadata_menu[metadata_type]['rangeboxes']['upper_default'] = \
            ['nan'] * len(instance.metadata_menu[metadata_type]['rangeboxes']['labels'])
        instance.metadata_menu[metadata_type]['rangeboxes']['apply_selected'] = []

        # reset checkbox labels            
        instance.metadata_menu[metadata_type]['navigation_buttons']['labels'] = \
            [metadata_var for metadata_var in instance.metadata_vars_to_read if
            (instance.standard_metadata[metadata_var]['metadata_type'] == metadata_type) &
            (instance.standard_metadata[metadata_var]['data_type'] == object)]
        
        # reset checkbox tooltips
        instance.metadata_menu[metadata_type]['navigation_buttons']['tooltips'] = \
            [instance.standard_metadata[metadata_var]['description'] for metadata_var in
            instance.metadata_menu[metadata_type]['navigation_buttons']['labels']]

        # reset checkboxes
        for metadata_var in instance.metadata_menu[metadata_type]['navigation_buttons']['labels']:
            # metadata variable already in dict?
            # then just reset lists
            if metadata_var in instance.metadata_menu[metadata_type]:
                instance.metadata_menu[metadata_type][metadata_var]['checkboxes']['labels'] = []
                instance.metadata_menu[metadata_type][metadata_var]['checkboxes']['keep_selected'] = []
                instance.metadata_menu[metadata_type][metadata_var]['checkboxes']['keep_default'] = []
                instance.metadata_menu[metadata_type][metadata_var]['checkboxes']['remove_selected'] = []
                instance.metadata_menu[metadata_type][metadata_var]['checkboxes']['remove_default'] = []
            # otherwise, create infrastructure to store metadata var information
            else:
                instance.metadata_menu[metadata_type][metadata_var] = {'window_title': metadata_var,
                                                                       'page_title': 'Filter stations by unique {} metadata'.format(metadata_var), 
                                                                       'checkboxes': {}}
                instance.metadata_menu[metadata_type][metadata_var]['checkboxes'] = {'labels': [], 
                                                                                     'keep_selected': [], 'keep_default': [],
                                                                                     'remove_selected': [], 'remove_default': []}
        
        # remove metadata type checkbox var dicts not in metadata_vars_to_read
        metadata_type_checkbox_vars = [metadata_type_var for metadata_type_var in instance.metadata_menu[metadata_type].keys() if metadata_type_var not in ['window_title', 'page_title', 'navigation_buttons', 'rangeboxes']]
        metadata_type_checkbox_vars_to_remove = [metadata_type_checkbox_var for metadata_type_checkbox_var in metadata_type_checkbox_vars if metadata_type_checkbox_var not in instance.metadata_vars_to_read]
        for metadata_type_checkbox_var_to_remove in metadata_type_checkbox_vars_to_remove:
            del instance.metadata_menu[metadata_type][metadata_type_checkbox_var_to_remove]


def update_qa(instance):
    """
    Update the quality assurance menu labels and mappings.

    Parameters
    ----------
    instance : object
        An object instance to be updated with quality assurance menu configurations.
    """

    # reset fields
    if instance.selected_network == 'actris/actris':
        qa_labels = {'Invalid Data Provider Flags - GHOST Decreed': 6, 
                     'Invalid Data Provider Flags - Network Decreed': 7}
    else:
        qa_labels = instance.standard_QA_name_to_QA_code

    qa_key = qa_labels.get
    qa_map_vars = qa_labels.values()
    instance.qa_menu['checkboxes']['labels'] = np.array(sorted(qa_labels, 
                                                               key=qa_key))
    instance.qa_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
    instance.qa_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
    instance.qa_menu['checkboxes']['map_vars'] = np.sort(list(qa_map_vars))


def update_coverage_fields(instance):
    """
    Update the data coverage menu labels and values based on the active data resolution and network type.

    Parameters
    ----------
    instance : object
        An object instance to be updated with coverage menu configurations.
    """

    # get previously set rangebox labels / values
    previous_mapped_labels = copy.deepcopy(instance.coverage_menu['rangeboxes']['map_vars'])
    previous_mapped_labels_old = copy.deepcopy(instance.coverage_menu['rangeboxes']['map_vars_old'])
    previous_lower = copy.deepcopy(instance.coverage_menu['rangeboxes']['current_lower'])
    
    # build coverage menu
    if (instance.reading_ghost) & (instance.ghost_features == 'max'):
        network_type = 'ghost'
    else:
        network_type = 'nonghost'

    # set resolution to get coverage information based on active resolution
    if instance.resolution in ['hourly', 'hourly_instantaneous']:
        resolution = 'hourly'
    elif instance.resolution in ['3hourly', '3hourly_instantaneous']:
        resolution = '3hourly'
    elif instance.resolution in ['6hourly', '6hourly_instantaneous']:
        resolution = '6hourly'
    elif instance.resolution == 'daily':
        resolution = 'daily'
    elif instance.resolution == 'monthly':
        resolution = 'monthly'
    elif instance.resolution == 'annual':
        resolution = 'annual'

    instance.coverage_menu['rangeboxes']['map_vars'] = instance.coverage_info[network_type][resolution]['map_vars']
    instance.coverage_menu['rangeboxes']['map_vars_old'] = instance.coverage_info[network_type][resolution]['map_vars_old']
    instance.coverage_menu['rangeboxes']['labels'] = instance.coverage_info[network_type][resolution]['labels']
    instance.coverage_menu['rangeboxes']['subtitles'] = instance.coverage_info[network_type][resolution]['subtitles']
    instance.coverage_menu['rangeboxes']['subtitle_inds'] = instance.coverage_info[network_type][resolution]['subtitle_inds']
    
    # initialise rangebox values --> for data coverage fields
    # the default is 0%, for max gap fields % the default is 100%
    instance.coverage_menu['rangeboxes']['current_lower'] = []
    for label_ii, (label_mapped, label_mapped_old) in enumerate(zip(instance.coverage_menu['rangeboxes']['map_vars'], instance.coverage_menu['rangeboxes']['map_vars_old'])):
        if 'max_gap' in label_mapped:
            instance.coverage_menu['rangeboxes']['current_lower'].append('100')
        else:
            instance.coverage_menu['rangeboxes']['current_lower'].append('0')

        # label previously existed?
        if (label_mapped in previous_mapped_labels) or (label_mapped_old in previous_mapped_labels_old):
            instance.coverage_menu['rangeboxes']['current_lower'][label_ii] = \
                previous_lower[previous_mapped_labels.index(label_mapped)]


def update_period_fields(instance):
    """
    Update the data period menu labels based on the active temporal resolution.

    Parameters
    ----------
    instance : object
        An object instance to be updated with period menu configurations.
    """

    # hourly/3hourly/6hourly temporal resolution?
    if 'hourly' in instance.resolution:
        instance.period_menu['checkboxes']['labels'] = ['Daytime', 'Nighttime', 'Weekday', 'Weekend',
                                                        'Spring', 'Summer', 'Autumn', 'Winter']
    # daily temporal resolution?
    elif instance.resolution == 'daily':
        instance.period_menu['checkboxes']['labels'] = ['Weekday', 'Weekend', 
                                                        'Spring', 'Summer', 'Autumn', 'Winter']
        # drop selected fields from higher temporal resolutions
        labels_to_remove = ['Daytime', 'Nighttime']
        for label in labels_to_remove:
            if label in instance.period_menu['checkboxes']['keep_selected']:
                instance.period_menu['checkboxes']['keep_selected'].remove(label)
            if label in instance.period_menu['checkboxes']['remove_selected']:
                instance.period_menu['checkboxes']['remove_selected'].remove(label)

    # monthly temporal resolution?
    elif instance.resolution == 'monthly':
        instance.period_menu['checkboxes']['labels'] = ['Spring', 'Summer', 'Autumn', 'Winter']
        # drop selected fields from higher temporal resolutions
        labels_to_remove = ['Daytime', 'Nighttime', 'Weekday', 'Weekend']
        for label in labels_to_remove:
            if label in instance.period_menu['checkboxes']['keep_selected']:
                instance.period_menu['checkboxes']['keep_selected'].remove(label)
            if label in instance.period_menu['checkboxes']['remove_selected']:
                instance.period_menu['checkboxes']['remove_selected'].remove(label)

    # annual temporal resolution?
    elif instance.resolution == 'annual':
        instance.period_menu['checkboxes']['labels'] = []
        # drop selected fields from higher temporal resolutions
        labels_to_remove = ['Daytime', 'Nighttime', 'Weekday', 'Weekend', 'Spring', 'Summer', 'Autumn', 'Winter']
        for label in labels_to_remove:
            if label in instance.period_menu['checkboxes']['keep_selected']:
                instance.period_menu['checkboxes']['keep_selected'].remove(label)
            if label in instance.period_menu['checkboxes']['remove_selected']:
                instance.period_menu['checkboxes']['remove_selected'].remove(label)

def update_metadata_fields(instance):
    """
    Update the metadata menu with unique categorical fields or numeric boundaries 
    derived from newly read data.
    
    Non-numeric metadata gets all the unique fields per metadata variable.
    Numeric metadata gets the minimum and maximum boundaries per metadata variable.
    If previously metadata settings for a field deviate from the default,
    then if the same field still exists then the settings (i.e. bounds or
    checkbox selection) are copied across, rather than setting to the default.

    Parameters
    ----------
    instance : object
        An object instance to be updated with metadata menu configurations.
    """

    metadata_menu = instance.metadata_menu
    metadata_vars = instance.metadata_vars_to_read
    standard_metadata = instance.standard_metadata
    metadata_in_memory = instance.metadata_in_memory
    networkspecies = instance.networkspecies

    # before doing anything check if have all current metadata variables in menu
    # if not, this is either because it is initialised empty, 
    # or because of changing between GHOST and non-GHOST data
    # if there is a difference, fill metadata menu with empty values
    reset_meta = False
    nav_labels = metadata_menu['navigation_buttons']['labels']

    for metadata_type in nav_labels:
        required_vars = [
            mv for mv in metadata_vars
            if standard_metadata[mv]['metadata_type'] == metadata_type
            and standard_metadata[mv]['data_type'] != object
        ]
        if len(required_vars) != len(metadata_menu[metadata_type]['rangeboxes']['labels']):
            reset_meta = True
            break

    # reinitialise metadata menu
    if reset_meta:
        init_metadata(instance)
        metadata_menu = instance.metadata_menu

    # normalise network/species to lists (avoid conversion inside loop)
    if isinstance(instance.network, str):
        instance.network = [instance.network]
    if isinstance(instance.species, str):
        instance.species = [instance.species]

    # main metadata update loop
    for meta_var in metadata_vars:
        md_info = standard_metadata[meta_var]
        metadata_type = md_info['metadata_type']
        dtype = md_info['data_type']

        # gather metadata across all networks/species
        chunks = []

        for netspec in networkspecies:
            arr = metadata_in_memory[netspec][meta_var]

            # apply valid station mask (non-GHOST only)
            if hasattr(instance, "valid_station_inds") and not instance.reading_ghost:
                if instance.temporal_colocation:
                    inds = instance.valid_station_inds_temporal_colocation[netspec][instance.observations_data_label]
                else:
                    inds = instance.valid_station_inds[netspec][instance.observations_data_label]
                chunks.append(arr[inds].ravel())
            else:
                chunks.append(arr.ravel())

        # one fast concatenate
        meta_vals = np.concatenate(chunks)

        # remove NaNs / "nan"
        if dtype == object:
            # object/string metadata
            mask = (meta_vals != 'nan') & (meta_vals != "") & (~pd.isna(meta_vals))
        else:
            # numeric metadata
            mask = ~np.isnan(meta_vals)

        meta_vals_clean = meta_vals[mask]

        # update metadata menu for categorical (object) fields
        if dtype == object:
            mm = metadata_menu[metadata_type][meta_var]['checkboxes']

            # get previous settings
            prev_fields = mm['labels']
            prev_set = set(prev_fields)
            prev_keep = set(mm['keep_selected'])
            prev_remove = set(mm['remove_selected'])

            # new unique sorted values
            new_fields = np.unique(meta_vals_clean)
            mm['labels'] = new_fields

            # rebuild keep/remove efficiently
            keep = []
            remove = []
            for f in new_fields:
                if f in prev_set:
                    if f in prev_keep:
                        keep.append(f)
                    if f in prev_remove:
                        remove.append(f)

            mm['keep_selected'] = keep
            mm['remove_selected'] = remove
            mm['keep_default'] = []
            mm['remove_default'] = []
            continue

        # update numeric metadata fields (rangeboxes)
        mm = metadata_menu[metadata_type]['rangeboxes']
        idx = mm['labels'].index(meta_var)

        prev_lower_def = mm['lower_default'][idx]
        prev_upper_def = mm['upper_default'][idx]
        prev_lower = mm['current_lower'][idx]
        prev_upper = mm['current_upper'][idx]

        if len(meta_vals_clean) > 0:
            mn = str(np.min(meta_vals_clean))
            mx = str(np.max(meta_vals_clean))

            # apply previous settings if stricter
            new_lower = mn
            if prev_lower not in ("nan", None) and prev_lower_def not in ("nan", None):
                if prev_lower > prev_lower_def:
                    new_lower = prev_lower

            new_upper = mx
            if prev_upper not in ("nan", None) and prev_upper_def not in ("nan", None):
                if prev_upper < prev_upper_def:
                    new_upper = prev_upper

            mm['current_lower'][idx] = new_lower
            mm['current_upper'][idx] = new_upper
            mm['lower_default'][idx] = mn
            mm['upper_default'][idx] = mx

        else:
            # No numeric metadata → set NAN
            mm['current_lower'][idx] = 'nan'
            mm['lower_default'][idx] = 'nan'
            mm['current_upper'][idx] = 'nan'
            mm['upper_default'][idx] = 'nan'

            # remove "apply selected" if invalid
            if meta_var in mm['apply_selected']:
                mm['apply_selected'].remove(meta_var)

def multispecies_conf(instance):
    """
    Configure multispecies filtering variables and widget states using settings loaded from a configuration file.

    Parameters
    ----------
    instance : object
        An object instance to be updated with multispecies menu configurations.
    """

    if hasattr(instance, 'filter_species'):
        filter_species = copy.deepcopy(instance.filter_species)
        for (networkspeci_ii, networkspeci), networkspeci_bounds in zip(enumerate(filter_species.keys()),
                                                                        filter_species.values()):

            for bounds in networkspeci_bounds:
                # update menu_current
                if ('networkspeci_' + str(networkspeci_ii)) not in instance.multispecies_menu['multispecies']['labels']:
                    instance.multispecies_menu['multispecies']['labels'].append('networkspeci_' + str(networkspeci_ii))

                # add values
                instance.multispecies_menu['multispecies']['current_lower'][networkspeci_ii] = bounds[0]
                instance.multispecies_menu['multispecies']['current_upper'][networkspeci_ii] = bounds[1]
                instance.multispecies_menu['multispecies']['current_filter_species_fill_value'][networkspeci_ii] = bounds[2]
                instance.multispecies_menu['multispecies']['apply_selected'][networkspeci_ii] = True

                # set initial selected config variables as set .conf files or defaults
                instance.selected_widget_network.update({networkspeci_ii: networkspeci.split('|')[0]})
                instance.selected_widget_matrix.update({networkspeci_ii: instance.parameter_dictionary[networkspeci.split('|')[1]]['matrix']})
                instance.selected_widget_species.update({networkspeci_ii: networkspeci.split('|')[1]})
                instance.selected_widget_lower.update({networkspeci_ii: bounds[0]})
                instance.selected_widget_upper.update({networkspeci_ii: bounds[1]})
                instance.selected_widget_filter_species_fill_value.update({networkspeci_ii: bounds[2]})
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

    for i, (label, label_old) in enumerate(zip(instance.coverage_menu['rangeboxes']['map_vars'], instance.coverage_menu['rangeboxes']['map_vars_old'])):
        if hasattr(instance, label):
            instance.coverage_menu['rangeboxes']['current_lower'][i] = str(getattr(instance, label))
        elif hasattr(instance, label_old):
            instance.coverage_menu['rangeboxes']['current_lower'][i] = str(getattr(instance, label_old))

def period_conf(instance):
    """
    Configure period filter variables and checkbox selections using settings loaded from a configuration file.

    Parameters
    ----------
    instance : object
        An object instance to be updated with period menu configurations.
    """

    from .configuration import split_options

    if hasattr(instance, 'period'):
        keeps, removes = split_options(instance, instance.period)
        instance.period_menu['checkboxes']['keep_selected'] = keeps
        instance.period_menu['checkboxes']['remove_selected'] = removes


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
        for i, label in enumerate(instance.metadata_menu[menu_type]['rangeboxes']['labels']):
            if hasattr(instance, label):
                current_lower.append(str(getattr(instance, label)[0]))
                current_upper.append(str(getattr(instance, label)[1]))
                apply_selected.append(label)
            else:
                current_lower.append(instance.metadata_menu[menu_type]['rangeboxes']['current_lower'][i])
                current_upper.append(instance.metadata_menu[menu_type]['rangeboxes']['current_upper'][i])
        instance.metadata_menu[menu_type]['rangeboxes']['current_lower'] = current_lower
        instance.metadata_menu[menu_type]['rangeboxes']['current_upper'] = current_upper
        instance.metadata_menu[menu_type]['rangeboxes']['apply_selected'] = apply_selected

        # and then treat the keep/remove
        for label in instance.metadata_menu[menu_type]['navigation_buttons']['labels']:
            if hasattr(instance, label):
                keeps, removes = split_options(instance, getattr(instance, label))
                instance.metadata_menu[menu_type][label]['checkboxes']['keep_selected'] = keeps
                instance.metadata_menu[menu_type][label]['checkboxes']['remove_selected'] = removes