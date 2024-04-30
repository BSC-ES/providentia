""" Functions to initialise and update fields (flags, period, metadata, etc.) """

import copy

import numpy as np
import pandas as pd


def init_flags(instance):
    """ Initialise internal structure to store selected flags.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'flag_menu'):
        instance.flag_menu = {'window_title':'FLAGS', 
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
    """ Initialise internal structure to store selected qa.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'qa_menu'):
        instance.qa_menu = {'window_title':'QA', 
                            'page_title':'Select standardised quality assurance flags to filter by', 
                            'checkboxes':{}}
        instance.qa_menu['select_buttons'] = ['all', 'clear', 'default']
    
    # reset fields
    instance.qa_menu['checkboxes']['labels'] = np.array(sorted(instance.standard_QA_name_to_QA_code, 
                                                               key=instance.standard_QA_name_to_QA_code.get))
    instance.qa_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
    instance.qa_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
    instance.qa_menu['checkboxes']['map_vars'] = np.sort(list(instance.standard_QA_name_to_QA_code.values()))
    

def init_experiments(instance):
    """ Initialise internal structure to store selected experiments.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'experiments_menu'):
        instance.experiments_menu = {'window_title': 'EXPERIMENTS', 
                                     'page_title': 'Select Experiment/s', 
                                     'checkboxes':{}}
        instance.experiments_menu['select_buttons'] = ['all', 'clear']
    
    # reset fields
    instance.experiments_menu['checkboxes']['labels'] = [] 
    instance.experiments_menu['checkboxes']['keep_default'] = [] 
    instance.experiments_menu['checkboxes']['keep_selected'] = [] 
    instance.experiments_menu['checkboxes']['map_vars'] = [] 


def init_multispecies(instance):
    """ Initialise internal structure to store multispecies fields.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'multispecies_menu'):
        instance.multispecies_menu = {'window_title': 'MULTISPECIES', 
                                      'page_title': 'Select Network/s and Specie/s to Filter by',
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


def init_representativity(instance):
    """ Initialise internal structure to store representativity fields.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'representativity_menu'):
        instance.representativity_menu = {'window_title': '% DATA REPRESENTATIVITY', 
                                          'page_title': 'Select Minimum Required % Data Representativity', 
                                          'rangeboxes':{}}
    
    # reset fields
    instance.representativity_menu['rangeboxes']['tooltips'] = []
    instance.representativity_menu['rangeboxes']['labels'] = [] 
    instance.representativity_menu['rangeboxes']['current_lower'] = []                                                  
    instance.representativity_menu['rangeboxes']['map_vars'] = []
    instance.representativity_menu['rangeboxes']['subtitles'] = []
    instance.representativity_menu['rangeboxes']['subtitle_inds'] = []


def init_period(instance):
    """ Initialise internal structure to store period fields.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'period_menu'):
        instance.period_menu = {'window_title': 'DATA PERIOD', 
                                'page_title': 'Select Data Periods', 
                                'checkboxes':{}}
    
    # reset fields
    instance.period_menu['checkboxes']['labels'] = []
    instance.period_menu['checkboxes']['keep_selected'] = []
    instance.period_menu['checkboxes']['remove_selected'] = []


def init_metadata(instance):
    """ Initialise internal structure to store metadata.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'metadata_menu'):
        instance.metadata_types = {'STATION POSITION': 'Filter stations by measurement position',
                                   'STATION CLASSIFICATIONS': 'Filter stations by station provided classifications',
                                   'STATION MISCELLANEOUS': 'Filter stations by miscellaneous station provided metadata',
                                   'GLOBALLY GRIDDED CLASSIFICATIONS': 'Filter stations by globally gridded classifications',
                                   'MEASUREMENT PROCESS INFORMATION': 'Filter stations by measurement process information'}
            
        instance.metadata_menu = {'window_title': 'METADATA', 
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


def update_representativity_fields(instance):
    """ Update the data representativity menu upon read.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """
 
    # get previously set rangebox labels / values
    previous_mapped_labels = copy.deepcopy(instance.representativity_menu['rangeboxes']['map_vars'])
    previous_lower = copy.deepcopy(instance.representativity_menu['rangeboxes']['current_lower'])

    # hourly temporal resolution?
    if (instance.resolution == 'hourly') or (instance.resolution == 'hourly_instantaneous'):
        # GHOST 
        if instance.reading_ghost:

            instance.representativity_menu['rangeboxes']['map_vars'] = ['hourly_native_representativity_percent',
                                                                        'daily_native_representativity_percent',
                                                                        'monthly_native_representativity_percent',
                                                                        'daily_representativity_percent',
                                                                        'monthly_representativity_percent',
                                                                        'all_representativity_percent']
            instance.representativity_menu['rangeboxes']['labels'] = ['Hourly',
                                                                      'Daily',
                                                                      'Monthly',
                                                                      'Daily',
                                                                      'Monthly',
                                                                      'All']

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Native', 'Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0, 3]

        # non-GHOST
        else:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['daily_representativity_percent',
                                                                        'monthly_representativity_percent',
                                                                        'all_representativity_percent']
            instance.representativity_menu['rangeboxes']['labels'] = ['Daily',
                                                                      'Monthly',
                                                                      'All']

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0]                                                                   

    # daily temporal resolution?
    elif (instance.resolution == 'daily') or (instance.resolution == '3hourly') or \
            (instance.resolution == '6hourly') or (instance.resolution == '3hourly_instantaneous') or \
            (instance.resolution == '6hourly_instantaneous'):
        # GHOST 
        if instance.reading_ghost:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['daily_native_representativity_percent',
                                                                        'monthly_native_representativity_percent',
                                                                        'monthly_representativity_percent',
                                                                        'all_representativity_percent']
            instance.representativity_menu['rangeboxes']['labels'] = ['Daily',
                                                                      'Monthly',
                                                                      'Monthly',
                                                                      'All']

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Native', 'Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0, 2]

        # non-GHOST
        else:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['monthly_representativity_percent',
                                                                        'all_representativity_percent']
            instance.representativity_menu['rangeboxes']['labels'] = ['Monthly',
                                                                      'All']

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0] 

    # monthly temporal resolution?
    elif instance.resolution == 'monthly':
        # GHOST 
        if instance.reading_ghost:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['monthly_native_representativity_percent',
                                                                        'all_representativity_percent']
            instance.representativity_menu['rangeboxes']['labels'] = ['Monthly',
                                                                      'All']
            
            instance.representativity_menu['rangeboxes']['subtitles'] = ['Native', 'Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0, 1]  

        # non-GHOST
        else:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['all_representativity_percent']
            instance.representativity_menu['rangeboxes']['labels'] = ['All']

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0] 
            

    # initialise rangebox values --> for data representativity fields
    # the default is 0%, for max gap fields % the default is 100%
    instance.representativity_menu['rangeboxes']['current_lower'] = []
    for label_ii, label_mapped in enumerate(instance.representativity_menu['rangeboxes']['map_vars']):
        if 'max_gap' in label_mapped:
            instance.representativity_menu['rangeboxes']['current_lower'].append('100')
        else:
            instance.representativity_menu['rangeboxes']['current_lower'].append('0')

        # label previously existed?
        if label_mapped in previous_mapped_labels:
            instance.representativity_menu['rangeboxes']['current_lower'][label_ii] = \
                previous_lower[previous_mapped_labels.index(label_mapped)]


def update_period_fields(instance):
    """ Update the data period menu upon read.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # hourly temporal resolution?
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


def update_metadata_fields(instance):
    """ Update the metadata menu object with metadata associated with newly read data
        for non-numeric metadata gets all the unique fields per metadata variable,
        and sets the available fields as such, and for numeric gets the minimum and maximum
        boundaries of each metadata variable.
        If previously metadata settings for a field deviate from the default,
        then if the same field still exists then the settings (i.e. bounds or
        checkbox selection) are copied across, rather than setting to the default.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # before doing anything check if have all current metadata variables in menu
    # if not, this is either because it is initialised empty, 
    # or because of changing between GHOST and non-GHOST data
    # if there is a difference, fill metadata menu with empty values
    reset_meta = False
    for metadata_type_ii, metadata_type in enumerate(instance.metadata_menu['navigation_buttons']['labels']):
        required_vars = [metadata_var for metadata_var in instance.metadata_vars_to_read
            if (instance.standard_metadata[metadata_var]['metadata_type'] == metadata_type)
                 & (instance.standard_metadata[metadata_var]['data_type'] != object)]
        if len(required_vars) != len(instance.metadata_menu[metadata_type]['rangeboxes']['labels']):
            reset_meta = True
            break

    # reinitialise metadata menu
    if reset_meta:
        init_metadata(instance)

    # update metadata menu
    for meta_var in instance.metadata_vars_to_read:
        
        # get all metadata values for field across all networks and species
        meta_var_field = []

        # transform to list if there is only one species / network
        if isinstance(instance.network, str):
            instance.network = [instance.network]
        if isinstance(instance.species, str):
            instance.species = [instance.species]

        for networkspeci in instance.networkspecies:
            meta_var_field.extend(instance.metadata_in_memory[networkspeci][meta_var].flatten())
        meta_var_field = np.array(meta_var_field)
        meta_var_field[np.where(meta_var_field == 'nan')[0]] = np.NaN

        # get metadata variable type/data type
        metadata_type = instance.standard_metadata[meta_var]['metadata_type']
        metadata_data_type = instance.standard_metadata[meta_var]['data_type']

        # remove NaNs from metadata values
        meta_var_field_nan_removed = meta_var_field[~pd.isnull(meta_var_field)]

        # update pop-up metadata menu object with read metadata values
        # for non-numeric metadata gets all the unique fields per metadata variable
        # and sets the available fields as such
        if metadata_data_type == object:
            # get previous fields in menu 
            previous_fields = copy.deepcopy(instance.metadata_menu[metadata_type][meta_var]['checkboxes']['labels'])
            previous_keep = copy.deepcopy(instance.metadata_menu[metadata_type][meta_var]['checkboxes']['keep_selected'])
            previous_remove = copy.deepcopy(instance.metadata_menu[metadata_type][meta_var]['checkboxes']['remove_selected'])

            # update new labels
            instance.metadata_menu[metadata_type][meta_var]['checkboxes']['labels'] = np.unique(meta_var_field_nan_removed)
            
            # if field previously existed, then copy across checkbox settings for field
            # else set initial checkboxes to be all blank
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
            previous_lower_default = copy.deepcopy(instance.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index])
            previous_upper_default = copy.deepcopy(instance.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index])
            previous_lower = copy.deepcopy(instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index])
            previous_upper = copy.deepcopy(instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index])

            # have some numeric values for existing metadata variable?
            if len(meta_var_field_nan_removed) > 0:
                min_val = str(np.min(meta_var_field_nan_removed))
                max_val = str(np.max(meta_var_field_nan_removed))

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
            # do not have update_representativity_fields some numeric values for metadata variable so set as 'nan'
            else:
                instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = 'nan'
                instance.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index] = 'nan'
                instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = 'nan'
                instance.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index] = 'nan'
                if meta_var in instance.metadata_menu[metadata_type]['rangeboxes']['apply_selected']:
                    instance.metadata_menu[metadata_type]['rangeboxes']['apply_selected'].remove(meta_var)


def multispecies_conf(instance):
    """ Function used when loading from a configuration file. 
        Sets defined multispecies filtering variables, rest of variables are set to default. 

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
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


def representativity_conf(instance):
    """ Function used when loading from a configuration file. 
        Sets defined representativity filter variables, rest of variables are set to default. 

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    for i, label in enumerate(instance.representativity_menu['rangeboxes']['map_vars']):
        if hasattr(instance, label):
            instance.representativity_menu['rangeboxes']['current_lower'][i] = str(getattr(instance, label))


def period_conf(instance):
    """ Function used when loading from a configuration file. 
        Sets defined period filter variables, rest of variables are set to default. 

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    from .configuration import split_options

    if hasattr(instance, 'period'):
        keeps, removes = split_options(instance, instance.period)
        instance.period_menu['checkboxes']['keep_selected'] = keeps
        instance.period_menu['checkboxes']['remove_selected'] = removes


def metadata_conf(instance):
    """ Function used when loading from a configuration file. 
        Sets defined metadata filter variables, rest of variables are set to default.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
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
