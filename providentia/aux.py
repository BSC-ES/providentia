""" Contains functions that are shared by the dashboard and the offline version. 
    Currently, it includes function related to the initialization and update
    of metadata, checks fields coming from conf files etc.
"""

from glob import glob
import os
import copy
import json
import itertools
import datetime
import sys
from netCDF4 import Dataset, chartostring
import numpy as np
import pandas as pd
import math
import pyproj
from scipy.spatial import cKDTree
import seaborn as sns


def multispecies_mapping(species):
    """ Map species special case str to multiple species names. """

    multi_species_map = {'vconcaerobin*':['vconcaerobin1','vconcaerobin2','vconcaerobin3','vconcaerobin4',
                         'vconcaerobin5','vconcaerobin6','vconcaerobin7','vconcaerobin8','vconcaerobin9',
                         'vconcaerobin10','vconcaerobin11','vconcaerobin12','vconcaerobin13','vconcaerobin14',
                         'vconcaerobin15','vconcaerobin16','vconcaerobin17','vconcaerobin18','vconcaerobin19',
                         'vconcaerobin20','vconcaerobin21','vconcaerobin22']}

    return multi_species_map[species]

def get_multispecies_aliases(networkspecies):
    """ Map networkspecies to networkspecies aliases.
        Also get label for alias.   
    """

    multispecies_aliases = {'vconcaerobin1': '0.05',
                            'vconcaerobin2': '0.066',
                            'vconcaerobin3': '0.086',
                            'vconcaerobin4': '0.113',
                            'vconcaerobin5': '0.148',
                            'vconcaerobin6': '0.194',
                            'vconcaerobin7': '0.255',
                            'vconcaerobin8': '0.335',
                            'vconcaerobin9': '0.439',
                            'vconcaerobin10': '0.576',
                            'vconcaerobin11': '0.756',
                            'vconcaerobin12': '0.992',
                            'vconcaerobin13': '1.302',
                            'vconcaerobin14': '1.708',
                            'vconcaerobin15': '2.241',
                            'vconcaerobin16': '2.940',
                            'vconcaerobin17': '3.857',
                            'vconcaerobin18': '5.061',
                            'vconcaerobin19': '6.641',
                            'vconcaerobin20': '8.713',
                            'vconcaerobin21': '11.432',
                            'vconcaerobin22': '15.00'
                            }

    multispecies_labels =  {'vconcaerobin1': 'Radius [µm]',
                            'vconcaerobin2': 'Radius [µm]',
                            'vconcaerobin3': 'Radius [µm]',
                            'vconcaerobin4': 'Radius [µm]',
                            'vconcaerobin5': 'Radius [µm]',
                            'vconcaerobin6': 'Radius [µm]',
                            'vconcaerobin7': 'Radius [µm]',
                            'vconcaerobin8': 'Radius [µm]',
                            'vconcaerobin9': 'Radius [µm]',
                            'vconcaerobin10': 'Radius [µm]',
                            'vconcaerobin11': 'Radius [µm]',
                            'vconcaerobin12': 'Radius [µm]',
                            'vconcaerobin13': 'Radius [µm]',
                            'vconcaerobin14': 'Radius [µm]',
                            'vconcaerobin15': 'Radius [µm]',
                            'vconcaerobin16': 'Radius [µm]',
                            'vconcaerobin17': 'Radius [µm]',
                            'vconcaerobin18': 'Radius [µm]',
                            'vconcaerobin19': 'Radius [µm]',
                            'vconcaerobin20': 'Radius [µm]',
                            'vconcaerobin21': 'Radius [µm]',
                            'vconcaerobin22': 'Radius [µm]'
                            }
    
    networkspecies_aliases = [multispecies_aliases[networkspeci] 
                              if networkspeci in multispecies_aliases else networkspeci 
                              for networkspeci in networkspecies]

    labels = np.unique([multispecies_labels[networkspeci] 
                        for networkspeci in networkspecies if networkspeci in multispecies_labels])
    if len(labels) == 1:
        unique_label = labels[0]
    else:
        unique_label = ''

    return networkspecies_aliases, unique_label

def exceedance_lim(species):
    """ Return the exceedance limit depending on the species input. 
        If species doesn't have a reported limit, returns np.NaN.

        :param species: name of species currently selected (e.g. sconco3)
        :type species: str
        :return: value of exceedance limit
        :rtype: int
    """

    exceedance_limits = {'sconco3': 90.21, 'sconcno2': 106.38}
    if species in exceedance_limits:
        return exceedance_limits[species]
    else:
        return np.NaN

def temporal_resolution_order_dict():
    """ Return temporal resolution order for menus as a dictionary.

        :return: temporal resolution order
        :rtype: dict
    """

    resolution_order_dict = {'hourly': 1, '3hourly': 2, '6hourly': 3, 'hourly_instantaneous': 4,
                             '3hourly_instantaneous': 5, '6hourly_instantaneous': 6,
                             'daily': 7, 'monthly': 8}

    return resolution_order_dict

def temp_axis_dict():
    """ Return temporal mapping as a dictionary used for the plots.

        :return: numbering of months/days
        :rtype: dict
    """

    map_dict = {'short': {'dayofweek': {0: 'M', 1: 'T', 2: 'W', 3: 'T', 4: 'F', 5: 'S', 6: 'S'},
                          'month': {1: 'J', 2: 'F', 3: 'M', 4: 'A', 5: 'M', 6: 'J',
                                    7: 'J', 8: 'A', 9: 'S', 10: 'O', 11: 'N', 12: 'D'}},
                'long': {'dayofweek': {0: 'Monday', 1: 'Tuesday', 2: 'Wednesday', 3: 'Thursday', 
                                       4: 'Friday', 5: 'Saturday', 6: 'Sunday'},
                         'month': {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June',
                                   7: 'July', 8: 'August', 9: 'September', 10: 'October', 11: 'November', 
                                   12: 'December'}}
                }

    return map_dict

def periodic_xticks():
    """ Return xticks for periodic subplots.

        :return dictionary of xticks per temporal resolution
        :rtype dict
    """

    return {'hour':np.arange(24, dtype=np.int), 
            'dayofweek':np.arange(7, dtype=np.int), 
            'month':np.arange(1, 13, dtype=np.int)}

def periodic_labels():
    """ Return axes labels for periodic subplots.

        :return: axes labels
        :rtype: dict
    """

    return {'hour':'H', 'dayofweek':'DoW', 'month':'M'}

def get_relevant_temporal_resolutions(resolution):        
    """ Get relevant temporal reolsutions for periodic plots, by selected temporal resolution.

        :param resolution: name of selected temporal resolution
        :type resolution: str
        :return: relevant temporal resolutions
        :rtype: list
    """

    if 'hourly' in resolution:
        relevant_temporal_resolutions = ['hour', 'dayofweek', 'month']
    elif resolution == 'daily':
        relevant_temporal_resolutions = ['dayofweek', 'month']
    elif resolution == 'monthly':
        relevant_temporal_resolutions = ['month']
    else:
        relevant_temporal_resolutions = []
        
    return relevant_temporal_resolutions

def get_land_polygon_resolution(selection):
    """ Get resolution of land polygons to plot on map.

        :param selection: name of selected temporal resolution
        :type resolution: str
        :return: selected land polygon resolution
        :rtype: list
    """
    
    land_polygon_resolutions = {'low': '110m','medium': '50m','high': '10m'}
    
    return land_polygon_resolutions[selection]

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
                                          'page_title': 'Select % Data Representativity Bounds', 
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
            & (instance.standard_metadata[metadata_var]['data_type'] != np.object)]

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
            (instance.standard_metadata[metadata_var]['data_type'] == np.object)]
        
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
                                                                        'hourly_native_max_gap_percent',
                                                                        'daily_native_representativity_percent',
                                                                        'daily_native_max_gap_percent',
                                                                        'monthly_native_representativity_percent',
                                                                        'monthly_native_max_gap_percent',
                                                                        'daily_representativity_percent',
                                                                        'daily_max_gap_percent',
                                                                        'monthly_representativity_percent',
                                                                        'monthly_max_gap_percent',
                                                                        'all_representativity_percent', 
                                                                        'all_max_gap_percent']
            
            #instance.representativity_menu['rangeboxes']['labels'] = ['Min Hourly Rep. %',
            #                                                          'Max Hourly Gap %',
            #                                                          'Min Daily Rep. %',
            #                                                          'Max Daily Gap %',
            #                                                          'Min Monthly Rep. %',
            #                                                          'Max Monthly Gap %',
            #                                                          'Min Daily Rep. %',
            #                                                          'Max Daily Gap %',
            #                                                          'Min Monthly Rep. %',
            #                                                          'Max Monthly Gap %',
            #                                                          'Min All Rep. %',
            #                                                          'Max All Gap %']
            instance.representativity_menu['rangeboxes']['labels'] = copy.deepcopy(instance.representativity_menu['rangeboxes']['map_vars'])

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Native', 'Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0, 6]

        # non-GHOST
        else:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['daily_representativity_percent',
                                                                        'daily_max_gap_percent',
                                                                        'monthly_representativity_percent',
                                                                        'monthly_max_gap_percent',
                                                                        'all_representativity_percent', 
                                                                        'all_max_gap_percent']

            #instance.representativity_menu['rangeboxes']['labels'] = ['Min Daily Rep. %',
            #                                                          'Max Daily Gap %',
            #                                                          'Min Monthly Rep. %',
            #                                                          'Max Monthly Gap %',
            #                                                          'Min All Rep. %',
            #                                                          'Max All Gap %']
            instance.representativity_menu['rangeboxes']['labels'] = copy.deepcopy(instance.representativity_menu['rangeboxes']['map_vars'])

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0]                                                                   

    # daily temporal resolution?
    elif (instance.resolution == 'daily') or (instance.resolution == '3hourly') or \
            (instance.resolution == '6hourly') or (instance.resolution == '3hourly_instantaneous') or \
            (instance.resolution == '6hourly_instantaneous'):
        # GHOST 
        if instance.reading_ghost:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['daily_native_representativity_percent',
                                                                        'daily_native_max_gap_percent',
                                                                        'monthly_native_representativity_percent',
                                                                        'monthly_native_max_gap_percent',
                                                                        'monthly_representativity_percent',
                                                                        'monthly_max_gap_percent',
                                                                        'all_representativity_percent', 
                                                                        'all_max_gap_percent']

            #instance.representativity_menu['rangeboxes']['labels'] = ['Min Daily Rep. %',
            #                                                          'Max Daily Gap %',
            #                                                          'Min Monthly Rep. %',
            #                                                          'Max Monthly Gap %',
            #                                                          'Min Monthly Rep. %',
            #                                                          'Max Monthly Gap %',
            #                                                          'Min All Rep. %',
            #                                                          'Max All Gap %']
            instance.representativity_menu['rangeboxes']['labels'] = copy.deepcopy(instance.representativity_menu['rangeboxes']['map_vars'])

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Native', 'Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0, 4]

        # non-GHOST
        else:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['monthly_representativity_percent',
                                                                        'monthly_max_gap_percent',
                                                                        'all_representativity_percent', 
                                                                        'all_max_gap_percent']

            #instance.representativity_menu['rangeboxes']['labels'] = ['Min Monthly Rep. %',
            #                                                          'Min All Rep. %',
            #                                                          'Max Monthly Gap %',
            #                                                          'Max All Gap %'] 
            instance.representativity_menu['rangeboxes']['labels'] = copy.deepcopy(instance.representativity_menu['rangeboxes']['map_vars'])

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0] 

    # monthly temporal resolution?
    elif instance.resolution == 'monthly':
        # GHOST 
        if instance.reading_ghost:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['monthly_native_representativity_percent',
                                                                        'monthly_native_max_gap_percent',
                                                                        'all_representativity_percent', 
                                                                        'all_max_gap_percent']

            #instance.representativity_menu['rangeboxes']['labels'] = ['Min Monthly Rep. %',
            #                                                          'Max Monthly Gap %',
            #                                                          'Min All Rep. %',
            #                                                          'Max All Gap %'] 
            instance.representativity_menu['rangeboxes']['labels'] = copy.deepcopy(instance.representativity_menu['rangeboxes']['map_vars'])

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Native', 'Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0, 2]  

        # non-GHOST
        else:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['all_representativity_percent', 
                                                                        'all_max_gap_percent']
            
            #instance.representativity_menu['rangeboxes']['labels'] = ['Min All Rep. %',
            #                                                          'Max All Gap %'] 
            instance.representativity_menu['rangeboxes']['labels'] = copy.deepcopy(instance.representativity_menu['rangeboxes']['map_vars'])

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
                 & (instance.standard_metadata[metadata_var]['data_type'] != np.object)]
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
        if metadata_data_type == np.object:
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

def valid_date(date_text):
    """ Determine if a date string is in the correct format. """

    try:
        datetime.datetime.strptime(str(date_text), '%Y%m%d')
        return True
    except Exception as e:
        return False

def check_for_ghost(network_name):
    """ Check whether the selected network comes from GHOST or not.
        All non-GHOST networks start with an asterisk at their name.
    """

    if '/' in network_name:
        return False
    else:
        return True

def load_conf(self, fpath=None):
    """ Load existing configurations from file
        for running offline Providentia.
    """

    from .configuration import read_conf

    if fpath is None:
        print("No configuration file found")
        sys.exit(1)

    # if DEFAULT is not present, then return
    if not os.path.isfile(fpath):
        print(("Error %s" % fpath))
        return

    self.sub_opts, self.all_sections, self.parent_section_names, self.subsection_names, self.filenames = read_conf(fpath)

def update_plotting_parameters(instance):
    """ Function that updates plotting parameters (colour and zorder) for data labels.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # generate a list of RGB tuples for number of experiments there are
    sns.reset_orig()
    clrs = sns.color_palette(instance.plot_characteristics_templates['general']['legend_color_palette'], n_colors=len(instance.data_labels)-1)

    # add colour and zorder for observations
    instance.plotting_params['observations']['colour'] = instance.plot_characteristics_templates['general']['obs_markerfacecolor']
    instance.plotting_params['observations']['zorder'] = instance.plot_characteristics_templates['general']['obs_zorder']

    # add colours and zorder for each experiment
    experiment_ind = 1
    for data_label in instance.data_labels:
        if data_label != 'observations':
            # define colour for experiment
            instance.plotting_params[data_label]['colour'] = clrs[experiment_ind-1]
            # define zorder for experiment (obs zorder + experiment_ind)
            instance.plotting_params[data_label]['zorder'] = \
                instance.plotting_params['observations']['zorder'] + experiment_ind
            # update count of experiments
            experiment_ind += 1

def get_ghost_observational_tree(instance):
    """ Create GHOST observational data tree as a nested dictionary,
        storing a list of start YYYYMM yearmonths per:
        network / resolution / matrix / speci

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :return: GHOST observational tree dictionary
        :rtype: dict
    """

    ghost_observation_data = {}

    # iterate through available networks
    for network in instance.available_networks:

        # check if directory for network exists
        # if not, continue
        if not os.path.exists('%s/%s/%s' % (instance.ghost_root, network, instance.ghost_version)):
            continue

        # write empty dictionary for network
        ghost_observation_data[network] = {}

        # iterate through available resolutions
        for resolution in instance.available_resolutions:

            # check if directory for resolution exists
            # if not, continue
            if not os.path.exists('%s/%s/%s/%s' % (instance.ghost_root, network, instance.ghost_version, resolution)):
                continue

            # write nested empty dictionary for resolution
            ghost_observation_data[network][resolution] = {}

            # get available species for network/resolution
            available_species = os.listdir('%s/%s/%s/%s' % (instance.ghost_root, network, instance.ghost_version, resolution))

            # iterate through available files per species
            for speci in available_species:

                # get all available netCDF files
                available_files = os.listdir(
                    '%s/%s/%s/%s/%s' % (instance.ghost_root, network, instance.ghost_version, resolution, speci))
                
                # coninue if have no files
                if len(available_files) == 0:
                    continue

                # get monthly start date (YYYYMM) of all files
                file_yearmonths = sorted([f.split('_')[-1][:6] for f in available_files if f != 'temporary'])
                
                # get matrix for current species
                matrix = instance.parameter_dictionary[speci]['matrix']
                if matrix not in ghost_observation_data[network][resolution]:
                    # write nested empty dictionary for matrix
                    ghost_observation_data[network][resolution][matrix] = {}

                # write nested dictionary for species, with associated file yearmonths
                ghost_observation_data[network][resolution][matrix][speci] = file_yearmonths

    return ghost_observation_data

def get_nonghost_observational_tree(instance, nonghost_observation_data_json):
    """ Fill non-GHOST observational data tree,
        storing a list of start YYYYMM yearmonths per:
        network / resolution / matrix / speci

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :param nonghost_observation_data_json: non-GHOST observational tree json
        :type nonghost_observation_data_json: json
        :return: non-GHOST observational tree dictionary
        :rtype: dict
    """

    nonghost_observation_data = {}

    # iterate through networks
    for network in nonghost_observation_data_json:

        # check if directory for network exists
        # if not, continue
        if not os.path.exists('%s/%s' % (instance.nonghost_root, network)):
            continue

        # write empty dictionary for network
        nonghost_observation_data[network] = {}

        # iterate through resolutions
        for resolution in nonghost_observation_data_json[network]:

            # check if directory for resolution exists
            # if not, continue
            if not os.path.exists('%s/%s/%s' % (instance.nonghost_root, network, resolution)):
                continue

            # write nested empty dictionary for resolution
            nonghost_observation_data[network][resolution] = {}

            # iterate through species
            for speci in nonghost_observation_data_json[network][resolution]:

                # get all available netCDF files 
                available_files = glob('%s/%s/%s/%s/%s_??????.nc' % (instance.nonghost_root, network, resolution, speci, speci))

                # coninue if have no files
                if len(available_files) == 0:
                    continue

                # get monthly start date (YYYYMM) of all files
                file_yearmonths = sorted([f.split('_')[-1][:6] for f in available_files])

                # get matrix for current species
                matrix = instance.parameter_dictionary[speci]['matrix']
                if matrix not in nonghost_observation_data[network][resolution]:
                    # write nested empty dictionary for matrix
                    nonghost_observation_data[network][resolution][matrix] = {}

                # write nested dictionary for species, with associated file yearmonths
                nonghost_observation_data[network][resolution][matrix][speci] = file_yearmonths
    
    return nonghost_observation_data

def get_experiment_tree(instance):
    """ Fill experiment data tree,
        storing a list of start YYYYMM yearmonths per:
        network / resolution / matrix / speci

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :return: experiment tree dictionary
        :rtype: dict
    """

    experiment_data = {}

    # get all different experiment names (from providentia-interpolation output dir)
    available_experiments = os.listdir('%s/%s' % (instance.exp_root, instance.ghost_version))

    # iterate through all observation data dict, and attempt to get recriprocal experiment files
    for network in instance.all_observation_data:
        for resolution in instance.all_observation_data[network]:
            for matrix in instance.all_observation_data[network][resolution]:
                for species in instance.all_observation_data[network][resolution][matrix]:
                    # iterate through available experiments
                    for experiment in available_experiments:
                        
                        # get folder where interpolated experiments are saved
                        if '/' not in network:
                            files_directory = '%s/%s/%s/%s/%s/%s' % (instance.exp_root, instance.ghost_version, 
                                                                     experiment, resolution, species, network)
                        else:
                            files_directory = '%s/%s/%s/%s/%s/%s' % (instance.exp_root, instance.ghost_version, 
                                                                      experiment, resolution, species,
                                                                      network.replace('/', '-'))
                        
                        # test if interpolated directory exists for experiment
                        # if it does not exit, continue
                        if not os.path.exists(files_directory):
                            continue
                        else:
                            # get all available netCDF files
                            available_files = os.listdir(files_directory)

                            # get monthly start date (YYYYMM) of all files
                            file_yearmonths = sorted([f.split('_')[-1][:6] for f in available_files])

                            # write nested dictionary for experiment, with associated file yearmonths
                            if len(file_yearmonths) > 0:
                                # if network/resolution/matrix/species/experiment not in dictionary yet, add it
                                if network not in experiment_data:
                                    experiment_data[network] = {}
                                if resolution not in experiment_data[network]:
                                    experiment_data[network][resolution] = {}
                                if matrix not in experiment_data[network][resolution]:
                                    experiment_data[network][resolution][matrix] = {}
                                if species not in experiment_data[network][resolution][matrix]:
                                    experiment_data[network][resolution][matrix][species] = {}
                                experiment_data[network][resolution][matrix][species][experiment] = file_yearmonths

    return experiment_data

def get_valid_obs_files_in_date_range(instance, start_date, end_date):
    """ Iterate through observational dictionary tree and return 
        a dictionary of available data in the selected daterange

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :param start_date: start date (e.g. "20201101")
        :type start_date: str
        :param end_date: end date (e.g. "20201101")
        :type end_date: str
        :return: available observational tree dictionary
        :rtype: dict
    """

    # create dictionary to store available observational data
    instance.available_observation_data = {}

    # check if start/end date are valid values, if not, return with no valid obs. files
    if (not valid_date(start_date)) or (not valid_date(end_date)):
        print(f'Warning: One of start date ({start_date}) or end date ({end_date}) are not valid.')
        return

    # check end date is > start date, if not, return with no valid obs. files
    if int(start_date) >= int(end_date):
        print(f'Warning: Start date ({start_date}) exceeds end date ({end_date}).')
        return

    # check start date and end date are both within if valid date range (19000101 - 20500101),
    # if not, return with no valid obs. files
    if (int(start_date) < 19000101) or (int(end_date) < 19000101) or (int(start_date) >= 20500101) or (int(end_date) >= 20500101):
        print(f'Warning: One of start date ({start_date}) or end date ({end_date}) are not valid.')
        return 

    # get start date on first of month
    start_date_firstdayofmonth = int(str(start_date)[:6] + '01')

    # iterate through observational dictionary
    for network in instance.all_observation_data:
        for resolution in instance.all_observation_data[network]:
            for matrix in instance.all_observation_data[network][resolution]:
                for speci in instance.all_observation_data[network][resolution][matrix]:
                    
                    # get file yearmonths
                    file_yearmonths = instance.all_observation_data[network][resolution][matrix][speci]

                    # get file yearmonths within date range
                    valid_file_yearmonths = sorted([ym for ym in file_yearmonths if
                                                    (int('{}01'.format(ym)) >= start_date_firstdayofmonth) & (int('{}01'.format(ym)) < int(end_date))])
                    
                    # add yearmonths to available observation data dict 
                    if len(valid_file_yearmonths) > 0:
                        if network not in instance.available_observation_data:
                            instance.available_observation_data[network] = {}
                        if resolution not in instance.available_observation_data[network]:
                            instance.available_observation_data[network][resolution] = {}
                        if matrix not in instance.available_observation_data[network][resolution]:
                            instance.available_observation_data[network][resolution][matrix] = {}
                        instance.available_observation_data[network][resolution][matrix][speci] = valid_file_yearmonths

def get_valid_experiments(instance, start_date, end_date, resolution, networks, species):
    """ Get valid experiments for daterange, and selected parameters.
        Update experiment pop-up with valid experiments.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :param start_date: start date (e.g. "20201101")
        :type start_date: str
        :param end_date: end date (e.g. "20201231")
        :type end_date: str
        :param resolution: resolution (e.g. "hourly")
        :type resolution: str
        :param networks: list of networks
        :type networks: list
        :param species: list of species
        :type species: list 
    """

    # get all different experiment names (from providentia-interpolation output dir)
    available_experiments = os.listdir('%s/%s' % (instance.exp_root, instance.ghost_version))

    # create dictionary to store available experiment data
    instance.available_experiment_data = {}

    #list for saving experiments to add to experiments pop-up 
    experiments_to_add = []

    # get start date on first of month
    start_date_firstdayofmonth = int(str(start_date)[:6] + '01')

    # iterate through networks and species
    for network, speci in zip(networks, species):

        # iterate through available experiments
        for experiment in available_experiments:

            # get folder where interpolated experiments are saved
            if '/' not in network:           
                files_directory = '%s/%s/%s/%s/%s/%s' % (instance.exp_root, instance.ghost_version, 
                                                         experiment, resolution, speci, network)
            else:
                files_directory = '%s/%s/%s/%s/%s/%s' % (instance.exp_root, instance.ghost_version, 
                                                          experiment, resolution, speci,
                                                          network.replace('/', '-'))
                
            # test if interpolated directory exists for experiment
            # if it does not exit, continue
            if not os.path.exists(files_directory):
                continue
            else:
                # get all available netCDF files
                available_files = os.listdir(files_directory)

            # get monthly start date (YYYYMM) of all files
            file_yearmonths = sorted([f.split('_')[-1][:6] for f in available_files])

            # write nested dictionary for experiment, with associated file yearmonths
            if len(file_yearmonths) > 0:

                # get file yearmonths within date range
                valid_file_yearmonths = sorted([ym for ym in file_yearmonths if 
                                                (int('{}01'.format(ym)) >= start_date_firstdayofmonth) & (int('{}01'.format(ym)) < int(end_date))])

                #if have valid files, then add experiment to pop-up menu, 
                # and add yearmonths to available experiment data
                if len(valid_file_yearmonths) > 0:
                    experiments_to_add.append(experiment)

                    if network not in instance.available_experiment_data:
                        instance.available_experiment_data[network] = {}
                    if resolution not in instance.available_experiment_data[network]:
                        instance.available_experiment_data[network][resolution] = {}
                    if speci not in instance.available_experiment_data[network][resolution]:
                        instance.available_experiment_data[network][resolution][speci] = {}
                    if experiment not in instance.available_experiment_data[network][resolution][speci]:
                        instance.available_experiment_data[network][resolution][speci][experiment] = valid_file_yearmonths

    # set list of experiment names to add on experiments pop-up
    if not instance.offline:
        experiments_to_add = np.array(sorted(experiments_to_add))
        instance.experiments_menu['checkboxes']['labels'] = experiments_to_add
        instance.experiments_menu['checkboxes']['map_vars'] = experiments_to_add

def get_basic_metadata(instance):     
    """ Get basic unique metadata across networkspecies wanting to read
        The basic fields are: station_reference, station_name, longitude, latitude, measurement_altitude, station_classification and area_classification

        If have multiple species, then spatially cocolocate across species 
        to get matching stations across stations.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :return: station_references per networkspecies, station_name per networkspecies, longitudes per networkspecies, latitudes per networkspecies, measurement altitudes per networkspecies, station_classifications per networkspecies, area_classifications per networkspecies, nonghost_units 
        :rtype: dict, dict, dict, dict, dict, dict, dict, dict
    """

    # define dictionaries for storing basic metadata across all species to read
    station_references = {}
    station_names = {}
    station_longitudes = {}
    station_latitudes = {}
    station_measurement_altitudes = {}
    station_classifications = {}
    area_classifications = {}
    nonghost_units = {}

    # iterate through network, speci pairs
    for networkspeci in (instance.networkspecies + instance.filter_networkspecies):
    
        # get indivudual network and species strings
        network = networkspeci.split('|')[0]
        speci = networkspeci.split('|')[1]

        # get species matrix
        matrix = instance.parameter_dictionary[speci]['matrix']

        # get file root
        # GHOST
        if instance.reading_ghost:
            file_root = '%s/%s/%s/%s/%s/%s_' % (instance.ghost_root, network,
                                                instance.ghost_version, instance.resolution,
                                                speci, speci)
        # non-GHOST
        else:
            file_root = '%s/%s/%s/%s/%s_' % (instance.nonghost_root, network,
                                             instance.resolution, speci, speci)

        # get relevant files
        relevant_files_before_filter = sorted([file_root+str(yyyymm)+'.nc' for yyyymm in instance.yearmonths])
        relevant_files = copy.deepcopy(relevant_files_before_filter)

        # drop files if they don't exist
        for file in relevant_files_before_filter:
            if not os.path.exists(file):
                relevant_files.remove(file)

        # if have 0 files to read for networkspeci, then drop networkspeci
        if len(relevant_files) == 0:
            if networkspeci in instance.networkspecies:
                instance.networkspecies.remove(networkspeci)
                print('Warning: There is no available observational data for the network|species: {}. Dropping.'.format(networkspeci))
            elif networkspeci in instance.filter_networkspecies:
                instance.filter_networkspecies.remove(networkspeci)
                del instance.filter_species[networkspeci]
                print('Warning: There is no available observational data for the filter network|species: {}. Dropping.'.format(networkspeci))
            continue
            
        # get station references, longitudes and latitudes for speci
        # GHOST
        if instance.reading_ghost:
            
            # define arrays for storing speci metadata
            speci_station_references = []
            speci_station_names = []
            speci_station_longitudes = []
            speci_station_latitudes = []
            speci_station_measurement_altitudes = []
            speci_station_classifications = []
            speci_area_classifications = []

            for relevant_file in relevant_files:
                ncdf_root = Dataset(relevant_file)
                speci_station_references = np.append(speci_station_references, ncdf_root['station_reference'][:])
                speci_station_names = np.append(speci_station_names, ncdf_root['station_name'][:])
                speci_station_longitudes = np.append(speci_station_longitudes, ncdf_root['longitude'][:])
                speci_station_latitudes = np.append(speci_station_latitudes, ncdf_root['latitude'][:])
                speci_station_measurement_altitudes = np.append(speci_station_measurement_altitudes, ncdf_root['measurement_altitude'][:])
                speci_station_classifications = np.append(speci_station_classifications, ncdf_root['station_classification'][:])
                speci_area_classifications = np.append(speci_area_classifications, ncdf_root['area_classification'][:])
                ncdf_root.close()

            speci_station_references, station_unique_indices = np.unique(speci_station_references, return_index=True)
            station_references[networkspeci] = speci_station_references
            station_names[networkspeci] = speci_station_names[station_unique_indices]
            station_longitudes[networkspeci] = speci_station_longitudes[station_unique_indices]
            station_latitudes[networkspeci] = speci_station_latitudes[station_unique_indices]
            station_measurement_altitudes[networkspeci] = speci_station_measurement_altitudes[station_unique_indices]
            station_classifications[networkspeci] = speci_station_classifications[station_unique_indices]
            area_classifications[networkspeci] = speci_area_classifications[station_unique_indices]
        
        # non-GHOST
        else:
        
            ncdf_root = Dataset(relevant_files[0])
            if 'station_reference' in ncdf_root.variables:
                station_reference_var = 'station_reference'
            elif 'station_code' in ncdf_root.variables:
                station_reference_var = 'station_code'
            elif 'station_name' in ncdf_root.variables:
                station_reference_var = 'station_name'
            else: 
                print('Error: {} cannot be read because it has no station_name.'.format(relevant_file))
                sys.exit()

            meta_shape = ncdf_root[station_reference_var].shape
            station_references[networkspeci] = ncdf_root[station_reference_var][:]
            meta_val_dtype = np.array([station_references[networkspeci][0]]).dtype
            if len(meta_shape) == 2:
                if meta_val_dtype == np.dtype(object):
                    station_references[networkspeci] = np.array([''.join(val) for val in station_references[networkspeci]])
                else:
                    station_references[networkspeci] = chartostring(station_references[networkspeci])

            # get indices of all non-NaN stations (can be NaN for some non-GHOST files)
            non_nan_station_indices = np.array([ref_ii for ref_ii, ref in enumerate(station_references[networkspeci]) if ref.lower() != 'nan'])
            station_references[networkspeci] = station_references[networkspeci][non_nan_station_indices]

            if "latitude" in ncdf_root.variables:
                station_longitudes[networkspeci] = ncdf_root['longitude'][non_nan_station_indices]
                station_latitudes[networkspeci] = ncdf_root['latitude'][non_nan_station_indices]
            else:
                station_longitudes[networkspeci] = ncdf_root['lon'][non_nan_station_indices]
                station_latitudes[networkspeci] = ncdf_root['lat'][non_nan_station_indices]

            if "station_name" in ncdf_root.variables:
                meta_shape = ncdf_root['station_name'].shape
                station_names[networkspeci] = ncdf_root['station_name'][non_nan_station_indices]
                meta_val_dtype = np.array([station_names[networkspeci][0]]).dtype
                if len(meta_shape) == 2:
                    if meta_val_dtype == np.dtype(object):
                        station_names[networkspeci] = np.array([''.join(val) for val in station_names[networkspeci]])
                    else:
                        station_names[networkspeci] = chartostring(station_names[networkspeci])

            if "station_classification" in ncdf_root.variables:
                meta_shape = ncdf_root['station_classification'].shape
                station_classifications[networkspeci] = ncdf_root['station_classification'][non_nan_station_indices]
                meta_val_dtype = np.array([station_classifications[networkspeci][0]]).dtype
                if len(meta_shape) == 2:
                    if meta_val_dtype == np.dtype(object):
                        station_classifications[networkspeci] = np.array([''.join(val) for val in station_classifications[networkspeci]])
                    else:
                        station_classifications[networkspeci] = chartostring(station_classifications[networkspeci])

            if "area_classification" in ncdf_root.variables:
                meta_shape = ncdf_root['area_classification'].shape
                area_classifications[networkspeci] = ncdf_root['area_classification'][non_nan_station_indices]
                meta_val_dtype = np.array([area_classifications[networkspeci][0]]).dtype
                if len(meta_shape) == 2:
                    if meta_val_dtype == np.dtype(object):
                        area_classifications[networkspeci] = np.array([''.join(val) for val in area_classifications[networkspeci]])
                    else:
                        area_classifications[networkspeci] = chartostring(area_classifications[networkspeci])

            # get non-GHOST measurement units
            nonghost_units[speci] = ncdf_root[speci].units

            ncdf_root.close()

    # if have more than 1 networkspecies (including filter networkspecies), and spatial_colocation is active,
    # then spatially colocate stations across species
    if (len((instance.networkspecies + instance.filter_networkspecies)) > 1) & (instance.spatial_colocation):
        # get intersecting station indices across species (handle both GHOST and non-GHOST cases)
        if instance.reading_ghost:
            intersecting_indices = spatial_colocation_ghost(station_longitudes, station_latitudes, station_measurement_altitudes)
        else:
            intersecting_indices = spatial_colocation_nonghost(station_references, station_longitudes, station_latitudes)
        
        # iterate through networkspecies specific intersecting indices, setting 
        for ns, ns_intersects in intersecting_indices.items():
            station_references[ns] = station_references[ns][ns_intersects]
            station_longitudes[ns] = station_longitudes[ns][ns_intersects]
            station_latitudes[ns] = station_latitudes[ns][ns_intersects]
            if ns in station_measurement_altitudes:
                station_measurement_altitudes[ns] = station_measurement_altitudes[ns][ns_intersects]
            if ns in station_names:
                station_names[ns] = station_names[ns][ns_intersects]
            if ns in station_classifications:
                station_classifications[ns] = station_classifications[ns][ns_intersects]
            if ns in area_classifications:
                area_classifications[ns] =  area_classifications[ns][ns_intersects] 

    return station_references, station_names, station_longitudes, station_latitudes, station_measurement_altitudes, station_classifications, area_classifications, nonghost_units

def spatial_colocation_nonghost(station_references, longitudes, latitudes):
    """ Given multiple species, return intersecting indices for matching stations across species (per network/species)
        for non-GHOST data.

        This is done by 
            1. Cross-checking the station references between species to get matching station_references
            2. Cross-checking matching longitude / latitude coordinates to a tolerance of 19m difference

        The tolerance is calculated by allowing for a tolerance of 19.053m in the 3 independent x,y,z dimensions, 
        as is done in GHOST to distinguish unique stations.
        Using Pythagoras in 3D √(11**2 +11**2 + 11**2) = 19.053.

        :param station_references: dictionary of station references per network/species
        :type station_references: dict
        :param longitudes: dictionary of longitudes per network/species
        :type longitudes: dict
        :param latitudes: dictionary of latitudes per network/species
        :type latitudes: dict
        :return: intersecting indices per network/species
        :rtype: dict
    """

    # get indices of intersection of station references across species
    #for networkspecies in station_references:
    #    station_references[networkspecies] = station_references[networkspecies].tolist()    
    intersecting_station_references = list(set.intersection(*map(set,list(station_references.values()))))
    intersecting_indices = {}
    for networkspecies in station_references:
        intersecting_indices[networkspecies] = np.array([list(station_references[networkspecies]).index(ref) 
                                                        for ref in intersecting_station_references], dtype=np.int)

    # set variable for first networkspecies
    firstnetworkspecies = list(intersecting_indices.keys())[0]

    # if have zero intersecting indices across species, then return with warning message
    if len(intersecting_indices[firstnetworkspecies]) == 0:
        print('Warning: No intersecting stations across networks/species')
        return intersecting_indices

    # if non-intersecting indices unaccounted for across species, 
    # then attempt to resolve them by matching longitudes / latitudes
    if len(intersecting_indices[firstnetworkspecies]) != len(station_references[firstnetworkspecies]):

        # set tolerance for matching longitudes and latitudes in metres
        tolerance = 19.053

        # get non-intersecting indices, longitudes and latitudes across speci
        non_intersecting_indices = {networkspecies: np.setdiff1d(np.arange(len(station_references[networkspecies])), intersecting_indices[networkspecies]) for networkspecies in station_references}
        non_intersecting_longitudes = {networkspecies: longitudes[networkspecies][non_intersecting_indices[networkspecies]] for networkspecies in longitudes}
        non_intersecting_latitudes = {networkspecies: latitudes[networkspecies][non_intersecting_indices[networkspecies]] for networkspecies in latitudes}

        # get non-intersecting station longitudes and latitudes for first speci
        firstnetworkspecies_longitudes = non_intersecting_longitudes[firstnetworkspecies]
        firstnetworkspecies_latitudes = non_intersecting_latitudes[firstnetworkspecies]

        # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
        # (Earth Centred, Earth Fixed) coordinates assuming WGS84 datum and ellipsoid, and that all heights equal zero
        # ECEF coordinates represent positions (in metres) as X, Y, Z coordinates, approximating the earth surface as an ellipsoid of revolution
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        firstnetworkspecies_x, firstnetworkspecies_y, firstnetworkspecies_z = pyproj.transform(lla, ecef, 
            firstnetworkspecies_longitudes, firstnetworkspecies_latitudes, 
            np.zeros(len(firstnetworkspecies_longitudes)), radians=False)
        
        # merge coordinates to 3D array
        firstnetworkspecies_xyz = np.column_stack((firstnetworkspecies_x, firstnetworkspecies_y, firstnetworkspecies_z))

        # iterate through all other speci, and get intersections (within tolerance) of longitudes and latitudes 
        # with first speci longitudes and latitudes
        pairwise_intersect_inds = {firstnetworkspecies:[]}
        for networkspecies in non_intersecting_longitudes:

            if networkspecies == firstnetworkspecies:
                continue

            nextnetworkspecies_longitudes = non_intersecting_longitudes[networkspecies]
            nextnetworkspecies_latitudes = non_intersecting_latitudes[networkspecies]

            # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
            nextnetworkspecies_x, nextnetworkspecies_y, nextnetworkspecies_z = pyproj.transform(lla, ecef, 
                nextnetworkspecies_longitudes, nextnetworkspecies_latitudes, 
                np.zeros(len(nextnetworkspecies_longitudes)), radians=False)
            
            # merge coordinates to 3D array
            nextnetworkspecies_xyz = np.column_stack((nextnetworkspecies_x, nextnetworkspecies_y, nextnetworkspecies_z))            

            # get all indices of next speci xyz coords, within tolerance of each first speci xyz coords
            idx = cKDTree(nextnetworkspecies_xyz).query_ball_point(firstnetworkspecies_xyz, tolerance)

            # get all indices where have non-duplicated and duplicated matched indices
            unique_idx, unique_idx_counts = np.unique(list(itertools.chain(*idx)), return_counts=True)
            nondup_idx = unique_idx[np.where(unique_idx_counts == 1)[0]]
            dup_idx = unique_idx[np.where(unique_idx_counts > 1)[0]]

            # gather list of indices for first speci and next speci of already resolved within tolerance indices
            fs_wtol_inds = []
            ns_wtol_inds = []

            # resolve all duplicated matched indices
            idx, unresolved_dup_idx, fs_wtol_inds, ns_wtol_inds = resolve_duplicate_spatial_colocation_matches(idx, 
                                                                    nondup_idx, dup_idx, 
                                                                    firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                                    fs_wtol_inds, ns_wtol_inds)

            # pass through again to resolve all unresovered duplicated indices after first pass
            if len(unresolved_dup_idx) > 0:
                idx, _, fs_wtol_inds, ns_wtol_inds = resolve_duplicate_spatial_colocation_matches(idx, 
                                                    nondup_idx, unresolved_dup_idx, 
                                                    firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                    fs_wtol_inds, ns_wtol_inds)

            # iterate though next speci within tolerance indices, per each first speci coord 
            for idx_ii, idx_l in enumerate(idx):
                # if position has already been resolved (through resolving duplicates) then continue
                if idx_ii in fs_wtol_inds:
                    continue
                
                # no matches, then append nothing
                elif len(idx_l) == 0:
                    continue

                # just 1 match, then append
                elif len(idx_l) == 1:
                    fs_wtol_inds.append(idx_ii) 
                    ns_wtol_inds.append(idx_l[0])

                # more than 1 match, then find which match has closest distance, then append
                else:
                    # find the dists between all relevant first speci xyz coords and next speci duplicate xyz coord
                    idx_l = np.array(idx_l)
                    dists = cKDTree([firstnetworkspecies_xyz[idx_ii]]).query(nextnetworkspecies_xyz[idx_l], k=1)[0]
                    ordered_idx_l = idx_l[np.argsort(dists)]
                    fs_wtol_inds.append(idx_ii) 
                    ns_wtol_inds.append(ordered_idx_l[0])

            # order matched indices for first species in ascending order, and order next speci indices in smae way
            if len(fs_wtol_inds) > 0:
                fs_wtol_inds, ns_wtol_inds = list(zip(*sorted(zip(fs_wtol_inds, ns_wtol_inds))))

                # set indices where first species differences are within tolerance, i.e. intersecting 
                pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)] = non_intersecting_indices[firstnetworkspecies][np.array(fs_wtol_inds)]
                pairwise_intersect_inds[firstnetworkspecies].extend(non_intersecting_indices[firstnetworkspecies][np.array(fs_wtol_inds)])
            
                # get indices where next species differences are within tolerance, i.e. intersecting 
                pairwise_intersect_inds[networkspecies] = non_intersecting_indices[networkspecies][np.array(ns_wtol_inds)]
            else:
                pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)] = np.array([], dtype=np.int)
                pairwise_intersect_inds[networkspecies] = np.array([], dtype=np.int)

        # get indices (for first networkspecies) where longitude and latitudes intersect across all species
        pairwise_intersect_inds_unique, counts = np.unique(pairwise_intersect_inds[firstnetworkspecies], return_counts=True)
        pairwise_intersect_inds[firstnetworkspecies] = pairwise_intersect_inds_unique[counts == (len(longitudes)-1)]

        if len(pairwise_intersect_inds[firstnetworkspecies]) > 0:
            # get specific intersect indices across all species, for rest of species
            for networkspecies in non_intersecting_longitudes:
                if networkspecies == firstnetworkspecies:
                    continue
                _, species_intersect_inds, _ = np.intersect1d(pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)], pairwise_intersect_inds[firstnetworkspecies], return_indices=True)
                pairwise_intersect_inds[networkspecies] = pairwise_intersect_inds[networkspecies][species_intersect_inds]

            # append newly found intersecting indices to previously found intersect inds
            for networkspecies in non_intersecting_longitudes:
                intersecting_indices[networkspecies] = np.array(np.append(intersecting_indices[networkspecies], pairwise_intersect_inds[networkspecies]), dtype=np.int)

    return intersecting_indices

def spatial_colocation_ghost(longitudes, latitudes, measurement_altitudes):
    """ Given multiple species, return intersecting indices for matching stations across species (per network/species)
        for GHOST data.

        This is done by cross-checking matching longitude / latitude / measurement altitudes coordinates 
        to a tolerance of 19.053m difference.
        This tolerance is calculated by allowing for a tolerance of 11m in the 3 independent x,y,z dimensions, 
        as is done in GHOST to distinguish unique stations.
        Using Pythagoras in 3D √(11**2 +11**2 + 11**2) = 19.053.

        A current limitation is that at one station there can be several measurement methods, which
        are represented as unique stations in GHOST. Currently, if these stations have the same measurement 
        position, simply the first of these stations will be preferentially chosen as a match.
        This could be better done by prioritising first by method, when have multiple matches.

        :param longitudes: dictionary of longitudes per network/species
        :type longitudes: dict
        :param latitudes: dictionary of latitudes per network/species
        :type latitudes: dict
        :param measurement_altitudes: dictionary of measurement altitudes per network/species
        :type measurement_altitudes: dict
        :return: intersecting indices per network/species
        :rtype: dict
    """

    # set tolerance for matching longitudes / latitudes / measurement_altitudes in metres
    tolerance = 19.053

    # set variable for first networkspecies
    firstnetworkspecies = list(longitudes.keys())[0]
    
    # get station coordinates for firstnetworkspecies
    firstnetworkspecies_longitudes = longitudes[firstnetworkspecies]
    firstnetworkspecies_latitudes = latitudes[firstnetworkspecies] 
    firstnetworkspecies_measurement_altitudes = measurement_altitudes[firstnetworkspecies]

    # convert longitudes / latitudes / measurement_altitudes in geographic coordinates to cartesian ECEF 
    # (Earth Centred, Earth Fixed) coordinates assuming WGS84 datum and ellipsoid, and that all heights equal zero
    # ECEF coordinates represent positions (in metres) as X, Y, Z coordinates, approximating the earth surface as an ellipsoid of revolution
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    firstnetworkspecies_x, firstnetworkspecies_y, firstnetworkspecies_z = pyproj.transform(lla, ecef, 
        firstnetworkspecies_longitudes, firstnetworkspecies_latitudes, firstnetworkspecies_measurement_altitudes,
        radians=False)

    # merge coordinates to 3D array
    firstnetworkspecies_xyz = np.column_stack((firstnetworkspecies_x, firstnetworkspecies_y, firstnetworkspecies_z))

    # iterate through all other speci, and get intersections (within tolerance) of 
    # longitudes / latitudes / measurement_altitudes, with first speci longitudes and latitudes
    pairwise_intersect_inds = {firstnetworkspecies:[]}
    for networkspecies in longitudes:

        if networkspecies == firstnetworkspecies:
            continue

        nextnetworkspecies_longitudes = longitudes[networkspecies]
        nextnetworkspecies_latitudes = latitudes[networkspecies]
        nextnetworkspecies_measurement_altitudes = measurement_altitudes[networkspecies]

        # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
        nextnetworkspecies_x, nextnetworkspecies_y, nextnetworkspecies_z = pyproj.transform(lla, 
            ecef, nextnetworkspecies_longitudes, nextnetworkspecies_latitudes, nextnetworkspecies_measurement_altitudes,
            radians=False)
        
        # merge coordinates to 3D array
        nextnetworkspecies_xyz = np.column_stack((nextnetworkspecies_x, nextnetworkspecies_y, nextnetworkspecies_z))            

        # get all indices of next speci xyz coords, within tolerance of each first speci xyz coords
        idx = cKDTree(nextnetworkspecies_xyz).query_ball_point(firstnetworkspecies_xyz, tolerance)

        # get all indices where have non-duplicated and duplicated matched indices
        unique_idx, unique_idx_counts = np.unique(list(itertools.chain(*idx)), return_counts=True)
        nondup_idx = unique_idx[np.where(unique_idx_counts == 1)[0]]
        dup_idx = unique_idx[np.where(unique_idx_counts > 1)[0]]

        # gather list of indices for first speci and next speci of already resolved within tolerance indices
        fs_wtol_inds = []
        ns_wtol_inds = []

        # resolve all duplicated matched indices
        idx, unresolved_dup_idx, fs_wtol_inds, ns_wtol_inds = resolve_duplicate_spatial_colocation_matches(idx, 
                                                                nondup_idx, dup_idx, 
                                                                firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                                fs_wtol_inds, ns_wtol_inds)

        # pass through again to resolve all unresovered duplicated indices after first pass
        if len(unresolved_dup_idx) > 0:
            idx, _, fs_wtol_inds, ns_wtol_inds = resolve_duplicate_spatial_colocation_matches(idx, 
                                                   nondup_idx, unresolved_dup_idx, 
                                                   firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                   fs_wtol_inds, ns_wtol_inds)

        # iterate though next speci within tolerance indices, per each first speci coord 
        for idx_ii, idx_l in enumerate(idx):
            # if position has already been resolved (through resolving duplicates) then continue
            if idx_ii in fs_wtol_inds:
                continue
            
            # no matches, then append nothing
            elif len(idx_l) == 0:
                continue

            # just 1 match, then append
            elif len(idx_l) == 1:
                fs_wtol_inds.append(idx_ii) 
                ns_wtol_inds.append(idx_l[0])

            # more than 1 match, then find which match has closest distance, then append
            else:
                # find the dists between all relevant first speci xyz coords and next speci duplicate xyz coord
                idx_l = np.array(idx_l)
                dists = cKDTree([firstnetworkspecies_xyz[idx_ii]]).query(nextnetworkspecies_xyz[idx_l], k=1)[0]
                ordered_idx_l = idx_l[np.argsort(dists)]
                fs_wtol_inds.append(idx_ii) 
                ns_wtol_inds.append(ordered_idx_l[0])

        # order matched indices for first species in ascending order, and order next speci indices in smae way
        if len(fs_wtol_inds) > 0:
            fs_wtol_inds, ns_wtol_inds = list(zip(*sorted(zip(fs_wtol_inds, ns_wtol_inds))))

            # set indices where first species differences are within tolerance, i.e. intersecting 
            pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)] = np.array(fs_wtol_inds)
            pairwise_intersect_inds[firstnetworkspecies].extend(np.array(fs_wtol_inds))
        
            # get indices where next species differences are within tolerance, i.e. intersecting 
            pairwise_intersect_inds[networkspecies] = np.array(ns_wtol_inds)
        else:
            pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)] = np.array([], dtype=np.int)
            pairwise_intersect_inds[networkspecies] = np.array([], dtype=np.int)

    # get indices (for first networkspecies) where longitude, latitudes and measurement_altitudes intersect across all species
    pairwise_intersect_inds_unique, counts = np.unique(pairwise_intersect_inds[firstnetworkspecies], return_counts=True)
    pairwise_intersect_inds[firstnetworkspecies] = pairwise_intersect_inds_unique[counts == (len(longitudes)-1)]

    # get specific intersect indices across all species, for rest of species
    intersecting_indices = {}
    for networkspecies in longitudes:
        if networkspecies == firstnetworkspecies:
            intersecting_indices[networkspecies] = np.array(pairwise_intersect_inds[networkspecies], dtype=np.int)
        else:
            _, species_intersect_inds, _ = np.intersect1d(pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)], pairwise_intersect_inds[firstnetworkspecies], return_indices=True)
            intersecting_indices[networkspecies] = np.array(pairwise_intersect_inds[networkspecies][species_intersect_inds], dtype=np.int)

    return intersecting_indices

def resolve_duplicate_spatial_colocation_matches(idx, nondup_idx, dup_idx, 
                                                 firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                 fs_wtol_inds, ns_wtol_inds):

    """ Function that resolves duplicate indices found during spatial colocation of 2 species.

        In spatial colocation it is neccessary to match each stations geographically within a specific
        tolerance, for 2 different species. 

        In some cases the stations within the tolerance can match for multiple stations. 
        In order to resolve this, for each duplicated station, it is set to match with the closest station to it
        in terms of 3D distance. This is done iteratively until there are no more duplicates.

        :param idx: per first networkspeci xyz coords, a list of indices of next networkspeci xyz coords within tolerance
        :type idx: array 
        :param nondup_idx: next networkspeci idx coords that are not duplicated across idx array 
        :type nondup_idx: array
        :param dup_idx: next networkspeci idx coords that are duplicated across idx array 
        :type dup_idx: array
        :param firstnetworkspecies_xyz: ECEF coordinates for first networkspeci stations
        :type firstnetworkspecies_xyz: array
        :param nextnetworkspecies_xyz: ECEF coordinates for next networkspeci stations
        :type nextnetworkspecies_xyz: array
        :param fs_wtol_inds: first networkspeci station indices within tolerance (i.e. have paired match)
        :type fs_wtol_inds: list
        :param ns_wtol_inds: next networkspeci station indices within tolerance (i.e. have paired match)
        :type ns_wtol_inds: list
        :return: idx, unresolved_dup_idx, fs_wtol_inds, ns_wtol_inds 
        :rtype: array, list, list, list
    """

    # resolve all duplicated matched indices by finding for which index the distance is closest
    unresolved_dup_idx = []
    for dup_index in dup_idx:

        # get all relevant first speci indices for which contain the duplicate next speci index
        relevant_idx = np.array([idx_ii for idx_ii, idx_l in enumerate(idx) if dup_index in idx_l])
        # find the dists between all relevant first speci xyz coords and next speci duplicate xyz coord
        dists = cKDTree([nextnetworkspecies_xyz[dup_index]]).query(firstnetworkspecies_xyz[relevant_idx], k=1)[0]
        # order the relevant first speci indices by dists (closest first)
        ordered_relevant_idx = relevant_idx[np.argsort(dists)]
        # iterate through ordered relevant first speci indices until find an ind for which can claim duplicate ind
        # keep going in order iteratively until have exhausted all options for duplicate ind
        for ordered_relevant_index_ii, ordered_relevant_index in enumerate(ordered_relevant_idx):
        
            # if current relevant first speci index has no competing indices, then append matched index,
            # remove duplicate index from other match lists, and then break out of iteration
            if len(idx[ordered_relevant_index]) == 1:
                fs_wtol_inds.append(ordered_relevant_index)
                ns_wtol_inds.append(dup_index)
                for next_ordered_relevant_index in ordered_relevant_idx[ordered_relevant_index_ii+1:]:
                    idx[next_ordered_relevant_index].remove(dup_index)
                break

            # otherwise, 
            # find the dists between the relevant first speci xyz coord and next speci duplicate xyz coords
            else:
                idx_l = np.array(idx[ordered_relevant_index])
                dists = cKDTree([firstnetworkspecies_xyz[ordered_relevant_index]]).query(nextnetworkspecies_xyz[idx_l], k=1)[0]
                ordered_idx_l = idx_l[np.argsort(dists)]

                # if the closest index is the dup index, or is another non-duplicate, then append it
                # if the closest index is the dup index, then remove duplicate index from other match lists 
                # and then break out of iteration
                if (ordered_idx_l[0] == dup_index) or (ordered_idx_l[0] in nondup_idx):
                    fs_wtol_inds.append(ordered_relevant_index)
                    ns_wtol_inds.append(ordered_idx_l[0])
                    if ordered_idx_l[0] == dup_index:
                        for next_ordered_relevant_index in ordered_relevant_idx[ordered_relevant_index_ii+1:]:
                            idx[next_ordered_relevant_index].remove(dup_index)
                        break

                # else, if the closest index is another duplicated index,
                # then cannot currently resolve position, so come back to this later
                else:
                    unresolved_dup_idx.append(dup_index)
                    break

    return idx, unresolved_dup_idx, fs_wtol_inds, ns_wtol_inds 

def show_message(read_instance, msg, msg_offline=None, from_conf=None):

    if read_instance.offline:
        if msg_offline is not None:
            print('Warning: ' + msg_offline)
        else:
            print('Warning: ' + msg)
    
    else:
        # there are some warnings that will only be shown if we launch the dashboard
        # using a configuration file (those in filter.py, read.py and configuration.py)
        if (from_conf is None) or (from_conf):
            from .dashboard_aux import MessageBox
            MessageBox(msg)
