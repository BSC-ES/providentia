""" Contains functions that are shared by the dashboard and the offline version. 
    Currently, it includes function related to the initialization and update
    of metadata, checks fields coming from conf files etc.
"""

from glob import glob
import os
import copy
import json
import datetime
import sys
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import math
import pyproj
from scipy.spatial import cKDTree
import seaborn as sns


def get_default_qa(instance, speci):
    """ Return the default values according to GHOST standards. 

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :return: QA flags' codes in list
        :rtype: list
    """

    if speci in instance.met_parameters:
        return sorted(instance.default_qa_met)
    else:
        return sorted(instance.default_qa_standard)

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
    instance.multispecies_menu['multispecies']['labels'] = ['networkspeci_0']
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
        for (networkspeci_ii, networkspeci), bounds in zip(enumerate(instance.filter_species.keys()),
                                                                     instance.filter_species.values()):
            
            # update menu_current
            if networkspeci_ii > 0:
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

            update_filter_species(instance, networkspeci_ii)

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
        :param data_label: label of data array
        :type data_label: text
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
        The basic fields are: station_reference, longitude, latitude, station_classification and area_classification

        If have multiple species, then spatially cocolocate across species 
        to get matching stations across stations.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
        :return: station_references per networkspecies, longitudes per networkspecies, latitudes per networkspecies, station_classifications per networkspecies, area_classifications per networkspecies 
        :rtype: dict, dict, dict, dict, dict
    """

    # define dictionaries for storing basic metadata across all species to read
    station_references = {}
    station_names = {}
    station_longitudes = {}
    station_latitudes = {}
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
            speci_station_classifications = []
            speci_area_classifications = []

            for relevant_file in relevant_files:
                ncdf_root = Dataset(relevant_file)
                speci_station_references = np.append(speci_station_references, ncdf_root['station_reference'][:])
                speci_station_names = np.append(speci_station_names, ncdf_root['station_name'][:])
                speci_station_longitudes = np.append(speci_station_longitudes, ncdf_root['longitude'][:])
                speci_station_latitudes = np.append(speci_station_latitudes, ncdf_root['latitude'][:])
                speci_station_classifications = np.append(speci_station_classifications, ncdf_root['station_classification'][:])
                speci_area_classifications = np.append(speci_area_classifications, ncdf_root['area_classification'][:])
                ncdf_root.close()

            speci_station_references, station_unique_indices = np.unique(speci_station_references, return_index=True)
            station_references[networkspeci] = speci_station_references
            station_names[networkspeci] = speci_station_names
            station_longitudes[networkspeci] = speci_station_longitudes[station_unique_indices]
            station_latitudes[networkspeci] = speci_station_latitudes[station_unique_indices]
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
            if ncdf_root[station_reference_var].dtype == np.str:
                station_references[networkspeci] = ncdf_root[station_reference_var][:]
            else:
                if ncdf_root[station_reference_var].dtype == np.dtype(object):
                    station_references[networkspeci] = np.array([''.join(val) for val in ncdf_root[station_reference_var][:]])
                else:
                    station_references[networkspeci] = np.array(
                        [st_name.tostring().decode('ascii').replace('\x00', '')
                        for st_name in ncdf_root[station_reference_var][:]], dtype=np.str)
            
            # get indices of all non-NaN stations (can be NaN for some non-GHOST files)
            non_nan_station_indices = np.array([ref_ii for ref_ii, ref in enumerate(station_references[networkspeci]) if ref.lower() != 'nan'])
            station_references[networkspeci] = station_references[networkspeci][non_nan_station_indices]

            if "station_name" in ncdf_root.variables:
                if ncdf_root['station_name'].dtype == np.str:
                    station_names[networkspeci] = ncdf_root['station_name'][non_nan_station_indices]
                else:
                    if ncdf_root['station_name'].dtype == np.dtype(object):
                        station_names[networkspeci] = np.array([''.join(val) for val in ncdf_root['station_name'][non_nan_station_indices]])
                    else:
                        station_names[networkspeci] = np.array(
                            [st_name.tostring().decode('ascii').replace('\x00', '')
                            for st_name in ncdf_root['station_name'][non_nan_station_indices]], dtype=np.str)

            if "latitude" in ncdf_root.variables:
                station_longitudes[networkspeci] = ncdf_root['longitude'][non_nan_station_indices]
                station_latitudes[networkspeci] = ncdf_root['latitude'][non_nan_station_indices]
            else:
                station_longitudes[networkspeci] = ncdf_root['lon'][non_nan_station_indices]
                station_latitudes[networkspeci] = ncdf_root['lat'][non_nan_station_indices]

            if "station_classification" in ncdf_root.variables:
                if ncdf_root['station_classification'].dtype == np.str:
                    station_classifications[networkspeci] = ncdf_root['station_classification'][non_nan_station_indices]
                else:
                    if ncdf_root['station_classification'].dtype == np.dtype(object):
                        station_classifications[networkspeci] = np.array([''.join(val) for val in ncdf_root['station_classification'][non_nan_station_indices]])
                    else:
                        station_classifications[networkspeci] = np.array(
                            [st_classification.tostring().decode('ascii').replace('\x00', '')
                            for st_classification in ncdf_root['station_classification'][non_nan_station_indices]], dtype=np.str)
            
            if "area_classification" in ncdf_root.variables:
                if ncdf_root['area_classification'].dtype == np.str:
                    area_classifications[networkspeci] = ncdf_root['area_classification'][non_nan_station_indices]
                else:
                    if ncdf_root['area_classification'].dtype == np.dtype(object):
                        area_classifications[networkspeci] = np.array([''.join(val) for val in ncdf_root['area_classification'][non_nan_station_indices]]) 
                    else:
                        area_classifications[networkspeci] = np.array(
                            [area_classification.tostring().decode('ascii').replace('\x00', '')
                            for area_classification in ncdf_root['area_classification'][non_nan_station_indices]], dtype=np.str)
            
            # get non-GHOST measurement units
            nonghost_units[speci] = ncdf_root[speci].units

            ncdf_root.close()

    # if have more than 1 networkspecies (including filter networkspecies), and spatial_colocation is active,
    # then spatially colocate stations across species
    if (len((instance.networkspecies + instance.filter_networkspecies)) > 1) & (instance.spatial_colocation):
        # get intersecting station indices across species
        intersecting_indices = spatial_colocation(instance.reading_ghost, station_references, station_longitudes, station_latitudes)
        
        # iterate through networkspecies specific intersecting indices, setting 
        for ns, ns_intersects in intersecting_indices.items():
            station_references[ns] = station_references[ns][ns_intersects]
            station_names[ns] = station_names[ns][ns_intersects]
            station_longitudes[ns] = station_longitudes[ns][ns_intersects]
            station_latitudes[ns] = station_latitudes[ns][ns_intersects]
            station_classifications[ns] = station_classifications[ns][ns_intersects]
            area_classifications[ns] =  area_classifications[ns][ns_intersects] 

    return station_references, station_names, station_longitudes, station_latitudes, station_classifications, area_classifications, nonghost_units

def spatial_colocation(reading_ghost, station_references, longitudes, latitudes):
    """ Given multiple species, return intersecting indices for matching stations across species (per network/species).

        This is done by 
            1. Cross-checking the station references between species to get matching station_references
            2. Cross-checking matching longitude/latitude coordinates to a tolerance of 20m difference

        :param reading_ghost: boolean informing if are using GHOST data or not
        :type reading_ghost: boolean
        :param station_references: dictionary of station references per network/species
        :type station_references: dict
        :param longitudes: dictionary of longitudes per network/species
        :type longitudes: dict
        :param latitudes: dictionary of latitudes per network/species
        :type latitudes: dict
        :return: intersecting indices per network/species
        :rtype: dict
    """

    # if are reading GHOST data remove method abbreviation from station_references
    if reading_ghost:
        station_references_no_method = {}
        for networkspecies in station_references:
            station_references_no_method[networkspecies] = ['_'.join(ref.split('_')[:-1]) 
                                                            for ref in station_references[networkspecies]]
    else:
        station_references_no_method = station_references

    # get indices of intersection of station references across species
    intersecting_station_references_no_method = list(set.intersection(*map(set,list(station_references_no_method.values()))))
    intersecting_indices = {}
    for networkspecies in station_references_no_method:
        intersecting_indices[networkspecies] = np.array(np.sort([station_references_no_method[networkspecies].index(ref) 
                                                        for ref in intersecting_station_references_no_method]), dtype=np.int)

    # set variable for first networkspecies
    firstnetworkspecies = list(intersecting_indices.keys())[0]

    # if have zero intersecting indices across species, then return with warning message
    if len(intersecting_indices[firstnetworkspecies]) == 0:
        print('Warning: No intersecting stations across networks/species')
        return intersecting_indices

    # if non-intersecting indices unaccounted for across species, 
    # then attempt to resolve them by matching longitudes / latitudes
    if len(intersecting_indices[firstnetworkspecies]) != len(station_references_no_method[firstnetworkspecies]):

        # set tolerance for matching longitudes and latitudes in metres
        tolerance = 20

        # get non-intersecting indices, longitudes and latitudes across speci
        non_intersecting_indices = {networkspecies: np.setdiff1d(np.arange(len(station_references[networkspecies])), intersecting_indices[networkspecies]) for networkspecies in station_references}
        non_intersecting_longitudes = {networkspecies: longitudes[networkspecies][non_intersecting_indices[networkspecies]] for networkspecies in longitudes}
        non_intersecting_latitudes = {networkspecies: latitudes[networkspecies][non_intersecting_indices[networkspecies]] for networkspecies in latitudes}

        # get non-intersecting station longitudes and latitudes for first speci
        speci_non_intersecting_longitudes = non_intersecting_longitudes[firstnetworkspecies]
        speci_non_intersecting_latitudes = non_intersecting_latitudes[firstnetworkspecies]

        # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
        # (Earth Centred, Earth Fixed) coordinates assuming WGS84 datum and ellipsoid, and that all heights equal zero
        # ECEF coordiantes represent positions (in metres) as X, Y, Z coordinates, approximating the earth surface as an ellipsoid of revolution
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        speci_non_intersecting_x, speci_non_intersecting_y, speci_non_intersecting_z = pyproj.transform(lla, ecef, 
            speci_non_intersecting_longitudes, speci_non_intersecting_latitudes, 
            np.zeros(len(speci_non_intersecting_longitudes)), radians=False)
        # merge coordinates to 2D array
        speci_non_intersecting_xy = np.column_stack((speci_non_intersecting_x, speci_non_intersecting_y))

        # iterate through all other speci, and get intersections (within tolerance) of longitudes and latitudes 
        # with first speci longitudes and latitudes
        pairwise_intersect_inds = {firstnetworkspecies:[]}
        for networkspecies in non_intersecting_longitudes:

            if networkspecies == firstnetworkspecies:
                continue

            next_speci_non_intersecting_longitudes = non_intersecting_longitudes[networkspecies]
            next_speci_non_intersecting_latitudes = non_intersecting_latitudes[networkspecies]

            # convert speci longitude and latitudes in geogroahic coordinates to cartesian ECEF 
            next_speci_non_intersecting_x, next_speci_non_intersecting_y, next_speci_non_intersecting_z = pyproj.transform(lla, 
                ecef, next_speci_non_intersecting_longitudes, next_speci_non_intersecting_latitudes, 
                np.zeros(len(next_speci_non_intersecting_longitudes)), radians=False)
            
            # merge coordinates to 2D array
            next_speci_non_intersecting_xy = np.column_stack((next_speci_non_intersecting_x, next_speci_non_intersecting_y))            

            # get closest differences between next speci lon,lat coords, with first speci lon lats coords (in dimensions of first species)
            dists_fs = cKDTree(next_speci_non_intersecting_xy).query(speci_non_intersecting_xy, k=1)[0]
            
            # get closest differences between first speci lon lats coords, with next speci lon,lat coords (in dimensions of next species)
            dists_ns = cKDTree(speci_non_intersecting_xy).query(next_speci_non_intersecting_xy, k=1)[0]

            # get indices where first species differences are within tolerance, i.e. intersecting 
            pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)] = non_intersecting_indices[firstnetworkspecies][np.where(dists_fs <= tolerance)[0]]
            pairwise_intersect_inds[firstnetworkspecies].extend(copy.deepcopy(pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)]))
           
            # get indices where next species differences are within tolerance, i.e. intersecting 
            pairwise_intersect_inds[networkspecies] = non_intersecting_indices[networkspecies][np.where(dists_ns <= tolerance)[0]]

        # get indices (for first networkspecies) where longitude and latitudes intersect across all species
        pairwise_intersect_inds_unique, counts = np.unique(pairwise_intersect_inds[firstnetworkspecies], return_counts=True)
        pairwise_intersect_inds[firstnetworkspecies] = pairwise_intersect_inds_unique[counts == (len(longitudes)-1)]

        # get specific intersect indices across all species, for rest of species
        for networkspecies in non_intersecting_longitudes:
            if networkspecies == firstnetworkspecies:
                continue
            _, species_intersect_inds, _ = np.intersect1d(pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)], pairwise_intersect_inds[firstnetworkspecies], return_indices=True)
            pairwise_intersect_inds[networkspecies] = pairwise_intersect_inds[networkspecies][species_intersect_inds]

        # append newly found intersecting indices to previously found intersect inds
        # while doing so also manke sure indices are sorted
        for networkspecies in non_intersecting_longitudes:
            intersecting_indices[networkspecies] = np.array(np.sort(np.append(intersecting_indices[networkspecies], pairwise_intersect_inds[networkspecies])), dtype=np.int)

    return intersecting_indices

def update_filter_species(instance, label_ii, add_filter_species=True):
    """ Function to update filter species after launching the dashboard with a configuration file or 
        by editing the fields in the multispecies filtering tab in the dashboard. 

        :param instance: Instance of class ProvidentiaMainWindow
        :type instance: object
        :param label_ii: Corresponding widget line in dashboard
        :type label_ii: int
        :param add_filter_species: boolean to indicate if networkspeci has to be added or removed
        :type add_filter_species: boolean
    """

    # get selected network, species and bounds
    network = instance.selected_widget_network[label_ii]
    speci = instance.selected_widget_species[label_ii]
    networkspeci = network + '|' + speci
    current_lower = instance.selected_widget_lower[label_ii]
    current_upper = instance.selected_widget_upper[label_ii]
    current_filter_species_fill_value = instance.selected_widget_filter_species_fill_value[label_ii]

    # if apply button is checked or filter_species in configuration file, add networkspecies in filter_species
    if add_filter_species:

        # do not add networkspeci to filter_species if filtered network is main
        if networkspeci != instance.networkspecies[0]:
            
            # add or update networkspeci
            # check selected lower and upper bounds and fill value are numbers or nan
            try:
                instance.filter_species[networkspeci] = [float(current_lower), float(current_upper),
                                                         float(current_filter_species_fill_value)]
                
            # if any of the fields are not numbers, return from function
            except ValueError:
                print("Warning: Data limit fields must be numeric")
                return

            # add if networkspeci is not already in qa_per_species
            if speci not in instance.qa_per_species:
                # get species in memory 
                species = copy.deepcopy(instance.species)
                filter_species = [val.split('|')[1] for val in list(copy.deepcopy(instance.filter_species).keys())]
                qa_species = species + filter_species
                
                # add
                qa_species.append(speci)
                instance.qa_per_species = {speci:get_default_qa(instance, speci) for speci in qa_species}

        # update le_minimum_value and le_maximum_value (data bounds) for main networkspeci
        else:
            instance.le_minimum_value.setText(current_lower)
            instance.le_maximum_value.setText(current_upper)
            instance.bounds_set_on_multispecies = True

    # if apply button is unchecked, remove networkspecies from filter_species
    else:
        # remove from filter_species
        if networkspeci in instance.filter_species.keys():
            del instance.filter_species[networkspeci]

        # remove from qa_per_species
        if speci in instance.qa_per_species:
            del instance.qa_per_species[speci]

    print('SELECTED SPECIES')
    print('- Network species', instance.networkspecies)
    print('- Filter species', instance.filter_species)
    