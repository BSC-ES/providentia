"""
Contains functions that are shared by the dashboard and the offline version. 
Currently, it includes function related to the initialization and update
of metadata, checks fields coming from conf files etc.
"""

from .configuration import split_options
from .configuration import parse_path
from .configuration import read_conf

from glob import glob
import os
import copy
import json
import datetime
import sys

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import pyproj
from scipy.spatial import cKDTree
import seaborn as sns

def which_bounds(instance, speci):
    """Returns lower/upper bounds for speci. 
    
    If there are bounds defined in a config file, take those values.
    Otherwise take default GHOST values

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    :param speci: The species current selected for evaluation (e.g. sconco3)
    :type speci: str
    :return: The lower and upper bound
    :rtype: np.float32
    """

    lower = np.float32(instance.parameter_dictionary[speci]['extreme_lower_limit'])
    upper = np.float32(instance.parameter_dictionary[speci]['extreme_upper_limit'])

    speci_index = instance.species.index(speci)

    if hasattr(instance, 'lower_bound'):
        lower_bound_split = [bound.strip() for bound in instance.lower_bound.split(',')]

        if len(lower_bound_split) > 1:
            if len(instance.species) != len(lower_bound_split):
                error = 'Error: "lower_bound" variable must be same length of number of species read.'
                sys.exit(error)
            lower = lower_bound_split[speci_index]
        else:
            lower = lower_bound_split[0]

    if hasattr(instance, 'upper_bound'):
        upper_bound_split = [bound.strip() for bound in instance.upper_bound.split(',')]

        if len(upper_bound_split) > 1:
            if len(instance.species) != len(upper_bound_split):
                error = 'Error: "upper_bound" variable must be same length of number of species read.'
                sys.exit(error)
            upper = upper_bound_split[speci_index]
        else:
            upper = upper_bound_split[0]

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

    if (return_defaults) or (not hasattr(instance, 'qa')):
        if instance.species in instance.met_parameters:
            return sorted(instance.default_qa_met)
        else:
            return sorted(instance.default_qa_standard)

    if hasattr(instance, 'qa'):
        # if conf has only 1 QA
        if isinstance(instance.qa, int):
            return [instance.qa]
        # empty string
        elif instance.qa == "":
            return []
        # if the QAs are written with their names
        elif isinstance(instance.qa, str):
            return sorted([instance.standard_QA_name_to_QA_code[q.strip()] for q in instance.qa.split(",")])
        # list of integer codes
        else:
            return sorted(list(instance.qa))

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
        # empty string
        elif instance.flags == "":
            return []
        # if the flags are written with their names
        elif isinstance(instance.flags, str):
            return [instance.standard_data_flag_name_to_data_flag_code[f.strip()]
                    for f in instance.flags.split(",")]
        # list of integer codes
        else:
            return list(instance.flags)
    else:
        return []

def multi_species_mapping(species):
    """Map species special case str to multiple species names"""

    #multi_species_map = {'vconcaerobin*':['vconcaerobin1','vconcaerobin2','vconcaerobin3','vconcaerobin4','vconcaerobin5','vconcaerobin6','vconcaerobin7','vconcaerobin8','vconcaerobin9','vconcaerobin10','vconcaerobin11','vconcaerobin12','vconcaerobin13','vconcaerobin14','vconcaerobin15','vconcaerobin16','vconcaerobin17','vconcaerobin18','vconcaerobin19','vconcaerobin20','vconcaerobin21','vconcaerobin22']}
    multi_species_map = {'vconcaerobin*':['vconcaerobin7','vconcaerobin8','vconcaerobin9','vconcaerobin10','vconcaerobin11','vconcaerobin12','vconcaerobin13','vconcaerobin14','vconcaerobin15','vconcaerobin16','vconcaerobin17','vconcaerobin18','vconcaerobin19','vconcaerobin20','vconcaerobin21','vconcaerobin22']}

    return multi_species_map[species]

def get_parameters(instance):
    """
    Handle parsing of required config parameters.

    Throw errors if parameters do not exist or are empty strings.

    If have multiple networks / species, parse these correctly.
    Also handling special case strings which map to binned species.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    # network
    if hasattr(instance, 'network'):
        # throw error if network is empty str
        if instance.network.strip() == '':
            error = 'Error: "network" field is empty in .conf file'
            sys.exit(error)
        # parse multiple networks
        elif ',' in instance.network:
            instance.network = [network.strip() for network in instance.network.split(',')]
        else:
            instance.network = [instance.network.strip()]
    else:
        error = 'Error: "network" field must be defined in .conf file'
        sys.exit(error)

    # species
    if hasattr(instance, 'species'):
        # throw error if species is empty str
        if instance.species.strip() == '':
            error = 'Error: "species" field is empty in .conf file'
            sys.exit(error)
        # parse multiple species
        elif ',' in instance.species:
            instance.species = [speci.strip() for speci in instance.species.split(',')]
        else:
            instance.species = [instance.species.strip()]
    # throw error if species is not defined
    else:
        error = 'Error: "species" field must be defined in .conf file'
        sys.exit(error)

    # if number of networks and species is not the same,
    # and len of one of network or species == 1,
    # then duplicate respestive network/species
    if len(instance.network) != len(instance.species):

        # 1 network?
        if len(instance.network) == 1:
            # duplicate network to match species len
            instance.network = instance.network * len(instance.species)

        # 1 species?
        elif len(instance.species) == 1:
            # duplicate species to match network len
            instance.species = instance.species * len(instance.network)

        # otherwise throw error
        else:
            error = 'Error: The number of networks and species is not the same.'
            sys.exit(error)

    # throw error if one of networks are non all GHOST or non-GHOST
    for network_ii, network in enumerate(instance.network):
        if network_ii == 0:
            previous_is_ghost = check_for_ghost(network)
        else:
            is_ghost = check_for_ghost(network)
            if is_ghost != previous_is_ghost:
                error = 'Error: Networks must be all GHOST or non-GHOST'
                sys.exit(error)
            previous_is_ghost = is_ghost

    # resolution
    if hasattr(instance, 'resolution'):
        # throw error if species if empty str
        if instance.resolution.strip() == '':
            error = 'Error: "resolution" field is empty in .conf file'
            sys.exit(error)
    # throw error if resolution is not defined
    else:
        error = 'Error: "resolution" field must be defined in .conf file'
        sys.exit(error)

    # start_date
    if hasattr(instance, 'start_date'):
        # throw error if start_date if empty str
        if str(instance.start_date).strip() == '':
            error = 'Error: "start_date" field is empty in .conf file'
            sys.exit(error)
    # throw error if start_date is not defined
    else:
        error = 'Error: "start_date" field must be defined in .conf file'
        sys.exit(error)

    # end_date
    if hasattr(instance, 'end_date'):
        # throw error if start_date if empty str
        if str(instance.end_date).strip() == '':
            error = 'Error: "end_date" field is empty in .conf file'
            sys.exit(error)
    # throw error if end_date is not defined
    else:
        error = 'Error: "end_date" field must be defined in .conf file'
        sys.exit(error)

    # map to multiple species if have * wildcard
    # also duplicate out associated network
    # remove any species for which there exists no data
    new_species = copy.deepcopy(instance.species)
    for speci_ii, speci in enumerate(instance.species): 
        if '*' in speci:
            mapped_species = multi_species_mapping(speci)
            del new_species[speci_ii]
            new_species[speci_ii:speci_ii] = mapped_species
            network_to_duplicate = instance.network[speci_ii]
            del instance.network[speci_ii]
            instance.network[speci_ii:speci_ii] = [network_to_duplicate]*len(mapped_species)
    instance.species = copy.deepcopy(new_species)

    # if are using dashboard then just take first network/species pair, as multivar not supported yet
    if (len(instance.network) > 1) & (len(instance.species) > 1) & (not instance.offline):
      instance.network = [instance.network[0]]
      instance.species = [instance.species[0]]
      print('Warning: Mutiple networks/species not supported for dashboard.\nFirst network / species taken.')

def get_experiments(instance):
    """If there are experiments coming from a config file,
    select those. Otherwise, return empty dict.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    :return: dict (exp:exp_legend_name)
    :rtype: dict
    """

    if hasattr(instance, 'experiments'):
        # empty string
        if instance.experiments.strip() == "":
            return {}
        # split experiments
        else:
            # have alternative experiment names for the legend, then parse them?
            if ('(' in instance.experiments) & (')' in instance.experiments):
                exps = [exp.strip() for exp in instance.experiments.split('(')[0].strip().split(",")]
                exps_legend = [exp_legend.strip() for exp_legend in instance.experiments.split('(')[1].split(')')[0].strip().split(",")]
            # otherwise set legend names as given experiment names in full
            else: 
                exps = [exp.strip() for exp in instance.experiments.split(",")]
                exps_legend = copy.deepcopy(exps)
            return {exp:exp_legend for exp,exp_legend in zip(exps,exps_legend)}
    else:
        return {}

def get_default_qa_codes(instance):
    """Retrieve default QA codes from GHOST_standards using the QA flags' names.

    A specific selection of qa are defined for met. parameters

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    :return: three lists which contain flags' codes
    :rtype: list
    """

    # get defaukt names from json files
    standard_qa_names = json.load(open(
        "providentia/conf/default_qa.json"))['standard']
    met_qa_names = json.load(open(
        "providentia/conf/default_qa.json"))['met']
    # get qa codes
    standard_qa = [instance.standard_QA_name_to_QA_code[qa_name]
                   for qa_name in standard_qa_names]
    met_qa = [instance.standard_QA_name_to_QA_code[qa_name]
              for qa_name in met_qa_names]

    return standard_qa, met_qa

def exceedance_lim(species):
    """Returns the exceedance limit depending on the species input. If
    species doesn't have a reported limit, returns np.NaN.

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
    """Returns temporal resolution order for menus as a dictionary

    :return: temporal resolution order
    :rtype: dict
    """

    resolution_order_dict = {'hourly': 1, '3hourly': 2, '6hourly': 3, 'hourly_instantaneous': 4,
                             '3hourly_instantaneous': 5, '6hourly_instantaneous': 6,
                             'daily': 7, 'monthly': 8}

    return resolution_order_dict

def temp_axis_dict():
    """Returns temporal mapping as a dictionary used for the plots.

    :return: numbering of months/days
    :rtype: dict
    """
    map_dict = {'dayofweek': {0: 'M', 1: 'T', 2: 'W', 3: 'T', 4: 'F', 5: 'S', 6: 'S'},
                'month': {1: 'J', 2: 'F', 3: 'M', 4: 'A', 5: 'M', 6: 'J',
                          7: 'J', 8: 'A', 9: 'S', 10: 'O', 11: 'N', 12: 'D'}
                }
    return map_dict

def periodic_xticks():
    """Returns xticks for periodic subplots.

    :return dictionary of xticks per temporal resolution
    :rtype dict
    """

    return {'hour':np.arange(24, dtype=np.int), 
            'dayofweek':np.arange(7, dtype=np.int), 
            'month':np.arange(1, 13, dtype=np.int)}

def periodic_labels():
    """Return axes labels for periodic subplots.

    :return: axes labels
    :rtype: dict
    """

    return {'hour':'H', 'dayofweek':'DoW', 'month':'M'}

def get_relevant_temporal_resolutions(resolution):        
    """Get relevant temporal reolsutions for periodic plots, by selected temporal resolution.

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
    return relevant_temporal_resolutions

def get_land_polygon_resolution(selection):
    """get resolution of land polygons to plot on map.

        :param selection: name of selected temporal resolution
    :type resolution: str
    :return: selected land polygon resolution
    :rtype: list
    """
    
    land_polygon_resolutions = {'low': '110m','medium': '50m','high': '10m'}
    return land_polygon_resolutions[selection]

def init_flags(instance):
    """Initialise internal structure to store selected flags.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'flag_menu'):
        instance.flag_menu = {'window_title':'FLAGS', 'page_title':'Select standardised data reporter provided flags to filter by', 'checkboxes':{}}
        instance.flag_menu['select_buttons'] = ['all', 'clear', 'default']
    # reset fields
    instance.flag_menu['checkboxes']['labels'] = np.array(sorted(instance.standard_data_flag_name_to_data_flag_code, key=instance.standard_data_flag_name_to_data_flag_code.get))
    instance.flag_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
    instance.flag_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
    instance.flag_menu['checkboxes']['map_vars'] = np.sort(list(instance.standard_data_flag_name_to_data_flag_code.values()))
    
def init_qa(instance):
    """Initialise internal structure to store selected qa.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'qa_menu'):
        instance.qa_menu = {'window_title':'QA', 'page_title':'Select standardised quality assurance flags to filter by', 'checkboxes':{}}
        instance.qa_menu['select_buttons'] = ['all', 'clear', 'default']
    # reset fields
    instance.qa_menu['checkboxes']['labels'] = np.array(sorted(instance.standard_QA_name_to_QA_code, key=instance.standard_QA_name_to_QA_code.get))
    instance.qa_menu['checkboxes']['remove_default'] = np.array([], dtype=np.uint8)
    instance.qa_menu['checkboxes']['remove_selected'] = np.array([], dtype=np.uint8)
    instance.qa_menu['checkboxes']['map_vars'] = np.sort(list(instance.standard_QA_name_to_QA_code.values()))
    
def init_experiments(instance):
    """Initialise internal structure to store selected experiments.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'experiments_menu'):
        instance.experiments_menu = {'window_title': 'EXPERIMENTS', 'page_title': 'Select Experiment/s', 'checkboxes':{}}
        instance.experiments_menu['select_buttons'] = ['all', 'clear']
    # reset fields
    instance.experiments_menu['checkboxes']['labels'] = [] 
    instance.experiments_menu['checkboxes']['keep_default'] = [] 
    instance.experiments_menu['checkboxes']['keep_selected'] = [] 
    instance.experiments_menu['checkboxes']['map_vars'] = [] 
    
def init_representativity(instance):
    """Initialise internal structure to store representativity fields.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'representativity_menu'):
        instance.representativity_menu = {'window_title': '% DATA REPRESENTATIVITY', 'page_title': 'Select % Data Representativity Bounds', 'rangeboxes':{}}
    # reset fields
    instance.representativity_menu['rangeboxes']['tooltips'] = []
    instance.representativity_menu['rangeboxes']['labels'] = [] 
    instance.representativity_menu['rangeboxes']['current_lower'] = []                                                  
    instance.representativity_menu['rangeboxes']['map_vars'] = []
    instance.representativity_menu['rangeboxes']['subtitles'] = []
    instance.representativity_menu['rangeboxes']['subtitle_inds'] = []

def init_period(instance):
    """Initialise internal structure to store period fields.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    # do not have object instance already?
    # if not, create it
    if not hasattr(instance, 'period_menu'):
        instance.period_menu = {'window_title': 'DATA PERIOD', 'page_title': 'Select Data Periods', 'checkboxes':{}}
    # reset fields
    instance.period_menu['checkboxes']['labels'] = []
    instance.period_menu['checkboxes']['keep_selected'] = []
    instance.period_menu['checkboxes']['remove_selected'] = []

def init_metadata(instance):
    """Initialise internal structure to store metadata.

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
            
        instance.metadata_menu = {'window_title': 'METADATA', 'page_title': 'Select metadata type to filter stations by',
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

        #reset rangebox labels    
        instance.metadata_menu[metadata_type]['rangeboxes']['labels'] = \
            [metadata_var for metadata_var in instance.metadata_vars_to_read
            if (instance.standard_metadata[metadata_var]['metadata_type'] == metadata_type)
            & (instance.standard_metadata[metadata_var]['data_type'] != np.object)]

        #reset rangebox tooltips
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

        #reset checkbox labels            
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
            #otherwise, create infracstructure to store metadata var information
            else:
                instance.metadata_menu[metadata_type][metadata_var] = {'window_title': metadata_var,
                                                                       'page_title': 'Filter stations by unique {} metadata'.format(metadata_var), 
                                                                       'checkboxes': {}}
                instance.metadata_menu[metadata_type][metadata_var]['checkboxes'] = {'labels': [], 
                                                                                     'keep_selected': [], 'keep_default': [],
                                                                                     'remove_selected': [], 'remove_default': []}
        
        #remove metadata type checkbox var dicts not in metadata_vars_to_read
        metadata_type_checkbox_vars = [metadata_type_var for metadata_type_var in instance.metadata_menu[metadata_type].keys() if metadata_type_var not in ['window_title', 'page_title', 'navigation_buttons', 'rangeboxes']]
        metadata_type_checkbox_vars_to_remove = [metadata_type_checkbox_var for metadata_type_checkbox_var in metadata_type_checkbox_vars if metadata_type_checkbox_var not in instance.metadata_vars_to_read]
        for metadata_type_checkbox_var_to_remove in metadata_type_checkbox_vars_to_remove:
            del instance.metadata_menu[metadata_type][metadata_type_checkbox_var_to_remove]

def update_representativity_fields(instance):
    """Update the data representativity menu upon read.

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
            
            instance.representativity_menu['rangeboxes']['labels'] = ['Min Hourly Rep. %',
                                                                      'Max Hourly Gap %',
                                                                      'Min Daily Rep. %',
                                                                      'Max Daily Gap %',
                                                                      'Min Monthly Rep. %',
                                                                      'Max Monthly Gap %',
                                                                      'Min Daily Rep. %',
                                                                      'Max Daily Gap %',
                                                                      'Min Monthly Rep. %',
                                                                      'Max Monthly Gap %',
                                                                      'Min All Rep. %',
                                                                      'Max All Gap %']

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

            instance.representativity_menu['rangeboxes']['labels'] = ['Daily',
                                                                      'Daily Max Gap',
                                                                      'Monthly',
                                                                      'Monthly Max Gap',
                                                                      'All',
                                                                      'All Max Gap']   

            instance.representativity_menu['rangeboxes']['labels'] = ['Min Daily Rep. %',
                                                                      'Max Daily Gap %',
                                                                      'Min Monthly Rep. %',
                                                                      'Max Monthly Gap %',
                                                                      'Min All Rep. %',
                                                                      'Max All Gap %']

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

            instance.representativity_menu['rangeboxes']['labels'] = ['Min Daily Rep. %',
                                                                      'Max Daily Gap %',
                                                                      'Min Monthly Rep. %',
                                                                      'Max Monthly Gap %',
                                                                      'Min Monthly Rep. %',
                                                                      'Max Monthly Gap %',
                                                                      'Min All Rep. %',
                                                                      'Max All Gap %']

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Native', 'Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0, 4]

        # non-GHOST
        else:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['monthly_representativity_percent',
                                                                        'monthly_max_gap_percent',
                                                                        'all_representativity_percent', 
                                                                        'all_max_gap_percent']

            instance.representativity_menu['rangeboxes']['labels'] = ['Min Monthly Rep. %',
                                                                      'Max Monthly Gap %',
                                                                      'Min All Rep. %',
                                                                      'Max All Gap %'] 

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

            instance.representativity_menu['rangeboxes']['labels'] = ['Min Monthly Rep. %',
                                                                      'Max Monthly Gap %',
                                                                      'Min All Rep. %',
                                                                      'Max All Gap %'] 

            instance.representativity_menu['rangeboxes']['subtitles'] = ['Native', 'Averaged']
            instance.representativity_menu['rangeboxes']['subtitle_inds'] = [0, 2]  

        # non-GHOST
        else:
            instance.representativity_menu['rangeboxes']['map_vars'] = ['all_representativity_percent', 
                                                                        'all_max_gap_percent']
            
            instance.representativity_menu['rangeboxes']['labels'] = ['Min All Rep. %',
                                                                      'Max All Gap %'] 

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
    """Update the data period menu upon read.

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
        
        #get all metadata values for field across all networks and species
        meta_var_field = []

        # transform to list if there is only one species / network
        if isinstance(instance.network, str):
            instance.network = [instance.network]
        if isinstance(instance.species, str):
            instance.species = [instance.species]

        for network, speci in zip(instance.network, instance.species):
            networkspeci = '{}|{}'.format(network, speci)
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
            # do not have some numeric values for metadata variable so set as 'nan'
            else:
                instance.metadata_menu[metadata_type]['rangeboxes']['current_lower'][meta_var_index] = 'nan'
                instance.metadata_menu[metadata_type]['rangeboxes']['lower_default'][meta_var_index] = 'nan'
                instance.metadata_menu[metadata_type]['rangeboxes']['current_upper'][meta_var_index] = 'nan'
                instance.metadata_menu[metadata_type]['rangeboxes']['upper_default'][meta_var_index] = 'nan'
                if meta_var in instance.metadata_menu[metadata_type]['rangeboxes']['apply_selected']:
                    instance.metadata_menu[metadata_type]['rangeboxes']['apply_selected'].remove(meta_var)

def representativity_conf(instance):
    """Function used when loading from a configuration file. 
    Sets defined representativity filter variables, rest of variables are set to default. 

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    for i, label in enumerate(instance.representativity_menu['rangeboxes']['map_vars']):
        if hasattr(instance, label):
            instance.representativity_menu['rangeboxes']['current_lower'][i] = str(getattr(instance, label))
        else:
            if 'max_gap' in label:
                instance.representativity_menu['rangeboxes']['current_lower'][i] = '100'
            else:
                instance.representativity_menu['rangeboxes']['current_lower'][i] = '0'

def period_conf(instance):
    """Function used when loading from a configuration file. 
    Sets defined period filter variables, rest of variables are set to default. 

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

    if hasattr(instance, 'period'):
        keeps, removes = split_options(instance.period)
        instance.period_menu['checkboxes']['keep_selected'] = keeps
        instance.period_menu['checkboxes']['remove_selected'] = removes
    else:
        instance.period_menu['checkboxes']['keep_selected'] =  []
        instance.period_menu['checkboxes']['remove_selected'] = []

def metadata_conf(instance):
    """Function used when loading from a configuration file. 
    Sets defined metadata filter variables, rest of variables are set to default.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object
    """

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
                current_lower.append('nan')
                current_upper.append('nan')
        instance.metadata_menu[menu_type]['rangeboxes']['current_lower'] = current_lower
        instance.metadata_menu[menu_type]['rangeboxes']['current_upper'] = current_upper
        instance.metadata_menu[menu_type]['rangeboxes']['apply_selected'] = apply_selected

        # and then treat the keep/remove
        for label in instance.metadata_menu[menu_type]['navigation_buttons']['labels']:
            if hasattr(instance, label):
                keeps, removes = split_options(getattr(instance, label))
                instance.metadata_menu[menu_type][label]['checkboxes']['keep_selected'] = keeps
                instance.metadata_menu[menu_type][label]['checkboxes']['remove_selected'] = removes
            else:
                instance.metadata_menu[menu_type][label]['checkboxes']['keep_selected'] = []
                instance.metadata_menu[menu_type][label]['checkboxes']['remove_selected'] = []

def valid_date(date_text):
    """Determines if a date string is in the correct format."""

    try:
        datetime.datetime.strptime(str(date_text), '%Y%m%d')
        return True
    except Exception as e:
        return False

def check_for_ghost(network_name):
    """
    Check whether the selected network comes from GHOST or not.
    All non-GHOST networks start with an asterisk at their name.
    """

    if '/' in network_name:
        return False
    else:
        return True


def load_conf(self, fpath=None):
    """Load existing configurations from file
    for running offline Providentia."""

    if fpath is None:
        print("No configuration file found")
        sys.exit(1)

    # if DEFAULT is not present, then return
    if not os.path.isfile(fpath):
        print(("Error %s" % fpath))
        return

    self.sub_opts, self.all_sections, self.parent_section_names, self.subsection_names, self.filenames = read_conf(fpath)

def update_plotting_parameters(instance):
    """
    Function that updates plotting parameters (colour and zorder) for data labels

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
    """
    Create GHOST observational data tree as a nested dictionary,
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
    """
    Fill non-GHOST observational data tree,
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

    #iterate through networks
    for network in nonghost_observation_data_json:

        # check if directory for network exists
        # if not, continue
        if not os.path.exists('%s/%s' % (instance.nonghost_root, network)):
            continue

        # write empty dictionary for network
        nonghost_observation_data[network] = {}

        #iterate through resolutions
        for resolution in nonghost_observation_data_json[network]:

            # check if directory for resolution exists
            # if not, continue
            if not os.path.exists('%s/%s/%s' % (instance.nonghost_root, network, resolution)):
                continue

            # write nested empty dictionary for resolution
            nonghost_observation_data[network][resolution] = {}

            #iterate through species
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
    """
    Fill experiment data tree,
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
                            files_directory = '%s/%s/%s/%s/%s/*%s' % (instance.exp_root, instance.ghost_version, 
                                                                      experiment, resolution, species,
                                                                      network.split('/')[0].upper())
                            
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
    """
    Iterate through observational dictionary tree and return 
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
        print('Warning: One of start date or end date are not valid.')
        return

    # check end date is > start date, if not, return with no valid obs. files
    if int(start_date) >= int(end_date):
        print('Warning: Start date exceeds end date.')
        return

    # check start date and end date are both within if valid date range (19000101 - 20500101),
    # if not, return with no valid obs. files
    if (int(start_date) < 19000101) or (int(end_date) < 19000101) or (int(start_date) >= 20500101) or (int(end_date) >= 20500101):
        print('Warning: One of start date or end date are not valid.')
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
    """
    Get valid experiments for daterange, and selected parameters.
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

    #iterate through networks and species
    for network, speci in zip(networks, species):

        #iterate through available experiments
        for experiment in available_experiments:

            # get folder where interpolated experiments are saved
            if '/' not in network:           
                files_directory = '%s/%s/%s/%s/%s/%s' % (instance.exp_root, instance.ghost_version, 
                                                         experiment, resolution, speci, network)
            else:
                files_directory = '%s/%s/%s/%s/%s/*%s' % (instance.exp_root, instance.ghost_version, 
                                                          experiment, resolution, speci,
                                                          network.split('/')[0].upper())
                
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

def spatial_colocation(reading_ghost, station_references, longitudes, latitudes):
    """ 
    Given multiple species, return intersecting station_references, longitudes and latitudes 
    across species.

    This is done by 
        1. Cross-checking the station references between species to get matching station_references
        2. Cross-checking matching longitude/latitude coordinates to a tolerance of 20m difference

    If longitude/latitude matching is needed, then take impose station reference of the first species 
    upon the rest of intersecting indices across species. 

    :param reading_ghost: boolean informing if are using GHOST data or not
    :type reading_ghost: boolean
    :param station_references: dictionary of station references per network/species
    :type station_references: dict
    :param longitudes: dictionary longitudes per network/species
    :type longitudes: dict
    :param latitudes: dictionary of latitudes per network/species
    :type latitudes: dict
    :return: intersecting station_references, intersecting longitudes, intersecting latitudes 
    :rtype: list, list, list
    """

    # if are reading ghost data remove method abbreviation from station_references
    if reading_ghost:
        station_references_no_method = {}
        for networkspecies in station_references:
            station_references_no_method[networkspecies] = ['_'.join(ref.split('_')[:-1]) for ref in station_references[networkspecies]]
    else:
        station_references_no_method = station_references

    # get indices of intersection of station references across species
    intersecting_station_references_no_method = list(set.intersection(*map(set,list(station_references_no_method.values()))))
    intersecting_indices = {}
    for networkspecies in station_references_no_method:
        intersecting_indices[networkspecies] = np.sort([station_references_no_method[networkspecies].index(ref) for ref in intersecting_station_references_no_method])

    # set intersecting station_references, longitudes and latitudes,
    # after matching station references across species
    firstnetworkspecies = list(intersecting_indices.keys())[0]
    if len(intersecting_indices[firstnetworkspecies]) > 0:
        intersecting_station_references = np.array(station_references_no_method[firstnetworkspecies])[intersecting_indices[firstnetworkspecies]]
        intersecting_longitudes = np.array(longitudes[firstnetworkspecies])[intersecting_indices[firstnetworkspecies]]
        intersecting_latitudes = np.array(latitudes[firstnetworkspecies])[intersecting_indices[firstnetworkspecies]]
    else:
        intersecting_station_references = np.array([])
        intersecting_longitudes = np.array([])
        intersecting_latitudes = np.array([])

    # if non-intersecting indices unaccounted for across species, 
    # then attempt to resolve them by matching longitudes / latitudes
    if len(intersecting_station_references) != len(station_references_no_method[firstnetworkspecies]):

        # set tolerance for matching longitudes and latitudes in metres
        tolerance = 20

        # get non-intersecting station references, longitudes and latitudes across speci
        non_intersecting_station_references = {networkspecies: (np.delete(station_references_no_method[networkspecies], intersecting_indices[networkspecies]) if len(intersecting_indices[networkspecies]) > 0 else np.array(station_references_no_method[networkspecies])) for networkspecies in station_references_no_method}
        non_intersecting_longitudes = {networkspecies: (np.delete(longitudes[networkspecies], intersecting_indices[networkspecies]) if len(intersecting_indices[networkspecies]) > 0 else np.array(longitudes[networkspecies])) for networkspecies in longitudes}
        non_intersecting_latitudes = {networkspecies: (np.delete(latitudes[networkspecies], intersecting_indices[networkspecies]) if len(intersecting_indices[networkspecies]) > 0 else np.array(latitudes[networkspecies])) for networkspecies in latitudes}

        # get non-intersecting station references, longitudes and latitudes for first speci
        speci_non_intersecting_station_references = non_intersecting_station_references[firstnetworkspecies]
        speci_non_intersecting_longitudes = non_intersecting_longitudes[firstnetworkspecies]
        speci_non_intersecting_latitudes = non_intersecting_latitudes[firstnetworkspecies]
        # convert speci longitude and latitudes in geogroahic coordinates to cartesian ECEF 
        # (Earth Centred, Earth Fixed) coordinates assuming WGS84 datum and ellipsoid, and that all heights equal zero
        # ECEF coordiantes represent positions (in metres) as X, Y, Z coordinates, approximating the earth surface as an ellipsoid of revolution
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        speci_non_intersecting_x, speci_non_intersecting_y, speci_non_intersecting_z = pyproj.transform(lla, ecef, speci_non_intersecting_longitudes, speci_non_intersecting_latitudes, np.zeros(len(speci_non_intersecting_longitudes)), radians=False)
        #merge coordinates to 2D array
        speci_non_intersecting_xy = np.column_stack((speci_non_intersecting_x, speci_non_intersecting_y))

        # iterate through all other speci, and get intersections (within tolerance) of longitudes and latitudes 
        # with first speci longitudes and latitudes
        pairwise_intersect_inds = []
        for networkspecies in non_intersecting_longitudes:

            if networkspecies == firstnetworkspecies:
                continue

            next_speci_non_intersecting_longitudes = non_intersecting_longitudes[networkspecies]
            next_speci_non_intersecting_latitudes = non_intersecting_latitudes[networkspecies]

            # convert speci longitude and latitudes in geogroahic coordinates to cartesian ECEF 
            next_speci_non_intersecting_x, next_speci_non_intersecting_y, next_speci_non_intersecting_z = pyproj.transform(lla, ecef, next_speci_non_intersecting_longitudes, next_speci_non_intersecting_latitudes, np.zeros(len(next_speci_non_intersecting_longitudes)), radians=False)
            # merge coordinates to 2D array
            next_speci_non_intersecting_xy = np.column_stack((next_speci_non_intersecting_x, next_speci_non_intersecting_y))            

            # get closest differences between next speci lon,lat coords, with first speci lon lats coords
            dists = cKDTree(next_speci_non_intersecting_xy).query(speci_non_intersecting_xy, k=1)[0]
            
            # get indices where differences are within tolerance, i.e. intersecting 
            pairwise_intersect_inds.extend(np.where(dists <= tolerance)[0])

        # get indices where longitude and latitudes intersect across all species
        pairwise_intersect_inds_unique, counts = np.unique(pairwise_intersect_inds, return_counts=True)
        intersecting_indices_lonlat = pairwise_intersect_inds_unique[counts == (len(longitudes)-1)]

        # append newly found intersecting station references, longitudes and latitudes
        if len(intersecting_indices_lonlat) > 0:
            intersecting_station_references = np.append(intersecting_station_references, speci_non_intersecting_station_references[intersecting_indices_lonlat])
            intersecting_longitudes = np.append(intersecting_longitudes, speci_non_intersecting_longitudes[intersecting_indices_lonlat])
            intersecting_latitudes = np.append(intersecting_latitudes, speci_non_intersecting_latitudes[intersecting_indices_lonlat])
        
    # sort arrays by station references (alphabetically)
    sorted_inds = intersecting_station_references.argsort()
    intersecting_station_references = intersecting_station_references[sorted_inds]
    intersecting_longitudes = intersecting_longitudes[sorted_inds]
    intersecting_latitudes = intersecting_latitudes[sorted_inds]

    return intersecting_station_references, intersecting_longitudes, intersecting_latitudes

def get_basic_metadata(instance, networks, species, resolution):     
    """
    Get basic unique metadata across networks / species wanting to read
    The basic fields are: station_reference, longitude, latitude

    If have multiple species, then spatially cocolocate across species 
    to get matching stations across stations.

    :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
    :type instance: object

    :param networks: list of networks
    :type networks: list
    :param species: list of species
    :type species: list
    :param resolution: selected temporal resolution
    :type resolution: str
    :return: station_references per network/species, longitudes per network/species, latitudes per network/species 
    :rtype: dict, dict, dict
    """

    # define dictionaries for storing basic metadata across all species to read
    station_references = {}
    station_longitudes = {}
    station_latitudes = {}

    # iterate through network, speci pair
    for network, speci in zip(networks, species):
    
        # get species matrix
        matrix = instance.parameter_dictionary[speci]['matrix']

        # get file root
        # GHOST
        if instance.reading_ghost:
            file_root = '%s/%s/%s/%s/%s/%s_' % (instance.ghost_root, network,
                                                instance.ghost_version, resolution,
                                                speci, speci)
        # non-GHOST
        else:
            file_root = '%s/%s/%s/%s/%s_' % (instance.nonghost_root, network,
                                                resolution, speci, speci)

        # get relevant files
        relevant_files = sorted([file_root+str(yyyymm)+'.nc' for yyyymm in instance.yearmonths])
    
        # get station references, longitudes and latitudes for speci
        # GHOST
        if instance.reading_ghost:
            
            # define arrays for storing speci metadata
            speci_station_references = []
            speci_station_longitudes = []
            speci_station_latitudes = []

            for relevant_file in relevant_files:
                ncdf_root = Dataset(relevant_file)
                speci_station_references = np.append(speci_station_references, ncdf_root['station_reference'][:])
                speci_station_longitudes = np.append(speci_station_longitudes, ncdf_root['longitude'][:])
                speci_station_latitudes = np.append(speci_station_latitudes, ncdf_root['latitude'][:])
                ncdf_root.close()

            speci_station_references, station_unique_indices = np.unique(speci_station_references, return_index=True)
            station_references['{}|{}'.format(network, speci)] = speci_station_references
            station_longitudes['{}|{}'.format(network, speci)] = speci_station_longitudes[station_unique_indices]
            station_latitudes['{}|{}'.format(network, speci)] = speci_station_latitudes[station_unique_indices]
        
        # non-GHOST
        else:
            
            ncdf_root = Dataset(relevant_files[0])
            station_references['{}|{}'.format(network, speci)] = np.array(
                [st_name.tostring().decode('ascii').replace('\x00', '')
                for st_name in ncdf_root['station_name'][:]], dtype=np.str)
            if "latitude" in ncdf_root.variables:
                station_longitudes['{}|{}'.format(network, speci)] = ncdf_root['longitude'][:]
                station_latitudes['{}|{}'.format(network, speci)] = ncdf_root['latitude'][:]
            else:
                station_longitudes['{}|{}'.format(network, speci)] = ncdf_root['lon'][:]
                station_latitudes['{}|{}'.format(network, speci)] = ncdf_root['lat'][:]
            ncdf_root.close()

    # if have more than 1 species to read, and spatial_colocation is active,
    # then spatially colocate stations across species
    if (len(species) > 1) & (instance.spatial_colocation):
        intersecting_station_references, intersecting_station_longitudes, intersecting_station_latitudes = \
            spatial_colocation(instance.reading_ghost, station_references, station_longitudes, station_latitudes)
        for networkspecies in station_references:
            station_references[networkspecies] = intersecting_station_references
            station_longitudes[networkspecies] = intersecting_station_longitudes
            station_latitudes[networkspecies] = intersecting_station_latitudes

    return station_references, station_longitudes, station_latitudes

def filter_by_species():
    pass