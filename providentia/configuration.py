""" Providentia Configuration Module """

import configparser
import copy
import json
import os
import sys
import re
import subprocess

from .aux import check_for_ghost, multi_species_mapping, get_default_qa

import numpy as np
import pandas as pd

MACHINE = os.environ.get('BSC_MACHINE', '')

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

def parse_path(dir, f):
    if os.path.isabs(f):
        return f
    else:
        return os.path.join(dir, f)

class ProvConfiguration:
    """ Configuration parameters definitions """

    def __init__(self, read_instance, **kwargs):
        
        self.read_instance = read_instance 

        # set variable defaults
        var_defaults = {
            'ghost_version': '1.4',
            'config_dir': os.path.join(os.environ['HOME'], '.providentia'),
            'cartopy_data_dir': '',
            'available_cpus': '',
            'n_cpus': '',
            'ghost_root': '',
            'nonghost_root': '',
            'exp_root': '',
            'offline': False,
            'available_resolutions': ['hourly', '3hourly', '6hourly', 'hourly_instantaneous',
                                      '3hourly_instantaneous', '6hourly_instantaneous',
                                      'daily', 'monthly'],
            'available_networks': ['AERONET_v3_lev1.5','AERONET_v3_lev2.0','CANADA_NAPS','CAPMoN','CHILE_SINCA',
                                   'EANET','EBAS','EEA_AIRBASE','EEA_AQ_eReporting','JAPAN_NIES','MEXICO_CDMX',
                                   'MITECO','NOAA_ISD','NOAA_ISD_EU','NOAA_ISD_IP','NOAA_ISD_NA',
                                   'SEARCH','UK_AIR','US_EPA_AQS','US_EPA_CASTNET','US_NADP_AMNet','US_NADP_AMoN',
                                   'WMO_WDCGG'], 
            'network': 'EBAS',
            'species': 'sconco3',
            'resolution': 'hourly',
            'start_date': '20180101',
            'end_date': '20180201',
            'experiments': {},
            'qa': None,
            'flags': None,
            'temporal_colocation': False,
            'spatial_colocation': True,
            'map_extent': None, 
            'filter_species': {},
            'lower_bound': None,
            'upper_bound': None,
            'report_type': 'standard',
            'report_summary': True,
            'report_stations': False,
            'report_title': 'Providentia Offline Report',
            'report_filename': 'PROVIDENTIA_Report',
            'active_dashboard_plots': ['timeseries', 'statsummary', 'distribution', 'periodic'],
            'plot_characteristics_filename': '',
            'fixed_section_vars':  ['ghost_version', 'config_dir', 'cartopy_data_dir', 'available_cpus', 'n_cpus',
                                    'ghost_root', 'nonghost_root', 'exp_root', 'offline',
                                    'available_resolutions', 'available_networks',
                                    'network', 'species', 'resolution', 'start_date', 'end_date', 'experiments', 
                                    'temporal_colocation', 'spatial_colocation', 'report_type', 'report_summary', 
                                    'report_stations', 'report_title', 'report_filename', 
                                    'plot_characteristics_filename']
        }

        # if variable is given by command line, set that value, otherwise set as default value 
        for k, val in var_defaults.items():
            val = kwargs.get(k, val)
            setattr(self.read_instance, k, self.parse_parameter(k, val))

    def parse_parameter(self, key, value):
        """ parse a parameter """

        # get available N CPUs
        if key == 'available_cpus':
            if (MACHINE == 'power') or (MACHINE == 'mn4'):
                bash_command = 'squeue -h -o "%C"'
                process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
                output, _ = process.communicate()
                return int(re.findall(r'\d+', str(output))[0])
            else:
                return int(os.cpu_count())

        elif key == 'cartopy_data_dir':
            # set cartopy data directory (needed on CTE-POWER/MN4/N3 as has no external
            # internet connection)

            if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3v2'):
                return '/gpfs/projects/bsc32/software/rhel/7.5/ppc64le/POWER9/software/Cartopy/0.17.0-foss-2018b-Python-3.7.0/lib/python3.7/site-packages/Cartopy-0.17.0-py3.7-linux-ppc64le.egg/cartopy/data'
            # on all other machines pull from internet

        elif key == 'n_cpus':
            # Define number of CPUs to process on (leave empty to automatically
            # utilise all available CPUs) NOTE: if this value is set higher than the
            # actual number of CPUs available, then the max number of CPUs is used.

            if (value == '') or (int(value) > self.read_instance.available_cpus):
                return self.read_instance.available_cpus

        elif key == 'ghost_root':
            # Define GHOST observational root data directory (if undefined it is
            # automatically taken from the BSC machine the tool is ran on)

            # set default if left undefined
            if value == '':
                # running on CTE-POWER/MN4/N3?
                if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3v2'):
                    return '/gpfs/projects/bsc32/AC_cache/obs/ghost'
                else:
                    # running on workstation?
                    return '/esarchive/obs/ghost'

        elif key == 'nonghost_root':
            # Define non-GHOST observational root data directory (if undefined it is
            # automatically taken from the BSC machine the tool is ran on)

            # set default if left undefined
            if value == '':
                # running on MN4?
                if (MACHINE == 'mn4'):
                    return '/gpfs/projects/bsc32/AC_cache/obs/nonghost'
                else:
                    return '/esarchive/obs'

        elif key == 'exp_root':
            # Define experiment root data directory
            # set experiment root data directory if left undefined
            if value == '':
                # not running on workstation?
                if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3v2'):
                    return '/gpfs/projects/bsc32/AC_cache/recon/exp_interp'
                else:
                    # running on workstation?
                    return '/esarchive/recon/prov_interp'

        elif key == 'ghost_version':
            # parse GHOST version

            # import GHOST standards
            sys.path.insert(1, os.path.join(CURRENT_PATH, 'dependencies/GHOST_standards/{}'.format(value)))
            from GHOST_standards import standard_parameters
            from GHOST_standards import get_standard_metadata
            from GHOST_standards import standard_data_flag_name_to_data_flag_code
            from GHOST_standards import standard_QA_name_to_QA_code

            # modify standard parameter dictionary to have BSC standard parameter names as keys (rather than GHOST)
            self.read_instance.parameter_dictionary = dict()
            for _, param_dict in standard_parameters.items():
                self.read_instance.parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict
            # get standard metadata dictionary
            self.read_instance.standard_metadata = get_standard_metadata({'standard_units': ''})
            # create list of GHOST metadata variables to read
            self.read_instance.ghost_metadata_vars_to_read = [key for key in self.read_instance.standard_metadata.keys() if
                                                              pd.isnull(self.read_instance.standard_metadata[key]['metadata_type']) == False]
            self.read_instance.ghost_metadata_dtype = [(key, self.read_instance.standard_metadata[key]['data_type']) 
                                                       for key in self.read_instance.ghost_metadata_vars_to_read]
            self.read_instance.standard_data_flag_name_to_data_flag_code = standard_data_flag_name_to_data_flag_code
            self.read_instance.standard_QA_name_to_QA_code = standard_QA_name_to_QA_code
            self.read_instance.met_parameters = ['dir10', 'spd10', 't2', 'rh2', 'sst', 'td2', 'pshltr',
                                                 'slp', 'acprec', 'acsnow', 'si', 'cldbot', 'vdist', 'cfracmax']

            return str(value)

        elif key == 'network':
            # parse network

            if isinstance(value, str):
                # throw error if network is empty str
                if value.strip() == '':
                    error = 'Error: "network" field is empty in .conf file'
                    sys.exit(error)
                # parse multiple networks
                elif ',' in value:
                    return [network.strip() for network in value.split(',')]
                else:
                    return [value.strip()]

        elif key == 'species':
            # parse species
            if isinstance(value, str):
                # throw error if species is empty str
                if value.strip() == '':
                    error = 'Error: "species" field is empty in .conf file'
                    sys.exit(error)
                # parse multiple species
                elif ',' in value:
                    return [speci.strip() for speci in value.split(',')]
                else:
                    return [value.strip()]

        elif key == 'qa':
            # parse qa

            # set default qa codes
            standard_qa_names = json.load(open("providentia/conf/default_qa.json"))['standard']
            met_qa_names = json.load(open("providentia/conf/default_qa.json"))['met']
            self.read_instance.default_qa_standard = [self.read_instance.standard_QA_name_to_QA_code[qa_name] 
                                                      for qa_name in standard_qa_names]
            self.read_instance.default_qa_met = [self.read_instance.standard_QA_name_to_QA_code[qa_name] 
                                                 for qa_name in met_qa_names]

            # if not None then set QA by that given
            if value is not None:
                # if conf has only 1 QA
                if isinstance(value, int):
                    return [value]
                # empty string
                elif value == "":
                    return []
                # if the QAs are written with their names
                elif isinstance(value, str):
                    return sorted([self.read_instance.standard_QA_name_to_QA_code[q.strip()] for q in value.split(",")])
                # list of integer codes
                else:
                    return sorted(list(value))
            # otherwise, set default QA per species (set later)
            else:
                # set qa to be empty dict (to be later filled)
                return {}

        elif key == 'flags':
            # parse flags

            # if not None then set flags by that given
            if value is not None:
                # if conf has only one flag
                if isinstance(value, int):
                    return [value]
                # empty string
                elif value == "":
                    return []
                # if the flags are written with their names
                elif isinstance(value, str):
                    return [self.read_instance.standard_data_flag_name_to_data_flag_code[f.strip()]
                            for f in value.split(",")]
                # list of integer codes
                else:
                    return sorted(list(value))
            # otherwise, set default (empty list)
            else:
                return []

        elif key == 'experiments':
            # parse experiments

            if isinstance(value, str):
                # empty string
                if value == "":
                    return {}
                # split experiments
                else:
                    # have alternative experiment names for the legend, then parse them?
                    if ('(' in value) & (')' in value):
                        exps = [exp.strip() for exp in value.split('(')[0].strip().split(",")]
                        exps_legend = [exp_legend.strip() for exp_legend in value.split('(')[1].split(')')[0].strip().split(",")]
                    # otherwise set legend names as given experiment names in full
                    else: 
                        exps = [exp.strip() for exp in value.split(",")]
                        exps_legend = copy.deepcopy(exps)
                    return {exp:exp_legend for exp,exp_legend in zip(exps,exps_legend)}

        elif key == 'map_extent':
            # parse map extent

            # map_extent empty?
            if not value:
                # if dashboard, if map extent not defined then fix to global default extent
                if not self.read_instance.offline:
                    return [-180, 180, -90, 90]
            # otherwise parse it
            else:
                if isinstance(value, str):
                    return [float(c.strip()) for c in value.split(',')]

        elif key == 'filter_species':
            # parse filter species

            # per networkspecies to filter by, save in dict as networkspecies:[lower_limit, upper_limit] 
            if isinstance(value, str):
                # strip all whitespace
                value_strip = "".join(value.split())
                # return empty dict if empty str
                if value_strip == '':
                    return {}
                else:
                    # split per networkspecies
                    networkspecies_split = value_strip.split('),')
                    # iterate through networkspecies, saving list of limits per networkspecies
                    filter_networkspecies_dict = {}
                    for networkspeci_split in networkspecies_split:
                        networkspeci_split_2 = networkspeci_split.split('(')
                        # get networkspeci
                        networkspeci = networkspeci_split_2[0].replace(':','|')
                        # get lower and upper limits
                        networkspeci_split_3 = networkspeci_split_2[1].split(',')
                        lower_limit = float(networkspeci_split_3[0])
                        upper_limit = float(networkspeci_split_3[1].replace(')',''))
                        # save limits per networkspecies
                        filter_networkspecies_dict[networkspeci] = [lower_limit, upper_limit]
                    return filter_networkspecies_dict

        elif key == 'lower_bound':
            #parse lower_bound
            
            # if not None then set lower_bound by that given
            # make sure it is a list of values
            if value is not None:
                if value == "":
                    return []
                elif isinstance(value, str):
                    return [np.float32(c.strip()) for c in value.split(',')]
                elif (isinstance(value, int)) or (isinstance(value, float)):
                    return [value]
                #otherwise must be a already a list of values
                else:
                    return value
            # lower_bound empty?
            # then set lower bound using GHOST extreme lower limit for all species in memory (set later)
            else:
                return {}

        elif key == 'upper_bound':
            # parse upper bound

            # if not None then set upper_bound by that given
            # make sure it is a list of values
            if value is not None:
                if value == "":
                    return []
                elif isinstance(value, str):
                    return [np.float32(c.strip()) for c in value.split(',')]
                elif (isinstance(value, int)) or (isinstance(value, float)):
                    return [value]
                #otherwise must be a already a list of values
                else:
                    return value
            # upper_bound empty?
            # then set upper bound using GHOST extreme upper limit for all species in memory (set later)
            else:
                return {}

        # if no special parsing treatment for variable, simply return value
        return value

    def check_validity(self):
        """check validity of set variables after parsing"""

        # check have network, otherwise throw error
        if not hasattr(self.read_instance, 'network'):
            error = 'Error: "network" field must be defined in .conf file'
            sys.exit(error)

        # check have species, otherwise throw error
        if not hasattr(self.read_instance, 'species'):
            error = 'Error: "species" field must be defined in .conf file'
            sys.exit(error)

        # if number of networks and species is not the same,
        # and len of one of network or species == 1,
        # then duplicate respestive network/species
        if len(self.read_instance.network) != len(self.read_instance.species):

            # 1 network?
            if len(self.read_instance.network) == 1:
                # duplicate network to match species len
                self.read_instance.network = self.read_instance.network * len(self.read_instance.species)

            # 1 species?
            elif len(self.read_instance.species) == 1:
                # duplicate species to match network len
                self.read_instance.species = self.read_instance.species * len(self.read_instance.network)

            # otherwise throw error
            else:
                error = 'Error: The number of networks and species is not the same.'
                sys.exit(error)

        # throw error if one of networks are non all GHOST or non-GHOST
        for network_ii, network in enumerate(self.read_instance.network):
            if network_ii == 0:
                previous_is_ghost = check_for_ghost(network)
            else:
                is_ghost = check_for_ghost(network)
                if is_ghost != previous_is_ghost:
                    error = 'Error: Networks must be all GHOST or non-GHOST'
                    sys.exit(error)
                previous_is_ghost = is_ghost

        # resolution
        if hasattr(self.read_instance, 'resolution'):
            # throw error if species if empty str
            if self.read_instance.resolution.strip() == '':
                error = 'Error: "resolution" field is empty in .conf file'
                sys.exit(error)
        # throw error if resolution is not defined
        else:
            error = 'Error: "resolution" field must be defined in .conf file'
            sys.exit(error)

        # start_date
        if hasattr(self.read_instance, 'start_date'):
            # throw error if start_date if empty str
            if str(self.read_instance.start_date).strip() == '':
                error = 'Error: "start_date" field is empty in .conf file'
                sys.exit(error)
        # throw error if start_date is not defined
        else:
            error = 'Error: "start_date" field must be defined in .conf file'
            sys.exit(error)

        # end_date
        if hasattr(self.read_instance, 'end_date'):
            # throw error if start_date if empty str
            if str(self.read_instance.end_date).strip() == '':
                error = 'Error: "end_date" field is empty in .conf file'
                sys.exit(error)
        # throw error if end_date is not defined
        else:
            error = 'Error: "end_date" field must be defined in .conf file'
            sys.exit(error)

        # if filter_species is active, and spatial_colocation is not active, then cannot filter by species
        # set filter_species to empty dict and advise user of this
        if (self.read_instance.filter_species) and (not self.read_instance.spatial_colocation):
            self.read_instance.filter_species = {}
            print('Warning: "spatial_colocation" must be set to True if wanting to use "filter_species" option.')

        # map to multiple species if have * wildcard
        # also duplicate out associated network
        # remove any species for which there exists no data
        new_species = copy.deepcopy(self.read_instance.species)
        for speci_ii, speci in enumerate(self.read_instance.species): 
            if '*' in speci:
                mapped_species = multi_species_mapping(speci)
                del new_species[speci_ii]
                new_species[speci_ii:speci_ii] = mapped_species
                network_to_duplicate = self.read_instance.network[speci_ii]
                del self.read_instance.network[speci_ii]
                self.read_instance.network[speci_ii:speci_ii] = [network_to_duplicate]*len(mapped_species)
        self.read_instance.species = copy.deepcopy(new_species)

        # create variable for all unique species (plus filter species)
        filter_species = []
        species_plus_filter_species = copy.deepcopy(self.read_instance.species)
        if self.read_instance.filter_species:
            for networkspeci in self.read_instance.filter_species:
                speci = networkspeci.split('|')[1]
                if speci not in self.read_instance.species:
                    filter_species.append(speci)
                    species_plus_filter_species.append(speci)
                    
        # set lower_bound and upper_bound as dicts with limits per species (including filter species)
        # if type is dict then set bound using GHOST extreme limits per species

        # lower_bound
        # type is dict, then set as default limit per species using GHOST limits
        if isinstance(self.read_instance.lower_bound, dict):
            self.read_instance.lower_bound = {speci:np.float32(self.read_instance.parameter_dictionary[speci]['extreme_lower_limit']) for speci in species_plus_filter_species}
        # otherwise set list values to dict, saving limits per species
        # if have just 1 limit apply for all read species, but if have multiple, set limits per species
        # throw error if have multiple lower bounds, but not equal to number of species to read  
        else:      
            lower_bound_dict = {}
            if len(self.read_instance.lower_bound) == 1:
                for speci in species_plus_filter_species:
                    lower_bound_dict[speci] = self.read_instance.lower_bound[0]
            elif len(self.read_instance.lower_bound) > 1:
                if len(self.read_instance.species) != len(self.read_instance.lower_bound):
                    error = 'Error: "lower_bound" variable must be same length as number of species read.'
                    sys.exit(error)
                else:
                    for speci_ii, speci in enumerate(self.read_instance.species):
                        lower_bound_dict[speci] = value[speci_ii] 
                    # add filter_species (using GHOST limits)
                    for speci in filter_species:
                        lower_bound_dict[speci] = np.float32(self.read_instance.parameter_dictionary[speci]['extreme_lower_limit'])
            self.read_instance.lower_bound = lower_bound_dict

        # upper_bound
        # type is dict, then set as default limit per species using GHOST limits
        if isinstance(self.read_instance.upper_bound, dict):
            self.read_instance.upper_bound = {speci:np.float32(self.read_instance.parameter_dictionary[speci]['extreme_upper_limit']) for speci in species_plus_filter_species}
        # otherwise set list values to dict, saving limits per species
        # if have just 1 limit apply for all read species, but if have multiple, set limits per species
        # throw error if have multiple upper bounds, but not equal to number of species to read  
        else:      
            upper_bound_dict = {}
            if len(self.read_instance.upper_bound) == 1:
                for speci in species_plus_filter_species:
                    upper_bound_dict[speci] = self.read_instance.upper_bound[0]
            elif len(self.read_instance.upper_bound) > 1:
                if len(self.read_instance.species) != len(self.read_instance.upper_bound):
                    error = 'Error: "upper_bound" variable must be same length as number of species read.'
                    sys.exit(error)
                else:
                    for speci_ii, speci in enumerate(self.read_instance.species):
                        upper_bound_dict[speci] = value[speci_ii] 
                    # add filter_species (using GHOST limits)
                    for speci in filter_species:
                        upper_bound_dict[speci] = np.float32(self.read_instance.parameter_dictionary[speci]['extreme_upper_limit'])
            self.read_instance.upper_bound = upper_bound_dict

        # create a variable to set qa per species (including filter species)
        if isinstance(self.read_instance.qa, dict):
            self.read_instance.qa_per_species = {speci:get_default_qa(self.read_instance, speci) for speci in species_plus_filter_species}
            # set qa to be first of qa per species pairs
            self.read_instance.qa = self.read_instance.qa_per_species[list(self.read_instance.qa_per_species.keys())[0]]
        else:
            self.read_instance.qa_per_species = {speci:self.read_instance.qa for speci in species_plus_filter_species}

        # if are using dashboard then just take first network/species pair, as multivar not supported yet
        if (len(self.read_instance.network) > 1) & (len(self.read_instance.species) > 1) & (not self.read_instance.offline):
            self.read_instance.network = [self.read_instance.network[0]]
            self.read_instance.species = [self.read_instance.species[0]]
            print('Warning: Multiple networks/species not supported for dashboard. First network / species taken.')

def read_conf(fpath=None):
    """Read configuration"""

    res = {}
    res_sub = {}
    config = {}
    all_sections = []
    all_sections_modified = []
    all_sections_commented = []
    repeated_subsections = []
    repeated_subsections_modified = {}
    subsections = []
    subsections_modified = []
    parent_sections = []
    filenames = []

    # get section names (e.g. [SECTIONA], [[Spain]]) and modified names (e.g. SECTIONA, SECTIONA-Spain)
    with open(fpath) as file:
        for line in file:
            if '[' in line and ']' in line and '[[' not in line and ']]' not in line:
                section = line.strip()
                #if first character is comment do not parse section
                if section[0] != '#':
                    section_modified = section.split('[')[1].split(']')[0]
                    if section_modified not in all_sections_modified:
                        parent_sections.append(section_modified)
                        all_sections_modified.append(section_modified)
                    else:
                        error = 'Error: It is not possible to have two sections with the same name.'
                        sys.exit(error)
            elif '[[' in line and ']]' in line:
                subsection = line.strip()
                #if first character is comment do not parse subsection
                if subsection[0] != '#':
                    subsection_modified = section_modified + '|' + line.split('[[')[1].split(']]')[0]
                    subsections.append(subsection)
                    subsections_modified.append(subsection_modified)
                    all_sections_modified.append(subsection_modified)

            if '[' in line and ']' in line:
                #if first character is comment then add section to list to avoid parsing
                if line.strip()[0] == '#':
                    all_sections_commented.append(line.strip())
                else:
                    all_sections.append(line.strip())

    # get repeated elements
    repetition_counts = {section:subsections.count(section) for section in subsections}
    for section, counts in repetition_counts.items():
        if counts > 1:
            repeated_subsections.append(section)
            repeated_subsections_modified[section] = [x for x in all_sections_modified 
                                                        if section.split('[[')[1].split(']]')[0] in x]

    # get attributes for each section and store in dict
    for (i, section), section_modified in zip(enumerate(all_sections), all_sections_modified):
        
        repetition = 0
        copy = False
        config[section_modified] = {}
        
        with open(fpath) as file:
            for line in file:
                line_strip = line.strip()
                # parsing all but last section 
                if section_modified != all_sections_modified[-1]:
                    # start of relevant section 
                    if line_strip == all_sections[i]:
                        if line_strip in repeated_subsections:
                            position = repeated_subsections_modified[section].index(section_modified)
                            if position == repetition:
                                copy = True
                            else:
                                copy = False
                            repetition += 1
                        else:
                            copy = True
                        continue
                    # start of next section, or commented section
                    elif (line_strip == all_sections[i+1]) or (line_strip in all_sections_commented):
                        copy = False
                        continue
                
                # parsing last section
                else:
                    # start of relevant section ?
                    if line_strip == all_sections[-1]:
                        if line_strip in repeated_subsections:
                            position = repeated_subsections_modified[section].index(section_modified)
                            if position == repetition:
                                copy = True
                            else:
                                copy = False
                            repetition += 1
                        else:
                            copy = True
                        continue

                    #start of commented section
                    elif line_strip in all_sections_commented:
                        copy = False
                        continue
                
                # set section attributes
                if copy:
                    if line_strip != '':
                        key = line_strip.split('=')[0].strip()
                        value = line_strip.split('=')[1].strip()
                        config[section_modified][key] = value

    # add section attributes to subsection if do not exist there (e.g. add SECTIONA values to SECTIONA-Spain)
    for section_modified in all_sections_modified:
        
        # determine if subsection or not
        if '|' in section_modified:
            is_subsection = True
            par_section = section_modified.split('|')[0]
            # add attributes from parent section
            for par_k, par_val in config[par_section].items():
                # if first character of key is comment character (#), do not parse this attribute
                if par_k[0] == '#':
                    continue
                try:
                    res_sub[par_k] = eval(par_val)
                except:
                    res_sub[par_k] = par_val
        else:
            is_subsection = False

        # add attributes from current section/subsection
        for k, val in config[section_modified].items():
            if not is_subsection:
                # store filename
                if k == 'report_filename':
                    filenames.append(val)

            # if first character of key is comment character (#), do not parse this attribute
            if k[0] == '#':
                continue
            # overwrite attributes from current subsection
            try:
                res_sub[k] = eval(val)
            except:
                res_sub[k] = val

        # store pairs into res variable
        res[section_modified] = res_sub

        # reset res variable
        res_sub = {}

    return res, all_sections_modified, parent_sections, subsections_modified, filenames
   
def write_conf(section, subsection, fpath, opts):
    """Write configurations on file. """

    config = configparser.RawConfigParser()

    # update configuration
    for section, section_name in zip(['section', 'subsection'], [section, subsection]):
        if opts[section]:
            config.add_section(section_name)
            for item in opts[section]:
                val = opts[section][item]
                config.set(section_name, item, val)

    # write configuration
    with open(fpath, 'w') as configfile:
        config.write(configfile)

def split_options(conf_string, separator="||"):
    """For the options in the configuration that define the keep and remove
    options. Returns the values in two lists, the keeps and removes"""
    
    keeps, removes = [], []

    if separator not in conf_string:
        if ("keep:" in conf_string) and ("remove:" not in conf_string):
            keep_start = conf_string.find("keep:")
            keeps = conf_string[keep_start+5:]
            keeps = keeps.split(",")
            keeps = [k.strip() for k in keeps]
        elif ("keep:" not in conf_string) and ("remove:" in conf_string):
            remove_start = conf_string.find("remove:")
            removes = conf_string[remove_start+7:]
            removes = removes.split(",")
            removes = [r.strip() for r in removes]
        elif ("keep:" in conf_string) and ("remove:" in conf_string):
            print('Warning: In order to define the keep and remove options, these must be separated by ||.')
    else:
        if "keep:" in conf_string:
            keep_start, keep_end = conf_string.find("keep:"), conf_string.find(separator)
            keeps = conf_string[keep_start+5:keep_end]
            keeps = keeps.split(",")
            keeps = [k.strip() for k in keeps]
        if "remove:" in conf_string:
            remove_start = conf_string.find("remove:")
            removes = conf_string[remove_start+7:]
            removes = removes.split(",")
            removes = [r.strip() for r in removes]
    
    return keeps, removes