""" Providentia Configuration Module """

import configparser
import os
import sys
import re
import subprocess

MACHINE = os.environ.get('BSC_MACHINE', '')

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

def parse_path(dir, f):
    # print "Opening data file", f
    if os.path.isabs(f):
        return f
    else:
        return os.path.join(dir, f)

class ProvConfiguration(object):
    """ Configuration parameters definitions """

    def __init__(self, **kwargs):
        self.config_dir = kwargs.get('config_dir',
                                     os.path.join(os.environ['HOME'], '.providentia'))
        self.ghost_version = kwargs.get('ghost_version', '1.4')
        self.cartopy_data_dir = kwargs.get('cartopy_data_dir', '')
        self.available_cpus = kwargs.get('available_cpus', '')
        self.n_cpus = kwargs.get('n_cpus', '')
        self.ghost_root = kwargs.get('ghost_root', '')
        self.nonghost_root = kwargs.get('nonghost_root', '')
        self.exp_root = kwargs.get('exp_root', '')
        self.offline = kwargs.get('offline', '')
        self.available_resolutions =\
            kwargs.get('available_resolutions',
                       ['hourly', '3hourly', '6hourly', 'hourly_instantaneous',
                       '3hourly_instantaneous', '6hourly_instantaneous',
                       'daily', 'monthly'])
        self.available_networks =\
            kwargs.get('available_networks',
                       ['AERONET_v3_lev1.5','AERONET_v3_lev2.0','CANADA_NAPS','CAPMoN','CHILE_SINCA',
                        'EANET','EBAS','EEA_AIRBASE','EEA_AQ_eReporting','JAPAN_NIES','MEXICO_CDMX',
                        'MITECO','NOAA_ISD','NOAA_ISD_EU','NOAA_ISD_IP','NOAA_ISD_NA',
                        'SEARCH','UK_AIR','US_EPA_AQS','US_EPA_CASTNET','US_NADP_AMNet','US_NADP_AMoN','WMO_WDCGG'])
        self.network = kwargs.get('network', '')
        self.species = kwargs.get('species', '')
        self.resolution = kwargs.get('resolution', '')
        self.start_date = kwargs.get('start_date', '')
        self.end_date = kwargs.get('end_date', '')
        self.position_1 = kwargs.get('position_1', 'map')
        self.position_2 = kwargs.get('position_2', 'timeseries')
        self.position_3 = kwargs.get('position_3', 'metadata')
        self.position_4 = kwargs.get('position_4', 'distribution')
        self.position_5 = kwargs.get('position_5', 'periodic')
        self.experiments = kwargs.get('experiments', '')
        self.temporal_colocation = kwargs.get('temporal_colocation', False)
        self.spatial_colocation = kwargs.get('spatial_colocation', True)
        self.filter_species = kwargs.get('filter_species', '')
        self.report_type = kwargs.get('report_type', 'standard')
        self.report_summary = kwargs.get('report_summary', True)
        self.report_stations = kwargs.get('report_stations', False)
        self.report_title = kwargs.get('report_title ', 'Report')
        self.report_filename = kwargs.get('report_filename', 'PROVIDENTIA_Report')
        self.map_extent = kwargs.get('map_extent', '-180, 180, -90, 90')
        self.plot_characteristics_filename = kwargs.get('plot_characteristics_filename', '')
        self.fixed_section_vars =  ['network', 'species', 'resolution', 'start_date', 'end_date', 'experiments', 
                                    'spatial_colocation', 'report_type', 'report_summary', 'report_stations',
                                    'report_title', 'report_filename', 'plot_characteristics_filename']

    def __setattr__(self, key, value):
        super(ProvConfiguration, self).__setattr__(key, self.parse_parameter(key, value))

    def parse_parameter(self, key, value):
        """ parse parameters """

        # get available N CPUs
        if key == 'available_cpus':
            if (MACHINE == 'power') or (MACHINE == 'mn4'):
                bash_command = 'squeue -h -o "%C"'
                process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
                output, _ = process.communicate()
                return int(re.findall(r'\d+', str(output))[0])

            value = int(os.cpu_count())

        elif key == 'cartopy_data_dir':
            # set cartopy data directory (needed on CTE-POWER/MN4/N3 as has no external
            # internet connection)

            if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3v2'):
                value = '/gpfs/projects/bsc32/software/rhel/7.5/ppc64le/POWER9/software/Cartopy/0.17.0-foss-2018b-Python-3.7.0/lib/python3.7/site-packages/Cartopy-0.17.0-py3.7-linux-ppc64le.egg/cartopy/data'
            # on all other machines pull from internet

        elif key == 'n_cpus':
            # Define number of CPUs to process on (leave empty to automatically
            # utilise all available CPUs) NOTE: if this value is set higher than the
            # actual number of CPUs available, then the max number of CPUs is used.

            if (value == '') or (int(value) > self.available_cpus):
                value = self.available_cpus

        elif key == 'ghost_root':
            # Define GHOST observational root data directory (if undefined it is
            # automatically taken from the BSC machine the tool is ran on)

            # set default if left undefined
            if value == '':
                # running on CTE-POWER/MN4/N3?
                if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3v2'):
                    value = '/gpfs/projects/bsc32/AC_cache/obs/ghost'
                else:
                    # running on workstation?
                    value = '/esarchive/obs/ghost'

        elif key == 'nonghost_root':
            # Define non-GHOST observational root data directory (if undefined it is
            # automatically taken from the BSC machine the tool is ran on)

            # set default if left undefined
            if value == '':
                # running on MN4?
                if (MACHINE == 'mn4'):
                    value = '/gpfs/projects/bsc32/AC_cache/obs/nonghost'
                else:
                    value = '/esarchive/obs'

        elif key == 'exp_root':
            # Define experiment root data directory
            # set experiment root data directory if left undefined
            if value == '':
                # not running on workstation?
                if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3v2'):
                    value = '/gpfs/projects/bsc32/AC_cache/recon/exp_interp'
                else:
                    # running on workstation?
                    value = '/esarchive/recon/prov_interp'

        return value

def read_conf(fpath=None):
    """Read configuration"""

    dconf_path = (os.path.join(CURRENT_PATH, 'conf/default.conf'))

    res = {}
    res_sub = {}

    if fpath == dconf_path:
        config = configparser.RawConfigParser(empty_lines_in_values=False)
        config.read(fpath) 

        for k, val in config.items('default'):
            try:
                res_sub[k] = eval(val)
            except:
                res_sub[k] = val

        res['default'] = res_sub
        all_sections_modified, parent_sections, subsections_modified, filenames = None, None, None, None
        
    else:
        config = {}
        all_sections = []
        all_sections_modified = []
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
                    section_modified = par_section = section.split('[')[1].split(']')[0]
                    if section_modified not in all_sections_modified:
                        parent_sections.append(section_modified)
                        all_sections_modified.append(section_modified)
                    else:
                        print('Error: It is not possible to have two sections with the same name.')
                        sys.exit()
                elif '[[' in line and ']]' in line:
                    subsection = line.strip()
                    subsection_modified = par_section + '|' + line.split('[[')[1].split(']]')[0]
                    subsections.append(subsection)
                    subsections_modified.append(subsection_modified)
                    all_sections_modified.append(subsection_modified)

                if '[' in line and ']' in line:
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
                    if section_modified != all_sections_modified[-1]:
                        if line.strip() == all_sections[i]:
                            if line.strip() in repeated_subsections:
                                position = repeated_subsections_modified[section].index(section_modified)
                                if position == repetition:
                                    copy = True
                                else:
                                    copy = False
                                repetition += 1
                            else:
                                copy = True
                            continue
                        elif line.strip() == all_sections[i+1]:
                            copy = False
                            continue
                    else:
                        if line.strip() == all_sections[-1]:
                            if line.strip() in repeated_subsections:
                                position = repeated_subsections_modified[section].index(section_modified)
                                if position == repetition:
                                    copy = True
                                else:
                                    copy = False
                                repetition += 1
                            else:
                                copy = True
                            continue
        
                    if copy:
                        if line.strip() != '':
                            key = line.split('=')[0].strip()
                            value = line.split('=')[1].strip()
                            config[section_modified][key] = value

        # regenerate sections (e.g. add SECTIONA values to SECTIONA-Spain)
        par_section = all_sections_modified[0]
        for section_modified in all_sections_modified:

            # get parent section
            if '|' in section_modified:
                is_subsection = True
            else:
                is_subsection = False
                par_section = section_modified

            # store key-values pairs
            for k, val in config[section_modified].items():
                if is_subsection:
                    # store pairs from parent section
                    for par_k, par_val in config[par_section].items():
                        try:
                            res_sub[par_k] = eval(par_val)
                        except:
                            res_sub[par_k] = par_val
                else:
                    # store filename
                    if k == 'report_filename':
                        filenames.append(val)
                        
                # store pairs from current section
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
            print('Warning!!!')
            print('In order to define the keep and remove options, these must be separated by ||.')
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