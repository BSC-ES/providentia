""" Providentia config module """
# -*- coding: utf-8 -*-
#
# Copyright (c) 2016 Barcelona Supercomputing Center
# @license: https://www.gnu.org/licenses/gpl-3.0.html
# @author: see AUTHORS file

import os
import sys
import configparser
import logging

from configargparse import ArgumentParser
import providentia
from . import prov_exceptions
import numpy as np

logging.basicConfig(level=logging.WARNING)
log = logging.getLogger(__name__)

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

class ProvArgumentParser(object):
    """ Argument Parser """

    def __init__(self):
        """
        Initialization of the arguments the parser can handle
        """

        try:
            self.parser = ArgumentParser(description='Main parser for Providentia.')
            self.parser.add_argument('-V', '--version', action='version',
                                     version=providentia.__version__,
                                     help="returns Providentia version number and exit")
            self.parser.add_argument('--config', #is_config_file=True,
                                     dest="config",
                                     help='specifies the config file to read'
                                     ) #required=False)
            self.parser.add_argument('--config_dir', #is_config_file=True,
                                     dest="config_dir",
                                     help='specifies the configuration directory where config files are'
                                     ) #required=False)
            self.parser.add_argument("--section",
                                     dest="section",
                                     help="config file section to read")
            # main options
            self.parser.add_argument("--ghost_version",
                                     dest="ghost_version",
                                     help="set GHOST version data to work with")
            self.parser.add_argument("--cartopy_data_dir",
                                     dest="cartopy_data_dir",
                                     help="set cartopy data directory")
            self.parser.add_argument("--n_cpus",
                                     dest="n_cpus",
                                     help="Define number of CPUs to process on")
            self.parser.add_argument("--obs_root",
                                     dest="obs_root",
                                     help="directory where is/are to be stored observations")
            self.parser.add_argument("--nonghost_root",
                                     dest="nonghost_root",
                                     help="directory where is/are to stored nonghost observations")
            self.parser.add_argument("--exp_root",
                                     dest="exp_root",
                                     help="set experiment root data directory")
            self.parser.add_argument("--offline",
                                     dest="offline",
                                     default=False,
                                     action='store_true',
                                     help="run Providentia offline")
#            self.parser.add_argument("--debug",
#                                     dest="debug",
#                                     help="debug (default=False)")
            self.parser.add_argument("--map_coastline_resolution",
                                     dest="map_coastline_resolution",
                                     help="define coastlines resolution")
            self.parser.add_argument("--available_networks",
                                     dest="available_networks",
                                     help="define available networks (default=['EBAS', 'EEA_AQ_eReporting'])")
            self.parser.add_argument("--network",
                                     dest="network",
                                     help="define network to load (e.g. 'EBAS', 'EEA_AQ_eReporting'")
            self.parser.add_argument("--resolution",
                                     dest="resolution",
                                     help="define data resolution (e.g. 'hourly', '3hourly', 'daily'")
            self.parser.add_argument("--matrix",
                                     dest="matrix",
                                     help="define species matrix (e.g. 'gas', 'aerosol'")
            self.parser.add_argument("--species",
                                     dest="species",
                                     help="define species to load (e.g. 'sconco3', 'pm10'")
            self.parser.add_argument("--start_date",
                                     dest="start_date",
                                     help="define start date in format as 20160101")
            self.parser.add_argument("--end_date",
                                     dest="end_date",
                                     help="define end date in format as 20170101")

        except Exception as error:
            log.error('Unhandled exception on Providentia: %s' % error, exc_info=True)

    #-----------------------------------------------------------------------
    # Parse arguments and preprocess
    #-----------------------------------------------------------------------
    def parse_args(self, args=None):
        """
        Parse arguments given to an executable
        :param args:
        """
        try:
            return self.parser.parse_args(args)
        except Exception as error:
            print(error)
            raise prov_exceptions.ProvArgumentParserException

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
                    if k == 'filename':
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
