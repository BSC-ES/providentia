""" Providentia config module """
# -*- coding: utf-8 -*-
#
# Copyright (c) 2016 Barcelona Supercomputing Center
# @license: https://www.gnu.org/licenses/gpl-3.0.html
# @author: see AUTHORS file

import os
import configparser
import logging

from configargparse import ArgumentParser
import providentia
from . import prov_exceptions

logging.basicConfig(level=logging.WARNING)
log = logging.getLogger(__name__)


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
#                                dest="debug",
#                                help="debug (default=False)")
            self.parser.add_argument("--sequential_colourmap",
                                     dest="sequential_colourmap",
                                     help="perceptually uniform sequential colourmap (default='viridis')")
            self.parser.add_argument("--sequential_colourmap_warm",
                                     dest="sequential_colourmap_warm",
                                     help="warm sequential colourmap (default='Reds')")
            self.parser.add_argument("--diverging_colourmap",
                                     dest="diverging_colourmap",
                                     help="diverging colourmap (default='bwr')")
            self.parser.add_argument("--unsel_station_markersize",
                                     dest="unsel_station_markersize",
                                     help="define marker sizes for unselected stations on map")
            self.parser.add_argument("--sel_station_markersize",
                                     dest="sel_station_markersize",
                                     help="define marker sizes for selected stations on map")
            self.parser.add_argument("--legend_markersize",
                                     dest="legend_markersize",
                                     help="define marker sizes on legend")
            self.parser.add_argument("--timeseries_markersize",
                                     dest="timeseries_markersize",
                                     help="define marker sizes on time series plot")
            self.parser.add_argument("--timeseries_n_slices",
                                     dest="timeseries_n_slices",
                                     help="define number of slices in data labels")
            self.parser.add_argument("--timeseries_last_step",
                                     dest="timeseries_last_step",
                                     help="define if data labels will include end date")
            self.parser.add_argument("--temp_aggregated_markersize",
                                     dest="temp_aggregated_markersize",
                                     help="define marker sizes on temporally aggregated plot")
            self.parser.add_argument("--temp_agg_expbias_markersize",
                                     dest="temp_agg_expbias_markersize",
                                     help="define marker sizes on temporally aggregated experiment bias plot")
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


def read_conf(section=None, fpath=None):
    """Read configuration"""

    config = configparser.RawConfigParser()
    config.read(fpath)
    #if no section defined, but just 1 section in file then set that as section   
    if (section is None) and (len(config.sections()) == 1):   
        section = config.sections()[0]
    #if section is undefined then cannot read
    if section is None:
        print('*** WARNING!!! CANNOT LOAD CONFIGURATION FILE AS NO SECTION DEFINED.')
        return None

    #convert numeric information appropriate types
    res = {}
    for k, val in config.items(section):
        try:
            res[k] = eval(val)
        except:
            res[k] = val
    return res


def read_offline_conf(fpath=None):
    """Read configuration files when running Providentia
    offline. When running offline, having a 'DEFAULT' section
    is mandatory. The 'DEFAULT' section contains options that
    are applied across all remaining sections."""

    config = configparser.RawConfigParser()
    config.read(fpath)

    # check if Default section has options
    if not config.defaults():
        return None

    #convert numeric information to appropriate types
    defaults = {}
    for k, val in config.items(config.default_section):
        try:
            defaults[k] = eval(val)
        except:
            defaults[k] = val

    # store remaining sections into dict
    res = {}
    for section in config.sections():
        res[section] = read_conf(section, fpath)
    return defaults, res


def write_conf(section, fpath, opts):
    """Write configurations on file. """

    config = configparser.RawConfigParser()

    # check if file exists
    if os.path.exists(fpath):
        config.read(fpath)

    # check if section exists
    if not config.has_section(section):
        config.add_section(section)

    # update configuration
    for item in opts:
        val = opts[item]
        config.set(section, item, val)

    # write configuration
    with open(fpath, 'w') as configfile:
        config.write(configfile)


def split_options(conf_string, separator="||"):
    """For the options in the configuration that define the keep and remove
    options. Returns the values in two lists, the keeps and removes"""
    keeps, removes = [], []
    if "keep:" in conf_string:
        keep_start, keep_end = conf_string.find("keep:"), conf_string.find(separator)
        keeps = conf_string[keep_start+5:keep_end]
        keeps = keeps.split(",")
        keeps = [k.strip() for k in keeps]

    if "remove:" in conf_string:
        remove_start = conf_string.find("remove:")
        removes = conf_string[remove_start+7:-1]
        removes = removes.split(",")
        removes = [r.strip() for r in removes]

    return keeps, removes
