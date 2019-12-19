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
from providentia import prov_exceptions


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
            self.parser.add_argument("--exp_root",
                                     dest="exp_root",
                                     help="set experiment root data directory")
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
            self.parser.add_argument("--time_series_markersize",
                                     dest="time_series_markersize",
                                     help="define marker sizes on time series plot")
            self.parser.add_argument("--temp_aggregated_markersize",
                                     dest="temp_aggregated_markersize",
                                     help="define marker sizes on temporally aggregated plot")
            self.parser.add_argument("--temp_agg_expbias_markersize",
                                     dest="temp_agg_expbias_markersize",
                                     help="define marker sizes on temporally aggregated experiment bias plot")

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
            #return self.parser.parse_args(args)
            return self.parser.parse_args(args)
        except Exception as error:
            print(error)
            raise prov_exceptions.ProvArgumentParserException


def read_conf(section=None, fpath=None):
    """ Read configuration """

    config = configparser.RawConfigParser()
    config.read(fpath)
    if section is None:
        return config.sections()

    res = {}
    for k, val in config.items(section):
        try:
            res[k] = eval(val)
        except:
            res[k] = val
    return res


def write_conf(section, fpath, opts):
    """ Write configurations on file. """

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
    with open(fpath, 'wb') as configfile:
        config.write(configfile)
