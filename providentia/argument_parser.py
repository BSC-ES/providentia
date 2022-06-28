""" Providentia argument parser module """
# -*- coding: utf-8 -*-
#
# Copyright (c) 2016 Barcelona Supercomputing Center
# @license: https://www.gnu.org/licenses/gpl-3.0.html
# @author: see AUTHORS file

import os
import logging

from configargparse import ArgumentParser
import providentia
from . import exceptions

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
            self.parser.add_argument('--config', 
                                     dest="config",
                                     help='specifies the config file to read') 
            self.parser.add_argument('--config_dir', 
                                     dest="config_dir",
                                     help='specifies the configuration directory where config files are') 
            self.parser.add_argument("--section",
                                     dest="section",
                                     help="config file section to read")
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
            self.parser.add_argument("--available_networks",
                                     dest="available_networks",
                                     help="define available networks (default=['EBAS', 'EEA_AQ_eReporting'])")
            self.parser.add_argument("--network",
                                     dest="selected_network",
                                     help="define network to load (e.g. 'EBAS', 'EEA_AQ_eReporting'")
            self.parser.add_argument("--resolution",
                                     dest="selected_resolution",
                                     help="define data resolution (e.g. 'hourly', '3hourly', 'daily'")
            self.parser.add_argument("--matrix",
                                     dest="selected_matrix",
                                     help="define species matrix (e.g. 'gas', 'aerosol'")
            self.parser.add_argument("--species",
                                     dest="selected_species",
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
            raise exceptions.ProvArgumentParserException