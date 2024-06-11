""" Providentia argument parser module """
# -*- coding: utf-8 -*-
#
# Copyright (c) 2016 Barcelona Supercomputing Center
# @license: https://www.gnu.org/licenses/gpl-3.0.html
# @author: see AUTHORS file

import logging
import os
import sys

from configargparse import ArgumentParser

import providentia

logging.basicConfig(level=logging.WARNING)
log = logging.getLogger(__name__)

class ProvArgumentParser(object):
    """ Class that handles the argument parser. """

    def __init__(self):
        """ Initialise the arguments the parser can handle. """

        try:
            self.parser = ArgumentParser(description='Main parser for Providentia.')
            self.parser.add_argument('-V', '--version', action='version',
                                     version=providentia.__version__,
                                     help="returns Providentia version number and exit")
            self.parser.add_argument('--debug',
                                     help="runs Providentia in debug mode, just reserving allocation")
            self.parser.add_argument('--interactive',
                                     help="runs Providentia interactive mode on Jupyter notebook")
            self.parser.add_argument('--conf', '--config', 
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
            self.parser.add_argument("--ghost_root",
                                     dest="ghost_root",
                                     help="root directory where GHOST observations are stored")
            self.parser.add_argument("--nonghost_root",
                                     dest="nonghost_root",
                                     help="root directory where non-GHOST observations are stored")
            self.parser.add_argument("--exp_root",
                                     dest="exp_root",
                                     help="set experiment root data directory")
            self.parser.add_argument("--generate_file_tree", '--file_tree', '--gft',
                                     dest="generate_file_tree",
                                     default=False,
                                     action='store_true',
                                     help="boolean switch whether to dynamically create observational filetrees")
            self.parser.add_argument("--offline",
                                     dest="offline",
                                     default=False,
                                     action='store_true',
                                     help="run Providentia offline mode")
            self.parser.add_argument("--download", "--dl",
                                     dest="download",
                                     default=False,
                                     action='store_true',
                                     help="run Providentia download mode")
            self.parser.add_argument("--network",
                                     dest="network",
                                     help="define network to load (e.g. 'EBAS', 'EEA_AQ_eReporting'")
            self.parser.add_argument("--species",
                                     dest="species",
                                     help="define species to load (e.g. 'sconco3', 'pm10'")
            self.parser.add_argument("--resolution",
                                     dest="resolution",
                                     help="define data resolution (e.g. 'hourly', '3hourly', 'daily'")
            self.parser.add_argument("--start_date",
                                     dest="start_date",
                                     help="define start date in format as 20160101")
            self.parser.add_argument("--end_date",
                                     dest="end_date",
                                     help="define end date in format as 20170101")
            self.parser.add_argument("--observations_data_label",
                                     dest="observations_data_label",
                                     help="alias for observations data label")
            self.parser.add_argument("--experiments",
                                     dest="experiments",
                                     help="experiments to read")
            self.parser.add_argument("--temporal_colocation",
                                     dest="temporal_colocation",
                                     help="activate temporal colocation betwen observations and experiments")
            self.parser.add_argument("--spatial_colocation",
                                     dest="spatial_colocation",
                                     help="activate spatial colocation between multiple read species")
            self.parser.add_argument("--filter_species",
                                     dest="filter_species",
                                     help="filter read species by other species within a given data range")
            self.parser.add_argument("--lower_bound",
                                     dest="lower_bound",
                                     help="filter out data below this set lower bound")
            self.parser.add_argument("--upper_bound",
                                     dest="upper_bound",
                                     help="filter out data above this set upper bound")
            self.parser.add_argument("--report_type",
                                     dest="report_type",
                                     help="define ")
            self.parser.add_argument("--report_summary",
                                     dest="report_summary",
                                     help="activate summary plots in offline report")
            self.parser.add_argument("--report_stations",
                                     dest="report_stations",
                                     help="activate station specific plots in offline report")
            self.parser.add_argument("--report_title",
                                     dest="report_title",
                                     help="offline report title")
            self.parser.add_argument("--report_filename",
                                     dest="report_filename",
                                     help="offline report filename")
            self.parser.add_argument("--map_extent",
                                     dest="map_extent",
                                     help="map extent for plots involving any maps")
            self.parser.add_argument("--active_dashboard_plots",
                                     dest="active_dashboard_plots",
                                     help="active plots on dashboard upon launch")                     
            self.parser.add_argument("--plot_characteristics_filename",
                                     dest="plot_characteristics_filename",
                                     help="set filename for plot characteristics")
            self.parser.add_argument("--harmonise_stations",
                                     dest="harmonise_stations",
                                     help="harmonise axes limits across stations for stations report")                     
            self.parser.add_argument("--harmonise_summary",
                                     dest="harmonise_summary",
                                     help="harmonise axes limits across subsections for summary report")
            self.parser.add_argument("--remove_extreme_stations",
                                     dest="remove_extreme_stations",
                                     help="remove extreme stations using defined statistic limits")


        except Exception as error:
            log.error('Unhandled exception on Providentia: %s' % error, exc_info=True)

    #-----------------------------------------------------------------------
    # Parse arguments and preprocess
    #-----------------------------------------------------------------------
    def parse_args(self, args=None):
        """ Parse arguments given to an executable :param args:. """

        try:
            return self.parser.parse_args(args)
        except Exception as error:
            print(error)
