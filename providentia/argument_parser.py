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
            self.parser.add_argument('--clean',
                                     action='store_true',
                                     help='cleans all output files')
            self.parser.add_argument("--section",
                                     dest="section",
                                     help="config file section to read")
            self.parser.add_argument("--ghost_version",
                                     dest="ghost_version",
                                     help="set GHOST version data to work with")
            self.parser.add_argument('--conf', '--config', 
                                     dest="config",
                                     help='specifies the config file to read') 
            self.parser.add_argument('--config_dir', 
                                     dest="config_dir",
                                     help='specifies the configuration directory where config files are') 
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
            self.parser.add_argument("--exp_to_interp_root",
                                     dest="exp_to_interp_root",
                                     help="set experiment to interpolate root data directory")
            self.parser.add_argument("--generate_file_tree", '--gft',
                                     dest="generate_file_tree",
                                     default=False,
                                     action='store_true',
                                     help="Whether to dynamically create observational filetrees")
            self.parser.add_argument("--disable_file_tree", '--dft',
                                     dest="disable_file_tree",
                                     default=False,
                                     action='store_true',
                                     help="Whether to disable dynamically creating observational filetrees")
            self.parser.add_argument("--report", "--reports", "--offline",
                                     dest="report",
                                     default=False,
                                     action='store_true',
                                     help="run Providentia report mode")
            self.parser.add_argument("--notebook", "--nb", "--jupyter",
                                     help="opens a Jupyter Notebook session")
            self.parser.add_argument("--download", "--dl",
                                     dest="download",
                                     default=False,
                                     action='store_true',
                                     help="run Providentia download mode")
            self.parser.add_argument('--interpolation','--interp','--interpolate', 
                                     action='store_true',
                                     help='runs Providentia Interpolation') 
            self.parser.add_argument("--slurm_job_id",
                                     dest="slurm_job_id",
                                     help="id of the interpolation sbatch job")                
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
            self.parser.add_argument("--domain",
                                     dest="domain",
                                     help="domain of the experiment")
            self.parser.add_argument("--ensemble_options",
                                     dest="ensemble_options",
                                     help="ensemble options of the experiment")
            self.parser.add_argument("--dataset",
                            dest="dataset",
                            help="CAMS dataset")
            self.parser.add_argument("--forecast_day",
                                     dest="forecast_day",
                                     help="day of the model forecast to analyse")
            self.parser.add_argument("--forecast",
                                     dest="forecast",
                                     help="indicates if data comes from forecast")
            self.parser.add_argument("--qa",
                                     dest="qa",
                                     help="list of qa flags (numbers or text) to use to filter data")
            self.parser.add_argument("--flags",
                                     dest="flags",
                                     help="list of network flags (numbers or text) to use to filter data")
            self.parser.add_argument("--add_qa",
                                     dest="add_qa",
                                     help="list of qa flags (numbers or text) to add to default qa flags, used to filter data")
            self.parser.add_argument("--subtract_qa",
                                     dest="subtract_qa",
                                     help="list of qa flags (numbers or text) to subtract from default qa flags, used to filter data")
            self.parser.add_argument("--add_flags",
                                     dest="add_flags",
                                     help="list of network flags (numbers or text) to add to default network flags, used to filter data")
            self.parser.add_argument("--subtract_flags",
                                     dest="subtract_flags",
                                     help="list of network flags (numbers or text) to subtract from default network flags, used to filter data")
            self.parser.add_argument("--temporal_colocation",
                                     dest="temporal_colocation",
                                     help="activate temporal colocation betwen observations and experiments")
            self.parser.add_argument("--spatial_colocation",
                                     dest="spatial_colocation",
                                     help="activate spatial colocation between multiple read species")
            self.parser.add_argument("--spatial_colocation_tolerance",
                                     dest="spatial_colocation_tolerance",
                                     help="set spatial colocation tolerance to match stations by longitudes/latitudes and/or measurement_altitudes (in metres)")
            self.parser.add_argument("--spatial_colocation_station_reference",
                                     dest="spatial_colocation_station_reference",
                                     help="use 'station_reference' variable for spatial colocation?")
            self.parser.add_argument("--spatial_colocation_station_name",
                                     dest="spatial_colocation_station_name",
                                     help="use 'station_name' variable for spatial colocation?")
            self.parser.add_argument("--spatial_colocation_longitude_latitude",
                                     dest="spatial_colocation_longitude_latitude",
                                     help="use 'longitude' and 'latitude' variables for spatial colocation?")
            self.parser.add_argument("--spatial_colocation_measurement_altitude",
                                     dest="spatial_colocation_measurement_altitude",
                                     help="use 'measurement_altitude' variable for spatial colocation?")
            self.parser.add_argument("--spatial_colocation_validation",
                                     dest="spatial_colocation_validation",
                                     help="validate spatial colocation intersections via position using 'spatial_colocation_tolerance'?")
            self.parser.add_argument("--spatial_colocation_validation_tolerance",
                                     dest="spatial_colocation_validation_tolerance",
                                     help="set spatial colocation validation tolerance to validate station reference/station name match of stations by longitude/latitude position (in metres)")
            self.parser.add_argument("--map_extent",
                                     dest="map_extent",
                                     help="map extent for plots involving any maps")
            self.parser.add_argument("--filter_species",
                                     dest="filter_species",
                                     help="filter read species by other species within a given data range")
            self.parser.add_argument("--calibration_factor",
                                     dest="calibration_factor",
                                     help="factor to calibrate data")
            self.parser.add_argument("--lower_bound",
                                     dest="lower_bound",
                                     help="filter out data below this set lower bound")
            self.parser.add_argument("--upper_bound",
                                     dest="upper_bound",
                                     help="filter out data above this set upper bound")
            self.parser.add_argument("--remove_extreme_stations",
                                     dest="remove_extreme_stations",
                                     help="remove extreme stations using defined statistic limits")
            self.parser.add_argument("--report_type",
                                     dest="report_type",
                                     help="define plot options")
            self.parser.add_argument("--report_summary",
                                     dest="report_summary",
                                     help="activate summary plots in reports")
            self.parser.add_argument("--report_stations",
                                     dest="report_stations",
                                     help="activate station specific plots in reports")
            self.parser.add_argument("--report_title",
                                     dest="report_title",
                                     help="report title")
            self.parser.add_argument("--report_filename",
                                     dest="report_filename",
                                     help="report filename")
            self.parser.add_argument("--harmonise_summary",
                                     dest="harmonise_summary",
                                     help="harmonise axes limits across subsections for summary report")
            self.parser.add_argument("--harmonise_stations",
                                     dest="harmonise_stations",
                                     help="harmonise axes limits across stations for stations report")                     
            self.parser.add_argument("--active_dashboard_plots",
                                     dest="active_dashboard_plots",
                                     help="active plots on dashboard upon launch")    
            self.parser.add_argument("--resampling_resolution",
                                     dest="resampling_resolution",
                                     help="temporal resolution to resample data")
            self.parser.add_argument("--statistic_mode",
                                     dest="statistic_mode",
                                     help="name of statistical mode")
            self.parser.add_argument("--statistic_aggregation",
                                     dest="statistic_aggregation",
                                     help="type of statistic aggregation")
            self.parser.add_argument("--periodic_statistic_mode",
                                     dest="periodic_statistic_mode",
                                     help="name of periodic statistical mode")
            self.parser.add_argument("--periodic_statistic_aggregation",
                                     dest="periodic_statistic_aggregation",
                                     help="type of periodic statistic aggregation")
            self.parser.add_argument("--timeseries_statistic_aggregation",
                                     dest="timeseries_statistic_aggregation",
                                     help="type of timeseries statistic aggregation")
            self.parser.add_argument("--plot_characteristics_filename",
                                     dest="plot_characteristics_filename",
                                     help="set filename for plot characteristics")
            self.parser.add_argument("--interp_n_neighbours",
                                     dest="interp_n_neighbours",
                                     help="number of N nearest neighbours used for interpolation")
            self.parser.add_argument("--interp_reverse_vertical_orientation",
                                     dest="interp_reverse_vertical_orientation",
                                     help="reverse vertical orientation of model levels for interpolation")
            self.parser.add_argument("--interp_chunk_size",
                                     dest="interp_chunk_size",
                                     help="minimum number of jobs to pack in each chunk for interpolation")
            self.parser.add_argument("--interp_job_array_limit",
                                     dest="interp_job_array_limit",
                                     help="maximum number of chunks in job array for interpolation")
            self.parser.add_argument("--interp_multiprocessing",
                                     dest="interp_multiprocessing",
                                     help="use multiprocessing instead of greasy to interpolate in HPC machines")
            self.parser.add_argument("--logfile",
                                     dest="logfile",
                                     action='store_true',
                                     help="redirect the output to a log file")
 
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
