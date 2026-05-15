""" Class for Providentia report mode """

import copy
import os
import sys

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as fa
import numpy as np
from packaging.version import Version
import pandas as pd
from pypdf import PdfReader, PdfWriter
import yaml
import re

from providentia.auxiliar import CURRENT_PATH, join, expand_plot_characteristics
from .configuration import load_conf
from .configuration import ProvConfiguration
from .fields_menus import (init_representativity, init_period, init_metadata,
                           update_representativity_fields, update_period_fields, update_metadata_fields,
                           representativity_conf, period_conf, metadata_conf)
from .filter import DataFilter
from .plotting import Plotting
from .plot_aux import get_taylor_diagram_ghelper, set_map_extent
from .plot_formatting import format_plot_options, format_axis, harmonise_xy_lims_paradigm, set_axis_label, set_axis_title
from .read import DataReader
from .read_aux import (generate_file_trees, get_possible_resampling_resolutions, 
                       get_periodic_nonrelevant_temporal_resolutions, get_periodic_relevant_temporal_resolutions, 
                       get_valid_models, get_valid_obs_files_in_date_range)
from .statistics import (calculate_statistic, get_fairmode_data,
                         generate_colourbar, get_selected_station_data, get_z_statistic_info)
from .warnings_prv import show_message

PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
fairmode_settings = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/fairmode.yaml')))


class Report:
    """ Class to create Providentia reports. """

    def __init__(self, **kwargs):
        """
        Initialise the Report instance.

        Parameters
        ----------
        **kwargs : dict
            Arbitrary keyword arguments used to configure the report 
            settings and override default configuration values.
        """
        
        # update self with command line arguments
        self.commandline_arguments = copy.deepcopy(kwargs)

        # make sure that we are not using Qt5 backend with matplotlib
        matplotlib.use('Agg')

        # load statistical yamls
        self.basic_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/basic_stats.yaml')))
        self.modbias_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/model_bias_stats.yaml')))

        # load representativity information
        self.representativity_info = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/internal/representativity.yaml')))

        # initialise default configuration variables
        # modified by commandline arguments, if given
        ProvConfiguration(self, **self.commandline_arguments)

        self.logger.info("Creating a Providentia Report...")

        # update variables from config file
        if self.config != '':  
            read_conf = False
            if os.path.exists(self.config):
                read_conf = True
            else: 
                if os.path.exists(join(self.config_dir, self.config)):
                    self.config = join(self.config_dir, self.config)
                    read_conf = True
            if read_conf:
                load_conf(self, self.config)
                self.from_conf = True
            else:
                error = 'Error: The path to the configuration file specified in the command line does not exist.'
                self.logger.error(error)
                sys.exit(1)
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            self.logger.error(error)
            sys.exit(1)

        # load report plot presets
        try:
            self.report_plots = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/report_plots.yaml')))
        except:
            error = "Error: Report plots file could not be read, check for common typos (e.g. missing double quotes or commas)."
            self.logger.error(error)
            sys.exit(1)

        # select section if passed through command line arguments
        if "section" in self.commandline_arguments.keys():
            if self.commandline_arguments["section"] in self.all_sections:
                index = self.parent_section_names.index(self.commandline_arguments["section"].split("·")[0])
                self.parent_section_names = [self.parent_section_names[index]]
            else:
                error = "Error: Section {} does not exist in configuration file.".format(self.commandline_arguments["section"])
                self.logger.error(error)
                sys.exit(1)

        # if no parent section names are found throw an error
        if len(self.parent_section_names) == 0:
            error = "Error: No sections were found in the configuration file, make sure to name them using square brackets."
            self.logger.error(error)
            sys.exit(1)

    def run(self):
        """Execute the Providentia report workflow for all configured sections."""

        # iterate through configuration sections
        for section_ind, section in enumerate(self.parent_section_names):
            self.logger.info('\nStarting to create PDF for {} section'.format(section))

            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # initialize plot characteristics
            self.plot_characteristics = dict()

            # reinitialise default configuration variables
            # modified by commandline arguments, if given
            provconf = ProvConfiguration(self, **self.commandline_arguments)

            # get GHOST version before reading configuration file section parameters
            current_ghost_version = self.ghost_version 

            # update self with section variables (if not passed via command line)
            for k, val in self.section_opts.items():
                if k not in self.commandline_arguments:
                    setattr(self, k, provconf.parse_parameter(k, val))

            # if first section or GHOST version has changed
            if (section_ind == 0) or (current_ghost_version != self.ghost_version):
                generate_file_trees(self)

                # initialise DataReader class
                self.datareader = DataReader(self)

            # check for self defined plot characteristics file
            if self.plot_characteristics_filename == '':
                self.plot_characteristics_filename = join(PROVIDENTIA_ROOT, 'settings/plot_characteristics.yaml')
            plot_characteristics = yaml.safe_load(open(self.plot_characteristics_filename))
            self.plot_characteristics_templates = expand_plot_characteristics(plot_characteristics, 'report')
            self.plot_characteristics = {}

            # initialise Plotting class
            self.plotting = Plotting(read_instance=self, canvas_instance=self)

            # add general plot characteristics to self (if not passed via command line)
            for k, val in self.plot_characteristics_templates['general'].items():
                if k not in self.commandline_arguments:
                    setattr(self, k, val)

            # now all variables have been parsed, check validity of those, throwing errors where necessary
            provconf.check_validity()

            # if some filename has not been provided through the configuration file use default names
            if 'report_filename' in self.section_opts:
                filename = self.section_opts['report_filename']
            else:
                # add a number next to the filename to avoid overwriting
                filename = f'{self.report_filename}_{section_ind}'

                msg = 'Warning: Report filename (report_filename) has not been defined in '
                msg += f"configuration file section '{section}'. Setting filename to be '{filename}'."
                self.logger.info(msg)

            # set some key configuration variables
            self.periodic_relevant_temporal_resolutions = get_periodic_relevant_temporal_resolutions(self.resolution)
            self.periodic_nonrelevant_temporal_resolutions = get_periodic_nonrelevant_temporal_resolutions(self.resolution)
            self.data_labels = [self.observations_data_label] + list(self.experiments.values())
            self.data_labels_raw = [self.observations_data_label] + list(self.experiments.keys())
            self.networkspecies = ['{}|{}'.format(network,speci) for network, speci in zip(self.network, self.species)]

            # get valid observations in date range
            get_valid_obs_files_in_date_range(self, self.start_date, self.end_date)

            # update available models for selected fields
            get_valid_models(self, self.start_date, self.end_date, self.resolution,
                             self.network, self.species)

            # read data
            self.datareader.read_setup(['reset'])
            
            # initialise previous QA, flags, filter species and calibration factor as section values
            self.previous_qa = copy.deepcopy(self.qa)
            self.previous_flags = copy.deepcopy(self.flags)
            self.previous_filter_species = copy.deepcopy(self.filter_species)
            self.previous_calibration_factor = copy.deepcopy(self.calibration_factor)

            # if no valid data has been found be to be read, then skip to next section
            if self.invalid_read:
                self.logger.info('No valid data for {} section'.format(section))
                continue
            
            # check if report type is valid
            if self.report_type not in self.report_plots.keys():
                error = 'Error: The report type {0} cannot be found in settings/report_plots.yaml. '.format(self.report_type)
                error += 'The available report types are {0}. Select one or create your own.'.format(list(self.report_plots.keys()))
                self.logger.error(error)
                sys.exit(1)

            # set plots that need to be made (summary and station specific)
            self.summary_plots_to_make = []
            self.station_plots_to_make = []
            if isinstance(self.report_plots[self.report_type], list):
                for plot_type in self.report_plots[self.report_type]:
                    # there can be no station specific plots for map plot type
                    if plot_type[:4] != 'map-':
                        self.station_plots_to_make.append(plot_type)
                    # there can be no summary specific contingency tables (if gerrity, then allow)
                    if (plot_type[:16] != 'contingencytable'):
                        self.summary_plots_to_make.append(plot_type)
                    else:
                        if ('gerrity' in plot_type.split('_')[1:]):
                            self.summary_plots_to_make.append(plot_type)
            elif isinstance(self.report_plots[self.report_type], dict):
                # get summary plots
                if 'summary' in self.report_plots[self.report_type].keys():
                    if not self.report_summary:
                        msg = 'report_summary is False, summary plots will not be created.'
                        show_message(self, msg)
                    else:
                        # there can be no summary specific contingency tables (if gerrity, then allow)
                        for plot_type in self.report_plots[self.report_type]['summary']:
                            if (plot_type[:16] != 'contingencytable'):
                                self.summary_plots_to_make.append(plot_type)
                            else:
                                if ('gerrity' in plot_type.split('_')[1:]):
                                    self.summary_plots_to_make.append(plot_type)
                # get station plots
                if 'station' in self.report_plots[self.report_type].keys():
                    if not self.report_stations:
                        msg = 'report_stations is False, station plots will not be created.'
                        show_message(self, msg)
                    else:
                        for plot_type in self.report_plots[self.report_type]['station']:
                            # there can be no station specific plots for map plot type
                            if plot_type[:4] != 'map-':
                                self.station_plots_to_make.append(plot_type)

            # TODO: For Taylor diagrams, remove this piece of code when we stop using Matplotlib 3.3
            if Version(matplotlib.__version__) < Version("3.8"):
                if plot_type[:6] == 'taylor':
                    if (plot_type in self.station_plots_to_make) or (plot_type in self.summary_plots_to_make):
                        error = 'It is not possible to create Taylor diagrams yet, please remove from settings/report_plots.yaml.'
                        sys.exit(error)

            # check if there are multispecies plots
            multispecies = False
            for plot_type in self.summary_plots_to_make:
                if 'multispecies' in plot_type:
                    multispecies = True
                    break
            for plot_type in self.station_plots_to_make:
                if 'multispecies' in plot_type:
                    multispecies = True
                    break
                
            if (multispecies) and (len(np.unique(list(self.measurement_units.values()))) > 1):
                msg = f"Units in the multispecies plots will be converted to 'multispecies_units' ({self.multispecies_units}) for consistency. "
                show_message(self, msg)
                if self.multispecies_units in [None, ""]:
                    error = f"Please specify the units in your configuration file by adding 'multispecies_units'. "
                    error += f"Units for each species are: {self.measurement_units}."
                    sys.exit(error)

            # set plot characteristics for all plot types (summary, station)
            self.plots_to_make = list(self.summary_plots_to_make)
            self.plots_to_make.extend(x for x in self.station_plots_to_make
                                      if x not in self.summary_plots_to_make)
            self.plotting.set_plot_characteristics(self.plots_to_make)
            
            # define dictionary to store plot figures per page
            self.plot_dictionary = {}

            # start making PDF
            self.start_pdf(filename)

            # remove section variables from memory 
            for k in self.section_opts:
                try:
                    vars(self).pop(k)
                except:
                    pass

            if section_ind != len(self.parent_section_names) - 1:
                self.logger.info('\n'+'='*70)

    def start_pdf(self, filename):
        """
        Initialises and builds the PDF document structure, managing plotting paradigms and page ordering.

        Parameters
        ----------
        filename : str
            The name or path for the generated PDF report.
        """
        
        # get path where reports will be saved
        if '/' in filename:
            if os.path.isdir(os.path.dirname(filename)):
                reports_path = filename
        else:
            reports_path = (join(PROVIDENTIA_ROOT, 'reports/')) + filename

        # create reports folder
        if not os.path.exists(os.path.dirname(reports_path)):
            if '/' in reports_path:
                self.logger.info('Path {0} does not exist and it will be created.'.format(os.path.dirname(reports_path)))
            os.makedirs(os.path.dirname(reports_path))

        # add termination .pdf to filenames
        if '.pdf' not in reports_path:
            reports_path += '.pdf'

        # create 'temp' reports path for first writing to before finalising (for compression) 
        reports_path_temp = '{}_temp.pdf'.format(reports_path.split('.pdf')[0])
        reports_doi_path_temp = '{}_doi_temp.pdf'.format(reports_path.split('.pdf')[0])

        # open new PDF file
        with PdfPages(reports_path_temp) as pdf:
            
            self.pdf = pdf

            # initialise dictionaries to store relevant page numebrs
            if self.report_summary:
                self.summary_pages = {}
            if self.report_stations:
                self.station_pages = {}

            # get subsection names
            self.child_subsection_names = [subsection_name for subsection_name in self.subsection_names 
                                           if self.section == subsection_name.split('·')[0]]
            
            # if subsection has been passed in command line arguments get subsection
            if ("section" in self.commandline_arguments.keys()) and ("·" in self.commandline_arguments["section"]):
                self.subsections = [self.commandline_arguments["section"]]
            else:
                if len(self.child_subsection_names) > 0:
                    self.subsections = self.child_subsection_names
                else:
                    self.subsections = [self.section]

            # make header
            self.plotting.set_plot_characteristics(['header'])
            self.plotting.make_header(self.pdf, self.plot_characteristics['header'])

            # create variables to keep track of minimum and maximum data ranges across subsections
            self.data_range_min_summary = {networkspeci:np.inf for networkspeci in self.networkspecies}
            self.data_range_min_station = {networkspeci:np.inf for networkspeci in self.networkspecies}
            self.data_range_max_summary = {networkspeci:0 for networkspeci in self.networkspecies}
            self.data_range_max_station = {networkspeci:0 for networkspeci in self.networkspecies}
            self.stddev_max_summary = {networkspeci:0 for networkspeci in self.networkspecies}
            self.stddev_max_station = {networkspeci:0 for networkspeci in self.networkspecies}

            # make all plots per subsection (for all plot types except distribution/taylor plots)
            summary_plots_to_make = [plot_type for plot_type in self.summary_plots_to_make 
                                     if ('distribution' not in plot_type) and ('taylor' not in plot_type)]
            station_plots_to_make = [plot_type for plot_type in self.station_plots_to_make 
                                     if ('distribution' not in plot_type) and ('taylor' not in plot_type)]
            self.make_plots_per_subsection(summary_plots_to_make, station_plots_to_make, 
                                           do_plot_geometry_setup=True)

            # make all plots per subsection
            # for distribution/taylor plot types --> done so to calculate data ranges 
            # across subsections first
            summary_plots_to_make = [plot_type for plot_type in self.summary_plots_to_make 
                                     if ('distribution' in plot_type) or ('taylor' in plot_type)]
            station_plots_to_make = [plot_type for plot_type in self.station_plots_to_make 
                                     if ('distribution' in plot_type) or ('taylor' in plot_type)]
            if (len(summary_plots_to_make) > 0) or (len(station_plots_to_make) > 0):
                self.make_plots_per_subsection(summary_plots_to_make, station_plots_to_make)
            
            # finalise formatting for plots
            # create colourbars
            # harmonise xy limits(not for map, heatmap or table, or when xlim and ylim defined)

            # remove header from plot characteristics dictionary
            if 'header' in list(self.plot_characteristics.keys()):
                del self.plot_characteristics['header']

            # set variables to inform when have formatted 1 set of networkspecies plots for stations
            formatted_networkspeci_plots = False        
            did_formatting = False

            # iterate through networks and species
            for networkspeci in self.networkspecies: 

                # iterate through plot types
                for plot_type in list(self.plot_characteristics.keys()):

                    # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                    plot_options = plot_type.split('_')[1:]

                    # if a multispecies plot is active then only format on first pass
                    if ('multispecies' in plot_options) & (formatted_networkspeci_plots):
                        continue

                    # get zstat information from plot_type
                    zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type=plot_type)

                    # get base plot type (without stat and options)
                    if zstat:
                        base_plot_type = plot_type.split('-')[0] 
                    else:
                        base_plot_type = plot_type.split('_')[0] 

                    # get relevant paradigm pages to harmonise axes limits for
                    paradigm_pages = {'summary':[], 'station':[]}
                    if (self.report_summary) & (self.report_stations):
                        if plot_type in self.summary_pages:
                            if networkspeci in self.summary_pages[plot_type]:
                                paradigm_pages['summary'] = self.summary_pages[plot_type][networkspeci]
                        # for multispecies plots of specific stations, spatial colocation needs to be on.
                        # if not, the plots will not exist so handle this
                        if ('multispecies' in plot_options) & (not self.spatial_colocation):
                            paradigm_pages['station'] = []
                        else:
                            if plot_type in self.station_pages:
                                if networkspeci in self.station_pages[plot_type]:
                                    for subsection in self.subsections:
                                        if subsection in self.station_pages[plot_type][networkspeci]:
                                            paradigm_pages['station'].extend(self.station_pages[plot_type][networkspeci][subsection])

                    elif self.report_summary:
                        if plot_type in self.summary_pages:
                            if networkspeci in self.summary_pages[plot_type]:
                                paradigm_pages['summary'] = self.summary_pages[plot_type][networkspeci]
                    
                    elif self.report_stations:
                        # for multispecies plots of specific stations, spatial colocation needs to be on.
                        # if not, the plots will not exist so handle this
                        if ('multispecies' in plot_options) & (not self.spatial_colocation):
                            paradigm_pages['station'] = []
                        else:
                            if plot_type in self.station_pages:
                                if networkspeci in self.station_pages[plot_type]:
                                    for subsection in self.subsections:
                                        if subsection in self.station_pages[plot_type][networkspeci]:
                                            paradigm_pages['station'].extend(self.station_pages[plot_type][networkspeci][subsection])

                    # iterate through paradigm pages
                    for plotting_paradigm, relevant_pages in paradigm_pages.items():                        
                        if len(relevant_pages) == 0:
                            continue

                        # get relevant axs and plot types per networkspeci / plot type
                        relevant_axs, relevant_data_labels = self.get_relevant_axs_per_networkspeci_plot_type(base_plot_type, 
                                                                                                              relevant_pages)

                        # if ax is not visible, remove ax from relevant_axs
                        # this means plot was not created
                        axs_to_remove = []
                        for ax in relevant_axs:
                            if not ax.get_visible():
                                axs_to_remove.append(ax)
                        relevant_axs = [ax for ax in relevant_axs if ax not in axs_to_remove]

                        # if have no relevant axs, continue to next paradigm
                        if len(relevant_axs) == 0:
                            continue

                        # get data ranges for plotting paradigm
                        if plotting_paradigm == 'summary':
                            data_range_min = self.data_range_min_summary[networkspeci]
                            data_range_max = self.data_range_max_summary[networkspeci]
                            stddev_max = self.stddev_max_summary[networkspeci]
                        elif plotting_paradigm == 'station':
                            data_range_min = self.data_range_min_station[networkspeci]
                            data_range_max = self.data_range_max_station[networkspeci]
                            stddev_max = self.stddev_max_station[networkspeci]

                        # generate colourbars for required plots in paradigm on each relevant page
                        if 'cb' in list(self.plot_characteristics[plot_type].keys()):
                            # get all cb_axs for plot_type across relevant pages
                            cb_axs = [self.plot_dictionary[relevant_page]['cb_ax'] for relevant_page in relevant_pages]
                            generate_colourbar(self, relevant_axs, cb_axs, zstat, self.plot_characteristics[plot_type], 
                                               networkspeci.split('|')[-1])

                        # harmonise xy limits for plot paradigm
                        harmonise = False
                        if (((plotting_paradigm == 'summary') and (self.harmonise_summary))
                            or ((plotting_paradigm == 'station') and (self.harmonise_stations))):
                            harmonise = True 
                        if (base_plot_type not in ['map', 'heatmap', 'table', 'taylor', 'fairmode-statsummary']): 
                            if base_plot_type == 'scatter':
                                harmonise_xy_lims_paradigm(self, self, relevant_axs, base_plot_type, 
                                                           self.plot_characteristics[plot_type], plot_options, 
                                                           relim=True, harmonise=harmonise)
                            else:
                                harmonise_xy_lims_paradigm(self, self, relevant_axs, base_plot_type,
                                                           self.plot_characteristics[plot_type], plot_options, 
                                                           relim=True, autoscale=True, harmonise=harmonise)
                        # update variable to reflect some formatting was performed
                        did_formatting = True

                # update variables to show if a networkspeci has been formatted
                if did_formatting:
                    formatted_networkspeci_plots = True

            # initialise arrays with pages that include multispecies plots
            self.summary_multispecies_pages = []
            self.station_multispecies_pages = []

            # save page figures
            valid_page = False
            real_page = 1
            for i, page in enumerate(self.plot_dictionary):
                # if page has no active data plotted, do not plot it
                n_page_plotted_labels = 0
                for ax_dict in self.plot_dictionary[page]['axs']:
                    n_page_plotted_labels += len(ax_dict['data_labels'])
                if n_page_plotted_labels > 0:
                    # the following variables will be used to set the order of multispecies plot
                    # in the report if there are
                    # get pages in PDF that contain multispecies plots
                    if 'multispecies' in self.plot_dictionary[page]['plot_type']:
                        if self.plot_dictionary[page]['paradigm'] == 'summary':
                            self.summary_multispecies_pages.append(real_page)
                        elif self.plot_dictionary[page]['paradigm'] == 'station':
                            self.station_multispecies_pages.append(real_page)

                    # save page where station plots start to be created
                    if ((self.plot_dictionary[page]['paradigm'] == 'station') and 
                        (not hasattr(self, 'paradigm_break_page'))):
                        self.paradigm_break_page = real_page

                    # save page figure                        
                    fig = self.plot_dictionary[page]['fig']
                    self.pdf.savefig(fig, dpi=self.dpi)
                    plt.close(fig)
                    real_page += 1
                    valid_page = True

            # if last page was reached that means there were plots to write (valid_page = True)
            if valid_page:
                self.logger.info(f'\nWriting PDF in {reports_path}')
            else:    
                self.logger.info('\n0 plots remain to write to PDF')

        # compress PDF using ghostscript if desired (at 300 DPI)
        if self.compression:
            self.logger.info('\nCompressing PDF\n')
            os.system("gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH -sOutputFile={} {}".format(reports_path,reports_path_temp))
            os.system("rm {}".format(reports_path_temp))
        else:
            os.system("mv {} {}".format(reports_path_temp, reports_path))

        # make page with DOIs
        doi_pdf = None
        if 'actris/actris' in self.network:
            self.plotting.set_plot_characteristics(['doi'])
            doi_pdf = self.plotting.make_doi_pdf(self.plot_characteristics['doi'], reports_doi_path_temp)

        # reorder pages
        # if only summary plots have been made, set paradigm break page to last page
        if not hasattr(self, 'paradigm_break_page'):
            pdf_file = PdfReader(open(reports_path, "rb"))
            self.paradigm_break_page = len(pdf_file.pages)
        
        self.reorder_pdf_pages(reports_path, reports_path, self.summary_multispecies_pages, 
                               self.station_multispecies_pages, self.paradigm_break_page, doi_pdf, 
                               reports_doi_path_temp)

    def setup_plot_geometry(self, plotting_paradigm, networkspeci, have_setup_multispecies):
        """
        Sets up the Matplotlib figure and axes geometry for either summary or station-specific plots.

        Parameters
        ----------
        plotting_paradigm : str
            The context of the plots, either 'summary' or 'station'.
        networkspeci : str
            The combined network and species identifier.
        have_setup_multispecies : bool
            Flag indicating if multi-species geometry has already been initialised.
        """

        # depending on plot type set plots to make
        if plotting_paradigm == 'summary':
            plots_to_make = copy.deepcopy(self.summary_plots_to_make)
        elif plotting_paradigm == 'station':
            plots_to_make = copy.deepcopy(self.station_plots_to_make)

        # iterate through plot types to make
        for plot_type in plots_to_make:

            # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
            plot_options = plot_type.split('_')[1:]
            
            # if making a multispecies plot per specific station, spatial colocation must be also active
            if 'multispecies' in plot_options:
                if plotting_paradigm == 'summary':
                    if have_setup_multispecies:
                        continue
                elif plotting_paradigm == 'station':
                    if not self.spatial_colocation:
                        msg = f'{plot_type} cannot be created per station '
                        msg += 'without activating the spatial colocation.'
                        show_message(self, msg)
                        continue
                    elif (have_setup_multispecies):
                        continue
            
            # create variables to store list of page numbers per plot type / networkspeci / subsection (if do not exist)
            if plotting_paradigm == 'summary':
                if plot_type not in self.summary_pages:
                    self.summary_pages[plot_type] = {}
                if networkspeci not in self.summary_pages[plot_type]:
                    self.summary_pages[plot_type][networkspeci] = []
            if plotting_paradigm == 'station':
                if plot_type not in self.station_pages:
                    self.station_pages[plot_type] = {}
                if networkspeci not in self.station_pages[plot_type]:
                    self.station_pages[plot_type][networkspeci] = {}
                if self.subsection not in self.station_pages[plot_type][networkspeci]:
                    self.station_pages[plot_type][networkspeci][self.subsection] = []

            # get zstat information from plot_type
            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type=plot_type)

            # get base plot type (without stat and options)
            if zstat:
                base_plot_type = plot_type.split('-')[0] 
            else:
                base_plot_type = plot_type.split('_')[0] 

            # get copy of plot characteristics for plot type
            plot_characteristics = copy.deepcopy(self.plot_characteristics[plot_type])
            plot_characteristics_vars = list(plot_characteristics.keys())

            # update page title depending on plot paradigm
            if plotting_paradigm == 'summary':
                if 'multispecies' in plot_options:
                    plot_characteristics['page_title']['t'] = '{} (Summary)\nmultispecies'.format(plot_characteristics['page_title']['t']) 
                else:
                    plot_characteristics['page_title']['t'] = '{} (Summary)\n{}'.format(plot_characteristics['page_title']['t'], 
                                                                                        networkspeci) 
            elif plotting_paradigm == 'station':
                if 'multispecies' in plot_options:
                    plot_characteristics['page_title']['t'] = '{} (Per Station)\n{}\nmultispecies'.format(plot_characteristics['page_title']['t'], 
                                                                                                          self.subsection) 
                else:
                    plot_characteristics['page_title']['t'] = '{} (Per Station)\n{}\n{}'.format(plot_characteristics['page_title']['t'], 
                                                                                                self.subsection, networkspeci) 

            # define number of plots per type
            n_plots_per_plot_type = False
            if base_plot_type == 'map':
                if 'obs' in plot_options:
                    n_plots_per_plot_type = len(self.subsections)
                elif z_statistic_sign == 'bias':
                    n_plots_per_plot_type = len(self.subsections) * \
                                            (len(self.data_labels) - 1)
                else:
                    n_plots_per_plot_type = len(self.subsections) * \
                                            len(self.data_labels)
            elif base_plot_type in ['heatmap', 'table']:
                if plotting_paradigm == 'summary':
                    n_plots_per_plot_type = 1
                elif plotting_paradigm == 'station':
                    n_plots_per_plot_type = self.n_stations
            else:
                if plotting_paradigm == 'summary':
                    if 'individual' in plot_options:
                        if (base_plot_type in ['scatter', 'fairmode-target', 'fairmode-statsummary']) or ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                            n_plots_per_plot_type = len(self.subsections) * \
                                                    (len(self.data_labels) - 1)

                        else:
                            n_plots_per_plot_type = len(self.subsections) * \
                                                    len(self.data_labels)
                    else:
                        n_plots_per_plot_type = len(self.subsections) 
                elif plotting_paradigm == 'station':
                    if 'individual' in plot_options:
                        if (base_plot_type in ['scatter', 'fairmode-target', 'fairmode-statsummary']) or ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                            n_plots_per_plot_type = self.n_stations * \
                                                    (len(self.data_labels) - 1) 
                        else:
                            n_plots_per_plot_type = self.n_stations * \
                                                    len(self.data_labels) 
                    else:
                        n_plots_per_plot_type = self.n_stations

            # get n pages per plot type
            n_pages_per_plot_type = int(np.ceil(n_plots_per_plot_type / (
                                    plot_characteristics['figure']['ncols'] * plot_characteristics['figure']['nrows'])))
            plot_ii_per_type = 0

            # iterate through n pages per plot type
            for page_n in range(self.n_total_pages, self.n_total_pages + n_pages_per_plot_type):

                if base_plot_type == 'map':
                    plot_characteristics['figure']['subplot_kw'] = {'projection': self.plotcrs}
                elif base_plot_type == 'taylor':
                    reference_stddev = 7.5
                    ghelper = get_taylor_diagram_ghelper(reference_stddev, plot_characteristics)
                    plot_characteristics['figure']['subplot_kw'] = {'axes_class': fa.FloatingAxes,
                                                                    'grid_helper': ghelper}
                fig, axs = plt.subplots(**plot_characteristics['figure'])

                # each page is handled as 1 figure
                # intialise page if not yet done
                if page_n not in self.plot_dictionary:
                    self.plot_dictionary[page_n] = {'fig': fig, 'plot_type': plot_type, 'axs': [], 
                                                    'paradigm': plotting_paradigm}

                # make page title?
                if 'page_title' in plot_characteristics_vars:
                    st = fig.suptitle(**plot_characteristics['page_title'])

                # iterate through axes (by row, then column)
                row_ii = -1
                col_ii = copy.deepcopy(plot_characteristics['figure']['ncols'])
                
                # flatten axis for iteration (when we have more than one axis per page)
                if not isinstance(axs, np.ndarray):
                    axs = np.array(axs)

                for ax_ii, ax in enumerate(axs.flatten()):

                    if col_ii == plot_characteristics['figure']['ncols']:
                        row_ii += 1
                        col_ii = 0
                    if row_ii == (plot_characteristics['figure']['nrows'] - 1):
                        last_row_on_page = True
                    else:
                        last_row_on_page = False

                    # force rasterized (bitmap) drawing in vector backend output.
                    ax.set_rasterized(True)

                    # keep iteratively plotting until have satisfied needed plots per type
                    if plot_ii_per_type < n_plots_per_plot_type:
                        
                        # determine if are on last valid row to plot
                        if (n_plots_per_plot_type - plot_ii_per_type) <= plot_characteristics['figure']['ncols']:
                            last_valid_row = True
                        else:
                            last_valid_row = False

                        # setup periodic plot type gridspec
                        if base_plot_type in ['periodic', 'periodic-violin']:
                            gs = gridspec.GridSpecFromSubplotSpec(100, 100, subplot_spec=ax.get_subplotspec())
                            grid_dict = dict()
                            grid_dict['hour'] = fig.add_subplot(gs[:46, :])
                            grid_dict['dayofweek'] = fig.add_subplot(gs[54:, 64:])
                            grid_dict['month'] = fig.add_subplot(gs[54:, :61])
                            self.plot_dictionary[page_n]['axs'].append({'handle':grid_dict, 'data_labels':[]})
                            ax.spines['top'].set_color('none')
                            ax.spines['bottom'].set_color('none')
                            ax.spines['left'].set_color('none')
                            ax.spines['right'].set_color('none')
                            ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
                            grid_dict['hour'].axis('off')
                            grid_dict['dayofweek'].axis('off')
                            grid_dict['month'].axis('off')
                            
                            # format axis
                            format_axis(self, self, grid_dict, base_plot_type, plot_characteristics,
                                        relevant_temporal_resolutions=self.periodic_relevant_temporal_resolutions)

                            # get references to periodic label annotations made, and then hide them
                            for relevant_temporal_resolution in self.periodic_relevant_temporal_resolutions:
                                annotations = [child for child in grid_dict[relevant_temporal_resolution].get_children() 
                                               if isinstance(child, matplotlib.text.Annotation)]
                                # hide annotations
                                for annotation in annotations:
                                    annotation.set_visible(False)
                        
                        # setup fairmode summary stats plot type gridspec
                        elif base_plot_type == 'fairmode-statsummary':
                            # get current species
                            speci = networkspeci.split('|')[1]

                            # get number of rows and columns
                            ncols = 4
                            nrows = 8 if speci in ["sconco3", "sconcno2", "pm10"] else 7
                            
                            # create gridspec and add it to a list
                            gs = gridspec.GridSpec(nrows, ncols, **plot_characteristics["gridspec_kw"])
                            grid_list = [fig.add_subplot(gs[i, j]) for i in range(nrows) for j in range(ncols)]
                            
                            # setup dictionary
                            self.plot_dictionary[page_n]['axs'].append({'handle':grid_list, 'data_labels':[]})

                            # format axis 
                            format_axis(self, self, grid_list, base_plot_type, plot_characteristics)
                        
                        # rest of plot types
                        else:
                            
                            # collect axes in dict
                            self.plot_dictionary[page_n]['axs'].append({'handle':ax, 'data_labels':[]})
                            
                            # format axis 
                            format_axis(self, self, ax, base_plot_type, plot_characteristics)

                    # turn off axes until some data is plottted
                    ax.axis('off')
                    ax.set_visible(False)
                    
                    # iterate plot and column counts
                    plot_ii_per_type += 1
                    col_ii += 1

                # set figure attributes
                # adjust subplots?
                if 'subplots_adjust' in plot_characteristics_vars:
                    if base_plot_type == 'contingencytable':
                        if 'gerrity' in plot_options:
                            fig.subplots_adjust(**plot_characteristics['subplots_adjust']['gerrity'])
                        else:
                            fig.subplots_adjust(**plot_characteristics['subplots_adjust']['contingency'])
                    else:
                        fig.subplots_adjust(**plot_characteristics['subplots_adjust'])

                # make legend?
                if 'legend' in plot_characteristics_vars:
                    if (base_plot_type in ['scatter', 'fairmode-target']) or ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                        set_obs = False
                    else:
                        set_obs = True
                    plot_characteristics['legend'] = self.plotting.make_legend_handles(plot_characteristics['legend'], 
                                                                                       set_obs=set_obs)
                    fig.legend(**plot_characteristics['legend']['plot'])

                # add colourbar axis to plot dictionary (if not already there)?
                if 'cb' in plot_characteristics_vars:
                    if 'cb_ax' not in self.plot_dictionary[page_n]:
                        self.plot_dictionary[page_n]['cb_ax'] = fig.add_axes(plot_characteristics['cb']['position'])
                        self.plot_dictionary[page_n]['cb_ax'].set_rasterized(True)
                            
                # add current page number
                if plotting_paradigm == 'summary':
                    self.summary_pages[plot_type][networkspeci].append(page_n)
                elif plotting_paradigm == 'station':
                    self.station_pages[plot_type][networkspeci][self.subsection].append(page_n)
                
                # add to total number of pages
                self.n_total_pages += 1

    def make_plots_per_subsection(self, summary_plots_to_make, station_plots_to_make, do_plot_geometry_setup=False):
        """
        Coordinates the iteration through report subsections, handling data filtering and triggering plot generation.

        Parameters
        ----------
        summary_plots_to_make : list
            List of plot types to be generated at the network-summary level.
        station_plots_to_make : list
            List of plot types to be generated for individual stations.
        do_plot_geometry_setup : bool, optional
            Flag to determine if Matplotlib figure geometry should be initialised (default is False).
        """

        # create variable to keep track if have setup summary plot geometry yet (done for all subsections at once)
        self.summary_plot_geometry_setup = False
        self.do_plot_geometry_setup = do_plot_geometry_setup

        # define dictionary to store stats from all subsections for heatmap and table plots
        self.stats_summary = {}
        self.stats_station = {}

        # define dictionary to store station references and dois per subsection
        self.station_reference_data = {}
        self.station_doi_data = {}

        # iterate through subsections
        for subsection_ind, subsection in enumerate(self.subsections):
            
            self.subsection_ind = subsection_ind
            self.subsection = subsection

            # create nested dictionary to store statistical information across all subsections
            if self.subsection not in self.stats_summary:
                self.stats_summary[self.subsection] = {}
            if self.subsection not in self.stats_station:
                self.stats_station[self.subsection] = {}

            # create nested dictionary to store DOI information across all subsections
            if self.subsection not in self.station_reference_data:
                self.station_reference_data[self.subsection] = {}
            if self.subsection not in self.station_doi_data:
                self.station_doi_data[self.subsection] = {}

            # initialise forecast models as None
            forecast_models = None

            # update the conf options for defined subsection
            if len(self.child_subsection_names) > 0:

                # we save the species to make sure that the species that do not have data after read
                # are not being set later in parse parameter
                read_species = copy.deepcopy(self.species)
                read_networkspecies = copy.deepcopy(self.networkspecies)

                # if have forecast active then save current model variable in memory as may need to set it if not re-reading data 
                if len(self.forecast) != 0:
                    forecast_models = copy.deepcopy(self.experiments)  

                # get subsection variables
                self.subsection_opts = self.sub_opts[self.subsection]

                # ensure all fixed section variables defined in subsection have same value as current section variables
                self.subsection_opts = {k: (self.section_opts[k] if k in self.fixed_section_vars else val) 
                                        for (k, val) in self.subsection_opts.items()}

                # reinitialise default configuration variables
                # modified by commandline arguments, if given
                provconf = ProvConfiguration(self, **self.commandline_arguments)

                # update subsection variables (if not passed via command line)
                for k, val in self.subsection_opts.items():
                    if k not in self.commandline_arguments:
                        value = provconf.parse_parameter(k, val, deactivate_warning=True)
                        setattr(self, k, value)

                # now all variables have been parsed, check validity of those, throwing errors where necessary
                provconf.check_validity(deactivate_warning=True)

                # keep only the species that were available when we read the parent section data
                self.species = read_species
                self.networkspecies = read_networkspecies

            # determine if need to re-read data (qa, flags, filter_species or calibration factor have changed)
            if (np.array_equal(self.qa, self.previous_qa) == False) or (
                np.array_equal(self.flags, self.previous_flags) == False) or (
                str(dict(sorted(self.filter_species.items()))) != str(dict(sorted(self.previous_filter_species.items())))) or (
                str(dict(sorted(self.calibration_factor.items()))) != str(dict(sorted(self.previous_calibration_factor.items())))):
                # reset data labels and models for forecast cases
                self.data_labels = self.original_data_labels
                self.data_labels_raw = self.original_data_labels_raw
                self.experiments = self.original_models
                # re-read data
                self.datareader.read_setup(['reset'])
            else:
                # if not re-reading data and forecast active, then overwrite models variable
                if forecast_models is not None:
                    self.experiments = forecast_models

                # update fields available for filtering
                init_representativity(self)
                update_representativity_fields(self)
                representativity_conf(self)
                init_period(self)
                update_period_fields(self)
                period_conf(self)
                init_metadata(self)

                # for non-GHOST delete valid station indices variables because we do not want to 
                # remove the stations with 0 valid measurements before the filter has been updated, 
                # this will happen later
                if hasattr(self, 'valid_station_inds') and (not self.reading_ghost):
                    delattr(self, 'valid_station_inds')
                    delattr(self, 'valid_station_inds_temporal_colocation')

                update_metadata_fields(self)
                metadata_conf(self)

            # set previous QA, flags, filter species and calibration factor as subsection
            self.previous_qa = copy.deepcopy(self.qa)
            self.previous_flags = copy.deepcopy(self.flags)
            self.previous_filter_species = copy.deepcopy(self.filter_species)
            self.previous_calibration_factor = copy.deepcopy(self.calibration_factor)

            # filter dataset for current subsection
            self.logger.info('\nFiltering data for {} subsection'.format(self.subsection))
            DataFilter(self)

            # for non-GHOST, we call update_metadata_fields after filtering to remove the stations that have
            # 0 valid measurements, to do this we need to have valid_station_inds, which is obtained 
            # after filtering
            if not self.reading_ghost:
                update_metadata_fields(self)

            # iterate through networks and species, creating plots
            self.n_total_pages = len(self.plot_dictionary)

            # make summary plots?
            if self.report_summary:

                # set variable to inform when have made 1 set of networkspecies plots for summary
                self.made_networkspeci_summary_plots = False

                # iterate through networkspecies
                for networkspeci in self.networkspecies:

                    # make summary plots
                    self.make_summary_plots(networkspeci, summary_plots_to_make)

                # update variable to keep track if have setup summary plot geometry yet for a subsection
                if self.made_networkspeci_summary_plots:
                    self.summary_plot_geometry_setup = True

            # make station specific plots?
            if self.report_stations and self.station_plots_to_make:      

                # set variable to inform when have made 1 set of networkspecies plots for stations
                self.made_networkspeci_station_plots = False        

                # iterate through networkspecies
                for networkspeci in self.networkspecies:

                    # make plots per station
                    self.make_station_plots(networkspeci, station_plots_to_make)

            # remove subsection variables from memory (if have subsections)
            # do not remove fixed section variables
            if (len(self.child_subsection_names) > 0):
                for k in self.subsection_opts:
                    if k not in self.fixed_section_vars:
                        try:
                            vars(self).pop(k)
                        except:
                            pass

    def reorder_pdf_pages(self, input_pdf, output_pdf, summary_multispecies_pages, 
                          station_multispecies_pages, paradigm_break_page, doi_pdf, reports_doi_path_temp):
        """
        Reorder PDF pages so that multispecies plots appear before other plots and DOI pages appear last.

        Parameters
        ----------
        read_instance : object
            Instance of class containing the PDF report.
        input_pdf : str
            Path to the original PDF file.
        output_pdf : str
            Path where the reordered PDF will be saved.
        summary_multispecies_pages : list
            Pages that contain summary multispecies plots.
        station_multispecies_pages : list
            Pages that contain station multispecies plots.
        paradigm_break_page : int
            Page index where station plots start.
        doi_pdf : str or None
            Path to DOI PDF to append at the end of the reordered PDF or None.
        reports_doi_path_temp : str
            Temporary path for the DOI PDF used during processing.

        Returns
        -------
        None
            Writes a reordered PDF file to `output_pdf`.
        """

        if (len(self.summary_multispecies_pages) > 0) or (len(self.station_multispecies_pages) > 0):
            self.logger.info('\nReordering pages')

        # Get pages
        summary_multispecies_pages = np.array(summary_multispecies_pages)
        station_multispecies_pages = np.array(station_multispecies_pages)

        # Get original order
        input_pdf_file = PdfReader(open(input_pdf, "rb"))
        all_pages = np.arange(len(input_pdf_file.pages))

        # Initialise page order
        page_order = copy.deepcopy(all_pages)

        # Move summary pages after page 0 (cover)
        if len(summary_multispecies_pages) > 0:
            summary_multispecies_pages = np.concatenate((np.array([0]), summary_multispecies_pages))
            summary_other_plots_pages = all_pages[~np.isin(all_pages, summary_multispecies_pages)]
            page_order = np.concatenate((
                summary_multispecies_pages, 
                summary_other_plots_pages)).tolist()
        
        # Move station pages after paradigm break page (when we start to see station plots)
        if len(station_multispecies_pages) > 0:
            station_pages = all_pages[paradigm_break_page:]
            station_other_plots_pages = station_pages[~np.isin(station_pages, station_multispecies_pages)]
            page_order = np.concatenate((
                page_order[:paradigm_break_page], 
                station_multispecies_pages, 
                station_other_plots_pages)).tolist()

        # Reorder pages
        output_pdf_file = PdfWriter()
        for page_number in page_order:
            output_pdf_file.add_page(input_pdf_file.pages[int(page_number)])

        # Add DOI pages at the end
        if doi_pdf is not None and os.path.exists(reports_doi_path_temp):
            input_doi_pdf = PdfReader(open(reports_doi_path_temp, "rb"))
            for page_number in range(len(input_doi_pdf.pages)):
                output_pdf_file.add_page(input_doi_pdf.pages[page_number])
            os.system("rm {}".format(reports_doi_path_temp))
        
        # Write the rearranged pages to a new PDF file
        self.logger.info(f'Writing {output_pdf}')
        with open(output_pdf, "wb") as outputStream:
            output_pdf_file.write(outputStream)

    def get_dois_from_station_reference(self, station_reference, networkspeci):
        """
        Get DOIs from station references. Only used for ACTRIS network downloads.

        Parameters
        ----------
        station_reference : str
            Current station reference
        networkspeci : str
            The combined network and species identifier.

        Returns
        -------
        list
            List of DOIs per month
        """

        current_station_dois = []
        for month_metadata in self.selected_station_metadata[networkspeci]:
            for metadata_vars in month_metadata:
                if metadata_vars[0] != self.current_station_reference:
                    continue
                for var in metadata_vars:
                    if isinstance(var, str) and re.compile(r'^https?://doi\.org/').match(var):
                        current_station_dois.append(var)

        return current_station_dois
                
    def make_summary_plots(self, networkspeci, summary_plots_to_make):
        """
        Generates all summary-level visualisations for a specific network-species combination within a subsection.

        Parameters
        ----------
        networkspeci : str
            The combined network and species identifier.
        summary_plots_to_make : list
            List of plot types to be generated at the network-summary level.
        """

        # setup plotting geometry for summary plots per networkspeci (for all subsections)
        if (not self.summary_plot_geometry_setup) & (self.do_plot_geometry_setup):
            self.setup_plot_geometry('summary', networkspeci, self.made_networkspeci_summary_plots)

        # get valid station inds for networkspeci 
        if self.temporal_colocation:
            self.relevant_station_inds = self.valid_station_inds_temporal_colocation[networkspeci][self.observations_data_label]
        else:
            self.relevant_station_inds = self.valid_station_inds[networkspeci][self.observations_data_label]  

        # get N stations for networkspeci
        self.n_stations = len(self.relevant_station_inds)
        
        # if have 0 relevant stations, continue to next networkspeci
        if self.n_stations == 0:
            self.logger.error('\nNo valid stations for {}, {}. Not making summmary plots'.format(networkspeci, self.subsection))
        else:
            self.logger.info('\nMaking {}, {} summary plots'.format(networkspeci, self.subsection)) 

        # create nested dictionary to store statistical information across all networkspecies
        if networkspeci not in self.stats_summary[self.subsection]:
            self.stats_summary[self.subsection][networkspeci] = {}

        # initialise the reference and DOI data for ACTRIS
        self.station_reference_data[self.subsection][networkspeci] = []
        self.station_doi_data[self.subsection][networkspeci] = []

        # get selected station data across all networkspecies for current section (if have not yet)
        if not self.made_networkspeci_summary_plots:
            
            # get selected station data
            get_selected_station_data(read_instance=self, canvas_instance=self, 
                                      networkspecies=self.networkspecies)

            # update data range min/maxes for summary paradigm 
            for ns in self.networkspecies:
                if self.selected_station_data_min[ns] < self.data_range_min_summary[ns]:
                    self.data_range_min_summary[ns] = copy.deepcopy(self.selected_station_data_min[ns]) 
                if self.selected_station_data_max[ns] > self.data_range_max_summary[ns]:
                    self.data_range_max_summary[ns] = copy.deepcopy(self.selected_station_data_max[ns])
                if self.selected_station_stddev_max[ns] > self.stddev_max_summary[ns]:
                    self.stddev_max_summary[ns] = copy.deepcopy(self.selected_station_stddev_max[ns])

        # if have no valid data across data labels (no observations or models), then set flag
        if not self.selected_station_data[networkspeci]: 
            have_nodata = True
        else:
            have_nodata = False

        # for actris, store current station reference and doi per month
        if networkspeci.split('|')[0] == 'actris/actris' and not have_nodata:
            for relevant_station_ind in self.relevant_station_inds:
                self.current_station_reference = self.station_references[networkspeci][relevant_station_ind]
                self.station_reference_data[self.subsection][networkspeci].append(self.current_station_reference)
                current_station_dois = self.get_dois_from_station_reference(self.current_station_reference, networkspeci)
                self.station_doi_data[self.subsection][networkspeci].append(current_station_dois)
            
        # iterate through plots to make
        for plot_type in summary_plots_to_make:

            # get zstat information from plot_type
            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type=plot_type)
            
            # get base plot type (without stat and options)
            if zstat:
                base_plot_type = plot_type.split('-')[0] 
            else:
                base_plot_type = plot_type.split('_')[0] 

            # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
            plot_options = plot_type.split('_')[1:]

            # make sure periodic, map, heatmap, taylor and table plots have a -[stat]
            if base_plot_type in ['periodic', 'map', 'heatmap', 'taylor', 'table'] and zstat is None:
                msg = f'{plot_type} plot needs a statistic -[stat].'
                show_message(self, msg)
                continue

            # check if plot type if a dataframe
            plot_type_df = self.get_plot_type_df(base_plot_type)

            # do not make plot, if plot type is not a dataframe or multispecies, and have no data
            if (not plot_type_df) & ('multispecies' not in plot_options) & (have_nodata):
                continue
            
            # update dictionary to store statistical information if plot type is a dataframe
            if plot_type_df:
                self.update_stats_tables('summary', base_plot_type, plot_type, zstat, networkspeci, plot_options)
                # do not make plot until last subsection
                if (self.subsection != self.subsections[-1]):
                    continue

            # for multispecies plots, continue in certain scenarios
            if ('multispecies' in plot_options):
                # for non-dataframe type plots, only plot in first subsection
                if (not plot_type_df) & (self.made_networkspeci_summary_plots):
                    continue
                # for dataframe plot types, only plot in last subsection and last networkspecies
                elif (plot_type_df) & (networkspeci != self.networkspecies[-1]):
                    continue

            # for timeseries chunking
            chunk_stat = None
            chunk_resolution = None
            if base_plot_type == 'timeseries':
                if zstat:
                    # get chunk statistic and resolution
                    chunk_stat = copy.deepcopy(zstat)
                    chunk_resolution = plot_type.split('-')[2].split('_')[0]

            self.logger.info('\nMaking summary {0}'.format(plot_type))

            if base_plot_type in ['fairmode-target', 'fairmode-statsummary']:
                # warning for fairmode plots if species aren't PM2.5, PM10, NO2 or O3
                speci = networkspeci.split('|')[1]
                if speci not in ['sconco3', 'sconcno2', 'pm10', 'pm2p5']:
                    msg = f'Fairmode plot cannot be created for {speci}.'
                    show_message(self, msg)
                    continue
                # warning for fairmode plots if resolution is not hourly or daily
                if ((speci in ['sconco3', 'sconcno2'] and self.active_resolution != 'hourly') 
                    or (speci in ['pm10', 'pm2p5'] and (self.active_resolution not in ['hourly', 'daily']))):
                    msg = 'Fairmode plot can only be created if the resolution is hourly (O3, NO2, PM2.5 and PM10) or daily (PM2.5 and PM10).'
                    show_message(self, msg)
                    continue

                # skip making plot if there is no valid data
                data, valid_station_idxs = get_fairmode_data(self, self, networkspeci, self.data_labels)
                if not any(valid_station_idxs):
                    self.logger.info(f'No data to generate FAIRMODE plot after filtering by coverage for {speci}.')
                    continue
            
            if base_plot_type == 'contingencytable':
                # warning for contingency tables if species aren't PM2.5, PM10, NO2, O3 or SO2
                speci = networkspeci.split('|')[1]
                if speci not in ['sconco3', 'sconcno2', 'pm10', 'pm2p5', 'sconcso2']:
                    msg = f'Contingency table cannot be created for {speci}.'
                    show_message(self, msg)
                    continue
                # warning for contingency tables if resolution is not hourly
                if self.active_resolution != 'hourly':
                    msg = 'Contingency table can only be created if the resolution is hourly.'
                    show_message(self, msg)
                    continue

            # make plot
            plot_indices = self.make_plot('summary', plot_type, plot_options, networkspeci)

            # do formatting for plot options
            relevant_axs, relevant_data_labels = self.get_relevant_axs_per_networkspeci_plot_type_page_ind(
                base_plot_type, plot_indices)
            format_plot_options(self, self, relevant_axs, relevant_data_labels, networkspeci, base_plot_type, 
                                plot_type, plot_options, map_extent=self.map_extent,
                                chunk_stat=chunk_stat, chunk_resolution=chunk_resolution)

        # update N total pages 
        self.n_total_pages = len(self.plot_dictionary)

        # update variable when summary plots have been made for a networkspecies
        self.made_networkspeci_summary_plots = True

    def make_station_plots(self, networkspeci, station_plots_to_make):
        """
        Generates all station-specific visualisations for a specific network-species combination.

        Parameters
        ----------
        networkspeci : str
            The combined network and species identifier.
        station_plots_to_make : list
            List of plot types to be generated for individual measurement stations.
        """

        # setup plotting geometry for station plots per networkspeci (for one subsection)
        if self.do_plot_geometry_setup:
            self.setup_plot_geometry('station', networkspeci, self.made_networkspeci_station_plots)

        # get valid station inds for networkspeci 
        if self.temporal_colocation:
            self.relevant_station_inds = self.valid_station_inds_temporal_colocation[networkspeci][self.observations_data_label]
        else:
            self.relevant_station_inds = self.valid_station_inds[networkspeci][self.observations_data_label]  

        # get N stations for networkspeci
        self.n_stations = len(self.relevant_station_inds)

        # if have 0 relevant stations, continue to next networkspeci
        if self.n_stations == 0:
            self.logger.info('No valid stations for {}, {}. Not making station plots'.format(networkspeci, self.subsection))
        else:
            self.logger.info('Making {}, {} station plots'.format(networkspeci, self.subsection)) 
        
        # create nested dictionary to store statistical information across all networkspecies
        if networkspeci not in self.stats_station:
            self.stats_station[self.subsection][networkspeci] = {}

        # initialise station ind as -1
        self.station_ind = -1

        # initialise or reset (if report_summary is True) the reference and DOI data for ACTRIS
        self.station_reference_data[self.subsection][networkspeci] = []
        self.station_doi_data[self.subsection][networkspeci] = []

        # iterate through stations
        for i, relevant_station_ind in enumerate(self.relevant_station_inds):
            
            # gather some information about current station
            self.station_ind += 1
            self.current_lon = round(self.station_longitudes[networkspeci][relevant_station_ind], 2)
            self.current_lat = round(self.station_latitudes[networkspeci][relevant_station_ind], 2)
            self.current_station_name = self.station_names[networkspeci][relevant_station_ind]
            self.current_station_reference = self.station_references[networkspeci][relevant_station_ind]

            # get selected station data for specific networkspeci and current section
            get_selected_station_data(read_instance=self, canvas_instance=self, 
                                      networkspecies=[networkspeci], 
                                      station_index=relevant_station_ind, 
                                      data_range_min=self.data_range_min_station, 
                                      data_range_max=self.data_range_max_station,
                                      stddev_max=self.stddev_max_station)

            # if have no valid data across data labels (no observations or models), then set flag
            if not self.selected_station_data[networkspeci]: 
                have_nodata = True
            else:
                have_nodata = False
            
            # update data range min/maxes for station paradigm
            if self.selected_station_data_min[networkspeci] < self.data_range_min_station[networkspeci]:
                self.data_range_min_station[networkspeci] = copy.deepcopy(self.selected_station_data_min[networkspeci]) 
            if self.selected_station_data_max[networkspeci] > self.data_range_max_station[networkspeci]:
                self.data_range_max_station[networkspeci] = copy.deepcopy(self.selected_station_data_max[networkspeci])
            if self.selected_station_stddev_max[networkspeci] > self.stddev_max_station[networkspeci]:
                self.stddev_max_station[networkspeci] = copy.deepcopy(self.selected_station_stddev_max[networkspeci])

            # for actris, store current station reference and doi per month
            if networkspeci.split('|')[0] == 'actris/actris' and not have_nodata:
                self.station_reference_data[self.subsection][networkspeci].append(self.current_station_reference)
                current_station_dois = self.get_dois_from_station_reference(self.current_station_reference, networkspeci)
                self.station_doi_data[self.subsection][networkspeci].append(current_station_dois)
            
            # iterate through plots to make
            for plot_type in station_plots_to_make:
                
                # get zstat information from plot_type
                zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type=plot_type)
                
                # get base plot type (without stat and options)
                if zstat:
                    base_plot_type = plot_type.split('-')[0] 
                else:
                    base_plot_type = plot_type.split('_')[0] 
                
                # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                plot_options = plot_type.split('_')[1:]

                # make sure periodic, map, heatmap, taylor and table plots have a -[stat]
                if base_plot_type in ['periodic', 'map', 'heatmap', 'taylor', 'table'] and zstat is None:
                    msg = f'{plot_type} plot needs a statistic -[stat].'
                    show_message(self, msg)
                    continue
                
                # check if plot type if a dataframe
                plot_type_df = self.get_plot_type_df(base_plot_type)

                # do not make plot, if plot type is not a dataframe or multispecies, and have no data
                if (not plot_type_df) & ('multispecies' not in plot_options) & (have_nodata):
                    continue

                # remove station measurement method code from references for GHOST data
                # if this has not already been removed
                if ('multispecies' in plot_options) & (self.reading_ghost) & ('_' in self.current_station_reference):
                    self.current_station_reference = '_'.join(self.current_station_reference.split("_")[:-1])

                # update dictionary to store statistical information if plot type is a dataframe
                if plot_type_df:
                    self.update_stats_tables('station', base_plot_type, plot_type, zstat, networkspeci, plot_options)

                # for multispecies plots, continue in certain scenarios
                if ('multispecies' in plot_options):
                    # spatial colocation needs to be turned on for multispecies station plots 
                    if not self.spatial_colocation:
                        continue
                    # for non-dataframe plot types, only plot in first subsection
                    elif (not plot_type_df) & (self.made_networkspeci_station_plots):
                        continue
                    # for dataframe plot types, only plot in last subsection and last networkspecies
                    elif (plot_type_df) & (networkspeci != self.networkspecies[-1]):
                        continue

                    # get selected station data for non-dataframe type plots across all networkspecies
                    if not plot_type_df:
                        get_selected_station_data(read_instance=self, canvas_instance=self, 
                                                  networkspecies=self.networkspecies, 
                                                  station_index=relevant_station_ind, 
                                                  data_range_min=self.data_range_min_station, 
                                                  data_range_max=self.data_range_max_station,
                                                  stddev_max=self.stddev_max_station)

                        # update data range min/maxes for station paradigm
                        for ns in self.networkspecies:
                            if self.selected_station_data_min[ns] < self.data_range_min_station[ns]:
                                self.data_range_min_station[ns] = copy.deepcopy(self.selected_station_data_min[ns]) 
                            if self.selected_station_data_max[ns] > self.data_range_max_station[ns]:
                                self.data_range_max_station[ns] = copy.deepcopy(self.selected_station_data_max[ns])
                            if self.selected_station_stddev_max[ns] > self.stddev_max_station[ns]:
                                self.stddev_max_station[ns] = copy.deepcopy(self.selected_station_stddev_max[ns])

                # for timeseries chunking
                chunk_stat = None
                chunk_resolution = None
                if base_plot_type == 'timeseries':
                    if zstat:
                        # get chunk statistic and resolution
                        chunk_stat = copy.deepcopy(zstat)
                        chunk_resolution = plot_type.split('-')[2].split('_')[0]

                self.logger.info('\nMaking station {2} for {3} ({0}/{1})'.format(i+1, 
                                                                    len(self.relevant_station_inds),
                                                                    plot_type, 
                                                                    self.current_station_name))   

                speci = networkspeci.split('|')[1]
                if base_plot_type in ['fairmode-target', 'fairmode-statsummary']:
                    # warning for fairmode plots if species aren't PM2.5, PM10, NO2 or O3
                    if speci not in ['sconco3', 'sconcno2', 'pm10', 'pm2p5']:
                        msg = f'Fairmode plot cannot be created for {speci} in {self.current_station_name}.'
                        show_message(self, msg)
                        continue
                    # warning for fairmode plots if resolution is not hourly or daily
                    if ((speci in ['sconco3', 'sconcno2'] and self.active_resolution != 'hourly') 
                        or (speci in ['pm10', 'pm2p5'] and (self.active_resolution not in ['hourly', 'daily']))):
                        msg = 'Fairmode plot can only be created if the resolution is hourly (O3, NO2, PM2.5 and PM10) or daily (PM2.5 and PM10).'
                        show_message(self, msg)
                        continue

                    # skip making plot if there is no valid data
                    data, valid_station_idxs = get_fairmode_data(self, self, networkspeci, self.data_labels)
                    if not any(valid_station_idxs):
                        self.logger.info(f'No data to generate FAIRMODE plot after filtering by coverage for {speci} in {self.current_station_name}.')
                        continue
                
                if base_plot_type == 'contingencytable':
                    # warning for contingency tables if species aren't PM2.5, PM10, NO2, O3 or SO2
                    if speci not in ['sconco3', 'sconcno2', 'pm10', 'pm2p5', 'sconcso2']:
                        msg = f'Contingency table cannot be created for {speci} in {self.current_station_name}.'
                        show_message(self, msg)
                        continue
                    # warning for contingency tables if resolution is not hourly
                    if self.active_resolution != 'hourly':
                        msg = 'Contingency table can only be created if the resolution is hourly.'
                        show_message(self, msg)
                        continue

                # make plot                
                plot_indices = self.make_plot('station', plot_type, plot_options, networkspeci)

                # do not format Taylor diagrams until last station
                if (base_plot_type == 'taylor'):
                    if relevant_station_ind != self.relevant_station_inds[-1]:
                        continue

                # do formatting for plot options
                relevant_axs, relevant_data_labels = self.get_relevant_axs_per_networkspeci_plot_type_page_ind(
                    base_plot_type, plot_indices)
                format_plot_options(self, self, relevant_axs, relevant_data_labels, networkspeci, base_plot_type, 
                                    plot_type, plot_options, map_extent=self.map_extent, 
                                    chunk_stat=chunk_stat, chunk_resolution=chunk_resolution)

        # update N total pages 
        self.n_total_pages = len(self.plot_dictionary)

        # update variable now station plots have been made for a networkspecies
        self.made_networkspeci_station_plots = True
    
    def get_plot_type_df(self, base_plot_type):   
        """
        Determines if a plot type requires data collation across subsections and networkspecies.

        Parameters
        ----------
        base_plot_type : str
            The primary identifier for the plot category (e.g. 'table', 'heatmap').

        Returns
        -------
        bool
            True if the plot type is a dataframe-based summary, False otherwise.
        """

        if base_plot_type in ['table', 'heatmap', 'statsummary']:
            plot_type_df = True
        else: 
            plot_type_df = False

        return plot_type_df

    def update_stats_tables(self, plotting_paradigm, base_plot_type, plot_type, zstat, networkspeci, plot_options):
        """
        Update statistical tables for plot types that collate data across various networkspecies/subsections.

        Parameters
        ----------
        plotting_paradigm : str
            The context of the plots, either 'summary' or 'station'.
        base_plot_type : str
            The primary identifier for the plot category (e.g. 'table', 'heatmap').
        plot_type : str
            The full plot type string including options and statistics.
        zstat : str
            The specific statistic to be calculated.
        networkspeci : str
            The combined network and species identifier.
        plot_options : list
            List of configuration options parsed from the plot type string.
        """

        # set statistics to calculate based on plot type
        if base_plot_type in ['heatmap', 'table']:
            stats = [zstat]
        elif base_plot_type == 'statsummary':
            if 'bias' in plot_options:
                stats = self.plot_characteristics[plot_type]['model_bias']
            else:
                stats = self.plot_characteristics[plot_type]['basic']

        # iterate through all data labels
        for data_label in self.data_labels:

            # create nested dictionary to store statistical information across all subsections
            for stat in stats:

                # get zstat information 
                zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=stat)

                # skip observations data label when plotting bias
                if (data_label == self.observations_data_label) & (z_statistic_sign == 'bias'):
                    continue

                if plotting_paradigm == 'summary':
                    if stat not in self.stats_summary[self.subsection][networkspeci]:
                        self.stats_summary[self.subsection][networkspeci][stat] = {}
                elif plotting_paradigm == 'station':
                    if self.current_station_reference not in self.stats_station[self.subsection][networkspeci]:
                        self.stats_station[self.subsection][networkspeci][self.current_station_reference] = {}
                    if stat not in self.stats_station[self.subsection][networkspeci][self.current_station_reference]:
                        self.stats_station[self.subsection][networkspeci][self.current_station_reference][stat] = {}
                
                # get stat for current data label
                if data_label in self.selected_station_data_labels[networkspeci]:
                    # if relevant stat is modbias stat, then ensure temporal colocation is active
                    if (base_plot_type == 'statsummary') and (stat in self.modbias_stats) and ((not self.temporal_colocation) or (len(self.data_labels) == 1)):
                        data_to_add = np.nan
                    # otherwise calculate statistic
                    else:
                        if z_statistic_sign == 'bias':
                            data_to_add = calculate_statistic(self, self, networkspeci, zstat, [self.observations_data_label], [data_label])
                        else:
                            data_to_add = calculate_statistic(self, self, networkspeci, zstat, [data_label], [])
                else:
                    data_to_add = np.nan
                
                # add data to dicts
                if plotting_paradigm == 'summary':
                    self.stats_summary[self.subsection][networkspeci][stat][data_label] = data_to_add
                elif plotting_paradigm == 'station':
                    self.stats_station[self.subsection][networkspeci][self.current_station_reference][stat][data_label] = data_to_add
            
        return None

    def make_plot(self, plotting_paradigm, plot_type, plot_options, networkspeci):
        """
        Directs the generation of specific plot types by coordinating data labels, axes, and plotting functions.

        Parameters
        ----------
        plotting_paradigm : str
            The context of the plots, either 'summary' or 'station'.
        plot_type : str
            The full plot type string including options and statistics.
        plot_options : list
            List of configuration options parsed from the plot type string.
        networkspeci : str
            The combined network and species identifier.

        Returns
        -------
        list
            A list of indices [page_number, page_ind] where plots were successfully created.
        """

        current_plot_ind = 0

        # create list to store index of saved plot information for plot_type
        # index is composed of nested list of [page_number, page_ind]
        plot_indices = []

        # get zstat information from plot_type
        zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type=plot_type)

        # get base plot type (without stat and options)
        if zstat:
            base_plot_type = plot_type.split('-')[0] 
        else:
            base_plot_type = plot_type.split('_')[0] 

        # set variable to know if we need to create plot for dataframes in last subsection
        plot_type_df = self.get_plot_type_df(base_plot_type)

        # get all data_labels for selected_station_data
        # if multispecies plot then ensure have data_labels from across all networkspecies
        if plot_type_df:

            if base_plot_type in ['heatmap', 'table']:                
                if z_statistic_sign == 'bias':
                    data_labels = list(self.experiments.values())
                else:
                    data_labels = [self.observations_data_label] + list(self.experiments.values())

            elif base_plot_type == 'statsummary':
                if 'bias' in plot_options:
                    data_labels = list(self.experiments.values())
                else:
                    data_labels = [self.observations_data_label] + list(self.experiments.values())

        elif 'multispecies' in plot_options: 
            data_labels = []
            for ns in self.selected_station_data:
                data_labels.extend(self.selected_station_data_labels[ns])
            _, idx = np.unique(data_labels, return_index=True)
            data_labels = np.array(data_labels)[np.sort(idx)].tolist()
        else:
            data_labels = copy.deepcopy(self.data_labels)

        # if have no valid data labels then return 
        if len(data_labels) == 0:
            msg = "Cannot make plot as have no valid data. Not making plot."
            show_message(self, msg)
            return plot_indices

        # if are making bias plot, and have no valid model data then cannot make plot type
        if (('bias' in plot_options) or (z_statistic_sign == 'bias')) & (len(data_labels) < 2) & (not plot_type_df):
            msg = "Cannot make a bias plot with no model data. Not making plot."
            show_message(self, msg)
            return plot_indices

        # do not make plot if bias and threshold plots are in plot options
        if ('bias' in plot_options) & ('threshold' in plot_options):
            msg = "Cannot make a bias plot showing threshold lines. Not making plot."
            show_message(self, msg)
            return plot_indices

        # do not make plot if hidedata is active but smooth is not in plot options
        if (base_plot_type == 'timeseries') and ('hidedata' in plot_options) and ('smooth' not in plot_options):
            msg = f"Cannot make {plot_type} because 'hidedata' plot option is set for "
            msg += "timeseries plot, but 'smooth' is not active. Not making plot."
            show_message(self, msg)
            return plot_indices
        
        # do not make plot if hidedata is active but regression is not in plot options
        if (base_plot_type == 'scatter') and ('hidedata' in plot_options) and ('regression' not in plot_options):
            msg = f"Cannot make {plot_type} because 'hidedata' plot option is set for "
            msg += "scatter lot, but 'regression' is not active. Not making plot."
            show_message(self, msg)
            return plot_indices
        
        # do not make Taylor diagram if statistic is not r or r2
        if (base_plot_type == 'taylor') and (zstat not in ['r', 'r2']):
            msg = f"Cannot make {plot_type} because statistic is not available or defined. "
            msg += "Choose between 'taylor-r' or 'taylor-r2'. Not making plot."
            show_message(self, msg)
            return plot_indices

        # do not make periodic plot if stat is MDA8
        if (base_plot_type == 'periodic') and (base_zstat == 'MDA8'):
            msg = f"Cannot make {plot_type} because MDA8 statistic is not available for periodic plots. Not making plot."
            show_message(self, msg)
            return plot_indices
        
        # do not make statsummary plot if stat is MDA8, and are making periodic statistic
        if (base_plot_type == 'statsummary') and (base_zstat == 'MDA8') and (z_statistic_period is not None):
            msg = f"Cannot make {plot_type} because MDA8 statistic is not available for periodic statistics. Not making plot."
            show_message(self, msg)
            return plot_indices

        if base_plot_type == 'contingencytable':
            if ((not self.temporal_colocation) 
                or ((self.temporal_colocation) and (len(data_labels) == 1 or len(data_labels) > 2))):
                # do not make contingencytable if temporal colocation if off
                if not self.temporal_colocation:
                    msg = f'Cannot make {plot_type} plots without activating the temporal colocation.'
                # do not make contingencytable if no model is loaded
                elif len(data_labels) == 1:
                    msg = f'Cannot make {plot_type} plots without loading models.'
                # do not make contingencytable if more than one model is loaded
                elif len(data_labels) > 2:
                    msg = f'Cannot make {plot_type} plots with more than 1 model.'
                show_message(self, msg)
                return plot_indices
                
        # get data labels without observations
        data_labels_sans_obs = copy.deepcopy(data_labels)
        if self.observations_data_label in data_labels_sans_obs:
            data_labels_sans_obs.remove(self.observations_data_label)

        # determine if have some data to plot
        plot_validity = False
        if ((base_plot_type in ['scatter', 'fairmode-target', 'fairmode-statsummary']) or ('bias' in plot_options) or 
            (z_statistic_sign == 'bias')):
            data_labels_to_test = copy.deepcopy(data_labels_sans_obs)
        else:
            data_labels_to_test = copy.deepcopy(data_labels)
        if plot_type_df:
            plot_validity = True
        elif 'multispecies' in plot_options:
            for ns in self.selected_station_data:
                for data_label in data_labels_to_test:
                    if data_label in self.selected_station_data_labels[ns]:
                        plot_validity = True
        else:
            for data_label in data_labels_to_test:
                if data_label in self.selected_station_data_labels[networkspeci]:
                    plot_validity = True
        if not plot_validity:
            msg = f'{plot_type} cannot be created because there is no available data.'
            show_message(self, msg)
            return plot_indices

        # get data ranges for plotting paradigm
        if plotting_paradigm == 'summary':
            data_range_min = self.data_range_min_summary[networkspeci]
            data_range_max = self.data_range_max_summary[networkspeci]
            stddev_max = self.stddev_max_summary[networkspeci]
        elif plotting_paradigm == 'station':
            data_range_min = self.data_range_min_station[networkspeci]
            data_range_max = self.data_range_max_station[networkspeci]
            stddev_max = self.stddev_max_station[networkspeci]

        # map plots (1 plot per data array/s (1 array if absolute plot,
        # 2 arrays if making bias plot), per subsection)
        if base_plot_type == 'map':
            
            # get necessary data labels to plot
            if z_statistic_sign == 'bias':
                z1 = [self.observations_data_label] * len(data_labels_sans_obs)
                z2 = data_labels_sans_obs
            elif 'obs' in plot_options:
                z1 = [self.observations_data_label]
                z2 = ['']
            else:
                z1 = data_labels
                z2 = ['']*len(data_labels)

            # iterate through relevant data labels making plots
            for z1_label, z2_label in zip(z1, z2):

                # skip map if we have no data for a specific label
                skip_map = False
                if (z2_label == '') and (z1_label not in self.selected_station_data_labels[networkspeci]):
                    skip_map = True
                elif ((z2_label != '') and ((z1_label not in self.selected_station_data_labels[networkspeci])
                                            or (z2_label not in self.selected_station_data_labels[networkspeci]))):
                    skip_map = True

                if skip_map:
                    if (z2_label == ''):
                        unavailable_label = '{}'.format(z1_label)
                    else:
                        unavailable_label = '{} - {}'.format(z2_label, z1_label)
                    msg = f'{plot_type} cannot be created because there is no available data of {unavailable_label}.'
                    show_message(self, msg)
                    return plot_indices
                    
                # get relevant page/axis to plot on
                axis_ind = (current_plot_ind * len(self.subsections)) + self.subsection_ind
                relevant_page, page_ind, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, 
                                                                                     plot_type, axis_ind)

                # set axis title
                if relevant_axis.get_title() == '':
                    if z2_label != '':
                        label = copy.deepcopy(z2_label)
                    else:
                        label = copy.deepcopy(z1_label) 

                    axis_title_label = '{}\n{} '.format(label, self.subsection)
                    if self.n_stations == 1:
                        axis_title_label += '(1 station)'
                    else:
                        axis_title_label += '({} stations)'.format(self.n_stations)
                    set_axis_title(self, relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

                # set map extent ? 
                if self.map_extent:
                    set_map_extent(self, relevant_axis, self.map_extent)

                # make map plot
                self.plotting.make_map(relevant_axis, networkspeci, self.plot_characteristics[plot_type], plot_options,
                                       zstat=zstat, labela=z1_label, labelb=z2_label)

                # save plot information for later formatting 
                if z2_label == '':
                    self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(z1_label)
                else:
                    self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(z2_label)
                plot_index = [relevant_page, page_ind]
                if plot_index not in plot_indices:
                    plot_indices.append(plot_index)

                # turn axis on
                relevant_axis.axis('on')
                relevant_axis.set_visible(True)

                # iterate number of plots made for current type of plot 
                current_plot_ind += 1     

        # other plots (1 plot per subsection with multiple data arrays for summary paradigm, 1 plot per subsection per station for station paradigm)
        elif base_plot_type not in ['heatmap', 'table', 'statsummary']:

            # if making individual plots, iterate through data labels one at a time, otherwise pass all data labels 
            # together
            if 'individual' in plot_options:
                iter_data_labels = copy.deepcopy(data_labels)
            else:
                iter_data_labels = [data_labels]
            
            # for timeseries chunking
            if base_plot_type == 'timeseries':
                if zstat:
                    # get chunk statistic and resolution
                    chunk_stat = copy.deepcopy(zstat)
                    chunk_resolution = plot_type.split('-')[2].split('_')[0]
                    
                    # get available chunk timeseries resolutions
                    available_timeseries_chunk_resolutions = get_possible_resampling_resolutions(self.active_resolution,
                                                                                                 daily_forecast=self.daily_forecast)

                    # show warning if it is not available
                    if chunk_resolution not in available_timeseries_chunk_resolutions:
                        msg = f'{plot_type} cannot be created because {chunk_resolution} '
                        msg += 'is not an available chunking resolution. '
                        if len(available_timeseries_chunk_resolutions) > 0:
                            msg += f'The available resolutions are: {available_timeseries_chunk_resolutions}'
                        show_message(self, msg)
                        return plot_indices

                    # show warning if chunk stat is MDA8 and active resolution is not hourly
                    if (chunk_stat == 'MDA8') and (self.active_resolution != 'hourly'):
                        msg = f'{plot_type} cannot be created because {chunk_stat} '
                        msg += 'is only available when active resolution is hourly.'
                        show_message(self, msg)
                        return plot_indices
                    
                    # show warning if chunk stat is MDA8 and chunk resolution is not daily
                    if (chunk_stat == 'MDA8') and (chunk_resolution != 'daily'):
                        msg = f'{plot_type} cannot be created because {chunk_stat} '
                        msg += 'is only available when chunk resolution is daily.'
                        show_message(self, msg)
                        return plot_indices

                else:
                    chunk_stat = None
                    chunk_resolution = None
            
            for data_labels in iter_data_labels:

                # skip observations
                if ((data_labels == self.observations_data_label) 
                    and ((base_plot_type in ['scatter', 'fairmode-target', 'fairmode-statsummary']) or ('bias' in plot_options) or 
                    (z_statistic_sign == 'bias'))):
                    continue

                # skip individual plots if we have no data for a specific label
                if 'individual' in plot_options:
                    if data_labels not in self.selected_station_data_labels[networkspeci]:
                        msg = f'{plot_type} cannot be created because there is no available data.'
                        show_message(self, msg)
                        return plot_indices

                if type(data_labels) != list:
                    data_labels = [data_labels]
                    
                # get relevant axis to plot on
                if plotting_paradigm == 'summary':
                    if 'individual' in plot_options:
                        if ((base_plot_type in ['scatter', 'fairmode-target', 'fairmode-statsummary']) or ('bias' in plot_options) or 
                            (z_statistic_sign == 'bias')):
                                axis_ind = (current_plot_ind + self.subsection_ind + (len(self.experiments) - 1) * self.subsection_ind)
                        else:
                            axis_ind = (current_plot_ind + self.subsection_ind + len(self.experiments) * self.subsection_ind)
                    else:
                        axis_ind = self.subsection_ind
                elif plotting_paradigm == 'station':
                    if 'individual' in plot_options:
                        if ((base_plot_type in ['scatter', 'fairmode-target', 'fairmode-statsummary']) or ('bias' in plot_options) or 
                            (z_statistic_sign == 'bias')):
                            axis_ind = (current_plot_ind + self.station_ind + (len(self.experiments) - 1) * self.station_ind)
                        else:
                            axis_ind = (current_plot_ind + self.station_ind + len(self.experiments) * self.station_ind)
                    else:
                        axis_ind = self.station_ind
                
                # get relevant axis
                relevant_page, page_ind, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, 
                                                                                     plot_type, axis_ind)

                # set axis title (only if not previously set)
                if isinstance(relevant_axis, dict):
                    for relevant_temporal_resolution, sub_ax in relevant_axis.items():
                        if relevant_temporal_resolution == 'hour':
                            axis_title = sub_ax.get_title()
                            break
                elif isinstance(relevant_axis, list):
                    axis_title = relevant_axis[0].get_title()
                else:
                    axis_title = relevant_axis.get_title()

                # set xlabel and ylabel (only if not previously set)
                if isinstance(relevant_axis, dict):
                    for relevant_temporal_resolution, sub_ax in relevant_axis.items():
                        if relevant_temporal_resolution in ['hour','month']:
                            axis_xlabel = 'NaN'
                            axis_ylabel = sub_ax.get_ylabel()
                            break
                elif isinstance(relevant_axis, list):
                    axis_xlabel = ""
                    axis_ylabel = ""
                else:
                    axis_xlabel = relevant_axis.get_xlabel()
                    axis_ylabel = relevant_axis.get_ylabel()

                # axis xlabel is empty?
                if (axis_xlabel == '') or ('[measurement_units]' in axis_xlabel):
                    if 'xlabel' in self.plot_characteristics[plot_type]:
                        xlabel = self.plot_characteristics[plot_type]['xlabel']['xlabel']
                        if '[measurement_units]' in xlabel:
                            xlabel = xlabel.replace('[measurement_units]', '[{}]'.format(self.measurement_units[networkspeci.split('|')[-1]]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                    else:
                        xlabel = ''
                    # set xlabel
                    if xlabel != '':
                        set_axis_label(relevant_axis, 'x', xlabel, self.plot_characteristics[plot_type])

                # axis ylabel is empty?
                if (axis_ylabel == '') or ('[measurement_units]' in axis_ylabel):
                    if base_plot_type in ['periodic'] or ((base_plot_type == 'timeseries') 
                                                          and (chunk_stat is not None) 
                                                          and (chunk_resolution is not None)):
                        if z_statistic_type == 'basic':
                            ylabel = self.basic_stats[base_zstat]['label']
                            ylabel_units = self.basic_stats[base_zstat]['units']
                        else:
                            ylabel = self.modbias_stats[base_zstat]['label']
                            ylabel_units = self.modbias_stats[base_zstat]['units']
                        if ylabel_units == '[measurement_units]':
                            ylabel_units = self.measurement_units[networkspeci.split('|')[-1]] 
                        if ylabel_units != '':
                            ylabel += ' [{}]'.format(ylabel_units)
                    else:
                        if 'ylabel' in self.plot_characteristics[plot_type]:
                            ylabel = self.plot_characteristics[plot_type]['ylabel']['ylabel']
                            if '[measurement_units]' in ylabel:
                                ylabel = ylabel.replace('[measurement_units]', '[{}]'.format(self.measurement_units[networkspeci.split('|')[-1]]))
                        else:
                            ylabel = ''
                    # set ylabel
                    if ylabel != '':
                        set_axis_label(relevant_axis, 'y', ylabel, self.plot_characteristics[plot_type])

                # get plotting function
                if base_plot_type in ['fairmode-target', 'fairmode-statsummary']:
                    func = getattr(self.plotting, 'make_{}'.format(base_plot_type.replace('-','_')))
                else:
                    func = getattr(self.plotting, 'make_{}'.format(base_plot_type.split('-')[0]))

                if base_plot_type == 'periodic':
                    func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         plot_options, zstat=zstat)    
                elif base_plot_type == 'distribution':
                    func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         plot_options, data_range_min=data_range_min, data_range_max=data_range_max) 
                elif base_plot_type == 'taylor':
                    func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         plot_options, zstat=zstat, stddev_max=stddev_max)
                elif base_plot_type == 'timeseries':
                    func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         plot_options, chunk_stat=chunk_stat, chunk_resolution=chunk_resolution)   
                else:
                    func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         plot_options) 

                # axis title is empty?
                # moved after plots have been created because in the case of the Taylor diagram 
                # the axis ghelper needs to be updated to be able to add the axis title
                if axis_title == '':
                    if plotting_paradigm == 'summary':
                        if self.n_stations == 1:
                            axis_title_label = '{} ({} station)'.format(self.subsection, self.n_stations)
                        else:
                            axis_title_label = '{} ({} stations)'.format(self.subsection, self.n_stations)
                    elif plotting_paradigm == 'station':
                        if base_plot_type == 'metadata':
                            axis_title_label = ''
                        else:
                            axis_title_label = '{}, {} ({:.{}f}, {:.{}f})'.format(self.current_station_reference,
                                                                                  self.current_station_name, 
                                                                                  self.current_lon,
                                                                                  self.plot_characteristics[plot_type]['round_decimal_places']['title'],
                                                                                  self.current_lat,
                                                                                  self.plot_characteristics[plot_type]['round_decimal_places']['title'])
                    
                    if base_plot_type in ['fairmode-target', 'fairmode-statsummary']:
                        speci = networkspeci.split('|')[1]
                        axis_title_label += '\n{}'.format(fairmode_settings[speci]['title'])

                    # set title
                    set_axis_title(self, relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

                # save plot information for later formatting
                self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].extend(data_labels)
                plot_index = [relevant_page, page_ind]
                if plot_index not in plot_indices:
                    plot_indices.append(plot_index)

                # turn axis/axes on
                if base_plot_type in ['periodic','periodic-violin']:
                    for relevant_temporal_resolution in self.periodic_relevant_temporal_resolutions:
                        relevant_axis[relevant_temporal_resolution].axis('on')
                        relevant_axis[relevant_temporal_resolution].set_visible(True)

                        # get references to periodic label annotations made, and then show them
                        annotations = [child for child in relevant_axis[relevant_temporal_resolution].get_children() 
                                        if isinstance(child, matplotlib.text.Annotation)]
                        # hide annotations
                        for annotation in annotations:
                            annotation.set_visible(True)
                
                elif base_plot_type == "fairmode-statsummary":
                    for i in range(len(relevant_axis)):
                        relevant_axis[i].axis('on')
                        relevant_axis[i].set_visible(True)
                else:
                    relevant_axis.axis('on')
                    relevant_axis.set_visible(True)

                # iterate number of plots made for current type of plot 
                current_plot_ind += 1     

        # make heatmap / table / statsummary plot
        # heatmap / table is one stat per networkspecies, subsections and data labels (columns are data labels, rows are networkspecies / subsections)
        # statsummary is multiple stats per networkspecies, subsections and data labels (columns are stats, rows are networkspecies / subsections / data labels)
        elif base_plot_type in ['heatmap', 'table', 'statsummary']:
            
            # get relevant axis to plot on
            axis_networkspeci = networkspeci
            if plotting_paradigm == 'summary':
                if base_plot_type == 'statsummary':
                    axis_ind = self.subsection_ind
                elif base_plot_type in ['heatmap', 'table']: 
                    axis_ind = 0
                if 'multispecies' in plot_options:
                    axis_networkspeci = list(self.summary_pages[plot_type].keys())[0]
            elif plotting_paradigm == 'station':
                axis_ind = self.station_ind
                if 'multispecies' in plot_options:
                    axis_networkspeci = list(self.station_pages[plot_type].keys())[0]
            relevant_page, page_ind, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, axis_networkspeci, 
                                                                                 plot_type, axis_ind)

            # convert stats_summary and stats_station dicts to dataframes
            if plotting_paradigm == 'summary':
                stats_to_plot = copy.deepcopy(self.stats_summary)
            elif plotting_paradigm == 'station':
                stats_to_plot = copy.deepcopy(self.stats_station)

            # get data for all subsections for summary
            if plotting_paradigm == 'summary':
                subsections = stats_to_plot.keys()
            # get data only for current subsection
            elif plotting_paradigm == 'station':
                subsections = [self.subsection]

            # get multiple networkspecies for multispecies (used in heatmaps, tables and statsummaries)
            if 'multispecies' in plot_options:
                networkspecies = self.networkspecies
            # get unique networkspeci
            else:
                networkspecies = [networkspeci]

            if base_plot_type in ['heatmap', 'table']:

                # create empty dataframe with networkspecies and subsections
                index = pd.MultiIndex.from_product([networkspecies, self.subsections],
                                                    names=["networkspecies", "subsections"])
                stats_df = pd.DataFrame(np.nan, index=index, columns=data_labels, dtype=np.float64)
                
                # convert stats_summary and stats_station dicts to dataframes
                for subsection in subsections:
                    for networkspeci in networkspecies:
                        stats_per_data_label = []
                        for data_label in data_labels:
                            # initialise stat with nan
                            stat_to_append = np.nan
                            # update stat
                            if networkspeci in stats_to_plot[subsection]:
                                if plotting_paradigm == 'summary':
                                    if zstat in stats_to_plot[subsection][networkspeci]:
                                        if data_label in stats_to_plot[subsection][networkspeci][zstat]:
                                            stat_to_append = stats_to_plot[subsection][networkspeci][zstat][data_label]
                                elif plotting_paradigm == 'station':
                                    if self.current_station_reference in stats_to_plot[subsection][networkspeci]:
                                        if zstat in stats_to_plot[subsection][networkspeci][self.current_station_reference]:
                                            if data_label in stats_to_plot[subsection][networkspeci][self.current_station_reference][zstat]:
                                                stat_to_append = stats_to_plot[subsection][networkspeci][self.current_station_reference][zstat][data_label]
                            stats_per_data_label.append(stat_to_append)

                        # get floats instead of arrays with 1 element each and save
                        stats_per_data_label = [stat_per_data_label[0] 
                                                if isinstance(stat_per_data_label, np.ndarray) 
                                                else stat_per_data_label 
                                                for stat_per_data_label in stats_per_data_label]
                        stats_df.loc[(networkspeci, subsection)] = stats_per_data_label
                
            elif base_plot_type == 'statsummary':

                # get stats
                if 'bias' in plot_options:
                    stats = self.plot_characteristics[plot_type]['model_bias']
                else:
                    stats = self.plot_characteristics[plot_type]['basic']

                # create empty dataframe with networkspecies and subsections
                index = pd.MultiIndex.from_product([self.networkspecies, self.subsections, data_labels],
                                                    names=["networkspecies", "subsections", "labels"])
                stats_df = pd.DataFrame(np.nan, index=index, columns=stats, dtype=np.float64)

                # convert stats_summary and stats_station dicts to dataframes
                for subsection in subsections:
                    for networkspeci in networkspecies:
                        for data_label in data_labels:
                            stats_per_data_label = []
                            for stat in stats:
                                # initialise stat as nan
                                stat_to_append = np.nan
                                # update stat
                                if networkspeci in stats_to_plot[subsection]:
                                    if plotting_paradigm == 'summary':
                                        if stat in stats_to_plot[subsection][networkspeci]:
                                            if data_label in stats_to_plot[subsection][networkspeci][stat]:
                                                stat_to_append = stats_to_plot[subsection][networkspeci][stat][data_label]
                                    elif plotting_paradigm == 'station':
                                        if self.current_station_reference in stats_to_plot[subsection][networkspeci]:
                                            if stat in stats_to_plot[subsection][networkspeci][self.current_station_reference]:
                                                if data_label in stats_to_plot[subsection][networkspeci][self.current_station_reference][stat]:
                                                    stat_to_append = stats_to_plot[subsection][networkspeci][self.current_station_reference][stat][data_label]
                                stats_per_data_label.append(stat_to_append)

                            # get floats instead of arrays with 1 element each and save
                            stats_per_data_label = [stat_per_data_label[0] 
                                                    if isinstance(stat_per_data_label, np.ndarray) 
                                                    else stat_per_data_label 
                                                    for stat_per_data_label in stats_per_data_label]
                            stats_df.loc[(networkspeci, subsection, data_label)] = stats_per_data_label

            # turn on relevant axis if dataframe has values or not all NaN
            if (len(stats_df.index) > 0) & (not stats_df.isnull().values.all()):
                
                if base_plot_type == 'statsummary':
                    # make statsummary
                    func = getattr(self.plotting, 'make_table')
                    func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         plot_options, stats, statsummary=True, subsection=self.subsection, 
                         plotting_paradigm=plotting_paradigm, stats_df=stats_df)
                else:
                    # make table/heatmap
                    func = getattr(self.plotting, 'make_{}'.format(base_plot_type))
                    func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         plot_options, zstat, subsection=self.subsection, plotting_paradigm=plotting_paradigm, 
                         stats_df=stats_df)
                
                # save plot information for later formatting
                self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].extend(data_labels)
                plot_index = [relevant_page, page_ind]
                if plot_index not in plot_indices:
                    plot_indices.append(plot_index)

                # turn axis on
                relevant_axis.axis('on')
                relevant_axis.set_visible(True)

            # if have no data then do not make plot
            else:
                msg = f'{plot_type} cannot be created because there is no available data.'
                show_message(self, msg)
                return plot_indices

            # set axis title
            if relevant_axis.get_title() == '':
                if plotting_paradigm == 'station':
                    axis_title_label = '{}, {} ({:.{}f}, {:.{}f})'.format(self.current_station_reference,
                                                                          self.current_station_name, 
                                                                          self.current_lon,
                                                                          self.plot_characteristics[plot_type]['round_decimal_places']['title'],
                                                                          self.current_lat,
                                                                          self.plot_characteristics[plot_type]['round_decimal_places']['title'])
                    set_axis_title(self, relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

        return plot_indices

    def get_relevant_page_axis(self, plotting_paradigm, networkspeci, plot_type, axis_ind):
        """
        Retrieves the specific page number, page index, and axis handle for a plot.

        Parameters
        ----------
        plotting_paradigm : str
            The context of the plots, either 'summary' or 'station'.
        networkspeci : str
            The combined network and species identifier.
        plot_type : str
            The full plot type string including options and statistics.
        axis_ind : int
            The sequential index of the plot being created for the given type.

        Returns
        -------
        int
            The global page number in the plot dictionary.
        int
            The index of the axis on that specific page.
        matplotlib.axes.Axes
            The handle to the Matplotlib axis (or dictionary of handles for complex grids).
        """

        # get axes associated with plot type
        if plotting_paradigm == 'summary':
            relevant_pages = self.summary_pages[plot_type][networkspeci]
        elif plotting_paradigm == 'station':
            relevant_pages = self.station_pages[plot_type][networkspeci][self.subsection]
            
        all_relevant_pages = []
        relevant_axes = []     
        page_inds = []
        for relevant_page in relevant_pages:
            relevant_axes.extend(self.plot_dictionary[relevant_page]['axs'])
            all_relevant_pages.extend([relevant_page]*len(self.plot_dictionary[relevant_page]['axs']))
            page_inds.extend(list(range(len(self.plot_dictionary[relevant_page]['axs']))))

        return all_relevant_pages[axis_ind], page_inds[axis_ind], relevant_axes[axis_ind]['handle']

    def get_relevant_axs_per_networkspeci_plot_type(self, base_plot_type, relevant_pages):
        """
        Retrieves all axis handles and associated data labels for a specific plot type across multiple pages.

        Parameters
        ----------
        base_plot_type : str
            The primary identifier for the plot category.
        relevant_pages : list
            List of page numbers to search for axes.

        Returns
        -------
        relevant_axs : list
            List of Matplotlib axis handles found on the specified pages.
        relevant_data_labels : list
            List of data label sets associated with each retrieved axis.
        """
       
        relevant_axs = []
        relevant_data_labels = []
        for relevant_page in relevant_pages:
            if base_plot_type in ['periodic', 'periodic-violin']:
                for ax in self.plot_dictionary[relevant_page]['axs']:
                    for relevant_temporal_resolution in self.periodic_relevant_temporal_resolutions:
                        relevant_axs.append(ax['handle'][relevant_temporal_resolution])
                        relevant_data_labels.append(ax['data_labels'])
            elif base_plot_type == 'fairmode-statsummary':
                for ax in self.plot_dictionary[relevant_page]['axs']:
                    for i in range(len(ax)):
                        relevant_axs.append(ax['handle'][i])
                    relevant_data_labels.append(ax['data_labels'])
            else:
                for ax in self.plot_dictionary[relevant_page]['axs']:
                    relevant_axs.append(ax['handle'])
                    relevant_data_labels.append(ax['data_labels'])

        return relevant_axs, relevant_data_labels
    
    def get_relevant_axs_per_networkspeci_plot_type_page_ind(self, base_plot_type, plot_indices):
        """
        Retrieves specific axis handles and data labels based on provided page and subplot indices.

        Parameters
        ----------
        base_plot_type : str
            The primary identifier for the plot category.
        plot_indices : list
            List of nested [page_number, page_ind] pairs identifying specific subplots.

        Returns
        -------
        relevant_axs : list
            List of Matplotlib axis handles corresponding to the requested indices.
        relevant_data_labels : list
            List of data label sets associated with each retrieved axis.
        """

        relevant_axs = []
        relevant_data_labels = []
        for relevant_page, page_ind in plot_indices:
            if base_plot_type in ['periodic', 'periodic-violin']:
                for relevant_temporal_resolution in self.periodic_relevant_temporal_resolutions:
                    relevant_axs.append(self.plot_dictionary[relevant_page]['axs'][page_ind]['handle'][relevant_temporal_resolution])
                    relevant_data_labels.append(self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'])
            else:
                relevant_axs.append(self.plot_dictionary[relevant_page]['axs'][page_ind]['handle'])
                relevant_data_labels.append(self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'])

        return relevant_axs, relevant_data_labels

def main(**kwargs):
    """
    Main entry point for generating offline reports.

    Parameters
    ----------
    **kwargs : dict
        Configuration arguments and command-line parameters passed to the Report initialiser.
    """
   
    report = Report(**kwargs)
    report.run()
