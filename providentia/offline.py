import os
import sys
import json
import copy

import datetime
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import mpl_toolkits.axisartist.floating_axes as fa
from matplotlib.projections import PolarAxes

from .read import DataReader
from .filter import DataFilter
from .plot import Plot
from .statistics import get_selected_station_data, calculate_statistic, generate_colourbar, get_z_statistic_info
from .configuration import ProvConfiguration
from providentia import aux

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))

class ProvidentiaOffline:
    """ Class to create Providentia offline reports. """

    # make sure that we are not using Qt5 backend with matplotlib
    matplotlib.use('Agg')

    def __init__(self, **kwargs):
        print("Starting Providentia offline...")

        # initialise default configuration variables
        # modified by commandline arguments, if given
        provconf = ProvConfiguration(self, **kwargs)

        # update self with command line arguments
        self.commandline_arguments = copy.deepcopy(kwargs)

        # update variables from config file
        if ('config' in kwargs) and (os.path.exists(kwargs['config'])):
            aux.load_conf(self, kwargs['config'])
            self.from_conf = True
        elif ('config' in kwargs) and (not os.path.exists(kwargs['config'])):     
            error = 'Error: The path to the configuration file specified in the command line does not exist.'
            sys.exit(error)
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            sys.exit(error)

        # load report plot presets
        self.report_plots = json.load(open(os.path.join(CURRENT_PATH, 'conf/report_plots.json')))

        # create dictionary of all available observational GHOST data
        self.all_observation_data = aux.get_ghost_observational_tree(self)

        # load dictionary with non-GHOST esarchive files to read
        nonghost_observation_data_json = json.load(open(os.path.join(CURRENT_PATH, 'conf/nonghost_files.json')))
        # merge to existing GHOST observational data dict if we have the path
        if self.nonghost_root is not None:
            nonghost_observation_data = aux.get_nonghost_observational_tree(self, nonghost_observation_data_json)
            self.all_observation_data = {**self.all_observation_data, **nonghost_observation_data}

        # initialise DataReader class
        self.datareader = DataReader(self)

        # iterate through configuration sections
        for section_ind, (filename, section) in enumerate(zip(self.filenames, self.parent_section_names)):

            print('Starting to create PDF for {} section'.format(section))

            # update for new section parameters
            self.report_filename = filename
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # initialize plot characteristics
            self.plot_characteristics = dict()

            # reinitialise default configuration variables
            # modified by commandline arguments, if given
            provconf = ProvConfiguration(self, **self.commandline_arguments)

            # update self with section variables
            for k, val in self.section_opts.items():
                setattr(self, k, provconf.parse_parameter(k, val))

            # check for self defined plot characteristics file
            if self.plot_characteristics_filename == '':
                self.plot_characteristics_filename = os.path.join(CURRENT_PATH, 'conf/plot_characteristics_offline.json')
            self.plot_characteristics_templates = json.load(open(self.plot_characteristics_filename))
            self.plot_characteristics = {}

            # error when using wrong custom plot characteristics path to launch dashboard
            if 'header' not in self.plot_characteristics_templates.keys():
                msg = 'It is not possible to use the dashboard plot characteristics path to generate offline reports. Consider adding another path to plot_characteristics_filename, as in: '
                msg += 'plot_characteristics_filename = dashboard:/path/plot_characteristics_dashboard.json, offline:/path/plot_characteristics_offline.json.'
                sys.exit(msg)

            # initialise Plot class
            self.plot = Plot(read_instance=self, canvas_instance=self)

            # add general plot characteristics to self
            for k, val in self.plot_characteristics_templates['general'].items():
                setattr(self, k, val)

            # now all variables have been parsed, check validity of those, throwing errors where necessary
            provconf.check_validity()

            # set some key configuration variables
            self.relevant_temporal_resolutions = aux.get_relevant_temporal_resolutions(self.resolution)
            self.nonrelevant_temporal_resolutions = aux.get_nonrelevant_temporal_resolutions(self.resolution)
            self.data_labels = ['observations'] + list(self.experiments.keys())
            self.networkspecies = ['{}|{}'.format(network,speci) for network, speci in zip(self.network, self.species)]

            # get valid observations in date range
            aux.get_valid_obs_files_in_date_range(self, self.start_date, self.end_date)

            # update available experiments for selected fields
            aux.get_valid_experiments(self, self.start_date, self.end_date, self.resolution,
                                      self.network, self.species)

            # read data
            self.datareader.read_setup(['reset'])
            
            # initialise previous QA, flags and filter species as section values
            self.previous_qa = copy.deepcopy(self.qa)
            self.previous_flags = copy.deepcopy(self.flags)
            self.previous_filter_species = copy.deepcopy(self.filter_species)

            # if no valid data has been found be to be read, then skip to next section
            if self.invalid_read:
                print('No valid data for {} section'.format(section))
                continue
            
            # check if report type is valid
            if self.report_type not in self.report_plots.keys():
                msg = 'Error: The report type {0} cannot be found in conf/report_plots.json. '.format(self.report_type)
                msg += 'The available report types are {0}. Select one or create your own.'.format(list(self.report_plots.keys()))
                sys.exit(msg)

            # set plots that need to be made (summary and station specific)
            self.summary_plots_to_make = []
            self.station_plots_to_make = []
            if isinstance(self.report_plots[self.report_type], list):
                self.summary_plots_to_make = self.report_plots[self.report_type]
                for plot_type in self.report_plots[self.report_type]:
                    # there can be no station specific plots for map plot type
                    if plot_type[:4] != 'map-':
                        self.station_plots_to_make.append(plot_type)
            elif isinstance(self.report_plots[self.report_type], dict):
                # get summary plots
                if 'summary' in self.report_plots[self.report_type].keys():
                    if not self.report_summary:
                        print('Warning: report_summary is False, summary plots will not be created.')
                    else:
                        self.summary_plots_to_make = self.report_plots[self.report_type]['summary']
                # get station plots
                if 'station' in self.report_plots[self.report_type].keys():
                    if not self.report_stations:
                        print('Warning: report_stations is False, station plots will not be created.')
                    else:
                        for plot_type in self.report_plots[self.report_type]['station']:
                            # there can be no station specific plots for map plot type
                            if plot_type[:4] != 'map-':
                                self.station_plots_to_make.append(plot_type)

            # TODO: For Taylor diagrams, remove this piece of code until Matplotlib 3.7.2 is available
            for diagram_name in ['taylor', 'taylor-r', 'taylor-r2']:
                if (diagram_name in self.station_plots_to_make) or (diagram_name in self.summary_plots_to_make):
                    error = 'It is not possible to create Taylor diagrams yet, please remove.'
                    sys.exit(error)

            # set plot characteristics for all plot types (summary, station)
            self.plots_to_make = list(self.summary_plots_to_make)
            self.plots_to_make.extend(x for x in self.station_plots_to_make
                                      if x not in self.summary_plots_to_make)
            self.plot.set_plot_characteristics(self.plots_to_make)
            
            # define dictionary to store plot figures per page
            self.plot_dictionary = {}

            # start making PDF
            self.start_pdf()

            # remove section variables from memory 
            for k in self.section_opts:
                try:
                    vars(self).pop(k)
                except:
                    pass

    def start_pdf(self):
        """ Create PDF document where plots will be stored. """

        # get path where reports will be saved
        if '/' in self.report_filename:
            if os.path.isdir(os.path.dirname(self.report_filename)):
                reports_path = self.report_filename
        else:
            reports_path = (os.path.join(CURRENT_PATH, '../reports/')) + self.report_filename

        # create reports folder
        if not os.path.exists(os.path.dirname(reports_path)):
            if '/' in self.report_filename:
                print('Path {0} does not exist and it will be created.'.format(os.path.dirname(self.report_filename)))
            os.makedirs(os.path.dirname(reports_path))

        # add termination .pdf to filenames
        if '.pdf' not in reports_path:
            reports_path += '.pdf'

        # open new PDF file
        with PdfPages(reports_path) as pdf:
            
            self.pdf = pdf

            # initialise dictionaries to store relevant page numebrs
            if self.report_summary:
                self.summary_pages = {}
            if self.report_stations:
                self.station_pages = {}

            # get subsection names
            self.child_subsection_names = [subsection_name for subsection_name in self.subsection_names 
                                           if self.section == subsection_name.split('·')[0]]
            if len(self.child_subsection_names) > 0:
                self.subsections = self.child_subsection_names
            else:
                self.subsections = [self.section]

            # make header
            self.plot.set_plot_characteristics(['header'])
            self.plot.make_header(self.pdf, self.plot_characteristics['header'])

            # create variables to keep track of minimum and maximum data ranges across subsections
            self.data_range_min_summary = {networkspeci:np.inf for networkspeci in self.networkspecies}
            self.data_range_min_station = {networkspeci:np.inf for networkspeci in self.networkspecies}
            self.data_range_max_summary = {networkspeci:0 for networkspeci in self.networkspecies}
            self.data_range_max_station = {networkspeci:0 for networkspeci in self.networkspecies}
            self.stddev_max_summary = {networkspeci:0 for networkspeci in self.networkspecies}
            self.stddev_max_station = {networkspeci:0 for networkspeci in self.networkspecies}

            # make all plots per subsection (for all plot types except distribution plot)
            summary_plots_to_make = [plot_type for plot_type in self.summary_plots_to_make 
                                     if ('distribution' not in plot_type) and ('taylor' not in plot_type)]
            station_plots_to_make = [plot_type for plot_type in self.station_plots_to_make 
                                     if ('distribution' not in plot_type) and ('taylor' not in plot_type)]
            self.make_plots_per_subsection(summary_plots_to_make, station_plots_to_make, 
                                           do_plot_geometry_setup=True)

            # make all plots per subsection
            # for distribution plot types --> done so to calclulate data ranges across subsections first
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
                        if base_plot_type not in ['map', 'heatmap', 'table', 'taylor']: 
                            if base_plot_type == 'scatter':
                                self.plot.harmonise_xy_lims_paradigm(relevant_axs, base_plot_type, 
                                                                     self.plot_characteristics[plot_type], plot_options, 
                                                                     relim=True)
                            elif base_plot_type != 'taylor':
                                self.plot.harmonise_xy_lims_paradigm(relevant_axs, base_plot_type,
                                                                     self.plot_characteristics[plot_type], plot_options, 
                                                                     relim=True, autoscale=True)
                        
                        # update variable to reflect some formatting was performed
                        did_formatting = True

                # update variables to show if a networkspeci has been formatted
                if did_formatting:
                    formatted_networkspeci_plots = True

            # save page figures
            valid_page = False
            for page in self.plot_dictionary:
                # if page has no active data plotted, do not plot it
                n_page_plotted_labels = 0
                for ax_dict in self.plot_dictionary[page]['axs']:
                    n_page_plotted_labels += len(ax_dict['data_labels'])
                if n_page_plotted_labels > 0:
                    if not valid_page:
                        print('\nWriting PDF\n')
                        valid_page = True
                    fig = self.plot_dictionary[page]['fig']
                    self.pdf.savefig(fig, dpi=self.dpi)
                    plt.close(fig)
            if not valid_page:
                print('\n0 plots remain to write to PDF\n')

    def setup_plot_geometry(self, plotting_paradigm, networkspeci, have_setup_multispecies):
        """ Setup plotting geometry for summary or station specific plots, per network/species. """

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
                        msg = f'Warning: {plot_type} cannot be created per station '
                        msg += 'without activating the spatial colocation.'
                        print(msg)
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
                        if (base_plot_type == 'scatter') or ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                            n_plots_per_plot_type = len(self.subsections) * \
                                                    (len(self.data_labels) - 1)
                        else:
                            n_plots_per_plot_type = len(self.subsections) * \
                                                    len(self.data_labels)
                    else:
                        n_plots_per_plot_type = len(self.subsections) 
                elif plotting_paradigm == 'station':
                    if 'individual' in plot_options:
                        if (base_plot_type == 'scatter') or ('bias' in plot_options) or (z_statistic_sign == 'bias'):
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
                    ghelper = self.plot.get_taylor_diagram_ghelper(reference_stddev, plot_characteristics)
                    plot_characteristics['figure']['subplot_kw'] = {'axes_class': fa.FloatingAxes,
                                                                    'grid_helper': ghelper}
                fig, axs = plt.subplots(**plot_characteristics['figure'])

                # each page is handled as 1 figure
                # intialise page if not yet done
                if page_n not in self.plot_dictionary:
                    self.plot_dictionary[page_n] = {'fig': fig, 'plot_type': plot_type, 'axs': []}

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
                            gs = gridspec.GridSpecFromSubplotSpec(100, 100, subplot_spec=ax)
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
                            self.plot.format_axis(grid_dict, base_plot_type, plot_characteristics, set_extent=False, 
                                                  relevant_temporal_resolutions=self.relevant_temporal_resolutions)

                            # get references to periodic label annotations made, and then hide them
                            for relevant_temporal_resolution in self.relevant_temporal_resolutions:
                                annotations = [child for child in grid_dict[relevant_temporal_resolution].get_children() 
                                               if isinstance(child, matplotlib.text.Annotation)]
                                # hide annotations
                                for annotation in annotations:
                                    annotation.set_visible(False)

                        # rest of plot types
                        else:
                            
                            # collect axes in dict
                            self.plot_dictionary[page_n]['axs'].append({'handle':ax, 'data_labels':[]})
                            
                            # format axis 
                            self.plot.format_axis(ax, base_plot_type, plot_characteristics, set_extent=False)

                    # turn off axes until some data is plottted
                    ax.axis('off')
                    ax.set_visible(False)
                    
                    # iterate plot and column counts
                    plot_ii_per_type += 1
                    col_ii += 1

                # set figure attributes
                # adjust subplots?
                if 'subplots_adjust' in plot_characteristics_vars:
                    fig.subplots_adjust(**plot_characteristics['subplots_adjust'])

                # make legend?
                if 'legend' in plot_characteristics_vars:
                    plot_characteristics['legend'] = self.plot.make_legend_handles(plot_characteristics['legend'])
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
        """ Function that calls making of all plots per subsection. """

        # create variable to keep track if have setup summary plot geometry yet (done for all subsections at once)
        self.summary_plot_geometry_setup = False
        self.do_plot_geometry_setup = do_plot_geometry_setup

        # define dictionary to store stats from all subsections for heatmap and table plots
        self.stats_summary = {}
        self.stats_station = {}

        # set default markersize from density
        if self.do_plot_geometry_setup:
            self.plot.map_markersize_from_density = False

        # iterate through subsections
        for subsection_ind, subsection in enumerate(self.subsections):

            self.subsection_ind = subsection_ind
            self.subsection = subsection

            # create nested dictionary to store statistical information across all subsections
            if self.subsection not in self.stats_summary:
                self.stats_summary[self.subsection] = {}
            if self.subsection not in self.stats_station:
                self.stats_station[self.subsection] = {}

            # update the conf options for defined subsection
            if len(self.child_subsection_names) > 0:
                # get subsection variables
                self.subsection_opts = self.sub_opts[self.subsection]

                # ensure all fixed section variables defined in subsection have same value as current section variables
                self.subsection_opts = {k: (self.section_opts[k] if k in self.fixed_section_vars else val) 
                                        for (k, val) in self.subsection_opts.items()}

                # reinitialise default configuration variables
                # modified by commandline arguments, if given
                provconf = ProvConfiguration(self, **self.commandline_arguments)

                # update subsection variables
                for k, val in self.subsection_opts.items():
                    setattr(self, k, provconf.parse_parameter(k, val))

                # now all variables have been parsed, check validity of those, throwing errors where necessary
                provconf.check_validity()

            # if have no experiments, force temporal colocation to be False
            if len(self.experiments) == 0:
                self.temporal_colocation = False    

            # determine if need to re-read data (qa, flags or filter_species have changed)
            if (self.qa != self.previous_qa) or (self.flags != self.previous_flags) or\
            (self.filter_species != self.previous_filter_species):
                # re-read data
                self.datareader.read_setup(['reset'])

            # update fields available for filtering
            aux.init_representativity(self)
            aux.update_representativity_fields(self)
            aux.representativity_conf(self)
            aux.init_period(self)
            aux.update_period_fields(self)
            aux.period_conf(self)
            aux.init_metadata(self)
            aux.update_metadata_fields(self)
            aux.metadata_conf(self)

            # set previous QA, flags and filter species as subsection
            self.previous_qa = copy.deepcopy(self.qa)
            self.previous_flags = copy.deepcopy(self.flags)
            self.previous_filter_species = copy.deepcopy(self.filter_species)
            
            # filter dataset for current subsection
            print('\nFiltering data for {} subsection'.format(self.subsection))
            DataFilter(self)
            
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
            if self.report_stations:      

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

    def make_summary_plots(self, networkspeci, summary_plots_to_make):
        """ Function which makes all of summary plots for a specific subsection/networkspeci. """

        # get valid station inds for networkspeci 
        if self.temporal_colocation and len(self.data_labels) > 1:
            self.relevant_station_inds = self.valid_station_inds_temporal_colocation[networkspeci]['observations']
        else:
            self.relevant_station_inds = self.valid_station_inds[networkspeci]['observations']  

        # get N stations for networkspeci
        self.n_stations = len(self.relevant_station_inds)

        # check if we have multispecies plots
        have_multispecies = [True if 'multispecies' in plot_type else False 
                                for plot_type in summary_plots_to_make]

        # if have 0 relevant stations, continue to next networkspeci
        if self.n_stations == 0:
            print('No valid stations for {}, {}. Not making summmary plots'.format(networkspeci, self.subsection))
            # do not stop if there is any multispecies plot and we are in last subsection
            # if last subsection has data for 0 stations, it would not create them
            if (have_multispecies) and (self.subsection == self.subsections[-1]):
                pass 
            else:
                return
        else:
            print('Making {}, {} summary plots'.format(networkspeci, self.subsection)) 

        # create nested dictionary to store statistical information across all networkspecies
        if networkspeci not in self.stats_summary[self.subsection]:
            self.stats_summary[self.subsection][networkspeci] = {}

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

        # if have no valid data across data labels (no observations or experiments), then continue to next networkspeci
        if not self.selected_station_data[networkspeci]: 
            # do not stop if there is any multispecies plot and we are in last subsection
            # if last subsection has data for 0 stations, it would not create them
            if (have_multispecies) and (self.subsection == self.subsections[-1]):
                pass 
            else:
                return
        
        # setup plotting geometry for summary plots per networkspeci (for all subsections)
        if (not self.summary_plot_geometry_setup) & (self.do_plot_geometry_setup):
            self.setup_plot_geometry('summary', networkspeci, self.made_networkspeci_summary_plots)

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

            # set variable to know if we need to create plot in last subsection
            plot_type_df = self.get_plot_type_df(base_plot_type)

            # update nested dictionary to store statistical information
            if plot_type_df:
                self.update_stats_tables('summary', base_plot_type, plot_type, zstat, networkspeci, plot_options)
                # do not make plot until last subsection (multispecies or not) if we need a dataframe
                if (self.subsection != self.subsections[-1]):
                    continue

            # do not make multispecies plots
            if ('multispecies' in plot_options):
                # unless we are in first instance (boxplot)
                if (not plot_type_df) and (self.made_networkspeci_summary_plots):
                    continue
                # until last subsection (table, heatmap, statsummary)
                elif (plot_type_df) and ((self.subsection != self.subsections[-1]) 
                                            or (networkspeci != self.networkspecies[-1])):
                    continue

            # make plot
            print('Making summary {0}'.format(plot_type))
            plot_indices = self.make_plot('summary', plot_type, plot_options, networkspeci)

            # do formatting
            relevant_axs, relevant_data_labels = self.get_relevant_axs_per_networkspeci_plot_type_page_ind(base_plot_type, 
                                                                                                            plot_indices)
            self.plot.do_formatting(relevant_axs, relevant_data_labels, networkspeci,
                                    base_plot_type, plot_type, plot_options, 'summary')

        # update N total pages 
        self.n_total_pages = len(self.plot_dictionary)

        # update variable when summary plots have been made for a networkspecies
        self.made_networkspeci_summary_plots = True

    def make_station_plots(self, networkspeci, station_plots_to_make):
        """ Function which makes all of station plots for a specific subsection/networkspeci. """

        # get valid station inds for networkspeci 
        if self.temporal_colocation and len(self.data_labels) > 1:
            self.relevant_station_inds = self.valid_station_inds_temporal_colocation[networkspeci]['observations']
        else:
            self.relevant_station_inds = self.valid_station_inds[networkspeci]['observations']  

        # get N stations for networkspeci
        self.n_stations = len(self.relevant_station_inds)

        # check if we have multispecies plots
        have_multispecies = [True if 'multispecies' in plot_type else False 
                                for plot_type in station_plots_to_make][0]

        # if have 0 relevant stations, continue to next networkspeci
        if self.n_stations == 0:
            print('No valid stations for {}, {}. Not making station plots'.format(networkspeci, self.subsection))
            # do not stop if there is any multispecies plot and we are in last subsection
            # if last subsection has data for 0 stations, it would not create them
            if (have_multispecies) and (networkspeci == self.networkspecies[-1]):
                pass 
            else:
                return
        else:
            print('Making {}, {} station plots'.format(networkspeci, self.subsection)) 
        
        # create nested dictionary to store statistical information across all networkspecies
        if networkspeci not in self.stats_station:
            self.stats_station[self.subsection][networkspeci] = {}

        # setup plotting geometry for station plots per networkspeci (for one subsection)
        if self.do_plot_geometry_setup:
            self.setup_plot_geometry('station', networkspeci, self.made_networkspeci_station_plots)
            
        print('setup geometry')

        # initialise station ind as -1
        self.station_ind = -1

        for i, relevant_station_ind in enumerate(self.relevant_station_inds):
            
            # gather some information about current station
            self.station_ind += 1
            self.current_lon = round(self.station_longitudes[networkspeci][relevant_station_ind], 2)
            self.current_lat = round(self.station_latitudes[networkspeci][relevant_station_ind], 2)
            self.current_station_name = self.station_names[networkspeci][relevant_station_ind]
            self.current_station_reference = self.station_references[networkspeci][relevant_station_ind]

            # get selected station data 
            get_selected_station_data(read_instance=self, canvas_instance=self, 
                                      networkspecies=[networkspeci], 
                                      station_index=relevant_station_ind, 
                                      data_range_min=self.data_range_min_station, 
                                      data_range_max=self.data_range_max_station,
                                      stddev_max=self.stddev_max_station)

            # if have no valid data across data labels (no observations or experiments), then continue to next station
            if not self.selected_station_data[networkspeci]:
                # do not stop if there is any multispecies plot and we are in last subsection
                # if last subsection has data for 0 stations, it would not create them
                if (have_multispecies) and (networkspeci == self.networkspecies[-1]):
                    pass 
                else:
                    continue
            
            # update data range min/maxes for station paradigm
            if self.selected_station_data_min[networkspeci] < self.data_range_min_station[networkspeci]:
                self.data_range_min_station[networkspeci] = copy.deepcopy(self.selected_station_data_min[networkspeci]) 
            if self.selected_station_data_max[networkspeci] > self.data_range_max_station[networkspeci]:
                self.data_range_max_station[networkspeci] = copy.deepcopy(self.selected_station_data_max[networkspeci])
            if self.selected_station_stddev_max[networkspeci] > self.stddev_max_station[networkspeci]:
                self.stddev_max_station[networkspeci] = copy.deepcopy(self.selected_station_stddev_max[networkspeci])

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

                # set variable to know if we need to create plot in last subsection
                plot_type_df = self.get_plot_type_df(base_plot_type)
                
                # remove station measurement method code from references for GHOST data
                # if this has not already been removed
                if (('multispecies' in plot_options) and (self.reading_ghost) 
                    and ('_' in self.current_station_reference)):
                    self.current_station_reference = '_'.join(self.current_station_reference.split("_")[:-1])

                # collect statistics in dictionaries
                if plot_type_df:

                    # do not update dictionary for multispecies plots if spatial colocation is turned off
                    if ('multispecies' in plot_options) and (not self.spatial_colocation):
                        continue
                    
                    # update nested dictionary to store statistical information
                    self.update_stats_tables('station', base_plot_type, plot_type, zstat, networkspeci, plot_options)
                    
                    # only plot last subsection and last networkspecies for multispecies
                    if ('multispecies' in plot_options) and (networkspeci != self.networkspecies[-1]):
                        continue
                    
                # multispecies (non dataframe multispecies plots)
                if ('multispecies' in plot_options) and (not plot_type_df):
                    # do plot if we are in first instance and spatial colocation is active
                    if (not self.made_networkspeci_station_plots) & (self.spatial_colocation):
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
                    else:
                        continue

                # make plot
                print('Making station {2} for {3} ({0}/{1})'.format(i+1, 
                                                                    len(self.relevant_station_inds),
                                                                    plot_type, 
                                                                    self.current_station_name))   
                
                plot_indices = self.make_plot('station', plot_type, plot_options, networkspeci)

                # do not format Taylor diagrams until last station
                if (base_plot_type == 'taylor'):
                    if relevant_station_ind != self.relevant_station_inds[-1]:
                        continue

                # do formatting
                relevant_axs, relevant_data_labels = self.get_relevant_axs_per_networkspeci_plot_type_page_ind(base_plot_type, 
                                                                                                               plot_indices)
                self.plot.do_formatting(relevant_axs, relevant_data_labels, networkspeci,
                                        base_plot_type, plot_type, plot_options, 'station')

        # update N total pages 
        self.n_total_pages = len(self.plot_dictionary)

        # update variable now station plots have been made for a networkspecies
        self.made_networkspeci_station_plots = True
    
    def get_plot_type_df(self, base_plot_type):   
        """ 
        Function that determines if plot type is one that involves collation of statistics across subsections
        and networkspecies.
        """

        if base_plot_type in ['table', 'heatmap', 'statsummary']:
            plot_type_df = True
        else: 
            plot_type_df = False

        return plot_type_df

    def update_stats_tables(self, plotting_paradigm, base_plot_type, plot_type, zstat, networkspeci, plot_options):
        """ Update statistical tables for plot types that collate data across various networkspecies/subsections. """ 

        # set statistics to calculate based on plot type
        if base_plot_type in ['heatmap', 'table']:
            stats = [zstat]
        elif base_plot_type == 'statsummary':
            if 'bias' in plot_options:
                stats = self.plot_characteristics[plot_type]['experiment_bias']
            else:
                stats = self.plot_characteristics[plot_type]['basic']

        # iterate through all data labels
        for data_label in self.data_labels:

            # create nested dictionary to store statistical information across all subsections
            for stat in stats:

                # get zstat information 
                zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=stat)

                # skip observations data label when plotting bias
                if (data_label == 'observations') & (z_statistic_sign == 'bias'):
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
                    # if relevant stat is expbias stat, then ensure temporal colocation is active
                    if (base_plot_type == 'statsummary') and (stat in expbias_stats) and (not self.temporal_colocation):
                        data_to_add = np.NaN
                    # otherwise calculate statistic
                    else:
                        if z_statistic_sign == 'bias':
                            data_to_add = calculate_statistic(self, self, networkspeci, zstat, ['observations'], [data_label])
                        else:
                            data_to_add = calculate_statistic(self, self, networkspeci, zstat, [data_label], [])
                else:
                    data_to_add = np.NaN

                # add data to dicts
                if plotting_paradigm == 'summary':
                    self.stats_summary[self.subsection][networkspeci][stat][data_label] = data_to_add
                elif plotting_paradigm == 'station':
                    self.stats_station[self.subsection][networkspeci][self.current_station_reference][stat][data_label] = data_to_add
            
        return None

    def make_plot(self, plotting_paradigm, plot_type, plot_options, networkspeci):
        """ Function that calls making of any type of plot. """

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

        # get all data_labels for selected_station_data
        # if multispecies plot then ensure have data_labels from across all networkspecies
        if 'multispecies' in plot_options: 
            data_labels = []
            for ns in self.selected_station_data:
                data_labels.extend(self.selected_station_data_labels[ns])
            _, idx = np.unique(data_labels, return_index=True)
            data_labels = np.array(data_labels)[np.sort(idx)].tolist()
        else:
            data_labels = copy.deepcopy(self.data_labels)

        # if have no valid data labels then return 
        if len(data_labels) == 0:
            return plot_indices

        # if are making bias plot, and have no valid experiment data then cannot make plot type
        if (('bias' in plot_options) or (z_statistic_sign == 'bias')) & (len(data_labels) < 2):
            return plot_indices

        # get data labels without observations
        data_labels_sans_obs = copy.deepcopy(data_labels)
        data_labels_sans_obs.remove('observations')

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
                z1 = ['observations'] * len(data_labels_sans_obs)
                z2 = data_labels_sans_obs
            else:
                z1 = data_labels
                z2 = ['']*len(data_labels)

            #iterate through relevant data labels making plots
            for z1_label, z2_label in zip(z1, z2):
                
                # get relevant page/axis to plot on
                axis_ind = (current_plot_ind * len(self.subsections)) + self.subsection_ind
                relevant_page, page_ind, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, 
                                                                                     plot_type, axis_ind)

                # set axis title
                if relevant_axis.get_title() == '':

                    if z1_label == 'observations':
                        label = copy.deepcopy(z1_label) 
                    else:
                        label = self.experiments[z1_label]
                    axis_title_label = '{}\n{} '.format(label, self.subsection)
                    if self.n_stations == 1:
                        axis_title_label += '(1 station)'
                    else:
                        axis_title_label += '({} stations)'.format(self.n_stations)
                    self.plot.set_axis_title(relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

                # set map extent ? 
                if self.map_extent:
                    self.plot.set_map_extent(relevant_axis)

                # calculate z statistic
                z_statistic, self.active_map_valid_station_inds = calculate_statistic(self, self, networkspeci,
                                                                                      zstat, z1_label, z2_label, 
                                                                                      map=True)

                # make map plot
                self.plot.make_map(relevant_axis, networkspeci, z_statistic, self.plot_characteristics[plot_type], 
                                   plot_options=plot_options)
                
                # save plot information for later formatting 
                if z2 == '':
                    self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(z1)
                else:
                    self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(z2)
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
            
            #if making indivdual plots, iterate through data labels one at a time, otherwise pass all data labels 
            #together
            if 'individual' in plot_options:
                iter_data_labels = copy.deepcopy(data_labels)
            else:
                iter_data_labels = [[data_labels]]
            for data_labels in iter_data_labels:

                if type(data_labels) != list:
                    data_labels = [data_labels]

                # get relevant axis to plot on
                if plotting_paradigm == 'summary':
                    if 'individual' in plot_options:
                        if ((base_plot_type == 'scatter') or ('bias' in plot_options) or 
                            (z_statistic_sign == 'bias')):
                                axis_ind = (current_plot_ind + self.subsection_ind + (len(self.experiments) - 1) * self.subsection_ind)
                        else:
                            axis_ind = (current_plot_ind + self.subsection_ind + len(self.experiments) * self.subsection_ind)
                    else:
                        axis_ind = self.subsection_ind
                    station_inds = copy.deepcopy(self.relevant_station_inds) 
                elif plotting_paradigm == 'station':
                    if 'individual' in plot_options:
                        if ((base_plot_type == 'scatter') or ('bias' in plot_options) or 
                            (z_statistic_sign == 'bias')):
                            axis_ind = (current_plot_ind + self.station_ind + (len(self.experiments) - 1) * self.station_ind)
                        else:
                            axis_ind = (current_plot_ind + self.station_ind + len(self.experiments) * self.station_ind)
                    else:
                        axis_ind = self.station_ind
                    station_inds = [self.relevant_station_inds[self.station_ind]]
                
                # get relevant axis
                relevant_page, page_ind, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, 
                                                                                     plot_type, axis_ind)

                # set axis title (only if not previously set)
                if isinstance(relevant_axis, dict):
                    for relevant_temporal_resolution, sub_ax in relevant_axis.items():
                        if relevant_temporal_resolution == 'hour':
                            axis_title = sub_ax.get_title()
                            break
                else:
                    axis_title = relevant_axis.get_title()

                # axis title is empty?
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
                                                                                  self.plot_characteristics[plot_type]['round_decimal_places'],
                                                                                  self.current_lat,
                                                                                  self.plot_characteristics[plot_type]['round_decimal_places'])
                            
                    # set title
                    self.plot.set_axis_title(relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

                # set xlabel and ylabel (only if not previously set)
                if isinstance(relevant_axis, dict):
                    for relevant_temporal_resolution, sub_ax in relevant_axis.items():
                        if relevant_temporal_resolution in ['hour','month']:
                            axis_xlabel = 'NaN'
                            axis_ylabel = sub_ax.get_ylabel()
                            break
                else:
                    axis_xlabel = relevant_axis.get_xlabel()
                    axis_ylabel = relevant_axis.get_ylabel()

                # axis xlabel is empty?
                if (axis_xlabel == '') or (axis_xlabel == 'measurement_units'):
                    if 'xlabel' in self.plot_characteristics[plot_type]:
                        if self.plot_characteristics[plot_type]['xlabel']['xlabel'] == 'measurement_units':
                            xlabel = '[{}]'.format(self.measurement_units[networkspeci.split('|')[-1]])
                        else:
                            xlabel = self.plot_characteristics[plot_type]['xlabel']['xlabel']
                    else:
                        xlabel = ''
                    # set xlabel
                    if xlabel != '':
                        self.plot.set_axis_label(relevant_axis, 'x', xlabel, self.plot_characteristics[plot_type])

                # axis ylabel is empty?
                if (axis_ylabel == '') or (axis_ylabel == 'measurement_units'):
                    if base_plot_type in ['periodic']:
                        if z_statistic_type == 'basic':
                            ylabel = basic_stats[base_zstat]['label']
                            ylabel_units = basic_stats[base_zstat]['units']
                        else:
                            ylabel = expbias_stats[base_zstat]['label']
                            ylabel_units = expbias_stats[base_zstat]['units']
                        if ylabel_units == 'measurement_units':
                            ylabel_units = self.measurement_units[networkspeci.split('|')[-1]] 
                        if ylabel_units != '':
                            ylabel = '[{}]'.format(ylabel_units)
                    else:
                        if 'ylabel' in self.plot_characteristics[plot_type]:
                            if self.plot_characteristics[plot_type]['ylabel']['ylabel'] == 'measurement_units':
                                ylabel = '[{}]'.format(self.measurement_units[networkspeci.split('|')[-1]])
                            else:
                                ylabel = self.plot_characteristics[plot_type]['ylabel']['ylabel']
                        else:
                            ylabel = ''
                    # set ylabel
                    if ylabel != '':
                        self.plot.set_axis_label(relevant_axis, 'y', ylabel, self.plot_characteristics[plot_type])

                # get plotting function                
                func = getattr(self.plot, 'make_{}'.format(base_plot_type.split('-')[0]))

                # determine if have some data to plot
                plot_validity = False
                if 'multispecies' in plot_options:
                    for ns in self.selected_station_data:
                        for data_label in data_labels:
                            if data_label in self.selected_station_data_labels[ns]:
                                plot_validity = True
                else:
                    for data_label in data_labels:
                        if data_label in self.selected_station_data_labels[networkspeci]:
                            plot_validity = True

                if plot_validity:
                    if base_plot_type == 'metadata':
                        func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                                plot_options=plot_options, station_inds=station_inds)
                    elif base_plot_type == 'periodic':
                        func(relevant_axis, networkspeci, data_labels, 
                                self.plot_characteristics[plot_type], zstat=zstat, plot_options=plot_options)
                    elif base_plot_type == 'distribution':
                        func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                            plot_options=plot_options, 
                            data_range_min=data_range_min, 
                            data_range_max=data_range_max) 
                    elif base_plot_type == 'taylor':
                        func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                            plot_options=plot_options, 
                            stddev_max=stddev_max)
                    else:
                        func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                            plot_options=plot_options) 
                    
                    # save plot information for later formatting
                    self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].extend(data_labels)
                    plot_index = [relevant_page, page_ind]
                    if plot_index not in plot_indices:
                        plot_indices.append(plot_index)

                    # turn axis/axes on
                    if base_plot_type in ['periodic','periodic-violin']:
                        for relevant_temporal_resolution in self.relevant_temporal_resolutions:
                            relevant_axis[relevant_temporal_resolution].axis('on')
                            relevant_axis[relevant_temporal_resolution].set_visible(True)

                            #get references to periodic label annotations made, and then show them
                            annotations = [child for child in relevant_axis[relevant_temporal_resolution].get_children() if isinstance(child, matplotlib.text.Annotation)]
                            # hide annotations
                            for annotation in annotations:
                                annotation.set_visible(True)
                    else:
                        relevant_axis.axis('on')
                        relevant_axis.set_visible(True)

                    # iterate number of plots made for current type of plot 
                    current_plot_ind += 1     

        # make plot heatmap / table / statsummary plot
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

            # get data labels (based on statistic type)
            if z_statistic_sign == 'bias':
                data_labels = list(self.experiments.keys())
            else:
                data_labels = ['observations'] + list(self.experiments.keys())

            if base_plot_type in ['heatmap', 'table']:
                
                # create empty dataframe with networkspecies and subsections
                index = pd.MultiIndex.from_product([networkspecies, self.subsections],
                                                    names=["networkspecies", "subsections"])
                stats_df = pd.DataFrame(np.nan, index=index, columns=data_labels, dtype=np.float)
                
                # convert stats_summary and stats_station dicts to dataframes
                for subsection in subsections:
                    for networkspeci in networkspecies:
                        stats_per_data_label = []
                        for data_label in data_labels:
                            # initialise stat with nan
                            stat_to_append = np.NaN
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
                        stats_df.loc[(networkspeci, subsection)] = stats_per_data_label
                
            elif base_plot_type == 'statsummary':

                # get stats
                if 'bias' in plot_options:
                    stats = self.plot_characteristics[plot_type]['experiment_bias']
                else:
                    stats = self.plot_characteristics[plot_type]['basic']

                # create empty dataframe with networkspecies and subsections
                index = pd.MultiIndex.from_product([self.networkspecies, self.subsections, data_labels],
                                                    names=["networkspecies", "subsections", "labels"])
                stats_df = pd.DataFrame(np.nan, index=index, columns=stats, dtype=np.float)

                # convert stats_summary and stats_station dicts to dataframes
                for subsection in subsections:
                    for networkspeci in networkspecies:
                        for data_label in data_labels:
                            if ('bias' in plot_options) and (data_label == 'observations'):
                                continue
                            stats_per_data_label = []
                            for stat in stats:
                                # initialise stat as nan
                                stat_to_append = np.NaN
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
                            stats_df.loc[(networkspeci, subsection, data_label)] = stats_per_data_label

            # turn on relevant axis if dataframe has values or not all NaN
            if (len(stats_df.index) > 0) & (not stats_df.isnull().values.all()):
                
                if base_plot_type == 'statsummary':
                    # make statsummary
                    func = getattr(self.plot, 'make_table')
                    func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         statsummary=True, plot_options=plot_options, subsection=self.subsection, 
                         plotting_paradigm=plotting_paradigm, stats_df=stats_df)
                else:
                    # make table/heatmap
                    func = getattr(self.plot, 'make_{}'.format(base_plot_type))
                    func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         plot_options=plot_options, subsection=self.subsection, 
                         plotting_paradigm=plotting_paradigm, stats_df=stats_df)
                
                # save plot information for later formatting
                self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].extend(data_labels)
                plot_index = [relevant_page, page_ind]
                if plot_index not in plot_indices:
                    plot_indices.append(plot_index)

                # turn axis on
                relevant_axis.axis('on')
                relevant_axis.set_visible(True)

            # set axis title
            if relevant_axis.get_title() == '':
                if plotting_paradigm == 'summary':
                    if (base_plot_type == 'statsummary') & ('multispecies' not in plot_options):
                        if self.n_stations == 1:
                            axis_title_label = '{} ({} station)'.format(self.subsection, self.n_stations)
                        else:
                            axis_title_label = '{} ({} stations)'.format(self.subsection, self.n_stations)
                    else:
                        axis_title_label = ''
                elif plotting_paradigm == 'station':
                    axis_title_label = '{}, {} ({:.{}f}, {:.{}f})'.format(self.current_station_reference,
                                                                          self.current_station_name, 
                                                                          self.current_lon,
                                                                          self.plot_characteristics[plot_type]['round_decimal_places'],
                                                                          self.current_lat,
                                                                          self.plot_characteristics[plot_type]['round_decimal_places'])
                self.plot.set_axis_title(relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

        return plot_indices

    def get_relevant_page_axis(self, plotting_paradigm, networkspeci, plot_type, axis_ind):
        """ Get relevant page and axis for current plot type/subsection/axis index. """

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
        """Get relevant axs per plot type"""
        
        relevant_axs = []
        relevant_data_labels = []
        for relevant_page in relevant_pages:
            if base_plot_type in ['periodic', 'periodic-violin']:
                for ax in self.plot_dictionary[relevant_page]['axs']:
                    for relevant_temporal_resolution in self.relevant_temporal_resolutions:
                        relevant_axs.append(ax['handle'][relevant_temporal_resolution])
                        relevant_data_labels.append(ax['data_labels'])
            else:
                for ax in self.plot_dictionary[relevant_page]['axs']:
                    relevant_axs.append(ax['handle'])
                    relevant_data_labels.append(ax['data_labels'])

        return relevant_axs, relevant_data_labels
    
    def get_relevant_axs_per_networkspeci_plot_type_page_ind(self, base_plot_type, plot_indices):
        """Get relevant axs per plot type and page ind"""

        relevant_axs = []
        relevant_data_labels = []
        for relevant_page, page_ind in plot_indices:
            if base_plot_type in ['periodic', 'periodic-violin']:
                for relevant_temporal_resolution in self.relevant_temporal_resolutions:
                    relevant_axs.append(self.plot_dictionary[relevant_page]['axs'][page_ind]['handle'][relevant_temporal_resolution])
                    relevant_data_labels.append(self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'])
            else:
                relevant_axs.append(self.plot_dictionary[relevant_page]['axs'][page_ind]['handle'])
                relevant_data_labels.append(self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'])

        return relevant_axs, relevant_data_labels

def main(**kwargs):
    """ Main function when running offine reports. """
   
    ProvidentiaOffline(**kwargs)
