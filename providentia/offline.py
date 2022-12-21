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

from .read import DataReader
from .filter import DataFilter
from .plot import Plot
from .statistics import to_pandas_dataframe
from .statistics import calculate_z_statistic
from .statistics import generate_colourbar
from .statistics import get_z_statistic_info
from .configuration import ProvConfiguration
from providentia import aux

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))

class ProvidentiaOffline:
    """ Class to create Providentia offline reports. """

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
        elif ('config' in kwargs) and (not os.path.exists(kwargs['config'])):     
            error = 'Error: The path to the configuration file specified in the command line does not exist.'
            sys.exit(error)
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            sys.exit(error)

        # load report plot presets
        self.report_plots = json.load(open(os.path.join(CURRENT_PATH, 'conf/report_plots.json')))

        # check for self defined plot characteristics file
        if self.plot_characteristics_filename == '':
            self.plot_characteristics_filename = os.path.join(CURRENT_PATH, 'conf/plot_characteristics_offline.json')
        self.plot_characteristics_templates = json.load(open(self.plot_characteristics_filename))
        self.plot_characteristics = {}

        # add general plot characteristics to self
        for k, val in self.plot_characteristics_templates['general'].items():
            setattr(self, k, val)

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

        # initialise Plot class
        self.plot = Plot(read_instance=self, canvas_instance=self)

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

            # now all variables have been parsed, check validity of those, throwing errors where necessary
            provconf.check_validity()

            # set some key configuration variables
            self.relevant_temporal_resolutions = aux.get_relevant_temporal_resolutions(self.resolution)
            self.data_labels = ['observations'] + list(self.experiments.keys())
            self.networkspecies = ['{}|{}'.format(network,speci) for network, speci in zip(self.network, self.species)]

            # get valid observations in date range
            aux.get_valid_obs_files_in_date_range(self, self.start_date, self.end_date)

            # update available experiments for selected fields
            aux.get_valid_experiments(self, self.start_date, self.end_date, self.resolution,
                                      self.network, self.species)

            # read data
            self.datareader.read_setup(['reset'])
            #initialise previous QA, flags and filter species as section values
            self.previous_qa = copy.deepcopy(self.qa)
            self.previous_flags = copy.deepcopy(self.flags)
            self.previous_filter_species = copy.deepcopy(self.filter_species)

            # if no valid data has been found be to be read, then skip to next section
            if self.invalid_read:
                print('No valid data for {} section'.format(section))
                continue

            # set plot characteristics
            try:
                plot_types = self.report_plots[self.report_type]
            except KeyError:
                msg = 'Error: The report type {0} cannot be found in conf/report_plots.json. '.format(self.report_type)
                msg += 'The available report types are {0}. Select one or create your own.'.format(list(self.report_plots.keys()))
                sys.exit(msg)
            self.plot.set_plot_characteristics(plot_types)

            # define dictionary to store plot figures per page
            self.plot_dictionary = {}

            # set plots that need to be made (summary and station specific)
            self.summary_plots_to_make = list(self.plot_characteristics.keys())
            self.station_plots_to_make = []
            for plot_type in self.summary_plots_to_make:
                # there can be no station specific plots for map plot type
                if plot_type[:4] != 'map-':
                    self.station_plots_to_make.append(plot_type)

            # start making PDF
            self.start_pdf()

            # remove section variables from memory (if not last section)
            if section_ind != (len(self.parent_section_names) - 1):
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
                                           if self.section == subsection_name.split('|')[0]]
            if len(self.child_subsection_names) > 0:
                self.subsections = self.child_subsection_names
            else:
                self.subsections = [self.section]

            # make header
            self.plot.set_plot_characteristics(['header'])
            self.plot.make_header(self.pdf, self.plot_characteristics['header'])

            # define dictionary to store stats from all subsections for heatmap and table plots
            self.subsection_stats_summary = {}
            self.subsection_stats_station = {}

            # create variables to keep track of minimum and maximum data ranges across subsections
            self.data_range_min_summary = {networkspeci:np.inf for networkspeci in self.networkspecies}
            self.data_range_min_station = {networkspeci:np.inf for networkspeci in self.networkspecies}
            self.data_range_max_summary = {networkspeci:0 for networkspeci in self.networkspecies}
            self.data_range_max_station = {networkspeci:0 for networkspeci in self.networkspecies}

            # make all plots per subsection (for all plot types except distribution plot)
            summary_plots_to_make_nondist = [plot_type for plot_type in self.summary_plots_to_make if 'distribution' not in plot_type]
            station_plots_to_make_nondist = [plot_type for plot_type in self.station_plots_to_make if 'distribution' not in plot_type]
            self.make_plots_per_subsection(summary_plots_to_make_nondist, station_plots_to_make_nondist, do_plot_geometry_setup=True)

            # make all plots per subsection (for only distribution plot types --> done so to calclulate data ranges across subsections first)
            summary_plots_to_make_dist = [plot_type for plot_type in self.summary_plots_to_make if 'distribution' in plot_type]
            station_plots_to_make_dist = [plot_type for plot_type in self.station_plots_to_make if 'distribution' in plot_type]
            if (len(summary_plots_to_make_dist) > 0) or (len(station_plots_to_make_dist) > 0):
                self.make_plots_per_subsection(summary_plots_to_make_dist, station_plots_to_make_dist)
            
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

                    # if a multispecies plot is wanted then this is only made on first instance of formatting a networkspeci 
                    if ('multispecies' in plot_options) & (formatted_networkspeci_plots):
                        continue

                    # get zstat information from plot_type
                    zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)

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

                        # if have no relevant axs, continue to next paradigm
                        if len(relevant_axs) == 0:
                            continue

                        # get data ranges for plotting paradigm
                        if plotting_paradigm == 'summary':
                            data_range_min = self.data_range_min_summary[networkspeci]
                            data_range_max = self.data_range_max_summary[networkspeci]
                        elif plotting_paradigm == 'station':
                            data_range_min = self.data_range_min_station[networkspeci]
                            data_range_max = self.data_range_max_station[networkspeci]

                        # generate colourbars for required plots in paradigm on each relevant page
                        if 'cb' in list(self.plot_characteristics[plot_type].keys()):
                            # get all cb_axs for plot_type across relevant pages
                            cb_axs = [self.plot_dictionary[relevant_page]['cb_ax'] for relevant_page in relevant_pages]
                            generate_colourbar(self, relevant_axs, cb_axs, zstat, self.plot_characteristics[plot_type], 
                                               networkspeci.split('|')[-1])

                        # harmonise xy limits for plot paradigm
                        if base_plot_type not in ['map', 'heatmap', 'table']: 
                            if base_plot_type == 'periodic-violin':
                                self.plot.harmonise_xy_lims_paradigm(relevant_axs, base_plot_type, 
                                                                     self.plot_characteristics[plot_type], plot_options, 
                                                                     ylim=[data_range_min, data_range_max],
                                                                     relim=True, autoscale_x=True)
                            elif base_plot_type == 'scatter':
                                self.plot.harmonise_xy_lims_paradigm(relevant_axs, base_plot_type, 
                                                                     self.plot_characteristics[plot_type], plot_options, 
                                                                     relim=True)
                            else:
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
            
            # if a multispecies plot is wanted then only make this on first instance of plotting a networkspecies
            # if making a multispecies plot per specific station, spatial colocation must be also active
            if 'multispecies' in plot_options:
                if plotting_paradigm == 'summary':
                    if have_setup_multispecies:
                        continue
                elif plotting_paradigm == 'station':
                    if (have_setup_multispecies) or (not self.spatial_colocation):
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
            zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)

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
                                annotations = [child for child in grid_dict[relevant_temporal_resolution].get_children() if isinstance(child, matplotlib.text.Annotation)]
                                # hide annotations
                                for annotation in annotations:
                                    annotation.set_visible(False)

                        # rest of plot types
                        else:
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
        summary_plot_geometry_setup = False

        # set default markersize from density
        if do_plot_geometry_setup:
            self.plot.map_markersize_from_density = False

        # iterate through subsections
        for subsection_ind, subsection in enumerate(self.subsections):

            self.subsection_ind = subsection_ind
            self.subsection = subsection

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

            #set previous QA, flags and filter species as subsection
            self.previous_qa = copy.deepcopy(self.qa)
            self.previous_flags = copy.deepcopy(self.flags)
            self.previous_filter_species = copy.deepcopy(self.filter_species)
            
            # filter dataset for current subsection
            print('\nFiltering data for {} subsection'.format(self.subsection))
            DataFilter(self)
            
            # iterate through networks and species, creating plots
            self.n_total_pages = 0

            # make summary plots?
            if self.report_summary:

                # set variable to inform when have made 1 set of networkspecies plots for summary
                made_networkspeci_summary_plots = False

                # determine if all plots are multispecies plots
                all_multispecies_plots  = np.all([True if 'multispecies' in plot_type.split('_')[1:] else False for plot_type in summary_plots_to_make])

                # iterate through networkspecies
                for networkspeci in self.networkspecies:

                    # if have 0 plot types that are non-multispecies, and have previously made networkspeci plots,
                    # then continue 
                    if (made_networkspeci_summary_plots) & (all_multispecies_plots): 
                        continue

                    # get valid station inds for networkspeci 
                    if self.temporal_colocation and len(self.data_labels) > 1:
                        self.relevant_station_inds = self.valid_station_inds_temporal_colocation[networkspeci]['observations']
                    else:
                        self.relevant_station_inds = self.valid_station_inds[networkspeci]['observations']  

                    # get N stations for networkspeci
                    self.n_stations = len(self.relevant_station_inds)

                    #if have 0 relevant stations, continue to next networkspeci
                    if self.n_stations == 0:
                        print('No valid stations for {}, {}. Not making summmary plots'.format(networkspeci, self.subsection))
                        continue
                    else:
                        print('Making {}, {} summary plots'.format(networkspeci, self.subsection)) 

                    if not made_networkspeci_summary_plots:
                        # get median timeseries across data from filtered data, and place it pandas dataframe
                        to_pandas_dataframe(read_instance=self, canvas_instance=self, 
                                            networkspecies=self.networkspecies)

                        # update data range min/maxes for summary paradigm 
                        for ns in self.networkspecies:
                            if self.selected_station_data_min[ns] < self.data_range_min_summary[ns]:
                                self.data_range_min_summary[ns] = copy.deepcopy(self.selected_station_data_min[ns]) 
                            if self.selected_station_data_max[ns] > self.data_range_max_summary[ns]:
                                self.data_range_max_summary[ns] = copy.deepcopy(self.selected_station_data_max[ns])
                    
                    # if have no valid data across data labels (no observations or experiments), then continue 
                    if not self.selected_station_data[networkspeci]:
                        continue

                    # setup plotting geometry for summary plots per networkspeci (for all subsections)
                    if (not summary_plot_geometry_setup) & (do_plot_geometry_setup):
                        self.setup_plot_geometry('summary', networkspeci, made_networkspeci_summary_plots)
                    
                    # iterate through plots to make
                    for plot_type in summary_plots_to_make:

                        # get zstat information from plot_type
                        zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)
                        
                        # get base plot type (without stat and options)
                        if zstat:
                            base_plot_type = plot_type.split('-')[0] 
                        else:
                            base_plot_type = plot_type.split('_')[0] 

                        #get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                        plot_options = plot_type.split('_')[1:]

                        # if a multispecies plot is wanted then only make this on first instance of plotting a networkspecies
                        if ('multispecies' in plot_options) & (made_networkspeci_summary_plots):
                            continue

                        # make plot
                        print('Making summary {0}'.format(plot_type))
                        plot_indices = self.make_plot('summary', plot_type, plot_options, networkspeci)
                        
                        # get relevant axs and plot types per networkspeci / plot type
                        relevant_axs, relevant_data_labels = self.get_relevant_axs_per_networkspeci_plot_type(base_plot_type, 
                                                                                                              [plot_indices[0][0]])

                        # do formatting
                        self.plot.do_formatting(relevant_axs, relevant_data_labels, networkspeci,
                                                base_plot_type, plot_type, plot_options, 'summary')

                    # update N total pages 
                    self.n_total_pages = len(self.plot_dictionary)

                    # update variable now summary plots have been made for a networkspecies
                    made_networkspeci_summary_plots = True

            # make station specific plots?
            if self.report_stations:      

               # set variable to inform when have made 1 set of networkspecies plots for stations
               made_networkspeci_station_plots = False        

               # determine if all plots are multispecies plots
               all_multispecies_plots = np.all([True if 'multispecies' in plot_type.split('_')[1:] else False for plot_type in station_plots_to_make])

               # iterate through networkspecies
               for networkspeci in self.networkspecies:

                    # if have 0 plot types that are non-multispecies, and have previously made networkspeci plots,
                    # then continue 
                    if (made_networkspeci_station_plots) & (all_multispecies_plots): 
                        continue

                    # get valid station inds for networkspeci 
                    if self.temporal_colocation and len(self.data_labels) > 1:
                        self.relevant_station_inds = self.valid_station_inds_temporal_colocation[networkspeci]['observations']
                    else:
                        self.relevant_station_inds = self.valid_station_inds[networkspeci]['observations']  

                    # get N stations for networkspeci
                    self.n_stations = len(self.relevant_station_inds)

                    #if have 0 relevant stations, continue to next networkspeci
                    if self.n_stations == 0:
                        print('No valid stations for {}, {}. Not making station plots'.format(networkspeci, self.subsection))
                        continue
                    else:
                        print('Making {}, {} station plots'.format(networkspeci, self.subsection)) 

                    # setup plotting geometry for station plots per networkspeci (for one subsection)
                    if do_plot_geometry_setup:
                        self.setup_plot_geometry('station', networkspeci, made_networkspeci_station_plots)
                        
                    # initialise station ind as -1
                    self.station_ind = -1

                    for i, relevant_station_ind in enumerate(self.relevant_station_inds):

                        # gather some information about current station
                        self.station_ind += 1
                        valid_station_references = self.metadata_in_memory[networkspeci]['station_reference'][relevant_station_ind, :]
                        self.current_station_reference = valid_station_references[pd.notnull(valid_station_references)][0]
                        valid_station_names = self.metadata_in_memory[networkspeci]['station_name'][relevant_station_ind, :]
                        self.current_station_name = valid_station_names[pd.notnull(valid_station_names)][0]
                        current_lons = self.metadata_in_memory[networkspeci]['longitude'][relevant_station_ind, :]
                        self.current_lon = round(current_lons[pd.notnull(current_lons)][0], 2)
                        current_lats = self.metadata_in_memory[networkspeci]['latitude'][relevant_station_ind, :]
                        self.current_lat = round(current_lats[pd.notnull(current_lats)][0], 2)
                        
                        # put station data in pandas dataframe
                        to_pandas_dataframe(read_instance=self, canvas_instance=self, 
                                            networkspecies=[networkspeci], 
                                            station_index=relevant_station_ind, 
                                            data_range_min=self.data_range_min_station, 
                                            data_range_max=self.data_range_max_station)

                        # if have no valid data across data labels (no observations or experiments), then continue 
                        if not self.selected_station_data[networkspeci]:
                            continue
                        
                        # update data range min/maxes for station paradigm
                        if self.selected_station_data_min[networkspeci] < self.data_range_min_station[networkspeci]:
                            self.data_range_min_station[networkspeci] = copy.deepcopy(self.selected_station_data_min[networkspeci]) 
                        if self.selected_station_data_max[networkspeci] > self.data_range_max_station[networkspeci]:
                            self.data_range_max_station[networkspeci] = copy.deepcopy(self.selected_station_data_max[networkspeci])

                        # iterate through plots to make
                        for plot_type in station_plots_to_make:

                            # get zstat information from plot_type
                            zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)
                            # get base plot type (without stat and options)
                            if zstat:
                                base_plot_type = plot_type.split('-')[0] 
                            else:
                                base_plot_type = plot_type.split('_')[0] 

                            # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                            plot_options = plot_type.split('_')[1:]

                            # if have multispecies option: 
                            #   put multiple species data in pandas dataframe if first networkspeci being plotted and spatial_colocation is active 
                            #   if not, then continue
                            if 'multispecies' in plot_options: 
                                if (not made_networkspeci_station_plots) & (self.spatial_colocation):
                                    to_pandas_dataframe(read_instance=self, canvas_instance=self, 
                                                        networkspecies=self.networkspecies, 
                                                        station_index=relevant_station_ind, 
                                                        data_range_min=self.data_range_min_station, 
                                                        data_range_max=self.data_range_max_station)
                                    
                                    # update data range min/maxes for station paradigm
                                    for ns in self.networkspecies:
                                        if self.selected_station_data_min[ns] < self.data_range_min_station[ns]:
                                            self.data_range_min_station[ns] = copy.deepcopy(self.selected_station_data_min[ns]) 
                                        if self.selected_station_data_max[ns] > self.data_range_max_station[ns]:
                                            self.data_range_max_station[ns] = copy.deepcopy(self.selected_station_data_max[ns])
                                else:
                                    continue
                            
                            # make plot
                            print('Making station {2} for {3} ({0}/{1})'.format(i+1, 
                                                                                len(self.relevant_station_inds),
                                                                                plot_type, 
                                                                                self.current_station_name))                                
                            plot_indices = self.make_plot('station', plot_type, plot_options, networkspeci)

                            # get relevant axs and plot types per networkspeci / plot type
                            relevant_axs, relevant_data_labels = self.get_relevant_axs_per_networkspeci_plot_type(base_plot_type, 
                                                                                                                  [plot_indices[0][0]])

                            # do formatting 
                            self.plot.do_formatting(relevant_axs, relevant_data_labels, networkspeci,
                                                    base_plot_type, plot_type, plot_options, 'station')

                    # update variable now station plots have been made for a networkspecies
                    made_networkspeci_station_plots = True

            # remove subsection variables from memory (if have one, and not last subsection)
            if (len(self.child_subsection_names) > 0) & (subsection_ind != (len(self.subsections) - 1)):
                for k in self.subsection_opts:
                    try:
                        vars(self).pop(k)
                    except:
                        pass

            # update variable to keep track if have setup summary plot geometry yet for a subsection
            if made_networkspeci_summary_plots:
                summary_plot_geometry_setup = True

    def make_plot(self, plotting_paradigm, plot_type, plot_options, networkspeci):
        """ Function that calls making of any type of plot. """

        self.current_plot_ind = 0

        # create list to store index of saved plot information for plot_type
        # index is composed of nested list of [page_number, page_ind]
        plot_indices = []

        # get zstat information from plot_type
        zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)

        # get base plot type (without stat and options)
        if zstat:
            base_plot_type = plot_type.split('-')[0] 
        else:
            base_plot_type = plot_type.split('_')[0] 

        # get all data_labels for selected_station_data
        # if multispecies plot then ensure have data_labels from across all 
        if 'multispecies' in plot_options: 
            all_data_labels = []
            for ns in self.selected_station_data:
                all_data_labels.extend(list(self.selected_station_data[ns].keys()))
            _, idx = np.unique(all_data_labels, return_index=True)
            all_data_labels = np.array(all_data_labels)[np.sort(idx)]
        else:
            all_data_labels = list(self.selected_station_data[networkspeci].keys())

        # if are making bias plot, and have no valid experiment data then cannot make plot type
        if ('bias' in plot_options) & (len(all_data_labels) < 2):
            return plot_indices

        # get data ranges for plotting paradigm
        if plotting_paradigm == 'summary':
            data_range_min = self.data_range_min_summary[networkspeci]
            data_range_max = self.data_range_max_summary[networkspeci]
        elif plotting_paradigm == 'station':
            data_range_min = self.data_range_min_station[networkspeci]
            data_range_max = self.data_range_max_station[networkspeci]

        # iterate through all data arrays
        first_data_label = True
        for data_label in all_data_labels:
            
            # set how experiment should be referred to in heatmap/table
            if data_label == 'observations':
                data_label_legend = copy.deepcopy(data_label)
            else:    
                data_label_legend = self.experiments[data_label]

            # do not show experiment data with 'obs' option set
            if ('obs' in plot_options) & (data_label != 'observations'):
                continue

            # skip observations data label when plotting bias
            if (z_statistic_sign == 'bias') & (data_label == 'observations'):
                continue
                
            # map plots (1 plot per data array/s (1 array if absolute plot,
            # 2 arrays if making bias plot), per subsection)
            if base_plot_type == 'map':
                
                # get necessary data labels to plot
                if z_statistic_sign == 'bias':
                    z1 = 'observations'
                    z2 = copy.deepcopy(data_label)
                else:
                    z1 = copy.deepcopy(data_label)
                    z2 = ''

                # get relevant page/axis to plot on
                axis_ind = (self.current_plot_ind * len(self.subsections)) + self.subsection_ind
                relevant_page, page_ind, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, 
                                                                                     plot_type, axis_ind)

                # set axis title
                if relevant_axis.get_title() == '':
                    if self.n_stations == 1:
                        axis_title_label = '{}\n{} ({} station)'.format(data_label, self.subsection, self.n_stations)
                    else:
                        axis_title_label = '{}\n{} ({} stations)'.format(data_label, self.subsection, self.n_stations)
                    self.plot.set_axis_title(relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

                # set map extent ? 
                if self.map_extent:
                    self.plot.set_map_extent(relevant_axis)

                # do not make map if there is no valid data for relevant data label/s
                if z1 not in all_data_labels:
                    continue

                if (z2 != '') & (z2 not in all_data_labels): 
                    continue

                # calculate z statistic
                self.z_statistic, active_map_valid_station_inds = calculate_z_statistic(self, z1, z2, zstat, networkspeci)
                self.active_map_valid_station_inds = active_map_valid_station_inds

                # make map plot
                self.plot.make_map(relevant_axis, networkspeci, self.z_statistic, self.plot_characteristics[plot_type], 
                                   plot_options=plot_options, first_data_label=first_data_label)
                
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

                first_data_label = False

            # heatmap and table
            elif base_plot_type in ['heatmap', 'table']:

                # create nested dictionary to store statistical information across all subsections
                if plotting_paradigm == 'summary':
                    if zstat not in self.subsection_stats_summary:
                        self.subsection_stats_summary[zstat] = {}
                    if data_label_legend not in self.subsection_stats_summary[zstat]:
                        self.subsection_stats_summary[zstat][data_label_legend] = {}
                elif plotting_paradigm == 'station':
                    if self.current_station_reference not in self.subsection_stats_station:
                        self.subsection_stats_station[self.current_station_reference] = {}
                    if zstat not in self.subsection_stats_station[self.current_station_reference]:
                        self.subsection_stats_station[self.current_station_reference][zstat] = {}
                    if data_label_legend not in self.subsection_stats_station[self.current_station_reference][zstat]:
                        self.subsection_stats_station[self.current_station_reference][zstat][data_label_legend] = {}
                
                # add stat for current data array (if has been calculated correctly, otherwise append NaNs)                                         
                if data_label in self.selected_station_data[networkspeci]:
                    if len(self.selected_station_data[networkspeci][data_label]['pandas_df']['data'].dropna()) > 0:
                        data_to_add = self.selected_station_data[networkspeci][data_label]['all'][zstat][0]
                    else:
                        data_to_add = np.NaN
                else:
                    data_to_add = np.NaN
                
                if plotting_paradigm == 'summary':
                    if self.subsection not in self.subsection_stats_summary[zstat][data_label_legend]:
                        self.subsection_stats_summary[zstat][data_label_legend][self.subsection] = data_to_add
                elif plotting_paradigm == 'station':
                    if self.subsection not in self.subsection_stats_station[self.current_station_reference][zstat][data_label_legend]:
                        self.subsection_stats_station[self.current_station_reference][zstat][data_label_legend][self.subsection] = data_to_add

            # other plots (1 plot per subsection with multiple data arrays for summary paradigm, 1 plot per subsection per station for station paradigm)
            elif base_plot_type != 'statsummary':
                # get relevant axis to plot on
                if plotting_paradigm == 'summary':
                    if 'individual' in plot_options:
                        if ((base_plot_type == 'scatter') or ('bias' in plot_options) or 
                            (z_statistic_sign == 'bias')):
                                axis_ind = (self.current_plot_ind + self.subsection_ind + (len(self.experiments) - 1) * self.subsection_ind)
                        else:
                            axis_ind = (self.current_plot_ind + self.subsection_ind + len(self.experiments) * self.subsection_ind)
                    else:
                        axis_ind = self.subsection_ind
                    station_inds = copy.deepcopy(self.relevant_station_inds) 
                elif plotting_paradigm == 'station':
                    if 'individual' in plot_options:
                        if ((base_plot_type == 'scatter') or ('bias' in plot_options) or 
                            (z_statistic_sign == 'bias')):
                            axis_ind = (self.current_plot_ind + self.station_ind + (len(self.experiments) - 1) * self.station_ind)
                        else:
                            axis_ind = (self.current_plot_ind + self.station_ind + len(self.experiments) * self.station_ind)
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

                # axis xlabel is empty or == ?
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

                # periodic plots
                if base_plot_type in ['periodic', 'periodic-violin']:
                    # skip observational array if plotting bias stat
                    if (z_statistic_sign == 'bias') & (data_label == 'observations'):
                        continue
                        
                    # determine if have some data to plot
                    plot_validity = False
                    if 'multispecies' in plot_options:
                        for ns in self.selected_station_data:
                            if data_label in self.selected_station_data[ns]:
                                if len(self.selected_station_data[ns][data_label]['pandas_df']['data'].dropna()) > 0:
                                    plot_validity = True
                    else:
                        if data_label in self.selected_station_data[networkspeci]:
                            if len(self.selected_station_data[networkspeci][data_label]['pandas_df']['data'].dropna()) > 0:
                                plot_validity = True

                    if plot_validity:
                        if base_plot_type == 'periodic':
                            self.plot.make_periodic(relevant_axis, networkspeci, data_label, 
                                                    self.plot_characteristics[plot_type], zstat=zstat, 
                                                    plot_options=plot_options, 
                                                    first_data_label=first_data_label)
                        elif base_plot_type == 'periodic-violin':
                            self.plot.make_periodic(relevant_axis, networkspeci, data_label, 
                                                    self.plot_characteristics[plot_type], 
                                                    plot_options=plot_options, 
                                                    first_data_label=first_data_label)
                        
                        # save plot information for later formatting 
                        self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(data_label)
                        plot_index = [relevant_page, page_ind]
                        if plot_index not in plot_indices:
                            plot_indices.append(plot_index)
                        first_data_label = False

                        # turn relevant axes on
                        for relevant_temporal_resolution in self.relevant_temporal_resolutions:
                            relevant_axis[relevant_temporal_resolution].axis('on')
                            relevant_axis[relevant_temporal_resolution].set_visible(True)

                            #get references to periodic label annotations made, and then show them
                            annotations = [child for child in relevant_axis[relevant_temporal_resolution].get_children() if isinstance(child, matplotlib.text.Annotation)]
                            # hide annotations
                            for annotation in annotations:
                                annotation.set_visible(True)

                # other plot types (except heatmap, table and statsummary) 
                else:
                    # skip observational array for bias/scatter plots
                    if data_label == 'observations' and ('bias' in plot_options or base_plot_type == 'scatter'):
                        continue
                    else:
                        func = getattr(self.plot, 'make_{}'.format(base_plot_type))
                    # determine if have some data to plot
                    plot_validity = False
                    if 'multispecies' in plot_options:
                        for ns in self.selected_station_data:
                            if data_label in self.selected_station_data[ns]:
                                if len(self.selected_station_data[ns][data_label]['pandas_df']['data'].dropna()) > 0:
                                    plot_validity = True
                    else:
                        if data_label in self.selected_station_data[networkspeci]:
                            if len(self.selected_station_data[networkspeci][data_label]['pandas_df']['data'].dropna()) > 0:
                                plot_validity = True

                    if plot_validity:
                        if base_plot_type == 'metadata':
                            func(relevant_axis, networkspeci, data_label, self.plot_characteristics[plot_type], 
                                plot_options=plot_options, first_data_label=first_data_label, station_inds=station_inds) 
                        elif base_plot_type == 'distribution':
                            func(relevant_axis, networkspeci, data_label, self.plot_characteristics[plot_type], 
                                plot_options=plot_options, first_data_label=first_data_label,
                                data_range_min=data_range_min, 
                                data_range_max=data_range_max) 
                        else:
                            func(relevant_axis, networkspeci, data_label, self.plot_characteristics[plot_type], 
                                plot_options=plot_options, first_data_label=first_data_label) 
                        
                        # save plot information for later formatting
                        self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(data_label)
                        plot_index = [relevant_page, page_ind]
                        if plot_index not in plot_indices:
                            plot_indices.append(plot_index)
                        first_data_label = False

                        # turn axis on
                        relevant_axis.axis('on')
                        relevant_axis.set_visible(True)

                first_data_label = False

            # iterate number of plots made for current type of plot 
            self.current_plot_ind += 1     

        # then make plot heatmap / table / statsummary plot
        if base_plot_type in ['heatmap', 'table', 'statsummary']:

            # get relevant axis to plot on
            create_plot = False
            if plotting_paradigm == 'summary':
                if base_plot_type == 'statsummary':
                    axis_ind = self.subsection_ind
                    create_plot = True
                else: 
                    axis_ind = 0
                    if self.subsection_ind == (len(self.subsections) - 1):
                        create_plot = True
            elif plotting_paradigm == 'station':
                axis_ind = self.station_ind
                create_plot = True
            
            # make plot ?
            if create_plot:
                relevant_page, page_ind, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, 
                                                                                     plot_type, axis_ind)
            
                if base_plot_type in ['heatmap', 'table']:
               
                    # convert subsection_stats dicts to dataframe, with subsection names as indices
                    if plotting_paradigm == 'summary':
                        stats_to_plot = copy.deepcopy(self.subsection_stats_summary)
                        for data_label in stats_to_plot[zstat].keys():
                            stats_per_data_label = []
                            for subsection, stat in stats_to_plot[zstat][data_label].items():
                                stats_per_data_label.append(stat)
                        stats_to_plot[zstat][data_label] = stats_per_data_label
                        stats_df = pd.DataFrame(data=stats_to_plot[zstat],
                                                index=self.subsections)
                    elif plotting_paradigm == 'station':
                        if self.current_station_reference not in self.subsection_stats_station:
                            stats_df = pd.DataFrame()
                        else:
                            stats_to_plot = copy.deepcopy(self.subsection_stats_station)
                            for data_label in stats_to_plot[self.current_station_reference][zstat].keys():
                                stats_per_data_label = []
                                for subsection, stat in stats_to_plot[self.current_station_reference][zstat][data_label].items():
                                    stats_per_data_label.append(stat)
                            stats_to_plot[self.current_station_reference][zstat][data_label] = stats_per_data_label
                            stats_df = pd.DataFrame(data=self.subsection_stats_station[self.current_station_reference][zstat],
                                                    index=[self.subsection])
        
                elif base_plot_type == 'statsummary':
                    # create structure to store data for statsummary plot
                    if 'bias' in plot_options:
                        relevant_zstats = self.plot_characteristics[plot_type]['experiment_bias']
                        relevant_data_labels = [data_label for data_label in all_data_labels if data_label != 'observations']
                    else:
                        relevant_zstats = self.plot_characteristics[plot_type]['basic']
                        relevant_data_labels = self.selected_station_data[networkspeci]
                        
                    stats_df = {relevant_zstat:[] for relevant_zstat in relevant_zstats}
                    for data_label in relevant_data_labels:
                        for relevant_zstat in relevant_zstats:                    
                            # if relevant stat is expbias stat, then ensure temporal colocation is active
                            # otherwise set value as NaN
                            if (relevant_zstat in expbias_stats) & (not self.temporal_colocation):
                                stat_val = np.NaN
                            else:
                                stat_val = self.selected_station_data[networkspeci][data_label]['all'][relevant_zstat][0]
                            stats_df[relevant_zstat].append(stat_val)
                    stats_df = pd.DataFrame(data=stats_df, index=relevant_data_labels)

                # set axis title
                if relevant_axis.get_title() == '':
                    if plotting_paradigm == 'summary':
                        if base_plot_type == 'statsummary':
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

                # turn off relevant axis if dataframe is empty or all NaN
                if (len(stats_df.index) > 0) & (not stats_df.isnull().values.all()):
                    # make plot
                    if base_plot_type == 'statsummary':
                        func = getattr(self.plot, 'make_table')
                        func(relevant_axis, stats_df, self.plot_characteristics[plot_type], plot_options=plot_options, 
                             statsummary=True)
                        # save plot information for later formatting
                        self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].extend(stats_df.index.tolist())

                    else:
                        func = getattr(self.plot, 'make_{}'.format(base_plot_type))
                        func(relevant_axis, stats_df, self.plot_characteristics[plot_type], plot_options=plot_options)
                        # save plot information for later formatting
                        self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].extend(stats_df.columns.tolist())

                    plot_index = [relevant_page, page_ind]
                    if plot_index not in plot_indices:
                        plot_indices.append(plot_index)

                    # turn axis on
                    relevant_axis.axis('on')
                    relevant_axis.set_visible(True)

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
        
def main(**kwargs):
    """ Main function when running offine reports. """
   
    ProvidentiaOffline(**kwargs)
