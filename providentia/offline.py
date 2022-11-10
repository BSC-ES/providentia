import os
import sys
import json
import copy

import numpy as np
import pandas as pd
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
    """Run Providentia offline reports"""

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
        # and merge to existing GHOST observational data dict if we have the path
        if self.nonghost_root is not None:
            nonghost_observation_data = aux.get_nonghost_observational_tree(self, nonghost_observation_data_json)
            self.all_observation_data = {**self.all_observation_data, **nonghost_observation_data}

        # initialise DataReader class
        self.datareader = DataReader(self)

        # initialise Plot class
        self.plot = Plot(read_instance=self, canvas_instance=self)

        # iterate through configuration sections
        for filename, section in zip(self.filenames, self.parent_section_names):

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
            vars(self).update({(k, provconf.parse_parameter(k, val)) for k, val in self.section_opts.items()})

            # now all variables have been parsed, check validity of those, throwing errors where necessary
            provconf.check_validity()

            # set some key configuration variables
            self.relevant_temporal_resolutions = aux.get_relevant_temporal_resolutions(self.resolution)
            self.data_labels = ['observations'] + list(self.experiments.keys())

            # get valid observations in date range
            aux.get_valid_obs_files_in_date_range(self, self.start_date, self.end_date)

            # update available experiments for selected fields
            aux.get_valid_experiments(self, self.start_date, self.end_date, self.resolution,
                                      self.network, self.species)

            # read data
            self.datareader.read_setup(['reset'])

            # initialise/reinitialise structures to store % data representativity, data periods and metadata
            aux.init_representativity(self)
            aux.init_period(self)
            aux.init_metadata(self)

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

            # remove section variables from memory
            for k in self.section_opts:
                try:
                    vars(self).pop(k)
                except:
                    pass

    def start_pdf(self):

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

            # get list of all networks and species strings
            self.networkspecies = ['{}|{}'.format(network,speci) for network, speci in zip(self.network, self.species)]

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
            self.make_plots_per_subsection(summary_plots_to_make_nondist, station_plots_to_make_nondist, plot_geometry_setup=True)

            # make all plots per subsection (for only distribution plot types --> done so to calclulate data ranges across subsections first)
            summary_plots_to_make_dist = [plot_type for plot_type in self.summary_plots_to_make if 'distribution' in plot_type]
            station_plots_to_make_dist = [plot_type for plot_type in self.station_plots_to_make if 'distribution' in plot_type]
            self.make_plots_per_subsection(summary_plots_to_make_dist, station_plots_to_make_dist)
            
            # do formatting to axes per networkspeci
            # create colourbars
            # harmonise xy limits(not for map, heatmap or table, or when xlim and ylim defined)
            # log axes
            # annotation
            # regression line
            # smooth line

            # remove header from plot characteristics dictionary
            if 'header' in list(self.plot_characteristics.keys()):
                del self.plot_characteristics['header']

            # iterate through networks and species
            for networkspeci_ii, networkspeci in enumerate(self.networkspecies): 

                # iterate through plot types
                for plot_type in list(self.plot_characteristics.keys()):

                    # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                    plot_options = plot_type.split('_')[1:]

                    # if a multispecies plot is wanted then this is only made when networkspeci_ii == 0
                    if ('multispecies' in plot_options) & (networkspeci_ii != 0):
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
                        paradigm_pages['summary'] = self.summary_pages[plot_type][networkspeci]
                        # for multispecies plots of specific stations, spatial colocation needs to be on.
                        # if not, the plots will not exist so handle this
                        if ('multispecies' in plot_options) & (not self.spatial_colocation):
                            paradigm_pages['station'] = []
                        else:
                            if plot_type in self.station_plots_to_make:
                                for subsection in self.subsections:
                                    paradigm_pages['station'].extend(self.station_pages[plot_type][networkspeci][subsection])

                    elif self.report_summary:
                        paradigm_pages['summary'] = self.summary_pages[plot_type][networkspeci]
                    
                    elif self.report_stations:
                        # for multispecies plots of specific stations, spatial colocation needs to be on.
                        # if not, the plots will not exist so handle this
                        if ('multispecies' in plot_options) & (not self.spatial_colocation):
                            paradigm_pages['station'] = []
                        else:
                            if plot_type in self.station_plots_to_make:
                                for subsection in self.subsections:
                                    paradigm_pages['station'].extend(self.station_pages[plot_type][networkspeci][subsection])

                    # iterate through paradigm pages
                    for plotting_paradigm, relevant_pages in paradigm_pages.items():                        
                        if len(relevant_pages) == 0:
                            continue
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

                        # iterate through all relevant axes for plot type in paradigm
                        for relevant_ax_ii, relevant_ax in enumerate(relevant_axs):

                            # log axes?
                            if 'logx' in plot_options:
                                log_validity = self.plot.log_validity(relevant_ax, 'logx')
                                if log_validity:
                                    self.plot.log_axes(relevant_ax, 'logx', base_plot_type, self.plot_characteristics[plot_type])
                                else:
                                    msg = "Warning: It is not possible to log the x-axis "
                                    msg += "in {0} with negative values.".format(plot_type)
                                    print(msg)
                            if 'logy' in plot_options:
                                log_validity = self.plot.log_validity(relevant_ax, 'logy')
                                if log_validity:
                                    self.plot.log_axes(relevant_ax, 'logy', base_plot_type, self.plot_characteristics[plot_type])
                                else:
                                    msg = "Warning: It is not possible to log the y-axis "
                                    msg += "in {0} with negative values.".format(plot_type)
                                    print(msg)

                            # annotation
                            if 'annotate' in plot_options:
                                if base_plot_type not in ['heatmap']:
                                    self.plot.annotation(relevant_ax, networkspeci, relevant_data_labels[relevant_ax_ii], 
                                                        base_plot_type, self.plot_characteristics[plot_type],
                                                        plot_options=plot_options)
                                    # annotate on first axis
                                    if base_plot_type in ['periodic', 'periodic-violin']:
                                        break

                            # regression line
                            if 'regression' in plot_options:
                                self.plot.linear_regression(relevant_ax, networkspeci, 
                                                            relevant_data_labels[relevant_ax_ii], 
                                                            base_plot_type,
                                                            self.plot_characteristics[plot_type], 
                                                            plot_options=plot_options)

                            # smooth line
                            if 'smooth' in plot_options:
                                self.plot.smooth(relevant_ax, networkspeci, relevant_data_labels[relevant_ax_ii], 
                                                 base_plot_type, self.plot_characteristics[plot_type], 
                                                 plot_options=plot_options)

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
            # save page figures
            print('\nWRITING PDF')
            for page in self.plot_dictionary:
                fig = self.plot_dictionary[page]['fig']
                self.pdf.savefig(fig, dpi=self.dpi)
                plt.close(fig)

    def setup_plot_geometry(self, plotting_paradigm, networkspeci, networkspeci_ii):
        """setup plotting geometry for summary or station specific plots, per network/species"""

        # depending on plot type set plots to make
        if plotting_paradigm == 'summary':
            self.plots_to_make = self.summary_plots_to_make
        elif plotting_paradigm == 'station':
            self.plots_to_make = self.station_plots_to_make

        # iterate through plot types to make
        for plot_type in self.plots_to_make:

            # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
            plot_options = plot_type.split('_')[1:]
            
            # if a multispecies plot is wanted then only make this when networkspeci_ii == 0
            # if making a multispecies plot per specific station, spatial colocation must be also active
            if 'multispecies' in plot_options:
                if plotting_paradigm == 'summary':
                    if networkspeci_ii != 0:
                        continue
                elif plotting_paradigm == 'station':
                    if (networkspeci_ii != 0) or (not self.spatial_colocation):
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
                    plot_characteristics['page_title']['t'] = '{} (Summary)\n{}'.format(plot_characteristics['page_title']['t'], networkspeci) 
            elif plotting_paradigm == 'station':
                if 'multispecies' in plot_options:
                    plot_characteristics['page_title']['t'] = '{} (Per Station)\n{}\nmultispecies'.format(plot_characteristics['page_title']['t'], self.subsection) 
                else:
                    plot_characteristics['page_title']['t'] = '{} (Per Station)\n{}\n{}'.format(plot_characteristics['page_title']['t'], self.subsection, networkspeci) 

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
                            grid_dict['hour'].set_axis_off()
                            grid_dict['dayofweek'].set_axis_off()
                            grid_dict['month'].set_axis_off()
                            
                            # turn on all axes that will be plotted on
                            for relevant_temporal_resolution in self.relevant_temporal_resolutions:
                                grid_dict[relevant_temporal_resolution].set_axis_on()
                            
                            # format axis
                            self.plot.format_axis(grid_dict, base_plot_type, plot_characteristics, set_extent=False, relevant_temporal_resolutions=self.relevant_temporal_resolutions)

                        # rest of plot types
                        else:
                            self.plot_dictionary[page_n]['axs'].append({'handle':ax, 'data_labels':[]})
                            
                            # format axis 
                            self.plot.format_axis(ax, base_plot_type, plot_characteristics, set_extent=False)

                    # no more plots to make on page? 
                    # then turn off unneeded axes
                    else:
                        ax.axis('off')
                        ax.set_visible(False)
                        ax.grid(False)

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

                # add colourbar axis to plot dictionary?
                if 'cb' in plot_characteristics_vars:
                    self.plot_dictionary[page_n]['cb_ax'] = fig.add_axes(plot_characteristics['cb']['position'])
                    self.plot_dictionary[page_n]['cb_ax'].set_rasterized(True)
                            
                # add current page number
                if plotting_paradigm == 'summary':
                    self.summary_pages[plot_type][networkspeci].append(page_n)
                elif plotting_paradigm == 'station':
                    self.station_pages[plot_type][networkspeci][self.subsection].append(page_n)
                
                # add to total number of pages
                self.n_total_pages += 1

    def make_plots_per_subsection(self, summary_plots_to_make, station_plots_to_make, plot_geometry_setup=False):
        """Function that calls making of all plots per subsection"""

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
                vars(self).update({(k, provconf.parse_parameter(k, val)) for k, val in self.subsection_opts.items()})

                # now all variables have been parsed, check validity of those, throwing errors where necessary
                provconf.check_validity()

            # if have no experiments, force temporal colocation to be False
            if len(self.experiments) == 0:
                self.temporal_colocation = False    

            # update fields available for filtering
            aux.update_representativity_fields(self)
            aux.representativity_conf(self)
            aux.update_period_fields(self)
            aux.period_conf(self)
            aux.update_metadata_fields(self)
            aux.metadata_conf(self)
            
            # filter dataset for current subsection
            print('\nFiltering data for {} subsection'.format(self.subsection))
            DataFilter(self)

            # iterate through all desired plots, making each one (summary or station specific plots)
            print('Making {} subsection plots'.format(self.subsection))  
            
            # iterate through networks and species, creating plots
            self.n_total_pages = 0

            for networkspeci_ii, networkspeci in enumerate(self.networkspecies):
                
                # get valid station inds for networkspeci 
                if self.temporal_colocation and len(self.data_labels) > 1:
                    self.relevant_station_inds = self.valid_station_inds_temporal_colocation[networkspeci]['observations']
                else:
                    self.relevant_station_inds = self.valid_station_inds[networkspeci]['observations']  

                # get N stations for networkspeci
                self.n_stations = len(self.relevant_station_inds)

                # make summary plots?
                if self.report_summary:

                    if networkspeci_ii == 0:
                        # get median timeseries across data from filtered data, and place it pandas dataframe
                        to_pandas_dataframe(read_instance=self, canvas_instance=self, 
                                            networkspecies=self.networkspecies)

                        # update data range min/maxes for summary paradigm
                        for ns in self.networkspecies:
                            if self.selected_station_data_min[ns] < self.data_range_min_summary[ns]:
                                self.data_range_min_summary[ns] = copy.deepcopy(self.selected_station_data_min[ns]) 
                            if self.selected_station_data_max[ns] > self.data_range_max_summary[ns]:
                                self.data_range_max_summary[ns] = copy.deepcopy(self.selected_station_data_max[ns])

                    if subsection_ind == 0:
                        # update plot characteristics
                        self.plot.set_plot_characteristics(summary_plots_to_make)
                    
                        # setup plotting geometry for summary plots per networkspeci (for all subsections)
                        if plot_geometry_setup:
                            self.setup_plot_geometry('summary', networkspeci, networkspeci_ii)
                    
                    # iterate through plots to make
                    for plot_type in summary_plots_to_make:
                        
                        #get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                        plot_options = plot_type.split('_')[1:]

                        # if a multispecies plot is wanted then only make this when networkspeci_ii == 0
                        if ('multispecies' in plot_options) & (networkspeci_ii != 0):
                            continue

                        # make plot
                        print('Making summary {0}'.format(plot_type))
                        self.make_plot('summary', plot_type, plot_options, networkspeci)

                    # update N total pages 
                    self.n_total_pages = len(self.plot_dictionary)

                # make station specific plots?
                if self.report_stations:               

                    # update plot characteristics
                    self.plot.set_plot_characteristics(station_plots_to_make)

                    #  setup plotting geometry for station plots per networkspeci (for one subsection)
                    if plot_geometry_setup:
                        self.setup_plot_geometry('station', networkspeci, networkspeci_ii)
                        
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
                        
                        # update data range min/maxes for station paradigm
                        for ns in self.networkspecies:
                            if self.selected_station_data_min[ns] < self.data_range_min_station[ns]:
                                self.data_range_min_station[ns] = copy.deepcopy(self.selected_station_data_min[ns]) 
                            if self.selected_station_data_max[ns] > self.data_range_max_station[ns]:
                                self.data_range_max_station[ns] = copy.deepcopy(self.selected_station_data_max[ns])

                        # iterate through plots to make
                        for plot_type in station_plots_to_make:

                            # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                            plot_options = plot_type.split('_')[1:]

                            # if have multispecies option: 
                            #   put multiple species data in pandas dataframe if networkspeci_ii == 0 and spatial_colocation is active 
                            #   if not, then continue
                            if 'multispecies' in plot_options: 
                                if (networkspeci_ii == 0) & (self.spatial_colocation):
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
                            self.make_plot('station', plot_type, plot_options, networkspeci)

            # remove subsection variables from memory (if have one)
            if len(self.child_subsection_names) > 0:
                for k in self.subsection_opts:
                    try:
                        vars(self).pop(k)
                    except:
                        pass

    def make_plot(self, plotting_paradigm, plot_type, plot_options, networkspeci):
        """Function that calls making of any type of plot"""

        self.current_plot_ind = 0
        
        # get zstat information from plot_type
        zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)

        # get base plot type (without stat and options)
        if zstat:
            base_plot_type = plot_type.split('-')[0] 
        else:
            base_plot_type = plot_type.split('_')[0] 

        # if are making bias plot, and have no valid experiment data then cannot make plot type
        if ('bias' in plot_options) & (len(self.selected_station_data[networkspeci]) < 2):
            return

        # get data ranges for plotting paradigm
        if plotting_paradigm == 'summary':
            data_range_min = self.data_range_min_summary[networkspeci]
            data_range_max = self.data_range_max_summary[networkspeci]
        elif plotting_paradigm == 'station':
            data_range_min = self.data_range_min_station[networkspeci]
            data_range_max = self.data_range_max_station[networkspeci]

        # iterate through all data arrays
        first_data_label = True 
        for n_data_label, data_label in enumerate(self.selected_station_data[networkspeci].keys()):

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
                relevant_page, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, plot_type, 
                                                                           axis_ind)

                # calculate number of created axis
                n_axes_plot_type = 0
                for page in self.plot_dictionary:
                    # get number of axes until current page
                    if page < relevant_page:
                        # avoid axes from other plot types
                        if plotting_paradigm == 'station':
                            if page >= self.station_pages[plot_type][networkspeci][self.subsection][0]:
                                for axis in self.plot_dictionary[page]['axs']:
                                    n_axes_plot_type += 1
                        elif plotting_paradigm == 'summary':
                            if page >= self.summary_pages[plot_type][networkspeci][0]:
                                for axis in self.plot_dictionary[page]['axs']:
                                    n_axes_plot_type += 1

                # get corresponding page index given an axis index
                if axis_ind >= len(self.plot_dictionary[relevant_page]['axs']):
                    page_ind = axis_ind - n_axes_plot_type
                else:
                    page_ind = axis_ind 
                
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

                # make map if there are data
                if not self.selected_station_data[networkspeci]:
                    relevant_axis.set_axis_off()
                    relevant_axis.set_visible(False)
                else:
                    # calculate z statistic
                    self.z_statistic, active_map_valid_station_inds = calculate_z_statistic(self, z1, z2, zstat, networkspeci)
                    self.active_map_valid_station_inds = active_map_valid_station_inds

                    # make map plot
                    self.plot.make_map(relevant_axis, networkspeci, self.z_statistic, self.plot_characteristics[plot_type], 
                                       plot_options=plot_options)
                    if z2 == '':
                        self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(z1)
                    else:
                        self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(z2)

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
                    if len(self.selected_station_data[networkspeci][data_label]['pandas_df']['data']) > 0:
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
                relevant_page, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, plot_type, axis_ind)

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
                            xlabel = self.measurement_units[networkspeci.split('|')[-1]]
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
                            ylabel = copy.deepcopy(ylabel_units)
                    else:
                        if 'ylabel' in self.plot_characteristics[plot_type]:
                            if self.plot_characteristics[plot_type]['ylabel']['ylabel'] == 'measurement_units':
                                ylabel = self.measurement_units[networkspeci.split('|')[-1]]
                            else:
                                ylabel = self.plot_characteristics[plot_type]['ylabel']['ylabel']
                        else:
                            ylabel = ''
                    # set ylabel
                    if ylabel != '':
                        self.plot.set_axis_label(relevant_axis, 'y', ylabel, self.plot_characteristics[plot_type])

                # make plot if there is data
                if not self.selected_station_data[networkspeci]:
                    # relevant axis is a dict of the different temporal aggregations in some cases (e.g. periodic plots)
                    if isinstance(relevant_axis, dict):
                        for temporal_aggregation_resolution, temporal_aggregation_relevant_axis in relevant_axis.items():
                            temporal_aggregation_relevant_axis.set_axis_off()
                    else:
                        relevant_axis.set_axis_off()
                else: 
                    # calculate number of created axis
                    n_axes_plot_type = 0
                    for page in self.plot_dictionary:
                        # get number of axes until current page
                        if page < relevant_page:
                            # avoid axes from other plot types
                            if plotting_paradigm == 'station':
                                if page >= self.station_pages[plot_type][networkspeci][self.subsection][0]:
                                    for axis in self.plot_dictionary[page]['axs']:
                                        n_axes_plot_type += 1
                            elif plotting_paradigm == 'summary':
                                if page >= self.summary_pages[plot_type][networkspeci][0]:
                                    for axis in self.plot_dictionary[page]['axs']:
                                        n_axes_plot_type += 1

                    # get corresponding page index given an axis index
                    if axis_ind >= len(self.plot_dictionary[relevant_page]['axs']):
                        page_ind = axis_ind - n_axes_plot_type
                    else:
                        page_ind = axis_ind 
                    
                    # periodic plots
                    if base_plot_type in ['periodic', 'periodic-violin']:
                        # skip observational array if plotting bias stat
                        if (z_statistic_sign == 'bias') & (data_label == 'observations'):
                            continue
                            
                        if data_label in self.selected_station_data[networkspeci]:
                            if len(self.selected_station_data[networkspeci][data_label]['pandas_df']['data']) > 0:
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
                                self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(data_label)
                                first_data_label = False

                    # other plot types (except heatmap, table and statsummary) 
                    else:
                        # skip observational array for bias/scatter plots
                        if data_label == 'observations' and ('bias' in plot_options or base_plot_type == 'scatter'):
                            continue
                        else:
                            func = getattr(self.plot, 'make_{}'.format(base_plot_type))
                        if data_label in self.selected_station_data[networkspeci]:
                            if len(self.selected_station_data[networkspeci][data_label]['pandas_df']['data']) > 0:
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
                                first_data_label = False
                                self.plot_dictionary[relevant_page]['axs'][page_ind]['data_labels'].append(data_label)        

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
                relevant_page, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, plot_type, 
                                                                       axis_ind)
            
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
                        relevant_data_labels = [data_label for data_label in self.selected_station_data[networkspeci] if data_label != 'observations']
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
                if (len(stats_df.index) == 0) or (stats_df.isnull().values.all()):
                    relevant_axis.set_axis_off()
                else:
                    # make plot
                    if base_plot_type == 'statsummary':
                        func = getattr(self.plot, 'make_table')
                        func(relevant_axis, stats_df, self.plot_characteristics[plot_type], plot_options=plot_options, 
                             statsummary=True)
                    else:
                        func = getattr(self.plot, 'make_{}'.format(base_plot_type))
                        func(relevant_axis, stats_df, self.plot_characteristics[plot_type], plot_options=plot_options)

    def get_relevant_page_axis(self, plotting_paradigm, networkspeci, plot_type, axis_ind):
        """get relevant page and axis for current plot type/subsection/axis index"""

        # get axes associated with plot type
        if plotting_paradigm == 'summary':
            relevant_pages = self.summary_pages[plot_type][networkspeci]
        elif plotting_paradigm == 'station':
            relevant_pages = self.station_pages[plot_type][networkspeci][self.subsection]

        all_relevant_pages = []
        relevant_axes = []         
        for relevant_page in relevant_pages:
            relevant_axes.extend(self.plot_dictionary[relevant_page]['axs'])
            all_relevant_pages.extend([relevant_page]*len(self.plot_dictionary[relevant_page]['axs']))

        return all_relevant_pages[axis_ind], relevant_axes[axis_ind]['handle']

def main(**kwargs):
    """Main function when running offine reports"""
    ProvidentiaOffline(**kwargs)
