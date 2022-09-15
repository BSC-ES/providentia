import os
import sys
import json
import copy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

from .read_aux import drop_nans
from .read import DataReader
from .filter import DataFilter
from .plot import Plot
from .statistics import to_pandas_dataframe
from .statistics import calculate_z_statistic
from .statistics import generate_colourbar
from .statistics import get_z_statistic_info
from .configuration import ProvConfiguration
from .init_standards import InitStandards
from providentia import aux

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
 
class ProvidentiaOffline(ProvConfiguration, InitStandards):
    """Run Providentia offline reports"""

    def __init__(self, **kwargs):
        ProvConfiguration.__init__(self, **kwargs)

        print("Starting Providentia offline...")

        # update from config file
        if 'config' in kwargs:
            aux.load_conf(self, kwargs['config'])
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            sys.exit(error)

        # update self from command line arguments
        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in kwargs.items()})

        # init GHOST standards
        InitStandards.__init__(self, ghost_root=self.ghost_root, ghost_version=self.ghost_version)
        
        # create dictionary of all available observational GHOST data
        self.all_observation_data = aux.get_ghost_observational_tree(self)

        # load dictionary with non-GHOST esarchive files to read
        nonghost_observation_data_json = json.load(open(os.path.join(CURRENT_PATH, 'conf/nonghost_files.json')))
        # and merge to existing GHOST observational data dict if we have the path
        if self.nonghost_root is not None:
            nonghost_observation_data = aux.get_nonghost_observational_tree(self, nonghost_observation_data_json)
            self.all_observation_data = {**self.all_observation_data, **nonghost_observation_data}

        # load necessary dictionaries
        self.report_plots = json.load(open(os.path.join(
            CURRENT_PATH, 'conf/report_plots.json')))

        # initialise DataReader class
        self.datareader = DataReader(self)

        # check for self defined plot characteristics file
        if self.plot_characteristics_filename == '':
            self.plot_characteristics_filename = os.path.join(
                CURRENT_PATH, 'conf/plot_characteristics_offline.json')
        self.plot_characteristics_templates = json.load(open(self.plot_characteristics_filename))
        self.plot_characteristics = {}

        # add general plot characteristics to self
        for k, val in self.plot_characteristics_templates['general'].items():
            setattr(self, k, val)

        # initialise Plot class
        self.plot = Plot(read_instance=self, canvas_instance=self)

        # iterate through configuration sections
        for section in self.parent_section_names:

            # remove old parameters
            other_sections = [value for value in self.all_sections if value != section]
            for other_section in other_sections:
                for k in self.sub_opts[other_section]:
                    try:
                        vars(self).pop(k)
                    except:
                        pass

            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[section]

            # update self with section variables
            vars(self).update({(k, self.parse_parameter(k, val)) for k, val in self.section_opts.items()})

            # update self from command line arguments
            vars(self).update({(k, self.parse_parameter(k, val)) for k, val in kwargs.items()})

            # get key configuration variables
            aux.get_parameters(self) 
            self.qa = aux.which_qa(self)
            self.flags = aux.which_flags(self)
            self.experiments = aux.get_experiments(self)
            self.relevant_temporal_resolutions = aux.get_relevant_temporal_resolutions(self.resolution)
            self.data_labels = ['observations'] + list(self.experiments.keys())

            # if have no experiments, force temporal colocation to be False
            if len(self.experiments) == 0:
                self.temporal_colocation = False
                self.defaults['temporal_colocation'] = False         
                for k in self.sub_opts.keys():
                    self.sub_opts[k]['temporal_colocation'] = False

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

            # initialise first page number to plot
            self.current_page_n = 1

            # initialise dictionaries to store relevant page numebrs
            if self.report_summary:
                self.summary_pages = {}
            if self.report_stations:
                self.station_pages = {}

            # get list of all networks and species strings
            networkspecies = ['{}-{}'.format(network,speci) for network, speci in zip(self.network, self.species)]

            # get subsections
            self.child_subsection_names = [subsection_name for subsection_name in self.subsection_names 
                                           if self.section == subsection_name.split('|')[0]]
            if len(self.child_subsection_names) > 0:
                self.subsections = self.child_subsection_names
            else:
                self.subsections = [self.section]

            # make header
            self.plot.set_plot_characteristics(['header'])
            self.plot.make_header(self.pdf, self.plot_characteristics['header'])
        
            # iterate through subsections
            for subsection_ind, subsection in enumerate(self.subsections):
 
                self.subsection_ind = subsection_ind
                self.subsection = subsection

                # update the conf options for this subsection
                if len(self.child_subsection_names) > 0:
                    for section in self.all_sections:
                        for k in self.sub_opts[section]:
                            try:
                                vars(self).pop(k)
                            except:
                                continue
                    vars(self).update({(k, self.parse_parameter(k, val)) for k, val in
                                      self.sub_opts[self.subsection].items()})

                # update fields available for filtering
                aux.update_representativity_fields(self)
                aux.representativity_conf(self)
                aux.update_period_fields(self)
                aux.period_conf(self)
                aux.update_metadata_fields(self)
                aux.metadata_conf(self)
                
                # filter dataset for current subsection
                print('\nFiltering Data for {} Subsection'.format(self.subsection))
                DataFilter(self)

                # iterate through all desired plots, making each one (summary or station specific plots)
                print('Making {} Subsection Plots'.format(self.subsection))  

                # iterate through networks and species
                for networkspeci_ii, networkspeci in enumerate(networkspecies):

                    # make summary plots?
                    if self.report_summary:

                        if networkspeci_ii == 0:
                            # get median timeseries across data from filtered data, and place it pandas dataframe
                            to_pandas_dataframe(read_instance=self, canvas_instance=self, networkspecies=networkspecies)

                        if subsection_ind == 0:
                            # update plot characteristics
                            self.plot.set_plot_characteristics(self.summary_plots_to_make, speci=networkspeci.split('-')[-1])
                        
                            # setup plotting geometry for summary plots per networkspeci (for all subsections)
                            self.setup_plot_geometry('summary', networkspeci, networkspeci_ii)

                        # iterate through plots to make
                        for plot_type in self.summary_plots_to_make:
                            
                            #get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                            plot_options = plot_type.split('_')[1:]

                            # if a multispecies plot is wanted then only make this when networkspeci_ii == 0
                            if ('multispecies' in plot_options) & (networkspeci_ii != 0):
                                continue

                            # make plot
                            self.make_plot('summary', plot_type, plot_options, networkspeci)
                    
                    # make station specific plots?
                    if self.report_stations:               

                        # update plot characteristics
                        self.plot.set_plot_characteristics(self.station_plots_to_make, speci=networkspeci.split('-')[-1])

                        # get valid station inds for networkspeci / subsection
                        if self.temporal_colocation:
                            valid_station_inds = self.valid_station_inds_temporal_colocation[networkspeci]['observations']
                        else:
                            valid_station_inds = self.valid_station_inds[networkspeci]['observations']

                        # initialise station ind as -1
                        self.station_ind = -1

                        # iterate through stations
                        for valid_station_ind in valid_station_inds:
                            
                            # gather some information about current station
                            self.station_ind += 1
                            self.current_station_reference = self.metadata_in_memory[networkspeci]['station_reference'][valid_station_ind, :][0]
                            self.current_station_name = self.metadata_in_memory[networkspeci]['station_name'][valid_station_ind, :][0]
                            self.current_lon = round(self.metadata_in_memory[networkspeci]['longitude'][valid_station_ind, :][0], 2)
                            self.current_lat = round(self.metadata_in_memory[networkspeci]['latitude'][valid_station_ind, :][0], 2)
                            
                            # put station data in pandas dataframe
                            # put multiple species data in pandas dataframe if a multispecies plot is wanted (spatial colocation must also be active)
                            have_multispecies_option = np.any(['multispecies' in plot_type.split('_')[1:] for plot_type in self.station_plots_to_make])
                            if (have_multispecies_option) & (networkspeci_ii == 0) & (self.spatial_colocation):
                                to_pandas_dataframe(read_instance=self, canvas_instance=self, networkspecies=networkspecies, station_index=valid_station_ind)
                            else:
                                to_pandas_dataframe(read_instance=self, canvas_instance=self, networkspecies=[networkspeci], station_index=valid_station_ind)
                            
                            #  setup plotting geometry for station plots per networkspeci (for one subsection)
                            self.setup_plot_geometry('station', networkspeci, networkspeci_ii, n_stations=len(valid_station_inds))

                            # iterate through plots to make
                            for plot_type in self.station_plots_to_make:

                                # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                                plot_options = plot_type.split('_')[1:]

                                # if a multispecies plot is wanted then only make this when networkspeci_ii == 0
                                # if making a multispecies plot per specific station, spatial colocation must be also active
                                if 'multispecies' in plot_options: 
                                    if (networkspeci_ii != 0) or (not self.spatial_colocation):
                                        continue

                                # make plot
                                self.make_plot('station', plot_type, plot_options, networkspeci)

            # do formatting to axes per networkspeci
            # create colourbars
            # harmonise xy limits(not for map, heatmap or table, or when xlim and ylim defined)
            # log axes
            # annotation
            # regression line
            # trend line

            # iterate through networks and species
            for networkspeci_ii, networkspeci in enumerate(networkspecies): 

                # update plot characteristics
                self.plot.set_plot_characteristics(self.plots_to_make, speci=networkspeci.split('-')[-1])

                # iterate through plot types
                for plot_type in self.plots_to_make:

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
                    if (self.report_summary) & (self.report_stations):
                        summary_pages = [self.summary_pages[plot_type][networkspeci]]
                        # for multispecies plots of specific stations, spatial colocation needs to be on.
                        # if not, the plots will not exist so handle this
                        if ('multispecies' in plot_options) & (not self.spatial_colocation) & (plot_type not in self.station_pages):
                            station_pages = []
                        else:
                            station_pages = [self.station_pages[plot_type][networkspeci][subsection] for subsection in self.subsections]
                        paradigm_pages = summary_pages + station_pages
                    elif self.report_summary:
                        paradigm_pages = [self.summary_pages[plot_type][networkspeci]]
                    elif self.report_stations:
                        # for multispecies plots of specific stations, spatial colocation needs to be on.
                        # if not, the plots will not exist so handle this
                        if ('multispecies' in plot_options) & (not self.spatial_colocation) & (plot_type not in self.station_pages):
                            paradigm_pages = []
                        else:
                            paradigm_pages = [self.station_pages[plot_type][networkspeci][subsection] for subsection in self.subsections]

                    # iterate through paradigm pages
                    for relevant_pages in paradigm_pages:

                        # get all relevant axes for plot_type/paradigm
                        if base_plot_type in ['periodic', 'periodic-violin']:
                            ax_types = ['hour', 'month', 'dayofweek']
                        else:
                            ax_types = ['']
                        
                        relevant_axs = []
                        relevant_data_labels = []
                        for ax_type in ax_types:
                            for relevant_page in relevant_pages:
                                if base_plot_type in ['periodic', 'periodic-violin']:
                                    for ax in self.plot_dictionary[relevant_page]['axs']:
                                        relevant_axs.append(ax['handle'][ax_type])
                                        relevant_data_labels.append(ax['data_labels'])
                                else:
                                    for ax in self.plot_dictionary[relevant_page]['axs']:
                                        relevant_axs.append(ax['handle'])
                                        relevant_data_labels.append(ax['data_labels'])

                        # generate colourbars for required plots in paradigm on each relevant page
                        if 'cb' in list(self.plot_characteristics[plot_type].keys()):
                            # get all cb_axs for plot_type across relevant pages
                            cb_axs = [self.plot_dictionary[relevant_page]['cb_ax'] for relevant_page in relevant_pages]
                            generate_colourbar(self, relevant_axs, cb_axs, zstat, self.plot_characteristics[plot_type], 
                                               networkspeci.split('-')[-1])

                        # harmonise xy limits for plot paradigm
                        if base_plot_type not in ['map', 'heatmap', 'table']: 
                            if base_plot_type in ['periodic', 'periodic-violin']:
                                self.plot.harmonise_xy_lims_paradigm(relevant_axs, base_plot_type, 
                                                                     self.plot_characteristics[plot_type], plot_options, 
                                                                     ylim=[self.selected_station_data_min[networkspeci], 
                                                                           self.selected_station_data_max[networkspeci]])
                            elif base_plot_type == 'scatter':
                                self.plot.harmonise_xy_lims_paradigm(relevant_axs, base_plot_type, 
                                                                     self.plot_characteristics[plot_type], plot_options, 
                                                                     relim=True)
                            else:
                                self.plot.harmonise_xy_lims_paradigm(relevant_axs, base_plot_type,
                                                                     self.plot_characteristics[plot_type], plot_options, 
                                                                     relim=True, autoscale=True)

                        # iterate through all relevant axes for plot type in paradigm
                        for relevant_ax_ii, relevant_ax in enumerate(relevant_axs):

                            # log axes?
                            if 'logx' in plot_options:
                                self.plot.log_axes(relevant_ax, 'logx')
                            if 'logy' in plot_options:
                                self.plot.log_axes(relevant_ax, 'logy')

                            # annotation
                            if 'annotate' in plot_options:
                                self.plot.annotation(relevant_ax, networkspeci, relevant_data_labels[relevant_ax_ii], 
                                                     self.plot_characteristics[plot_type], plot_type, 
                                                     plot_options=plot_options)
                                # annotate in first axis
                                if base_plot_type in ['periodic', 'periodic-violin']:
                                    break

                            # regression line
                            if 'regression' in plot_options:
                                self.plot.linear_regression(relevant_ax, networkspeci, 
                                                            relevant_data_labels[relevant_ax_ii], 
                                                            self.plot_characteristics[plot_type], 
                                                            plot_options=plot_options)

                            # trend line
                            if 'trend' in plot_options:
                                self.plot.trend(relevant_ax, networkspeci, relevant_data_labels[relevant_ax_ii], 
                                                self.plot_characteristics[plot_type], plot_options=plot_options)
    
            # save page figures
            print('WRITING PDF')
            for figure_n in self.plot_dictionary:
                fig = self.plot_dictionary[figure_n]['fig']
                self.pdf.savefig(fig, dpi=self.dpi)
                plt.close(fig)

    def setup_plot_geometry(self, plotting_paradigm, networkspeci, networkspeci_ii, n_stations=0):
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
                plot_characteristics['page_title']['t'] = '{} (Summary)\n{}'.format(plot_characteristics['page_title']['t'], networkspeci) 
            elif plotting_paradigm == 'station':
                plot_characteristics['page_title']['t'] = '{} (Per Station)\n{}\n{}'.format(plot_characteristics['page_title']['t'], self.subsection, networkspeci) 

            # update markersize in plot characteristics (timeseries and scatter plots)
            if (base_plot_type == 'timeseries') or (base_plot_type == 'scatter'):
                self.plot.get_markersize(networkspeci, self.plot_characteristics[plot_type])

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
            elif base_plot_type in ['heatmap','table']:
                if plotting_paradigm == 'summary':
                    n_plots_per_plot_type = 1
                elif plotting_paradigm == 'station':
                    n_plots_per_plot_type = n_stations
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
                            n_plots_per_plot_type = n_stations * \
                                                    (len(self.data_labels) - 1) 
                        else:
                            n_plots_per_plot_type = n_stations * \
                                                    len(self.data_labels) 
                    else:
                        n_plots_per_plot_type = n_stations 

            # get n pages per plot type
            n_pages_per_plot_type = int(np.ceil(n_plots_per_plot_type / (
                                    plot_characteristics['figure']['ncols'] * plot_characteristics['figure']['nrows'])))
            plot_ii_per_type = 0

            # iterate through n pages per plot
            for page_n in range(n_pages_per_plot_type):
                if base_plot_type == 'map':
                    plot_characteristics['figure']['subplot_kw'] = {'projection': self.plotcrs}
                fig, axs = plt.subplots(**plot_characteristics['figure'])
                # each page is handled as 1 figure, with multiple axes
                plot_reference = ''
                self.plot_dictionary[self.current_page_n] = {'fig': fig, 'plot_type': plot_type, 'axs': []}

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
                        if base_plot_type in ['periodic','periodic-violin']:
                            gs = gridspec.GridSpecFromSubplotSpec(20, 20, subplot_spec=ax)
                            grid_dict = dict()
                            grid_dict['hour'] = fig.add_subplot(gs[:9, :])
                            grid_dict['month'] = fig.add_subplot(gs[11:, :11])
                            grid_dict['dayofweek'] = fig.add_subplot(gs[11:, 13:])
                            self.plot_dictionary[self.current_page_n]['axs'].append({'handle':grid_dict, 'data_labels':[]})
                            ax.spines['top'].set_color('none')
                            ax.spines['bottom'].set_color('none')
                            ax.spines['left'].set_color('none')
                            ax.spines['right'].set_color('none')
                            ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
                            grid_dict['hour'].set_axis_off()
                            grid_dict['month'].set_axis_off()
                            grid_dict['dayofweek'].set_axis_off()

                            # turn on all axes that will be plotted on, and add yaxis grid to each axis, and change axis label tick sizes
                            for relevant_temporal_resolution_ii, relevant_temporal_resolution in enumerate(self.relevant_temporal_resolutions):
                                # turn axis on
                                grid_dict[relevant_temporal_resolution].set_axis_on()
                                # format axis
                                self.plot.format_axis(grid_dict[relevant_temporal_resolution], base_plot_type, plot_characteristics, relevant_temporal_resolution=relevant_temporal_resolution, col_ii=col_ii, last_valid_row=last_valid_row, last_row_on_page=last_row_on_page)

                        # rest of plot types
                        else:
                            self.plot_dictionary[self.current_page_n]['axs'].append({'handle':ax, 'data_labels':[]})
                            # format axis 
                            self.plot.format_axis(ax, base_plot_type, plot_characteristics, col_ii=col_ii, last_valid_row=last_valid_row, last_row_on_page=last_row_on_page)

                    # no more plots to make on page? 
                    # then turn off unneeded axes
                    else:
                        ax.set_axis_off()

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
                    self.plot_dictionary[self.current_page_n]['cb_ax'] = fig.add_axes(plot_characteristics['cb']['position'])
                    self.plot_dictionary[self.current_page_n]['cb_ax'].set_rasterized(True)
                            
                # add current page number
                if plotting_paradigm == 'summary':
                    self.summary_pages[plot_type][networkspeci].append(self.current_page_n)
                elif plotting_paradigm == 'station':
                    self.station_pages[plot_type][networkspeci][self.subsection].append(self.current_page_n)
                
                self.current_page_n += 1

    def make_plot(self, plotting_paradigm, plot_type, plot_options, networkspeci):
        """Function that calls making of any type of plot"""

        # get zstat information from plot_type
        zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)

        # get base plot type (without stat and options)
        if zstat:
            base_plot_type = plot_type.split('-')[0] 
        else:
            base_plot_type = plot_type.split('_')[0] 

        # count how many plots are made per plot type
        current_plot_ind = 0

        # if are making bias plot, and have no valid experiment data then cannot make plot type
        if ('bias' in plot_options) & (len(self.selected_station_data[networkspeci]) < 2):
            return

        # iterate through all data arrays 
        for n_data_label, data_label in enumerate(self.data_labels):

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
                axis_ind = (current_plot_ind * len(self.subsections)) + self.subsection_ind
                relevant_page, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, plot_type, axis_ind)
                
                # set axis title
                # TODO: Fix this
                """
                if relevant_axis.get_title() == '':
                    axis_title_label = '{}\n{}'.format(data_label, self.subsection)
                    self.plot.set_axis_title(relevant_axis, axis_title_label, self.plot_characteristics[plot_type])
                """

                # make map if there are data
                if not self.selected_station_data[networkspeci]:
                    relevant_axis.set_axis_off()
                    relevant_axis.set_visible(False)
                else:
                    # calculate z statistic
                    z_statistic, active_map_valid_station_inds = calculate_z_statistic(self, z1, z2, zstat, networkspeci)
                    self.active_map_valid_station_inds = active_map_valid_station_inds

                    # make map plot
                    self.plot.make_map(relevant_axis, networkspeci, z_statistic, self.plot_characteristics[plot_type], plot_options=plot_options)
                    if z2 == '':
                        self.plot_dictionary[relevant_page]['axs'][axis_ind]['data_labels'].append(z1)
                    else:
                        self.plot_dictionary[relevant_page]['axs'][axis_ind]['data_labels'].append(z2)

            # heatmap and table
            elif base_plot_type in ['heatmap','table']:

                # first subsection?
                # then create nested dictionary to store statistical information across all subsections
                if self.subsection_ind == 0:
                    if plotting_paradigm == 'summary':
                        if zstat not in self.subsection_stats_summary:
                            self.subsection_stats_summary[zstat] = {}
                        self.subsection_stats_summary[zstat][data_label_legend] = []
                    elif plotting_paradigm == 'station':
                        if self.current_station_reference not in self.subsection_stats_station:
                            self.subsection_stats_station[self.current_station_reference] = {}
                        if zstat not in self.subsection_stats_station[self.current_station_reference]:
                            self.subsection_stats_station[self.current_station_reference][zstat] = {}
                        self.subsection_stats_station[self.current_station_reference][zstat][data_label_legend] = []

                # add stat for current data array (if has been calculated correctly, otherwise append NaNs)                                         
                if data_label in self.selected_station_data[networkspeci]:
                    if len(self.selected_station_data[networkspeci][data_label]['pandas_df']['data']) > 0:
                        data_to_add = self.selected_station_data[networkspeci][data_label]['all'][zstat][0]
                    else:
                        data_to_add = np.NaN
                else:
                    data_to_add = np.NaN

                if plotting_paradigm == 'summary':
                    self.subsection_stats_summary[zstat][data_label_legend].append(data_to_add)
                elif plotting_paradigm == 'station':
                    self.subsection_stats_station[self.current_station_reference][zstat][data_label_legend].append(data_to_add)

            # other plots (1 plot per subsection with multiple data arrays for summary paradigm, 1 plot per subsection per station for station paradigm)
            else:

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
                elif plotting_paradigm == 'station':
                    if 'individual' in plot_options:
                        if ((base_plot_type == 'scatter') or ('bias' in plot_options) or 
                            (z_statistic_sign == 'bias')):
                            axis_ind = (current_plot_ind + self.station_ind + (len(self.experiments) - 1) * self.station_ind)
                        else:
                            axis_ind = (current_plot_ind + self.station_ind + len(self.experiments) * self.station_ind)
                    else:
                        axis_ind = self.station_ind

                # get relevant axis
                relevant_page, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, plot_type, 
                                                                           axis_ind)
                
                # set axis title
                if isinstance(relevant_axis, dict):
                    for sub_ax in relevant_axis.values():
                        axis_title = sub_ax.get_title()
                        break
                else:
                    axis_title = relevant_axis.get_title()
                if axis_title == '':
                    if plotting_paradigm == 'summary':
                        axis_title_label = self.subsection
                    elif plotting_paradigm == 'station':
                        axis_title_label = '{} ({}, {})'.format(self.current_station_name, self.current_lon, self.current_lat)
                    self.plot.set_axis_title(relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

                # make plot if there is data
                if not self.selected_station_data[networkspeci]:
                    # relevant axis is a dict of the different temporal aggregations in some cases (e.g. periodic plots)
                    if isinstance(relevant_axis, dict):
                        for temporal_aggregation_resolution, temporal_aggregation_relevant_axis in relevant_axis.items():
                            temporal_aggregation_relevant_axis.set_axis_off()
                    else:
                        relevant_axis.set_axis_off()
                else:
                    # periodic plots
                    if base_plot_type in ['periodic','periodic-violin']:
                        # skip observational array if plotting bias stat
                        if (z_statistic_sign == 'bias') & (data_label == 'observations'):
                            continue
                        if data_label in self.selected_station_data[networkspeci]:
                            if len(self.selected_station_data[networkspeci][data_label]['pandas_df']['data']) > 0:
                                if base_plot_type == 'periodic':
                                    self.plot.make_periodic(relevant_axis, networkspeci, data_label, self.plot_characteristics[plot_type], zstat=zstat, plot_options=plot_options)
                                elif base_plot_type == 'periodic-violin':
                                    self.plot.make_periodic(relevant_axis, networkspeci, data_label, self.plot_characteristics[plot_type], plot_options=plot_options)
                                self.plot_dictionary[relevant_page]['axs'][axis_ind]['data_labels'].append(data_label)

                    # other plot types (except heatmap and table) 
                    else:
                        # skip observational array for bias/scatter plots
                        if data_label == 'observations' and ('bias' in plot_options or base_plot_type == 'scatter'):
                            continue
                        else:
                            func = getattr(self.plot, 'make_{}'.format(base_plot_type))
                        
                        if data_label in self.selected_station_data[networkspeci]:
                            if len(self.selected_station_data[networkspeci][data_label]['pandas_df']['data']) > 0:
                                func(relevant_axis, networkspeci, data_label, self.plot_characteristics[plot_type], plot_options=plot_options) 
                                self.plot_dictionary[relevant_page]['axs'][axis_ind]['data_labels'].append(data_label)              

            # iterate number of plots made for current type of plot 
            current_plot_ind += 1

        #last subsection?
        #then make plot heatmap/table plot
        if (base_plot_type in ['heatmap','table']) & (self.subsection_ind == (len(self.subsections) - 1)):
            
            # get relevant axis to plot on
            if plotting_paradigm == 'summary':
                axis_ind = 0
            elif plotting_paradigm == 'station':
                axis_ind = self.station_ind
            relevant_page, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, networkspeci, plot_type, axis_ind)

            #convert subsection_stats dicts to dataframe, with subsection names as indices
            if plotting_paradigm == 'summary':
                stats_df = pd.DataFrame(data=self.subsection_stats_summary[zstat],index=self.subsections)
            elif plotting_paradigm == 'station':
                if self.current_station_reference not in self.subsection_stats_station:
                    stats_df = pd.DataFrame()
                else:
                    stats_df = pd.DataFrame(data=self.subsection_stats_station[self.current_station_reference][zstat],index=self.subsections)

            # set axis title
            if relevant_axis.get_title() == '':
                if plotting_paradigm == 'summary':
                    axis_title_label = ''
                elif plotting_paradigm == 'station':
                    axis_title_label = '{} ({}, {})'.format(self.current_station_name, self.current_lon, self.current_lat)
                self.plot.set_axis_title(relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

            #turn off relevant axis if dataframe is empty or all NaN
            if (len(stats_df.index) == 0) or (stats_df.isnull().values.all()):
                relevant_axis.set_axis_off()
            else:
                #round dataframe
                stats_df = stats_df.round(self.plot_characteristics[plot_type]['round_decimal_places'])

                #make plot
                func = getattr(self.plot, 'make_{}'.format(base_plot_type))
                func(relevant_axis, stats_df, self.plot_characteristics[plot_type], plot_options=plot_options)
                #self.plot_dictionary[relevant_page]['axs'][axis_ind]['data_labels'].append() 

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
