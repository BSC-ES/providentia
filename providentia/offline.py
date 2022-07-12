import os
import sys
import json
import copy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from providentia import aux

from .read_aux import drop_nans
from .read import DataReader
from .filter import DataFilter
from .plot import Plot
from .statistics import to_pandas_dataframe
from .statistics import generate_colourbar
from .statistics import get_z_statistic_info
from .configuration import ProvConfiguration
from .init_standards import InitStandards

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
 
class ProvidentiaOffline(ProvConfiguration, InitStandards):
    """Run Providentia offline reports"""

    def __init__(self, read_type='parallel', **kwargs):
        
        print("Starting Providentia offline")

        #add passed arguments to self
        self.read_type = read_type

        #load configuration 
        ProvConfiguration.__init__(self, **kwargs)
        if 'config' in kwargs:
            self.load_conf(kwargs['config'])
        else:
            print("No configuration file found. The path to the config file must be added as an argument.")
            sys.exit(1)
        # update self from command line arguments
        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in kwargs.items()})

        # init GHOST standards
        InitStandards.__init__(self, obs_root=self.obs_root, ghost_version=self.ghost_version)
        
        # load necessary dictionaries
        self.report_plots = json.load(open(os.path.join(
            CURRENT_PATH, 'conf/report_plots.json')))
        #check for self defined plot characteristics file
        if hasattr(self, 'plot_characteristics_filename'):
            self.plot_characteristics_templates = json.load(self.plot_characteristics_filename)
        else:
            self.plot_characteristics_templates = json.load(open(os.path.join(
                CURRENT_PATH, 'conf/plot_characteristics_offline.json')))
        self.plot_characteristics = {}

        #add general plot characteristics to self
        for k, val in self.plot_characteristics_templates['general'].items():
            setattr(self, k, val)
        
        #get key config variables
        self.station_subset_names = list(self.sub_opts.keys())
        self.active_qa = aux.which_qa(self)
        self.active_flags = aux.which_flags(self)
        self.reading_nonghost = aux.check_for_ghost(self.active_network)
        self.experiments_legend = aux.get_experiments(self)
        self.relevant_temporal_resolutions = aux.get_relevant_temporal_resolutions(self.active_resolution)
                
        #if have no experiments, force temporal colocation to be False
        if len(self.experiments_legend) == 0:
            self.temporal_colocation = False
            self.defaults['temporal_colocation'] = False         
            for k in self.sub_opts.keys():
                self.sub_opts[k]['temporal_colocation'] = False

        # get all netCDF monthly files per species
        if not self.reading_nonghost:
            species_files = os.listdir('%s/%s/%s/%s/%s' % (self.obs_root, self.active_network,
                                                           self.ghost_version, self.active_resolution,
                                                           self.active_species))
        else:
            species_files = os.listdir('%s/%s/%s/%s/%s' % (self.nonghost_root, self.active_network.lower()[1:],
                                                           self.active_matrix, self.active_resolution,
                                                           self.active_species))

        # get monthly start date (YYYYMM) of all species files
        species_files_yearmonths = \
            [int(f.split('_')[-1][:6] + '01') for f in species_files if f != 'temporary']

        # initialise structure to store all obs
        self.all_observation_data = {self.active_network: {
            self.active_resolution: {self.active_matrix: {
                self.active_species: species_files_yearmonths}}}}

        # initialise struxture to store metadata
        self.metadata_types, self.metadata_menu = aux.init_metadata(self)

        # initialise DataReader class
        self.datareader = DataReader(self)
        # read data
        self.datareader.get_valid_obs_files_in_date_range(self.active_start_date, self.active_end_date)
        self.datareader.get_valid_experiment_files_in_date_range(self.active_start_date, self.active_end_date,
                                                                 self.active_resolution, self.active_network,
                                                                 self.active_species)
        self.datareader.read_setup(self.active_resolution, self.active_start_date, self.active_end_date,
                                   self.active_network, self.active_species, self.active_matrix)

        # create data dictionaries to fill
        self.datareader.reset_data_in_memory()
        self.metadata_inds_to_fill = np.arange(len(self.relevant_yearmonths))

        # read observations
        self.datareader.read_data('observations', self.active_start_date, self.active_end_date,
                                  self.active_network, self.active_resolution,
                                  self.active_species, self.active_matrix)
        # read selected experiments (iterate through)
        for exp in self.experiments_legend.keys():
            self.datareader.read_data(exp, self.active_start_date, self.active_end_date, self.active_network,
                                      self.active_resolution, self.active_species, self.active_matrix)

        # update dictionary of plotting parameters (colour and zorder etc.) for each data array in memory
        self.datareader.update_plotting_parameters()

        # initialise Plot class
        self.plot = Plot(read_instance=self, canvas_instance=self)

        #set plot characteristics
        plot_types = self.report_plots[self.report_type]
        plot_types.extend(['header'])
        self.plot.set_plot_characteristics(plot_types)

        # define dictionary to store plot figures per page
        self.plot_dictionary = {}

        #set plots that need to be made (summary and station specific)
        self.summary_plots_to_make = list(self.plot_characteristics.keys())
        self.station_plots_to_make = []
        for plot_type in self.summary_plots_to_make:
            #there can be no station specific plots for map plot type
            if plot_type[:4] != 'map-':
                self.station_plots_to_make.append(plot_type)

        #start making PDF
        self.start_pdf()

    def start_pdf(self):

        reports_path = (os.path.join(CURRENT_PATH, '../reports/'))

        if hasattr(self, 'report_filename'):
            if '/' in self.report_filename:
                filename = self.report_filename
            else:
                filename = '{}/{}.pdf'.format(reports_path,self.report_filename)             

            filename = self.report_filename + '.pdf'
        else:
            filename = "Providentia_offline_report.pdf"

        # open new PDF file
        with PdfPages(filename) as pdf:
            self.pdf = pdf

            #initialise first page number to plot
            self.current_page_n = 1

            #make header
            self.plot.make_header(self.plot_characteristics['header'])

            #setup plotting geometry for summary plots if required 
            if self.report_summary:   
                self.setup_plot_geometry('summary')
        
            #set station index to be -1 if station plots required
            if self.report_stations:   
                self.station_ind = -1

            # iterate through station_subset_names
            for station_subset_ind, station_subset in enumerate(self.station_subset_names):
                self.station_subset_ind = station_subset_ind
                self.station_subset = station_subset

                # update the conf options for this subset
                if station_subset_ind != 0:
                    for k in self.sub_opts[prv_station]:
                        vars(self).pop(k)
                vars(self).update({(k, self.parse_parameter(k, val)) for k, val in
                                   self.sub_opts[station_subset].items()})
                prv_station = station_subset
                aux.update_metadata_fields(self)
                aux.meta_from_conf(self)

                # create and update the representativity options
                self.representativity_menu = aux.representativity_fields(self, self.active_resolution)
                aux.representativity_conf(self)
                self.minimum_value, self.maximum_value = aux.which_bounds(self, self.active_species)

                print('\nFiltering Data for {} Subset'.format(station_subset))
                # filter dataset for current station_subset
                DataFilter(self)

                #iterate through all desired plots, making each one (summary or station specific plots)
                print('Making {} Subset Plots'.format(station_subset))  

                # get number of stations
                n_stations = len(self.datareader.plotting_params['observations']['valid_station_inds'])              

                # make summary plots?
                if self.report_summary:
                    # get median timeseries across data from filtered data, and place it pandas dataframe
                    to_pandas_dataframe(read_instance=self, canvas_instance=self)
                    #iterate through plots and make each one for subset group
                    for plot_type in self.summary_plots_to_make:
                        self.make_plot('summary', plot_type)
                    
                # make station specific plots?
                if self.report_stations:
                    #setup plotting geometry for station plots
                    self.setup_plot_geometry('station', n_stations=n_stations)
                    for actual_station_ind in self.datareader.plotting_params['observations']['valid_station_inds']:
                        #gather some information about current station
                        self.station_ind += 1
                        self.current_station_reference = self.datareader.metadata_in_memory['station_reference'][actual_station_ind, :][0]
                        self.current_station_name = self.datareader.metadata_in_memory['station_name'][actual_station_ind, :][0]
                        self.current_lon = round(self.datareader.metadata_in_memory['longitude'][actual_station_ind, :][0], 2)
                        self.current_lat = round(self.datareader.metadata_in_memory['latitude'][actual_station_ind, :][0], 2)
                        #put station data in pandas dataframe
                        to_pandas_dataframe(read_instance=self, canvas_instance=self, station_index=actual_station_ind)
                        #iterate through plots and make each one for subset group
                        for plot_type in self.station_plots_to_make:
                            self.make_plot('station', plot_type)

            #do formatting to axes
            #create colourbars
            #harmonise xy limits(not for map, heatmap or table, or when xlim and ylim defined)
            #log axes
            #annotation
            #regression line
            #trend line
            for plot_type in self.station_plots_to_make:

                #get zstat information from plot_type
                zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)

                #get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                plot_options = plot_type.split('_')[1:]

                #get base plot type (without stat and options)
                if zstat:
                    base_plot_type = plot_type.split('-')[0] 
                else:
                    base_plot_type = plot_type.split('_')[0] 

                #get relevant paradigm pages to harmonise axes limits for
                if (self.report_summary) & (self.report_stations):
                    paradigm_pages = [self.plot_characteristics[plot_type]['summary_pages'],
                                          self.plot_characteristics[plot_type]['station_pages']]
                elif self.report_summary:
                    paradigm_pages = [self.plot_characteristics[plot_type]['summary_pages']]
                elif self.report_stations:
                    paradigm_pages = [self.plot_characteristics[plot_type]['station_pages']]

                #iterate through paradigm pages
                for relevant_pages in paradigm_pages:

                    #get all relevant axes for plot_type/paradigm
                    if base_plot_type in ['periodic','periodic-violin']:
                        ax_types = ['hour','month','dayofweek']
                    else:
                        ax_types = ['']
                    for ax_type in ax_types:
                        relevant_axs = []
                        relevant_data_labels = []
                        for relevant_page in relevant_pages:
                            if base_plot_type in ['periodic','periodic-violin']:
                                for axs in self.plot_dictionary[relevant_page]['axs']:
                                    relevant_axs.append(axs['handle'][ax_type])
                                    relevant_data_labels.append(axs['data_labels'])
                            else:
                                relevant_axs.extend(self.plot_dictionary[relevant_page]['axs']['handle'])
                                relevant_data_labels.append(self.plot_dictionary[relevant_page]['axs']['data_labels'])

                    #generate colourbars for required plots in paradigm on each relevant page
                    if 'cb' in self.plot_characteristics[plot_type]:
                        #get all cb_axs for plot_type across relevant pages
                        cb_axs = [self.plot_dictionary[relevant_page]['cb_ax'] for relevant_page in relevant_pages]
                        generate_colourbar(self, relevant_axs, cb_axs, zstat, self.plot_characteristics[plot_type])

                    #harmonise xy limits for plot paradigm
                    if base_plot_type not in ['map','heatmap','table']: 
                        self.plot.harmonise_xy_lims_paradigm(relevant_axs, base_plot_type, self.plot_characteristics[plot_type], plot_options)

                    #iterate through all relevant axes for plot type in paradigm
                    for relevant_ax_ii, relevant_ax in enumerate(relevant_axs):

                        #log axes?
                        if 'logx' in plot_options:
                            self.plot.log_axes(relevant_ax, 'logx')
                        if 'logy' in plot_options:
                            self.plot.log_axes(relevant_ax, 'logy')

                        #annotation
                        self.plot.annotation(relevant_ax, relevant_data_labels[relevant_ax_ii], self.plot_characteristics[plot_type], plot_options=plot_options)

                        #regression line
                        self.plot.linear_regression(relevant_ax, relevant_data_labels[relevant_ax_ii], self.plot_characteristics[plot_type], plot_options=plot_options)

                        #trend line
                        self.plot.trend(relevant_ax, relevant_data_labels[relevant_ax_ii], self.plot_characteristics[plot_type], plot_options=plot_options)

            # save page figures
            print('WRITING PDF')
            for figure_n in self.plot_dictionary.keys():
                fig = self.plot_dictionary[figure_n]['fig']
                self.pdf.savefig(fig, dpi=self.dpi)
                plt.close(fig)

    def setup_plot_geometry(self, plotting_paradigm, n_stations=0):
        """setup plotting geometry for summary or station specific plots"""

        #depending on plot type set plots to make
        if plotting_paradigm == 'summary':
            plots_to_make = self.summary_plots_to_make
        elif plotting_paradigm == 'station':
            plots_to_make = self.station_plots_to_make
        
        # iterate through plot types to make
        for plot_type in plots_to_make:

            #add list to append page numbers per plotting paradigm
            if plotting_paradigm == 'summary':
                self.plot_characteristics[plot_type]['summary_pages'] = []
            elif plotting_paradigm == 'station':
                self.plot_characteristics[plot_type]['station_pages'] = []

            # get plot characteristics
            plot_characteristics = copy.deepcopy(self.plot_characteristics[plot_type])
            plot_characteristics_vars = list(plot_characteristics.keys())

            # update page title depending on plot paradigm
            if plotting_paradigm == 'summary':
                plot_characteristics['page_title']['t'] = '{} (Summary)'.format(plot_characteristics['page_title']['t']) 
            elif plotting_paradigm == 'station':
                plot_characteristics['page_title']['t'] = '{} (Per Station)'.format(plot_characteristics['page_title']['t']) 

            #get zstat information from plot_type
            zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)

            #get options defined to configure plot (e.g. bias, individual, annotate, etc.)
            plot_options = plot_type.split('_')[1:]

            #get base plot type (without stat and options)
            if zstat:
                base_plot_type = plot_type.split('-')[0] 
            else:
                base_plot_type = plot_type.split('_')[0] 

            # define number of plots per type
            n_plots_per_plot_type = False
            if base_plot_type == 'map':
                if 'obs' in plot_options:
                    n_plots_per_plot_type = len(self.station_subset_names)
                elif z_statistic_sign == 'bias':
                    n_plots_per_plot_type = len(self.station_subset_names) * \
                                            (len(list(self.datareader.data_in_memory.keys())) - 1)
                else:
                    n_plots_per_plot_type = len(self.station_subset_names) * \
                                            len(list(self.datareader.data_in_memory.keys()))
            elif base_plot_type in ['heatmap','table']:
                if plotting_paradigm == 'summary':
                    n_plots_per_plot_type = 1
                elif plotting_paradigm == 'station':
                    n_plots_per_plot_type = copy.deepcopy(n_stations)
            else:
                if plotting_paradigm == 'summary':
                    if 'individual' in plot_options:
                        if (base_plot_type == 'scatter') or ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                            n_plots_per_plot_type = len(self.station_subset_names) * \
                                                    (len(list(self.datareader.data_in_memory.keys())) - 1)
                        else:
                            n_plots_per_plot_type = len(self.station_subset_names) * \
                                                    len(list(self.datareader.data_in_memory.keys()))
                    else:
                        n_plots_per_plot_type = len(self.station_subset_names)
                elif plotting_paradigm == 'station':
                    if 'individual' in plot_options:
                        if (base_plot_type == 'scatter') or ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                            n_plots_per_plot_type = n_stations * \
                                                    (len(list(self.datareader.data_in_memory.keys())) - 1)
                        else:
                            n_plots_per_plot_type = n_stations * \
                                                    len(list(self.datareader.data_in_memory.keys()))
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
                #each page is handled as 1 figure, with multiple axes
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

                    #Force rasterized (bitmap) drawing in vector backend output.
                    ax.set_rasterized(True)

                    #keep iteratively plotting until have satisfied needed plots per type
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
                                #turn axis on
                                grid_dict[relevant_temporal_resolution].set_axis_on()
                                #format axis
                                self.plot.format_axis(grid_dict[relevant_temporal_resolution], base_plot_type, plot_characteristics, relevant_temporal_resolution=relevant_temporal_resolution, relevant_temporal_resolution_ii=relevant_temporal_resolution_ii, col_ii=col_ii, last_valid_row=last_valid_row, last_row_on_page=last_row_on_page)

                        #rest of plot types
                        else:
                            self.plot_dictionary[self.current_page_n]['axs'].append({'handle':ax, 'data_labels':[]})
                            #format axis 
                            self.plot.format_axis(ax, base_plot_type, plot_characteristics, col_ii=col_ii, last_valid_row=last_valid_row, last_row_on_page=last_row_on_page)

                    #no more plots to make on page? 
                    #then turn off unneeded axes
                    else:
                        ax.set_axis_off()
                    
                    plot_ii_per_type += 1
                    col_ii += 1

                #set figure attributes
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
                    self.plot_characteristics[plot_type]['summary_pages'].append(self.current_page_n)
                elif plotting_paradigm == 'station':
                    self.plot_characteristics[plot_type]['station_pages'].append(self.current_page_n)
                
                self.current_page_n += 1

    def make_plot(self, plotting_paradigm, plot_type):
        """Function that calls making of any type of plot"""

        #get zstat information from plot_type
        zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)

        #get options defined to configure plot (e.g. bias, individual, annotate, etc.)
        plot_options = plot_type.split('_')[1:]

        #get base plot type (without stat and options)
        if zstat:
            base_plot_type = plot_type.split('-')[0] 
        else:
            base_plot_type = plot_type.split('_')[0] 

        # count how many plots are made per plot type
        current_plot_ind = 0

        # iterate through all data arrays 
        original_data_array_labels = list(self.datareader.data_in_memory.keys())
        for n_data_label, original_data_label in enumerate(original_data_array_labels):

            #set how experiment should be referred to in heatmap/table
            if original_data_label == 'observations':
                data_label_legend = copy.deepcopy(original_data_label)
            else:    
                data_label_legend = self.experiments_legend[original_data_label]

            #get data label to plot
            data_label = aux.get_data_label(self.temporal_colocation, original_data_label)

            # do not show experiment data with 'obs' option set
            if ('obs' in plot_options) & (original_data_label != 'observations'):
                continue

            #skip observations data label when plotting bias
            if (z_statistic_sign == 'bias') & (original_data_label == 'observations'):
                continue
                
            # map plots (1 plot per data array/s (1 array if absolute plot,
            # 2 arrays if making bias plot), per subset)
            if base_plot_type == 'map':
                
                 # get necessary data labels to plot
                if z_statistic_sign == 'bias':
                    z1 = aux.get_data_label(self.temporal_colocation, 'observations')
                    z2 = copy.deepcopy(data_label)
                else:
                    z1 = copy.deepcopy(data_label)
                    z2 = ''

                # get relevant page/axis to plot on
                axis_ind = (current_plot_ind * len(self.station_subset_names)) + self.station_subset_ind
                relevant_page, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, plot_type, axis_ind)
                
                # set axis title
                if self.relevant_axis.get_title() == '':
                    axis_title_label = '{}\n{}'.format(original_data_label, self.station_subset)
                    self.plot.set_axis_title(relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

                # make map if there are data
                if not self.selected_station_data:
                    relevant_axis.set_axis_off()
                    relevant_axis.set_visible(False)
                else:
                    # calculate z statistic
                    z_statistic, active_map_valid_station_inds = calculate_z_statistic(self, z1, z2, zstat, self.temporal_colocation)
                    self.active_map_valid_station_inds = active_map_valid_station_inds

                    # make map plot
                    self.plot.make_map(relevant_axis, z_statistic, self.plot_characteristics[plot_type], plot_options=plot_options)
                    if z2 == '':
                        self.plot_dictionary[relevant_page]['axs'][axis_ind]['data_labels'].append(z1)
                    else:
                        self.plot_dictionary[relevant_page]['axs'][axis_ind]['data_labels'].append(z2)

            #heatmap and table
            elif base_plot_type in ['heatmap','table']:

                #first subset?
                #then create nested dictionary to store statistical information across all subsets 
                if self.station_subset_ind == 0:
                    if plotting_paradigm == 'summary':
                        if zstat not in self.subset_stats_summary:
                            self.subset_stats_summary[zstat] = {}
                        self.subset_stats_summary[zstat][data_label_legend] = []
                    elif plotting_paradigm == 'station':
                        if self.current_station_reference not in self.subset_stats_station:
                            self.subset_stats_station[self.current_station_reference] = {}
                        if zstat not in self.subset_stats_station[self.current_station_reference]:
                            self.subset_stats_station[self.current_station_reference][zstat] = {}
                        self.subset_stats_station[self.current_station_reference][zstat][data_label_legend] = []

                #add stat for current data array (if has been calculated correctly, otherwise append NaNs)                                         
                if data_label in list(self.selected_station_data.keys()):
                    if len(self.selected_station_data[data_label]['pandas_df']['data']) > 0:
                        data_to_add = self.selected_station_data[data_label]['all'][zstat][0]
                    else:
                        data_to_add = np.NaN
                else:
                    data_to_add = np.NaN

                if plotting_paradigm == 'summary':
                    self.subset_stats_summary[zstat][data_label_legend].append(data_to_add)
                elif plotting_paradigm == 'station':
                    self.subset_stats_station[self.current_station_reference][zstat][data_label_legend].append(data_to_add)

            # other plots (1 plot per subset with multiple data arrays for summary paradigm, 1 plot per subset per station for station paradigm)
            else:

                # get relevant axis to plot on
                if plotting_paradigm == 'summary':
                    if 'individual' in plot_options:
                        if ((base_plot_type == 'scatter') or ('bias' in plot_options) or 
                            (z_statistic_sign == 'bias')):
                                axis_ind = (current_plot_ind + self.station_subset_ind + (len(self.experiments_legend) - 1) * self.station_subset_ind)
                        else:
                            axis_ind = (current_plot_ind + self.station_subset_ind + len(self.experiments_legend) * self.station_subset_ind)
                    else:
                        axis_ind = self.station_subset_ind
                elif plotting_paradigm == 'station':
                    if 'individual' in plot_options:
                        if ((base_plot_type == 'scatter') or ('bias' in plot_options) or 
                            (z_statistic_sign == 'bias')):
                            axis_ind = (current_plot_ind + self.station_ind + (len(self.experiments_legend) - 1) * self.station_ind)
                        else:
                            axis_ind = (current_plot_ind + self.station_ind + len(self.experiments_legend) * self.station_ind)
                    else:
                        axis_ind = self.station_ind

                #get relevant axis
                relevant_page, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, plot_type, axis_ind)
                
                # set axis title
                if self.relevant_axis.get_title() == '':
                    if plotting_paradigm == 'summary':
                        axis_title_label = self.station_subset
                    elif plotting_paradigm == 'station':
                        axis_title_label = '{}, {} ({}, {})'.format(self.station_subset, self.current_station_name, self.current_lon, self.current_lat)
                    self.plot.set_axis_title(relevant_axis, axis_title_label, self.plot_characteristics[plot_type])

                # make plot if there is data
                if not self.selected_station_data:
                    # relevant axis is a dict of the different temporal aggregations in some cases (e.g. periodic plots)
                    if isinstance(relevant_axis, dict):
                        for temporal_aggregation_resolution, temporal_aggregation_relevant_axis in relevant_axis.items():
                            temporal_aggregation_relevant_axis.set_axis_off()
                    else:
                        relevant_axis.set_axis_off()
                else:
                    # Periodic plots
                    if base_plot_type in ['periodic','periodic-violin']:
                        # skip observational array if plotting bias stat
                        if (z_statistic_sign == 'bias') & (original_data_label == 'observations'):
                            continue
                        if data_label in list(self.selected_station_data.keys()):
                            if len(self.selected_station_data[data_label]['pandas_df']['data']) > 0:
                                if base_plot_type == 'periodic':
                                    self.make_periodic(relevant_axis, data_label, self.plot_characteristics[plot_type], zstat=zstat, plot_options=plot_options)
                                elif base_plot_type == 'periodic-violin':
                                    self.make_periodic(relevant_axis, data_label, self.plot_characteristics[plot_type], plot_options=plot_options)
                                self.plot_dictionary[relevant_page]['axs'][axis_ind]['data_labels'].append(data_label)

                    # Other plot types (distribution, timeseries, scatter)
                    else:
                        # skip observational array for bias/scatter plots
                        if original_data_label == 'observations' and ('bias' in plot_options or base_plot_type == 'scatter'):
                            continue
                        else:
                            func = getattr(self.plot, 'make_{}'.format(base_plot_type))
                        
                        if data_label in list(self.selected_station_data.keys()):
                            if len(self.selected_station_data[data_label]['pandas_df']['data']) > 0:
                                func(relevant_axis, data_label, self.plot_characteristics[plot_type], plot_options=plot_options) 
                                self.plot_dictionary[relevant_page]['axs'][axis_ind]['data_labels'].append(data_label)              

            # iterate number of plots made for current type of plot 
            current_plot_ind += 1

        #last subset?
        #then make plot heatmap/table plot
        if (base_plot_type in ['heatmap','table']) & (self.station_subset_ind == (len(self.station_subset_names) - 1)):
            
            # get relevant axis to plot on
            if plotting_paradigm == 'summary':
                axis_ind = 0
            elif plotting_paradigm == 'station':
                axis_ind = self.station_ind
            relevant_page, relevant_axis = self.get_relevant_page_axis(plotting_paradigm, plot_type, axis_ind)

            #convert subset_stats dicts to dataframe, with subset names as indices
            if plotting_paradigm == 'summary':
                stats_df = pd.DataFrame(data=self.subset_stats_summary[zstat],index=self.station_subset_names)
            elif plotting_paradigm == 'station':
                if self.current_station_reference not in self.subset_stats_station:
                    stats_df = pd.DataFrame()
                else:
                    stats_df = pd.DataFrame(data=self.subset_stats_station[self.current_station_reference][zstat],index=self.station_subset_names)

            # set axis title
            if self.relevant_axis.get_title() == '':
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

    def get_relevant_page_axis(self, plotting_paradigm, plot_type, axis_ind):
        """get relevant page and axis for current plot type/subset/axis index"""

        # get axes associated with plot type
        if plotting_paradigm == 'summary':
            relevant_pages = self.plot_characteristics[plot_type]['summary_pages']
        elif plotting_paradigm == 'station':
            relevant_pages = self.plot_characteristics[plot_type]['station_pages']

        all_relevant_pages = []
        relevant_axes = []         
        for relevant_page in relevant_pages:
            relevant_axes.extend(self.plot_dictionary[relevant_page]['axs'])
            all_relevant_pages.extend([relevant_page]*len(self.plot_dictionary[relevant_page]['axs']))

        return all_relevant_pages[axis_ind], relevant_axes[axis_ind]['handle']

def main(**kwargs):
    """Main function when running offine reports"""
    ProvidentiaOffline(**kwargs)
