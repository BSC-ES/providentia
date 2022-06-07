from ast import If
import os
from select import POLLOUT
import sys
import json
import copy

import numpy as np
import pandas as pd
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from providentia import aux
import scipy.stats as st
import seaborn as sns
from datetime import datetime

from .reading import drop_nans
from .calculate import Stats
from .calculate import ExpBias
from .prov_read import DataReader
from .filter import DataFilter
from .configuration import ProvConfiguration
from .init_standards import InitStandards
from .config import read_conf

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class ProvidentiaOffline(ProvConfiguration, InitStandards):
    """Run Providentia offline reports"""

    def __init__(self, read_type='parallel', **kwargs):
        
        print("Starting Providentia offline")
        # portrait/landscape page figsize
        self.portrait_figsize = (8.27, 11.69)
        self.landscape_figsize = (11.69, 8.27)
        self.dpi = 200
        ProvConfiguration.__init__(self, **kwargs)
        self.read_type = read_type

        # define projections for map plot and actual geographic coordinates
        self.datacrs = ccrs.PlateCarree()
        self.plotcrs = ccrs.Robinson()
        land_polygon_resolutions = {'low': '110m',
                                    'medium': '50m',
                                    'high': '10m'}
        self.feature = cfeature.NaturalEarthFeature('physical', 'land',
                                                    land_polygon_resolutions[
                                                    self.map_coastline_resolution],
                                                    facecolor='0.85')

        # update from config file
        if 'config' in kwargs:
            self.load_conf(kwargs['config'])
        else:
            print("No configuration file found. The path to the config file must be added as an argument.")
            sys.exit(1)

        # update from command line
        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in kwargs.items()})
        
        # init GHOST standards
        InitStandards.__init__(self, obs_root=self.obs_root,
                               ghost_version=self.ghost_version)
        cartopy.config['pre_existing_data_dir'] = self.cartopy_data_dir

        # load necessary dictionaries
        self.basic_stats_dict = json.load(open(os.path.join(
            CURRENT_PATH, 'conf/basic_stats_dict.json')))
        self.expbias_dict = json.load(open(os.path.join(
            CURRENT_PATH, 'conf/experiment_bias_stats_dict.json')))
        self.plots_per_report_type = json.load(open(os.path.join(
            CURRENT_PATH, 'conf/plots_per_report_type.json')))

        # get config parameters as variables
        for filename, section in zip(self.filenames, self.parent_section_names):
            
            self.filename = filename
            self.section = section
            self.section_opts = self.sub_opts[section]
            vars(self).update({(k, self.parse_parameter(k, val)) for k, val in self.section_opts.items()})

            #self.bounding_box = {'longitude': {'min': -12, 'max': 34}, 'latitude': {'min': 30, 'max': 46}}
            self.active_qa = aux.which_qa(self)
            self.active_flags = aux.which_flags(self)
            self.reading_nonghost = aux.check_for_ghost(self.selected_network)
            self.experiments_to_read = aux.get_experiments(self)
            #if have no experiments, force temporal colocation to be False
            if len(self.experiments_to_read) == 0:
                self.temporal_colocation = False    
                for k in self.sub_opts.keys():
                    self.sub_opts[k]['temporal_colocation'] = False

            # get all netCDF monthly files per species
            if not self.reading_nonghost:
                species_files = os.listdir('%s/%s/%s/%s/%s' % (self.obs_root, self.selected_network,
                                                            self.ghost_version, self.selected_resolution,
                                                            self.selected_species))
            else:
                species_files = os.listdir('%s/%s/%s/%s/%s' % (self.nonghost_root, self.selected_network.lower()[1:],
                                                            self.selected_matrix, self.selected_resolution,
                                                            self.selected_species))

            # get monthly start date (YYYYMM) of all species files
            species_files_yearmonths = \
                [int(f.split('_')[-1][:6] + '01') for f in species_files if f != 'temporary']

            # initialize structure to store all obs
            self.all_observation_data = {self.selected_network: {
                self.selected_resolution: {self.selected_matrix: {
                    self.selected_species: species_files_yearmonths}}}}

            self.metadata_types, self.metadata_menu = aux.init_metadata(self)

            # initialize DataReader
            self.datareader = DataReader(self)
            # read
            self.datareader.get_valid_obs_files_in_date_range(self.start_date, self.end_date)
            self.datareader.get_valid_experiment_files_in_date_range()
            self.datareader.read_setup(self.selected_resolution, self.start_date, self.end_date,
                                    self.selected_network, self.selected_species, self.selected_matrix)

            # create data dictionaries to fill
            self.datareader.reset_data_in_memory()
            self.metadata_inds_to_fill = np.arange(len(self.relevant_yearmonths))

            # read observations
            self.datareader.read_data('observations', self.start_date, self.end_date,
                                      self.selected_network, self.selected_resolution,
                                      self.selected_species, self.selected_matrix)
            # read selected experiments (iterate through)
            for exp in self.experiments_to_read:
                self.datareader.read_data(exp, self.selected_start_date, self.selected_end_date, self.selected_network,
                                          self.selected_resolution, self.selected_species, self.selected_matrix)

            # update dictionary of plotting parameters (colour and zorder etc.) for each data array
            self.datareader.update_plotting_parameters()

            #load characteristics_per_plot_type dictionary (dependent on read memory)
            sys.path.insert(1, os.path.join(CURRENT_PATH, 'conf'))
            from characteristics_per_plot_type import get_characteristics_per_plot_type, get_characteristics_per_plot_type_templates
            self.characteristics_per_plot_type = get_characteristics_per_plot_type(self)
            self.characteristics_per_plot_type_templates = get_characteristics_per_plot_type_templates(self)

            # update characteristics_per_plot_type
            plot_types_to_remove = []
            for plot_type in self.plots_per_report_type[self.report_type]:
                
                # get base zstat
                if '-' in plot_type:
                    if '_' in plot_type:
                        base_zstat = plot_type.split('-')[1].split('_')[0]
                    else:
                        base_zstat = plot_type.split('-')[1]
                else:
                    base_zstat = None

                # add new keys to make maps and periodic plots (with or without bias)
                if '-' in plot_type and 'violin' not in plot_type:

                    self.stats_dict = {**self.basic_stats_dict, **self.expbias_dict}
                    
                    if base_zstat in self.stats_dict:
                        # add new keys from templates
                        if plot_type[:4] == 'map-':
                            if '_individual' not in plot_type:
                                self.characteristics_per_plot_type[plot_type] = copy.deepcopy(self.characteristics_per_plot_type_templates['map'])
                            else:
                                # it is not possible to create individual plots for maps
                                print(f'Warning: {plot_type} could not be created.')
                                plot_types_to_remove.append(plot_type)
                                continue
                        elif 'periodic-' in plot_type:
                            self.characteristics_per_plot_type[plot_type] = copy.deepcopy(self.characteristics_per_plot_type_templates['periodic'])
                        elif 'heatmap-' in plot_type:
                            self.characteristics_per_plot_type[plot_type] = copy.deepcopy(self.characteristics_per_plot_type['heatmap'])
                        
                        # remove periodic plots from list when there is no observations data to show
                        if (('_obs' in plot_type and '_bias' in plot_type) or 
                            ('_obs' in plot_type and '_bias' not in plot_type and base_zstat in list(self.expbias_dict.keys()))):
                            print(f'Warning: {plot_type} could not be created.')
                            plot_types_to_remove.append(plot_type)
                            continue

                        # change title for new keys 
                        if '_bias' in plot_type:
                            page_title = '{} bias'.format(self.stats_dict[base_zstat]['label'])
                            self.characteristics_per_plot_type[plot_type]['page_title']['t'] = page_title
                        else:
                            page_title = self.stats_dict[base_zstat]['label']
                            self.characteristics_per_plot_type[plot_type]['page_title']['t'] = page_title

                # add new keys to make individual plots
                elif '_individual' in plot_type:
                    if any(valid_plot in plot_type for valid_plot in ('periodic-violin', 'timeseries', 'distribution', 'scatter')):
                        self.characteristics_per_plot_type[plot_type] = copy.deepcopy(self.characteristics_per_plot_type[plot_type.split('_')[0]])
                    elif 'periodic-' in plot_type and 'violin' not in plot_type:
                        continue
                    else:
                        # it is not possible to create individual plots for heatmaps and maps
                        print(f'Warning: {plot_type} could not be created.')
                        plot_types_to_remove.append(plot_type)

                # add new keys to make plots with annotations
                elif '_annotate' in plot_type:
                    if 'heatmap-' not in plot_type:
                        self.characteristics_per_plot_type[plot_type] = copy.deepcopy(self.characteristics_per_plot_type[plot_type.split('_')[0]])
                    else:
                        # it is not possible to add annotations to heatmaps
                        print(f'Warning: {plot_type} could not be created.')
                        plot_types_to_remove.append(plot_type)

                # add new keys to make plots only showing observations data
                elif '_obs' in plot_type:
                    if not any(invalid_plot in plot_type for invalid_plot in ('heatmap', 'scatter', '_bias')):
                        self.characteristics_per_plot_type[plot_type] = copy.deepcopy(self.characteristics_per_plot_type[plot_type.split('_')[0]])
                    else:
                        # it is not possible to show only observational data in heatmaps, scatter plots, bias plots and plots with exp stats
                        print(f'Warning: {plot_type} could not be created.')
                        plot_types_to_remove.append(plot_type)

                # add new keys to make plots which show experiment bias (i.e. exp - obs) (e.g. distribution_bias)
                elif '_bias' in plot_type:
                    self.characteristics_per_plot_type[plot_type] = copy.deepcopy(self.characteristics_per_plot_type[plot_type.split('_')[0]])

                #add new keys to make plots with trends
                elif '_trend' in plot_type:
                    #trend only currently possible for timeseries
                    if 'timeseries' in plot_type:
                        self.characteristics_per_plot_type[plot_type] = copy.deepcopy(self.characteristics_per_plot_type[plot_type.split('_')[0]])
                    else:
                        print(f'Warning: {plot_type} could not be created.')
                        plot_types_to_remove.append(plot_type)

            self.plots_per_report_type[self.report_type] = [plot_type for plot_type in self.plots_per_report_type[self.report_type] if plot_type not in plot_types_to_remove]

            self.exceedance_limit = aux.exceedance_lim(self.selected_species)
            self.temporal_axis_mapping_dict = aux.temp_axis_dict()

            self.start_pdf()

    def load_conf(self, fpath=None):
        """Load existing configurations from file
        for running offline Providentia."""

        if fpath is None:
            print("No configuration file found")
            sys.exit(1)

        # if DEFAULT is not present, then return
        if not os.path.isfile(fpath):
            print(("Error %s" % fpath))
            return

        self.sub_opts, self.all_sections, self.parent_section_names, self.subsection_names, self.filenames = read_conf(fpath)

    def start_pdf(self):

        # get path where reports will be saved
        reports_path = (os.path.join(CURRENT_PATH, '../reports/'))

        # open new PDF file
        with PdfPages(reports_path + self.filename + '.pdf') as pdf:
            
            self.pdf = pdf
            self.make_header()

            # define dictionary to store plot figures
            self.plot_dictionary = {}

            # iterate through subsets
            self.child_subsection_names = [subsection_name for subsection_name in self.subsection_names 
                                           if self.section == subsection_name.split('|')[0]]
            if len(self.child_subsection_names) > 0:
                self.station_subset_names = self.child_subsection_names
            else:
                self.station_subset_names = self.section

            if isinstance(self.station_subset_names, str):
                self.station_subset_names = [self.station_subset_names]

            #setup plotting framework
            self.setup_plot_framework()

            #setup plotting geometry for summary plots if required 
            if self.summary_plots:   
                self.setup_plot_geometry('summary')
        
            #set station index to be -1 if station plots required
            if self.station_plots:   
                self.station_ind = -1

            for station_subset_ind, station_subset in enumerate(self.station_subset_names):

                self.station_subset_ind = station_subset_ind
                self.station_subset = station_subset

                # update the conf options for this subset
                if len(self.child_subsection_names) > 0:
                    for section in self.all_sections:
                        for k in self.sub_opts[section]:
                            try:
                                vars(self).pop(k)
                            except:
                                continue
                    vars(self).update({(k, self.parse_parameter(k, val)) for k, val in
                                      self.sub_opts[station_subset].items()})
                aux.update_metadata_fields(self)
                aux.meta_from_conf(self)

                # create and update the representativity options
                self.representativity_menu = aux.representativity_fields(self, self.selected_resolution)
                aux.representativity_conf(self)
                self.minimum_value, self.maximum_value = aux.which_bounds(self, self.selected_species)

                print('\nFiltering Data for {} Subset'.format(station_subset))
                # filter dataset for current station_subset
                DataFilter(self)

                #iterate through all desired plots, making each one (summary or station specific plots)
                print('Making {} Subset Plots'.format(station_subset))  

                # get number of stations
                n_stations = len(self.datareader.plotting_params['observations']['valid_station_inds'])              

                # make summary plots?
                if self.summary_plots:
                    # get median timeseries across data from filtered data, and place it pandas dataframe
                    self.to_pandas_dataframe()
                    #iterate through plots and make each one for subset group
                    for plot_type in self.summary_plots_to_make:
                        self.make_plot('summary', plot_type, n_stations)
                    #generate colourbars for map plots 
                    relevant_axs = []
                    cb_axs = []
                    previous_plot_type = ''
                    for figure_n in self.plot_dictionary.keys():
                        plot_type = self.plot_dictionary[figure_n]['plot_type']
                        if (previous_plot_type != plot_type) & (len(relevant_axs) > 0):
                            self.generate_colourbar(cb_axs, relevant_axs, previous_plot_type)
                            relevant_axs = []
                            cb_axs = []
                        if 'cb' in list(self.characteristics_per_plot_type[plot_type].keys()):
                            relevant_axs.extend(self.plot_dictionary[figure_n]['axs'])
                            cb_axs.append(self.plot_dictionary[figure_n]['cb_ax'])
                        previous_plot_type = copy.deepcopy(plot_type)
                    if len(relevant_axs) > 0:
                        self.generate_colourbar(cb_axs, relevant_axs, previous_plot_type)
                    
                # make station specific plots?
                if self.station_plots:
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
                        self.to_pandas_dataframe(station_index=actual_station_ind)
                        #iterate through plots and make each one for subset group
                        for plot_type in self.station_plots_to_make:
                            self.make_plot('station', plot_type, n_stations)

            #harmonise xy limits (not for maps and heatmaps)
            for plot_type in self.station_plots_to_make:
                if 'map-' not in plot_type:
                    self.harmonise_xy_lims(plot_type)

            # save page figures
            print('WRITING PDF')
            for figure_n in self.plot_dictionary.keys():
                fig = self.plot_dictionary[figure_n]['fig']
                self.pdf.savefig(fig, dpi=self.dpi)
                plt.close(fig)

    def make_header(self):
        
        # set tile
        page = plt.figure(figsize=self.portrait_figsize)
        if hasattr(self, 'report_title'):
            txt = self.report_title
        else:
            txt = 'Providentia Offline Report'
        page.text(0.5, 0.9, txt, transform=page.transFigure,
                  size=20, ha="center", va='top', wrap=True)

        txt = 'Network = {}\nTemporal Resolution = {}\n' \
              'Species = {}\nDate Range = {} - {}\nExperiments = {}\n'\
            .format(self.selected_network,
                    self.selected_resolution,
                    self.selected_species,
                    self.start_date,
                    self.end_date, self.experiments_to_read)

        # txt += 'Station Subset/s = {}\n'.format(self.station_subset_names)
        if hasattr(self, 'bounding_box'):
            txt += 'Bounding Box = [{}E:{}E, {}N:{}N]\n'.\
                format(self.bounding_box['longitude']['min'],
                       self.bounding_box['longitude']['max'],
                       self.bounding_box['latitude']['min'],
                       self.bounding_box['latitude']['max'])

        page.text(0.5, 0.82, txt, transform=page.transFigure,
                  weight='light', size=15, ha="center", va='top', wrap=True)

        self.pdf.savefig(page, dpi=self.dpi, backend='pgf')
        plt.close(page)

    def setup_plot_framework(self):
        """setup plotting framework for PDF"""

        #set plots that need to be made (summary and station specific)
        self.summary_plots_to_make = self.plots_per_report_type[self.report_type]
        self.station_plots_to_make = []

        # initialize array with plots to be removed
        plot_types_to_remove = []
        
        for plot_type in self.summary_plots_to_make: 
            
            #for station specific plots, there can be no specific plots for map/heatmap
            if 'map-' not in plot_type:
                self.station_plots_to_make.append(plot_type)
    
            # do not create scatter plot if the temporal colocation is not active
            if ('scatter' in plot_type) and (self.temporal_colocation == False):
                plot_types_to_remove.append(plot_type)
                print(f'Warning: {plot_type} could not be created.')

        #if no experiments are defined, remove all bias statistic plots and notify user of this
        if len(list(self.datareader.data_in_memory.keys())) == 1:
            print('*** WARNING!!! NO EXPERIMENTS DEFINED, SO NO BIAS PLOTS WILL BE MADE')
            for plot_type in self.summary_plots_to_make:
                if 'heatmap-' in plot_type:
                    plot_types_to_remove.append(plot_type)
                elif '_bias' in plot_type:
                    plot_types_to_remove.append(plot_type) 
                # get zstat
                elif '-' in plot_type:
                    if ('_individual' in plot_type) or ('_annotate' in plot_type) or ('_obs' in plot_type):
                        if '_bias' not in plot_type:
                            zstat = plot_type.split('_')[0].split('-')[1]
                        else:
                            zstat = plot_type.split('_')[0].split('-')[1] + '_bias'
                    else:
                        zstat = plot_type.split('-')[1]
                    # get zstat type (basic or expbias) and sign (absolute or bias)
                    z_statistic_type, z_statistic_sign, base_zstat = get_z_statistic_type_sign(self.basic_stats_dict, zstat) 
                    #if zstat sign is bias but have 0 models, remove plot
                    if z_statistic_sign == 'bias':
                        plot_types_to_remove.append(plot_type)

        self.summary_plots_to_make = [plot_type for plot_type in self.summary_plots_to_make if plot_type not in plot_types_to_remove]
        self.station_plots_to_make = [plot_type for plot_type in self.station_plots_to_make if plot_type not in plot_types_to_remove]

        #initialise first page number to plot
        self.current_page_n = 1

    def setup_plot_geometry(self, plotting_paradigm, n_stations=0):
        """setup plotting geometry for summary or station specific plots"""

        #depending on plot type set plots to make
        if plotting_paradigm == 'summary':
            plots_to_make = self.summary_plots_to_make
        elif plotting_paradigm == 'station':
            plots_to_make = self.station_plots_to_make
        
        # iterate through plot types to make
        for plot_type in plots_to_make:

            # get plot characteristics
            plot_characteristics = copy.deepcopy(self.characteristics_per_plot_type[plot_type])
            plot_characteristics_vars = list(plot_characteristics.keys())
            
            # update page title depending on plot paradigm
            if plotting_paradigm == 'summary':
                plot_characteristics['page_title']['t'] = '{} (Summary)'.format(plot_characteristics['page_title']['t']) 
                if 'bias_title' in list(plot_characteristics.keys()):
                    plot_characteristics['bias_title']['t'] = '{} (Summary)'.format(plot_characteristics['bias_title']['t']) 
            elif plotting_paradigm == 'station':
                plot_characteristics['page_title']['t'] = '{} (Per Station)'.format(plot_characteristics['page_title']['t']) 
                if 'bias_title' in list(plot_characteristics.keys()):
                    plot_characteristics['bias_title']['t'] = '{} (Summary)'.format(plot_characteristics['bias_title']['t']) 

            # get zstat
            if '-' in plot_type:
                if ('_individual' in plot_type) or ('_annotate' in plot_type) or ('_obs' in plot_type):
                    if '_bias' not in plot_type:
                        zstat = plot_type.split('_')[0].split('-')[1]
                    else:
                        zstat = plot_type.split('_')[0].split('-')[1] + '_bias'
                else:
                    zstat = plot_type.split('-')[1]
                # get zstat type (basic or expbias) and sign (absolute or bias)
                z_statistic_type, z_statistic_sign, base_zstat = get_z_statistic_type_sign(self.basic_stats_dict, zstat)
            else:
                zstat, base_zstat = None, None

            if 'periodic-' in plot_type and 'violin' not in plot_type:
                if z_statistic_type == 'basic':
                    if base_zstat not in ['Data %', 'Exceedances']:
                        plot_characteristics['axis_ylabel']['ylabel'] = self.datareader.measurement_units 
                    else:
                        plot_characteristics['axis_ylabel']['ylabel'] = self.basic_stats_dict[base_zstat]['label']
                else:
                    plot_characteristics['axis_ylabel']['ylabel'] = self.expbias_dict[base_zstat]['label']

            # define number of plots per type
            n_plots_per_plot_type = False
            if 'heatmap-' in plot_type:
                n_plots_per_plot_type = len(plot_type.split('-')) - 1
                plot_type = 'heatmap'
            elif plot_type[:4] == 'map-':
                if '_obs' in plot_type:
                    n_plots_per_plot_type = len(self.station_subset_names)
                elif z_statistic_sign == 'bias':
                    n_plots_per_plot_type = len(self.station_subset_names) * \
                                            (len(list(self.datareader.data_in_memory.keys())) - 1)
                else:
                    n_plots_per_plot_type = len(self.station_subset_names) * \
                                            len(list(self.datareader.data_in_memory.keys()))
            else:
                if plotting_paradigm == 'summary':
                    if '_individual' in plot_type:
                        if (('scatter' in plot_type) or ('_bias' in plot_type) or 
                            ('-' in plot_type and base_zstat in list(self.expbias_dict.keys()))):
                            n_plots_per_plot_type = len(self.station_subset_names) * \
                                                    (len(list(self.datareader.data_in_memory.keys())) - 1)
                        else:
                            n_plots_per_plot_type = len(self.station_subset_names) * \
                                                    len(list(self.datareader.data_in_memory.keys()))
                    else:
                        n_plots_per_plot_type = len(self.station_subset_names)
                elif plotting_paradigm == 'station':
                    if '_individual' in plot_type:
                        if (('scatter' in plot_type) or ('_bias' in plot_type) or 
                            ('-' in plot_type and base_zstat in list(self.expbias_dict.keys()))):
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
                fig, axs = plt.subplots(**plot_characteristics['figure'])
                #each page is handled as 1 figure, with multiple axes
                plot_reference = ''
                self.plot_dictionary[self.current_page_n] = {'fig': fig, 'axs': [], 'plot_type': plot_type}

                # make page title?
                if 'page_title' in plot_characteristics_vars:
                    if '_bias' in plot_type and 'bias_title' in list(plot_characteristics.keys()):
                        st = fig.suptitle(**plot_characteristics['bias_title'])
                    else:
                        st = fig.suptitle(**plot_characteristics['page_title'])

                # iterate through axes (by row, then column)
                row_ii = -1
                col_ii = copy.deepcopy(plot_characteristics['figure']['ncols'])
               
                for ax in axs.flatten():
                    if col_ii == plot_characteristics['figure']['ncols']:
                        row_ii += 1
                        col_ii = 0
                    if row_ii == (plot_characteristics['figure']['nrows'] - 1):
                        last_row_on_page = True
                    else:
                        last_row_on_page = False

                    ax.set_rasterized(True)
                    if plot_ii_per_type < n_plots_per_plot_type:

                        # determine if are on last valid row to plot
                        if (n_plots_per_plot_type - plot_ii_per_type) <= plot_characteristics['figure']['ncols']:
                            last_valid_row = True
                        else:
                            last_valid_row = False

                        # setup periodic plot type gridspec
                        if 'periodic-' in plot_type:
                            gs = gridspec.GridSpecFromSubplotSpec(20, 20, subplot_spec=ax)
                            grid_dict = dict()
                            grid_dict['hour'] = fig.add_subplot(gs[:9, :])
                            grid_dict['month'] = fig.add_subplot(gs[11:, :11])
                            grid_dict['dayofweek'] = fig.add_subplot(gs[11:, 13:])
                            self.plot_dictionary[self.current_page_n]['axs'].append(grid_dict)
                            ax.spines['top'].set_color('none')
                            ax.spines['bottom'].set_color('none')
                            ax.spines['left'].set_color('none')
                            ax.spines['right'].set_color('none')
                            ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
                        else:
                            self.plot_dictionary[self.current_page_n]['axs'].append(ax)

                        # make axis title?
                        if 'axis_title' in plot_characteristics_vars:
                            ax.set_title(**plot_characteristics['axis_title'])

                        # make axis xlabel (only on last row on page/last valid row of visible axes)?
                        if 'axis_xlabel' in plot_characteristics_vars:
                            if last_valid_row or last_row_on_page:
                                ax.set_xlabel(**plot_characteristics['axis_xlabel'])

                        # make axis ylabel (only on leftmost column of visible axes)?
                        if ('axis_ylabel' in plot_characteristics_vars) and (col_ii == 0):
                            if ('periodic' in plot_type) and (self.selected_resolution != 'hourly'):
                                pass
                            else:
                                ax.set_ylabel(**plot_characteristics['axis_ylabel'])

                        # if are sharing xticks, and not on last row on page/last
                        # valid row, then ensure current axis xticks are hidden
                        if ('xtick_share' in plot_characteristics_vars) and (not last_valid_row) and (not last_row_on_page):
                            plt.setp(ax.get_xticklabels(), visible=False)

                        # if are sharing yticks, and not on left column, then ensure current axis yticks are hidden
                        if ('ytick_share' in plot_characteristics_vars) and (col_ii != 0):
                            plt.setp(ax.get_yticklabels(), visible=False)

                        # add gridlines?
                        if 'grid' in plot_characteristics_vars:
                            ax.grid(**plot_characteristics['grid'])

                        # add coastlines/gridlines to map?
                        if plot_type[:4] == 'map-':
                            ax.add_feature(self.feature)
                            ax.gridlines(linestyle='-', alpha=0.4)
                            if hasattr(self, 'bounding_box'):
                                if 'longitude' in self.bounding_box.keys():
                                    if 'min' in self.bounding_box['longitude']:
                                        min_lon, null = self.plotcrs.transform_point(
                                            self.bounding_box['longitude']['min'], 0,
                                            src_crs=self.datacrs)
                                        ax.set_xlim(left=min_lon)
                                    if 'max' in self.bounding_box['longitude']:
                                        max_lon, null = self.plotcrs.transform_point(
                                            self.bounding_box['longitude']['max'], 0,
                                            src_crs=self.datacrs)
                                        ax.set_xlim(right=max_lon)
                                if 'latitude' in self.bounding_box.keys():
                                    if 'min' in self.bounding_box['latitude']:
                                        null, min_lat = self.plotcrs.transform_point(0, self.bounding_box[
                                            'latitude']['min'], src_crs=self.datacrs)
                                        ax.set_ylim(bottom=min_lat)
                                    if 'max' in self.bounding_box['latitude']:
                                        null, max_lat = self.plotcrs.transform_point(0, self.bounding_box[
                                            'latitude']['max'], src_crs=self.datacrs)
                                        ax.set_ylim(top=max_lat)
                    else:
                        ax.set_visible(False)
                    
                    plot_ii_per_type += 1
                    col_ii += 1
                    
                # tight layout?
                if 'tightlayout' in plot_characteristics_vars:
                    fig.tight_layout()

                # adjust subplots?
                if 'subplots_adjust' in plot_characteristics_vars:
                    fig.subplots_adjust(**plot_characteristics['subplots_adjust'])

                # make legend?
                if 'legend' in plot_characteristics_vars:
                    leg_dict = plot_characteristics['legend']
                    leg_dict['handles'] = self.make_legend_elements()
                    fig.legend(**leg_dict)

                # add colourbar axis to plot dictionary?
                if 'cb' in plot_characteristics_vars:
                    self.plot_dictionary[self.current_page_n]['cb_ax'] = fig.add_axes(plot_characteristics['cb']['position'])
                    self.plot_dictionary[self.current_page_n]['cb_ax'].set_rasterized(True)
                            
                # add current page number
                if plotting_paradigm == 'summary':
                    self.characteristics_per_plot_type[plot_type]['summary_pages'].append(self.current_page_n)
                elif plotting_paradigm == 'station':
                    self.characteristics_per_plot_type[plot_type]['station_pages'].append(self.current_page_n)
                
                self.current_page_n += 1

    def make_legend_elements(self):
        """Function that makes legend elements"""

        # create legend elements
        # add observations element
        legend_elements = [Line2D([0], [0], marker='o', color='white',
                                  markerfacecolor=self.datareader.plotting_params['observations']['colour'],
                                  markersize=self.legend_markersize, label='observations')]
        # add element for each experiment
        for experiment_ind, experiment in enumerate(sorted(list(self.datareader.data_in_memory.keys()))):
            if experiment != 'observations':
                # add experiment element
                legend_elements.append(Line2D([0], [0], marker='o', color='white',
                                              markerfacecolor=self.datareader.plotting_params[experiment]['colour'],
                                              markersize=self.legend_markersize,
                                              label=experiment))
        return legend_elements

    def to_pandas_dataframe(self, station_index=False):
        """Function that takes data in memory puts it in a pandas dataframe.
        For summary plots this involves take the median timeseries across the timeseries.
        For station plots it is just the station in question.
        Also temporally aggregate selected data dataframes (by hour, day of week, month),
        and if have some experiment data associated with selected stations, calculate
        temporally aggregated basic statistic differences and bias statistics between
        observations and experiment data arrays.
        """

        # create new dictionary to store selection station data by data array
        self.selected_station_data = {}
       
        # iterate through data arrays in data in memory filtered dictionary
        for data_label in self.data_in_memory_filtered.keys():

            # if colocation is not active, do not convert colocated data arrays to pandas data frames
            if not self.temporal_colocation:
                if 'colocated' in data_label:
                    continue
            # else, if colocation is active, do not convert non-colocated data arrays to pandas data frames
            elif self.temporal_colocation:
                if 'colocated' not in data_label:
                    continue

            # get data for selected stations
            if station_index:
                data_array = self.data_in_memory_filtered[data_label][self.selected_species][station_index, :]
            else:
                data_array = self.data_in_memory_filtered[data_label][self.selected_species][
                             self.datareader.plotting_params[data_label]['valid_station_inds'], :]

            # if data array has no valid data for selected stations, do not create a pandas dataframe
            # data array has valid data and?
            if data_array.size and np.isnan(data_array).all() == False:

                # add nested dictionary for data array name to selection station data dictionary
                self.selected_station_data[data_label] = {}
                # if making summary plots, take cross station median of selected data for data array, and place it in a pandas
                # dataframe -->  add to selected station data dictionary
                if station_index:
                    self.selected_station_data[data_label]['pandas_df'] = pd.DataFrame(data_array, index=self.time_array, columns=['data'])
                else:
                    self.selected_station_data[data_label]['pandas_df'] = pd.DataFrame(np.nanmedian(data_array, axis=0),
                                                                                       index=self.time_array,
                                                                                       columns=['data'])
        
                # get min/max across all selected station data
                selected_station_data_all = [self.selected_station_data[data_label]['pandas_df']['data'] for data_label in self.selected_station_data.keys()]
                selected_station_data_all_flat = [item for items in selected_station_data_all for item in items]
                if len(selected_station_data_all_flat) > 0:
                    self.selected_station_data_min = np.nanmin(selected_station_data_all_flat)
                    self.selected_station_data_max = np.nanmax(selected_station_data_all_flat)
                
                # temporally aggregate selected data dataframes (by hour, day of week, month)
                self.pandas_temporal_aggregation()

                # if have some experiment data associated with selected stations, calculate
                # temporally aggregated basic statistic differences and bias statistics between
                # observations and experiment data arrays
                if len(self.experiments_to_read) > 0:
                    self.calculate_temporally_aggregated_experiment_bias_statistics()

    def pandas_temporal_aggregation(self):
        """Function that aggregates pandas dataframe data, for all data arrays,
        into desired temporal groupings also calculates all defined basic
        statistics for each individual temporal grouping
        """

        # define statistics to calculate (all basic statistics)
        statistics_to_calculate = list(self.basic_stats_dict.keys())

        # define all temporal aggregation resolutions that will be used to aggregate data
        # (variable by temporal resolution of data in memory)
        if 'hourly' in self.selected_resolution:
            self.temporal_aggregation_resolutions = ['hour', 'dayofweek', 'month', 'all']
        elif self.selected_resolution == 'daily':
            self.temporal_aggregation_resolutions = ['dayofweek', 'month', 'all']
        elif self.selected_resolution == 'monthly':
            self.temporal_aggregation_resolutions = ['month', 'all']

        # iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            # define all possible xticks for temporal resolution
            if temporal_aggregation_resolution == 'hour':
                all_xticks = np.arange(24, dtype=np.int)
            elif temporal_aggregation_resolution == 'dayofweek':
                all_xticks = np.arange(7, dtype=np.int)
            elif temporal_aggregation_resolution == 'month':
                all_xticks = np.arange(1, 13, dtype=np.int)

            # iterate through data arrays names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):

                # create nested dictionary inside selected station data dictionary for storing
                # aggregated data by data array label and temporal aggregation resolution
                self.selected_station_data[data_label][temporal_aggregation_resolution] = {}

                # if temporal group is all then simply take data in pandas df as is
                if temporal_aggregation_resolution == 'all':
                    full_grouped_data = [self.selected_station_data[data_label]['pandas_df']['data'].dropna()]
                else:
                    # else, aggregate data array into desired temporal groups (dropping NaNs)
                    grouped_data = [g['data'].dropna() for n, g in
                                    self.selected_station_data[data_label]['pandas_df'].groupby(
                                        getattr(self.selected_station_data[data_label]['pandas_df'].index,
                                                temporal_aggregation_resolution))]
                    # drop groups which have no data
                    grouped_data = [group for group in grouped_data if len(group) > 0]
                    # get xticks for groups which have valid data (i.e. which hours/days/months have valid data)
                    valid_xticks = [getattr(group.index, temporal_aggregation_resolution)[0] for group in grouped_data]
                    # create array of the size of the full range of each aggregation periods,
                    # initialised with empty lists per element (i.e. 24 for hourly aggregation)
                    full_grouped_data = [[] for _ in range(len(all_xticks))]
                    # place valid grouped data in correct positions within full array
                    full_group_indices_to_place = np.array([np.where(all_xticks == valid_xtick)[0][0]
                                                            for valid_xtick in valid_xticks], dtype=np.int)
                    for grouped_data_ii, full_group_index_to_place in enumerate(full_group_indices_to_place):
                        full_grouped_data[full_group_index_to_place] = grouped_data[grouped_data_ii]
                    # add valid xticks for group to selected data dictionary
                    # (i.e. the group xtick indexes which have valid data)
                    self.selected_station_data[data_label][
                        temporal_aggregation_resolution]['valid_xticks'] = valid_xticks

                # add full grouped data to selected data dictionary
                self.selected_station_data[data_label][
                    temporal_aggregation_resolution]['grouped_data'] = full_grouped_data

                # calculate basic statistics in each group and add them to selected station data dictionary
                for stat in statistics_to_calculate:
                    # get specific statistic dictionary (containing necessary
                    # information for calculation of selected statistic)
                    stats_dict = self.basic_stats_dict[stat]
                    # load default statistic arguments for passing to statistical function
                    function_arguments = stats_dict['arguments']
                    # if stat is exceedances then add threshold value (if available)  # TODO: what is this?
                    if stat == 'Exceedances':
                        function_arguments['threshold'] = self.exceedance_limit
                    # create empty array for storing calculated statistic by group
                    stat_output_by_group = []
                    # iterate through grouped data
                    for group in full_grouped_data:
                        # only calculate statistic if have valid data in group
                        if len(group) > 0:
                            # calculate statistic (appending to all group statistic output array)
                            stat_output_by_group = \
                                np.append(stat_output_by_group,
                                          getattr(Stats, stats_dict['function'])(group, **function_arguments))
                        # if no valid data in group, append NaN
                        else:
                            stat_output_by_group = np.append(stat_output_by_group, np.NaN)
                    # save statistical output by group to selected station data dictionary
                    self.selected_station_data[data_label][temporal_aggregation_resolution][stat] = stat_output_by_group

    def calculate_temporally_aggregated_experiment_bias_statistics(self):
        """Function that calculates temporally aggregated basic statistic
        differences and bias statistics between observations and experiment
        data arrays
        """

        # define all basic statistics that will be subtracted
        # (each experiment - observations) for each temporal aggregation
        basic_statistics = list(self.basic_stats_dict.keys())
        # define all experiment bias statistics that will be calculated
        # between each experiment and observations for each temporal aggregation
        bias_statistics = list(self.expbias_dict.keys())

        # iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:
            # iterate through data arrays names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):
                # make sure the data array is an experimental one
                if data_label.split('_')[0] != 'observations':
                    # get relevant aggregated observational statistics dictionary (i.e. colocated or not)
                    if not self.temporal_colocation:
                        relevant_aggregated_observations_dict = \
                            self.selected_station_data['observations'][temporal_aggregation_resolution]
                    else:
                        exp = data_label.split('_colocatedto_')[0]
                        relevant_aggregated_observations_dict = \
                            self.selected_station_data[
                                'observations_colocatedto_%s' % exp][temporal_aggregation_resolution]

                    # calculate temporally aggregated basic statistic differences between experiment and observations
                    # iterate through basic statistics
                    for basic_stat in basic_statistics:
                        # create empty array for storing calculated experiment-observations
                        # difference statistic by group
                        stat_diff_by_group = []
                        # iterate through aggregated index times
                        for group_ii in range(len(relevant_aggregated_observations_dict[basic_stat])):
                            # get observational and experiment group aggregated statistic
                            group_obs_stat = relevant_aggregated_observations_dict[basic_stat][group_ii]
                            group_exp_stat = self.selected_station_data[data_label][
                                temporal_aggregation_resolution][basic_stat][group_ii]

                            # take difference between observations and experiment statistics, if both values finite
                            if (np.isfinite(group_obs_stat)) & (np.isfinite(group_exp_stat)):
                                # calculate difference statistic (experiment - observations)
                                stat_diff_by_group = np.append(stat_diff_by_group, group_exp_stat-group_obs_stat)
                            # else, if one (or both) of observations/experiment statistics are NaN, append NaN
                            else:
                                stat_diff_by_group = np.append(stat_diff_by_group, np.NaN)
                        # save statistical difference output by group to selected station data dictionary
                        self.selected_station_data[data_label][temporal_aggregation_resolution][
                            '%s_bias' % (basic_stat)] = stat_diff_by_group

                    # if colocation is active, calculate temporally aggregated experiment bias
                    # statistical differences between experiment and observations
                    if self.temporal_colocation:
                        # iterate through bias statistics
                        for bias_stat in bias_statistics:
                            # get specific statistic dictionary (containing necessary information
                            # for calculation of selected statistic)
                            stats_dict = self.expbias_dict[bias_stat]
                            # load default statistic arguments for passing to statistical function
                            function_arguments = stats_dict['arguments']
                            # create empty array for storing calculated statistic by group
                            stat_output_by_group = []
                            # iterate through experimental grouped data
                            for group_ii, exp_group in enumerate(self.selected_station_data[data_label]
                                                                 [temporal_aggregation_resolution]['grouped_data']):
                                # add aggregated observational and experiment group data as
                                # arguments to pass to statistical function
                                function_arguments['obs'] = relevant_aggregated_observations_dict['grouped_data'][
                                    group_ii]
                                function_arguments['exp'] = exp_group

                                # calculate experiment bias statistic between observations and
                                # experiment group data
                                try:
                                    # calculate experiment bias statistic
                                    stat_output_by_group = \
                                        np.append(stat_output_by_group,
                                                  getattr(ExpBias, stats_dict['function'])(**function_arguments))
                                # if can not calculate statistic, append NaN
                                except Exception as e:
                                    stat_output_by_group = np.append(stat_output_by_group, np.NaN)
                            # save experiment bias statistic by group to selected station data dictionary
                            self.selected_station_data[data_label][temporal_aggregation_resolution][
                                '%s' % (bias_stat)] = stat_output_by_group

    def make_annotation(self, relevant_axis, data_label, plot_type, base_zstat):
        """Function to add annotation in plot showing the stats"""
        
        data_label_annotate = data_label.split('_colocated')[0]
        stats = self.characteristics_per_plot_type[plot_type]['annotate_stats']

        if '_annotate' in plot_type and stats:
            # reset string in box only if plots are aggregated
            if '_individual' not in plot_type and plot_type[:4] != 'map-':
                if len(list(self.datareader.data_in_memory.keys())) == 2:
                    self.str_to_annotate = []
                    self.colors = []
                else:
                    # for more than one experiment, initialize str only once, at first iteration
                    if ((data_label_annotate == list(self.datareader.data_in_memory.keys())[0]) or
                        ((('scatter' in plot_type) or ('_bias' in plot_type) or ('-' in plot_type and base_zstat in list(self.expbias_dict.keys()))) and 
                           data_label_annotate == list(self.datareader.data_in_memory.keys())[1])):
                        self.str_to_annotate = []
                        self.colors = []
            else:
                self.str_to_annotate = []
                self.colors = []

            # get stats
            stats_annotate = []
            for zstat in stats:
                if zstat in list(self.selected_station_data[data_label]['all']):
                    stats_annotate.append(zstat + ': ' + str(round(self.selected_station_data[data_label]['all'][zstat][0], 2)))

            # show number of stations if defined
            if self.characteristics_per_plot_type[plot_type]['annotate_text']['n_stations'] == True and not self.str_to_annotate:
                self.colors.append('black')
                if '_individual' in plot_type:
                    self.str_to_annotate.append('Stations: 1')
                else:
                    self.str_to_annotate.append('Stations: ' + str(len(self.datareader.plotting_params['observations']['valid_station_inds'])))

            # get colors
            self.colors.append(self.datareader.plotting_params[data_label]['colour'])

            # generate annotation
            self.str_to_annotate.append(', '.join(stats_annotate))

            # add annotation to plot
            # see loc options at https://matplotlib.org/3.1.0/api/offsetbox_api.html
            if (('_individual' in plot_type) or (plot_type[:4] == 'map-') or
                ('_individual' not in plot_type and data_label_annotate == list(self.datareader.data_in_memory.keys())[-1])):
                lines = [TextArea(line, textprops=dict(color=color, size=self.characteristics_per_plot_type[plot_type]['annotate_text']['fontsize'])) 
                                  for line, color in zip(self.str_to_annotate, self.colors)]
                bbox = AnchoredOffsetbox(child=VPacker(children=lines, align="left", pad=0, sep=1),
                                         loc=self.characteristics_per_plot_type[plot_type]['annotate_text']['loc'],
                                         bbox_transform=relevant_axis.transAxes)
                bbox.patch.set(**self.characteristics_per_plot_type[plot_type]['annotate_bbox']) 
                relevant_axis.add_artist(bbox)

        elif '_annotate' in plot_type and not stats_annotate:
            print(f'The annotation statistics have not been defined for {plot_type} in characteristics_per_plot_type.py')

    def make_plot(self, plotting_paradigm, plot_type, n_stations=0):
        """Function that calls making of any type of plot"""

        # table?
        if 'table' in plot_type:
            self.make_table()

        # heatmap plot?
        elif 'heatmap-' in plot_type:
            heatmap_types = plot_type.split('-')[1:]
            if self.station_subset_ind == 0:
                self.heatmap_dict = {}
                for heatmap_type in heatmap_types:
                    self.heatmap_dict[heatmap_type] = {}
                    for original_data_label in self.datareader.data_in_memory.keys():
                        if original_data_label != 'observations':
                            # get experiment name
                            # exp_name = get_exp_name(original_data_label)  # does not work if expid is same and the domain is different
                            exp_name = original_data_label
                            self.heatmap_dict[heatmap_type][exp_name] = []
            for heatmap_type in heatmap_types:
                for original_data_label in self.datareader.data_in_memory.keys():
                    if original_data_label != 'observations':
                        exp_name = original_data_label
                        if self.temporal_colocation:
                            data_label = '{}_colocatedto_observations'.format(original_data_label)
                        else:
                            data_label = original_data_label
                        if data_label in list(self.selected_station_data.keys()):
                            if len(self.selected_station_data[data_label]['pandas_df']['data']) > 0:
                                if heatmap_type in list(self.basic_stats_dict.keys()):
                                    self.heatmap_dict[heatmap_type][exp_name].append(
                                        self.selected_station_data[data_label]['all'][
                                            '{}_bias'.format(heatmap_type)][0])
                                else:
                                    self.heatmap_dict[heatmap_type][exp_name].append(
                                        self.selected_station_data[data_label]['all'][heatmap_type][0])
                            else:
                                self.heatmap_dict[heatmap_type][exp_name].append(np.NaN)
                        else:
                            self.heatmap_dict[heatmap_type][exp_name].append(np.NaN)

            if self.station_subset_ind == (len(self.station_subset_names) - 1):
                for heatmap_type_ii, heatmap_type in enumerate(heatmap_types):
                    relevant_axis = self.get_relevant_axis(plotting_paradigm, 'heatmap', heatmap_type_ii)
                    # make heatmap if there are data
                    if not self.selected_station_data:
                        relevant_axis.set_axis_off()
                    else:
                        axis_title_characteristics = copy.deepcopy(self.characteristics_per_plot_type['heatmap']['axis_title'])
                        axis_title_characteristics['label'] = heatmap_type
                        relevant_axis.set_title(**axis_title_characteristics)
                        if self.heatmap_dict:
                            heatmap_df = pd.DataFrame(data=self.heatmap_dict[heatmap_type],
                                                    index=self.station_subset_names)
                            self.make_heatmap(relevant_axis, heatmap_type, heatmap_df)

        # other plot type?
        else:
            # count how many plots are made per plot type
            current_plot_ind = 0

            # iterate through all data arrays
            original_data_array_labels = list(self.datareader.data_in_memory.keys())
            
            for n_data_label, original_data_label in enumerate(original_data_array_labels):

                # do not show experiment data with _obs configuration (if possible)
                if '_obs' in plot_type and original_data_label != 'observations':
                    continue

                # get zstat
                if '-' in plot_type:
                    if ('_individual' in plot_type) or ('_annotate' in plot_type) or ('_obs' in plot_type):
                        if '_bias' not in plot_type:
                            zstat = plot_type.split('_')[0].split('-')[1]
                        else:
                            zstat = plot_type.split('_')[0].split('-')[1] + '_bias'
                    else:
                        zstat = plot_type.split('-')[1]
                    # get zstat type (basic or expbias) and sign (absolute or bias)
                    z_statistic_type, z_statistic_sign, base_zstat = get_z_statistic_type_sign(self.basic_stats_dict, zstat)
                else:
                    zstat, base_zstat = None, None
                    
                # map plots (1 plot per data array/s (1 array if absolute plot,
                # 2 arrays if making bias plot), per subset)
                if plot_type[:4] == 'map-':
                    # get necessary data arrays
                    if '_obs' in plot_type:
                        if self.temporal_colocation:
                            z1 = 'observations_colocatedto_experiments'
                        else:
                            z1 = 'observations'
                        z2 = ''
                    elif z_statistic_sign == 'bias':
                        if original_data_label == 'observations':
                            continue
                        if self.temporal_colocation:
                            z1 = 'observations_colocatedto_{}'.format(original_data_label)
                            z2 = '{}_colocatedto_observations'.format(original_data_label)
                        else:
                            z1 = 'observations'
                            z2 = original_data_label
                    else:
                        if original_data_label == 'observations':
                            if self.temporal_colocation:
                                z1 = 'observations_colocatedto_experiments'
                            else:
                                z1 = 'observations'
                        else:
                            if self.temporal_colocation:
                                z1 = '{}_colocatedto_observations'.format(original_data_label)
                            else:
                                z1 = original_data_label
                        z2 = ''

                    # get data label to create annotations
                    if not z2:
                        data_label = z1
                    else:
                        data_label = z2

                    # get relevant axis to plot on
                    relevant_axis = self.get_relevant_axis(plotting_paradigm, plot_type, (current_plot_ind * len(
                                                           self.station_subset_names)) + self.station_subset_ind)
                    
                    # make map if there are data
                    if not self.selected_station_data:
                        relevant_axis.set_axis_off()
                    else:
                        # make map plot
                        n_stations = self.make_map(relevant_axis, z1, z2, zstat)
                        
                        # set axis title
                        axis_title_characteristics = copy.deepcopy(self.characteristics_per_plot_type[plot_type]['axis_title'])
                        if plot_type == 'map-Data %_obs':
                            axis_title_characteristics['label'] = '{} {}\n({} Stations)'.format(original_data_label, self.station_subset, n_stations)
                        else:
                            axis_title_characteristics['label'] = '{}\n{}'.format(original_data_label, self.station_subset)
                        relevant_axis.set_title(**axis_title_characteristics)

                        # add annotation
                        self.make_annotation(relevant_axis, data_label, plot_type, base_zstat)

                # other plots (1 plot per subset with multiple data arrays for summary paradigm, 1 plot per subset per station for station paradigm)
                else:

                    # get necessary data array to plot
                    if original_data_label == 'observations':
                        if self.temporal_colocation:
                            data_label = 'observations_colocatedto_experiments'
                        else:
                            data_label = 'observations'
                    else:
                        if self.temporal_colocation:
                            data_label = '{}_colocatedto_observations'.format(original_data_label)
                        else:
                            data_label = original_data_label
                    
                    # get relevant axis to plot on
                    if plotting_paradigm == 'summary':
                        if '_individual' in plot_type:
                            if (('scatter' in plot_type) or ('_bias' in plot_type) or 
                                ('-' in plot_type and base_zstat in list(self.expbias_dict.keys()))):
                                relevant_axis = self.get_relevant_axis(plotting_paradigm, plot_type, 
                                                                       (current_plot_ind + self.station_subset_ind + 
                                                                       (len(self.experiments_to_read) - 1) * self.station_subset_ind))
                            else:
                                relevant_axis = self.get_relevant_axis(plotting_paradigm, plot_type, 
                                                                       (current_plot_ind + self.station_subset_ind + 
                                                                       len(self.experiments_to_read) * self.station_subset_ind))
                        else:
                            relevant_axis = self.get_relevant_axis(plotting_paradigm, plot_type, self.station_subset_ind)
                    elif plotting_paradigm == 'station':
                        if '_individual' in plot_type:
                            if (('scatter' in plot_type) or ('_bias' in plot_type) or 
                                ('-' in plot_type and base_zstat in list(self.expbias_dict.keys()))):
                                relevant_axis = self.get_relevant_axis(plotting_paradigm, plot_type, 
                                                                       (current_plot_ind + self.station_ind + 
                                                                       (len(self.experiments_to_read) - 1) * self.station_ind))
                            else:
                                relevant_axis = self.get_relevant_axis(plotting_paradigm, plot_type, 
                                                                       (current_plot_ind + self.station_ind + 
                                                                       len(self.experiments_to_read) * self.station_ind))
                        else:
                            relevant_axis = self.get_relevant_axis(plotting_paradigm, plot_type, self.station_ind)
                    
                    # make plot if there are data
                    if not self.selected_station_data:
                        # relevant axis is a dict of the different temporal aggregations in some cases (e.g. periodic plots)
                        if isinstance(relevant_axis, dict):
                            for temporal_aggregation_resolution, temporal_aggregation_relevant_axis in relevant_axis.items():
                                temporal_aggregation_relevant_axis.set_axis_off()
                        else:
                            relevant_axis.set_axis_off()
                    else:
                        # Periodic plots
                        if 'periodic-' in plot_type:
                            if 'violin' in plot_type:
                                self.make_periodic(relevant_axis, data_label, plot_type, zstat = None)
                            else:
                                # skip observational array if plotting bias stat
                                if (z_statistic_sign == 'bias') & (original_data_label == 'observations'):
                                    continue
                                func = getattr(self, 'make_periodic')
                                if data_label in list(self.selected_station_data.keys()):
                                    if len(self.selected_station_data[data_label]['pandas_df']['data']) > 0:
                                        func(relevant_axis, data_label, plot_type, zstat)
                                        
                            # set axis title
                            axis_title_characteristics = copy.deepcopy(self.characteristics_per_plot_type[plot_type]['axis_title'])
                            if plotting_paradigm == 'summary':
                                axis_title_characteristics['label'] = self.station_subset
                            elif plotting_paradigm == 'station':
                                axis_title_characteristics['label'] = '{}, {} ({}, {})'.format(self.station_subset, self.current_station_name, self.current_lon, self.current_lat)
                            relevant_axis['hour'].set_title(**axis_title_characteristics, y = 1.15)

                            # add annotation
                            self.make_annotation(relevant_axis['hour'], data_label, plot_type, base_zstat)

                        # Other plot types (distribution, timeseries, scatter)
                        else:
                            # skip observational array
                            if original_data_label == 'observations' and ('_bias' in plot_type or 'scatter' in plot_type):
                                continue
                            else:
                                func = getattr(self, 'make_{}'.format(plot_type.split('_')[0]))
                                
                            if data_label in list(self.selected_station_data.keys()):
                                if len(self.selected_station_data[data_label]['pandas_df']['data']) > 0:
                                    # make standard plot, without bias
                                    if '_bias' in plot_type:
                                        func(relevant_axis, data_label, plot_type, bias=True)
                                    # make plot without bias
                                    else:
                                        func(relevant_axis, data_label, plot_type)                           

                            # set axis title
                            axis_title_characteristics = copy.deepcopy(self.characteristics_per_plot_type[plot_type]['axis_title'])
                            if plotting_paradigm == 'summary':
                                axis_title_characteristics['label'] = self.station_subset
                            elif plotting_paradigm == 'station':
                                axis_title_characteristics['label'] = '{}, {} ({}, {})'.format(self.station_subset, self.current_station_name, self.current_lon, self.current_lat)
                            relevant_axis.set_title(**axis_title_characteristics)

                            # add annotation
                            self.make_annotation(relevant_axis, data_label, plot_type, base_zstat)

                # iterate number of plots made for current type of plot 
                # skip when there is no data and we work only with experiment data to avoid empty plots
                if ((not self.selected_station_data) and 
                    (original_data_label == original_data_array_labels[-2]) and
                    (('scatter' in plot_type) or ('_bias' in plot_type) or 
                    ('-' in plot_type and base_zstat in list(self.expbias_dict.keys())))):
                    pass
                else:
                    current_plot_ind += 1

    def make_map(self, relevant_axis, z1, z2, zstat):
        """plot map"""

        # calculate z statistic
        z_statistic, active_map_valid_station_inds = self.calculate_z_statistic(z1, z2, zstat)
        
        # get colourmap
        _, _, _, z_colourmap = self.generate_colourbar_detail(zstat, 0, 0)
        
        # plot new station points on map - coloured by currently active z statisitic
        relevant_axis.scatter(self.datareader.station_longitudes[active_map_valid_station_inds],
                              self.datareader.station_latitudes[active_map_valid_station_inds],
                              s=self.unsel_station_markersize, c=z_statistic,
                              cmap=z_colourmap, zorder=3, transform=self.datacrs,
                              linewidth=0.0, alpha=None)

        return len(active_map_valid_station_inds)

    def get_relevant_axis(self, plotting_paradigm, plot_type, current_plot_ind):
        """get relevant axis for current plot type/subset/index"""

        # get axes associated with plot type
        if plotting_paradigm == 'summary':
            relevant_pages = self.characteristics_per_plot_type[plot_type]['summary_pages']
        elif plotting_paradigm == 'station':
            relevant_pages = self.characteristics_per_plot_type[plot_type]['station_pages']

        relevant_axes = [] 
        for relevant_page in relevant_pages:
            relevant_axes.extend(self.plot_dictionary[relevant_page]['axs'])

        return relevant_axes[current_plot_ind]

    def calculate_z_statistic(self, z1, z2, zstat):
        """Function that calculates selected z statistic for map"""

        # get zstat type (basic or expbias) and sign (absolute or bias)
        z_statistic_type, z_statistic_sign, base_zstat = get_z_statistic_type_sign(self.basic_stats_dict, zstat)

        # get dictionary containing necessary information for calculation of selected statistic
        if z_statistic_type == 'basic':
            stats_dict = self.basic_stats_dict[base_zstat]
        else:
            stats_dict = self.expbias_dict[base_zstat]

        # read selected data arrays (after subsetting arrays by intersection of z1/z2 valid station indices)
        # and calculate desired Z statistic (after removing NaNs from arrays)

        # get active map valid station indices (i.e. the indices of the stations data to plot on the map)
        # if only have z1, valid map indices are those simply for the z1 array
        if z2 == '':
            active_map_valid_station_inds = \
                self.datareader.plotting_params[z1]['valid_station_inds']
        else:
            # if have z2 array, get intersection of z1 and z2 valid station indices
            active_map_valid_station_inds = \
                np.intersect1d(self.datareader.plotting_params[z1]['valid_station_inds'],
                               self.datareader.plotting_params[z2]['valid_station_inds'])

        # read z1 data
        z1_array_data = \
            self.data_in_memory_filtered[z1][self.selected_species][active_map_valid_station_inds,:]
        # drop NaNs and reshape to object list of station data arrays (if not checking data %)
        if base_zstat != 'Data %':
            z1_array_data = drop_nans(z1_array_data)
        else:
            z1_array_data.tolist()

        # create empty array to store z statistic
        z_statistic = np.empty(len(z1_array_data))

        # if have no z2 data, calculate 'absolute' basic statistic
        if z2 == '':

            # load default selected z statistic arguments for passing to statistical function
            function_arguments = stats_dict['arguments']

            # iterate through stations calculating statistic
            for z_ii in range(len(z_statistic)):

                # and calculate its statistics
                z_statistic[z_ii] = \
                    getattr(Stats, stats_dict['function'])(z1_array_data[z_ii], **function_arguments)

        # else, read z2 data then calculate 'difference' statistic
        else:

            # read z2 data
            z2_array_data = \
                self.data_in_memory_filtered[z2][self.selected_species][active_map_valid_station_inds,:]
            # drop NaNs and reshape to object list of station data arrays (if not checking data %)
            if base_zstat != 'Data %':
                z2_array_data = drop_nans(z2_array_data)
            else:
                z2_array_data = z2_array_data.tolist()
            # is the difference statistic basic (i.e. mean)?
            if z_statistic_type == 'basic':

                # load default selected z statistic arguments and make separate arguments
                # dictionaries for z1/z2 calculations (as doing 2 separate calculations for z1/z2 and subtracting)
                function_arguments_z1 = stats_dict['arguments']
                function_arguments_z2 = copy.deepcopy(function_arguments_z1)

                # iterate through stations calculating statistic
                for z_ii in range(len(z_statistic)):

                    # calculate statistics for z1/z2 arrays and subtract z2-z1
                    z_statistic[z_ii] = \
                        getattr(Stats, stats_dict['function'])(z2_array_data[z_ii], **function_arguments_z2) - \
                        getattr(Stats, stats_dict['function'])(z1_array_data[z_ii], **function_arguments_z1)

            # else, is the difference statistic an experiment bias statistic (i.e. r)?
            elif z_statistic_type == 'expbias':

                # load default selected z statistic arguments for passing to statistical function
                function_arguments = stats_dict['arguments']

                # iterate through stations calculating statistic
                for z_ii in range(len(z_statistic)):

                    # set station z1/z2 arrays as arguments in argument dictionary
                    function_arguments['obs'] = z1_array_data[z_ii]
                    function_arguments['exp'] = z2_array_data[z_ii]

                    # calculate statistic
                    z_statistic[z_ii] = getattr(ExpBias, stats_dict['function'])(**function_arguments)

        # if any station z statistics come out as NaN/inf, remove respective
        # stations from active map valid station indices
        # also cut z_statistic to remove invalid NaNs/infs
        valid_z_statistic_boolean = np.isfinite(z_statistic)
        active_map_valid_station_inds = active_map_valid_station_inds[valid_z_statistic_boolean]
        z_statistic = z_statistic[valid_z_statistic_boolean]

        return z_statistic, active_map_valid_station_inds

    def generate_colourbar_detail(self, zstat, plotted_min, plotted_max):

        # get zstat type (basic or expbias) and sign (absolute or bias)
        z_statistic_type, z_statistic_sign, base_zstat = get_z_statistic_type_sign(self.basic_stats_dict, zstat)

        # get dictionary containing necessary information for calculation of selected statistic
        # also set label for colourbar
        if z_statistic_type == 'basic':
            stats_dict = self.basic_stats_dict[base_zstat]
            if base_zstat not in ['Data %','Exceedances']:
                label_units = ' ({})'.format(self.datareader.measurement_units)
            else:
                label_units = ''
        else:
            stats_dict = self.expbias_dict[base_zstat]
            label_units = ''

        # generate z colourbar label
        if z_statistic_sign == 'absolute':
            z_label = '{} {}'.format(stats_dict['label'], label_units)
        else:
            if z_statistic_type == 'basic':
                z_label = '{} bias {}'.format(stats_dict['label'], label_units)
            else:
                z_label = '{} {}'.format(stats_dict['label'], label_units)

        # set colourbar for z statistic
        # first check if have defined colourbar for z statistic, if so use that
        if 'colourbar' in list(stats_dict.keys()):
            z_colourmap = getattr(self, stats_dict['colourbar'])
        # else, set appropriate colourmap for the type of statistic
        else:
            # statistic is 'absolute'? if so use sequential colourbar
            if z_statistic_sign == 'absolute':
                z_colourmap = self.sequential_colourmap
            # the statistic is 'bias'? if so use diverging colourbar
            else:
                z_colourmap = self.diverging_colourmap

        # check if have defined vmin
        have_defined_vmin = False
        if 'vmin' in list(stats_dict.keys()):
            # if z statistic type is 'basic' and taking a difference, then DO NOT use defined vmin
            if (z_statistic_type == 'basic') & (z_statistic_sign == 'absolute'):
                have_defined_vmin = True
            elif z_statistic_type == 'expbias':
                have_defined_vmin = True
        # have defined zmin?
        if have_defined_vmin:
            z_vmin = stats_dict['vmin']
        # else, take vmin as minimum range value of calculated statistic
        else:
            z_vmin = plotted_min

        # check if have defined vmax
        have_defined_vmax = False
        if 'vmax' in list(stats_dict.keys()):
            # if z statistic type is 'basic' and taking a difference, then DO NOT use defined vmax
            if (z_statistic_type == 'basic') & (z_statistic_sign == 'absolute'):
                have_defined_vmax = True
            elif z_statistic_type == 'expbias':
                have_defined_vmax = True
        # have defined zmax?
        if have_defined_vmax:
            z_vmax = stats_dict['vmax']
        # else, take vmax as maximum range value of calculated statistic
        else:
            z_vmax = plotted_max

        # if z statistic is a bias stat, and do not have one of vmin/vmax pre-defined,
        # force vmin/vmax to be symmetrical across 0
        if (z_statistic_sign == 'bias') & (not have_defined_vmin) & (not have_defined_vmax):
            limit_stat = np.max(np.abs([z_vmin, z_vmax]))
            z_vmin = -limit_stat
            z_vmax = limit_stat

        return z_vmin, z_vmax, z_label, z_colourmap

    def generate_colourbar(self, cb_axs, axs, plot_type):

        # get axes z min/maxes
        ax_min = []
        ax_max = []
        for ax in axs:
            if len(ax.collections) > 0:
                try:
                    ax_min.append(np.nanmin(ax.collections[0].get_array()))
                    ax_max.append(np.nanmax(ax.collections[0].get_array()))
                except:
                    continue

        plotted_min = np.nanmin(ax_min)
        plotted_max = np.nanmax(ax_max)

        #get zstat 
        if '-' in plot_type:
            if ('_individual' in plot_type) or ('_annotate' in plot_type) or ('_obs' in plot_type):
                if '_bias' not in plot_type:
                    zstat = plot_type.split('_')[0].split('-')[1]
                else:
                    zstat = plot_type.split('_')[0].split('-')[1] + '_bias'
            else:
                zstat = plot_type.split('-')[1]
        else:
            zstat, base_zstat = None, None

        # get colourbar limits/label
        z_vmin, z_vmax, z_label, _ = self.generate_colourbar_detail(zstat, plotted_min, plotted_max)

        # generate colourbar tick array
        tick_array = np.linspace(z_vmin, z_vmax, 5, endpoint=True)

        # plot colourbar
        norm = matplotlib.colors.Normalize(vmin=z_vmin, vmax=z_vmax)
        for cb_ax in cb_axs:
            cb = matplotlib.colorbar.ColorbarBase(cb_ax, cmap=axs[0].collections[0].get_cmap(),
                                                  norm=norm, orientation='horizontal',ticks=tick_array)
            # plot colorbar label
            self.characteristics_per_plot_type[plot_type]['cb_xlabel']['xlabel'] = z_label
            cb.ax.set_xlabel(self.characteristics_per_plot_type[plot_type]['cb_xlabel']['xlabel'],
                             size=self.characteristics_per_plot_type[plot_type]['cb_xlabel']['fontsize'])
            # set cb ticks
            cb.ax.tick_params(labelsize=self.characteristics_per_plot_type[plot_type]['cb_xticks']['labelsize'])

        # update figure axes (to take account of new colourbar limits)
        for ax in axs:
            if len(ax.collections) > 0:
                ax.collections[0].set_clim(vmin=z_vmin,vmax=z_vmax)

    def harmonise_xy_lims(self, plot_type):

        if plot_type[:4] == 'map-':
            if hasattr(self, 'bounding_box'):
                return

        for relevant_pages in (self.characteristics_per_plot_type[plot_type]['summary_pages'],
                               self.characteristics_per_plot_type[plot_type]['station_pages']):

            # initialize arrays to save lower and upper limits in all axes
            all_xlim_lower = []
            all_xlim_upper = []
            all_ylim_lower = []
            all_ylim_upper = []

            if 'periodic-' in plot_type:
                ax_types = ['hour','month','dayofweek']
            else:
                ax_types = ['']
            for ax_type in ax_types:
                relevant_axes = []
                for relevant_page in relevant_pages:
                    if 'periodic-' in plot_type:
                        for axs in self.plot_dictionary[relevant_page]['axs']:
                            relevant_axes.append(axs[ax_type])
                    else:
                        relevant_axes.extend(self.plot_dictionary[relevant_page]['axs'])
                
                # get lower and upper limits in all axes and save in array
                for ax in relevant_axes:
                    if ax.lines:
                        xlim_lower, xlim_upper = ax.get_xlim()
                        ylim_lower, ylim_upper = ax.get_ylim()
                        all_xlim_lower.append(xlim_lower)
                        all_xlim_upper.append(xlim_upper)
                        all_ylim_lower.append(ylim_lower)
                        all_ylim_upper.append(ylim_upper)

                # get minimum and maximum from all axes and set limits
                for ax in relevant_axes:
                    if ax.lines:
                        if '_bias' in plot_type:
                            # if there is bias center plots y limits around 0
                            if np.abs(np.max(all_ylim_upper)) >= np.abs(np.min(all_ylim_lower)):
                                ylim_min = -np.abs(np.max(all_ylim_upper))
                                ylim_max = np.abs(np.max(all_ylim_upper))
                            elif np.abs(np.max(all_ylim_upper)) < np.abs(np.min(all_ylim_lower)):
                                ylim_min = -np.abs(np.min(all_ylim_lower))
                                ylim_max = np.abs(np.min(all_ylim_lower))
                            ax.set_ylim(ylim_min, ylim_max)
                        elif 'scatter' in plot_type:
                            lim_min = min(xlim_lower, ylim_lower)
                            lim_max = max(xlim_upper, ylim_upper)
                            ax.set_xlim([lim_min, lim_max])
                            ax.set_ylim([lim_min, lim_max])
                            ax.set_aspect('equal', adjustable='box')
                        else:
                            if 'periodic-' not in plot_type:
                                ax.set_xlim(np.min(all_xlim_lower), np.max(all_xlim_upper))
                            ax.set_ylim(np.min(all_ylim_lower), np.max(all_ylim_upper))

    def make_periodic(self, relevant_axis, data_label, plot_type, zstat):
        """Function that makes the temporally aggregated experiment bias statistic plots"""

        # define dictionaries defining the relevant axis, axis titles, x axis ticks
        # (and an empty nested dictionary for storing plot objects) for different temporal aggregation resolutions
        hour_aggregation_dict = {
            'ax': relevant_axis['hour'], 'title': 'H', 'xticks': np.arange(24, dtype=np.int), 'plots': {}}
        dayofweek_aggregation_dict = {
            'ax': relevant_axis['dayofweek'], 'title': 'DoW', 'xticks': np.arange(7, dtype=np.int), 'plots': {}}
        month_aggregation_dict = {
            'ax': relevant_axis['month'], 'title': 'M', 'xticks': np.arange(1, 13, dtype=np.int), 'plots': {}}

        # based on the temporal resolution of the data, combine the relevant temporal aggregation dictionaries
        aggregation_dict = {'hour': hour_aggregation_dict, 'dayofweek': dayofweek_aggregation_dict, 'month': month_aggregation_dict}
        aggregation_dict_to_remove = {}
        if self.selected_resolution == 'daily':
            aggregation_dict_to_remove = {'hour': hour_aggregation_dict}
            hour_aggregation_dict['title'] = ''
        if self.selected_resolution == 'monthly':
            aggregation_dict_to_remove = {'dayofweek': dayofweek_aggregation_dict, 'hour': hour_aggregation_dict}

        # turn on all axes that will be plotted on, and add yaxis grid to each axis, and change axis label tick sizes
        for temporal_aggregation_resolution in list(aggregation_dict.keys()):
            # add yaxis grid
            aggregation_dict[temporal_aggregation_resolution]['ax'].set_axisbelow(True)
            aggregation_dict[temporal_aggregation_resolution]['ax'].yaxis.grid(color='lightgrey', alpha=0.8, zorder=0)
            # add axis aggregation resolution label
            aggregation_dict[temporal_aggregation_resolution]['ax'].annotate(
                aggregation_dict[temporal_aggregation_resolution]['title'], (0, 1),
                xytext=(2, -2), xycoords='axes fraction', textcoords='offset points', fontsize=8, ha='left', va='top')
            # change axis tick labels
            aggregation_dict[temporal_aggregation_resolution]['ax'].tick_params(labelsize=8)

        # remove plots if there is no aggregation dict for them
        if aggregation_dict_to_remove:
            for temporal_aggregation_resolution in list(aggregation_dict_to_remove.keys()):
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_axis_off()
        
        if 'violin' in plot_type:

            # iterate through all defined temporal aggregation resolutions
            for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:
                
                if temporal_aggregation_resolution != 'all':

                    # if colocation is active, for plotting observational aggregated data,
                    # use the 'observations_colocatedto_experiments' array
                    if self.temporal_colocation:
                        if data_label.split('_')[0] == 'observations':
                            if data_label != 'observations_colocatedto_experiments':
                                continue
                        # print only relevant data, otherwise we get double lines for multiple exps
                        if data_label.split('_')[-1] not in ("observations", "experiments"):
                            continue

                    # get grouped data for current temporal aggregation resolution
                    grouped_data = self.selected_station_data[data_label][temporal_aggregation_resolution]['grouped_data']
                    # drop any groups which have no data
                    grouped_data = [group for group in grouped_data if len(group) > 0]

                    # make violin plot --> add plotted object to aggregation dictionary
                    aggregation_dict[temporal_aggregation_resolution]['plots'][data_label] = \
                        aggregation_dict[temporal_aggregation_resolution]['ax'].violinplot(grouped_data, positions=self.selected_station_data[data_label][temporal_aggregation_resolution]['valid_xticks'], points=100, widths=0.85, showmeans=False, showmedians=False, showextrema=False)

                    # set x axis limits
                    aggregation_dict[temporal_aggregation_resolution]['ax'].set_xlim(np.min(aggregation_dict[temporal_aggregation_resolution]['xticks'])-0.5, np.max(aggregation_dict[temporal_aggregation_resolution]['xticks'])+0.5)
                    # set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
                    if temporal_aggregation_resolution == 'hour':
                        aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(aggregation_dict[temporal_aggregation_resolution]['xticks'][::3])
                    else:
                        aggregation_dict[temporal_aggregation_resolution]['ax'].\
                            set_xticks(aggregation_dict[temporal_aggregation_resolution]['xticks'])
                        aggregation_dict[temporal_aggregation_resolution]['ax'].\
                            set_xticklabels([self.temporal_axis_mapping_dict[temporal_aggregation_resolution][xtick]
                                            for xtick in aggregation_dict[temporal_aggregation_resolution]['xticks']])
                    
                    # get plotted object per data array
                    violin_plot = aggregation_dict[temporal_aggregation_resolution]['plots'][data_label]

                    # get xticks (all valid aggregated time indexes) and medians to plot
                    xticks = aggregation_dict[temporal_aggregation_resolution]['xticks']
                    medians = self.selected_station_data[data_label][temporal_aggregation_resolution]['p50']
                    # split arrays if there are any temporal gaps to avoid
                    # line drawn being interpolated across missing values
                    inds_to_split = np.where(np.diff(xticks) > 1)[0]
                    if len(inds_to_split) == 0:
                        aggregation_dict[temporal_aggregation_resolution]['ax'].plot(
                            xticks, medians, marker='o',
                            color=self.datareader.plotting_params[data_label]['colour'],
                            markersize=self.temp_agg_markersize, linewidth=0.5)
                    else:
                        inds_to_split += 1
                        start_ind = 0
                        for end_ind in inds_to_split:
                            aggregation_dict[temporal_aggregation_resolution]['ax'].plot(
                                xticks[start_ind:end_ind], medians[start_ind:end_ind],
                                marker='o', color=self.datareader.plotting_params[data_label]['colour'],
                                markersize=self.temp_agg_markersize, linewidth=0.5)
                            start_ind = end_ind
                        aggregation_dict[temporal_aggregation_resolution]['ax'].plot(
                            xticks[start_ind:], medians[start_ind:], marker='o',
                            color=self.datareader.plotting_params[data_label]['colour'],
                            markersize=self.temp_agg_markersize, linewidth=0.5)

                    # update plotted objects with necessary colour and alpha
                    for patch in violin_plot['bodies']:
                        
                        patch.set_facecolor(self.datareader.plotting_params[data_label]['colour'])
                        if data_label.split('_')[0] == 'observations':
                            patch.set_alpha(0.7)
                        else:
                            patch.set_alpha(0.5)
                        # if have at least 1 experiment data array, split the violin plot across the horizontal
                        # (observations on left, experiment violin_plots on right)
                        if len(list(self.datareader.data_in_memory.keys())) > 1 and '_individual' not in plot_type:
                            m = np.mean(patch.get_paths()[0].vertices[:, 0])
                            # observations on left
                            if data_label.split('_')[0] == 'observations':
                                patch.get_paths()[0].vertices[:, 0] = np.clip(
                                    patch.get_paths()[0].vertices[:, 0], -np.inf, m)
                            # experiments on right
                            else:
                                patch.get_paths()[0].vertices[:, 0] = np.clip(
                                    patch.get_paths()[0].vertices[:, 0], m, np.inf)
         
                    # plot title (with units)
                    # if selected data resolution is 'hourly', plot the title on off the hourly aggregation axis
                    if 'hourly' in self.selected_resolution:
                        aggregation_dict['hour']['ax'].set_title('Temporal Distributions (%s)' % self.datareader.measurement_units,
                                                                                                fontsize=8.0, loc='left')
                    # otherwise, plot the units on the monthly aggregation axis
                    else:
                        aggregation_dict['month']['ax'].set_title('Temporal Distributions (%s)' % self.datareader.measurement_units,
                                                                                                fontsize=8.0, loc='left')

        else:

            # Make experiment bias plots for active bias statistic for all relevant temporal aggregation resolutions
            
            # iterate through all defined temporal aggregation resolutions
            for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:
                if temporal_aggregation_resolution != 'all':
                    #aggregation_dict[temporal_aggregation_resolution]['plots'][data_label] = \
                    aggregation_dict[temporal_aggregation_resolution]['ax'].plot(
                        aggregation_dict[temporal_aggregation_resolution]['xticks'],
                        self.selected_station_data[data_label][temporal_aggregation_resolution][zstat],
                        color=self.datareader.plotting_params[data_label]['colour'],
                        marker='o',
                        markersize=self.temp_agg_expbias_markersize, linewidth=0.5)

                    # set x axis limits
                    aggregation_dict[temporal_aggregation_resolution]['ax'].set_xlim(
                        np.min(aggregation_dict[temporal_aggregation_resolution]['xticks'])-0.5,
                        np.max(aggregation_dict[temporal_aggregation_resolution]['xticks'])+0.5)

                    # set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
                    if temporal_aggregation_resolution == 'hour':
                        aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(
                            aggregation_dict[temporal_aggregation_resolution]['xticks'][::3])
                    else:
                        aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(
                            aggregation_dict[temporal_aggregation_resolution]['xticks'])
                        aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticklabels(
                            [self.temporal_axis_mapping_dict[temporal_aggregation_resolution][xtick] for xtick
                                in aggregation_dict[temporal_aggregation_resolution]['xticks']])

                    # set xticks
                    aggregation_dict[temporal_aggregation_resolution]['ax'].xaxis.set_tick_params(labelsize=self.characteristics_per_plot_type[plot_type]['xticks']['labelsize'])
                    
                    # set yticks
                    aggregation_dict[temporal_aggregation_resolution]['ax'].yaxis.set_tick_params(labelsize=self.characteristics_per_plot_type[plot_type]['yticks']['labelsize'])
                    
                    # get zstat type (basic or expbias) and sign (absolute or bias)
                    z_statistic_type, z_statistic_sign, base_zstat = get_z_statistic_type_sign(self.basic_stats_dict, zstat)

                    if z_statistic_sign == 'bias':
                        # get value/s of minimum bias for statistic
                        if z_statistic_type == 'basic':
                            minimum_bias = self.basic_stats_dict[base_zstat]['minimum_bias']
                        else:
                            minimum_bias = self.expbias_dict[base_zstat]['minimum_bias']
                        # plot horizontal line/s across x axis at value/s of minimum experiment bias
                        for mb in minimum_bias:
                            aggregation_dict[temporal_aggregation_resolution]['ax'].axhline(y=mb, linestyle='--', linewidth=1.0,
                                                                                            color='black', zorder=0)

    def make_timeseries(self, relevant_axis, data_label, plot_type, bias=False):
        
        #bias plot?
        if bias:
            #observations and experiment must be colocated for this
            ts_obs = self.selected_station_data['observations_colocatedto_experiments']['pandas_df']
            original_data_label = data_label.split('_colocatedto')[0]
            ts_model = self.selected_station_data['{}_colocatedto_observations'.format(original_data_label)]['pandas_df'] 
            ts = ts_model - ts_obs
        else:
            ts = self.selected_station_data[data_label]['pandas_df']

        #configure size of plots if have very few points
        if (data_label.split('_')[0] == 'observations') or (bias):
            if len(ts.dropna().index) < 200:
                self.timeseries_markersize = 4.0 

        # make time series plot
        relevant_axis.plot(ts.dropna(),
                           color=self.datareader.plotting_params[data_label]['colour'],
                           marker='o', markeredgecolor=None, mew=0,
                           markersize=self.timeseries_markersize, linestyle='None')

        # plot trend line?
        if '_trend' in plot_type:
            relevant_axis.plot(ts.rolling(self.characteristics_per_plot_type['timeseries']['trend']['window'], min_periods=self.characteristics_per_plot_type['timeseries']['trend']['min_points'], center=True).mean().dropna(),
                               color=self.datareader.plotting_params[data_label]['colour'],
                               linewidth=1.0,
                               zorder=self.datareader.plotting_params[data_label]['zorder']+len(list(self.datareader.plotting_params.keys())))

        # get user-defined characteristics for xticks
        define = self.characteristics_per_plot_type['timeseries']['xticks']['define']
        n_slices = self.characteristics_per_plot_type['timeseries']['xticks']['n_slices']
        add_last_step = self.characteristics_per_plot_type['timeseries']['xticks']['last_step']

        # get steps in days or months
        label = list(self.selected_station_data.keys())[-1]
        steps = self.selected_station_data[label]['pandas_df'].dropna().index.values
        last_step = self.selected_station_data[label]['pandas_df'].dropna().index.values[-1]

        # get start and end dates
        timeseries_start_date = pd.to_datetime(steps[0])
        timeseries_end_date = pd.to_datetime(steps[-1])

        # get number of months and days
        n_months = (timeseries_end_date.year - timeseries_start_date.year) * 12 + (timeseries_end_date.month - timeseries_start_date.month)
        n_days = (timeseries_end_date - timeseries_start_date).days

        # get months that are complete
        months = pd.date_range(timeseries_start_date, timeseries_end_date, freq='MS')
        if months.size > 0 and (months[-1] != steps[-1]):
            months = months[:-1]
            n_months -= 1

        # define time slices
        if define:
            if n_months >= 3:
                steps = months
            elif n_days < 7:
                relevant_axis.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M'))
            slices = int(np.ceil(len(steps) / n_slices))

            # use default axes if the number of timesteps is lower than the number of slices
            if slices >= 1:
                xticks = steps[0::slices]
            else:
                xticks = relevant_axis.xaxis.get_ticks()
            
            # add last step to xticks
            if add_last_step and (xticks[-1] != last_step):
                xticks = np.append(xticks, last_step)
            
            # set xticks
            relevant_axis.xaxis.set_ticks(xticks) 

        # set xticks size and rotation
        relevant_axis.xaxis.set_tick_params(labelsize=self.characteristics_per_plot_type['timeseries']['xticks']['labelsize'], 
                                            rotation=self.characteristics_per_plot_type['timeseries']['xticks']['rotation'])
        # set yticks size
        relevant_axis.yaxis.set_tick_params(labelsize=self.characteristics_per_plot_type['timeseries']['yticks']['labelsize'])


    def make_distribution(self, relevant_axis, data_label, plot_type, bias=False):
        
        # make distribution plot
        minmax_diff = self.selected_station_data_max - self.selected_station_data_min
        if pd.isnull(self.parameter_dictionary[self.selected_species]['minimum_resolution']):
            n_samples = 2000
        else:
            n_samples = int(np.around(minmax_diff/(self.parameter_dictionary[self.selected_species]['minimum_resolution']/4.0),0))
            if n_samples < 2000:
                n_samples = 2000
        x_grid = np.linspace(self.selected_station_data_min, self.selected_station_data_max, n_samples, endpoint=True)

        # setup bias plot
        if bias:
            if self.temporal_colocation:
                PDF_obs = st.gaussian_kde(self.selected_station_data['observations_colocatedto_experiments']['pandas_df']['data'].dropna())
            else:
                PDF_obs = st.gaussian_kde(self.selected_station_data['observations']['pandas_df']['data'].dropna())
            PDF_model = st.gaussian_kde(self.selected_station_data[data_label]['pandas_df']['data'].dropna())
            PDF_sampled = PDF_model(x_grid) - PDF_obs(x_grid)
            # plot horizontal line across x axis at 0
            relevant_axis.axhline(y=0.0, linestyle='--', linewidth=1.0, color='black', zorder=0)

        # setup standard plot
        else:
            PDF = st.gaussian_kde(self.selected_station_data[data_label]['pandas_df']['data'].dropna())
            PDF_sampled = PDF(x_grid)

        # make plot
        relevant_axis.plot(x_grid, PDF_sampled, linewidth=1, color=self.datareader.plotting_params[data_label]['colour'])

        # make xaxis logged for certain species
        if self.selected_species in ['sconcno', 'sconcno2']:
            relevant_axis.set_xlim(0.01, 10)

        # set xticks
        relevant_axis.xaxis.set_tick_params(labelsize=self.characteristics_per_plot_type['distribution']['xticks']['labelsize'])
        # set yticks
        relevant_axis.yaxis.set_tick_params(labelsize=self.characteristics_per_plot_type['distribution']['yticks']['labelsize'])

    def make_heatmap(self, relevant_axis, base_zstat, heatmap_df):
        
        # get zstat type (basic or expbias)
        z_statistic_type = get_z_statistic_type(self.basic_stats_dict, base_zstat)

        # if z_statistic_type == 'basic', affix '_bias' to base_zstat
        if z_statistic_type == 'basic':
            zstat = '{}_bias'.format(base_zstat)
        else:
            zstat = copy.deepcopy(base_zstat)

        # get colourbar limits/colourmap
        z_vmin, z_vmax, _, z_colourmap = self.generate_colourbar_detail(zstat, heatmap_df.min().min(), heatmap_df.max().max())

        # get colorbar label
        if z_statistic_type == 'basic':
            if base_zstat not in ['Data %','Exceedances']:
                self.characteristics_per_plot_type['heatmap']['cb_xlabel']['xlabel'] = self.datareader.measurement_units
            else:
                self.characteristics_per_plot_type['heatmap']['cb_xlabel']['xlabel'] = self.basic_stats_dict[base_zstat]['label']
        else:
            self.characteristics_per_plot_type['heatmap']['cb_xlabel']['xlabel'] = self.expbias_dict[base_zstat]['label']

        # plot heatmap
        ax = sns.heatmap(heatmap_df, vmin=z_vmin, vmax=z_vmax, cmap=z_colourmap, xticklabels=1, yticklabels=1, 
                         square=True, annot=self.characteristics_per_plot_type['heatmap']['annot'], 
                         cbar_kws = {'use_gridspec':False,'orientation':'horizontal'}, ax=relevant_axis)

        # axis cuts off due to bug in matplotlib 3.1.1 - hack fix. Remove in Future!
        bottom, top = relevant_axis.get_ylim()
        relevant_axis.set_ylim(bottom + 0.5, top - 0.5)

        # set xticks (rotating ticks)
        relevant_axis.xaxis.set_tick_params(labelsize=self.characteristics_per_plot_type['heatmap']['xticks']['labelsize'], 
                                            rotation=self.characteristics_per_plot_type['heatmap']['xticks']['rotation'])
        # set yticks (rotating ticks)
        relevant_axis.yaxis.set_tick_params(labelsize=self.characteristics_per_plot_type['heatmap']['yticks']['labelsize'], 
                                            rotation=self.characteristics_per_plot_type['heatmap']['yticks']['rotation'])
        # set colorbar label
        cb = ax.collections[0].colorbar
        cb.ax.set_xlabel(self.characteristics_per_plot_type['heatmap']['cb_xlabel']['xlabel'], 
                         size=self.characteristics_per_plot_type['heatmap']['cb_xlabel']['fontsize'])
        # set cb ticks
        cb.ax.tick_params(labelsize=self.characteristics_per_plot_type['heatmap']['cb_xticks']['labelsize'])

    def make_scatter(self, relevant_axis, data_label, plot_type):

        # Get observations data
        observations_data = self.selected_station_data['observations_colocatedto_experiments']['pandas_df'].dropna()['data']

        # Get experiment data
        experiment_data = self.selected_station_data[data_label]['pandas_df'].dropna()['data']

        # Create scatter plot
        relevant_axis.plot(observations_data, experiment_data, linestyle = 'None',
                           marker=self.characteristics_per_plot_type['scatter']['markers']['character'],
                           markersize=self.characteristics_per_plot_type['scatter']['markers']['size'],
                           color=self.datareader.plotting_params[data_label]['colour'])
        
        # Set ticks
        relevant_axis.xaxis.set_tick_params(labelsize=self.characteristics_per_plot_type['scatter']['xticks']['labelsize'])
        relevant_axis.yaxis.set_tick_params(labelsize=self.characteristics_per_plot_type['scatter']['yticks']['labelsize'])

        # Set labels
        relevant_axis.set_xlabel(xlabel=self.characteristics_per_plot_type['scatter']['axis_xlabel']['xlabel'], 
                                 fontsize=self.characteristics_per_plot_type['scatter']['axis_xlabel']['fontsize'])
        relevant_axis.set_ylabel(ylabel=self.characteristics_per_plot_type['scatter']['axis_ylabel']['ylabel'], 
                                 fontsize=self.characteristics_per_plot_type['scatter']['axis_ylabel']['fontsize'])

        # Add line 1:1
        line = Line2D([0, 1], [0, 1], color='lightgrey', linewidth=1, linestyle='--')
        line.set_transform(relevant_axis.transAxes)
        relevant_axis.add_line(line)

    def make_table(self):
        pass

def get_z_statistic_type(stats_dict, zstat):
    """Function that checks if the z statistic is basic or expbias statistic"""
    # check if the chosen statistic is a basic statistic
    if zstat in stats_dict.keys():
        return 'basic'
    # if not a basic statistic, it must be an experiment bias statistic
    else:
        return 'expbias'

def get_z_statistic_sign(zstat, zstat_type):
    """Function that checks if the z statistic is an absolute or bias statistic"""
    # statistic is bias?
    if ('_bias' in zstat) or (zstat_type == 'expbias'):
        return 'bias'
    # statistic is bias?
    else:
        return 'absolute'

def get_z_statistic_type_sign(stats_dict, zstat):
    """Function that gets z statistic type and sign"""                    
    # get base name name of zstat (dropping _bias suffix)
    base_zstat = zstat.split('_bias')[0]
    # get zstat type (basic or expbias)
    z_statistic_type = get_z_statistic_type(stats_dict, base_zstat)
    # get zstat sign (absolute or bias)
    z_statistic_sign = get_z_statistic_sign(zstat, z_statistic_type)
    return z_statistic_type, z_statistic_sign, base_zstat
    
def get_exp_name(fullname):
    """ Return experiments id/name from its GHOST identifier in /gpfs.
    For example, a3p1-regional-000 returns a3p1. In the case of cams,
    e.g. cams61_emep_ph3-eu-000 it returns emep."""
    short = fullname[:fullname.find('-')]
    if 'cams' in short:
        short = fullname.split('_')[1].upper()
    return short

def main_offline(**kwargs):
    """Main function when running offine reports"""
    ProvidentiaOffline(**kwargs)

