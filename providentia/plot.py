import copy
from datetime import datetime
import json
import os

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib 
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st
import seaborn as sns

from .statistics import get_z_statistic_info
from providentia import aux

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))

class Plot:
    """
    "Class that makes plots and handles plot configuration options when defined"
    """

    def __init__(self, read_instance=None, canvas_instance=None):
        self.read_instance = read_instance
        self.canvas_instance = canvas_instance

        # set cartopy data directory
        cartopy.config['pre_existing_data_dir'] = self.read_instance.cartopy_data_dir

        # define projections for map plot and actual geographic coordinates
        self.canvas_instance.datacrs = ccrs.PlateCarree()
        proj_class = getattr(ccrs, self.canvas_instance.plot_characteristics_templates['map']['projection']) 
        self.canvas_instance.plotcrs = proj_class()
        #setup land polygons on map
        self.canvas_instance.feature = cfeature.NaturalEarthFeature('physical', 'land',
            aux.get_land_polygon_resolution(self.canvas_instance.plot_characteristics_templates['map']['map_coastline_resolution']), 
              **self.canvas_instance.plot_characteristics_templates['map']['land_polygon'])

        #set miscellaneous vars
        self.canvas_instance.temporal_axis_mapping_dict = aux.temp_axis_dict()
        self.canvas_instance.periodic_xticks = aux.periodic_xticks()
        self.canvas_instance.periodic_labels = aux.periodic_labels()

    def set_plot_characteristics(self, plot_types):
        """Iterate through all plots to make, and determine if they can and cannot be made.
           Update plot characteristics associated with specific plot types due to plot options. 

            :param plot_types: plot types to create 
            :type plot_types: list  
        """

        # add all valid defined plots to self.plot_characteristics
        for plot_type in plot_types:

            #get options defined to configure plot (e.g. bias, individual, annotate, etc.)
            plot_options = plot_type.split('_')[1:]

            #get zstat information from plot_type (if available)
            zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(plot_type=plot_type)

            # remove plots where setting 'obs' and 'bias' options together
            if ('obs' in plot_options) & ('bias' in plot_options): 
                print(f"Warning: {plot_type} could not be created as 'obs' and 'bias' options set together")
                continue

            #if no experiments are defined, remove all bias plots 
            if ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                if len(list(self.read_instance.datareader.data_in_memory.keys())) == 1:
                    print(f'Warning: No experiments defined, so {plot_type} bias plot cannot be created')
                    continue

            # add new keys to make plots with stats (map, periodic, heatmap, table)
            if zstat:

                #get base plot type (without stat and options)
                base_plot_type = plot_type.split('-')[0] 
                #combine basic and expbias stats dicts together
                stats_dict = {**basic_stats, **expbias_stats}
                
                #check all defined plot options are allowed for current plot type, if not, remove plot from list to be created
                if not all(plot_option in self.canvas_instance.plot_characteristics_templates[base_plot_type]['plot_options'] for plot_option in plot_options):
                    print(f'Warning: {plot_type} could not be created as some plot options are not valid')
                    continue
                #check desired statistic is defined in stats dict
                if base_zstat not in stats_dict:
                    print(f"Warning: {plot_type} could not be created as {base_zstat} not defined in Providentia's statistical library.")
                    continue
                # remove plots where setting 'obs', but z_statistic_sign is 'bias'
                elif ('obs' in plot_options) & (z_statistic_sign == 'bias'):
                    print(f"Warning: {plot_type} could not be created as are plotting a bias statistic but 'obs' option is set")
                    continue
                #otherwise add information for plot type from base plot type template
                else:
                    self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[base_plot_type])

                    # set page title
                    if 'bias' in plot_options:
                        self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'] = '{} {} bias'.format(self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'],stats_dict[base_zstat]['label'])
                    else:
                        self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'] = '{} {}'.format(self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'], stats_dict[base_zstat]['label'])

                    #set ylabel for periodic plots
                    if base_plot_type == 'periodic':
                        if z_statistic_type == 'basic':
                            if base_zstat not in ['Data %', 'Exceedances']:
                                self.canvas_instance.plot_characteristics[plot_type]['ylabel']['ylabel'] = self.read_instance.datareader.measurement_units 
                            else:
                                self.canvas_instance.plot_characteristics[plot_type]['ylabel']['ylabel'] = basic_stats[base_zstat]['label']
                        else:
                            self.canvas_instance.plot_characteristics[plot_type]['ylabel']['ylabel'] = expbias_stats[base_zstat]['label']

                    # get colorbar label for heatmap
                    if 'heatmap' == base_plot_type:
                        if z_statistic_type == 'basic':
                            if base_zstat not in ['Data %', 'Exceedances']:
                                self.canvas_instance.plot_characteristics[plot_type]['cb_xlabel']['xlabel'] = self.read_instance.datareader.measurement_units
                            else:
                                self.canvas_instance.plot_characteristics[plot_type]['cb_xlabel']['xlabel'] = basic_stats[base_zstat]['label']
                        else:
                            self.canvas_instance.plot_characteristics[plot_type]['cb_xlabel']['xlabel'] = expbias_stats[base_zstat]['label']

                    #define dictionary to store stats from all subsets for heatmap and table plots
                    if base_plot_type in ['heatmap','table']:
                        self.canvas_instance.subset_stats_summary = {}
                        self.canvas_instance.subset_stats_station = {}

            #add new keys for plots without stats
            else:

                #get base plot type (without options)
                base_plot_type = plot_type.split('_')[0] 

                #check all defined plot options are allowed for current plot type, if not, remove plot from list to be created
                if not all(plot_option in self.canvas_instance.plot_characteristics_templates[base_plot_type]['plot_options'] for plot_option in plot_options):
                    print(f'Warning: {plot_type} could not be created as some plot options are not valid')
                    continue
                # do not create scatter plot if the temporal colocation is not active
                elif ('scatter' == base_plot_type) & (self.canvas_instance.temporal_colocation == False):
                    print(f'Warning: {plot_type} could not be created as temporal colocation is not active')
                    continue
                #do not create timeseries bias plot if the temporal colocation is not active
                elif ('timeseries' == base_plot_type) & ('bias' in plot_options) & (self.canvas_instance.temporal_colocation == False):
                    print(f'Warning: {plot_type} could not be created as temporal colocation is not active')
                    continue
                #otherwise add information for plot type for base plot type 
                else:
                    self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[base_plot_type])

                    # change page title if have 'bias' option
                    if 'bias' in plot_options:
                        self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'] = '{} bias'.format(self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'])

            #add figsize for plot orientation
            if 'orientation' in self.canvas_instance.plot_characteristics[plot_type]:
                if self.canvas_instance.plot_characteristics[plot_type]['orientation'] == 'landscape':
                    self.canvas_instance.plot_characteristics[plot_type]['figure']['figsize'] = self.canvas_instance.landscape_figsize
                elif self.canvas_instance.plot_characteristics[plot_type]['orientation'] == 'portrait':
                    self.canvas_instance.plot_characteristics[plot_type]['figure']['figsize'] = self.canvas_instance.portrait_figsize

            #add measurement units to xlabel/ylabel if needed
            if 'xlabel' in self.canvas_instance.plot_characteristics[plot_type]:
                if self.canvas_instance.plot_characteristics[plot_type]['xlabel']['xlabel'] == 'measurement_units':
                    self.canvas_instance.plot_characteristics[plot_type]['xlabel']['xlabel'] = self.read_instance.datareader.measurement_units
            if 'ylabel' in self.canvas_instance.plot_characteristics[plot_type]:
                if self.canvas_instance.plot_characteristics[plot_type]['ylabel']['ylabel'] == 'measurement_units':
                    self.canvas_instance.plot_characteristics[plot_type]['ylabel']['ylabel'] = self.read_instance.datareader.measurement_units


    def format_axis(self, ax, base_plot_type, plot_characteristics, relevant_temporal_resolution='hour', relevant_temporal_resolution_ii=0, col_ii=0, last_valid_row=True, last_row_on_page=True):
        """Format a plotting axis.
        
        :param ax: axis object
        :type ax: object
        :param base_plot_type: plot to make, without statistical information
        :type base_plot_type: str  
        :param plot_characteristics: plot characteristics
        :type plot_characteristics: dict
        :param relevant_temporal_resolution: the relevant temporal resolution of axis (for periodic plots) 
        :type relevant_temporal_resolution: str 
        :param relevant_temporal_resolution_ii: the relevant temporal resolution axis index (for periodic plots) 
        :type relevant_temporal_resolution_ii: int
        :param col_ii: column index (for offline report)
        :type col_ii: int
        :param last_valid_row: boolean informing if last valid row to plot on (for offline report)
        :type last_valid_row: boolean
        :param last_row_on_page: boolean informing if last valid row on page (for offline report)
        :type last_row_on_page: boolean
        """

        # get plot characteristics vars
        plot_characteristics_vars = list(plot_characteristics.keys())

        #format axis

        #set axis ticks and gridlines below all artists
        ax.set_axisbelow(True)

        # make axis title?
        if 'axis_title' in plot_characteristics_vars:
            ax.set_title(**plot_characteristics['axis_title'])

        # make axis xlabel (only on last row on page/last valid row of visible axes)?
        if 'xlabel' in plot_characteristics_vars:
            if last_valid_row or last_row_on_page:
                ax.set_xlabel(**plot_characteristics['xlabel'])

        # make axis ylabel (only on leftmost column of visible axes, and for first resolution ind)?
        if ('ylabel' in plot_characteristics_vars) & (col_ii == 0) & (relevant_temporal_resolution_ii == 0):
            ax.set_ylabel(**plot_characteristics['ylabel'])

        #set xtick params ?
        if 'xtick_params' in plot_characteristics_vars:
            ax.xaxis.set_tick_params(**plot_characteristics['xtick_params'])

        #set ytick params ?
        if 'ytick_params' in plot_characteristics_vars:
            ax.yaxis.set_tick_params(**plot_characteristics['ytick_params'])

        # if are sharing xticks, and not on last row on page/last
        # valid row, then ensure current axis xticks are hidden
        if ('xtick_share' in plot_characteristics_vars) and (not last_valid_row) and (not last_row_on_page):
            plt.setp(ax.get_xticklabels(), visible=False)

        # if are sharing yticks, and not on left column, then ensure current axis yticks are hidden
        if ('ytick_share' in plot_characteristics_vars) and (col_ii != 0):
            plt.setp(ax.get_yticklabels(), visible=False)

        # set xlim?
        if 'xlim' in plot_characteristics_vars:
            ax.set_xlim(**plot_characteristics_vars['xlim'])

        # set ylim? 
        if 'ylim' in plot_characteristics_vars:
            ax.set_ylim(**plot_characteristics_vars['ylim'])

        # add gridlines (x and y)?
        if 'grid' in plot_characteristics_vars:
            ax.grid(**plot_characteristics['grid'])

        # add x gridlines?
        if 'xgrid' in plot_characteristics_vars:
            ax.xaxis.grid(**plot_characteristics['xgrid'])

        # add y gridlines?
        if 'ygrid' in plot_characteristics_vars:
            ax.yaxis.grid(**plot_characteristics['ygrid'])

        # remove spines?
        if 'remove_spines' in plot_characteristics_vars:
            for side in plot_characteristics['remove_spines']:
                ax.spines[side].set_visible(False)

        # make plot aspect ratio is equal
        # (ensure ticks and ticklabels are same also)
        if 'equal_aspect' in plot_characteristics_vars:
            ax.set_aspect('equal', adjustable='box')
            xticklocs = ax.get_xticks()
            xticklabels = ax.get_xticklabels()
            yticklocs = ax.get_yticks()
            yticklabels = ax.get_yticklabels()
            if len(xticklocs) < len(yticklocs):
                ax.set_yticks(xticklocs)
                ax.set_yticklabels(xticklabels)
            elif len(yticklocs) < len(xticklocs):
                ax.set_xticks(yticklocs)
                ax.set_xticklabels(yticklabels)

        #handle formatting specific to plot types
        if base_plot_type in ['periodic','periodic-violin']:

            # add axis resolution label 
            ax.annotate(self.canvas_instance.periodic_labels[relevant_temporal_resolution], **plot_characteristics['label'])
            # set x axis limits
            ax.set_xlim(np.min(self.canvas_instance.periodic_xticks[relevant_temporal_resolution])-0.5,np.max(self.canvas_instance.periodic_xticks[relevant_temporal_resolution])+0.5)

            # set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
            if relevant_temporal_resolution == 'hour':
                plot_characteristics['xticks'] = self.canvas_instance.periodic_xticks[relevant_temporal_resolution][::3]
                ax.set_xticks(plot_characteristics['xticks'])
            else:
                plot_characteristics['xticks'] = self.canvas_instance.periodic_xticks[relevant_temporal_resolution]
                ax.set_xticks(plot_characteristics['xticks'])
                ax.set_xticklabels([self.canvas_instance.temporal_axis_mapping_dict[relevant_temporal_resolution][xtick] for xtick
                                                                            in self.canvas_instance.periodic_xticks[relevant_temporal_resolution]])
        
        elif base_plot_type == 'map':

            #add land polygons
            ax.add_feature(self.canvas_instance.feature)

            #add gridlines ?
            if 'gridlines' in plot_characteristics_vars:
                ax.gridlines(crs=self.canvas_instance.datacrs, **plot_characteristics['gridlines'])

            #set map_extent
            if hasattr(self.read_instance, 'map_extent'):
                map_extent = self.read_instance.map_extent
            else:
                map_extent = plot_characteristics['map_extent']
                self.read_instance.map_extent = map_extent
            ax.set_extent(map_extent, crs=self.canvas_instance.datacrs)

    def make_legend_handles(self, plot_characteristics_legend, plot_options=[]):
        """Make legend element handles
        
        :param plot_characteristics_legend: plot characteristics for relevant legend
        :type plot_characteristics_legend: dict
        :return: plot_characteristics_legend with handles updated
        :rtype: dict
        :param plot_options: list of options to configure plot  
        :type plot_options: list        
        """

        # create legend elements
        # add observations element
        legend_elements = [Line2D([0], [0], 
                                  marker=plot_characteristics_legend['plot']['handles']['marker'], 
                                  color=plot_characteristics_legend['plot']['handles']['color'],
                                  markerfacecolor=self.read_instance.datareader.plotting_params['observations']['colour'],
                                  markersize=plot_characteristics_legend['plot']['handles']['markersize'], 
                                  label=plot_characteristics_legend['plot']['handles']['obs_label'])]
                                  
        # add element for each experiment
        for experiment_ind, experiment in enumerate(sorted(list(self.read_instance.datareader.data_in_memory.keys()))):
            if experiment != 'observations':
                # add experiment element
                legend_elements.append(Line2D([0], [0], 
                                              marker=plot_characteristics_legend['plot']['handles']['marker'],  
                                              color=plot_characteristics_legend['plot']['handles']['color'],
                                              markerfacecolor=self.read_instance.datareader.plotting_params[experiment]['colour'],
                                              markersize=plot_characteristics_legend['plot']['handles']['markersize'],
                                              label=self.read_instance.experiments_legend[experiment]))
        
        plot_characteristics_legend['plot']['handles'] = legend_elements
        return plot_characteristics_legend

    def make_experiment_domain_polygons(self, plot_options=[]):
        """Make experiment domain polygons
        
        :return: grid_edge_polygons
        :rtype: list
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """

        grid_edge_polygons = []

        # iterate through read experiments and plot grid domain edges on map
        for experiment in sorted(list(self.read_instance.datareader.data_in_memory.keys())):
            if experiment != 'observations':
                # create matplotlib polygon object from experiment grid edge map projection coordinates
                grid_edge_outline_poly = \
                    Polygon(np.vstack((self.read_instance.datareader.plotting_params[experiment]['grid_edge_longitude'],
                                       self.read_instance.datareader.plotting_params[experiment]['grid_edge_latitude'])).T,
                                       edgecolor=self.read_instance.datareader.plotting_params[experiment]['colour'],
                                       transform=self.canvas_instance.datacrs,
                                       **self.canvas_instance.experiment_domain_polygon)
                # append polygon
                grid_edge_polygons.append(grid_edge_outline_poly)
            
        return grid_edge_polygons

    def make_metadata(self, relevant_axis, plot_characteristics, plot_options=[]):
        """Make metadata summary plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param plot_characteristics: plot characteristics 
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """

        print('MAKE METADATA')

        # get some details of the station metadata axis --> to set limit for wrapping text
        # get axis bounding box
        ax_bbox = relevant_axis.get_window_extent().transformed(self.canvas_instance.figure.dpi_scale_trans.inverted())
        # get axis dimensions in inches
        ax_width_inches = ax_bbox.width
        # get axis dimensions in pixels
        ax_width_px = ax_width_inches * self.canvas_instance.figure.dpi

        # initialise string to plot on axis
        str_to_plot = ''

        # setup some variables for handling if just one or multiple stations selected
        if len(self.canvas_instance.relative_selected_station_inds) == 1:
            var_str_name = 'name_one' 
            str_to_plot += '{}'.format(self.read_instance.station_references[self.canvas_instance.relative_selected_station_inds][0])
        else:
            var_str_name = 'name_multiple' 
            str_to_plot += '{} Stations'.format(len(self.canvas_instance.relative_selected_station_inds))
        #iterate n vars per line and add spacing
        current_n_vars_per_line = 1
        str_to_plot += (' '*plot_characteristics['var_spacing'])

        #non-GHOST (add longitude and latitude)
        if self.read_instance.reading_nonghost:
            if len(self.canvas_instance.relative_selected_station_inds) == 1:
                #lon
                str_to_plot += '{}: {}'.format(plot_characteristics['non-ghost_vars']['longitude']['name_one'],
                    round(self.read_instance.datareader.station_latitudes[self.canvas_instance.relative_selected_station_inds][0], plot_characteristics['non-ghost_vars']['longitude']['dp']))
                #spacing
                str_to_plot += (' '*plot_characteristics['var_spacing'])
                #lat
                str_to_plot += '{}: {}'.format(plot_characteristics['non-ghost_vars']['latitude']['name_one'],
                    round(self.read_instance.datareader.station_latitudes[self.canvas_instance.relative_selected_station_inds][0], plot_characteristics['non-ghost_vars']['latitude']['dp']))

        #GHOST
        else:
            #iterate through defined variables and add them
            for ghost_var, ghost_var_dict in plot_characteristics['ghost_vars'].items():
                
                #check var str name exists in ghost var dict, if not move to next var
                if var_str_name not in ghost_var_dict:
                    continue

                #round decimal places if float
                if 'dp' in ghost_var_dict:
                    str_to_plot += '{}: {}'.format(ghost_var_dict[var_str_name],
                        round(np.nanmedian(self.read_instance.datareader.metadata_in_memory[ghost_var][
                            self.canvas_instance.relative_selected_station_inds].astype(np.float32)), ghost_var_dict['dp']))
                #if str then get unique elements or percentage dependent on n uniques
                else:
                    # gather all selected station metadata for current meta variable
                    all_current_meta = self.read_instance.datareader.metadata_in_memory[ghost_var][
                        self.canvas_instance.relative_selected_station_inds].flatten().astype(np.str)

                    # get counts of all unique metadata elements across selected stations
                    unique_meta, meta_counts = np.unique(all_current_meta, return_counts=True)
                    # get number of unique metadata elements across selected stations
                    n_unique_meta = len(unique_meta)

                    # 1 unique metadata element? then return it
                    if n_unique_meta == 1:
                        str_to_plot += '{}: {}'.format(ghost_var_dict[var_str_name], unique_meta[0])
                    # if have > N unique metadata elements, just return count of the elements across the selected stations
                    elif n_unique_meta > plot_characteristics['max_uniques']:
                        str_to_plot += '{}: {} uniques\n'.format(ghost_var_dict[var_str_name], n_unique_meta)
                    # otherwise, get percentage of unique metadata elements across selected stations
                    else:
                        meta_pc = (100. / len(all_current_meta)) * meta_counts
                        meta_pc = ['{:.1f}%'.format(meta) for meta in meta_pc]
                        # create string for variable to plot
                        str_to_plot += '{}: {}\n'.format(ghost_var_dict[var_str_name], ', '.join(
                            [':'.join([str(var), pc]) for var, pc in zip(unique_meta, meta_pc)]))

                #add units
                if 'units' in ghost_var_dict:
                    if ghost_var_dict['units'] != '':
                        str_to_plot += ' {}'.format(ghost_var_dict['units'])

                #iterate current_n_vars_per_line
                current_n_vars_per_line += 1

                #if are on limit of vars allowed per line then break to new line
                if current_n_vars_per_line == plot_characteristics['max_vars_per_row']:
                    str_to_plot += '\n'
                    current_n_vars_per_line = 0
                #otherwise, add spacing between variables on line
                else:
                    str_to_plot += (' '*plot_characteristics['var_spacing'])

        # plot string to axis
        plot_txt = relevant_axis.text(0.0, 1.0, str_to_plot, transform=relevant_axis.transAxes, **plot_characteristics['plot'])

        # modify limit to wrap text as axis width in pixels  --> hack as matplotlib
        # automatically sets limit as figure width
        plot_txt._get_wrap_line_width = lambda: ax_width_px

    def make_header(self, plot_characteristics, plot_options=[]):
        """Make header
        
        :param plot_characteristics: plot characteristics 
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """

        # set header title
        page = plt.figure(**plot_characteristics['figure'])
        if hasattr(self.read_instance, 'report_title'):
            txt = self.read_instance.report_title
        else:
            txt = 'Providentia Offline Report'
        plot_characteristics['page_title']['s'] = txt   
        plot_characteristics['page_title']['transform'] = page.transFigure
        page.text(**plot_characteristics['page_title'])

        #set header main text
        txt = 'Network = {}\nTemporal Resolution = {}\n' \
              'Species = {}\nDate Range = {} - {}\nExperiments = {}\n' \
              'Subsets = {}\n' \
            .format(self.read_instance.active_network,
                    self.read_instance.active_resolution,
                    self.read_instance.active_species,
                    self.read_instance.active_start_date,
                    self.read_instance.active_end_date, 
                    list(self.read_instance.experiments_legend.values()),
                    self.read_instance.station_subset_names)
        plot_characteristics['page_text']['s'] = txt   
        plot_characteristics['page_text']['transform'] = page.transFigure
        page.text(**plot_characteristics['page_text'])

        self.pdf.savefig(page, dpi=self.canvas_instance.dpi)
        plt.close(page)

    def make_map(self, relevant_axis, z_statistic, plot_characteristics, plot_options=[]):
        """make map plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param z_statistic: calculated z statistic to plot
        :type z_statistic: np.array
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """
      
        # plot new station points on map - coloured by currently active z statisitic
        relevant_axis.scatter(self.read_instance.datareader.station_longitudes[self.canvas_instance.active_map_valid_station_inds], 
                              self.read_instance.datareader.station_latitudes[self.canvas_instance.active_map_valid_station_inds], 
                              c=z_statistic, transform=self.canvas_instance.datacrs,
                              **plot_characteristics['plot'])

    def make_timeseries(self, relevant_axis, data_label, plot_characteristics, plot_options=[]):
        """make timeseries plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param data_label: name of data array to plot
        :type data_label: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """

        print('MAKE TIMESERIES')
        
        #bias plot?
        if 'bias' in plot_options:
            #observations and experiment must be colocated for this
            ts_obs = self.canvas_instance.selected_station_data['observations_colocatedto_experiments']['pandas_df']
            original_data_label = data_label.split('_colocatedto')[0]
            ts_model = self.canvas_instance.selected_station_data['{}_colocatedto_observations'.format(original_data_label)]['pandas_df'] 
            ts = ts_model - ts_obs
        else:
            ts = self.canvas_instance.selected_station_data[data_label]['pandas_df']

        #configure size of plots if have very few points 
        if len(ts.dropna().index) < plot_characteristics['markersize_npoints_threshold']:
            markersize = plot_characteristics['markersize']['few_points'] 
        else:
            markersize = plot_characteristics['markersize']['standard'] 

        # make time series plot
        relevant_axis.plot(ts.dropna(), 
                           color=self.read_instance.datareader.plotting_params[data_label]['colour'], 
                           markersize=markersize,
                           **plot_characteristics['plot'])



        # get user-defined characteristics for xticks
        n_slices = plot_characteristics['xtick_alteration']['n_slices']
        last_step = plot_characteristics['xtick_alteration']['last_step']
        define = plot_characteristics['xtick_alteration']['define']

        # get number of months and days
        timeseries_start_date = datetime.strptime(str(self.read_instance.active_start_date), '%Y%m%d')
        timeseries_end_date = datetime.strptime(str(self.read_instance.active_end_date), '%Y%m%d')
        n_months = (timeseries_end_date.year - timeseries_start_date.year) * 12 + (timeseries_end_date.month - timeseries_start_date.month)
        n_days = (timeseries_end_date - timeseries_start_date).days
        
        # get months and days (steps)
        label = list(self.canvas_instance.selected_station_data.keys())[-1]
        steps = self.canvas_instance.selected_station_data[label]['pandas_df'].dropna().index.values
        months = pd.date_range(timeseries_start_date, timeseries_end_date, freq='MS')
        # drop last month if it is not complete
        if months[-1] != steps[-1]:
            months = months[:-1]

        # define time slices
        if define:
            if n_months >= 3:
                steps = months
            elif n_days < 7:
                relevant_axis.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M'))
            slices = int(np.ceil(len(steps) / n_slices))

            # use default axes if slices is below 1 (do not define)
            if slices >= 1:
                xticks = steps[0::slices]
                if last_step and (steps[-1] not in xticks):
                    xticks = np.append(xticks, steps[-1])
                relevant_axis.xaxis.set_ticks(xticks) 

    def make_periodic(self, relevant_axis, data_label, plot_characteristics, zstat=None, plot_options=[]):
        """Make period or period-violin plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param data_label: name of data array to plot
        :type data_label: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param zstat: name of statistic
        :type zstat: str
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """

        print('MAKE PERIODIC')

        # iterate through all relevant temporal aggregation resolutions
        for relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions:

            #get subplot axis
            relevant_subplot_axis = relevant_axis[relevant_temporal_resolution]

            #violin plot type?
            if not zstat:

                # get grouped data for current temporal aggregation resolution
                grouped_data = self.canvas_instance.selected_station_data[data_label][relevant_temporal_resolution]['grouped_data']
                # drop any groups which have no data
                grouped_data = [group for group in grouped_data if len(group) > 0]

                # make violin plot 
                violin_plot = relevant_subplot_axis.violinplot(grouped_data, positions=self.canvas_instance.selected_station_data[data_label][relevant_temporal_resolution]['valid_xticks'], **plot_characteristics['plot_violin'])

                #plot p50
                xticks = self.canvas_instance.periodic_xticks[relevant_temporal_resolution]
                medians = self.canvas_instance.selected_station_data[data_label][relevant_temporal_resolution]['p50']
                median_zorder = self.read_instance.datareader.plotting_params[data_label]['zorder']+len(list(self.read_instance.datareader.plotting_params.keys()))
                # split arrays if there are any temporal gaps to avoid
                # line drawn being interpolated across missing values
                inds_to_split = np.where(np.diff(xticks) > 1)[0]
                if len(inds_to_split) == 0:
                    relevant_subplot_axis.plot(xticks, medians, color=self.read_instance.datareader.plotting_params[data_label]['colour'], zorder=median_zorder, **plot_characteristics['plot_p50'])
                else:
                    inds_to_split += 1
                    start_ind = 0
                    for end_ind in inds_to_split:
                        relevant_subplot_axis.plot(xticks[start_ind:end_ind], medians[start_ind:end_ind], color=self.read_instance.datareader.plotting_params[data_label]['colour'], zorder=median_zorder, **plot_characteristics['plot_p50'])
                        start_ind = end_ind
                    relevant_subplot_axis.plot(xticks[start_ind:], medians[start_ind:], color=self.read_instance.datareader.plotting_params[data_label]['colour'], zorder=median_zorder, **plot_characteristics['plot_p50'])

                # update plotted objects with necessary colour and alpha
                for patch in violin_plot['bodies']:
                    patch.set_facecolor(self.read_instance.datareader.plotting_params[data_label]['colour'])
                    if data_label.split('_')[0] == 'observations':
                        patch.set_alpha(plot_characteristics['patch']['alpha_obs'])
                    else:
                        patch.set_alpha(plot_characteristics['patch']['alpha_exp'])
                    # if have at least 1 experiment data array, split the violin plot across the horizontal
                    # (observations on left, experiment violin_plots on right)
                    if (len(list(self.canvas_instance.selected_station_data.keys())) > 1) & ('individual' not in plot_options):
                        m = np.mean(patch.get_paths()[0].vertices[:, 0])
                        # observations on left
                        if data_label.split('_')[0] == 'observations':
                            patch.get_paths()[0].vertices[:, 0] = np.clip(patch.get_paths()[0].vertices[:, 0], -np.inf, m)
                        # experiments on right
                        else:
                            patch.get_paths()[0].vertices[:, 0] = np.clip(patch.get_paths()[0].vertices[:, 0], m, np.inf)

            #standard periodic plot type
            else:
                # make plot
                relevant_subplot_axis.plot(self.canvas_instance.periodic_xticks[relevant_temporal_resolution],self.canvas_instance.selected_station_data[data_label][relevant_temporal_resolution][zstat],
                                           color=self.read_instance.datareader.plotting_params[data_label]['colour'], zorder=self.read_instance.datareader.plotting_params[data_label]['zorder'],
                                           **plot_characteristics['plot'])
                    
                # plot horizontal line/s across x axis at value/s of minimum experiment bias
                zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat)
                if z_statistic_sign == 'bias':
                    # get value/s of minimum bias for statistic
                    if z_statistic_type == 'basic':
                        minimum_bias = basic_stats[base_zstat]['minimum_bias']
                    else:
                        minimum_bias = expbias_stats[base_zstat]['minimum_bias']
                    for mb in minimum_bias:
                        relevant_subplot_axis.axhline(y=mb, **plot_characteristics['bias_line'])

    def make_distribution(self, relevant_axis, data_label, plot_characteristics, plot_options=[]):
        """Make distribution plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param data_label: name of data array to plot
        :type data_label: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """

        # make distribution plot
        minmax_diff = self.canvas_instance.selected_station_data_max - self.canvas_instance.selected_station_data_min
        if pd.isnull(self.read_instance.parameter_dictionary[self.read_instance.active_species]['minimum_resolution']):
            n_samples = plot_characteristics['pdf_min_samples']
        else:
            n_samples = int(np.around(minmax_diff/(self.read_instance.parameter_dictionary[self.active_species]['minimum_resolution']/4.0),0))
            if n_samples < plot_characteristics['pdf_min_samples']:
                n_samples = plot_characteristics['pdf_min_samples']
        x_grid = np.linspace(self.canvas_instance.selected_station_data_min, self.canvas_instance.selected_station_data_max, n_samples, endpoint=True)

        # setup bias plot
        if 'bias' in plot_options:
            if self.canvas_instance.temporal_colocation:
                PDF_obs = st.gaussian_kde(self.canvas_instance.selected_station_data['observations_colocatedto_experiments']['pandas_df']['data'].dropna())
            else:
                PDF_obs = st.gaussian_kde(self.canvas_instance.selected_station_data['observations']['pandas_df']['data'].dropna())
            PDF_model = st.gaussian_kde(self.canvas_instance.selected_station_data[data_label]['pandas_df']['data'].dropna())
            PDF_sampled = PDF_model(x_grid) - PDF_obs(x_grid)
            # plot horizontal line across x axis at 0
            relevant_axis.axhline(**plot_characteristics['bias_line'])

        # setup standard plot
        else:
            PDF = st.gaussian_kde(self.canvas_instance.selected_station_data[data_label]['pandas_df']['data'].dropna())
            PDF_sampled = PDF(x_grid)

        # make plot
        relevant_axis.plot(x_grid, PDF_sampled, color=self.read_instance.datareader.plotting_params[data_label]['colour'], **plot_characteristics['plot'])

    def make_scatter(self, relevant_axis, data_label, plot_characteristics, plot_options=[]):
        """Make scatter plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param data_label: name of data array to plot
        :type data_label: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """

        # get observations data
        observations_data = self.canvas_instance.selected_station_data['observations_colocatedto_experiments']['pandas_df'].dropna()['data']

        # get experiment data
        experiment_data = self.canvas_instance.selected_station_data[data_label]['pandas_df'].dropna()['data']

        # create scatter plot
        relevant_axis.plot(observations_data, experiment_data, 
                           color=self.read_instance.datareader.plotting_params[data_label]['colour']
                           **plot_characteristics['plot'])

        # add 1:1 line
        relevant_axis.plot([0, 1], [0, 1], transform=relevant_axis.transAxes, **plot_characteristics['1:1_line'])
        # add 1:2 line
        relevant_axis.plot([0, 1], [0, 0.5], transform=relevant_axis.transAxes, **plot_characteristics['1:2_line'])
        # add 2:1 line
        relevant_axis.plot([0, 0.5], [0, 1], transform=relevant_axis.transAxes, **plot_characteristics['2:1_line'])

    def make_heatmap(self, relevant_axis, stat_df, plot_characteristics, plot_options=[]):
        """Make heatmap plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param stat_df: dataframe of calculated statistical information
        :type stat_df: object
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """
        # determine if want to add annotations or not from plot_options
        if 'annotate' in plot_options:
            annotate = True
        else:
            annotate = False

        # plot heatmap
        ax = sns.heatmap(stat_df, 
                         ax=relevant_axis, annot=annotate,
                         **plot_characteristics['plot'])

        # axis cuts off due to bug in matplotlib 3.1.1 - hack fix. Remove in Future!
        bottom, top = relevant_axis.get_ylim()
        relevant_axis.set_ylim(bottom + 0.5, top - 0.5)

    def make_table(self, relevant_axis, stat_df, plot_characteristics, plot_options=[]):
        """Make table plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param stat_df: dataframe of calculated statistical information
        :type stat_df: object
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plot  
        :type plot_options: list
        """

        # turn off axis to make table
        relevant_axis.axis('off')

        # make table
        table = relevant_axis.table(cellText=stat_df.values, colLabels=stat_df.columns, rowLabels=stat_df.index, loc='center')
        #table.set_fontsize(18)

    def log_axes(self, relevant_axis, log_ax):
        """Log plot axes

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param log_ax: which axis to log
        :type log_ax: str
        """

        if log_ax == 'logx':
            relevant_axis.set_xscale('log')
        if log_ax == 'logy':
            relevant_axis.set_yscale('log')

    def linear_regression(self, relevant_axis, data_labels, plot_characteristics, plot_options=[]):
        """Add linear regression to plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param data_labels: names of plotted data arrays (names potentially modified by colocation)  
        :type data_labels: list
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plots
        :type plot_options: list
        """
    
        #get observations data
        observations_data = self.canvas_instance.selected_station_data['observations_colocatedto_experiments']['pandas_df'].dropna()['data']

        #iterate through experiment data, making regression line to obseervations
        for data_label in data_labels:
            if data_label.split('_')[0] != 'observations':
                experiment_data = self.canvas_instance.selected_station_data[data_label]['pandas_df'].dropna()['data']

                m, b = np.polyfit(observations_data, experiment_data, deg=1)
                relevant_axis.plot(observations_data, m*observations_data+b, 
                                   color=self.read_instance.datareader.plotting_params[data_label]['colour'],
                                   zorder=self.read_instance.datareader.plotting_params[data_label]['zorder']+len(list(self.read_instance.datareader.plotting_params.keys())),
                                   transform=relevant_axis.transAxes,
                                   **plot_characteristics['regression'])

    def trend(self, relevant_axis, data_labels, plot_characteristics, plot_options=[]):
        """Add trendline to plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param data_labels: names of plotted data arrays (names potentially modified by colocation)  
        :type data_labels: list
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plots
        :type plot_options: list
        """

        #iterate through plotted data arrays making trendline
        for data_label in data_labels:

            #bias plot?
            if 'bias' in plot_options:
                #observations and experiment must be colocated for this
                ts_obs = self.canvas_instance.selected_station_data['observations_colocatedto_experiments']['pandas_df']
                original_data_label = data_label.split('_colocatedto')[0]
                ts_model = self.canvas_instance.selected_station_data['{}_colocatedto_observations'.format(original_data_label)]['pandas_df'] 
                ts = ts_model - ts_obs
            #normal plot?
            else:
                ts = self.canvas_instance.selected_station_data[data_label]['pandas_df']
            #make trendline
            relevant_axis.plot(ts.rolling(plot_characteristics['trend']['window'], min_periods=plot_characteristics['trend']['min_points'], center=True).mean().dropna(),
                               color=self.read_instance.datareader.plotting_params[data_label]['colour'],
                               zorder=self.read_instance.datareader.plotting_params[data_label]['zorder']+len(list(self.read_instance.datareader.plotting_params.keys())),
                               **plot_characteristics['trend']['trend_format'])


    def annotation(self, relevant_axis, data_labels, plot_characteristics, plot_options=[]):
        """Add statistical annotations to plot

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param data_labels: names of plotted data arrays (names potentially modified by colocation)  
        :type data_labels: list
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plots
        :type plot_options: list
        """
        
        #get stats wished to be annotated
        stats = plot_characteristics['annotate_stats']

        #if no stats defined, then return
        if len(stats) == 0:
            print(f'No annotation statistics have not been defined for {plot_type} in plot_characteristics_offline.py')
            return

        #initialise list of strs to annotate, and colours of annotations
        str_to_annotate = []
        colours = []

        #iterate through plotted data labels
        for data_label in data_labels:

            # get stats
            stats_annotate = []
            for zstat in stats:
                if zstat in list(self.canvas_instance.selected_station_data[data_label]['all']):
                    stats_annotate.append(zstat + ': ' + str(round(self.canvas_instance.selected_station_data[data_label]['all'][zstat][0], plot_characteristics['annotate_text']['round_decimal_places'])))

            # show number of stations if defined
            if plot_characteristics['annotate_text']['n_stations']:
                colours.append('black')
                if 'individual' in plot_options:
                    str_to_annotate.append('Stations: 1')
                else:
                    str_to_annotate.append('Stations: ' + str(len(self.read_instance.datareader.plotting_params['observations']['valid_station_inds'])))

            # get colors
            colours.append(self.read_instance.datareader.plotting_params[data_label]['colour'])

            # generate annotation
            str_to_annotate.append(', '.join(stats_annotate))

        # add annotation to plot
        # see loc options at https://matplotlib.org/3.1.0/api/offsetbox_api.html
        lines = [TextArea(line, textprops=dict(color=colour, size=plot_characteristics['annotate_text']['fontsize'])) 
                            for line, colour in zip(str_to_annotate, colours)]
        bbox = AnchoredOffsetbox(child=VPacker(children=lines, align="left", pad=0, sep=1),
                                    loc=plot_characteristics['annotate_text']['loc'],
                                    bbox_transform=relevant_axis.transAxes)
        bbox.patch.set(**plot_characteristics['annotate_bbox']) 
        relevant_axis.add_artist(bbox)

    def harmonise_xy_lims_paradigm(self, relevant_axs, base_plot_type, plot_characteristics, plot_options):
        """Harmonises xy limits across paradigm of plot type, unless axis limits have been defined
        
        :param relevant_axs: relevant axes
        :type relevant_axs: list
        :param base_plot_type: plot type name
        :type base_plot_type: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plots
        :type plot_options: list
        """

        # initialise arrays to save lower and upper limits in all axes
        all_xlim_lower = []
        all_xlim_upper = []
        all_ylim_lower = []
        all_ylim_upper = []

        #initialise variables for setting axis limits
        xlim_min = None
        xlim_max = None
        ylim_min = None
        ylim_max = None

        # get lower and upper limits across all relevant axes
        for ax in relevant_axs:
            if ax.lines:
                xlim_lower, xlim_upper = ax.get_xlim()
                ylim_lower, ylim_upper = ax.get_ylim()
                all_xlim_lower.append(xlim_lower)
                all_xlim_upper.append(xlim_upper)
                all_ylim_lower.append(ylim_lower)
                all_ylim_upper.append(ylim_upper)

        # get minimum and maximum from all axes and set limits
        for ax in relevant_axs:
            if ax.lines:

                #get xlims
                if base_plot_type not in ['periodic','periodic-violin']:
                    xlim_min = np.min(all_xlim_lower)
                    xlim_max = np.max(all_xlim_upper)

                #get ylims
                if 'bias' in plot_options:
                    # if there is bias center plots y limits around 0
                    if np.abs(np.max(all_ylim_upper)) >= np.abs(np.min(all_ylim_lower)):
                        ylim_min = -np.abs(np.max(all_ylim_upper))
                        ylim_max = np.abs(np.max(all_ylim_upper))
                    elif np.abs(np.max(all_ylim_upper)) < np.abs(np.min(all_ylim_lower)):
                        ylim_min = -np.abs(np.min(all_ylim_lower))
                        ylim_max = np.abs(np.min(all_ylim_lower))
                else:
                    ylim_min = np.min(all_ylim_lower) 
                    ylim_max = np.max(all_ylim_upper)
            
            #set xlim
            if (xlim_min) & (xlim_max) & ('xlim' not in plot_characteristics):
                ax.set_xlim(xlim_min, xlim_max)

            #set ylim
            if (ylim_min) & (ylim_max) & ('ylim' not in plot_characteristics):
                ax.set_ylim(ylim_min, ylim_max)

    def set_axis_title(self, relevant_axis, title, plot_characteristics):
        """Set title of plot axis

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param title: axis title
        :type title: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        """    

        axis_title_characteristics = copy.deepcopy(plot_characteristics['axis_title'])
        axis_title_characteristics['label'] = title

        relevant_axis.set_title(**axis_title_characteristics)

    def set_axis_label(self, relevant_axis, label_ax, label, plot_characteristics):
        """Set label of plot axis

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param label_ax: which axis to set label of
        :type label_ax: str
        :param label: axis label
        :type label: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        """

        #return if label is empty str
        if label == '':
            return

        #get appropriate axis for plotting label for plots with multiple sub-axes
        if type(relevant_axis) == dict:
            for sub_ax in relevant_axis.values():
                relevant_axis = sub_ax
                break

        #set label
        if label_ax == 'x':
            axis_label_characteristics = copy.deepcopy(plot_characteristics['xlabel'])
            axis_label_characteristics['xlabel'] = label
            relevant_axis.set_xlabel(**axis_label_characteristics) 
        elif label_ax == 'y':
            axis_label_characteristics = copy.deepcopy(plot_characteristics['ylabel'])
            axis_label_characteristics['ylabel'] = label
            relevant_axis.set_ylabel(**axis_label_characteristics) 