import copy
from datetime import datetime
import json
import os
import time

import math
import pyproj
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from KDEpy import FFTKDE
import matplotlib as mpl 
from matplotlib.dates import num2date
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf
import matplotlib.style as mplstyle
import numpy as np
import pandas as pd
import scipy.stats as st
import seaborn as sns
from itertools import groupby

from .statistics import get_z_statistic_info, calculate_statistic, boxplot_inner_fences
from .aux import (get_land_polygon_resolution, temp_axis_dict, periodic_xticks, periodic_labels,
                  get_multispecies_aliases, show_message, kde_fft)
from .read_aux import drop_nans

# speed up transformations in cartopy
pyproj.set_use_global_context()

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))

class Plot:
    """ Class that makes plots and handles plot configuration options when defined. """

    def __init__(self, read_instance=None, canvas_instance=None):
        
        self.read_instance = read_instance
        self.canvas_instance = canvas_instance

        # set cartopy data directory
        cartopy.config['pre_existing_data_dir'] = self.read_instance.cartopy_data_dir

        # define projections for map plot and actual geographic coordinates
        self.canvas_instance.datacrs = ccrs.PlateCarree()
        proj_class = getattr(ccrs, self.canvas_instance.plot_characteristics_templates['map']['projection']) 
        self.canvas_instance.plotcrs = proj_class()
        
        # set miscellaneous vars
        self.canvas_instance.temporal_axis_mapping_dict = temp_axis_dict()
        self.canvas_instance.periodic_xticks = periodic_xticks()
        self.canvas_instance.periodic_labels = periodic_labels()

    def set_plot_characteristics(self, plot_types, zstat=False):
        """ Iterate through all plots to make, and determine if they can and cannot be made.
            Update plot characteristics associated with specific plot types due to plot options. 

            :param plot_types: plot types to create 
            :type plot_types: list  
            :param zstat: z statistic str 
            :type zstat: str
        """

        # add all valid defined plots to plot_characteristics
        for plot_type in plot_types:
            
            # initialize condition
            valid_plot_type = True

            # do not create empty plots
            if plot_type == 'None':
                continue

            # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
            plot_options = plot_type.split('_')[1:]

            # get zstat information from plot_type
            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type)
            
            # remove plots where setting 'obs' and 'bias' options together
            if ('obs' in plot_options) & ('bias' in plot_options): 
                if self.read_instance.offline:
                    print(f"Warning: {plot_type} cannot not be created as 'obs' and 'bias' options set together.")
                    valid_plot_type = False

            # if no experiments are defined, remove all bias plots 
            if ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                if len(self.read_instance.data_labels) == 1:
                    if self.read_instance.offline:
                        print(f'Warning: No experiments defined, so {plot_type} bias plot cannot be created.')
                        valid_plot_type = False

            # if are making an experiment bias plot, and temporal_colocation is off, then remove plot
            if (z_statistic_type == 'expbias') & (not self.read_instance.temporal_colocation):
                if self.read_instance.offline:
                    print(f'Warning: To calculate the experiment bias stat {zstat}, temporal_colocation must be set to True, so {plot_type} plot cannot be created.')
                    valid_plot_type = False

            # add new keys to make plots with stats (map, periodic, heatmap, table)
            if zstat:

                # get base plot type (without stat and options)
                base_plot_type = plot_type.split('-')[0] 
                # combine basic and expbias stats dicts together
                stats_dict = {**basic_stats, **expbias_stats}
                
                # check all defined plot options are allowed for current plot type
                if not all(plot_option in self.canvas_instance.plot_characteristics_templates[base_plot_type]['plot_options'] for plot_option in plot_options):
                    if self.read_instance.offline:
                        print(f'Warning: {plot_type} cannot be created as some plot options are not valid.')
                        valid_plot_type = False
                # check desired statistic is defined in stats dict
                if base_zstat not in stats_dict:
                    if self.read_instance.offline:
                        print(f"Warning: {plot_type} cannot be created as {base_zstat} not defined in Providentia's statistical library.")
                        valid_plot_type = False
                # remove plots where setting 'obs', but z_statistic_sign is 'bias'
                elif ('obs' in plot_options) & (z_statistic_sign == 'bias'):
                    if self.read_instance.offline:
                        print(f"Warning: {plot_type} cannot be created as are plotting a bias statistic but 'obs' option is set.")
                        valid_plot_type = False

                if not valid_plot_type:
                    if plot_type in self.read_instance.summary_plots_to_make:
                        self.read_instance.summary_plots_to_make.remove(plot_type)
                    if plot_type in self.read_instance.station_plots_to_make:
                        self.read_instance.station_plots_to_make.remove(plot_type)
                    continue

                # add information for plot type from base plot type template 
                self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[base_plot_type])

                # set page title 
                if 'page_title' in self.canvas_instance.plot_characteristics[plot_type]:
                    if 't' in self.canvas_instance.plot_characteristics[plot_type]['page_title'].keys():
                        if 'bias' in plot_options:
                            self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'] = '{} {} bias'.format(self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'],
                                                                                                                          stats_dict[base_zstat]['label'])
                        else:
                            self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'] = '{} {}'.format(self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'], 
                                                                                                                     stats_dict[base_zstat]['label'])
            # add new keys for plots without stats
            else:

                # get base plot type (without options)
                base_plot_type = plot_type.split('_')[0] 

                # check all defined plot options are allowed for current plot type
                if not all(plot_option in self.canvas_instance.plot_characteristics_templates[base_plot_type]['plot_options'] for plot_option in plot_options):
                    if self.read_instance.offline:
                        print(f'Warning: {plot_type} cannot be created as some plot options are not valid.')
                        valid_plot_type = False
                # warning for scatter plot if the temporal colocation is not active
                elif ('scatter' == base_plot_type) & (not self.read_instance.temporal_colocation):
                    if self.read_instance.offline:
                        print(f'Warning: {plot_type} cannot be created as temporal colocation is not active.')
                        valid_plot_type = False
                # warning for timeseries bias plot if the temporal colocation is not active
                elif ('timeseries' == base_plot_type) & ('bias' in plot_options) & (not self.read_instance.temporal_colocation):
                    if self.read_instance.offline:
                        print(f'Warning: {plot_type} cannot be created as temporal colocation is not active.')
                        valid_plot_type = False

                # break loop if the plot type is not valid and remove plot type from lists
                if not valid_plot_type:
                    if plot_type in self.read_instance.summary_plots_to_make:
                        self.read_instance.summary_plots_to_make.remove(plot_type)
                    if plot_type in self.read_instance.station_plots_to_make:
                        self.read_instance.station_plots_to_make.remove(plot_type)
                    continue

                # add information for plot type for base plot type 
                self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[base_plot_type])

                # change page title if have 'bias' option
                if 'page_title' in self.canvas_instance.plot_characteristics[plot_type]:
                    if 'bias' in plot_options:
                        self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'] = '{} bias'.format(self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'])

            # add figsize for plot orientation
            if 'orientation' in self.canvas_instance.plot_characteristics[plot_type]:
                if self.canvas_instance.plot_characteristics[plot_type]['orientation'] == 'landscape':
                    self.canvas_instance.plot_characteristics[plot_type]['figure']['figsize'] = self.canvas_instance.landscape_figsize
                elif self.canvas_instance.plot_characteristics[plot_type]['orientation'] == 'portrait':
                    self.canvas_instance.plot_characteristics[plot_type]['figure']['figsize'] = self.canvas_instance.portrait_figsize

    def format_axis(self, ax, base_plot_type, plot_characteristics, 
                    col_ii=0, last_valid_row=True, last_row_on_page=True, set_extent=True,
                    relevant_temporal_resolutions=['hour','dayofweek','month']):
        """ Format a plotting axis.
        
            :param ax: axis object
            :type ax: object
            :param base_plot_type: plot to make, without statistical information
            :type base_plot_type: str  
            :param plot_characteristics: plot characteristics
            :type plot_characteristics: dict
            :param col_ii: column index (for offline report)
            :type col_ii: int
            :param last_valid_row: boolean informing if last valid row to plot on (for offline report)
            :type last_valid_row: boolean
            :param last_row_on_page: boolean informing if last valid row on page (for offline report)
            :type last_row_on_page: boolean
            :param map_extent: boolean informing if wanting to set map_extent or not
            :type map_extent: boolean
            :param relevant_temporal_resolutions: list of relevant temporal resolutions
            :type relevant_temporal_resolutions: list
        """

        # get plot characteristics vars
        plot_characteristics_vars = list(plot_characteristics.keys())

        # get appropriate axes for nested axes
        axs_to_format = []
        temporal_resolutions_per_ax = []
        if isinstance(ax, dict):
            for relevant_temporal_resolution, sub_ax in ax.items():
                if relevant_temporal_resolution in relevant_temporal_resolutions:
                    axs_to_format.append(sub_ax)
                    temporal_resolutions_per_ax.append(relevant_temporal_resolution)
        else:
            axs_to_format.append(ax)
            temporal_resolutions_per_ax.append('')

        # iterate though relevant axes (and relevant temporal resolutions for periodic plots)
        for ax_to_format, relevant_temporal_resolution in zip(axs_to_format, temporal_resolutions_per_ax): 

            # set axis ticks and gridlines below all artists
            ax_to_format.set_axisbelow(True)

            # make axis xlabel?
            if 'xlabel' in plot_characteristics_vars:
                ax_to_format.set_xlabel(**plot_characteristics['xlabel'])

            # make axis ylabel (only on leftmost column of visible axes)?
            if 'ylabel' in plot_characteristics_vars:
                ax_to_format.set_ylabel(**plot_characteristics['ylabel'])

            # set xtick params ?
            if 'xtick_params' in plot_characteristics_vars:
                ax_to_format.xaxis.set_tick_params(**plot_characteristics['xtick_params'])

            # set ytick params ?
            if 'ytick_params' in plot_characteristics_vars:
                ax_to_format.yaxis.set_tick_params(**plot_characteristics['ytick_params'])

            # if are sharing xticks, and not on last row on page/last
            # valid row, then ensure current axis xticks are hidden
            if ('xtick_share' in plot_characteristics_vars) and (not last_valid_row) and (not last_row_on_page):
                plt.setp(ax_to_format.get_xticklabels(), visible=False)

            # if are sharing yticks, and not on left column, then ensure current axis yticks are hidden
            if ('ytick_share' in plot_characteristics_vars) and (col_ii != 0):
                plt.setp(ax_to_format.get_yticklabels(), visible=False)

            # set xlim?
            if 'xlim' in plot_characteristics_vars:
                ax_to_format.set_xlim(**plot_characteristics_vars['xlim'])

            # set ylim? 
            if 'ylim' in plot_characteristics_vars:
                ax_to_format.set_ylim(**plot_characteristics_vars['ylim'])

            # add gridlines (x and y)?
            if 'grid' in plot_characteristics_vars:
                ax_to_format.grid(**plot_characteristics['grid'])

            # add x gridlines?
            if 'xgrid' in plot_characteristics_vars:
                ax_to_format.xaxis.grid(**plot_characteristics['xgrid'])

            # add y gridlines?
            if 'ygrid' in plot_characteristics_vars:
                ax_to_format.yaxis.grid(**plot_characteristics['ygrid'])

            # set x axis decimal places?
            if 'round_decimal_places' in plot_characteristics_vars:
                if 'x' in plot_characteristics['round_decimal_places']:
                    ax_to_format.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.{}f'.format(plot_characteristics['round_decimal_places']['x'])))

            # set y axis decimal places?
            if 'round_decimal_places' in plot_characteristics_vars:
                if 'y' in plot_characteristics['round_decimal_places']:
                    ax_to_format.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.{}f'.format(plot_characteristics['round_decimal_places']['y'])))

            # remove spines?
            if 'remove_spines' in plot_characteristics_vars:
                for side in plot_characteristics['remove_spines']:
                    ax_to_format.spines[side].set_visible(False)

                for side in list(set(['top', 'bottom', 'right', 'left']).symmetric_difference(plot_characteristics['remove_spines'])):
                    ax_to_format.spines[side].set_visible(True)

            # handle formatting specific to plot types
            if base_plot_type in ['periodic','periodic-violin']:

                # add axis resolution label 
                ax_to_format.annotate(self.canvas_instance.periodic_labels[relevant_temporal_resolution], **plot_characteristics['label'])

                # set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
                if relevant_temporal_resolution == 'hour':
                    plot_characteristics['xticks'] = self.canvas_instance.periodic_xticks[relevant_temporal_resolution][::3]
                    ax_to_format.set_xticks(plot_characteristics['xticks'])
                else:
                    plot_characteristics['xticks'] = self.canvas_instance.periodic_xticks[relevant_temporal_resolution]
                    ax_to_format.set_xticks(plot_characteristics['xticks'])
                    ax_to_format.set_xticklabels([self.canvas_instance.temporal_axis_mapping_dict['short'][relevant_temporal_resolution][xtick] for xtick
                                                  in self.canvas_instance.periodic_xticks[relevant_temporal_resolution]])
            # map specific formatting
            elif base_plot_type == 'map':

                # add land polygons
                feature = cfeature.NaturalEarthFeature(category='physical', name='land',
                                                       scale=get_land_polygon_resolution(self.canvas_instance.plot_characteristics_templates['map']['map_coastline_resolution']), 
                                                       **self.canvas_instance.plot_characteristics_templates['map']['land_polygon'])
                ax_to_format.add_feature(feature)

                # add gridlines ?
                if 'gridlines' in plot_characteristics_vars:
                    ax_to_format.gridlines(crs=self.canvas_instance.datacrs, 
                                           **plot_characteristics['gridlines'])

                # set map extent (if wanted)
                if set_extent:
                    self.set_map_extent(ax_to_format)

    def set_map_extent(self, ax):
        """ Set map extent, done set_xlim and set_ylim rather than set_extent 
            to avoid axis cutting off slightly (https://github.com/SciTools/cartopy/issues/697).
        """

        mlon = np.mean(self.read_instance.map_extent[:2])
        mlat = np.mean(self.read_instance.map_extent[2:])
        xtrm_data = np.array([[self.read_instance.map_extent[0], mlat], [mlon, self.read_instance.map_extent[2]], 
                              [self.read_instance.map_extent[1], mlat], [mlon, self.read_instance.map_extent[3]]])
        proj_to_data = self.canvas_instance.datacrs._as_mpl_transform(ax) - ax.transData
        xtrm = proj_to_data.transform(xtrm_data)
        ax.set_xlim(xtrm[:,0].min(), xtrm[:,0].max())
        ax.set_ylim(xtrm[:,1].min(), xtrm[:,1].max())

    def get_map_extent(self, ax):
        """ Get map extent from xlim and ylim. """ 

        # get plot extent
        coords = np.array(ax.get_extent())
        current_xlim = coords[0:2]
        current_ylim = coords[2:4]

        # calculate means
        mlon = np.mean(current_xlim)
        mlat = np.mean(current_ylim)

        # get coordinates
        xcoords = np.array([current_xlim[0], mlon, current_xlim[1], mlon])
        ycoords = np.array([mlat, current_ylim[0], mlat, current_ylim[1]])

        # transform coordinates to projected data
        transformed_coords = self.canvas_instance.datacrs.transform_points(self.canvas_instance.plotcrs, 
                                                                           xcoords, ycoords)[:, :2]
    
        # keep longitudes between -180 and 180
        lon_change = False
        if (np.isnan(transformed_coords[0, 0])) or (transformed_coords[0, 0] == -179.99999999999932):
            transformed_coords[0, 0] = -180
            lon_change = True
        if (np.isnan(transformed_coords[2, 0])) or (transformed_coords[2, 0] == 179.99999999999932):
            transformed_coords[2, 0] = 180  
            lon_change = True 

        # keep latitudes between -90 and 90
        lat_change = False
        if (np.isnan(transformed_coords[1, 1])) or (transformed_coords[1, 1] == -89.99999999999966):
            transformed_coords[1, 1] = -90
            lat_change = True  
        if (np.isnan(transformed_coords[3, 1])) or (transformed_coords[3, 1] == 89.99999999999966):
            transformed_coords[3, 1] = 90
            lat_change = True  

        # recalculate means
        if lon_change or lat_change:
            # recalculate longitude means
            mlon = np.mean(np.array([transformed_coords[0, 0], transformed_coords[2, 0]]))
            transformed_coords[1, 0] = mlon
            transformed_coords[3, 0] = mlon

            # recalculate latitude means
            mlat = np.mean(np.array([transformed_coords[1, 1], transformed_coords[3, 1]]))
            transformed_coords[0, 1] = mlat
            transformed_coords[2, 1] = mlat

        # get map extent
        map_extent = [transformed_coords[:,0].min(), transformed_coords[:,0].max(),
                      transformed_coords[:,1].min(), transformed_coords[:,1].max()]

        return map_extent

    def set_equal_axes(self, ax, plot_options):
        """ Set equal aspect and limits (useful for scatter plots). """

        # set equal aspect
        # only if no axis is in log scale
        if 'logy' in plot_options or 'logx' in plot_options:
            ax.set_aspect(aspect='auto')
        else:
            ax.set_aspect(aspect='equal', adjustable='box')

        if len(ax.lines) == 0:
            return None

        # get min and max values for ticks
        for i, line in enumerate(ax.lines):
            line_xdata = line.get_xdata()
            line_ydata = line.get_ydata()
            if list(line_xdata) != [0, 1] and list(line_xdata) != [0, 0.5]:
                break
            
        xtickmin = np.nanmin(line_xdata)
        xtickmax = np.nanmax(line_xdata)
        ytickmin = np.nanmin(line_ydata)
        ytickmax = np.nanmax(line_ydata)

        for line in ax.lines[i:]:
            if np.nanmin(line.get_xdata()) < xtickmin:
                xtickmin = np.nanmin(line.get_xdata())
            if np.nanmax(line.get_xdata()) > xtickmax:
                xtickmax = np.nanmax(line.get_xdata())
            if np.nanmin(line.get_ydata()) < ytickmin:
                ytickmin = np.nanmin(line.get_ydata())
            if np.nanmax(line.get_ydata()) > ytickmax:
                ytickmax = np.nanmax(line.get_ydata())

        # compare min and max across axes
        if xtickmin < ytickmin:
            tickmin = xtickmin
        else:
            tickmin = ytickmin
        if xtickmax > ytickmax:
            tickmax = xtickmax
        else:
            tickmax = ytickmax

        # set equal ticks
        ax.set_xlim(tickmin, tickmax)
        ax.set_ylim(tickmin, tickmax)

        return None

    def make_legend_handles(self, plot_characteristics_legend, plot_options=[]):
        """ Make legend element handles.
        
            :param plot_characteristics_legend: plot characteristics for relevant legend
            :type plot_characteristics_legend: dict
            :param plot_options: list of options to configure plot  
            :type plot_options: list    
            :return: plot_characteristics_legend with handles updated
            :rtype: dict
        """

        # create legend elements
        # add observations element
        legend_elements = [Line2D([0], [0], 
                                  marker=plot_characteristics_legend['handles']['marker'], 
                                  color=plot_characteristics_legend['handles']['color'],
                                  markerfacecolor=self.read_instance.plotting_params['observations']['colour'],
                                  markersize=plot_characteristics_legend['handles']['markersize'], 
                                  label=plot_characteristics_legend['handles']['obs_label'])]
                                  
        # add element for each experiment
        for experiment in self.read_instance.data_labels:
            if experiment != 'observations':
                # add experiment element
                legend_elements.append(Line2D([0], [0], 
                                              marker=plot_characteristics_legend['handles']['marker'],  
                                              color=plot_characteristics_legend['handles']['color'],
                                              markerfacecolor=self.read_instance.plotting_params[experiment]['colour'],
                                              markersize=plot_characteristics_legend['handles']['markersize'],
                                              label=self.read_instance.experiments[experiment]))
        
        plot_characteristics_legend['plot']['handles'] = legend_elements
        
        return plot_characteristics_legend

    def make_experiment_domain_polygons(self, plot_options=[]):
        """ Make experiment domain polygons.
        
            :param plot_options: list of options to configure plot  
            :type plot_options: list
            :return: grid_edge_polygons
            :rtype: list
        """

        grid_edge_polygons = []

        # iterate through read experiments and plot grid domain edges on map
        for experiment in self.read_instance.data_labels:
            if experiment != 'observations':
                # create matplotlib polygon object from experiment grid edge map projection coordinates
                grid_edge_outline_poly = \
                    Polygon(np.vstack((self.read_instance.plotting_params[experiment]['grid_edge_longitude'],
                                       self.read_instance.plotting_params[experiment]['grid_edge_latitude'])).T,
                                       edgecolor=self.read_instance.plotting_params[experiment]['colour'],
                                       transform=self.canvas_instance.datacrs,
                                       **self.canvas_instance.experiment_domain_polygon)
                # append polygon
                grid_edge_polygons.append(grid_edge_outline_poly)
            
        return grid_edge_polygons

    def make_header(self, pdf, plot_characteristics, plot_options=[]):
        """ Make header.
        
            :param plot_characteristics: plot characteristics 
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plot  
            :type plot_options: list
        """

        # set header title
        page = plt.figure(**plot_characteristics['figure'])
        plot_characteristics['page_title']['s'] = self.read_instance.report_title 
        plot_characteristics['page_title']['transform'] = page.transFigure
        page.text(**plot_characteristics['page_title'])

        # if len of network or species uniques is 1, set that instead of long list of duplicates
        _, idx = np.unique(self.read_instance.network, return_index=True)
        network_to_write = np.array(self.read_instance.network)[np.sort(idx)]
        _, idx = np.unique(self.read_instance.species, return_index=True)
        species_to_write = np.array(self.read_instance.species)[np.sort(idx)]

        # set header main text
        txt = 'Network : {}\n' \
              'Species : {}\n' \
              'Temporal Resolution : {}\n' \
              'Date Range : {} - {}\n' \
              'Experiments : {}\n' \
              'Temporal Colocation : {}\n' \
              'Spatial Colocation : {}\n' \
              .format(network_to_write,
                      species_to_write,
                      self.read_instance.resolution,
                      self.read_instance.start_date,
                      self.read_instance.end_date,
                      list(self.read_instance.experiments.values()),
                      self.read_instance.temporal_colocation,
                      self.read_instance.spatial_colocation 
                      )

        # add filter species to header if have it set
        if self.read_instance.filter_species: 
            txt += 'Filter Species : {}\n'.format(self.read_instance.filter_species)

        # add calibration factor to header if have it set
        if hasattr(self.read_instance, 'calibration_factor'):
            if self.read_instance.calibration_factor: 
                txt += 'Calibration : {}\n'.format(self.read_instance.calibration_factor)

        # add subsections to header
        txt += 'Subsections : {}\n'.format(self.read_instance.subsections)

        plot_characteristics['page_text']['s'] = txt   
        plot_characteristics['page_text']['transform'] = page.transFigure
        page.text(**plot_characteristics['page_text'])

        pdf.savefig(page, dpi=self.canvas_instance.dpi)
        plt.close(page)

    def make_metadata(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options=[],
                      station_inds=[]):
        """ Make metadata summary plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: name of data arrays to plot
            :type data_labels: str
            :param plot_characteristics: plot characteristics 
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plot  
            :type plot_options: list
            :param station_inds: list of relevant station indices
            :type station_inds: list
        """

        # initialise string to plot on axis
        str_to_plot = ''
        
        # get relevant station indices
        if not self.read_instance.offline:
            station_inds = self.canvas_instance.relative_selected_station_inds

        # setup some variables for handling if just one or multiple stations selected
        if len(station_inds) == 1:
            var_str_name = 'name_one' 
            str_to_plot += 'Station: {}'.format(self.read_instance.station_references[networkspeci][station_inds][0])
        else:
            var_str_name = 'name_multiple' 
            str_to_plot += '{} Stations'.format(len(station_inds))
        
        # iterate n vars per line and add spacing
        current_n_vars_per_line = 1

        # set up read for GHOST and non-GHOST data types 
        if self.read_instance.reading_ghost:
            char_var = 'ghost_vars'
        else:
            char_var = 'non-ghost_vars'
    
        # iterate through defined variables and add them
        for ghost_var, ghost_var_dict in plot_characteristics[char_var].items():

            # check var str name exists in ghost var dict, if not move to next var
            if var_str_name not in ghost_var_dict:
                continue

            # if are on limit of vars allowed per line then break to new line
            if current_n_vars_per_line == plot_characteristics['max_vars_per_row']:
                str_to_plot += '\n'
                current_n_vars_per_line = 0
            # otherwise, add spacing between variables on line
            else:
                str_to_plot += (' '*plot_characteristics['var_spacing'])

            # round decimal places if float
            if 'dp' in ghost_var_dict:
                str_to_plot += '{}: {:.{}f}'.format(ghost_var_dict[var_str_name],
                                                    np.nanmedian(self.read_instance.metadata_in_memory[networkspeci][ghost_var][
                                                    station_inds].astype(np.float32)), 
                                                    ghost_var_dict['dp'])

            # if str then get unique elements or percentage dependent on n uniques
            else:
                # gather all selected station metadata for current meta variable
                all_current_meta = self.read_instance.metadata_in_memory[networkspeci][ghost_var][
                    station_inds].flatten().astype(np.str)

                # get counts of all unique metadata elements across selected stations
                unique_meta, meta_counts = np.unique(all_current_meta, return_counts=True)
                # get number of unique metadata elements across selected stations
                n_unique_meta = len(unique_meta)

                # 1 unique metadata element? then return it
                if n_unique_meta == 1:
                    str_to_plot += '{}: {}'.format(ghost_var_dict[var_str_name], unique_meta[0])
                # if have > N unique metadata elements, just return count of the elements across the selected stations
                elif n_unique_meta > plot_characteristics['max_uniques']:
                    str_to_plot += '{}: {} uniques'.format(ghost_var_dict[var_str_name], n_unique_meta)
                # otherwise, get percentage of unique metadata elements across selected stations
                else:
                    meta_pc = (100. / len(all_current_meta)) * meta_counts
                    meta_pc = ['{:.1f}%'.format(meta) for meta in meta_pc]
                    # create string for variable to plot
                    str_to_plot += '{}: {}'.format(ghost_var_dict[var_str_name], ', '.join(
                        [':'.join([str(var), pc]) for var, pc in zip(unique_meta, meta_pc)]))

            # add units
            if 'units' in ghost_var_dict:
                if ghost_var_dict['units'] != '':
                    str_to_plot += ' {}'.format(ghost_var_dict['units'])

            # iterate current_n_vars_per_line
            current_n_vars_per_line += 1

        # plot string to axis
        plot_txt = relevant_axis.text(0.0, 1.0, str_to_plot, transform=relevant_axis.transAxes, **plot_characteristics['plot'])

        # modify limit to wrap text as axis width in pixels
        if not self.read_instance.offline:
            # get axis bounding box
            ax_bbox = relevant_axis.get_window_extent().transformed(self.canvas_instance.figure.dpi_scale_trans.inverted())
            
            # get axis dimensions in inches
            ax_width_inches = ax_bbox.width

            # get axis dimensions in pixels
            ax_width_px = ax_width_inches * self.canvas_instance.figure.dpi

        else:
            # get axis dimensions in pixels
            ax_width_px = relevant_axis.bbox.width * plot_characteristics['figure']['nrows']

        # automatically sets limit as figure width
        plot_txt._get_wrap_line_width = lambda: ax_width_px

        # track plot elements if using dashboard 
        if not self.read_instance.offline:
            self.track_plot_elements('observations', 'metadata', 'plot', plot_txt, bias=False)

    def make_map(self, relevant_axis, networkspeci, z_statistic, plot_characteristics, plot_options=[]):
        """ Make map plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param z_statistic: calculated z statistic to plot
            :type z_statistic: np.array
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plot  
            :type plot_options: list
        """

        # get marker size (for offline)
        if self.read_instance.offline:
            self.get_markersize(relevant_axis, 'map', networkspeci, plot_characteristics)

        # plot new station points on map - coloured by currently active z statisitic
        self.stations_scatter = relevant_axis.scatter(self.read_instance.station_longitudes[networkspeci][self.canvas_instance.active_map_valid_station_inds], 
                                                      self.read_instance.station_latitudes[networkspeci][self.canvas_instance.active_map_valid_station_inds], 
                                                      c=z_statistic, transform=self.canvas_instance.datacrs,
                                                      **plot_characteristics['plot'])
        
        # track plot elements if using dashboard 
        if not self.read_instance.offline:
            self.track_plot_elements('observations', 'map', 'plot', [self.stations_scatter], bias=False)

    def make_timeseries(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options=[]):
        """ Make timeseries plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: name of data arrays to plot
            :type data_labels: list
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plot  
            :type plot_options: list
        """

        # if 'obs' in plot_options, set data labels to just 'observations'
        if 'obs' in plot_options:
            data_labels = ['observations']

         # plot horizontal line across x axis at 0 if bias plot
        if 'bias' in plot_options:
            bias =  True
            bias_line = relevant_axis.axhline(**plot_characteristics['bias_line'])
            # track plot elements if using dashboard 
            if not self.read_instance.offline:
                self.track_plot_elements('ALL', 'timeseries', 'bias_line', [bias_line], bias=bias)
        else:
            bias = False

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # iterate through data labels
        for data_label in cut_data_labels:

            # bias plot?
            if bias:

                # skip if 'observations' data label
                if data_label == 'observations':
                    continue

                bias = True
                ts_obs = self.canvas_instance.selected_station_data[networkspeci]['timeseries']['observations']
                ts_model = self.canvas_instance.selected_station_data[networkspeci]['timeseries'][data_label] 
                ts = ts_model - ts_obs

            else:
                bias = False
                ts = self.canvas_instance.selected_station_data[networkspeci]['timeseries'][data_label]

            # get marker size (for offline)
            if self.read_instance.offline:
                self.get_markersize(relevant_axis, 'timeseries', networkspeci, plot_characteristics, data=ts)

            # make timeseries plot
            self.timeseries_plot = relevant_axis.plot(ts, 
                                                      color=self.read_instance.plotting_params[data_label]['colour'], 
                                                      **plot_characteristics['plot'])

            # update maximum smooth value
            if not self.read_instance.offline:
                if self.canvas_instance.timeseries_smooth_window_sl.value() != (len(ts)*2 - 1):
                    self.canvas_instance.timeseries_smooth_window_sl.setMaximum(len(ts)*2 - 1)

            # track plot elements if using dashboard 
            if not self.read_instance.offline:
                self.track_plot_elements(data_label, 'timeseries', 'plot', self.timeseries_plot, bias=bias)

    def make_periodic(self, relevant_axis, networkspeci, data_labels, plot_characteristics, 
                      zstat=None, plot_options=[]):
        """ Make period or period-violin plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: name of data arrays to plot
            :type data_labels: list
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param zstat: name of statistic
            :type zstat: str
            :param plot_options: list of options to configure plot  
            :type plot_options: list
        """

        # if 'obs' in plot_options, set data labels to just 'observations'
        if 'obs' in plot_options:
            data_labels = ['observations']

        # determine if 'bias' in plot_options
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # get zstat information
        if zstat:
            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat)

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # hide non-relevant resolution axes
        for nonrelevant_temporal_resolution in self.read_instance.nonrelevant_temporal_resolutions:
            # get subplot axis
            relevant_sub_ax = relevant_axis[nonrelevant_temporal_resolution]
            relevant_sub_ax.axis('off')
            relevant_sub_ax.set_visible(False)

        # iterate through all relevant temporal aggregation resolutions
        for relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions:

            # get subplot axis
            relevant_sub_ax = relevant_axis[relevant_temporal_resolution]

            # un-hide relevant resolution axis
            relevant_sub_ax.axis('on')
            relevant_sub_ax.set_visible(True)

            # set xticks
            xticks = self.canvas_instance.selected_station_data[networkspeci][relevant_temporal_resolution]['valid_xticks']

            # violin plot type?
            if not zstat:

                # calculate medians
                medians = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, 'Median', 
                                              data_labels, [], period=relevant_temporal_resolution)

                # calculate PDF for data label
                period_x_grid, PDFs_sampled = self.make_distribution(relevant_axis, networkspeci, data_labels, 
                                                                     plot_characteristics, 
                                                                     violin_resolution=relevant_temporal_resolution)

                # iterate through data labels and plot violins
                for data_label_ii, data_label in enumerate(cut_data_labels): 

                    # list to save all violins per data label
                    violins = []

                    # get median zorder
                    median_zorder = self.read_instance.plotting_params[data_label]['zorder']+len(cut_data_labels)
                    
                    # get alpha
                    if data_label == 'observations':
                        alpha = plot_characteristics['patch']['alpha_obs']
                    else:
                        alpha = plot_characteristics['patch']['alpha_exp']

                    # make plot of median
                    median_plots = relevant_sub_ax.plot(xticks, medians[:, data_label_ii], 
                                                        color=self.read_instance.plotting_params[data_label]['colour'], 
                                                        zorder=median_zorder, 
                                                        **plot_characteristics['plot']['median'])

                    # make violin plot
                    for period_ii in range(len(xticks)):

                        # get x_grid for period
                        x_grid = period_x_grid[period_ii]

                        if relevant_temporal_resolution == 'month':
                            period_pos = period_ii + 1
                        else:
                            period_pos = period_ii
                        PDF_sampled = PDFs_sampled[data_label_ii, period_ii, :]
                        if not np.all(np.isnan(PDF_sampled)):
                            
                            PDF_sampled = 0.5 * plot_characteristics['plot']['violin']['widths'] * PDF_sampled / PDF_sampled.max()
                            self.violin_plot = relevant_sub_ax.fill_betweenx(x_grid, -PDF_sampled + period_pos, PDF_sampled + period_pos,
                                                                             facecolor=self.read_instance.plotting_params[data_label]['colour'], 
                                                                             alpha=alpha)

                            # if have at least 1 valid experiment data array, split the violin plot across the horizontal
                            # (observations on left, experiment violin_plots on right)
                            if (len(cut_data_labels) > 1) & ('individual' not in plot_options):
                                m = np.mean(self.violin_plot.get_paths()[0].vertices[:, 0])
                                # observations on left
                                if data_label == 'observations':
                                    self.violin_plot.get_paths()[0].vertices[:, 0] = np.clip(self.violin_plot.get_paths()[0].vertices[:, 0], -np.inf, m)
                                # experiments on right
                                else:
                                    self.violin_plot.get_paths()[0].vertices[:, 0] = np.clip(self.violin_plot.get_paths()[0].vertices[:, 0], m, np.inf)

                            # save violin
                            violins.append(self.violin_plot)
                    
                    # add hidden poins for data limits per data label to allow for limit harmonisation
                    if len(violins) > 0:
                        limit_plot = relevant_sub_ax.plot([period_pos,period_pos], 
                                                          [np.min(period_x_grid),np.max(period_x_grid)],
                                                           alpha=0.0)
                        violins.extend(limit_plot)

                    # track plot elements if using dashboard 
                    if not self.read_instance.offline:
                        self.track_plot_elements(data_label, 'periodic-violin', 'violin_plot_{}'.format(relevant_temporal_resolution), violins, bias=False)
                        self.track_plot_elements(data_label, 'periodic-violin', 'Median_plot_{}'.format(relevant_temporal_resolution), median_plots, bias=False)

            # periodic plot type
            else:

                # plot horizontal line/s across x axis at value/s of minimum experiment bias (if bias statistic)
                if z_statistic_sign == 'bias':
                    # get value/s of minimum bias for statistic
                    if z_statistic_type == 'basic':
                        minimum_bias = basic_stats[base_zstat]['minimum_bias']
                    else:
                        minimum_bias = expbias_stats[base_zstat]['minimum_bias']
                    bias_lines = []
                    for mb in minimum_bias:
                        bias_lines += [relevant_sub_ax.axhline(y=mb, **plot_characteristics['bias_line'])]
                    # track plot elements if using dashboard 
                    if not self.read_instance.offline:
                        self.track_plot_elements('ALL', 'periodic', 'bias_line_{}'.format(relevant_temporal_resolution), 
                                                 bias_lines, bias=bias)

                # calculate statistic in each of periodic groups per data label
                if z_statistic_sign == 'bias':
                    if 'observations' in cut_data_labels:
                        cut_data_labels.remove('observations')
                    stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstat, 
                                                    ['observations']*len(cut_data_labels), cut_data_labels, 
                                                    period=relevant_temporal_resolution)
                else:
                    stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, 
                                                     zstat, data_labels, [], period=relevant_temporal_resolution)

                # iterate through data labels and plot violins
                for data_label_ii, data_label in enumerate(cut_data_labels): 

                    # skip observational array if bias stat
                    if (z_statistic_sign == 'bias') & (data_label == 'observations'):
                        continue

                    # make plot
                    self.periodic_plots = relevant_sub_ax.plot(xticks, stats_calc[:, data_label_ii], 
                                                               color=self.read_instance.plotting_params[data_label]['colour'], 
                                                               zorder=self.read_instance.plotting_params[data_label]['zorder'], 
                                                               **plot_characteristics['plot'])

                    # track plot elements if using dashboard 
                    if not self.read_instance.offline:
                        self.track_plot_elements(data_label, 'periodic', 'plot_{}'.format(relevant_temporal_resolution), self.periodic_plots, bias=bias)

    def make_distribution(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options=[],
                          data_range_min=None, data_range_max=None, violin_resolution=None):
        """ Make distribution plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: name of data arrays to plot
            :type data_labels: list
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plot  
            :type plot_options: list
            :param data_range_min: minimum data range of distribution plot grid 
            :type data_range_min: float
            :param data_range_max: maximum data range of distribution plot grid 
            :type data_range_max: float
            :param violin_resolution: if are calculating distribution for violin plot, this is set to temporal resolution of groupings
            :type violin_resolution: int
        """

        # determine if 'bias' in plot_options
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # if 'obs' in plot_options, set data labels to just 'observations'
        if 'obs' in plot_options:
            data_labels = ['observations']

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # set data ranges for distribution plot grid if not set explicitly
        if not data_range_min:
            data_range_min = self.canvas_instance.selected_station_data_min[networkspeci]

        if not data_range_max:
            data_range_max = self.canvas_instance.selected_station_data_max[networkspeci]

        # set xgrid for calculating distribution
        # if calculating for period n_samples is set to pdf_min_samples
        # otherwise it is inferred from data (if above minimum value)
        if violin_resolution:
            n_samples = plot_characteristics['pdf_min_samples']
        else:
            minmax_diff = data_range_max - data_range_min
            if pd.isnull(self.read_instance.parameter_dictionary[networkspeci.split('|')[1]]['minimum_resolution']):
                n_samples = plot_characteristics['pdf_min_samples']
            else:
                n_samples = int(np.around(minmax_diff/(self.read_instance.parameter_dictionary[networkspeci.split('|')[1]]['minimum_resolution']/100.0),0))
                if n_samples < plot_characteristics['pdf_min_samples']:
                    n_samples = plot_characteristics['pdf_min_samples']

        # round n_samples to next next power of 2 (for fft optimisation)
        n_samples = 2 ** np.ceil(np.log2(n_samples)) 
        # set x_grid
        x_grid = np.linspace(data_range_min,data_range_max,int(n_samples))

        # plot horizontal line across x axis at 0 if bias plot
        # also remove observations from cut_data_labels
        if bias:
            bias_line = [relevant_axis.axhline(**plot_characteristics['bias_line'])]
            # track plot elements if using dashboard 
            if not self.read_instance.offline:
                self.track_plot_elements('ALL', 'distribution', 'bias_line', bias_line, bias=bias)
            if 'observations' in cut_data_labels:
                cut_data_labels.remove('observations')

        # if violin plot setup arrays for saving data
        if violin_resolution:
            PDFs_sampled = np.full((len(cut_data_labels), len(self.canvas_instance.periodic_xticks[violin_resolution]), int(n_samples)), np.NaN, dtype=np.float32)

        # iterate through data labels
        for data_label_ii, data_label in enumerate(cut_data_labels):
        
            PDF_sampled_calculated = False

            # setup bias plot
            if bias:

                # calculate obs PDF on first pass
                if data_label_ii == 0:
                    kde_data_obs = drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index('observations'),0,:])
                    
                    #filter out data outside data range bounds
                    #kde_data_obs = kde_data_obs[(kde_data_obs > data_range_min) & (kde_data_obs < data_range_max)]
                    
                    # check if all values are equal in the dataframe
                    if kde_data_obs.size == 0:
                        print('Warning: The kernel density cannot be calculated because there are no valid observational values.')
                        return
                    elif np.all(kde_data_obs == kde_data_obs[0]):
                        print('Warning: The kernel density cannot be calculated because all observational values are equal.')
                        return
                    else:
                        PDF_obs_sampled = kde_fft(kde_data_obs, xgrid=x_grid)
                        #PDF_fit = FFTKDE(kernel='gaussian', bw='scott').fit(kde_data_obs)
                        #PDF_obs_sampled = PDF_fit.evaluate(x_grid)


                # calculate model PDF
                kde_data_model = drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:])
                
                #filter out data outside data range bounds
                #kde_data_model = kde_data_model[(kde_data_model > data_range_min) & (kde_data_model < data_range_max)]
                
                # check if all values are equal in the dataframe
                if kde_data_model.size == 0:
                    print('Warning: The kernel density cannot be calculated because there are no valid values for {} experiment.'.format(data_label))
                    continue
                elif np.all(kde_data_model == kde_data_model[0]):
                    print('Warning: The kernel density cannot be calculated because all values for {} experiment are equal.'.format(data_label))
                    continue
                # calculate PDF
                PDF_model_sampled = kde_fft(kde_data_model, xgrid=x_grid)
                #PDF_fit = FFTKDE(kernel='gaussian', bw='scott').fit(kde_data_model)
                #PDF_model_sampled = PDF_fit.evaluate(x_grid)

                # calculate difference
                PDF_sampled = PDF_model_sampled - PDF_obs_sampled
                PDF_sampled_calculated = True

            # setup standard plot
            else:

                # if first data label and calculating distributions for violin plot,
                # calculate the x_grid / data ranges per period  
                # use min for min data range and upper inner Tukey fence for max data range
                if (violin_resolution != None) & (data_label_ii == 0):
                    period_data_range_min = []
                    period_data_range_max = []
                    period_x_grid = []
                    # iterate through periods
                    for group in self.canvas_instance.selected_station_data[networkspeci][violin_resolution]['active_mode']:
                        lower_inner_fence, upper_inner_fence = boxplot_inner_fences(group)
                        min_data = np.nanmin(group)
                        period_data_range_min.append(min_data)
                        period_data_range_max.append(upper_inner_fence)
                        period_x_grid.append(np.linspace(min_data,upper_inner_fence,int(n_samples)))

                # get data (flattened and drop NaNs)
                if violin_resolution:
                    kde_data_grouped = [drop_nans(group[valid_data_labels.index(data_label)].flatten()) for group in self.canvas_instance.selected_station_data[networkspeci][violin_resolution]['active_mode']]
                else:    
                    kde_data_grouped = [drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:])]

                # iterate through kde data groups
                for period_ii, kde_data in enumerate(kde_data_grouped):

                    # get relevant data ranges / x_grid, for violin period distribution calculation
                    if violin_resolution:
                        data_range_min = period_data_range_min[period_ii]
                        data_range_max = period_data_range_max[period_ii]
                        x_grid = period_x_grid[period_ii]

                    #filter out data outside data range bounds
                    #kde_data = kde_data[(kde_data > data_range_min) & (kde_data < data_range_max)]

                    # check if all values are equal in the dataframe
                    if kde_data.size == 0:
                        if not violin_resolution:
                            print('Warning: The kernel density cannot be calculated because there are no valid values for {}.'.format(data_label))
                        continue
                    elif np.all(kde_data == kde_data[0]):
                        if not violin_resolution:
                            print('Warning: The kernel density cannot be calculated because all {} values are equal.'.format(data_label))
                        continue
                    else:
                        PDF_sampled = kde_fft(kde_data, xgrid=x_grid)
                        #PDF_fit = FFTKDE(kernel='gaussian', bw='scott').fit(kde_data)
                        #PDF_sampled = PDF_fit.evaluate(x_grid)
                        # save PDF for violin plot
                        if violin_resolution:
                            PDFs_sampled[data_label_ii, period_ii, :] = PDF_sampled
                        else:
                            PDF_sampled_calculated = True

            # make plot (if not making PDF for violin plot)
            if PDF_sampled_calculated:
                # make plot
                self.distribution_plot = relevant_axis.plot(x_grid, PDF_sampled, 
                                                            color=self.read_instance.plotting_params[data_label]['colour'], 
                                                            **plot_characteristics['plot'])

                # track plot elements if using dashboard 
                if not self.read_instance.offline:
                    self.track_plot_elements(data_label, 'distribution', 'plot', self.distribution_plot, bias=bias)

        # if have made PDFs for violin plot then return it
        if violin_resolution:
            return period_x_grid, PDFs_sampled

    def make_scatter(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options=[]):
        """ Make scatter plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: name of data arrays to plot
            :type data_labels: list
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plot  
            :type plot_options: list
        """

        # if 'obs' in plot_options, set data labels to just 'observations'
        if 'obs' in plot_options:
            data_labels = ['observations']

        # add 1:1 line (if in plot_characteristics)
        if '1:1_line' in plot_characteristics:
            relevant_axis.plot([0, 1], [0, 1], transform=relevant_axis.transAxes, 
                            **plot_characteristics['1:1_line'])
        # add 1:2 line (if in plot_characteristics)
        if '1:2_line' in plot_characteristics:
            relevant_axis.plot([0, 1], [0, 0.5], transform=relevant_axis.transAxes, 
                            **plot_characteristics['1:2_line'])     
        # add 2:1 line (if in plot_characteristics)
        if '2:1_line' in plot_characteristics:
            relevant_axis.plot([0, 0.5], [0, 1], transform=relevant_axis.transAxes, 
                            **plot_characteristics['2:1_line'])

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get observations data (flattened)
        observations_data = self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index('observations'),0,:]

        # determine if number of points per data array exceeds max limit,
        # if so subset arrays
        subset = False
        data_array_size = observations_data.size
        if data_array_size > plot_characteristics['max_points']:
            subset = True
            inds_subset = np.random.choice(data_array_size, size=plot_characteristics['max_points'], replace=False)
            observations_data = observations_data[inds_subset]

        # iterate through data labels
        for data_label in cut_data_labels:

            # continue for observations data label
            if data_label == 'observations':
                continue

            # get experiment data (flattened)
            experiment_data = self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:]

            # subset data if neccessary
            if subset:
                experiment_data = experiment_data[inds_subset]

            # get marker size (for offline)
            if self.read_instance.offline:
                self.get_markersize(relevant_axis, 'scatter', networkspeci, plot_characteristics, data=observations_data)

            # create scatter plot
            self.scatter_plot = relevant_axis.plot(observations_data, experiment_data, 
                                                   color=self.read_instance.plotting_params[data_label]['colour'],
                                                   **plot_characteristics['plot'])

            # track plot elements if using dashboard 
            if not self.read_instance.offline:
                self.track_plot_elements(data_label, 'scatter', 'plot', self.scatter_plot, bias=False)

    def make_boxplot(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options=[]):
        """ Make boxplot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: name of data arrays to plot
            :type data_labels: list
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plot  
            :type plot_options: list
        """

        # if 'obs' in plot_options, set data labels to just 'observations'
        if 'obs' in plot_options:
            data_labels = ['observations']

        # if multispecies in plot options then make plot for all networkspecies
        if 'multispecies' in plot_options:
            networkspecies = self.read_instance.networkspecies
        else:
            networkspecies = [networkspeci]

        # iterate through networkspecies
        ns_current = 0
        for ns in networkspecies:

            # get valid data labels for networkspeci
            valid_data_labels = self.canvas_instance.selected_station_data_labels[ns]

            # cut data_labels for those in valid data labels
            cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

            # only proceed if have some data labels to plot
            if len(cut_data_labels) > 0:

                # get data label width and spacing
                if ('individual' in plot_options) or ('obs' in plot_options) or (len(self.read_instance.networkspecies) == 1):
                    widths = plot_characteristics['group_widths']['singlespecies']
                else:
                    available_width = plot_characteristics['group_widths']['multispecies']
                    remainder_width = 1.0 - available_width
                    start_point = -0.5 + (remainder_width / 2.0)
                    widths = available_width / (len(cut_data_labels) + 0.15)
                    spacing = (available_width - (widths * len(cut_data_labels))) / (len(cut_data_labels) - 1)

                # get plot positions
                if ('individual' in plot_options) or ('obs' in plot_options):
                    positions = [ns_current]
                elif (len(self.read_instance.networkspecies) == 1):
                    positions = np.arange(len(cut_data_labels))
                else:
                    positions = [((start_point + (widths/2.0)) + (spacing * data_label_ii) + (widths * data_label_ii)) + ns_current 
                                for data_label_ii in range(len(cut_data_labels))]  
                
                # iterate ns_current
                ns_current += 1

                # iterate through cut data labels making plot
                for data_label_ii, data_label in enumerate(cut_data_labels):

                    # get data (flattened and drop NaNs)
                    data_array = drop_nans(self.canvas_instance.selected_station_data[ns]['flat'][valid_data_labels.index(data_label),0,:])

                    # make boxplot
                    boxplot = relevant_axis.boxplot(x=data_array, positions=[positions[data_label_ii]], widths=widths, 
                                                    **plot_characteristics['plot'])

                    # set box colour
                    for element in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
                        plt.setp(boxplot[element], color=self.read_instance.plotting_params[data_label]['colour'])
                    # set fill colour to be white
                    for patch in boxplot['boxes']:
                        patch.set(facecolor='white')

                    # track plot elements if using dashboard 
                    if not self.read_instance.offline:
                        self.track_plot_elements(data_label, 'boxplot', 'plot', boxplot, bias=False)

        # set xticklabels 
        # labels for multispecies plot
        xtick_params = copy.deepcopy(plot_characteristics['xtick_params'])
        xticklabel_params = copy.deepcopy(plot_characteristics['xticklabel_params'])
        if ('multispecies' in plot_options) & (len(self.read_instance.networkspecies) > 1):
            xticks = np.arange(len(self.read_instance.networkspecies))
            #if all networks or species are same, drop them from xtick label
            if len(np.unique(self.read_instance.network)) == 1:
                xtick_labels = copy.deepcopy(self.read_instance.species)
            elif len(np.unique(self.read_instance.species)) == 1:
                xtick_labels = copy.deepcopy(self.read_instance.network)
            else:
                xtick_labels = copy.deepcopy(self.read_instance.networkspecies)
            # get aliases for multispecies (if have any)
            xtick_labels, xlabel = get_multispecies_aliases(xtick_labels)
            # set xlabel if xlabels have changed due to alias, and have one to set
            if xlabel != '':
                plot_characteristics['xlabel']['xlabel'] = xlabel
                relevant_axis.set_xlabel(**plot_characteristics['xlabel'])
        # labels for standard plot
        else:
            data_labels_to_plot = copy.deepcopy(cut_data_labels)
            for data_label_ii, data_label in enumerate(cut_data_labels):
                if data_label == 'observations':
                    if 'legend' in plot_characteristics:
                        obs_label = plot_characteristics['legend']['handles']['obs_label']
                    elif 'legend' in self.canvas_instance.plot_characteristics_templates.keys():
                        obs_label = self.canvas_instance.plot_characteristics_templates['legend']['handles']['obs_label']
                    else:
                        obs_label = 'Observations'
                    data_labels_to_plot[data_label_ii] = obs_label
                else:
                    data_labels_to_plot[data_label_ii] = self.read_instance.experiments[data_label]
            xticks = positions
            if ('individual' in plot_options) or ('obs' in plot_options):
                xtick_labels = [data_labels_to_plot[cut_data_labels.index(data_label)]]
                relevant_axis.xaxis.set_tick_params()
            else:
                xtick_labels = data_labels_to_plot

        #modify xticks to be horizontal as just have 1 label
        if len(xtick_labels) == 1:
            xtick_params['rotation'] = 0
            xticklabel_params = {}

        # set xticks / xticklabels
        relevant_axis.set_xticks(xticks)
        relevant_axis.xaxis.set_tick_params(**xtick_params)
        relevant_axis.set_xticklabels(xtick_labels, **xticklabel_params)

    def make_heatmap(self, relevant_axis, networkspeci, data_labels, plot_characteristics, 
                     zstats=None, plot_options=[], subsection=None, plotting_paradigm=None, stats_df=None):
        """ Make heatmap plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: name of data arrays to plot
            :type data_labels: list
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param zstats: name of statistics
            :type zstats: list
            :param plot_options: list of options to configure plot  
            :type plot_options: list
            :param subsection: str of currently active subsection
            :type subsection: str
            :param plotting_paradigm: plotting paradigm (summary or station in offline reports)
            :type plotting_paradigm: str
            :param stats_df: dataframe of previously calculated statistics, default is None
            :param stats_df: pandas dataframe
        """

        # determine if want to add annotations or not from plot_options
        if 'annotate' in plot_options:
            annotate = True
        else:
            annotate = False

        # bias plot?
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # if statistical dataframe is not provided then create it
        if not isinstance(stats_df, pd.DataFrame):

            # get valid data labels for networkspeci
            valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

            # cut data_labels for those in valid data labels
            cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

            # calculate statistics
            if bias:
                if 'observations' in cut_data_labels:
                    cut_data_labels.remove('observations')
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstats, 
                                                ['observations']*len(cut_data_labels), cut_data_labels)
            else:
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstats, 
                                                 cut_data_labels, [])

            # create stats dataframe
            stats_df = pd.DataFrame(data=stats_calc, 
                                    index=cut_data_labels,
                                    dtype=np.float)

        # get observations label
        if 'legend' in plot_characteristics:
            obs_label = plot_characteristics['legend']['handles']['obs_label'] 
        elif 'legend' in self.canvas_instance.plot_characteristics_templates.keys():
            obs_label = self.canvas_instance.plot_characteristics_templates['legend']['handles']['obs_label'] 
        else:
            obs_label = 'Observations'

        # round dataframe
        stats_df = stats_df.round(plot_characteristics['round_decimal_places']['table'])

        # get subsections
        subsections = list(np.unique(stats_df.index.get_level_values(1)))

        # get relevant data
        if plotting_paradigm == 'station':
            stats_df = stats_df.iloc[stats_df.index.get_level_values('subsections') == subsection]
        if 'multispecies' not in plot_options:
            stats_df = stats_df.iloc[stats_df.index.get_level_values('networkspecies') == networkspeci]
        else:
            # replace subsection name by networkspecies if there is only one
            if (len(subsections) == 1) or (plotting_paradigm == 'station'):
                stats_df = stats_df.droplevel(level='subsections')

        # plot heatmap
        heatmap = sns.heatmap(stats_df, 
                              ax=relevant_axis, 
                              annot=annotate,
                              **plot_characteristics['plot'])

        # remove networkspecies-subsections label from y-axis
        relevant_axis.set_ylabel("")

        # if there is only one subsection or station data
        if (plotting_paradigm == 'station') or (len(subsections) == 1):
            # for multispecies, remove network names from labels
            if ('multispecies' in plot_options) and (not plot_characteristics['multispecies']['network_names']):
                if not plot_characteristics['multispecies']['network_names']:
                    yticklabels = [networkspeci.split('|')[1] if '|' in networkspeci else networkspeci 
                                   for networkspeci in stats_df.index]
            # for non multispecies, remove subsection names from labels
            elif ('multispecies' not in plot_options) and (not plot_characteristics['parent_section_names']):
                yticklabels = []
                for subsection_label in stats_df.index.get_level_values(1):
                    if "·" in subsection_label:
                        subsection_label = subsection_label.split('·')[1]
                    yticklabels.append(subsection_label)
            # keep original labels
            else:
                yticklabels = stats_df.index
        # if there is summary data for more than one subsection
        elif (plotting_paradigm == 'summary') and (len(subsections) > 1):
            # remove parent names from subsections
            if not plot_characteristics['parent_section_names']:
                yticklabels = []
                for subsection_label in stats_df.index.get_level_values(1):
                    if "·" in subsection_label:
                        subsection_label = subsection_label.split('·')[1]
                    yticklabels.append(subsection_label)
            # keep original labels
            else:
                yticklabels = stats_df.index.get_level_values(1)
        relevant_axis.set_yticklabels(yticklabels, **plot_characteristics['yticklabels'])

        # update observation and experiments labels
        xticklabels = stats_df.columns
        xticklabels = [obs_label if item == 'observations' else item if item == '' else self.read_instance.experiments[item] for item in xticklabels]
        relevant_axis.set_xticklabels(xticklabels, **plot_characteristics['xticklabels'])

        # format for multispecies
        if 'multispecies' in plot_options:
            # if we have more than one subsection and we are plotting summaries
            if (len(subsections) > 1) and (plotting_paradigm == 'summary'):
                # add horizontal lines to separate networkspecies
                networkspecies = list(stats_df.index.get_level_values(0)[::(len(subsections))])
                y_separators = []
                for networkspeci_ii in range(len(networkspecies)+1):
                    y_separators.append(len(subsections)*networkspeci_ii)
                relevant_axis.hlines(y=y_separators, xmin=plot_characteristics['multispecies']['xmin'], 
                                     xmax=0, clip_on=False, **plot_characteristics['multispecies']['hlines'])

                # annotate networkspecies names on the left
                for networkspeci, y_separator in zip(networkspecies, y_separators[:-1]):
                    y_position = (len(subsections)/2) + y_separator
                    # remove network names from networkspecies
                    if not plot_characteristics['multispecies']['network_names']:
                        networkspeci_label = networkspeci.split('|')[1]
                    else:
                        networkspeci_label = networkspeci
                    relevant_axis.annotate(s=networkspeci_label, annotation_clip=False,
                                           xy=(plot_characteristics['multispecies']['xmin'], y_position), 
                                           **plot_characteristics['multispecies']['yticklabels'])
        
        # format for non multispecies
        else:
            # axis cuts off due to bug in matplotlib 3.1.1 - hack fix. Remove in Future!
            # if len(stats_df.index) > 1:
            #     bottom, top = relevant_axis.get_ylim()
            #     relevant_axis.set_ylim(bottom + 0.5, top - 0.5)

            # vertically align yticklabels due to bug again in matplotlib - hack fix. Remove in Future!
            for tick in relevant_axis.get_yticklabels():
                tick.set_verticalalignment("center")

        # track plot elements if using dashboard 
        if not self.read_instance.offline:
            self.track_plot_elements('observations', 'heatmap', 'plot', heatmap, bias=bias)

    def merge_cells(self, table, cells):
        """
        Merge cells in matplotlib.Table. Reference: https://stackoverflow.com/a/53819765/12684122.

        :param table: table
        :type table: matplotlib.Table
        :param cells: cells to merge in table coordinates, e.g. [(0,1), (0,0), (0,2)]
        :type cells: list
        """

        cells_array = [np.asarray(c) for c in cells]
        h = np.array([cells_array[i+1][0] - cells_array[i][0] for i in range(len(cells_array) - 1)])
        v = np.array([cells_array[i+1][1] - cells_array[i][1] for i in range(len(cells_array) - 1)])

        # if it's a horizontal merge, all values for `h` are 0
        if not np.any(h):
            # sort by horizontal coord
            cells = np.array(sorted(list(cells), key=lambda v: v[1]))
            edges = ['BTL'] + ['BT' for i in range(len(cells) - 2)] + ['BTR']
        elif not np.any(v):
            cells = np.array(sorted(list(cells), key=lambda h: h[0]))
            edges = ['TRL'] + ['RL' for i in range(len(cells) - 2)] + ['BRL']
        else:
            raise ValueError("Only horizontal and vertical merges allowed")

        for cell, e in zip(cells, edges):
            table[cell[0], cell[1]].visible_edges = e
            
        txts = [table[cell[0], cell[1]].get_text() for cell in cells]
        tpos = [np.array(t.get_position()) for t in txts]

        # transpose the text of the left cell
        trans = (tpos[-1] - tpos[0])/2

        # didn't had to check for ha because I only want ha='center'
        txts[0].set_transform(mpl.transforms.Affine2D().translate(*trans))
        for txt in txts[1:]:
            txt.set_visible(False)

    def make_table(self, relevant_axis, networkspeci, data_labels, plot_characteristics,
                   zstats=None, statsummary=False, plot_options=[],
                   subsection=None, plotting_paradigm=None, stats_df=None):
        """ Make table plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: name of data arrays to plot
            :type data_labels: list
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param zstats: name of statistics
            :type zstats: list
            :param statsummary: boolean indiciating if making alternative statistical summary table plot  
            :type statsummary: boolean
            :param plot_options: list of options to configure plot  
            :type plot_options: list
            :param subsection: str of currently active subsection
            :type subsection: str
            :param plotting_paradigm: plotting paradigm (summary or station in offline reports)
            :type plotting_paradigm: str
            :param stats_df: dataframe of previously calculated statistics, default is None
            :param stats_df: pandas dataframe
        """

        # turn off axis to make table
        relevant_axis.axis('off')

        # bias plot?
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # if statistical dataframe is not provided then create it
        if not isinstance(stats_df, pd.DataFrame):

            # get valid data labels for networkspeci
            valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

            # cut data_labels for those in valid data labels
            cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

            # calculate statistics
            if bias:
                if 'observations' in cut_data_labels:
                    cut_data_labels.remove('observations')
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstats, 
                                                ['observations']*len(cut_data_labels), cut_data_labels)
            else:
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstats, 
                                                 cut_data_labels, [])

            # create stats dataframe
            stats_df = pd.DataFrame(data=stats_calc, 
                                    index=cut_data_labels,
                                    dtype=np.float)

        # rename columns to save space
        columns = {}
        for column in stats_df.columns:
            new_colum_name = copy.deepcopy(column)
            # rename columns to replace Diurnal, Weekly and Monthly by D, W, M
            if 'diurnal' in column:
                new_colum_name = column.replace('diurnal', 'D')
            elif 'weekly' in column:
                new_colum_name = column.replace('weekly', 'W')
            elif 'monthly' in column:
                new_colum_name = column.replace('monthly', 'M')
            # remove _bias from columns
            if '_bias' in new_colum_name:
                new_colum_name = new_colum_name.replace('_bias', '')   
            columns[column] = new_colum_name
        stats_df = stats_df.rename(columns=columns)
  
        # get column and row labels
        col_labels = stats_df.columns.tolist()
        row_labels = stats_df.index.tolist()

        # get observations label
        if 'legend' in plot_characteristics:
            obs_label = plot_characteristics['legend']['handles']['obs_label'] 
        elif 'legend' in self.canvas_instance.plot_characteristics_templates.keys():
            obs_label = self.canvas_instance.plot_characteristics_templates['legend']['handles']['obs_label'] 
        else:
            obs_label = 'Observations'

        # round dataframe
        stats_df = stats_df.round(plot_characteristics['round_decimal_places']['table'])

        # offline reports
        if self.read_instance.offline:
            
            # get relevant data
            if 'multispecies' not in plot_options:
                stats_df = stats_df.iloc[stats_df.index.get_level_values('networkspecies') == networkspeci]
            if plotting_paradigm == 'station':
                stats_df = stats_df.iloc[stats_df.index.get_level_values('subsections') == subsection]

            # get labels
            networkspecies = list(stats_df.index.get_level_values('networkspecies'))
            subsections = list(stats_df.index.get_level_values('subsections'))
            if statsummary:
                data_labels = list(stats_df.index.get_level_values('labels'))
                stats = list(stats_df.columns)
            else:
                data_labels = list(stats_df.columns)

            # reset index after filtering
            stats_df = stats_df.reset_index()

            # hide subsections from station plots or if there is only 1 section
            if self.read_instance.offline:
                if plotting_paradigm == 'station' or len(np.unique(subsections)) == 1:
                    stats_df = stats_df.drop(columns='subsections')
        
            # hide networkspecies from plots that are not multispecies
            if 'multispecies' not in plot_options:
                stats_df = stats_df.drop(columns='networkspecies')
        
            # remove parent names from subsections
            if ('subsections' in stats_df.columns) and (not plot_characteristics['parent_section_names']):
                stats_df['subsections'] = [subsection_label.split('·')[1] if '·' in subsection_label 
                                           else subsection_label for subsection_label in subsections]

            # remove network names from networkspecies
            if (('multispecies' in plot_options) and ('networkspecies' in stats_df.columns) and 
                (not plot_characteristics['multispecies']['network_names'])):
                stats_df['networkspecies'] = [networkspeci_label.split('|')[1] for networkspeci_label in networkspecies]
            
            # get number of "empty" cells (without stats) and 
            # column labels (hide networkspecies, subsections and data labels)
            if statsummary:
                empty_cells = len(stats_df.columns) - len(stats)
                col_labels = ['']*empty_cells + stats
                stats_df['labels'] = [obs_label if item == 'observations' else item if item == '' else self.read_instance.experiments[item] for item in stats_df['labels']]
            else:
                empty_cells = len(stats_df.columns) - len(data_labels)
                col_labels = ['']*empty_cells + data_labels
                col_labels = [obs_label if item == 'observations' else item if item == '' else self.read_instance.experiments[item] for item in col_labels]

        # dashboard
        else:
            # there is only statsummary
            if statsummary:

                # get labels
                data_labels = list(stats_df.index)
                stats = list(stats_df.columns)
                
                # reset index
                stats_df = stats_df.reset_index()
                
                # update observation and experiments labels
                stats_df['index'] = [obs_label if item == 'observations' else item if item == '' else self.read_instance.experiments[item] for item in stats_df['index']]

                # get number of "empty" cells (without stats)
                empty_cells = 1
                col_labels = ['']*empty_cells + stats

        # set cell colors
        if statsummary:
            if 'cell_colours' in plot_characteristics:
                if plot_characteristics['cell_colours']:
                    cell_colours = [[]] * (stats_df.shape[1])
                    for col in range(stats_df.shape[1]):
                        # custom colors for data labels cells
                        if col == (empty_cells-1):
                            for data_label in data_labels:
                                # observations in white
                                if data_label == 'observations':
                                    color = 'white'
                                # experiments in legend colors
                                else:
                                    color = self.read_instance.plotting_params[data_label]['colour']
                                cell_colours[col].append(color)
                        # white for other cells
                        else:
                            cell_colours[col] = ['white'] * stats_df.shape[0]
                    plot_characteristics['plot']['cellColours'] = np.array(cell_colours).T
        else:
            if 'col_colours' in plot_characteristics:
                if plot_characteristics['col_colours']:
                    col_colours = []
                    for data_label in data_labels:
                        # observations in white
                        if data_label == 'observations':
                            color = 'white'
                        # experiments in legend colors
                        else:
                            color = self.read_instance.plotting_params[data_label]['colour']
                        col_colours.extend([color])
                    plot_characteristics['plot']['colColours'] = ['white']*empty_cells + col_colours

        # make table
        table = relevant_axis.table(cellText=stats_df.values, 
                                    colLabels=col_labels,
                                    **plot_characteristics['plot'])

        # merge cells in networkspecies and subsections columns (if any)
        if self.read_instance.offline:
            column_ii = 0
            for column, rows in zip(['networkspecies', 'subsections'], (networkspecies, subsections)):
                if column in stats_df.columns:
                    # count consecutive duplicates
                    count_dups = [sum(1 for _ in group) for _, group in groupby(rows)]
                    
                    # merge cells that have consecutive duplicates
                    current_row = 0
                    for count_ii, count in enumerate(count_dups):
                        cells_to_merge = [(current_row + i, column_ii) for i in range(1, count+1)]
                        self.merge_cells(table, cells_to_merge)
                        current_row += count
                    column_ii += 1

        # adjust cell height
        if 'cell_height' in plot_characteristics:
            table.scale(1, plot_characteristics['cell_height'])

        # adjust fontsize
        if 'fontsize' in plot_characteristics:
            table.auto_set_font_size(False)
            table.set_fontsize(plot_characteristics['fontsize'])
            table.auto_set_column_width(np.arange(-1, len(col_labels)+1))

        # track plot elements if using dashboard 
        if not self.read_instance.offline:
            if statsummary:
                self.track_plot_elements('observations', 'statsummary', 'plot', [table], bias=bias)
            else:
                self.track_plot_elements('observations', 'table', 'plot', [table], bias=bias)
    
    def get_taylor_diagram_ghelper_info(self, reference_stddev, plot_characteristics, extend=False):
        """ Make Taylor diagram plot axis extremes and labels. 
            
            :param reference_stddev: reference standard deviation to set the limits
            :type reference_stddev: float
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param extend: extend to negative correlations
            :type undo: boolean
        """

        tmin = 0

        # diagram extended to negative correlations
        if extend:
            tmax = np.pi
        # diagram limited to positive correlations
        else:
            tmax = np.pi/2

        # get standard deviation axis extent
        srange = plot_characteristics['srange']
        smin = srange[0] * reference_stddev
        smax = srange[1] * reference_stddev

        # correlation labels
        if extend:
            rlocs = np.array(plot_characteristics['rlocs']['all'])
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            rlocs = np.array(plot_characteristics['rlocs']['positive'])

        # convert correlation values into polar angles
        tlocs = np.arccos(rlocs)
        gl1 = gf.FixedLocator(tlocs)
        tf1 = gf.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        return tmin, tmax, smin, smax, gl1, tf1
    
    def get_taylor_diagram_ghelper(self, reference_stddev, plot_characteristics, extend=False):
        """ Make Taylor diagram plot grid helper. 

            :param reference_stddev: reference standard deviation to set the limits
            :type reference_stddev: float
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param extend: extend to negative correlations
            :type undo: boolean
        """

        # get axis extremes
        tmin, tmax, smin, smax, gl1, tf1 = self.get_taylor_diagram_ghelper_info(reference_stddev, plot_characteristics, 
                                                                                extend)

        # get grid helper
        ghelper = fa.GridHelperCurveLinear(PolarAxes.PolarTransform(),
                                           extremes=(tmin, tmax, smin, smax),
                                           grid_locator1=gl1, tick_formatter1=tf1)

        return ghelper

    def make_taylor(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options=[], 
                    stddev_max=None):
        """ Make Taylor diagram plot.
            Reference: https://gist.github.com/ycopin/3342888.

            See explanation of calculations here:
            https://waterprogramming.wordpress.com/2020/12/22/taylor-diagram/

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: name of data arrays to plot
            :type data_labels: list
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plot  
            :type plot_options: list
            :param stddev_max: maximum standard deviation
            :type stddev_max: float
        """

        # get polar axis here
        self.taylor_polar_relevant_axis = relevant_axis.get_aux_axes(PolarAxes.PolarTransform())

        # calculate statistics
        stats_dict = {}

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get data labels without observations
        obs_index = cut_data_labels.index('observations')
        data_labels_sans_obs = copy.deepcopy(cut_data_labels)
        if 'observations' in data_labels_sans_obs:
            data_labels_sans_obs.remove('observations')

        # standard deviation - absolute calculations of observations and models
        stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, 'StdDev', 
                                         data_labels, [])
        stats_dict['StdDev'] = stats_calc   

        # correlation - between observations and model
        corr_stat = plot_characteristics['taylor']['corr_stat']
        stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, corr_stat, 
                                         ['observations']*len(data_labels_sans_obs), data_labels_sans_obs)
        stats_calc = np.insert(stats_calc, obs_index, np.NaN)
        stats_dict['StdDev'] = corr_stat 

        # get maximum stddev in dataframe (if not defined)
        if not stddev_max:
            stddev_max = np.nanmax(stats_dict["StdDev"])

        # create stats dataframe
        stats_df = pd.DataFrame(data=stats_dict, 
                                index=data_labels,
                                dtype=np.float)

        # get labels 
        rlabel = stats_df.columns[0]
        xylabel = stats_df.columns[1]

        # check if we need to extend the axis to negative correlations
        extend = False
        if np.nanmin(stats_df[rlabel]) < 0:
            extend = True

        # update axis extremes and labels
        tmin, tmax, smin, smax, gl1, tf1 = self.get_taylor_diagram_ghelper_info(stddev_max, 
                                                                                plot_characteristics,
                                                                                extend)
        relevant_axis.get_grid_helper().update_grid_finder(
            extreme_finder=fa.ExtremeFinderFixed((tmin, tmax, smin, smax)),
            grid_locator1=gl1, tick_formatter1=tf1)

        # update axis position and size in dashboard
        if not self.read_instance.offline:

            # find Taylor plot position in layout
            for plot_position in range(2, 6):
                plot_type = getattr(self.read_instance, 'position_{}'.format(plot_position))
                if plot_type == 'taylor':
                    break
            
            # changing the extend reduces the size of the plot and changes its start position
            if extend:
                old_position = relevant_axis.get_position().bounds
                if plot_position == 2:
                    new_position = (0.60, 0.42, 0.288, 0.594)
                elif plot_position == 3:
                    new_position = (0.03, 0, 0.256, 0.56)
                elif plot_position == 4:
                    new_position = (0.37, 0, 0.256, 0.56)
                elif plot_position == 5:
                    new_position = (0.70, 0, 0.256, 0.56)
            else:
                if plot_position == 2:
                    new_position = (0.64, 0.57, 0.16, 0.33)
                elif plot_position == 3:
                    new_position = (0.08, 0.08, 0.16, 0.35)
                elif plot_position == 4:
                    new_position = (0.41, 0.08, 0.16, 0.35)
                elif plot_position == 5:
                    new_position = (0.69, 0.08, 0.16, 0.35)
            relevant_axis.set_position(new_position)

        # clear axis, add grid and adjust limits 
        # as suggested by the Matpotlib devs in https://github.com/matplotlib/matplotlib/issues/25426
        relevant_axis.clear()
        relevant_axis.grid(**plot_characteristics['grid'])
        relevant_axis.adjust_axes_lim()

        # create and hide annotation at (0, 0)
        if not self.read_instance.offline:
            self.canvas_instance.create_taylor_annotation()

        # adjust top axis (curve)
        relevant_axis.axis['top'].set_axis_direction('bottom')
        relevant_axis.axis['top'].toggle(ticklabels=True, label=True)
        relevant_axis.axis['top'].major_ticklabels.set_axis_direction('top')
        relevant_axis.axis['top'].major_ticklabels.set(**plot_characteristics['rtick_params'])
        relevant_axis.axis['top'].label.set_text(rlabel)
        relevant_axis.axis['top'].label.set_fontsize(plot_characteristics['rlabel']['fontsize'])
        relevant_axis.axis['top'].label.set_axis_direction('top')

        # adjust right axis (y axis)
        relevant_axis.axis['right'].set_axis_direction('top')
        relevant_axis.axis['right'].toggle(ticklabels=True)
        relevant_axis.axis['right'].major_ticklabels.set_axis_direction('bottom' if extend else 'left')
        relevant_axis.axis['right'].major_ticklabels.set(**plot_characteristics['xytick_params'])

        # adjust left axis (x axis)
        relevant_axis.axis['left'].set_axis_direction('bottom')
        relevant_axis.axis['left'].major_ticklabels.set(**plot_characteristics['xytick_params'])

        # hide bottom axis ticks and tick labels
        relevant_axis.axis['bottom'].major_ticklabels.set_visible(False)
        relevant_axis.axis['bottom'].major_ticks.set_visible(False)

        # show label
        if extend:
            relevant_axis.axis['bottom'].label.set_text(xylabel)
            relevant_axis.axis['bottom'].label.set_pad(10)
            relevant_axis.axis['bottom'].label.set_fontsize(plot_characteristics['xylabel']['fontsize'])
        else:
            relevant_axis.axis['left'].label.set_text(xylabel)
            relevant_axis.axis['left'].label.set_fontsize(plot_characteristics['xylabel']['fontsize'])  

        # add contours around observations standard deviation
        reference_stddev = stats_df[xylabel][obs_index]
        num_levels = plot_characteristics['contours']['levels']['number']
        rs, ts = np.meshgrid(np.linspace(smin, smax), np.linspace(0, tmax))
        rms = np.sqrt(reference_stddev**2 + rs**2 - 2*reference_stddev*rs*np.cos(ts))
        contours = self.taylor_polar_relevant_axis.contour(ts, rs, rms, num_levels,
            **plot_characteristics['contours']['style']['general'])

        # add contour labels
        self.taylor_polar_relevant_axis.clabel(contours, contours.levels, inline=True, fmt = '%r', fontsize=6)

        # add reference contour of observational standard deviation
        ref_x = np.linspace(0, tmax)
        ref_y = np.zeros_like(ref_x) + reference_stddev
        self.taylor_polar_relevant_axis.plot(ref_x, ref_y, **plot_characteristics['contours']['style']['obs'])

        # add models
        for data_label, stddev, corr_stat in zip(stats_df.index, 
                                                 stats_df[xylabel], 
                                                 stats_df[rlabel]):
            if data_label == 'observations':
                continue
            self.taylor_plot = self.taylor_polar_relevant_axis.plot(np.arccos(corr_stat), stddev,
                                                                    **plot_characteristics['plot'],
                                                                    mfc=self.read_instance.plotting_params[data_label]['colour'], 
                                                                    mec=self.read_instance.plotting_params[data_label]['colour'],
                                                                    label=data_label) 

            # track plot elements if using dashboard 
            if not self.read_instance.offline:
                self.track_plot_elements(data_label, 'taylor', 'plot', self.taylor_plot, bias=False)
          
    def log_axes(self, relevant_axis, log_ax, plot_characteristics, undo=False):
        """ Log plot axes.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param log_ax: which axis to log
            :type log_ax: str
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param undo: unlog plot axes
            :type undo: boolean
        """

        if not undo:
            if log_ax == 'logx':
                relevant_axis.set_xscale('log')
                
            if log_ax == 'logy':
                relevant_axis.set_yscale('log')

        else:
            if log_ax == 'logx':
                relevant_axis.set_xscale('linear')
           
            if log_ax == 'logy':
                relevant_axis.set_yscale('linear')
            
    def linear_regression(self, relevant_axis, networkspeci, data_labels, base_plot_type, plot_characteristics, 
                          plot_options=[]):
        """ Add linear regression to plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: names of plotted data arrays  
            :type data_labels: list
            :param base_plot_type: plot type, without statistical information
            :type base_plot_type: str
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plots
            :type plot_options: list
        """
    
        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get observations data (flattened and drop NaNs)
        observations_data = drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index('observations'),0,:])

        # determine if number of points per data array exceeds max limit,
        # if so subset arrays
        subset = False
        data_array_size = observations_data.size
        if data_array_size > plot_characteristics['max_points']:
            subset = True
            inds_subset = np.random.choice(data_array_size, size=plot_characteristics['max_points'], replace=False)
            observations_data = observations_data[inds_subset]

        # iterate through experiment data, making regression line to observations
        for data_label in cut_data_labels:
            if data_label != 'observations':
                # get experiement data (flattened and drop NaNs)
                experiment_data = drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:])
                # subset data if neccessary
                if subset:
                    experiment_data = experiment_data[inds_subset]
                m, b = np.polyfit(observations_data, experiment_data, deg=1)
                regression_line = relevant_axis.plot(observations_data, m*observations_data+b, 
                                                        color=self.read_instance.plotting_params[data_label]['colour'],
                                                        zorder=self.read_instance.plotting_params[data_label]['zorder']+len(cut_data_labels),
                                                        **plot_characteristics['regression'])
                
                # track plot elements if using dashboard 
                if not self.read_instance.offline:
                    self.track_plot_elements(data_label, base_plot_type, 'regression', regression_line, bias=False)

    def smooth(self, relevant_axis, networkspeci, data_labels, base_plot_type, plot_characteristics, plot_options=[]):
        """ Add smooth line to plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: names of plotted data arrays   
            :type data_labels: list
            :param base_plot_type: plot type, without statistical information
            :type base_plot_type: str
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plots
            :type plot_options: list
        """

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # iterate through plotted data arrays making smooth line
        for data_label in cut_data_labels:

            # bias plot?
            if 'bias' in plot_options:
                # skip to next data label if making bias, and data label == 'observations'
                if data_label == 'observations':
                    continue
                ts_obs = self.canvas_instance.selected_station_data[networkspeci]['timeseries']['observations']
                ts_model = self.canvas_instance.selected_station_data[networkspeci]['timeseries'][data_label] 
                ts = ts_model - ts_obs
                bias = True
            # normal plot?
            else:
                ts = self.canvas_instance.selected_station_data[networkspeci]['timeseries'][data_label]
                bias = False

            # make smooth line
            smooth_line_data = ts.rolling(plot_characteristics['smooth']['window'], 
                                     min_periods=plot_characteristics['smooth']['min_points'], 
                                     center=True).mean()
            smooth_line = relevant_axis.plot(smooth_line_data,
                                             color=self.read_instance.plotting_params[data_label]['colour'],
                                             zorder=self.read_instance.plotting_params[data_label]['zorder']+len(cut_data_labels),
                                             **plot_characteristics['smooth']['format'])

            # track plot elements if using dashboard 
            if not self.read_instance.offline:
                self.track_plot_elements(data_label, base_plot_type, 'smooth', smooth_line, bias=bias)

    def annotation(self, relevant_axis, networkspeci, data_labels, base_plot_type, plot_characteristics,
                   plot_options=[], plotting_paradigm=None):
        """ Add statistical annotations to plot.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param data_labels: names of plotted data arrays 
            :type data_labels: list
            :param base_plot_type: plot type, without statistical information
            :type base_plot_type: str
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plots
            :type plot_options: list
            :param plotting_paradigm: plotting paradigm (summary or station in offline reports)
            :type plotting_paradigm: str
        """

        # get stats wished to be annotated
        stats = plot_characteristics['annotate_stats']
        
        # if no stats defined, then return
        if len(stats) == 0:
            msg_dashboard = 'No annotation statistics are defined for {} in plot_characteristics_dashboard.json.'.format(base_plot_type)
            msg_offline = 'No annotation statistics are defined for {} in plot_characteristics_offline.json.'.format(base_plot_type)
            show_message(self.read_instance, msg=msg_dashboard, msg_offline=msg_offline)
            return

        # initialise list of strs to annotate, and colours of annotations
        str_to_annotate = []
        colours = []

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # bias plot?  
        if 'bias' in plot_options:
            bias = True
            if 'observations' in cut_data_labels:
                cut_data_labels.remove('observations')
        else:
            bias = False

        # avoid plotting stats for observations data for scatter plots
        if base_plot_type == 'scatter':
            if 'observations' in cut_data_labels:
                cut_data_labels.remove('observations')

        # generate annotation str to plot

        # show number of stations if defined
        if plot_characteristics['annotate_text']['n_stations']:
            colours.append('black')
            if self.read_instance.offline:
                if plotting_paradigm == 'station':
                    str_to_annotate.append('Stations: 1')
                else:
                    str_to_annotate.append('Stations: ' + str(len(self.read_instance.station_inds)))
            else:
                str_to_annotate.append('Stations: ' + str(len(self.read_instance.station_inds)))

        # generate annotation line by line (one line per data label, for all stats)
        for data_label_ii, data_label in enumerate(cut_data_labels):

            # get colour for data label
            colours.append(self.read_instance.plotting_params[data_label]['colour'])

            # iterate through stats to calculate
            stats_annotate = []
            for zstat in stats:

                # get zstat information
                zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat)

                # calculate stats
                if (bias) or (z_statistic_sign == 'bias') :
                    stat_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstat, 
                                                    ['observations'], [data_label])
                else:
                    stat_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstat, 
                                                    [data_label], [])

                # format annotation line
                stats_annotate.append("{0}: {1:.{2}f}".format(zstat, stat_calc[0],
                                                              plot_characteristics['annotate_text']['round_decimal_places'])) 

            # append annotation line
            if (plot_characteristics['annotate_text']['exp_labels']):
                if data_label == 'observations':
                    if 'legend' in plot_characteristics:
                        str_to_append = plot_characteristics['legend']['handles']['obs_label'] + ' | ' + ', '.join(stats_annotate)
                    elif 'legend' in self.canvas_instance.plot_characteristics_templates.keys():
                        str_to_append = self.canvas_instance.plot_characteristics_templates['legend']['handles']['obs_label'] + ' | ' + ', '.join(stats_annotate)
                    else:
                        str_to_append = 'Observations | ' + ', '.join(stats_annotate)
                else:
                    str_to_append = self.read_instance.experiments[data_label] + ' | ' + ', '.join(stats_annotate)
            else:
                str_to_append = ', '.join(stats_annotate)
            str_to_annotate.append(str_to_append)

        if plot_characteristics['annotate_text']['color'] != "":
            colours = [plot_characteristics['annotate_text']['color']]*len(cut_data_labels)

        # add annotation to plot
        # see loc options at https://matplotlib.org/3.1.0/api/offsetbox_api.html
        lines = [TextArea(line, textprops=dict(color=colour, 
                                               size=plot_characteristics['annotate_text']['fontsize'])) 
                    for line, colour in zip(str_to_annotate, colours)]
        bbox = AnchoredOffsetbox(child=VPacker(children=lines, align='left', pad=0, sep=1),
                                 loc=plot_characteristics['annotate_text']['loc'],
                                 bbox_transform=relevant_axis.transAxes)
        bbox.zorder = plot_characteristics['annotate_bbox']['zorder']
        bbox.patch.set(**plot_characteristics['annotate_bbox'])
        relevant_axis.add_artist(bbox)

        # track plot elements if using dashboard 
        if not self.read_instance.offline:
            self.track_plot_elements('ALL', base_plot_type, 'annotate', [bbox], bias=bias)

    def get_no_margin_lim(self, ax, lim):
        """ Get true limits of plot area. """

        # xlim
        if lim == 'xlim':
            xlim = ax.get_xlim()
            xwidth = xlim[1] - xlim[0]
            lower_lim = xlim[0] + (0.5 * ax.margins()[0]) / (0.5 + ax.margins()[0]) * xwidth
            upper_lim = xlim[1] - (0.5 * ax.margins()[0]) / (0.5 + ax.margins()[0]) * xwidth

        # ylim
        if lim == 'ylim':
            ylim = ax.get_ylim()
            ywidth = ylim[1] - ylim[0]
            lower_lim = ylim[0] + (0.5 * ax.margins()[1]) / (0.5 + ax.margins()[1]) * ywidth
            upper_lim = ylim[1] - (0.5 * ax.margins()[1]) / (0.5 + ax.margins()[1]) * ywidth

        return lower_lim, upper_lim

    def log_validity(self, relevant_axis, log_ax):
        """ Determine if log operation for a given axes is valid (no values <= 0).
        
            :param relevant_axis: relevant axes
            :type relevant_axis: list
            :param log_ax: which axis to log
            :type log_ax: str
            :return: validity to log axis
            :rtype: boolean
        """

        if log_ax == 'logx':
            lower_lim, _ = self.get_no_margin_lim(relevant_axis, 'xlim')
            if round(lower_lim, 2) >= 0:
                validity = True
            else:
                validity = False
        
        if log_ax == 'logy':
            lower_lim, _ = self.get_no_margin_lim(relevant_axis, 'ylim')
            if round(lower_lim, 2) >= 0:
                validity = True
            else:
                validity = False

        return validity

    def do_formatting(self, relevant_axs, relevant_data_labels, networkspeci,
                      base_plot_type, plot_type, plot_options, plotting_paradigm):
        """ Function that handles formatting of a plot axis,
            based on given plot options.

            :param relevant_axs: relevant axes
            :type relevant_axs: list
            :param relevant_data_labels: names of plotted data arrays 
            :type relevant_data_labels: list
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param base_plot_type: plot type, without statistical information
            :type base_plot_type: str
            :param plot_type: plot type
            :type plot_type: str
            :param plot_options: list of options to configure plots
            :type plot_options: list
            :param plotting_paradigm: plotting paradigm (summary or station in offline reports)
            :type plotting_paradigm: str
        """

        for relevant_ax_ii, relevant_ax in enumerate(relevant_axs):

            # log axes?
            if 'logx' in plot_options:            
                log_validity = self.log_validity(relevant_ax, 'logx')
                if log_validity:
                    self.log_axes(relevant_ax, 'logx', self.canvas_instance.plot_characteristics[plot_type])
                else:
                    msg = "Warning: It is not possible to log the x-axis "
                    msg += "in {0} with negative values.".format(plot_type)
                    print(msg)

            if 'logy' in plot_options:
                log_validity = self.log_validity(relevant_ax, 'logy')
                if log_validity:
                    self.log_axes(relevant_ax, 'logy', self.canvas_instance.plot_characteristics[plot_type])
                else:
                    msg = "Warning: It is not possible to log the y-axis "
                    msg += "in {0} with negative values.".format(plot_type)
                    print(msg)

            # annotation
            if 'annotate' in plot_options:
                if base_plot_type not in ['heatmap']:
                    self.annotation(relevant_ax, networkspeci, 
                                    relevant_data_labels[relevant_ax_ii], base_plot_type, 
                                    self.canvas_instance.plot_characteristics[plot_type],
                                    plot_options=plot_options,
                                    plotting_paradigm=plotting_paradigm)
                    # annotate on first axis
                    if base_plot_type in ['periodic', 'periodic-violin']:
                        break
            
            # regression line
            if 'regression' in plot_options:
                self.linear_regression(relevant_ax, networkspeci, 
                                       relevant_data_labels[relevant_ax_ii], base_plot_type, 
                                       self.canvas_instance.plot_characteristics[plot_type], 
                                       plot_options=plot_options)

            # smooth line
            if 'smooth' in plot_options:
                self.smooth(relevant_ax, networkspeci,
                            relevant_data_labels[relevant_ax_ii], base_plot_type, 
                            self.canvas_instance.plot_characteristics[plot_type], 
                            plot_options=plot_options)

    def track_plot_elements(self, data_label, base_plot_type, element_type, plot_object, bias=False):
        """ Function that tracks plotted lines and collections
            that will be removed/added when picking up legend elements on dashboard.

            :param data_label: name of data array to plot
            :type data_label: str
            :param base_plot_type: plot type, without statistical information
            :type base_plot_type: str
            :param element_type: type of element
            :type element_type: str
            :param plot_object: plotted element object
            :type plot_object: object
            :param bias: boolean stating if plot is a bias plot
            :type bias: boolean
        """

        # set variable name to access plot elements (absolute or bias versions)
        if not bias:
            plot_element_varname = 'absolute'
        else:
            plot_element_varname = 'bias'

        # add dictionary for plot_type elements if does not yet exist
        if base_plot_type not in self.canvas_instance.plot_elements:
            self.canvas_instance.plot_elements[base_plot_type] =  {'active': 'absolute', 'absolute': {}} 

        # add plot_element_varname if does not yet exist
        if plot_element_varname not in self.canvas_instance.plot_elements[base_plot_type]:
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname] = {}

        # add dictionary for data label if does not yet exist
        if data_label not in self.canvas_instance.plot_elements[base_plot_type][plot_element_varname]:
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label] = {}

        # add list for element type if does not yet exist
        if element_type not in self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label]:
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] = []
        # if does exist already then remove plot element return, as element has already been plotted
        else:
            return 

        # track plot elements
        # periodic plot specific elements 
        if (base_plot_type in ['periodic', 'periodic-violin']) & (data_label != 'ALL'):
            # add list of collections
            if isinstance(plot_object, dict):
                self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['bodies']
            # add list of lines
            elif isinstance(plot_object, list):
                self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object
        # boxplot plot specific elements
        elif (base_plot_type == 'boxplot') & (data_label != 'ALL'):
            # add lines for all boxplot elements
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['boxes']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['medians']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['whiskers']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['caps']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['fliers']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['means']
        # do not save elements for plot objects that can not be made invisisble
        elif (base_plot_type in ['metadata', 'map']) & (data_label != 'ALL'):
            pass
        # all other plot elements
        else:
            # add list of lines
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object
            
        # set element visibility
        if (data_label not in self.canvas_instance.plot_elements['data_labels_active']) & (data_label != 'ALL'):
            for element in self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type]:
                element.set_visible(False)

    def harmonise_xy_lims_paradigm(self, relevant_axs, base_plot_type, plot_characteristics, plot_options, 
                                   xlim=None,  ylim=None, relim=False, autoscale=False, autoscale_x=False, 
                                   autoscale_y=False, bias_centre=False):
        """ Harmonise xy limits across paradigm of plot type, unless axis limits have been defined.
        
            :param relevant_axs: relevant axes
            :type relevant_axs: list
            :param base_plot_type: plot type, without statistical information
            :type base_plot_type: str
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: list of options to configure plots
            :type plot_options: list
            :param xlim: xlimits to set
            :type xlim: list
            :param ylim: ylimits to set
            :type ylim: list
            :param relim: turn on relimiting of axes limits (when updating plotted data on axis)
            :type relim: boolean
            :param autoscale: Autoscale the axis view to the data (both x and y axes)
            :type autoscale: boolean
            :param autoscale_x: Autoscale the x axis view to the data
            :type autoscale_x: boolean
            :param autoscale_y: Autoscale the x axis view to the data
            :type autoscale_y: boolean
            :param bias_centre: centre bias plots at 0 on the y axis
            :type bias_centre: boolean   
        """
        
        # initialise arrays to save lower and upper limits in all axes
        all_xlim_lower = []
        all_xlim_upper = []
        all_ylim_lower = []
        all_ylim_upper = []

        # initialise variables for setting axis limits
        xlim_min = None
        xlim_max = None
        ylim_min = None
        ylim_max = None
        
        if not isinstance(relevant_axs, list):
            # if changes only apply to one axis, put it in list
            if not isinstance(relevant_axs, dict):
                relevant_axs = [relevant_axs]
            # transform dictionaries into lists
            else:
                relevant_axs = [relevant_axs[relevant_temporal_resolution] for 
                                relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions]

        # get mapped resolution per axis for periodic plots
        if base_plot_type in ['periodic', 'periodic-violin']:
            mapped_resolutions = self.read_instance.relevant_temporal_resolutions*(int(len(relevant_axs)/len(self.read_instance.relevant_temporal_resolutions)))

        # remove any axes from relevant_axs which are not active (only for offline)
        if self.read_instance.offline:
            relevant_axs_active = []
            mapped_resolutions_active = []
            for ax_ii, ax in enumerate(relevant_axs):
                if ax.axison:
                    relevant_axs_active.append(ax)
                    if base_plot_type in ['periodic', 'periodic-violin']:
                        mapped_resolutions_active.append(mapped_resolutions[ax_ii])
        else:
            relevant_axs_active = relevant_axs
            if base_plot_type in ['periodic', 'periodic-violin']:
                mapped_resolutions_active = mapped_resolutions

        # get lower and upper limits across all relevant axes
        for ax in relevant_axs_active:
            if 'equal_aspect' in plot_characteristics:
                self.set_equal_axes(ax, plot_options)
            else:
                ax.set_aspect('auto')
            if relim:
                ax.relim(visible_only=True)
            if autoscale:
                ax.autoscale(tight=False)
            if autoscale_x:
                ax.autoscale(axis='x', tight=False)
            if autoscale_y:
                ax.autoscale(axis='y', tight=False)
            if xlim is None and ('xlim' not in plot_characteristics):
                if base_plot_type not in ['periodic', 'periodic-violin', 'timeseries']:
                    xlim_lower, xlim_upper = ax.get_xlim()
                elif base_plot_type == 'timeseries':
                    xlim_lower, xlim_upper = self.get_no_margin_lim(ax, 'xlim')
                    try:
                        xlim_lower = num2date(xlim_lower)
                        xlim_upper = num2date(xlim_upper)
                    except ValueError:
                        continue
                if base_plot_type not in ['periodic', 'periodic-violin']:
                    all_xlim_lower.append(xlim_lower)
                    all_xlim_upper.append(xlim_upper)
            if ylim is None and ('ylim' not in plot_characteristics):
                ylim_lower, ylim_upper = ax.get_ylim()
                all_ylim_lower.append(ylim_lower)
                all_ylim_upper.append(ylim_upper)

        # get minimum and maximum from all axes and set limits
        for ax in relevant_axs_active:
            # get xlim
            if xlim is None and ('xlim' not in plot_characteristics):
                if base_plot_type not in ['periodic', 'periodic-violin']:
                    xlim_min = np.min(all_xlim_lower)
                    xlim_max = np.max(all_xlim_upper)
                    xlim = xlim_min, xlim_max
            elif 'xlim' in plot_characteristics:
                xlim = plot_characteristics['xlim']

            # set xlim
            if xlim is not None:
                if (base_plot_type not in ['timeseries', 'boxplot']) and ('equal_aspect' not in plot_characteristics):
                    ax.set_xlim(xlim)

            # get ylim
            if ylim is None and ('ylim' not in plot_characteristics):
                ylim_min = np.min(all_ylim_lower) 
                ylim_max = np.max(all_ylim_upper)
                # if have bias_centre option, centre around zero
                if ('bias' in plot_options) & (bias_centre):                    
                    if np.abs(np.max(all_ylim_upper)) >= np.abs(np.min(all_ylim_lower)):
                        ylim_min = -np.abs(np.max(all_ylim_upper))
                        ylim_max = np.abs(np.max(all_ylim_upper))
                    elif np.abs(np.max(all_ylim_upper)) < np.abs(np.min(all_ylim_lower)):
                        ylim_min = -np.abs(np.min(all_ylim_lower))
                        ylim_max = np.abs(np.min(all_ylim_lower))
                ylim = ylim_min, ylim_max
            elif 'ylim' in plot_characteristics:
                ylim = plot_characteristics['ylim']

            # set ylim
            if (ylim is not None) and ('equal_aspect' not in plot_characteristics):
                ax.set_ylim(ylim)
            
        # get minimum and maximum from all axes and set limits for periodic plots
        if base_plot_type in ['periodic', 'periodic-violin']:
            if xlim is None and ('xlim' not in plot_characteristics):
                for temporal_resolution, sub_ax in zip(mapped_resolutions_active, relevant_axs_active):
                    # adjust plot x axis to have correct margin on edges
                    xlim_lower, xlim_upper = sub_ax.get_xlim()
                    first_valid_x = self.canvas_instance.periodic_xticks[temporal_resolution][(np.abs(self.canvas_instance.periodic_xticks[temporal_resolution] - xlim_lower)).argmin()]
                    last_valid_x = self.canvas_instance.periodic_xticks[temporal_resolution][(np.abs(self.canvas_instance.periodic_xticks[temporal_resolution] - xlim_upper)).argmin()]
                    if temporal_resolution == 'hour':
                        xlim_lower = first_valid_x - 0.65
                        xlim_upper = last_valid_x + 0.65
                    elif temporal_resolution == 'dayofweek':
                        xlim_lower = first_valid_x - 0.55
                        xlim_upper = last_valid_x + 0.55
                    elif temporal_resolution == 'month':
                        xlim_lower = first_valid_x - 0.55
                        xlim_upper = last_valid_x + 0.55
                    xlim = xlim_lower, xlim_upper
                    sub_ax.set_xlim(xlim)
            elif 'xlim' in plot_characteristics:
                xlim = plot_characteristics['xlim']
                for temporal_resolution, sub_ax in zip(mapped_resolutions_active, relevant_axs_active):
                    sub_ax.set_xlim(xlim)

        # get minimum and maximum from all axes and set limits for timeseries
        elif base_plot_type == 'timeseries':
            if plot_characteristics['xtick_alteration']['define']:
            
                # get steps for all data labels
                steps = pd.date_range(xlim[0], xlim[1], freq=self.read_instance.active_frequency_code)

                # get number of months and days
                n_months = (12*(xlim[1].year - xlim[0].year) + (xlim[1].month - xlim[0].month))
                n_days = (xlim[1] - xlim[0]).days

                # get months that are complete
                months_start = pd.date_range(xlim[0], xlim[1], freq='MS')
                months_end = pd.date_range(xlim[0], xlim[1], freq='M')
                if months_start.size > 1:
                    if (xlim[1] - months_end[-1]).days >= 1:
                        months = months_start[:-1]
                    else:
                        months = months_start
                else:
                    months = months_start

                # show hours if number of days is less than 7
                if n_days < 7:
                    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d %H:%M'))
                else:
                    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d'))

                # define time slices
                if n_months >= 3:
                    steps = months
                slices = int(np.ceil(len(steps) / int(plot_characteristics['xtick_alteration']['n_slices'])))

                # use default axes if the number of timesteps is lower than the number of slices
                if slices >= 1:
                    xticks = steps[0::slices]
                else:
                    xticks = ax.xaxis.get_ticks()

                # transform to numpy.datetime64
                if not isinstance(xticks[0], np.datetime64):
                    xticks = [x.to_datetime64() for x in xticks]
                if not isinstance(xlim[1], np.datetime64):
                    xlim = xlim[0], np.datetime64(xlim[1])

                # add last step to xticks
                if plot_characteristics['xtick_alteration']['last_step'] and (xticks[-1] != xlim[1]):
                    xticks = np.append(xticks, xlim[1])

                # set modified xticks
                for ax in relevant_axs_active:
                    ax.xaxis.set_ticks(xticks)

    def set_axis_title(self, relevant_axis, title, plot_characteristics, 
                       relevant_temporal_resolutions=['hour']):
        """ Set title of plot axis.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param title: axis title
            :type title: str
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param relevant_temporal_resolutions: relevant temporal resolutions to plot title on   
            :type relevant_temporal_resolutions: list
        """    

        # return if title is empty str
        if title == '':
            return

        # get appropriate axis for plotting label for plots with multiple sub-axes
        axs_to_set_title = []
        if isinstance(relevant_axis, dict):
            for relevant_temporal_resolution, sub_ax in relevant_axis.items():
                if relevant_temporal_resolution in relevant_temporal_resolutions:
                    axs_to_set_title.append(sub_ax)
                    break
        else:
            axs_to_set_title.append(relevant_axis)

        # set title for appropriate axes
        axis_title_characteristics = copy.deepcopy(plot_characteristics['axis_title'])
        axis_title_characteristics['label'] = title
        for relevant_axis in axs_to_set_title:
            relevant_axis.set_title(**axis_title_characteristics)

    def set_axis_label(self, relevant_axis, label_ax, label, plot_characteristics, 
                       relevant_temporal_resolutions=['hour', 'month']):
        """ Set label of plot axis.

            :param relevant_axis: axis to plot on 
            :type relevant_axis: object
            :param label_ax: which axis to set label of
            :type label_ax: str
            :param label: axis label
            :type label: str
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
        """

        # return if label is empty str
        if label == '':
            return

        # get appropriate axis for plotting label for plots with multiple sub-axes (hour and month axes)
        axs_to_set_label = []
        if isinstance(relevant_axis, dict):
            for relevant_temporal_resolution, sub_ax in relevant_axis.items():
                if relevant_temporal_resolution in relevant_temporal_resolutions:
                    axs_to_set_label.append(sub_ax)
                # remove day of week axis label if setting ylabel
                if (relevant_temporal_resolution == 'dayofweek') & (label_ax == 'y'):                           
                    sub_ax.yaxis.set_tick_params(which='both', labelleft=False)
                    sub_ax.set_ylabel('')
        else:
            axs_to_set_label.append(relevant_axis)

        # set label for appropriate axes
        for relevant_axis in axs_to_set_label:
            if label_ax == 'x':
                axis_label_characteristics = copy.deepcopy(plot_characteristics['xlabel'])
                axis_label_characteristics['xlabel'] = label
                relevant_axis.set_xlabel(**axis_label_characteristics)
            elif label_ax == 'y':
                axis_label_characteristics = copy.deepcopy(plot_characteristics['ylabel'])
                axis_label_characteristics['ylabel'] = label
                relevant_axis.set_ylabel(**axis_label_characteristics)

    def get_markersize(self, relevant_axis, base_plot_type, networkspeci, plot_characteristics, data=None):
        """ Set markersize for plot.
        
            :param base_plot_type: plot type, without statistical information
            :type base_plot_type: str
            :param networkspeci: str of currently active network and species 
            :type networkspeci: str
            :param plot_characteristics: plot characteristics  
            :type plot_characteristics: dict
            :param data: data array to be plotted
            :type data: numpy array
        """

        if base_plot_type in ['timeseries', 'scatter']:
            
            if plot_characteristics['plot']['markersize'] == '':

                # get minimum number of non-NaN data points for plot type across data labels
                min_points = np.count_nonzero(~np.isnan(data))

                # configure size of plots if have very few points
                if min_points < plot_characteristics['markersize_npoints_threshold']:
                    markersize = plot_characteristics['markersize']['few_points'] 
                else:
                    markersize = plot_characteristics['markersize']['standard'] 

                # add to plot_characteristics json
                plot_characteristics['plot']['markersize'] = markersize
        
        elif base_plot_type == 'map':

            if (plot_characteristics['plot']['s'] == '') or (self.map_markersize_from_density):

                # calculate marker size considering points density
                n_points = len(self.read_instance.station_longitudes[networkspeci][self.canvas_instance.active_map_valid_station_inds])
                
                # calculate figure area and density
                # divide area by 1000 so the function below makes sense
                area = (relevant_axis.bbox.width * relevant_axis.bbox.height) / 1000
                density = n_points / area

                # marker size is calculated using an exponential equation
                # the maximum size is 40 (very low densities)
                # see https://earth.bsc.es/gitlab/ac/Providentia/-/issues/210
                plot_characteristics['plot']['s'] = 1.2**(-density)*40
                self.map_markersize_from_density = True
                
            else:

                self.map_markersize_from_density = False
