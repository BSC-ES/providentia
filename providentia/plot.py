""" Functions to generate plots """

import copy
import json
import os
import sys
import yaml

import cartopy
import cartopy.crs as ccrs
from KDEpy import FFTKDE
from itertools import groupby
import matplotlib
import matplotlib.image as mpimg
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.projections import PolarAxes
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as fa
import numpy as np
from packaging.version import Version
import pandas as pd
import pyproj
import seaborn as sns

from .calculate import ExpBias, Stats
from .dashboard_interactivity import HoverAnnotation
from .statistics import boxplot_inner_fences, calculate_statistic, get_z_statistic_info, get_z_statistic_type
from .read_aux import drop_nans
from .plot_aux import (create_chunked_timeseries, get_multispecies_aliases, 
                       get_taylor_diagram_ghelper_info, kde_fft, merge_cells, periodic_labels, 
                       periodic_xticks, round_decimal_places, temp_axis_dict)


# speed up transformations in cartopy
pyproj.set_use_global_context()

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])

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

    def set_plot_characteristics(self, plot_types, zstat=False, data_labels=None, format={}):
        """ Iterate through all plots to make, and determine if they can and cannot be made.
            Update plot characteristics associated with specific plot types due to plot options. 

            :param plot_types: Plot types to create 
            :type plot_types: list  
            :param zstat: Statistic 
            :type zstat: str
            :param data_labels: Data arrays to plot
            :type data_labels: list 
            :param format: Dictionary to overwrite default formatting 
            :type format: dict
        """

        # if data_labels are not defined, take all in memory
        if data_labels is None:
            data_labels = copy.deepcopy(self.read_instance.data_labels)

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
            
            # check if plot type is correct for offline and interactive modes
            if (self.read_instance.offline) or (self.read_instance.interactive):

                # remove plots where setting 'obs' and 'bias' options together
                if ('obs' in plot_options) & ('bias' in plot_options): 
                    print(f"Warning: {plot_type} cannot not be created as 'obs' and 'bias' options set together.")
                    valid_plot_type = False

                # if no experiments are defined, remove all bias plots 
                elif ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                    if len(data_labels) == 1:
                        print(f'Warning: No experiments defined, so {plot_type} bias plot cannot be created.')
                        valid_plot_type = False

                # if are making an experiment bias plot, and temporal_colocation is off, then remove plot
                elif (z_statistic_type == 'expbias') & (not self.read_instance.temporal_colocation):
                    print(f'Warning: To calculate the experiment bias stat {zstat}, temporal_colocation must be set to True, so {plot_type} plot cannot be created.')
                    valid_plot_type = False

            # add new keys to make plots with stats (map, periodic, heatmap, table)
            if zstat:

                # get base plot type (without stat and options)
                base_plot_type = plot_type.split('-')[0] 
                # combine basic and expbias stats dicts together
                stats_dict = {**self.read_instance.basic_stats, **self.read_instance.expbias_stats}
                
                # check if plot type is correct for offline and interactive modes
                if (self.read_instance.offline) or (self.read_instance.interactive):

                    # check all defined plot options are allowed for current plot type
                    invalid_plot_options = [plot_option for plot_option in plot_options if plot_option not in self.canvas_instance.plot_characteristics_templates[base_plot_type]['plot_options']]
                    if len(invalid_plot_options) > 0:
                        print(f'Warning: {plot_type} cannot be created as {invalid_plot_options} plot options are not valid.')
                        valid_plot_type = False
                    
                    # check desired statistic is defined in stats dict
                    elif base_zstat not in stats_dict:
                        print(f"Warning: {plot_type} cannot be created as {base_zstat} not defined in Providentia's statistical library.")
                        valid_plot_type = False
                    
                    # remove plots where setting 'obs', but z_statistic_sign is 'bias'
                    elif ('obs' in plot_options) & (z_statistic_sign == 'bias'):
                        print(f"Warning: {plot_type} cannot be created as are plotting a bias statistic but 'obs' option is set.")
                        valid_plot_type = False

                if not valid_plot_type:
                    if self.read_instance.offline:
                        if plot_type in self.read_instance.summary_plots_to_make:
                            self.read_instance.summary_plots_to_make.remove(plot_type)
                        if plot_type in self.read_instance.station_plots_to_make:
                            self.read_instance.station_plots_to_make.remove(plot_type)
                    elif self.read_instance.interactive:
                        return valid_plot_type
                    continue

                # add information for plot type 
                # first try get it from custom plot charcteristics, and then from base plot type template 
                if plot_type in self.canvas_instance.plot_characteristics_templates:
                    self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[plot_type])
                else:
                    try:
                        self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[base_plot_type])
                    except KeyError:
                        error = f'Error: Plot type {plot_type} is not available. Remove from settings/report_plots.yaml'
                        sys.exit(error)

                # overwrite default plot characteristics with custom formatting
                for format_var in format:
                    self.canvas_instance.plot_characteristics[plot_type][format_var] = format[format_var]

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

                # check if plot type is correct for offline and interactive modes
                if (self.read_instance.offline) or (self.read_instance.interactive):

                    # check all defined plot options are allowed for current plot type
                    invalid_plot_options = [plot_option for plot_option in plot_options if plot_option not in self.canvas_instance.plot_characteristics_templates[base_plot_type]['plot_options']]
                    if len(invalid_plot_options) > 0:
                        print(f'Warning: {plot_type} cannot be created as {invalid_plot_options} plot options are not valid.')
                        valid_plot_type = False

                    # warning for scatter, taylor and fairmode plots if the temporal colocation is not active 
                    elif (base_plot_type in ['scatter', 'taylor', 'fairmode-target']) & (not self.read_instance.temporal_colocation):
                        print(f'Warning: {plot_type} cannot be created as temporal colocation is not active.')
                        valid_plot_type = False

                    # warning for scatter, taylor and fairmode plots if have no experiments
                    elif (base_plot_type in ['scatter', 'taylor', 'fairmode-target']) & (len(data_labels) == 1):
                        print(f'Warning: No experiments defined, so {plot_type} cannot be created.')
                        valid_plot_type = False

                    # warning for timeseries bias plot if the temporal colocation is not active
                    elif ('timeseries' == base_plot_type) & ('bias' in plot_options) & (not self.read_instance.temporal_colocation):
                        print(f'Warning: {plot_type} cannot be created as temporal colocation is not active.')
                        valid_plot_type = False

                    # warning for timeseries bias plot if have no experiments
                    elif ('timeseries' == base_plot_type) & ('bias' in plot_options) & (len(data_labels) == 1):
                        print(f'Warning: No experiments defined, so {plot_type} cannot be created.')
                        valid_plot_type = False

                # break loop if the plot type is not valid and remove plot type from lists
                if not valid_plot_type:
                    if self.read_instance.offline:
                        if plot_type in self.read_instance.summary_plots_to_make:
                            self.read_instance.summary_plots_to_make.remove(plot_type)
                        if plot_type in self.read_instance.station_plots_to_make:
                            self.read_instance.station_plots_to_make.remove(plot_type)
                    elif self.read_instance.interactive:
                        return valid_plot_type
                    continue

                # add information for plot type for base plot type 
                if plot_type in self.canvas_instance.plot_characteristics_templates:
                    self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[plot_type])
                else:
                    self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[base_plot_type])

                # overwrite default plot characteristics with custom formatting
                for format_var in format:
                    self.canvas_instance.plot_characteristics[plot_type][format_var] = format[format_var]

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

        # return valid plot type if interactive mode
        if self.read_instance.interactive:
            return valid_plot_type

    def make_legend_handles(self, plot_characteristics_legend, data_labels=None, set_obs=True):
        """ Make legend element handles.
        
            :param plot_characteristics_legend: Plot characteristics for relevant legend
            :type plot_characteristics_legend: dict
            :param data_labels: Data arrays to plot
            :type data_labels: list 
            :param set_obs: Indicates if to set observations in legend or not  
            :type set_obs: boolean 
            :return: plot_characteristics_legend with handles updated
            :rtype: dict
        """

        # if data_labels are not defined, take all in memory
        if data_labels is None:
            data_labels = copy.deepcopy(self.read_instance.data_labels)

        # create legend elements
        legend_elements = []

        # add observations element, if available, and set_obs == True
        if (self.read_instance.observations_data_label in data_labels) and (set_obs):
            legend_elements.append(Line2D([0], [0], 
                                marker=plot_characteristics_legend['handles']['marker'], 
                                color=plot_characteristics_legend['handles']['color'],
                                markerfacecolor=self.read_instance.plotting_params[self.read_instance.observations_data_label]['colour'],
                                markersize=plot_characteristics_legend['handles']['markersize'], 
                                label=self.read_instance.observations_data_label))
                                  
        # add element for each experiment
        for experiment in data_labels:
            if experiment != self.read_instance.observations_data_label:
                # add experiment element
                legend_elements.append(Line2D([0], [0], 
                                              marker=plot_characteristics_legend['handles']['marker'],  
                                              color=plot_characteristics_legend['handles']['color'],
                                              markerfacecolor=self.read_instance.plotting_params[experiment]['colour'],
                                              markersize=plot_characteristics_legend['handles']['markersize'],
                                              label=experiment))
        
        plot_characteristics_legend['plot']['handles'] = legend_elements
        
        return plot_characteristics_legend

    def make_experiment_domain_polygons(self, data_labels=None):
        """ Make experiment domain polygons.
            
            :param data_labels: Data arrays to plot
            :type data_labels: list 
            :return: grid_edge_polygons
            :rtype: list
        """

        # if data_labels are not defined, or are just observations, then take all labels in memory
        if data_labels is None:
            data_labels = copy.deepcopy(self.read_instance.data_labels)
        elif len(data_labels) == 1:
            if data_labels[0] == self.read_instance.observations_data_label:
                data_labels = copy.deepcopy(self.read_instance.data_labels)

        grid_edge_polygons = []
        
        # iterate through read experiments and plot grid domain edges on map
        for experiment in data_labels:
            if experiment != self.read_instance.observations_data_label:
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

    def make_header(self, pdf, plot_characteristics):
        """ Make header.
        
            :param plot_characteristics: Plot characteristics 
            :type plot_characteristics: dict
        """

        # create header page
        page = plt.figure(**plot_characteristics['figure'])

        # get logo
        logo_path = os.path.join(PROVIDENTIA_ROOT, 'assets/logoline.png')
        logo = mpimg.imread(logo_path)

        # place logo on top center
        page_width, page_height = page.get_size_inches()
        xo = (page_width * page.dpi) / 2 + 500
        yo = page_height * page.dpi * 2 + 700
        page.figimage(logo, xo=int(xo), yo=int(yo))

        # set header title
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

    def make_metadata(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options):
        """ Make metadata summary plot.

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param data_labels: Data arrays to plot
            :type data_labels: str
            :param plot_characteristics: Plot characteristics 
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
        """

        # initialise string to plot on axis
        str_to_plot = ''
        
        # get number of selected stations
        n_stations = self.canvas_instance.selected_station_metadata[networkspeci].shape[0]

        # set first line of metadata print to be either selected stations or number of stations 
        if n_stations == 1:
            var_str_name = 'name_one' 
            station_references = self.canvas_instance.selected_station_metadata[networkspeci][
                                 'station_reference'].flatten()
            station_reference = station_references[~pd.isnull(station_references)].astype(str)[0]
            str_to_plot += 'Station: {}'.format(station_reference)
        else:
            var_str_name = 'name_multiple' 
            str_to_plot += '{} Stations'.format(n_stations)
        
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
            
            # check if var is in metadata
            if ghost_var not in self.read_instance.metadata_vars_to_read:
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
                                                    np.nanmedian(self.canvas_instance.selected_station_metadata[
                                                    networkspeci][ghost_var].astype(np.float32)), 
                                                    ghost_var_dict['dp'])

            # if str then get unique elements or percentage dependent on n uniques
            else:
                # gather all selected station metadata for current meta variable
                all_current_meta = self.canvas_instance.selected_station_metadata[networkspeci][
                                   ghost_var].flatten()
                # remove NaNs
                all_current_meta = all_current_meta[~pd.isnull(all_current_meta)].astype(str)

                # get counts of all unique metadata elements across selected stations
                unique_meta, meta_counts = np.unique(all_current_meta, return_counts=True)
                # get number of unique metadata elements across selected stations
                n_unique_meta = len(unique_meta)

                # 0 unique metadata elements, then must be all NaN, so set as NaN
                if n_unique_meta == 0:
                    str_to_plot += '{}: {}'.format(ghost_var_dict[var_str_name], 'nan')
                # 1 unique metadata element? then return it
                elif n_unique_meta == 1:
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
        if self.read_instance.offline: 
            
            # get axis dimensions in pixels
            ax_width_px = relevant_axis.bbox.width * plot_characteristics['figure']['nrows']
            
        elif self.read_instance.interactive:

            # get axis dimensions in pixels
            ax_width_px = relevant_axis.bbox.width

        else:
            # get axis bounding box
            ax_bbox = relevant_axis.get_window_extent().transformed(self.canvas_instance.figure.dpi_scale_trans.inverted())
            
            # get axis dimensions in inches
            ax_width_inches = ax_bbox.width

            # get axis dimensions in pixels
            ax_width_px = ax_width_inches * self.canvas_instance.figure.dpi

        # automatically sets limit as figure width
        plot_txt._get_wrap_line_width = lambda: ax_width_px

        # track plot elements if using dashboard 
        if (not self.read_instance.offline) and (not self.read_instance.interactive):
            self.track_plot_elements(self.read_instance.observations_data_label, 'metadata', 'plot', [plot_txt], bias=False)

    def make_map(self, relevant_axis, networkspeci, plot_characteristics, plot_options, zstat=None, labela='', 
                 labelb=''):
        """ Make map plot.

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
            :param zstat: Statistic to plot
            :type zstat: str
            :param labela: Label of first dataset
            :type labela: str
            :param labelb: Label of second dataset (if defined then a bias plot is made)
            :type labelb: str
        """

        # calculate statistic
        z_statistic, active_map_valid_station_inds = calculate_statistic(self.read_instance, self.canvas_instance, 
                                                                         networkspeci, zstat, [labela], [labelb], 
                                                                         map=True)

        # get marker size (for offline and interactive)
        if (self.read_instance.offline) or (self.read_instance.interactive):
            self.get_markersize(relevant_axis, 'map', networkspeci, plot_characteristics, 
                                active_map_valid_station_inds=active_map_valid_station_inds)
        # if using dashboard make z_statistic and active_map_valid_station_inds class variables
        else:
            self.canvas_instance.z_statistic = z_statistic
            self.canvas_instance.active_map_valid_station_inds = active_map_valid_station_inds

        # plot new station points on map - coloured by currently active z statisitic
        self.stations_scatter = relevant_axis.scatter(self.read_instance.station_longitudes[networkspeci][active_map_valid_station_inds], 
                                                      self.read_instance.station_latitudes[networkspeci][active_map_valid_station_inds], 
                                                      c=z_statistic, transform=self.canvas_instance.datacrs,
                                                      **plot_characteristics['plot'])
        
        # track plot elements if using dashboard 
        if (not self.read_instance.offline) and (not self.read_instance.interactive):
            self.track_plot_elements(self.read_instance.observations_data_label, 'map', 'plot', [self.stations_scatter], bias=False)

    def make_timeseries(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options, 
                        chunk_stat=None, chunk_resolution=None):
        """ Make timeseries plot.

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param data_labels: Data arrays to plot
            :type data_labels: list
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
            :param chunk_stat: Chunk statistic
            :type chunk_stat: str
            :param chunk_resolution: Chunk resolution
            :type chunk_resolution: str
        """

        # skip making timeseries (points) for offline and interactive mode
        # we do not apply this in the dashboard to avoid being unable to see the points on certain changes
        if ((self.read_instance.offline) or (self.read_instance.interactive)) and ('hidedata' in plot_options):
            return

        # if 'obs' in plot_options, set data labels to just observations data label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

        # get bias and set if bias line will be added
        if 'bias' in plot_options:
            bias =  True
            add_bias_line = True
        else:
            bias = False
            add_bias_line = False

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get chunking stat and resolution in dashboard
        if (not self.read_instance.offline) and (not self.read_instance.interactive):
            chunk_stat = self.canvas_instance.timeseries_chunk_stat.currentText()
            chunk_resolution = self.canvas_instance.timeseries_chunk_resolution.currentText()
            chunk_stat = None if chunk_stat == 'None' else chunk_stat
            chunk_resolution = None if chunk_resolution == 'None' else chunk_resolution
        
        # chunk timeseries
        if (chunk_stat is not None) and (chunk_resolution is not None):
            timeseries_data = create_chunked_timeseries(self.read_instance, self.canvas_instance, chunk_stat, 
                                                        chunk_resolution, networkspeci, cut_data_labels, bias)

            # if it is a bias chunk statistic, add bias line
            z_statistic_type = get_z_statistic_type(chunk_stat)
            if z_statistic_type == 'expbias':
                add_bias_line = True
        # normal timeseries
        else:
            timeseries_data = self.canvas_instance.selected_station_data[networkspeci]["timeseries"]
        
        # plot horizontal line across x axis at 0 if bias plot
        if add_bias_line:
            if 'bias' not in plot_options:
               plot_characteristics['bias_line']['y'] = self.read_instance.expbias_stats[chunk_stat]['minimum_bias']
            bias_line = relevant_axis.axhline(**plot_characteristics['bias_line'])
            # track plot elements if using dashboard 
            if (not self.read_instance.offline) and (not self.read_instance.interactive):
                self.track_plot_elements('ALL', 'timeseries', 'bias_line', [bias_line], bias=bias)

        # iterate through data labels
        for data_label in cut_data_labels:

            # bias plot?
            if bias:

                # skip if data label is for observations
                if data_label == self.read_instance.observations_data_label:
                    continue
                
                # chunk bias timeseries
                if (chunk_stat is not None) and (chunk_resolution is not None):
                    ts = timeseries_data[data_label] 
                # normal bias timeseries
                else:
                    ts_obs = timeseries_data[self.read_instance.observations_data_label]
                    ts_model = timeseries_data[data_label] 
                    ts = ts_model - ts_obs

            else:
                ts = timeseries_data[data_label]
            
            # get marker size (for offline and interactive)
            if (self.read_instance.offline) or (self.read_instance.interactive):
                self.get_markersize(relevant_axis, 'timeseries', networkspeci, plot_characteristics, data=ts)

            # make timeseries plot
            self.timeseries_plot = relevant_axis.plot(ts, 
                                                      color=self.read_instance.plotting_params[data_label]['colour'], 
                                                      **plot_characteristics['plot'])

            # update maximum smooth value
            if (not self.read_instance.offline) and (not self.read_instance.interactive):
                self.canvas_instance.timeseries_smooth_window_sl.setMaximum(len(ts))
                # To get straight line
                # if self.canvas_instance.timeseries_smooth_window_sl.value() != (len(ts)*2 - 1):
                #     self.canvas_instance.timeseries_smooth_window_sl.setMaximum(int(len(ts)*2 - 1))

            # track plot elements if using dashboard 
            if (not self.read_instance.offline) and (not self.read_instance.interactive):
                self.track_plot_elements(data_label, 'timeseries', 'plot', self.timeseries_plot, bias=bias)

    def make_periodic(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options, zstat=None):
        """ Make period or period-violin plot.

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param data_labels: Data arrays to plot
            :type data_labels: list
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
            :param zstat: Statistic
            :type zstat: str

        """

        # if 'obs' in plot_options, set data labels to just observations data label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

        # determine if 'bias' in plot_options
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # get zstat information
        if zstat is not None:
            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat)

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]
        
        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get number of experiments in data  labels
        if self.read_instance.observations_data_label in cut_data_labels:
            n_exps = len(cut_data_labels) - 1
        else:
            n_exps = len(cut_data_labels) 

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
            if zstat is None:

                # calculate medians
                medians = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, 'Median', 
                                              cut_data_labels, [], period=relevant_temporal_resolution)

                # calculate PDF for data label
                period_x_grid, PDFs_sampled = self.make_distribution(relevant_axis, networkspeci, data_labels, 
                                                                     plot_characteristics, plot_options,
                                                                     violin_resolution=relevant_temporal_resolution)

                # iterate through data labels and plot violins
                for data_label_ii, data_label in enumerate(cut_data_labels): 

                    # list to save all violins per data label
                    violins = []

                    # get median zorder
                    median_zorder = self.read_instance.plotting_params[data_label]['zorder']+len(cut_data_labels)
                    
                    # get alpha and violin fill information 
                    if data_label == self.read_instance.observations_data_label:
                        alpha = plot_characteristics['violin_alphas']['alpha_obs']
                        violin_fill = plot_characteristics['violin_fill_obs']
                    else:
                        alpha = plot_characteristics['violin_alphas']['alpha_exp']
                        if (n_exps == 1) or ('individual' in plot_options):
                            violin_fill = plot_characteristics['violin_fill_1model']
                        else:
                            violin_fill = plot_characteristics['violin_fill_2+models']

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
                            
                            PDF_sampled = 0.5 * plot_characteristics['violin_widths'] * PDF_sampled / PDF_sampled.max()

                            # make violin plot (filled or unfilled)
                            if violin_fill:
                                self.violin_plot = relevant_sub_ax.fill_betweenx(x_grid, -PDF_sampled + period_pos, PDF_sampled + period_pos,
                                                                                 facecolor=self.read_instance.plotting_params[data_label]['colour'], 
                                                                                 alpha=alpha,
                                                                                 **plot_characteristics['plot']['violin'])
                            else:
                                self.violin_plot = relevant_sub_ax.fill_betweenx(x_grid, -PDF_sampled + period_pos, PDF_sampled + period_pos,
                                                                                 facecolor='None', edgecolor=self.read_instance.plotting_params[data_label]['colour'], 
                                                                                 **plot_characteristics['plot']['violin'])

                            # if have more than 1 valid data array (both obs and model), 
                            # split the violin plot across the horizontal
                            # (observations on left, experiment violin_plots on right)
                            if ((n_exps > 0) and (self.read_instance.observations_data_label in cut_data_labels) and
                               ('individual' not in plot_options)):
                                m = np.mean(self.violin_plot.get_paths()[0].vertices[:, 0])
                                # observations on left
                                if data_label == self.read_instance.observations_data_label:
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
                    if (not self.read_instance.offline) and (not self.read_instance.interactive):
                        self.track_plot_elements(data_label, 'periodic-violin', 'violin_plot_{}'.format(relevant_temporal_resolution), violins, bias=False)
                        self.track_plot_elements(data_label, 'periodic-violin', 'Median_plot_{}'.format(relevant_temporal_resolution), median_plots, bias=False)

            # periodic plot type
            else:

                # plot horizontal line/s across x axis at value/s of minimum experiment bias (if bias statistic)
                if z_statistic_sign == 'bias':
                    # get value/s of minimum bias for statistic
                    if z_statistic_type == 'basic':
                        minimum_bias = self.read_instance.basic_stats[base_zstat]['minimum_bias']
                    else:
                        minimum_bias = self.read_instance.expbias_stats[base_zstat]['minimum_bias']
                    bias_lines = []
                    for mb in minimum_bias:
                        bias_lines += [relevant_sub_ax.axhline(y=mb, **plot_characteristics['bias_line'])]
                    # track plot elements if using dashboard 
                    if (not self.read_instance.offline) and (not self.read_instance.interactive):
                        self.track_plot_elements('ALL', 'periodic', 'bias_line_{}'.format(relevant_temporal_resolution), 
                                                 bias_lines, bias=bias)

                # calculate statistic in each of periodic groups per data label
                if z_statistic_sign == 'bias':
                    if self.read_instance.observations_data_label in cut_data_labels:
                        cut_data_labels.remove(self.read_instance.observations_data_label)
                    stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstat, 
                                                    [self.read_instance.observations_data_label]*len(cut_data_labels), 
                                                    cut_data_labels, period=relevant_temporal_resolution)
                else:
                    stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, 
                                                     zstat, cut_data_labels, [], period=relevant_temporal_resolution)

                # iterate through data labels and plot violins
                for data_label_ii, data_label in enumerate(cut_data_labels): 

                    # skip observational array if bias stat
                    if (z_statistic_sign == 'bias') & (data_label == self.read_instance.observations_data_label):
                        continue

                    # make plot
                    self.periodic_plots = relevant_sub_ax.plot(xticks, stats_calc[:, data_label_ii], 
                                                               color=self.read_instance.plotting_params[data_label]['colour'], 
                                                               zorder=self.read_instance.plotting_params[data_label]['zorder'], 
                                                               **plot_characteristics['plot'])

                    # track plot elements if using dashboard 
                    if (not self.read_instance.offline) and (not self.read_instance.interactive):
                        self.track_plot_elements(data_label, 'periodic', 'plot_{}'.format(relevant_temporal_resolution), self.periodic_plots, bias=bias)

    def make_distribution(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options,
                          data_range_min=None, data_range_max=None, violin_resolution=None):
        """ Make distribution plot.

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param data_labels: Data arrays to plot
            :type data_labels: list
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
            :param data_range_min: Minimum data range of distribution plot grid 
            :type data_range_min: float
            :param data_range_max: Maximum data range of distribution plot grid 
            :type data_range_max: float
            :param violin_resolution: If are calculating distribution for violin plot, this is set to temporal resolution of groupings
            :type violin_resolution: int
        """

        # determine if 'bias' in plot_options
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # if 'obs' in plot_options, set data labels to just observations label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # set data ranges for distribution plot grid if not set explicitly
        if data_range_min is None:
            data_range_min = self.canvas_instance.selected_station_data_min[networkspeci]
        
        if data_range_max is None:
            data_range_max = self.canvas_instance.selected_station_data_max[networkspeci]

        # set xgrid for calculating distribution
        # if calculating for period n_samples is set to pdf_min_samples
        # otherwise it is inferred from data (if above minimum value)
        if violin_resolution is not None:
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
            if (not self.read_instance.offline) and (not self.read_instance.interactive):
                self.track_plot_elements('ALL', 'distribution', 'bias_line', bias_line, bias=bias)
            if self.read_instance.observations_data_label in cut_data_labels:
                cut_data_labels.remove(self.read_instance.observations_data_label)

        # if violin plot setup arrays for saving data
        if violin_resolution is not None:
            PDFs_sampled = np.full((len(cut_data_labels), len(self.canvas_instance.periodic_xticks[violin_resolution]), int(n_samples)), np.NaN, dtype=np.float32)

        # iterate through data labels
        for data_label_ii, data_label in enumerate(cut_data_labels):
        
            PDF_sampled_calculated = False

            # setup bias plot
            if bias:

                # calculate obs PDF on first pass
                if data_label_ii == 0:
                    kde_data_obs = drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(self.read_instance.observations_data_label),0,:])
                    
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

                        if PDF_obs_sampled is None:
                            msg = 'Warning: The kernel bandwidth is 0 for {}. '.format(data_label)
                            msg += 'The distribution plot will be created and not include data for this label. '
                            msg += 'To change the bandwith, we recommend increasing the number of '
                            msg += 'pdf_min_samples in the plot characteristics settings files.'
                            print(msg)
                            continue

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
                if PDF_model_sampled is None:
                    msg = 'Warning: The kernel bandwidth is 0 for {}. '.format(data_label)
                    msg += 'The distribution plot will be created and not include data for this label. '
                    msg += 'To change the bandwith, we recommend increasing the number of '
                    msg += 'pdf_min_samples in the plot characteristics settings files.'
                    print(msg)
                    continue
    
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
                if (violin_resolution is not None) & (data_label_ii == 0):
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
                if violin_resolution is not None:
                    kde_data_grouped = [drop_nans(group[valid_data_labels.index(data_label)].flatten()) for group in self.canvas_instance.selected_station_data[networkspeci][violin_resolution]['active_mode']]
                else:    
                    kde_data_grouped = [drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:])]

                # iterate through kde data groups
                for period_ii, kde_data in enumerate(kde_data_grouped):

                    # get relevant data ranges / x_grid, for violin period distribution calculation
                    if violin_resolution is not None:
                        data_range_min = period_data_range_min[period_ii]
                        data_range_max = period_data_range_max[period_ii]
                        x_grid = period_x_grid[period_ii]
                    
                    #filter out data outside data range bounds
                    #kde_data = kde_data[(kde_data > data_range_min) & (kde_data < data_range_max)]

                    # check if all values are equal in the dataframe
                    if kde_data.size == 0:
                        if violin_resolution is None:
                            print('Warning: The kernel density cannot be calculated because there are no valid values for {}.'.format(data_label))
                        continue
                    elif np.all(kde_data == kde_data[0]):
                        if violin_resolution is None:
                            print('Warning: The kernel density cannot be calculated because all {} values are equal.'.format(data_label))
                        continue
                    else:
                        PDF_sampled = kde_fft(kde_data, xgrid=x_grid)
                        if PDF_sampled is None:
                            msg = 'Warning: The kernel bandwidth is 0 for {}. '.format(data_label)
                            msg += 'The distribution plot will be created and not include data for this label. '
                            msg += 'To change the bandwith, we recommend increasing the number of '
                            msg += 'pdf_min_samples in the plot characteristics settings files.'
                            print(msg)
                            continue
                        
                        #PDF_fit = FFTKDE(kernel='gaussian', bw='scott').fit(kde_data)
                        #PDF_sampled = PDF_fit.evaluate(x_grid)
                        # save PDF for violin plot
                        if violin_resolution is not None:
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
                if (not self.read_instance.offline) and (not self.read_instance.interactive):
                    self.track_plot_elements(data_label, 'distribution', 'plot', self.distribution_plot, bias=bias)

        # if have made PDFs for violin plot then return it
        if violin_resolution is not None:
            return period_x_grid, PDFs_sampled

    def make_scatter(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options):
        """ Make scatter plot.

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param data_labels: Data arrays to plot
            :type data_labels: list
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
        """
        
        # if 'obs' in plot_options, set data labels to just observations data label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

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

        # skip making scatter for offline and interactive mode
        # we do not apply this in the dashboard to avoid being unable to see the points on certain changes
        if ((self.read_instance.offline) or (self.read_instance.interactive)) and ('hidedata' in plot_options):
            return
        
        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get observations data (flattened)
        observations_data = self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(self.read_instance.observations_data_label),0,:]

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
            if data_label == self.read_instance.observations_data_label:
                continue

            # get experiment data (flattened)
            experiment_data = self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:]

            # subset data if neccessary
            if subset:
                experiment_data = experiment_data[inds_subset]

            # get marker size (for offline and interactive)
            if (self.read_instance.offline) or (self.read_instance.interactive):
                self.get_markersize(relevant_axis, 'scatter', networkspeci, plot_characteristics, data=observations_data)

            # create scatter plot
            self.scatter_plot = relevant_axis.plot(observations_data, experiment_data, 
                                                   color=self.read_instance.plotting_params[data_label]['colour'],
                                                   **plot_characteristics['plot'])

            # track plot elements if using dashboard 
            if (not self.read_instance.offline) and (not self.read_instance.interactive):
                self.track_plot_elements(data_label, 'scatter', 'plot', self.scatter_plot, bias=False)

    def make_boxplot(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options):
        """ Make boxplot.

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param data_labels: Data arrays to plot
            :type data_labels: list
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
        """

        # if 'obs' in plot_options, set data labels to just observations data label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

        # if multispecies in plot options then make plot for all networkspecies
        if 'multispecies' in plot_options:
            networkspecies = self.read_instance.networkspecies
        else:
            networkspecies = [networkspeci]

        # if norm in plot options, then get factor for normalisation (per observations and experiment)
        if ('normmedian' in plot_options) or ('normmean' in plot_options) or ('normmin' in plot_options) or ('normmax' in plot_options) or ('normsum' in plot_options):

            norm_factor = {}

            # iterate through networkspecies
            for ns in networkspecies:
            
                # get valid data labels for networkspeci
                valid_data_labels = self.canvas_instance.selected_station_data_labels[ns]
                
                # cut data_labels for those in valid data labels
                cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

                # iterate through cut data labels making plot
                for data_label_ii, data_label in enumerate(cut_data_labels):

                    if data_label not in norm_factor:
                        norm_factor[data_label] = []

                    # get data (flattened and drop NaNs)
                    data_array = drop_nans(self.canvas_instance.selected_station_data[ns]['flat'][valid_data_labels.index(data_label),0,:])

                    #calculate values for normalisation
                    if 'normmedian' in plot_options:
                        norm_factor[data_label].append(np.nanmedian(data_array))
                    elif 'normmean' in plot_options:
                        norm_factor[data_label].append(np.nanmean(data_array))
                    elif 'normmin' in plot_options:
                        norm_factor[data_label].append(np.nanmin(data_array))    
                    elif 'normmax' in plot_options:
                        norm_factor[data_label].append(np.nanmax(data_array)) 
                    elif 'normsum' in plot_options:
                        norm_factor[data_label].append(np.nansum(data_array)) 

            # get ultimate factor per observations / experiments for normalisation
            for data_label in norm_factor:
                norm_factor[data_label] = sum(np.diff(bins)*norm_factor[data_label])
            
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
                if (('individual' in plot_options) or ('obs' in plot_options) 
                    or (len(self.read_instance.networkspecies) == 1) or (len(cut_data_labels) == 1)):
                    widths = plot_characteristics['group_widths']['singlespecies']
                else:
                    available_width = plot_characteristics['group_widths']['multispecies']
                    remainder_width = 1.0 - available_width
                    start_point = -0.5 + (remainder_width / 2.0)
                    widths = available_width / (len(cut_data_labels) + 0.15)
                    spacing = (available_width - (widths * len(cut_data_labels))) / (len(cut_data_labels) - 1)

                # get plot positions
                if ('individual' in plot_options) or ('obs' in plot_options) or (len(cut_data_labels) == 1):
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

                    # normalise data array
                    if ('normmedian' in plot_options) or ('normmean' in plot_options) or ('normmin' in plot_options) or ('normmax' in plot_options) or ('normsum' in plot_options):
                        data_array = data_array * norm_factor[data_label]

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
                    if (not self.read_instance.offline) and (not self.read_instance.interactive):
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
            xticks = positions
            xtick_labels = copy.deepcopy(cut_data_labels)

        #modify xticks to be horizontal as just have 1 label
        if len(xtick_labels) == 1:
            xtick_params['rotation'] = 0
            xticklabel_params = {}

        # set xticks / xticklabels
        relevant_axis.set_xticks(xticks)
        relevant_axis.xaxis.set_tick_params(**xtick_params)
        relevant_axis.set_xticklabels(xtick_labels, **xticklabel_params)

    def make_heatmap(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options,
                     zstat=None, subsection=None, plotting_paradigm=None, stats_df=None):
        """ Make heatmap plot.

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param data_labels: Data arrays to plot
            :type data_labels: list
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
            :param zstat: Statistic
            :type zstat: str
            :param subsection: Currently active subsection
            :type subsection: str
            :param plotting_paradigm: Plotting paradigm (summary or station in offline reports)
            :type plotting_paradigm: str
            :param stats_df: Dataframe of previously calculated statistics, default is None
            :param stats_df: pandas dataframe
        """

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
                if self.read_instance.observations_data_label in cut_data_labels:
                    cut_data_labels.remove(self.read_instance.observations_data_label)
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, [zstat], 
                                                [self.read_instance.observations_data_label]*len(cut_data_labels), 
                                                cut_data_labels)
            else:
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, [zstat], 
                                                 cut_data_labels, [])

            # create stats dataframe
            stats_df = pd.DataFrame(data=stats_calc, 
                                    index=cut_data_labels,
                                    dtype=np.float64)

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

        # determine if want to add annotations or not from plot_options
        if 'annotate' in plot_options:
            # get rounded labels
            decimal_places = plot_characteristics['round_decimal_places']['table']
            if Version(pd.__version__) >= Version("2.1.0"):
                annotate = stats_df.map(lambda x: round_decimal_places(x, decimal_places))
            else:
                annotate = stats_df.applymap(lambda x: round_decimal_places(x, decimal_places))
        else:
            annotate = False

        # plot heatmap
        heatmap = sns.heatmap(stats_df, 
                              ax=relevant_axis, 
                              annot=annotate,
                              fmt='',
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

        # set xticklables
        relevant_axis.set_xticklabels(stats_df.columns, **plot_characteristics['xticklabels'])

        # axis cuts off due to bug in matplotlib 3.1.1 - hack fix. Remove in Future!
        if len(stats_df.index) > 1:
            bottom, top = relevant_axis.get_ylim()
            relevant_axis.set_ylim(bottom + 0.5, top - 0.5)

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
                    
                    if Version(matplotlib.__version__) >= Version("3.3"):
                        relevant_axis.annotate(text=networkspeci_label, annotation_clip=False,
                                               xy=(plot_characteristics['multispecies']['xmin'], y_position), 
                                               **plot_characteristics['multispecies']['yticklabels'])
                    else:
                        relevant_axis.annotate(s=networkspeci_label, annotation_clip=False,
                                               xy=(plot_characteristics['multispecies']['xmin'], y_position), 
                                               **plot_characteristics['multispecies']['yticklabels'])        
        
        # format for non multispecies
        else:

            # vertically align yticklabels due to bug again in matplotlib - hack fix. Remove in Future!
            for tick in relevant_axis.get_yticklabels():
                tick.set_verticalalignment("center")

        # track plot elements if using dashboard 
        if (not self.read_instance.offline) and (not self.read_instance.interactive):
            self.track_plot_elements(self.read_instance.observations_data_label, 'heatmap', 'plot', heatmap, bias=bias)

    def make_table(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options,
                   zstats=None, statsummary=False, subsection=None, plotting_paradigm=None, stats_df=None):
        """ Make table plot.

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param data_labels: Data arrays to plot
            :type data_labels: list
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
            :param zstats: Statistics
            :type zstats: list
            :param statsummary: To indicate if making alternative statistical summary table plot  
            :type statsummary: boolean
            :param subsection: Currently active subsection
            :type subsection: str
            :param plotting_paradigm: Plotting paradigm (summary or station in offline reports)
            :type plotting_paradigm: str
            :param stats_df: Dataframe of previously calculated statistics, default is None
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
                if self.read_instance.observations_data_label in cut_data_labels:
                    cut_data_labels.remove(self.read_instance.observations_data_label)
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstats, 
                                                [self.read_instance.observations_data_label]*len(cut_data_labels), 
                                                cut_data_labels)
            else:
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstats, 
                                                 cut_data_labels, [])

            # create stats dataframe
            stats_df = pd.DataFrame(data=stats_calc, 
                                    index=cut_data_labels,
                                    dtype=np.float64)
        
        # when we have 1 stat in the statsummary, the column name is 0
        # we need to rename it to the stat name
        if (len(stats_df.columns) == 1) and (stats_df.columns[0] == 0):
            stats_df = stats_df.rename(columns={stats_df.columns[0]: zstats[0]})
        
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

        # round dataframe
        decimal_places = plot_characteristics['round_decimal_places']['table']
        if Version(pd.__version__) >= Version("2.1.0"):
            stats_df = stats_df.map(lambda x: round_decimal_places(x, decimal_places))
        else: 
            stats_df = stats_df.applymap(lambda x: round_decimal_places(x, decimal_places))

        # offline reports
        if (self.read_instance.offline) or (self.read_instance.interactive):
            
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
            if (self.read_instance.offline) or (self.read_instance.interactive):
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
            else:
                empty_cells = len(stats_df.columns) - len(data_labels)
                col_labels = ['']*empty_cells + data_labels

        # dashboard
        else:
            # there is only statsummary
            if statsummary:

                # get labels
                data_labels = list(stats_df.index)
                stats = list(stats_df.columns)
                
                # reset index
                stats_df = stats_df.reset_index()

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
                                if data_label == self.read_instance.observations_data_label:
                                    color = 'white'
                                # experiments in legend colors
                                else:
                                    color = self.read_instance.plotting_params[data_label]['colour']
                                cell_colours[col].append(color)
                        # white for other cells
                        else:
                            cell_colours[col] = ['white'] * stats_df.shape[0]
                    if stats_df.shape == (1, 1):
                        plot_characteristics['plot']['cellColours'] = np.array(cell_colours, dtype=object)
                    else:
                        plot_characteristics['plot']['cellColours'] = np.array(cell_colours, dtype=object).T
        else:
            if 'col_colours' in plot_characteristics:
                if plot_characteristics['col_colours']:
                    col_colours = []
                    for data_label in data_labels:
                        # observations in white
                        if data_label == self.read_instance.observations_data_label:
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
        if (self.read_instance.offline) or (self.read_instance.interactive):
            column_ii = 0
            for column, rows in zip(['networkspecies', 'subsections'], (networkspecies, subsections)):
                if column in stats_df.columns:
                    # count consecutive duplicates
                    count_dups = [sum(1 for _ in group) for _, group in groupby(rows)]
                    
                    # merge cells that have consecutive duplicates
                    current_row = 0
                    for count_ii, count in enumerate(count_dups):
                        cells_to_merge = [(current_row + i, column_ii) for i in range(1, count+1)]
                        merge_cells(table, cells_to_merge)
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
        if (not self.read_instance.offline) and (not self.read_instance.interactive):
            if statsummary:
                self.track_plot_elements(self.read_instance.observations_data_label, 'statsummary', 'plot', [table], bias=bias)
            else:
                self.track_plot_elements(self.read_instance.observations_data_label, 'table', 'plot', [table], bias=bias)
    
    def make_taylor(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options, zstat,
                    stddev_max=None):
        """ Make Taylor diagram plot.
            Reference: https://gist.github.com/ycopin/3342888.

            See explanation of calculations here:
            https://waterprogramming.wordpress.com/2020/12/22/taylor-diagram/

            :param relevant_axis: Axis to plot on 
            :type relevant_axis: object
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param data_labels: Data arrays to plot
            :type data_labels: list
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param plot_options: Options to configure plot  
            :type plot_options: list
            :param zstat: Statistic
            :type zstat: str
            :param stddev_max: Maximum standard deviation
            :type stddev_max: float
        """

        if (self.read_instance.offline) or (self.read_instance.interactive):
            self.taylor_polar_relevant_axis = relevant_axis.get_aux_axes(
                PolarAxes.PolarTransform(apply_theta_transforms=False))

        # calculate statistics
        stats_dict = {}

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get data labels without observations
        obs_index = cut_data_labels.index(self.read_instance.observations_data_label)
        data_labels_sans_obs = copy.deepcopy(cut_data_labels)
        if self.read_instance.observations_data_label in data_labels_sans_obs:
            data_labels_sans_obs.remove(self.read_instance.observations_data_label)

        # standard deviation - absolute calculations of observations and models
        stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, 'StdDev', 
                                         data_labels, [])
        stats_dict['StdDev'] = stats_calc   

        # correlation - between observations and model
        stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstat, 
                                         [self.read_instance.observations_data_label]*len(data_labels_sans_obs), 
                                         data_labels_sans_obs)
        stats_calc = np.insert(stats_calc, obs_index, np.NaN)
        stats_dict[zstat] = stats_calc
        
        # get maximum stddev in dataframe (if not defined)
        if stddev_max is None:
            stddev_max = np.nanmax(stats_dict["StdDev"])

        # create stats dataframe
        stats_df = pd.DataFrame(data=stats_dict, 
                                index=data_labels,
                                dtype=np.float64)

        # get labels
        rlabel = stats_df.columns[1]
        xylabel = stats_df.columns[0]

        # check if we need to extend the axis to negative correlations
        extend = False
        if np.nanmin(stats_df[rlabel]) < 0:
            extend = True

        # update axis extremes and labels
        tmin, tmax, smin, smax, gl1, tf1 = get_taylor_diagram_ghelper_info(stddev_max, 
                                                                           plot_characteristics,
                                                                           extend)
        relevant_axis.get_grid_helper().update_grid_finder(
            extreme_finder=fa.ExtremeFinderFixed((tmin, tmax, smin, smax)),
            grid_locator1=gl1, tick_formatter1=tf1)

        # update axis position and size in dashboard
        if (not self.read_instance.offline) and (not self.read_instance.interactive):

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
        reference_stddev = stats_df[xylabel].iloc[obs_index]
        num_levels = plot_characteristics['contours']['levels']['number']
        rs, ts = np.meshgrid(np.linspace(smin, smax), np.linspace(0, tmax))
        rms = np.sqrt(reference_stddev**2 + rs**2 - 2*reference_stddev*rs*np.cos(ts))
        contours = self.taylor_polar_relevant_axis.contour(ts, rs, rms, num_levels,
            **plot_characteristics['contours']['style']['general'])

        # add contour labels
        self.taylor_polar_relevant_axis.clabel(contours, contours.levels, inline=True, fmt = "%.2f", fontsize=6)

        # add reference contour of observational standard deviation
        ref_x = np.linspace(0, tmax)
        ref_y = np.zeros_like(ref_x) + reference_stddev
        self.taylor_polar_relevant_axis.plot(ref_x, ref_y, **plot_characteristics['contours']['style']['obs'])

        # set up annotations
        if (not self.read_instance.offline) and (not self.read_instance.interactive):
            # create annotation on hover
            annotation = HoverAnnotation(self.canvas_instance, 
                                        'taylor', 
                                        self.taylor_polar_relevant_axis,
                                        plot_characteristics, 
                                        add_vline=False)
            self.canvas_instance.annotations['taylor'] = annotation.annotation
            self.canvas_instance.annotations_lock['taylor'] = False

            # connect axis to hover function
            self.canvas_instance.figure.canvas.mpl_connect('motion_notify_event', 
                lambda event: annotation.hover_annotation(event, 'taylor'))

        # add models
        for data_label, stddev, corr_stat in zip(stats_df.index, 
                                                 stats_df[xylabel], 
                                                 stats_df[rlabel]):
            if data_label == self.read_instance.observations_data_label:
                continue
            self.taylor_plot = self.taylor_polar_relevant_axis.plot(np.arccos(corr_stat), stddev,
                                                                    **plot_characteristics['plot'],
                                                                    mfc=self.read_instance.plotting_params[data_label]['colour'], 
                                                                    mec=self.read_instance.plotting_params[data_label]['colour'],
                                                                    label=data_label) 

            # track plot elements if using dashboard 
            if (not self.read_instance.offline) and (not self.read_instance.interactive):
                self.track_plot_elements(data_label, 'taylor', 'plot', self.taylor_plot, bias=False)

    def make_fairmode_target(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options):
        
        # add target
        main_circle = plt.Circle(**plot_characteristics['auxiliar']['circle']['main'])
        relevant_axis.add_patch(main_circle)

        # add a black circle with radius 1 (continuous line)
        big_circle = plt.Circle(**plot_characteristics['auxiliar']['circle']['big'])
        relevant_axis.add_patch(big_circle)

        # add a black circle with radius 0.5 (dotted line)
        small_circle = plt.Circle(**plot_characteristics['auxiliar']['circle']['small'])
        relevant_axis.add_patch(small_circle)

        # add text in the sides
        relevant_axis.text(**plot_characteristics['auxiliar']['sides']['top'], transform=relevant_axis.transAxes)
        relevant_axis.text(**plot_characteristics['auxiliar']['sides']['bottom'], transform=relevant_axis.transAxes)
        relevant_axis.text(**plot_characteristics['auxiliar']['sides']['left'], transform=relevant_axis.transAxes)
        relevant_axis.text(**plot_characteristics['auxiliar']['sides']['right'], transform=relevant_axis.transAxes)

        # add diagonal lines (y = x and y = -x)
        xmin = np.min(plot_characteristics['xticks']['ticks'])
        xmax = np.max(plot_characteristics['xticks']['ticks'])
        ymin = np.min(plot_characteristics['yticks']['ticks'])
        ymax = np.max(plot_characteristics['yticks']['ticks'])
        relevant_axis.plot([xmin, xmax], [ymin, ymax], **plot_characteristics['auxiliar']['crosses']['increasing'])
        relevant_axis.plot([xmin, xmax], [ymax, ymin], **plot_characteristics['auxiliar']['crosses']['decreasing'])

        # get data for fairmode
        fairmode_settings = yaml.safe_load(
            open(os.path.join(PROVIDENTIA_ROOT, 
                              'settings/fairmode.yaml')))

        speci = networkspeci.split('|')[1]
        speci_settings = fairmode_settings[speci]
        u_95r_RV = speci_settings['u_95r_RV']
        RV = speci_settings['RV']
        alpha = speci_settings['alpha']
        beta = speci_settings['beta']
        coverage = speci_settings['coverage']
        
        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get observations data
        observations_data = self.canvas_instance.selected_station_data[networkspeci]['per_station'][0,:,:]

        # TODO: Resample to daily for PM, calculate MDA8 for ozone, add flag to make sure resolution is hourly
        # TODO: Make sure days with less than 75% coverage are nan (avoid NO2)

        # drop stations that have a coverage of less than coverage
        obs_representativity = Stats.calculate_data_avail_fraction(observations_data)
        valid_station_idxs = obs_representativity >= coverage
        observations_data = observations_data[valid_station_idxs, :]

        # get station references without nans
        # if for all timesteps it is nan, set nan as station reference
        station_references = self.canvas_instance.selected_station_metadata[networkspeci][
                                'station_reference'][valid_station_idxs, :]
        valid_station_references = []
        for station_reference in station_references:
            first_valid_station_reference = next((
                reference for reference in station_reference 
                if reference == reference), 'nan')
            valid_station_references.append(first_valid_station_reference)

        n_stations = len(valid_station_references)

        # get station area classification without nans
        # if for all timesteps it is nan, set nan as station reference
        station_area_classifications = self.canvas_instance.selected_station_metadata[networkspeci][
                                        'area_classification'][valid_station_idxs, :]
        valid_station_area_classifications = []
        for station_area_classification in station_area_classifications:
            first_valid_station_area_classification = next((
                classification for classification in station_area_classification 
                if classification == classification), 'nan')
            valid_station_area_classifications.append(first_valid_station_area_classification)

        # initialise annotation text
        if 'annotate' in plot_options:
            annotate_text = f"α={alpha}\n"
            annotate_text += f"β={beta}\n"
            annotate_text += f"RV={RV}\n"
            annotate_text += f"U₉₅,ᵣᴿⱽ={u_95r_RV}\n\n\n"
            annotate_text += f"{n_stations} stations with\ncoverage above {coverage}%\n\n\n"
        
        # iterate through data labels
        for data_label in cut_data_labels:

            # continue for observations data label
            if data_label == self.read_instance.observations_data_label:
                continue

            # get experiment data
            experiment_data = self.canvas_instance.selected_station_data[networkspeci]['per_station'][valid_data_labels.index(data_label),valid_station_idxs,:]

            # calculate MQI for the current station
            x_points = []
            y_points = []
            stations = []
            bad_stations = []
            area_classifications = []

            # get FAIRMODE statistics per station
            mqi_array = np.full(n_stations, np.nan)
         
            for station_idx, (station, area_classification) in enumerate(
                zip(valid_station_references, valid_station_area_classifications)):

                st_observations_data = observations_data[station_idx, :]
                st_experiment_data = experiment_data[station_idx, :]
                
                x, y, mqi = ExpBias.calculate_fairmode_target_stats(st_observations_data, st_experiment_data, 
                                                                    u_95r_RV, RV, alpha, beta)

                if station != np.nan:
                    x_points.append(x)
                    y_points.append(y)
                    mqi_array[station_idx] = mqi
                    stations.append(station)
                    area_classifications.append(area_classification)

            # plot data
            for x, y, station, mqi, area_classification in (zip(
                x_points, y_points, stations, mqi_array, area_classifications)):
                if area_classification not in plot_characteristics['area_classification']['markers']:
                    marker = 'h'
                else:
                    marker = plot_characteristics['area_classification']['markers'][area_classification]
                self.fairmode_plot = relevant_axis.plot(x, y, markeredgecolor=self.read_instance.plotting_params[data_label]['colour'], 
                                                        marker=marker, **plot_characteristics['plot'])
                if mqi > 1:
                    bad_stations.append(station)

                # track plot elements if using dashboard 
                if (not self.read_instance.offline) and (not self.read_instance.interactive):
                    self.track_plot_elements(data_label, 'fairmode-target', 'plot', self.fairmode_plot, bias=False)

            if 'annotate' in plot_options:

                # calculate MQI90
                mqi_sorted = sorted(mqi_array[~np.isnan(mqi_array)])
                i_90 = int(0.9 * len(mqi_sorted)) - 1
                MQI90 = mqi_sorted[i_90]
                MQI90_formatted = f"{MQI90:.2f}"

                # add MQI90
                annotate_text += f"{data_label}\n"
                annotate_text += f"MQI₉₀ = {MQI90_formatted}\n"

                # add bad stations
                if len(bad_stations) == 1:
                    stations_name = 'station'
                else:
                    stations_name = 'stations'
                annotate_text += f'{len(bad_stations)} {stations_name} with MQI > 1'
                if len(bad_stations) > 0:
                    annotate_text += f': {bad_stations}'
                if data_label != cut_data_labels[-1]:
                    annotate_text += '\n\n'

        # add annotation text
        if 'annotate' in plot_options:
            relevant_axis.text(2.2, 0, annotate_text, **plot_characteristics['annotate_text'])
            
        # update axis labels
        relevant_axis.set_xticks(**plot_characteristics['xticks'])
        relevant_axis.set_yticks(**plot_characteristics['yticks'])

        # create legend
        legend_elements = []
        for area_classification in np.unique(area_classifications):
            if area_classification not in plot_characteristics['area_classification']['markers']:
                marker = 'h'
            else:
                marker = plot_characteristics['area_classification']['markers'][area_classification]
            legend_element = mlines.Line2D([], [], marker=marker, 
                                           label=area_classification, 
                                           **plot_characteristics['area_classification']['plot'])
            legend_elements.append(legend_element)
        relevant_axis.legend(handles=legend_elements, 
                             **plot_characteristics['area_classification']['legend'])

    def track_plot_elements(self, data_label, base_plot_type, element_type, plot_object, bias=False):
        """ Function that tracks plotted lines and collections
            that will be removed/added when picking up legend elements on dashboard.

            :param data_label: Data array to plot
            :type data_label: str
            :param base_plot_type: Plot type, without statistical information
            :type base_plot_type: str
            :param element_type: Element type
            :type element_type: str
            :param plot_object: Plotted element object
            :type plot_object: object
            :param bias: Indicates if plot is a bias plot
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
        # do not save elements for plot objecplot type, without statistical information
        # all other plot elements
        else:
            # add list of lines
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object
            
        # set element visibility
        if (data_label not in self.canvas_instance.plot_elements['data_labels_active']) & (data_label != 'ALL'):
            for element in self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type]:
                element.set_visible(False)

    def get_markersize(self, relevant_axis, base_plot_type, networkspeci, plot_characteristics, 
                       data=None, active_map_valid_station_inds=[]):
        """ Set markersize for plot.
        
            :param base_plot_type: Plot type, without statistical information
            :type base_plot_type: str
            :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
            :type networkspeci: str
            :param plot_characteristics: Plot characteristics  
            :type plot_characteristics: dict
            :param data: Data array to be plotted
            :type data: numpy array
            :param active_map_valid_station_inds: Valid map indices to plot
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

                # add to plot_characteristics yaml
                plot_characteristics['plot']['markersize'] = markersize
        
        elif base_plot_type == 'map':

            if plot_characteristics['plot']['s'] == '':

                # calculate marker size considering points density
                n_points = len(self.read_instance.station_longitudes[networkspeci][active_map_valid_station_inds])
                
                # calculate figure area and density
                # divide area by 1000 so the function below makes sense
                area = (relevant_axis.bbox.width * relevant_axis.bbox.height) / 1000
                density = n_points / area

                # marker size is calculated using an exponential equation
                # the maximum size is 40 (very low densities)
                # see https://earth.bsc.es/gitlab/ac/Providentia/-/issues/210
                plot_characteristics['plot']['s'] = 1.2**(-density)*40
