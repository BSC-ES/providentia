""" Class to create interactive session"""

import copy
import datetime
import json
import os
import sys
import yaml

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.projections import PolarAxes
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as fa
import numpy as np
import pandas as pd

from .configuration import load_conf
from .configuration import ProvConfiguration
from .fields_menus import (init_metadata, init_period, init_representativity, metadata_conf,
                           update_metadata_fields, update_period_fields, update_representativity_fields,
                           period_conf, representativity_conf)
from .filter import DataFilter
from .plot import Plot
from .plot_aux import get_taylor_diagram_ghelper
from .plot_formatting import (format_plot_options, format_axis, set_axis_label, set_axis_title, 
                              harmonise_xy_lims_paradigm)
from .read import DataReader
from .read_aux import (generate_file_trees, get_lower_resolutions, 
                       get_nonrelevant_temporal_resolutions, get_relevant_temporal_resolutions, 
                       get_valid_experiments, get_valid_obs_files_in_date_range)
from .statistics import (calculate_statistic, generate_colourbar, generate_colourbar_detail, 
                         get_fairmode_data, get_selected_station_data, get_z_statistic_info)
from .writing import export_configuration, export_data_npz, export_netcdf

from providentia.auxiliar import CURRENT_PATH, join, expand_plot_characteristics

PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
fairmode_settings = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/fairmode.yaml')))

# do not print deprecated warnings
import warnings
warnings.filterwarnings("ignore")

# determine if using jupyter notebook or not
try:
    __IPYTHON__
    jupyter_session = True
except NameError:
    jupyter_session = False

class Interactive:
    """ Class to create interactive Providentia session."""

    def __init__(self, **kwargs):

        self.kwargs = kwargs

        # set mode as interactive
        self.kwargs['interactive'] = True

        # load statistical yamls
        self.basic_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/basic_stats.yaml')))
        self.expbias_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/experiment_bias_stats.yaml')))

        # load representativity information
        self.representativity_info = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/internal/representativity.yaml')))

        # set configuration variables, as well as any other defined variables
        valid_config = self.set_config(**self.kwargs)
        # if configuration read was not valid then return here
        if not valid_config:
            return

        # generate file trees if needed
        generate_file_trees(self)
        
        # initialise DataReader class
        self.datareader = DataReader(self)

        # check for self defined plot characteristics file
        if self.plot_characteristics_filename == '':
            if self.tests:
                mode = 'tests'
            else:
                mode = 'interactive'
            self.plot_characteristics_filename = join(PROVIDENTIA_ROOT, 'settings/plot_characteristics.yaml')
        plot_characteristics = yaml.safe_load(open(self.plot_characteristics_filename))
        self.plot_characteristics_templates = expand_plot_characteristics(plot_characteristics, mode)

        # initialise Plot class
        self.plot = Plot(read_instance=self, canvas_instance=self)

        # add general plot characteristics to self
        for k, val in self.plot_characteristics_templates['general'].items():
            setattr(self, k, val)

        # set some key configuration variables
        self.relevant_temporal_resolutions = get_relevant_temporal_resolutions(self.resolution)
        self.nonrelevant_temporal_resolutions = get_nonrelevant_temporal_resolutions(self.resolution)
        self.data_labels = [self.observations_data_label] + list(self.experiments.values())
        self.data_labels_raw = [self.observations_data_label] + list(self.experiments.keys())
        self.networkspecies = ['{}|{}'.format(network,speci) for network, speci in zip(self.network, self.species)]

        # get valid observations in date range
        get_valid_obs_files_in_date_range(self, self.start_date, self.end_date)

        # update available experiments for selected fields
        get_valid_experiments(self, self.start_date, self.end_date, self.resolution,
                              self.network, self.species)

        # read data
        self.read()  
        if self.invalid_read:
            print('No valid data to read')
            return

        # filter
        self.reset_filter(initialise=True)

        # set variable to know if data is in intial state or not
        self.initialised = True

    def read(self):
        """ Wrapper method to read data. """

        print('Reading data')

        # read data
        self.datareader.read_setup(['reset'])

    def filter(self):
        """ Wrapper method to filter data. """

        print('Filtering data')

        # filter data
        DataFilter(self)

        # get selected station data
        get_selected_station_data(read_instance=self, canvas_instance=self, networkspecies=self.networkspecies)

    def reset_filter(self, initialise=False):
        """ Wrapper method to reset filter data.

        :param initialise: Indicates whether to reset data to initial state when class was initialised 
        :type initialise: boolean, optional
        """

        print(f'Resetting filter for {self.subsection}')
   
        # initialise structures to store fields        
        init_representativity(self)
        init_period(self)
        init_metadata(self)

        # update available fields
        update_representativity_fields(self)
        update_period_fields(self)

        # for non-GHOST delete valid station indices variables because we do not want to 
        # remove the stations with 0 valid measurements before the filter has been updated, 
        # this will happen later
        if hasattr(self, 'valid_station_inds') and (not self.reading_ghost):
            delattr(self, 'valid_station_inds')
            delattr(self, 'valid_station_inds_temporal_colocation')

        update_metadata_fields(self)
        
        # apply set fields at initalisation for filtering
        if initialise:
            representativity_conf(self)
            period_conf(self)
            metadata_conf(self)
            
        # re-filter 
        self.filter()

        # set variable to know if data is in intial state or not
        if initialise:
            self.initialised = True
        else:
            self.initialised = False

    def make_plot(self, plot, data_labels=None, labela='', labelb='', title=None, xlabel=None, ylabel=None, 
                  cb=True, legend=True, set_obs_legend=True, map_extent=None, annotate=False, bias=False, 
                  domain=False, hidedata=False, logx=False, logy=False, multispecies=False, regression=False, 
                  smooth=False, threshold=False, plot_options=None, save=False, return_plot=False, format=None,
                  width=None, height=None):
        """ Wrapper method to make a Providentia plot.

        :param plot: Plot type
        :type plot: str
        :param data_labels: Data arrays to plot, defaults to None
        :type data_labels: list, optional
        :param labela: Label of first dataset, defaults to ''
        :type labela: str, optional
        :param labelb: Label of second dataset (if defined then a bias plot is made), defaults to ''
        :type labelb: str, optional
        :param title: Axes title, defaults to None
        :type title: str, optional
        :param xlabel: Label on x axes, defaults to None
        :type xlabel: str, optional
        :param ylabel: Label on y axes, defaults to None
        :type ylabel: str, optional
        :param cb: Indicates if colorbar appears on plot, defaults to True
        :type cb: bool, optional
        :param legend: Indicates if legend appears on plot, defaults to True
        :type legend: bool, optional
        :param set_obs_legend: Indicates if observations appear on legend, defaults to True
        :type set_obs_legend: bool, optional
        :param map_extent: Map extent, defaults to None
        :type map_extent: list, optional
        :param annotate: Indicates if there are annotations, defaults to False
        :type annotate: bool, optional
        :param bias: Indicates if data is biased, defaults to False
        :type bias: bool, optional
        :param domain: Indicates if domain shows in maps, defaults to False
        :type domain: bool, optional
        :param hidedata: Indicates if data points are hidden in plot, defaults to False
        :type hidedata: bool, optional
        :param logx: Indicates if the scale of the x axis is log, defaults to False
        :type logx: bool, optional
        :param logy: Indicates if the scale of the y axis is log, defaults to False
        :type logy: bool, optional
        :param multispecies: Indicates if plot has multispecies, defaults to False
        :type multispecies: bool, optional
        :param regression: Indicates if scatter plot has regression line/s, defaults to False
        :type regression: bool, optional
        :param smooth: Indicates if timeseries has smooth line/s, defaults to False
        :type smooth: bool, optional
        :param threshold: Indicates if plot has threshold line/s, defaults to False
        :type threshold: bool, optional
        :param plot_options: List with plot options, defaults to None
        :type plot_options: list, optional
        :param save: Indicates if you want to save the figure, defaults to False
        :type save: bool, optional
        :param return_plot: Indicates if you want to get the axes, defaults to False
        :type return_plot: bool, optional
        :param format: Format to overwrite the plots format taken from plot characteristics 
        :type format: dict, optional
        :param width: Figure width
        :type format: int, float, optional
        :param height: Figure height
        :type format: int, float, optional
        :return: matplotlib.axes._axes.Axes
        :rtype: Plot axes
        :return: List of statistics (used in statsummary)
        :rtype: list
        """

        # close any previously open figures
        plt.close()

        # define default argument mutables
        if data_labels is None:
            data_labels = []
        if plot_options is None:
            plot_options = []
        if format is None:
            format = {}

        # if any of plot options are given via keywords, put them in a list (with other passed plot options)
        if annotate:
            if 'annotate' not in plot_options:
                plot_options.append('annotate')
            # if passed argument is a list, then use that for stat list (if valid)
            if type(annotate) == list:
                annotation_stats = copy.deepcopy(annotate)
        if bias:
            if 'bias' not in plot_options:
                plot_options.append('bias')
        if domain:
            if 'domain' not in plot_options:
                plot_options.append('domain')
        if hidedata:
            if 'hidedata' not in plot_options:
                plot_options.append('hidedata')
        if logx:
            if 'logx' not in plot_options:
                plot_options.append('logx')
        if logy:
            if 'logy' not in plot_options:
                plot_options.append('logy')
        if multispecies:
            if 'multispecies' not in plot_options:
                plot_options.append('multispecies')
        if regression:
            if 'regression' not in plot_options:
                plot_options.append('regression')
        if smooth:
            if 'smooth' not in plot_options:
                plot_options.append('smooth')
            # if passed argument is a str/int/float, then use that for smoothing window
            if type(smooth) == str:
                try:
                    smooth = int(smooth)
                except:
                    pass
            if (type(smooth) == int) or (type(smooth) == float):
                smooth_window = int(smooth)
        if threshold:
            if 'threshold' not in plot_options:
                plot_options.append('threshold')

        # get base plot type (no plot options), and plot type (with plot options)
        base_plot_type = copy.deepcopy(plot)
        if len(plot_options) > 0:
            plot_type = '{}_{}'.format(base_plot_type, '_'.join(plot_options))
        else:
            plot_type = copy.deepcopy(plot)

        # get zstat for required plots
        base_plot_type_split = base_plot_type.split('-')
        if (len(base_plot_type_split) > 1) & (base_plot_type not in ['periodic-violin', 'fairmode-target']):
            base_plot_type = base_plot_type_split[0]
            zstat = base_plot_type_split[1]
        else:
            zstat = None

        # get networkspeci to plot (for non-multispecies plots), taking first one preferentially
        networkspeci = self.networkspecies[0]
        speci = networkspeci.split('|')[-1]

        # for timeseries chunking
        chunk_stat = None
        chunk_resolution = None
        if base_plot_type == 'timeseries':
            if zstat is not None:
                # get chunk statistic and resolution
                chunk_stat = copy.deepcopy(zstat)
                chunk_resolution = plot_type.split('-')[2].split('_')[0]
                
                # get zstat information 
                zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=chunk_stat)
                
                # check if chunk resolution is available
                if self.resampling_resolution is None:
                    available_timeseries_chunk_resolutions = list(get_lower_resolutions(self.resolution))
                else:
                    available_timeseries_chunk_resolutions = list(get_lower_resolutions(self.resampling_resolution))

                # show warning if it is not
                if chunk_resolution not in available_timeseries_chunk_resolutions:
                    msg = f'Warning: {plot_type} cannot be created because {chunk_resolution} '
                    msg += 'is not an available chunking resolution.'
                    if len(available_timeseries_chunk_resolutions) > 0:
                        msg += f'The available resolutions are: {available_timeseries_chunk_resolutions}'
                    print(msg)
                    return

                # show warning if chunk stat is NStations and mode is not Temporal|Spatial
                if chunk_stat == 'NStations' and self.statistic_mode != 'Temporal|Spatial':
                    msg = f'Warning: {plot_type} cannot be created because {chunk_stat} '
                    msg += 'it is only available when Temporal|Spatial mode is active.'
                    print(msg)
                    return
                    
        # get zstat information 
        zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type=plot_type) 

        # if only 1 label passed for map plot, and stat is a bias statistic then throw error
        if (base_plot_type == 'map') & (z_statistic_sign == 'bias') & (labelb == ''):
            print("Warning: Plotting a bias statistic, and only 1 label is set. Not making plot.")
            return
        
        # if bias and threshold plots are in plot options throw error
        if ('bias' in plot_options) & ('threshold' in plot_options):
            print("Warning: Cannot make a bias plot showing threshold lines. Not making plot.")
            return

        # do not make plot if hidedata is active but smooth is not in plot options
        if (base_plot_type == 'timeseries') and ('hidedata' in plot_options) and ('smooth' not in plot_options):
            msg = f"Warning: Cannot make {plot_type} because 'hidedata' plot option is set for "
            msg += "timeseries plot, but 'smooth' is not active. Not making plot."
            print(msg)
            return
        
        # do not make plot if hidedata is active but regression is not in plot options
        if (base_plot_type == 'scatter') and ('hidedata' in plot_options) and ('regression' not in plot_options):
            msg = f"Warning: Cannot make {plot_type} because 'hidedata' plot option is set for "
            msg += "scatter lot, but 'regression' is not active. Not making plot."
            print(msg)
            return
        
        # do not make Taylor diagram if statistic is not r or r2
        if (base_plot_type == 'taylor') and (zstat not in ['r', 'r2']):
            msg = f"Warning: Cannot make {plot_type} because statistic is not available or defined. "
            msg += "Choose between 'taylor-r' or 'taylor-r2'. Not making plot."
            print(msg)
            return
        
        # do not make FAIRMODE target plot if species not in list or resolution not hourly
        if base_plot_type == 'fairmode-target':
            if speci not in ['sconco3', 'sconcno2', 'pm10', 'pm2p5']:
                print(f'Warning: Fairmode target plot cannot be created for {speci}.')
                return
            if ((speci in ['sconco3', 'sconcno2'] and self.resolution != 'hourly') 
                or (speci in ['pm10', 'pm2p5'] and (self.resolution not in ['hourly', 'daily']))):
                print('Warning: Fairmode target plot can only be created if the resolution is hourly (O3, NO2, PM2.5 and PM10) or daily (PM2.5 and PM10).')
                return
            
            # skip making plot if there is no valid data
            data, valid_station_idxs = get_fairmode_data(self, self, networkspeci, self.resolution, self.data_labels)
            if not any(valid_station_idxs):
                print(f'No data after filtering by coverage for {speci}.')
                return

        # get data labels for plot
        if len(data_labels) == 0:
            data_labels = copy.deepcopy(self.data_labels)
        # if any passed data labels are not available then pass warning
        else:
            invalid_data_labels = [data_label for data_label in data_labels if data_label not in self.data_labels]
            data_labels = [data_label for data_label in data_labels if data_label in self.data_labels]
            if len(data_labels) == 0:
                print("Warning: None of the passed data labels are available. Not making plot.")
                return
            elif len(invalid_data_labels) > 0:
                print("Warning: Passed data labels {} are not available.".format(invalid_data_labels))

        # set plot characteristics
        self.plot_characteristics = dict()
        valid_plot = self.plot.set_plot_characteristics([plot_type], format=format, data_labels=data_labels)

        # if after setting plot charateristics it has been determined plot is not valid, then return
        if not valid_plot:
            return

        # adjust plot option attributes if passed
        if ('annotation_stats' in locals()) & ('annotate_stats' in self.plot_characteristics[plot_type]):
            self.plot_characteristics[plot_type]['annotate_stats'] = annotation_stats
        if ('smooth_window' in locals()) & ('smooth' in self.plot_characteristics[plot_type]):
            self.plot_characteristics[plot_type]['smooth']['window'] = smooth_window

        # if map extent passed not passed as argument, set it as value from .conf in memory
        if (not map_extent) and (self.map_extent):
            map_extent = copy.deepcopy(self.map_extent)

        # create figure
        if (width is not None) and (height is not None):
            fig = plt.figure(figsize=(width, height))
        else:
            print("Warning: Width and/or height have not been passed. The default values will be set.")
            fig = plt.figure(figsize=self.plot_characteristics[plot_type]['figsize'])

        # create axes
        if base_plot_type == 'map':
            ax = fig.add_subplot(111, projection=self.plotcrs)
        elif base_plot_type == 'taylor':            
            reference_stddev = 7.5
            ghelper = get_taylor_diagram_ghelper(reference_stddev, self.plot_characteristics[plot_type])
            ax = fig.add_subplot(111, axes_class=fa.FloatingAxes, grid_helper=ghelper)
        else:
            ax = fig.add_subplot(111)

        if base_plot_type in ['periodic', 'periodic-violin']:
            gs = gridspec.GridSpecFromSubplotSpec(100, 100, subplot_spec=ax.get_subplotspec())
            grid_dict = dict()
            grid_dict['hour'] = fig.add_subplot(gs[:46, :])
            grid_dict['dayofweek'] = fig.add_subplot(gs[54:, 64:])
            grid_dict['month'] = fig.add_subplot(gs[54:, :62])
            ax.spines['top'].set_color('none')
            ax.spines['bottom'].set_color('none')
            ax.spines['left'].set_color('none')
            ax.spines['right'].set_color('none')
            ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
            relevant_ax = grid_dict
        else:
            relevant_ax = ax

        # adjust margings and subplot spacing if defined
        if 'subplots_adjust' in self.plot_characteristics[plot_type]:
            fig.subplots_adjust(**self.plot_characteristics[plot_type]['subplots_adjust'])

        # get plotting function
        if base_plot_type == 'statsummary':
            func = getattr(self.plot, 'make_table')
        elif base_plot_type == 'fairmode-target':
            func = getattr(self.plot, 'make_fairmode_target')
        elif base_plot_type != 'legend':
            func = getattr(self.plot, 'make_{}'.format(base_plot_type.split('-')[0]))
         
        # set boolean on whether to plot obs in legend or not, and relevant data labels (data labels plotted)
        if (base_plot_type == 'scatter') or ('bias' in plot_options) or (z_statistic_sign == 'bias'):
            set_obs_legend = False
            relevant_data_labels = list(self.experiments.values())
        else:
            relevant_data_labels = copy.deepcopy(data_labels)

        # if multispecies is active then use all networkspecies, otherwise take first
        if 'multispecies' in plot_options:
            networkspecies = copy.deepcopy(self.networkspecies)
        # take first defined networkspeci
        else:
            networkspecies = [networkspeci]
            if len(self.networkspecies) > 1:
                print("Warning: More than 1 network or species defined, can only plot for 1 pair. Taking {}.".format(networkspeci))
            
        # legend plot (on its own axis)
        if base_plot_type == 'legend':
            legend_handles = self.make_legend(plot_type, data_labels=data_labels, set_obs=set_obs_legend)
            relevant_ax.legend(**legend_handles)

        # map plot
        elif base_plot_type == 'map':
            # get map data labels to plot
            # if no specific labels defined then take first data label and give warning
            if (labela == '') & (labelb == ''):
                labela = data_labels[0]
                print("Warning: No specific data labels set, plotting first available data label: {}.".format(labela))
            # labelb defined but labela for some reason, set labela to be labelb, and labelb empty str 
            elif (labela == ''):
                labela = copy.deepcopy(labelb)
                labelb = ''
            # set map title
            if z_statistic_sign == 'absolute':
                map_title = '{}'.format(labela)
            elif z_statistic_sign == 'bias':
                map_title = '{}'.format(labelb)

            func(relevant_ax, networkspeci, self.plot_characteristics[plot_type], plot_options, zstat=zstat, 
                 labela=labela, labelb=labelb)
        # periodic plot
        elif base_plot_type == 'periodic':
            func(grid_dict, networkspeci, data_labels, self.plot_characteristics[plot_type], plot_options, zstat=zstat)
        # make statsummary plot
        elif base_plot_type == 'statsummary':
            
            # get stats to plot
            if 'bias' in plot_options:
                stats_to_plot = self.plot_characteristics[plot_type]['experiment_bias']
            else:
                stats_to_plot = self.plot_characteristics[plot_type]['basic']

            # create empty dataframe with networkspecies and subsections
            index = pd.MultiIndex.from_product([self.networkspecies, self.subsections, relevant_data_labels],
                                                names=["networkspecies", "subsections", "labels"])
            stats_df = pd.DataFrame(np.nan, index=index, columns=stats_to_plot, dtype=np.float64)
            
            # fill dataframe
            is_initial = copy.deepcopy(self.initialised)
            kwargs = copy.deepcopy(self.kwargs)
            # save current subsection 
            orig_ss = copy.deepcopy(self.subsection)
            for ss in self.subsections:
                kwargs['subsection'] = ss
                self.set_config(**kwargs)
                # filter data
                self.reset_filter(initialise=True)
                for ns in networkspecies:
                    for dl in relevant_data_labels:
                        stats_per_data_label = []
                        for stp in stats_to_plot:
                            # get zstat information 
                            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=stp)
                            # calculate statistic
                            if dl in self.selected_station_data_labels[ns]:
                                # if relevant stat is expbias stat, then ensure temporal colocation is active
                                if (base_plot_type == 'statsummary') and (stp in self.expbias_stats) and ((not self.temporal_colocation) or (len(self.data_labels) == 1)):
                                    stats_per_data_label.append(np.NaN)
                                # otherwise calculate statistic
                                else:
                                    if z_statistic_sign == 'bias':
                                        stats_per_data_label.append(calculate_statistic(self, self, ns, zstat, [self.observations_data_label], [dl]))
                                    else:
                                        stats_per_data_label.append(calculate_statistic(self, self, ns, zstat, [dl], []))
                            else:
                                stats_per_data_label.append(np.NaN)

                        # get floats instead of arrays with 1 element each and save
                        stats_per_data_label = [stat_per_data_label[0] 
                                                if isinstance(stat_per_data_label, np.ndarray) 
                                                else stat_per_data_label 
                                                for stat_per_data_label in stats_per_data_label]
                        
                        # put data in dataframe
                        stats_df.loc[(ns, ss, dl)] = stats_per_data_label
                    
                # remove subsection variables from memory (if have subsections)
                # do not remove fixed section variables
                for k in self.subsection_opts:
                    if k not in self.fixed_section_vars:
                        try:
                            vars(self).pop(k)
                        except:
                            pass

            # make plot
            func(relevant_ax, networkspeci, relevant_data_labels, self.plot_characteristics[plot_type], plot_options,
                 statsummary=True, plotting_paradigm='summary', stats_df=stats_df)     

            # re-filter for original subsection
            kwargs['subsection'] = orig_ss
            self.set_config(**kwargs)
            if is_initial:
                self.reset_filter(initialise=True)
            else:
                self.reset_filter()

        # make heatmap / table plot
        elif base_plot_type in ['heatmap','table']:  

            # create empty dataframe with networkspecies and subsections
            index = pd.MultiIndex.from_product([networkspecies, self.subsections],
                                                names=["networkspecies", "subsections"])
            stats_df = pd.DataFrame(np.nan, index=index, columns=relevant_data_labels, dtype=np.float64)
            
            # fill dataframe
            is_initial = copy.deepcopy(self.initialised)
            kwargs = copy.deepcopy(self.kwargs)
            # save current subsection 
            orig_ss = copy.deepcopy(self.subsection)
            for ss in self.subsections:
                kwargs['subsection'] = ss
                self.set_config(**kwargs)
                # filter data
                self.reset_filter(initialise=True)
                for ns in networkspecies:
                    stat_per_data_labels = []
                    for dl in relevant_data_labels:
                        # calculate statistic
                        if dl in self.selected_station_data_labels[ns]:
                            if z_statistic_sign == 'bias':
                                stat_per_data_labels.append(calculate_statistic(self, self, ns, zstat, [self.observations_data_label], [dl]))
                            else:
                                stat_per_data_labels.append(calculate_statistic(self, self, ns, zstat, [dl], []))
                        else:
                            stat_per_data_labels.append(np.NaN)

                    # get floats instead of arrays with 1 element each and save
                    stat_per_data_labels = [stat_per_data_label[0] 
                                            if isinstance(stat_per_data_label, np.ndarray) 
                                            else stat_per_data_label 
                                            for stat_per_data_label in stat_per_data_labels]

                    # put data in dataframe
                    stats_df.loc[(ns, ss)] = stat_per_data_labels

                # remove subsection variables from memory (if have subsections)
                # do not remove fixed section variables
                for k in self.subsection_opts:
                    if k not in self.fixed_section_vars:
                        try:
                            vars(self).pop(k)
                        except:
                            pass

            # make plot
            func(relevant_ax, networkspeci, relevant_data_labels, 
                 self.plot_characteristics[plot_type], plot_options, plotting_paradigm='summary', 
                 stats_df=stats_df)

            # re-filter for original subsection
            kwargs['subsection'] = orig_ss
            self.set_config(**kwargs)
            if is_initial:
                self.reset_filter(initialise=True)
            else:
                self.reset_filter()

        # make timeseries plot
        elif base_plot_type == 'timeseries':
            func(relevant_ax, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                 plot_options, chunk_stat=chunk_stat, chunk_resolution=chunk_resolution)
        
        # make taylor diagram plot
        elif base_plot_type == 'taylor':
            stddev_max = self.selected_station_stddev_max[networkspeci]
            func(relevant_ax, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                 plot_options, zstat=zstat, stddev_max=stddev_max)
            
        # other plots
        elif base_plot_type != 'legend': 
            func(relevant_ax, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                 plot_options)

        # get number of total available stations, and individual station information if just have 1 station
        if self.temporal_colocation:
            station_inds = self.valid_station_inds_temporal_colocation[networkspeci][self.observations_data_label]
        else:
            station_inds = self.valid_station_inds[networkspeci][self.observations_data_label]  
        n_stations = len(station_inds)
        if n_stations == 1:
            station_ind = station_inds[0]
            current_lon = round(self.station_longitudes[networkspeci][station_ind], 2)
            current_lat = round(self.station_latitudes[networkspeci][station_ind], 2)
            current_station_name = self.station_names[networkspeci][station_ind]
            current_station_reference = self.station_references[networkspeci][station_ind]
        elif n_stations == 0:
            print('No valid stations for {} in {} subsection. Not making {} plot'.format(networkspeci, self.subsection, plot_type))
            return
            
        # set title, xlabel and ylabel for plots

        # set xlabel / ylabel
        if base_plot_type == 'periodic' or ((base_plot_type == 'timeseries') 
                                            and (chunk_stat is not None) 
                                            and (chunk_resolution is not None)):
            
            if not ylabel:         
                if z_statistic_type == 'basic':
                    ylabel = self.basic_stats[base_zstat]['label']
                    ylabel_units = self.basic_stats[base_zstat]['units']
                else:
                    ylabel = self.expbias_stats[base_zstat]['label']
                    ylabel_units = self.expbias_stats[base_zstat]['units']
                if ylabel_units == '[measurement_units]':
                    ylabel_units = self.measurement_units[speci] 
                if ylabel_units != '':
                    ylabel += ' [{}]'.format(ylabel_units)

        elif base_plot_type not in ['legend', 'metadata', 'map', 'heatmap', 'table', 'statsummary', 'taylor']:
            if not xlabel:
                if 'xlabel' in self.plot_characteristics[plot_type]:
                    xlabel = self.plot_characteristics[plot_type]['xlabel']['xlabel']
                    if '[measurement_units]' in xlabel:
                        xlabel = xlabel.replace('[measurement_units]', '[{}]'.format(self.measurement_units[speci]))

            if not ylabel:
                if 'ylabel' in self.plot_characteristics[plot_type]:
                    ylabel = self.plot_characteristics[plot_type]['ylabel']['ylabel']
                    if '[measurement_units]' in ylabel:
                        ylabel = ylabel.replace('[measurement_units]', '[{}]'.format(self.measurement_units[speci]))

        # set title
        if not title:
            if zstat is not None:
                if 'axis_title' in self.plot_characteristics[plot_type]:
                    title = self.plot_characteristics[plot_type]['axis_title']['label']
                    if title == '':
                        stat_label = generate_colourbar_detail(self, zstat, 0, 1, self.plot_characteristics[plot_type], 
                                                               speci, only_label=True)
                        if '[' in stat_label:
                            stat_label = stat_label.split('[')[0].strip()
                        if n_stations == 1:
                            title = '{} for {}, {} ({:.{}f}, {:.{}f})'.format(stat_label, current_station_reference,
                                                                              current_station_name, 
                                                                              current_lon,
                                                                              self.plot_characteristics[plot_type]['round_decimal_places']['title'],
                                                                              current_lat,
                                                                              self.plot_characteristics[plot_type]['round_decimal_places']['title'])
                            if base_plot_type == 'map':
                                title = '{} {}'.format(map_title, title)

                        else:
                            if zstat == 'NStations':
                                title = '{}'.format(stat_label)
                            else:
                                title = '{} at {} stations'.format(stat_label, n_stations)
                            if base_plot_type == 'map':
                                title = '{} {}'.format(map_title, title)

            elif base_plot_type not in ['legend', 'metadata']:
                if 'axis_title' in self.plot_characteristics[plot_type]:
                    title = self.plot_characteristics[plot_type]['axis_title']['label']
                    if title == '':
                        if n_stations == 1:
                            title = '{}, {} ({:.{}f}, {:.{}f})'.format(current_station_reference,
                                                                       current_station_name, 
                                                                       current_lon,
                                                                       self.plot_characteristics[plot_type]['round_decimal_places']['title'],
                                                                       current_lat,
                                                                       self.plot_characteristics[plot_type]['round_decimal_places']['title'])
                        else:
                            title = '{} stations'.format(n_stations)

                    if base_plot_type == 'fairmode-target':
                        speci = networkspeci.split('|')[1]
                        title += '\n{}'.format(fairmode_settings[speci]['title'])

        # overwrite passed xlabels and ylabels
        if title:
            set_axis_title(self, relevant_ax, title, self.plot_characteristics[plot_type])
        if xlabel:
            set_axis_label(relevant_ax, 'x', xlabel, self.plot_characteristics[plot_type])
        if ylabel:
            set_axis_label(relevant_ax, 'y', ylabel, self.plot_characteristics[plot_type])

        # format plot axis/axes
        format_axis(self, self, relevant_ax, base_plot_type, self.plot_characteristics[plot_type], 
                    map_extent=map_extent)

        # format plot options
        format_plot_options(self, self, relevant_ax, [relevant_data_labels], networkspeci, 
                            base_plot_type, plot_type, plot_options, map_extent=map_extent, 
                            chunk_stat=chunk_stat, chunk_resolution=chunk_resolution)                         

        # handle harmonisation of axes
        if base_plot_type == 'scatter':
            harmonise_xy_lims_paradigm(self, self, relevant_ax, base_plot_type, 
                                        self.plot_characteristics[plot_type], plot_options, relim=True)
        elif base_plot_type not in ['legend', 'metadata', 'map', 'taylor', 'fairmode-target']:
            harmonise_xy_lims_paradigm(self, self, relevant_ax, base_plot_type, 
                                        self.plot_characteristics[plot_type], plot_options, relim=True, 
                                        autoscale=True)

        # make legend (embedded on plot axis)
        if (legend) & (base_plot_type != 'legend'):
            if 'legend' in self.plot_characteristics[plot_type]:
            
                # only make map legend in 'domain' plot option is a active 
                # also remove observations from legend
                valid_legend = True
                if base_plot_type == 'map':
                    if 'domain' in plot_options:
                        set_obs_legend = False
                    else:
                        valid_legend = False

                if valid_legend:
                    legend_handles = self.make_legend(plot_type, data_labels=data_labels, set_obs=set_obs_legend)
                    if base_plot_type in ['periodic', 'periodic-violin']:
                        try:
                            ax_to_plot = self.plot_characteristics[plot_type]['legend']['handles']['ax']
                        except:
                            print("Warning: axis to plot legend on not defined for plot type in plot_characteristics.yaml, or passed via 'format' argument.\nTaking first available axis.")
                            ax_to_plot = self.relevant_temporal_resolutions[0]
                        if ax_to_plot not in self.relevant_temporal_resolutions:
                            print("Warning: defined axis to plot legend on not available for data resolution of read data.\nInstead, taking first available axis.")
                            ax_to_plot = self.relevant_temporal_resolutions[0]
                        relevant_ax[ax_to_plot].legend(**legend_handles)
                    else:
                        if base_plot_type == 'fairmode-target':
                            print("Warning: Data labels legend cannot be plotted, create single legend using make_plot function.")
                        else:
                            relevant_ax.legend(**legend_handles)

        # make colourbar (embedded on plot axis)
        if 'cb' in self.plot_characteristics[plot_type]:
            cb_ax = self.make_colourbar(fig, relevant_ax, zstat, speci, plot_type)
            # hide colourbar if requested, we still need to create it to get the correct cmap / bounds in the maps
            if not cb:
                cb_ax.set_visible(False)

        # if save is passed then save plot and return
        if save:
            # if save is boolean then auto generate fname
            if type(save) == bool:
                figure_fname = join(PROVIDENTIA_ROOT, 'plots/{}.png'.format(plot_type))
            else:
                figure_fname = copy.deepcopy(save)
            print('Saving {} figure to {}'.format(plot_type, figure_fname))
            # save figure
            plt.savefig(figure_fname, bbox_inches='tight')
            return None
        # elif return_plot is passed then return plot axis/axes
        elif return_plot:
            return fig
        # otherwise show plot
        else:
            plt.show()

    def make_colourbar(self, fig, plot_ax, stat, speci, plot_type):
        """ Wrapper method to make colourbar.

        :param fig: Figure
        :type fig: matplotlib.figure
        :param plot_ax: Axis
        :type plot_ax: matplotlib.axes
        :param stat: Statistic
        :type stat: str
        :param speci: Species
        :type speci: str
        :param plot_type: Plot type
        :type plot_type: str
        """

        # create cb axis
        cb_ax = fig.add_axes(self.plot_characteristics[plot_type]['cb']['position'])
        cb_ax.set_rasterized(True)

        # generate colourbar
        generate_colourbar(self, [plot_ax], [cb_ax], stat, self.plot_characteristics[plot_type], speci)

        return cb_ax
    
    def make_legend(self, plot_type, data_labels=None, set_obs=True):
        """ Wrapper method to make legend.

        :param plot_type: Plot type
        :type plot_type: str
        :param data_labels: Data arrays to plot, defaults to None
        :type data_labels: list, optional
        :param set_obs: Indicates if observations appear on legend, defaults to True
        :type set_obs: bool, optional
        :return: Legend
        :rtype: dict
        """
        
        if plot_type == 'legend':
            legend_characteristics = self.plot_characteristics['legend']
        elif 'legend' in self.plot_characteristics[plot_type]:
            legend_characteristics = self.plot_characteristics[plot_type]['legend']
        else:
            print("Warning: 'legend' not defined for plot type in plot_characteristics.yaml")
            return

        legend_handles = self.plot.make_legend_handles(legend_characteristics, data_labels=data_labels, set_obs=set_obs)
        
        return legend_handles['plot']

    def calculate_stat(self, stat, labela='', labelb='', per_station=False):
        """ Wrapper method to calculate statistic/s.

        :param stat: Statistic
        :type stat: str
        :param labela: Label of first dataset, defaults to ''
        :type labela: str, optional
        :param labelb: Label of second dataset (if defined then a bias plot is made), defaults to ''
        :type labelb: str, optional
        :param per_station: Indicates if the station data is per station or for all stations, defaults to False
        :type per_station: bool, optional
        :return: Statistic value
        :rtype: np.ndarray
        """

        # if no specific labels defined then take first data label and give warning
        if (labela == '') & (labelb == ''):
            labela = self.data_labels[0]
            print("Warning: No specific data labels set, plotting first available data label: {}.".format(labela))
        # labelb defined but labela for some reason, set labela to be labelb, and labelb empty str 
        elif (labela == ''):
            labela = copy.deepcopy(labelb)
            labelb = ''

        # get zstat information 
        zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=stat) 

        # if only 1 label passed and stat is a bias statistic then throw error
        if (z_statistic_sign == 'bias') & (labelb == ''):
              print("Warning: Calculating a bias statistic, and only 1 label is set. Cannot calculate statistic.")
              return

        # get networkspeci to calculate for
        networkspeci = self.networkspecies[0]
        if len(self.networkspecies) > 1:
            print("Warning: More than 1 network or species defined, can only calculate for 1. Taking {}.".format(networkspeci))

        if per_station:
            stat = calculate_statistic(self, self, networkspeci, stat, [labela], [labelb], per_station=True)
        else:
            stat = calculate_statistic(self, self, networkspeci, stat, [labela], [labelb])
        
        return stat

    def set_config(self, **kwargs):
        """ Wrapper method to set configuration variables. """

        # initialise default configuration variables
        # modified by passed arguments, if given
        provconf = ProvConfiguration(self, **kwargs)

        # for any passed arguments not in default Providentia variables, now set them to self
        for kwarg in kwargs:
            if kwarg not in provconf.var_defaults:
                setattr(self, kwarg, kwargs[kwarg])

        # update variables to set from config file
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
                error = 'Error: The path to the configuration file passed as an argument does not exist.'
                print(error)
                return False
        else:
            error = "Error: The configuration file must be given as an argument: e.g. 'config=...'"
            print(error)
            return False

        # parse section
        # if section name provided, try and use that
        # otherwise take first defined section name
        self.sections = copy.deepcopy(self.parent_section_names)

        #check if configuration file has a section title
        if len(self.sections) == 0:
            error = "Error: No sections were found in the configuration file, make sure to name them using square brackets."
            print(error)
            return False
    
        have_section = False
        if hasattr(self, 'section'): 
            # check that section actually exists
            if self.section in self.sections:
                have_section = True
            else:
                print("Warning: Defined section {} does not exist in configuration file.".format(self.section))
        if not have_section:
            self.section = self.sections[0]
            if len(self.sections) > 1:
                print("Warning: Taking first defined section ({}) to be read.".format(self.section))

        # update self with section variables (if not overwritten by kwargs)
        self.section_opts = self.sub_opts[self.section]
        for k, val in self.section_opts.items():
            if k not in kwargs:
                setattr(self, k, provconf.parse_parameter(k, val))

        # parse subsection
        # if subsection name is provided, try and use that
        # otherwise take first defined subsection name
        # if have no subsections, section is set as subsection name
        have_subsection = False
        # get subsection names
        self.subsections = [subsection_name for subsection_name in self.subsection_names 
                            if self.section == subsection_name.split('·')[0]]
        self.subsections_reduced = [subsection_name.split('·')[1] for subsection_name in self.subsections]

        if hasattr(self, 'subsection'): 
            # check that subsection actually exists
            if self.subsection in self.subsections:
                have_subsection = True
            elif self.subsection in self.subsections_reduced:
                have_subsection = True
                self.subsection = self.subsections[self.subsections_reduced.index(self.subsection)]
            elif self.subsection == self.section:
                have_subsection = True
            else:
                print("Warning: Defined subsection {} does not exist in configuration file.".format(self.subsection))

        if len(self.subsections) > 0:
            if not have_subsection:
                self.subsection = self.subsections[0]
                have_subsection = True
                if len(self.subsections) > 1:
                    print("Warning: Taking first defined subsection ({}) to be read.".format(self.subsection))
        else:
            self.subsections = [self.section]
            self.subsection = self.subsections[0]

        # update self with subsection variables (if have subsection) 
        if have_subsection:   
            # get subsection variables
            self.subsection_opts = self.sub_opts[self.subsection]
            # ensure all fixed section variables defined in subsection have same value as current section variables
            self.subsection_opts = {k: (self.section_opts[k] if k in self.fixed_section_vars else val) 
                                    for (k, val) in self.subsection_opts.items()}
            # update subsection variables
            for k, val in self.subsection_opts.items():
                if k not in kwargs:
                    setattr(self, k, provconf.parse_parameter(k, val))

        # now all variables have been parsed, check validity of those, throwing errors where necessary
        provconf.check_validity()

        return True

    def print_config(self, conf=None, config=None):
        """ Print selected config file to console.

        :param conf: Configuration file name, defaults to None
        :type conf: str, optional
        :param config: Configuration file name, defaults to None
        :type config: str, optional
        """

        # if conf or config not None, then print that file
        if conf:
            pass
        elif config:
            conf = copy.deepcopy(config) 
        # otherwise take it to be file previously loaded
        else:
            conf = copy.deepcopy(self.config)

        # check file exists
        if not os.path.isfile(conf):
            print("Warning: The passed .conf file: '{}' does not exist.".format(conf))
        # otherwise, print conf
        else:
            with open(conf, "r") as f:
                print(f.read())

    def select_station(self, station):
        """ Wrapper method to select specific station/s.

        :param station: Station reference
        :type station: str
        """
        
        if type(station) == 'str':
            stations_to_keep = [station]
        else:
            stations_to_keep = station
        self.metadata_menu['STATION MISCELLANEOUS']['station_reference']['checkboxes']['keep_selected'] = stations_to_keep

        # filter for station/s    
        self.filter()

    def apply_filter(self, field, limit=None, keep=None, remove=None, lower=None, upper=None):
        """ Wrapper method to select specific station/s.

        :param field: field to filter by
        :type field: str
        :param limit: limit for filtering representativity fields
        :type limit: str, optional
        :param keep: data to keep
        :type keep: str, optional
        :param remove: data to remove
        :type remove: str, optional
        :param lower: lower bound to retain data
        :type lower: str, optional
        :param upper: upper bound to retain data
        :type upper: str, optional
        """

        # variable to know if to filter or not
        do_filter = False

        # make sure keep and remove arguments are lists
        if keep is not None:
            if type(keep) == str:
                keep = [keep]
        if remove is not None:
            if type(remove) == str:
                remove = [remove]

        # field is a represenativity field?
        if field in self.representativity_menu['rangeboxes']['map_vars']:
            do_filter = True

            field_index = self.representativity_menu['rangeboxes']['map_vars'].index(field)
            # ensure limit is set for field
            if limit is not None:
                self.representativity_menu['rangeboxes']['current_lower'][field_index] = limit
            else:
                print("Warning: When filtering by representativity field: {}, 'limit' must be passed as an argument.".format(field)) 
                return

        # field is a period field?
        elif field == 'period':
            do_filter = True

            # if neither keep or remove are defined, filtering cannot be done
            if (keep is None) and (remove is None): 
                print("Warning: When filtering by a period field, 'keep' or 'remove' must be passed as arguments.")
                return

            if keep is not None:
                new_keep = copy.deepcopy(self.period_menu['checkboxes']['keep_selected'])
                for item in keep:
                    if item not in new_keep:
                        new_keep.append(item)
                self.period_menu['checkboxes']['keep_selected'] = new_keep
            if remove is not None:
                new_remove = copy.deepcopy(self.period_menu['checkboxes']['remove_selected'])
                for item in remove:
                    if item not in new_remove:
                        new_remove.append(item)
                self.period_menu['checkboxes']['remove_selected'] = new_remove

        # fields is a metadata field?
        else:
            for menu_type in self.metadata_types:
                if field in self.metadata_menu[menu_type]['rangeboxes']['labels']:
                    do_filter = True

                    # if neither lower or upper are defined, filtering cannot be done
                    if (lower is None) and (upper is None): 
                        print("Warning: When filtering by a numeric metadata field, 'lower' or 'upper' must be passed as arguments.")
                        return

                    field_index = self.metadata_menu[menu_type]['rangeboxes']['labels'].index(field)
                    if lower is not None:
                        self.metadata_menu[menu_type]['rangeboxes']['current_lower'][field_index] = lower
                    if upper is not None:
                        self.metadata_menu[menu_type]['rangeboxes']['current_upper'][field_index] = upper           
                    if field not in self.metadata_menu[menu_type]['rangeboxes']['apply_selected']:
                        self.metadata_menu[menu_type]['rangeboxes']['apply_selected'].append(field)
                    break

                elif field in self.metadata_menu[menu_type]['navigation_buttons']['labels']:
                    do_filter = True

                    # if neither keep or remove are defined, filtering cannot be done
                    if (keep is None) and (remove is None): 
                        print("Warning: When filtering by a text period field, 'keep' or 'remove' must be passed as arguments.")
                        return

                    if keep is not None:
                        new_keep = copy.deepcopy(self.metadata_menu[menu_type][field]['checkboxes']['keep_selected'])
                        for item in keep:
                            if item not in new_keep:
                                new_keep.append(item)
                        self.metadata_menu[menu_type][field]['checkboxes']['keep_selected'] = new_keep
                    if remove is not None:
                        new_remove = copy.deepcopy(self.metadata_menu[menu_type][field]['checkboxes']['remove_selected'])
                        for item in remove:
                            if item not in new_remove:
                                new_remove.append(item)
                        self.metadata_menu[menu_type][field]['checkboxes']['remove_selected'] = new_remove
                    break

        # do filtering?
        if do_filter: 
            self.filter()
        # otherwise set warning that field was not found
        else:
            print('Warning: {} not available for filtering.'.format(field)) 
        
    def save(self, fname='', format='nc'):
        """ Wrapper method to save current data/ metadata in memory.

        :param fname: File name, defaults to ''
        :type fname: str, optional
        :param format: File format, defaults to 'nc'
        :type format: str, optional
        """

        # set fname if not provided
        if fname == '':
            date_str = datetime.datetime.today().strftime('%Y%m%d_%H%M')
            fname = join(PROVIDENTIA_ROOT, 'saved_data/PRV_{}'.format(date_str))

        if format in ['conf','config','.conf']:
            fname = '{}.conf'.format(fname)
            export_configuration(self, fname)

        elif format in ['netCDF', 'netcdf', 'netCDF4', 'netcdf4', 'nc', '.nc']:
            fname = '{}.nc'.format(fname)
            export_netcdf(self, fname)

        elif format in ['npz','.npz','np','.np','npy','.npy','numpy']:
            fname = '{}.npz'.format(fname)
            export_data_npz(self, fname)

        print('Data saved to {}'.format(fname))

    def get_data(self, format='nc'):
        """ Wrapper method return data / metadata in specific format.

        :param format: File format, defaults to 'nc'
        :type format: str, optional
        :return: Data
        :rtype: numpy.ndarray
        """

        # for non-ghost, update networkspecies name 
        # (e.g. ghost_btx/ghost_btx|sconcc6h6 -> ghost_btx-ghost_btx|sconcc6h6)
        # this will avoid permission denied errors
        networkspeci = self.networkspecies[0]
        if '/' in networkspeci:
            networkspeci = networkspeci.replace('/', '-')

        # set temporary fname for writing
        temporary_fname = join(PROVIDENTIA_ROOT, 'saved_data/temp_{}'.format(networkspeci))
        
        # check if temporary fname already exists 
        if os.path.isfile(temporary_fname):
            # if so, keep iterating until find fname is new
            invalid_fname = True
            dup_count = 2
            while invalid_fname:
                temporary_fname = join(PROVIDENTIA_ROOT, 'saved_data/temp_{}_{}'.format(networkspeci, dup_count))
                if os.path.isfile(temporary_fname):
                    dup_count += 1
                else:
                    invalid_fname = False
        
        if format in ['netCDF', 'netcdf', 'netCDF4', 'netcdf4', 'nc', '.nc']:
            data = export_netcdf(self, temporary_fname, set_in_memory=True)

        elif format in ['xr', '.xr', 'xarr', 'xarray','Xarray']:
            data = export_netcdf(self, temporary_fname, set_in_memory=True, xarray=True)

        elif format in ['npz','.npz','np','.np','npy','.npy','numpy']:
            data = export_data_npz(self, temporary_fname, set_in_memory=True)

        return data

    def get_var(self, var=''):
        """ Wrapper method to return specific data / metadata variable.

        :param var: Variable name, defaults to ''
        :type var: str, optional
        :return: Data
        :rtype: numpy.ndarray
        """

        # if variable is undefined then print warning
        if var == '':
            print("Warning: Variable to read is undefined.")
            return 
        else:
            data = self.get_data(format='nc')
            if var not in data.variables.keys():
                # try adding networkspeci to variable, if just have 1 networkspecies
                if len(self.networkspecies) == 1:
                    test_var = '{}_{}'.format(self.networkspecies[0], var)
                    if test_var in data.variables.keys():
                        var_data = data[test_var][:]
                        return var_data
                print("Warning: Variable '{}' is not defined".format(var))
            else:
                var_data = data[var][:]
                return var_data
            