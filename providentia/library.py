""" Class for Providentia library mode """

import copy
import datetime
import os
import random
import sys

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as fa
import numpy as np
import pandas as pd
import yaml

from providentia.auxiliar import CURRENT_PATH, join, expand_plot_characteristics, Tee
from .configuration import load_conf
from .configuration import ProvConfiguration
from .fields_menus import (init_metadata, init_period, init_representativity, metadata_conf,
                           update_metadata_fields, update_period_fields, update_representativity_fields,
                           period_conf, representativity_conf)
from .filter import DataFilter
from .plotting import Plotting
from .plot_aux import get_taylor_diagram_ghelper
from .plot_formatting import (format_plot_options, format_axis, set_axis_label, set_axis_title, 
                              harmonise_xy_lims_paradigm)
from .read import DataReader
from .read_aux import (generate_file_trees, get_possible_resampling_resolutions, 
                       get_periodic_nonrelevant_temporal_resolutions, get_periodic_relevant_temporal_resolutions, 
                       get_valid_models, get_valid_obs_files_in_date_range)
from .statistics import (calculate_statistic, generate_colourbar, generate_colourbar_detail, 
                         get_fairmode_data, get_selected_station_data, get_z_statistic_info)
from .warnings_prv import show_message
from .writing import export_configuration, export_data_npz, export_netcdf


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

class Providentia:
    """Class for Providentia Library mode"""

    def __init__(self, config, **kwargs):
        """
        Initialize the library mode environment

        Parameters
        ----------
        config : str
            Path to the configuration file.
        **kwargs : dict
            Optional command-line arguments that override default configuration values.
        """
        
        # set config to self
        self.config = config

        # set kwargs to self
        self.kwargs = kwargs

        # make configuration file visible to other modes
        self.kwargs['config'] = config

        # update kwargs to detect library mode
        self.kwargs['library'] = True

        # load statistical yamls
        self.basic_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/basic_stats.yaml')))
        self.modbias_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/model_bias_stats.yaml')))

        # load representativity information
        self.representativity_info = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/internal/representativity.yaml')))

        # set configuration variables, as well as any other defined variables
        self.valid_config = self.set_config(**self.kwargs)

        # initialise DataReader class
        self.datareader = DataReader(self)

        # check for self defined plot characteristics file
        if self.tests:
            mode = 'tests'
        else:
            mode = 'library'
        if self.plot_characteristics_filename == '':
            self.plot_characteristics_filename = join(PROVIDENTIA_ROOT, 'settings/plot_characteristics.yaml')
        plot_characteristics = yaml.safe_load(open(self.plot_characteristics_filename))
        self.plot_characteristics_templates = expand_plot_characteristics(plot_characteristics, mode)

        # initialise Plotting class
        self.plotting = Plotting(read_instance=self, canvas_instance=self)

        # add general plot characteristics to self
        for k, val in self.plot_characteristics_templates['general'].items():
            if k not in kwargs:
                setattr(self, k, val)

        # set some key configuration variables
        self.periodic_relevant_temporal_resolutions = get_periodic_relevant_temporal_resolutions(self.resolution)
        self.periodic_nonrelevant_temporal_resolutions = get_periodic_nonrelevant_temporal_resolutions(self.resolution)
        self.data_labels = [self.observations_data_label] + list(self.experiments.values())
        self.data_labels_raw = [self.observations_data_label] + list(self.experiments.keys())
        self.networkspecies = ['{}|{}'.format(network,speci) for network, speci in zip(self.network, self.species)]

    def read(self):
        """Wrapper method to read data."""

        self.logger.info('Reading data')

        # read data
        self.datareader.read_setup(['reset'])

    def apply_filter(self):
        """Method to apply filters to data."""

        self.logger.info('Filtering data')

        # filter data
        DataFilter(self)

        # get selected station data
        get_selected_station_data(read_instance=self, canvas_instance=self, networkspecies=self.networkspecies)

    def filter(self, field, limit=None, keep=None, remove=None, lower=None, upper=None):
        """
        Wrapper method to filter data.

        Parameters
        ----------
        field : str
            Field to filter by.
        limit : str, optional
            Limit for filtering representativity fields.
        keep : str or list of str, optional
            Data to keep.
        remove : str or list of str, optional
            Data to remove.
        lower : str, optional
            Lower bound to retain data.
        upper : str, optional
            Upper bound to retain data.
        """

        # check have valid conf and have loaded data
        valid_config = self.have_valid_config()
        if not valid_config:
            return
        loaded_data = self.loaded_data()
        if not loaded_data:
            return

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
                msg = "When filtering by representativity field: {}, 'limit' must be passed as an argument.".format(field)
                show_message(self, msg)
                return

        # field is a period field?
        elif field == 'period':
            do_filter = True

            # if neither keep or remove are defined, filtering cannot be done
            if (keep is None) and (remove is None): 
                msg = "When filtering by a period field, 'keep' or 'remove' must be passed as arguments."
                show_message(self, msg)
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
                        msg = "When filtering by a numeric metadata field, 'lower' or 'upper' must be passed as arguments."
                        show_message(self, msg)
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
                        msg = "When filtering by a text period field, 'keep' or 'remove' must be passed as arguments."
                        show_message(self, msg)
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
            self.apply_filter()
        # otherwise set warning that field was not found
        else:
            msg = '{} not available for filtering.'.format(field)
            show_message(self, msg)
        
    def filter_station(self, station):
        """ 
        Wrapper method to filter specific station/s.

        Parameters
        ----------
        station : str
            Station reference
        """
        
        # check have valid conf and have loaded data
        valid_config = self.have_valid_config()
        if not valid_config:
            return
        loaded_data = self.loaded_data()
        if not loaded_data:
            return

        if type(station) == 'str':
            stations_to_keep = [station]
        else:
            stations_to_keep = station
        self.metadata_menu['STATION MISCELLANEOUS']['station_reference']['checkboxes']['keep_selected'] = stations_to_keep

        # filter for station/s    
        self.apply_filter()

    def reset(self, initialise=False):
        """ 
        Wrapper method to reset filter data.

        Parameters
        ----------
        initialise : bool, optional
            Indicates whether to reset data to initial state when class was initialised 
        """

        # check have valid conf and have loaded data
        valid_config = self.have_valid_config()
        if not valid_config:
            return
        loaded_data = self.loaded_data()
        if not loaded_data:
            return

        if initialise:
            self.logger.info(f'Resetting data filters to when class was initialised, loading {self.subsection} subsection filters.')
        else:
            self.logger.info(f'Resetting all data filters.')

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
        self.apply_filter()

        # for non-GHOST, we call update_metadata_fields after filtering to remove the stations that have
        # 0 valid measurements, to do this we need to have valid_station_inds, which is obtained 
        # after filtering
        if not self.reading_ghost:
            update_metadata_fields(self)

        # set variable to know if data is in intial state or not
        if initialise:
            self.initialised = True
        else:
            self.initialised = False

    def plot(self, plot, data_labels=None, labela='', labelb='', title=None, xlabel=None, ylabel=None, cb=True, 
             legend=True, set_obs_legend=True, map_extent=None, annotate=False, bias=False, domain=False, 
             hidedata=False, logx=False, logy=False, multispecies=False, regression=False, smooth=False, 
             threshold=False, gerrity=False, plot_options=None, save=False, return_plot=False, format=None, 
             width=None, height=None):
        """ 
        Wrapper method to make a Providentia plot.

        Parameters
        ----------
        plot : str
            Plot type.
        data_labels : list of str, optional
            Data arrays to plot, defaults to None.
        labela : str, optional
            Label of first dataset, defaults to ''.
        labelb : str, optional
            Label of second dataset (if defined then a bias plot is made), defaults to ''.
        title : str, optional
            Axes title, defaults to None.
        xlabel : str, optional
            Label on x axes, defaults to None.
        ylabel : str, optional
            Label on y axes, defaults to None.
        cb : bool, optional
            Indicates if colorbar appears on plot, defaults to True.
        legend : bool, optional
            Indicates if legend appears on plot, defaults to True.
        set_obs_legend : bool, optional
            Indicates if observations appear on legend, defaults to True.
        map_extent : list, optional
            Map extent, defaults to None.
        annotate : bool or list, optional
            Indicates if there are annotations or a list of statistics to annotate, defaults to False.
        bias : bool, optional
            Indicates if data is biased, defaults to False.
        domain : bool, optional
            Indicates if domain shows in maps, defaults to False.
        hidedata : bool, optional
            Indicates if data points are hidden in plot, defaults to False.
        logx : bool, optional
            Indicates if the scale of the x axis is log, defaults to False.
        logy : bool, optional
            Indicates if the scale of the y axis is log, defaults to False.
        multispecies : bool, optional
            Indicates if plot has multispecies, defaults to False.
        regression : bool, optional
            Indicates if scatter plot has regression line/s, defaults to False.
        smooth : bool, int, float, optional
            Indicates if timeseries has smooth line/s or the smoothing window, defaults to False.
        threshold : bool, optional
            Indicates if plot has threshold line/s, defaults to False.
        gerrity : bool, optional
            Indicates if plot shows Gerrity scores per station
        plot_options : list, optional
            List with plot options, defaults to None.
        save : bool or str, optional
            Indicates if you want to save the figure, defaults to False.
        return_plot : bool, optional
            Indicates if you want to get the figure object, defaults to False.
        format : dict, optional
            Format to overwrite the plot characteristics, defaults to None.
        width : int or float, optional
            Figure width, defaults to None.
        height : int or float, optional
            Figure height, defaults to None.

        Returns
        -------
        matplotlib.figure.Figure or None
            Returns the figure object if `return_plot` is True. Otherwise, displays the plot or saves it to file.
        """

        # check have valid conf and have loaded data
        valid_config = self.have_valid_config()
        if not valid_config:
            return
        loaded_data = self.loaded_data()
        if not loaded_data:
            return

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
        if gerrity:
            if 'gerrity' not in plot_options:
                plot_options.append('gerrity')          

        # get base plot type (no plot options), and plot type (with plot options)
        base_plot_type = copy.deepcopy(plot)
        if len(plot_options) > 0:
            plot_type = '{}_{}'.format(base_plot_type, '_'.join(plot_options))
        else:
            plot_type = copy.deepcopy(plot)

        # get zstat for required plots
        base_plot_type_split = base_plot_type.split('-')
        if (len(base_plot_type_split) > 1) & (base_plot_type not in ['periodic-violin', 'fairmode-target', 
                                                                     'fairmode-statsummary']):
            base_plot_type = base_plot_type_split[0]
            zstat = base_plot_type_split[1]
        else:
            zstat = None

        # get networkspeci to plot (for non-multispecies plots), taking first one preferentially
        if len(self.networkspecies) > 0:
            networkspeci = self.networkspecies[0]
        else:
            msg = 'There are no available species.'
            show_message(self, msg)
            return
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
                
                # get available chunk timeseries resolutions
                available_timeseries_chunk_resolutions = get_possible_resampling_resolutions(self.active_resolution,
                                                                                            daily_forecast=self.daily_forecast)

                # show warning if it is not
                if chunk_resolution not in available_timeseries_chunk_resolutions:
                    msg = f'{plot_type} cannot be created because {chunk_resolution} '
                    msg += 'is not an available chunking resolution. '
                    if len(available_timeseries_chunk_resolutions) > 0:
                        msg += f'The available resolutions are: {available_timeseries_chunk_resolutions}'
                    show_message(self, msg)
                    return

                # show warning if chunk stat is MDA8 and active resolution is not hourly
                if (chunk_stat == 'MDA8') and (self.active_resolution != 'hourly'):
                    msg = f'{plot_type} cannot be created because {chunk_stat} '
                    msg += 'is only available when active resolution is hourly.'
                    show_message(self, msg)
                    return
                
                # show warning if chunk stat is MDA8 and chunk resolution is not daily
                if (chunk_stat == 'MDA8') and (chunk_resolution != 'daily'):
                    msg = f'{plot_type} cannot be created because {chunk_stat} '
                    msg += 'is only available when chunk resolution is daily.'
                    show_message(self, msg)
                    return

        # get zstat information 
        zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type=plot_type) 

        # if only 1 label passed for map plot, and stat is a bias statistic then throw error
        if (base_plot_type == 'map') & (z_statistic_sign == 'bias') & (labelb == ''):
            msg = "Plotting a bias statistic, and only 1 label is set. Not making plot."
            show_message(self, msg)
            return
        
        # if bias and threshold plots are in plot options throw error
        if ('bias' in plot_options) & ('threshold' in plot_options):
            msg = "Cannot make a bias plot showing threshold lines. Not making plot."
            show_message(self, msg)
            return

        # do not make plot if hidedata is active but smooth is not in plot options
        if (base_plot_type == 'timeseries') and ('hidedata' in plot_options) and ('smooth' not in plot_options):
            msg = f"Cannot make {plot_type} because 'hidedata' plot option is set for "
            msg += "timeseries plot, but 'smooth' is not active. Not making plot."
            show_message(self, msg)
            return
        
        # do not make plot if hidedata is active but regression is not in plot options
        if (base_plot_type == 'scatter') and ('hidedata' in plot_options) and ('regression' not in plot_options):
            msg = f"Cannot make {plot_type} because 'hidedata' plot option is set for "
            msg += "scatter lot, but 'regression' is not active. Not making plot."
            show_message(self, msg)
            return
        
        # do not make Taylor diagram if statistic is not r or r2
        if (base_plot_type == 'taylor') and (zstat not in ['r', 'r2']):
            msg = f"Cannot make {plot_type} because statistic is not available or defined. "
            msg += "Choose between 'taylor-r' or 'taylor-r2'. Not making plot."
            show_message(self, msg)
            return
        
        # do not make periodic plot if stat is MDA8
        if (base_plot_type == 'periodic') and (base_zstat == 'MDA8'):
            msg = f"Cannot make {plot_type} because MDA8 statistic is not available for periodic plots. Not making plot."
            show_message(self, msg)
            return
        
        # do not make statsummary plot if stat is MDA8, and are making periodic statistic
        if (base_plot_type == 'statsummary') and (base_zstat == 'MDA8') and (z_statistic_period is not None):
            msg = f"Cannot make {plot_type} because MDA8 statistic is not available for periodic statistics. Not making plot."
            show_message(self, msg)
            return

        if base_plot_type in ['fairmode-target','fairmode-statsummary']:
            # warning for fairmode plots if species aren't PM2.5, PM10, NO2 or O3
            if speci not in ['sconco3', 'sconcno2', 'pm10', 'pm2p5']:
                msg = f'Fairmode plot cannot be created for {speci}.'
                show_message(self, msg)
                return
            # warning for fairmode plots if resolution is not hourly or daily
            if ((speci in ['sconco3', 'sconcno2'] and self.active_resolution != 'hourly') 
                or (speci in ['pm10', 'pm2p5'] and (self.active_resolution not in ['hourly', 'daily']))):
                msg = 'Fairmode plot can only be created if the resolution is hourly (O3, NO2, PM2.5 and PM10) or daily (PM2.5 and PM10).'
                show_message(self, msg)
                return
            
            # skip making plot if there is no valid data
            data, valid_station_idxs = get_fairmode_data(self, self, networkspeci, self.data_labels)
            if not any(valid_station_idxs):
                self.logger.info(f'No data after filtering by coverage for {speci}.')
                return
            
        if base_plot_type == 'contingencytable':
            # warning for contingency table if species aren't PM2.5, PM10, NO2, O3, or SO2
            if speci not in ['sconco3', 'sconcno2', 'pm10', 'pm2p5', 'sconcso2']:
                msg = f'Contingency table cannot be created for {speci}.'
                show_message(self, msg)
                return

        # get data labels for plot
        if len(data_labels) == 0:
            data_labels = copy.deepcopy(self.data_labels)
        # if any passed data labels are not available then pass warning
        else:
            invalid_data_labels = [data_label for data_label in data_labels if data_label not in self.data_labels]
            data_labels = [data_label for data_label in data_labels if data_label in self.data_labels]
            if len(data_labels) == 0:
                msg = "None of the passed data labels are available. Not making plot."
                show_message(self, msg)
                return
            elif len(invalid_data_labels) > 0:
                msg = "Passed data labels {} are not available.".format(invalid_data_labels)
                show_message(self, msg)

        # set plot characteristics
        self.plot_characteristics = dict()
        valid_plot = self.plotting.set_plot_characteristics([plot_type], format=format, data_labels=data_labels)

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
            msg = "Width and/or height have not been passed. The default values will be set."
            show_message(self, msg)
            fig = plt.figure(figsize=self.plot_characteristics[plot_type]['figsize'])

        # create axes
        main_gs = gridspec.GridSpec(2, 1, **self.plot_characteristics[plot_type]['main_gs'])
        if base_plot_type == 'map':
            ax = fig.add_subplot(main_gs[0], projection=self.plotcrs)
        elif base_plot_type == 'taylor':            
            reference_stddev = 7.5
            ghelper = get_taylor_diagram_ghelper(reference_stddev, self.plot_characteristics[plot_type])
            ax = fig.add_subplot(main_gs[0], axes_class=fa.FloatingAxes, grid_helper=ghelper)
        elif base_plot_type == "fairmode-statsummary":
            # get current species
            speci = networkspeci.split('|')[1]

            # get number of rows and columns
            ncols = 4
            nrows = 8 if speci in ["sconco3", "sconcno2", "pm10"] else 7

            # create gridspec and add it to a list
            gs =  gridspec.GridSpecFromSubplotSpec(nrows, ncols, subplot_spec=main_gs[0],
                                                **self.plot_characteristics["fairmode-statsummary"]["gridspec_kw"])
            ax = [fig.add_subplot(gs[i, j]) for i in range(nrows) for j in range(ncols)]
        else:
            ax = fig.add_subplot(main_gs[0])

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
            func = getattr(self.plotting, 'make_table')
        elif base_plot_type in ['fairmode-target', 'fairmode-statsummary']:
            func = getattr(self.plotting, 'make_{}'.format(base_plot_type.replace('-','_')))
        elif base_plot_type != 'legend':
            func = getattr(self.plotting, 'make_{}'.format(base_plot_type.split('-')[0]))
        
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
                msg = "More than 1 network or species defined, can only plot for 1 pair. Taking {}.".format(networkspeci)
                show_message(self, msg)
            
        # legend plot (on its own axis)
        if base_plot_type == 'legend':
            legend_handles = self.legend(plot_type, data_labels=data_labels, set_obs=set_obs_legend)
            relevant_ax.legend(**legend_handles)

        # map plot
        elif base_plot_type == 'map':
            # get map data labels to plot
            # if no specific labels defined then take first data label and give warning
            if (labela == '') & (labelb == ''):
                labela = data_labels[0]
                msg = "No specific data labels set, plotting first available data label: {}.".format(labela)
                show_message(self, msg)
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
                stats_to_plot = self.plot_characteristics[plot_type]['model_bias']
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
                self.reset(initialise=True)
                for ns in networkspecies:
                    for dl in relevant_data_labels:
                        stats_per_data_label = []
                        for stp in stats_to_plot:
                            # get zstat information 
                            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=stp)
                            # calculate statistic
                            if dl in self.selected_station_data_labels[ns]:
                                # if relevant stat is modbias stat, then ensure temporal colocation is active
                                if (base_plot_type == 'statsummary') and (stp in self.modbias_stats) and ((not self.temporal_colocation) or (len(self.data_labels) == 1)):
                                    stats_per_data_label.append(np.nan)
                                # otherwise calculate statistic
                                else:
                                    if z_statistic_sign == 'bias':
                                        stats_per_data_label.append(calculate_statistic(self, self, ns, zstat, [self.observations_data_label], [dl]))
                                    else:
                                        stats_per_data_label.append(calculate_statistic(self, self, ns, zstat, [dl], []))
                            else:
                                stats_per_data_label.append(np.nan)

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
                self.reset(initialise=True)
            else:
                self.reset()

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
                self.reset(initialise=True)
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
                            stat_per_data_labels.append(np.nan)

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
                self.reset(initialise=True)
            else:
                self.reset()

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

        # get relevant station inds
        if self.temporal_colocation:
            inds_array = self.valid_station_inds_temporal_colocation
        else:
            inds_array = self.valid_station_inds

        if labela != '':
            labela_station_inds = inds_array[networkspeci][labela]
            if labelb == '':
                station_inds = copy.deepcopy(labela_station_inds)
            else:
                labelb_station_inds = inds_array[networkspeci][labelb]
                station_inds = np.unique([labela_station_inds,labelb_station_inds])
        else:
            label_station_inds = []
            if data_labels:
                dls = copy.deepcopy(data_labels)
            else:
                dls = copy.deepcopy(self.data_labels)
            for dl in dls:
                if ('bias' in plot_options) & (dl == self.observations_data_label):
                    continue
                label_station_inds.extend(inds_array[networkspeci][dl])
            station_inds = np.unique(label_station_inds)

        # get number of total available stations, and individual station information if just have 1 station
        n_stations = len(station_inds)
        if n_stations == 1:
            station_ind = station_inds[0]
            current_lon = round(self.station_longitudes[networkspeci][station_ind], 2)
            current_lat = round(self.station_latitudes[networkspeci][station_ind], 2)
            current_station_name = self.station_names[networkspeci][station_ind]
            current_station_reference = self.station_references[networkspeci][station_ind]
        elif n_stations == 0:
            self.logger.info('No valid stations for {} in {} subsection. Not making {} plot'.format(networkspeci, self.subsection, plot_type))
            return
        
        # set xlabel / ylabel
        if base_plot_type == 'periodic' or ((base_plot_type == 'timeseries') 
                                            and (chunk_stat is not None) 
                                            and (chunk_resolution is not None)):
            
            if not ylabel:         
                if z_statistic_type == 'basic':
                    ylabel = self.basic_stats[base_zstat]['label']
                    ylabel_units = self.basic_stats[base_zstat]['units']
                else:
                    ylabel = self.modbias_stats[base_zstat]['label']
                    ylabel_units = self.modbias_stats[base_zstat]['units']
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

            if (zstat is not None) & (base_plot_type not in ['statsummary']):
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

                    if base_plot_type in ['fairmode-target','fairmode-statsummary']:
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
        elif base_plot_type not in ['legend', 'metadata', 'map', 'taylor', 'fairmode-statsummary']:
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
                    legend_handles = self.legend(plot_type, data_labels=data_labels, set_obs=set_obs_legend)
                    legend_ax = fig.add_subplot(main_gs[1])
                    legend_ax.axis("off")
                    legend_ax.legend(**legend_handles)

        # make colourbar (embedded on plot axis)
        if 'cb' in self.plot_characteristics[plot_type]:
            cb_ax = self.colourbar(fig, relevant_ax, zstat, speci, plot_type)
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
            self.logger.info('Saving {} figure to {}'.format(plot_type, figure_fname))
            # save figure
            plt.savefig(figure_fname, bbox_inches='tight')
            return None
        # elif return_plot is passed then return plot axis/axes
        elif return_plot:
            return fig
        # otherwise show plot
        else:
            plt.show()

    def colourbar(self, fig, plot_ax, stat, speci, plot_type):
        """
        Wrapper method to make a colourbar.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure object to which the colourbar will be added.
        plot_ax : matplotlib.axes.Axes
            Axis for the main plot.
        stat : str
            Statistic to display in the colourbar.
        speci : str
            Species for which the colourbar applies.
        plot_type : str
            Plot type, used to fetch plot characteristics.

        Returns
        -------
        matplotlib.axes.Axes
            Axis object of the created colourbar.
        """

        # create cb axis
        cb_ax = fig.add_axes(self.plot_characteristics[plot_type]['cb']['position'])
        cb_ax.set_rasterized(True)

        # generate colourbar
        generate_colourbar(self, [plot_ax], [cb_ax], stat, self.plot_characteristics[plot_type], speci)

        return cb_ax
    
    def legend(self, plot_type, data_labels=None, set_obs=True):
        """
        Wrapper method to make legend.

        Parameters
        ----------
        plot_type : str
            Plot type for which the legend will be generated.
        data_labels : list, optional
            Data arrays to plot, by default None.
        set_obs : bool, optional
            Indicates if observations appear on the legend, by default True.

        Returns
        -------
        dict
            Legend handles for the plot.
        """
        
        if plot_type == 'legend':
            legend_characteristics = self.plot_characteristics['legend']
        elif 'legend' in self.plot_characteristics[plot_type]:
            legend_characteristics = self.plot_characteristics[plot_type]['legend']
        else:
            msg = "'legend' not defined for plot type in plot_characteristics.yaml"
            show_message(self, msg)
            return

        legend_handles = self.plotting.make_legend_handles(legend_characteristics, data_labels=data_labels, set_obs=set_obs)
        
        return legend_handles['plot']

    def statistic(self, stat, labela='', labelb='', per_station=False, period=None, chunk=None):
        """
        Wrapper method to calculate statistic/s.

        Parameters
        ----------
        stat : str
            Statistic to calculate.
        labela : str, optional
            Label of first dataset, by default ''.
        labelb : str, optional
            Label of second dataset (if defined then a bias statistic is calculated), by default ''.
        per_station : bool, optional
            Indicates if the statistic should be calculated per station or integrated for all stations, by default False.
        period : str, optional
            Period to group data into for calculation of statistics. Current options: 'hour', 'dayofweek', 'month'.
        chunk : str, optional
            Chunked temporal window to calculate statistics for. Current options: 'daily', 'monthly', 'annual'.

        Returns
        -------
        np.ndarray
            Calculated statistic values.
        """

        # check have valid conf and have loaded data
        valid_config = self.have_valid_config()
        if not valid_config:
            return
        loaded_data = self.loaded_data()
        if not loaded_data:
            return

        # if no specific labels defined then take first data label and give warning
        if (labela == '') & (labelb == ''):
            labela = self.data_labels[0]
            msg = "No specific data labels set, plotting first available data label: {}.".format(labela)
            show_message(self, msg)
        # labelb defined but labela for some reason, set labela to be labelb, and labelb empty str 
        elif (labela == ''):
            labela = copy.deepcopy(labelb)
            labelb = ''

        # get zstat information 
        zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=stat) 

        # combine basic and modbias stats dicts together
        stats_dict = {**self.basic_stats, **self.modbias_stats}

        # check desired statistic is defined in stats dict
        if base_zstat not in stats_dict:
            msg = f"{base_zstat} not defined in Providentia's statistical library. Cannot calculate statistic."
            show_message(self, msg)
            return

        # if only 1 label passed and stat is a bias statistic then throw error
        elif (z_statistic_sign == 'bias') & (labelb == ''):
            msg = "Calculating a bias statistic, and only 1 label is set. Cannot calculate statistic."
            show_message(self, msg)
            return

        # if calculating bias stat but temporal_colocation is not active, then throw error
        elif (z_statistic_type == 'modbias') & (not self.temporal_colocation):
            msg = f'To calculate the model bias stat {zstat}, temporal_colocation must be set to True. Cannot calculate statistic.'
            show_message(self, msg)
            return

        # throw error if both period and chunk are given
        elif (period is not None) & (chunk is not None):
            msg = f"Cannot calculate statistic when both 'period' and 'chunk' are given."
            show_message(self, msg)
            return

        # get networkspeci to calculate for
        networkspeci = self.networkspecies[0]
        if len(self.networkspecies) > 1:
            msg = "More than 1 network or species defined, can only calculate for 1. Taking {}.".format(networkspeci)
            show_message(self, msg)

        # calculate statistic
        stat = calculate_statistic(self, self, networkspeci, stat, [labela], [labelb], per_station=per_station, 
                                period=period, chunk_resolution=chunk)
        
        return stat

    def set_config(self, **kwargs):
        """
        Wrapper method to set configuration variables.

        Parameters
        ----------
        **kwargs : dict
            Optional command-line arguments that override default configuration values.

        Returns
        -------
        bool
            True if configuration is successfully set, False if there is an error.
        """

        # if have forecast active then save current model variable in memory as will want to set it again 
        # after resetting configuration variables, as are not re-reading 
        forecast_models = None
        if hasattr(self, 'forecast'): 
            if len(self.forecast) != 0:
                if self.experiments:
                    forecast_models = copy.deepcopy(self.experiments)
            
        # initialise default configuration variables
        # modified by passed arguments, if given
        self.provconf = ProvConfiguration(self, **kwargs)

        # for any passed arguments not in default Providentia variables, now set them to self
        for kwarg in kwargs:
            if kwarg not in self.provconf.var_defaults:
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
                self.logger.info(error)
                return False
        else:
            error = "Error: The configuration file must be given as an argument: e.g. 'config=...'"
            self.logger.info(error)
            return False

        # parse section
        # if section name provided, try and use that
        # otherwise take first defined section name
        self.sections = copy.deepcopy(self.parent_section_names)

        # check if configuration file has a section title
        if len(self.sections) == 0:
            error = "Error: No sections were found in the configuration file, make sure to name them using square brackets."
            self.logger.info(error)
            return False
    
        self.have_section = False
        if hasattr(self, 'section'): 
            # check that section actually exists
            if self.section in self.all_sections:
                self.have_section = True
            else:
                msg = 'Error: The section specified in the command line ({0}) does not exist.'.format(self.section)
                msg += '\nTip: For subsections, add the name of the parent section followed by an interpunct (·) '
                msg += 'before the subsection name (e.g. SECTIONA·Spain). Available: {0}'.format(self.all_sections)
                show_message(self, msg)

        if not self.have_section:
            self.section = self.sections[0]
            if len(self.sections) > 1:
                msg = "Taking first defined section ({}).".format(self.section)
                show_message(self, msg)

        # update self with section variables (if not overwritten by kwargs)
        self.section_opts = self.sub_opts[self.section]
        for k, val in self.section_opts.items():
            if k not in kwargs:
                setattr(self, k, self.provconf.parse_parameter(k, val))

        # parse subsection (if not already parsed from section)
        # if subsection name is provided, try and use that
        # otherwise take first defined subsection name
        # if have no subsections, section is set as subsection name
        have_subsection = False
        # get subsection names
        self.subsections = [subsection_name for subsection_name in self.subsection_names 
                            if self.section == subsection_name.split('·')[0]]
        self.subsections_reduced = [subsection_name.split('·')[1] for subsection_name in self.subsections]

        #give warning if have previously defined subsection in section name, but it is defined again 
        if (hasattr(self, 'subsection')) & ('·' in self.section): 
            msg = "Defined subsection {} is not taken into account as it is already passed in section {}.".format(self.subsection, self.section)
            show_message(self, msg)

        # check that subsection actually exists if defined
        elif (hasattr(self, 'subsection')) & ('·' not in self.section): 
            
            if self.subsection in self.subsections:
                have_subsection = True
            elif self.subsection in self.subsections_reduced:
                have_subsection = True
                self.subsection = self.subsections[self.subsections_reduced.index(self.subsection)]
            elif self.subsection == self.section:
                have_subsection = True
            else:
                msg = "Defined subsection {} does not exist in configuration file.".format(self.subsection)
                show_message(self, msg)

        # reduce multiple subsections to first one if none defined
        if len(self.subsections) > 0:
            if not have_subsection:
                self.subsection = self.subsections[0]
                have_subsection = True
                if len(self.subsections) > 1:
                    msg = "Using the first defined subsection ({}).".format(self.subsection)
                    show_message(self, msg)
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
            
            # update subsection variables (if not overwritten by kwargs)
            for k, val in self.subsection_opts.items():
                if k not in kwargs:
                    setattr(self, k, self.provconf.parse_parameter(k, val))

        # now all variables have been parsed, check validity of those, throwing errors where necessary
        self.provconf.check_validity()

        # overwrite models variable if forecast active
        if forecast_models is not None:
            self.experiments = forecast_models

        return True

    def print_config(self, conf=None, config=None):
        """ 
        Print selected configuration file to the console.

        Parameters
        ----------
        conf : str, optional
            Configuration file name, defaults to None.
        config : str, optional
            Alternative Configuration file name, defaults to None.
            If both conf and config are None, the currently loaded configuration file is printed.
        """

        # check have valid conf
        valid_config = self.have_valid_config()
        if not valid_config:
            return

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
            msg = "The passed .conf file: '{}' does not exist.".format(conf)
            show_message(self, msg)
        # otherwise, print conf
        else:
            with open(conf, "r") as f:
                self.logger.info(f.read())

    def __str__(self):
        """
        Set the active configuration as the string representation.

        Returns
        -------
        str
            Empty string after printing the active configuration to the console.
        """

        self.print_config(conf=self.config)
        return ''

    def save(self, fname='', format='nc'):
        """ 
        Wrapper method to save current data/ metadata in memory.

        Parameters
        ----------
        fname : str, optional
            File name, defaults to ''.
        format : str, optional
            File format, defaults to 'nc'.
        """

        # check have valid conf and have loaded data
        valid_config = self.have_valid_config()
        if not valid_config:
            return
        loaded_data = self.loaded_data()
        if not loaded_data:
            return

        # set fname if not provided
        if fname == '':
            date_str = datetime.datetime.today().strftime('%Y%m%d_%H%M')
            fname = join(PROVIDENTIA_ROOT, 'saved_data/PRV_{}'.format(date_str))

        # remove extension if provided in fname
        if '.' in fname:
            fname = fname.split('.')[0]

        if format in ['conf','config','.conf']:
            fname = '{}.conf'.format(fname)
            export_configuration(self, fname)

        elif format in ['netCDF', 'netcdf', 'netCDF4', 'netcdf4', 'nc', '.nc']:
            fname = '{}.nc'.format(fname)
            export_netcdf(self, fname)

        elif format in ['npz','.npz','np','.np','npy','.npy','numpy']:
            fname = '{}.npz'.format(fname)
            export_data_npz(self, fname)

        self.logger.info('Data saved to {}'.format(fname))

    def data(self, format='nc'):
        """ 
        Wrapper method return data / metadata in specific format.

        Parameters
        ----------
        format : str, optional
            File format, defaults to 'nc'.

        Returns
        -------
        data : numpy.ndarray or xarray.Dataset
            Data returned in the requested format.
        """

        # check have valid conf and have loaded data
        valid_config = self.have_valid_config()
        if not valid_config:
            return
        loaded_data = self.loaded_data()
        if not loaded_data:
            return

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

    def variable(self, var=''):
        """
        Wrapper method to return a specific data/metadata variable.

        Parameters
        ----------
        var : str, optional
            Variable name to read, defaults to ''.

        Returns
        -------
        numpy.ndarray
            Data of the requested variable. Returns None if 
            variable is undefined or not available.
        """

        # check have valid conf and have loaded data
        valid_config = self.have_valid_config()
        if not valid_config:
            return
        loaded_data = self.loaded_data()
        if not loaded_data:
            return

        # if variable is undefined then print warning
        if var == '':
            msg = "Variable to read is undefined."
            show_message(self, msg)
            return 
        else:
            data = self.data(format='nc')
            if var not in data.variables.keys():
                # try adding networkspeci to variable, if just have 1 networkspecies
                if len(self.networkspecies) == 1:
                    test_var = '{}_{}'.format(self.networkspecies[0], var)
                    if test_var in data.variables.keys():
                        var_data = data[test_var][:]
                        return var_data
                msg = "Variable '{}' is not defined".format(var)
                show_message(self, msg)
            else:
                var_data = data[var][:]
                return var_data

    def load(self):
        """
        Wrapper method to load data into memory, apply initial filtering, 
        and mark data as initialized.
        """

        # check have valid conf
        valid_config = self.have_valid_config()
        if not valid_config:
            return

        # update filetrees
        generate_file_trees(self)

        # get valid observations in date range
        get_valid_obs_files_in_date_range(self, self.start_date, self.end_date)

        # update available models for selected fields
        get_valid_models(self, self.start_date, self.end_date, self.resolution,
                            self.network, self.species)
        
        # reset configuration variables in case new data has been downloaded
        self.valid_config = self.set_config(**self.kwargs)
        self.data_labels = [self.observations_data_label] + list(self.experiments.values())
        self.data_labels_raw = [self.observations_data_label] + list(self.experiments.keys())
        self.networkspecies = ['{}|{}'.format(network,speci) for network, speci in zip(self.network, self.species)]
            
        # read data
        self.read()  
        if self.invalid_read:
            self.logger.info('No valid data to read')
            return

        # filter
        self.reset(initialise=True)

        # set variable to know if data is in intial state or not
        self.initialised = True
    
    def loaded_data(self):
        """
        Helper method for determining if data has been read
        
        Returns
        -------
        bool
            True if data has been loaded, False otherwise.
        """

        if not hasattr(self, 'data_in_memory'):
            self.logger.info('Error: Data has not been loaded. Use the load() method')
            return False
        else:
            return True
        
    def have_valid_config(self):
        """
        Helper method for determining if a valid .conf file has been read
        
        Returns
        -------
        bool
            True if a valid configuration file has been read, False otherwise.
        """

        if not self.valid_config:
            self.logger.info("Error: A valid configuration file has not been read. Please reinitialise your Providentia object with a valid file: prv.Providentia('filename.conf')")
            return False
        else:
            return True

    def download(self, **kwargs):
        """
        Wrapper function for initialising Download class
        
        Parameters
        ----------
        **kwargs : dict
            Optional command-line arguments that override default configuration values.
        """

        # check have valid conf
        valid_config = self.have_valid_config()
        if not valid_config:
            return

        from .download import main

        kwargs['download'] = True
        parent_kwargs = copy.deepcopy(self.kwargs) 
        parent_kwargs.update(kwargs)
        main(**parent_kwargs)

    def dashboard(self, **kwargs):
        """ 
        Wrapper function for initialising Dashboard class        
        
        Parameters
        ----------
        **kwargs : dict
            Optional command-line arguments that override default configuration values.
        """

        from .dashboard import main
        
        kwargs['dashboard'] = True
        parent_kwargs = copy.deepcopy(self.kwargs) 
        parent_kwargs.update(kwargs)
        main(**parent_kwargs)

    def interpolate(self, **kwargs):
        """
        Wrapper function for initialising Interpolation class        
        
        Parameters
        ----------
        **kwargs : dict
            Optional command-line arguments that override default configuration values.
        """
        
        # check have valid conf
        valid_config = self.have_valid_config()
        if not valid_config:
            return

        from .interpolation import experiment_interpolation_submission as interpolation
        
        kwargs['interpolation'] = True
        
        unique_id = f"{random.randint(0, 999999):06d}"
        kwargs['slurm_job_id'] = unique_id

        log_path = join(
            PROVIDENTIA_ROOT,
            "logs",
            "interpolation",
            "management_logs",
            f"{unique_id}.out"
        )

        # save original stdout
        orig_stdout = sys.stdout
        with open(log_path, "w") as f:
            sys.stdout = Tee(orig_stdout, f)
            try:
                # do interpolation
                parent_kwargs = copy.deepcopy(self.kwargs) 
                parent_kwargs.update(kwargs)
                interpolation.main(**parent_kwargs)
            finally:
                # reset stdout
                sys.stdout = orig_stdout

        # reuse logger
        self.provconf.switch_logging()

    def report(self, **kwargs):
        """
        Wrapper function for initialising Report class        
        
        Parameters
        ----------
        **kwargs : dict
            Optional command-line arguments that override default configuration values.
        """

        # check have valid conf
        valid_config = self.have_valid_config()
        if not valid_config:
            return

        from .report import main
        
        kwargs['report'] = True
        parent_kwargs = copy.deepcopy(self.kwargs) 
        parent_kwargs.update(kwargs)
        main(**parent_kwargs)