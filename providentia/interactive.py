""" Class to create interactive session"""

import copy
import datetime
import json
import os
import sys

import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.projections import PolarAxes
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as fa
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QApplication, QWidget

from .configuration import load_conf
from .configuration import ProvConfiguration
from .fields_menus import (init_metadata, init_period, init_representativity, metadata_conf,
                           update_metadata_fields, update_period_fields, update_representativity_fields,
                           period_conf, representativity_conf)
from .filter import DataFilter
from .plot import Plot
from .plot_formatting import (format_plot_options, format_axis, set_axis_label, set_axis_title, 
                              harmonise_xy_lims_paradigm)
from .read import DataReader
from .read_aux import (get_ghost_observational_tree, get_nonghost_observational_tree, 
                       get_nonrelevant_temporal_resolutions, get_relevant_temporal_resolutions, 
                       get_valid_experiments, get_valid_obs_files_in_date_range)
from .statistics import (calculate_statistic, generate_colourbar, generate_colourbar_detail, 
                         get_selected_station_data, get_z_statistic_info)
from .writing import export_configuration, export_data_npz, export_netcdf

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])

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

        #set mode as interactive
        kwargs['interactive'] = True

        # load statistical jsons
        self.basic_stats = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/basic_stats.json')))
        self.expbias_stats = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/experiment_bias_stats.json')))

        #set configuration variables, as well as any other defined variables
        valid_config = self.set_config(**kwargs)
        # if configuration read was not valid then return here
        if not valid_config:
            return

        # create dictionary of all available observational GHOST data
        self.all_observation_data = get_ghost_observational_tree(self)

        # load dictionary with non-GHOST esarchive files to read
        nonghost_observation_data_json = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/nonghost_files.json')))
        # merge to existing GHOST observational data dict if we have the path
        if self.nonghost_root is not None:
            nonghost_observation_data = get_nonghost_observational_tree(self, nonghost_observation_data_json)
            self.all_observation_data = {**self.all_observation_data, **nonghost_observation_data}

        # initialise DataReader class
        self.datareader = DataReader(self)

        # check for self defined plot characteristics file
        if self.plot_characteristics_filename == '':
            self.plot_characteristics_filename = os.path.join(PROVIDENTIA_ROOT, 'settings/plot_characteristics_interactive.json')
        self.plot_characteristics_templates = json.load(open(self.plot_characteristics_filename))

        # initialise Plot class
        self.plot = Plot(read_instance=self, canvas_instance=self)

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

        # if have no experiments, force temporal colocation to be False
        if len(self.experiments) == 0:
            self.temporal_colocation = False  

        # read data
        self.read()  

        # filter
        self.reset_filter()

    def read(self):
        """wrapper method to read data"""

        print('Reading data')

        # read data
        self.datareader.read_setup(['reset'])

        if self.invalid_read:
            print('No valid data to read')
            return

    def filter(self):
        """wrapper method to filter data"""

        print('Filtering data')

        # filter data
        DataFilter(self)

        # get selected station data
        get_selected_station_data(read_instance=self, canvas_instance=self, networkspecies=self.networkspecies)

    def reset_filter(self):
        """wrapper method to reset filter data"""

        print('Resetting filter')

        # reset representativity fields        
        init_representativity(self)
        update_representativity_fields(self)
        representativity_conf(self)

        # reset period fields 
        init_period(self)
        update_period_fields(self)
        period_conf(self)

        # reset metadata
        init_metadata(self)
        update_metadata_fields(self)
        metadata_conf(self)

        #re-filter 
        self.filter()

    def make_plot(self, plot, data_labels=[], labela='', labelb='', title=None, xlabel=None, ylabel=None, 
                  cb=True, legend=True, map_extent=None, plot_options=[], return_plot=False, set_obs_legend=True, 
                  format={}):
        """wrapper method to make a Providentia plot"""

        # get base plot type (no plot options), and plot type (with plot options)
        base_plot_type = copy.deepcopy(plot)
        if len(plot_options) > 0:
            plot_type = '{}_{}'.format(base_plot_type, '_'.join(plot_options))
        else:
            plot_type = copy.deepcopy(plot)
        # get zstat for required plots
        base_plot_type_split = base_plot_type.split('-')
        if (len(base_plot_type_split) > 1) & (base_plot_type != 'periodic-violin'):
            base_plot_type = base_plot_type_split[0]
            zstat = base_plot_type_split[1]
        else:
            zstat = None

        # get zstat information 
        zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type=plot_type) 

        # get data labels for plot
        if len(data_labels) == 0:
            data_labels = copy.deepcopy(self.data_labels)
        # if any passed data labels are not available then pass warning and return
        else:
            invalid_data_labels = [data_label for data_label in data_labels if data_label not in self.data_labels]
            if len(invalid_data_labels) > 0:
                print("Warning: Passed data labels {} are not available. Not making plot.".format(invalid_data_labels))
                return

        # set plot characteristics
        self.plot_characteristics = dict()
        valid_plot = self.plot.set_plot_characteristics([plot_type], format=format, data_labels=data_labels)

        # if after setting plot charateristics it has been determined plot is not valid, then return
        if not valid_plot:
            return

        # if map extent passed not passed as argument, set it as value from .conf in memory
        if (not map_extent) and (self.map_extent):
            map_extent = copy.deepcopy(self.map_extent)

        # create figure and axis/axes for plot
        # differentiation in approach between jupyter and standard call for creation of figure geometry
        if jupyter_session:
            fig = plt.figure(figsize=self.plot_characteristics[plot_type]['figsize'])
        else:
            app = QApplication(sys.argv)
            desktop = app.desktop()
            screenRect = desktop.screenGeometry()
            width = screenRect.width()
            height = screenRect.height()
            dpi = plt.rcParams['figure.dpi']
            px = 1.0/dpi
            fig = plt.figure(figsize=(width*px,height*px))

        #create axes
        if base_plot_type == 'map':
            ax = fig.add_subplot(111, projection=self.plotcrs)
        else:
            ax = fig.add_subplot(111)

        if base_plot_type in ['periodic', 'periodic-violin']:
            gs = gridspec.GridSpecFromSubplotSpec(100, 100, subplot_spec=ax)
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
        elif base_plot_type != 'legend':
            func = getattr(self.plot, 'make_{}'.format(base_plot_type.split('-')[0]))
         
        # set boolean on whether to plot obs in legend or not
        if (base_plot_type == 'scatter') or ('bias' in plot_options) or (z_statistic_sign == 'bias'):
            set_obs_legend = False

        # get networkspeci to plot
        networkspeci = self.networkspecies[0]
        speci = networkspeci.split('|')[-1]
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
                print("Warning: No specific data labels set, plotting first available data label: {}.".format(z1_label))
            # labelb defined but labela for some reason, set labela to be labelb, and labelb empty str 
            elif (labela == ''):
                labela = copy.deepcopy(labelb)
                labelb = ''
            func(relevant_ax, networkspeci, self.plot_characteristics[plot_type], zstat=zstat, 
                 labela=labela, labelb=labelb, plot_options=plot_options)
        # periodic plot
        elif base_plot_type == 'periodic':
            func(grid_dict, networkspeci, data_labels, self.plot_characteristics[plot_type], zstat=zstat, 
                 plot_options=plot_options)
        # make statsummary plot
        elif base_plot_type == 'statsummary':
            relevant_zstats = self.active_statsummary_stats['basic']
            func(relevant_ax, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                 zstats=relevant_zstats, statsummary=True, plot_options=plot_options) 

            func(relevant_axis, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                         statsummary=True, plot_options=plot_options, subsection=self.subsection, 
                         plotting_paradigm=plotting_paradigm, stats_df=stats_df)     

        # make heatmap / table plot
        elif base_plot_type in ['heatmap','table']:                
            func(relevant_ax, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                 plot_options=plot_options, subsection=self.subsection, 
                 plotting_paradigm=plotting_paradigm, stats_df=stats_df)
                           
        # other plots
        else: 
            func(relevant_ax, networkspeci, data_labels, self.plot_characteristics[plot_type], 
                 plot_options=plot_options)

        # get number of total available stations, and individual station information if just have 1 station
        if (self.temporal_colocation) and (len(data_labels) > 1):
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

        # set title, xlabel and ylabel for plots

        # set xlabel / ylabel
        if base_plot_type == 'periodic':
            if not ylabel:
                
                if 'ylabel' in self.plot_characteristics[plot_type]:
                    ylabel = self.plot_characteristics[plot_type]['ylabel']['ylabel']
                else:
                    ylabel = ''

                if ylabel == '':                            
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
            if zstat:
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
                        else:
                            title = '{} at {} stations'.format(stat_label, n_stations)

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

        # overwrite passed xlabels and ylabels
        if title:
            #self.plot_characteristics[plot_type]['axis_title']['label'] = title
            set_axis_title(self, relevant_ax, title, self.plot_characteristics[plot_type])
        if xlabel:
            #self.plot_characteristics[plot_type]['xlabel']['xlabel'] = xlabel
            set_axis_label(relevant_ax, 'x', xlabel, self.plot_characteristics[plot_type])
        if ylabel:
            #self.plot_characteristics[plot_type]['ylabel']['ylabel'] = ylabel
            set_axis_label(relevant_ax, 'y', ylabel, self.plot_characteristics[plot_type])

        # format plot axis/axes
        format_axis(self, self, relevant_ax, base_plot_type, self.plot_characteristics[plot_type], 
                    map_extent=map_extent)

        # format plot options
        format_plot_options(self, self, relevant_ax, [data_labels], networkspeci, base_plot_type, plot_type, 
                            plot_options)                         

        # handle harmonisation of axes
        if base_plot_type not in ['legend', 'metadata', 'map', 'taylor']:
            if base_plot_type == 'scatter':
                harmonise_xy_lims_paradigm(self, self, relevant_ax, base_plot_type, 
                                           self.plot_characteristics[plot_type], plot_options, relim=True)
            else:
                harmonise_xy_lims_paradigm(self, self, relevant_ax, base_plot_type, 
                                           self.plot_characteristics[plot_type], plot_options, relim=True, 
                                           autoscale=True)

        # make legend (embedded on plot axis)
        if (legend) & (base_plot_type != 'legend'):
            if 'legend' in self.plot_characteristics[plot_type]:
                legend_handles = self.make_legend(plot_type, data_labels=data_labels, set_obs=set_obs_legend)
                if base_plot_type in ['periodic', 'periodic-violin']:
                    try:
                        ax_to_plot = self.plot_characteristics[plot_type]['legend']['handles']['ax']
                    except:
                        print("Warning: axis to plot legend on not defined for plot type in plot_characteristics_interactive.json, or passed via 'format' argument.\n\Taking first available axis.")
                        ax_to_plot = self.relevant_temporal_resolutions[0]
                    if ax_to_plot not in self.relevant_temporal_resolutions:
                        print("Warning: defined axis to plot legend on not available for data resolution of read data.\nInstead, taking first available axis.")
                        ax_to_plot = self.relevant_temporal_resolutions[0]
                    relevant_ax[ax_to_plot].legend(**legend_handles)
                else:
                    relevant_ax.legend(**legend_handles)

        # make colourbar (embedded on plot axis)
        if cb:
            if 'cb' in self.plot_characteristics[plot_type]:
                self.make_colourbar(fig, relevant_ax, zstat, speci, plot_type)

        # if return_plot is passed then return plot axis/axes
        if return_plot:
            return relevant_ax
        # otherwise show plot
        else:
            plt.show()

    def make_colourbar(self, fig, plot_ax, stat, speci, plot_type):
        """wrapper method to make colourbar"""

        # create cb axis
        cb_ax = fig.add_axes(self.plot_characteristics[plot_type]['cb']['position'])
        cb_ax.set_rasterized(True)

        # generate colourbar
        generate_colourbar(self, [plot_ax], [cb_ax], stat, self.plot_characteristics[plot_type], speci)

    def make_legend(self, plot_type, data_labels=None, set_obs=True):
        """wrapper method to make legend"""
        
        if plot_type == 'legend':
            legend_characteristics = self.plot_characteristics['legend']
        elif 'legend' in self.plot_characteristics[plot_type]:
            legend_characteristics = self.plot_characteristics[plot_type]['legend']
        else:
            print("Warning: 'legend' not defined for plot type in plot_characteristics_interactive.json")
            return

        legend_handles = self.plot.make_legend_handles(legend_characteristics, data_labels=data_labels, set_obs=set_obs)
        return legend_handles['plot']

    def calculate_statistic(self, stat='', labela='', labelb='', per_station=False):
        """wrapper method to calculate statistic/s"""

        # if no specific labels defined then take first data label and give warning
        if (labela == '') & (labelb == ''):
            labela = data_labels[0]
            print("Warning: No specific data labels set, plotting first available data label: {}.".format(z1_label))
        # labelb defined but labela for some reason, set labela to be labelb, and labelb empty str 
        elif (labela == ''):
            labela = copy.deepcopy(labelb)
            labelb = ''

        if per_station:
            stat, active_map_valid_station_inds = calculate_statistic(self, self, networkspeci, stat, [labela], [labelb], map=True)
        else:
            stat = calculate_statistic(self, self, networkspeci, stat, [labela], [labelb])
        return stat

    def set_config(self, **kwargs):
        """wrapper method to set configuration variables"""

        # initialise default configuration variables
        # modified by passed arguments, if given
        provconf = ProvConfiguration(self, **kwargs)

        # update variables from config file
        if self.config != '':
            read_conf = False 
            if os.path.exists(self.config):
                read_conf = True
            else:
                if os.path.exists(os.path.join(self.config_dir, self.config)):
                    self.config = os.path.join(self.config_dir, self.config)
                    read_conf = True
            if read_conf:
                load_conf(self, self.config)
                self.from_conf = True
            else:
                error = 'Error: The path to the configuration file passed as an argument does not exist.'
                return False
        else:
            error = "Error: The configuration file must be given as an argument: e.g. 'config=...'"
            return False

        # parse section
        # if section name provided, try and use that
        # otherwise take first defined section name
        have_section = False
        if hasattr(self, 'section'): 
            # check that section actually exists
            if self.section in self.parent_section_names:
                have_section = True
            else:
                print("Warning: Defined section {} does not exist in configuration file.".format(self.section))
        if not have_section:
            self.section = self.parent_section_names[0]
            if len(self.parent_section_names) > 1:
                print("Warning: Taking first defined section ({}) to be read.".format(self.section))

        # update self with section variables
        self.section_opts = self.sub_opts[self.section]
        for k, val in self.section_opts.items():
            setattr(self, k, provconf.parse_parameter(k, val))

        # parse subsection
        # if subsection name is provided, try and use that
        # otherwise take first defined subsection name
        # if have no subsections, section is set as subsection name
        have_subsection = False
        # get subsection names
        self.child_subsection_names = [subsection_name for subsection_name in self.subsection_names 
                                        if self.section == subsection_name.split('·')[0]]
        if hasattr(self, 'subsection'): 
            # check that subsection actually exists
            if self.subsection in self.child_subsection_names:
                have_subsection = True
            else:
                print("Warning: Defined subsection {} does not exist in configuration file.".format(self.subsection))

        if len(self.child_subsection_names) > 0:
            if not have_subsection:
                self.subsection = self.child_subsection_names[0]
                have_subsection = True
                if len(self.child_subsection_names) > 1:
                    print("Warning: Taking first defined subsection ({}) to be read.".format(self.subsection))
        else:
            self.subsection = [self.section]

        # update self with subsection variables (if have subsection) 
        if have_subsection:   
            # get subsection variables
            self.subsection_opts = self.sub_opts[self.subsection]
            # ensure all fixed section variables defined in subsection have same value as current section variables
            self.subsection_opts = {k: (self.section_opts[k] if k in self.fixed_section_vars else val) 
                                    for (k, val) in self.subsection_opts.items()}
            # update subsection variables
            for k, val in self.subsection_opts.items():
                setattr(self, k, provconf.parse_parameter(k, val))

        # now all variables have been parsed, check validity of those, throwing errors where necessary
        provconf.check_validity()

        return True

    def print_config(self, conf=None, config=None):
        """print selected config file to console"""

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
        """wrapper method to select specific station/s"""
        
        if type(station) == 'str':
            stations_to_keep = [station]
        else:
            stations_to_keep = station
        self.metadata_menu['STATION MISCELLANEOUS']['station_reference']['checkboxes']['keep_selected'] = stations_to_keep

        #filter for station/s    
        self.filter()

    def save(self, fname='', format='nc'):
        """wrapper method to save current data/ metadata in memory"""

        # set fname if not provided
        if fname == '':
            date_str = datetime.datetime.today().strftime('%Y%m%d_%H%M')
            fname = os.path.join(PROVIDENTIA_ROOT, 'saved_data/PRV_{}'.format(date_str))

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
        """wrapper method return data / metadata in specific format"""

        # set temporary fname for writing
        temporary_fname = os.path.join(PROVIDENTIA_ROOT, 'saved_data/temp')

        if format in ['netCDF', 'netcdf', 'netCDF4', 'netcdf4', 'nc', '.nc']:
            data = export_netcdf(self, temporary_fname, set_in_memory=True)

        elif format in ['xr', '.xr', 'xarr', 'xarray','Xarray']:
            data = export_netcdf(self, temporary_fname, set_in_memory=True, xarray=True)

        elif format in ['npz','.npz','np','.np','npy','.npy','numpy']:
            data = export_data_npz(self, temporary_fname, set_in_memory=True)

        return data

    def get_var(self, var=''):
        """wrapper method to return specific data / metadata variable"""

        # if variable is undefined then print warning
        if var == '':
            print("Warning: Variable to read is undefined.")
            return 
        else:
            data = self.get_data(format='npz')
            var_data = data[var]
            return var_data