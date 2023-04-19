from .filter import DataFilter
from .statistics import *
from .plot import Plot

import os
import copy
import json
import sys
import time
import datetime
import math
from weakref import WeakKeyDictionary
import time

import matplotlib
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg \
        as FigureCanvas
from matplotlib.backend_bases import FigureCanvasBase
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.path import Path
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec
import matplotlib.style as mplstyle
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
from PyQt5 import QtCore, QtGui, QtWidgets 
from .dashboard_aux import set_formatting, ComboBox, CheckableComboBox, LassoSelector
from .aux import get_relevant_temporal_resolutions, show_message

# make sure that we are using Qt5 backend with matplotlib
matplotlib.use('Qt5Agg')
register_matplotlib_converters()

# use matplotlib fast style: https://matplotlib.org/stable/users/explain/performance.html
mplstyle.use('fast')

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))
formatting_dict = json.load(open(os.path.join(CURRENT_PATH, 'conf/stylesheet.json')))

class MPLCanvas(FigureCanvas):
    """ Class that handles the creation and updates of
        a matplotlib canvas, and associated subplots.
    """

    def __init__(self, read_instance):
        """ Initialise the MPL canvas. """

        # create figure and canvas objects
        self.figure = Figure(dpi=100)
        FigureCanvas.__init__(self, self.figure)

        # add passed arguments to self
        self.read_instance = read_instance

        # get characteristics per plot type
        self.plot_characteristics_templates = self.read_instance.plot_characteristics_templates
        self.plot_characteristics = {}

        # initialise plot elements
        self.plot_elements = {}

        # add general plot characteristics to self
        for k, val in self.plot_characteristics_templates['general'].items():
            setattr(self, k, val)

        # initialise some key vars
        self.filter_data = None

        # initialise Plot class
        self.plot = Plot(read_instance=self.read_instance, canvas_instance=self)

        # setup gridding of canvas
        self.gridspec = gridspec.GridSpec(self.gridspec_nrows, self.gridspec_ncols)
        self.gridspec.update(**self.gridspec_format)

        # define all possible plots
        self.all_plots = ['legend', 'map', 'timeseries', 'periodic-violin', 'periodic', 
                          'metadata', 'distribution', 'scatter', 'statsummary', 'boxplot']

        # define all possible plots in layout options
        self.layout_options = ['None', 'boxplot', 'distribution', 'metadata', 'periodic', 
                               'periodic-violin', 'scatter', 'statsummary', 'timeseries']

        # parse active dashboard plot string        
        if isinstance(self.read_instance.active_dashboard_plots, str):
            self.read_instance.active_dashboard_plots = [c.strip() for c in self.read_instance.active_dashboard_plots.split(',')]
        
        # stop running if plot type in active_dashboard_plots does not exist
        for plot_type in self.read_instance.active_dashboard_plots:
            if plot_type not in self.all_plots:
                error = "Error: Plot type {0} is not an option. ".format(plot_type)
                error += "The available plots are: {0}.".format(self.all_plots[2:])
                sys.exit(error)

        # initialize layout positions
        self.read_instance.position_1 = 'map'
        self.read_instance.position_2 = self.read_instance.active_dashboard_plots[0]
        self.read_instance.position_3 = self.read_instance.active_dashboard_plots[1]
        self.read_instance.position_4 = self.read_instance.active_dashboard_plots[2]
        self.read_instance.position_5 = self.read_instance.active_dashboard_plots[3]

        # update plot characteristics for all plots
        for plot_type in self.all_plots:
            # for plot types with zstat, initialise with default zstat (Mean)
            if plot_type in ['map', 'periodic']:
                self.plot.set_plot_characteristics([plot_type], zstat='Mean')
            else:
                self.plot.set_plot_characteristics([plot_type])

        # create map, colorbar and legend plot axes
        self.plot_axes = {}
        self.plot_axes['map'] = self.figure.add_subplot(self.gridspec.new_subplotspec((2, 0), 
                                                        rowspan=44, colspan=42), projection=self.plotcrs)
        self.plot_axes['cb'] = self.figure.add_axes([0.0455, 0.536, 0.3794, 0.02])
        self.plot_axes['legend'] = self.figure.add_subplot(self.gridspec.new_subplotspec((0, 47), 
                                                           rowspan=8, colspan=53))

        # add settings menus
        self.generate_interactive_elements()

        # create rest of plot axes (default: timeseries, statsummary, distribution, periodic)
        for position, plot_type in enumerate(self.read_instance.active_dashboard_plots):
            self.read_instance.update_plot_axis(self, position + 2, plot_type)
        
        # update layout fields
        self.read_instance.update_layout_fields(self)

        # initialise variable of valid station indices plotted on map as empty list
        self.active_map_valid_station_inds = np.array([], dtype=np.int)

        # setup blocker for picker events
        self.figure.canvas.mpl_connect('axes_enter_event', self.picker_block_func)

        # setup legend line selection
        self.figure.canvas.mpl_connect('pick_event', self.legend_picker_func)

        # setup interactive lasso on map
        self.lasso_left = LassoSelector(self, self.plot_axes['map'], onselect=self.onlassoselect_left,
                                        useblit=False, props=self.lasso, button=[1])
        self.lasso_right = LassoSelector(self, self.plot_axes['map'], onselect=self.onlassoselect_right,
                                         useblit=False, props=self.lasso, button=[3])

        # setup station annotations
        self.create_station_annotation()
        self.map_annotation_disconnect = False
        self.map_annotation_event = self.figure.canvas.mpl_connect('motion_notify_event', self.hover_map_annotation)

        # setup zoom on scroll wheel on map
        self.lock_zoom = False
        self.figure.canvas.mpl_connect('scroll_event', self.zoom_map_func)

        # format axes for map, legend and active_dashboard_plots
        for plot_type in ['map', 'legend'] + self.read_instance.active_dashboard_plots:
            self.plot.format_axis(self.plot_axes[plot_type], plot_type, self.plot_characteristics[plot_type])

        # hide all plot axes
        for plot_type, ax in self.plot_axes.items():
            self.remove_axis_elements(ax, plot_type)

    def update_MPL_canvas(self):
        """ Function that updates MPL canvas upon clicking
            the 'READ' button, and when colocating data.
        """

        # reset relative index lists of selected station on map as empty lists
        self.previous_relative_selected_station_inds = np.array([], dtype=np.int)
        self.relative_selected_station_inds = np.array([], dtype=np.int)
        self.absolute_selected_station_inds = np.array([], dtype=np.int)

        # reset plot_elements
        self.plot_elements = {}
        self.plot_elements['data_labels_active'] = []
        for data_label in self.read_instance.data_labels:
            self.plot_elements['data_labels_active'].append(data_label)

        # update plotted map z statistic
        self.update_map_z_statistic()

        # plot domain edges on map and legend if have valid data
        if len(self.active_map_valid_station_inds) > 0:

            # plot experiment grid domain edges on map
            self.update_experiment_domain_edges()

            # update legend
            self.update_legend()

        # draw changes
        self.figure.canvas.draw()

        return None

    def reset_ax_navigation_toolbar_stack(self, ax):
        """ Function which resets the navigation toolbar stack
            for a given axis with the current view limits.
        """

        # get appropriate axes for nested axes
        axs_to_reset = []
        if isinstance(ax, dict):
            for relevant_temporal_resolution, sub_ax in ax.items():
                if relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions:
                    axs_to_reset.append(sub_ax)
        else:
            axs_to_reset.append(ax)

        # check if have axes dictionaries in stack list
        for ax_to_reset in axs_to_reset:

            if len(self.read_instance.navi_toolbar._nav_stack) == 0:
                # if don't have an axes dictionary in stack list, create one with current
                # axis in dictionary with current view limits
                self.read_instance.navi_toolbar._nav_stack.push(
                    WeakKeyDictionary({ax_to_reset: (ax_to_reset._get_view(), (ax_to_reset.get_position(True).frozen(), ax_to_reset.get_position().frozen()))}))

            # if have existing axes dictionaries in stack list, iterate through stack list
            # removing given axis from all stack list dictionaries
            else:
                for axes_dict in self.read_instance.navi_toolbar._nav_stack:
                    if ax_to_reset in axes_dict.keyrefs():
                        axes_dict.pop(ax_to_reset)

                # now add axis to first dictionary in stack, with the current view limits
                self.read_instance.navi_toolbar._nav_stack[0][ax_to_reset] = \
                    (ax_to_reset._get_view(), (ax_to_reset.get_position(True).frozen(), ax_to_reset.get_position().frozen()))

        return None

    def handle_data_filter_update(self):
        """ Function which handles updates of data filtering. """

        # return if nothing has been loaded yet
        if not hasattr(self.read_instance, 'data_in_memory'):
            return

        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        if self.filter_data is None:
            self.filter_data = DataFilter(self.read_instance)
        else:
            self.filter_data.filter_all()
            self.update_active_map()

        QtWidgets.QApplication.restoreOverrideCursor()

        return None

    def handle_resampling_update(self):
        """ Function which handles updates of resampling. """
        
        if not self.read_instance.block_MPL_canvas_updates:

            # if there are experiments selected, temporal colocation must be on
            if ((self.read_instance.temporal_colocation == False) 
                and (len(self.read_instance.data_labels) > 1)
                and (self.read_instance.cb_resampling_resolution.currentText() != 'None')):
                    msg = 'It is not possible to resample the data without first activating the temporal colocation.'
                    show_message(self.read_instance, msg)
                    self.read_instance.block_config_bar_handling_updates = True
                    self.read_instance.cb_resampling_resolution.setCurrentText('None')
                    self.read_instance.block_config_bar_handling_updates = False
                    return

            # activate or deactivate resampling
            self.read_instance.resampling_resolution = self.read_instance.cb_resampling_resolution.currentText()
            if self.read_instance.resampling_resolution == 'None':
                self.read_instance.resampling = False
            else:
                self.read_instance.resampling = True

            # if have selected stations on map, then now remake plots
            if hasattr(self, 'relative_selected_station_inds'):
                if len(self.relative_selected_station_inds) > 0:
                    
                    # update mouse cursor to a waiting cursor
                    QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

                    # update associated plots with selected stations
                    self.update_associated_active_dashboard_plots()

                    # draw changes
                    self.figure.canvas.draw()

                    # restore mouse cursor to normal
                    QtWidgets.QApplication.restoreOverrideCursor()

        return None

    def update_active_map(self):
        """ Function that updates plotted map z statistic and updates associated plots. """

        if not self.read_instance.block_MPL_canvas_updates:

            # update plotted map z statistic
            self.update_map_z_statistic()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds,
                                  self.relative_selected_station_inds):
                self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw()

        return None

    def handle_aggregation_update(self):
        """ Function that handles the update of the MPL canvas
            when we change the way we aggregate the stations data upon clicking APPLY button.
        """

        if not self.read_instance.block_MPL_canvas_updates:

            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # get statistic
            self.read_instance.selected_aggregation_statistic = self.read_instance.cb_aggregation_statistic.currentText()
            
            # update statistic in memory
            self.read_instance.aggregation_statistic = self.read_instance.selected_aggregation_statistic 

            # update associated plots with selected stations
            self.update_associated_active_dashboard_plots()
            
            # draw changes
            self.figure.canvas.draw()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()

        return None
    
    def handle_temporal_colocate_update(self):
        """ Function that handles the update of the MPL canvas
            with colocated data upon checking of the temporal colocate checkbox.
        """

        if not self.read_instance.block_MPL_canvas_updates:
            
            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # if only have 1 data array in memory (i.e. observations), no colocation is possible,
            # therefore set colocation to be False, and return
            if len(self.read_instance.data_labels) == 1:
                self.read_instance.temporal_colocation = False
                return

            # else, if have loaded experiment data, check if colocate checkbox is checked or unchecked
            check_state = self.read_instance.ch_colocate.checkState()

            # update variable to inform plotting functions whether to use colocated data/or not
            if check_state == QtCore.Qt.Checked:
                self.read_instance.temporal_colocation = True
            else:
                self.read_instance.temporal_colocation = False

            # update map z statistic/ periodic statistic comboboxes (without updating canvas)
            self.read_instance.block_MPL_canvas_updates = True
            self.handle_map_z_statistic_update()
            self.handle_periodic_statistic_update()
            self.read_instance.block_MPL_canvas_updates = False

            # update layout fields
            self.read_instance.update_layout_fields(self)

            # update plotted map z statistic
            self.update_map_z_statistic()

            # update associated plots with selected stations
            self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()

        return None

    def unselect_map_checkboxes(self):
        """ Function to uncheck All, Intersect and Extent checkboxes without updating canvas. """

        self.read_instance.block_MPL_canvas_updates = True

        # if select all stations checkbox is checked then uncheck it
        if self.read_instance.ch_select_all.checkState() == QtCore.Qt.Checked:
            self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)

        # if select intersect stations checkbox is checked then uncheck it
        elif self.read_instance.ch_intersect.checkState() == QtCore.Qt.Checked:
            self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)

        # if select extent stations checkbox is checked then uncheck it
        elif self.read_instance.ch_extent.checkState() == QtCore.Qt.Checked:
            self.read_instance.ch_extent.setCheckState(QtCore.Qt.Unchecked)

        self.read_instance.block_MPL_canvas_updates = False

    def update_map_z_statistic(self):
        """ Function that updates plotted z statistic on map, with colourbar. """

        # remove axis elements from map/cb
        self.remove_axis_elements(self.plot_axes['map'], 'map')
        self.remove_axis_elements(self.plot_axes['cb'], 'cb')

        # get zstat name from combobox 
        base_zstat = self.map_z_stat.currentText()
        zstat = get_z_statistic_comboboxes(base_zstat, second_data_label=self.map_z2.currentText())

        # calculate map z statistic (for selected z statistic) --> updating active map valid station indices
        self.z_statistic, active_map_valid_station_inds = calculate_z_statistic(self.read_instance, 
                                                                                self.map_z1.currentText(), 
                                                                                self.map_z2.currentText(), 
                                                                                zstat, 
                                                                                self.read_instance.networkspeci)
        self.active_map_valid_station_inds = active_map_valid_station_inds 
        
        # update absolute selected plotted station indices with respect to new active map valid station indices
        self.absolute_selected_station_inds = np.array(
            [np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind in
             self.relative_selected_station_inds if selected_ind in self.active_map_valid_station_inds],
            dtype=np.int)

        # if have no valid active map indices, reset absolute/relative
        # selected station indices to be empty lists
        # also uncheck select all/intersect/extent checkboxes
        if len(self.active_map_valid_station_inds) == 0:

            # unselect all/intersect/extent checkboxes
            self.read_instance.block_MPL_canvas_updates = True
            self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
            self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
            self.read_instance.ch_extent.setCheckState(QtCore.Qt.Unchecked)
            self.read_instance.block_MPL_canvas_updates = False

            # clear previously selected relative/absolute station indices
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)
            self.relative_selected_station_inds = np.array([], dtype=np.int)
            self.absolute_selected_station_inds = np.array([], dtype=np.int)
            self.absolute_non_selected_station_inds = np.array([], dtype=np.int)

            # plot map with 0 stations
            self.plot.make_map(self.plot_axes['map'], self.read_instance.networkspeci, self.z_statistic, 
                               self.plot_characteristics['map'])

            # activate map
            self.activate_axis(self.plot_axes['map'], 'map')

        # otherwise plot valid active map stations on map
        else:
            # if any of the currently selected stations are not in the current active map
            # valid station indices --> unselect selected stations (and associated plots)
            # also uncheck select all/intersect/extent checkboxes
            if not np.all(np.in1d(self.relative_selected_station_inds, self.active_map_valid_station_inds)):
                
                # unselect all/intersect/extent checkboxes
                self.unselect_map_checkboxes()
                
                # reset relative/absolute selected station indices to be empty lists
                self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)
                self.relative_selected_station_inds = np.array([], dtype=np.int)
                self.absolute_selected_station_inds = np.array([], dtype=np.int)

            # get absolute non-selected station inds
            self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                                 self.absolute_selected_station_inds))[0]

            # plot new station points on map - coloured by currently active z statisitic, setting up plot picker
            self.plot.make_map(self.plot_axes['map'], self.read_instance.networkspeci, self.z_statistic, self.plot_characteristics['map'])

            # create 2D numpy array of plotted station coordinates
            self.map_points_coordinates = np.vstack((self.read_instance.station_longitudes[self.read_instance.networkspeci][self.active_map_valid_station_inds], 
                                                     self.read_instance.station_latitudes[self.read_instance.networkspeci][self.active_map_valid_station_inds])).T

            # generate colourbar
            generate_colourbar(self.read_instance, [self.plot_axes['map']], [self.plot_axes['cb']], zstat, 
                               self.plot_characteristics['map'], self.read_instance.species[0])

            # activate map/cb axes
            self.activate_axis(self.plot_axes['map'], 'map')
            self.activate_axis(self.plot_axes['cb'], 'cb')

        # update plot options
        self.update_plot_options(plot_types=['map'])

        # re-draw (needed to update plotted colours before update_map_station_selection)
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

        # update map selection appropriately for z statistic
        self.update_map_station_selection()

        return None

    def update_map_station_selection(self):
        """ Function that updates the visual selection of stations on map. """

        # update map title
        if len(self.relative_selected_station_inds) == 1:
            axis_title_label = '{} Selected'.format(
                self.read_instance.station_references[self.read_instance.networkspeci][self.relative_selected_station_inds[0]])
        else:
            axis_title_label = '{} Selected Stations of {} Available'.format(
                len(self.relative_selected_station_inds), len(self.active_map_valid_station_inds))
        self.plot.set_axis_title(self.plot_axes['map'], axis_title_label, self.plot_characteristics['map'])
        self.plot_characteristics['map']['axis_title']['label'] = axis_title_label

        # reset alphas and marker sizes of stations (if have some stations on map)
        if len(self.active_map_valid_station_inds) > 0:
            
            # set markersize of all stations (initally assuming zero stations are selected)
            markersizes = np.full(len(self.active_map_valid_station_inds), self.plot_characteristics['map']['marker_zero_stations_selected']['s'])
            
            for collection in self.plot_axes['map'].collections:
                if isinstance(collection, matplotlib.collections.PathCollection):
                    
                    rgba_tuples = collection.get_facecolor()
                    
                    # set alpha of all stations (initally assuming zero stations are selected)
                    rgba_tuples[:, -1] = self.plot_characteristics['map']['marker_zero_stations_selected']['alpha']
                    
                    # have selected stations?
                    if len(self.relative_selected_station_inds) > 0:
                        
                        # update markersize and alphas of non-selected stations
                        markersizes[self.absolute_non_selected_station_inds] = self.plot_characteristics['map']['marker_unselected']['s']
                        rgba_tuples[self.absolute_non_selected_station_inds, -1] = self.plot_characteristics['map']['marker_unselected']['alpha']
                        
                        # update markersize and alphas of selected stations
                        markersizes[self.absolute_selected_station_inds] = self.plot_characteristics['map']['marker_selected']['s']
                        rgba_tuples[self.absolute_selected_station_inds, -1] = self.plot_characteristics['map']['marker_selected']['alpha']
                    
                    # set new markersizes and alphas
                    collection.set_sizes(markersizes)
                    collection.set_facecolor(rgba_tuples)
        
        # redraw points
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def update_periodic_statistic(self, zstat):

        # update periodic plot/s if have some stations selected on map
        if len(self.relative_selected_station_inds) > 0:

            # clear and turn off all relevant axes before updating
            self.remove_axis_elements(self.plot_axes['periodic'], 'periodic')

            # get zstat information 
            zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat)
            
            # set plot options as empty
            plot_options = []

            # set new ylabel
            if z_statistic_type == 'basic':
                ylabel = basic_stats[base_zstat]['label']
                ylabel_units = basic_stats[base_zstat]['units']
            else:
                ylabel = expbias_stats[base_zstat]['label']
                ylabel_units = expbias_stats[base_zstat]['units']
            if ylabel_units == 'measurement_units':
                ylabel_units = '[{}]'.format(self.read_instance.measurement_units[self.read_instance.species[0]]) 
            if ylabel_units != '':
                ylabel = copy.deepcopy(ylabel_units)

            # update periodic plot
            first_data_label = True
            for data_label in self.selected_station_data[self.read_instance.networkspeci]:
                # skip observational array if bias stat
                if (z_statistic_sign == 'bias') & (data_label == 'observations'):
                    continue
                self.plot.make_periodic(self.plot_axes['periodic'], self.read_instance.networkspeci, data_label, 
                                        self.plot_characteristics['periodic'], zstat=zstat, 
                                        first_data_label=first_data_label)
                first_data_label = False

            # reset axes limits (harmonising across subplots for periodic plots) 
            self.plot.harmonise_xy_lims_paradigm(self.plot_axes['periodic'], 'periodic', 
                                                    self.plot_characteristics['periodic'], 
                                                    plot_options, relim=True, autoscale=True)

            # set ylabel
            self.plot.set_axis_label(self.plot_axes['periodic'], 'y', ylabel, self.plot_characteristics['periodic'])

            # activate axis
            self.activate_axis(self.plot_axes['periodic'], 'periodic')

            # reset navigation toolbar stack for plot
            self.reset_ax_navigation_toolbar_stack(self.plot_axes['periodic'])

            # update plot options
            self.update_plot_options(plot_types=['periodic'])

            # draw changes
            self.figure.canvas.draw()
            self.figure.canvas.flush_events()

    def update_associated_active_dashboard_plot(self, plot_type):
        """ Function that updates a plot associated with selected stations on map. """
    
        if hasattr(self, 'relative_selected_station_inds'):
            if len(self.relative_selected_station_inds) > 0:

                # get relevant axis
                ax = self.plot_axes[plot_type]

                # get options defined to configure plot 
                plot_options = []

                # get plotting function for specific plot
                if plot_type == 'statsummary':
                    func = getattr(self.plot, 'make_table')
                else:
                    func = getattr(self.plot, 'make_{}'.format(plot_type.split('-')[0]))

                # set ylabel for periodic plot
                if plot_type in ['periodic']:
                    # get currently selected periodic statistic name
                    zstat = self.periodic_stat.currentText()
                    # get zstat information 
                    zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat) 
                    # set new ylabel
                    if z_statistic_type == 'basic':
                        ylabel = basic_stats[base_zstat]['label']
                        ylabel_units = basic_stats[base_zstat]['units']
                    else:
                        ylabel = expbias_stats[base_zstat]['label']
                        ylabel_units = expbias_stats[base_zstat]['units']
                    if ylabel_units == 'measurement_units':
                        ylabel_units = self.read_instance.measurement_units[self.read_instance.species[0]] 
                    if ylabel_units != '':
                        ylabel = '[{}]'.format(ylabel_units)
                    xlabel = ''

                # create structure to store data for statsummary plot
                elif plot_type in ['statsummary']:
                    # get list of statistics to create lists for
                    relevant_zstats = self.plot_characteristics[plot_type]['basic']
                    stats_df = {relevant_zstat:[] for relevant_zstat in relevant_zstats}
                    xlabel = ''
                    ylabel = ''
                # setup xlabel / ylabel for other plot_types
                else:    
                    # set new xlabel
                    if 'xlabel' in self.plot_characteristics[plot_type]:
                        if self.plot_characteristics[plot_type]['xlabel']['xlabel'] == 'measurement_units':
                            xlabel = '[{}]'.format(self.read_instance.measurement_units[self.read_instance.species[0]])
                        else:
                            xlabel = self.plot_characteristics[plot_type]['xlabel']['xlabel']
                    else:
                        xlabel = ''

                    # set new ylabel
                    if 'ylabel' in self.plot_characteristics[plot_type]:
                        if self.plot_characteristics[plot_type]['ylabel']['ylabel'] == 'measurement_units':
                            ylabel = '[{}]'.format(self.read_instance.measurement_units[self.read_instance.species[0]])
                        else:
                            ylabel = self.plot_characteristics[plot_type]['ylabel']['ylabel']
                    else:
                        ylabel = ''

                # iterate through data array names in selected station data dictionary
                first_data_label = True
                for data_label in self.selected_station_data[self.read_instance.networkspeci]:

                    # call function to update plot
                    # periodic plot
                    if plot_type in ['periodic']:
                        # skip observational array if bias stat and data array is observations
                        if (z_statistic_sign == 'bias') & (data_label == 'observations'):
                            continue
                        func(ax, self.read_instance.networkspeci, data_label, self.plot_characteristics[plot_type], 
                             zstat=zstat, plot_options=plot_options, first_data_label=first_data_label)
                        first_data_label = False
                    # gather data for statsummary plot
                    elif plot_type in ['statsummary']:
                        for relevant_zstat in relevant_zstats:
                            stats_df[relevant_zstat].append(self.selected_station_data[self.read_instance.networkspeci][data_label]['all'][relevant_zstat][0])
                    # other plots
                    else: 
                        if plot_type == 'metadata':
                            if data_label != 'observations':
                                continue
                        if plot_type == 'scatter':
                            if data_label == 'observations':
                                continue
                        func(ax, self.read_instance.networkspeci, data_label, self.plot_characteristics[plot_type], 
                             plot_options=plot_options, first_data_label=first_data_label)
                        first_data_label = False

                # make statsummary plot
                if plot_type in ['statsummary']:
                    stats_df = pd.DataFrame(data=stats_df, index=self.selected_station_data[self.read_instance.networkspeci])
                    func(ax, stats_df, self.plot_characteristics[plot_type], plot_options=plot_options, statsummary=True)

                # reset axes limits (harmonising across subplots for periodic plots) 
                if plot_type == 'periodic-violin':
                    self.plot.harmonise_xy_lims_paradigm(ax, plot_type, self.plot_characteristics[plot_type], 
                                                         plot_options, ylim=[self.selected_station_data_min[self.read_instance.networkspeci], 
                                                                             self.selected_station_data_max[self.read_instance.networkspeci]],
                                                         relim=True, autoscale_x=True)
                elif plot_type == 'scatter':
                        self.plot.harmonise_xy_lims_paradigm(ax, plot_type, self.plot_characteristics[plot_type], 
                                                             plot_options, relim=True)
                else:
                    self.plot.harmonise_xy_lims_paradigm(ax, plot_type, self.plot_characteristics[plot_type], 
                                                         plot_options, relim=True, autoscale=True)

                # set xlabel
                self.plot.set_axis_label(ax, 'x', xlabel, self.plot_characteristics[plot_type])

                # set ylabel
                self.plot.set_axis_label(ax, 'y', ylabel, self.plot_characteristics[plot_type])

                # activate axis
                self.activate_axis(ax, plot_type)

                # reset navigation toolbar stack for plot
                self.reset_ax_navigation_toolbar_stack(ax)

                # update plot options
                self.update_plot_options(plot_types=[plot_type])

    def update_associated_active_dashboard_plots(self):
        """ Function that updates all plots associated with selected stations on map. """

        # clear all previously plotted artists from selected station plots and hide axes 
        for plot_type in self.read_instance.active_dashboard_plots:
            self.remove_axis_elements(self.plot_axes[plot_type], plot_type)

        # skip plots that need active temporal colocation and experiments data
        self.plots_to_skip = []
        for plot_type in ['scatter']:
            if plot_type in self.read_instance.active_dashboard_plots:
                if ((not self.read_instance.temporal_colocation) 
                    or ((self.read_instance.temporal_colocation) and (len(self.read_instance.experiments) == 0))): 
                    # show warning
                    if (not self.read_instance.temporal_colocation):
                        msg = f'It is not possible to make {plot_type} plots without activating the temporal colocation.'
                    else:
                        msg = f'It is not possible to make {plot_type} plots without loading experiments.'
                    show_message(self.read_instance, msg)

                    # remove plot type from options in layout options
                    self.plots_to_skip.append(plot_type)
                else:
                    # add plot type to options in layout options
                    if plot_type in self.plots_to_skip:
                        self.plots_to_skip.remove(plot_type)
            
            # update plot options in buttons
            self.read_instance.update_layout_fields(self)

        # if have selected stations on map, then now remake plots
        if hasattr(self, 'relative_selected_station_inds'):
            if len(self.relative_selected_station_inds) > 0:
                
                # put selected data for each data array into pandas dataframe
                to_pandas_dataframe(read_instance=self.read_instance, canvas_instance=self, 
                                    networkspecies=[self.read_instance.networkspeci])

                # iterate through active_dashboard_plots
                for plot_type in self.read_instance.active_dashboard_plots:
                    if plot_type not in self.plots_to_skip:
                        # if there are no temporal resolutions (only yearly), skip periodic plots
                        if ((plot_type in ['periodic', 'periodic-violin']) and 
                            (not self.read_instance.relevant_temporal_resolutions)):
                            msg = 'It is not possible to make periodic plots using annual resolution data.'
                            show_message(self.read_instance, msg)
                            continue
                    
                        # if temporal colocation is turned off, skip scatter plot
                        if (plot_type == 'scatter') and (not self.read_instance.temporal_colocation):
                            msg = 'It is not possible to make scatter plots without activating the temporal colocation.'
                            show_message(self.read_instance, msg)
                            continue

                        # update plot
                        self.update_associated_active_dashboard_plot(plot_type)
            
            # update map plot options
            self.update_plot_options(plot_types=['map'])

    def update_experiment_domain_edges(self):
        """ Function that plots grid domain edges of experiments in memory. """

        # remove grid domain polygon if previously plotted
        inds_to_remove = []
        for pol_ii, pol in enumerate(self.plot_axes['map'].patches):
            if isinstance(pol, matplotlib.patches.Polygon):
                inds_to_remove.append(pol_ii)
        self.plot_axes['map'].patches = list(np.delete(np.array(self.plot_axes['map'].patches), inds_to_remove))

        # create grid edge polygons for experiments in memory
        grid_edge_polygons = self.plot.make_experiment_domain_polygons()

        # plot grid edge polygons on map
        for grid_edge_polygon in grid_edge_polygons:
            self.plot_axes['map'].add_patch(grid_edge_polygon)

    def update_legend(self):
        """ Function that updates legend. """

        # create legend element handles
        legend_plot_characteristics = self.plot.make_legend_handles(copy.deepcopy(self.plot_characteristics['legend']))

        # plot legend
        self.legend = self.plot_axes['legend'].legend(**legend_plot_characteristics['plot'], 
                                                      prop=legend_plot_characteristics['prop'])

        # un-hide legend
        self.activate_axis(self.plot_axes['legend'], 'legend')

        # setup element picker in legend
        for legend_label in self.legend.texts:
            legend_label.set_picker(True)
        
        return None

    def handle_map_z_statistic_update(self):
        """ Function which handles update of map z statistic upon interaction with map comboboxes. """

        if not self.read_instance.block_config_bar_handling_updates:

            # update map z statistic comboboxes
            # set variable that blocks configuration bar handling updates until all
            # changes to the z statistic comboboxes are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected items
            selected_z_stat = self.map_z_stat.currentText()
            selected_z1_array = self.map_z1.currentText()
            selected_z2_array = self.map_z2.currentText()

            # if selected_z_stat and selected_z1_array are empty strings it is
            # because they being initialised for the first time
            # force them to be 'observations' and first basic z statistic respectively
            if selected_z_stat == '':
                selected_z_stat = self.read_instance.basic_z_stats[0]
            if hasattr(self.read_instance, 'map_z'):
                if self.read_instance.map_z in self.read_instance.basic_z_stats:
                    selected_z_stat = self.read_instance.map_z
            if selected_z1_array == '':
                selected_z1_array = 'observations'

            # update z statistic field to all basic stats if colocation not-active OR z2
            # array not selected, else select basic+bias stats
            if (not self.read_instance.temporal_colocation) or (selected_z2_array == ''):
                z_stat_items = copy.deepcopy(self.read_instance.basic_z_stats)
            else:
                z_stat_items = copy.deepcopy(self.read_instance.basic_and_bias_z_stats)

            # remove selected z1/z2 items from opposite z2/z1 comboboxes (if have value
            # selected, i.e. z2 array not empty string)
            if selected_z2_array != '':
                z1_items = np.delete(self.read_instance.z1_arrays,
                                     np.where(self.read_instance.z1_arrays == selected_z2_array)[0])
            else:
                z1_items = self.read_instance.z1_arrays
            z2_items = np.delete(self.read_instance.z2_arrays,
                                 np.where(self.read_instance.z2_arrays == selected_z1_array)[0])

            # update all comboboxes (clear, then add items)
            self.map_z_stat.clear()
            self.map_z1.clear()
            self.map_z2.clear()
            self.map_z_stat.addItems(z_stat_items)
            self.map_z1.addItems(z1_items)
            self.map_z2.addItems(z2_items)

            # maintain currently selected z statistic (if exists in new item list)
            if selected_z_stat in z_stat_items:
                self.map_z_stat.setCurrentText(selected_z_stat)

            # maintain currently selected z1/z2 arrays
            self.map_z1.setCurrentText(selected_z1_array)
            self.map_z2.setCurrentText(selected_z2_array)

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted map z statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_map_z_statistic()

        return None

    def handle_periodic_statistic_update(self):
        """ Function that handles update of plotted periodic statistic
            upon interaction with periodic statistic combobox.
        """

        if not self.read_instance.block_config_bar_handling_updates:

            # update periodic statistic comboboxes

            # set variable that blocks configuration bar handling updates until all changes
            # to the periodic statistic combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected statistic
            zstat = self.periodic_stat.currentText()

            # update periodic statistics, to all basic stats
            # if colocation not-active, and basic+bias stats if colocation active
            if not self.read_instance.temporal_colocation:
                available_periodic_stats = copy.deepcopy(self.read_instance.basic_z_stats)
            else:
                available_periodic_stats = copy.deepcopy(self.read_instance.basic_and_bias_z_stats)

            # if base_zstat is empty string, it is because fields are being initialised for the first time
            if zstat == '':
                # set periodic stat to be first available stat
                zstat = available_periodic_stats[0]

            # update periodic statistic combobox (clear, then add items)
            self.periodic_stat.clear()
            self.periodic_stat.addItems(available_periodic_stats)

            # maintain currently selected periodic statistic (if exists in new item list)
            if zstat in available_periodic_stats:
                self.periodic_stat.setCurrentText(zstat)

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted periodic statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_periodic_statistic(zstat)

        return None

    def remove_axis_elements(self, ax, plot_type):
        """ Remove all plotted axis elements and hide axis. """
       
        # get appropriate axes for nested axes
        axs_to_remove = []
        if isinstance(ax, dict):
            for relevant_temporal_resolution, sub_ax in ax.items():
                axs_to_remove.append(sub_ax)
        else:
            axs_to_remove.append(ax)

        # iterate through axes
        for ax_to_remove in axs_to_remove:

            # remove all plotted axis elements
            if plot_type == 'legend':
                leg = ax_to_remove.get_legend()
                if leg:
                    leg.remove()

            elif plot_type == 'map':
                inds_to_remove = []
                for artist_ii, artist in enumerate(ax_to_remove.artists):         
                    if type(artist) == AnchoredOffsetbox:
                        inds_to_remove.append(artist_ii)
                ax_to_remove.artists = list(np.delete(np.array(ax_to_remove.artists), inds_to_remove))

                inds_to_remove = []
                for col_ii, col in enumerate(ax_to_remove.collections): 
                    # TODO: Put line collection back into place when we turn on the auto_update in gridlines
                    # if ((isinstance(col, matplotlib.collections.PathCollection))
                    #     or (isinstance(col, matplotlib.collections.LineCollection))):
                    if (isinstance(col, matplotlib.collections.PathCollection)):
                        inds_to_remove.append(col_ii)
                ax_to_remove.collections = list(np.delete(np.array(ax_to_remove.collections), inds_to_remove))

            elif plot_type == 'cb':
                ax_to_remove.artists = []
                ax_to_remove.collections = [] 

            elif plot_type == 'timeseries':
                ax_to_remove.lines = [self.timeseries_vline]
                ax_to_remove.artists = []
            
            elif plot_type == 'periodic':
                ax_to_remove.lines = []
                ax_to_remove.artists = []

            elif plot_type == 'periodic-violin':
                ax_to_remove.lines = []
                ax_to_remove.artists = []
                ax_to_remove.collections = []

            elif plot_type == 'metadata':
                ax_to_remove.texts = []

            elif plot_type == 'distribution':
                ax_to_remove.lines = [self.distribution_vline]
                ax_to_remove.artists = []

            elif plot_type == 'scatter':
                ax_to_remove.lines = []
                ax_to_remove.artists = []

            elif plot_type == 'statsummary':
                ax_to_remove.tables = []

            elif plot_type == 'boxplot':
                ax_to_remove.lines = []
                ax_to_remove.artists = []

            # hide axis
            ax_to_remove.axis('off')
            ax_to_remove.set_visible(False)

        # append vertical lines to periodic and periodic violin plots
        if plot_type == 'periodic':
            for resolution in self.plot_axes['periodic'].keys():
                ax_to_remove.lines.append(self.periodic_vline[resolution])
        
        elif plot_type == 'periodic-violin':
            for resolution in self.plot_axes['periodic-violin'].keys():
                ax_to_remove.lines.append(self.periodic_violin_vline[resolution])
                
        # hide plot buttons
        if plot_type == 'map':
            self.map_menu_button.hide()
            self.map_save_button.hide()

        elif plot_type == 'timeseries':
            self.timeseries_menu_button.hide()
            self.timeseries_save_button.hide()
        
        elif plot_type == 'periodic':
            self.periodic_menu_button.hide()
            self.periodic_save_button.hide()

        elif plot_type == 'periodic-violin':
            self.periodic_violin_menu_button.hide()
            self.periodic_violin_save_button.hide()

        elif plot_type == 'metadata':
            self.metadata_menu_button.hide()
            self.metadata_save_button.hide()

        elif plot_type == 'distribution':
            self.distribution_menu_button.hide()
            self.distribution_save_button.hide()

        elif plot_type == 'scatter':
            self.scatter_menu_button.hide()
            self.scatter_save_button.hide()

        elif plot_type == 'statsummary':
            self.statsummary_menu_button.hide()
            self.statsummary_save_button.hide()

        elif plot_type == 'boxplot':
            self.boxplot_menu_button.hide()
            self.boxplot_save_button.hide()

        # hide interactive elements
        if plot_type in self.interactive_elements:
            for element in self.interactive_elements[plot_type]['elements']:
                if isinstance(element, dict):
                    for sub_element in element.values():
                        sub_element.hide()
                else:
                    element.hide()

        # remove tracked plot elements
        if plot_type in self.plot_elements:
            self.plot_elements[plot_type]['absolute'] = {}
            if 'bias' in self.plot_elements[plot_type]:
                del self.plot_elements[plot_type]['bias'] 

        # hide layout options buttons when there is no selected data
        if hasattr(self, 'relative_selected_station_inds'):
            if len(self.relative_selected_station_inds) == 0:
                for position in range(2, 6):
                    cb_position = getattr(self.read_instance, 'cb_position_{}'.format(position))
                    cb_position.hide()

        return None

    def activate_axis(self, ax, plot_type):
        """ Un-hide axis. """

        # get appropriate axes for nested axes
        axs_to_activate = []
        if isinstance(ax, dict):
            for relevant_temporal_resolution, sub_ax in ax.items():
                if relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions:
                    axs_to_activate.append(sub_ax)
        else:
            axs_to_activate.append(ax)

        # iterate through axes
        for ax_to_activate in axs_to_activate:

            # un-hide axis
            ax_to_activate.axis('on')
            ax_to_activate.set_visible(True)

        # show plot buttons
        if plot_type == 'map':
            self.map_menu_button.show()
            self.map_save_button.show()

        elif plot_type == 'timeseries':
            self.timeseries_menu_button.show()
            self.timeseries_save_button.show()

        elif plot_type == 'periodic':
            self.periodic_menu_button.show()
            self.periodic_save_button.show()

        elif plot_type == 'periodic-violin':
            self.periodic_violin_menu_button.show()
            self.periodic_violin_save_button.show()

        elif plot_type == 'metadata':
            # hide for now because there are no options
            # self.metadata_menu_button.show()
            self.metadata_save_button.show()

        elif plot_type == 'distribution':
            self.distribution_menu_button.show()
            self.distribution_save_button.show()
    
        elif plot_type == 'scatter':
            self.scatter_menu_button.show()
            self.scatter_save_button.show()

        elif plot_type == 'statsummary':
            self.statsummary_menu_button.show()
            self.statsummary_save_button.show()

        elif plot_type == 'boxplot':
            self.boxplot_menu_button.show()
            self.boxplot_save_button.show()
        
        # show layout options buttons when there are selected stations
        if len(self.relative_selected_station_inds) > 0:
            for position in range(2, 6):
                cb_position = getattr(self.read_instance, 'cb_position_{}'.format(position))
                cb_position.show()

        return None

    def update_plot_options(self, plot_types):
        """ Uncheck checked boxes in plot configuration options under menus and
            reapply check with new data. This can be done for all currently active plot types, 
            or just one specific type.
        """

        for plot_type in plot_types:
            all_plot_options = self.plot_characteristics[plot_type]['plot_options']
            checked_options = self.read_instance.current_plot_options[plot_type]
            if plot_type == 'periodic-violin':
                plot_type = 'periodic_violin'
            cb_options = getattr(self, plot_type + '_options')
            for checked_option in checked_options:
                index = all_plot_options.index(checked_option)
                cb_options.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                cb_options.model().item(index).setCheckState(QtCore.Qt.Checked)

        return None

    def select_all_stations(self):
        """ Function that selects/unselects all plotted stations
            (and associated plots) upon ticking of checkbox.
        """

        if not self.read_instance.block_MPL_canvas_updates:

            # check if checkbox to select all stations is checked or unchecked
            check_state = self.read_instance.ch_select_all.checkState()

            # show warning and uncheck box
            if not hasattr(self, 'relative_selected_station_inds'):
                if check_state == QtCore.Qt.Checked:
                    msg = 'Data must be read into memory before selecting the data.'
                    show_message(self.read_instance, msg)
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False
                    return

            # make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

            # if checkbox is checked, select all plotted stations
            if check_state == QtCore.Qt.Checked:
                self.relative_selected_station_inds = copy.deepcopy(self.active_map_valid_station_inds)

                # if select intersect stations checkbox is checked then uncheck it (without updating canvas)
                if self.read_instance.ch_intersect.checkState() == QtCore.Qt.Checked:
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False

                # if select extent stations checkbox is checked then uncheck it (without updating canvas)
                if self.read_instance.ch_extent.checkState() == QtCore.Qt.Checked:
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_extent.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False

            # if checkbox is unchecked then unselect all plotted stations
            elif check_state == QtCore.Qt.Unchecked:
                self.relative_selected_station_inds = np.array([], dtype=np.int)

            # update absolute selected station indices (indices relative to plotted stations on map)
            self.absolute_selected_station_inds = np.arange(len(self.relative_selected_station_inds), dtype=np.int)

            # get absolute non-selected station inds
            self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                                 self.absolute_selected_station_inds))[0]

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
                self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw()

        return None

    def select_intersect_stations(self):
        """ Function that selects/unselects intersection of
            stations and all experiment domains (and associated plots)
            upon ticking of checkbox.
        """

        if not self.read_instance.block_MPL_canvas_updates:

            # check if checkbox to select intersection of stations is checked or unchecked
            check_state = self.read_instance.ch_intersect.checkState()

            # show warning and uncheck box
            if not hasattr(self, 'relative_selected_station_inds'):
                if check_state == QtCore.Qt.Checked:
                    msg = 'Data must be read into memory before selecting the data.'
                    show_message(self.read_instance, msg)
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False
                    return

            # make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

            # if checkbox is unchecked then unselect all plotted stations
            if check_state == QtCore.Qt.Unchecked:
                self.relative_selected_station_inds = np.array([], dtype=np.int)
                self.absolute_selected_station_inds = np.array([], dtype=np.int)
                self.absolute_non_selected_station_inds = np.arange(len(self.relative_selected_station_inds),
                                                                    dtype=np.int)

            # else, if checkbox is checked then select all stations which intersect with all loaded experiment domains
            elif check_state == QtCore.Qt.Checked:

                # if have only observations loaded into memory, select all plotted stations
                if len(self.read_instance.data_labels) == 1:
                    self.relative_selected_station_inds = copy.deepcopy(self.active_map_valid_station_inds)
                    self.absolute_selected_station_inds = np.arange(len(self.relative_selected_station_inds),
                                                                    dtype=np.int)
                    self.absolute_non_selected_station_inds = np.array([], dtype=np.int)
                # else, define list of lists to get intersection between (active_map_valid_station_inds, 
                # and valid station indices associated with each loaded experiment array)
                else:
                    intersect_lists = [self.active_map_valid_station_inds]
                    for data_label in self.read_instance.data_labels:
                        if data_label != 'observations':
                            if self.read_instance.temporal_colocation:
                                valid_station_inds = self.read_instance.valid_station_inds_temporal_colocation[self.read_instance.networkspeci][data_label]
                            else:
                                valid_station_inds = self.read_instance.valid_station_inds[self.read_instance.networkspeci][data_label]
                            intersect_lists.append(valid_station_inds)
                    
                    # get intersect between active map valid station indices and valid station indices
                    # associated with each loaded experiment array --> relative selected station indcies
                    self.relative_selected_station_inds = np.sort(list(set.intersection(*map(set, intersect_lists))))
                    
                    # get absolute selected station indices (indices relative to plotted stations on map)
                    self.absolute_selected_station_inds = \
                        np.array([np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind
                                  in self.relative_selected_station_inds], dtype=np.int)
                    
                    #get absolute non-selected station inds
                    self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                                         self.absolute_selected_station_inds))[0]

                # if select all stations checkbox is checked then uncheck it (without updating canvas)
                if self.read_instance.ch_select_all.checkState() == QtCore.Qt.Checked:
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False

                 # if select extent checkbox is checked then uncheck it (without updating canvas)
                if self.read_instance.ch_extent.checkState() == QtCore.Qt.Checked:
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_extent.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
                self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw()

        return None

    def select_extent_stations(self):
        """ Function that selects/unselects the
            stations for the current map extent (and associated plots)
            upon ticking of checkbox.
        """

        if not self.read_instance.block_MPL_canvas_updates:
            
            # check if checkbox to select extent of stations is checked or unchecked
            check_state = self.read_instance.ch_extent.checkState()

            # show warning and uncheck box
            if not hasattr(self, 'relative_selected_station_inds'):
                if check_state == QtCore.Qt.Checked:
                    msg = 'Data must be read into memory before selecting the data.'
                    show_message(self.read_instance, msg)
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_extent.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False
                    return
        
            # get map extent (in data coords)
            self.read_instance.map_extent = self.plot.get_map_extent(self.plot_axes['map'])

            # make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

            # if checkbox is checked, select all plotted stations
            if check_state == QtCore.Qt.Checked:

                # make copy of current full array relative selected stations indices, before selecting new ones
                self.relative_selected_station_inds = copy.deepcopy(self.active_map_valid_station_inds)

                # get inds of stations inside map extent
                extent_station_inds = []
                for ind, lon, lat in zip(self.relative_selected_station_inds, 
                                         self.read_instance.station_longitudes[self.read_instance.networkspeci][self.relative_selected_station_inds],
                                         self.read_instance.station_latitudes[self.read_instance.networkspeci][self.relative_selected_station_inds]):

                    if ((lon >= self.read_instance.map_extent[0]) and (lon <= self.read_instance.map_extent[1]) 
                        and (lat >= self.read_instance.map_extent[2]) and (lat <= self.read_instance.map_extent[3])):
                        extent_station_inds.append(ind)
                                           
                self.relative_selected_station_inds = extent_station_inds

                # if select intersect stations checkbox is checked then uncheck it (without updating canvas)
                if self.read_instance.ch_intersect.checkState() == QtCore.Qt.Checked:
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False

                # if select all stations checkbox is checked then uncheck it (without updating canvas)
                if self.read_instance.ch_select_all.checkState() == QtCore.Qt.Checked:
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False

            # if checkbox is unchecked then unselect all plotted stations
            elif check_state == QtCore.Qt.Unchecked:
                self.relative_selected_station_inds = np.array([], dtype=np.int)

            # get absolute selected station indices (indices relative to plotted stations on map)
            self.absolute_selected_station_inds = \
                np.array([np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind
                            in self.relative_selected_station_inds], dtype=np.int)

            # get absolute unselected station indices
            self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                                 self.absolute_selected_station_inds))[0]

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
                self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw()

        return None

    def onlassoselect_left(self, verts):
        """ Function that handles station selection upon lasso selection with left click.

            Operation:
            Select all stations within lasso boundaries.
            If a click is made rather than using lasso, then select nearest station within tolerance.

            If no station is found with left click, all stations are unselected.
        """

        # check if have any plotted stations on map, if not, return
        if len(self.active_map_valid_station_inds) == 0:
            return

        # unselect all/intersect/extent checkboxes
        self.unselect_map_checkboxes()

        # make copy of current full array absolute abd relative selected stations indices, before selecting new ones
        self.previous_absolute_selected_station_inds = copy.deepcopy(self.absolute_selected_station_inds)
        self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

        # get coordinates of drawn lasso
        lasso_path = Path(verts)
        lasso_path_vertices = lasso_path.vertices

        # transform lasso coordinates from projected to standard geographic coordinates
        lasso_path.vertices = self.datacrs.transform_points(self.plotcrs, lasso_path_vertices[:, 0], lasso_path_vertices[:, 1])[:, :2]

        # get absolute selected indices of stations on map (the station coordinates contained within lasso)
        self.absolute_selected_station_inds = np.nonzero(lasso_path.contains_points(self.map_points_coordinates))[0]

        # if have no valid selected indices, add a small tolerance (variable by visible map extent) to try get a match 
        if len(self.absolute_selected_station_inds) == 0:
            # take first selected point coordinates and get matches of stations within tolerance 
            self.read_instance.map_extent = self.plot.get_map_extent(self.plot_axes['map'])
            tolerance = np.average([self.read_instance.map_extent[1]-self.read_instance.map_extent[0],
                                    self.read_instance.map_extent[3]-self.read_instance.map_extent[2]]) / 100.0
            point_coordinates = lasso_path.vertices[0:1,:]
            sub_abs_vals = np.abs(self.map_points_coordinates[None,:,:] - point_coordinates[:,None,:])
            self.absolute_selected_station_inds = np.arange(len(self.active_map_valid_station_inds))[np.all(np.any(sub_abs_vals<=tolerance,axis=0),axis=1)]
            
            # if more than 1 point selected, limit this to be just nearest point
            if len(self.absolute_selected_station_inds) > 1:
                self.absolute_selected_station_inds = np.array([self.absolute_selected_station_inds[np.argmin(np.sum(sub_abs_vals[0,self.absolute_selected_station_inds,:],axis=1))]], dtype=np.int)

        # get absolute non-selected station inds
        self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                             self.absolute_selected_station_inds))[0]

        # get selected station indices with respect to all available stations
        self.relative_selected_station_inds = self.map_selected_station_inds_to_all_available_inds(self.absolute_selected_station_inds)

        # hide lasso after selection
        self.lasso_left.set_visible(False)

        # if selected stations have changed from previous selected, update station selection and associated plots
        if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
            self.update_map_station_selection()
            self.update_associated_active_dashboard_plots()

        # draw changes
        self.figure.canvas.draw()

        return None

    def onlassoselect_right(self, verts):
        """ Function that handles station selection upon lasso selection with right click.
        
            Operation:
            Unselect station/s (if station/s currently selected), 
            or Select station/s (if station/s currently unselected).

            If no station is found with right click, nothing happens.
        """

        # check if have any plotted stations on map, if not, return
        if len(self.active_map_valid_station_inds) == 0:
            return

        # unselect all/intersect/extent checkboxes
        self.unselect_map_checkboxes()

        # make copy of current full array relative selected stations indices, before selecting new ones
        previous_absolute_selected_station_inds = copy.deepcopy(self.absolute_selected_station_inds)
        previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

        # get coordinates of drawn lasso
        lasso_path = Path(verts)
        lasso_path_vertices = lasso_path.vertices
        # transform lasso coordinates from projected to standard geographic coordinates
        lasso_path.vertices = self.datacrs.transform_points(self.plotcrs, lasso_path_vertices[:, 0], lasso_path_vertices[:, 1])[:, :2]
        # get absolute selected indices of stations on map (the station coordinates contained within lasso)
        absolute_selected_station_inds = np.nonzero(lasso_path.contains_points(self.map_points_coordinates))[0]

        # if have no valid selected indices, add a small tolerance (variable by visible map extent) to try get a match 
        if len(absolute_selected_station_inds) == 0:
            # take first selected point coordinates and get matches of stations within tolerance
            self.read_instance.map_extent = self.plot.get_map_extent(self.plot_axes['map'])
            tolerance = np.average([self.read_instance.map_extent[1]-self.read_instance.map_extent[0],self.read_instance.map_extent[3]-self.read_instance.map_extent[2]]) / 100.0
            point_coordinates = lasso_path.vertices[0:1,:]
            sub_abs_vals = np.abs(self.map_points_coordinates[None,:,:] - point_coordinates[:,None,:])
            absolute_selected_station_inds = np.arange(len(self.active_map_valid_station_inds))[np.all(np.any(sub_abs_vals<=tolerance,axis=0),axis=1)]
            # if more than 1 point selected, limit this to be just nearest point
            if len(absolute_selected_station_inds) > 1:
                absolute_selected_station_inds = np.array([absolute_selected_station_inds[np.argmin(np.sum(sub_abs_vals[0,absolute_selected_station_inds,:],axis=1))]], dtype=np.int)

        # if have zero stations selected then return, doing nothing to selection
        if len(absolute_selected_station_inds) == 0:
            return 

        # update absolute indices of selected stations
        # remove stations that were previously selected, and add stations that were not previously selected
        self.absolute_selected_station_inds = np.setxor1d(previous_absolute_selected_station_inds, absolute_selected_station_inds)
        self.previous_absolute_selected_station_inds = previous_absolute_selected_station_inds

        # get absolute non-selected station inds
        self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                             self.absolute_selected_station_inds))[0]

        # get new relative selected indices with respect to all available stations
        self.previous_relative_selected_station_inds = previous_relative_selected_station_inds
        self.relative_selected_station_inds = self.map_selected_station_inds_to_all_available_inds(self.absolute_selected_station_inds)

        # hide lasso after selection
        self.lasso_right.set_visible(False)

        # if selected stations have changed from previous selected, update station selection and associated plots
        if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
            self.update_map_station_selection()
            self.update_associated_active_dashboard_plots()

        # draw changes
        self.figure.canvas.draw()

        return None

    def map_selected_station_inds_to_all_available_inds(self, selected_map_inds):
        """ Take the indices of selected stations on the map
            (potentially a subset of all available stations), and returns the indices
            of the stations inside the full loaded data arrays.
        """

        # index the array of indices of stations plotted on the map (indexed with respect to
        # all available stations), with the absolute indices of the subset of plotted selected stations
        return self.active_map_valid_station_inds[selected_map_inds]

    def generate_interactive_elements(self):
        """ Function to create settings menus for each plot and their elements."""

        self.interactive_elements = {}

        # LAYOUT OPTIONS #
        # add position 2 plot selector
        self.read_instance.cb_position_2 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.read_instance.cb_position_2.AdjustToContents
        self.read_instance.cb_position_2.setToolTip('Select plot type in top right position')
        self.read_instance.cb_position_2.currentTextChanged.connect(self.read_instance.handle_layout_update)
        self.read_instance.cb_position_2.hide()

        # add position 3 plot selector
        self.read_instance.cb_position_3 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.read_instance.cb_position_3.AdjustToContents
        self.read_instance.cb_position_3.setToolTip('Select plot type in bottom left position')
        self.read_instance.cb_position_3.currentTextChanged.connect(self.read_instance.handle_layout_update)
        self.read_instance.cb_position_3.hide()

        # add position 4 plot selector
        self.read_instance.cb_position_4 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.read_instance.cb_position_4.AdjustToContents
        self.read_instance.cb_position_4.setToolTip('Select plot type in bottom centre position')
        self.read_instance.cb_position_4.currentTextChanged.connect(self.read_instance.handle_layout_update)
        self.read_instance.cb_position_4.hide()

        # add position 5 plot selector
        self.read_instance.cb_position_5 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.read_instance.cb_position_5.AdjustToContents
        self.read_instance.cb_position_5.setToolTip('Select plot type in bottom right position')
        self.read_instance.cb_position_5.currentTextChanged.connect(self.read_instance.handle_layout_update)
        self.read_instance.cb_position_5.hide()

        # MAP SETTINGS MENU #
        # add button to map to show and hide settings menu
        self.map_menu_button = set_formatting(QtWidgets.QPushButton(self), 
                                              formatting_dict['settings_icon'])
        self.map_menu_button.setObjectName('map_menu_button')
        self.map_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.map_menu_button.setIconSize(QtCore.QSize(31, 37))
        self.map_menu_button.hide()

        # add white container
        self.map_container = set_formatting(QtWidgets.QWidget(self), 
                                            formatting_dict['settings_container'])
        self.map_container.setGeometry(self.map_menu_button.geometry().x()-230,
                                       self.map_menu_button.geometry().y()+25, 
                                       250, 430)
        self.map_container.hide()

        # add settings label
        self.map_settings_label = set_formatting(QtWidgets.QLabel('SETTINGS', self), 
                                                 formatting_dict['settings_label'])
        self.map_settings_label.setGeometry(self.map_menu_button.geometry().x()-220, 
                                            self.map_menu_button.geometry().y()+30, 
                                            230, 20)
        self.map_settings_label.hide()

        # add map stat label ('Statistic') to layout
        self.map_z_stat_label = QtWidgets.QLabel('Statistic', self)
        self.map_z_stat_label.setGeometry(self.map_menu_button.geometry().x()-220, 
                                          self.map_menu_button.geometry().y()+50, 
                                          230, 20)
        self.map_z_stat_label.hide()

        # add map stat combobox
        self.map_z_stat = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.map_z_stat.AdjustToContents
        self.map_z_stat.setGeometry(self.map_menu_button.geometry().x()-220, 
                                    self.map_menu_button.geometry().y()+75, 
                                    110, 20)
        self.map_z_stat.hide()

        # add map dataset 1 label ('Dataset 1') to layout
        self.map_z1_label = QtWidgets.QLabel('Dataset 1', self)
        self.map_z1_label.setGeometry(self.map_menu_button.geometry().x()-220, 
                                      self.map_menu_button.geometry().y()+100, 
                                      230, 20)
        self.map_z1_label.hide()

        # add map dataset 1 combobox
        self.map_z1 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.map_z1.AdjustToContents
        self.map_z1.setGeometry(self.map_menu_button.geometry().x()-220, 
                                self.map_menu_button.geometry().y()+125, 
                                110, 20)
        self.map_z1.hide()

        # add map dataset 2 label ('Dataset 2') to layout
        self.map_z2_label = QtWidgets.QLabel('Dataset 2', self)
        self.map_z2_label.setGeometry(self.map_menu_button.geometry().x()-100, 
                                      self.map_menu_button.geometry().y()+100, 
                                      230, 20)
        self.map_z2_label.hide()

        # add map dataset 2 combobox
        self.map_z2 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.map_z2.AdjustToContents
        self.map_z2.setGeometry(self.map_menu_button.geometry().x()-100, 
                                self.map_menu_button.geometry().y()+125, 
                                110, 20)
        self.map_z2.hide()

        # add map general text for unselected stations ('Unselected stations')
        self.map_unsel_label = QtWidgets.QLabel("Unselected stations", self)
        self.map_unsel_label.setGeometry(self.map_menu_button.geometry().x()-220,
                                         self.map_menu_button.geometry().y()+150, 
                                         230, 20)
        self.map_unsel_label.hide()

        # add map markersize slider name ('Size') to layout
        self.map_markersize_unsel_sl_label = QtWidgets.QLabel('Size', self)
        self.map_markersize_unsel_sl_label.setStyleSheet("QLabel { font-style: italic; }")
        self.map_markersize_unsel_sl_label.setGeometry(self.map_menu_button.geometry().x()-220,
                                                       self.map_menu_button.geometry().y()+175, 
                                                       230, 20)
        self.map_markersize_unsel_sl_label.hide()

        # add map markersize unselected stations slider
        self.map_markersize_unsel_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.map_markersize_unsel_sl.setObjectName('map_markersize_unsel_sl')
        self.map_markersize_unsel_sl.setMinimum(0)
        self.map_markersize_unsel_sl.setMaximum(80)
        self.map_markersize_unsel_sl.setValue(self.plot_characteristics['map']['marker_unselected']['s'])
        self.map_markersize_unsel_sl.setTickInterval(2)
        self.map_markersize_unsel_sl.setTracking(False)
        self.map_markersize_unsel_sl.setGeometry(self.map_menu_button.geometry().x()-220, 
                                                 self.map_menu_button.geometry().y()+200, 
                                                 230, 20)
        self.map_markersize_unsel_sl.hide()

        # add map opacity slider name ('Opacity') to layout
        self.map_opacity_unsel_sl_label = QtWidgets.QLabel('Opacity', self)
        self.map_opacity_unsel_sl_label.setStyleSheet("QLabel { font-style: italic; }")
        self.map_opacity_unsel_sl_label.setGeometry(self.map_menu_button.geometry().x()-220, 
                                                    self.map_menu_button.geometry().y()+225, 
                                                    230, 20)
        self.map_opacity_unsel_sl_label.hide()

        # add map opacity unselected stations slider
        self.map_opacity_unsel_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.map_opacity_unsel_sl.setObjectName('map_opacity_unsel_sl')
        self.map_opacity_unsel_sl.setMinimum(0)
        self.map_opacity_unsel_sl.setMaximum(10)
        self.map_opacity_unsel_sl.setValue(self.plot_characteristics['map']['marker_unselected']['alpha']*10)
        self.map_opacity_unsel_sl.setTickInterval(1)
        self.map_opacity_unsel_sl.setTracking(False)
        self.map_opacity_unsel_sl.setGeometry(self.map_menu_button.geometry().x()-220, 
                                              self.map_menu_button.geometry().y()+250, 
                                              230, 20)
        self.map_opacity_unsel_sl.hide()

        # add map general text for selected stations ('Selected stations')
        self.map_sel_label = QtWidgets.QLabel("Selected stations", self)
        self.map_sel_label.setGeometry(self.map_menu_button.geometry().x()-220, 
                                       self.map_menu_button.geometry().y()+275, 
                                       230, 20)
        self.map_sel_label.hide()

        # add map markersize slider name ('Size') to layout
        self.map_markersize_sel_sl_label = QtWidgets.QLabel('Size', self)
        self.map_markersize_sel_sl_label.setStyleSheet("QLabel { font-style: italic; }")
        self.map_markersize_sel_sl_label.setGeometry(self.map_menu_button.geometry().x()-220, 
                                                     self.map_menu_button.geometry().y()+300, 
                                                     230, 20)
        self.map_markersize_sel_sl_label.hide()

        # add map markersize selected stations slider
        self.map_markersize_sel_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.map_markersize_sel_sl.setObjectName('map_markersize_sel_sl')
        self.map_markersize_sel_sl.setMinimum(0)
        self.map_markersize_sel_sl.setMaximum(80)
        self.map_markersize_sel_sl.setValue(self.plot_characteristics['map']['marker_selected']['s'])
        self.map_markersize_sel_sl.setTickInterval(2)
        self.map_markersize_sel_sl.setTracking(False)
        self.map_markersize_sel_sl.setGeometry(self.map_menu_button.geometry().x()-220, 
                                               self.map_menu_button.geometry().y()+325, 
                                               230, 20)
        self.map_markersize_sel_sl.hide()

        # add map opacity slider name ('Opacity') to layout
        self.map_opacity_sel_sl_label = QtWidgets.QLabel('Opacity', self)
        self.map_opacity_sel_sl_label.setStyleSheet("QLabel { font-style: italic; }")
        self.map_opacity_sel_sl_label.setGeometry(self.map_menu_button.geometry().x()-220, 
                                                  self.map_menu_button.geometry().y()+350, 
                                                  230, 20)
        self.map_opacity_sel_sl_label.hide()

        # add map opacity selected stations slider
        self.map_opacity_sel_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.map_opacity_sel_sl.setObjectName('map_opacity_sel_sl')
        self.map_opacity_sel_sl.setMinimum(0)
        self.map_opacity_sel_sl.setMaximum(10)
        self.map_opacity_sel_sl.setValue(self.plot_characteristics['map']['marker_selected']['alpha']*10)
        self.map_opacity_sel_sl.setTickInterval(1)
        self.map_opacity_sel_sl.setTracking(False)
        self.map_opacity_sel_sl.setGeometry(self.map_menu_button.geometry().x()-220, 
                                            self.map_menu_button.geometry().y()+375, 
                                            230, 20)
        self.map_opacity_sel_sl.hide()

        # add map options name ('Options') to layout
        self.map_options_label = QtWidgets.QLabel("Options", self)
        self.map_options_label.setGeometry(self.map_menu_button.geometry().x()-220,
                                           self.map_menu_button.geometry().y()+400, 
                                           230, 20)
        self.map_options_label.hide()

        # add map options checkboxes
        self.map_options = CheckableComboBox(self)
        self.map_options.setObjectName('map_options')
        self.map_options.addItems(self.plot_characteristics['map']['plot_options'])        
        self.map_options.setGeometry(self.map_menu_button.geometry().x()-220, 
                                     self.map_menu_button.geometry().y()+425, 
                                     230, 20)
        self.map_options.currentTextChanged.connect(self.update_plot_option)
        self.map_options.hide()

        # add map figure save button
        self.map_save_button = set_formatting(QtWidgets.QPushButton(self), 
                                              formatting_dict['settings_icon'])
        self.map_save_button.setObjectName('map_save_button')
        self.map_save_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_fig_icon.png")))
        self.map_save_button.setIconSize(QtCore.QSize(20, 20))
        self.map_save_button.hide()

        # set show/hide actions
        self.map_elements = [self.map_container, self.map_settings_label, 
                             self.map_z_stat_label, self.map_z_stat, 
                             self.map_z1_label, self.map_z1, 
                             self.map_z2_label, self.map_z2, 
                             self.map_unsel_label, 
                             self.map_markersize_unsel_sl_label, self.map_markersize_unsel_sl, 
                             self.map_opacity_unsel_sl_label, self.map_opacity_unsel_sl, 
                             self.map_sel_label, self.map_markersize_sel_sl, 
                             self.map_markersize_sel_sl_label, self.map_opacity_sel_sl_label, 
                             self.map_opacity_sel_sl, self.map_options_label, 
                             self.map_options
                             ]
        self.interactive_elements['map'] = {'button': self.map_menu_button, 
                                            'hidden': True,
                                            'elements': self.map_elements,
                                            'markersize_sl': [self.map_markersize_unsel_sl, 
                                                              self.map_markersize_sel_sl],
                                            'opacity_sl': [self.map_opacity_unsel_sl, 
                                                           self.map_opacity_sel_sl],
                                            'linewidth_sl': []
                                            }
        self.map_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.map_markersize_unsel_sl.valueChanged.connect(self.update_markersize_func)
        self.map_markersize_sel_sl.valueChanged.connect(self.update_markersize_func)
        self.map_opacity_unsel_sl.valueChanged.connect(self.update_opacity_func)
        self.map_opacity_sel_sl.valueChanged.connect(self.update_opacity_func)
        self.map_save_button.clicked.connect(self.save_axis_figure_func)
        self.map_z_stat.currentTextChanged.connect(self.handle_map_z_statistic_update)
        self.map_z1.currentTextChanged.connect(self.handle_map_z_statistic_update)
        self.map_z2.currentTextChanged.connect(self.handle_map_z_statistic_update)

        # TIMESERIES PLOT SETTINGS MENU #
        # add button to timeseries to show and hide settings menu
        self.timeseries_menu_button = set_formatting(QtWidgets.QPushButton(self), 
                                                     formatting_dict['settings_icon'])
        self.timeseries_menu_button.setObjectName('timeseries_menu_button')
        self.timeseries_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.timeseries_menu_button.setIconSize(QtCore.QSize(31, 37))
        self.timeseries_menu_button.hide()

        # add white container
        self.timeseries_container = set_formatting(QtWidgets.QWidget(self), 
                                                   formatting_dict['settings_container'])
        self.timeseries_container.setGeometry(self.timeseries_menu_button.geometry().x()-230,
                                              self.timeseries_menu_button.geometry().y()+25, 
                                              250, 180)
        self.timeseries_container.hide()

        # add settings label
        self.timeseries_settings_label = set_formatting(QtWidgets.QLabel('SETTINGS', self), 
                                                        formatting_dict['settings_label'])
        self.timeseries_settings_label.setGeometry(self.timeseries_menu_button.geometry().x()-220, 
                                                   self.timeseries_menu_button.geometry().y()+30, 
                                                   230, 20)
        self.timeseries_settings_label.hide()

        # add timeseries markersize slider name ('Size') to layout
        self.timeseries_markersize_sl_label = QtWidgets.QLabel('Size', self)
        self.timeseries_markersize_sl_label.setGeometry(self.timeseries_menu_button.geometry().x()-220,
                                                        self.timeseries_menu_button.geometry().y()+50, 
                                                        230, 20)
        self.timeseries_markersize_sl_label.hide()

        # add timeseries markersize slider
        self.timeseries_markersize_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.timeseries_markersize_sl.setObjectName('timeseries_markersize_sl')
        self.timeseries_markersize_sl.setMinimum(0)
        self.timeseries_markersize_sl.setMaximum(self.plot_characteristics['timeseries']['plot']['markersize']*10)
        self.timeseries_markersize_sl.setValue(self.plot_characteristics['timeseries']['plot']['markersize'])
        self.timeseries_markersize_sl.setTickInterval(2)
        self.timeseries_markersize_sl.setTracking(False)
        self.timeseries_markersize_sl.setGeometry(self.timeseries_menu_button.geometry().x()-220, 
                                                  self.timeseries_menu_button.geometry().y()+75, 
                                                  230, 20)
        self.timeseries_markersize_sl.hide()

        # add timeseries smooth slider name ('Smooth') to layout
        self.timeseries_smooth_sl_label = QtWidgets.QLabel('Smooth', self)
        self.timeseries_smooth_sl_label.setGeometry(self.timeseries_menu_button.geometry().x()-220,
                                                    self.timeseries_menu_button.geometry().y()+100, 
                                                    230, 20)
        self.timeseries_smooth_sl_label.hide()

        # add timeseries smooth slider
        self.timeseries_smooth_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.timeseries_smooth_sl.setObjectName('timeseries_smooth_sl')
        self.timeseries_smooth_sl.setMinimum(0)
        self.timeseries_smooth_sl.setValue(0)
        self.timeseries_smooth_sl.setTickInterval(2)
        self.timeseries_smooth_sl.setTracking(False)
        self.timeseries_smooth_sl.setGeometry(self.timeseries_menu_button.geometry().x()-220, 
                                              self.timeseries_menu_button.geometry().y()+125, 
                                              230, 20)
        self.timeseries_smooth_sl.hide()

        # add timeseries plot options name ('Options') to layout
        self.timeseries_options_label = QtWidgets.QLabel("Options", self)
        self.timeseries_options_label.setGeometry(self.timeseries_menu_button.geometry().x()-220,
                                                  self.timeseries_menu_button.geometry().y()+150, 
                                                  230, 20)
        self.timeseries_options_label.hide()

        # add timeseries plot options checkboxes
        self.timeseries_options = CheckableComboBox(self)
        self.timeseries_options.setObjectName('timeseries_options')
        self.timeseries_options.addItems(self.plot_characteristics['timeseries']['plot_options'])        
        self.timeseries_options.setGeometry(self.timeseries_menu_button.geometry().x()-220, 
                                            self.timeseries_menu_button.geometry().y()+175, 
                                            230, 20)
        self.timeseries_options.currentTextChanged.connect(self.update_plot_option)
        self.timeseries_options.hide()

        # add timeseries figure save button
        self.timeseries_save_button = set_formatting(QtWidgets.QPushButton(self), 
                                                     formatting_dict['settings_icon'])
        self.timeseries_save_button.setObjectName('timeseries_save_button')
        self.timeseries_save_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_fig_icon.png")))
        self.timeseries_save_button.setIconSize(QtCore.QSize(20, 20))
        self.timeseries_save_button.hide()

        # set show/hide actions
        self.timeseries_elements = [self.timeseries_container, self.timeseries_settings_label, 
                                    self.timeseries_markersize_sl_label, self.timeseries_markersize_sl,
                                    self.timeseries_smooth_sl_label, self.timeseries_smooth_sl,
                                    self.timeseries_options_label, self.timeseries_options]
        self.interactive_elements['timeseries'] = {'button': self.timeseries_menu_button, 
                                                   'hidden': True,
                                                   'elements': self.timeseries_elements,
                                                   'markersize_sl': [self.timeseries_markersize_sl],
                                                   'opacity_sl': [],
                                                   'linewidth_sl': [],
                                                   'smooth_sl': [self.timeseries_smooth_sl]
                                                   }
        self.timeseries_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.timeseries_markersize_sl.valueChanged.connect(self.update_markersize_func)
        self.timeseries_smooth_sl.valueChanged.connect(self.update_smooth_func)
        self.timeseries_save_button.clicked.connect(self.save_axis_figure_func)

        # PERIODIC PLOT SETTINGS MENU #
        # add button to periodic plot to show and hide settings menu
        self.periodic_menu_button = set_formatting(QtWidgets.QPushButton(self), 
                                                   formatting_dict['settings_icon'])
        self.periodic_menu_button.setObjectName('periodic_menu_button')
        self.periodic_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.periodic_menu_button.setIconSize(QtCore.QSize(31, 37))
        self.periodic_menu_button.hide()

        # add white container
        self.periodic_container = set_formatting(QtWidgets.QWidget(self), 
                                                 formatting_dict['settings_container'])
        self.periodic_container.setGeometry(self.periodic_menu_button.geometry().x()-230, 
                                            self.periodic_menu_button.geometry().y()+25, 
                                            250, 230)
        self.periodic_container.hide()

        # add settings label
        self.periodic_settings_label = set_formatting(QtWidgets.QLabel('SETTINGS', self), 
                                                      formatting_dict['settings_label'])
        self.periodic_settings_label.setGeometry(self.periodic_menu_button.geometry().x()-220, 
                                                 self.periodic_menu_button.geometry().y()+30, 
                                                 230, 20)
        self.periodic_settings_label.hide()

        # add periodic stat label ('Statistic') to layout
        self.periodic_stat_label = QtWidgets.QLabel("Statistic", self)
        self.periodic_stat_label.setGeometry(self.periodic_menu_button.geometry().x()-220, 
                                             self.periodic_menu_button.geometry().y()+50, 
                                             230, 20)
        self.periodic_stat_label.hide()
        
        # add periodic stat combobox
        self.periodic_stat = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.periodic_stat.AdjustToContents
        self.periodic_stat.setGeometry(self.periodic_menu_button.geometry().x()-220, 
                                       self.periodic_menu_button.geometry().y()+75, 
                                       110, 20)
        self.periodic_stat.hide()

        # add periodic markersize slider name ('Size') to layout
        self.periodic_markersize_sl_label = QtWidgets.QLabel('Size', self)
        self.periodic_markersize_sl_label.setGeometry(self.periodic_menu_button.geometry().x()-220, 
                                                      self.periodic_menu_button.geometry().y()+100, 
                                                      230, 20)
        self.periodic_markersize_sl_label.hide()

        # add periodic markersize slider
        self.periodic_markersize_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.periodic_markersize_sl.setObjectName('periodic_markersize_sl')
        self.periodic_markersize_sl.setMinimum(0)
        self.periodic_markersize_sl.setMaximum(self.plot_characteristics['periodic']['plot']['markersize']*10)
        self.periodic_markersize_sl.setValue(self.plot_characteristics['periodic']['plot']['markersize'])
        self.periodic_markersize_sl.setTickInterval(2)
        self.periodic_markersize_sl.setTracking(False)
        self.periodic_markersize_sl.setGeometry(self.periodic_menu_button.geometry().x()-220, 
                                                self.periodic_menu_button.geometry().y()+125, 
                                                230, 20)
        self.periodic_markersize_sl.hide()

        # add periodic line width slider name ('Line width') to layout
        self.periodic_linewidth_sl_label = QtWidgets.QLabel("Line width", self)
        self.periodic_linewidth_sl_label.setGeometry(self.periodic_menu_button.geometry().x()-220, 
                                                     self.periodic_menu_button.geometry().y()+150, 
                                                     230, 20)
        self.periodic_linewidth_sl_label.hide()

        # add periodic line width slider
        self.periodic_linewidth_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.periodic_linewidth_sl.setObjectName('periodic_linewidth_sl')
        self.periodic_linewidth_sl.setMinimum(0)
        self.periodic_linewidth_sl.setMaximum(self.plot_characteristics['periodic']['plot']['linewidth']*100)
        self.periodic_linewidth_sl.setValue(self.plot_characteristics['periodic']['plot']['linewidth']*10)
        self.periodic_linewidth_sl.setTickInterval(2)
        self.periodic_linewidth_sl.setTracking(False)
        self.periodic_linewidth_sl.setGeometry(self.periodic_menu_button.geometry().x()-220, 
                                               self.periodic_menu_button.geometry().y()+175, 
                                               230, 20)
        self.periodic_linewidth_sl.hide()

        # add periodic plot options name ('Options') to layout
        self.periodic_options_label = QtWidgets.QLabel("Options", self)
        self.periodic_options_label.setGeometry(self.periodic_menu_button.geometry().x()-220,
                                                self.periodic_menu_button.geometry().y()+200, 
                                                230, 20)
        self.periodic_options_label.hide()

        # add periodic plot options checkboxes
        self.periodic_options = CheckableComboBox(self)
        self.periodic_options.setObjectName('periodic_options')
        self.periodic_options.addItems(self.plot_characteristics['periodic']['plot_options'])        
        self.periodic_options.setGeometry(self.periodic_menu_button.geometry().x()-220, 
                                          self.periodic_menu_button.geometry().y()+225, 
                                          230, 20)
        self.periodic_options.currentTextChanged.connect(self.update_plot_option)
        self.periodic_options.hide()

        # add periodic figure save button
        self.periodic_save_button = set_formatting(QtWidgets.QPushButton(self), 
                                                   formatting_dict['settings_icon'])
        self.periodic_save_button.setObjectName('periodic_save_button')
        self.periodic_save_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_fig_icon.png")))
        self.periodic_save_button.setIconSize(QtCore.QSize(20, 20))
        self.periodic_save_button.hide()

        # set show/hide actions
        self.periodic_elements = [self.periodic_container, self.periodic_settings_label, 
                                  self.periodic_stat_label, self.periodic_stat, 
                                  self.periodic_markersize_sl_label, self.periodic_markersize_sl,
                                  self.periodic_linewidth_sl_label, self.periodic_linewidth_sl,
                                  self.periodic_options_label, self.periodic_options
                                 ]
        self.interactive_elements['periodic'] = {'button': self.periodic_menu_button, 
                                                 'hidden': True,
                                                 'elements': self.periodic_elements,
                                                 'markersize_sl': [self.periodic_markersize_sl],
                                                 'opacity_sl': [],
                                                 'linewidth_sl': [self.periodic_linewidth_sl]
                                                 }
        self.periodic_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.periodic_markersize_sl.valueChanged.connect(self.update_markersize_func)
        self.periodic_linewidth_sl.valueChanged.connect(self.update_linewidth_func)
        self.periodic_stat.currentTextChanged.connect(self.handle_periodic_statistic_update)
        self.periodic_save_button.clicked.connect(self.save_axis_figure_func)

        # PERIODIC VIOLIN PLOT SETTINGS MENU #
        # add button to periodic violin plot to show and hide settings menu
        self.periodic_violin_menu_button = set_formatting(QtWidgets.QPushButton(self), 
                                                          formatting_dict['settings_icon'])
        self.periodic_violin_menu_button.setObjectName('periodic_violin_menu_button')
        self.periodic_violin_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.periodic_violin_menu_button.setIconSize(QtCore.QSize(31, 37))
        self.periodic_violin_menu_button.hide()

        # add white container
        self.periodic_violin_container = set_formatting(QtWidgets.QWidget(self), 
                                                        formatting_dict['settings_container'])
        self.periodic_violin_container.setGeometry(self.periodic_violin_menu_button.geometry().x()-230, 
                                                   self.periodic_violin_menu_button.geometry().y()+25, 
                                                   250, 180)
        self.periodic_violin_container.hide()

        # add settings label
        self.periodic_violin_settings_label = set_formatting(QtWidgets.QLabel('SETTINGS', self), 
                                                             formatting_dict['settings_label'])
        self.periodic_violin_settings_label.setGeometry(self.periodic_violin_menu_button.geometry().x()-220,
                                                        self.periodic_violin_menu_button.geometry().y()+30, 
                                                        230, 20)
        self.periodic_violin_settings_label.hide()

        # add periodic violin markersize slider name ('Size') to layout
        self.periodic_violin_markersize_sl_label = QtWidgets.QLabel('Size', self)
        self.periodic_violin_markersize_sl_label.setGeometry(self.periodic_violin_menu_button.geometry().x()-220,
                                                             self.periodic_violin_menu_button.geometry().y()+50, 
                                                             230, 20)
        self.periodic_violin_markersize_sl_label.hide()

        # add periodic violin markersize slider
        self.periodic_violin_markersize_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.periodic_violin_markersize_sl.setObjectName('periodic_violin_markersize_sl')
        self.periodic_violin_markersize_sl.setMinimum(0)
        self.periodic_violin_markersize_sl.setMaximum(self.plot_characteristics['periodic-violin']['plot']['p50']['markersize']*10)
        self.periodic_violin_markersize_sl.setValue(self.plot_characteristics['periodic-violin']['plot']['p50']['markersize'])
        self.periodic_violin_markersize_sl.setTickInterval(2)
        self.periodic_violin_markersize_sl.setTracking(False)
        self.periodic_violin_markersize_sl.setGeometry(self.periodic_violin_menu_button.geometry().x()-220,
                                                       self.periodic_menu_button.geometry().y()+75, 
                                                       230, 20)
        self.periodic_violin_markersize_sl.hide()

        # add periodic violin line width slider name ('Line width') to layout
        self.periodic_violin_linewidth_sl_label = QtWidgets.QLabel("Line width", self)
        self.periodic_violin_linewidth_sl_label.setGeometry(self.periodic_violin_menu_button.geometry().x()-220, 
                                                            self.periodic_violin_menu_button.geometry().y()+100, 
                                                            230, 20)
        self.periodic_violin_linewidth_sl_label.hide()

        # add periodic violin line width slider
        self.periodic_violin_linewidth_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.periodic_violin_linewidth_sl.setObjectName('periodic_violin_linewidth_sl')
        self.periodic_violin_linewidth_sl.setMinimum(0)
        self.periodic_violin_linewidth_sl.setMaximum(self.plot_characteristics['periodic-violin']['plot']['p50']['linewidth']*100)
        self.periodic_violin_linewidth_sl.setValue(self.plot_characteristics['periodic-violin']['plot']['p50']['linewidth']*10)
        self.periodic_violin_linewidth_sl.setTickInterval(2)
        self.periodic_violin_linewidth_sl.setTracking(False)
        self.periodic_violin_linewidth_sl.setGeometry(self.periodic_violin_menu_button.geometry().x()-220, 
                                                      self.periodic_violin_menu_button.geometry().y()+125, 
                                                      230, 20)
        self.periodic_violin_linewidth_sl.hide()

        # add periodic violin plot options name ('Options') to layout
        self.periodic_violin_options_label = QtWidgets.QLabel("Options", self)
        self.periodic_violin_options_label.setGeometry(self.periodic_violin_menu_button.geometry().x()-220,
                                                       self.periodic_violin_menu_button.geometry().y()+150, 
                                                       230, 20)
        self.periodic_violin_options_label.hide()

        # add periodic violin plot options checkboxes
        self.periodic_violin_options = CheckableComboBox(self)
        self.periodic_violin_options.setObjectName('periodic_violin_options')
        self.periodic_violin_options.addItems(self.plot_characteristics['periodic-violin']['plot_options'])        
        self.periodic_violin_options.setGeometry(self.periodic_violin_menu_button.geometry().x()-220, 
                                                 self.periodic_violin_menu_button.geometry().y()+175, 
                                                 230, 20)
        self.periodic_violin_options.currentTextChanged.connect(self.update_plot_option)
        self.periodic_violin_options.hide()

        # add periodic violin figure save button
        self.periodic_violin_save_button = set_formatting(QtWidgets.QPushButton(self), 
                                                          formatting_dict['settings_icon'])
        self.periodic_violin_save_button.setObjectName('periodic_violin_save_button')
        self.periodic_violin_save_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_fig_icon.png")))
        self.periodic_violin_save_button.setIconSize(QtCore.QSize(20, 20))
        self.periodic_violin_save_button.hide()

        # set show/hide actions
        self.periodic_violin_elements = [self.periodic_violin_container, self.periodic_violin_settings_label, 
                                         self.periodic_violin_markersize_sl_label, self.periodic_violin_markersize_sl,
                                         self.periodic_violin_linewidth_sl_label, self.periodic_violin_linewidth_sl,
                                         self.periodic_violin_options_label, self.periodic_violin_options
                                        ]
        self.interactive_elements['periodic-violin'] = {'button': self.periodic_violin_menu_button, 
                                                        'hidden': True,
                                                        'elements': self.periodic_violin_elements,
                                                        'markersize_sl': [self.periodic_violin_markersize_sl],
                                                        'opacity_sl': [],
                                                        'linewidth_sl': [self.periodic_violin_linewidth_sl]
                                                        }
        self.periodic_violin_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.periodic_violin_markersize_sl.valueChanged.connect(self.update_markersize_func)
        self.periodic_violin_linewidth_sl.valueChanged.connect(self.update_linewidth_func)
        self.periodic_violin_save_button.clicked.connect(self.save_axis_figure_func)

        # METADATA PLOT SETTINGS MENU #
        # add button to metadata plot to show and hide settings menu
        self.metadata_menu_button = set_formatting(QtWidgets.QPushButton(self), 
                                                   formatting_dict['settings_icon'])
        self.metadata_menu_button.setObjectName('metadata_menu_button')
        self.metadata_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.metadata_menu_button.setIconSize(QtCore.QSize(31, 37))
        self.metadata_menu_button.hide()

        # add white container
        self.metadata_container = set_formatting(QtWidgets.QWidget(self), 
                                                 formatting_dict['settings_container'])
        self.metadata_container.setGeometry(self.metadata_menu_button.geometry().x()-230,
                                            self.metadata_menu_button.geometry().y()+25, 
                                            250, 75)
        self.metadata_container.hide()

        # add settings label
        self.metadata_settings_label = set_formatting(QtWidgets.QLabel('SETTINGS', self), 
                                                      formatting_dict['settings_label'])
        self.metadata_settings_label.setGeometry(self.metadata_menu_button.geometry().x()-220, 
                                                 self.metadata_menu_button.geometry().y()+30, 
                                                 230, 20)
        self.metadata_settings_label.hide()

        # add metadata plot options name ('Options') to layout
        self.metadata_options_label = QtWidgets.QLabel("Options", self)
        self.metadata_options_label.setGeometry(self.metadata_menu_button.geometry().x()-220,
                                                self.metadata_menu_button.geometry().y()+50, 
                                                230, 20)
        self.metadata_options_label.hide()

        # add metadata plot options checkboxes
        self.metadata_options = CheckableComboBox(self)
        self.metadata_options.setObjectName('metadata_options')
        self.metadata_options.addItems(self.plot_characteristics['metadata']['plot_options'])        
        self.metadata_options.setGeometry(self.metadata_menu_button.geometry().x()-220, 
                                          self.metadata_menu_button.geometry().y()+75, 
                                          230, 20)
        self.metadata_options.currentTextChanged.connect(self.update_plot_option)
        self.metadata_options.hide()

        # add metadata figure save button
        self.metadata_save_button = set_formatting(QtWidgets.QPushButton(self), 
                                                   formatting_dict['settings_icon'])
        self.metadata_save_button.setObjectName('metadata_save_button')
        self.metadata_save_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_fig_icon.png")))
        self.metadata_save_button.setIconSize(QtCore.QSize(20, 20))
        self.metadata_save_button.hide()

        # set show/hide actions
        self.metadata_elements = [self.metadata_container, self.metadata_settings_label, 
                                  self.metadata_options_label, self.metadata_options]
        self.interactive_elements['metadata'] = {'button': self.metadata_menu_button, 
                                                 'hidden': True,
                                                 'elements': self.metadata_elements,
                                                 'markersize_sl': [],
                                                 'opacity_sl': [],
                                                 'linewidth_sl': []
                                                 }
        self.metadata_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.metadata_save_button.clicked.connect(self.save_axis_figure_func)

        # DISTRIBUTION PLOT SETTINGS MENU #
        # add button to distribution plot to show and hide settings menu
        self.distribution_menu_button = set_formatting(QtWidgets.QPushButton(self), 
                                                       formatting_dict['settings_icon'])
        self.distribution_menu_button.setObjectName('distribution_menu_button')
        self.distribution_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.distribution_menu_button.setIconSize(QtCore.QSize(31, 37))
        self.distribution_menu_button.hide()

        # add white container
        self.distribution_container = set_formatting(QtWidgets.QWidget(self), 
                                                     formatting_dict['settings_container'])
        self.distribution_container.setGeometry(self.distribution_menu_button.geometry().x()-230,
                                                self.distribution_menu_button.geometry().y()+25, 
                                                250, 130)
        self.distribution_container.hide()

        # add settings label
        self.distribution_settings_label = set_formatting(QtWidgets.QLabel('SETTINGS', self), 
                                                          formatting_dict['settings_label'])
        self.distribution_settings_label.setGeometry(self.distribution_menu_button.geometry().x()-220, 
                                                     self.distribution_menu_button.geometry().y()+30, 
                                                     230, 20)
        self.distribution_settings_label.hide()

        # add distribution plot line width slider name ('Line width') to layout
        self.distribution_linewidth_sl_label = QtWidgets.QLabel("Line width", self)
        self.distribution_linewidth_sl_label.setGeometry(self.distribution_menu_button.geometry().x()-220,
                                                         self.distribution_menu_button.geometry().y()+50, 
                                                         230, 20)
        self.distribution_linewidth_sl_label.hide()

        # add distribution plot line width slider
        self.distribution_linewidth_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.distribution_linewidth_sl.setObjectName('distribution_linewidth_sl')
        self.distribution_linewidth_sl.setMinimum(0)
        self.distribution_linewidth_sl.setMaximum(self.plot_characteristics['distribution']['plot']['linewidth']*100)
        self.distribution_linewidth_sl.setValue(self.plot_characteristics['distribution']['plot']['linewidth']*10)
        self.distribution_linewidth_sl.setTickInterval(2)
        self.distribution_linewidth_sl.setTracking(False)
        self.distribution_linewidth_sl.setGeometry(self.distribution_menu_button.geometry().x()-220, 
                                                   self.distribution_menu_button.geometry().y()+75, 
                                                   230, 20)
        self.distribution_linewidth_sl.hide()

        # add distribution plot options name ('Options') to layout
        self.distribution_options_label = QtWidgets.QLabel("Options", self)
        self.distribution_options_label.setGeometry(self.distribution_menu_button.geometry().x()-220,
                                                    self.distribution_menu_button.geometry().y()+100, 
                                                    230, 20)
        self.distribution_options_label.hide()

        # add distribution plot options checkboxes
        self.distribution_options = CheckableComboBox(self)
        self.distribution_options.setObjectName('distribution_options')
        self.distribution_options.addItems(self.plot_characteristics['distribution']['plot_options'])        
        self.distribution_options.setGeometry(self.distribution_menu_button.geometry().x()-220, 
                                              self.distribution_menu_button.geometry().y()+125, 
                                              230, 20)
        self.distribution_options.currentTextChanged.connect(self.update_plot_option)
        self.distribution_options.hide()

        # add distribution figure save button
        self.distribution_save_button = set_formatting(QtWidgets.QPushButton(self), 
                                                       formatting_dict['settings_icon'])
        self.distribution_save_button.setObjectName('distribution_save_button')
        self.distribution_save_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_fig_icon.png")))
        self.distribution_save_button.setIconSize(QtCore.QSize(20, 20))
        self.distribution_save_button.hide()

        # set show/hide actions
        self.distribution_elements = [self.distribution_container, self.distribution_settings_label, 
                                      self.distribution_linewidth_sl_label, self.distribution_linewidth_sl,
                                      self.distribution_options_label, self.distribution_options]
        self.interactive_elements['distribution'] = {'button': self.distribution_menu_button, 
                                                     'hidden': True,
                                                     'elements': self.distribution_elements,
                                                     'markersize_sl': [],
                                                     'opacity_sl': [],
                                                     'linewidth_sl': [self.distribution_linewidth_sl]
                                                    }
        self.distribution_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.distribution_linewidth_sl.valueChanged.connect(self.update_linewidth_func)
        self.distribution_save_button.clicked.connect(self.save_axis_figure_func)

        # SCATTER PLOT SETTINGS MENU #
        # add button to scatter plot to show and hide settings menu
        self.scatter_menu_button = set_formatting(QtWidgets.QPushButton(self), 
                                                  formatting_dict['settings_icon'])
        self.scatter_menu_button.setObjectName('scatter_menu_button')
        self.scatter_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.scatter_menu_button.setIconSize(QtCore.QSize(31, 37))
        self.scatter_menu_button.hide()

        # add white container
        self.scatter_container = set_formatting(QtWidgets.QWidget(self),
                                                formatting_dict['settings_container'])
        self.scatter_container.setGeometry(self.scatter_menu_button.geometry().x()-230,
                                           self.scatter_menu_button.geometry().y()+25, 
                                           250, 130)
        self.scatter_container.hide()

        # add settings label
        self.scatter_settings_label = set_formatting(QtWidgets.QLabel('SETTINGS', self), 
                                                     formatting_dict['settings_label'])
        self.scatter_settings_label.setGeometry(self.scatter_menu_button.geometry().x()-220, 
                                                self.scatter_menu_button.geometry().y()+30, 
                                                230, 20)
        self.scatter_settings_label.hide()

        # add scatter plot markersize slider name ('Size') to layout
        self.scatter_markersize_sl_label = QtWidgets.QLabel('Size', self)
        self.scatter_markersize_sl_label.setGeometry(self.scatter_menu_button.geometry().x()-220,
                                                    self.scatter_menu_button.geometry().y()+50, 
                                                    230, 20)
        self.scatter_markersize_sl_label.hide()

        # add scatter plot markersize slider
        self.scatter_markersize_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.scatter_markersize_sl.setObjectName('scatter_markersize_sl')
        self.scatter_markersize_sl.setMinimum(0)
        self.scatter_markersize_sl.setMaximum(self.plot_characteristics['scatter']['plot']['markersize']*100)
        self.scatter_markersize_sl.setValue(self.plot_characteristics['scatter']['plot']['markersize']*10)
        self.scatter_markersize_sl.setTickInterval(2)
        self.scatter_markersize_sl.setTracking(False)
        self.scatter_markersize_sl.setGeometry(self.scatter_menu_button.geometry().x()-220, 
                                              self.scatter_menu_button.geometry().y()+75, 
                                              230, 20)
        self.scatter_markersize_sl.hide()

        # add scatter plot options name ('Options') to layout
        self.scatter_options_label = QtWidgets.QLabel("Options", self)
        self.scatter_options_label.setGeometry(self.scatter_menu_button.geometry().x()-220,
                                               self.scatter_menu_button.geometry().y()+100, 
                                               230, 20)
        self.scatter_options_label.hide()

        # add scatter plot options checkboxes
        self.scatter_options = CheckableComboBox(self)
        self.scatter_options.setObjectName('scatter_options')
        self.scatter_options.addItems(self.plot_characteristics['scatter']['plot_options'])        
        self.scatter_options.setGeometry(self.scatter_menu_button.geometry().x()-220, 
                                         self.scatter_menu_button.geometry().y()+125, 
                                         230, 20)
        self.scatter_options.currentTextChanged.connect(self.update_plot_option)
        self.scatter_options.hide()

        # add scatter figure save button
        self.scatter_save_button = set_formatting(QtWidgets.QPushButton(self), formatting_dict['settings_icon'])
        self.scatter_save_button.setObjectName('scatter_save_button')
        self.scatter_save_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_fig_icon.png")))
        self.scatter_save_button.setIconSize(QtCore.QSize(20, 20))
        self.scatter_save_button.hide()

        # set show/hide actions
        self.scatter_elements = [self.scatter_container, self.scatter_settings_label, 
                                 self.scatter_markersize_sl_label, self.scatter_markersize_sl,
                                 self.scatter_options_label, self.scatter_options]
        self.interactive_elements['scatter'] = {'button': self.scatter_menu_button, 
                                                'hidden': True,
                                                'elements': self.scatter_elements,
                                                'markersize_sl': [self.scatter_markersize_sl],
                                                'opacity_sl': [],
                                                'linewidth_sl': []
                                               }
        self.scatter_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.scatter_markersize_sl.valueChanged.connect(self.update_markersize_func)
        self.scatter_save_button.clicked.connect(self.save_axis_figure_func)

        # STATSUMMARY PLOT SETTINGS MENU #
        # add button to statsummary plot to show and hide settings menu
        self.statsummary_menu_button = set_formatting(QtWidgets.QPushButton(self), formatting_dict['settings_icon'])
        self.statsummary_menu_button.setObjectName('statsummary_menu_button')
        self.statsummary_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.statsummary_menu_button.setIconSize(QtCore.QSize(31, 37))
        self.statsummary_menu_button.hide()

        # add white container
        self.statsummary_container = set_formatting(QtWidgets.QWidget(self), formatting_dict['settings_container'])
        self.statsummary_container.setGeometry(self.statsummary_menu_button.geometry().x()-230,
                                               self.statsummary_menu_button.geometry().y()+25, 
                                               250, 80)
        self.statsummary_container.hide()

        # add settings label
        self.statsummary_settings_label = set_formatting(QtWidgets.QLabel('SETTINGS', self), 
                                                                          formatting_dict['settings_label'])
        self.statsummary_settings_label.setGeometry(self.statsummary_menu_button.geometry().x()-220, 
                                                    self.statsummary_menu_button.geometry().y()+30, 
                                                    230, 20)
        self.statsummary_settings_label.hide()

        # add statsummary plot options name ('Options') to layout
        self.statsummary_options_label = QtWidgets.QLabel("Options", self)
        self.statsummary_options_label.setGeometry(self.statsummary_menu_button.geometry().x()-220,
                                                   self.statsummary_menu_button.geometry().y()+50, 
                                                   230, 20)
        self.statsummary_options_label.hide()

        # add statsummary plot options checkboxes
        self.statsummary_options = CheckableComboBox(self)
        self.statsummary_options.setObjectName('statsummary_options')
        self.statsummary_options.addItems(self.plot_characteristics['statsummary']['plot_options'])        
        self.statsummary_options.setGeometry(self.statsummary_menu_button.geometry().x()-220, 
                                             self.statsummary_menu_button.geometry().y()+75, 
                                             230, 20)
        self.statsummary_options.currentTextChanged.connect(self.update_plot_option)
        self.statsummary_options.hide()

        # add statsummary figure save button
        self.statsummary_save_button = set_formatting(QtWidgets.QPushButton(self), formatting_dict['settings_icon'])
        self.statsummary_save_button.setObjectName('statsummary_save_button')
        self.statsummary_save_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_fig_icon.png")))
        self.statsummary_save_button.setIconSize(QtCore.QSize(20, 20))
        self.statsummary_save_button.hide()

        # set show/hide actions
        self.statsummary_elements = [self.statsummary_container, self.statsummary_settings_label, 
                                     self.statsummary_options_label, self.statsummary_options]
        self.interactive_elements['statsummary'] = {'button': self.statsummary_menu_button, 
                                                    'hidden': True,
                                                    'elements': self.statsummary_elements,
                                                    'markersize_sl': [],
                                                    'opacity_sl': [],
                                                    'linewidth_sl': []
                                                    }
        self.statsummary_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.statsummary_save_button.clicked.connect(self.save_axis_figure_func)

        # BOXPLOT PLOT SETTINGS MENU #
        # add button to boxplot to show and hide settings menu
        self.boxplot_menu_button = set_formatting(QtWidgets.QPushButton(self), formatting_dict['settings_icon'])
        self.boxplot_menu_button.setObjectName('boxplot_menu_button')
        self.boxplot_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.boxplot_menu_button.setIconSize(QtCore.QSize(31, 37))
        self.boxplot_menu_button.hide()

        # add white container
        self.boxplot_container = set_formatting(QtWidgets.QWidget(self), formatting_dict['settings_container'])
        self.boxplot_container.setGeometry(self.boxplot_menu_button.geometry().x()-230,
                                           self.boxplot_menu_button.geometry().y()+25, 
                                           250, 80)
        self.boxplot_container.hide()

        # add settings label
        self.boxplot_settings_label = set_formatting(QtWidgets.QLabel('SETTINGS', self), 
                                                     formatting_dict['settings_label'])
        self.boxplot_settings_label.setGeometry(self.boxplot_menu_button.geometry().x()-220, 
                                                self.boxplot_menu_button.geometry().y()+30, 
                                                230, 20)
        self.boxplot_settings_label.hide()

        # add boxplot options name ('Options') to layout
        self.boxplot_options_label = QtWidgets.QLabel("Options", self)
        self.boxplot_options_label.setGeometry(self.boxplot_menu_button.geometry().x()-220,
                                               self.boxplot_menu_button.geometry().y()+50, 
                                               230, 20)
        self.boxplot_options_label.hide()

        # add boxplot options checkboxes
        self.boxplot_options = CheckableComboBox(self)
        self.boxplot_options.setObjectName('boxplot_options')
        self.boxplot_options.addItems(self.plot_characteristics['boxplot']['plot_options'])        
        self.boxplot_options.setGeometry(self.boxplot_menu_button.geometry().x()-220, 
                                         self.boxplot_menu_button.geometry().y()+75, 
                                         230, 20)
        self.boxplot_options.currentTextChanged.connect(self.update_plot_option)
        self.boxplot_options.hide()

        # add boxplot figure save button
        self.boxplot_save_button = set_formatting(QtWidgets.QPushButton(self), formatting_dict['settings_icon'])
        self.boxplot_save_button.setObjectName('boxplot_save_button')
        self.boxplot_save_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/save_fig_icon.png")))
        self.boxplot_save_button.setIconSize(QtCore.QSize(20, 20))
        self.boxplot_save_button.hide()

        # set show/hide actions
        self.boxplot_elements = [self.boxplot_container, self.boxplot_settings_label, 
                                 self.boxplot_options_label, self.boxplot_options]
        self.interactive_elements['boxplot'] = {'button': self.boxplot_menu_button, 
                                                'hidden': True,
                                                'elements': self.boxplot_elements,
                                                'markersize_sl': [],
                                                'opacity_sl': [],
                                                'linewidth_sl': []
                                               }
        self.boxplot_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.boxplot_save_button.clicked.connect(self.save_axis_figure_func)

        # create array with buttons and elements to edit when the canvas is resized or the plots are changed
        self.menu_buttons = [self.map_menu_button, self.timeseries_menu_button,
                             self.periodic_menu_button, self.periodic_violin_menu_button,
                             self.metadata_menu_button, self.distribution_menu_button, 
                             self.scatter_menu_button, self.statsummary_menu_button, 
                             self.boxplot_menu_button]

        self.save_buttons = [self.map_save_button, self.timeseries_save_button,
                             self.periodic_save_button, self.periodic_violin_save_button,
                             self.metadata_save_button, self.distribution_save_button, 
                             self.scatter_save_button, self.statsummary_save_button,
                             self.boxplot_save_button]

        self.elements = [self.map_elements, self.timeseries_elements, 
                         self.periodic_elements, self.periodic_violin_elements,
                         self.metadata_elements, self.distribution_elements, 
                         self.scatter_elements, self.statsummary_elements, 
                         self.boxplot_elements]

        # make sure white containers are above buttons
        for element in self.elements:
            for sub_element in element:
                if isinstance(sub_element, dict):
                    for val in sub_element.values():
                        val.raise_()
                else:
                    sub_element.raise_()

        return None

    def interactive_elements_button_func(self):
        """ Function to show and hide elements in setting menus. """

        event_source = self.sender()
        for key, val in self.interactive_elements.items():
            if event_source == self.interactive_elements[key]['button']:
                hidden = self.interactive_elements[key]['hidden']
                elements = self.interactive_elements[key]['elements']
                break

        if hidden:
            for element in elements: 
                if isinstance(element, dict):
                    for sub_element in element.values():
                        sub_element.show()
                else:
                    element.show()

            self.interactive_elements[key]['hidden'] = False
        else:
            for element in elements: 
                if isinstance(element, dict):
                    for sub_element in element.values():
                        sub_element.hide()
                else:
                    element.hide()
            self.interactive_elements[key]['hidden'] = True

        return None

    def update_markersize_func(self):
        """ Function to handle the update of the markers size. """

        event_source = self.sender()
        source_object = event_source.objectName()
        if '_unsel' in source_object:
            loc = 0
        elif '_sel' in source_object:
            loc = 1
        else:
            loc = 0
        for key, val in self.interactive_elements.items():
            if event_source in self.interactive_elements[key]['markersize_sl']:
                markersize = self.interactive_elements[key]['markersize_sl'][loc].value()
                break
        self.update_markersize(self.plot_axes[key], key, markersize, event_source)

        return None

    def update_opacity_func(self):
        """ Function to handle the update of the markers opacity. """

        event_source = self.sender()
        source_object = event_source.objectName()
        if '_unsel' in source_object:
            loc = 0
        elif '_sel' in source_object:
            loc = 1
        else:
            loc = 0
        for key, val in self.interactive_elements.items():
            if event_source in self.interactive_elements[key]['opacity_sl']:
                opacity = self.interactive_elements[key]['opacity_sl'][loc].value()/10
                break
        self.update_opacity(self.plot_axes[key], key, opacity, event_source)

        return None

    def update_linewidth_func(self):
        """ Function to handle the update of the lines widths. """

        event_source = self.sender()
        for key, val in self.interactive_elements.items():
            if event_source in self.interactive_elements[key]['linewidth_sl']:
                linewidth = self.interactive_elements[key]['linewidth_sl'][0].value()/10
                break
        self.update_linewidth(self.plot_axes[key], key, linewidth)

        return None

    def update_smooth_func(self):
        
        # get source
        event_source = self.sender()
        plot_type = event_source.objectName().split('_smooth')[0]
        for element in self.interactive_elements[plot_type]['smooth_sl']:
            smooth_window = element.value()
            break

        self.update_smooth_window(self.plot_axes[plot_type], plot_type, smooth_window, 
                                  self.read_instance.current_plot_options[plot_type])

        return None

    def update_plot_option(self):
        """ Function to handle the update of the plot options. """
        
        if not self.read_instance.block_MPL_canvas_updates:

            # get source
            event_source = self.sender()
            plot_type_alt = event_source.objectName().split('_options')[0]

            # correct perodic-violin name
            if plot_type_alt == 'periodic_violin':
                plot_type = 'periodic-violin'
            else:
                plot_type = copy.deepcopy(plot_type_alt)

            # an option is selected or there are options in previous to undo?
            if event_source.currentData() or self.read_instance.previous_plot_options[plot_type]:

                # get plot options (previous and currently selected)
                self.read_instance.current_plot_options[plot_type] = copy.deepcopy(event_source.currentData())
                all_plot_options = self.plot_characteristics[plot_type]['plot_options']
                
                for option in all_plot_options:
                    
                    # get index to raise errors and uncheck options
                    index = all_plot_options.index(option)

                    # if do not have selected station_station_data in memory, then no data has been read
                    # so return
                    if not hasattr(self, 'selected_station_data'):
                        msg = 'Select at least one station in the plot to apply options.'
                        show_message(self.read_instance, msg)
                        self.read_instance.block_MPL_canvas_updates = True
                        event_source.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                        self.read_instance.block_MPL_canvas_updates = False
                        return

                    # return from function if selected_station_data has not been updated for new species yet.
                    if self.read_instance.networkspeci not in self.selected_station_data:
                        return

                    # undo plot options that were selected before but not now
                    if ((option in self.read_instance.previous_plot_options[plot_type]) 
                        and (option not in self.read_instance.current_plot_options[plot_type])):
                        undo = True
                    # do plot options that were not selected before
                    elif ((option not in self.read_instance.previous_plot_options[plot_type]) 
                        and (option in self.read_instance.current_plot_options[plot_type])):
                        undo = False
                    # do nothing if options were selected before and now
                    elif ((option in self.read_instance.previous_plot_options[plot_type]) 
                        and (option in self.read_instance.current_plot_options[plot_type])):
                        continue
                    # do nothing if options were never selected
                    elif ((option not in self.read_instance.previous_plot_options[plot_type]) 
                        and (option not in self.read_instance.current_plot_options[plot_type])):
                        continue

                    # if plot type not in plot_elements, then return
                    if plot_type not in self.plot_elements:
                        return

                    # if no selected stations then remove all plot_elements for active plot_options,
                    # and then return
                    if (len(self.relative_selected_station_inds) == 0):
                        for active_type in self.plot_elements[plot_type]:
                            if active_type != 'active':
                                for data_label in self.plot_elements[plot_type][active_type]:
                                    for plot_option in self.read_instance.current_plot_options[plot_type]:
                                        if plot_option in self.plot_elements[plot_type][active_type][data_label]:
                                            for plot_element in self.plot_elements[plot_type][active_type][data_label][plot_option]:
                                                plot_element.remove()
                                            del self.plot_elements[plot_type][active_type][data_label][plot_option]
                        return 

                    # remove current option elements (both absolute and bias)
                    for active_type in self.plot_elements[plot_type]:
                        if active_type != 'active':
                            for data_label in self.plot_elements[plot_type][active_type]:
                                if option in self.plot_elements[plot_type][active_type][data_label]:
                                    for plot_element in self.plot_elements[plot_type][active_type][data_label][option]:
                                        plot_element.remove()
                                    del self.plot_elements[plot_type][active_type][data_label][option]

                    # if option is 'bias', then remove current option elements (both absolute and bias)
                    if option == 'bias':
                        for active_type in self.plot_elements[plot_type]:
                            if active_type != 'active':
                                for data_label in self.plot_elements[plot_type][active_type]:
                                    for plot_option in self.read_instance.current_plot_options[plot_type]:
                                        if plot_option in self.plot_elements[plot_type][active_type][data_label]:
                                            for plot_element in self.plot_elements[plot_type][active_type][data_label][plot_option]:
                                                plot_element.remove()
                                            del self.plot_elements[plot_type][active_type][data_label][plot_option] 

                    # options 'logy' and 'logx' 
                    # only plot if axis has all positive values
                    if (option == 'logy') or (option == 'logx'):
                        if isinstance(self.plot_axes[plot_type], dict):
                            for temporal_resolution, sub_ax in self.plot_axes[plot_type].items():
                                log_validity = self.plot.log_validity(sub_ax, option)
                                if log_validity:
                                    self.plot.log_axes(sub_ax, option, self.plot_characteristics[plot_type], undo=undo)
                                else:
                                    msg = "It is not possible to log the {0}-axis ".format(option[-1])
                                    msg += "in {0} with negative values.".format(plot_type)
                                    show_message(self.read_instance, msg)
                                    self.read_instance.block_MPL_canvas_updates = True
                                    event_source.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                                    self.read_instance.block_MPL_canvas_updates = False
                                    return None
                        else:
                            log_validity = self.plot.log_validity(self.plot_axes[plot_type], option)
                            if log_validity:
                                self.plot.log_axes(self.plot_axes[plot_type], option, self.plot_characteristics[plot_type], 
                                                undo=undo)
                            else:
                                msg = "It is not possible to log the {0}-axis ".format(option[-1])
                                msg += "in {0} with negative values.".format(plot_type)
                                show_message(self.read_instance, msg)
                                self.read_instance.block_MPL_canvas_updates = True
                                event_source.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                                self.read_instance.block_MPL_canvas_updates = False
                                return None

                    # option 'annotate'
                    # only plot if have selected stations (for map annotations)
                    elif option == 'annotate':
                        if not undo:
                            if isinstance(self.plot_axes[plot_type], dict):
                                for relevant_temporal_resolution, sub_ax in self.plot_axes[plot_type].items():
                                    if relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions:
                                        self.plot.annotation(sub_ax,
                                                            self.read_instance.networkspeci,
                                                            list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                                            plot_type,
                                                            self.plot_characteristics[plot_type], 
                                                            plot_options=self.read_instance.current_plot_options[plot_type])
                                        break
                            else:
                                self.plot.annotation(self.plot_axes[plot_type], 
                                                    self.read_instance.networkspeci,
                                                    list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                                    plot_type,
                                                    self.plot_characteristics[plot_type], 
                                                    plot_options=self.read_instance.current_plot_options[plot_type])

                    # option 'smooth'
                    elif option == 'smooth':
                        if not undo:
                            self.plot.smooth(self.plot_axes[plot_type], 
                                             self.read_instance.networkspeci,
                                             list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                             plot_type,
                                             self.plot_characteristics[plot_type], 
                                             plot_options=self.read_instance.current_plot_options[plot_type])
                          
                    # option 'regression'
                    elif option == 'regression':
                        if not undo:
                            self.plot.linear_regression(self.plot_axes[plot_type], 
                                                        self.read_instance.networkspeci,
                                                        list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                                        plot_type,
                                                        self.plot_characteristics[plot_type],  
                                                        plot_options=self.read_instance.current_plot_options[plot_type])

                    # option 'bias'
                    elif option == 'bias':

                        # firstly if just 1 data label then cannot make bias plot 
                        if len(list(self.selected_station_data[self.read_instance.networkspeci].keys())) == 1:
                            msg = 'It is not possible to make a bias plot with just observations loaded.'
                            show_message(self.read_instance, msg)
                            self.read_instance.block_MPL_canvas_updates = True
                            event_source.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                            self.read_instance.block_MPL_canvas_updates = False
                            self.plot_elements[plot_type]['active'] = 'absolute'
                            self.read_instance.current_plot_options[plot_type].remove('bias')

                            # create other active plot option elements for now absolute plot (if do not already exist)
                            self.redraw_active_options(list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                                       plot_type, 'absolute', self.read_instance.current_plot_options[plot_type])

                        # if bias option is enabled then first check if bias elements stored
                        elif not undo:

                            # update active (bias)
                            self.plot_elements[plot_type]['active'] = 'bias' 

                            # handle some special cases for periodic plot
                            if plot_type in ['periodic']:

                                # get currently selected periodic statistic name
                                base_zstat = self.periodic_stat.currentText()
                                zstat = get_z_statistic_comboboxes(base_zstat, second_data_label='model')

                                # get zstat information 
                                zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat) 

                                # if get_z_statistic_type == 'expbias' then return as bias already plotted
                                if z_statistic_type == 'expbias':
                                    self.read_instance.block_MPL_canvas_updates = True
                                    event_source.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                                    self.read_instance.block_MPL_canvas_updates = False
                                    self.plot_elements[plot_type]['active'] = 'absolute'

                            # handle some special cases for statsummary plot
                            elif plot_type in ['statsummary']:

                                # create structure to store data for statsummary plot
                                relevant_zstats = self.plot_characteristics[plot_type]['experiment_bias']
                                stats_df = {relevant_zstat:[] for relevant_zstat in relevant_zstats}

                            # get plotting function for specific plot
                            if plot_type == 'statsummary':
                                func = getattr(self.plot, 'make_table')
                            else:
                                func = getattr(self.plot, 'make_{}'.format(plot_type.split('-')[0]))

                            # iterate through valid data labels 
                            first_data_label = True
                            for data_label in list(self.selected_station_data[self.read_instance.networkspeci].keys()) + ['ALL']:

                                # hide absolute plot elements
                                if data_label in self.plot_elements[plot_type]['absolute']:
                                    for element_type in self.plot_elements[plot_type]['absolute'][data_label]:
                                        for element in self.plot_elements[plot_type]['absolute'][data_label][element_type]:
                                            element.set_visible(False)

                                # if have bias elements pre-stored then simply show bias elements (if data label on legend is active)
                                set_bias = False
                                if 'bias' in self.plot_elements[plot_type]:
                                    if data_label in self.plot_elements[plot_type]['bias']:
                                        set_bias = True
                                        if (data_label in self.plot_elements['data_labels_active']) or (data_label == 'ALL'):
                                            for element_type in self.plot_elements[plot_type]['bias'][data_label]:
                                                for element in self.plot_elements[plot_type]['bias'][data_label][element_type]:
                                                    element.set_visible(True)

                                # if do not already have bias elements, then make them (tracking plot elements also) 
                                if not set_bias:
                                    # skip if data array is observations / ALL
                                    if data_label in ['observations', 'ALL']:
                                        continue

                                    # gather data for statsummary plot
                                    if plot_type in ['statsummary']:
                                        for relevant_zstat in relevant_zstats:
                                            # if relevant stat is expbias stat, then ensure temporal colocation is active
                                            # otherwise set value as NaN
                                            if (relevant_zstat in expbias_stats) & (not self.read_instance.temporal_colocation):
                                                stat_val = np.NaN
                                            else:
                                                stat_val = self.selected_station_data[self.read_instance.networkspeci][data_label]['all'][relevant_zstat][0]
                                            stats_df[relevant_zstat].append(stat_val)

                                    # other plot types
                                    else:
                                        # call plotting function
                                        if plot_type in ['periodic']:
                                            func(self.plot_axes[plot_type], self.read_instance.networkspeci, data_label, 
                                                self.plot_characteristics[plot_type], zstat=zstat, plot_options=self.read_instance.current_plot_options[plot_type], 
                                                first_data_label=first_data_label)
                                        else: 
                                            func(self.plot_axes[plot_type], self.read_instance.networkspeci, data_label, 
                                                self.plot_characteristics[plot_type], plot_options=self.read_instance.current_plot_options[plot_type], 
                                                first_data_label=first_data_label)
                                        first_data_label = False

                            # make statsummary bias plot (if not previously made)
                            if (plot_type in ['statsummary']) & ('bias' not in self.plot_elements[plot_type]):
                                if len(stats_df[list(stats_df.keys())[0]]) > 0:
                                    index = [data_label for data_label in self.selected_station_data[self.read_instance.networkspeci] if data_label != 'observations']
                                    stats_df = pd.DataFrame(data=stats_df,index=index)
                                    func(self.plot_axes[plot_type], stats_df, self.plot_characteristics[plot_type], 
                                         plot_options=self.read_instance.current_plot_options[plot_type], statsummary=True)

                            # create other active plot option elements for bias plot (if do not already exist)
                            self.redraw_active_options(list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                                       plot_type, 'bias', self.read_instance.current_plot_options[plot_type])

                        # if bias option is not enabled then hide bias plot elements and show absolute plots again
                        else:

                            # update active (absolute)
                            self.plot_elements[plot_type]['active'] = 'absolute' 

                            # iterate through valid data labels
                            for data_label in list(self.selected_station_data[self.read_instance.networkspeci].keys()) + ['ALL']:

                                # hide bias plot elements 
                                if 'bias' in self.plot_elements[plot_type]:
                                    if data_label in self.plot_elements[plot_type]['bias']:
                                        for element_type in self.plot_elements[plot_type]['bias'][data_label]:
                                            for element in self.plot_elements[plot_type]['bias'][data_label][element_type]:
                                                element.set_visible(False)

                                # show absolute plot elements (if data label on legend is active)
                                if data_label in self.plot_elements[plot_type]['absolute']:
                                    if (data_label in self.plot_elements['data_labels_active']) or (data_label == 'ALL'):
                                        for element_type in self.plot_elements[plot_type]['absolute'][data_label]:
                                            for element in self.plot_elements[plot_type]['absolute'][data_label][element_type]:
                                                element.set_visible(True)

                            # create other active plot option elements for absolute plot (if do not already exist)
                            self.redraw_active_options(list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                                       plot_type, 'absolute', self.read_instance.current_plot_options[plot_type])

                    # reset axes limits (harmonising across subplots for periodic plots) 
                    if plot_type != 'map':
                        if plot_type == 'periodic-violin':
                            self.plot.harmonise_xy_lims_paradigm(self.plot_axes[plot_type], plot_type, 
                                                                 self.plot_characteristics[plot_type], 
                                                                 self.read_instance.current_plot_options[plot_type], 
                                                                 ylim=[self.selected_station_data_min[self.read_instance.networkspeci], 
                                                                       self.selected_station_data_max[self.read_instance.networkspeci]],
                                                                 relim=True, autoscale_x=True)
                        elif plot_type == 'scatter':
                            self.plot.harmonise_xy_lims_paradigm(self.plot_axes[plot_type], plot_type, 
                                                                 self.plot_characteristics[plot_type], 
                                                                 self.read_instance.current_plot_options[plot_type], 
                                                                 relim=True)
                        else:
                            self.plot.harmonise_xy_lims_paradigm(self.plot_axes[plot_type], plot_type, 
                                                                 self.plot_characteristics[plot_type], 
                                                                 self.read_instance.current_plot_options[plot_type], 
                                                                 relim=True, autoscale=True)                       

                # save current plot options as previous
                self.read_instance.previous_plot_options[plot_type] = self.read_instance.current_plot_options[plot_type]

                # draw changes
                self.figure.canvas.draw()

        return None

    def redraw_active_options(self, data_labels, plot_type, active, plot_options):
        """ Redraw active plot option elements when moving between absolute and bias plots,
            if do not already exist.
        """

        # if 'bias' is active, remove 'observations' from data_labels
        data_labels_alt = copy.deepcopy(data_labels)
        if active == 'bias':
            data_labels_alt.remove('observations')

        # iterate through plot_options
        for plot_option in plot_options:

            if plot_option == 'annotate':
                if isinstance(self.plot_axes[plot_type], dict):
                    for relevant_temporal_resolution, sub_ax in self.plot_axes[plot_type].items():
                        if relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions:
                            self.plot.annotation(sub_ax,
                                                 self.read_instance.networkspeci,
                                                 data_labels, 
                                                 plot_type,
                                                 self.plot_characteristics[plot_type], 
                                                 plot_options=plot_options)
                            break
                else:
                    self.plot.annotation(self.plot_axes[plot_type], 
                                         self.read_instance.networkspeci,
                                         data_labels, 
                                         plot_type,
                                         self.plot_characteristics[plot_type],
                                         plot_options=plot_options)

            elif plot_option == 'smooth':
                self.plot.smooth(self.plot_axes[plot_type], 
                                 self.read_instance.networkspeci,
                                 data_labels_alt,
                                 plot_type,
                                 self.plot_characteristics[plot_type], 
                                 plot_options=plot_options)
            
            elif plot_option == 'regression':
                self.plot.linear_regression(self.plot_axes[plot_type], 
                                            self.read_instance.networkspeci,
                                            data_labels_alt, 
                                            plot_type,
                                            self.plot_characteristics[plot_type],  
                                            plot_options=plot_options)

    def update_markersize(self, ax, plot_type, markersize, event_source):
        """ Update markers size for each plot type. """
        
        # set markersize
        if plot_type in ['timeseries', 'periodic', 'scatter', 'periodic-violin']:

            if isinstance(ax, dict):
                for sub_ax in ax.values():
                    for line in sub_ax.lines:
                        line.set_markersize(markersize)
            else:
                for line in ax.lines:
                    line.set_markersize(markersize)

            # update characteristics per plot type
            # this is made to keep the changes when selecting stations with lasso
            if plot_type in ['timeseries', 'periodic', 'scatter']:
                self.plot_characteristics[plot_type]['plot']['markersize'] = markersize
            elif plot_type == 'periodic-violin':
                self.plot_characteristics[plot_type]['plot']['p50']['markersize'] = markersize

        elif plot_type == 'map':

            markersizes = self.plot_axes['map'].collections[-1].get_sizes()

            # zero selected and unselected stations
            if event_source == self.interactive_elements[plot_type]['markersize_sl'][0]:

                # actually have zero selected stations currently?
                # if so, update active markersizes
                if len(self.absolute_selected_station_inds) == 0:
                    markersizes[:] = markersize
                    for collection in self.plot_axes['map'].collections:
                        if isinstance(collection, matplotlib.collections.PathCollection):
                            collection.set_sizes(markersizes)

                # actually have selected stations currently?
                # if so, update active opacities
                elif (len(self.absolute_non_selected_station_inds) > 0) & (len(self.absolute_selected_station_inds) > 0):
                    markersizes[self.absolute_non_selected_station_inds] = markersize
                    for collection in self.plot_axes['map'].collections:
                        if isinstance(collection, matplotlib.collections.PathCollection):
                            collection.set_sizes(markersizes)

                # update characteristics per plot type for unselected stations
                self.plot_characteristics['map']['marker_unselected']['s'] = markersize

                # update characteristics per plot type for zero selected stations
                self.plot_characteristics['map']['marker_zero_stations_selected']['s'] = markersize

            # selected stations
            elif event_source == self.interactive_elements[plot_type]['markersize_sl'][1]:

                # actually have selected stations currently?
                # if so, update active opacities
                if len(self.absolute_selected_station_inds) > 0:
                    markersizes[self.absolute_selected_station_inds] = markersize
                    for collection in self.plot_axes['map'].collections:
                        if isinstance(collection, matplotlib.collections.PathCollection):
                            collection.set_sizes(markersizes)

                # update characteristics per plot type
                self.plot_characteristics['map']['marker_selected']['s'] = markersize

        # redraw points
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

        return None

    def update_opacity(self, ax, plot_type, opacity, event_source):
        """ Update markers opacity for each plot type. """

        # set opacity
        if plot_type == 'map':

            opacities = self.plot_axes['map'].collections[-1].get_facecolor()

            # zero selected and unselected stations
            if event_source == self.interactive_elements[plot_type]['opacity_sl'][0]:

                # actually have zero selected stations currently?
                # if so, update active opacities
                if len(self.absolute_selected_station_inds) == 0:
                    opacities[:, -1] = opacity
                    for collection in self.plot_axes['map'].collections:
                        if isinstance(collection, matplotlib.collections.PathCollection):
                            collection.set_facecolor(opacities)

                # actually have selected stations currently?
                # if so, update active opacities
                elif (len(self.absolute_non_selected_station_inds) > 0) & (len(self.absolute_selected_station_inds) > 0):
                    opacities[self.absolute_non_selected_station_inds, -1] = opacity
                    for collection in self.plot_axes['map'].collections:
                        if isinstance(collection, matplotlib.collections.PathCollection):
                            collection.set_facecolor(opacities)

                # update characteristics per plot type for unselected stations
                self.plot_characteristics['map']['marker_zero_stations_selected']['alpha'] = opacity

                # update characteristics per plot type for zero selected stations
                self.plot_characteristics['map']['marker_unselected']['alpha'] = opacity

            # selected stations
            elif event_source == self.interactive_elements[plot_type]['opacity_sl'][1]:

                # actually have selected stations currently?
                # if so, update active opacities
                if len(self.absolute_selected_station_inds) > 0:
                    opacities[self.absolute_selected_station_inds, -1] = opacity
                    for collection in self.plot_axes['map'].collections:
                        if isinstance(collection, matplotlib.collections.PathCollection):
                            collection.set_facecolor(opacities)

                # update characteristics per plot type
                self.plot_characteristics['map']['marker_selected']['alpha'] = opacity

        # redraw points
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

        return None

    def update_linewidth(self, ax, plot_type, linewidth):
        """ Update line widths for each plot type. """
        
        # set linewidth
        if isinstance(ax, dict):
            for sub_ax in ax.values():
                for line in sub_ax.lines:
                    line.set_linewidth(linewidth)
        else:
            for line in ax.lines:
                line.set_linewidth(linewidth)

        # update characteristics per plot type
        if plot_type == 'periodic-violin':
            self.plot_characteristics[plot_type]['plot']['p50']['linewidth'] = linewidth
        else:
            self.plot_characteristics[plot_type]['plot']['linewidth'] = linewidth

        # redraw points
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

        return None

    def update_smooth_window(self, ax, plot_type, smooth_window, plot_options):

        # update characteristics per plot type
        self.plot_characteristics[plot_type]['smooth']['window'] = smooth_window
        
        # get index of smooth in plot options
        all_plot_options = self.plot_characteristics[plot_type]['plot_options']
        index = all_plot_options.index('smooth')

        # remove smooth plot option
        self.timeseries_options.model().item(index).setCheckState(QtCore.Qt.Unchecked)

        # create smooth lines
        if smooth_window > 0:
            
            # add smooth plot option
            self.timeseries_options.model().item(index).setCheckState(QtCore.Qt.Checked)

        # redraw points
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

        return None

    def update_violin_widths(self, ax, plot_type, width):
        """ Update violin widths for violin plots. """
        
        # set violin widths
        if plot_type == 'periodic-violin':
            if isinstance(ax, dict):
                for sub_ax in ax.values():
                    widths = np.full(len(self.active_map_valid_station_inds), width)
                    for collection in sub_ax.collections:
                        # TODO: Change widths - Try to find a function like set_widths
                        pass
                       
            else:
                for collection in ax.collections:
                    # TODO: Change widths - Try to find a function like set_widths
                    pass

        # update characteristics per plot type
        self.plot_characteristics[plot_type]['plot']['violin']['widths'] = widths

        # redraw points
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

        return None

    def save_axis_figure_dialog(self, plot_type, relevant_temporal_resolution=None):
        """ Function to create the dialog box to save each plot figure. """

        default_filename = '{0}-{1}-{2}-{3}-{4}-{5}-{6}.png'.format(self.read_instance.network[0],
                                                                    self.read_instance.species[0],
                                                                    self.read_instance.resolution, 
                                                                    self.read_instance.start_date, 
                                                                    self.read_instance.end_date,
                                                                    plot_type, str(relevant_temporal_resolution))

        if relevant_temporal_resolution == None:
            default_filename = default_filename.split('-None')[0] + '.png'

        options = QtWidgets.QFileDialog.Options()
        options |=  QtWidgets.QFileDialog.DontUseNativeDialog
        figure_path, _ =  QtWidgets.QFileDialog.getSaveFileName(self.read_instance, "Choose folder to save figure", 
                                                                default_filename, "All Files (*);;Figures (*.png)", 
                                                                options=options)
        if figure_path:
            return figure_path

    def save_axis_figure_func(self):
        """ Function to save each plot figure. """
        
        # get option and plot names
        event_source = self.sender()
        plot_type = event_source.objectName().split('_save')[0]
        if plot_type == 'periodic_violin':
            plot_type = 'periodic-violin'
        
        # set extent expansion
        for i, position in enumerate([self.read_instance.position_1, 
                                      self.read_instance.position_2, 
                                      self.read_instance.position_3, 
                                      self.read_instance.position_4, 
                                      self.read_instance.position_5]):
            if plot_type == position:
                if i + 1 == 1:
                    expand_x, expand_y = 1.4, 1.4
                elif i + 1 == 2:
                    if plot_type == 'scatter':
                        expand_x, expand_y = 1.30, 1.25
                    else:
                        expand_x, expand_y = 1.2, 1.3
                else:
                    if plot_type == 'scatter':
                        expand_x, expand_y = 1.25, 1.2
                    else:
                        expand_x, expand_y = 1.2, 1.2
                continue

        # remove titles
        for key in self.read_instance.active_dashboard_plots:
            if isinstance(self.plot_axes[key], dict):
                for relevant_temporal_resolution, sub_ax in self.plot_axes[key].items():
                    if relevant_temporal_resolution in ['hour']:
                        sub_ax.set_title(label='', 
                                         fontsize=self.plot_characteristics[key]['axis_title']['fontsize'],
                                         loc=self.plot_characteristics[key]['axis_title']['loc'])
            else:
                self.plot_axes[key].set_title(label='',
                                              fontsize=self.plot_characteristics[key]['axis_title']['fontsize'],
                                              loc=self.plot_characteristics[key]['axis_title']['loc'])
       
        # hide colourbar
        if plot_type != 'map':
            self.plot_axes['cb'].set_visible(False)

        # draw changes
        self.figure.canvas.draw()

        # make screenshot and save
        if isinstance(self.plot_axes[plot_type], dict):
            for relevant_temporal_resolution, sub_ax in self.plot_axes[plot_type].items():
                extent = sub_ax.get_window_extent().transformed(self.figure.dpi_scale_trans.inverted())
                if relevant_temporal_resolution == 'hour':
                    expand_x, expand_y = 1.15, 1.3
                elif relevant_temporal_resolution == 'dayofweek':
                    expand_x, expand_y = 1.25, 1.3
                elif relevant_temporal_resolution == 'month':
                    expand_x, expand_y = 1.15, 1.3
                
                # get folder where figure will be saved
                figure_path = self.save_axis_figure_dialog(plot_type, relevant_temporal_resolution)
                
                # save figure
                if figure_path is not None:
                    self.figure.savefig(figure_path, bbox_inches=extent.expanded(expand_x, expand_y))
        else:
            extent = self.plot_axes[plot_type].get_window_extent().transformed(self.figure.dpi_scale_trans.inverted())
            
            # get folder where figure will be saved
            figure_path = self.save_axis_figure_dialog(plot_type)

            # save figure
            if figure_path is not None:
                self.figure.savefig(figure_path, bbox_inches=extent.expanded(expand_x, expand_y))

        # add titles
        for key in self.read_instance.active_dashboard_plots:
            if isinstance(self.plot_axes[key], dict):
                for relevant_temporal_resolution, sub_ax in self.plot_axes[key].items():
                    if relevant_temporal_resolution in ['hour']:
                        sub_ax.set_title(**self.plot_characteristics[key]['axis_title'])
            else:
                self.plot_axes[key].set_title(**self.plot_characteristics[key]['axis_title'])

        # show colourbar
        self.plot_axes['cb'].set_visible(True)

        # draw changes
        self.figure.canvas.draw()

        return None

    def create_station_annotation(self):
        """ Create annotation at (0, 0) that will be updated later. """

        # get axes transformation
        transform = self.datacrs._as_mpl_transform(self.plot_axes['map'])

        # in the newest version of matplotlib, s corresponds to text
        self.station_annotation = self.plot_axes['map'].annotate(s='', xy=(0, 0), xycoords=transform,
                                                                 **self.plot_characteristics['map']['stations_annotate'],
                                                                 bbox={**self.plot_characteristics['map']['stations_annotate_bbox']},
                                                                 arrowprops={**self.plot_characteristics['map']['stations_annotate_arrowprops']})

        self.station_annotation.set_visible(False)

        return None

    def update_station_annotation(self, annotation_index):
        """ Update annotation for each station that is hovered. """

        # retrieve stations references and coordinates
        station_names = self.read_instance.metadata_in_memory[self.read_instance.networkspeci]['station_name'][self.active_map_valid_station_inds[annotation_index['ind'][0]]]
        station_reference = self.read_instance.station_references[self.read_instance.networkspeci][self.active_map_valid_station_inds[annotation_index['ind'][0]]]
        station_location = self.plot.stations_scatter.get_offsets()[annotation_index['ind'][0]]
        station_value = self.z_statistic[annotation_index['ind'][0]]

        # update location
        self.station_annotation.xy = station_location

        # update bbox position
        self.read_instance.map_extent = self.plot.get_map_extent(self.plot_axes['map'])
        lat_min = self.read_instance.map_extent[2]
        lat_max = self.read_instance.map_extent[3]
        if station_location[1] > ((lat_max + lat_min) / 2):
            self.station_annotation.set_y(-10)
            self.station_annotation.set_va('top')
        else:
            self.station_annotation.set_y(10)
            self.station_annotation.set_va('bottom')

        # create annotation text
        text_label = ('Station: {0}\n').format(station_names[pd.notnull(station_names)][0])
        text_label += ('Reference: {0}\n').format(station_reference)
        text_label += ('Longitude: {0:.2f}\n').format(station_location[0])
        text_label += ('Latitude: {0:.2f}\n').format(station_location[1])
        text_label += ('{0}: {1:.2f}').format(self.map_z_stat.currentText(), station_value)
        self.station_annotation.set_text(text_label)

        return None
        
    def hover_map_annotation(self, event):
        """ Show or hide annotation for each station that is hovered. """

        if (not self.lock_zoom) & (event.inaxes == self.plot_axes['map']):

            # activate hover over map
            if (hasattr(self.plot, 'stations_scatter')):
                    
                is_contained, annotation_index = self.plot.stations_scatter.contains(event)
                
                if is_contained:
                    # update annotation if hovered
                    self.update_station_annotation(annotation_index)
                    self.station_annotation.set_visible(True)
                else:
                    # hide annotation if not hovered
                    if self.station_annotation.get_visible():
                        self.station_annotation.set_visible(False)

                # redraw points
                self.figure.canvas.draw()

        return None

    def create_timeseries_annotation(self):
        """ Create annotation at (0, 0) that will be updated later. """

        # in the newest version of matplotlib, s corresponds to text
        self.timeseries_annotation = self.plot_axes['timeseries'].annotate(s='', xy=(0, 0), xycoords='data',
                                                                           **self.plot_characteristics['timeseries']['marker_annotate'],
                                                                           bbox={**self.plot_characteristics['timeseries']['marker_annotate_bbox']},
                                                                           arrowprops={**self.plot_characteristics['timeseries']['marker_annotate_arrowprops']})
        self.timeseries_annotation.set_visible(False)

        return None

    def create_timeseries_annotation_vline(self):
        """ Create annotation vertical line at (0, 0) that will be updated later. """

        # add vertical line
        self.timeseries_vline = self.plot_axes['timeseries'].axvline(0, **self.plot_characteristics['timeseries']['marker_annotate_vline'])
        self.timeseries_vline.set_visible(False)

        return None

    def update_timeseries_annotation(self, annotation_index):
        """ Update annotation for each timeseries point that is hovered. """
        
        for data_label in self.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.timeseries_annotate_data_label:
                
                # skip observations for bias plot
                if self.plot_elements['timeseries']['active'] == 'bias' and data_label == 'observations':
                    continue
                
                # do not annotate if plot is cleared
                if data_label not in self.plot_elements['timeseries'][self.plot_elements['timeseries']['active']].keys():
                    continue

                # retrieve time and concentration
                line = self.plot_elements['timeseries'][self.plot_elements['timeseries']['active']][data_label]['plot'][0]
                time = line.get_xdata()[annotation_index['ind'][0]]
                concentration = line.get_ydata()[annotation_index['ind'][0]]

                # update location
                self.timeseries_annotation.xy = (time, concentration)
                self.timeseries_vline.set_xdata(time)

                # update bbox position
                time_middle = line.get_xdata()[math.floor((len(line.get_xdata()) - 1)/2)]
                if time > time_middle:
                    self.timeseries_annotation.set_x(-10)
                    self.timeseries_annotation.set_ha('right')
                else:
                    self.timeseries_annotation.set_x(10)
                    self.timeseries_annotation.set_ha('left')

                # create annotation text
                text_label = ('Time: {0}').format(time.astype('datetime64[us]').astype(datetime.datetime).strftime("%m/%d/%Y %H:%M:%S"))
        
        for data_label in self.plot_elements['data_labels_active']:
            
            # skip observations for bias plot
            if self.plot_elements['timeseries']['active'] == 'bias' and data_label == 'observations':
                continue

            # do not annotate if plot is cleared
            if data_label not in self.plot_elements['timeseries'][self.plot_elements['timeseries']['active']].keys():
                continue

            # retrieve concentration
            line = self.plot_elements['timeseries'][self.plot_elements['timeseries']['active']][data_label]['plot'][0]
            concentration = line.get_ydata()[np.where(line.get_xdata() == time)[0]]

            # for all labels if there is data
            if len(concentration) >= 1:
                if data_label != 'observations':
                    exp_alias = self.read_instance.experiments[data_label]
                    text_label += ('\n{0}: {1:.2f}').format(exp_alias, concentration[0])
                else:
                    text_label += ('\n{0}: {1:.2f}').format(self.plot_characteristics['legend']['handles']['obs_label'], 
                                                            concentration[0])

        self.timeseries_annotation.set_text(text_label)

        return None

    def hover_timeseries_annotation(self, event):
        """ Show or hide annotation for each point that is hovered in the timeseries plot. """

        # activate hover over timeseries
        if ('timeseries' in self.read_instance.active_dashboard_plots):
            if event.inaxes == self.plot_axes['timeseries']:
                if ((hasattr(self.plot, 'timeseries_plot')) and ('timeseries' in self.plot_elements)
                    and (self.lock_timeseries_annotation == False)):

                    # lock annotation
                    self.lock_timeseries_annotation = True
                    is_contained = False

                    for data_label in self.plot_elements['data_labels_active']:

                        # skip observations for bias plot
                        if self.plot_elements['timeseries']['active'] == 'bias' and data_label == 'observations':
                            continue

                        # do not annotate if plot is cleared
                        if data_label not in self.plot_elements['timeseries'][self.plot_elements['timeseries']['active']].keys():
                            continue

                        line = self.plot_elements['timeseries'][self.plot_elements['timeseries']['active']][data_label]['plot'][0]
                        is_contained, annotation_index = line.contains(event)
                        if is_contained:
                            self.timeseries_annotate_data_label = data_label
                            break
                    
                    if is_contained:
                        # update annotation if hovered
                        self.update_timeseries_annotation(annotation_index)
                        self.timeseries_annotation.set_visible(True)
                        self.timeseries_vline.set_visible(True)
                    else:
                        # hide annotation if not hovered
                        if self.timeseries_annotation.get_visible():
                            self.timeseries_annotation.set_visible(False)
                            self.timeseries_vline.set_visible(False)
                            
                    # redraw points
                    self.figure.canvas.draw()
                    self.figure.canvas.flush_events()
                        
                    # unlock annotation 
                    self.lock_timeseries_annotation = False

        return None

    def create_scatter_annotation(self):
        """ Create annotation at (0, 0) that will be updated later. """

        # in the newest version of matplotlib, s corresponds to text
        self.scatter_annotation = self.plot_axes['scatter'].annotate(s='', xy=(0, 0), xycoords='data',
                                                                     **self.plot_characteristics['scatter']['marker_annotate'],
                                                                     bbox={**self.plot_characteristics['scatter']['marker_annotate_bbox']},
                                                                     arrowprops={**self.plot_characteristics['scatter']['marker_annotate_arrowprops']})
        self.scatter_annotation.set_visible(False)

        return None

    def update_scatter_annotation(self, annotation_index):

        for data_label in self.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.scatter_annotate_data_label:
                
                # do not annotate if plot is cleared
                if data_label not in self.plot_elements['scatter'][self.plot_elements['scatter']['active']].keys():
                    continue

                # retrieve time and concentration
                line = self.plot_elements['scatter'][self.plot_elements['scatter']['active']][data_label]['plot'][0]
                concentration_x = line.get_xdata()[annotation_index['ind'][0]]
                concentration_y = line.get_ydata()[annotation_index['ind'][0]]

                # update location
                self.scatter_annotation.xy = (concentration_x, concentration_y)

                # update bbox position
                concentration_x_middle = line.get_xdata()[math.floor((len(line.get_xdata()) - 1)/2)]
                if concentration_x > concentration_x_middle:
                    self.scatter_annotation.set_x(-10)
                    self.scatter_annotation.set_ha('right')
                else:
                    self.scatter_annotation.set_x(10)
                    self.scatter_annotation.set_ha('left')

                # create annotation text
                # observations label
                text_label = ('{0}: {1:.2f}').format(self.plot_characteristics['legend']['handles']['obs_label'], 
                                                       concentration_x)
                # experiment label
                exp_alias = self.read_instance.experiments[data_label]
                text_label += ('\n{0}: {1:.2f}').format(exp_alias, concentration_y)
    
        self.scatter_annotation.set_text(text_label)

        return None

    def hover_scatter_annotation(self, event):
        """ Show or hide annotation for each point that is hovered in the scatter plot. """
        
        # activate hover over scatter
        if ('scatter' in self.read_instance.active_dashboard_plots):
            if event.inaxes == self.plot_axes['scatter']:
                if ((hasattr(self.plot, 'scatter_plot')) and ('scatter' in self.plot_elements)
                    and (self.lock_scatter_annotation == False)):

                    # lock annotation
                    self.lock_scatter_annotation = True
                    is_contained = False

                    for data_label in self.plot_elements['data_labels_active']:

                        # do not annotate if plot is cleared
                        if data_label not in self.plot_elements['scatter'][self.plot_elements['scatter']['active']].keys():
                            continue

                        line = self.plot_elements['scatter'][self.plot_elements['scatter']['active']][data_label]['plot'][0]
                        is_contained, annotation_index = line.contains(event)
                        if is_contained:
                            self.scatter_annotate_data_label = data_label
                            break
                    
                    if is_contained:
                        # update annotation if hovered
                        self.update_scatter_annotation(annotation_index)
                        self.scatter_annotation.set_visible(True)
                    else:
                        # hide annotation if not hovered
                        if self.scatter_annotation.get_visible():
                            self.scatter_annotation.set_visible(False)
                            
                    # redraw points
                    self.figure.canvas.draw()
                    self.figure.canvas.flush_events()
                        
                    # unlock annotation 
                    self.lock_scatter_annotation = False

        return None

    def create_distribution_annotation(self):
        """ Create annotation at (0, 0) that will be updated later. """

        # in the newest version of matplotlib, s corresponds to text
        self.distribution_annotation = self.plot_axes['distribution'].annotate(s='', xy=(0, 0), xycoords='data',
                                                                               **self.plot_characteristics['distribution']['marker_annotate'],
                                                                               bbox={**self.plot_characteristics['distribution']['marker_annotate_bbox']},
                                                                               arrowprops={**self.plot_characteristics['distribution']['marker_annotate_arrowprops']})
        self.distribution_annotation.set_visible(False)

        return None

    def create_distribution_annotation_vline(self):
        """ Create annotation vertical line at (0, 0) that will be updated later. """

        # add vertical line
        self.distribution_vline = self.plot_axes['distribution'].axvline(0, **self.plot_characteristics['distribution']['marker_annotate_vline'])
        self.distribution_vline.set_visible(False)

        return None

    def update_distribution_annotation(self, annotation_index):
        """ Update annotation for each distribution point that is hovered. """
        
        for data_label in self.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.distribution_annotate_data_label:
                
                # skip observations for bias plot
                if self.plot_elements['distribution']['active'] == 'bias' and data_label == 'observations':
                    continue
                
                # do not annotate if plot is cleared
                if data_label not in self.plot_elements['distribution'][self.plot_elements['distribution']['active']].keys():
                    continue

                # retrieve time and concentration
                line = self.plot_elements['distribution'][self.plot_elements['distribution']['active']][data_label]['plot'][0]
                concentration = line.get_xdata()[annotation_index['ind'][0]]
                density = line.get_ydata()[annotation_index['ind'][0]]

                # update location
                self.distribution_annotation.xy = (concentration, density)
                self.distribution_vline.set_xdata(concentration)

                # update bbox position
                concentration_middle = line.get_xdata()[math.floor((len(line.get_xdata()) - 1)/2)]
                if concentration > concentration_middle:
                    self.distribution_annotation.set_x(-10)
                    self.distribution_annotation.set_ha('right')
                else:
                    self.distribution_annotation.set_x(10)
                    self.distribution_annotation.set_ha('left')

                # create annotation text
                text_label = ('{0}: {1:.2f}').format(self.read_instance.species[0], concentration)
        
        for data_label in self.plot_elements['data_labels_active']:
            
            # skip observations for bias plot
            if self.plot_elements['distribution']['active'] == 'bias' and data_label == 'observations':
                continue

            # do not annotate if plot is cleared
            if data_label not in self.plot_elements['distribution'][self.plot_elements['distribution']['active']].keys():
                continue

            # retrieve density
            line = self.plot_elements['distribution'][self.plot_elements['distribution']['active']][data_label]['plot'][0]
            density = line.get_ydata()[np.where(line.get_xdata() == concentration)[0]]

            # for all labels if there is data
            if len(density) >= 1:
                if data_label != 'observations':
                    exp_alias = self.read_instance.experiments[data_label]
                    text_label += ('\n{0}: {1:.2f}').format(exp_alias, density[0])
                else:
                    text_label += ('\n{0}: {1:.2f}').format(self.plot_characteristics['legend']['handles']['obs_label'], 
                                                            density[0])

        self.distribution_annotation.set_text(text_label)

        return None

    def hover_distribution_annotation(self, event):
        """ Show or hide annotation for each point that is hovered in the distribution plot. """

        # activate hover over distribution
        if ('distribution' in self.read_instance.active_dashboard_plots):
            if event.inaxes == self.plot_axes['distribution']:
                if ((hasattr(self.plot, 'distribution_plot')) and ('distribution' in self.plot_elements)
                    and (self.lock_distribution_annotation == False)):

                    # lock annotation
                    self.lock_distribution_annotation = True
                    is_contained = False

                    for data_label in self.plot_elements['data_labels_active']:

                        # skip observations for bias plot
                        if self.plot_elements['distribution']['active'] == 'bias' and data_label == 'observations':
                            continue

                        # do not annotate if plot is cleared
                        if data_label not in self.plot_elements['distribution'][self.plot_elements['distribution']['active']].keys():
                            continue

                        line = self.plot_elements['distribution'][self.plot_elements['distribution']['active']][data_label]['plot'][0]
                        is_contained, annotation_index = line.contains(event)
                        if is_contained:
                            self.distribution_annotate_data_label = data_label
                            break
                    
                    if is_contained:
                        # update annotation if hovered
                        self.update_distribution_annotation(annotation_index)
                        self.distribution_annotation.set_visible(True)
                        self.distribution_vline.set_visible(True)
                    else:
                        # hide annotation if not hovered
                        if self.distribution_annotation.get_visible():
                            self.distribution_annotation.set_visible(False)
                            self.distribution_vline.set_visible(False)
                            
                    # redraw points
                    self.figure.canvas.draw()
                    self.figure.canvas.flush_events()
                        
                    # unlock annotation 
                    self.lock_distribution_annotation = False

        return None

    def create_periodic_annotation(self):
        """ Create annotation at (0, 0) that will be updated later. """

        # in the newest version of matplotlib, s corresponds to text
        self.periodic_annotation = dict()
        for resolution in self.plot_axes['periodic'].keys():
            self.periodic_annotation[resolution] = self.plot_axes['periodic'][resolution].annotate(s='', xy=(0, 0), xycoords='data',
                                                                                                   **self.plot_characteristics['periodic']['marker_annotate'],
                                                                                                   bbox={**self.plot_characteristics['periodic']['marker_annotate_bbox']},
                                                                                                   arrowprops={**self.plot_characteristics['periodic']['marker_annotate_arrowprops']})
            self.periodic_annotation[resolution].set_visible(False)

        return None

    def create_periodic_annotation_vline(self):
        """ Create annotation vertical line at (0, 0) that will be updated later. """

        # add vertical line
        self.periodic_vline = dict()
        for resolution in self.plot_axes['periodic'].keys():
            self.periodic_vline[resolution] = self.plot_axes['periodic'][resolution].axvline(0, **self.plot_characteristics['periodic']['marker_annotate_vline'])
            self.periodic_vline[resolution].set_visible(False)

        return None

    def update_periodic_annotation(self, annotation_index, resolution):
        """ Update annotation for each periodic point that is hovered. """
        
        for data_label in self.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.periodic_annotate_data_label:
                
                # skip observations for bias plot
                if self.plot_elements['periodic']['active'] == 'bias' and data_label == 'observations':
                    continue
                
                # do not annotate if plot is cleared
                if data_label not in self.plot_elements['periodic'][self.plot_elements['periodic']['active']].keys():
                    continue

                # retrieve time and concentration
                line = self.plot_elements['periodic'][self.plot_elements['periodic']['active']][data_label]['plot_' + resolution][0]
                time = line.get_xdata()[annotation_index['ind'][0]]
                concentration = line.get_ydata()[annotation_index['ind'][0]]

                # update location
                self.periodic_annotation[resolution].xy = (time, concentration)
                self.periodic_vline[resolution].set_xdata(time)

                # update bbox position
                time_middle = line.get_xdata()[math.floor((len(line.get_xdata()) - 1)/2)]
                if time > time_middle:
                    self.periodic_annotation[resolution].set_x(-10)
                    self.periodic_annotation[resolution].set_ha('right')
                else:
                    self.periodic_annotation[resolution].set_x(10)
                    self.periodic_annotation[resolution].set_ha('left')

                # create annotation text
                if resolution == 'hour':
                    resolution_text = 'Hour'
                    time_text = time
                else:
                    time_options = [self.temporal_axis_mapping_dict['long'][resolution][xtick] 
                                    for xtick in self.periodic_xticks[resolution]]
                    time_text = time_options[time-1]
                    if resolution == 'dayofweek':
                        resolution_text = 'Day'
                    elif resolution == 'month':
                        resolution_text = 'Month'
                text_label = ('{0}: {1}').format(resolution_text, time_text)
        
        for data_label in self.plot_elements['data_labels_active']:
            
            # skip observations for bias plot
            if self.plot_elements['periodic']['active'] == 'bias' and data_label == 'observations':
                continue

            # do not annotate if plot is cleared
            if data_label not in self.plot_elements['periodic'][self.plot_elements['periodic']['active']].keys():
                continue

            # retrieve concentration
            line = self.plot_elements['periodic'][self.plot_elements['periodic']['active']][data_label]['plot_' + resolution][0]
            concentration = line.get_ydata()[np.where(line.get_xdata() == time)[0]]
            
            # for all labels if there is data
            if len(concentration) >= 1:
                if data_label != 'observations':
                    exp_alias = self.read_instance.experiments[data_label]
                    text_label += ('\n{0}: {1:.2f}').format(exp_alias, concentration[0])
                else:
                    text_label += ('\n{0}: {1:.2f}').format(self.plot_characteristics['legend']['handles']['obs_label'], 
                                                            concentration[0])
        
        self.periodic_annotation[resolution].set_text(text_label)

        return None

    def hover_periodic_annotation(self, event):
        """ Show or hide annotation for each point that is hovered in the periodic plot. """
        
        # activate hover over periodic
        if ('periodic' in self.read_instance.active_dashboard_plots):
            if hasattr(self.read_instance, 'relevant_temporal_resolutions'):
                for resolution in self.read_instance.relevant_temporal_resolutions:
                    if event.inaxes == self.plot_axes['periodic'][resolution]:
                        if ((hasattr(self.plot, 'periodic_plot')) and ('periodic' in self.plot_elements)
                            and (self.lock_periodic_annotation[resolution] == False)):

                            # lock annotation
                            self.lock_periodic_annotation[resolution] = True
                            is_contained = False

                            for data_label in self.plot_elements['data_labels_active']:

                                # skip observations for bias plot
                                if self.plot_elements['periodic']['active'] == 'bias' and data_label == 'observations':
                                    continue

                                # do not annotate if plot is cleared
                                if data_label not in self.plot_elements['periodic'][self.plot_elements['periodic']['active']].keys():
                                    continue
                                
                                line = self.plot_elements['periodic'][self.plot_elements['periodic']['active']][data_label]['plot_' + resolution][0]
                                is_contained, annotation_index = line.contains(event)
                                if is_contained:
                                    self.periodic_annotate_data_label = data_label
                                    break
                            
                            if is_contained:
                                # update annotation if hovered
                                self.update_periodic_annotation(annotation_index, resolution)
                                self.periodic_annotation[resolution].set_visible(True)
                                self.periodic_vline[resolution].set_visible(True)
                            else:
                                # hide annotation if not hovered
                                if self.periodic_annotation[resolution].get_visible():
                                    self.periodic_annotation[resolution].set_visible(False)
                                    self.periodic_vline[resolution].set_visible(False)
                                    
                            # redraw points
                            self.figure.canvas.draw()
                            self.figure.canvas.flush_events()
                                
                            # unlock annotation 
                            self.lock_periodic_annotation[resolution] = False

        return None

    def create_periodic_violin_annotation(self):
        """ Create annotation at (0, 0) that will be updated later. """

        # in the newest version of matplotlib, s corresponds to text
        self.periodic_violin_annotation = dict()
        for resolution in self.plot_axes['periodic-violin'].keys():
            self.periodic_violin_annotation[resolution] = self.plot_axes['periodic-violin'][resolution].annotate(s='', xy=(0, 0), xycoords='data',
                                                                                                                 **self.plot_characteristics['periodic-violin']['marker_annotate'],
                                                                                                                 bbox={**self.plot_characteristics['periodic-violin']['marker_annotate_bbox']},
                                                                                                                 arrowprops={**self.plot_characteristics['periodic-violin']['marker_annotate_arrowprops']})
            self.periodic_violin_annotation[resolution].set_visible(False)

        return None

    def create_periodic_violin_annotation_vline(self):
        """ Create annotation vertical line at (0, 0) that will be updated later. """

        # add vertical line
        self.periodic_violin_vline = dict()
        for resolution in self.plot_axes['periodic-violin'].keys():
            self.periodic_violin_vline[resolution] = self.plot_axes['periodic-violin'][resolution].axvline(0, **self.plot_characteristics['periodic-violin']['marker_annotate_vline'])
            self.periodic_violin_vline[resolution].set_visible(False)

        return None

    def update_periodic_violin_annotation(self, annotation_index, resolution):
        """ Update annotation for each periodic violin point that is hovered. """
        
        for data_label in self.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.periodic_violin_annotate_data_label:
                
                # skip observations for bias plot
                if self.plot_elements['periodic-violin']['active'] == 'bias' and data_label == 'observations':
                    continue
                
                # do not annotate if plot is cleared
                if data_label not in self.plot_elements['periodic-violin'][self.plot_elements['periodic-violin']['active']].keys():
                    continue

                # retrieve time and concentration
                line = self.plot_elements['periodic-violin'][self.plot_elements['periodic-violin']['active']][data_label]['p50_plot_' + resolution][0]
                time = line.get_xdata()[annotation_index['ind'][0]]
                concentration = line.get_ydata()[annotation_index['ind'][0]]

                # update location
                self.periodic_violin_annotation[resolution].xy = (time, concentration)
                self.periodic_violin_vline[resolution].set_xdata(time)

                # update bbox position
                time_middle = line.get_xdata()[math.floor((len(line.get_xdata()) - 1)/2)]
                if time > time_middle:
                    self.periodic_violin_annotation[resolution].set_x(-10)
                    self.periodic_violin_annotation[resolution].set_ha('right')
                else:
                    self.periodic_violin_annotation[resolution].set_x(10)
                    self.periodic_violin_annotation[resolution].set_ha('left')

                # create annotation text
                if resolution == 'hour':
                    resolution_text = 'Hour'
                    time_text = time
                else:
                    time_options = [self.temporal_axis_mapping_dict['long'][resolution][xtick] 
                                    for xtick in self.periodic_xticks[resolution]]
                    time_text = time_options[time-1]
                    if resolution == 'dayofweek':
                        resolution_text = 'Day'
                    elif resolution == 'month':
                        resolution_text = 'Month'
                text_label = ('{0}: {1}').format(resolution_text, time_text)
        
        for data_label in self.plot_elements['data_labels_active']:
            
            # skip observations for bias plot
            if self.plot_elements['periodic-violin']['active'] == 'bias' and data_label == 'observations':
                continue

            # do not annotate if plot is cleared
            if data_label not in self.plot_elements['periodic-violin'][self.plot_elements['periodic-violin']['active']].keys():
                continue

            # retrieve concentration
            line = self.plot_elements['periodic-violin'][self.plot_elements['periodic-violin']['active']][data_label]['p50_plot_' + resolution][0]
            concentration = line.get_ydata()[np.where(line.get_xdata() == time)[0]]
            
            # for all labels if there is data
            if len(concentration) >= 1:
                if data_label != 'observations':
                    exp_alias = self.read_instance.experiments[data_label]
                    text_label += ('\n{0}: {1:.2f}').format(exp_alias, concentration[0])
                else:
                    text_label += ('\n{0}: {1:.2f}').format(self.plot_characteristics['legend']['handles']['obs_label'], 
                                                            concentration[0])
        
        self.periodic_violin_annotation[resolution].set_text(text_label)

        return None

    def hover_periodic_violin_annotation(self, event):
        """ Show or hide annotation for each point that is hovered in the periodic violin plot. """
        
        # activate hover over periodic violin
        if ('periodic-violin' in self.read_instance.active_dashboard_plots):
            if hasattr(self.read_instance, 'relevant_temporal_resolutions'):
                for resolution in self.read_instance.relevant_temporal_resolutions:
                    if event.inaxes == self.plot_axes['periodic-violin'][resolution]:
                        if ((hasattr(self.plot, 'violin_plot')) and ('periodic-violin' in self.plot_elements)
                            and (self.lock_periodic_violin_annotation[resolution] == False)):
                            
                            # lock annotation
                            self.lock_periodic_violin_annotation[resolution] = True
                            is_contained = False

                            for data_label in self.plot_elements['data_labels_active']:

                                # skip observations for bias plot
                                if self.plot_elements['periodic-violin']['active'] == 'bias' and data_label == 'observations':
                                    continue

                                # do not annotate if plot is cleared
                                if data_label not in self.plot_elements['periodic-violin'][self.plot_elements['periodic-violin']['active']].keys():
                                    continue
                                
                                line = self.plot_elements['periodic-violin'][self.plot_elements['periodic-violin']['active']][data_label]['p50_plot_' + resolution][0]
                                is_contained, annotation_index = line.contains(event)
                                if is_contained:
                                    self.periodic_violin_annotate_data_label = data_label
                                    break
                            
                            if is_contained:
                                # update annotation if hovered
                                self.update_periodic_violin_annotation(annotation_index, resolution)
                                self.periodic_violin_annotation[resolution].set_visible(True)
                                self.periodic_violin_vline[resolution].set_visible(True)
                            else:
                                # hide annotation if not hovered
                                if self.periodic_violin_annotation[resolution].get_visible():
                                    self.periodic_violin_annotation[resolution].set_visible(False)
                                    self.periodic_violin_vline[resolution].set_visible(False)
                                    
                            # redraw points
                            self.figure.canvas.draw()
                            self.figure.canvas.flush_events()
                                
                            # unlock annotation 
                            self.lock_periodic_violin_annotation[resolution] = False

        return None

    def zoom_map_func(self, event):
        """ Function to handle zoom on map using scroll wheel. """

        if event.inaxes == self.plot_axes['map']:
            
            if self.lock_zoom == False:

                # lock zoom
                self.lock_zoom = True

                # get the current x and y limits
                current_xlim = self.plot_axes['map'].get_xlim()
                current_ylim = self.plot_axes['map'].get_ylim()

                # get position of cursor
                xdata = event.xdata
                ydata = event.ydata
                base_scale = self.plot_characteristics['map']['base_scale']

                if event.button == 'up':
                    # deal with zoom in
                    scale_factor = base_scale
                elif event.button == 'down':
                    # deal with zoom out
                    scale_factor = 1/base_scale
                else:
                    # exceptions
                    scale_factor = 1
                
                if event.button == 'up' or event.button == 'down':
                    
                    # set new limits
                    self.plot_axes['map'].set_xlim([xdata - (xdata - current_xlim[0]) / scale_factor, 
                                                    xdata + (current_xlim[1] - xdata) / scale_factor])
                    self.plot_axes['map'].set_ylim([ydata - (ydata - current_ylim[0]) / scale_factor, 
                                                    ydata + (current_ylim[1] - ydata) / scale_factor])
                    
                    # save map extent (in data coords)
                    self.read_instance.map_extent = self.plot.get_map_extent(self.plot_axes['map'])
                    
                    # redraw points
                    self.figure.canvas.draw()
                    self.figure.canvas.flush_events()
                
                    # update buttons (previous-forward) history
                    self.read_instance.navi_toolbar.push_current()
                    self.read_instance.navi_toolbar.set_history_buttons()

                    # unlock zoom
                    self.lock_zoom = False

        return None

    def picker_block_func(self, event):
        """ Block or unblock the station and legend pick functions
            to avoid interferences.
        """

        if event.inaxes == self.plot_axes['legend']:
            # unblock legend picker in legend
            self.lock_legend_pick = False
        
        else:
            # block legend picker
            self.lock_legend_pick = True

        return None

    def legend_picker_func(self, event):
        """ Function to handle legend picker. """

        if self.lock_legend_pick == False:
            if self.plot_elements:
                # lock legend pick
                self.lock_legend_pick = True
            
                # get legend label information
                legend_label = event.artist
                data_label = legend_label.get_text()

                # transform legend label into data labels
                for exp_label, exp_alias in self.read_instance.experiments.items():
                    if data_label == exp_alias:
                        data_label = exp_label
                        continue
                if data_label == self.plot_characteristics['legend']['handles']['obs_label']:
                    data_label = 'observations'

                if data_label not in self.plot_elements['data_labels_active']:
                    visible = True
                    # put observations label always first in pop-ups on hover
                    if data_label == 'observations':
                        self.plot_elements['data_labels_active'].insert(0, data_label)
                    # put experiment labels in the same order as in the legend
                    else:
                        self.plot_elements['data_labels_active'].insert(list(self.read_instance.experiments.keys()).index(data_label)+1, 
                                                                        data_label)
                else:
                    visible = False
                    self.plot_elements['data_labels_active'].remove(data_label)

                # iterate through plot types stored in plot_elements (if have selected stations)
                if len(self.relative_selected_station_inds) > 0:
                    for plot_type_alt in self.plot_elements:  

                        # correct perodic-violin name for plot_options
                        if plot_type_alt == 'periodic-violin':
                            plot_type = 'periodic_violin'
                        else:
                            plot_type = copy.deepcopy(plot_type_alt)

                        if plot_type not in ['data_labels_active', 'metadata', 'map', 'heatmap', 
                                             'table', 'statsummary']:

                            # get currently selected options for plot
                            plot_options = self.read_instance.current_plot_options[plot_type]
                        
                            # get active (absolute / bias)
                            active = self.plot_elements[plot_type]['active']

                            # change visibility of plot elements (if data label in plot elements dictionary)
                            if data_label in self.plot_elements[plot_type][active]:
                                for element_type in self.plot_elements[plot_type][active][data_label]:
                                    for plot_element in self.plot_elements[plot_type][active][data_label][element_type]:
                                        if visible:
                                            plot_element.set_visible(True)
                                        else:
                                            plot_element.set_visible(False)

                            # reset axes limits (harmonising across subplots for periodic plots) 
                            if plot_type == 'periodic-violin':
                                self.plot.harmonise_xy_lims_paradigm(self.plot_axes[plot_type], plot_type, 
                                                                     self.plot_characteristics[plot_type], 
                                                                     plot_options, 
                                                                     ylim=[self.selected_station_data_min[self.read_instance.networkspeci], 
                                                                           self.selected_station_data_max[self.read_instance.networkspeci]],
                                                                     relim=True, autoscale_x=True)
                            elif plot_type == 'scatter':
                                self.plot.harmonise_xy_lims_paradigm(self.plot_axes[plot_type], plot_type, 
                                                                     self.plot_characteristics[plot_type], 
                                                                     plot_options, relim=True)
                            else: 
                                self.plot.harmonise_xy_lims_paradigm(self.plot_axes[plot_type], plot_type, 
                                                                     self.plot_characteristics[plot_type], 
                                                                     plot_options, relim=True, autoscale=True)

                # change font weight of label
                legend_label._fontproperties = self.legend.get_texts()[0]._fontproperties.copy()
                if visible:
                    legend_label.set_fontweight('bold')
                else:
                    legend_label.set_fontweight('regular')

                # redraw points
                self.figure.canvas.draw()
                self.figure.canvas.flush_events()
                
                # unlock legend pick 
                self.lock_legend_pick = False
    
        return None