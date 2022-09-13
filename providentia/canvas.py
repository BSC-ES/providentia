from .read_aux import drop_nans
from .filter import DataFilter
from .statistics import to_pandas_dataframe
from .statistics import calculate_z_statistic
from .statistics import generate_colourbar
from .statistics import get_z_statistic_info
from .statistics import get_z_statistic_comboboxes
from .plot import Plot

import copy
import json
import os
import sys
import time
from weakref import WeakKeyDictionary

import matplotlib
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg \
        as FigureCanvas

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector, Slider
import matplotlib.gridspec as gridspec
from pandas.plotting import register_matplotlib_converters
from PyQt5 import QtCore, QtGui, QtWidgets 

import numpy as np
import pandas as pd

# make sure that we are using Qt5 backend with matplotlib
matplotlib.use('Qt5Agg')
register_matplotlib_converters()

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))

class MPLCanvas(FigureCanvas):
    """Class that handles the creation and updates of
    a matplotlib canvas, and associated subplots
    """

    def __init__(self, read_instance):
        """Initialise the MPL canvas"""

        # create figure and canvas objects
        self.figure = Figure(dpi=100)
        FigureCanvas.__init__(self, self.figure)

        # add passed arguments to self
        self.read_instance = read_instance

        # get characteristics per plot type
        self.plot_characteristics_templates = self.read_instance.plot_characteristics_templates
        self.plot_characteristics = {}

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

        # create plot axes on grid
        self.plot_axes = {}
        self.plot_axes['map'] = self.figure.add_subplot(self.gridspec.new_subplotspec((2, 0), rowspan=44, colspan=42), projection=self.plotcrs)
        self.plot_axes['legend'] = self.figure.add_subplot(self.gridspec.new_subplotspec((0, 47), rowspan=8,  colspan=53))
        self.plot_axes['timeseries'] = self.figure.add_subplot(self.gridspec.new_subplotspec((12, 50), rowspan=34, colspan=50))
        self.plot_axes['periodic-violin'] = {}
        self.plot_axes['periodic-violin']['hour'] = self.figure.add_subplot(self.gridspec.new_subplotspec((56, 70), rowspan=20, colspan=30))
        self.plot_axes['periodic-violin']['dayofweek'] = self.figure.add_subplot(self.gridspec.new_subplotspec((80, 90), rowspan=20, colspan=10))
        self.plot_axes['periodic-violin']['month'] = self.figure.add_subplot(self.gridspec.new_subplotspec((80, 70), rowspan=20, colspan=18))
        self.plot_axes['distribution'] = self.figure.add_subplot(self.gridspec.new_subplotspec((56, 35), rowspan=44, colspan=30))
        self.plot_axes['metadata'] = self.figure.add_subplot(self.gridspec.new_subplotspec((56, 0),  rowspan=44, colspan=30))
        self.plot_axes['cb'] = self.figure.add_axes([0.0455, 0.536, 0.3794, 0.02])

        # define plots to update upon selecting stations
        self.selected_station_plots = ['timeseries', 'periodic-violin', 'distribution', 'metadata']

        # define all possible plots
        self.all_plots = ['legend', 'map', 'timeseries', 'periodic-violin', 'periodic', 
                          'metadata', 'distribution', 'scatter']

        # update plot characteristics for all plots
        for plot_type in self.all_plots:
            # for plot types with zstat, initialise with default zstat (Mean)
            if plot_type in ['map', 'periodic']:
                self.plot.set_plot_characteristics([plot_type], zstat='Mean')
            else:
                self.plot.set_plot_characteristics([plot_type])

        # initialise variable of valid station indices plotted on map as empty list
        self.active_map_valid_station_inds = np.array([], dtype=np.int)

        # add settings menus
        self.generate_interactive_elements()

        # setup blocker for picker events
        self.figure.canvas.mpl_connect('axes_enter_event', self.picker_block_func)

        # setup legend line selection
        self.figure.canvas.mpl_connect('pick_event', self.legend_picker_func)

        # setup interactive lasso on map
        self.lasso_left = LassoSelector(self.plot_axes['map'], onselect=self.onlassoselect_left,
                                   useblit=True, lineprops=self.lasso, button=[1])
        self.lasso_right = LassoSelector(self.plot_axes['map'], onselect=self.onlassoselect_right,
                                   useblit=True, lineprops=self.lasso, button=[3])

        # setup station annotations
        self.create_station_annotation()
        self.annotation_visible = False
        self.lock_annotation = False
        self.station_annotation_event = self.figure.canvas.mpl_connect('motion_notify_event', self.hover_station_annotation)

        # setup zoom on scroll wheel on map
        self.lock_zoom = False
        self.figure.canvas.mpl_connect('scroll_event', self.zoom_map_func)

        # format axes for map and legend
        for plot_type in ['map', 'legend']:
            ax = self.plot_axes[plot_type]
            if type(ax) == dict:
                for relevant_temporal_resolution, sub_ax in ax.items():
                    self.plot.format_axis(sub_ax, plot_type, self.plot_characteristics[plot_type], 
                                          relevant_temporal_resolution=relevant_temporal_resolution, 
                                          col_ii=-1)
            else:
                self.plot.format_axis(ax, plot_type, self.plot_characteristics[plot_type])

        # hide all plot axes
        for plot_type, ax in self.plot_axes.items():
            if type(ax) == dict:
                for relevant_temporal_resolution, sub_ax in ax.items():
                    self.remove_axis_elements(sub_ax, plot_type)
            else:
                self.remove_axis_elements(ax, plot_type)

    def update_MPL_canvas(self):
        """Function that updates MPL canvas upon clicking
        the 'READ' button, and when colocating data
        """

        print('UPDATE CANVAS')

        # clear and then hide all axes 
        for plot_type, ax in self.plot_axes.items():
            if type(ax) == dict:
                for sub_ax in ax.values():
                    self.remove_axis_elements(sub_ax, plot_type)
            else:
                self.remove_axis_elements(ax, plot_type)

        # reset relative index lists of selected station on map as empty lists
        self.previous_relative_selected_station_inds = np.array([], dtype=np.int)
        self.relative_selected_station_inds = np.array([], dtype=np.int)
        self.absolute_selected_station_inds = np.array([], dtype=np.int)

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
        """Function which resets the navigation toolbar stack
        for a given axis with the current view limits"""

        # check if have axes dictionaries in stack list
        if len(self.read_instance.navi_toolbar._nav_stack) == 0:
            # if don't have an axes dictionary in stack list, create one with current
            # axis in dictionary with current view limits
            self.read_instance.navi_toolbar._nav_stack.push(
                WeakKeyDictionary({ax: (ax._get_view(), (ax.get_position(True).frozen(), ax.get_position().frozen()))}))

        # if have existing axes dictionaries in stack list, iterate through stack list
        # removing given axis from all stack list dictionaries
        else:
            for axes_dict in self.read_instance.navi_toolbar._nav_stack:
                if ax in axes_dict.keyrefs():
                    axes_dict.pop(ax)

            # now add axis to first dictionary in stack, with the current view limits
            self.read_instance.navi_toolbar._nav_stack[0][ax] = \
                (ax._get_view(), (ax.get_position(True).frozen(), ax.get_position().frozen()))

        return None

    def handle_data_filter_update(self):
        """Function which handles updates of data filtering"""

        print('UPDATE DATA FILTER')

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

    def update_active_map(self):
        """Function that updates plotted map z statistic and updates associated plots"""

        if not self.read_instance.block_MPL_canvas_updates:

            print('UPDATE ACTIVE MAP')

            # update plotted map z statistic
            self.update_map_z_statistic()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds,
                                  self.relative_selected_station_inds):
                self.update_associated_selected_station_plots()

            # draw changes
            self.figure.canvas.draw()

        return None

    def handle_temporal_colocate_update(self):
        """Function that handles the update of the MPL canvas
        with colocated data upon checking of the temporal colocate checkbox"""

        if not self.read_instance.block_MPL_canvas_updates:

            print('TEMPORAL COLOCATE UPDATE')

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

            # update z statistic/experiment bias comboboxes (without updating canvas)
            self.read_instance.block_MPL_canvas_updates = True
            self.handle_map_z_statistic_update()
            self.handle_experiment_bias_update()
            self.read_instance.block_MPL_canvas_updates = False

            # update plotted map z statistic
            self.update_map_z_statistic()

            # update layout field options 
            self.read_instance.update_layout_fields()

            # update associated plots with selected stations
            self.update_associated_selected_station_plots()

            # draw changes
            self.figure.canvas.draw()

        return None

    def update_map_z_statistic(self):
        """Function that updates plotted z statistic on map, with colourbar"""

        print('UPDATE MAP Z')

        # remove axis elements from map/cb
        self.remove_axis_elements(self.plot_axes['map'], 'map')
        self.remove_axis_elements(self.plot_axes['cb'], 'cb')

        # get zstat name from combobox 
        base_zstat = self.read_instance.cb_z_stat.currentText()
        zstat = get_z_statistic_comboboxes(base_zstat, second_data_label=self.read_instance.cb_z2.currentText())

        # calculate map z statistic (for selected z statistic) --> updating active map valid station indices
        z_statistic, active_map_valid_station_inds = calculate_z_statistic(self.read_instance, self.read_instance.cb_z1.currentText(), self.read_instance.cb_z2.currentText(), zstat, self.read_instance.networkspeci)
        self.active_map_valid_station_inds = active_map_valid_station_inds 
        
        # update absolute selected plotted station indices with respect to new active map valid station indices
        self.absolute_selected_station_inds = np.array(
            [np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind in
             self.relative_selected_station_inds if selected_ind in self.active_map_valid_station_inds],
            dtype=np.int)

        # if have no valid active map indices, reset absolute/relative
        # selected station indices to be empty lists
        # also uncheck select all/intersect checkboxes
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

        # otherwise plot valid active map stations on map
        else:
            # if any of the currently selected stations are not in the current active map
            # valid station indices --> unselect selected stations (and associated plots)
            # also uncheck select all/intersect checkboxes
            if not np.all(np.in1d(self.relative_selected_station_inds, self.active_map_valid_station_inds)):
                # unselect all/intersect checkboxes
                self.read_instance.block_MPL_canvas_updates = True
                self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
                self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
                self.read_instance.block_MPL_canvas_updates = False
                # reset relative/absolute selected station indices to be empty lists
                self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)
                self.relative_selected_station_inds = np.array([], dtype=np.int)
                self.absolute_selected_station_inds = np.array([], dtype=np.int)
                self.absolute_non_selected_station_inds = np.array([], dtype=np.int)

            # plot new station points on map - coloured by currently active z statisitic, setting up plot picker
            self.plot.make_map(self.plot_axes['map'], self.read_instance.networkspeci, z_statistic, self.plot_characteristics['map'])

            # create 2D numpy array of plotted station coordinates
            self.map_points_coordinates = np.vstack((self.read_instance.station_longitudes[self.read_instance.networkspeci][self.active_map_valid_station_inds], 
                                                     self.read_instance.station_latitudes[self.read_instance.networkspeci][self.active_map_valid_station_inds])).T

            # generate colourbar
            generate_colourbar(self.read_instance, [self.plot_axes['map']], [self.plot_axes['cb']], zstat, self.plot_characteristics['map'], self.read_instance.species[0])

            # activate map/cb axes
            self.activate_axis(self.plot_axes['map'], 'map')
            self.activate_axis(self.plot_axes['cb'], 'cb')

        # update plot options
        self.update_plot_options()

        # re-draw (needed to update plotted colours before update_map_station_selection)
        self.figure.canvas.draw()

        # update map selection appropriately for z statistic
        self.update_map_station_selection()

        return None

    def update_map_station_selection(self):
        """Function that updates the visual selection of stations on map"""

        print('UPDATE STATION SELECTION')

        # update map title
        if len(self.relative_selected_station_inds) == 1:
            axis_title_label = '{} Selected'.format(
                self.read_instance.station_references[self.read_instance.networkspeci][self.relative_selected_station_inds][0])
        else:
            axis_title_label = '{} Selected Stations of {} Available'.format(
                len(self.relative_selected_station_inds), len(self.active_map_valid_station_inds))
        self.plot.set_axis_title(self.plot_axes['map'], axis_title_label, self.plot_characteristics['map'])

        # reset alphas and marker sizes of all plotted stations (if have some stations on map)
        if len(self.active_map_valid_station_inds) > 0:
            markersizes = np.full(len(self.active_map_valid_station_inds), self.plot_characteristics['map']['marker_unselected']['s'])
            for collection in self.plot_axes['map'].collections:
                if isinstance(collection, matplotlib.collections.PathCollection):
                    rgba_tuples = collection.get_facecolor()
                    rgba_tuples[:, -1] = self.plot_characteristics['map']['marker_selected']['alpha']
                    collection.set_facecolor(rgba_tuples)
                    collection.set_sizes(markersizes)

            # if have some selected stations, update map plot with station selection
            # (reducing alpha of non-selected stations, and increasing marker size of selected stations)
            if len(self.relative_selected_station_inds) > 0:

                self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                                     self.absolute_selected_station_inds))[0]
                
                # decrease alpha of non-selected stations
                if len(self.absolute_non_selected_station_inds) > 0:
                    rgba_tuples[self.absolute_non_selected_station_inds, -1] = self.plot_characteristics['map']['marker_unselected']['alpha']
                    for collection in self.plot_axes['map'].collections:
                        if isinstance(collection, matplotlib.collections.PathCollection):
                            collection.set_facecolor(rgba_tuples)

                # increase marker size of selected stations
                markersizes[self.absolute_selected_station_inds] = self.plot_characteristics['map']['marker_selected']['s']
                for collection in self.plot_axes['map'].collections:
                    if isinstance(collection, matplotlib.collections.PathCollection):
                        collection.set_sizes(markersizes)

        # redraw points
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def update_associated_selected_station_plot(self, plot_type):
        """Function that updates a plot associated with selected stations on map"""
    
        print('UPDATE SELECTED STATION PLOT')

        if hasattr(self, 'relative_selected_station_inds'):
            if len(self.relative_selected_station_inds) > 0:

                # get relevant axis
                ax = self.plot_axes[plot_type]

                # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                plot_options = plot_type.split('_')[1:]

                # format axis
                if type(ax) == dict:
                    for relevant_temporal_resolution, sub_ax in ax.items():
                        self.plot.format_axis(sub_ax, plot_type, self.plot_characteristics[plot_type], 
                                            relevant_temporal_resolution=relevant_temporal_resolution, 
                                            col_ii=-1)
                else:
                    self.plot.format_axis(ax, plot_type, self.plot_characteristics[plot_type])

                # get plotting function for specific plot
                func = getattr(self.plot, 'make_{}'.format(plot_type.split('-')[0]))

                # get zstat for plots with statistics, and set new axis title and ylabel
                if plot_type in ['periodic']:
                    # get currently selected experiment bias statistic name
                    base_zstat = self.read_instance.cb_experiment_bias_stat.currentText()
                    # if have no zstat, it is because no experiments are loaded so cannot make periodic bias plot
                    if not base_zstat:
                        return 
                    zstat = get_z_statistic_comboboxes(base_zstat, second_data_label='model')
                    # get zstat information 
                    zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat) 
                    # add bias to plot options if stat is bias
                    if z_statistic_sign == 'bias':
                        plot_options.append('bias')  
                    #set new axis title
                    axis_title = '{}-{}'.format(plot_type, zstat) 
                    # set new ylabel
                    if z_statistic_type == 'basic':
                        if base_zstat not in ['Data%', 'Exceedances']:
                            ylabel = self.read_instance.measurement_units[self.read_instance.species[0]] 
                        else:
                            ylabel = basic_stats[base_zstat]['label']
                    else:
                        ylabel = expbias_stats[base_zstat]['label']
                else:
                    #set new axis title
                    axis_title = plot_type 
                    
                    # set new ylabel
                    if 'ylabel' in self.plot_characteristics[plot_type]:
                        if self.plot_characteristics[plot_type]['ylabel']['ylabel'] == 'measurement_units':
                            ylabel = self.read_instance.measurement_units[self.read_instance.species[0]]
                        else:
                            ylabel = self.plot_characteristics[plot_type]['ylabel']['ylabel']
                    else:
                        ylabel = ''

                    # set new xlabel
                    if 'xlabel' in self.plot_characteristics[plot_type]:
                        if self.plot_characteristics[plot_type]['xlabel']['xlabel'] == 'measurement_units':
                            xlabel = self.read_instance.measurement_units[self.read_instance.species[0]]
                        else:
                            xlabel = self.plot_characteristics[plot_type]['xlabel']['xlabel']
                    else:
                        xlabel = ''

                # if are making bias plot, and have no valid experiment data then cannot make plot type
                if ('bias' in plot_options) & (len(self.selected_station_data[self.read_instance.networkspeci]) < 2):
                    return

                # iterate through data array names in selected station data dictionary
                valid_data_labels = []
                for data_label in self.selected_station_data[self.read_instance.networkspeci]:

                    # call function to update plot
                    if plot_type in ['periodic']:
                        # skip observational array if bias stat and data array is observations
                        if (z_statistic_sign == 'bias') & (data_label == 'observations'):
                            continue
                        func(ax, self.read_instance.networkspeci, data_label, self.plot_characteristics[plot_type], zstat=zstat, plot_options=plot_options)
                    else: 
                        if plot_type == 'metadata':
                            if data_label != 'observations':
                                continue
                        if plot_type == 'scatter':
                            if data_label == 'observations':
                                continue
                        func(ax, self.read_instance.networkspeci, data_label, self.plot_characteristics[plot_type], plot_options=plot_options)
                    valid_data_labels.append(data_label)
                
                # format axes reset axes limits (harmonise across subplots for periodic plots), reset navigation toolbar stack, and set axis title / ylabel
                if type(ax) == dict:
                    relevant_axs = [ax[relevant_temporal_resolution] for relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions]
                    if plot_type == 'periodic-violin':
                        self.plot.harmonise_xy_lims_paradigm(relevant_axs, plot_type, self.plot_characteristics[plot_type], 
                                                             plot_options, ylim=[self.selected_station_data_min[self.read_instance.networkspeci], 
                                                             self.selected_station_data_max[self.read_instance.networkspeci]])
                    else:
                        self.plot.harmonise_xy_lims_paradigm(relevant_axs, plot_type, self.plot_characteristics[plot_type], 
                                                             plot_options, autoscale_y=True)

                    set_title = False
                    for relevant_temporal_resolution, sub_ax in ax.items():
                        #do not show axis if temporal resolution is not relevant
                        if relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions:
                            # set ylabel (if sub_ax in first column) and title (if sub_ax in first column and not previously set)   
                            if relevant_temporal_resolution in ['hour','month']:
                                self.plot.set_axis_label(sub_ax, 'y', ylabel, self.plot_characteristics[plot_type])
                                if not set_title:
                                    self.plot.set_axis_title(sub_ax, axis_title, self.plot_characteristics[plot_type])
                                    set_title = True
                            self.activate_axis(sub_ax, plot_type)
                            self.reset_ax_navigation_toolbar_stack(sub_ax)
                else:
                    if plot_type == 'scatter':
                        self.plot.harmonise_xy_lims_paradigm([ax], plot_type, self.plot_characteristics[plot_type], 
                                                             plot_options, relim=True)
                    else: 
                        self.plot.harmonise_xy_lims_paradigm([ax], plot_type, self.plot_characteristics[plot_type], 
                                                             plot_options, relim=True, autoscale=True)
                    self.plot.set_axis_title(ax, axis_title, self.plot_characteristics[plot_type])
                    self.plot.set_axis_label(ax, 'x', xlabel, self.plot_characteristics[plot_type])
                    self.plot.set_axis_label(ax, 'y', ylabel, self.plot_characteristics[plot_type])
                    self.activate_axis(ax, plot_type)
                    self.reset_ax_navigation_toolbar_stack(ax)

                # configure legend picker per plot type / data_label
                self.configure_legend_picker(plot_type, valid_data_labels)

    def update_associated_selected_station_plots(self):
        """Function that updates all plots associated with selected stations on map"""

        print('UPDATE SELECTED STATION PLOTS')

        # clear all previously plotted artists from selected station plots and hide axes 
        for plot_type in self.selected_station_plots:
            ax = self.plot_axes[plot_type]
            if type(ax) == dict:
                for sub_ax in ax.values():
                    self.remove_axis_elements(sub_ax, plot_type)
            else:
                self.remove_axis_elements(ax, plot_type)

        # remove legend picker lines
        for legend_label in self.lined:
            relevant_plots = list(self.lined[legend_label]['lines_per_plot'].keys())
            for plot_type in relevant_plots:
                del self.lined[legend_label]['lines_per_plot'][plot_type]

        # if have selected stations on map, then now remake plots
        if hasattr(self, 'relative_selected_station_inds'):
            if len(self.relative_selected_station_inds) > 0:

                # put selected data for each data array into pandas dataframes
                to_pandas_dataframe(read_instance=self.read_instance, canvas_instance=self, networkspecies=[self.read_instance.networkspeci])

                # iterate through selected_station_plots
                for plot_type in self.selected_station_plots:

                    # update plot
                    self.update_associated_selected_station_plot(plot_type)

    def configure_legend_picker(self, plot_type, legend_labels):
        """ Function that creates a dictionary with the lines and collections
            that will be removed when picking up each legend element
        """

        # get plot lines from timeseries and periodic violin plots
        ax_plot_lines = {}
        ax_plot_collections = {}
        if plot_type in ['timeseries', 'periodic-violin', 'periodic', 'distribution', 'scatter']:
            
                for legend_label in legend_labels:
                    self.lined[legend_label]['lines_per_plot'][plot_type] = []
                
                ax = self.plot_axes[plot_type]
                if type(ax) == dict:
                    for key, sub_ax in ax.items(): 
                        if sub_ax.get_visible():

                            for i, legend_label in enumerate(legend_labels):

                                visible = self.lined[legend_label]['visible']

                                ax_plot_lines[key] = {}
                                ax_plot_collections[key] = {}
                                last_break = 0
                                break_sum = 0

                                # remove [0, 0] lines (periodic plots)
                                sub_ax_lines = [line for line in sub_ax.lines if line.get_ydata() != [0, 0]]
                                
                                # get line per legend element
                                ax_plot_lines[key][legend_label] = sub_ax_lines[i]
                                
                                # are there collections? (skip periodic plots)
                                if sub_ax.collections:
                                
                                    # get corresponding patches according to breaks
                                    # some patches do not exist, so breaks are not always 0-24, 24-48, etc...
                                    # and change depending on the legend element
                                    if i != 0:
                                        last_break = len(sub_ax_lines[i-1].get_ydata()[~np.isnan(sub_ax_lines[i-1].get_ydata())])
                                    curr_break = len(sub_ax_lines[i].get_ydata()[~np.isnan(sub_ax_lines[i].get_ydata())])
                                    break_sum += last_break

                                    # get patches per legend elements
                                    ax_plot_collections[key][legend_label] = sub_ax.collections[break_sum:break_sum+curr_break]

                                # transform single elements to list
                                if not isinstance(ax_plot_lines[key][legend_label], list):
                                    ax_plot_lines[key][legend_label] = [ax_plot_lines[key][legend_label]]

                                # add lines to self.lined
                                self.lined[legend_label]['lines_per_plot'][plot_type] += ax_plot_lines[key][legend_label]

                                # there are no collections for periodic bias plots (only lines)
                                if sub_ax.collections:

                                    # transform single elements to list
                                    if not isinstance(ax_plot_collections[key][legend_label], list):
                                        ax_plot_collections[key][legend_label] = [ax_plot_collections[key][legend_label]]
                                    
                                    # add collections to self.lined
                                    self.lined[legend_label]['lines_per_plot'][plot_type] += ax_plot_collections[key][legend_label]

                                # set line visibility
                                if not visible:
                                    for line in self.lined[legend_label]['lines_per_plot'][plot_type]:
                                        line.set_visible(False)

                else:
                    if ax.get_visible():

                        for i, legend_label in enumerate(legend_labels):

                            visible = self.lined[legend_label]['visible']

                            ax_plot_lines[legend_label] = ax.lines[i]

                            # transform single elements to list
                            if not isinstance(ax_plot_lines[legend_label], list):
                                ax_plot_lines[legend_label] = [ax_plot_lines[legend_label]]

                            # add lines to self.lined
                            self.lined[legend_label]['lines_per_plot'][plot_type] += ax_plot_lines[legend_label]

                            # set line visibility
                            if not visible:
                                for line in self.lined[legend_label]['lines_per_plot'][plot_type]:
                                    line.set_visible(False)

        return None

    def update_experiment_domain_edges(self):
        """Function that plots grid domain edges of experiments in memory"""

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
        """Function that updates legend"""

        # create legend element handles
        legend_plot_characteristics = self.plot.make_legend_handles(copy.deepcopy(self.plot_characteristics['legend']))

        # plot legend
        self.legend = self.plot_axes['legend'].legend(**legend_plot_characteristics['plot'])

        # un-hide legend
        self.activate_axis(self.plot_axes['legend'], 'legend')

        # setup element picker in legend
        self.lined = {}
        for data_label, legend_label in zip(self.read_instance.data_labels, self.legend.texts):
            legend_label.set_picker(True)
            self.lined[data_label] = {'visible':True, 'lines_per_plot':{}}

        return None

    def handle_map_z_statistic_update(self):
        """Define function which handles update of map z statistic upon interaction with map comboboxes"""

        if not self.read_instance.block_config_bar_handling_updates:

            print('MAP Z COMBOBOX UPDATE')

            # update map z statistic comboboxes
            # set variable that blocks configuration bar handling updates until all
            # changes to the z statistic comboboxes are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected items
            selected_z_stat = self.read_instance.cb_z_stat.currentText()
            selected_z1_array = self.read_instance.cb_z1.currentText()
            selected_z2_array = self.read_instance.cb_z2.currentText()

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
            self.read_instance.cb_z_stat.clear()
            self.read_instance.cb_z1.clear()
            self.read_instance.cb_z2.clear()
            self.read_instance.cb_z_stat.addItems(z_stat_items)
            self.read_instance.cb_z1.addItems(z1_items)
            self.read_instance.cb_z2.addItems(z2_items)
            # maintain currently selected z statistic (if exists in new item list)
            if selected_z_stat in z_stat_items:
                self.read_instance.cb_z_stat.setCurrentText(selected_z_stat)
            # maintain currently selected z1/z2 arrays
            self.read_instance.cb_z1.setCurrentText(selected_z1_array)
            self.read_instance.cb_z2.setCurrentText(selected_z2_array)

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            if not self.read_instance.block_MPL_canvas_updates:

                # update plotted map z statistic
                self.update_map_z_statistic()

        return None

    def handle_experiment_bias_update(self):
        """Define function that handles update of plotted experiment bias statistics"""

        if not self.read_instance.block_config_bar_handling_updates:

            print('UPDATE EXP BIAS')
            # update experiment bias comboboxes

            # set variable that blocks configuration bar handling updates until all changes
            # to the experiment bias comboboxes are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected items
            selected_experiment_bias_type = self.read_instance.cb_experiment_bias_type.currentText()
            base_zstat = self.read_instance.cb_experiment_bias_stat.currentText()

            # update experiment bias statistics, to all basic stats
            # if colocation not-active, and basic+bias stats if colocation active
            if not self.read_instance.temporal_colocation:
                available_experiment_bias_stats = copy.deepcopy(self.read_instance.basic_z_stats)
            else:
                available_experiment_bias_stats = copy.deepcopy(self.read_instance.basic_and_bias_z_stats)

            # if base_zstat is empty string, it is because fields are being initialised for the first time
            if base_zstat == '':
                # set experiment bias stat to be first available stat
                base_zstat = available_experiment_bias_stats[0]
                if hasattr(self.read_instance, 'exp_bias_stat'):
                    if self.read_instance.exp_bias_stat in available_experiment_bias_stats:
                        base_zstat = self.read_instance.exp_bias_stat

            # update statistic combobox (clear, then add items)
            self.read_instance.cb_experiment_bias_stat.clear()
            self.read_instance.cb_experiment_bias_stat.addItems(available_experiment_bias_stats)

            # maintain currently selected bias statistic (if exists in new item list)
            if base_zstat in available_experiment_bias_stats:
                self.read_instance.cb_experiment_bias_stat.setCurrentText(base_zstat)

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            if not self.read_instance.block_MPL_canvas_updates:

                # update experiment bias plot/s if have some stations selected on map
                if len(self.relative_selected_station_inds) > 0:

                    # clear and turn off all relevant axes before updating
                    for sub_ax in self.plot_axes['periodic'].values():
                        self.remove_axis_elements(sub_ax, 'periodic')

                    # remove legend picker lines
                    for legend_label in self.lined:
                        if 'periodic' in self.lined[legend_label]['lines_per_plot']:
                            del self.lined[legend_label]['lines_per_plot']['periodic']

                    # get currently selected experiment bias statistic name
                    zstat = get_z_statistic_comboboxes(base_zstat, second_data_label='model')
                    # get zstat information 
                    zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat)
                    # add bias to plot options if neccessary
                    plot_options = []
                    if '_bias' in zstat:
                        plot_options.append('bias')  

                    # set new axis title
                    axis_title = 'periodic-{}'.format(zstat)
                    # set new ylabel
                    if z_statistic_type == 'basic':
                        if base_zstat not in ['Data%', 'Exceedances']:
                            ylabel = self.read_instance.measurement_units[self.read_instance.species[0]] 
                        else:
                            ylabel = basic_stats[base_zstat]['label']
                    else:
                        ylabel = expbias_stats[base_zstat]['label']

                    # if experiment bias type == 'Station Median' --> update plotted experiment bias plots
                    if selected_experiment_bias_type == 'Station Median':

                        valid_data_labels = []
                        for data_label in self.selected_station_data[self.read_instance.networkspeci]:
                            # skip observational array if bias stat
                            if (z_statistic_sign == 'bias') & (data_label == 'observations'):
                                continue
                            self.plot.make_periodic(self.plot_axes['periodic'], self.read_instance.networkspeci, data_label, self.plot_characteristics['periodic'], zstat=zstat)
                            valid_data_labels.append(data_label)

                    # harmonise axes limits across subplots 
                    relevant_axs = [self.plot_axes['periodic'][relevant_temporal_resolution] for relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions]
                    self.plot.harmonise_xy_lims_paradigm(relevant_axs, 'periodic', self.plot_characteristics['periodic'], 
                                                         plot_options, autoscale_y=True)
                    set_title = False
                    # un-hide axes, and reset navigation toolbar stack, and set axis title and ylabel
                    for relevant_temporal_resolution, sub_ax in self.plot_axes['periodic'].items():
                        # do not show axis if temporal resolution is not relevant
                        if relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions:
                            # set ylabel (if sub_ax in first column) and title (if sub_ax in first column and not previously set)   
                            if relevant_temporal_resolution in ['hour','month']:
                                self.plot.set_axis_label(sub_ax, 'y', ylabel, self.plot_characteristics['periodic'])
                                if not set_title:
                                    self.plot.set_axis_title(sub_ax, axis_title, self.plot_characteristics['periodic'])
                                    set_title = True
                            self.activate_axis(sub_ax, 'periodic')
                            self.reset_ax_navigation_toolbar_stack(sub_ax)

                    # configure legend picker per plot type / data_label
                    self.configure_legend_picker('periodic', valid_data_labels)

                    # update plot options
                    self.update_plot_options()

                    # draw changes
                    self.figure.canvas.draw()

        return None

    def remove_axis_elements(self, ax, plot_type):
        """ Removes all plotted axis elements,
            and hide axis.
        """

        #remove all plotted axis elements
        if plot_type == 'legend':
            leg = ax.get_legend()
            if leg:
                leg.remove()

        elif plot_type == 'map':
            inds_to_remove = []
            for artist_ii, artist in enumerate(ax.artists):         
                if type(artist) == AnchoredOffsetbox:
                    inds_to_remove.append(artist_ii)
            ax.artists = list(np.delete(np.array(ax.artists), inds_to_remove))

            inds_to_remove = []
            for col_ii, col in enumerate(ax.collections):            
                if isinstance(col, matplotlib.collections.PathCollection):
                    inds_to_remove.append(col_ii)
            ax.collections = list(np.delete(np.array(ax.collections), inds_to_remove))
            self.map_menu_button.hide()

        elif plot_type == 'cb':
            ax.artists = []
            ax.collections = [] 

        elif plot_type == 'timeseries':
            ax.lines = []
            ax.artists = []
            self.timeseries_menu_button.hide()
        
        elif plot_type == 'periodic':
            ax.lines = []
            ax.artists = []
            self.periodic_menu_button.hide()

        elif plot_type == 'periodic-violin':
            ax.lines = []
            ax.artists = []
            ax.collections = []
            self.periodic_violin_menu_button.hide()

        elif plot_type == 'metadata':
            ax.texts = []

        elif plot_type == 'distribution':
            ax.lines = []
            ax.artists = []
            self.distribution_menu_button.hide()

        elif plot_type == 'scatter':
            ax.lines = []
            ax.artists = []
            self.scatter_menu_button.hide()

        if plot_type in self.interactive_elements:
            for element in self.interactive_elements[plot_type]['elements']:
                if isinstance(element, dict):
                    for sub_element in element.values():
                        sub_element.hide()
                else:
                    element.hide()

        # hide axis
        ax.axis('off')
        ax.set_visible(False)
        ax.grid(False)
        if plot_type != 'map':
            ax.clear()

        return None

    def activate_axis(self, ax, plot_type):
        """Un-hide axis"""

        # un-hide axis
        ax.axis('on')
        ax.set_visible(True)

        if plot_type != 'legend' and plot_type != 'cb' and plot_type != 'metadata':
            ax.grid(True)

        if plot_type == 'map':
            self.map_menu_button.show()

        elif plot_type == 'timeseries':
            self.timeseries_menu_button.show()

        elif plot_type == 'periodic':
            self.periodic_menu_button.show()

        elif plot_type == 'periodic-violin':
            self.periodic_violin_menu_button.show()

        elif plot_type == 'distribution':
            self.distribution_menu_button.show()
       
        elif plot_type == 'scatter':
            self.scatter_menu_button.show()

        return None

    def update_plot_options(self):
        """ Uncheck checked boxes in plot configuration options under menus and
            reapply check with new data
        """
        
        for option_box in self.options:
            if option_box.isChecked():
                option_box.setCheckState(QtCore.Qt.Unchecked)
                option_box.setCheckState(QtCore.Qt.Checked)

        return None

    def select_all_stations(self):
        """Define function that selects/unselects all plotted stations
        (and associated plots) upon ticking of checkbox"""

        print('SELECT ALL')

        if not self.read_instance.block_MPL_canvas_updates:

            # make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

            # check if checkbox to select all stations is checked or unchecked
            check_state = self.read_instance.ch_select_all.checkState()

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

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
                self.update_associated_selected_station_plots()

            # update plot options
            self.update_plot_options()

            # draw changes
            self.figure.canvas.draw()

        return None

    def select_intersect_stations(self):
        """Define function that selects/unselects intersection of
        stations and all experiment domains (and associated plots)
        upon ticking of checkbox
        """

        print('SELECT INTERSECT')

        if not self.read_instance.block_MPL_canvas_updates:

            # make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

            # check if checkbox to select intersection of stations is checked or unchecked
            check_state = self.read_instance.ch_intersect.checkState()

            # if checkbox is unchecked then unselect all plotted stations
            if check_state == QtCore.Qt.Unchecked:
                self.relative_selected_station_inds = np.array([], dtype=np.int)
                self.absolute_selected_station_inds = np.array([], dtype=np.int)

            # else, if checkbox is checked then select all stations which intersect with all loaded experiment domains
            elif check_state == QtCore.Qt.Checked:

                # if have only observations loaded into memory, select all plotted stations
                if len(self.read_instance.data_labels) == 1:
                    self.relative_selected_station_inds = copy.deepcopy(self.active_map_valid_station_inds)
                    self.absolute_selected_station_inds = np.arange(len(self.relative_selected_station_inds),
                                                                    dtype=np.int)
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
                self.update_associated_selected_station_plots()

            # update plot options
            self.update_plot_options()

            # draw changes
            self.figure.canvas.draw()

        return None

    def select_extent_stations(self):
        """Define function that selects/unselects the
        stations for the current map extent (and associated plots)
        upon ticking of checkbox
        """

        if not self.read_instance.block_MPL_canvas_updates:

            # make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

            # check if checkbox to select extent of stations is checked or unchecked
            check_state = self.read_instance.ch_extent.checkState()

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

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
                self.update_associated_selected_station_plots()

            # update plot options
            self.update_plot_options()

            # draw changes
            self.figure.canvas.draw()

        return None

    def onlassoselect_left(self, verts):
        """Function that handles station selection upon lasso selection with left click.

           Operation:
           Select all stations within lasso boundaries.
           If a click is made rather than using lasso, then select nearest station within tolerance.

           If no station is found with left click, all stations are unselected.
        """

        print('ON LASSO LEFT')

        # check if have any plotted stations on map, if not, return
        if len(self.active_map_valid_station_inds) == 0:
            return

        if self.lock_station_pick == False:

            # lock stations pick
            self.lock_station_pick = True

        # unselect all/intersect checkboxes
        self.read_instance.block_MPL_canvas_updates = True
        self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.ch_extent.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.block_MPL_canvas_updates = False

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
            #take first selected point coordinates and get matches of stations within tolerance 
            current_map_extent = self.plot_axes['map'].get_extent(crs=self.datacrs)
            tolerance = np.average([current_map_extent[1]-current_map_extent[0],current_map_extent[3]-current_map_extent[2]]) / 100.0
            point_coordinates = lasso_path.vertices[0:1,:]
            sub_abs_vals = np.abs(self.map_points_coordinates[None,:,:] - point_coordinates[:,None,:])
            self.absolute_selected_station_inds = np.arange(len(self.active_map_valid_station_inds))[np.all(np.any(sub_abs_vals<=tolerance,axis=0),axis=1)]
            # if more than 1 point selected, limit this to be just nearest point
            if len(self.absolute_selected_station_inds) > 1:
                self.absolute_selected_station_inds = np.array([self.absolute_selected_station_inds[np.argmin(np.sum(sub_abs_vals[0,self.absolute_selected_station_inds,:],axis=1))]], dtype=np.int)

        # get selected station indices with respect to all available stations
        self.relative_selected_station_inds = self.map_selected_station_inds_to_all_available_inds(self.absolute_selected_station_inds)

        # hide lasso after selection
        self.lasso_left.set_visible(False)

        # if selected stations have changed from previous selected, update station selection and associated plots
        if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
            self.update_map_station_selection()
            self.update_associated_selected_station_plots()

        # update plot options
        self.update_plot_options()

        # draw changes
        self.figure.canvas.draw()

        # unlock stations pick 
        self.lock_station_pick = False

        return None

    def onlassoselect_right(self, verts):
        """Function that handles station selection upon lasso selection with right click.
        
           Operation:
           Unselect station/s (if station/s currently selected), 
           or Select station/s (if station/s currently unselected).

           If no station is found with right click, nothing happens.
        """

        print('ON LASSO RIGHT')

        # check if have any plotted stations on map, if not, return
        if len(self.active_map_valid_station_inds) == 0:
            return

        # unselect all/intersect checkboxes
        self.read_instance.block_MPL_canvas_updates = True
        self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.ch_extent.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.block_MPL_canvas_updates = False

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
            #take first selected point coordinates and get matches of stations within tolerance 
            current_map_extent = self.plot_axes['map'].get_extent(crs=self.datacrs)
            tolerance = np.average([current_map_extent[1]-current_map_extent[0],current_map_extent[3]-current_map_extent[2]]) / 100.0
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

        # get new relative selected indices with respect to all available stations
        self.previous_relative_selected_station_inds = previous_relative_selected_station_inds
        self.relative_selected_station_inds = self.map_selected_station_inds_to_all_available_inds(self.absolute_selected_station_inds)

        # hide lasso after selection
        self.lasso_right.set_visible(False)

        # if selected stations have changed from previous selected, update station selection and associated plots
        if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
            self.update_map_station_selection()
            self.update_associated_selected_station_plots()

        # update plot options
        self.update_plot_options()

        # draw changes
        self.figure.canvas.draw()

        return None

    def map_selected_station_inds_to_all_available_inds(self, selected_map_inds):
        """ Takes the indices of selected stations on the map
            (potentially a subset of all available stations), and returns the indices
            of the stations inside the full loaded data arrays.
        """

        # index the array of indices of stations plotted on the map (indexed with respect to
        # all available stations), with the absolute indices of the subset of plotted selected stations
        return self.active_map_valid_station_inds[selected_map_inds]

    def generate_interactive_elements(self):

        """ Function to create settings menus for each plot and their elements."""

        self.interactive_elements = {}
        self.options = []

        # MAP SETTINGS MENU #
        # add button to map to show and hide settings menu
        self.map_menu_button = QtWidgets.QPushButton(self)
        self.map_menu_button.setObjectName('map_menu_button')
        self.map_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.map_menu_button.setIconSize(QtCore.QSize(31, 35))
        self.map_menu_button.setStyleSheet("QPushButton { border: None;} "
                                           "QPushButton:hover { border: None; }")
        self.map_menu_button.hide()

        # add white container
        self.map_container = QtWidgets.QWidget(self)
        self.map_container.setStyleSheet("QWidget { padding: 5px; background-color: white; "
                                         "border: 1px solid; border-color: lightgrey; border-radius: 5px; }")
        self.map_container.setGeometry(self.map_menu_button.geometry().x()-210,
                                       self.map_menu_button.geometry().y()+35, 
                                       235, 350)
        self.map_container.hide()

        # add settings label
        self.map_settings_label = QtWidgets.QLabel("Settings", self)
        self.map_settings_label.setStyleSheet("QLabel { font-weight: bold; }")
        self.map_settings_label.setGeometry(self.map_menu_button.geometry().x()-200, 
                                            self.map_menu_button.geometry().y()+35, 
                                            200, 20)
        self.map_settings_label.hide()

        # add map general text for unselected stations ('Unselected stations')
        self.map_unsel_label = QtWidgets.QLabel("Unselected stations", self)
        self.map_unsel_label.setStyleSheet("QLabel { font-style: italic; }")
        self.map_unsel_label.setGeometry(self.map_menu_button.geometry().x()-200,
                                         self.map_menu_button.geometry().y()+60, 
                                         200, 20)
        self.map_unsel_label.hide()

        # add map markersize slider name ('Size') to layout
        self.map_markersize_unsel_sl_label = QtWidgets.QLabel("Size", self)
        self.map_markersize_unsel_sl_label.setGeometry(self.map_menu_button.geometry().x()-200,
                                                       self.map_menu_button.geometry().y()+85, 
                                                       200, 20)
        self.map_markersize_unsel_sl_label.hide()

        # add map markersize unselected stations slider
        self.map_markersize_unsel_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.map_markersize_unsel_sl.setObjectName('map_markersize_unsel_sl')
        self.map_markersize_unsel_sl.setMinimum(0)
        self.map_markersize_unsel_sl.setMaximum(self.plot_characteristics['map']['marker_unselected']['s']*10)
        self.map_markersize_unsel_sl.setValue(self.plot_characteristics['map']['marker_unselected']['s'])
        self.map_markersize_unsel_sl.setTickInterval(5)
        self.map_markersize_unsel_sl.setGeometry(self.map_menu_button.geometry().x()-200, 
                                                 self.map_menu_button.geometry().y()+110, 
                                                 200, 20)
        self.map_markersize_unsel_sl.hide()

        # add map opacity slider name ('Opacity') to layout
        self.map_opacity_unsel_sl_label = QtWidgets.QLabel("Opacity", self)
        self.map_opacity_unsel_sl_label.setGeometry(self.map_menu_button.geometry().x()-200, 
                                                    self.map_menu_button.geometry().y()+135, 
                                                    200, 20)
        self.map_opacity_unsel_sl_label.hide()

        # add map opacity unselected stations slider
        self.map_opacity_unsel_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.map_opacity_unsel_sl.setObjectName('map_opacity_unsel_sl')
        self.map_opacity_unsel_sl.setMinimum(0)
        self.map_opacity_unsel_sl.setMaximum(10)
        self.map_opacity_unsel_sl.setValue(self.plot_characteristics['map']['marker_unselected']['alpha']*10)
        self.map_opacity_unsel_sl.setTickInterval(1)
        self.map_opacity_unsel_sl.setGeometry(self.map_menu_button.geometry().x()-200, 
                                              self.map_menu_button.geometry().y()+160, 
                                              200, 20)
        self.map_opacity_unsel_sl.hide()

        # add map general text for selected stations ('Selected stations')
        self.map_sel_label = QtWidgets.QLabel("Selected stations", self)
        self.map_sel_label.setStyleSheet("QLabel { font-style: italic; }")
        self.map_sel_label.setGeometry(self.map_menu_button.geometry().x()-200, 
                                       self.map_menu_button.geometry().y()+185, 
                                       200, 20)
        self.map_sel_label.hide()

        # add map markersize slider name ('Size') to layout
        self.map_markersize_sel_sl_label = QtWidgets.QLabel("Size", self)
        self.map_markersize_sel_sl_label.setGeometry(self.map_menu_button.geometry().x()-200, 
                                                     self.map_menu_button.geometry().y()+210, 
                                                     200, 20)
        self.map_markersize_sel_sl_label.hide()

        # add map markersize selected stations slider
        self.map_markersize_sel_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.map_markersize_sel_sl.setObjectName('map_markersize_sel_sl')
        self.map_markersize_sel_sl.setMinimum(0)
        self.map_markersize_sel_sl.setMaximum(self.plot_characteristics['map']['marker_selected']['s']*10)
        self.map_markersize_sel_sl.setValue(self.plot_characteristics['map']['marker_unselected']['s'])
        self.map_markersize_sel_sl.setTickInterval(5)
        self.map_markersize_sel_sl.setGeometry(self.map_menu_button.geometry().x()-200, 
                                               self.map_menu_button.geometry().y()+235, 
                                               200, 20)
        self.map_markersize_sel_sl.hide()

        # add map opacity slider name ('Opacity') to layout
        self.map_opacity_sel_sl_label = QtWidgets.QLabel("Opacity", self)
        self.map_opacity_sel_sl_label.setGeometry(self.map_menu_button.geometry().x()-200, 
                                                  self.map_menu_button.geometry().y()+260, 
                                                  200, 20)
        self.map_opacity_sel_sl_label.hide()

        # add map opacity selected stations slider
        self.map_opacity_sel_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.map_opacity_sel_sl.setObjectName('map_opacity_sel_sl')
        self.map_opacity_sel_sl.setMinimum(0)
        self.map_opacity_sel_sl.setMaximum(10)
        self.map_opacity_sel_sl.setValue(self.plot_characteristics['map']['marker_selected']['alpha']*10)
        self.map_opacity_sel_sl.setTickInterval(1)
        self.map_opacity_sel_sl.setGeometry(self.map_menu_button.geometry().x()-200, 
                                            self.map_menu_button.geometry().y()+285, 
                                            200, 20)
        self.map_opacity_sel_sl.hide()

        # add map options name ('Options') to layout
        self.map_options_label = QtWidgets.QLabel("Options", self)
        self.map_options_label.setGeometry(self.map_menu_button.geometry().x()-200,
                                           self.map_menu_button.geometry().y()+310, 
                                           200, 20)
        self.map_options_label.hide()

        # add map options checkboxes
        self.map_options = dict()
        y_move = 0
        for option in self.plot_characteristics['map']['plot_options']:

            self.map_options[option] = QtWidgets.QCheckBox(option, self)
            self.map_options[option].setObjectName('map_option_' + option)
            self.map_options[option].setGeometry(self.map_menu_button.geometry().x()-200, 
                                                 self.map_menu_button.geometry().y()+345+y_move, 
                                                 200, 20)
            self.map_options[option].stateChanged.connect(self.update_plot_option)
            self.map_options[option].hide()
            y_move += 25 

            self.options.append(self.map_options[option])

        # set show/hide actions
        self.map_elements = [self.map_container, self.map_settings_label, 
                             self.map_unsel_label, self.map_markersize_unsel_sl_label, 
                             self.map_markersize_unsel_sl, self.map_opacity_unsel_sl_label, 
                             self.map_opacity_unsel_sl, self.map_sel_label, 
                             self.map_markersize_sel_sl, self.map_markersize_sel_sl_label,
                             self.map_opacity_sel_sl_label, self.map_opacity_sel_sl,
                             self.map_options_label, self.map_options
                             ]
        self.interactive_elements['map'] = {'button': self.map_menu_button, 
                                            'hidden': True,
                                            'elements': self.map_elements,
                                            'markersize_sl': [self.map_markersize_unsel_sl, 
                                                              self.map_markersize_sel_sl],
                                            'opacity_sl': [self.map_opacity_unsel_sl, 
                                                           self.map_opacity_sel_sl],
                                            'linewidth_sl': [],
                                            'widths_sl': []
                                            }
        self.map_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.map_markersize_sel_sl.valueChanged.connect(self.update_markersize_func)
        self.map_markersize_unsel_sl.valueChanged.connect(self.update_markersize_func)
        self.map_opacity_sel_sl.valueChanged.connect(self.update_opacity_func)
        self.map_opacity_unsel_sl.valueChanged.connect(self.update_opacity_func)
        
        # TIMESERIES PLOT SETTINGS MENU #
        # add button to timeseries to show and hide settings menu
        self.timeseries_menu_button = QtWidgets.QPushButton(self)
        self.timeseries_menu_button.setObjectName('timeseries_menu_button')
        self.timeseries_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.timeseries_menu_button.setIconSize(QtCore.QSize(31, 35))
        self.timeseries_menu_button.setStyleSheet("QPushButton { border: None; }"
                                                  "QPushButton:hover { border: None; }")
        self.timeseries_menu_button.hide()

        # add white container
        self.timeseries_container = QtWidgets.QWidget(self)
        self.timeseries_container.setStyleSheet("QWidget { padding: 5px; background-color: white; "
                                                "border: 1px solid; border-color: lightgrey; border-radius: 5px; }")
        self.timeseries_container.setGeometry(self.timeseries_menu_button.geometry().x()-210,
                                              self.timeseries_menu_button.geometry().y()+25, 
                                              235, 150)
        self.timeseries_container.hide()

        # add settings label
        self.timeseries_settings_label = QtWidgets.QLabel("Settings", self)
        self.timeseries_settings_label.setStyleSheet("QLabel { font-weight: bold; }")
        self.timeseries_settings_label.setGeometry(self.timeseries_menu_button.geometry().x()-200, 
                                                   self.timeseries_menu_button.geometry().y()+25, 
                                                   200, 20)
        self.timeseries_settings_label.hide()

        # add timeseries markersize slider name ('Markersize') to layout
        self.timeseries_markersize_sl_label = QtWidgets.QLabel("Markersize", self)
        self.timeseries_markersize_sl_label.setGeometry(self.timeseries_menu_button.geometry().x()-200,
                                                        self.timeseries_menu_button.geometry().y()+50, 
                                                        200, 20)
        self.timeseries_markersize_sl_label.hide()

        # add timeseries markersize slider
        self.timeseries_markersize_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.timeseries_markersize_sl.setObjectName('timeseries_markersize_sl')
        self.timeseries_markersize_sl.setMinimum(0)
        self.timeseries_markersize_sl.setMaximum(self.plot_characteristics['timeseries']['plot']['markersize']*10)
        self.timeseries_markersize_sl.setValue(self.plot_characteristics['timeseries']['plot']['markersize'])
        self.timeseries_markersize_sl.setTickInterval(2)
        self.timeseries_markersize_sl.setGeometry(self.timeseries_menu_button.geometry().x()-200, 
                                                  self.timeseries_menu_button.geometry().y()+75, 
                                                  200, 20)
        self.timeseries_markersize_sl.hide()

        # add timeseries plot options name ('Options') to layout
        self.timeseries_options_label = QtWidgets.QLabel("Options", self)
        self.timeseries_options_label.setGeometry(self.timeseries_menu_button.geometry().x()-200,
                                                  self.timeseries_menu_button.geometry().y()+100, 
                                                  200, 20)
        self.timeseries_options_label.hide()

        # add timeseries plot options checkboxes
        self.timeseries_options = dict()
        x_move, y_move = 0, 0
        for option in self.plot_characteristics['timeseries']['plot_options']:
            
            self.timeseries_options[option] = QtWidgets.QCheckBox(option, self)
            self.timeseries_options[option].setObjectName('timeseries_option_' + option)
            self.timeseries_options[option].setGeometry(self.timeseries_menu_button.geometry().x()-200+x_move, 
                                                        self.timeseries_menu_button.geometry().y()+125+y_move, 
                                                        200, 20)
            self.timeseries_options[option].stateChanged.connect(self.update_plot_option)
            self.timeseries_options[option].hide()

            # create new column after 2 elements
            if y_move == 25:
                x_move = 100
                y_move = 0
            else:
                y_move += 25 
            
            self.options.append(self.timeseries_options[option])

        # set show/hide actions
        self.timeseries_elements = [self.timeseries_container, self.timeseries_settings_label, 
                                    self.timeseries_markersize_sl_label, self.timeseries_markersize_sl,
                                    self.timeseries_options_label, self.timeseries_options]
        self.interactive_elements['timeseries'] = {'button': self.timeseries_menu_button, 
                                                   'hidden': True,
                                                   'elements': self.timeseries_elements,
                                                   'markersize_sl': [self.timeseries_markersize_sl],
                                                   'opacity_sl': [],
                                                   'linewidth_sl': [],
                                                   'widths_sl': []
                                                   }
        self.timeseries_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.timeseries_markersize_sl.valueChanged.connect(self.update_markersize_func)

        # PERIODIC PLOT SETTINGS MENU #
        # add button to periodic plot to show and hide settings menu
        self.periodic_menu_button = QtWidgets.QPushButton(self)
        self.periodic_menu_button.setObjectName('periodic_menu_button')
        self.periodic_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.periodic_menu_button.setIconSize(QtCore.QSize(31, 35))
        self.periodic_menu_button.setStyleSheet("QPushButton { border: None; } "
                                                "QPushButton:hover { border: None; }")
        self.periodic_menu_button.hide()

        # add white container
        self.periodic_container = QtWidgets.QWidget(self)
        self.periodic_container.setStyleSheet("QWidget { padding: 5px; background-color: white; "
                                              "border: 1px solid; border-color: lightgrey; border-radius: 5px; }")
        self.periodic_container.setGeometry(self.periodic_menu_button.geometry().x()-220, 
                                            self.periodic_menu_button.geometry().y()+25, 
                                            235, 200)
        self.periodic_container.hide()

        # add settings label
        self.periodic_settings_label = QtWidgets.QLabel("Settings", self)
        self.periodic_settings_label.setStyleSheet("QLabel { font-weight: bold; }")
        self.periodic_settings_label.setGeometry(self.periodic_menu_button.geometry().x()-210, 
                                                 self.periodic_menu_button.geometry().y()+25, 
                                                 200, 20)
        self.periodic_settings_label.hide()

        # add periodic markersize slider name ('Size') to layout
        self.periodic_markersize_sl_label = QtWidgets.QLabel("Size", self)
        self.periodic_markersize_sl_label.setGeometry(self.periodic_menu_button.geometry().x()-210, 
                                                      self.periodic_menu_button.geometry().y()+50, 
                                                      200, 20)
        self.periodic_markersize_sl_label.hide()

        # add periodic markersize slider
        self.periodic_markersize_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.periodic_markersize_sl.setObjectName('periodic_markersize_sl')
        self.periodic_markersize_sl.setMinimum(0)
        self.periodic_markersize_sl.setMaximum(self.plot_characteristics['periodic']['plot']['markersize']*10)
        self.periodic_markersize_sl.setValue(self.plot_characteristics['periodic']['plot']['markersize'])
        self.periodic_markersize_sl.setTickInterval(2)
        self.periodic_markersize_sl.setGeometry(self.periodic_menu_button.geometry().x()-210, 
                                                self.periodic_menu_button.geometry().y()+75, 
                                                200, 20)
        self.periodic_markersize_sl.hide()

        # add periodic line width slider name ('Line width') to layout
        self.periodic_linewidth_sl_label = QtWidgets.QLabel("Line width", self)
        self.periodic_linewidth_sl_label.setGeometry(self.periodic_menu_button.geometry().x()-210, 
                                                     self.periodic_menu_button.geometry().y()+100, 
                                                     200, 20)
        self.periodic_linewidth_sl_label.hide()

        # add periodic line width slider
        self.periodic_linewidth_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.periodic_linewidth_sl.setObjectName('periodic_linewidth_sl')
        self.periodic_linewidth_sl.setMinimum(0)
        self.periodic_linewidth_sl.setMaximum(self.plot_characteristics['periodic']['plot']['linewidth']*100)
        self.periodic_linewidth_sl.setValue(self.plot_characteristics['periodic']['plot']['linewidth']*10)
        self.periodic_linewidth_sl.setTickInterval(2)
        self.periodic_linewidth_sl.setGeometry(self.periodic_menu_button.geometry().x()-210, 
                                                self.periodic_menu_button.geometry().y()+125, 
                                                200, 20)
        self.periodic_linewidth_sl.hide()

        # add periodic plot options name ('Options') to layout
        self.periodic_options_label = QtWidgets.QLabel("Options", self)
        self.periodic_options_label.setGeometry(self.periodic_menu_button.geometry().x()-210,
                                                self.periodic_menu_button.geometry().y()+150, 
                                                200, 20)
        self.periodic_options_label.hide()

        # add periodic plot options checkboxes
        self.periodic_options = dict()
        x_move, y_move = 0, 0
        for option in self.plot_characteristics['periodic']['plot_options']:
            
            self.periodic_options[option] = QtWidgets.QCheckBox(option, self)
            self.periodic_options[option].setObjectName('periodic_option_' + option)
            self.periodic_options[option].setGeometry(self.periodic_menu_button.geometry().x()-210+x_move, 
                                                      self.periodic_menu_button.geometry().y()+175+y_move, 
                                                      200, 20)
            self.periodic_options[option].stateChanged.connect(self.update_plot_option)
            self.periodic_options[option].hide()

            # create new column after 2 elements
            if y_move == 25:
                x_move = 100
                y_move = 0
            else:
                y_move += 25 
            
            self.options.append(self.periodic_options[option])

        # set show/hide actions
        self.periodic_elements = [self.periodic_container, self.periodic_settings_label, 
                                  self.periodic_markersize_sl_label, self.periodic_markersize_sl,
                                  self.periodic_linewidth_sl_label, self.periodic_linewidth_sl,
                                  self.periodic_options_label, self.periodic_options
                                 ]
        self.interactive_elements['periodic'] = {'button': self.periodic_menu_button, 
                                                 'hidden': True,
                                                 'elements': self.periodic_elements,
                                                 'markersize_sl': [self.periodic_markersize_sl],
                                                 'opacity_sl': [],
                                                 'linewidth_sl': [self.periodic_linewidth_sl],
                                                 'widths_sl': []
                                                 }
        self.periodic_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.periodic_markersize_sl.valueChanged.connect(self.update_markersize_func)
        self.periodic_linewidth_sl.valueChanged.connect(self.update_linewidth_func)

        # PERIODIC VIOLIN PLOT SETTINGS MENU #
        # add button to periodic violin plot to show and hide settings menu
        self.periodic_violin_menu_button = QtWidgets.QPushButton(self)
        self.periodic_violin_menu_button.setObjectName('periodic_violin_menu_button')
        self.periodic_violin_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.periodic_violin_menu_button.setIconSize(QtCore.QSize(31, 35))
        self.periodic_violin_menu_button.setStyleSheet("QPushButton { border: None; } "
                                                       "QPushButton:hover { border: None; }")
        self.periodic_violin_menu_button.hide()

        # add white container
        self.periodic_violin_container = QtWidgets.QWidget(self)
        self.periodic_violin_container.setStyleSheet("QWidget { padding: 5px; background-color: white; "
                                                     "border: 1px solid; border-color: lightgrey; border-radius: 5px; }")
        self.periodic_violin_container.setGeometry(self.periodic_violin_menu_button.geometry().x()-210, 
                                                   self.periodic_violin_menu_button.geometry().y()+25, 
                                                   235, 200)
        self.periodic_violin_container.hide()

        # add settings label
        self.periodic_violin_settings_label = QtWidgets.QLabel("Settings", self)
        self.periodic_violin_settings_label.setStyleSheet("QLabel { font-weight: bold; }")
        self.periodic_violin_settings_label.setGeometry(self.periodic_violin_menu_button.geometry().x()-200,
                                                        self.periodic_violin_menu_button.geometry().y()+25, 
                                                        200, 20)
        self.periodic_violin_settings_label.hide()

        # add periodic violin markersize slider name ('Size') to layout
        self.periodic_violin_markersize_sl_label = QtWidgets.QLabel("Size", self)
        self.periodic_violin_markersize_sl_label.setGeometry(self.periodic_violin_menu_button.geometry().x()-200,
                                                             self.periodic_violin_menu_button.geometry().y()+50, 
                                                             200, 20)
        self.periodic_violin_markersize_sl_label.hide()

        # add periodic violin markersize slider
        self.periodic_violin_markersize_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.periodic_violin_markersize_sl.setObjectName('periodic_violin_markersize_sl')
        self.periodic_violin_markersize_sl.setMinimum(0)
        self.periodic_violin_markersize_sl.setMaximum(self.plot_characteristics['periodic-violin']['plot']['p50']['markersize']*10)
        self.periodic_violin_markersize_sl.setValue(self.plot_characteristics['periodic-violin']['plot']['p50']['markersize'])
        self.periodic_violin_markersize_sl.setTickInterval(2)
        self.periodic_violin_markersize_sl.setGeometry(self.periodic_violin_menu_button.geometry().x()-200,
                                                       self.periodic_menu_button.geometry().y()+75, 
                                                       200, 20)
        self.periodic_violin_markersize_sl.hide()

        # add periodic violin line width slider name ('Line width') to layout
        self.periodic_violin_linewidth_sl_label = QtWidgets.QLabel("Line width", self)
        self.periodic_violin_linewidth_sl_label.setGeometry(self.periodic_violin_menu_button.geometry().x()-200, 
                                                            self.periodic_violin_menu_button.geometry().y()+100, 
                                                            200, 20)
        self.periodic_violin_linewidth_sl_label.hide()

        # add periodic violin line width slider
        self.periodic_violin_linewidth_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.periodic_violin_linewidth_sl.setObjectName('periodic_violin_linewidth_sl')
        self.periodic_violin_linewidth_sl.setMinimum(0)
        self.periodic_violin_linewidth_sl.setMaximum(self.plot_characteristics['periodic-violin']['plot']['p50']['linewidth']*100)
        self.periodic_violin_linewidth_sl.setValue(self.plot_characteristics['periodic-violin']['plot']['p50']['linewidth']*10)
        self.periodic_violin_linewidth_sl.setTickInterval(2)
        self.periodic_violin_linewidth_sl.setGeometry(self.periodic_violin_menu_button.geometry().x()-200, 
                                                      self.periodic_violin_menu_button.geometry().y()+125, 
                                                      200, 20)
        self.periodic_violin_linewidth_sl.hide()

        # add periodic violin plot options name ('Options') to layout
        self.periodic_violin_options_label = QtWidgets.QLabel("Options", self)
        self.periodic_violin_options_label.setGeometry(self.periodic_violin_menu_button.geometry().x()-200,
                                                       self.periodic_violin_menu_button.geometry().y()+150, 
                                                       200, 20)
        self.periodic_violin_options_label.hide()

        # add periodic violin plot options checkboxes
        self.periodic_violin_options = dict()
        x_move, y_move = 0, 0
        for option in self.plot_characteristics['periodic-violin']['plot_options']:
            
            self.periodic_violin_options[option] = QtWidgets.QCheckBox(option, self)
            self.periodic_violin_options[option].setObjectName('periodic_violin_option_' + option)
            self.periodic_violin_options[option].setGeometry(self.periodic_violin_menu_button.geometry().x()-200+x_move, 
                                                             self.periodic_violin_menu_button.geometry().y()+175+y_move, 
                                                             200, 20)
            self.periodic_violin_options[option].stateChanged.connect(self.update_plot_option)
            self.periodic_violin_options[option].hide()

            # create new column after 2 elements
            if y_move == 25:
                x_move = 100
                y_move = 0
            else:
                y_move += 25 

            self.options.append(self.periodic_violin_options[option])

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
                                                        'linewidth_sl': [self.periodic_violin_linewidth_sl],
                                                        'widths_sl': []
                                                        }
        self.periodic_violin_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.periodic_violin_markersize_sl.valueChanged.connect(self.update_markersize_func)
        self.periodic_violin_linewidth_sl.valueChanged.connect(self.update_linewidth_func)

        # DISTRIBUTION PLOT SETTINGS MENU #
        # add button to distribution plot to show and hide settings menu
        self.distribution_menu_button = QtWidgets.QPushButton(self)
        self.distribution_menu_button.setObjectName('distribution_menu_button')
        self.distribution_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.distribution_menu_button.setIconSize(QtCore.QSize(31, 35))
        self.distribution_menu_button.setStyleSheet("QPushButton { border: None; }"
                                                    "QPushButton:hover { border: None; }")
        self.distribution_menu_button.hide()

        # add white container
        self.distribution_container = QtWidgets.QWidget(self)
        self.distribution_container.setStyleSheet("QWidget { padding: 5px; background-color: white; "
                                                  "border: 1px solid; border-color: lightgrey; border-radius: 5px; }")
        self.distribution_container.setGeometry(self.distribution_menu_button.geometry().x()-210,
                                                self.distribution_menu_button.geometry().y()+25, 
                                                235, 150)
        self.distribution_container.hide()

        # add settings label
        self.distribution_settings_label = QtWidgets.QLabel("Settings", self)
        self.distribution_settings_label.setStyleSheet("QLabel { font-weight: bold; }")
        self.distribution_settings_label.setGeometry(self.distribution_menu_button.geometry().x()-200, 
                                                     self.distribution_menu_button.geometry().y()+25, 
                                                     200, 20)
        self.distribution_settings_label.hide()

        # add distribution plot line width slider name ('Line width') to layout
        self.distribution_linewidth_sl_label = QtWidgets.QLabel("Line width", self)
        self.distribution_linewidth_sl_label.setGeometry(self.distribution_menu_button.geometry().x()-200,
                                                         self.distribution_menu_button.geometry().y()+50, 
                                                         200, 20)
        self.distribution_linewidth_sl_label.hide()

        # add distribution plot line width slider
        self.distribution_linewidth_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.distribution_linewidth_sl.setObjectName('distribution_linewidth_sl')
        self.distribution_linewidth_sl.setMinimum(0)
        self.distribution_linewidth_sl.setMaximum(self.plot_characteristics['distribution']['plot']['linewidth']*100)
        self.distribution_linewidth_sl.setValue(self.plot_characteristics['distribution']['plot']['linewidth']*10)
        self.distribution_linewidth_sl.setTickInterval(2)
        self.distribution_linewidth_sl.setGeometry(self.distribution_menu_button.geometry().x()-200, 
                                                   self.distribution_menu_button.geometry().y()+75, 
                                                   200, 20)
        self.distribution_linewidth_sl.hide()

        # add distribution plot options name ('Options') to layout
        self.distribution_options_label = QtWidgets.QLabel("Options", self)
        self.distribution_options_label.setGeometry(self.distribution_menu_button.geometry().x()-200,
                                                    self.distribution_menu_button.geometry().y()+100, 
                                                    200, 20)
        self.distribution_options_label.hide()

        # add distribution plot options checkboxes
        self.distribution_options = dict()
        x_move, y_move = 0, 0
        for option in self.plot_characteristics['distribution']['plot_options']:
            
            self.distribution_options[option] = QtWidgets.QCheckBox(option, self)
            self.distribution_options[option].setObjectName('distribution_option_' + option)
            self.distribution_options[option].setGeometry(self.distribution_menu_button.geometry().x()-200+x_move, 
                                                          self.distribution_menu_button.geometry().y()+125+y_move, 
                                                          200, 20)
            self.distribution_options[option].stateChanged.connect(self.update_plot_option)
            self.distribution_options[option].hide()

            # create new column after 2 elements
            if y_move == 25:
                x_move = 100
                y_move = 0
            else:
                y_move += 25 
            
            self.options.append(self.distribution_options[option])

        # set show/hide actions
        self.distribution_elements = [self.distribution_container, self.distribution_settings_label, 
                                      self.distribution_linewidth_sl_label, self.distribution_linewidth_sl,
                                      self.distribution_options_label, self.distribution_options]
        self.interactive_elements['distribution'] = {'button': self.distribution_menu_button, 
                                                     'hidden': True,
                                                     'elements': self.distribution_elements,
                                                     'markersize_sl': [],
                                                     'opacity_sl': [],
                                                     'linewidth_sl': [self.distribution_linewidth_sl],
                                                     'widths_sl': []
                                                    }
        self.distribution_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.distribution_linewidth_sl.valueChanged.connect(self.update_linewidth_func)

        # SCATTER PLOT SETTINGS MENU #
        # add button scatter plot to show and hide settings menu
        self.scatter_menu_button = QtWidgets.QPushButton(self)
        self.scatter_menu_button.setObjectName('scatter_menu_button')
        self.scatter_menu_button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, "resources/menu_icon.png")))
        self.scatter_menu_button.setIconSize(QtCore.QSize(31, 35))
        self.scatter_menu_button.setStyleSheet("QPushButton { border: None; }"
                                               "QPushButton:hover { border: None; }")
        self.scatter_menu_button.hide()

        # add white container
        self.scatter_container = QtWidgets.QWidget(self)
        self.scatter_container.setStyleSheet("QWidget { padding: 5px; background-color: white; "
                                             "border: 1px solid; border-color: lightgrey; border-radius: 5px; }")
        self.scatter_container.setGeometry(self.scatter_menu_button.geometry().x()-210,
                                           self.scatter_menu_button.geometry().y()+25, 
                                           235, 150)
        self.scatter_container.hide()

        # add settings label
        self.scatter_settings_label = QtWidgets.QLabel("Settings", self)
        self.scatter_settings_label.setStyleSheet("QLabel { font-weight: bold; }")
        self.scatter_settings_label.setGeometry(self.scatter_menu_button.geometry().x()-200, 
                                                self.scatter_menu_button.geometry().y()+25, 
                                                200, 20)
        self.scatter_settings_label.hide()

        # add scatter plot markersize slider name ('Size') to layout
        self.scatter_markersize_sl_label = QtWidgets.QLabel("Size", self)
        self.scatter_markersize_sl_label.setGeometry(self.scatter_menu_button.geometry().x()-200,
                                                    self.scatter_menu_button.geometry().y()+50, 
                                                    200, 20)
        self.scatter_markersize_sl_label.hide()

        # add scatter plot markersize slider
        self.scatter_markersize_sl = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.scatter_markersize_sl.setObjectName('scatter_markersize_sl')
        self.scatter_markersize_sl.setMinimum(0)
        self.scatter_markersize_sl.setMaximum(self.plot_characteristics['scatter']['plot']['markersize']*100)
        self.scatter_markersize_sl.setValue(self.plot_characteristics['scatter']['plot']['markersize']*10)
        self.scatter_markersize_sl.setTickInterval(2)
        self.scatter_markersize_sl.setGeometry(self.scatter_menu_button.geometry().x()-200, 
                                              self.scatter_menu_button.geometry().y()+75, 
                                              200, 20)
        self.scatter_markersize_sl.hide()

        # add scatter plot options name ('Options') to layout
        self.scatter_options_label = QtWidgets.QLabel("Options", self)
        self.scatter_options_label.setGeometry(self.scatter_menu_button.geometry().x()-200,
                                               self.scatter_menu_button.geometry().y()+100, 
                                               200, 20)
        self.scatter_options_label.hide()

        # add scatter plot options checkboxes
        self.scatter_options = dict()
        x_move, y_move = 0, 0
        for option in self.plot_characteristics['scatter']['plot_options']:
            
            self.scatter_options[option] = QtWidgets.QCheckBox(option, self)
            self.scatter_options[option].setObjectName('scatter_option_' + option)
            self.scatter_options[option].setGeometry(self.scatter_menu_button.geometry().x()-200+x_move, 
                                                     self.scatter_menu_button.geometry().y()+125+y_move, 
                                                     200, 20)
            self.scatter_options[option].stateChanged.connect(self.update_plot_option)
            self.scatter_options[option].hide()

            # create new column after 2 elements
            if y_move == 25:
                x_move = 100
                y_move = 0
            else:
                y_move += 25 
            
            self.options.append(self.scatter_options[option])

        # set show/hide actions
        self.scatter_elements = [self.scatter_container, self.scatter_settings_label, 
                                 self.scatter_markersize_sl_label, self.scatter_markersize_sl,
                                 self.scatter_options_label, self.scatter_options]
        self.interactive_elements['scatter'] = {'button': self.scatter_menu_button, 
                                                'hidden': True,
                                                'elements': self.scatter_elements,
                                                'markersize_sl': [self.scatter_markersize_sl],
                                                'opacity_sl': [],
                                                'linewidth_sl': [],
                                                'widths_sl': []
                                               }
        self.scatter_menu_button.clicked.connect(self.interactive_elements_button_func)
        self.scatter_markersize_sl.valueChanged.connect(self.update_markersize_func)

        # create array with buttons and elements to edit when the canvas is resized or the plots are changed
        self.buttons = [self.map_menu_button, self.timeseries_menu_button,
                        self.periodic_menu_button, self.periodic_violin_menu_button,
                        self.distribution_menu_button, self.scatter_menu_button]

        self.elements = [self.map_elements, self.timeseries_elements, 
                         self.periodic_elements, self.periodic_violin_elements,
                         self.distribution_elements, self.scatter_elements]

        return None

    def interactive_elements_button_func(self):
        """Function to show and hide elements in setting menus"""

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
        """Function to handle the update of the markers size."""

        event_source = self.sender()
        loc = [1 if '_sel' in event_source.objectName() else 0][0]
        for key, val in self.interactive_elements.items():
            if event_source in self.interactive_elements[key]['markersize_sl']:
                markersize = self.interactive_elements[key]['markersize_sl'][loc].value()
                break
        self.update_markersize(self.plot_axes[key], key, markersize, event_source)

        return None

    def update_opacity_func(self):
        """Function to handle the update of the markers opacity."""

        event_source = self.sender()
        loc = [1 if '_sel' in event_source.objectName() else 0][0]
        for key, val in self.interactive_elements.items():
            if event_source in self.interactive_elements[key]['opacity_sl']:
                opacity = self.interactive_elements[key]['opacity_sl'][loc].value()/10
                break
        self.update_opacity(self.plot_axes[key], key, opacity, event_source)

        return None

    def update_linewidth_func(self):
        """Function to handle the update of the lines widths."""

        event_source = self.sender()
        for key, val in self.interactive_elements.items():
            if event_source in self.interactive_elements[key]['linewidth_sl']:
                linewidth = self.interactive_elements[key]['linewidth_sl'][0].value()/10
                break
        self.update_linewidth(self.plot_axes[key], key, linewidth)

        return None

    def update_violin_widths_func(self):
        """Function to handle the update of the violin widths."""

        event_source = self.sender()
        for key, val in self.interactive_elements.items():
            if event_source in self.interactive_elements[key]['widths_sl']:
                widths = self.interactive_elements[key]['widths_sl'][0].value()/10
                break
        self.update_violin_widths(self.plot_axes[key], key, widths)

        return None

    def update_plot_option(self):
        """ Function to handle the update of the plot options. """

        # get option and plot names
        event_source = self.sender()
        option = event_source.objectName().split('option_')[1]
        key = event_source.objectName().split('_option')[0]
        if key == 'periodic_violin':
            key = 'periodic-violin'

        # is box checked?
        check_state = event_source.isChecked()
        
        # options 'logy' and 'logx'
        if option == 'logy' or option == 'logx':
            # log y axis if box is checked
            if check_state:
                if isinstance(self.plot_axes[key], dict):
                    for sub_ax in self.plot_axes[key].values():
                        self.plot.log_axes(sub_ax,
                                           option, 
                                           event_source,
                                           self.plot_characteristics[key])
                else:
                    self.plot.log_axes(self.plot_axes[key], 
                                       option, 
                                       event_source,
                                       self.plot_characteristics[key])
            # undo log y axis if box is unchecked
            elif not check_state:
                if isinstance(self.plot_axes[key], dict):
                    for sub_ax in self.plot_axes[key].values():
                        self.plot.log_axes(sub_ax, 
                                           option,
                                           event_source, 
                                           self.plot_characteristics[key],
                                           undo=True)
                else:
                    self.plot.log_axes(self.plot_axes[key], 
                                       option, 
                                       event_source,
                                       self.plot_characteristics[key],
                                       undo=True)

        # option 'annotate'
        if option == 'annotate':
            # add annotation if box is checked
            if check_state:
                if isinstance(self.plot_axes[key], dict):
                    for sub_ax in self.plot_axes[key].values():
                        self.plot.annotation(sub_ax,
                                             self.read_instance.networkspeci,
                                             list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                             self.plot_characteristics[key], 
                                             key,
                                             plot_options=[])
                        break
                else:
                    self.plot.annotation(self.plot_axes[key], 
                                         self.read_instance.networkspeci,
                                         list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                         self.plot_characteristics[key], 
                                         key,
                                         plot_options=[])
            # remove annotation if box is unchecked
            elif not check_state:
                if isinstance(self.plot_axes[key], dict):
                    for sub_ax in self.plot_axes[key].values():
                        self.plot.annotation(sub_ax, 
                                             self.read_instance.networkspeci,
                                             list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                             self.plot_characteristics[key], 
                                             key,
                                             plot_options=[], 
                                             undo=True)
                else:
                    self.plot.annotation(self.plot_axes[key], 
                                         self.read_instance.networkspeci,
                                         list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                         self.plot_characteristics[key], 
                                         key,
                                         plot_options=[], 
                                         undo=True)
        
        # option 'trend'
        if option == 'trend':
            # add trend line if box is checked
            if check_state:
                self.plot.trend(self.plot_axes[key], 
                                self.read_instance.networkspeci,
                                list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                self.plot_characteristics[key], 
                                plot_options=[])
            # remove trend line if box is unchecked
            elif not check_state:
                self.plot.trend(self.plot_axes[key], 
                                self.read_instance.networkspeci,
                                list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                self.plot_characteristics[key], 
                                plot_options=[], 
                                undo=True)
        
        if option == 'regression':
            # add regression line if box is checked
            if check_state:
                self.plot.linear_regression(self.plot_axes[key], 
                                            self.read_instance.networkspeci,
                                            list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                            self.plot_characteristics[key],  
                                            plot_options=[])
            # remove regression line if box is unchecked
            elif not check_state:
                self.plot.linear_regression(self.plot_axes[key], 
                                            self.read_instance.networkspeci,
                                            list(self.selected_station_data[self.read_instance.networkspeci].keys()), 
                                            self.plot_characteristics[key],  
                                            plot_options=[],
                                            undo=True)

        # draw changes
        self.figure.canvas.draw()

        return None

    def update_markersize(self, ax, plot_type, markersize, event_source):
        """ Update markers size for each plot type. """
        
        # set markersize
        if plot_type in ['timeseries', 'periodic', 'scatter']:

            if type(ax) == dict:
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

            if event_source == self.interactive_elements[plot_type]['markersize_sl'][0]:

                # set markersize for unselected stations
                markersizes[self.absolute_non_selected_station_inds] = markersize
                for collection in self.plot_axes['map'].collections:
                    if isinstance(collection, matplotlib.collections.PathCollection):
                        collection.set_sizes(markersizes)

                # update characteristics per plot type
                self.plot_characteristics['map']['marker_unselected']['s'] = markersize

            elif event_source == self.interactive_elements[plot_type]['markersize_sl'][1]:

                # set markersize for selected stations
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
            if event_source == self.interactive_elements[plot_type]['opacity_sl'][0]:

                # set opacity for unselected stations
                opacities[self.absolute_non_selected_station_inds, -1] = opacity
                for collection in self.plot_axes['map'].collections:
                    if isinstance(collection, matplotlib.collections.PathCollection):
                        collection.set_facecolor(opacities)

                # update characteristics per plot type
                self.plot_characteristics['map']['marker_unselected']['alpha'] = opacity

            elif event_source == self.interactive_elements[plot_type]['opacity_sl'][1]:

                # set opacity for selected stations
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
        if type(ax) == dict:
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

    def update_violin_widths(self, ax, plot_type, width):
        """ Update violin widths for violin plots. """
        
        # set violin widths
        if plot_type == 'periodic-violin':
            if type(ax) == dict:
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

    def create_station_annotation(self):
        """ Create annotation at (0, 0) that will be updated later. """

        # get axes transformation
        transform = self.datacrs._as_mpl_transform(self.plot_axes['map'])

        # in the newest version of matplotlib, s corresponds to text
        self.station_annotation = self.plot_axes['map'].annotate(s='', xy=(0, 0), xycoords=transform,
                                                                  **self.plot_characteristics['map']['stations_annotate'],
                                                                  bbox={**self.plot_characteristics['map']['stations_annotate_bbox']},
                                                                  arrowprops={**self.plot_characteristics['map']['stations_annotate_arrowprops']})

        return None

    def update_station_annotation(self, annotation_index):
        """ Update annotation for each station that is hovered. """

        # retrieve stations references and coordinates
        station_reference = self.read_instance.station_references[self.read_instance.networkspeci][self.active_map_valid_station_inds][annotation_index['ind'][0]]
        station_location = self.plot.stations_scatter.get_offsets()[annotation_index['ind'][0]]

        # update location
        self.station_annotation.xy = station_location

        # update bbox position
        lat_min = self.plot_axes['map'].get_extent(crs=self.datacrs)[2]
        lat_max = self.plot_axes['map'].get_extent(crs=self.datacrs)[3]
        if station_location[1] > ((lat_max + lat_min) / 2):
            self.station_annotation.set_y(-10)
            self.station_annotation.set_va('top')
        else:
            self.station_annotation.set_y(10)
            self.station_annotation.set_va('bottom')

        # create annotation text
        text_label = ('Station: {0}\n').format(station_reference)
        text_label += ('Lon: {0:.2f}\n').format(station_location[0])
        text_label += ('Lat: {0:.2f}').format(station_location[1])
        self.station_annotation.set_text(text_label)

        return None
        
    def hover_station_annotation(self, event):
        """ Show or hide annotation for each station that is hovered. """

        # activate hover over map if any
        if event.inaxes == self.plot_axes['map']:
            if (hasattr(self.plot, 'stations_scatter')) and (self.lock_annotation == False):
                # lock annotation
                self.lock_annotation = True

                if event.button == None:
                    
                    is_contained, annotation_index = self.plot.stations_scatter.contains(event)
                    
                    if is_contained:
                        # update annotation if hovered
                        self.update_station_annotation(annotation_index)
                        self.station_annotation.set_visible(True)
                        self.annotation_visible = True
                    else:
                        # hide annotation if not hovered
                        if self.annotation_visible:
                            self.station_annotation.set_visible(False)
                            self.annotation_visible = False

                    # redraw points
                    self.figure.canvas.draw()
                    self.figure.canvas.flush_events()
                
                # unlock annotation 
                self.lock_annotation = False

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
                
                if scale_factor != 1:
                
                    # set new limits
                    self.plot_axes['map'].set_xlim([xdata - (xdata - current_xlim[0]) / scale_factor, 
                                                    xdata + (current_xlim[1] - xdata) / scale_factor])
                    self.plot_axes['map'].set_ylim([ydata - (ydata - current_ylim[0]) / scale_factor, 
                                                    ydata + (current_ylim[1] - ydata) / scale_factor])

                    # save map extent
                    self.read_instance.map_extent = self.plot_axes['map'].get_extent(crs=self.datacrs)

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

        if event.inaxes == self.plot_axes['map']:
            # block legend picker inside map
            self.lock_station_pick = False
            self.lock_legend_pick = True
        
        elif event.inaxes == self.plot_axes['legend']:
            # block stations picker inside legend
            self.lock_station_pick = True
            self.lock_legend_pick = False
        
        else:
            # block stations picker and legend picker
            self.lock_station_pick = True
            self.lock_legend_pick = True

        return None

    def legend_picker_func(self, event):
        """ Function to handle legend picker. """

        if self.lock_legend_pick == False:
            if self.lined:
                # lock legend pick
                self.lock_legend_pick = True
            
                # get plot lines
                legend_label = event.artist
                legend_label_text = legend_label.get_text().lower()
                visible = copy.deepcopy(self.lined[legend_label_text]['visible'])

                # iterate through plot types
                for plot_type in self.lined[legend_label_text]['lines_per_plot']:

                    print(len(self.lined[legend_label_text]['lines_per_plot'][plot_type]))
                    plot_lines = self.lined[legend_label_text]['lines_per_plot'][plot_type]
                
                    # change visibility of lines
                    for plot_line in plot_lines:
                        if visible:
                            plot_line.set_visible(False)
                        else:
                            plot_line.set_visible(True)

                # change font weight of label
                legend_label._fontproperties = self.legend.get_texts()[0]._fontproperties.copy()
                if visible:
                    legend_label.set_fontweight('regular')
                    self.lined[legend_label_text]['visible'] = False
                else:
                    legend_label.set_fontweight('bold')
                    self.lined[legend_label_text]['visible'] = True

                # redraw points
                self.figure.canvas.draw()
                self.figure.canvas.flush_events()
                
                # unlock legend pick 
                self.lock_legend_pick = False
    
        return None