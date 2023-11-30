""" Class to generate canvas """

import copy
import json
import os
import sys

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg \
        as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.path import Path
import matplotlib.style as mplstyle
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
from PyQt5 import QtCore, QtWidgets 
from weakref import WeakKeyDictionary

from .canvas_menus import SettingsMenu
from .dashboard_elements import ComboBox
from .dashboard_elements import set_formatting
from .dashboard_interactivity import LassoSelector, HoverAnnotation
from .dashboard_interactivity import legend_picker_func, picker_block_func, zoom_map_func
from .filter import DataFilter
from .plot import Plot
from .plot_aux import get_map_extent
from .plot_formatting import format_axis, harmonise_xy_lims_paradigm, set_axis_label, set_axis_title
from .plot_options import annotation, linear_regression, log_axes, log_validity, smooth
from .statistics import *
from .warnings import show_message

# make sure that we are using Qt5 backend with matplotlib
matplotlib.use('Qt5Agg')
register_matplotlib_converters()

# use matplotlib fast style: https://matplotlib.org/stable/users/explain/performance.html
mplstyle.use('fast')

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
basic_stats = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/experiment_bias_stats.json')))
formatting_dict = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/stylesheet.json')))
settings_dict = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/canvas_menus.json')))


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

        # add general plot characteristics to self
        for k, val in self.plot_characteristics_templates['general'].items():
            setattr(self, k, val)

        # initialise some key vars
        self.filter_data = None

        # initialise Plot class
        self.plot = Plot(read_instance=self.read_instance, canvas_instance=self)

        # initialise annotations
        self.annotations = dict()
        self.annotations_lock = dict()
        self.annotations_vline = dict()
        for plot_type in ['periodic', 'periodic-violin']:
            self.annotations[plot_type] = dict()
            self.annotations_lock[plot_type] = dict()
            self.annotations_vline[plot_type] = dict()

        # setup gridding of canvas
        self.gridspec = gridspec.GridSpec(self.gridspec_nrows, self.gridspec_ncols)
        self.gridspec.update(**self.gridspec_format)

        # define all possible plots
        self.all_plots = ['legend', 'map', 'timeseries', 'periodic-violin', 'periodic', 
                          'metadata', 'distribution', 'scatter', 'statsummary', 'boxplot',
                          'taylor']

        # define all possible plots in layout options
        self.layout_options = ['None', 'boxplot', 'distribution', 'metadata', 'periodic', 
                               'periodic-violin', 'scatter', 'statsummary', 'timeseries', 'taylor']

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

        # initialise plot elements
        self.plot_elements = {}
        self.current_plot_options = {}
        self.previous_plot_options = {}

        # initialise statsummary dict
        self.read_instance.current_statsummary_stats = {}
        self.read_instance.current_statsummary_stats['basic'] = {}
        self.read_instance.current_statsummary_stats['expbias'] = {}
        for periodic_cycle in ['None', 'Diurnal', 'Weekly', 'Monthly']:
            self.read_instance.current_statsummary_stats['basic'][periodic_cycle] = []
            self.read_instance.current_statsummary_stats['expbias'][periodic_cycle] = []

        # update plot characteristics for all plots, and initialise plot options per plot type 
        for plot_type in self.all_plots:
            # for plot types with zstat, initialise with default zstat (Mean)
            if plot_type in ['map', 'periodic']:
                self.plot.set_plot_characteristics([plot_type], zstat='Mean')
            else:
                self.plot.set_plot_characteristics([plot_type])
            self.current_plot_options[plot_type] = []
            self.previous_plot_options[plot_type] = []

        # create map, colorbar and legend plot axes
        self.plot_axes = {}
        self.plot_axes['map'] = self.figure.add_subplot(self.gridspec.new_subplotspec((2, 0), 
                                                        rowspan=44, colspan=42), projection=self.plotcrs)
        self.plot_axes['cb'] = self.figure.add_axes([0.0255, 0.536, 0.3794, 0.02])
        self.plot_axes['legend'] = self.figure.add_subplot(self.gridspec.new_subplotspec((0, 47), 
                                                           rowspan=8, colspan=53))

        # add settings menus
        self.generate_interactive_elements()

        # create rest of plot axes (default: timeseries, statsummary, distribution, periodic)
        # also show plot type buttons
        self.annotation_elements = []
        for position, plot_type in enumerate(self.read_instance.active_dashboard_plots):
            
            # update plot axis
            self.read_instance.update_plot_axis(self, position + 2, plot_type)

            # gather menu, save buttons and elements for plot type 
            for menu_button, save_button in zip(self.menu_buttons, self.save_buttons):
                menu_plot_type = menu_button.objectName().split('_menu')[0]
                if plot_type == 'periodic-violin':
                    plot_type = 'periodic_violin'
                # proceed once have objects for plot type
                if plot_type == menu_plot_type:
                    menu_button.show()
                    save_button.show()

        # show map buttons
        self.map_menu.buttons['settings_button'].show()
        self.map_menu.buttons['save_button'].show()

        # update layout fields
        self.read_instance.update_layout_fields(self)

        # initialise variable of valid station indices plotted on map as empty list
        self.active_map_valid_station_inds = np.array([], dtype=np.int64)

        # setup blocker for picker events
        self.figure.canvas.mpl_connect('axes_enter_event', lambda event: picker_block_func(self, event))

        # setup legend line selection
        self.figure.canvas.mpl_connect('pick_event', lambda event: legend_picker_func(self, event))

        # setup interactive lasso on map
        self.lasso_left = LassoSelector(self, self.plot_axes['map'], onselect=self.onlassoselect_left,
                                        useblit=False, props=self.lasso, button=[1])
        self.lasso_right = LassoSelector(self, self.plot_axes['map'], onselect=self.onlassoselect_right,
                                         useblit=False, props=self.lasso, button=[3])

        # setup station annotations
        self.annotations['map'] = HoverAnnotation(self, 'map', self.plot_axes['map'], self.plot_characteristics['map'])
        self.map_annotation_disconnect = False
        self.map_annotation_event = self.figure.canvas.mpl_connect('motion_notify_event', 
            self.annotations['map'].hover_map_annotation)

        # setup zoom on scroll wheel on map
        self.lock_zoom = False
        self.figure.canvas.mpl_connect('scroll_event', zoom_map_func)

        # format axes for map, legend and active_dashboard_plots
        for plot_type in ['map', 'legend'] + self.read_instance.active_dashboard_plots:
            format_axis(self, self.read_instance, self.plot_axes[plot_type], plot_type, 
                        self.plot_characteristics[plot_type], set_extent=True)

        # create covers to hide parts of canvas when updating / plotting
        self.canvas_cover = set_formatting(QtWidgets.QWidget(self), formatting_dict['canvas_cover'])
        self.top_right_canvas_cover = set_formatting(QtWidgets.QWidget(self), formatting_dict['canvas_cover'])
        self.top_right_canvas_cover.hide() 
        self.lower_canvas_cover = set_formatting(QtWidgets.QWidget(self), formatting_dict['canvas_cover'])
        self.lower_canvas_cover.hide()
        # place partial canvas covers below map elements 
        for element in self.map_elements:
            element.raise_()

    def update_MPL_canvas(self):
        """ Function that updates MPL canvas upon clicking
            the 'READ' button, and when colocating data.
        """

        # reset relative index lists of selected station on map as empty lists
        self.previous_relative_selected_station_inds = np.array([], dtype=np.int64)
        self.relative_selected_station_inds = np.array([], dtype=np.int64)
        self.absolute_selected_station_inds = np.array([], dtype=np.int64)

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

        # uncover map, but hide plotting axes
        self.canvas_cover.hide()
        self.top_right_canvas_cover.show() 
        self.lower_canvas_cover.show()

        # draw changes
        self.figure.canvas.draw_idle()

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
                    self.figure.canvas.draw_idle()

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
            self.figure.canvas.draw_idle()

        return None
    
    def handle_statistic_mode_update(self):
        """ Function that handles the update of the MPL canvas
            when we change the statistical calculation mode
        """

        if not self.read_instance.block_MPL_canvas_updates:

            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # get statistic
            self.read_instance.selected_statistic_mode = self.read_instance.cb_statistic_mode.currentText()
            
            # update timeseries statistic if it is different than the selected statistic aggregation
            if ((self.read_instance.selected_statistic_mode == 'Spatial|Temporal')
                and (self.timeseries_stat.currentText() != self.read_instance.selected_statistic_aggregation)):
                self.read_instance.block_MPL_canvas_updates = True
                self.timeseries_stat.setCurrentText(self.read_instance.selected_statistic_aggregation)
                self.read_instance.block_MPL_canvas_updates = False

            # update statistic in memory
            self.read_instance.statistic_mode = self.read_instance.selected_statistic_mode 

            # update associated plots with selected stations
            self.update_associated_active_dashboard_plots()
            
            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()   

        return None

    def handle_statistic_aggregation_update(self):
        """ Function that handles the update of the MPL canvas
            when we change the aggregation statistic
        """

        if not self.read_instance.block_MPL_canvas_updates:

            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # get statistic
            self.read_instance.selected_statistic_aggregation = self.read_instance.cb_statistic_aggregation.currentText()
            
            # update timeseries statistic if it is different than the selected statistic aggregation
            if ((self.read_instance.selected_statistic_mode == 'Spatial|Temporal')
                and (self.timeseries_stat.currentText() != self.read_instance.selected_statistic_aggregation)):
                self.read_instance.block_MPL_canvas_updates = True
                self.timeseries_stat.setCurrentText(self.read_instance.selected_statistic_aggregation)
                self.read_instance.block_MPL_canvas_updates = False

            # update statistic in memory
            self.read_instance.statistic_aggregation = self.read_instance.selected_statistic_aggregation 

            # update associated plots with selected stations
            self.update_associated_active_dashboard_plots()
            
            # draw changes
            self.figure.canvas.draw_idle()

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

            # if only have < 2 data arrays in memory, no colocation is possible,
            # therefore set colocation to be False, and return
            invalid = False
            if self.read_instance.data_labels == None:
                invalid = True
            elif len(self.read_instance.data_labels) == 1:
                invalid = True
            if invalid:
                msg = 'Load experiments before activating the temporal colocation'
                show_message(self.read_instance, msg)
                self.read_instance.temporal_colocation = False
                self.read_instance.block_MPL_canvas_updates = True
                self.read_instance.ch_colocate.setCheckState(QtCore.Qt.Unchecked)
                self.read_instance.block_MPL_canvas_updates = False
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
            self.handle_timeseries_statistic_update()
            self.handle_periodic_statistic_update()
            self.handle_statsummary_statistics_update()
            self.handle_statsummary_cycle_update()
            self.handle_statsummary_periodic_aggregation_update()
            self.handle_statsummary_periodic_mode_update()

            # self.handle_taylor_correlation_statistic_update()
            self.read_instance.block_MPL_canvas_updates = False

            # update layout fields
            self.read_instance.update_layout_fields(self)

            # update plotted map z statistic
            self.update_map_z_statistic()

            # update associated plots with selected stations
            self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw_idle()

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
        if self.map_z2.currentText() == '':
            zstat = get_z_statistic_comboboxes(base_zstat)
        else:
            zstat = get_z_statistic_comboboxes(base_zstat, bias=True)

        # plot map for zstat --> updating active map valid station indices and setting up plot picker
        self.plot.make_map(self.plot_axes['map'], self.read_instance.networkspeci, self.plot_characteristics['map'], 
                           zstat=zstat, labela=self.map_z1.currentText(), labelb=self.map_z2.currentText())
        
        # update absolute selected plotted station indices with respect to new active map valid station indices
        self.absolute_selected_station_inds = np.array(
            [np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind in
             self.relative_selected_station_inds if selected_ind in self.active_map_valid_station_inds],
            dtype=np.int64)

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
            self.relative_selected_station_inds = np.array([], dtype=np.int64)
            self.absolute_selected_station_inds = np.array([], dtype=np.int64)
            self.absolute_non_selected_station_inds = np.array([], dtype=np.int64)

        # else, if any of the currently selected stations are not in the current active map
        # valid station indices --> unselect selected stations (and associated plots)
        # also uncheck select all/intersect/extent checkboxes
        else:

            if not np.all(np.in1d(self.relative_selected_station_inds, self.active_map_valid_station_inds)):
                
                # unselect all/intersect/extent checkboxes
                self.unselect_map_checkboxes()
                
                # reset relative/absolute selected station indices to be empty lists
                self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)
                self.relative_selected_station_inds = np.array([], dtype=np.int64)
                self.absolute_selected_station_inds = np.array([], dtype=np.int64)

            # get absolute non-selected station inds
            self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                                 self.absolute_selected_station_inds))[0]

            # create 2D numpy array of plotted station coordinates
            self.map_points_coordinates = np.vstack((self.read_instance.station_longitudes[self.read_instance.networkspeci][self.active_map_valid_station_inds], 
                                                     self.read_instance.station_latitudes[self.read_instance.networkspeci][self.active_map_valid_station_inds])).T

            # generate colourbar
            generate_colourbar(self.read_instance, [self.plot_axes['map']], [self.plot_axes['cb']], zstat, 
                               self.plot_characteristics['map'], self.read_instance.species[0])

        # update plot options
        self.update_plot_options(plot_types=['map'])

        # redraw plot (needed to update plotted colours before update_map_station_selection)
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
        set_axis_title(self.read_instance, self.plot_axes['map'], axis_title_label, self.plot_characteristics['map'])
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
        
        # redraw plot
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def update_associated_active_dashboard_plot(self, plot_type):
        """ Function that updates a plot associated with selected stations on map. """

        if hasattr(self, 'relative_selected_station_inds'):
            if len(self.relative_selected_station_inds) > 0:

                # clear all previously plotted artists for plot type
                self.remove_axis_elements(self.plot_axes[plot_type], plot_type)

                # get relevant axis
                ax = self.plot_axes[plot_type]

                # get options defined to configure plot 
                plot_options = copy.deepcopy(self.current_plot_options[plot_type])

                # get plotting function for specific plot
                if plot_type == 'statsummary':
                    func = getattr(self.plot, 'make_table')
                else:
                    func = getattr(self.plot, 'make_{}'.format(plot_type.split('-')[0]))

                # set ylabel for periodic plot
                if plot_type == 'periodic':
                    # get currently selected periodic statistic name
                    base_zstat = self.periodic_stat.currentText()
                    if 'bias' in plot_options:
                        zstat = get_z_statistic_comboboxes(base_zstat, bias=True)
                    else:
                        zstat = get_z_statistic_comboboxes(base_zstat)
                    
                    # get zstat information 
                    zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat) 
                    
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
                elif plot_type == 'statsummary':
                    xlabel = ''
                    ylabel = ''
                
                # create structure to store data for Taylor diagram
                elif plot_type == 'taylor':
                    # get r or r2 as correlation statistic
                    corr_stat = self.plot_characteristics[plot_type]['corr_stat']
                    relevant_zstats = [corr_stat, "StdDev"]
                    stats_df = {relevant_zstat:[] for relevant_zstat in relevant_zstats}

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

                # call function to update plot
                # periodic plot
                if plot_type == 'periodic':
                    func(ax, self.read_instance.networkspeci, self.read_instance.data_labels, 
                         self.plot_characteristics[plot_type], zstat=zstat, plot_options=plot_options)
                # make statsummary plot
                elif plot_type == 'statsummary':
                    if 'bias' in plot_options:
                        relevant_zstats = self.active_statsummary_stats['expbias']
                    else:
                        relevant_zstats = self.active_statsummary_stats['basic']
                    
                    func(ax, self.read_instance.networkspeci, self.read_instance.data_labels, 
                         self.plot_characteristics[plot_type], 
                         zstats=relevant_zstats, statsummary=True, plot_options=plot_options)                
                # other plots
                else:
                    func(ax, self.read_instance.networkspeci, self.read_instance.data_labels, 
                         self.plot_characteristics[plot_type], plot_options=plot_options)

                # reset axes limits (harmonising across subplots for periodic plots) 
                if plot_type == 'scatter':
                    harmonise_xy_lims_paradigm(self, self.read_instance, ax, plot_type, 
                                               self.plot_characteristics[plot_type], plot_options, relim=True)
                elif plot_type != 'taylor':
                    harmonise_xy_lims_paradigm(self, self.read_instance, ax, plot_type, 
                                               self.plot_characteristics[plot_type], plot_options, relim=True, autoscale=True)

                # skip setting axes labels for Taylor diagram
                if plot_type != 'taylor':
                    # set xlabel
                    set_axis_label(ax, 'x', xlabel, self.plot_characteristics[plot_type])
                    # set ylabel
                    set_axis_label(ax, 'y', ylabel, self.plot_characteristics[plot_type])

                # reset navigation toolbar stack for plot
                self.reset_ax_navigation_toolbar_stack(ax)

                # update plot options
                self.update_plot_options(plot_types=[plot_type])

    def update_associated_active_dashboard_plots(self):
        """ Function that updates all plots associated with selected stations on map. """

        #start = time.time()

        # update dashboard plots
        if hasattr(self, 'relative_selected_station_inds'):
            # have no selected stations, so clear all previously plotted artists from selected station plots
            # cover plotting axes also
            if len(self.relative_selected_station_inds) == 0:      
                for plot_type in self.read_instance.active_dashboard_plots:
                    self.remove_axis_elements(self.plot_axes[plot_type], plot_type)
                self.top_right_canvas_cover.show() 
                self.lower_canvas_cover.show()

            elif len(self.relative_selected_station_inds) > 0:
                
                # get selected station data
                get_selected_station_data(read_instance=self.read_instance, canvas_instance=self, 
                                          networkspecies=[self.read_instance.networkspeci])

                # iterate through active_dashboard_plots
                for plot_type_ii, plot_type in enumerate(self.read_instance.active_dashboard_plots):

                    #plot_start = time.time()

                    # if there are no temporal resolutions (only yearly), skip periodic plots
                    if ((plot_type in ['periodic', 'periodic-violin']) and 
                        (not self.read_instance.relevant_temporal_resolutions)):
                        msg = 'It is not possible to make periodic plots using annual resolution data.'
                        show_message(self.read_instance, msg)
                        self.read_instance.handle_layout_update('None', sender=plot_type_ii+2) 
                        continue
                    
                    # if temporal colocation is turned off or there are no experiments, skip scatter plot
                    if plot_type in ['scatter', 'taylor']:
                        if ((not self.read_instance.temporal_colocation) 
                            or ((self.read_instance.temporal_colocation) and (len(self.read_instance.experiments) == 0))):
                            if (not self.read_instance.temporal_colocation):
                                msg = f'It is not possible to make {plot_type} plots without activating the temporal colocation.'
                            else:
                                msg = f'It is not possible to make {plot_type} plots without loading experiments.'
                            show_message(self.read_instance, msg)
                            self.read_instance.handle_layout_update('None', sender=plot_type_ii+2)
                            continue

                    # update plot
                    self.update_associated_active_dashboard_plot(plot_type)

                    #print('{}: {}'.format(plot_type, time.time()-plot_start))

                # un-hide plotting axes
                self.top_right_canvas_cover.hide() 
                self.lower_canvas_cover.hide()

            # update map plot options
            self.update_plot_options(plot_types=['map'])

            #print('TOTAL CANVAS UPDATE: {}'.format(time.time()-start))

    def update_experiment_domain_edges(self):
        """ Function that plots grid domain edges of experiments in memory. """

        # remove grid domain polygon if previously plotted
        self.remove_axis_objects(self.plot_axes['map'].patches, types_to_remove=[matplotlib.patches.Polygon])

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
            # force them to be observations label and first basic z statistic respectively
            if selected_z_stat == '':
                selected_z_stat = self.read_instance.basic_z_stats[0]
            if hasattr(self.read_instance, 'map_z'):
                if self.read_instance.map_z in self.read_instance.basic_z_stats:
                    selected_z_stat = self.read_instance.map_z
            if selected_z1_array == '':
                selected_z1_array = copy.deepcopy(self.read_instance.observations_data_label)

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

    def handle_timeseries_statistic_update(self):
        """ Function that handles update of plotted timeseries statistic
            upon interaction with timeseries statistic combobox.
        """

        if not self.read_instance.block_config_bar_handling_updates:

            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # set variable that blocks configuration bar handling updates until all changes
            # to the timeseries statistic combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected statistic
            zstat = self.timeseries_stat.currentText()

            # update timeseries statistics
            available_timeseries_stats = ['Mean', 'Median', 'p1', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95', 'p99']

            # if base_zstat is empty string, it is because fields are being initialised for the first time
            if zstat == '':
                # set timeseries stat to be first available stat
                zstat = available_timeseries_stats[0]

            # update timeseries statistic combobox (clear, then add items)
            self.timeseries_stat.clear()
            self.timeseries_stat.addItems(available_timeseries_stats)

            # maintain currently selected timeseries statistic (if exists in new item list)
            if zstat in available_timeseries_stats:
                self.timeseries_stat.setCurrentText(zstat)

            # update aggregation statistic if it is different than timeseries statistic
            if ((self.read_instance.statistic_mode == 'Spatial|Temporal')
                and (zstat != self.read_instance.cb_statistic_aggregation.currentText())):
                self.read_instance.block_MPL_canvas_updates = True
                self.read_instance.cb_statistic_aggregation.setCurrentText(zstat)
                self.read_instance.block_MPL_canvas_updates = False

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted timeseries statistic
            if not self.read_instance.block_MPL_canvas_updates:
                # update selected data on all active plots
                if self.read_instance.statistic_mode == 'Spatial|Temporal':
                    self.update_associated_active_dashboard_plots()

                # update selected data on timeseries plot
                elif self.read_instance.statistic_mode in ['Temporal|Spatial', 'Flattened']:
                    if len(self.station_inds[self.read_instance.networkspeci]) >= 1:
                        # update timeseries data
                        timeseries_stat = self.timeseries_stat.currentText()
                        aggregated_data = aggregation(self.read_instance.data_array, timeseries_stat, axis=1)
                        self.selected_station_data[self.read_instance.networkspeci]['timeseries'] = pd.DataFrame(
                            aggregated_data.T, 
                            columns=self.selected_station_data_labels[self.read_instance.networkspeci], 
                            index=self.read_instance.time_index)

                        # update plot                                                                         
                        self.update_associated_active_dashboard_plot('timeseries')

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()

        return None

    def handle_periodic_statistic_update(self):
        """ Function that handles update of plotted periodic statistic
            upon interaction with periodic statistic combobox.
        """

        if not self.read_instance.block_config_bar_handling_updates:

            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

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
                self.update_associated_active_dashboard_plot('periodic')

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()

        return None

    def handle_taylor_correlation_statistic_update(self):
        
        if not self.read_instance.block_config_bar_handling_updates:

            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # update taylor diagram correlation statistic combobox
            # set variable that blocks configuration bar handling updates until all
            # changes to the statistic combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected items
            corr_stat = self.taylor_corr_stat.currentText()

            # get available stats
            available_corr_stats = ['r', 'r2']

            # if correlation stat is empty string, it is because fields are being initialised for the first time
            if corr_stat == '':
                # set stat to be the one in plot characteristics
                corr_stat = self.plot_characteristics['taylor']['corr_stat']

            # update statistic combobox (clear, then add items)
            self.taylor_corr_stat.clear()
            self.taylor_corr_stat.addItems(available_corr_stats)

            # maintain currently selected statistic
            self.taylor_corr_stat.setCurrentText(corr_stat)

            # update dictionary
            self.plot_characteristics['taylor']['corr_stat'] = corr_stat

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted taylor diagram statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot('taylor')

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()

        return None

    def get_active_statsummary_stats(self, statistic_type):
        """ Get active statistics from dictionary of statsummary statistics in list """
        
        active_statsummary_stats = copy.deepcopy(self.read_instance.current_statsummary_stats[statistic_type])
        active_statsummary_stats = [stat for sublist in 
                                    list(active_statsummary_stats.values()) 
                                    for stat in sublist]

        return active_statsummary_stats
    
    def check_statsummary_stats(self):
        """ Function that checks the statistics in the statsummary statistic combobox. """

        # get stats to check for the selected periodic cycle
        periodic_cycle = self.statsummary_cycle.currentText()
        if periodic_cycle == '':
            periodic_cycle = 'None'
        plot_options = self.current_plot_options['statsummary']
        statistic_type = 'basic' if 'bias' not in plot_options else 'expbias'
        if 'bias' in plot_options:
            items = list(copy.deepcopy(self.read_instance.basic_and_bias_z_stats))
        else:
            items = list(copy.deepcopy(self.read_instance.basic_z_stats))
        if periodic_cycle != 'None':
            items = [stat + '-' + periodic_cycle.lower() for stat in items]
        self.statsummary_stat.clear()
        self.statsummary_stat.addItems(items)
        checked_options = copy.deepcopy(self.read_instance.current_statsummary_stats[statistic_type][periodic_cycle])
        checked_options = [option.split('_bias')[0] if '_bias' in option else option 
                           for option in checked_options]
        checked_options_in_items = list(set(checked_options) & set(items))
        
        # check stats in combobox
        if checked_options_in_items:
            for checked_option in checked_options_in_items:
                index = items.index(checked_option)
                self.statsummary_stat.model().item(index).setCheckState(QtCore.Qt.Checked)
        # leave empty
        else:
            self.statsummary_stat.lineEdit().setText('')

    def handle_statsummary_statistics_update(self):
        """ Function that handles update of plotted statsummary statistics
            upon interaction with statistic comboboxes.
        """

        if not self.read_instance.block_config_bar_handling_updates:
            
            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # set variable that blocks configuration bar handling updates until all changes
            # to the statsummary statistics combobox are made
            self.read_instance.block_config_bar_handling_updates = True
            
            # get all possible stats
            plot_options = self.current_plot_options['statsummary']
            statistic_type = 'basic' if 'bias' not in plot_options else 'expbias'
            
            # initialise stats
            if not hasattr(self, 'active_statsummary_stats'):
                
                # get initial stats from plot characteristics
                periodic_cycle = 'None'
                self.read_instance.current_statsummary_stats['basic']['None'] = self.plot_characteristics['statsummary']['basic']
                self.read_instance.current_statsummary_stats['expbias']['None'] = self.plot_characteristics['statsummary']['experiment_bias']
                self.active_statsummary_stats = {'basic': self.get_active_statsummary_stats('basic'),
                                                 'expbias': self.get_active_statsummary_stats('expbias')}

                # check stats for the selected periodic cycle
                self.check_statsummary_stats()

            # get stats from selection
            else:
                # save previous stats in list
                previous_active_statsummary_stats = copy.deepcopy(self.active_statsummary_stats[statistic_type])

                # remove bias from options to get correct active stats
                if statistic_type == 'expbias':
                    previous_active_statsummary_stats = [option.split('_bias')[0] if '_bias' in option else option 
                                                         for option in previous_active_statsummary_stats]
                
                # update stats
                periodic_cycle = self.statsummary_cycle.currentText()
                self.read_instance.current_statsummary_stats[statistic_type][periodic_cycle] = copy.deepcopy(self.statsummary_stat.currentData())

                # save current stats in list
                current_active_statsummary_stats = copy.deepcopy(self.get_active_statsummary_stats(statistic_type))

                # check stats for the selected periodic cycle
                self.check_statsummary_stats()

                # get active stats
                current_not_previous = list(set(current_active_statsummary_stats).difference(previous_active_statsummary_stats))
                previous_not_current = list(set(previous_active_statsummary_stats).difference(current_active_statsummary_stats))
                self.active_statsummary_stats[statistic_type] = copy.deepcopy(previous_active_statsummary_stats)
                if current_not_previous:
                    for stat in current_not_previous:
                        self.active_statsummary_stats[statistic_type].append(stat)
                if previous_not_current:
                    for stat in previous_not_current:
                        self.active_statsummary_stats[statistic_type].remove(stat)

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted statsummary statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot('statsummary')

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()

    def handle_statsummary_cycle_update(self):

        if not self.read_instance.block_config_bar_handling_updates:
            
            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # set variable that blocks configuration bar handling updates until all changes
            # to the statsummary periodic cycle combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected cycle
            periodic_cycle = self.statsummary_cycle.currentText()
            
            # update periodi cycles
            available_cycles = ['None', 'Diurnal', 'Weekly', 'Monthly']

            # if cycle is empty string, it is because fields are being initialised for the first time
            if periodic_cycle == '':
                # set statsummary cycle to be None
                periodic_cycle = available_cycles[0]

            # update statsummary cycle combobox (clear, then add items)
            self.statsummary_cycle.clear()
            self.statsummary_cycle.addItems(available_cycles)

            # maintain currently selected periodic cycle (if exists in new item list)
            if periodic_cycle in available_cycles:
                self.statsummary_cycle.setCurrentText(periodic_cycle)

            # check stats for the selected periodic cycle
            self.check_statsummary_stats()

            self.read_instance.block_config_bar_handling_updates = False

            # update plotted statsummary statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot('statsummary')

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()

        return None
    
    def handle_statsummary_periodic_aggregation_update(self):
        """ Function that handles update of plotted statsummary periodic aggregation statistic
            upon interaction with combobox.
        """

        if not self.read_instance.block_config_bar_handling_updates:

            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # set variable that blocks configuration bar handling updates until all changes
            # to the statsummary periodic aggregation combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected statistic
            stat = self.statsummary_periodic_aggregation.currentText()

            # update periodic aggregation
            available_stats = ['Mean', 'Median', 'p1', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95', 'p99']

            # if stat is empty string, it is because fields are being initialised for the first time
            if stat == '':
                # set periodic stat to be first available stat
                stat = available_stats[0]

            # update periodic aggregation combobox (clear, then add items)
            self.statsummary_periodic_aggregation.clear()
            self.statsummary_periodic_aggregation.addItems(available_stats)

            # maintain currently selected statistic (if exists in new item list)
            if stat in available_stats:
                self.statsummary_periodic_aggregation.setCurrentText(stat)
                self.read_instance.periodic_statistic_aggregation = copy.deepcopy(stat)

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted statsummary statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot('statsummary')

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()

        return None
    
    def handle_statsummary_periodic_mode_update(self):
        """ Function that handles update of plotted statsummary periodic aggregation mode
            upon interaction with combobox.
        """

        if not self.read_instance.block_config_bar_handling_updates:

            # update mouse cursor to a waiting cursor
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

            # set variable that blocks configuration bar handling updates until all changes
            # to the statsummary periodic mode combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected mode
            mode = self.statsummary_periodic_mode.currentText()

            # update periodic mode
            available_modes = ['Independent', 'Cycle']

            # if mode is empty string, it is because fields are being initialised for the first time
            if mode == '':
                # set mode to be first available stat
                mode = available_modes[0]

            # update periodic mode combobox (clear, then add items)
            self.statsummary_periodic_mode.clear()
            self.statsummary_periodic_mode.addItems(available_modes)

            # maintain currently selected mode (if exists in new item list)
            if mode in available_modes:
                self.statsummary_periodic_mode.setCurrentText(mode)
                self.read_instance.periodic_statistic_mode = copy.deepcopy(mode)

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted statsummary statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot('statsummary')

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            QtWidgets.QApplication.restoreOverrideCursor()

        return None
    
    def remove_axis_objects(self, ax_elements, elements_to_skip=[], types_to_remove=[]):
        """ Remove objects (artists, lines, collections, patches) from axis. """
        
        # put elements from dicts in lists for periodic plots
        if isinstance(elements_to_skip, dict):
            elements_to_skip = list(elements_to_skip.values())
        
        # it is not possible to remove the elements directly with index from a FeatureArtist
        # first we need to find the indices of the elements to be removed and remove them 
        # by index one by one (sorting is needed)
        inds_to_remove = []
        for element_ii, element in enumerate(ax_elements):
            if element not in elements_to_skip:
                if len(types_to_remove) > 0:
                    if isinstance(element, tuple(types_to_remove)):
                        inds_to_remove.append(element_ii)
                else:
                    inds_to_remove.append(element_ii)

        # remove
        for element_ii in sorted(inds_to_remove, reverse=True):
            ax_elements[element_ii].remove()

        return None
    
    def remove_axis_elements(self, ax, plot_type):
        """ Remove all plotted axis elements."""
       
        # get appropriate axes for nested axes
        axs_to_remove = []
        if isinstance(ax, dict):
            for relevant_temporal_resolution, sub_ax in ax.items():
                axs_to_remove.append(sub_ax)
        else:
            if plot_type == 'taylor':
                axs_to_remove.append(self.plot.taylor_polar_relevant_axis)
            axs_to_remove.append(ax)

        # iterate through axes
        for ax_to_remove in axs_to_remove:

            # remove all plotted axis elements
            if plot_type == 'legend':
                leg = ax_to_remove.get_legend()
                if leg:
                    leg.remove()

            elif plot_type == 'map':
                self.remove_axis_objects(ax_to_remove.artists, types_to_remove=[AnchoredOffsetbox])
                self.remove_axis_objects(ax_to_remove.collections, types_to_remove=[matplotlib.collections.PathCollection])
                # # TODO: Put line collection back into place when we turn on the auto_update in gridlines
                # self.remove_axis_objects(ax_to_remove.collections, types_to_remove=[matplotlib.collections.PathCollection],
                #                                                                     matplotlib.collections.LineCollection])

            elif plot_type == 'cb':
                for objects in [ax_to_remove.artists, ax_to_remove.collections]:
                    self.remove_axis_objects(objects)

            elif plot_type == 'timeseries':
                for objects in [ax_to_remove.lines, ax_to_remove.artists]:
                    self.remove_axis_objects(objects, elements_to_skip=[self.annotations_vline['timeseries']])

            elif plot_type == 'periodic':
                for objects in [ax_to_remove.lines, ax_to_remove.artists]:
                    self.remove_axis_objects(objects, elements_to_skip=self.annotations_vline['periodic'].values())

            elif plot_type == 'periodic-violin':
                for objects in [ax_to_remove.lines, ax_to_remove.artists, ax_to_remove.collections]:
                    self.remove_axis_objects(objects, elements_to_skip=self.annotations_vline['periodic-violin'].values())

            elif plot_type == 'metadata':
                self.remove_axis_objects(ax_to_remove.texts)

            elif plot_type == 'distribution':
                for objects in [ax_to_remove.lines, ax_to_remove.artists]:
                    self.remove_axis_objects(objects, elements_to_skip=[self.annotations_vline['distribution']])

            elif plot_type == 'statsummary':
                self.remove_axis_objects(ax_to_remove.tables)

            elif plot_type in ['taylor', 'boxplot', 'scatter']:
                for objects in [ax_to_remove.lines, ax_to_remove.artists]:
                    self.remove_axis_objects(objects)

        # remove tracked plot elements
        if plot_type in self.plot_elements:
            self.plot_elements[plot_type]['absolute'] = {}
            if 'bias' in self.plot_elements[plot_type]:
                del self.plot_elements[plot_type]['bias']

        # hide annotation boxes and lines
        for element in self.annotation_elements:
            if isinstance(element, dict):
                for val in element.values():
                    val.set_visible(False)
            else:
                element.set_visible(False)

        return None

    def update_plot_options(self, plot_types):
        """ Uncheck checked boxes in plot configuration options under menus and
            reapply check with new data. This can be done for all currently active plot types, 
            or just one specific type.
        """

        if not isinstance(plot_types, list):
            plot_types = [plot_types]

        for plot_type in plot_types:
            all_plot_options = self.plot_characteristics[plot_type]['plot_options']
            checked_options = self.current_plot_options[plot_type]
            if plot_type == 'periodic-violin':
                plot_type = 'periodic_violin'
            cb_options = getattr(self, plot_type + '_options')

            # uncheck all options
            for option in all_plot_options:
                index = all_plot_options.index(option)
                self.read_instance.block_MPL_canvas_updates = True
                cb_options.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                self.read_instance.block_MPL_canvas_updates = False
            
            # check selected options
            for checked_option_ii, checked_option in enumerate(checked_options):
                index = all_plot_options.index(checked_option)
                if checked_option_ii < (len(checked_options)-1):
                    self.read_instance.block_MPL_canvas_updates = True
                    cb_options.model().item(index).setCheckState(QtCore.Qt.Checked)
                    self.read_instance.block_MPL_canvas_updates = False
                else:
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
                self.relative_selected_station_inds = np.array([], dtype=np.int64)

            # update absolute selected station indices (indices relative to plotted stations on map)
            self.absolute_selected_station_inds = np.arange(len(self.relative_selected_station_inds), dtype=np.int64)

            # get absolute non-selected station inds
            self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                                 self.absolute_selected_station_inds))[0]

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
                self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw_idle()

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
                self.relative_selected_station_inds = np.array([], dtype=np.int64)
                self.absolute_selected_station_inds = np.array([], dtype=np.int64)
                self.absolute_non_selected_station_inds = np.arange(len(self.relative_selected_station_inds),
                                                                    dtype=np.int64)

            # else, if checkbox is checked then select all stations which intersect with all loaded experiment domains
            elif check_state == QtCore.Qt.Checked:

                # if have only observations loaded into memory, select all plotted stations
                if len(self.read_instance.data_labels) == 1:
                    self.relative_selected_station_inds = copy.deepcopy(self.active_map_valid_station_inds)
                    self.absolute_selected_station_inds = np.arange(len(self.relative_selected_station_inds),
                                                                    dtype=np.int64)
                    self.absolute_non_selected_station_inds = np.array([], dtype=np.int64)
                # else, define list of lists to get intersection between (active_map_valid_station_inds, 
                # and valid station indices associated with each loaded experiment array)
                else:
                    intersect_lists = [self.active_map_valid_station_inds]
                    for data_label in self.read_instance.data_labels:
                        if data_label != self.read_instance.observations_data_label:
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
                                  in self.relative_selected_station_inds], dtype=np.int64)
                    
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
            self.figure.canvas.draw_idle()

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
            self.read_instance.map_extent = get_map_extent(self)

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
                self.relative_selected_station_inds = np.array([], dtype=np.int64)

            # get absolute selected station indices (indices relative to plotted stations on map)
            self.absolute_selected_station_inds = \
                np.array([np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind
                            in self.relative_selected_station_inds], dtype=np.int64)

            # get absolute unselected station indices
            self.absolute_non_selected_station_inds = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                                 self.absolute_selected_station_inds))[0]

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
                self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw_idle()

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
            self.read_instance.map_extent = get_map_extent(self)
            tolerance = np.average([self.read_instance.map_extent[1]-self.read_instance.map_extent[0],
                                    self.read_instance.map_extent[3]-self.read_instance.map_extent[2]]) / 100.0
            point_coordinates = lasso_path.vertices[0:1,:]
            sub_abs_vals = np.abs(self.map_points_coordinates[None,:,:] - point_coordinates[:,None,:])
            self.absolute_selected_station_inds = np.arange(len(self.active_map_valid_station_inds))[np.all(np.any(sub_abs_vals<=tolerance,axis=0),axis=1)]
            
            # if more than 1 point selected, limit this to be just nearest point
            if len(self.absolute_selected_station_inds) > 1:
                self.absolute_selected_station_inds = np.array([self.absolute_selected_station_inds[np.argmin(np.sum(sub_abs_vals[0,self.absolute_selected_station_inds,:],axis=1))]], dtype=np.int64)

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
        self.figure.canvas.draw_idle()

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
            self.read_instance.map_extent = get_map_extent(self)
            tolerance = np.average([self.read_instance.map_extent[1]-self.read_instance.map_extent[0],self.read_instance.map_extent[3]-self.read_instance.map_extent[2]]) / 100.0
            point_coordinates = lasso_path.vertices[0:1,:]
            sub_abs_vals = np.abs(self.map_points_coordinates[None,:,:] - point_coordinates[:,None,:])
            absolute_selected_station_inds = np.arange(len(self.active_map_valid_station_inds))[np.all(np.any(sub_abs_vals<=tolerance,axis=0),axis=1)]
            # if more than 1 point selected, limit this to be just nearest point
            if len(absolute_selected_station_inds) > 1:
                absolute_selected_station_inds = np.array([absolute_selected_station_inds[np.argmin(np.sum(sub_abs_vals[0,absolute_selected_station_inds,:],axis=1))]], dtype=np.int64)

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
        self.figure.canvas.draw_idle()

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
        self.read_instance.cb_position_2.setToolTip('Select plot type in top right position')
        self.read_instance.cb_position_2.currentTextChanged.connect(self.read_instance.handle_layout_update)
        # self.read_instance.cb_position_2.hide()

        # add position 3 plot selector
        self.read_instance.cb_position_3 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.read_instance.cb_position_3.setToolTip('Select plot type in bottom left position')
        self.read_instance.cb_position_3.currentTextChanged.connect(self.read_instance.handle_layout_update)
        # self.read_instance.cb_position_3.hide()

        # add position 4 plot selector
        self.read_instance.cb_position_4 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.read_instance.cb_position_4.setToolTip('Select plot type in bottom centre position')
        self.read_instance.cb_position_4.currentTextChanged.connect(self.read_instance.handle_layout_update)
        # self.read_instance.cb_position_4.hide()

        # add position 5 plot selector
        self.read_instance.cb_position_5 = set_formatting(ComboBox(self), formatting_dict['combobox_menu'])
        self.read_instance.cb_position_5.setToolTip('Select plot type in bottom right position')
        self.read_instance.cb_position_5.currentTextChanged.connect(self.read_instance.handle_layout_update)
        # self.read_instance.cb_position_5.hide()

        # MAP SETTINGS MENU #
        # create map settings menu
        self.map_menu = SettingsMenu(plot_type='map', canvas_instance=self)
        self.map_options = self.map_menu.checkable_comboboxes['options']
        self.map_elements = self.map_menu.get_elements()

        # get stats
        self.map_z_stat = self.map_menu.comboboxes['z_stat']
        self.map_z1 = self.map_menu.comboboxes['z1']
        self.map_z2 = self.map_menu.comboboxes['z2']

        # get sliders and update values
        self.map_markersize_unsel_sl = self.map_menu.sliders['markersize_unsel_sl']
        self.map_markersize_unsel_sl.setValue(self.plot_characteristics['map']['marker_unselected']['s'])
        self.map_opacity_unsel_sl = self.map_menu.sliders['opacity_unsel_sl']
        self.map_opacity_unsel_sl.setValue(self.plot_characteristics['map']['marker_unselected']['alpha']*10)
        self.map_markersize_sel_sl = self.map_menu.sliders['markersize_sel_sl']
        self.map_markersize_sel_sl.setValue(self.plot_characteristics['map']['marker_selected']['s'])
        self.map_opacity_sel_sl = self.map_menu.sliders['opacity_sel_sl']
        self.map_opacity_sel_sl.setValue(self.plot_characteristics['map']['marker_selected']['alpha']*10)
        
        # get map interactive dictionary
        self.interactive_elements['map'] = {'hidden': True,
                                            'markersize_sl': [self.map_markersize_unsel_sl, 
                                                              self.map_markersize_sel_sl],
                                            'opacity_sl': [self.map_opacity_unsel_sl, 
                                                           self.map_opacity_sel_sl]
                                            }

        # TIMESERIES PLOT SETTINGS MENU #
        # create timeseries settings menu
        self.timeseries_menu = SettingsMenu(plot_type='timeseries', canvas_instance=self)
        self.timeseries_options = self.timeseries_menu.checkable_comboboxes['options']
        self.timeseries_elements = self.timeseries_menu.get_elements()

        # get stats
        self.timeseries_stat = self.timeseries_menu.comboboxes['stat']

        # get sliders and update values
        self.timeseries_markersize_sl = self.timeseries_menu.sliders['markersize_sl']
        self.timeseries_markersize_sl.setMaximum(self.plot_characteristics['timeseries']['plot']['markersize']*10)
        self.timeseries_markersize_sl.setValue(self.plot_characteristics['timeseries']['plot']['markersize'])
        self.timeseries_smooth_linewidth_sl = self.timeseries_menu.sliders['smooth_linewidth_sl']
        self.timeseries_smooth_linewidth_sl.setMaximum(self.plot_characteristics['timeseries']['smooth']['format']['linewidth']*100)
        self.timeseries_smooth_linewidth_sl.setValue(self.plot_characteristics['timeseries']['smooth']['format']['linewidth']*10)
        self.timeseries_smooth_window_sl = self.timeseries_menu.sliders['smooth_window_sl']

        # get timeseries interactive dictionary
        self.interactive_elements['timeseries'] = {'hidden': True,
                                                   'markersize_sl': [self.timeseries_markersize_sl],
                                                   'linewidth_sl': [self.timeseries_smooth_linewidth_sl],
                                                   'smooth_window_sl': [self.timeseries_smooth_window_sl]
                                                   }
        
        # PERIODIC PLOT SETTINGS MENU #
        # create periodic settings menu
        self.periodic_menu = SettingsMenu(plot_type='periodic', canvas_instance=self)
        self.periodic_options = self.periodic_menu.checkable_comboboxes['options']
        self.periodic_elements = self.periodic_menu.get_elements()

        # get stats
        self.periodic_stat = self.periodic_menu.comboboxes['stat']
        
        # get sliders and update values
        self.periodic_markersize_sl = self.periodic_menu.sliders['markersize_sl']
        self.periodic_markersize_sl.setMaximum(self.plot_characteristics['periodic']['plot']['markersize']*10)
        self.periodic_markersize_sl.setValue(self.plot_characteristics['periodic']['plot']['markersize'])
        self.periodic_linewidth_sl = self.periodic_menu.sliders['linewidth_sl']
        self.periodic_linewidth_sl.setMaximum(self.plot_characteristics['periodic']['plot']['linewidth']*100)
        self.periodic_linewidth_sl.setValue(self.plot_characteristics['periodic']['plot']['linewidth']*10)

        # get periodic interactive dictionary
        self.interactive_elements['periodic'] = {'hidden': True,
                                                 'markersize_sl': [self.periodic_markersize_sl],
                                                 'linewidth_sl': [self.periodic_linewidth_sl]
                                                 }

        # PERIODIC VIOLIN PLOT SETTINGS MENU #
        # create periodic violin settings menu
        self.periodic_violin_menu = SettingsMenu(plot_type='periodic_violin', canvas_instance=self)
        self.periodic_violin_options = self.periodic_violin_menu.checkable_comboboxes['options']
        self.periodic_violin_elements = self.periodic_violin_menu.get_elements()

        # get sliders and update values
        self.periodic_violin_markersize_sl = self.periodic_violin_menu.sliders['markersize_sl']
        self.periodic_violin_markersize_sl.setMaximum(self.plot_characteristics['periodic-violin']['plot']['median']['markersize']*10)
        self.periodic_violin_markersize_sl.setValue(self.plot_characteristics['periodic-violin']['plot']['median']['markersize'])
        self.periodic_violin_linewidth_sl = self.periodic_violin_menu.sliders['linewidth_sl']
        self.periodic_violin_linewidth_sl.setMaximum(self.plot_characteristics['periodic-violin']['plot']['median']['linewidth']*100)
        self.periodic_violin_linewidth_sl.setValue(self.plot_characteristics['periodic-violin']['plot']['median']['linewidth']*10)

        # get periodic violin interactive dictionary
        self.interactive_elements['periodic_violin'] = {'hidden': True,
                                                        'markersize_sl': [self.periodic_violin_markersize_sl],
                                                        'linewidth_sl': [self.periodic_violin_linewidth_sl]
                                                        }

        # METADATA PLOT SETTINGS MENU #
        # create metadata settings menu
        self.metadata_menu = SettingsMenu(plot_type='metadata', canvas_instance=self)
        self.metadata_options = self.metadata_menu.checkable_comboboxes['options']
        self.metadata_elements = self.metadata_menu.get_elements()

        # get metadata interactive dictionary
        self.interactive_elements['metadata'] = {'hidden': True
                                                 }

        # DISTRIBUTION PLOT SETTINGS MENU #
        # create distribution settings menu
        self.distribution_menu = SettingsMenu(plot_type='distribution', canvas_instance=self)
        self.distribution_options = self.distribution_menu.checkable_comboboxes['options']
        self.distribution_elements = self.distribution_menu.get_elements()

        # get sliders and update values
        self.distribution_linewidth_sl = self.distribution_menu.sliders['linewidth_sl']
        self.distribution_linewidth_sl.setMaximum(self.plot_characteristics['distribution']['plot']['linewidth']*100)
        self.distribution_linewidth_sl.setValue(self.plot_characteristics['distribution']['plot']['linewidth']*10)

        # get distribution interactive dictionary
        self.interactive_elements['distribution'] = {'hidden': True,
                                                     'linewidth_sl': [self.distribution_linewidth_sl]
                                                    }

        # SCATTER PLOT SETTINGS MENU #
        # create scatter settings menu
        self.scatter_menu = SettingsMenu(plot_type='scatter', canvas_instance=self)
        self.scatter_options = self.scatter_menu.checkable_comboboxes['options']
        self.scatter_elements = self.scatter_menu.get_elements()

        # get sliders and update values
        self.scatter_markersize_sl = self.scatter_menu.sliders['markersize_sl']
        self.scatter_markersize_sl.setMaximum(self.plot_characteristics['scatter']['plot']['markersize']*10)
        self.scatter_markersize_sl.setValue(self.plot_characteristics['scatter']['plot']['markersize'])
        self.scatter_regression_linewidth_sl = self.scatter_menu.sliders['regression_linewidth_sl']
        self.scatter_regression_linewidth_sl.setMaximum(self.plot_characteristics['scatter']['regression']['linewidth']*100)
        self.scatter_regression_linewidth_sl.setValue(self.plot_characteristics['scatter']['regression']['linewidth']*10)

        # get scatter interactive dictionary
        self.interactive_elements['scatter'] = {'hidden': True,
                                                'markersize_sl': [self.scatter_markersize_sl],
                                                'linewidth_sl': [self.scatter_regression_linewidth_sl]
                                               }

        # STATSUMMARY PLOT SETTINGS MENU #
        # create statsummary settings menu
        self.statsummary_menu = SettingsMenu(plot_type='statsummary', canvas_instance=self)
        self.statsummary_options = self.statsummary_menu.checkable_comboboxes['options']
        self.statsummary_elements = self.statsummary_menu.get_elements()

        # get stats and add items to cycle
        self.statsummary_stat = self.statsummary_menu.checkable_comboboxes['stat']
        self.statsummary_cycle = self.statsummary_menu.comboboxes['cycle']
        self.statsummary_cycle.addItems(['None', 'Diurnal', 'Weekly', 'Monthly'])
        self.statsummary_periodic_mode = self.statsummary_menu.comboboxes['periodic_mode']
        self.statsummary_periodic_aggregation = self.statsummary_menu.comboboxes['periodic_aggregation']

        # get statsummary interactive dictionary
        self.interactive_elements['statsummary'] = {'hidden': True
                                                    }

        # BOXPLOT PLOT SETTINGS MENU #
        # create boxplot settings menu
        self.boxplot_menu = SettingsMenu(plot_type='boxplot', canvas_instance=self)
        self.boxplot_options = self.boxplot_menu.checkable_comboboxes['options']
        self.boxplot_elements = self.boxplot_menu.get_elements()

        # get boxplot interactive dictionary
        self.interactive_elements['boxplot'] = {'hidden': True
                                               }

        # TAYLOR DIAGRAM SETTINGS MENU #
        # create taylor diagram settings menu
        self.taylor_menu = SettingsMenu(plot_type='taylor', canvas_instance=self)
        self.taylor_options = self.taylor_menu.checkable_comboboxes['options']
        self.taylor_elements = self.taylor_menu.get_elements()

        # get stat
        self.taylor_corr_stat = self.taylor_menu.comboboxes['corr_stat']

        # get sliders and update values
        self.taylor_markersize_sl = self.taylor_menu.sliders['markersize_sl']
        self.taylor_markersize_sl.setMaximum(self.plot_characteristics['taylor']['plot']['markersize']*10)
        self.taylor_markersize_sl.setValue(self.plot_characteristics['taylor']['plot']['markersize'])

        # get statsummary interactive dictionary
        self.interactive_elements['taylor'] = {'hidden': True,
                                               'markersize_sl': [self.taylor_markersize_sl]
                                               }

        # create array with buttons and elements to edit when the canvas is resized or the plots are changed
        self.menu_buttons = []
        self.save_buttons = []
        self.elements = []
        for plot_type in settings_dict.keys():
            if plot_type == 'periodic-violin':
                plot_type = 'periodic_violin'
            self.menu_buttons.append(getattr(self, plot_type + '_menu').buttons['settings_button'])
            self.save_buttons.append(getattr(self, plot_type + '_menu').buttons['save_button'])
            self.elements.append(getattr(self, plot_type + '_elements'))

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
            if event_source == getattr(self, key + '_menu').buttons['settings_button']:
                hidden = self.interactive_elements[key]['hidden']
                elements = getattr(self, key + '_elements')
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
            if 'markersize_sl' in self.interactive_elements[key]:
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
            if 'opacity_sl' in self.interactive_elements[key]:
                if event_source in self.interactive_elements[key]['opacity_sl']:
                    opacity = self.interactive_elements[key]['opacity_sl'][loc].value()/10
                    break
        self.update_opacity(self.plot_axes[key], key, opacity, event_source)

        return None

    def update_linewidth_func(self):
        """ Function to handle the update of the lines widths. """

        event_source = self.sender()
        for key, val in self.interactive_elements.items():
            if 'linewidth_sl' in self.interactive_elements[key]:
                if event_source in self.interactive_elements[key]['linewidth_sl']:
                    linewidth = self.interactive_elements[key]['linewidth_sl'][0].value()/10
                    break
        self.update_linewidth(self.plot_axes[key], key, linewidth)

        return None

    def update_smooth_window_func(self):
        
        # get source
        event_source = self.sender()
        plot_type = event_source.objectName().split('_smooth')[0]
        for element in self.interactive_elements[plot_type]['smooth_window_sl']:
            smooth_window = element.value()
            break

        self.update_smooth_window(self.plot_axes[plot_type], plot_type, smooth_window, 
                                  self.current_plot_options[plot_type])

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
            if event_source.currentData() or self.previous_plot_options[plot_type]:

                self.current_plot_options[plot_type] = copy.deepcopy(event_source.currentData())
                all_plot_options = event_source.currentData(all=True)

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
                    if ((option in self.previous_plot_options[plot_type]) 
                        and (option not in self.current_plot_options[plot_type])):
                        undo = True
                    # make plot option if currently selected
                    elif option in self.current_plot_options[plot_type]:
                        undo = False
                    # do nothing if options were never selected
                    elif ((option not in self.previous_plot_options[plot_type]) 
                        and (option not in self.current_plot_options[plot_type])): 
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
                                    for plot_option in self.current_plot_options[plot_type]:
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

                    # if option is 'bias', then remove all other current option elements (both absolute and bias)
                    if option == 'bias':
                        for active_type in self.plot_elements[plot_type]:
                            if active_type != 'active':
                                for data_label in self.plot_elements[plot_type][active_type]:
                                    for plot_option in self.current_plot_options[plot_type]:
                                        if plot_option in self.plot_elements[plot_type][active_type][data_label]:
                                            for plot_element in self.plot_elements[plot_type][active_type][data_label][plot_option]:
                                                plot_element.remove()
                                            del self.plot_elements[plot_type][active_type][data_label][plot_option] 

                    # options 'logy' and 'logx' 
                    # only plot if axis has all positive values
                    if (option == 'logy') or (option == 'logx'):
                        if isinstance(self.plot_axes[plot_type], dict):
                            for temporal_resolution, sub_ax in self.plot_axes[plot_type].items():
                                log_valid = log_validity(sub_ax, option)
                                if log_valid:
                                    log_axes(sub_ax, option, self.plot_characteristics[plot_type], undo=undo)
                                else:
                                    msg = "It is not possible to log the {0}-axis ".format(option[-1])
                                    msg += "in {0} with negative values.".format(plot_type)
                                    show_message(self.read_instance, msg)
                                    self.read_instance.block_MPL_canvas_updates = True
                                    event_source.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                                    self.read_instance.block_MPL_canvas_updates = False
                                    return None
                        else:
                            log_valid = log_validity(self.plot_axes[plot_type], option)
                            if log_valid:
                                log_axes(self.plot_axes[plot_type], option, self.plot_characteristics[plot_type], undo=undo)
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
                                        annotation(self, 
                                                   self.read_instance, 
                                                   sub_ax,
                                                   self.read_instance.networkspeci,
                                                   self.read_instance.data_labels, 
                                                   plot_type,
                                                   self.plot_characteristics[plot_type], 
                                                   plot_options=self.current_plot_options[plot_type])
                                        break
                            else:
                                annotation(self, 
                                           self.read_instance, 
                                           self.plot_axes[plot_type], 
                                           self.read_instance.networkspeci,
                                           self.read_instance.data_labels, 
                                           plot_type,
                                           self.plot_characteristics[plot_type], 
                                           plot_options=self.current_plot_options[plot_type])

                    # option 'smooth'
                    elif option == 'smooth':
                        if not undo:
                            smooth(self, 
                                   self.read_instance, 
                                   self.plot_axes[plot_type], 
                                   self.read_instance.networkspeci,
                                   self.read_instance.data_labels, 
                                   plot_type,
                                   self.plot_characteristics[plot_type], 
                                   plot_options=self.current_plot_options[plot_type])
                          
                    # option 'regression'
                    elif option == 'regression':
                        if not undo:
                            linear_regression(self, 
                                              self.read_instance, 
                                              self.plot_axes[plot_type], 
                                              self.read_instance.networkspeci,
                                              self.read_instance.data_labels, 
                                              plot_type,
                                              self.plot_characteristics[plot_type],  
                                              plot_options=self.current_plot_options[plot_type])

                    # option 'bias'
                    elif option == 'bias':

                        # firstly if just 1 data label then cannot make bias plot 
                        if len(self.read_instance.data_labels) == 1:
                            msg = 'It is not possible to make a bias plot with just observations loaded.'
                            show_message(self.read_instance, msg)
                            self.read_instance.block_MPL_canvas_updates = True
                            event_source.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                            self.read_instance.block_MPL_canvas_updates = False
                            self.plot_elements[plot_type]['active'] = 'absolute'
                            self.current_plot_options[plot_type].remove('bias')

                            # create other active plot option elements for now absolute plot (if do not already exist)
                            self.redraw_active_options(self.read_instance.data_labels, plot_type, 
                                                       'absolute', self.current_plot_options[plot_type])

                        # if bias option is enabled then first check if bias elements stored
                        elif not undo:

                            # update active (bias)
                            self.plot_elements[plot_type]['active'] = 'bias' 

                            # handle some special case for periodic plot
                            if plot_type == 'periodic':

                                # get currently selected periodic statistic name
                                base_zstat = self.periodic_stat.currentText()
                                zstat = get_z_statistic_comboboxes(base_zstat, bias=True)

                                # get zstat information 
                                zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat) 

                                # if get_z_statistic_type == 'expbias' then return as bias already plotted
                                if z_statistic_type == 'expbias':
                                    self.read_instance.block_MPL_canvas_updates = True
                                    event_source.model().item(index).setCheckState(QtCore.Qt.Unchecked)
                                    self.read_instance.block_MPL_canvas_updates = False
                                    self.plot_elements[plot_type]['active'] = 'absolute'

                            # iterate through valid data labels 
                            bias_labels_to_plot = []
                            for data_label in self.read_instance.data_labels + ['ALL']:

                                # hide absolute plot elements
                                if data_label in self.plot_elements[plot_type]['absolute']:
                                    for element_type in self.plot_elements[plot_type]['absolute'][data_label]:
                                        for element in self.plot_elements[plot_type]['absolute'][data_label][element_type]:
                                            element.set_visible(False)

                                # if have bias elements pre-stored then simply show bias elements (if data label on legend is active)
                                if 'bias' in self.plot_elements[plot_type]:
                                    if data_label in self.plot_elements[plot_type]['bias']:
                                        if (data_label in self.plot_elements['data_labels_active']) or (data_label == 'ALL'):
                                            for element_type in self.plot_elements[plot_type]['bias'][data_label]:
                                                for element in self.plot_elements[plot_type]['bias'][data_label][element_type]:
                                                    element.set_visible(True)
                                    else:
                                        if plot_type == 'statsummary':
                                            if data_label == self.read_instance.observations_data_label:
                                                bias_labels_to_plot.append(data_label) 
                                        else:
                                            if data_label not in [self.read_instance.observations_data_label, 'ALL']:
                                                bias_labels_to_plot.append(data_label) 
                                else:
                                    if plot_type == 'statsummary':
                                        if data_label == self.read_instance.observations_data_label:
                                            bias_labels_to_plot.append(data_label) 
                                    else:
                                        if data_label not in [self.read_instance.observations_data_label, 'ALL']:
                                            bias_labels_to_plot.append(data_label) 

                            # if do not already have bias elements, then make them (tracking plot elements also) 
                            if bias_labels_to_plot:

                                # get plotting function for specific plot
                                if plot_type == 'statsummary':
                                    func = getattr(self.plot, 'make_table')
                                else:
                                    func = getattr(self.plot, 'make_{}'.format(plot_type.split('-')[0]))

                                # call function to update plot
                                # periodic plot
                                if plot_type =='periodic':
                                    func(self.plot_axes[plot_type], self.read_instance.networkspeci, 
                                         bias_labels_to_plot, self.plot_characteristics[plot_type], zstat=zstat, 
                                         plot_options=self.current_plot_options[plot_type])
                                # make statsummary plot
                                elif plot_type == 'statsummary':
                                    relevant_zstats = self.active_statsummary_stats['expbias']
                                    func(self.plot_axes[plot_type], self.read_instance.networkspeci, 
                                         self.read_instance.data_labels, self.plot_characteristics[plot_type], 
                                         zstats=relevant_zstats, statsummary=True, 
                                         plot_options=self.current_plot_options[plot_type])                
                                # other plots
                                else: 
                                    func(self.plot_axes[plot_type], self.read_instance.networkspeci, 
                                         bias_labels_to_plot, self.plot_characteristics[plot_type], 
                                         plot_options=self.current_plot_options[plot_type])

                            # create other active plot option elements for bias plot (if do not already exist)
                            self.redraw_active_options(self.read_instance.data_labels, plot_type, 
                                                       'bias', self.current_plot_options[plot_type])

                        # if bias option is not enabled then hide bias plot elements and show absolute plots again
                        else:
                            
                            # update active (absolute)
                            self.plot_elements[plot_type]['active'] = 'absolute' 

                            # handle some special case for periodic plot
                            if plot_type == 'periodic':

                                # get currently selected periodic statistic name
                                base_zstat = self.periodic_stat.currentText()
                                zstat = get_z_statistic_comboboxes(base_zstat, bias=False)

                                # get zstat information 
                                zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat) 

                            # iterate through valid data labels
                            absolute_labels_to_plot = []
                            for data_label in list(self.read_instance.data_labels) + ['ALL']:

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
                                else:
                                    if plot_type == 'statsummary':
                                        if data_label == self.read_instance.observations_data_label:
                                            absolute_labels_to_plot.append(data_label) 
                                    else:
                                        if data_label != 'ALL':
                                            absolute_labels_to_plot.append(data_label) 

                            # if do not already have absolute elements, then make them (tracking plot elements also) 
                            if absolute_labels_to_plot:

                                # get plotting function for specific plot
                                if plot_type == 'statsummary':
                                    func = getattr(self.plot, 'make_table')
                                else:
                                    func = getattr(self.plot, 'make_{}'.format(plot_type.split('-')[0]))

                                # call function to update plot
                                # periodic plot
                                if plot_type =='periodic':
                                    func(self.plot_axes[plot_type], self.read_instance.networkspeci, 
                                         absolute_labels_to_plot, self.plot_characteristics[plot_type], zstat=zstat, 
                                         plot_options=self.current_plot_options[plot_type])
                                # make statsummary plot
                                elif plot_type == 'statsummary':
                                    relevant_zstats = self.active_statsummary_stats['basic']
                                    func(self.plot_axes[plot_type], self.read_instance.networkspeci, 
                                         self.read_instance.data_labels, self.plot_characteristics[plot_type], 
                                         zstats=relevant_zstats, statsummary=True, 
                                         plot_options=self.current_plot_options[plot_type])                
                                # other plots
                                else: 
                                    func(self.plot_axes[plot_type], self.read_instance.networkspeci, 
                                         absolute_labels_to_plot, self.plot_characteristics[plot_type], 
                                         plot_options=self.current_plot_options[plot_type])

                            # create other active plot option elements for absolute plot (if do not already exist)
                            self.redraw_active_options(self.read_instance.data_labels, 
                                                       plot_type, 'absolute', self.current_plot_options[plot_type])

                    # check stats for the selected periodic cycle
                    if plot_type in ['statsummary']:
                        self.read_instance.block_config_bar_handling_updates = True
                        self.check_statsummary_stats()
                        self.read_instance.block_config_bar_handling_updates = False

                    # reset axes limits (harmonising across subplots for periodic plots) 
                    if plot_type != 'map':
                        if plot_type == 'scatter':
                            harmonise_xy_lims_paradigm(self, self.read_instance, self.plot_axes[plot_type], plot_type, 
                                                       self.plot_characteristics[plot_type], 
                                                       self.current_plot_options[plot_type], relim=True)
                        elif plot_type != 'taylor':
                            harmonise_xy_lims_paradigm(self, self.read_instance, self.plot_axes[plot_type], plot_type, 
                                                       self.plot_characteristics[plot_type], 
                                                       self.current_plot_options[plot_type], relim=True, autoscale=True)                       

                # save current plot options as previous
                self.previous_plot_options[plot_type] = self.current_plot_options[plot_type]

                # draw changes
                self.figure.canvas.draw_idle()

        return None

    def redraw_active_options(self, data_labels, plot_type, active, plot_options):
        """ Redraw active plot option elements when moving between absolute and bias plots,
            if do not already exist.
        """

        # if 'bias' is active, remove observations data label from data labels
        data_labels_alt = copy.deepcopy(data_labels)
        if active == 'bias':
            data_labels_alt.remove(self.read_instance.observations_data_label)

        # iterate through plot_options
        for plot_option in plot_options:

            if plot_option == 'annotate':
                if isinstance(self.plot_axes[plot_type], dict):
                    for relevant_temporal_resolution, sub_ax in self.plot_axes[plot_type].items():
                        if relevant_temporal_resolution in self.read_instance.relevant_temporal_resolutions:
                            annotation(self, self.read_instance, sub_ax, self.read_instance.networkspeci, data_labels, 
                                       plot_type, self.plot_characteristics[plot_type], plot_options=plot_options)
                            break
                else:
                    annotation(self, self.read_instance, self.plot_axes[plot_type], self.read_instance.networkspeci,
                               data_labels, plot_type, self.plot_characteristics[plot_type], plot_options=plot_options)

            elif plot_option == 'smooth':
                smooth(self, self.read_instance, self.plot_axes[plot_type], self.read_instance.networkspeci,
                       data_labels_alt, plot_type, self.plot_characteristics[plot_type], plot_options=plot_options)
            
            elif plot_option == 'regression':
                linear_regression(self, self.read_instance, self.plot_axes[plot_type], self.read_instance.networkspeci,
                                  data_labels_alt, plot_type, self.plot_characteristics[plot_type],  
                                  plot_options=plot_options)

    def update_markersize(self, ax, plot_type, markersize, event_source):
        """ Update markers size for each plot type. """
        
        # set markersize
        if plot_type in ['timeseries', 'periodic', 'scatter', 'periodic-violin', 'taylor']:
            
            if isinstance(ax, dict):
                for sub_ax in ax.values():
                    for line in sub_ax.lines:
                        line.set_markersize(markersize)
            else:
                if plot_type == 'taylor':
                    for line in self.plot.taylor_polar_relevant_axis.lines:
                        line.set_markersize(markersize)
                else:
                    for line in ax.lines:
                        line.set_markersize(markersize)

            # update characteristics per plot type
            # this is made to keep the changes when selecting stations with lasso
            if plot_type in ['timeseries', 'periodic', 'scatter', 'taylor']:
                self.plot_characteristics[plot_type]['plot']['markersize'] = markersize
            elif plot_type == 'periodic-violin':
                self.plot_characteristics[plot_type]['plot']['median']['markersize'] = markersize

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

        # draw changes
        self.figure.canvas.draw_idle()

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
        self.figure.canvas.draw_idle()

        return None

    def update_linewidth(self, ax, plot_type, linewidth):
        """ Update line widths for each plot type. """
        
        # set linewidth
        if isinstance(ax, dict):
            for sub_ax in ax.values():
                for line in sub_ax.lines:
                    if line not in self.annotation_elements:
                        line.set_linewidth(linewidth)
        else:
            for line in ax.lines:
                if line not in self.annotation_elements:
                    if (((plot_type == 'scatter') and ((list(line.get_xdata()) == [0, 0.5])
                        or (list(line.get_xdata()) == [0, 1])))):
                        continue
                    else:
                        line.set_linewidth(linewidth)
                else:
                    line.set_linewidth(linewidth)

        # update characteristics per plot type
        if plot_type == 'periodic-violin':
            self.plot_characteristics[plot_type]['plot']['median']['linewidth'] = linewidth
        elif plot_type == 'timeseries':
            self.plot_characteristics[plot_type]['smooth']['format']['linewidth'] = linewidth
        elif plot_type == 'regression':
            self.plot_characteristics[plot_type]['regression']['linewidth'] = linewidth
        else:
            self.plot_characteristics[plot_type]['plot']['linewidth'] = linewidth

        # draw changes
        self.figure.canvas.draw_idle()

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

        # draw changes
        self.figure.canvas.draw_idle()

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

        # draw changes
        self.figure.canvas.draw_idle()

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
        self.figure.canvas.draw_idle()

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
        self.figure.canvas.draw_idle()

        return None
