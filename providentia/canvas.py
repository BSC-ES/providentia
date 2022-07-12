from .read_aux import drop_nans
from .filter import DataFilter
from .statistics import to_pandas_dataframe
from .statistics import calculate_z_statistic
from .statistics import generate_colourbar
from .statistics import get_z_statistic_info
from .statistics import get_z_statistic_comboboxes
from .plot import Plot
from providentia import aux

import copy
import json
import os
from weakref import WeakKeyDictionary

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg \
        as FigureCanvas

from matplotlib.figure import Figure
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector
import matplotlib.gridspec as gridspec
from pandas.plotting import register_matplotlib_converters
from PyQt5 import QtCore
from PyQt5 import QtWidgets

import numpy as np
import pandas as pd

# Make sure that we are using Qt5 backend with matplotlib
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

        #create figure and canvas objects
        self.figure = Figure(dpi=100)
        FigureCanvas.__init__(self, self.figure)

        #add passed arguments to self
        self.read_instance = read_instance

        # load characteristics per plot type
        #check for self defined plot characteristics file
        if hasattr(self, 'plot_characteristics_filename'):
            self.plot_characteristics_templates = json.load(self.plot_characteristics_filename)
        else:
            self.plot_characteristics_templates = json.load(open(os.path.join(
                CURRENT_PATH, 'conf/plot_characteristics_dashboard.json')))
        self.plot_characteristics = {}

        #add general plot characteristics to self
        for k, val in self.plot_characteristics_templates['general'].items():
            setattr(self, k, val)

        #initialise some key vars
        self.temporal_colocation = False
        self.filter_data = None

        # initialise Plot class
        self.plot = Plot(read_instance=self.read_instance, canvas_instance=self)

        # setup gridding of canvas
        self.gridspec = gridspec.GridSpec(self.gridspec_nrows, self.gridspec_ncols)
        self.gridspec.update(**self.gridspec_format)

        #create plot axes on grid
        self.plot_axes = {}
        self.plot_axes['map'] = self.figure.add_subplot(self.gridspec.new_subplotspec((0, 0), rowspan=44, colspan=45), projection=self.plotcrs)
        self.plot_axes['legend'] = self.figure.add_subplot(self.gridspec.new_subplotspec((0, 47), rowspan=8,  colspan=53))
        self.plot_axes['timeseries'] = self.figure.add_subplot(self.gridspec.new_subplotspec((12, 50), rowspan=35, colspan=50))
        self.plot_axes['periodic-violin'] = {}
        self.plot_axes['periodic-violin']['hour'] = self.figure.add_subplot(self.gridspec.new_subplotspec((56, 70), rowspan=20, colspan=30))
        self.plot_axes['periodic-violin']['dayofweek'] = self.figure.add_subplot(self.gridspec.new_subplotspec((80, 91), rowspan=20, colspan=9))
        self.plot_axes['periodic-violin']['month'] = self.figure.add_subplot(self.gridspec.new_subplotspec((80, 70), rowspan=20, colspan=18))
        self.plot_axes['periodic'] = {}
        self.plot_axes['periodic']['hour'] = self.figure.add_subplot(self.gridspec.new_subplotspec((56, 35), rowspan=20, colspan=30))
        self.plot_axes['periodic']['dayofweek'] = self.figure.add_subplot(self.gridspec.new_subplotspec((80, 56), rowspan=20, colspan=9))
        self.plot_axes['periodic']['month'] = self.figure.add_subplot(self.gridspec.new_subplotspec((80, 35), rowspan=20, colspan=18))
        self.plot_axes['metadata'] = self.figure.add_subplot(self.gridspec.new_subplotspec((56, 0),  rowspan=44, colspan=30))
        self.plot_axes['cb'] = self.figure.add_axes([0.0553, 0.53, 0.35, 0.02])

        # setup interactive picker/lasso on map
        self.figure.canvas.mpl_connect('pick_event', self.on_click)
        self.lasso = LassoSelector(self.plot_axes['map'], onselect=self.onlassoselect,
                                    useblit=True, lineprops=self.lasso)
        # initialise variable that informs whether to use picker/lasso for updates
        self.map_already_updated = False
        # initialise variable of valid station indices plotted on map as empty list
        self.active_map_valid_station_inds = np.array([], dtype=np.int)

        #define plots to update upon selecting stations
        self.selected_station_plots = ['timeseries','periodic-violin','periodic','metadata']

        #set plot characteristics for base plot types
        self.plot.set_plot_characteristics(['legend','map'])
        #update plot characteristics for selected_station_plots
        self.plot.set_plot_characteristics(self.selected_station_plots)

        #Format and then hide all axes
        for plot_type, ax in self.plot_axes.items():
            if type(ax) == dict:
                for relevant_temporal_resolution, sub_ax in ax.items():
                    if relevant_temporal_resolution in ['hour','month']:
                        col_ii = 0
                    else:
                        col_ii = 1
                    self.plot.format_axis(sub_ax, plot_type, self.plot_characteristics[plot_type], relevant_temporal_resolution=relevant_temporal_resolution, col_ii=col_ii)
                    self.remove_axis_elements(sub_ax, plot_type)
            else:
                if plot_type != 'cb':
                    self.plot.format_axis(ax, plot_type, self.plot_characteristics[plot_type])
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

        # set map extents
        self.plot_axes['map'].set_extent(self.read_instance.map_extent, crs=self.datacrs)
        #set navigation toolbar stack
        self.reset_ax_navigation_toolbar_stack(self.plot_axes['map'])

        # update plotted map z statistic
        self.update_map_z_statistic()

        # plot experiment grid domain edges on map
        self.update_experiment_domain_edges()

        # update legend
        self.update_legend()

        # draw changes
        self.figure.canvas.draw()

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

    def handle_data_filter_update(self):
        """Function which handles updates data filtering by
        selected lower/upper limit bounds, selected metadata, period codes, 
        and selected minimum data availability %
        """

        print('UPDATE DATA FILTER')

        # return if nothing has been loaded yet
        if not hasattr(self.read_instance.datareader, 'data_in_memory'):
            return

        check_state = self.read_instance.ch_colocate.checkState()
        # update variable to inform plotting functions whether to use colocated data/or not
        if check_state == QtCore.Qt.Checked:
            self.temporal_colocation = True
        else:
            self.temporal_colocation = False

        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        if self.filter_data is None:
            self.filter_data = DataFilter(self.read_instance)
        else:
            self.filter_data.filter_all()
            self.update_active_map()
        QtWidgets.QApplication.restoreOverrideCursor()

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

    def handle_colocate_update(self):
        """Function that handles the update of the MPL canvas
        with colocated data upon checking of the colocate checkbox"""

        if not self.read_instance.block_MPL_canvas_updates:

            print('COLOCATE UPDATE')

            # if only have 1 data array in memory (i.e. observations), no colocation is possible,
            # therefore set colocation_active to be False, and return
            if len(list(self.read_instance.data_in_memory_filtered.keys())) == 1:
                self.temporal_colocation = False
                return

            # else, if have loaded experiment data, check if colocate checkbox is checked or unchecked
            check_state = self.read_instance.ch_colocate.checkState()
            # update variable to inform plotting functions whether to use colocated data/or not
            if check_state == QtCore.Qt.Checked:
                self.temporal_colocation = True
            else:
                self.temporal_colocation = False

            # update z statistic/experiment bias comboboxes (without updating canvas)
            self.read_instance.block_MPL_canvas_updates = True
            self.handle_map_z_statistic_update()
            self.handle_experiment_bias_update()
            self.read_instance.block_MPL_canvas_updates = False

            # update plotted map z statistic
            self.update_map_z_statistic()

            # update associated plots with selected stations
            self.update_associated_selected_station_plots()

            # draw changes
            self.figure.canvas.draw()

    def update_map_z_statistic(self):
        """Function that updates plotted z statistic on map, with colourbar"""

        print('UPDATE MAP Z')

        # remove axis elements from map/cb
        self.remove_axis_elements(self.plot_axes['map'], 'map')
        self.remove_axis_elements(self.plot_axes['cb'], 'cb')

        # get necessary data arrays
        z1 = aux.get_data_label(self.temporal_colocation, self.read_instance.cb_z1.currentText())
        z2 = aux.get_data_label(self.temporal_colocation, self.read_instance.cb_z2.currentText())

        #get zstat name from combobox 
        base_zstat = self.read_instance.cb_z_stat.currentText()
        zstat = get_z_statistic_comboboxes(base_zstat, second_data_label=z2)

        # calculate map z statistic (for selected z statistic) --> updating active map valid station indices
        z_statistic, active_map_valid_station_inds = calculate_z_statistic(self.read_instance, z1, z2, zstat, self.temporal_colocation)
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
            # unselect all/intersect checkboxes
            self.read_instance.block_MPL_canvas_updates = True
            self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
            self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
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

            # plot new station points on map - coloured by currently active z statisitic, setting up plot picker
            self.plot.make_map(self.plot_axes['map'], z_statistic, self.plot_characteristics['map'])

            # create 2D numpy array of plotted station coordinates
            self.map_points_coordinates = np.vstack((self.read_instance.datareader.station_longitudes[self.active_map_valid_station_inds], self.read_instance.datareader.station_latitudes[self.active_map_valid_station_inds])).T

            #generate colourbar
            generate_colourbar(self.read_instance, [self.plot_axes['map']], [self.plot_axes['cb']], zstat, self.plot_characteristics['map'])

            #activate map/cb axes
            self.activate_axis(self.plot_axes['map'], 'map')
            self.activate_axis(self.plot_axes['cb'], 'cb')

        #re-draw (needed to update plotted colours before update_map_station_selection)
        self.figure.canvas.draw()

        # update map selection appropriately for z statistic
        self.update_map_station_selection()

    def update_map_station_selection(self):
        """Function that updates the visual selection of stations on map"""

        print('UPDATE STATION SELECTION')

        # update map title
        if len(self.relative_selected_station_inds) == 1:
            axis_title_label = '{} Selected'.format(
                self.read_instance.station_references[self.relative_selected_station_inds][0])
        else:
            axis_title_label = '{} Stations Selected of {} Available'.format(
                len(self.relative_selected_station_inds), len(self.active_map_valid_station_inds))
        self.plot.set_axis_title(self.plot_axes['map'], axis_title_label, self.plot_characteristics['map'])

        # reset alphas and marker sizes of all plotted stations (if have some stations on map)
        if len(self.active_map_valid_station_inds) > 0:
            # reset alphas
            rgba_tuples = self.plot_axes['map'].collections[-1].get_facecolor()
            rgba_tuples[:, -1] = self.plot_characteristics['map']['marker_selected']['alpha']
            self.plot_axes['map'].collections[-1].set_facecolor(rgba_tuples)
            # reset marker sizes
            marker_sizes = np.full(len(self.active_map_valid_station_inds), self.plot_characteristics['map']['marker_unselected']['s'])
            self.plot_axes['map'].collections[-1].set_sizes(marker_sizes)

            # if have some selected stations, update map plot with station selection
            # (reducing alpha of non-selected stations, and increasing marker size of selected stations)
            if len(self.relative_selected_station_inds) > 0:

                absolute_non_selected_stations = np.nonzero(~np.in1d(range(len(self.active_map_valid_station_inds)),
                                                                     self.absolute_selected_station_inds))[0]
                
                # decrease alpha of non-selected stations
                if len(absolute_non_selected_stations) > 0:
                    rgba_tuples[absolute_non_selected_stations, -1] = self.plot_characteristics['map']['marker_unselected']['alpha']
                    self.plot_axes['map'].collections[-1].set_facecolor(rgba_tuples)

                # increase marker size of selected stations
                marker_sizes[self.absolute_selected_station_inds] = self.plot_characteristics['map']['marker_selected']['s']
                self.plot_axes['map'].collections[-1].set_sizes(marker_sizes)
        
        #redraw points
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def update_associated_selected_station_plots(self):
        """Function that updates all plots associated with selected stations on map"""

        print('UPDATE SELECTED STATION PLOTS')

        # clear all previously plotted artists from selected station plots and hide axes 
        for plot_type in self.selected_station_plots:
            ax  = self.plot_axes[plot_type]
            if type(ax) == dict:
                for sub_ax in ax.values():
                    self.remove_axis_elements(sub_ax, plot_type)
            else:
                self.remove_axis_elements(ax, plot_type)

        # if have selected stations on map, then now remake plots
        if len(self.relative_selected_station_inds) > 0:

            # put selected data for each data array into pandas dataframes
            to_pandas_dataframe(read_instance=self.read_instance, canvas_instance=self)

            #iterate through selected_station_plots
            for plot_type in self.selected_station_plots:

                #get relevant axis
                ax = self.plot_axes[plot_type]

                #get options defined to configure plot (e.g. bias, individual, annotate, etc.)
                plot_options = plot_type.split('_')[1:]

                #get plotting function for specific plot
                func = getattr(self.plot, 'make_{}'.format(plot_type.split('-')[0]))

                #update ylabel for plot
                if plot_type in ['periodic']:
                    # get currently selected experiment bias statistic name
                    base_zstat = self.read_instance.cb_experiment_bias_stat.currentText()
                    #if have no zstat, it is because no experiments are loaded so cannot make plot type
                    if not base_zstat:
                        continue
                    zstat = get_z_statistic_comboboxes(base_zstat, second_data_label='model')
                    # get zstat information 
                    zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat)
                    #add bias to plot options if stat is bias
                    if z_statistic_sign == 'bias':
                        plot_options.append('bias')      
                    if z_statistic_type == 'basic':
                        if base_zstat not in ['Data %', 'Exceedances']:
                            ax_label = self.read_instance.datareader.measurement_units 
                        else:
                            ax_label = basic_stats[base_zstat]['label']
                    else:
                        ax_label = expbias_stats[base_zstat]['label']
                elif plot_type == 'metadata':
                    ax_label = ''
                else:
                    ax_label = self.read_instance.datareader.measurement_units
                self.plot.set_axis_label(ax, 'y', ax_label, self.plot_characteristics[plot_type])

                #add title for plot
                if plot_type in ['periodic']:
                    ax_title = '{}-{}'.format(plot_type, zstat)
                else:
                    ax_title = plot_type
                self.plot.set_axis_title(ax, ax_title, self.plot_characteristics[plot_type])

                # iterate through data array names in selected station data dictionary
                for data_label in list(self.selected_station_data.keys()):

                    #get original data label
                    original_data_label = data_label.split('_')[0]

                    #call function to update plot
                    if plot_type in ['periodic']:
                        # skip observational array if bias stat and data array is observations
                        if (z_statistic_sign == 'bias') & (original_data_label == 'observations'):
                            continue
                        func(ax, data_label, self.plot_characteristics[plot_type], zstat=zstat, plot_options=plot_options)
                    elif plot_type == 'metadata':
                        if original_data_label != 'observations':
                            continue
                        func(ax, self.plot_characteristics[plot_type], plot_options=plot_options)
                    else:
                        func(ax, data_label, self.plot_characteristics[plot_type], plot_options=plot_options)

                # un-hide axes, reset axes limits (harmonise across subplots for perioduc plots), and reset navigation toolbar stack
                if type(ax) == dict:
                    if plot_type == 'periodic-violin':
                        self.plot.harmonise_xy_lims_paradigm(list(ax.values()), plot_type, self.plot_characteristics[plot_type], plot_options, ylim=[self.selected_station_data_min, self.selected_station_data_max])
                    else:
                        self.plot.harmonise_xy_lims_paradigm(list(ax.values()), plot_type, self.plot_characteristics[plot_type], plot_options, relim=True)
                    for sub_ax in ax.values():
                        self.activate_axis(sub_ax, plot_type)
                        self.reset_ax_navigation_toolbar_stack(sub_ax)
                else:    
                    ax.relim()
                    ax.autoscale()
                    self.activate_axis(ax, plot_type)
                    self.reset_ax_navigation_toolbar_stack(ax)


    def update_experiment_domain_edges(self):
        """Function that plots grid domain edges of experiments in memory"""

        #remove grid domain polygon if previously plotted
        for pol in self.plot_axes['map'].patches:  
            if type(pol) == matplotlib.patches.Polygon:
                pol.remove()

        #create grid edge polygons for experiments in memory
        grid_edge_polygons = self.plot.make_experiment_domain_polygons()

        # plot grid edge polygons on map
        for grid_edge_polygon in grid_edge_polygons:
            self.plot_axes['map'].add_patch(grid_edge_polygon)

    def update_legend(self):
        """Function that updates legend"""

        # create legend element handles
        legend_plot_characteristics = self.plot.make_legend_handles(copy.deepcopy(self.plot_characteristics['legend']))

        # plot legend
        self.plot_axes['legend'].legend(**legend_plot_characteristics['plot'])

        #un-hide legend
        self.activate_axis(self.plot_axes['legend'], 'legend')

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
            if (not self.temporal_colocation) or (selected_z2_array == ''):
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

    def handle_experiment_bias_update(self):
        """Define function that handles update of plotted experiment bias statistics"""

        if not self.read_instance.block_config_bar_handling_updates:

            print('UPDATE EXP BIAS')

            # if no experiment data loaded, do not update
            if len(self.read_instance.experiment_bias_types) > 0:

                # update experiment bias comboboxes

                # set variable that blocks configuration bar handling updates until all changes
                # to the experiment bias comboboxes are made
                self.read_instance.block_config_bar_handling_updates = True

                # get currently selected items
                selected_experiment_bias_type = self.read_instance.cb_experiment_bias_type.currentText()
                base_zstat = self.read_instance.cb_experiment_bias_stat.currentText()

                # update experiment bias statistics (used for Aggregated field), to all basic stats
                # if colocation not-active, and basic+bias stats if colocation active
                if not self.temporal_colocation:
                    available_experiment_bias_stats = copy.deepcopy(self.read_instance.basic_z_stats)
                else:
                    available_experiment_bias_stats = copy.deepcopy(self.read_instance.basic_and_bias_z_stats)

                # if selected bias type is empty string, it is because fields are being initialised for the first time
                if selected_experiment_bias_type == '':
                    # set experiment bias type to be first available type
                    selected_experiment_bias_type = self.read_instance.experiment_bias_types[0]
                    # set experiment bias stat to be first available stat
                    base_zstat = available_experiment_bias_stats[0]
                    if hasattr(self.read_instance, 'exp_bias_stat'):
                        if self.read_instance.exp_bias_stat in available_experiment_bias_stats:
                            base_zstat = self.read_instance.exp_bias_stat

                # if selected bias type is 'Rank', then there are no stat options so force the
                # available items and selected stat to be empty
                if selected_experiment_bias_type == 'Rank':
                    available_experiment_bias_stats = []
                    base_zstat = ''

                # update all comboboxes (clear, then add items)
                self.read_instance.cb_experiment_bias_type.clear()
                self.read_instance.cb_experiment_bias_stat.clear()
                self.read_instance.cb_experiment_bias_type.addItems(self.read_instance.experiment_bias_types)
                self.read_instance.cb_experiment_bias_stat.addItems(available_experiment_bias_stats)

                # update selected values
                self.read_instance.cb_experiment_bias_type.setCurrentText(selected_experiment_bias_type)
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

                        # get currently selected experiment bias statistic name
                        zstat = get_z_statistic_comboboxes(base_zstat, second_data_label='model')
                        # get zstat information 
                        zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat)
                        #add bias to plot options if stat is bias
                        plot_options = []
                        if z_statistic_sign == 'bias':
                            plot_options.append('bias')  
                        #update ylabel text
                        if z_statistic_type == 'basic':
                            if base_zstat not in ['Data %', 'Exceedances']:
                                ax_label = self.read_instance.datareader.measurement_units 
                            else:
                                ax_label = basic_stats[base_zstat]['label']
                        else:
                            ax_label = expbias_stats[base_zstat]['label']
                        self.plot.set_axis_label(self.plot_axes['periodic'], 'y', ax_label, self.plot_characteristics['periodic'])

                        #update axis title
                        self.plot.set_axis_title(self.plot_axes['periodic'], 'periodic-{}'.format(zstat), self.plot_characteristics['periodic'])

                        # if experiment bias type == 'Aggregated' --> update plotted experiment bias plots
                        if selected_experiment_bias_type == 'Aggregated':

                            for data_label in list(self.selected_station_data.keys()):
                                # skip observational array if bias stat
                                if (z_statistic_sign == 'bias') & (data_label.split('_')[0] == 'observations'):
                                    continue
                                self.plot.make_periodic(self.plot_axes['periodic'], data_label, self.plot_characteristics['periodic'], zstat=zstat)

                        # un-hide axes, harmonise axes limits across subplots, and reset navigation toolbar stack
                        self.plot.harmonise_xy_lims_paradigm(list(self.plot_axes['periodic'].values()), 'periodic', self.plot_characteristics['periodic'], plot_options, relim=True)
                        for sub_ax in self.plot_axes['periodic'].values():
                            self.activate_axis(sub_ax, 'periodic')
                            self.reset_ax_navigation_toolbar_stack(sub_ax)

                        # draw changes
                        self.figure.canvas.draw()

    def select_all_stations(self):
        """Define function that selects/unselects all plotted stations
        (and associated plots) upon ticking of checkbox"""

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

            # draw changes
            self.figure.canvas.draw()

    def select_intersect_stations(self):
        """Define function that selects/unselects intersection of
        stations and all experiment domains (and associated plots)
        upon ticking of checkbox
        """

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
                if len(list(self.read_instance.datareader.data_in_memory.keys())) == 1:
                    self.relative_selected_station_inds = copy.deepcopy(self.active_map_valid_station_inds)
                    self.absolute_selected_station_inds = np.arange(len(self.relative_selected_station_inds),
                                                                    dtype=np.int)
                # else, define list of lists to get intersection between (active_map_
                # valid_station_inds, and valid station indices associated with each loaded experiment array)
                else:
                    intersect_lists = [self.active_map_valid_station_inds]
                    for data_label in list(self.read_instance.datareader.data_in_memory.keys()):
                        if data_label != 'observations':
                            intersect_lists.append(
                                self.read_instance.datareader.plotting_params[data_label]['valid_station_inds'])
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

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
                self.update_associated_selected_station_plots()

            # draw changes
            self.figure.canvas.draw()

    # define functions that handle interactive station selection on map
    # the selection methods are individual station selection, or multiple selection with lasso

    def on_click(self, event):
        """Function that handles single station selection upon mouse click."""

        # update variable to inform lasso handler that map as already been updated (to not redraw map)
        # the on_click function is only called when a station index has been selected
        # the variable will be reset by lasso handler (which is always called after on_click)
        self.map_already_updated = True

        # if the mouse click is not recognised as left(1)/right(3) then return from function
        if event.mouseevent.button not in [1, 3]:
            return

        # make copy of current full array relative selected stations indices, before selecting new ones
        self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)
        self.previous_absolute_selected_station_inds = copy.deepcopy(self.absolute_selected_station_inds)

        # get absolute selected index of station on map
        absolute_selected_station_inds = np.array([event.ind[0]], dtype=np.int)

        # get new relative selected index with respect to all available stations
        relative_selected_station_inds = self.map_selected_station_inds_to_all_available_inds(
            absolute_selected_station_inds)

        # if left click (code of 1) --> select station
        if event.mouseevent.button == 1:
            self.absolute_selected_station_inds = absolute_selected_station_inds
            self.relative_selected_station_inds = relative_selected_station_inds
        # if right click (code of 3) --> unselect station (if currently selected),
        # select station (if currently unselected)
        elif event.mouseevent.button == 3:
            relative_index = np.where(self.previous_relative_selected_station_inds == relative_selected_station_inds)[0]
            if len(relative_index) > 0:
                self.relative_selected_station_inds = np.delete(self.relative_selected_station_inds, relative_index)
            else:
                self.relative_selected_station_inds = np.append(self.relative_selected_station_inds,
                                                                relative_selected_station_inds)

            absolute_index = np.where(self.previous_absolute_selected_station_inds == absolute_selected_station_inds)[0]
            if len(absolute_index) > 0:
                self.absolute_selected_station_inds = np.delete(self.absolute_selected_station_inds, absolute_index)
            else:
                self.absolute_selected_station_inds = np.append(self.absolute_selected_station_inds,
                                                                absolute_selected_station_inds)

        # update map station selection
        self.update_map_station_selection()

        # if selected stations have changed from previous selected, update associated plots
        if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
            self.update_associated_selected_station_plots()

        # draw changes
        self.figure.canvas.draw()

    def onlassoselect(self, verts):
        """Function that handles multiple
        station selection upon lasso drawing
        """
        # unselect all/intersect checkboxes
        self.read_instance.block_MPL_canvas_updates = True
        self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
        self.read_instance.block_MPL_canvas_updates = False

        # check if have any plotted stations on map, if not, return
        if len(self.active_map_valid_station_inds) == 0:
            return

        # check if map as already been processed by on_click mouse click handling function, if so, return
        # the on_click function will always be called before lasso handler
        if self.map_already_updated:
            self.map_already_updated = False
            return

        # make copy of current full array relative selected stations indices, before selecting new ones
        self.previous_relative_selected_station_inds = copy.deepcopy(self.relative_selected_station_inds)

        # get coordinates of drawn lasso
        lasso_path = Path(verts)
        lasso_path_vertices = lasso_path.vertices
        # transform lasso coordinates from projected to standard geographic coordinates
        lasso_path.vertices = \
            self.datacrs.transform_points(self.plotcrs, lasso_path_vertices[:, 0], lasso_path_vertices[:, 1])[:, :2]
        # get absolute selected indices of stations on map (the station coordinates contained within lasso)
        self.absolute_selected_station_inds = np.nonzero(lasso_path.contains_points(self.map_points_coordinates))[0]

        # get selected station indices with respect to all available stations
        self.relative_selected_station_inds = \
            self.map_selected_station_inds_to_all_available_inds(self.absolute_selected_station_inds)

        # update map station selection
        self.update_map_station_selection()

        # hide lasso after selection
        self.lasso.set_visible(False)

        # if selected stations have changed from previous selected, update associated plots
        if not np.array_equal(self.previous_relative_selected_station_inds, self.relative_selected_station_inds):
            self.update_associated_selected_station_plots()

        # draw changes
        self.figure.canvas.draw()

    def map_selected_station_inds_to_all_available_inds(self, selected_map_inds):
        """Takes the indices of selected stations on the map
        (potentially a subset of all available stations), and returns the indices
        of the stations inside the full loaded data arrays
        """

        # index the array of indices of stations plotted on the map (indexed with respect to
        # all available stations), with the absolute indices of the subset of plotted selected stations
        return self.active_map_valid_station_inds[selected_map_inds]

    def remove_axis_elements(self, ax, plot_type):
        """Removes all plotted axis elements,
           and hide axis.
        """

        #remove all plotted axis elements
        if plot_type == 'legend':
            leg = ax.get_legend()
            if leg:
                leg.remove()

        elif plot_type == 'map':
            for col in ax.collections:            
                if type(col) == matplotlib.collections.PathCollection:
                    col.remove()

        elif plot_type == 'cb':
            ax.artists = []
            ax.collections = [] 

        elif plot_type == 'timeseries':
            ax.lines = []
        
        elif plot_type == 'periodic':
            ax.lines = []

        elif plot_type == 'periodic-violin':
            ax.lines = []
            ax.collections = []

        elif plot_type == 'metadata':
            ax.texts = []
        
        #hide axis
        ax.axis('off')
        ax.set_visible(False)

    def activate_axis(self, ax, plot_type):
        """Un-hide axis"""

        #un-hide axis
        ax.axis('on')
        ax.set_visible(True)







        
