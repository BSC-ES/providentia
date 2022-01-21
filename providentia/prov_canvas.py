""" MPL Canvas module """
from .calculate import Stats
from .calculate import ExpBias
from .reading import drop_nans
from .filter import DataFilter
from providentia import aux

import copy
from weakref import WeakKeyDictionary

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg \
        as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector
from matplotlib.gridspec import GridSpec
from pandas.plotting import register_matplotlib_converters
from PyQt5 import QtCore
from PyQt5 import QtWidgets

import numpy as np
import pandas as pd

# Make sure that we are using Qt5 backend with matplotlib
matplotlib.use('Qt5Agg')
register_matplotlib_converters()


class MPLCanvas(FigureCanvas):
    """Class that handles the creation and updates of
    a matplotlib canvas, and associated subplots
    """

    def __init__(self, read_instance, parent=None):
        """Initialise the MPL canvas"""

        self.figure = Figure(dpi=100)
        self.canvas = FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)
        self.read_instance = read_instance
        self.filter_data = None

        # set cartopy data directory
        cartopy.config['pre_existing_data_dir'] = self.read_instance.cartopy_data_dir

        # setup gridding of canvas
        self.gridspec = GridSpec(100, 100)
        self.gridspec.update(left=0.01, right=0.99, top=0.96, bottom=0.04, wspace=0.00, hspace=0.00)

        # map_ax =              gridspec.new_subplotspec((0, 0),   rowspan=45, colspan=45)
        legend_ax = self.gridspec.new_subplotspec((0, 46),  rowspan=8,  colspan=54)
        ts_ax = self.gridspec.new_subplotspec((12, 54), rowspan=34, colspan=46)
        violin_hours_ax = self.gridspec.new_subplotspec((57, 70), rowspan=19, colspan=30)
        violin_months_ax = self.gridspec.new_subplotspec((81, 70), rowspan=19, colspan=18)
        violin_days_ax = self.gridspec.new_subplotspec((81, 91), rowspan=19, colspan=9)
        exp_bias_hours_ax = self.gridspec.new_subplotspec((57, 35), rowspan=19, colspan=30)
        exp_bias_months_ax = self.gridspec.new_subplotspec((81, 35), rowspan=19, colspan=18)
        exp_bias_days_ax = self.gridspec.new_subplotspec((81, 56), rowspan=19, colspan=9)
        station_metadata_ax = self.gridspec.new_subplotspec((55, 0),  rowspan=45, colspan=30)

        # create subplot axes on grid
        # self.map_ax =              self.figure.add_subplot(map_ax)
        self.legend_ax = self.figure.add_subplot(legend_ax)
        self.ts_ax = self.figure.add_subplot(ts_ax)
        self.violin_hours_ax = self.figure.add_subplot(violin_hours_ax)
        self.violin_months_ax = self.figure.add_subplot(violin_months_ax)
        self.violin_days_ax = self.figure.add_subplot(violin_days_ax)
        self.exp_bias_hours_ax = self.figure.add_subplot(exp_bias_hours_ax)
        self.exp_bias_months_ax = self.figure.add_subplot(exp_bias_months_ax)
        self.exp_bias_days_ax = self.figure.add_subplot(exp_bias_days_ax)
        self.station_metadata_ax = self.figure.add_subplot(station_metadata_ax)

        # map colorbar ax
        self.cbar_ax = self.figure.add_axes([0.02, 0.525, 0.25, 0.0175])

        # Set map variables
        # create variable to create map axis on first data read
        self.map_initialised = False

        # define projections for map plot and actual geographic coordinates
        self.datacrs = ccrs.PlateCarree()
        self.plotcrs = ccrs.Robinson()

        # Turning off specific spines of time series axis
        self.ts_ax.spines["top"].set_visible(False)
        self.ts_ax.spines["right"].set_visible(False)

        # Hide all axes
        self.cbar_ax.axis('off')
        self.legend_ax.axis('off')
        self.ts_ax.axis('off')
        self.violin_hours_ax.axis('off')
        self.violin_months_ax.axis('off')
        self.violin_days_ax.axis('off')
        self.exp_bias_hours_ax.axis('off')
        self.exp_bias_months_ax.axis('off')
        self.exp_bias_days_ax.axis('off')
        self.station_metadata_ax.axis('off')

        # load basic statistics dictionary from conf, set in self
        # for rest of functions to see
        # json.load(open('providentia/conf/basic_stats_dict.json'))
        self.bstats_dict = self.read_instance.basic_stats_dict

        # same for experiment bias statistics
        # json.load(open('providentia/conf/experiment_bias_stats_dict.json'))
        self.expbias_dict = self.read_instance.expbias_dict

        # Define dictionary for mapping days of week/months as integers to
        # equivalent strings for writing on axes
        self.temporal_axis_mapping_dict = aux.temp_axis_dict()

    def update_MPL_canvas(self):
        """Function that updates MPL canvas upon clicking
        the 'READ' button, and when colocating data
        """

        # clear all axes (except map)
        self.legend_ax.cla()
        self.ts_ax.cla()
        self.violin_hours_ax.cla()
        self.violin_months_ax.cla()
        self.violin_days_ax.cla()
        self.exp_bias_hours_ax.cla()
        self.exp_bias_months_ax.cla()
        self.exp_bias_days_ax.cla()
        self.station_metadata_ax.cla()

        # hide all axes (except map)
        self.legend_ax.axis('off')
        self.ts_ax.axis('off')
        self.violin_hours_ax.axis('off')
        self.violin_months_ax.axis('off')
        self.violin_days_ax.axis('off')
        self.exp_bias_hours_ax.axis('off')
        self.exp_bias_months_ax.axis('off')
        self.exp_bias_days_ax.axis('off')
        self.station_metadata_ax.axis('off')

        # add map axis if map not yet initialised (i.e. first time loading some data)
        if not self.map_initialised:

            # create map axis
            map_ax = self.gridspec.new_subplotspec((0, 0), rowspan=45, colspan=45)
            self.map_ax = self.figure.add_subplot(map_ax, projection=self.plotcrs)
            # set map extents
            self.map_ax.set_global()

            # add coastlines and gridelines
            land_polygon_resolutions = {'low': '110m', 'medium': '50m', 'high': '10m'}
            feature = cfeature.NaturalEarthFeature('physical', 'land',
                                                   land_polygon_resolutions[self.read_instance.map_coastline_resolution],
                                                   facecolor='0.85')
            self.map_ax.add_feature(feature)
            self.map_ax.gridlines(linestyle='-', alpha=0.4)

            # reset the navigation toolbar stack for the map axis with the current view limits
            self.reset_ax_navigation_toolbar_stack(self.map_ax)

            # setup interactive picker/lasso on map
            self.figure.canvas.mpl_connect('pick_event', self.on_click)
            self.lasso = LassoSelector(self.map_ax, onselect=self.onlassoselect,
                                       useblit=True, lineprops=dict(alpha=0.5, color='hotpink', linewidth=1))
            # initialise variable that informs whether to use picker/lasso for updates
            self.map_already_updated = False

            # initialise variable of valid station indices plotted on map as empty list
            self.active_map_valid_station_inds = np.array([], dtype=np.int)

            # update variable to indicate map is now initialised
            self.map_initialised = True

        # define all temporal aggregation resolutions that will be used to aggregate data
        # (variable by temporal resolution of data in memory)
        if 'hourly' in self.read_instance.active_resolution:
            self.temporal_aggregation_resolutions = ['hour', 'dayofweek', 'month']
        elif self.read_instance.active_resolution == 'daily':
            self.temporal_aggregation_resolutions = ['dayofweek', 'month']
        elif self.read_instance.active_resolution == 'monthly':
            self.temporal_aggregation_resolutions = ['month']

        # reset relative index lists of selected station on map as empty lists
        self.previous_relative_selected_station_inds = np.array([], dtype=np.int)
        self.relative_selected_station_inds = np.array([], dtype=np.int)
        self.absolute_selected_station_inds = np.array([], dtype=np.int)

        # calculate map z statistic (for selected z statistic) --> updating active map valid station indices
        self.calculate_z_statistic()

        # update plotted map z statistic
        self.update_map_z_statistic()

        # plot experiment grid domain edges on map
        self.update_experiment_grid_domain_edges()

        # update legend
        self.update_legend()

        # draw changes
        self.draw()

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
        selected lower/upper limit bounds, selected measurement
        methods and selected minimum data availability %"""

        # return if nothing has been loaded yet
        if not hasattr(self.read_instance.datareader, 'data_in_memory'):
            return

        check_state = self.read_instance.ch_colocate.checkState()
        # update variable to inform plotting functions whether to use colocated data/or not
        if check_state == QtCore.Qt.Checked:
            self.colocate_active = True
        else:
            self.colocate_active = False

        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        if self.filter_data is None:
            self.filter_data = DataFilter(self.read_instance)
        else:
            self.filter_data.filter_all()
            self.update_active_map()
        QtWidgets.QApplication.restoreOverrideCursor()

    def update_active_map(self):
        if not self.read_instance.block_MPL_canvas_updates:
            # calculate map z statistic (for selected z statistic) --> updating active map valid station indices
            self.calculate_z_statistic()

            # update plotted map z statistic
            self.update_map_z_statistic()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(self.previous_relative_selected_station_inds,
                                  self.relative_selected_station_inds):
                self.update_associated_selected_station_plots()

            # draw changes
            self.draw()

    def handle_colocate_update(self):
        """Function that handles the update of the MPL canvas
        with colocated data upon checking of the colocate checkbox"""

        if not self.read_instance.block_MPL_canvas_updates:

            # if only have 1 data array in memory (i.e. observations), no colocation is possible,
            # therefore set colocation_active to be False, and return
            if len(list(self.read_instance.data_in_memory_filtered.keys())) == 1:
                self.colocate_active = False
                return

            # else, if have loaded experiment data, check if colocate checkbox is checked or unchecked
            check_state = self.read_instance.ch_colocate.checkState()
            # update variable to inform plotting functions whether to use colocated data/or not
            if check_state == QtCore.Qt.Checked:
                self.colocate_active = True
            else:
                self.colocate_active = False

            # update z statistic/experiment bias comboboxes (without updating canvas)
            self.read_instance.block_MPL_canvas_updates = True
            self.handle_map_z_statistic_update()
            self.handle_experiment_bias_update()
            self.read_instance.block_MPL_canvas_updates = False

            # update plotted map z statistic (if necessary)
            # calculate map z statistic (for selected z statistic) --> updating active map valid station indices
            self.calculate_z_statistic()

            # update plotted map z statistic
            self.update_map_z_statistic()

            # update associated plots with selected stations
            self.update_associated_selected_station_plots()

            # draw changes
            self.draw()

    def update_map_z_statistic(self):
        """Function that updates plotted z statistic on map, with colourbar"""

        # clear previously plotted station points
        try:
            self.map_points.remove()
        except:
            pass

        # clear/hide current colourbar axis
        self.cbar_ax.cla()
        self.cbar_ax.axis('off')

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
            self.map_points = self.map_ax.scatter(self.read_instance.datareader.station_longitudes[self.active_map_valid_station_inds],
                                                  self.read_instance.datareader.station_latitudes[self.active_map_valid_station_inds],
                                                  s=self.read_instance.unsel_station_markersize, c=self.z_statistic,
                                                  vmin=self.z_vmin, vmax=self.z_vmax, cmap=self.z_colourmap, picker=1,
                                                  zorder=3, transform=self.datacrs, linewidth=0.0, alpha=None)
            # create 2D numpy array of plotted station coordinates
            self.map_points_coordinates = np.vstack((self.read_instance.datareader.station_longitudes[self.active_map_valid_station_inds], self.read_instance.datareader.station_latitudes[self.active_map_valid_station_inds])).T

            # create colour normalisation instance
            colour_norm = matplotlib.colors.Normalize(vmin=self.z_vmin, vmax=self.z_vmax)
            # normalise z statistic to colourmap range
            norm_z_statisitic = colour_norm(self.z_statistic)
            # get list of rgba tuples to set map point colours
            self.rgba_tuples = matplotlib.cm.get_cmap(self.z_colourmap)(norm_z_statisitic)

            # colourbar
            # generate colourbar tick array
            tick_array = np.linspace(self.z_vmin, self.z_vmax, 5, endpoint=True)
            # plot colourbar
            cb = self.figure.colorbar(self.map_points, orientation='horizontal', cax=self.cbar_ax, label='',
                                      ticks=tick_array, extend='max')
            self.cbar_ax.tick_params(labelsize=8.0)
            # plot colourbar label
            self.cbar_ax.yaxis.set_label_position("right")
            self.cbar_ax.set_ylabel(self.z_label, fontsize=8.0, rotation=0, ha='left', va='top')
            # turn colourbar axis on
            self.cbar_ax.axis('on')

            # call update of map drawing (this is a hack to force map plot object to be updated
            # correctly --> only done when draw is called)
            self.draw()

        # update map selection appropriately for z statistic
        self.update_map_station_selection()

    def update_map_station_selection(self):
        """Function that updates the visual selection of stations on map"""

        # update map title
        if len(self.relative_selected_station_inds) == 1:
            self.map_ax.set_title('%s Selected' % (
                self.read_instance.station_references[self.relative_selected_station_inds][0]), fontsize=8.5, pad=3)
        else:
            self.map_ax.set_title('%s Stations Selected of %s Available' % (len(
                self.relative_selected_station_inds), len(self.active_map_valid_station_inds)), fontsize=8.5, pad=3)

        # reset alphas of all plotted stations (if have some stations on map)
        if len(self.active_map_valid_station_inds) > 0:
            self.rgba_tuples[:, -1] = 1.0
            marker_sizes = np.full(len(self.z_statistic), self.read_instance.unsel_station_markersize)
            self.map_points.set_facecolor(self.rgba_tuples)
            # reset marker sizes of all plotted stations
            self.map_points.set_sizes(marker_sizes)

            # if have some selected stations, update map plot with station selection
            # (reducing alpha of non-selected stations, and increasing marker size of selected stations)
            if len(self.relative_selected_station_inds) > 0:

                # decrease alpha of non-selected stations
                absolute_non_selected_stations = np.nonzero(~np.in1d(range(self.z_statistic.shape[0]),
                                                                     self.absolute_selected_station_inds))[0]
                if len(absolute_non_selected_stations) > 0:
                    self.rgba_tuples[absolute_non_selected_stations, -1] = 0.25
                    self.map_points.set_facecolor(self.rgba_tuples)

                # increase marker size of selected stations
                marker_sizes[self.absolute_selected_station_inds] = self.read_instance.sel_station_markersize
                self.map_points.set_sizes(marker_sizes)

    def update_associated_selected_station_plots(self):
        """Function that updates all plots associated with selected stations on map"""

        # clear/hide relevant axes
        # clear axes
        self.ts_ax.cla()
        self.violin_hours_ax.cla()
        self.violin_months_ax.cla()
        self.violin_days_ax.cla()
        self.exp_bias_hours_ax.cla()
        self.exp_bias_months_ax.cla()
        self.exp_bias_days_ax.cla()
        self.station_metadata_ax.cla()

        # hide axes
        self.ts_ax.axis('off')
        self.violin_hours_ax.axis('off')
        self.violin_months_ax.axis('off')
        self.violin_days_ax.axis('off')
        self.exp_bias_hours_ax.axis('off')
        self.exp_bias_months_ax.axis('off')
        self.exp_bias_days_ax.axis('off')
        self.station_metadata_ax.axis('off')

        # if have selected stations on map, then now remake plots
        if len(self.relative_selected_station_inds) > 0:

            # put selected data for each data array into pandas dataframes
            self.to_pandas_dataframe()

            # temporally aggregate selected data dataframes (by hour, day of week, month)
            self.pandas_temporal_aggregation()

            # if have some experiment data associated with selected stations, calculate
            # temporally aggregated basic statistic differences and bias statistics between
            # observations and experiment data arrays
            if len(list(self.selected_station_data.keys())) > 1:
                self.calculate_temporally_aggregated_experiment_bias_statistics()

            # update time series plot
            self.update_time_series_plot()

            # update violin plots for temporally aggregated data
            self.update_violin_plots()

            # if have some experiment data associated with selected stations,
            # update experiment bias plots for temporally aggregated statistical
            # differences/biases between observations and experiments
            if len(list(self.selected_station_data.keys())) > 1:
                self.update_experiment_bias_aggregated_plots()

            # update plotted station selected metadata
            self.update_selected_station_metadata()

    def update_experiment_grid_domain_edges(self):
        """Function that plots grid domain edges of experiments in memory"""

        # iterate through previously plotted experiment domain edge polygons, clearing them
        try:
            for grid_edge_polygon in self.grid_edge_polygons:
                grid_edge_polygon.remove()
        except:
            pass

        # create new array for storing plotted experiment domain edge polygons
        self.grid_edge_polygons = []

        # iterate through read experiments and plot grid domain edges on map
        for experiment in sorted(list(self.read_instance.datareader.data_in_memory.keys())):
            if experiment != 'observations':
                # compute map projection coordinates for each pair of
                # longitude/latitude experiment grid edge coordinates
                # exp_x,exp_y = self.bm(self.read_instance.datareader.data_in_memory[experiment]['grid_edge_longitude'],
                # self.read_instance.datareader.data_in_memory[experiment]['grid_edge_latitude'])
                # create matplotlib polygon object from experiment grid edge map projection coordinates
                grid_edge_outline_poly = \
                    Polygon(np.vstack((self.read_instance.datareader.plotting_params[experiment]['grid_edge_longitude'],
                                       self.read_instance.datareader.plotting_params[experiment]['grid_edge_latitude'])).T,
                                       edgecolor=self.read_instance.datareader.plotting_params[experiment]['colour'],
                                       linewidth=1, linestyle='--', fill=False, zorder=2, transform=self.datacrs)
                # plot grid edge polygon on map
                self.grid_edge_polygons.append(self.map_ax.add_patch(grid_edge_outline_poly))

    def update_legend(self):
        """Function that updates legend"""

        # create legend elements
        # add observations element
        legend_elements = [Line2D([0], [0], marker='o', color='white',
                                  markerfacecolor=self.read_instance.datareader.plotting_params['observations']['colour'],
                                  markersize=self.read_instance.legend_markersize, label='observations')]
        # add element for each experiment
        for experiment_ind, experiment in enumerate(sorted(list(self.read_instance.datareader.data_in_memory.keys()))):
            if experiment != 'observations':
                # add experiment element
                legend_elements.append(Line2D([0], [0], marker='o', color='white',
                                              markerfacecolor=self.read_instance.datareader.plotting_params[experiment]['colour'],
                                              markersize=self.read_instance.legend_markersize, label=experiment))

        # plot legend
        self.legend_ax.legend(handles=legend_elements, loc='best', mode='expand', ncol=4, fontsize=9.0)

    def to_pandas_dataframe(self):
        """Function that takes selected station data within data arrays and puts it into a pandas dataframe"""

        # create new dictionary to store selection station data by data array
        self.selected_station_data = {}

        # iterate through data arrays in data in memory filtered dictionary
        for data_label in list(self.read_instance.data_in_memory_filtered.keys()):

            # if colocation is not active, do not convert colocated data arrays to pandas data frames
            if not self.colocate_active:
                if 'colocated' in data_label:
                    continue
            # else, if colocation is active, do not convert non-colocated data arrays to pandas data frames
            elif self.colocate_active:
                if 'colocated' not in data_label:
                    continue

            # observational arrays
            if data_label.split('_')[0] == 'observations':
                # get data for selected stations
                data_array = self.read_instance.data_in_memory_filtered[data_label][
                                 self.read_instance.active_species][self.relative_selected_station_inds,:]
            # experiment arrays
            else:
                # get intersect between selected station indices and valid available indices for experiment data array
                valid_selected_station_indices = np.intersect1d(self.relative_selected_station_inds, self.read_instance.datareader.plotting_params[data_label]['valid_station_inds'])
                # get data for valid selected stations
                data_array = \
                    self.read_instance.data_in_memory_filtered[data_label][self.read_instance.active_species][valid_selected_station_indices,:]

            # if data array has no valid data for selected stations, do not create a pandas dataframe
            # data array has valid data?
            if data_array.size:
                # add nested dictionary for data array name to selection station data dictionary
                self.selected_station_data[data_label] = {}
                # take cross station median of selected data for data array, and place it in a pandas
                # dataframe -->  add to selected station data dictionary
                self.selected_station_data[data_label]['pandas_df'] = pd.DataFrame(np.nanmedian(data_array, axis=0),
                                                                                   index=self.read_instance.time_array,
                                                                                   columns=['data'])

    def pandas_temporal_aggregation(self):
        """Function that aggregates pandas dataframe data, for all data arrays,
        into desired temporal groupings also calculates all defined basic
        statistics for each individual temporal grouping
        """

        # define statistics to calculate (all basic statistics)
        statistics_to_calculate = list(self.bstats_dict.keys())

        # iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            # define all possible xticks for temporal resolution
            if temporal_aggregation_resolution == 'hour':
                all_xticks = np.arange(24, dtype=np.int)
            elif temporal_aggregation_resolution == 'dayofweek':
                all_xticks = np.arange(7, dtype=np.int)
            elif temporal_aggregation_resolution == 'month':
                all_xticks = np.arange(1, 13, dtype=np.int)

            # iterate through data arrays names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):

                # create nested dictionary inside selected station data dictionary for storing
                # aggregated data by data array label and temporal aggregation resolution
                self.selected_station_data[data_label][temporal_aggregation_resolution] = {}

                # aggregate data array into desired temporal groups (dropping NaNs)
                grouped_data = [g['data'].dropna() for n, g in self.selected_station_data[data_label]['pandas_df'].groupby(getattr(self.selected_station_data[data_label]['pandas_df'].index, temporal_aggregation_resolution))]
                # drop groups which have no data
                grouped_data = [group for group in grouped_data if len(group) > 0]
                # get xticks for groups which have valid data (i.e. which hours/days/months have valid data)
                valid_xticks = [getattr(group.index, temporal_aggregation_resolution)[0] for group in grouped_data]
                # create array of the size of the full range of each aggregation periods,
                # initialised with empty lists per element (i.e. 24 for hourly aggregation)
                full_grouped_data = [[] for _ in range(len(all_xticks))]
                # place valid grouped data in correct positions within full array
                full_group_indices_to_place = np.array([np.where(all_xticks == valid_xtick)[0][0]
                                                        for valid_xtick in valid_xticks], dtype=np.int)
                for grouped_data_ii, full_group_index_to_place in enumerate(full_group_indices_to_place):
                    full_grouped_data[full_group_index_to_place] = grouped_data[grouped_data_ii]

                # add full grouped data to selected data dictionary
                self.selected_station_data[data_label][temporal_aggregation_resolution]['grouped_data'] = full_grouped_data
                # add valid xticks for group to selected data dictionary
                # (i.e. the group xtick indexes which have valid data)
                self.selected_station_data[data_label][temporal_aggregation_resolution]['valid_xticks'] = valid_xticks

                # calculate basic statistics in each group and add them to selected station data dictionary
                for stat in statistics_to_calculate:
                    # get specific statistic dictionary (containing necessary
                    # information for calculation of selected statistic)
                    stats_dict = self.bstats_dict[stat]
                    # load default statistic arguments for passing to statistical function
                    function_arguments = stats_dict['arguments']
                    # create empty array for storing calculated statistic by group
                    stat_output_by_group = []
                    # iterate through grouped data
                    for group in full_grouped_data:
                        # only calculate statistic if have valid data in group
                        if len(group) > 0:
                            # calculate statistic (appending to all group statistic output array)
                            stat_output_by_group = np.append(
                                stat_output_by_group, getattr(Stats, stats_dict['function'])(group, **function_arguments))
                        # if no valid data in group, append NaN
                        else:
                            stat_output_by_group = np.append(stat_output_by_group, np.NaN)
                    # save statistical output by group to selected station data dictionary
                    self.selected_station_data[data_label][temporal_aggregation_resolution][stat] = stat_output_by_group

    def calculate_temporally_aggregated_experiment_bias_statistics(self):
        """Function that calculates temporally aggregated basic statistic
        differences and bias statistics between observations and experiment
        data arrays
        """

        # define all basic statistics that will be subtracted
        # (each experiment - observations) for each temporal aggregation
        basic_statistics = list(self.bstats_dict.keys())
        # define all experiment bias statistics that will be calculated
        # between each experiment and observations for each temporal aggregation
        bias_statistics = list(self.expbias_dict.keys())

        # iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            # iterate through data arrays names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):

                # make sure the data array is an experimental one
                if data_label.split('_')[0] != 'observations':

                    # get relevant aggregated observational statistics dictionary (i.e. colocated or not)
                    if not self.colocate_active:
                        relevant_aggregated_observations_dict = \
                            self.selected_station_data['observations'][temporal_aggregation_resolution]
                    else:
                        exp = data_label.split('_colocatedto_')[0]
                        relevant_aggregated_observations_dict = \
                            self.selected_station_data[
                                'observations_colocatedto_%s' % exp][temporal_aggregation_resolution]

                    # calculate temporally aggregated basic statistic differences between experiment and observations
                    # iterate through basic statistics
                    for basic_stat in basic_statistics:

                        # create empty array for storing calculated experiment-observations
                        # difference statistic by group
                        stat_diff_by_group = []

                        # iterate through aggregated index times
                        for group_ii in range(len(relevant_aggregated_observations_dict[basic_stat])):

                            # get observational and experiment group aggregated statistic
                            group_obs_stat = relevant_aggregated_observations_dict[basic_stat][group_ii]
                            group_exp_stat = self.selected_station_data[data_label][
                                temporal_aggregation_resolution][basic_stat][group_ii]

                            # take difference between observations and experiment statistics, if both values finite
                            if (np.isfinite(group_obs_stat)) & (np.isfinite(group_exp_stat)):
                                # calculate difference statistic (experiment - observations)
                                stat_diff_by_group = np.append(stat_diff_by_group, group_exp_stat-group_obs_stat)
                            # else, if one (or both) of observations/experiment statistics are NaN, append NaN
                            else:
                                stat_diff_by_group = np.append(stat_diff_by_group, np.NaN)
                        # save statistical difference output by group to selected station data dictionary
                        self.selected_station_data[data_label][temporal_aggregation_resolution][
                            '%s_bias' % (basic_stat)] = stat_diff_by_group

                    # if colocation is active, calculate temporally aggregated experiment bias
                    # statistical differences between experiment and observations
                    if self.colocate_active:

                        # iterate through bias statistics
                        for bias_stat in bias_statistics:

                            # get specific statistic dictionary (containing necessary information
                            # for calculation of selected statistic)
                            stats_dict = self.expbias_dict[bias_stat]
                            # load default statistic arguments for passing to statistical function
                            function_arguments = stats_dict['arguments']
                            # create empty array for storing calculated statistic by group
                            stat_output_by_group = []
                            # iterate through experimental grouped data
                            for group_ii, exp_group in enumerate(self.selected_station_data[data_label]
                                                                 [temporal_aggregation_resolution]['grouped_data']):
                                # add aggregated observational and experiment group data as
                                # arguments to pass to statistical function
                                function_arguments['obs'] = relevant_aggregated_observations_dict['grouped_data'][
                                    group_ii]
                                function_arguments['exp'] = exp_group

                                # calculate experiment bias statistic between observations and
                                # experiment group data 
                                try:
                                    # calculate experiment bias statistic
                                    stat_output_by_group = \
                                        np.append(stat_output_by_group,
                                                  getattr(ExpBias, stats_dict['function'])(**function_arguments))
                                # if can not calculate statistic, append NaN
                                except:
                                    stat_output_by_group = np.append(stat_output_by_group, np.NaN)
                            # save experiment bias statistic by group to selected station data dictionary
                            self.selected_station_data[data_label][temporal_aggregation_resolution][
                                '%s_bias' % (bias_stat)] = stat_output_by_group

    def update_time_series_plot(self):
        """Function that updates time series plot upon selection of station/s"""

        # turn axis on
        self.ts_ax.axis('on')

        # iterate through data array names in selected station data dictionary
        for data_label in list(self.selected_station_data.keys()):

            # if colocation is active, for observations use the 'observations_colocatedto_experiments',
            if self.colocate_active:
                if data_label.split('_')[0] == 'observations':
                    if data_label != 'observations_colocatedto_experiments':
                        continue
                # print only relevant data, otherwise we get double lines for multiple exps
                if data_label.split('_')[-1] not in ("observations", "experiments"):
                    continue

            # plot time series data
            self.data_array_ts = \
                self.ts_ax.plot(self.selected_station_data[data_label]['pandas_df'].dropna(),
                                color=self.read_instance.datareader.plotting_params[data_label]['colour'],
                                marker='o', markeredgecolor=None, mew=0,
                                markersize=self.read_instance.time_series_markersize,
                                linestyle='None',
                                zorder=self.read_instance.datareader.plotting_params[data_label]['zorder'])

        # set axes labels
        if self.read_instance.datareader.measurement_units == 'unitless':
            self.ts_ax.set_ylabel('{}'.format(self.read_instance.datareader.measurement_units), fontsize=8.0)
        else:
            self.ts_ax.set_ylabel('{} ({})'.format(self.read_instance.parameter_dictionary[self.read_instance.active_species]['axis_label'],self.read_instance.datareader.measurement_units), fontsize=8.0)

        # plot grid
        self.ts_ax.grid(color='lightgrey', alpha=0.8)

        # set axis tick label sizes
        self.ts_ax.tick_params(labelsize=8.0)

        # as are re-plotting on time series axis, reset the navigation
        # toolbar stack dictionaries entries associated with time series axis
        self.reset_ax_navigation_toolbar_stack(self.ts_ax)

    def update_violin_plots(self):
        """function that updates violin plots of temporally aggregated data upon selection of station/s"""

        # define dictionaries defining the relevant axis, axis titles, x axis ticks
        # (and an empty nested dictionary for storing plot objects) for different temporal aggregation resolutions
        hour_aggregation_dict = {
            'ax': self.violin_hours_ax, 'title': 'H', 'xticks': np.arange(24, dtype=np.int), 'plots': {}}
        dayofweek_aggregation_dict = {
            'ax': self.violin_days_ax, 'title': 'DoW', 'xticks': np.arange(7, dtype=np.int), 'plots': {}}
        month_aggregation_dict = {
            'ax': self.violin_months_ax, 'title': 'M', 'xticks': np.arange(1, 13, dtype=np.int), 'plots': {}}

        # based on the temporal resolution of the data, combine the relevant temporal aggregation dictionaries
        if 'hourly' in self.read_instance.active_resolution:
            aggregation_dict = {'hour': hour_aggregation_dict, 'dayofweek': dayofweek_aggregation_dict,
                                'month': month_aggregation_dict}
        elif self.read_instance.active_resolution == 'daily':
            aggregation_dict = {'dayofweek': dayofweek_aggregation_dict, 'month': month_aggregation_dict}
        elif self.read_instance.active_resolution == 'monthly':
            aggregation_dict = {'month': month_aggregation_dict}

        # turn on all axes that will be plotted on, and add yaxis grid to each axis, and change axis label tick sizes
        for temporal_aggregation_resolution in list(aggregation_dict.keys()):
            # turn on axis
            aggregation_dict[temporal_aggregation_resolution]['ax'].axis('on')
            # add yaxis grid
            aggregation_dict[temporal_aggregation_resolution]['ax'].yaxis.grid(color='lightgrey', alpha=0.8)
            # add axis aggregation resolution label
            aggregation_dict[temporal_aggregation_resolution]['ax'].\
                annotate(aggregation_dict[temporal_aggregation_resolution]['title'], (0, 1),
                         xytext=(2, -2), xycoords='axes fraction', textcoords='offset points',
                         fontsize=9.0, ha='left', va='top')
            # change axis tick labels
            aggregation_dict[temporal_aggregation_resolution]['ax'].tick_params(labelsize=8.0)

        # now, make violin plots for each temporally aggregated data array,
        # for all relevant temporal aggregation resolutions

        # iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            # create arrays for storing all calculated aggregated p5 and p95
            # (across all data arrays) for later limiting ylim
            all_p5 = []
            all_p95 = []

            # iterate through data array names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):

                # if colocation is active, for plotting observational aggregated data,
                # use the 'observations_colocatedto_experiments' array
                if self.colocate_active:
                    if data_label.split('_')[0] == 'observations':
                        if data_label != 'observations_colocatedto_experiments':
                            continue
                    # print only relevant data, otherwise we get double lines for multiple exps
                    if data_label.split('_')[-1] not in ("observations", "experiments"):
                        continue

                # get grouped data for current temporal aggregation resolution
                grouped_data = self.selected_station_data[data_label][temporal_aggregation_resolution]['grouped_data']
                # drop any groups which have no data
                grouped_data = [group for group in grouped_data if len(group) > 0]

                # make violin plot --> add plotted object to aggregation dictionary
                aggregation_dict[temporal_aggregation_resolution]['plots'][data_label] = \
                    aggregation_dict[temporal_aggregation_resolution]['ax'].violinplot(grouped_data, positions=self.selected_station_data[data_label][temporal_aggregation_resolution]['valid_xticks'], points=100, widths=0.85, showmeans=False, showmedians=False, showextrema=False)

                # append aggregated p5/p95 for data array, to all_p5/all_p95 arrays
                all_p5 = np.append(
                    all_p5, self.selected_station_data[data_label][temporal_aggregation_resolution]['p5'])
                all_p95 = np.append(
                    all_p95, self.selected_station_data[data_label][temporal_aggregation_resolution]['p95'])

            # set x axis limits
            aggregation_dict[temporal_aggregation_resolution]['ax'].set_xlim(np.min(aggregation_dict[temporal_aggregation_resolution]['xticks'])-0.5, np.max(aggregation_dict[temporal_aggregation_resolution]['xticks'])+0.5)
            # set y axis limits (use minimum p5 and maximum p95 across aggregations, across all data arrays)
            aggregation_dict[temporal_aggregation_resolution]['ax'].set_ylim(np.nanmin(all_p5), np.nanmax(all_p95))
            # set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
            if temporal_aggregation_resolution == 'hour':
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(aggregation_dict[temporal_aggregation_resolution]['xticks'][::3])
            else:
                aggregation_dict[temporal_aggregation_resolution]['ax'].\
                    set_xticks(aggregation_dict[temporal_aggregation_resolution]['xticks'])
                aggregation_dict[temporal_aggregation_resolution]['ax'].\
                    set_xticklabels([self.temporal_axis_mapping_dict[temporal_aggregation_resolution][xtick]
                                     for xtick in aggregation_dict[temporal_aggregation_resolution]['xticks']])

        # format violin plots
        # iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            # iterate through data arrays have plotted data
            for data_label in list(aggregation_dict[temporal_aggregation_resolution]['plots'].keys()):

                # get plotted object per data array
                violin_plot = aggregation_dict[temporal_aggregation_resolution]['plots'][data_label]

                # update plotted objects with necessary colour, zorder and alpha
                for patch in violin_plot['bodies']:
                    patch.set_facecolor(self.read_instance.datareader.plotting_params[data_label]['colour'])
                    patch.set_zorder(self.read_instance.datareader.plotting_params[data_label]['zorder'])
                    if data_label.split('_')[0] == 'observations':
                        patch.set_alpha(0.7)
                    else:
                        patch.set_alpha(0.5)
                    # if have at least 1 experiment data array, split the violin plot across the horizontal
                    # (observations on left, experiments on right)
                    if len(list(aggregation_dict[temporal_aggregation_resolution]['plots'].keys())) > 1:
                        m = np.mean(patch.get_paths()[0].vertices[:, 0])
                        # observations on left
                        if data_label.split('_')[0] == 'observations':
                            patch.get_paths()[0].vertices[:, 0] = np.clip(
                                patch.get_paths()[0].vertices[:, 0], -np.inf, m)
                        # experiments on right
                        else:
                            patch.get_paths()[0].vertices[:, 0] = np.clip(
                                patch.get_paths()[0].vertices[:, 0], m, np.inf)

                # overplot time series of medians over boxes in necessary color

                # generate zorder to overplot medians in same order as violin plots are ordered, but on top of them
                median_zorder = (self.read_instance.datareader.plotting_params['observations']['zorder']+len(
                    list(aggregation_dict[temporal_aggregation_resolution]['plots'].keys())) - 1) + \
                                self.read_instance.datareader.plotting_params[data_label]['zorder']

                # get xticks (all valid aggregated time indexes) and medians to plot
                xticks = aggregation_dict[temporal_aggregation_resolution]['xticks']
                medians = self.selected_station_data[data_label][temporal_aggregation_resolution]['p50']
                # split arrays if there are any temporal gaps to avoid
                # line drawn being interpolated across missing values
                inds_to_split = np.where(np.diff(xticks) > 1)[0]
                if len(inds_to_split) == 0:
                    aggregation_dict[temporal_aggregation_resolution]['ax'].plot(
                        xticks, medians, marker='o',
                        color=self.read_instance.datareader.plotting_params[data_label]['colour'],
                        markersize=self.read_instance.temp_agg_markersize, linewidth=0.5, zorder=median_zorder)
                else:
                    inds_to_split += 1
                    start_ind = 0
                    for end_ind in inds_to_split:
                        aggregation_dict[temporal_aggregation_resolution]['ax'].plot(
                            xticks[start_ind:end_ind], medians[start_ind:end_ind],
                            marker='o', color=self.read_instance.datareader.plotting_params[data_label]['colour'],
                            markersize=self.read_instance.temp_agg_markersize, linewidth=0.5, zorder=median_zorder)
                        start_ind = end_ind
                    aggregation_dict[temporal_aggregation_resolution]['ax'].plot(
                        xticks[start_ind:], medians[start_ind:], marker='o',
                        color=self.read_instance.datareader.plotting_params[data_label]['colour'],
                        markersize=self.read_instance.temp_agg_markersize, linewidth=0.5, zorder=median_zorder)

        # plot title (with units)
        # if selected data resolution is 'hourly', plot the title on off the hourly aggregation axis
        if 'hourly' in self.read_instance.active_resolution:
            self.violin_hours_ax.set_title('Temporal Distributions (%s)' % self.read_instance.datareader.measurement_units,
                                           fontsize=8.0, loc='left')
        # otherwise, plot the units on the monthly aggregation axis
        else:
            self.violin_months_ax.set_title('Temporal Distributions (%s)' % self.read_instance.datareader.measurement_units,
                                            fontsize=8.0, loc='left')

        # as are re-plotting on violin plot axes, reset the navigation toolbar stack
        # dictionaries entries associated with each of the axes
        for temporal_aggregation_resolution in list(aggregation_dict.keys()):
            self.reset_ax_navigation_toolbar_stack(aggregation_dict[temporal_aggregation_resolution]['ax'])

    def update_experiment_bias_aggregated_plots(self):
        """Function that updates the temporally aggregated experiment bias statistic plots"""

        # get currently selected experiment bias statistic
        selected_stat = self.read_instance.cb_experiment_bias_stat.currentText()
        # get name of array to retrieve pre-calculated bias for selected statistic
        selected_experiment_bias_stat = '%s_bias' % selected_stat

        # define dictionaries defining the relevant axis, axis titles, x axis ticks
        # (and an empty nested dictionary for storing plot objects) for different temporal aggregation resolutions
        hour_aggregation_dict = {
            'ax': self.exp_bias_hours_ax, 'title': 'H', 'xticks': np.arange(24, dtype=np.int), 'plots': {}}
        dayofweek_aggregation_dict = {
            'ax': self.exp_bias_days_ax, 'title': 'DoW', 'xticks': np.arange(7, dtype=np.int), 'plots': {}}
        month_aggregation_dict = {
            'ax': self.exp_bias_months_ax, 'title': 'M',   'xticks': np.arange(1, 13, dtype=np.int), 'plots': {}}

        # based on the temporal resolution of the data, combine the relevant temporal aggregation dictionaries
        if 'hourly' in self.read_instance.active_resolution:
            aggregation_dict = {
                'hour': hour_aggregation_dict, 'dayofweek': dayofweek_aggregation_dict, 'month': month_aggregation_dict}
        elif self.read_instance.active_resolution == 'daily':
            aggregation_dict = {'dayofweek': dayofweek_aggregation_dict, 'month': month_aggregation_dict}
        elif self.read_instance.active_resolution == 'monthly':
            aggregation_dict = {'month': month_aggregation_dict}

        # turn on all axes that will be plotted on, and add yaxis grid to each axis, and change axis label tick sizes
        for temporal_aggregation_resolution in list(aggregation_dict.keys()):
            # turn on axis
            aggregation_dict[temporal_aggregation_resolution]['ax'].axis('on')
            # add yaxis grid
            aggregation_dict[temporal_aggregation_resolution]['ax'].yaxis.grid(color='lightgrey',alpha=0.8)
            # add axis aggregation resolution label
            aggregation_dict[temporal_aggregation_resolution]['ax'].annotate(
                aggregation_dict[temporal_aggregation_resolution]['title'], (0, 1),
                xytext=(2, -2), xycoords='axes fraction', textcoords='offset points', fontsize=9, ha='left', va='top')
            # change axis tick labels
            aggregation_dict[temporal_aggregation_resolution]['ax'].tick_params(labelsize=8.0)

        # plot title (with units)
        # create title string to plot (based on type of statistic plotting)
        if selected_stat not in self.read_instance.basic_z_stats:
            stats_dict = self.expbias_dict[selected_stat]
            plot_title = 'Experiment %s' % (stats_dict['label'])
        else:
            stats_dict = self.bstats_dict[selected_stat]
            if selected_stat != 'Data %':
                title_units = ' (%s)' % self.read_instance.datareader.measurement_units
            else:
                title_units = ''
            plot_title = 'Experiment %s bias%s' % (stats_dict['label'], title_units)

        # if selected data resolution is 'hourly', plot the title on off the hourly aggregation axis
        if self.read_instance.active_resolution == 'hourly':
            self.exp_bias_hours_ax.set_title(plot_title, fontsize=8.0, loc='left')
        # otherwise, plot the units on the monthly aggregation axis
        else:
            self.exp_bias_months_ax.set_title(plot_title, fontsize=8.0, loc='left')

        # get value/s of minimum bias for statistic
        minimum_bias = stats_dict['minimum_bias']

        # ------------------------------------------------------------------------------------------------# 
        # now, make experiment bias plots for active bias statistic for all relevant temporal aggregation resolutions

        # iterate through all defined temporal aggregation resolutions
        for temporal_aggregation_resolution in self.temporal_aggregation_resolutions:

            # iterate through data array names in selected station data dictionary
            for data_label in list(self.selected_station_data.keys()):
                # if colocation is selected, print only relevant labels
                # to avoid double lines being printed
                if self.colocate_active:
                    if data_label.split('_')[-1] != 'observations':
                        continue
                # if data array is observational, continue to next experiment data array
                if data_label.split('_')[0] == 'observations':
                    continue
                # else, make temporally aggregated plot for currently active experiment bias statistic
                else:
                    aggregation_dict[temporal_aggregation_resolution]['plots'][data_label] = \
                        aggregation_dict[temporal_aggregation_resolution]['ax'].plot(
                            aggregation_dict[temporal_aggregation_resolution]['xticks'],
                            self.selected_station_data[data_label][temporal_aggregation_resolution]
                            [selected_experiment_bias_stat],
                            color=self.read_instance.datareader.plotting_params[data_label]['colour'],
                            marker='o', zorder=self.read_instance.datareader.plotting_params[data_label]['zorder'],
                            markersize=self.read_instance.temp_agg_expbias_markersize, linewidth=0.5)

            # set x axis limits
            aggregation_dict[temporal_aggregation_resolution]['ax'].set_xlim(
                np.min(aggregation_dict[temporal_aggregation_resolution]['xticks'])-0.5,
                np.max(aggregation_dict[temporal_aggregation_resolution]['xticks'])+0.5)
            # set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
            if temporal_aggregation_resolution == 'hour':
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(
                    aggregation_dict[temporal_aggregation_resolution]['xticks'][::3])
            else:
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticks(
                    aggregation_dict[temporal_aggregation_resolution]['xticks'])
                aggregation_dict[temporal_aggregation_resolution]['ax'].set_xticklabels(
                    [self.temporal_axis_mapping_dict[temporal_aggregation_resolution][xtick] for xtick
                     in aggregation_dict[temporal_aggregation_resolution]['xticks']])

            # plot horizontal line/s across x axis at value/s of minimum experiment bias
            for mb in minimum_bias:
                aggregation_dict[temporal_aggregation_resolution]['ax'].axhline(y=mb, linestyle='--', linewidth=1.0,
                                                                                color='black', zorder=0)

        # as are re-plotting on experiment bias axes,
        # reset the navigation toolbar stack dictionaries entries associated with each of the axes
        for temporal_aggregation_resolution in list(aggregation_dict.keys()):
            self.reset_ax_navigation_toolbar_stack(aggregation_dict[temporal_aggregation_resolution]['ax'])

    def update_selected_station_metadata(self):
        """Function which updates the plotted metadata
        detail of selected stations on the map
        """
        # get some details of the station metadata axis --> to set limit for wrapping text
        # get axis bounding box
        ax_bbox = self.station_metadata_ax.get_window_extent().transformed(self.figure.dpi_scale_trans.inverted())
        # get axis dimensions in inches
        ax_width_inches = ax_bbox.width
        # get axis dimensions in pixels
        ax_width_px = ax_width_inches * self.figure.dpi

        # initialise string to plot on axis
        str_to_plot = ''

        # define dictionary with metadata variable names paired with metadata variable names to plot
        metadata_variable_naming = {'station_reference': 'Station Reference',
                                    'station_name': 'Station Name',
                                    'latitude': 'Latitude',
                                    'longitude': 'Longitude',
                                    'measurement_altitude': 'Measurement Altitude',
                                    'country': 'Country',
                                    'network': 'Network',
                                    'area_classification': 'Area Classification',
                                    'station_classification': 'Station Classification',
                                    'main_emission_source': 'Main Emission Source',
                                    'land_use': 'Land Use',
                                    'terrain': 'Terrain',
                                    'measurement_scale': 'Measurement Scale',
                                    'representative_radius': 'Representative Radius',
                                    'ESDAC_Iwahashi_landform_classification': 'Iwahashi Landform',
                                    'ESDAC_Meybeck_landform_classification': 'Meybeck Landform',
                                    'Joly-Peuch_classification_code': 'Joly-Peuch Class',
                                    'Koppen-Geiger_classification': 'Koppen-Geiger Class',
                                    'MODIS_MCD12C1_v6_IGBP_land_use': 'MODIS Land Use',
                                    'WMO_region': 'WMO Region',
                                    'WWF_TEOW_terrestrial_ecoregion': 'WWF TEOW Terrestrial Ecoregion',
                                    'WWF_TEOW_biogeographical_realm': 'WWF TEOW Biogeographical Realm',
                                    'WWF_TEOW_biome': 'WWF TEOW Biome',
                                    'GHSL_settlement_model_classification': 'GHSL Settlement Class',
                                    'GSFC_coastline_proximity': 'To Coast',
                                    'primary_sampling_type': 'Sampling Instrument Type',
                                    'sample_preparation_types': 'Sample Preparation',
                                    'measurement_methodology': 'Measurement Method',
                                    'measuring_instrument_name': 'Measuring Instrument',
                                    'measuring_instrument_sampling_type': 'Measuring Instrument Sampling'
                                    }

        # is just 1 station selected?
        if len(self.relative_selected_station_inds) == 1:

            # get station reference of selected station
            selected_station_reference = self.read_instance.station_references[self.relative_selected_station_inds][0]

            if self.read_instance.reading_nonghost:
                str_to_plot += "%s   " % (selected_station_reference)
                str_to_plot += "Latitude: {:.4f}   ".format(
                    self.read_instance.datareader.station_latitudes[self.relative_selected_station_inds][0])
                str_to_plot += "Longitude: {:.4f}\n".format(
                    self.read_instance.datareader.station_longitudes[self.relative_selected_station_inds][0])

            else:
                # add station reference, latitude, longitude, measurement altitude, GSFC coastline proximity,
                # GPW population density and NOAA-DMSP-OLS nighttime stable lights
                str_to_plot += "%s   " % (selected_station_reference)
                str_to_plot += "Latitude: {:.4f}   ".format(
                    self.read_instance.datareader.station_latitudes[self.relative_selected_station_inds][0])
                str_to_plot += "Longitude: {:.4f}\n".format(
                    self.read_instance.datareader.station_longitudes[self.relative_selected_station_inds][0])
                str_to_plot += "Measurement Altitude: {:.2f}m   ".format(np.nanmedian(
                    self.read_instance.datareader.metadata_in_memory['measurement_altitude'][
                        self.relative_selected_station_inds].astype(np.float32)))
                str_to_plot += "To Coast: {:.2f}km\n".format(np.nanmedian(
                    self.read_instance.datareader.metadata_in_memory['GSFC_coastline_proximity'][
                        self.relative_selected_station_inds].astype(np.float32)))
                str_to_plot += "Population Density: {:.1f} people/km–2\n".format(np.nanmedian(
                    self.read_instance.datareader.metadata_in_memory['GHSL_population_density'][
                        self.relative_selected_station_inds].astype(np.float32)))
                str_to_plot += "Nighttime Lights: {:.1f}\n".format(np.nanmedian(
                    self.read_instance.datareader.metadata_in_memory['NOAA-DMSP-OLS_v4_nighttime_stable_lights'][
                        self.relative_selected_station_inds].astype(np.float32)))

                # define other metadata variables to plot, in order to plot
                # (plotting all unique associated metadata values)
                metadata_vars_to_plot = ['station_name', 'country', 'area_classification',
                                         'station_classification', 'terrain',
                                         'land_use', 'MODIS_MCD12C1_v6_IGBP_land_use',
                                         'GHSL_settlement_model_classification',
                                         'measurement_methodology', 'measuring_instrument_name']

                # iterate through metadata variables
                for meta_var in metadata_vars_to_plot:

                    # gather all selected station metadata for current meta variable
                    all_current_meta = self.read_instance.datareader.metadata_in_memory[meta_var][self.relative_selected_station_inds].flatten().astype(np.str)

                    # get counts of all unique metadata elements for selected station
                    unique_meta, meta_counts = np.unique(all_current_meta, return_counts=True)
                    # get number of unique metadata elements across selected stations
                    n_unique_meta = len(unique_meta)

                    # 1 unique metadata element?
                    if n_unique_meta == 1:
                        meta_string = '{}: {}\n'.format(metadata_variable_naming[meta_var], unique_meta[0])
                    # elif have > 2 unique metadata elements, just return count of the elements for the selected station
                    elif n_unique_meta > 2:
                        meta_string = '{}: {} unique elements\n'.format(metadata_variable_naming[meta_var], n_unique_meta)
                    # otherwise, get percentage of unique metadata elements across selected stations
                    else:
                        meta_pc = (100./len(all_current_meta))*meta_counts
                        meta_pc = ['{:.1f}%'.format(meta) for meta in meta_pc]
                        #create string for variable to plot
                        meta_string = '{}: {}\n'.format(metadata_variable_naming[meta_var], ', '.join([':'.join([str(var),pc]) for var, pc in zip(unique_meta, meta_pc)]))

                    #add meta string to str_to_plot
                    str_to_plot += meta_string

        # more than 1 station selected?
        else:
            if self.read_instance.reading_nonghost:
                str_to_plot += "%s Stations Selected\n" % (len(self.relative_selected_station_inds))
            else:
                # TODO looks like this variable is not used
                # get station references of all selected stations
                selected_station_references = self.read_instance.station_references[self.relative_selected_station_inds]

                # add N stations selected, in N countries
                str_to_plot += "%s Stations Selected\n" % (len(self.relative_selected_station_inds))
                # add median measurement altitude
                str_to_plot += "Median Measurement Altitude: {:.2f}m   ".format(np.nanmedian(
                    self.read_instance.datareader.metadata_in_memory['measurement_altitude'][
                        self.relative_selected_station_inds].astype(np.float32)))
                # add median GSFC coastline proximity
                str_to_plot += "Median To Coast: {:.2f}km\n".format(np.nanmedian(
                    self.read_instance.datareader.metadata_in_memory['GSFC_coastline_proximity'][
                        self.relative_selected_station_inds].astype(np.float32)))
                # add median GPW population density
                str_to_plot += "Median Population Density: {:.1f} people/km–2\n".format(np.nanmedian(
                    self.read_instance.datareader.metadata_in_memory['GHSL_population_density'][
                        self.relative_selected_station_inds].astype(np.float32)))
                # add median NOAA-DMSP-OLS nighttime lights
                str_to_plot += "Median Nighttime Lights: {:.1f}\n".format(np.nanmedian(
                    self.read_instance.datareader.metadata_in_memory['NOAA-DMSP-OLS_v4_nighttime_stable_lights'][
                        self.relative_selected_station_inds].astype(np.float32)))

                # get percentage of element occurrences across selected stations, for certain metadata variables
                metadata_vars_get_pc = ['country', 'area_classification', 'station_classification',
                                        'terrain', 'land_use', 'MODIS_MCD12C1_v6_IGBP_land_use',
                                        'GHSL_settlement_model_classification',
                                        'measurement_methodology', 'measuring_instrument_name']

                # iterate through metadata variables
                for meta_var in metadata_vars_get_pc:

                    # gather all selected station metadata for current meta variable
                    all_current_meta = self.read_instance.datareader.metadata_in_memory[meta_var][
                        self.relative_selected_station_inds].flatten().astype(np.str)

                    # get counts of all unique metadata elements across selected stations
                    unique_meta, meta_counts = np.unique(all_current_meta, return_counts=True)
                    # get number of unique metadata elements across selected stations
                    n_unique_meta = len(unique_meta)

                    # if have > 4 unique metadata elements, just return count
                    # of the elements across the selected stations
                    if n_unique_meta > 4:
                        meta_string = '{}: {} unique elements\n'.format(metadata_variable_naming[meta_var], n_unique_meta)
                    # otherwise, get percentage of unique metadata elements across selected stations
                    else:
                        meta_pc = (100. / len(all_current_meta)) * meta_counts
                        meta_pc = ['{:.1f}%'.format(meta) for meta in meta_pc]
                        # create string for variable to plot
                        meta_string = '{}: {}\n'.format(metadata_variable_naming[meta_var], ', '.join(
                            [':'.join([str(var), pc]) for var, pc in zip(unique_meta, meta_pc)]))

                    # add meta string to str_to_plot
                    str_to_plot += meta_string

        # plot string to axis
        plot_txt = self.station_metadata_ax.text(0.0, 1.0, str_to_plot, ha='left', va='top', fontsize=8.0,
                                                 transform=self.station_metadata_ax.transAxes, wrap=True,
                                                 linespacing=1.5)
        # modify limit to wrap text as axis width in pixels  --> hack as matplotlib
        # automatically sets limit as figure width
        plot_txt._get_wrap_line_width = lambda: ax_width_px

    def calculate_z_statistic(self):
        """Function that calculates selected z statistic for map"""

        # get relevant observational array (dependent on colocation)
        if not self.colocate_active:
            obs_array = self.read_instance.datareader.plotting_params['observations']['valid_station_inds']
        else:
            obs_array = \
                self.read_instance.datareader.plotting_params['observations_colocatedto_experiments']['valid_station_inds']

        # before doing anything check if have any valid station data for observations,
        # if not update active map valid station indices to be empty list and return
        if len(obs_array) == 0:
            self.active_map_valid_station_inds = np.array([], dtype=np.int)
            return

        # get selected z1/z2 data arrays
        z1_selected_name = self.read_instance.cb_z1.currentText()
        z2_selected_name = self.read_instance.cb_z2.currentText()
        # check if have a selected z2 array (i.e. z2 name is not empty string)
        if z2_selected_name == '':
            have_z2 = False
        else:
            have_z2 = True

        # get dictionary containing necessary information for calculation of selected statistic
        # check if the chosen statistic is a basic statistic
        z_statistic_name = self.read_instance.cb_z_stat.currentText()
        if z_statistic_name in list(self.bstats_dict.keys()):
            z_statistic_type = 'basic'
            stats_dict = self.bstats_dict[z_statistic_name]
            # set label units for statistic
            if z_statistic_name != 'Data %':
                label_units = ' ({})'.format(self.read_instance.datareader.measurement_units)
            else:
                label_units = ''
        # if not a basic statistic, it must be an experiment bias statistic
        else:
            z_statistic_type = 'bias'
            stats_dict = self.expbias_dict[z_statistic_name]
            label_units = ''

        # set colourbar for z statistic
        # first check if have defined colourbar for z statistic, if so use that
        if 'colourbar' in list(stats_dict.keys()):
            self.z_colourmap = getattr(self.read_instance, stats_dict['colourbar'])
        # else, set appropriate colourmap for the type of statistic
        else:
            # if only have selected z1 array, the statistic is 'absolute', so use sequential colourbar
            if not have_z2:
                self.z_colourmap = self.read_instance.sequential_colourmap
            # if have selected z1 and z2 arrays, the statistic is 'difference', so use diverging colourbar
            else:
                self.z_colourmap = self.read_instance.diverging_colourmap

        # generate z colourbar label
        if not have_z2:
            self.z_label = '{}\n{} {}'.format(z1_selected_name, stats_dict['label'], label_units)
        else:
            self.z_label = '{} - {}\n{} {}'.format(z2_selected_name, z1_selected_name, stats_dict['label'], label_units)

        # if colocation is active, set appropriate z1/z2 arrays to read to get colocated data arrays
        if self.colocate_active:
            # don't have z2 array?
            if not have_z2:
                if z1_selected_name == 'observations':
                    z1_array_to_read = 'observations_colocatedto_experiments'
                else:
                    z1_array_to_read = '{}_colocatedto_observations'.format(z1_selected_name)
            # have z2 array?
            elif have_z2:
                if z1_selected_name == 'observations':
                    z1_array_to_read = 'observations_colocatedto_{}'.format(z2_selected_name)
                else:
                    if z2_selected_name == 'observations':
                        z1_array_to_read = '{}_colocatedto_observations'.format(z1_selected_name)
                    else:   
                        z1_array_to_read = '{}_colocatedto_{}'.format(z1_selected_name, z2_selected_name) 

                if z2_selected_name == 'observations':
                    z2_array_to_read = 'observations_colocatedto_{}'.format(z1_selected_name)
                else:
                    if z1_selected_name == 'observations':
                        z2_array_to_read = '{}_colocatedto_observations'.format(z2_selected_name)
                    else:
                        z2_array_to_read = '{}_colocatedto_{}'.format(z2_selected_name, z1_selected_name)

        # else, simply use selected z1/z2 array names to read uncolocated data arrays
        else:
            z1_array_to_read = copy.deepcopy(z1_selected_name)
            z2_array_to_read = copy.deepcopy(z2_selected_name)

        # read selected data arrays (after subsetting arrays by intersection of z1/z2 valid station indices)
        # and calculate desired Z statistic (after removing NaNs from arrays)

        # get active map valid station indices (i.e. the indices of the stations data to plot on the map)
        # if only have z1, valid map indices are those simply for the z1 array
        if not have_z2:
            self.active_map_valid_station_inds = \
                self.read_instance.datareader.plotting_params[z1_array_to_read]['valid_station_inds']
        else:
            # if have z2 array, get intersection of z1 and z2 valid station indices
            self.active_map_valid_station_inds = \
                np.intersect1d(self.read_instance.datareader.plotting_params[z1_array_to_read]['valid_station_inds'],
                               self.read_instance.datareader.plotting_params[z2_array_to_read]['valid_station_inds'])

        # update absolute selected plotted station indices with respect to new active map valid station indices
        self.absolute_selected_station_inds = np.array(
            [np.where(self.active_map_valid_station_inds == selected_ind)[0][0] for selected_ind in
             self.relative_selected_station_inds if selected_ind in self.active_map_valid_station_inds],
            dtype=np.int)

        # read z1 data
        z1_array_data = \
            self.read_instance.data_in_memory_filtered[z1_array_to_read][self.read_instance.active_species][self.active_map_valid_station_inds,:]
        # drop NaNs and reshape to object list of station data arrays (if not checking data %)
        if z_statistic_name != 'Data %':
            z1_array_data = drop_nans(z1_array_data)
        else:
            z1_array_data.tolist()

        # create empty array to store z statistic
        self.z_statistic = np.empty(len(z1_array_data))

        # if have no z2 data, calculate 'absolute' basic statistic
        if not have_z2:

            # load default selected z statistic arguments for passing to statistical function
            function_arguments = stats_dict['arguments']

            # iterate through stations calculating statistic
            for z_ii in range(len(self.z_statistic)):

                # and calculate its statistics
                self.z_statistic[z_ii] = \
                    getattr(Stats, stats_dict['function'])(z1_array_data[z_ii], **function_arguments)

        # else, read z2 data then calculate 'difference' statistic
        else:
            # read z2 data
            z2_array_data = \
                self.read_instance.data_in_memory_filtered[z2_array_to_read][self.read_instance.active_species][self.active_map_valid_station_inds,:]
            # drop NaNs and reshape to object list of station data arrays (if not checking data %)
            if z_statistic_name != 'Data %':
                z2_array_data = drop_nans(z2_array_data)
            else:
                z2_array_data = z2_array_data.tolist()
            # is the difference statistic basic (i.e. mean)?
            if z_statistic_type == 'basic':

                # load default selected z statistic arguments and make separate arguments
                # dictionaries for z1/z2 calculations (as doing 2 separate calculations for z1/z2 and subtracting)
                function_arguments_z1 = stats_dict['arguments']
                function_arguments_z2 = copy.deepcopy(function_arguments_z1)

                # iterate through stations calculating statistic
                for z_ii in range(len(self.z_statistic)):

                    # calculate statistics for z1/z2 arrays and subtract z2-z1
                    self.z_statistic[z_ii] = \
                        getattr(Stats, stats_dict['function'])(z2_array_data[z_ii], **function_arguments_z2) - \
                        getattr(Stats, stats_dict['function'])(z1_array_data[z_ii], **function_arguments_z1)

            # else, is the difference statistic an experiment bias statistic (i.e. r)?
            elif z_statistic_type == 'bias':

                # load default selected z statistic arguments for passing to statistical function
                function_arguments = stats_dict['arguments']

                # iterate through stations calculating statistic
                for z_ii in range(len(self.z_statistic)):

                    # set station z1/z2 arrays as arguments in argument dictionary
                    function_arguments['obs'] = z1_array_data[z_ii]
                    function_arguments['exp'] = z2_array_data[z_ii]

                    # calculate statistic
                    self.z_statistic[z_ii] = getattr(ExpBias, stats_dict['function'])(**function_arguments)

        # if any station z statistics come out as NaN/inf, remove respective
        # stations from active map valid station indices
        # also cut z_statistic to remove invalid NaNs/infs

        valid_z_statistic_boolean = np.isfinite(self.z_statistic) 
        self.active_map_valid_station_inds = self.active_map_valid_station_inds[valid_z_statistic_boolean]
        self.z_statistic = self.z_statistic[valid_z_statistic_boolean]

        # calculate z vmin/vmax
        # if have defined z vmin OR vmax for chosen statistic, use them by priority
        # (if not taking a difference with a 'basic' statistic)

        # check if have defined vmin
        have_defined_vmin = False
        if 'vmin' in list(stats_dict.keys()):
            # if z statistic type is 'basic' and taking a difference, then DO NOT use defined vmin
            if (z_statistic_type == 'basic') & (not have_z2):
                have_defined_vmin = True
            elif z_statistic_type == 'bias':
                have_defined_vmin = True
        # have defined zmin?
        if have_defined_vmin:
            self.z_vmin = stats_dict['vmin']
        # else, take vmin as minimum range value of calculated statistic
        else:
            self.z_vmin = np.nanmin(self.z_statistic)

        # check if have defined vmax
        have_defined_vmax = False
        if 'vmax' in list(stats_dict.keys()):
            # if z statistic type is 'basic' and taking a difference, then DO NOT use defined vmax
            if (z_statistic_type == 'basic') & (not have_z2):
                have_defined_vmax = True
            elif z_statistic_type == 'bias':
                have_defined_vmax = True
        # have defined zmax?
        if have_defined_vmax:
            self.z_vmax = stats_dict['vmax']
        # else, take vmax as maximum range value of calculated statistic
        else:
            self.z_vmax = np.nanmax(self.z_statistic)

        # if z statistic is a 'difference', and do not have one of vmin/vmax pre-defined,
        # force vmin/vmax to be symmetrical across 0
        if have_z2 & (not have_defined_vmin) & (not have_defined_vmax):
            limit_stat = np.max(np.abs([self.z_vmin, self.z_vmax]))
            self.z_vmin = -limit_stat
            self.z_vmax = limit_stat

    def handle_map_z_statistic_update(self):
        """Define function which handles update of map z statistic"""

        if not self.read_instance.block_config_bar_handling_updates:

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
            if (not self.colocate_active) or (selected_z2_array == ''):
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

                # calculate map z statistic (for selected z statistic) --> updating active map valid station indices
                self.calculate_z_statistic()

                # update plotted map z statistic
                self.update_map_z_statistic()

                # draw changes
                self.draw()

    def handle_experiment_bias_update(self):
        """Define function that handles update of plotted experiment bias statistics"""

        if not self.read_instance.block_config_bar_handling_updates:

            # if no experiment data loaded, do not update
            if len(self.read_instance.experiment_bias_types) > 0:

                # update experiment bias comboboxes

                # set variable that blocks configuration bar handling updates until all changes
                # to the experiment bias comboboxes are made
                self.read_instance.block_config_bar_handling_updates = True

                # get currently selected items
                selected_experiment_bias_type = self.read_instance.cb_experiment_bias_type.currentText()
                selected_experiment_bias_stat = self.read_instance.cb_experiment_bias_stat.currentText()

                # update experiment bias statistics (used for Aggregated field), to all basic stats
                # if colocation not-active, and basic+bias stats if colocation active
                if not self.colocate_active:
                    available_experiment_bias_stats = copy.deepcopy(self.read_instance.basic_z_stats)
                else:
                    available_experiment_bias_stats = copy.deepcopy(self.read_instance.basic_and_bias_z_stats)

                # if selected bias type is empty string, it is because fields are being initialised for the first time
                if selected_experiment_bias_type == '':
                    # set experiment bias type to be first available type
                    selected_experiment_bias_type = self.read_instance.experiment_bias_types[0]
                    # set experiment bias stat to be first available stat
                    selected_experiment_bias_stat = available_experiment_bias_stats[0]
                    if hasattr(self.read_instance, 'exp_bias_stat'):
                        if self.read_instance.exp_bias_stat in available_experiment_bias_stats:
                            selected_experiment_bias_stat = self.read_instance.exp_bias_stat

                # if selected bias type is 'Rank', then there are no stat options so force the
                # available items and selected stat to be empty
                if selected_experiment_bias_type == 'Rank':
                    available_experiment_bias_stats = []
                    selected_experiment_bias_stat = ''

                # update all comboboxes (clear, then add items)
                self.read_instance.cb_experiment_bias_type.clear()
                self.read_instance.cb_experiment_bias_stat.clear()
                self.read_instance.cb_experiment_bias_type.addItems(self.read_instance.experiment_bias_types)
                self.read_instance.cb_experiment_bias_stat.addItems(available_experiment_bias_stats)

                # update selected values
                self.read_instance.cb_experiment_bias_type.setCurrentText(selected_experiment_bias_type)
                # maintain currently selected bias statistic (if exists in new item list)
                if selected_experiment_bias_stat in available_experiment_bias_stats:
                    self.read_instance.cb_experiment_bias_stat.setCurrentText(selected_experiment_bias_stat)

                # allow handling updates to the configuration bar again
                self.read_instance.block_config_bar_handling_updates = False

                if not self.read_instance.block_MPL_canvas_updates:

                    # update experiment bias plot/s if have some stations selected on map
                    if len(self.relative_selected_station_inds) > 0:

                        # clear and turn off all relevant axes before updating
                        self.exp_bias_hours_ax.cla()
                        self.exp_bias_months_ax.cla()
                        self.exp_bias_days_ax.cla()
                        self.exp_bias_hours_ax.axis('off')
                        self.exp_bias_months_ax.axis('off')
                        self.exp_bias_days_ax.axis('off')

                        # if experiment bias type == 'Aggregated' --> update plotted experiment bias plots
                        if selected_experiment_bias_type == 'Aggregated':
                            self.update_experiment_bias_aggregated_plots()

                        # draw changes
                        self.draw()

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
            self.draw()

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
            self.draw()

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
            # if len(self.previous_relative_selected_station_inds) > 0:
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
        self.draw()

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
        self.draw()

    def map_selected_station_inds_to_all_available_inds(self, selected_map_inds):
        """Takes the indices of selected stations on the map
        (potentially a subset of all available stations), and returns the indices
        of the stations inside the full loaded data arrays
        """

        # index the array of indices of stations plotted on the map (indexed with respect to
        # all available stations), with the absolute indices of the subset of plotted selected stations
        return self.active_map_valid_station_inds[selected_map_inds]

