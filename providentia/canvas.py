""" Class for Dashboard matplotlib canvas """

import copy
import functools
import inspect
import datetime
import math
import sys
import yaml
from weakref import WeakKeyDictionary

import cartopy.crs as ccrs
import matplotlib
from matplotlib.backend_bases import MouseButton
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.path import Path
import matplotlib.style as mplstyle
import numpy as np
from packaging.version import Version
import pandas as pd
from pandas.plotting import register_matplotlib_converters
from PyQt5 import QtCore, QtGui, QtWidgets

from providentia.auxiliar import (
    CURRENT_PATH,
    join,
    COLOUR_PRESET_CUSTOM,
    get_colour_presets,
    get_map_colours,
    get_role_colourmap,
)
from .canvas_menus import SettingsMenu, set_slider_enabled
from .dashboard_elements import ComboBox
from .dashboard_elements import set_formatting, set_cursor, unset_cursor
from .dashboard_elements import populate_colourmap_combobox, select_colourmap
from .dashboard_elements import (
    populate_projection_combobox,
    populate_colour_combobox,
    LAND_COLOUR_OPTIONS,
    OCEAN_COLOUR_OPTIONS,
)
from .dashboard_interactivity import HoverAnnotation
from .dashboard_interactivity import (
    legend_picker_func,
    picker_block_func,
    zoom_map_func,
)
from .fields_menus import update_metadata_fields
from .filter import DataFilter
from .plotting import Plotting
from .plot_aux import (
    get_map_extent,
    get_map_marker_size,
    download_plot_data_to_csv,
)
from .plot_formatting import (
    format_axis,
    fit_boxplot_xticklabels,
    harmonise_xy_lims_paradigm,
    log_validity,
    set_axis_label,
    set_axis_title,
    draw_map_features,
    remove_map_features,
    draw_map_gridlines,
)
from .plot_options import annotation, linear_regression, log_axes, smooth, threshold
from .read_aux import get_possible_resampling_resolutions, get_frequency_code
from .statistics import (
    get_z_statistic_comboboxes,
    generate_colourbar,
    get_selected_station_data,
    get_z_statistic_type,
    get_z_statistic_info,
    resolve_colourmap,
    get_colourmap_role,
)
from .warnings_prv import show_message

# make sure that we are using Qt5 backend with matplotlib
matplotlib.use("Qt5Agg")
register_matplotlib_converters()

# use matplotlib fast style: https://matplotlib.org/stable/users/explain/performance.html
mplstyle.use("fast")

PROVIDENTIA_ROOT = "/".join(CURRENT_PATH.split("/")[:-1])
settings_dict = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings/internal/canvas_menus.yaml"))
)

# tuning constants for automatic map marker opacity - see
# Canvas.apply_automatic_marker_style(). Size comes from the shared
# get_map_marker_size() in plot_aux.py, so it matches report and library.
MAP_AUTO_SIZING_MIN_OPACITY = 0.4
MAP_AUTO_SIZING_MAX_OPACITY = 1.0
# zoom ratio (current view's linear scale vs the projection's full global
# extent) at which marker size/opacity reach their maximum - beyond this,
# they're clamped rather than continuing to grow
MAP_AUTO_SIZING_REFERENCE_ZOOM = 15.0
# selected stations render this much bigger than unselected ones - large,
# deliberately, since size/opacity alone (no edge/outline) now carry the
# whole "selected vs unselected" distinction
MAP_AUTO_SIZING_SELECTED_SIZE_BOOST = 2.9
# once there's an active selection, unselected stations dim to this
# fraction of their normal opacity, so the selection reads clearly
MAP_AUTO_SIZING_UNSELECTED_OPACITY_DIM = 0.28


def restores_settings_guard(method):
    """
    Decorator which guarantees a settings handler leaves the dashboard usable,
    however it exits.

    Each handler raises block_config_bar_handling_updates while it works, so the
    controls it changes don't re-trigger each other, and every handler returns
    immediately when it is already set. Lowering it as the last statement rather
    than in a finally meant any exception in between left it raised for the rest
    of the session, silently disabling the whole settings menu. It is restored
    to whatever it was on entry, not simply cleared, so a handler called from
    inside another cannot lower the outer one's guard on its way out.

    Parameters
    ----------
    method : function
        Settings handler to wrap

    Returns
    -------
    function
        Wrapped handler
    """

    # how many arguments the handler itself takes, so that anything Qt adds
    # beyond them can be dropped. A signal hands its slot the value that
    # changed, and PyQt drops what the slot has no room for - but it reads the
    # slot's signature to know that, and through this wrapper every handler
    # looks as though it takes anything, so the trimming is done here instead
    parameters = list(inspect.signature(method).parameters.values())
    takes_anything = any(
        parameter.kind is parameter.VAR_POSITIONAL for parameter in parameters
    )
    accepted = (
        len(
            [
                parameter
                for parameter in parameters
                if parameter.kind
                in (parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD)
            ]
        )
        - 1
    )

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        if not takes_anything:
            args = args[:accepted]
        previous = self.read_instance.block_config_bar_handling_updates
        try:
            return method(self, *args, **kwargs)
        finally:
            self.read_instance.block_config_bar_handling_updates = previous
            # only restores if this call is the one that set it
            unset_cursor(self.read_instance.cursor_function, method.__name__)

    return wrapper


class Canvas(FigureCanvas):
    """
    Class that handles the creation and updates of
    a matplotlib canvas, and associated subplots
    """

    def __init__(self, read_instance):
        """
        Initialise the MPL canvas

        Parameters
        ----------
        read_instance : object
            Instance of class Dashboard or Report
        """

        # create figure and canvas objects
        self.figure = Figure(dpi=100)
        FigureCanvas.__init__(self, self.figure)

        # add passed arguments to self
        self.read_instance = read_instance

        # get characteristics per plot type
        self.plot_characteristics_templates = (
            self.read_instance.plot_characteristics_templates
        )
        self.plot_characteristics = {}

        # add general plot characteristics to self (if not passed via command line)
        for k, val in self.plot_characteristics_templates["general"].items():
            if k not in self.read_instance.commandline_arguments:
                setattr(self, k, val)

        # initialise some key vars
        self.filter_data = None

        # initialise Plotting class
        self.plotting = Plotting(read_instance=self.read_instance, canvas_instance=self)

        # setup gridding of canvas
        self.gridspec = gridspec.GridSpec(self.gridspec_nrows, self.gridspec_ncols)
        self.gridspec.update(**self.gridspec_format)

        # define all possible plots
        self.all_plots = [
            "legend",
            "map",
            "timeseries",
            "periodic-violin",
            "periodic",
            "metadata",
            "distribution",
            "histogram",
            "scatter",
            "statsummary",
            "boxplot",
            "taylor",
            "fairmode-target",
            "fairmode-statsummary",
            "contingencytable",
        ]

        # define all possible plots in layout options
        self.layout_options = [
            "None",
            "boxplot",
            "distribution",
            "histogram",
            "metadata",
            "periodic",
            "periodic-violin",
            "scatter",
            "statsummary",
            "timeseries",
            "taylor",
            "fairmode-target",
            "fairmode-statsummary",
            "contingencytable",
        ]

        # stop running if plot type in active_dashboard_plots does not exist
        for plot_type in self.read_instance.active_dashboard_plots:
            if plot_type not in self.all_plots + ["None"]:
                error = "Error: Plot type {0} is not an option. ".format(plot_type)
                error += "The available plots are: {0}.".format(self.all_plots[2:])
                self.read_instance.logger.error(error)
                sys.exit(1)

        # initialize layout positions
        self.read_instance.position_1 = "map"
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
        self.read_instance.current_statsummary_stats["basic"] = {}
        self.read_instance.current_statsummary_stats["modbias"] = {}
        for periodic_cycle in ["None", "Diurnal", "Weekly", "Monthly"]:
            self.read_instance.current_statsummary_stats["basic"][periodic_cycle] = []
            self.read_instance.current_statsummary_stats["modbias"][periodic_cycle] = []

        # update plot characteristics for all plots, and initialise plot options per plot type
        for plot_type in self.all_plots:
            # for plot types with zstat, initialise with default zstat (Mean)
            if plot_type in ["map", "periodic"]:
                self.plotting.set_plot_characteristics(
                    [plot_type], zstat="Mean", data_labels=["dummy"]
                )
            else:
                self.plotting.set_plot_characteristics(
                    [plot_type], data_labels=["dummy"]
                )
            self.current_plot_options[plot_type] = []
            self.previous_plot_options[plot_type] = []

        # create map, colorbar and legend plot axes
        self.plot_axes = {}
        self.plot_axes["map"] = self.figure.add_subplot(
            self.gridspec.new_subplotspec((2, 0), rowspan=44, colspan=42),
            projection=self.plotcrs,
        )
        self.plot_axes["cb"] = self.figure.add_axes([0.0255, 0.536, 0.3794, 0.02])
        self.plot_axes["legend"] = self.figure.add_subplot(
            self.gridspec.new_subplotspec((0, 47), rowspan=8, colspan=53)
        )

        # add settings menus
        self.generate_interactive_elements()

        # create rest of plot axes (default: timeseries, statsummary, distribution, periodic)
        # also show plot type buttons
        for position, plot_type in enumerate(self.read_instance.active_dashboard_plots):
            # update plot axis
            self.read_instance.update_plot_axis(self, position + 2, plot_type)

            # gather menu, save buttons and elements for plot type
            for menu_button, save_button, save_data_button in zip(
                self.menu_buttons, self.save_buttons, self.save_data_buttons
            ):
                menu_plot_type = menu_button.objectName().split("_menu")[0]
                if plot_type in [
                    "periodic_violin",
                    "fairmode_target",
                    "fairmode_statsummary",
                ]:
                    plot_type = plot_type.replace("_", "-")
                # proceed once have objects for plot type
                if plot_type == menu_plot_type:
                    menu_button.show()
                    save_button.show()
                    save_data_button.show()

        # show map buttons
        self.map_menu.buttons["settings_button"].show()
        self.map_menu.buttons["save_button"].show()
        self.map_menu.buttons["save_data_button"].show()

        # update layout fields
        self.read_instance.update_layout_fields(self)

        # initialise variable of valid station indices plotted on map as empty list
        self.active_map_valid_station_inds = np.array([], dtype=np.int32)

        # setup blocker for picker events
        self.axes_enter_event = self.figure.canvas.mpl_connect(
            "axes_enter_event", lambda event: picker_block_func(self, event)
        )

        # setup legend line selection
        self.legend_pick = self.figure.canvas.mpl_connect(
            "pick_event", lambda event: legend_picker_func(self, event)
        )

        # setup picker for station selection (left and right click)
        self.lasso_active = False
        self.station_pick = self.figure.canvas.mpl_connect(
            "button_press_event", self.station_select
        )

        # setup canvas annotations
        self.canvas_annotation = HoverAnnotation(self)
        self.canvas_annotation_event = self.figure.canvas.mpl_connect(
            "motion_notify_event", self.canvas_annotation.hover_annotation
        )

        # setup zoom on scroll wheel on map
        self.zoom_map_event = self.figure.canvas.mpl_connect(
            "scroll_event", lambda event: zoom_map_func(self, event)
        )

        # format axes for map, legend and active_dashboard_plots
        for plot_type in ["map", "legend"] + self.read_instance.active_dashboard_plots:
            if plot_type != "None":
                format_axis(
                    self.read_instance,
                    self,
                    self.plot_axes[plot_type],
                    plot_type,
                    self.plot_characteristics[plot_type],
                    map_extent=self.read_instance.map_extent,
                )

        # create covers to hide parts of canvas when updating / plotting
        self.canvas_cover = set_formatting(
            QtWidgets.QWidget(self), self.read_instance.formatting_dict["canvas_cover"]
        )
        self.top_right_canvas_cover = set_formatting(
            QtWidgets.QWidget(self), self.read_instance.formatting_dict["canvas_cover"]
        )
        self.top_right_canvas_cover.hide()
        self.lower_canvas_cover = set_formatting(
            QtWidgets.QWidget(self), self.read_instance.formatting_dict["canvas_cover"]
        )
        self.lower_canvas_cover.hide()
        # place partial canvas covers below map elements
        for element in self.map_elements:
            element.raise_()

    def update_MPL_canvas(self):
        """
        Function that updates MPL canvas upon clicking
        the 'READ' button, and when colocating data
        """

        # reset relative index lists of selected station on map as empty lists
        self.previous_relative_selected_station_inds = np.array([], dtype=np.int32)
        self.relative_selected_station_inds = np.array([], dtype=np.int32)
        self.absolute_selected_station_inds = np.array([], dtype=np.int32)

        # reset plot_elements
        self.plot_elements = {}
        self.plot_elements["data_labels_active"] = []
        for data_label in self.read_instance.data_labels:
            self.plot_elements["data_labels_active"].append(data_label)

        # add map domain plot option if on first read
        if (self.read_instance.first_read) & (
            "domain" not in self.current_plot_options["map"]
        ):
            self.current_plot_options["map"].append("domain")

        # update legend
        self.update_legend()

        # update plotted map z statistic
        self.update_map_z_statistic()

        # uncover map, but hide plotting axes
        self.canvas_cover.hide()
        self.cover_plot_axes()

        # draw changes
        self.figure.canvas.draw_idle()

        return None

    def reset_ax_navigation_toolbar_stack(self, ax):
        """
        Function which resets the navigation toolbar stack
        for a given axis with the current view limits

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes
        """

        # get appropriate axes for nested axes
        axs_to_reset = []
        if isinstance(ax, dict):
            for relevant_temporal_resolution, sub_ax in ax.items():
                if (
                    relevant_temporal_resolution
                    in self.read_instance.periodic_relevant_temporal_resolutions
                ):
                    axs_to_reset.append(sub_ax)
        elif isinstance(ax, list):
            axs_to_reset = copy.copy(ax)
        else:
            axs_to_reset.append(ax)

        # check if have axes dictionaries in stack list
        for ax_to_reset in axs_to_reset:
            if len(self.read_instance.navi_toolbar._nav_stack) == 0:
                # if don't have an axes dictionary in stack list, create one with current
                # axis in dictionary with current view limits
                if ax_to_reset.get_figure():
                    self.read_instance.navi_toolbar._nav_stack.push(
                        WeakKeyDictionary(
                            {
                                ax_to_reset: (
                                    ax_to_reset._get_view(),
                                    (
                                        ax_to_reset.get_position(True).frozen(),
                                        ax_to_reset.get_position().frozen(),
                                    ),
                                )
                            }
                        )
                    )

            # if have existing axes dictionaries in stack list, iterate through stack list
            # removing given axis from all stack list dictionaries
            else:
                for axes_dict in self.read_instance.navi_toolbar._nav_stack:
                    if ax_to_reset in axes_dict.keyrefs():
                        axes_dict.pop(ax_to_reset)

                # now add axis to first dictionary in stack, with the current view limits
                if ax_to_reset.get_figure():
                    self.read_instance.navi_toolbar._nav_stack[0][ax_to_reset] = (
                        ax_to_reset._get_view(),
                        (
                            ax_to_reset.get_position(True).frozen(),
                            ax_to_reset.get_position().frozen(),
                        ),
                    )

        return None

    def handle_data_filter_update(self):
        """
        Function which handles updates of data filtering
        """

        # return if nothing has been loaded yet
        if not hasattr(self.read_instance, "data_in_memory"):
            return None

        # update mouse cursor to a waiting cursor
        self.read_instance.cursor_function = set_cursor(
            self.read_instance.cursor_function, "handle_data_filter_update"
        )

        # filter data
        # if filter class not yet intialised, then do so
        if self.filter_data is None:
            self.filter_data = DataFilter(self.read_instance)
        # if it is, update filters, and update map and associated plots
        else:
            self.filter_data.filter_all()
            self.update_active_map()

        # restore mouse cursor to normal
        unset_cursor(self.read_instance.cursor_function, "handle_data_filter_update")

        return None

    def handle_resampling_update(self):
        """
        Function which handles updates of resampling
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # return if nothing has been loaded yet
            if not hasattr(self.read_instance, "data_in_memory"):
                return None

            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_resampling_update"
            )

            # update resampling statistics
            self.update_resampling_statistics()

            # deactivate MPL canvas updates while updating things
            original_block_MPL_canvas_updates = copy.deepcopy(
                self.read_instance.block_MPL_canvas_updates
            )
            self.read_instance.block_MPL_canvas_updates = True

            # disable MDA8 stat where neccessary
            self.handle_statsummary_statistics_update()
            self.handle_periodic_statistic_update()
            self.update_timeseries_chunk_statistics()

            # # restore block_MPL_canvas_updates
            self.read_instance.block_MPL_canvas_updates = (
                original_block_MPL_canvas_updates
            )

            # update layout fields
            self.read_instance.update_layout_fields(self)

            # update plots?
            if not self.read_instance.block_MPL_canvas_updates:
                # update plotted map z statistic
                self.update_map_z_statistic()

                # if have selected stations on map, then now remake plots
                if hasattr(self, "relative_selected_station_inds"):
                    if len(self.relative_selected_station_inds) > 0:
                        # update associated plots with selected stations
                        self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(self.read_instance.cursor_function, "handle_resampling_update")

        return None

    @restores_settings_guard
    def update_resampling_statistics(self):
        """
        Update resampling statistics
        """

        # turn off handling updates to configuration bar
        self.read_instance.block_config_bar_handling_updates = True

        # get available resampling resolutions, removing base resolution
        available_resampling_resolutions = get_possible_resampling_resolutions(
            self.read_instance.resolution,
            daily_forecast=self.read_instance.daily_forecast,
        )
        available_resampling_resolutions.remove(
            self.read_instance.resolution.split("_instantaneous")[0]
        )

        # remove resolutions if resampled data would be less than 2 timesteps
        resampling_resolutions = copy.deepcopy(available_resampling_resolutions)
        for resampling_resolution in resampling_resolutions:
            # get active frequency code
            active_frequency_code = get_frequency_code(resampling_resolution)

            # get test time array of new resolution
            start_date = str(self.read_instance.start_date)
            end_date = str(self.read_instance.end_date)
            time_array = pd.date_range(
                start=datetime.datetime(
                    int(start_date[:4]), int(start_date[4:6]), int(start_date[6:8])
                ),
                end=datetime.datetime(
                    int(end_date[:4]), int(end_date[4:6]), int(end_date[6:8])
                ),
                freq=active_frequency_code,
            )[:-1]

            # remove resolution when the data consists only of less than 2 timesteps
            if len(time_array) < 2:
                available_resampling_resolutions.remove(resampling_resolution)

        # get current selected resampling resolution
        current_resampling_resolution = (
            self.read_instance.cb_resampling_resolution.currentText()
        )

        # update resampling resolution field
        available_resampling_resolutions = [
            "None",
        ] + available_resampling_resolutions
        self.read_instance.cb_resampling_resolution.clear()
        self.read_instance.cb_resampling_resolution.addItems(
            available_resampling_resolutions
        )

        # if currently selected resampling resolution in the available resampling options, set it as active again
        # otherwise it will be None
        if current_resampling_resolution in available_resampling_resolutions:
            self.read_instance.cb_resampling_resolution.setCurrentText(
                current_resampling_resolution
            )
        # show message stating that are setting resampling resolution to None
        else:
            msg = f"The resampling resolution will be set to 'None' as the active resampling resolution ({current_resampling_resolution}) is not allowed."
            show_message(self.read_instance, msg)

        # update resampling resolution
        self.read_instance.resampling_resolution = (
            self.read_instance.cb_resampling_resolution.currentText()
        )

        # set active resolution, resampling_resolution when set, otherwise resolution
        if self.read_instance.resampling_resolution != "None":
            self.read_instance.active_resolution = (
                self.read_instance.resampling_resolution
            )
        else:
            self.read_instance.active_resolution = self.read_instance.resolution

        # allow handling updates to the canvas and configuration bar again
        self.read_instance.block_config_bar_handling_updates = False

    def update_active_map(self):
        """
        Function that updates plotted map z statistic and updates associated plots
        """

        if not self.read_instance.block_MPL_canvas_updates:
            # make copy of current full array relative selected stations indices
            self.previous_relative_selected_station_inds = copy.deepcopy(
                self.relative_selected_station_inds
            )

            # update plotted map z statistic
            self.update_map_z_statistic()

            # update associated plots
            self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw_idle()

        return None

    @restores_settings_guard
    def handle_statistic_mode_update(self):
        """
        Function that handles the update of the MPL canvas
        when we change the statistical calculation mode
        """

        if not self.read_instance.block_MPL_canvas_updates:
            # set variable that blocks configuration bar handling updates until all
            # changes to the statistic mode are made
            self.read_instance.block_config_bar_handling_updates = True

            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_statistic_mode_update"
            )

            # update mode
            self.read_instance.selected_statistic_mode = (
                self.read_instance.cb_statistic_mode.currentText()
            )
            self.read_instance.statistic_mode = (
                self.read_instance.selected_statistic_mode
            )

            # update aggregation statistic
            self.update_aggregation_statistic()

            # update timeseries aggregation statistic
            self.update_timeseries_aggregation_statistic()

            # update chunk statistic
            self.update_timeseries_chunk_statistics()

            self.read_instance.block_config_bar_handling_updates = False

            # update associated plots with selected stations
            self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function, "handle_statistic_mode_update"
            )

        return None

    @restores_settings_guard
    def handle_statistic_aggregation_update(self):
        """
        Function that handles the update of the MPL canvas
        when we change the aggregation statistic
        """

        if not self.read_instance.block_MPL_canvas_updates:
            # set variable that blocks configuration bar handling updates until all
            # changes to the statistics are made
            self.read_instance.block_config_bar_handling_updates = True

            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function,
                "handle_statistic_aggregation_update",
            )

            # update aggregation statistic
            self.update_aggregation_statistic()

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update associated plots with selected stations
            self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function,
                "handle_statistic_aggregation_update",
            )

        return None

    def handle_temporal_colocate_update(self):
        """
        Function that handles the update of the MPL canvas
        with colocated data upon checking of the temporal colocate checkbox
        """

        # update mouse cursor to a waiting cursor
        self.read_instance.cursor_function = set_cursor(
            self.read_instance.cursor_function, "handle_temporal_colocate_update"
        )

        # else, if have loaded model data, check if colocate checkbox is checked or unchecked
        check_state = self.read_instance.ch_colocate.checkState()

        # update variable to inform plotting functions whether to use colocated data/or not
        # turn colocation on
        if check_state == QtCore.Qt.Checked:
            self.read_instance.temporal_colocation = True
            # need to update plots?
            if len(self.read_instance.data_labels) < 2:
                if self.read_instance.temporal_colocation_active:
                    update_plots = True
                else:
                    update_plots = False
                self.read_instance.temporal_colocation_active = False
            else:
                if self.read_instance.temporal_colocation_active:
                    update_plots = False
                else:
                    update_plots = True
                self.read_instance.temporal_colocation_active = True
        # turn colocation off
        else:
            self.read_instance.temporal_colocation = False
            # need to update plots?
            if self.read_instance.temporal_colocation_active:
                update_plots = True
            else:
                update_plots = False
            self.read_instance.temporal_colocation_active = False

        # update metadata fields
        if self.read_instance.station_cap:
            update_metadata_fields(self.read_instance, cap=True)
        else:
            update_metadata_fields(self.read_instance)

        # update layout fields
        self.read_instance.update_layout_fields(self)

        # update canvas plots (if neccessary)
        if (update_plots) or (self.read_instance.performing_read):
            self.read_instance.block_MPL_canvas_updates = True
            # update plot statistics
            self.handle_map_z_statistic_update()
            self.handle_timeseries_chunk_statistic_update()
            self.handle_periodic_statistic_update()
            self.handle_statsummary_statistics_update()
            self.handle_statsummary_cycle_update()
            self.handle_statsummary_periodic_aggregation_update()
            self.handle_statsummary_periodic_mode_update()
            if self.read_instance.temporal_colocation_active:
                self.handle_taylor_correlation_statistic_update()
                self.handle_fairmode_target_classification_update()
            self.read_instance.block_MPL_canvas_updates = False

            # if not performing read then update plots
            if not self.read_instance.performing_read:
                # update plotted map z statistic
                self.update_map_z_statistic()

                # update associated plots with selected stations
                self.update_associated_active_dashboard_plots()

                # draw changes
                self.figure.canvas.draw_idle()

        # restore mouse cursor to normal
        unset_cursor(
            self.read_instance.cursor_function, "handle_temporal_colocate_update"
        )

        return None

    def unselect_map_checkboxes(self):
        """
        Function to uncheck All, Intersect and Extent checkboxes without updating canvas
        """

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
        """
        Function that updates plotted z statistic on map, with colourbar
        """

        # remove axis elements from map/cb
        self.remove_axis_elements(self.plot_axes["map"], "map")
        self.remove_axis_elements(self.plot_axes["cb"], "cb")

        # check if labels that have set for map exist in current data labels
        # if not then reset map plot
        labela = self.map_z1.currentText()
        labelb = self.map_z2.currentText()
        if labela not in self.read_instance.data_labels:
            self.map_z1.setCurrentText(self.read_instance.observations_data_label)
            self.map_z2.setCurrentTextText("")
        elif (labelb not in self.read_instance.data_labels) & (labelb != ""):
            self.map_z1.setCurrentText(self.read_instance.observations_data_label)
            self.map_z2.setCurrentText("")

        # get zstat name from combobox
        base_zstat = self.map_z_stat.currentText()
        if self.map_z2.currentText() == "":
            zstat = get_z_statistic_comboboxes(base_zstat)
        else:
            zstat = get_z_statistic_comboboxes(base_zstat, bias=True)

        # ensure label that have in memory still exists

        # plot map for zstat --> updating active map valid station indices and setting up plot picker
        self.plotting.make_map(
            self.plot_axes["map"],
            self.read_instance.networkspeci,
            self.plot_characteristics["map"],
            self.current_plot_options["map"],
            zstat=zstat,
            labela=self.map_z1.currentText(),
            labelb=self.map_z2.currentText(),
        )

        # update absolute selected plotted station indices with respect to new active map valid station indices
        self.absolute_selected_station_inds = np.array(
            [
                np.where(self.active_map_valid_station_inds == selected_ind)[0][0]
                for selected_ind in self.relative_selected_station_inds
                if selected_ind in self.active_map_valid_station_inds
            ],
            dtype=np.int32,
        )

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
            self.previous_relative_selected_station_inds = copy.deepcopy(
                self.relative_selected_station_inds
            )
            self.relative_selected_station_inds = np.array([], dtype=np.int32)
            self.absolute_selected_station_inds = np.array([], dtype=np.int32)
            self.absolute_non_selected_station_inds = np.array([], dtype=np.int32)

        # else, if any of the currently selected stations are not in the current active map
        # valid station indices --> unselect selected stations (and associated plots)
        # also uncheck select all/intersect/extent checkboxes
        else:
            if not np.all(
                np.in1d(
                    self.relative_selected_station_inds,
                    self.active_map_valid_station_inds,
                )
            ):
                # unselect all/intersect/extent checkboxes
                self.unselect_map_checkboxes()

                # reset relative/absolute selected station indices to be empty lists
                self.previous_relative_selected_station_inds = copy.deepcopy(
                    self.relative_selected_station_inds
                )
                self.relative_selected_station_inds = np.array([], dtype=np.int32)
                self.absolute_selected_station_inds = np.array([], dtype=np.int32)

            # get absolute non-selected station inds
            self.absolute_non_selected_station_inds = np.nonzero(
                ~np.in1d(
                    range(len(self.active_map_valid_station_inds)),
                    self.absolute_selected_station_inds,
                )
            )[0]

            # create 2D numpy array of plotted station coordinates
            self.map_points_coordinates = np.vstack(
                (
                    self.read_instance.station_longitudes[
                        self.read_instance.networkspeci
                    ][self.active_map_valid_station_inds],
                    self.read_instance.station_latitudes[
                        self.read_instance.networkspeci
                    ][self.active_map_valid_station_inds],
                )
            ).T

            # generate colourbar
            resolved_vmin, resolved_vmax = generate_colourbar(
                self.read_instance,
                [self.plot_axes["map"]],
                [self.plot_axes["cb"]],
                zstat,
                self.plot_characteristics["map"],
                self.read_instance.species[0],
                cmap_override=getattr(self.read_instance, "map_colourmap_override", None),
                vmin_override=getattr(self.read_instance, "map_vmin_override", None),
                vmax_override=getattr(self.read_instance, "map_vmax_override", None),
                discrete_override=getattr(self.read_instance, "map_discrete_override", None),
                n_discrete_override=getattr(
                    self.read_instance, "map_n_discrete_override", None
                ),
                n_ticks_override=getattr(self.read_instance, "map_n_ticks_override", None),
            )

            # show the colourbar's actual resolved limits in the settings
            # menu's limit fields, so editing one starts from the real
            # current number rather than a blank field. Guarded, as this can
            # run before the map settings menu exists
            if hasattr(self, "map_cb_min") and (resolved_vmin is not None):
                self.map_cb_min.setText(f"{resolved_vmin:.4g}")
            if hasattr(self, "map_cb_max") and (resolved_vmax is not None):
                self.map_cb_max.setText(f"{resolved_vmax:.4g}")

        # update plot options
        self.update_plot_options(plot_types=["map"])

        # resolve each collection's per-point face colours from its scalar
        # mappable, which update_map_station_selection() reads to apply the
        # selected/unselected styling on top of. A full draw+flush here
        # painted that intermediate state to screen (the map appearing to
        # flash selected then unselect) and pumped the event loop
        # mid-handler, swallowing the click that triggered it
        for collection in self.plot_axes["map"].collections:
            if isinstance(collection, matplotlib.collections.PathCollection):
                collection.update_scalarmappable()

        # update map selection appropriately for z statistic
        self.update_map_station_selection()

        return None

    def update_map_station_selection(self):
        """
        Function that updates the visual selection of stations on map
        """

        # recompute automatic marker size/opacity (a no-op if automatic
        # sizing is off) before anything below reads the marker styles -
        # every full map redraw and selection change routes through here
        self.apply_automatic_marker_style()

        # update map title
        if len(self.relative_selected_station_inds) == 1:
            axis_title_label = "{} Selected".format(
                self.read_instance.station_references[self.read_instance.networkspeci][
                    self.relative_selected_station_inds[0]
                ]
            )
        else:
            axis_title_label = "{} Selected Stations of {} Available".format(
                len(self.relative_selected_station_inds),
                len(self.active_map_valid_station_inds),
            )
        set_axis_title(
            self.read_instance,
            self.plot_axes["map"],
            axis_title_label,
            self.plot_characteristics["map"],
        )
        self.plot_characteristics["map"]["axis_title"]["label"] = axis_title_label

        # reset alphas and marker sizes of stations (if have some stations on map)
        if len(self.active_map_valid_station_inds) > 0:
            # set markersize of all stations (initally assuming zero stations are selected)
            markersizes = np.full(
                len(self.active_map_valid_station_inds),
                self.plot_characteristics["map"]["marker_zero_stations_selected"]["s"],
            )

            for collection in self.plot_axes["map"].collections:
                if isinstance(collection, matplotlib.collections.PathCollection):
                    if Version(matplotlib.__version__) < Version("3.4"):
                        opacities = collection.get_facecolor()

                        # set alpha of all stations (initally assuming zero stations are selected)
                        opacities[:, -1] = self.plot_characteristics["map"][
                            "marker_zero_stations_selected"
                        ]["alpha"]

                        # have selected stations?
                        if len(self.relative_selected_station_inds) > 0:
                            # update markersize and alphas of non-selected stations
                            markersizes[
                                self.absolute_non_selected_station_inds
                            ] = self.plot_characteristics["map"]["marker_unselected"][
                                "s"
                            ]
                            opacities[
                                self.absolute_non_selected_station_inds, -1
                            ] = self.plot_characteristics["map"]["marker_unselected"][
                                "alpha"
                            ]

                            # update markersize and alphas of selected stations
                            markersizes[
                                self.absolute_selected_station_inds
                            ] = self.plot_characteristics["map"]["marker_selected"]["s"]
                            opacities[
                                self.absolute_selected_station_inds, -1
                            ] = self.plot_characteristics["map"]["marker_selected"][
                                "alpha"
                            ]

                        # set new markersizes and alphas
                        collection.set_sizes(markersizes)
                        collection.set_facecolor(opacities)

                    else:
                        opacities = collection.get_facecolor()[:, -1]

                        # set alpha of all stations (initally assuming zero stations are selected)
                        opacities[:] = self.plot_characteristics["map"][
                            "marker_zero_stations_selected"
                        ]["alpha"]

                        # have selected stations?
                        if len(self.relative_selected_station_inds) > 0:
                            # update markersize and alphas of non-selected stations
                            markersizes[
                                self.absolute_non_selected_station_inds
                            ] = self.plot_characteristics["map"]["marker_unselected"][
                                "s"
                            ]
                            opacities[
                                self.absolute_non_selected_station_inds
                            ] = self.plot_characteristics["map"]["marker_unselected"][
                                "alpha"
                            ]

                            # update markersize and alphas of selected stations
                            markersizes[
                                self.absolute_selected_station_inds
                            ] = self.plot_characteristics["map"]["marker_selected"]["s"]
                            opacities[
                                self.absolute_selected_station_inds
                            ] = self.plot_characteristics["map"]["marker_selected"][
                                "alpha"
                            ]

                        # set new markersizes and alphas
                        collection.set_sizes(markersizes)
                        collection.set_alpha(opacities)

        # redraw plot
        self.figure.canvas.draw()

        # repaint(), not flush_events() - the latter is processEvents() with
        # nothing excluded, so it delivers queued user input too. Running
        # inside a settings handler, that delivered the mouse release
        # mid-handler and the button's clicked signal never fired
        self.figure.canvas.repaint()

    def update_associated_active_dashboard_plot(self, plot_type):
        """
        Function that updates a plot associated with selected stations on map

        Parameters
        ----------
        plot_type : str
            Plot type
        """

        if hasattr(self, "relative_selected_station_inds"):
            if len(self.relative_selected_station_inds) > 0:
                # get numeric position of plot type in dashboard
                plot_type_position = self.get_plot_type_position(plot_type)

                # if there are no temporal resolutions (only yearly), skip periodic plots
                if (plot_type in ["periodic", "periodic-violin"]) and (
                    self.read_instance.active_resolution == "annual"
                ):
                    msg = "It is not possible to make periodic plots using annual resolution data."
                    show_message(self.read_instance, msg)
                    self.read_instance.handle_layout_update(
                        "None", sender=plot_type_position
                    )
                    return

                # if temporal colocation is turned off or there are no models, skip some plots
                if plot_type in [
                    "scatter",
                    "taylor",
                    "fairmode-target",
                    "fairmode-statsummary",
                    "contingencytable",
                ]:
                    if (not self.read_instance.temporal_colocation) or (
                        (self.read_instance.temporal_colocation)
                        and (len(self.read_instance.data_labels) == 1)
                    ):
                        if not self.read_instance.temporal_colocation:
                            msg = f"It is not possible to make {plot_type} plots without activating the temporal colocation."
                        else:
                            msg = f"It is not possible to make {plot_type} plots without loading models."
                        show_message(self.read_instance, msg)
                        self.read_instance.handle_layout_update(
                            "None", sender=plot_type_position
                        )
                        return

                speci = self.read_instance.networkspeci.split("|")[1]
                if plot_type == "contingencytable":
                    # if we have more than one model, skip contingency table
                    if len(self.read_instance.data_labels) > 2:
                        msg = f"It is not possible to make {plot_type} plots with more than 1 model."
                        show_message(self.read_instance, msg)
                        self.read_instance.handle_layout_update(
                            "None", sender=plot_type_position
                        )
                        return
                    # if do not have correct species or resolution, cannot make contingency table
                    if speci not in [
                        "sconco3",
                        "sconcno2",
                        "pm10",
                        "pm2p5",
                        "sconcso2",
                    ]:
                        msg = f"Contingency table cannot be created for {speci}."
                        show_message(self.read_instance, msg)
                        self.read_instance.handle_layout_update(
                            "None", sender=plot_type_position
                        )
                        return
                    if self.read_instance.active_resolution != "hourly":
                        msg = "Contingency table can only be created if the resolution is hourly."
                        show_message(self.read_instance, msg)
                        self.read_instance.handle_layout_update(
                            "None", sender=plot_type_position
                        )
                        return

                # if do not have correct resolution or species, cannot make fairmode plots
                if plot_type in ["fairmode-target", "fairmode-statsummary"]:
                    if speci not in ["sconco3", "sconcno2", "pm10", "pm2p5"]:
                        msg = f"Fairmode plot cannot be created for {speci}."
                        show_message(self.read_instance, msg)
                        self.read_instance.handle_layout_update(
                            "None", sender=plot_type_position
                        )
                        return
                    if (
                        speci in ["sconco3", "sconcno2"]
                        and self.read_instance.active_resolution != "hourly"
                    ) or (
                        speci in ["pm10", "pm2p5"]
                        and (
                            self.read_instance.active_resolution
                            not in ["hourly", "daily"]
                        )
                    ):
                        msg = "Fairmode plot can only be created if the resolution is hourly (O3, NO2, PM2.5 and PM10) or daily (PM2.5 and PM10)."
                        show_message(self.read_instance, msg)
                        self.read_instance.handle_layout_update(
                            "None", sender=plot_type_position
                        )
                        return

                # clear all previously plotted artists for plot type
                self.remove_axis_elements(self.plot_axes[plot_type], plot_type)

                # get relevant axis
                ax = self.plot_axes[plot_type]

                # get options defined to configure plot
                plot_options = copy.deepcopy(self.current_plot_options[plot_type])

                # get plotting function for specific plot
                if plot_type == "statsummary":
                    func = getattr(self.plotting, "make_table")
                elif plot_type in ["fairmode-target", "fairmode-statsummary"]:
                    func = getattr(
                        self.plotting, "make_{}".format(plot_type.replace("-", "_"))
                    )
                else:
                    func = getattr(
                        self.plotting, "make_{}".format(plot_type.split("-")[0])
                    )

                # get timeseries chunking info
                if plot_type == "timeseries":
                    chunk_stat = self.timeseries_chunk_stat.currentText()
                    chunk_resolution = self.timeseries_chunk_resolution.currentText()

                # set ylabel for periodic plot
                if plot_type == "periodic" or (
                    (plot_type == "timeseries")
                    and (chunk_stat != "None")
                    and (chunk_resolution != "None")
                ):
                    # get information on periodic stat
                    if plot_type == "periodic":
                        base_zstat = self.periodic_stat.currentText()

                        if "bias" in plot_options:
                            zstat = get_z_statistic_comboboxes(base_zstat, bias=True)
                        else:
                            zstat = get_z_statistic_comboboxes(base_zstat)

                        # get zstat information
                        (
                            zstat,
                            base_zstat,
                            z_statistic_type,
                            z_statistic_sign,
                            z_statistic_period,
                        ) = get_z_statistic_info(zstat=zstat)

                    # get information on timeseries stat
                    elif plot_type == "timeseries":
                        # get zstat information
                        (
                            zstat,
                            base_zstat,
                            z_statistic_type,
                            z_statistic_sign,
                            z_statistic_period,
                        ) = get_z_statistic_info(zstat=chunk_stat)

                    # set new ylabel
                    if z_statistic_type == "basic":
                        ylabel = self.read_instance.basic_stats[base_zstat]["label"]
                        ylabel_units = self.read_instance.basic_stats[base_zstat][
                            "units"
                        ]
                    else:
                        ylabel = self.read_instance.modbias_stats[base_zstat]["label"]
                        ylabel_units = self.read_instance.modbias_stats[base_zstat][
                            "units"
                        ]
                    if ylabel_units == "[measurement_units]":
                        ylabel_units = self.read_instance.measurement_units[
                            self.read_instance.species[0]
                        ]
                    if ylabel_units != "":
                        ylabel += " [{}]".format(ylabel_units)
                    xlabel = ""

                    # if statistic type is 'modbias' and 'bias' in plot options, remove bias from plot options
                    if (z_statistic_type == "modbias") and ("bias" in plot_options):
                        bias_index = self.plot_characteristics[plot_type][
                            "plot_options"
                        ].index("bias")
                        plot_options.remove("bias")
                        self.plot_elements[plot_type]["active"] = "absolute"
                        self.read_instance.block_MPL_canvas_updates = True
                        self.periodic_options.model().item(bias_index).setCheckState(
                            QtCore.Qt.Unchecked
                        )
                        self.read_instance.block_MPL_canvas_updates = False

                # create structure to store data for Taylor diagram
                elif plot_type == "taylor":
                    # get r or r2 as correlation statistic
                    corr_stat = self.plot_characteristics[plot_type]["corr_stat"]
                    relevant_zstats = [corr_stat, "StdDev"]

                # setup xlabel / ylabel for other plot_types
                else:
                    # set new xlabel
                    if "xlabel" in self.plot_characteristics[plot_type]:
                        xlabel = self.plot_characteristics[plot_type]["xlabel"][
                            "xlabel"
                        ]
                        if "[measurement_units]" in xlabel:
                            xlabel = xlabel.replace(
                                "[measurement_units]",
                                "[{}]".format(
                                    self.read_instance.measurement_units[
                                        self.read_instance.species[0]
                                    ]
                                ),
                            )
                    else:
                        xlabel = ""

                    # set new ylabel
                    if "ylabel" in self.plot_characteristics[plot_type]:
                        ylabel = self.plot_characteristics[plot_type]["ylabel"][
                            "ylabel"
                        ]
                        if "[measurement_units]" in ylabel:
                            ylabel = ylabel.replace(
                                "[measurement_units]",
                                "[{}]".format(
                                    self.read_instance.measurement_units[
                                        self.read_instance.species[0]
                                    ]
                                ),
                            )
                    else:
                        ylabel = ""

                # call function to update plot
                # periodic plot
                if plot_type == "periodic":
                    func(
                        ax,
                        self.read_instance.networkspeci,
                        self.read_instance.data_labels,
                        self.plot_characteristics[plot_type],
                        plot_options,
                        zstat=zstat,
                    )
                # make statsummary plot
                elif plot_type == "statsummary":
                    if "bias" in plot_options:
                        relevant_zstats = self.active_statsummary_stats["modbias"]
                    else:
                        relevant_zstats = self.active_statsummary_stats["basic"]

                    func(
                        ax,
                        self.read_instance.networkspeci,
                        self.read_instance.data_labels,
                        self.plot_characteristics[plot_type],
                        plot_options,
                        zstats=relevant_zstats,
                        statsummary=True,
                    )
                # make taylor diagram
                elif plot_type == "taylor":
                    corr_stat = self.plot_characteristics["taylor"]["corr_stat"]
                    func(
                        ax,
                        self.read_instance.networkspeci,
                        self.read_instance.data_labels,
                        self.plot_characteristics[plot_type],
                        plot_options,
                        corr_stat,
                    )
                # other plots
                else:
                    func(
                        ax,
                        self.read_instance.networkspeci,
                        self.read_instance.data_labels,
                        self.plot_characteristics[plot_type],
                        plot_options,
                    )

                # reset axes limits (harmonising across subplots for periodic plots)
                if plot_type not in ["map", "taylor", "fairmode-statsummary"]:
                    if plot_type == "scatter":
                        harmonise_xy_lims_paradigm(
                            self.read_instance,
                            self,
                            ax,
                            plot_type,
                            self.plot_characteristics[plot_type],
                            plot_options,
                            relim=True,
                        )
                    else:
                        harmonise_xy_lims_paradigm(
                            self.read_instance,
                            self,
                            ax,
                            plot_type,
                            self.plot_characteristics[plot_type],
                            plot_options,
                            relim=True,
                            autoscale=True,
                        )

                # set axes labels
                if plot_type not in ["taylor", "contingencytable", "statsummary"]:
                    # set xlabel
                    set_axis_label(
                        ax, "x", xlabel, self.plot_characteristics[plot_type]
                    )
                    # set ylabel
                    set_axis_label(
                        ax, "y", ylabel, self.plot_characteristics[plot_type]
                    )

                # reset navigation toolbar stack for plot
                self.reset_ax_navigation_toolbar_stack(ax)

                # update plot options, except for plots with no options in dashboard
                if plot_type not in ["metadata", "fairmode-statsummary"]:
                    self.update_plot_options(plot_types=[plot_type])

    def get_plot_type_position(self, plot_type):
        """
        Function that returns numeric position of plot type within the dashboard

        Parameters
        ----------
        plot_type : str
            Plot type
        """

        position_vars = [
            self.read_instance.position_2,
            self.read_instance.position_3,
            self.read_instance.position_4,
            self.read_instance.position_5,
        ]
        positions = [2, 3, 4, 5]
        for position, position_var in zip(positions, position_vars):
            if plot_type == position_var:
                return position

    def update_associated_active_dashboard_plots(self):
        """
        Function that updates all plots associated with selected stations on map
        """

        # update dashboard plots
        if hasattr(self, "relative_selected_station_inds"):
            # have no selected stations, so clear all previously plotted artists from selected station plots
            # cover plotting axes also
            active_plots = self.read_instance.active_dashboard_plots
            if len(self.relative_selected_station_inds) == 0:
                for plot_type in active_plots:
                    if plot_type != "None":
                        self.remove_axis_elements(self.plot_axes[plot_type], plot_type)
                self.cover_plot_axes()

            elif len(self.relative_selected_station_inds) > 0:
                # get selected station data
                get_selected_station_data(
                    read_instance=self.read_instance,
                    canvas_instance=self,
                    networkspecies=[self.read_instance.networkspeci],
                )

                # iterate through active_dashboard_plots
                for plot_type in active_plots:
                    # update plot
                    if plot_type != "None":
                        self.update_associated_active_dashboard_plot(plot_type)

                # un-hide plotting axes
                self.top_right_canvas_cover.hide()
                self.lower_canvas_cover.hide()

            # update map plot options
            self.update_plot_options(plot_types=["map"])

    def update_model_domain_edges(self):
        """
        Function that plots grid domain edges of models in memory
        """

        # remove grid domain polygon if previously plotted
        self.remove_axis_objects(
            self.plot_axes["map"].patches, types_to_remove=[matplotlib.patches.Polygon]
        )

        # create grid edge polygons for models in memory
        grid_edge_polygons = self.plotting.make_model_domain_polygons()

        # plot grid edge polygons on map
        for grid_edge_polygon in grid_edge_polygons:
            self.plot_axes["map"].add_patch(grid_edge_polygon)

    def update_legend(self):
        """
        Function that updates legend
        """

        # create legend element handles
        legend_plot_characteristics = self.plotting.make_legend_handles(
            copy.deepcopy(self.plot_characteristics["legend"])
        )

        # plot legend
        self.legend = self.plot_axes["legend"].legend(
            **legend_plot_characteristics["plot"],
            prop=legend_plot_characteristics["prop"],
        )

        # setup element picker in legend, and clip legend text to axis bounds.
        # gid carries the real data label through, independent of the display
        # text - matplotlib doesn't preserve a gid set on the handles passed
        # into legend(), so it is redone here on the legend's own text
        for legend_label, data_label in zip(
            self.legend.texts, legend_plot_characteristics["data_labels_ordered"]
        ):
            legend_label.set_gid(data_label)
            legend_label.set_picker(True)

            # a label whose data has been clicked off is drawn in regular
            # weight rather than bold (see _toggle_legend_visibility()), and
            # this builds the legend afresh - so without putting that back,
            # anything hidden returns looking as though it were showing,
            # which a rename did every time it rebuilt the legend
            active_labels = self.plot_elements.get("data_labels_active")
            if (active_labels is not None) and (data_label not in active_labels):
                legend_label.set_fontweight("regular")

        return None

    @restores_settings_guard
    def handle_map_z_statistic_update(self):
        """
        Function which handles update of map z statistic upon interaction with map comboboxes
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_z_statistic_update"
            )
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)

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
            if selected_z_stat == "":
                selected_z_stat = self.read_instance.basic_z_stats[0]
            if hasattr(self.read_instance, "map_z"):
                if self.read_instance.map_z in self.read_instance.basic_z_stats:
                    selected_z_stat = self.read_instance.map_z

            # if z1 is initialised or labels have changed (by setting data label through a conf
            # for instance, we will reset z1)
            if (selected_z1_array == "") or (
                selected_z1_array not in self.read_instance.data_labels
            ):
                selected_z1_array = copy.deepcopy(
                    self.read_instance.observations_data_label
                )

            # update z statistic field to all basic stats if colocation not-active OR z2
            # array not selected, else select basic+bias stats
            if (
                (not self.read_instance.temporal_colocation)
                or (selected_z2_array == "")
                or (len(self.read_instance.data_labels) == 1)
            ):
                z_stat_items = copy.deepcopy(self.read_instance.basic_z_stats)
            else:
                z_stat_items = copy.deepcopy(self.read_instance.basic_and_bias_z_stats)

            # remove nonsensical available map stats
            nonsensical_map_stats = ["NStations", "NUniqueStations", "MDA8"]
            for nonsensical_map_stat in nonsensical_map_stats:
                if nonsensical_map_stat in z_stat_items:
                    z_stat_items = z_stat_items[z_stat_items != nonsensical_map_stat]

            # remove selected z1/z2 items from opposite z2/z1 comboboxes (if have value
            # selected, i.e. z2 array not empty string)
            if selected_z2_array != "":
                z1_items = np.delete(
                    self.read_instance.z1_arrays,
                    np.where(self.read_instance.z1_arrays == selected_z2_array)[0],
                )
            else:
                z1_items = self.read_instance.z1_arrays
            z2_items = np.delete(
                self.read_instance.z2_arrays,
                np.where(self.read_instance.z2_arrays == selected_z1_array)[0],
            )

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

            # a manually-typed colourbar limit rarely makes sense carried
            # over to a different statistic (a concentration range typed for
            # "Mean" would mis-scale "Bias"), so clear it back to auto
            self.map_cb_min.clear()
            self.map_cb_max.clear()
            self.read_instance.map_vmin_override = None
            self.read_instance.map_vmax_override = None
            # the colourmap suited to one statistic rarely suits the next, so
            # fall back to the statistic's own - see sync_map_colourmap()
            self.read_instance.map_colourmap_override = None
            self.sync_map_colourmap()
            # with that choice dropped the basemap may match its preset again
            self.sync_map_colour_preset()

            # update plotted map z statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_map_z_statistic()

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_z_statistic_update"
            )

        return None

    @restores_settings_guard
    def handle_map_colourmap_update(self):
        """
        Function which handles update of the map colourbar's colourmap upon
        interaction with the map colourmap combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_colourmap_update"
            )
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            self.read_instance.map_colourmap_override = self.map_colourmap.currentText()
            self.sync_map_colour_preset()

            # update plotted map z statistic (re-generates the colourbar too)
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_map_z_statistic()

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_colourmap_update"
            )

        return None

    @restores_settings_guard
    def handle_map_colourmap_scale_update(self):
        """
        Function which handles update of the map colourbar's discrete/
        continuous scale upon interaction with the map colourmap scale
        combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_colourmap_scale_update"
            )
            # force the cursor change to paint before the work below restores
            # it, as the set+unset can otherwise happen within one event loop
            # pass and never show. ExcludeUserInputEvents so this pump cannot
            # process a fresh click mid-handler
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            is_discrete = self.map_colourmap_scale.currentText() == "Discrete"
            self.read_instance.map_discrete_override = is_discrete
            self.sync_map_n_sections_visibility()

            # update plotted map z statistic (re-generates the colourbar too)
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_map_z_statistic()

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function, "handle_map_colourmap_scale_update"
            )

        return None

    @restores_settings_guard
    def handle_map_n_sections_update(self):
        """
        Function which handles update of the map colourbar's number of
        discrete sections upon interaction with the map sections combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_n_sections_update"
            )
            # see handle_map_colourmap_scale_update() for why this is here
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            self.read_instance.map_n_discrete_override = self._map_count_field_value(
                self.map_n_sections,
                self.read_instance.map_n_discrete_override,
                "chunks",
            )

            # update plotted map z statistic (re-generates the colourbar too)
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_map_z_statistic()

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function, "handle_map_n_sections_update"
            )

        return None

    @restores_settings_guard
    def handle_map_n_ticks_update(self):
        """
        Function which handles update of the number of tick labels shown
        on the map colourbar upon interaction with the map labels combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_n_ticks_update"
            )
            # see handle_map_colourmap_scale_update() for why this is here
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            self.read_instance.map_n_ticks_override = self._map_count_field_value(
                self.map_n_ticks, self.read_instance.map_n_ticks_override, "labels"
            )

            # update plotted map z statistic (re-generates the colourbar too)
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_map_z_statistic()

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function, "handle_map_n_ticks_update"
            )

        return None

    @restores_settings_guard
    def handle_map_colour_preset_update(self):
        """
        Function which applies a ready-made land/ocean/colourmap
        combination (see settings/colourmaps.yaml) upon interaction with the map
        colour preset combobox - one selection instead of setting the
        three individually.

        Selecting "Custom" does nothing: it isn't a preset, it's what the
        box falls back to showing once any of the three has been changed
        on its own, so the selector never claims a preset that no longer
        matches what's on screen (see sync_map_colour_preset()).
        """

        if not self.read_instance.block_config_bar_handling_updates:
            preset_name = self.map_colour_preset.currentText()
            preset = get_colour_presets().get(preset_name)
            if preset is None:
                return None

            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_colour_preset_update"
            )
            # see handle_map_colourmap_scale_update() for why this is here
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            # held across all three controls, so their own handlers no-op
            # and the map is redrawn once at the end rather than per control
            self.read_instance.block_config_bar_handling_updates = True

            map_template = self.plot_characteristics_templates["map"]
            self.set_map_colour_preset(preset_name)
            map_template["land_polygon"]["facecolor"] = preset["land"]
            map_template["ocean_polygon"]["facecolor"] = preset["ocean"]
            populate_colour_combobox(
                self.map_land_colour, LAND_COLOUR_OPTIONS, current=preset["land"]
            )
            populate_colour_combobox(
                self.map_ocean_colour, OCEAN_COLOUR_OPTIONS, current=preset["ocean"]
            )
            # a preset gives a colourmap per statistic type rather than one
            # colourmap, so clear any explicit choice and let the statistic's
            # own role pick from the preset
            self.read_instance.map_colourmap_override = None
            self.sync_map_colourmap()

            if not self.read_instance.block_MPL_canvas_updates:
                # the colourmap change goes through the z statistic
                # redraw, which rebuilds the colourbar; the land/ocean
                # colours are cartopy features and need their own refresh
                self.refresh_map_features()
                self.refresh_map_gridlines()
                self.update_map_z_statistic()

            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_colour_preset_update"
            )

        return None

    def set_map_colour_preset(self, preset_name):
        """
        Function which records the active map colour preset.

        Written to the map's plot characteristics as well as the template it
        was copied from: the two are separate dicts (see Plotting.make_plot(),
        which deep copies the template once), and the colourbar resolves its
        colourmap from the copy while the basemap features are drawn from the
        template.

        Parameters
        ----------
        preset_name : str
            Name of the preset, or an empty string for a custom combination
        """

        self.plot_characteristics_templates["map"]["colour_preset"] = preset_name
        if "map" in self.plot_characteristics:
            self.plot_characteristics["map"]["colour_preset"] = preset_name

        return None

    def sync_map_colourmap(self):
        """
        Function which points the colourmap selector at the colourmap the
        current statistic actually resolves to, so the box always reports what
        is on screen.

        The colourmap a statistic gets depends on what it measures - an error
        that is never negative reads differently to a signed bias - so it is
        resolved per statistic rather than held fixed. Only an explicit choice
        in the selector overrides it, and that is cleared whenever the
        statistic changes, same as the colourbar limits.
        """

        resolved = getattr(self.read_instance, "map_colourmap_override", None)
        if not resolved and not self.map_z_stat.currentText():
            # nothing plotted yet, so there is no statistic to resolve against
            return None

        if not resolved:
            # compose the statistic exactly as update_map_z_statistic() does:
            # with a second dataset selected the map shows the bias form, which
            # needs a different colourmap to the absolute one of the same name
            zstat = get_z_statistic_comboboxes(
                self.map_z_stat.currentText(),
                bias=self.map_z2.currentText() != "",
            )
            resolved = resolve_colourmap(
                self.read_instance,
                zstat,
                self.plot_characteristics_templates["map"],
                self.read_instance.networkspeci.split("|")[-1],
            )

        if resolved:
            # showing the resolved colourmap must not read as choosing it, so
            # the settings guard is held while it is set - not blockSignals(),
            # which on this editable combobox also stops Qt updating the line
            # edit the closed box actually displays, leaving the old name on
            # screen (see ComboBox in dashboard_elements.py). Restored to what
            # it was rather than cleared, as callers already hold it
            previous = self.read_instance.block_config_bar_handling_updates
            self.read_instance.block_config_bar_handling_updates = True
            select_colourmap(self.map_colourmap, resolved)
            self.read_instance.block_config_bar_handling_updates = previous

        return None

    def sync_map_colour_preset(self):
        """
        Function which points the colour preset selector at whichever preset
        the current land and ocean colours match, or at "Custom" when they
        match none of them.

        Called after either colour changes individually, so the box stops
        naming a preset the moment the basemap stops being that preset. The
        colourmap takes no part in the match: a preset now carries one
        colourmap per statistic type rather than a single one, and a choice
        made in the colourmap selector only lasts until the statistic changes.
        """

        # compare colours as resolved RGB, not as the strings they happen to
        # be written as - the config stores the default land colour as "0.85"
        # where the preset spells it "#D9D9D9", so a string comparison never
        # matched and the selector opened on "Custom"
        def same_colour(first, second):
            try:
                return matplotlib.colors.to_hex(
                    first
                ).lower() == matplotlib.colors.to_hex(second).lower()
            except ValueError:
                return first == second

        # the resolved colours, not the raw ones - land/ocean are left empty in
        # the config when they come from the preset, and comparing those empty
        # values against every preset would always fall through to "Custom"
        map_template = self.plot_characteristics_templates["map"]
        current_land, current_ocean = get_map_colours(map_template)
        # a preset is its basemap and its colourmaps together, so a colourmap
        # chosen by hand takes the selector off the preset just as a changed
        # land or ocean colour does. A preset carries one colourmap per
        # statistic type, so the comparison is against the one for the
        # statistic on screen. With no choice made there is nothing to compare
        # and the basemap alone decides
        presets = get_colour_presets()
        chosen_colourmap = getattr(self.read_instance, "map_colourmap_override", None)
        role = None
        if chosen_colourmap:
            role = get_colourmap_role(
                get_z_statistic_comboboxes(
                    self.map_z_stat.currentText(),
                    bias=self.map_z2.currentText() != "",
                )
            )
        matched = next(
            (
                name
                for name, preset in presets.items()
                if same_colour(preset["land"], current_land)
                and same_colour(preset["ocean"], current_ocean)
                and (role is None or preset.get(role) == chosen_colourmap)
            ),
            COLOUR_PRESET_CUSTOM,
        )
        if matched == COLOUR_PRESET_CUSTOM:
            # the basemap no longer matches a preset exactly, so write the
            # resolved colours back - nothing is left to be filled in later.
            # colour_preset is deliberately left naming the preset it came
            # from: it still says where the colourmaps come from, and a
            # basemap tweaked from a dark preset needs to keep that preset's
            # colourmaps rather than fall back to the light ones in defaults
            map_template["land_polygon"]["facecolor"] = current_land
            map_template["ocean_polygon"]["facecolor"] = current_ocean
        else:
            self.set_map_colour_preset(matched)
        # setCurrentIndex(), not setCurrentText() - see the comment on
        # map_colourmap_scale in generate_interactive_elements()
        names = list(presets) + [COLOUR_PRESET_CUSTOM]
        self.map_colour_preset.setCurrentIndex(names.index(matched))

        return None

    def handle_map_panel_reset(self):
        """
        Function which returns every control in the map settings menu's
        "Map" sub-panel (projection, land/ocean colour, map
        resolution, country borders, gridlines) to the state the dashboard
        started up in, upon clicking the reset control beside its title.

        block_config_bar_handling_updates is held for the whole run rather
        than left to the individual handlers, so setting six controls
        redraws the map once at the end instead of once per control.
        """

        # the projection is captured before the defaults are applied, so
        # the redraw below can tell whether it actually changed
        projection_before = self.map_projection.currentText()

        def redraw():
            if self.map_projection.currentText() != projection_before:
                # a different projection needs the whole axes rebuilding,
                # which redraws the features and gridlines with it
                self.rebuild_map_axes(self.map_projection.currentText())
            else:
                self.refresh_map_features()
                self.refresh_map_gridlines()
            # the colourmap is reset here too, so the colourbar and the
            # points' colours need rebuilding either way
            self.update_map_z_statistic()

        self._reset_map_subpanel(
            "handle_map_panel_reset", self._apply_map_panel_defaults, redraw
        )

        return None

    def handle_map_colourbar_panel_reset(self):
        """
        Function which returns every control in the map settings menu's
        "Colourbar" sub-panel (limits, colourmap, scale, number of
        labels/chunks) to the state the dashboard started up in, upon
        clicking the reset control beside its title. See
        handle_map_panel_reset().
        """

        self._reset_map_subpanel(
            "handle_map_colourbar_panel_reset",
            self._apply_map_colourbar_defaults,
            self.update_map_z_statistic,
        )

        return None

    def handle_map_points_panel_reset(self):
        """
        Function which returns every control in the map settings menu's
        "Points" sub-panel (automatic sizing, and the unselected/selected
        marker size and opacity sliders) to the state the dashboard
        started up in, upon clicking the reset control beside its title.
        See handle_map_panel_reset().
        """

        def redraw():
            if self.read_instance.map_auto_marker_sizing:
                # hands the sliders back to the automatic computation, so
                # they end up mirroring it rather than the config numbers
                self.apply_automatic_marker_style()
            self.update_map_station_selection()

        self._reset_map_subpanel(
            "handle_map_points_panel_reset", self._apply_map_points_defaults, redraw
        )

        return None

    @restores_settings_guard
    def _reset_map_subpanel(self, name, apply_defaults, redraw):
        """
        Shared implementation behind the three sub-panel reset controls:
        restore that panel's startup state, then redraw once.

        Parameters
        ----------
        name : str
            Handler name, for the busy-cursor ownership bookkeeping.
        apply_defaults : callable
            Sets this panel's controls (and the read_instance overrides
            behind them) back to their captured startup values.
        redraw : callable
            The single redraw to run once everything is set.
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # see handle_map_cb_limits_reset() - apply what the boxes hold
            # before undoing it
            self.apply_map_cb_fields()

            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, name
            )
            # see handle_map_colourmap_scale_update() for why this is here
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            # held across every control this sets, so their own change
            # handlers all no-op and the redraw happens once, below
            self.read_instance.block_config_bar_handling_updates = True

            apply_defaults()

            if not self.read_instance.block_MPL_canvas_updates:
                redraw()

            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(self.read_instance.cursor_function, name)

        return None

    def _apply_map_panel_defaults(self):
        """Restore the "Map" sub-panel's startup control state."""

        defaults = self.map_panel_startup_defaults
        self.map_projection.setCurrentIndex(defaults["projection"])
        self.map_land_colour.setCurrentIndex(defaults["land_colour"])
        self.map_ocean_colour.setCurrentIndex(defaults["ocean_colour"])
        self.map_resolution.setCurrentIndex(defaults["resolution"])
        self.map_borders.setChecked(defaults["borders"])
        self.map_gridlines.setChecked(defaults["gridlines"])

        # the controls above are the menu's view of these template values;
        # the map itself reads them from here, so both have to be put back
        map_template = self.plot_characteristics_templates["map"]
        map_template["projection"] = self.map_projection.currentText()
        map_template["land_polygon"]["facecolor"] = LAND_COLOUR_OPTIONS[
            self.map_land_colour.currentText()
        ]
        map_template["ocean_polygon"]["facecolor"] = OCEAN_COLOUR_OPTIONS[
            self.map_ocean_colour.currentText()
        ]
        map_template["map_resolution"] = self.map_resolution.currentText()
        map_template["borders"]["visible"] = defaults["borders"]
        self.plot_characteristics["map"]["gridlines"]["visible"] = defaults["gridlines"]
        # resolved once the basemap is back, so the colourmap follows the
        # preset the restored land/ocean colours belong to
        self.read_instance.map_colourmap_override = None
        self.sync_map_colour_preset()
        self.sync_map_colourmap()

        return None

    def _apply_map_colourbar_defaults(self):
        """Restore the "Colourbar" sub-panel's startup control state."""

        defaults = self.map_colourbar_startup_defaults
        self.map_colourmap_scale.setCurrentIndex(defaults["colourmap_scale"])
        self.map_n_ticks.setText(defaults["n_ticks"])
        self.map_n_sections.setText(defaults["n_sections"])
        # limits go back to automatic - they have no startup value of
        # their own, being filled in from whatever each redraw resolves
        self.map_cb_min.clear()
        self.map_cb_max.clear()

        # the colourmap goes back to following the statistic, the same way the
        # limits go back to automatic - a colourmap pinned here would outlast
        # the reset and keep overriding every statistic that followed
        self.read_instance.map_colourmap_override = None
        self.read_instance.map_discrete_override = (
            self.map_colourmap_scale.currentText() == "Discrete"
        )
        self.read_instance.map_n_discrete_override = (
            int(defaults["n_sections"]) if defaults["n_sections"] else None
        )
        self.read_instance.map_n_ticks_override = (
            int(defaults["n_ticks"]) if defaults["n_ticks"] else None
        )
        self.read_instance.map_vmin_override = None
        self.read_instance.map_vmax_override = None
        self.sync_map_n_sections_visibility()
        self.sync_map_colourmap()
        # the basemap has not changed, but the colourmap has, so the
        # combination may no longer be the preset the selector is naming
        self.sync_map_colour_preset()

        return None

    def _apply_map_points_defaults(self):
        """
        Restore the "Points" sub-panel's startup control state.

        While automatic sizing is on, the four sliders are not user state
        at all - apply_automatic_marker_style() rewrites them on every
        redraw so they mirror the size and opacity it computed for the
        current zoom and selection. Forcing them back to the captured
        startup numbers therefore moved every slider even when nothing had
        been changed, which is what made reset look like it was undoing
        edits that were never made. So each control is only written when
        it actually differs, and when automatic sizing is left on the
        sliders are handed straight back to it to re-mirror, rather than
        being pinned to config values it is about to overwrite anyway.
        """

        defaults = self.map_points_startup_defaults
        automatic = defaults["auto_sizing"]

        if self.map_auto_sizing.isChecked() != automatic:
            self.map_auto_sizing.setChecked(automatic)
        self.read_instance.map_auto_marker_sizing = automatic

        if not automatic:
            # only meaningful as user state when the sliders are the ones
            # actually driving the look
            for slider, value in (
                (self.map_markersize_unsel_sl, defaults["markersize_unsel"]),
                (self.map_opacity_unsel_sl, defaults["opacity_unsel"]),
                (self.map_markersize_sel_sl, defaults["markersize_sel"]),
                (self.map_opacity_sel_sl, defaults["opacity_sel"]),
            ):
                if slider.value() != value:
                    slider.setValue(value)

            map_characteristics = self.plot_characteristics["map"]
            map_characteristics["marker_unselected"]["s"] = defaults["markersize_unsel"]
            map_characteristics["marker_unselected"]["alpha"] = (
                defaults["opacity_unsel"] / 10
            )
            map_characteristics["marker_selected"]["s"] = defaults["markersize_sel"]
            map_characteristics["marker_selected"]["alpha"] = (
                defaults["opacity_sel"] / 10
            )

        self.sync_map_sizing_sliders_enabled()

        return None

    def _map_count_field_value(self, lineedit, current, what):
        """
        Read one of the colourbar's count fields, holding it to a minimum
        of one.

        A colourbar with no chunks or no labels isn't a lesser version of
        one, it's a broken one, so an empty or zero entry is refused and
        the field put back to the value it had rather than being taken as
        "automatic" - which is what an empty field used to mean here, and
        made it impossible to tell a deliberate reset from a typo.

        Parameters
        ----------
        lineedit : MenuLineEdit
            The field being read.
        current : int or None
            The value currently in force, restored if the entry is refused.
        what : str
            Name of the quantity, for the message shown when refusing.

        Returns
        -------
        int or None
            The value to apply.
        """

        text = lineedit.text().strip()
        value = int(text) if text.isdigit() else 0
        if value < 1:
            lineedit.setText(str(current) if current else "")
            msg = f"The number of colourbar {what} must be at least 1."
            show_message(self.read_instance, msg)
            return current
        return value

    @restores_settings_guard
    def handle_map_cb_limits_reset(self):
        """
        Function which clears both map colourbar limit fields back to
        automatic, upon clicking the small reset control beside them.

        Equivalent to emptying both fields by hand - the override is
        dropped and each limit falls back to whatever would otherwise be
        resolved (the plotted data's range, or a per-statistic/
        configuration-file default) - but as one click rather than two
        edits, since going back to automatic is much the most common
        thing to want after trying a manual limit. The fields are
        repopulated with the newly resolved numbers by
        update_map_z_statistic(), so they don't stay blank.
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # whatever the boxes currently hold is applied first, so a
            # reset always undoes a known state rather than racing the
            # edit that prompted it - see apply_map_cb_fields()
            self.apply_map_cb_fields()

            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_cb_limits_reset"
            )
            # see handle_map_colourmap_scale_update() for why this is here
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            self.read_instance.map_vmin_override = None
            self.read_instance.map_vmax_override = None
            self.map_cb_min.clear()
            self.map_cb_max.clear()

            # update plotted map z statistic (re-generates the colourbar too)
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_map_z_statistic()

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function, "handle_map_cb_limits_reset"
            )

        return None

    @restores_settings_guard
    def handle_map_cb_limits_update(self):
        """
        Function which handles update of the map colourbar's min/max limits
        upon editing the colourbar limit fields. Each field is kept
        populated with the actual resolved limit (see
        update_map_z_statistic()), so editing one starts from the real
        current number; clearing a field back to blank falls back to
        whichever limit would otherwise be resolved (the plotted data's
        range, or any per-statistic/configuration-file default) - see
        generate_colourbar_detail() in statistics.py.
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_cb_limits_update"
            )
            # see handle_map_colourmap_scale_update() for why this is here
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            for lineedit, override_attr in (
                (self.map_cb_min, "map_vmin_override"),
                (self.map_cb_max, "map_vmax_override"),
            ):
                text = lineedit.text().strip()
                if text == "":
                    setattr(self.read_instance, override_attr, None)
                    continue
                try:
                    value = float(text)
                except ValueError:
                    value = float("nan")
                if not np.isfinite(value):
                    setattr(self.read_instance, override_attr, None)
                    lineedit.clear()
                    msg = (
                        "Colourbar limits must be numbers. "
                        f"'{text}' will be set to automatic."
                    )
                    show_message(self.read_instance, msg)
                    continue
                setattr(self.read_instance, override_attr, value)

            # limits outside the plotted data range are allowed, as showing a
            # wider range than the data occupies is legitimate. A minimum
            # above the maximum is not - matplotlib raises out of the colour
            # normalisation - so only that case is refused here
            vmin = getattr(self.read_instance, "map_vmin_override", None)
            vmax = getattr(self.read_instance, "map_vmax_override", None)
            if (vmin is not None) and (vmax is not None) and (vmin > vmax):
                self.read_instance.map_vmin_override = None
                self.read_instance.map_vmax_override = None
                self.map_cb_min.clear()
                self.map_cb_max.clear()
                msg = (
                    "The colourbar minimum cannot be greater than the "
                    "maximum. Both limits will be set to automatic."
                )
                show_message(self.read_instance, msg)

            # update plotted map z statistic (re-generates the colourbar too)
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_map_z_statistic()

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function, "handle_map_cb_limits_update"
            )

        return None

    @restores_settings_guard
    def handle_map_projection_update(self):
        """
        Function which handles update of the map's projection upon
        interaction with the map projection combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_projection_update"
            )
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            if not self.read_instance.block_MPL_canvas_updates:
                self.rebuild_map_axes(self.map_projection.currentText())

            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_projection_update"
            )

        return None

    @restores_settings_guard
    def handle_map_land_colour_update(self):
        """
        Function which handles update of the map's land colour upon
        interaction with the map land colour combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_land_colour_update"
            )
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            selected_label = self.map_land_colour.currentText()
            self.plot_characteristics_templates["map"]["land_polygon"][
                "facecolor"
            ] = LAND_COLOUR_OPTIONS[selected_label]
            self.sync_map_colour_preset()
            if not self.read_instance.block_MPL_canvas_updates:
                self.refresh_map_features()
                # borders and gridlines both take their colour from the
                # basemap's brightness (see map_feature_ink()), so a
                # land/ocean/resolution change has to redraw the
                # gridlines as well or they keep the old contrast
                self.refresh_map_gridlines()

            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_land_colour_update"
            )

        return None

    @restores_settings_guard
    def handle_map_ocean_colour_update(self):
        """
        Function which handles update of the map's ocean colour upon
        interaction with the map ocean colour combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_ocean_colour_update"
            )
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            selected_label = self.map_ocean_colour.currentText()
            self.plot_characteristics_templates["map"]["ocean_polygon"][
                "facecolor"
            ] = OCEAN_COLOUR_OPTIONS[selected_label]
            self.sync_map_colour_preset()
            if not self.read_instance.block_MPL_canvas_updates:
                self.refresh_map_features()
                # borders and gridlines both take their colour from the
                # basemap's brightness (see map_feature_ink()), so a
                # land/ocean/resolution change has to redraw the
                # gridlines as well or they keep the old contrast
                self.refresh_map_gridlines()

            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_ocean_colour_update"
            )

        return None

    @restores_settings_guard
    def handle_map_borders_update(self):
        """
        Function which handles update of country border visibility upon
        interaction with the map borders checkbox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_borders_update"
            )
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            self.plot_characteristics_templates["map"]["borders"][
                "visible"
            ] = self.map_borders.isChecked()
            if not self.read_instance.block_MPL_canvas_updates:
                self.refresh_map_features()
                # borders and gridlines both take their colour from the
                # basemap's brightness (see map_feature_ink()), so a
                # land/ocean/resolution change has to redraw the
                # gridlines as well or they keep the old contrast
                self.refresh_map_gridlines()

            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_borders_update"
            )

        return None

    @restores_settings_guard
    def handle_map_gridlines_update(self):
        """
        Function which handles update of gridline visibility upon
        interaction with the map gridlines checkbox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_gridlines_update"
            )
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            self.plot_characteristics["map"]["gridlines"][
                "visible"
            ] = self.map_gridlines.isChecked()
            if not self.read_instance.block_MPL_canvas_updates:
                self.refresh_map_gridlines()

            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_gridlines_update"
            )

        return None

    @restores_settings_guard
    def handle_map_resolution_update(self):
        """
        Function which handles update of the map (coastline/land polygon) resolution
        upon interaction with the map resolution combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_resolution_update"
            )
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            self.plot_characteristics_templates["map"][
                "map_resolution"
            ] = self.map_resolution.currentText()
            if not self.read_instance.block_MPL_canvas_updates:
                self.refresh_map_features()
                # borders and gridlines both take their colour from the
                # basemap's brightness (see map_feature_ink()), so a
                # land/ocean/resolution change has to redraw the
                # gridlines as well or they keep the old contrast
                self.refresh_map_gridlines()

            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_resolution_update"
            )

        return None

    def layout_map_nav_buttons(self):
        """
        Function which sizes each of the map settings menu's three
        sub-panel nav buttons ("Map", "Points", "Colourbar") to its own
        text and centres them as a row across the panel.

        Sized to the text rather than stretched to the full panel width,
        so they read as three small buttons side by side instead of three
        stacked bars. Measured from the widget's own font rather than
        hard-coded in canvas_menus.yaml, since the same label is a
        different number of pixels wide on each platform's stylesheet -
        a fixed width would either clip the text somewhere or leave one
        button visibly padded.

        All three go on one line when they fit, spread so that the outer
        two sit against the panel's edges - level with the controls above
        them, which span the same width. When they don't fit (a wider font
        than the stylesheet's own, say), "Colourbar" - the widest, and the
        odd one out of the pair that naturally belong together - drops to a
        second line centred beneath the other two, and everything below the
        row shifts down to make space for it.

        The button text is fixed, never gaining an "open" marker: the row
        is sized to its labels, so changing them would re-measure and
        re-centre every button underneath the pointer on each click, and
        could tip the row between one line and two as panels are opened
        and closed.
        """

        settings_button = self.map_menu.buttons["settings_button"]
        panel_left = settings_button.x() - 220
        panel_width = 230
        gap = 6
        row_y = settings_button.y() + 210
        row_height = 20

        buttons = [
            self.map_panel_button,
            self.map_points_panel_button,
            self.map_colourbar_panel_button,
        ]
        widths = []
        for button in buttons:
            # the text's own width plus room for the button's border and
            # a little breathing space either side
            text_width = button.fontMetrics().boundingRect(button.text()).width()
            width = text_width + 16
            widths.append(width)
            button.setFixedWidth(width)
            button.resize(width, row_height)

        def place(row_buttons, row_widths, y):
            total = sum(row_widths) + gap * (len(row_widths) - 1)
            x = panel_left + max(0, (panel_width - total) // 2)
            for button, width in zip(row_buttons, row_widths):
                button.move(x, y)
                x += width + gap

        def spread(row_buttons, row_widths, y):
            """Sit the outer buttons against the panel's edges, gaps even."""

            step = (panel_width - sum(row_widths)) // (len(row_widths) - 1)
            x = panel_left
            for button_ii, (button, width) in enumerate(zip(row_buttons, row_widths)):
                # the last is placed against the right edge rather than at
                # whatever the gaps have added up to, so rounding cannot
                # leave it a pixel or two short of the controls above
                if button_ii == len(row_buttons) - 1:
                    x = panel_left + panel_width - width
                button.move(x, y)
                x += width + step

        wrapped = sum(widths) + gap * 2 > panel_width
        if not wrapped:
            spread(buttons, widths, row_y)
        else:
            place(buttons[:2], widths[:2], row_y)
            place(buttons[2:], widths[2:], row_y + row_height + 5)

        # the row is the last thing in the panel, so a wrapped second line
        # only needs the panel itself to grow to keep it inside
        container = self.map_menu.containers["container"]
        container.resize(container.width(), 220 + (row_height + 5 if wrapped else 0))

        return None

    def handle_map_subpanel_toggle(self, name):
        """
        Function which shows or hides one of the map settings menu's
        sub-panels ("Map", "Points" or "Colourbar") upon clicking its nav
        button. Pure UI show/hide - no data changes, so unlike the other
        map handlers this doesn't touch block_config_bar_handling_updates
        or the busy cursor (matches interactive_elements_button_func,
        which the main Settings toggle itself uses for the same reason).

        A sub-panel's widgets are built as siblings of (not nested
        inside) the map settings menu's own widgets, same as every other
        settings menu control, so - unlike a real child widget, which
        would inherit its parent's stacking automatically - each one needs
        raising above the rest of the canvas explicitly on show, or it's
        technically visible but painted behind other menu content and
        never actually seen.

        All three panels share the same on-screen spot (see
        canvas_menus.yaml), so opening one closes whichever other was
        open.

        Parameters
        ----------
        name : str
            Which sub-panel to toggle - a key of self.map_subpanels.
        """

        expand = self.map_open_subpanel != name

        # only one panel can occupy the shared spot at a time
        self.reset_map_subpanels()

        if expand:
            _button, elements, _label = self.map_subpanels[name]
            for element in elements:
                element.show()
                element.raise_()

            if name == "colourbar":
                # the loop above just unconditionally showed the chunks
                # field along with everything else - re-apply the "only
                # meaningful for a discrete scale" rule on top of that,
                # rather than special-casing it out of the element list
                self.sync_map_n_sections_visibility()
            elif name == "points":
                # re-apply the enabled/disabled (manual vs automatic)
                # state of the sliders - element.show() above only
                # affects visibility, not whether they're interactive
                self.sync_map_sizing_sliders_enabled()

            self.map_open_subpanel = name

        return None

    def handle_map_panel_toggle(self):
        """
        Function which toggles the map settings menu's "Map" sub-panel
        (projection, land/ocean colour, map resolution, country
        borders, gridlines) - see handle_map_subpanel_toggle().
        """

        self.handle_map_subpanel_toggle("map")

        return None

    def handle_map_points_panel_toggle(self):
        """
        Function which toggles the map settings menu's "Points" sub-panel
        (automatic sizing, plus the unselected/selected marker size and
        opacity sliders) - see handle_map_subpanel_toggle().
        """

        self.handle_map_subpanel_toggle("points")

        return None

    def handle_map_colourbar_panel_toggle(self):
        """
        Function which toggles the map settings menu's "Colourbar"
        sub-panel (limits, colourmap, scale, number of labels/chunks) -
        see handle_map_subpanel_toggle().
        """

        self.handle_map_subpanel_toggle("colourbar")

        return None

    def sync_map_n_sections_visibility(self):
        """
        Function which shows or hides the map colourbar's "№ Chunks" field
        (and its label) depending on whether the colourmap scale is
        currently discrete - it has no effect for a continuous scale, so
        it's hidden entirely rather than left visible but inert. Called
        both when the scale itself changes and when the Colourbar panel is
        (re-)expanded, since showing the whole panel would otherwise
        unconditionally reveal it regardless of scale.
        """

        is_discrete = self.map_colourmap_scale.currentText() == "Discrete"
        self.map_n_sections_label.setVisible(is_discrete)
        self.map_n_sections.setVisible(is_discrete)

        return None

    def apply_map_cb_fields(self):
        """
        Function which reads the four colourbar fields and applies exactly
        what they currently contain, with no validation messages and no
        redraw.

        Called by the reset controls before they reset, so the values in
        the boxes are always what gets undone - whether or not the fields
        were ever "committed" by Enter, focus or a click. The reset
        handlers no longer depend on any of that machinery having run.

        Deliberately silent: a warning box opened here would be a modal
        dialog raised in the middle of handling the reset click, which
        stops the click reaching the button at all. Anything invalid is
        simply left as automatic, and the reset that follows replaces it a
        moment later anyway.
        """

        for lineedit, override_attr in (
            (self.map_cb_min, "map_vmin_override"),
            (self.map_cb_max, "map_vmax_override"),
        ):
            text = lineedit.text().strip()
            try:
                value = float(text) if text else None
            except ValueError:
                value = None
            if (value is not None) and (not np.isfinite(value)):
                value = None
            setattr(self.read_instance, override_attr, value)
            lineedit.mark_committed()

        vmin = self.read_instance.map_vmin_override
        vmax = self.read_instance.map_vmax_override
        if (vmin is not None) and (vmax is not None) and (vmin > vmax):
            # matplotlib refuses to normalise this - see
            # handle_map_cb_limits_update()
            self.read_instance.map_vmin_override = None
            self.read_instance.map_vmax_override = None

        for lineedit, override_attr in (
            (self.map_n_ticks, "map_n_ticks_override"),
            (self.map_n_sections, "map_n_discrete_override"),
        ):
            text = lineedit.text().strip()
            if text.isdigit() and int(text) >= 1:
                setattr(self.read_instance, override_attr, int(text))
            lineedit.mark_committed()

        return None

    def commit_map_pending_edits(self):
        """
        Function which applies any map settings field that has been typed
        into but not yet confirmed - so a value takes effect whether it
        was finished with Enter, by clicking away, by pressing another
        control, or by closing the panel or the whole settings menu.

        Each field's handler is invoked directly here rather than by
        emitting the widget's own "committed" signal. Routing it through
        the signal meant the value's application depended on Qt delivering
        that signal, on SettingsMenu.connect()'s dispatch gate, and on the
        ordering between the two - and a value typed and then abandoned by
        clicking straight onto a reset repeatedly failed to take effect
        somewhere along that chain. Calling the handler is the same work
        with none of the indirection, and cannot be delivered late or out
        of order.
        """

        for lineedit, handler in (
            (self.map_cb_min, self.handle_map_cb_limits_update),
            (self.map_cb_max, self.handle_map_cb_limits_update),
            (self.map_n_ticks, self.handle_map_n_ticks_update),
            (self.map_n_sections, self.handle_map_n_sections_update),
        ):
            if not lineedit.has_pending_edit():
                continue
            # marked before the handler runs, so a handler that writes
            # back to the field (the limits are repopulated with whatever
            # they resolve to) can't leave it looking pending again
            lineedit.mark_committed()
            handler()

        return None

    def reset_map_subpanels(self):
        """
        Function which hides every one of the map settings menu's
        sub-panels - connected as an extra listener on the map menu's own
        settings button (see generate_interactive_elements()), alongside
        interactive_elements_button_func, so they always start hidden
        whenever Settings is opened or closed, rather than staying open
        (or left visible with the rest of the menu hidden) from however
        they were previously left.

        Unconditional, not a check-then-toggle: self.map_open_subpanel
        only tracks which nav button has been clicked, but
        interactive_elements_button_func (the other listener on the same
        settings_button click, which always runs first) shows every
        element in self.map_elements wholesale - including every
        sub-panel's, which need to be in that list so dashboard.py's
        layout handling keeps repositioning them (see the comment above
        self.map_panel_elements). So the widgets can already be visible
        here even while the flag still says nothing is open, and a
        guard on that flag would wrongly skip re-hiding them.
        """

        # a value typed into one of the colourbar fields and left there -
        # no Enter, no click away - is still a real edit the moment the
        # panel goes away, so commit it before hiding rather than
        # discarding it silently
        self.commit_map_pending_edits()

        for _name, (_button, elements, _label) in self.map_subpanels.items():
            for element in elements:
                element.hide()
        self.map_open_subpanel = None

        return None

    def sync_map_sizing_sliders_enabled(self):
        """
        Function which enables or disables (greys out, non-interactive)
        the map's manual size/opacity sliders depending on whether
        automatic sizing is on - they have no effect while it is, since
        apply_automatic_marker_style() overwrites whatever they're set to
        on every redraw/zoom/selection change. They stay visible either
        way (just disabled) rather than being hidden, and
        apply_automatic_marker_style() keeps them showing its own current
        values while automatic is on, so switching back to manual has an
        obvious, always-in-place starting point. Called both when the
        automatic-sizing checkbox itself changes and when the sizing
        panel is (re-)expanded.
        """

        manual = not self.map_auto_sizing.isChecked()
        for slider in (
            self.map_markersize_unsel_sl,
            self.map_opacity_unsel_sl,
            self.map_markersize_sel_sl,
            self.map_opacity_sel_sl,
        ):
            set_slider_enabled(slider, manual)

        return None

    @restores_settings_guard
    def handle_map_auto_sizing_update(self):
        """
        Function which handles toggling automatic map point size/opacity
        upon interaction with the map sizing panel's checkbox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_map_auto_sizing_update"
            )
            QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
            self.read_instance.block_config_bar_handling_updates = True

            self.read_instance.map_auto_marker_sizing = self.map_auto_sizing.isChecked()
            self.sync_map_sizing_sliders_enabled()

            if self.read_instance.map_auto_marker_sizing:
                self.apply_automatic_marker_style()

            self.read_instance.block_config_bar_handling_updates = False

            unset_cursor(
                self.read_instance.cursor_function, "handle_map_auto_sizing_update"
            )

        return None

    def apply_automatic_marker_style(self):
        """
        Automatically set the map's unselected/selected marker size and
        opacity as the map is navigated, so points stay identifiable - big
        and well spaced enough to read individually when zoomed in close,
        small enough not to overplot into an unreadable blob when zoomed
        out to (near) the full globe - without the user having to keep
        adjusting the manual sliders. Selected stations get a size/opacity
        boost over unselected ones, which additionally dim back a little
        once there's an active selection, so the selection reads clearly
        against the rest.

        Size comes from get_map_marker_size() in plot_aux.py, shared with
        report and library so the same map reads the same way in every
        mode. Opacity is dashboard only, and tracks the zoom ratio.

        A no-op if automatic sizing is off (read_instance.map_auto_marker_sizing)
        - safe to call unconditionally from every place the map's zoom or
        selection can change, rather than needing each call site to check
        the setting itself.

        Called on: every map redraw (see update_map_station_selection(),
        which every full map rebuild already routes through), scroll-
        wheel zoom (zoom_map_func() in dashboard_interactivity.py),
        toolbar box-zoom and the "world" reset-to-global-view button
        (toolbar.py), and turning automatic sizing on.

        Applies directly to every station on the map, rebuilt from
        scratch each call rather than a partial in-place update -
        deliberately not going through update_markersize()/update_opacity()
        (which the manual sliders use), since those assume the collection
        already has a per-point sizes/colours array to update selected
        indices into. That's not true right after make_map() creates a
        fresh collection (freshly plotted points start with a single
        scalar size/colour until something gives every point its own),
        and update_map_station_selection() - which this runs ahead of, on
        every map rebuild - is exactly the case where that happens.
        """

        if not getattr(self.read_instance, "map_auto_marker_sizing", False):
            return None

        ax = self.plot_axes["map"]

        # current view's extent, against the projection's own global
        # extent as a "fully zoomed out" reference - both in the same
        # (projected) coordinate space, so this ratio is meaningful
        # regardless of which projection is currently selected
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        current_area = abs(xlim[1] - xlim[0]) * abs(ylim[1] - ylim[0])
        x0, x1 = ax.projection.x_limits
        y0, y1 = ax.projection.y_limits
        global_area = abs(x1 - x0) * abs(y1 - y0)

        if (current_area <= 0) or (global_area <= 0):
            zoom_factor = 1.0
        else:
            # sqrt: area shrinks with the *square* of linear zoom, and
            # marker size/opacity should track linear zoom, not area
            zoom_factor = math.sqrt(global_area / current_area)
        zoom_factor = np.clip(zoom_factor, 1.0, MAP_AUTO_SIZING_REFERENCE_ZOOM)

        # size from the density of the stations currently in view, shared with
        # report and library (get_map_marker_size()), so the same map reads the
        # same way in every mode. Density already carries the zoom - fewer
        # stations left in view means bigger markers - as well as how crowded
        # the network is and how large the panel is, none of which a zoom ratio
        # on its own can see
        networkspeci = self.read_instance.networkspeci
        station_inds = getattr(self, "active_map_valid_station_inds", [])
        unsel_size = get_map_marker_size(
            ax,
            self.datacrs,
            self.read_instance.station_longitudes[networkspeci][station_inds],
            self.read_instance.station_latitudes[networkspeci][station_inds],
        )

        # opacity still tracks the zoom ratio directly, scaled between the
        # fully-zoomed-out and fully-zoomed-in references - dashboard only,
        # where panning and zooming makes overplotting come and go
        zoom_progress = (zoom_factor - 1.0) / (MAP_AUTO_SIZING_REFERENCE_ZOOM - 1.0)
        unsel_opacity = MAP_AUTO_SIZING_MIN_OPACITY + (
            MAP_AUTO_SIZING_MAX_OPACITY - MAP_AUTO_SIZING_MIN_OPACITY
        ) * zoom_progress

        selected = getattr(
            self, "absolute_selected_station_inds", np.array([], dtype=np.int32)
        )
        any_selected = len(selected) > 0
        sel_size = unsel_size * MAP_AUTO_SIZING_SELECTED_SIZE_BOOST
        sel_opacity = unsel_opacity
        if any_selected:
            unsel_opacity = unsel_opacity * MAP_AUTO_SIZING_UNSELECTED_OPACITY_DIM

        # keep the underlying config in sync too - read by anything resolving
        # marker style from it directly, and by the manual sliders if
        # automatic sizing is switched back off
        for key, size, opacity in (
            ("marker_unselected", unsel_size, unsel_opacity),
            ("marker_selected", sel_size, sel_opacity),
            ("marker_zero_stations_selected", unsel_size, unsel_opacity),
        ):
            self.plot_characteristics["map"][key]["s"] = size
            self.plot_characteristics["map"][key]["alpha"] = opacity

        # apply directly to every station currently on the map, rebuilt
        # from scratch each call (see the docstring for why)
        n_stations = len(getattr(self, "active_map_valid_station_inds", []))
        if n_stations > 0:
            sizes = np.full(n_stations, unsel_size)
            alphas = np.full(n_stations, unsel_opacity)
            if any_selected:
                sizes[selected] = sel_size
                alphas[selected] = sel_opacity
            for collection in self.plot_axes["map"].collections:
                if isinstance(collection, matplotlib.collections.PathCollection):
                    collection.set_sizes(sizes)
                    if Version(matplotlib.__version__) < Version("3.4"):
                        colours = collection.get_facecolor()
                        if colours.shape[0] == n_stations:
                            colours[:, -1] = alphas
                            collection.set_facecolor(colours)
                    else:
                        collection.set_alpha(alphas)

        # keep the manual sliders in sync with the computed values, so
        # switching to manual mode starts from the current look rather
        # than wherever they were last left
        for slider, value in (
            (self.map_markersize_unsel_sl, unsel_size),
            (self.map_markersize_sel_sl, sel_size),
        ):
            slider.blockSignals(True)
            slider.setValue(int(round(value)))
            slider.blockSignals(False)
        for slider, value in (
            (self.map_opacity_unsel_sl, unsel_opacity),
            (self.map_opacity_sel_sl, sel_opacity),
        ):
            slider.blockSignals(True)
            slider.setValue(int(round(value * 10)))
            slider.blockSignals(False)

        self.figure.canvas.draw_idle()

        return None

    def rebuild_map_axes(self, projection_name):
        """
        Recreate the map axes with a different cartopy projection.

        Cartopy locks a GeoAxes' projection at creation - there's no in-place
        way to change it - so this removes the current map axes and creates
        a fresh one at the same gridspec position, then redraws everything
        onto it (map features, station data, colourbar, and any toggled
        options like domain edges). Resets the view to the new projection's
        global extent, since a remembered pixel extent from the old
        projection doesn't carry meaning in a new one.

        Parameters
        ----------
        projection_name : str
            A valid cartopy.crs projection name (see get_valid_projections()
            in dashboard_elements.py).
        """

        self.plot_characteristics_templates["map"]["projection"] = projection_name
        self.plotcrs = getattr(ccrs, projection_name)()

        self.figure.delaxes(self.plot_axes["map"])
        self.plot_axes["map"] = self.figure.add_subplot(
            self.gridspec.new_subplotspec((2, 0), rowspan=44, colspan=42),
            projection=self.plotcrs,
        )

        # re-apply the map's one-time dressing (features/gridlines) to the
        # new axes; map_extent left at its default (False) resets to a
        # global view rather than reusing the old projection's pixel extent
        format_axis(
            self.read_instance,
            self,
            self.plot_axes["map"],
            "map",
            self.plot_characteristics["map"],
        )
        self.read_instance.map_extent = get_map_extent(self)

        # redraw station data, colourbar, and any toggled options (domain
        # edges/annotations) onto the new axes
        self.update_map_z_statistic()

        # the toolbar's back/forward view history still references the
        # removed axes - reset it rather than leave it holding a dead one
        self.read_instance.navi_toolbar.update()

    def refresh_map_features(self):
        """
        Redraws the map's ocean/land/country-border cartopy features after
        the user changes their colour, visibility, or the map
        resolution from the map settings menu. Does nothing if the map isn't
        using the "providentia" background - a custom background image or
        cartopy's shaded relief doesn't have these features to refresh.
        """

        if self.plot_characteristics["map"]["background"] != "providentia":
            return

        remove_map_features(self.map_feature_artists)
        self.map_feature_artists = draw_map_features(self, self.plot_axes["map"])
        # draw(), not draw_idle() - every caller wraps this in the Providentia
        # busy cursor, and a deferred repaint would land after that cursor
        # had already been restored
        self.figure.canvas.draw()

    def refresh_map_gridlines(self):
        """
        Redraws the map's gridlines after the user toggles them on/off from
        the map settings menu. Applies regardless of map background, unlike
        refresh_map_features().
        """

        if getattr(self, "map_gridliner", None) is not None:
            self.map_gridliner.remove()
        self.map_gridliner = draw_map_gridlines(
            self, self.plot_axes["map"], self.plot_characteristics["map"]["gridlines"]
        )
        # draw(), not draw_idle() - see refresh_map_features()
        self.figure.canvas.draw()

    def handle_timeseries_aggregation_statistic_update(self):
        """
        Function that handles update of timeseries aggregation statistic
        upon interaction with timeseries aggregation statistic combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function,
                "handle_timeseries_aggregation_statistic_update",
            )

            # update timeseries aggregation statistic
            self.update_timeseries_aggregation_statistic()

            # update plotted timeseries statistic
            if not self.read_instance.block_MPL_canvas_updates:
                # update selected data on timeseries plot
                # get selected station data
                get_selected_station_data(
                    read_instance=self.read_instance,
                    canvas_instance=self,
                    networkspecies=[self.read_instance.networkspeci],
                )

                # update plot
                self.update_associated_active_dashboard_plot("timeseries")

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function,
                "handle_timeseries_aggregation_statistic_update",
            )

        return None

    def handle_timeseries_chunk_statistic_update(self):
        """
        Function that handles update of plotted timeseries chunk statistic
        upon interaction with timeseries chunk statistic/resolution comboboxes
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function,
                "handle_timeseries_chunk_statistic_update",
            )

            # update chunk statistic / resolution
            self.update_timeseries_chunk_statistics()

            # get new chunk stat and chunk resolution
            new_chunk_stat = self.timeseries_chunk_stat.currentText()
            new_chunk_resolution = self.timeseries_chunk_resolution.currentText()

            # update plotted timeseries chunk statistic
            if not self.read_instance.block_MPL_canvas_updates:
                # work out if need to update timseries plot or not
                # calculate previous number not None
                if (self.read_instance.previous_chunk_stat == "None") & (
                    self.read_instance.previous_chunk_resolution == "None"
                ):
                    previous_set = 0
                elif (self.read_instance.previous_chunk_stat != "None") & (
                    self.read_instance.previous_chunk_resolution != "None"
                ):
                    previous_set = 2
                else:
                    previous_set = 1

                # calculate new number not None
                if (new_chunk_stat == "None") & (new_chunk_resolution == "None"):
                    new_set = 0
                elif (new_chunk_stat != "None") & (new_chunk_resolution != "None"):
                    new_set = 2
                else:
                    new_set = 1

                # work out if just one of chunk stat or resolution have changed
                if (self.read_instance.previous_chunk_stat != new_chunk_stat) or (
                    self.read_instance.previous_chunk_resolution != new_chunk_resolution
                ):
                    have_changed = True
                else:
                    have_changed = False

                do_update = False
                # if both fields active and at least one has changed from previous then update plot
                if (new_set == 2) & (have_changed):
                    do_update = True
                # elif previous fields were both active and at least one is not, then update plot
                elif (previous_set == 2) & (new_set != 2):
                    do_update = True

                # set previous chunk stat and resolution
                self.read_instance.previous_chunk_stat = new_chunk_stat
                self.read_instance.previous_chunk_resolution = new_chunk_resolution

                # update timeseries plot
                if do_update:
                    # get selected station data
                    get_selected_station_data(
                        read_instance=self.read_instance,
                        canvas_instance=self,
                        networkspecies=[self.read_instance.networkspeci],
                    )

                    # update plot
                    self.update_associated_active_dashboard_plot("timeseries")

                    # draw changes
                    self.figure.canvas.draw_idle()

            # set previous chunk stat and resolution
            else:
                self.read_instance.previous_chunk_stat = new_chunk_stat
                self.read_instance.previous_chunk_resolution = new_chunk_resolution

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function,
                "handle_timeseries_chunk_statistic_update",
            )

        return None

    @restores_settings_guard
    def update_timeseries_chunk_statistics(self):
        """
        Update timeseries chunk statistic and aggregation statistic
        """

        # turn off handling updates to configuration bar
        self.read_instance.block_config_bar_handling_updates = True

        # get currently selected statistic
        chunk_stat = self.timeseries_chunk_stat.currentText()

        # get currently selected resolution
        chunk_resolution = self.timeseries_chunk_resolution.currentText()

        # update timeseries chunk statistics
        if (not self.read_instance.temporal_colocation) or (
            len(self.read_instance.data_labels) == 1
        ):
            available_timeseries_chunk_stats = [
                "None",
            ] + list(copy.deepcopy(self.read_instance.basic_z_stats))
        else:
            available_timeseries_chunk_stats = [
                "None",
            ] + list(copy.deepcopy(self.read_instance.basic_and_bias_z_stats))

        # update available timeseries chunk resolutions
        available_timeseries_chunk_resolutions = [
            "None",
        ] + get_possible_resampling_resolutions(
            self.read_instance.active_resolution,
            daily_forecast=self.read_instance.daily_forecast,
        )

        # if active resolution is not hourly, then MDA8 cannot be available as stat
        if self.read_instance.active_resolution != "hourly":
            available_timeseries_chunk_stats.remove("MDA8")
            if chunk_stat == "MDA8":
                chunk_stat = "None"
                chunk_resolution = "None"
                msg = "The timeseries chunk statistic and chunk resolution will be set to 'None' as MDA8 can only be calculated when the active resolution is hourly."
                show_message(self.read_instance, msg)

        # if chunk stat is MDA8, the chunk resolution has to be None or daily (if available)
        if chunk_stat == "MDA8":
            available_timeseries_chunk_resolutions = ["None", "daily"]
            if chunk_resolution not in available_timeseries_chunk_resolutions:
                chunk_resolution = "None"
                msg = "The timeseries chunk resolution will be set to 'None' as MDA8 can only be calculated for a daily chunk resolution."
                show_message(self.read_instance, msg)

        # if zstat is empty string, it is because fields are being initialised for the first time
        if chunk_stat == "":
            # set timeseries chunk statistic to be None
            chunk_stat = available_timeseries_chunk_stats[0]

        # if resolution is empty string, it is because fields are being initialised for the first time
        if chunk_resolution == "":
            # set timeseries resolution to be None
            chunk_resolution = available_timeseries_chunk_resolutions[0]

        # update timeseries chunk statistic combobox (clear, then add items)
        self.timeseries_chunk_stat.clear()
        self.timeseries_chunk_stat.addItems(available_timeseries_chunk_stats)

        # update timeseries chunk resolution combobox (clear, then add items)
        self.timeseries_chunk_resolution.clear()
        self.timeseries_chunk_resolution.addItems(
            available_timeseries_chunk_resolutions
        )

        # maintain currently selected timeseries resolution (if exists in new item list)
        if chunk_resolution in available_timeseries_chunk_resolutions:
            self.timeseries_chunk_resolution.setCurrentText(chunk_resolution)
        # if does not exist then chunk resolution is None, so make chunk stat that also
        # also throw error stating why it happened
        else:
            chunk_stat = "None"
            msg = "Timeseries chunk resolution and statistic will be set to 'None' "
            msg += (
                f"as the active chunk resolution ({chunk_resolution}) is not valid for "
            )
            msg += "the active resolution or resampling resolution. "
            show_message(self.read_instance, msg)

        # maintain currently selected timeseries chunk statistic (if exists in new item list)
        if chunk_stat in available_timeseries_chunk_stats:
            self.timeseries_chunk_stat.setCurrentText(chunk_stat)
        # if does not exist then chunk stat is None, so make chunk resolution that also
        else:
            self.timeseries_chunk_resolution.setCurrentText("None")

        # allow handling updates to the configuration bar again
        self.read_instance.block_config_bar_handling_updates = False

    @restores_settings_guard
    def handle_periodic_statistic_update(self):
        """
        Function that handles update of plotted periodic statistic
        upon interaction with periodic statistic combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_periodic_statistic_update"
            )

            # set variable that blocks configuration bar handling updates until all changes
            # to the periodic statistic combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected statistic
            zstat = self.periodic_stat.currentText()

            # update periodic statistics, to all basic stats
            # if colocation not-active, and basic+bias stats if colocation active
            if (not self.read_instance.temporal_colocation) or (
                len(self.read_instance.data_labels) == 1
            ):
                available_periodic_stats = copy.deepcopy(
                    self.read_instance.basic_z_stats
                )
            else:
                available_periodic_stats = copy.deepcopy(
                    self.read_instance.basic_and_bias_z_stats
                )

            # remove MDA8 from available stats
            if "MDA8" in available_periodic_stats:
                available_periodic_stats = np.delete(
                    available_periodic_stats,
                    np.where(available_periodic_stats == "MDA8")[0],
                )

            # if base_zstat is empty string, it is because fields are being initialised for the first time
            if zstat == "":
                # set periodic stat to be first available stat
                zstat = available_periodic_stats[0]

            # update periodic statistic combobox (clear, then add items)
            self.periodic_stat.clear()
            self.periodic_stat.addItems(available_periodic_stats)

            # maintain currently selected periodic statistic (if exists in new item list)
            if zstat in available_periodic_stats:
                self.periodic_stat.setCurrentText(zstat)
            elif zstat == "MDA8":
                msg = f"Periodic statistic is being reset to {self.periodic_stat.currentText()}. MDA8 can only be calculated when the active resolution is hourly."
                show_message(self.read_instance, msg)

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted periodic statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot("periodic")

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function, "handle_periodic_statistic_update"
            )

        return None

    @restores_settings_guard
    def handle_taylor_correlation_statistic_update(self):
        """
        Function that handles update of correlation statistic
        upon interaction with Taylor diagram statistic combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function,
                "handle_taylor_correlation_statistic_update",
            )

            # set variable that blocks configuration bar handling updates until all
            # changes to the statistic combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected items
            corr_stat = self.taylor_corr_stat.currentText()

            # get available stats
            available_corr_stats = ["r", "r2"]

            # if correlation stat is empty string, it is because fields are being initialised for the first time
            if corr_stat == "":
                # set stat to be the one in plot characteristics
                corr_stat = available_corr_stats[0]

            # update statistic combobox (clear, then add items)
            self.taylor_corr_stat.clear()
            self.taylor_corr_stat.addItems(available_corr_stats)

            # maintain currently selected statistic
            self.taylor_corr_stat.setCurrentText(corr_stat)

            # update dictionary
            self.plot_characteristics["taylor"]["corr_stat"] = corr_stat

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted taylor diagram statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot("taylor")

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function,
                "handle_taylor_correlation_statistic_update",
            )

        return None

    def get_active_statsummary_stats(self, statistic_type):
        """
        Get active statistics from dictionary of statsummary statistics in list

        Parameters
        ----------
        statistic_type : str
            Statistic type
        """

        active_statsummary_stats = copy.deepcopy(
            self.read_instance.current_statsummary_stats[statistic_type]
        )
        active_statsummary_stats = [
            stat
            for sublist in list(active_statsummary_stats.values())
            for stat in sublist
        ]

        return active_statsummary_stats

    def check_statsummary_stats(self):
        """
        Function that checks the statistics in the statsummary statistic combobox
        """

        # get stats to check for the selected periodic cycle
        periodic_cycle = self.statsummary_cycle.currentText()
        if periodic_cycle == "":
            periodic_cycle = "None"
        plot_options = self.current_plot_options["statsummary"]
        statistic_type = "basic" if "bias" not in plot_options else "modbias"
        if "bias" in plot_options:
            items = list(copy.deepcopy(self.read_instance.basic_and_bias_z_stats))
        else:
            items = list(copy.deepcopy(self.read_instance.basic_z_stats))
        # remove MDA8 as option if periodic cycle is not None
        if ("MDA8" in items) & (periodic_cycle != "None"):
            items.remove("MDA8")
        if periodic_cycle != "None":
            items = [stat + "-" + periodic_cycle.lower() for stat in items]
        self.statsummary_stat.clear()
        self.statsummary_stat.addItems(items)
        checked_options = copy.deepcopy(
            self.read_instance.current_statsummary_stats[statistic_type][periodic_cycle]
        )
        checked_options = [
            option.split("_bias")[0] if "_bias" in option else option
            for option in checked_options
        ]
        checked_options_in_items = list(set(checked_options) & set(items))

        # check stats in combobox
        if checked_options_in_items:
            for checked_option in checked_options_in_items:
                index = items.index(checked_option)
                self.statsummary_stat.model().item(index).setCheckState(
                    QtCore.Qt.Checked
                )
        # leave empty
        else:
            self.statsummary_stat.lineEdit().setText("")

    @restores_settings_guard
    def handle_statsummary_statistics_update(self):
        """
        Function that handles update of plotted statsummary statistics
        upon interaction with statistic comboboxes
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function,
                "handle_statsummary_statistics_update",
            )

            # set variable that blocks configuration bar handling updates until all changes
            # to the statsummary statistics combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get all possible stats
            plot_options = self.current_plot_options["statsummary"]
            statistic_type = "basic" if "bias" not in plot_options else "modbias"

            # initialise stats
            if not hasattr(self, "active_statsummary_stats"):
                # get initial stats from plot characteristics
                periodic_cycle = "None"
                self.read_instance.current_statsummary_stats["basic"][
                    "None"
                ] = self.plot_characteristics["statsummary"]["basic"]
                self.read_instance.current_statsummary_stats["modbias"][
                    "None"
                ] = self.plot_characteristics["statsummary"]["model_bias"]
                self.active_statsummary_stats = {
                    "basic": self.get_active_statsummary_stats("basic"),
                    "modbias": self.get_active_statsummary_stats("modbias"),
                }

                # check stats for the selected stats
                self.check_statsummary_stats()

            # get stats from selection
            else:
                # save previous stats in list
                previous_active_statsummary_stats = copy.deepcopy(
                    self.active_statsummary_stats[statistic_type]
                )

                # remove bias from options to get correct active stats
                if statistic_type == "modbias":
                    previous_active_statsummary_stats = [
                        option.split("_bias")[0] if "_bias" in option else option
                        for option in previous_active_statsummary_stats
                    ]

                # update stats
                periodic_cycle = self.statsummary_cycle.currentText()
                self.read_instance.current_statsummary_stats[statistic_type][
                    periodic_cycle
                ] = copy.deepcopy(self.statsummary_stat.currentData())

                # save current stats in list
                current_active_statsummary_stats = copy.deepcopy(
                    self.get_active_statsummary_stats(statistic_type)
                )

                # check stats for the selected periodic cycle
                self.check_statsummary_stats()

                # get active stats
                current_not_previous = list(
                    set(current_active_statsummary_stats).difference(
                        previous_active_statsummary_stats
                    )
                )
                previous_not_current = list(
                    set(previous_active_statsummary_stats).difference(
                        current_active_statsummary_stats
                    )
                )
                self.active_statsummary_stats[statistic_type] = copy.deepcopy(
                    previous_active_statsummary_stats
                )
                if current_not_previous:
                    for stat in current_not_previous:
                        self.active_statsummary_stats[statistic_type].append(stat)
                if previous_not_current:
                    for stat in previous_not_current:
                        self.active_statsummary_stats[statistic_type].remove(stat)

            # remove MDA8 if active resolution is not hourly
            if self.read_instance.active_resolution != "hourly":
                basic_stats_to_remove = []
                bias_stats_to_remove = []
                for stat in self.active_statsummary_stats["basic"]:
                    if "MDA8" in stat:
                        basic_stats_to_remove.append(stat)
                for stat in self.active_statsummary_stats["modbias"]:
                    if "MDA8" in stat:
                        bias_stats_to_remove.append(stat)

                if (len(basic_stats_to_remove) > 0) or (len(bias_stats_to_remove) > 0):
                    for stat in basic_stats_to_remove:
                        self.active_statsummary_stats["basic"].remove(stat)
                    for stat in bias_stats_to_remove:
                        self.active_statsummary_stats["modbias"].remove(stat)

                    for stat_type in self.read_instance.current_statsummary_stats:
                        for resolution in self.read_instance.current_statsummary_stats[
                            stat_type
                        ]:
                            stats_to_remove = []
                            for stat in self.read_instance.current_statsummary_stats[
                                stat_type
                            ][resolution]:
                                if "MDA8" in stat:
                                    stats_to_remove.append(stat)
                            for stat in stats_to_remove:
                                self.read_instance.current_statsummary_stats[stat_type][
                                    resolution
                                ].remove(stat)

                    msg = "Removing all MDA8 statistics from statsummary plot as active resolution is not hourly."
                    show_message(self.read_instance, msg)

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted statsummary statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot("statsummary")

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function,
                "handle_statsummary_statistics_update",
            )

    @restores_settings_guard
    def handle_statsummary_cycle_update(self):
        """
        Function that handles update of statsummary periodic cycle
        upon interaction with statistic comboboxes
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "handle_statsummary_cycle_update"
            )

            # set variable that blocks configuration bar handling updates until all changes
            # to the statsummary periodic cycle combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # get currently selected cycle
            periodic_cycle = self.statsummary_cycle.currentText()

            # update periodi cycles
            available_cycles = ["None", "Diurnal", "Weekly", "Monthly"]

            # if cycle is empty string, it is because fields are being initialised for the first time
            if periodic_cycle == "":
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
                self.update_associated_active_dashboard_plot("statsummary")

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function, "handle_statsummary_cycle_update"
            )

        return None

    @restores_settings_guard
    def handle_statsummary_periodic_aggregation_update(self):
        """
        Function that handles update of plotted statsummary periodic aggregation statistic
        upon interaction with combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function,
                "handle_statsummary_periodic_aggregation_update",
            )

            # set variable that blocks configuration bar handling updates until all changes
            # to the statsummary periodic aggregation combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # update statistic
            self.read_instance.selected_periodic_statistic_aggregation = (
                self.statsummary_periodic_aggregation.currentText()
            )
            self.read_instance.periodic_statistic_aggregation = (
                self.read_instance.selected_periodic_statistic_aggregation
            )

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted statsummary statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot("statsummary")

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function,
                "handle_statsummary_periodic_aggregation_update",
            )

        return None

    @restores_settings_guard
    def handle_statsummary_periodic_mode_update(self):
        """
        Function that handles update of plotted statsummary periodic aggregation mode
        upon interaction with combobox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function,
                "handle_statsummary_periodic_mode_update",
            )

            # set variable that blocks configuration bar handling updates until all changes
            # to the statsummary periodic mode combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # update statistic mode
            self.read_instance.selected_periodic_statistic_mode = (
                self.statsummary_periodic_mode.currentText()
            )
            self.read_instance.periodic_statistic_mode = (
                self.read_instance.selected_periodic_statistic_mode
            )

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update plotted statsummary statistic
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot("statsummary")

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function,
                "handle_statsummary_periodic_mode_update",
            )

        return None

    @restores_settings_guard
    def handle_fairmode_target_classification_update(self):
        """
        Function that handles update of station classification on FAIRMODE target plot
        upon interaction with temporal colocation checkbox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # update mouse cursor to a waiting cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function,
                "handle_fairmode_target_classification_update",
            )

            # set variable that blocks configuration bar handling updates until all changes
            # to the classification combobox are made
            self.read_instance.block_config_bar_handling_updates = True

            # update classification type
            self.plot_characteristics["fairmode-target"]["markers"][
                "type"
            ] = self.fairmode_target_classification.currentText()

            # allow handling updates to the configuration bar again
            self.read_instance.block_config_bar_handling_updates = False

            # update FAIRMODE target plot
            if not self.read_instance.block_MPL_canvas_updates:
                self.update_associated_active_dashboard_plot("fairmode-target")

            # draw changes
            self.figure.canvas.draw_idle()

            # restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function,
                "handle_fairmode_target_classification_update",
            )

        return None

    def remove_axis_objects(
        self, ax_elements, elements_to_skip=None, types_to_remove=None
    ):
        """
        Remove objects (artists, lines, collections, patches) from axis
        """

        # define default argument mutables
        if elements_to_skip is None:
            elements_to_skip = []
        if types_to_remove is None:
            types_to_remove = []

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
        """
        Remove all plotted axis elements
        """

        # get appropriate axes for nested axes
        axs_to_remove = []
        if isinstance(ax, dict):
            axs_to_remove = list(ax.values())
        elif isinstance(ax, list):
            axs_to_remove = ax
        else:
            if plot_type == "taylor":
                axs_to_remove.append(self.plotting.taylor_polar_relevant_axis)
            axs_to_remove.append(ax)

        # iterate through axes
        for ax_to_remove in axs_to_remove:
            # remove all plotted axis elements
            if plot_type == "legend":
                leg = ax_to_remove.get_legend()
                if leg:
                    leg.remove()

            elif plot_type == "map":
                self.remove_axis_objects(
                    ax_to_remove.artists, types_to_remove=[AnchoredOffsetbox]
                )
                self.remove_axis_objects(
                    ax_to_remove.collections,
                    types_to_remove=[matplotlib.collections.PathCollection],
                )
                # # TODO: Put line collection back into place when we turn on the auto_update in gridlines
                # self.remove_axis_objects(ax_to_remove.collections, types_to_remove=[matplotlib.collections.PathCollection],
                #                                                                     matplotlib.collections.LineCollection])

            elif plot_type == "cb":
                for objects in [ax_to_remove.artists, ax_to_remove.collections]:
                    self.remove_axis_objects(objects)

            elif plot_type == "timeseries":
                for objects in [ax_to_remove.lines, ax_to_remove.artists]:
                    self.remove_axis_objects(objects)

            elif plot_type == "periodic":
                for objects in [ax_to_remove.lines, ax_to_remove.artists]:
                    self.remove_axis_objects(objects)

            elif plot_type == "periodic-violin":
                for objects in [
                    ax_to_remove.lines,
                    ax_to_remove.artists,
                    ax_to_remove.collections,
                ]:
                    self.remove_axis_objects(objects)

            elif plot_type == "metadata":
                self.remove_axis_objects(ax_to_remove.texts)

            elif plot_type in ["distribution", "histogram"]:
                for objects in [ax_to_remove.lines, ax_to_remove.artists]:
                    self.remove_axis_objects(objects)

            elif plot_type in ["statsummary", "contingencytable"]:
                self.remove_axis_objects(ax_to_remove.tables)

            elif plot_type in ["taylor", "scatter"]:
                for objects in [ax_to_remove.lines, ax_to_remove.artists]:
                    self.remove_axis_objects(objects)

            elif plot_type in ["boxplot"]:
                for objects in [
                    ax_to_remove.artists,
                    ax_to_remove.patches,
                    ax_to_remove.lines,
                ]:
                    self.remove_axis_objects(objects)

            elif plot_type == "fairmode-target":
                for objects in [
                    ax_to_remove.lines,
                    ax_to_remove.artists,
                    ax_to_remove.patches,
                    ax_to_remove.texts,
                ]:
                    self.remove_axis_objects(objects)

            elif plot_type == "fairmode-statsummary":
                self.remove_axis_objects(ax_to_remove.lines)

        # remove tracked plot elements
        if plot_type in self.plot_elements:
            self.plot_elements[plot_type]["absolute"] = {}
            if "bias" in self.plot_elements[plot_type]:
                del self.plot_elements[plot_type]["bias"]

        return None

    def update_plot_options(self, plot_types):
        """
        Uncheck checked boxes in plot configuration options under menus and
        reapply check with new data. This can be done for all currently active plot types,
        or just one specific type

        Parameters
        ----------
        plot_type : str
            Plot type
        """

        if not isinstance(plot_types, list):
            plot_types = [plot_types]

        for plot_type in plot_types:
            all_plot_options = self.plot_characteristics[plot_type]["plot_options"]
            checked_options = self.current_plot_options[plot_type]
            if plot_type in [
                "periodic-violin",
                "fairmode-target",
                "fairmode-statsummary",
            ]:
                plot_type = plot_type.replace("-", "_")
            cb_options = getattr(self, plot_type + "_options")

            if plot_type == "contingencytable":
                # if more than one station is selected, select gerrity as it is already showing
                n_stations = self.selected_station_data[
                    self.read_instance.networkspeci
                ]["per_station"].shape[1]

            # uncheck all options
            for option in all_plot_options:
                index = all_plot_options.index(option)
                if plot_type == "contingencytable" and n_stations > 1:
                    state = QtCore.Qt.Checked
                else:
                    state = QtCore.Qt.Unchecked
                self.read_instance.block_MPL_canvas_updates = True
                cb_options.model().item(index).setCheckState(state)
                self.read_instance.block_MPL_canvas_updates = False

            # check selected options
            for checked_option_ii, checked_option in enumerate(checked_options):
                index = all_plot_options.index(checked_option)
                if checked_option_ii < (len(checked_options) - 1):
                    self.read_instance.block_MPL_canvas_updates = True
                    cb_options.model().item(index).setCheckState(QtCore.Qt.Checked)
                    self.read_instance.block_MPL_canvas_updates = False
                else:
                    cb_options.model().item(index).setCheckState(QtCore.Qt.Checked)

        return None

    def select_all_stations(self):
        """
        Function that selects/unselects all plotted stations
        (and associated plots) upon ticking of checkbox
        """

        if not self.read_instance.block_MPL_canvas_updates:
            # check if checkbox to select all stations is checked or unchecked
            check_state = self.read_instance.ch_select_all.checkState()

            # show warning and uncheck box
            if not hasattr(self, "relative_selected_station_inds"):
                if check_state == QtCore.Qt.Checked:
                    msg = "Data must be read into memory before selecting the data."
                    show_message(self.read_instance, msg)
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_select_all.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False
                    return

            # set mouse cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "select_all_stations"
            )

            # make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(
                self.relative_selected_station_inds
            )

            # if checkbox is checked, select all plotted stations
            if check_state == QtCore.Qt.Checked:
                self.relative_selected_station_inds = copy.deepcopy(
                    self.active_map_valid_station_inds
                )

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
                self.relative_selected_station_inds = np.array([], dtype=np.int32)

            # update absolute selected station indices (indices relative to plotted stations on map)
            self.absolute_selected_station_inds = np.arange(
                len(self.relative_selected_station_inds), dtype=np.int32
            )

            # get absolute non-selected station inds
            self.absolute_non_selected_station_inds = np.nonzero(
                ~np.in1d(
                    range(len(self.active_map_valid_station_inds)),
                    self.absolute_selected_station_inds,
                )
            )[0]

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(
                self.previous_relative_selected_station_inds,
                self.relative_selected_station_inds,
            ):
                self.update_associated_active_dashboard_plots()

                # draw changes
                self.figure.canvas.draw_idle()

            # Restore mouse cursor to normal
            unset_cursor(self.read_instance.cursor_function, "select_all_stations")

        return None

    def select_intersect_stations(self):
        """
        Function that selects/unselects intersection of
        stations and all model domains (and associated plots)
        upon ticking of checkbox
        """

        if not self.read_instance.block_MPL_canvas_updates:
            # check if checkbox to select intersection of stations is checked or unchecked
            check_state = self.read_instance.ch_intersect.checkState()

            # show warning and uncheck box
            if not hasattr(self, "relative_selected_station_inds"):
                if check_state == QtCore.Qt.Checked:
                    msg = "Data must be read into memory before selecting the data."
                    show_message(self.read_instance, msg)
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_intersect.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False
                    return

            # set mouse cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "select_intersect_stations"
            )

            # make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(
                self.relative_selected_station_inds
            )

            # if checkbox is unchecked then unselect all plotted stations
            if check_state == QtCore.Qt.Unchecked:
                self.relative_selected_station_inds = np.array([], dtype=np.int32)
                self.absolute_selected_station_inds = np.array([], dtype=np.int32)
                self.absolute_non_selected_station_inds = np.arange(
                    len(self.relative_selected_station_inds), dtype=np.int32
                )

            # else, if checkbox is checked then select all stations which intersect with all loaded model domains
            elif check_state == QtCore.Qt.Checked:
                # if have only observations loaded into memory, select all plotted stations
                if len(self.read_instance.data_labels) == 1:
                    self.relative_selected_station_inds = copy.deepcopy(
                        self.active_map_valid_station_inds
                    )
                    self.absolute_selected_station_inds = np.arange(
                        len(self.relative_selected_station_inds), dtype=np.int32
                    )
                    self.absolute_non_selected_station_inds = np.array(
                        [], dtype=np.int32
                    )
                # else, define list of lists to get intersection between (active_map_valid_station_inds,
                # and valid station indices associated with each loaded model array)
                else:
                    intersect_lists = [self.active_map_valid_station_inds]
                    for data_label in self.read_instance.data_labels:
                        if data_label != self.read_instance.observations_data_label:
                            if self.read_instance.temporal_colocation:
                                valid_station_inds = self.read_instance.valid_station_inds_temporal_colocation[
                                    self.read_instance.networkspeci
                                ][
                                    data_label
                                ]
                            else:
                                valid_station_inds = (
                                    self.read_instance.valid_station_inds[
                                        self.read_instance.networkspeci
                                    ][data_label]
                                )
                            intersect_lists.append(valid_station_inds)

                    # get intersect between active map valid station indices and valid station indices
                    # associated with each loaded model array --> relative selected station indcies
                    self.relative_selected_station_inds = np.sort(
                        list(set.intersection(*map(set, intersect_lists)))
                    )

                    # get absolute selected station indices (indices relative to plotted stations on map)
                    self.absolute_selected_station_inds = np.array(
                        [
                            np.where(
                                self.active_map_valid_station_inds == selected_ind
                            )[0][0]
                            for selected_ind in self.relative_selected_station_inds
                        ],
                        dtype=np.int32,
                    )

                    # get absolute non-selected station inds
                    self.absolute_non_selected_station_inds = np.nonzero(
                        ~np.in1d(
                            range(len(self.active_map_valid_station_inds)),
                            self.absolute_selected_station_inds,
                        )
                    )[0]

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
            if not np.array_equal(
                self.previous_relative_selected_station_inds,
                self.relative_selected_station_inds,
            ):
                self.update_associated_active_dashboard_plots()

                # draw changes
                self.figure.canvas.draw_idle()

            # Restore mouse cursor to normal
            unset_cursor(
                self.read_instance.cursor_function, "select_intersect_stations"
            )

        return None

    def select_extent_stations(self):
        """
        Function that selects/unselects the
        stations for the current map extent (and associated plots)
        upon ticking of checkbox
        """

        if not self.read_instance.block_MPL_canvas_updates:
            # check if checkbox to select extent of stations is checked or unchecked
            check_state = self.read_instance.ch_extent.checkState()

            # show warning and uncheck box
            if not hasattr(self, "relative_selected_station_inds"):
                if check_state == QtCore.Qt.Checked:
                    msg = "Data must be read into memory before selecting the data."
                    show_message(self.read_instance, msg)
                    self.read_instance.block_MPL_canvas_updates = True
                    self.read_instance.ch_extent.setCheckState(QtCore.Qt.Unchecked)
                    self.read_instance.block_MPL_canvas_updates = False
                    return

            # set mouse cursor
            self.read_instance.cursor_function = set_cursor(
                self.read_instance.cursor_function, "select_extent_stations"
            )

            # get map extent (in data coords)
            self.read_instance.map_extent = get_map_extent(self)

            # make copy of current full array relative selected stations indices, before selecting new ones
            self.previous_relative_selected_station_inds = copy.deepcopy(
                self.relative_selected_station_inds
            )

            # if checkbox is checked, select all plotted stations
            if check_state == QtCore.Qt.Checked:
                # make copy of current full array relative selected stations indices, before selecting new ones
                self.relative_selected_station_inds = copy.deepcopy(
                    self.active_map_valid_station_inds
                )

                # get inds of stations inside map extent
                extent_station_inds = []
                for ind, lon, lat in zip(
                    self.relative_selected_station_inds,
                    self.read_instance.station_longitudes[
                        self.read_instance.networkspeci
                    ][self.relative_selected_station_inds],
                    self.read_instance.station_latitudes[
                        self.read_instance.networkspeci
                    ][self.relative_selected_station_inds],
                ):
                    if (
                        (lon >= self.read_instance.map_extent[0])
                        and (lon <= self.read_instance.map_extent[1])
                        and (lat >= self.read_instance.map_extent[2])
                        and (lat <= self.read_instance.map_extent[3])
                    ):
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
                self.relative_selected_station_inds = np.array([], dtype=np.int32)

            # get absolute selected station indices (indices relative to plotted stations on map)
            self.absolute_selected_station_inds = np.array(
                [
                    np.where(self.active_map_valid_station_inds == selected_ind)[0][0]
                    for selected_ind in self.relative_selected_station_inds
                ],
                dtype=np.int32,
            )

            # get absolute unselected station indices
            self.absolute_non_selected_station_inds = np.nonzero(
                ~np.in1d(
                    range(len(self.active_map_valid_station_inds)),
                    self.absolute_selected_station_inds,
                )
            )[0]

            # update map station selection
            self.update_map_station_selection()

            # if selected stations have changed from previous selected, update associated plots
            if not np.array_equal(
                self.previous_relative_selected_station_inds,
                self.relative_selected_station_inds,
            ):
                self.update_associated_active_dashboard_plots()

                # draw changes
                self.figure.canvas.draw_idle()

            # Restore mouse cursor to normal
            unset_cursor(self.read_instance.cursor_function, "select_extent_stations")

        return None

    def station_select(self, event):
        """
        Select station on map

        Parameters
        ----------
        event : matplotlib.backend_bases.Event
            Event
        """

        # return if not on map axis
        if event.inaxes != self.plot_axes["map"]:
            return

        # return if lasso active
        if self.lasso_active:
            return

        # check if have any plotted stations on map, if not, return
        if len(self.active_map_valid_station_inds) == 0:
            return

        # if canvas drawing is locked, then return if not owner
        if self.figure.canvas.widgetlock.locked():
            if not self.figure.canvas.widgetlock.isowner(self.station_pick):
                return
        # else, lock canvas drawing
        else:
            self.figure.canvas.widgetlock(self.station_pick)

        # set mouse cursor
        self.read_instance.cursor_function = set_cursor(
            self.read_instance.cursor_function, "station_select"
        )

        # unselect all/intersect/extent checkboxes
        self.unselect_map_checkboxes()

        # make copy of current full array absolute abd relative selected stations indices, before selecting new ones
        previous_absolute_selected_station_inds = copy.deepcopy(
            self.absolute_selected_station_inds
        )
        previous_relative_selected_station_inds = copy.deepcopy(
            self.relative_selected_station_inds
        )

        # get coordinates of selected point
        verts = [(event.xdata, event.ydata)]
        lasso_path = Path(verts)
        lasso_path_vertices = lasso_path.vertices

        # transform lasso coordinates from projected to standard geographic coordinates
        lasso_path.vertices = self.datacrs.transform_points(
            self.plotcrs, lasso_path_vertices[:, 0], lasso_path_vertices[:, 1]
        )[:, :2]

        # get absolute selected indices of stations on map
        absolute_selected_station_inds = np.nonzero(
            lasso_path.contains_points(self.map_points_coordinates)
        )[0]

        # if have no valid selected indices, add a small tolerance (variable by visible map extent) to try get a match
        if len(absolute_selected_station_inds) == 0:
            # take first selected point coordinates and get matches of stations within tolerance
            self.read_instance.map_extent = get_map_extent(self)
            tolerance = (
                np.average(
                    [
                        self.read_instance.map_extent[1]
                        - self.read_instance.map_extent[0],
                        self.read_instance.map_extent[3]
                        - self.read_instance.map_extent[2],
                    ]
                )
                / 100.0
            )
            point_coordinates = lasso_path.vertices[0:1, :]
            sub_abs_vals = np.abs(
                self.map_points_coordinates[None, :, :] - point_coordinates[:, None, :]
            )
            absolute_selected_station_inds = np.arange(
                len(self.active_map_valid_station_inds)
            )[np.all(np.any(sub_abs_vals <= tolerance, axis=0), axis=1)]

            # if more than 1 point selected, limit this to be just nearest point
            if len(absolute_selected_station_inds) > 1:
                absolute_selected_station_inds = np.array(
                    [
                        absolute_selected_station_inds[
                            np.argmin(
                                np.sum(
                                    sub_abs_vals[0, absolute_selected_station_inds, :],
                                    axis=1,
                                )
                            )
                        ]
                    ],
                    dtype=np.int32,
                )

        # handle left click event
        if event.button is MouseButton.LEFT:
            # set absolute selected inds to self
            self.absolute_selected_station_inds = absolute_selected_station_inds

        # handle right click event
        elif event.button is MouseButton.RIGHT:
            # if have zero stations selected then return, doing nothing to selection
            if len(absolute_selected_station_inds) == 0:
                # restore mouse cursor to normal
                unset_cursor(self.read_instance.cursor_function, "station_select")
                return

            # update absolute indices of selected stations
            # remove stations that were previously selected, and add stations that were not previously selected
            self.absolute_selected_station_inds = np.setxor1d(
                previous_absolute_selected_station_inds, absolute_selected_station_inds
            )

        # update previous selected absolute and relative inds
        self.previous_absolute_selected_station_inds = (
            previous_absolute_selected_station_inds
        )
        self.previous_relative_selected_station_inds = (
            previous_relative_selected_station_inds
        )

        # get absolute non-selected station inds
        self.absolute_non_selected_station_inds = np.nonzero(
            ~np.in1d(
                range(len(self.active_map_valid_station_inds)),
                self.absolute_selected_station_inds,
            )
        )[0]

        # get new relative selected indices with respect to all available stations
        self.relative_selected_station_inds = (
            self.map_selected_station_inds_to_all_available_inds(
                self.absolute_selected_station_inds
            )
        )

        # if selected stations have changed from previous selected, update station selection and associated plots
        if not np.array_equal(
            self.previous_relative_selected_station_inds,
            self.relative_selected_station_inds,
        ):
            self.update_map_station_selection()
            self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw_idle()

        # unlock canvas drawing
        if self.figure.canvas.widgetlock.isowner(self.station_pick):
            self.figure.canvas.widgetlock.release(self.station_pick)

        # restore mouse cursor to normal
        unset_cursor(self.read_instance.cursor_function, "station_select")

    def onlassoselect(self, verts):
        """
        Function that handles station selection upon lasso selection with left click.

        Operation:
        Select all stations within lasso boundaries.
        If a click is made rather than using lasso, then select nearest station within tolerance.

        If no station is found with left click, all stations are unselected.

        Parameters
        ----------
        verts : array
            Path vertices
        """

        # check if have any plotted stations on map, if not, return
        if len(self.active_map_valid_station_inds) == 0:
            return

        # set mouse cursor
        self.read_instance.cursor_function = set_cursor(
            self.read_instance.cursor_function, "onlassoselect"
        )

        # unselect all/intersect/extent checkboxes
        self.unselect_map_checkboxes()

        # make copy of current full array absolute abd relative selected stations indices, before selecting new ones
        self.previous_absolute_selected_station_inds = copy.deepcopy(
            self.absolute_selected_station_inds
        )
        self.previous_relative_selected_station_inds = copy.deepcopy(
            self.relative_selected_station_inds
        )

        # get coordinates of drawn lasso
        lasso_path = Path(verts)
        lasso_path_vertices = lasso_path.vertices

        # transform lasso coordinates from projected to standard geographic coordinates
        lasso_path.vertices = self.datacrs.transform_points(
            self.plotcrs, lasso_path_vertices[:, 0], lasso_path_vertices[:, 1]
        )[:, :2]

        # get absolute selected indices of stations on map (the station coordinates contained within lasso)
        self.absolute_selected_station_inds = np.nonzero(
            lasso_path.contains_points(self.map_points_coordinates)
        )[0]

        # if have no valid selected indices, add a small tolerance (variable by visible map extent) to try get a match
        if len(self.absolute_selected_station_inds) == 0:
            # take first selected point coordinates and get matches of stations within tolerance
            self.read_instance.map_extent = get_map_extent(self)
            tolerance = (
                np.average(
                    [
                        self.read_instance.map_extent[1]
                        - self.read_instance.map_extent[0],
                        self.read_instance.map_extent[3]
                        - self.read_instance.map_extent[2],
                    ]
                )
                / 100.0
            )
            point_coordinates = lasso_path.vertices[0:1, :]
            sub_abs_vals = np.abs(
                self.map_points_coordinates[None, :, :] - point_coordinates[:, None, :]
            )
            self.absolute_selected_station_inds = np.arange(
                len(self.active_map_valid_station_inds)
            )[np.all(np.any(sub_abs_vals <= tolerance, axis=0), axis=1)]
            # if more than 1 point selected, limit this to be just nearest point
            if len(self.absolute_selected_station_inds) > 1:
                self.absolute_selected_station_inds = np.array(
                    [
                        self.absolute_selected_station_inds[
                            np.argmin(
                                np.sum(
                                    sub_abs_vals[
                                        0, self.absolute_selected_station_inds, :
                                    ],
                                    axis=1,
                                )
                            )
                        ]
                    ],
                    dtype=np.int32,
                )

        # get absolute non-selected station inds
        self.absolute_non_selected_station_inds = np.nonzero(
            ~np.in1d(
                range(len(self.active_map_valid_station_inds)),
                self.absolute_selected_station_inds,
            )
        )[0]

        # get selected station indices with respect to all available stations
        self.relative_selected_station_inds = (
            self.map_selected_station_inds_to_all_available_inds(
                self.absolute_selected_station_inds
            )
        )

        # hide lasso after selection
        self.lasso_event.set_visible(False)

        # if selected stations have changed from previous selected, update station selection and associated plots
        if not np.array_equal(
            self.previous_relative_selected_station_inds,
            self.relative_selected_station_inds,
        ):
            self.update_map_station_selection()
            self.update_associated_active_dashboard_plots()

            # draw changes
            self.figure.canvas.draw_idle()

        # restore mouse cursor to normal
        unset_cursor(self.read_instance.cursor_function, "onlassoselect")

        return None

    def map_selected_station_inds_to_all_available_inds(self, selected_map_inds):
        """
        Take the indices of selected stations on the map
        (potentially a subset of all available stations), and returns the indices
        of the stations inside the full loaded data arrays

        Parameters
        ----------
        selected_map_inds : numpy.array
            Selected station indices
        """

        # index the array of indices of stations plotted on the map (indexed with respect to
        # all available stations), with the absolute indices of the subset of plotted selected stations
        return self.active_map_valid_station_inds[selected_map_inds]

    def generate_interactive_elements(self):
        """
        Function to create settings menus for each plot and their elements
        """

        self.interactive_elements = {}

        # LAYOUT OPTIONS #
        # add position 2 plot selector
        self.read_instance.cb_position_2 = set_formatting(
            ComboBox(self), self.read_instance.formatting_dict["menu_combobox"]
        )
        self.read_instance.cb_position_2.setToolTip(
            "Select plot type in top right position"
        )
        self.read_instance.cb_position_2.currentTextChanged.connect(
            self.read_instance.handle_layout_update
        )

        # add position 3 plot selector
        self.read_instance.cb_position_3 = set_formatting(
            ComboBox(self), self.read_instance.formatting_dict["menu_combobox"]
        )
        self.read_instance.cb_position_3.setToolTip(
            "Select plot type in bottom left position"
        )
        self.read_instance.cb_position_3.currentTextChanged.connect(
            self.read_instance.handle_layout_update
        )

        # add position 4 plot selector
        self.read_instance.cb_position_4 = set_formatting(
            ComboBox(self), self.read_instance.formatting_dict["menu_combobox"]
        )
        self.read_instance.cb_position_4.setToolTip(
            "Select plot type in bottom centre position"
        )
        self.read_instance.cb_position_4.currentTextChanged.connect(
            self.read_instance.handle_layout_update
        )

        # add position 5 plot selector
        self.read_instance.cb_position_5 = set_formatting(
            ComboBox(self), self.read_instance.formatting_dict["menu_combobox"]
        )
        self.read_instance.cb_position_5.setToolTip(
            "Select plot type in bottom right position"
        )
        self.read_instance.cb_position_5.currentTextChanged.connect(
            self.read_instance.handle_layout_update
        )

        # MAP SETTINGS MENU #
        # create map settings menu
        self.map_menu = SettingsMenu(plot_type="map", canvas_instance=self)
        self.map_options = self.map_menu.checkable_comboboxes["options"]

        # get stats
        self.map_z_stat = self.map_menu.comboboxes["z_stat"]
        self.map_z1 = self.map_menu.comboboxes["z1"]
        self.map_z2 = self.map_menu.comboboxes["z2"]

        # get colourbar limit fields - populated with the actual resolved
        # limits after each redraw (see update_map_z_statistic()), not a
        # generic "auto" placeholder; no override until the user edits one
        # away from that - see handle_map_cb_limits_update()
        self.map_cb_min = self.map_menu.lineedits["cb_min"]
        self.map_cb_max = self.map_menu.lineedits["cb_max"]
        self.read_instance.map_vmin_override = None
        self.read_instance.map_vmax_override = None

        # MAP "POINTS" SUB-MENU #
        # a separate floating panel for the marker size/opacity controls and
        # automatic sizing, kept out of the main panel so opening Settings
        # isn't a wall of controls
        self.map_sizing_container = self.map_menu.containers["sizing_container"]

        # get sliders and update values
        self.map_markersize_unsel_sl = self.map_menu.sliders["markersize_unsel_sl"]
        self.map_markersize_unsel_sl.setValue(
            int(self.plot_characteristics["map"]["marker_unselected"]["s"])
        )
        self.map_opacity_unsel_sl = self.map_menu.sliders["opacity_unsel_sl"]
        self.map_opacity_unsel_sl.setValue(
            int(self.plot_characteristics["map"]["marker_unselected"]["alpha"] * 10)
        )
        self.map_markersize_sel_sl = self.map_menu.sliders["markersize_sel_sl"]
        self.map_markersize_sel_sl.setValue(
            int(self.plot_characteristics["map"]["marker_selected"]["s"])
        )
        self.map_opacity_sel_sl = self.map_menu.sliders["opacity_sel_sl"]
        self.map_opacity_sel_sl.setValue(
            int(self.plot_characteristics["map"]["marker_selected"]["alpha"] * 10)
        )

        # get map interactive dictionary
        self.interactive_elements["map"] = {
            "hidden": True,
            "markersize_sl": [self.map_markersize_unsel_sl, self.map_markersize_sel_sl],
            "opacity_sl": [self.map_opacity_unsel_sl, self.map_opacity_sel_sl],
        }

        # automatic size/opacity, on by default - keeps points well spaced as
        # the map is zoomed and boosts selected stations over unselected ones
        # (see apply_automatic_marker_style()). The manual sliders stay
        # functional as an override, just hidden while automatic is on
        self.map_auto_sizing = self.map_menu.checkboxes["auto_sizing"]
        automatic = bool(
            self.plot_characteristics["map"].get("marker_automatic", True)
        )
        self.map_auto_sizing.setChecked(automatic)
        self.read_instance.map_auto_marker_sizing = automatic

        self.map_sizing_elements = [
            # container first, so raising this list in order (see
            # handle_map_sizing_toggle()) puts it behind everything else
            # drawn on top of it
            self.map_sizing_container,
            self.map_menu.labels["sizing_title"],
            self.map_auto_sizing,
            self.map_menu.labels["sizing_unsel_label"],
            self.map_menu.labels["markersize_unsel_sl_label"],
            self.map_markersize_unsel_sl,
            self.map_menu.labels["opacity_unsel_sl_label"],
            self.map_opacity_unsel_sl,
            self.map_menu.labels["sizing_sel_label"],
            self.map_menu.labels["markersize_sel_sl_label"],
            self.map_markersize_sel_sl,
            self.map_menu.labels["opacity_sel_sl_label"],
            self.map_opacity_sel_sl,
        ]
        self.sync_map_sizing_sliders_enabled()

        # MAP / COLOURBAR SUB-MENUS #
        # the cosmetic and rarely-touched controls live in their own floating
        # panels rather than the main one. All three share the same on-screen
        # spot, so only one is ever open at a time
        self.map_panel_button = self.map_menu.buttons["map_panel_button"]
        self.map_points_panel_button = self.map_menu.buttons["points_panel_button"]
        self.map_colourbar_panel_button = self.map_menu.buttons[
            "colourbar_panel_button"
        ]
        self.map_panel_container = self.map_menu.containers["map_panel_container"]
        self.map_colourbar_panel_container = self.map_menu.containers[
            "colourbar_panel_container"
        ]

        # get colourmap selector - dashboard only. It opens on whatever the
        # first statistic resolves to rather than a fixed colourmap, and
        # re-resolves whenever the statistic changes; only an explicit choice
        # here overrides that. See sync_map_colourmap()
        self.map_colourmap = self.map_menu.comboboxes["colourmap"]
        populate_colourmap_combobox(
            self.map_colourmap,
            current=get_role_colourmap(
                "sequential",
                self.plot_characteristics_templates["map"].get("colour_preset"),
            ),
        )
        self.read_instance.map_colourmap_override = None

        # get colourmap scale (continuous/discrete) and number-of-chunks
        # controls, preseeded from plot_characteristics.yaml's map.cb.n_discrete.
        # The chunks field only means anything for a discrete scale, so it is
        # hidden entirely for continuous
        self.map_colourmap_scale = self.map_menu.comboboxes["colourmap_scale"]
        self.map_n_sections_label = self.map_menu.labels["n_sections_label"]
        self.map_n_sections = self.map_menu.lineedits["n_sections"]
        self.map_colourmap_scale.addItems(["Continuous", "Discrete"])
        default_n_discrete = self.plot_characteristics_templates["map"]["cb"].get(
            "n_discrete"
        )
        is_discrete = bool(default_n_discrete)
        # setCurrentIndex(), not setCurrentText() - on this editable ComboBox
        # the latter updates the line edit's text without moving
        # currentIndex, so the selection silently reverts to index 0 the next
        # time something reads it
        self.map_colourmap_scale.setCurrentIndex(1 if is_discrete else 0)
        # free-text integer fields rather than a fixed dropdown of counts,
        # so any number can be asked for; seeded with the same default the
        # dropdown used to preselect. setText() (not placeholder) so the
        # value is really there to be read back and edited, not just hinted
        self.map_n_sections.setValidator(QtGui.QIntValidator(1, 256, self.map_n_sections))
        if is_discrete:
            self.map_n_sections.setText(str(int(default_n_discrete)))
        self.read_instance.map_discrete_override = is_discrete
        self.read_instance.map_n_discrete_override = (
            int(default_n_discrete) if is_discrete else None
        )

        # get labels (tick count) control - dashboard-only, same as above.
        # Seeded from plot_characteristics.yaml's map.cb.n_ticks (7 by
        # default), the value every basic/bias statistic falls back to
        # unless it defines its own.
        self.map_n_ticks = self.map_menu.lineedits["n_ticks"]
        default_n_ticks = self.plot_characteristics_templates["map"]["cb"].get(
            "n_ticks"
        )
        self.map_n_ticks.setValidator(QtGui.QIntValidator(1, 256, self.map_n_ticks))
        if default_n_ticks:
            self.map_n_ticks.setText(str(int(default_n_ticks)))
        self.read_instance.map_n_ticks_override = (
            int(default_n_ticks) if default_n_ticks else None
        )

        # colourbar limit fields (moved here from the main panel) plus the
        # small reset control that clears both back to automatic
        self.map_cb_reset = self.map_menu.buttons["cb_reset"]

        # get map projection selector, preselected to whatever
        # plot_characteristics.yaml's map.projection currently resolves to
        # (Robinson by default)
        self.map_projection = self.map_menu.comboboxes["projection"]
        populate_projection_combobox(
            self.map_projection,
            current=self.plot_characteristics_templates["map"]["projection"],
        )

        # get map detail controls (land/ocean colour, country borders,
        # map resolution) - all read from plot_characteristics_templates
        # (see draw_map_features()), same as projection above
        # seeded from the resolved colours, not the raw ones - land/ocean are
        # left empty in the config when they come from the preset, and an empty
        # value is not one of the options, so the box would fall back to its
        # first entry and name a colour the map is not drawn in
        map_template = self.plot_characteristics_templates["map"]
        land_colour, ocean_colour = get_map_colours(map_template)
        self.map_land_colour = self.map_menu.comboboxes["land_colour"]
        populate_colour_combobox(
            self.map_land_colour, LAND_COLOUR_OPTIONS, current=land_colour
        )
        self.map_ocean_colour = self.map_menu.comboboxes["ocean_colour"]
        populate_colour_combobox(
            self.map_ocean_colour, OCEAN_COLOUR_OPTIONS, current=ocean_colour
        )
        # colour preset selector, sitting above the individual land/ocean
        # controls it drives - see handle_map_colour_preset_update()
        self.map_colour_preset = self.map_menu.comboboxes["colour_preset"]
        self.map_colour_preset.addItems(
            list(get_colour_presets()) + [COLOUR_PRESET_CUSTOM]
        )
        self.sync_map_colour_preset()

        self.map_borders = self.map_menu.checkboxes["borders"]
        self.map_borders.setChecked(map_template["borders"]["visible"])

        # gridlines read from plot_characteristics (not _templates) - see
        # draw_map_gridlines() for why this differs from the four above
        self.map_gridlines = self.map_menu.checkboxes["gridlines"]
        self.map_gridlines.setChecked(
            self.plot_characteristics["map"]["gridlines"]["visible"]
        )
        self.map_resolution = self.map_menu.comboboxes["resolution"]
        resolution_options = ["low", "medium", "high"]
        self.map_resolution.addItems(resolution_options)
        # setCurrentIndex(), not setCurrentText() - see the comment above
        # on map_colourmap_scale
        self.map_resolution.setCurrentIndex(
            resolution_options.index(map_template["map_resolution"])
        )

        # elements belonging to each sub-panel, listed separately so
        # handle_map_subpanel_toggle() can show/hide just one panel's worth,
        # but also included in self.map_elements below - dashboard.py only
        # repositions what is in that list, and anything left out never moves
        # again from its initial position. reset_map_subpanels() re-hides
        # these right after the settings button shows map_elements wholesale.
        # Container first in each list, so raising in order puts it behind
        self.map_panel_elements = [
            self.map_panel_container,
            self.map_menu.labels["map_panel_title"],
            self.map_menu.labels["projection_label"],
            self.map_projection,
            self.map_menu.labels["colour_preset_label"],
            self.map_colour_preset,
            self.map_menu.labels["land_colour_label"],
            self.map_land_colour,
            self.map_menu.labels["ocean_colour_label"],
            self.map_ocean_colour,
            self.map_menu.labels["resolution_label"],
            self.map_resolution,
            self.map_borders,
            self.map_gridlines,
        ]
        self.map_colourbar_panel_elements = [
            self.map_colourbar_panel_container,
            self.map_menu.labels["colourbar_panel_title"],
            self.map_menu.labels["cb_limits_label"],
            self.map_cb_min,
            self.map_cb_max,
            self.map_cb_reset,
            self.map_menu.labels["colourmap_label"],
            self.map_colourmap,
            self.map_menu.labels["colourmap_scale_label"],
            self.map_colourmap_scale,
            self.map_menu.labels["n_ticks_label"],
            self.map_n_ticks,
            self.map_n_sections_label,
            self.map_n_sections,
        ]

        # each sub-panel's own reset control, sitting beside its title -
        # added to that panel's element list so it shows and hides with it
        self.map_panel_reset = self.map_menu.buttons["map_panel_reset"]
        self.map_colourbar_panel_reset = self.map_menu.buttons[
            "colourbar_panel_reset"
        ]
        self.map_points_panel_reset = self.map_menu.buttons["sizing_reset"]
        self.map_panel_elements.append(self.map_panel_reset)
        self.map_colourbar_panel_elements.append(self.map_colourbar_panel_reset)
        self.map_sizing_elements.append(self.map_points_panel_reset)

        # the state every panel's reset control returns to: whatever each
        # control holds at the end of construction, which is what the
        # dashboard opens with. Captured rather than re-derived from the yaml
        # at reset time, so the two can't drift apart
        self.map_panel_startup_defaults = {
            "projection": self.map_projection.currentIndex(),
            "land_colour": self.map_land_colour.currentIndex(),
            "ocean_colour": self.map_ocean_colour.currentIndex(),
            "resolution": self.map_resolution.currentIndex(),
            "borders": self.map_borders.isChecked(),
            "gridlines": self.map_gridlines.isChecked(),
        }
        # the colourmap is not captured here: it has no fixed startup value,
        # being resolved from whichever statistic is on screen, so both resets
        # clear the override and let it resolve again
        self.map_colourbar_startup_defaults = {
            "colourmap_scale": self.map_colourmap_scale.currentIndex(),
            "n_ticks": self.map_n_ticks.text(),
            "n_sections": self.map_n_sections.text(),
        }
        self.map_points_startup_defaults = {
            "auto_sizing": self.map_auto_sizing.isChecked(),
            "markersize_unsel": self.map_markersize_unsel_sl.value(),
            "opacity_unsel": self.map_opacity_unsel_sl.value(),
            "markersize_sel": self.map_markersize_sel_sl.value(),
            "opacity_sel": self.map_opacity_sel_sl.value(),
        }

        # every sub-panel, keyed by the nav button that opens it, so the
        # three toggles/resets are one shared implementation rather than
        # three near-identical copies
        self.map_subpanels = {
            "map": (self.map_panel_button, self.map_panel_elements, "Map"),
            "points": (self.map_points_panel_button, self.map_sizing_elements, "Points"),
            "colourbar": (
                self.map_colourbar_panel_button,
                self.map_colourbar_panel_elements,
                "Colourbar",
            ),
        }
        self.map_open_subpanel = None
        # deliberately not calling sync_map_n_sections_visibility() here -
        # sub-panel elements all start hidden, and it calls setVisible(True)
        # for a Discrete scale (the config default), which would show
        # n_sections on startup before its panel is ever opened

        # size each nav button to its own text and centre the row
        self.layout_map_nav_buttons()

        # everything built for "map", main panel and sub-panels alike - see
        # above for why the sub-panels stay in this list. get_elements() never
        # included buttons, so the nav buttons and the resets are added
        # explicitly to be repositioned and hidden with the rest of the menu
        self.map_elements = self.map_menu.get_elements() + [
            self.map_panel_button,
            self.map_points_panel_button,
            self.map_colourbar_panel_button,
            self.map_cb_reset,
            self.map_panel_reset,
            self.map_colourbar_panel_reset,
            self.map_points_panel_reset,
        ]

        # hide every sub-panel whenever Settings is opened or closed, so none
        # is left open, or visible with the rest of the menu hidden, from a
        # previous session
        self.map_menu.buttons["settings_button"].clicked.connect(
            self.reset_map_subpanels
        )

        # TIMESERIES PLOT SETTINGS MENU #
        # create timeseries settings menu
        self.timeseries_menu = SettingsMenu(
            plot_type="timeseries", canvas_instance=self
        )
        self.timeseries_options = self.timeseries_menu.checkable_comboboxes["options"]
        self.timeseries_elements = self.timeseries_menu.get_elements()

        # get aggregation stat, chunk stat and chunk resolution
        self.timeseries_stat = self.timeseries_menu.comboboxes["stat"]
        self.timeseries_chunk_stat = self.timeseries_menu.comboboxes["chunk_stat"]
        self.timeseries_chunk_resolution = self.timeseries_menu.comboboxes[
            "chunk_resolution"
        ]

        # get sliders and update values
        self.timeseries_markersize_sl = self.timeseries_menu.sliders["markersize_sl"]
        self.timeseries_markersize_sl.setMaximum(
            int(self.plot_characteristics["timeseries"]["plot"]["markersize"] * 10)
        )
        self.timeseries_markersize_sl.setValue(
            int(self.plot_characteristics["timeseries"]["plot"]["markersize"])
        )
        self.timeseries_smooth_linewidth_sl = self.timeseries_menu.sliders[
            "smooth_linewidth_sl"
        ]
        self.timeseries_smooth_linewidth_sl.setMaximum(
            int(
                self.plot_characteristics["timeseries"]["smooth"]["format"]["linewidth"]
                * 100
            )
        )
        self.timeseries_smooth_linewidth_sl.setValue(
            int(
                self.plot_characteristics["timeseries"]["smooth"]["format"]["linewidth"]
                * 10
            )
        )
        self.timeseries_smooth_window_sl = self.timeseries_menu.sliders[
            "smooth_window_sl"
        ]
        # the smooth line starts off, its window at zero, and is turned on by
        # setting a window (see update_smooth_window()) - the same way the
        # scatter plot's regression line is by its width. Held to what the
        # slider shows, as the window the plot characteristics carry is what a
        # smooth line is drawn with once one is asked for, not a line already
        # on screen
        self.timeseries_smooth_window_sl.setValue(0)
        self.plot_characteristics["timeseries"]["smooth"]["window"] = 0
        self.plot_characteristics_templates["timeseries"]["smooth"]["window"] = 0
        self.timeseries_smooth_min_points_sl = self.timeseries_menu.sliders[
            "smooth_min_points_sl"
        ]
        self.timeseries_smooth_min_points_sl.setValue(
            int(
                self.plot_characteristics["timeseries"]["smooth"][
                    "min_points_percentage"
                ]
            )
        )

        # get timeseries interactive dictionary
        self.interactive_elements["timeseries"] = {
            "hidden": True,
            "markersize_sl": [self.timeseries_markersize_sl],
            "linewidth_sl": [self.timeseries_smooth_linewidth_sl],
            "smooth_window_sl": [self.timeseries_smooth_window_sl],
            "smooth_min_points_sl": [self.timeseries_smooth_min_points_sl],
        }

        # PERIODIC PLOT SETTINGS MENU #
        # create periodic settings menu
        self.periodic_menu = SettingsMenu(plot_type="periodic", canvas_instance=self)
        self.periodic_options = self.periodic_menu.checkable_comboboxes["options"]
        self.periodic_elements = self.periodic_menu.get_elements()

        # get stats
        self.periodic_stat = self.periodic_menu.comboboxes["stat"]

        # get sliders and update values
        self.periodic_markersize_sl = self.periodic_menu.sliders["markersize_sl"]
        self.periodic_markersize_sl.setMaximum(
            int(self.plot_characteristics["periodic"]["plot"]["markersize"] * 10)
        )
        self.periodic_markersize_sl.setValue(
            int(self.plot_characteristics["periodic"]["plot"]["markersize"])
        )
        self.periodic_linewidth_sl = self.periodic_menu.sliders["linewidth_sl"]
        self.periodic_linewidth_sl.setMaximum(
            int(self.plot_characteristics["periodic"]["plot"]["linewidth"] * 100)
        )
        self.periodic_linewidth_sl.setValue(
            int(self.plot_characteristics["periodic"]["plot"]["linewidth"] * 10)
        )

        # get periodic interactive dictionary
        self.interactive_elements["periodic"] = {
            "hidden": True,
            "markersize_sl": [self.periodic_markersize_sl],
            "linewidth_sl": [self.periodic_linewidth_sl],
        }

        # PERIODIC VIOLIN PLOT SETTINGS MENU #
        # create periodic violin settings menu
        self.periodic_violin_menu = SettingsMenu(
            plot_type="periodic_violin", canvas_instance=self
        )
        self.periodic_violin_options = self.periodic_violin_menu.checkable_comboboxes[
            "options"
        ]
        self.periodic_violin_elements = self.periodic_violin_menu.get_elements()

        # get sliders and update values
        self.periodic_violin_markersize_sl = self.periodic_violin_menu.sliders[
            "markersize_sl"
        ]
        self.periodic_violin_markersize_sl.setMaximum(
            int(
                self.plot_characteristics["periodic-violin"]["plot"]["median"][
                    "markersize"
                ]
                * 10
            )
        )
        self.periodic_violin_markersize_sl.setValue(
            int(
                self.plot_characteristics["periodic-violin"]["plot"]["median"][
                    "markersize"
                ]
            )
        )
        self.periodic_violin_linewidth_sl = self.periodic_violin_menu.sliders[
            "linewidth_sl"
        ]
        self.periodic_violin_linewidth_sl.setMaximum(
            int(
                self.plot_characteristics["periodic-violin"]["plot"]["median"][
                    "linewidth"
                ]
                * 100
            )
        )
        self.periodic_violin_linewidth_sl.setValue(
            int(
                self.plot_characteristics["periodic-violin"]["plot"]["median"][
                    "linewidth"
                ]
                * 10
            )
        )

        # get periodic violin interactive dictionary
        self.interactive_elements["periodic_violin"] = {
            "hidden": True,
            "markersize_sl": [self.periodic_violin_markersize_sl],
            "linewidth_sl": [self.periodic_violin_linewidth_sl],
        }

        # METADATA PLOT SETTINGS MENU #
        # create metadata settings menu
        self.metadata_menu = SettingsMenu(plot_type="metadata", canvas_instance=self)
        self.metadata_elements = self.metadata_menu.get_elements()

        # get metadata interactive dictionary
        self.interactive_elements["metadata"] = {"hidden": True}

        # DISTRIBUTION PLOT SETTINGS MENU #
        # create distribution settings menu
        self.distribution_menu = SettingsMenu(
            plot_type="distribution", canvas_instance=self
        )
        self.distribution_options = self.distribution_menu.checkable_comboboxes[
            "options"
        ]
        self.distribution_elements = self.distribution_menu.get_elements()

        # get sliders and update values
        self.distribution_linewidth_sl = self.distribution_menu.sliders["linewidth_sl"]
        self.distribution_linewidth_sl.setMaximum(
            int(self.plot_characteristics["distribution"]["plot"]["linewidth"] * 100)
        )
        self.distribution_linewidth_sl.setValue(
            self.plot_characteristics["distribution"]["plot"]["linewidth"] * 10
        )

        # get distribution interactive dictionary
        self.interactive_elements["distribution"] = {
            "hidden": True,
            "linewidth_sl": [self.distribution_linewidth_sl],
        }

        # HISTOGRAM PLOT SETTINGS MENU #
        # create histogram settings menu
        self.histogram_menu = SettingsMenu(plot_type="histogram", canvas_instance=self)
        self.histogram_options = self.histogram_menu.checkable_comboboxes["options"]
        self.histogram_elements = self.histogram_menu.get_elements()

        # get sliders and update values
        self.histogram_linewidth_sl = self.histogram_menu.sliders["linewidth_sl"]
        self.histogram_linewidth_sl.setMaximum(
            int(self.plot_characteristics["histogram"]["plot"]["linewidth"] * 100)
        )
        self.histogram_linewidth_sl.setValue(
            self.plot_characteristics["histogram"]["plot"]["linewidth"] * 10
        )

        # the bin count slider spans the range the automatic count is held to,
        # and is set to whatever count is drawn until it is moved - see
        # sync_histogram_bins_slider()
        self.histogram_bins_sl = self.histogram_menu.sliders["bins_sl"]
        self.histogram_bins_sl.setMinimum(
            self.plot_characteristics["histogram"]["min_bins"]
        )
        self.histogram_bins_sl.setMaximum(
            self.plot_characteristics["histogram"]["max_bins"]
        )

        # automatic bin count, on by default - worked out from the data every
        # time the histogram is drawn (see get_histogram_bin_edges()). The
        # slider stays showing the count in use while automatic is on, just
        # disabled, so taking manual control starts from what is on screen
        self.histogram_auto_bins = self.histogram_menu.checkboxes["auto_bins"]
        self.histogram_auto_bins.setChecked(
            self.plot_characteristics["histogram"]["bins"] == "auto"
        )
        self.sync_histogram_bins_slider_enabled()

        # get histogram interactive dictionary
        self.interactive_elements["histogram"] = {
            "hidden": True,
            "linewidth_sl": [self.histogram_linewidth_sl],
            "auto_bins": [self.histogram_auto_bins],
            "bins_sl": [self.histogram_bins_sl],
        }

        # SCATTER PLOT SETTINGS MENU #
        # create scatter settings menu
        self.scatter_menu = SettingsMenu(plot_type="scatter", canvas_instance=self)
        self.scatter_options = self.scatter_menu.checkable_comboboxes["options"]
        self.scatter_elements = self.scatter_menu.get_elements()

        # get sliders and update values
        self.scatter_markersize_sl = self.scatter_menu.sliders["markersize_sl"]
        self.scatter_markersize_sl.setMaximum(
            int(self.plot_characteristics["scatter"]["plot"]["markersize"] * 10)
        )
        self.scatter_markersize_sl.setValue(
            int(self.plot_characteristics["scatter"]["plot"]["markersize"])
        )
        self.scatter_regression_linewidth_sl = self.scatter_menu.sliders[
            "regression_linewidth_sl"
        ]
        self.scatter_regression_linewidth_sl.setMaximum(
            int(self.plot_characteristics["scatter"]["regression"]["linewidth"] * 100)
        )
        # the regression line starts off, its width at zero, and is turned on by
        # setting a width - as the timeseries smooth line is by its window (see
        # update_regression_linewidth())
        self.scatter_regression_linewidth_sl.setValue(0)
        self.plot_characteristics["scatter"]["regression"]["linewidth"] = 0.0
        self.plot_characteristics_templates["scatter"]["regression"]["linewidth"] = 0.0

        # get scatter interactive dictionary
        self.interactive_elements["scatter"] = {
            "hidden": True,
            "markersize_sl": [self.scatter_markersize_sl],
            "linewidth_sl": [self.scatter_regression_linewidth_sl],
        }

        # FAIRMODE TARGET PLOT SETTINGS MENU #
        # create fairmode target settings menu
        self.fairmode_target_menu = SettingsMenu(
            plot_type="fairmode_target", canvas_instance=self
        )
        self.fairmode_target_options = self.fairmode_target_menu.checkable_comboboxes[
            "options"
        ]
        self.fairmode_target_elements = self.fairmode_target_menu.get_elements()
        self.fairmode_target_classification = self.fairmode_target_menu.comboboxes[
            "classification"
        ]
        self.fairmode_target_classification.addItems(["Area", "Station"])

        # get sliders and update values
        self.fairmode_target_markersize_sl = self.fairmode_target_menu.sliders[
            "markersize_sl"
        ]
        self.fairmode_target_markersize_sl.setMaximum(
            int(self.plot_characteristics["fairmode-target"]["plot"]["markersize"] * 10)
        )
        self.fairmode_target_markersize_sl.setValue(
            int(self.plot_characteristics["fairmode-target"]["plot"]["markersize"])
        )

        # get fairmode target interactive dictionary
        self.interactive_elements["fairmode_target"] = {
            "hidden": True,
            "markersize_sl": [self.fairmode_target_markersize_sl],
        }

        # FAIRMODE STATSUMMARY PLOT SETTINGS MENU #
        # create fairmode statsummary settings menu
        self.fairmode_statsummary_menu = SettingsMenu(
            plot_type="fairmode_statsummary", canvas_instance=self
        )
        self.fairmode_statsummary_elements = (
            self.fairmode_statsummary_menu.get_elements()
        )

        # get sliders and update values
        self.fairmode_statsummary_markersize_sl = (
            self.fairmode_statsummary_menu.sliders["markersize_sl"]
        )
        self.fairmode_statsummary_markersize_sl.setMaximum(
            int(self.plot_characteristics["fairmode-target"]["plot"]["markersize"] * 2)
        )
        self.fairmode_statsummary_markersize_sl.setValue(
            int(self.plot_characteristics["fairmode-target"]["plot"]["markersize"])
        )

        # get fairmode statsummary interactive dictionary
        self.interactive_elements["fairmode_statsummary"] = {
            "hidden": True,
            "markersize_sl": [self.fairmode_statsummary_markersize_sl],
        }

        # STATSUMMARY PLOT SETTINGS MENU #
        # create statsummary settings menu
        self.statsummary_menu = SettingsMenu(
            plot_type="statsummary", canvas_instance=self
        )
        self.statsummary_options = self.statsummary_menu.checkable_comboboxes["options"]
        self.statsummary_elements = self.statsummary_menu.get_elements()

        # get stats and add items to cycle
        self.statsummary_stat = self.statsummary_menu.checkable_comboboxes["stat"]
        self.statsummary_cycle = self.statsummary_menu.comboboxes["cycle"]
        self.statsummary_cycle.addItems(["None", "Diurnal", "Weekly", "Monthly"])
        self.statsummary_periodic_mode = self.statsummary_menu.comboboxes[
            "periodic_mode"
        ]
        self.statsummary_periodic_aggregation = self.statsummary_menu.comboboxes[
            "periodic_aggregation"
        ]

        # get statsummary interactive dictionary
        self.interactive_elements["statsummary"] = {"hidden": True}

        # BOXPLOT PLOT SETTINGS MENU #
        # create boxplot settings menu
        self.boxplot_menu = SettingsMenu(plot_type="boxplot", canvas_instance=self)
        self.boxplot_options = self.boxplot_menu.checkable_comboboxes["options"]
        self.boxplot_elements = self.boxplot_menu.get_elements()

        # whether the category labels fit along the x-axis (horizontal, or
        # rotated) is worked out fresh every time the boxplot is drawn (see
        # Plotting.make_boxplot() / fit_boxplot_xticklabels()); this
        # checkbox just reports what was decided, and lets that decision be
        # overridden for what is currently on screen without redoing the
        # whole plot (see handle_boxplot_xlabels_update())
        self.boxplot_xlabels = self.boxplot_menu.checkboxes["xlabels"]
        self.boxplot_xtick_cache = None

        # get boxplot interactive dictionary
        self.interactive_elements["boxplot"] = {"hidden": True}

        # TAYLOR DIAGRAM SETTINGS MENU #
        # create taylor diagram settings menu
        self.taylor_menu = SettingsMenu(plot_type="taylor", canvas_instance=self)
        self.taylor_options = self.taylor_menu.checkable_comboboxes["options"]
        self.taylor_elements = self.taylor_menu.get_elements()

        # get stat
        self.taylor_corr_stat = self.taylor_menu.comboboxes["corr_stat"]

        # get sliders and update values
        self.taylor_markersize_sl = self.taylor_menu.sliders["markersize_sl"]
        self.taylor_markersize_sl.setMaximum(
            int(self.plot_characteristics["taylor"]["plot"]["markersize"] * 10)
        )
        self.taylor_markersize_sl.setValue(
            int(self.plot_characteristics["taylor"]["plot"]["markersize"])
        )

        # get statsummary interactive dictionary
        self.interactive_elements["taylor"] = {
            "hidden": True,
            "markersize_sl": [self.taylor_markersize_sl],
        }

        # CONTINGENCY TABLE SETTINGS MENU #
        # create contingency table settings menu
        self.contingencytable_menu = SettingsMenu(
            plot_type="contingencytable", canvas_instance=self
        )
        self.contingencytable_options = self.contingencytable_menu.checkable_comboboxes[
            "options"
        ]
        self.contingencytable_elements = self.contingencytable_menu.get_elements()

        # get contingency table interactive dictionary
        self.interactive_elements["contingencytable"] = {"hidden": True}

        # create array with buttons and elements to edit when the canvas is resized or the plots are changed
        self.menu_buttons = []
        self.save_buttons = []
        self.save_data_buttons = []
        self.elements = []
        for plot_type in settings_dict.keys():
            if plot_type in [
                "periodic-violin",
                "fairmode-target",
                "fairmode-statsummary",
            ]:
                plot_type = plot_type.replace("-", "_")

            self.menu_buttons.append(
                getattr(self, plot_type + "_menu").buttons["settings_button"]
            )
            self.save_buttons.append(
                getattr(self, plot_type + "_menu").buttons["save_button"]
            )
            self.save_data_buttons.append(
                getattr(self, plot_type + "_menu").buttons["save_data_button"]
            )
            self.elements.append(getattr(self, plot_type + "_elements"))

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
        """
        Function to show and hide elements in setting menus
        """

        event_source = self.sender()
        for key, val in self.interactive_elements.items():
            if event_source == getattr(self, key + "_menu").buttons["settings_button"]:
                hidden = self.interactive_elements[key]["hidden"]
                elements = getattr(self, key + "_elements")
                break

        if hidden:
            # raised as well as shown, so the menu just opened is drawn over any
            # other left open that it overlaps, rather than wherever the menus'
            # fixed order puts it. In the order the elements are listed, which
            # puts the menu's own panel beneath its controls
            for element in elements:
                if isinstance(element, dict):
                    for sub_element in element.values():
                        sub_element.show()
                        sub_element.raise_()
                else:
                    element.show()
                    element.raise_()

            self.interactive_elements[key]["hidden"] = False
        else:
            self.close_settings_menus([key])

        return None

    def close_settings_menus(self, keys=None):
        """
        Function which closes settings menus that are open - e.g. those of the
        plots covered while their data is re-read, or while no stations are
        selected, which would otherwise be left open over plots no longer
        drawn.

        Parameters
        ----------
        keys : list, optional
            Menus to close, as keys of self.interactive_elements (default is
            None, i.e. every menu but the map's, whose plot is never covered
            on its own)
        """

        if keys is None:
            keys = [key for key in self.interactive_elements if key != "map"]

        for key in keys:
            if self.interactive_elements[key]["hidden"]:
                continue
            for element in getattr(self, key + "_elements"):
                if isinstance(element, dict):
                    for sub_element in element.values():
                        sub_element.hide()
                else:
                    element.hide()
            self.interactive_elements[key]["hidden"] = True

        return None

    def cover_plot_axes(self):
        """
        Function which covers every plot but the map's, closing their settings
        menus - see close_settings_menus().
        """

        self.close_settings_menus()
        self.top_right_canvas_cover.show()
        self.lower_canvas_cover.show()

        return None

    def update_markersize_func(self):
        """
        Function to handle the update of the markers size
        """

        event_source = self.sender()
        source_object = event_source.objectName()
        if "_unsel" in source_object:
            loc = 0
        elif "_sel" in source_object:
            loc = 1
        else:
            loc = 0

        for key, val in self.interactive_elements.items():
            if "markersize_sl" in self.interactive_elements[key]:
                if event_source in self.interactive_elements[key]["markersize_sl"]:
                    markersize = self.interactive_elements[key]["markersize_sl"][
                        loc
                    ].value()
                    break

        # correct perodic-violin and fairmode plots names
        if key in ["periodic_violin", "fairmode_target", "fairmode_statsummary"]:
            key = key.replace("_", "-")

        # busy cursor around the redraw, same as every other settings
        # handler - re-styling every point on a densely populated map is
        # not instant, and without this the slider was one of the few
        # controls that left the plain arrow cursor up while it worked
        self.read_instance.cursor_function = set_cursor(
            self.read_instance.cursor_function, "update_markersize_func"
        )
        # see handle_map_colourmap_scale_update() for why this is here
        QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)

        self.update_markersize(self.plot_axes[key], key, markersize, event_source)

        # repaint synchronously, inside the busy-cursor scope: the update
        # above only schedules a deferred draw, which would otherwise land
        # after the cursor had been restored and show the plain pointer
        # (or matplotlib's own generic wait cursor) for the slow part
        self.figure.canvas.draw()

        unset_cursor(self.read_instance.cursor_function, "update_markersize_func")

        return None

    def update_opacity_func(self):
        """
        Function to handle the update of the markers opacity
        """

        event_source = self.sender()
        source_object = event_source.objectName()
        if "_unsel" in source_object:
            loc = 0
        elif "_sel" in source_object:
            loc = 1
        else:
            loc = 0
        for key, val in self.interactive_elements.items():
            if "opacity_sl" in self.interactive_elements[key]:
                if event_source in self.interactive_elements[key]["opacity_sl"]:
                    opacity = (
                        self.interactive_elements[key]["opacity_sl"][loc].value() / 10
                    )
                    break
        # see update_markersize_func() for why the busy cursor is here
        self.read_instance.cursor_function = set_cursor(
            self.read_instance.cursor_function, "update_opacity_func"
        )
        QtCore.QCoreApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)

        self.update_opacity(self.plot_axes[key], key, opacity, event_source)

        # repaint synchronously, inside the busy-cursor scope: the update
        # above only schedules a deferred draw, which would otherwise land
        # after the cursor had been restored and show the plain pointer
        # (or matplotlib's own generic wait cursor) for the slow part
        self.figure.canvas.draw()

        unset_cursor(self.read_instance.cursor_function, "update_opacity_func")

        return None

    def update_linewidth_func(self):
        """
        Function to handle the update of the lines widths
        """

        # get source
        event_source = self.sender()
        for key, val in self.interactive_elements.items():
            if "linewidth_sl" in self.interactive_elements[key]:
                if event_source in self.interactive_elements[key]["linewidth_sl"]:
                    linewidth = (
                        self.interactive_elements[key]["linewidth_sl"][0].value() / 10
                    )
                    break

        # correct perodic-violin name
        if key == "periodic_violin":
            key = "periodic-violin"

        # the scatter plot's slider sets the width of its regression line,
        # which it also turns on and off - see update_regression_linewidth()
        if key == "scatter":
            self.update_regression_linewidth(linewidth)
            return None

        self.update_linewidth(self.plot_axes[key], key, linewidth)

        return None

    def update_regression_linewidth(self, linewidth):
        """
        Function to handle the update of the scatter plot's regression line
        width, which turns the line on and off with it: there is nothing to see
        at a width of zero, so the regression plot option follows the slider
        the same way the timeseries smooth line follows its window (see
        update_smooth_window()).

        Parameters
        ----------
        linewidth : float
            Width to draw the regression line with
        """

        # written to the template as well as the copy being drawn from, as the
        # two are separate dicts (see Plotting.set_plot_characteristics())
        self.plot_characteristics["scatter"]["regression"]["linewidth"] = linewidth
        self.plot_characteristics_templates["scatter"]["regression"][
            "linewidth"
        ] = linewidth

        # get index of regression in plot options
        all_plot_options = self.plot_characteristics["scatter"]["plot_options"]
        index = all_plot_options.index("regression")

        # remove regression plot option, so the line is drawn again at the new
        # width rather than left as it was
        self.scatter_options.model().item(index).setCheckState(QtCore.Qt.Unchecked)

        # create regression line
        if linewidth > 0:
            # add regression plot option
            self.scatter_options.model().item(index).setCheckState(QtCore.Qt.Checked)

        # draw changes
        self.figure.canvas.draw_idle()

        return None

    def sync_histogram_bins_slider_enabled(self):
        """
        Function which enables the bin count slider only while the count is
        being set by hand, as an automatic count overwrites whatever it is set
        to every time the histogram is drawn.

        The slider stays in place rather than being hidden, showing the count
        the automatic rule arrived at, so that taking manual control has an
        obvious starting point.
        """

        set_slider_enabled(
            self.histogram_bins_sl, not self.histogram_auto_bins.isChecked()
        )

        return None

    def handle_histogram_auto_bins_update(self):
        """
        Function which handles toggling the automatic bin count upon
        interaction with the histogram settings menu's checkbox
        """

        if not self.read_instance.block_config_bar_handling_updates:
            # going back to automatic hands the count back to the data;
            # coming off it holds the histogram at the count on screen
            if self.histogram_auto_bins.isChecked():
                n_bins = "auto"
            else:
                n_bins = self.histogram_bins_sl.value()

            self.sync_histogram_bins_slider_enabled()
            self.set_histogram_bins(n_bins)

        return None

    def set_histogram_bins(self, n_bins):
        """
        Function which sets the number of bins the histogram is drawn with and
        remakes it.

        Parameters
        ----------
        n_bins : int or str
            Number of bins, or "auto" to work it out from the data
        """

        # written to the template as well as the copy being drawn from, as the
        # two are separate dicts (see Plotting.make_plot())
        self.plot_characteristics["histogram"]["bins"] = n_bins
        self.plot_characteristics_templates["histogram"]["bins"] = n_bins

        # remake the plot, as the bins decide what is counted rather than how
        # what has been counted is drawn
        self.update_associated_active_dashboard_plot("histogram")
        self.figure.canvas.draw_idle()

        return None

    def update_histogram_bins_func(self):
        """
        Function to handle the update of the number of histogram bins
        """

        # a number set here stands in for the automatic count until it is
        # changed again, the same way it would if written into the plot
        # characteristics
        n_bins = self.histogram_bins_sl.value()
        self.set_histogram_bins(n_bins)

        return None

    def sync_histogram_bins_slider(self, n_bins):
        """
        Function which points the bin count slider at the number of bins the
        histogram has just been drawn with, so it always reports what is on
        screen rather than a number nothing was drawn with.

        Parameters
        ----------
        n_bins : int
            Number of bins the histogram was drawn with
        """

        # showing the count must not read as setting it, or every redraw would
        # pin the bins to whatever the last one worked out
        self.histogram_bins_sl.blockSignals(True)
        self.histogram_bins_sl.setValue(int(n_bins))
        self.histogram_bins_sl.blockSignals(False)

        return None

    def handle_boxplot_xlabels_update(self):
        """
        Function which handles toggling the boxplot's x-axis labels upon
        interaction with the boxplot settings menu's checkbox.

        Applies the click directly to what is already on screen (showing at
        whatever rotation last fitted, or the steepest tried if none did, or
        hiding outright) rather than remaking the whole plot, as nothing
        about the boxplot itself has changed
        """

        if not self.read_instance.block_config_bar_handling_updates:
            if self.boxplot_xtick_cache is not None:
                self.read_instance.cursor_function = set_cursor(
                    self.read_instance.cursor_function,
                    "handle_boxplot_xlabels_update",
                )
                QtCore.QCoreApplication.processEvents(
                    QtCore.QEventLoop.ExcludeUserInputEvents
                )

                fit_boxplot_xticklabels(
                    self.plot_axes["boxplot"],
                    forced=self.boxplot_xlabels.isChecked(),
                    **self.boxplot_xtick_cache,
                )
                self.figure.canvas.draw()

                unset_cursor(
                    self.read_instance.cursor_function,
                    "handle_boxplot_xlabels_update",
                )

        return None

    def sync_boxplot_xlabels_checkbox(self, shown):
        """
        Function which points the boxplot's "X-axis labels" checkbox at
        whether the labels were actually drawn on this pass, so it always
        reports what is on screen rather than a choice nothing was drawn
        with.

        Parameters
        ----------
        shown : bool
            Whether the boxplot's x-axis category labels are currently shown
        """

        # showing the outcome must not read as the user having set it, or
        # the next click would toggle away from whatever the automatic fit
        # just decided rather than overriding it
        self.boxplot_xlabels.blockSignals(True)
        self.boxplot_xlabels.setChecked(shown)
        self.boxplot_xlabels.blockSignals(False)

        return None

    def update_smooth_window_func(self):
        """
        Function to handle the update of the smooth window
        """

        # get source
        event_source = self.sender()
        plot_type = event_source.objectName().split("_smooth")[0]
        for element in self.interactive_elements[plot_type]["smooth_window_sl"]:
            smooth_window = element.value()
            break

        self.update_smooth_window(plot_type, smooth_window)

        return None

    def update_smooth_min_points_func(self):
        """
        Function to handle the update of the smooth minimum points
        """

        # get source
        event_source = self.sender()
        plot_type = event_source.objectName().split("_smooth")[0]
        for element in self.interactive_elements[plot_type]["smooth_min_points_sl"]:
            smooth_min_points = element.value()
            break

        self.update_smooth_min_points(plot_type, smooth_min_points)

        return None

    @restores_settings_guard
    def update_plot_option(self):
        """
        Function to handle the update of the plot options
        """

        if not self.read_instance.block_MPL_canvas_updates:
            # get source
            event_source = self.sender()
            plot_type_alt = event_source.objectName().split("_options")[0]

            # correct perodic-violin name
            if plot_type_alt in [
                "periodic_violin",
                "fairmode_target",
                "fairmode_statsummary",
            ]:
                plot_type = plot_type_alt.replace("_", "-")
            else:
                plot_type = copy.deepcopy(plot_type_alt)

            # force Taylor diagram to show bias statistics
            if "taylor" in plot_type:
                z_statistic_sign = "bias"
            # define z_statistic_sign to be absolute (will be overwritten if bias)
            else:
                z_statistic_sign = "absolute"

            # an option is selected or there are options in previous to undo?
            if event_source.currentData() or self.previous_plot_options[plot_type]:
                self.current_plot_options[plot_type] = copy.deepcopy(
                    event_source.currentData()
                )
                orig_plot_options = event_source.currentData(all=True)
                mod_plot_options = copy.deepcopy(orig_plot_options)

                # ensure bias option is handled first (to set z_statistic_sign)
                if "bias" in mod_plot_options:
                    mod_plot_options.remove("bias")
                    mod_plot_options.insert(0, "bias")

                # disable bias when threshold is active and viceversa
                if (
                    ("bias" in self.previous_plot_options[plot_type])
                    and ("bias" in self.current_plot_options[plot_type])
                    and ("threshold" in self.current_plot_options[plot_type])
                ):
                    msg = "Bias will be deactivated to show threshold lines"
                    show_message(self.read_instance, msg)
                    bias_index = orig_plot_options.index("bias")
                    self.update_option_on_combobox(event_source, bias_index)
                    self.current_plot_options[plot_type].remove("bias")
                if (
                    ("threshold" in self.previous_plot_options[plot_type])
                    and ("bias" in self.current_plot_options[plot_type])
                    and ("threshold" in self.current_plot_options[plot_type])
                    and (len(self.read_instance.data_labels) > 1)
                ):
                    msg = "Thresholds will be deactivated to show bias plots"
                    show_message(self.read_instance, msg)
                    threshold_index = orig_plot_options.index("threshold")
                    self.update_option_on_combobox(event_source, threshold_index)
                    self.current_plot_options[plot_type].remove("threshold")

                # apply smooth/regression if hide data is checked for timeseries
                if "hidedata" in self.current_plot_options[plot_type]:
                    # get option to check
                    if (plot_type == "timeseries") and (
                        "smooth" not in self.current_plot_options[plot_type]
                    ):
                        self.current_plot_options[plot_type].append("smooth")
                        index_to_check = orig_plot_options.index("smooth")
                    elif (plot_type == "scatter") and (
                        "regression" not in self.current_plot_options[plot_type]
                    ):
                        self.current_plot_options[plot_type].append("regression")
                        index_to_check = orig_plot_options.index("regression")

                    if "index_to_check" in locals():
                        # check option in combobox
                        self.update_option_on_combobox(
                            event_source, index_to_check, uncheck=False
                        )

                        # ensure hidedata option is handled second (to show smooth/regression after hiding data)
                        mod_plot_options.remove("hidedata")
                        mod_plot_options.insert(1, "hidedata")

                for option in mod_plot_options:
                    # get index to raise errors and uncheck options (in original plot options order)
                    index = orig_plot_options.index(option)

                    # if any option other than domain is selected
                    if (
                        len(
                            [
                                option
                                for option in self.current_plot_options[plot_type]
                                if option != "domain"
                            ]
                        )
                        >= 1
                    ):
                        # return if do not have selected station_station_data in memory, then no data plotted yet
                        if not hasattr(self, "selected_station_data"):
                            msg = "Select at least one station in the plot to apply options."
                            show_message(self.read_instance, msg)
                            self.update_option_on_combobox(event_source, index)
                            return None

                        # return from function if selected_station_data has not been updated for new species yet
                        if (
                            self.read_instance.networkspeci
                            not in self.selected_station_data
                        ):
                            return None

                    # undo plot options that were selected before but not now
                    if (option in self.previous_plot_options[plot_type]) and (
                        option not in self.current_plot_options[plot_type]
                    ):
                        undo = True
                    # make plot option if currently selected
                    elif option in self.current_plot_options[plot_type]:
                        undo = False
                    # do nothing if options were never selected
                    elif (option not in self.previous_plot_options[plot_type]) and (
                        option not in self.current_plot_options[plot_type]
                    ):
                        continue

                    # if plot type not in plot_elements, then return
                    if plot_type not in self.plot_elements:
                        return None

                    # if no selected stations then remove all plot_elements for active plot_options,
                    # and then return
                    if len(self.relative_selected_station_inds) == 0:
                        for active_type in self.plot_elements[plot_type]:
                            if active_type != "active":
                                for data_label in self.plot_elements[plot_type][
                                    active_type
                                ]:
                                    for plot_option in self.current_plot_options[
                                        plot_type
                                    ]:
                                        # do not remove domain even if there are no selected stations
                                        if (
                                            plot_option
                                            in self.plot_elements[plot_type][
                                                active_type
                                            ][data_label]
                                        ) and (plot_option != "domain"):
                                            for plot_element in self.plot_elements[
                                                plot_type
                                            ][active_type][data_label][plot_option]:
                                                plot_element.remove()
                                            del self.plot_elements[plot_type][
                                                active_type
                                            ][data_label][plot_option]
                        # do not skip applying domain to map
                        if (plot_type != "map") and (plot_option != "domain"):
                            return None

                    # remove current option elements (both absolute and bias)
                    for active_type in self.plot_elements[plot_type]:
                        if active_type != "active":
                            for data_label in self.plot_elements[plot_type][
                                active_type
                            ]:
                                if (
                                    option
                                    in self.plot_elements[plot_type][active_type][
                                        data_label
                                    ]
                                ):
                                    for plot_element in self.plot_elements[plot_type][
                                        active_type
                                    ][data_label][option]:
                                        plot_element.remove()
                                    del self.plot_elements[plot_type][active_type][
                                        data_label
                                    ][option]

                    # if option is 'bias', then remove all other current option elements (both absolute and bias)
                    if option == "bias":
                        for active_type in self.plot_elements[plot_type]:
                            if active_type != "active":
                                for data_label in self.plot_elements[plot_type][
                                    active_type
                                ]:
                                    for plot_option in self.current_plot_options[
                                        plot_type
                                    ]:
                                        if (
                                            plot_option
                                            in self.plot_elements[plot_type][
                                                active_type
                                            ][data_label]
                                        ):
                                            for plot_element in self.plot_elements[
                                                plot_type
                                            ][active_type][data_label][plot_option]:
                                                plot_element.remove()
                                            del self.plot_elements[plot_type][
                                                active_type
                                            ][data_label][plot_option]

                    # options 'logy' and 'logx'
                    # only plot if axis has all positive values
                    if (option == "logy") or (option == "logx"):
                        if isinstance(self.plot_axes[plot_type], dict):
                            for temporal_resolution, sub_ax in self.plot_axes[
                                plot_type
                            ].items():
                                log_valid = log_validity(sub_ax, option)
                                if log_valid:
                                    log_axes(
                                        sub_ax,
                                        option,
                                        self.plot_characteristics[plot_type],
                                        undo=undo,
                                    )
                                else:
                                    msg = "It is not possible to log the {0}-axis ".format(
                                        option[-1]
                                    )
                                    msg += "in {0} with negative values.".format(
                                        plot_type
                                    )
                                    show_message(self.read_instance, msg)
                                    self.update_option_on_combobox(event_source, index)
                                    return None
                        else:
                            log_valid = log_validity(self.plot_axes[plot_type], option)
                            if log_valid:
                                log_axes(
                                    self.plot_axes[plot_type],
                                    option,
                                    self.plot_characteristics[plot_type],
                                    undo=undo,
                                )
                            else:
                                msg = "It is not possible to log the {0}-axis ".format(
                                    option[-1]
                                )
                                msg += "in {0} with negative values.".format(plot_type)
                                show_message(self.read_instance, msg)
                                self.update_option_on_combobox(event_source, index)
                                return None

                    # option 'annotate'
                    # only plot if have selected stations (for map annotations)
                    elif option == "annotate":
                        if not undo:
                            if isinstance(self.plot_axes[plot_type], dict):
                                for (
                                    relevant_temporal_resolution,
                                    sub_ax,
                                ) in self.plot_axes[plot_type].items():
                                    if (
                                        relevant_temporal_resolution
                                        in self.read_instance.periodic_relevant_temporal_resolutions
                                    ):
                                        annotation(
                                            self.read_instance,
                                            self,
                                            sub_ax,
                                            self.read_instance.networkspeci,
                                            self.read_instance.data_labels,
                                            plot_type,
                                            self.plot_characteristics[plot_type],
                                            self.current_plot_options[plot_type],
                                            plot_z_statistic_sign=z_statistic_sign,
                                        )
                                        break
                            else:
                                annotation(
                                    self.read_instance,
                                    self,
                                    self.plot_axes[plot_type],
                                    self.read_instance.networkspeci,
                                    self.read_instance.data_labels,
                                    plot_type,
                                    self.plot_characteristics[plot_type],
                                    self.current_plot_options[plot_type],
                                    plot_z_statistic_sign=z_statistic_sign,
                                )

                    # option 'smooth'
                    elif option == "smooth":
                        if not undo:
                            # uncheck option in combobox if window is 0
                            if (
                                self.plot_characteristics[plot_type]["smooth"]["window"]
                                == 0
                            ):
                                msg = "It is not possible to show the smooth line "
                                msg += "if window is 0, increase it in advance."
                                show_message(self.read_instance, msg)
                                self.update_option_on_combobox(event_source, index)
                                return None
                            smooth(
                                self.read_instance,
                                self,
                                self.plot_axes[plot_type],
                                self.read_instance.networkspeci,
                                self.read_instance.data_labels,
                                plot_type,
                                self.plot_characteristics[plot_type],
                                self.current_plot_options[plot_type],
                            )

                    # option 'hidedata'
                    elif option == "hidedata":
                        # uncheck option in combobox if smooth window is 0
                        if plot_type == "timeseries":
                            if (
                                self.plot_characteristics[plot_type]["smooth"]["window"]
                                == 0
                            ):
                                msg = "It is not possible to show the smooth line "
                                msg += "if window is 0, increase it in advance."
                                show_message(self.read_instance, msg)
                                self.update_option_on_combobox(event_source, index)
                                # Deactivate also smooth
                                smooth_index = orig_plot_options.index("smooth")
                                self.update_option_on_combobox(
                                    event_source, smooth_index
                                )
                                return None
                        active_type = (
                            "bias"
                            if "bias" in self.current_plot_options[plot_type]
                            else "absolute"
                        )
                        for data_label in self.plot_elements[plot_type][active_type]:
                            if (
                                "plot"
                                in self.plot_elements[plot_type][active_type][
                                    data_label
                                ]
                            ):
                                for element in self.plot_elements[plot_type][
                                    active_type
                                ][data_label]["plot"]:
                                    if not undo:
                                        element.set_visible(False)
                                    else:
                                        element.set_visible(True)

                    # option 'domain'
                    elif option == "domain":
                        if not undo:
                            # plot model grid domain edges on map
                            self.update_model_domain_edges()
                        else:
                            # remove grid domain polygon if previously plotted
                            self.remove_axis_objects(
                                self.plot_axes["map"].patches,
                                types_to_remove=[matplotlib.patches.Polygon],
                            )

                    # option 'regression'
                    elif option == "regression":
                        if not undo:
                            # uncheck option in combobox if line width is 0
                            if (
                                self.plot_characteristics[plot_type]["regression"][
                                    "linewidth"
                                ]
                                == 0
                            ):
                                msg = "It is not possible to show the regression line "
                                msg += "if line width is 0, increase it in advance."
                                show_message(self.read_instance, msg)
                                self.update_option_on_combobox(event_source, index)
                                return None
                            linear_regression(
                                self.read_instance,
                                self,
                                self.plot_axes[plot_type],
                                self.read_instance.networkspeci,
                                self.read_instance.data_labels,
                                plot_type,
                                self.plot_characteristics[plot_type],
                                self.current_plot_options[plot_type],
                            )

                    # option 'gerrity'
                    elif option == "gerrity":
                        # clear all previously plotted artists for plot type
                        self.remove_axis_elements(self.plot_axes[plot_type], plot_type)

                        # make plot again considering plot option
                        func = getattr(self.plotting, "make_contingencytable")
                        func(
                            self.plot_axes[plot_type],
                            self.read_instance.networkspeci,
                            self.read_instance.data_labels,
                            self.plot_characteristics[plot_type],
                            self.current_plot_options[plot_type],
                        )

                    # option 'threshold'
                    elif option == "threshold":
                        if not undo:
                            if isinstance(self.plot_axes[plot_type], dict):
                                for (
                                    relevant_temporal_resolution,
                                    sub_ax,
                                ) in self.plot_axes[plot_type].items():
                                    if (
                                        relevant_temporal_resolution
                                        in self.read_instance.periodic_relevant_temporal_resolutions
                                    ):
                                        threshold(
                                            self.read_instance,
                                            self,
                                            sub_ax,
                                            self.read_instance.networkspeci,
                                            plot_type,
                                            self.plot_characteristics[plot_type],
                                        )
                            else:
                                threshold(
                                    self.read_instance,
                                    self,
                                    self.plot_axes[plot_type],
                                    self.read_instance.networkspeci,
                                    plot_type,
                                    self.plot_characteristics[plot_type],
                                )

                    # option 'bias'
                    elif option == "bias":
                        # firstly if just 1 data label then cannot make bias plot
                        if len(self.read_instance.data_labels) == 1:
                            msg = "It is not possible to make a bias plot with just observations loaded."
                            show_message(self.read_instance, msg)
                            self.update_option_on_combobox(event_source, index)
                            self.plot_elements[plot_type]["active"] = "absolute"
                            self.current_plot_options[plot_type].remove("bias")

                            # create other active plot option elements for now absolute plot (if do not already exist)
                            self.redraw_active_options(
                                self.read_instance.data_labels,
                                plot_type,
                                "absolute",
                                self.current_plot_options[plot_type],
                                z_statistic_sign=z_statistic_sign,
                            )
                            return None

                        # if bias option is enabled then first check if bias elements stored
                        elif not undo:
                            # update active (bias)
                            self.plot_elements[plot_type]["active"] = "bias"

                            # handle some special case for periodic plot
                            if plot_type == "periodic":
                                # get currently selected periodic statistic name
                                base_zstat = self.periodic_stat.currentText()
                                zstat = get_z_statistic_comboboxes(
                                    base_zstat, bias=True
                                )

                                # get zstat information
                                (
                                    zstat,
                                    base_zstat,
                                    z_statistic_type,
                                    z_statistic_sign,
                                    z_statistic_period,
                                ) = get_z_statistic_info(zstat=zstat)

                                # if get_z_statistic_type == 'modbias' then return as bias already plotted
                                if z_statistic_type == "modbias":
                                    self.update_option_on_combobox(event_source, index)
                                    self.plot_elements[plot_type]["active"] = "absolute"
                                    self.current_plot_options[
                                        plot_type
                                    ] = copy.deepcopy(
                                        self.previous_plot_options[plot_type]
                                    )
                                    return None

                            if plot_type == "timeseries":
                                # get currently selected chunk statistic and resolution name
                                chunk_stat = self.timeseries_chunk_stat.currentText()
                                chunk_resolution = (
                                    self.timeseries_chunk_resolution.currentText()
                                )

                                # if get_z_statistic_type == 'modbias' then return as bias already plotted
                                z_statistic_type = get_z_statistic_type(chunk_stat)
                                if z_statistic_type == "modbias":
                                    # chunk timeseries is active?
                                    if (chunk_stat != "None") and (
                                        chunk_resolution != "None"
                                    ):
                                        self.update_option_on_combobox(
                                            event_source, index
                                        )
                                        self.plot_elements[plot_type][
                                            "active"
                                        ] = "absolute"
                                        self.current_plot_options[
                                            plot_type
                                        ] = copy.deepcopy(
                                            self.previous_plot_options[plot_type]
                                        )
                                        return None

                            # iterate through valid data labels
                            bias_labels_to_plot = []
                            for data_label in self.read_instance.data_labels + ["ALL"]:
                                # hide absolute plot elements
                                if (
                                    data_label
                                    in self.plot_elements[plot_type]["absolute"]
                                ):
                                    for element_type in self.plot_elements[plot_type][
                                        "absolute"
                                    ][data_label]:
                                        for element in self.plot_elements[plot_type][
                                            "absolute"
                                        ][data_label][element_type]:
                                            element.set_visible(False)

                                # if have bias elements pre-stored then simply show bias elements (if data label on legend is active)
                                if "bias" in self.plot_elements[plot_type]:
                                    if (
                                        data_label
                                        in self.plot_elements[plot_type]["bias"]
                                    ):
                                        if (
                                            data_label
                                            in self.plot_elements["data_labels_active"]
                                        ) or (data_label == "ALL"):
                                            for element_type in self.plot_elements[
                                                plot_type
                                            ]["bias"][data_label]:
                                                for element in self.plot_elements[
                                                    plot_type
                                                ]["bias"][data_label][element_type]:
                                                    element.set_visible(True)
                                    else:
                                        if plot_type == "statsummary":
                                            if (
                                                data_label
                                                == self.read_instance.observations_data_label
                                            ):
                                                bias_labels_to_plot.append(data_label)
                                        else:
                                            if data_label not in [
                                                self.read_instance.observations_data_label,
                                                "ALL",
                                            ]:
                                                bias_labels_to_plot.append(data_label)
                                else:
                                    if plot_type == "statsummary":
                                        if (
                                            data_label
                                            == self.read_instance.observations_data_label
                                        ):
                                            bias_labels_to_plot.append(data_label)
                                    else:
                                        if data_label not in [
                                            self.read_instance.observations_data_label,
                                            "ALL",
                                        ]:
                                            bias_labels_to_plot.append(data_label)

                            # if do not already have bias elements, then make them (tracking plot elements also)
                            if bias_labels_to_plot:
                                # get plotting function for specific plot
                                if plot_type == "statsummary":
                                    func = getattr(self.plotting, "make_table")
                                elif plot_type in [
                                    "fairmode-target",
                                    "fairmode-statsummary",
                                ]:
                                    func = getattr(
                                        self.plotting,
                                        "make_{}".format(plot_type.replace("-", "_")),
                                    )
                                else:
                                    func = getattr(
                                        self.plotting,
                                        "make_{}".format(plot_type.split("-")[0]),
                                    )

                                # call function to update plot
                                # periodic plot
                                if plot_type == "periodic":
                                    func(
                                        self.plot_axes[plot_type],
                                        self.read_instance.networkspeci,
                                        bias_labels_to_plot,
                                        self.plot_characteristics[plot_type],
                                        self.current_plot_options[plot_type],
                                        zstat=zstat,
                                    )
                                # make statsummary plot
                                elif plot_type == "statsummary":
                                    relevant_zstats = self.active_statsummary_stats[
                                        "modbias"
                                    ]
                                    func(
                                        self.plot_axes[plot_type],
                                        self.read_instance.networkspeci,
                                        self.read_instance.data_labels,
                                        self.plot_characteristics[plot_type],
                                        self.current_plot_options[plot_type],
                                        zstats=relevant_zstats,
                                        statsummary=True,
                                    )
                                # other plots
                                else:
                                    func(
                                        self.plot_axes[plot_type],
                                        self.read_instance.networkspeci,
                                        bias_labels_to_plot,
                                        self.plot_characteristics[plot_type],
                                        self.current_plot_options[plot_type],
                                    )

                            # create other active plot option elements for bias plot (if do not already exist)
                            self.redraw_active_options(
                                self.read_instance.data_labels,
                                plot_type,
                                "bias",
                                self.current_plot_options[plot_type],
                                z_statistic_sign=z_statistic_sign,
                            )

                        # if bias option is not enabled then hide bias plot elements and show absolute plots again
                        else:
                            # update active (absolute)
                            self.plot_elements[plot_type]["active"] = "absolute"

                            # handle some special case for periodic plot
                            if plot_type == "periodic":
                                # get currently selected periodic statistic name
                                base_zstat = self.periodic_stat.currentText()
                                zstat = get_z_statistic_comboboxes(
                                    base_zstat, bias=False
                                )

                                # get zstat information
                                (
                                    zstat,
                                    base_zstat,
                                    z_statistic_type,
                                    z_statistic_sign,
                                    z_statistic_period,
                                ) = get_z_statistic_info(zstat=zstat)

                            # iterate through valid data labels
                            absolute_labels_to_plot = []
                            for data_label in list(self.read_instance.data_labels) + [
                                "ALL"
                            ]:
                                # hide bias plot elements
                                if "bias" in self.plot_elements[plot_type]:
                                    if (
                                        data_label
                                        in self.plot_elements[plot_type]["bias"]
                                    ):
                                        for element_type in self.plot_elements[
                                            plot_type
                                        ]["bias"][data_label]:
                                            for element in self.plot_elements[
                                                plot_type
                                            ]["bias"][data_label][element_type]:
                                                element.set_visible(False)

                                # show absolute plot elements (if data label on legend is active)
                                if (
                                    data_label
                                    in self.plot_elements[plot_type]["absolute"]
                                ):
                                    if (
                                        data_label
                                        in self.plot_elements["data_labels_active"]
                                    ) or (data_label == "ALL"):
                                        for element_type in self.plot_elements[
                                            plot_type
                                        ]["absolute"][data_label]:
                                            for element in self.plot_elements[
                                                plot_type
                                            ]["absolute"][data_label][element_type]:
                                                element.set_visible(True)
                                else:
                                    if plot_type == "statsummary":
                                        if (
                                            data_label
                                            == self.read_instance.observations_data_label
                                        ):
                                            absolute_labels_to_plot.append(data_label)
                                    else:
                                        if data_label != "ALL":
                                            absolute_labels_to_plot.append(data_label)

                            # if do not already have absolute elements, then make them (tracking plot elements also)
                            if absolute_labels_to_plot:
                                # get plotting function for specific plot
                                if plot_type == "statsummary":
                                    func = getattr(self.plotting, "make_table")
                                elif plot_type in [
                                    "fairmode-target",
                                    "fairmode-statsummary",
                                ]:
                                    func = getattr(
                                        self.plotting,
                                        "make_{}".format(plot_type.replace("-", "_")),
                                    )
                                else:
                                    func = getattr(
                                        self.plotting,
                                        "make_{}".format(plot_type.split("-")[0]),
                                    )

                                # call function to update plot
                                # periodic plot
                                if plot_type == "periodic":
                                    func(
                                        self.plot_axes[plot_type],
                                        self.read_instance.networkspeci,
                                        absolute_labels_to_plot,
                                        self.plot_characteristics[plot_type],
                                        self.current_plot_options[plot_type],
                                        zstat=zstat,
                                    )
                                # make statsummary plot
                                elif plot_type == "statsummary":
                                    relevant_zstats = self.active_statsummary_stats[
                                        "basic"
                                    ]
                                    func(
                                        self.plot_axes[plot_type],
                                        self.read_instance.networkspeci,
                                        self.read_instance.data_labels,
                                        self.plot_characteristics[plot_type],
                                        self.current_plot_options[plot_type],
                                        zstats=relevant_zstats,
                                        statsummary=True,
                                    )
                                # other plots
                                else:
                                    func(
                                        self.plot_axes[plot_type],
                                        self.read_instance.networkspeci,
                                        absolute_labels_to_plot,
                                        self.plot_characteristics[plot_type],
                                        self.current_plot_options[plot_type],
                                    )

                            # create other active plot option elements for absolute plot (if do not already exist)
                            self.redraw_active_options(
                                self.read_instance.data_labels,
                                plot_type,
                                "absolute",
                                self.current_plot_options[plot_type],
                                z_statistic_sign=z_statistic_sign,
                            )

                    # check stats for the selected periodic cycle
                    if plot_type in ["statsummary"]:
                        self.read_instance.block_config_bar_handling_updates = True
                        self.check_statsummary_stats()
                        self.read_instance.block_config_bar_handling_updates = False

                    # reset axes limits (harmonising across subplots for periodic plots)
                    if plot_type not in ["map", "taylor", "fairmode-statsummary"]:
                        if plot_type == "scatter":
                            harmonise_xy_lims_paradigm(
                                self.read_instance,
                                self,
                                self.plot_axes[plot_type],
                                plot_type,
                                self.plot_characteristics[plot_type],
                                self.current_plot_options[plot_type],
                                relim=True,
                            )
                        else:
                            harmonise_xy_lims_paradigm(
                                self.read_instance,
                                self,
                                self.plot_axes[plot_type],
                                plot_type,
                                self.plot_characteristics[plot_type],
                                self.current_plot_options[plot_type],
                                relim=True,
                                autoscale=True,
                            )

                # save current plot options as previous
                self.previous_plot_options[plot_type] = self.current_plot_options[
                    plot_type
                ]

                # draw changes
                self.figure.canvas.draw_idle()

        return None

    def redraw_active_options(
        self, data_labels, plot_type, active, plot_options, z_statistic_sign="absolute"
    ):
        """
        Redraw active plot option elements when moving between absolute and bias plots,
        if do not already exist

        Parameters
        ----------
        data_labels : list
            Data labels
        plot_type : str
            Plot type
        active : str
            'bias' if bias is active, else 'absolute'
        plot_options : list
            Plot options
        z_statistic_sign : str
            Statistic sign
        """

        # if 'bias' is active, remove observations data label from data labels
        data_labels_alt = copy.deepcopy(data_labels)
        if active == "bias":
            data_labels_alt.remove(self.read_instance.observations_data_label)

        # iterate through plot_options
        for plot_option in plot_options:
            if plot_option == "annotate":
                if isinstance(self.plot_axes[plot_type], dict):
                    for relevant_temporal_resolution, sub_ax in self.plot_axes[
                        plot_type
                    ].items():
                        if (
                            relevant_temporal_resolution
                            in self.read_instance.periodic_relevant_temporal_resolutions
                        ):
                            annotation(
                                self.read_instance,
                                self,
                                sub_ax,
                                self.read_instance.networkspeci,
                                data_labels,
                                plot_type,
                                self.plot_characteristics[plot_type],
                                plot_options,
                                plot_z_statistic_sign=z_statistic_sign,
                            )
                            break
                else:
                    annotation(
                        self.read_instance,
                        self,
                        self.plot_axes[plot_type],
                        self.read_instance.networkspeci,
                        data_labels,
                        plot_type,
                        self.plot_characteristics[plot_type],
                        plot_options,
                        plot_z_statistic_sign=z_statistic_sign,
                    )

            elif plot_option == "smooth":
                smooth(
                    self.read_instance,
                    self,
                    self.plot_axes[plot_type],
                    self.read_instance.networkspeci,
                    data_labels_alt,
                    plot_type,
                    self.plot_characteristics[plot_type],
                    plot_options,
                )

            elif plot_option == "threshold":
                if isinstance(self.plot_axes[plot_type], dict):
                    for relevant_temporal_resolution, sub_ax in self.plot_axes[
                        plot_type
                    ].items():
                        if (
                            relevant_temporal_resolution
                            in self.read_instance.periodic_relevant_temporal_resolutions
                        ):
                            threshold(
                                self.read_instance,
                                self,
                                sub_ax,
                                self.read_instance.networkspeci,
                                plot_type,
                                self.plot_characteristics[plot_type],
                            )
                else:
                    threshold(
                        self.read_instance,
                        self,
                        self.plot_axes[plot_type],
                        self.read_instance.networkspeci,
                        plot_type,
                        self.plot_characteristics[plot_type],
                    )

            elif plot_option == "regression":
                linear_regression(
                    self.read_instance,
                    self,
                    self.plot_axes[plot_type],
                    self.read_instance.networkspeci,
                    data_labels_alt,
                    plot_type,
                    self.plot_characteristics[plot_type],
                    plot_options,
                )

    def update_markersize(self, ax, plot_type, markersize, event_source):
        """
        Update markers size for each plot type

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes
        plot_type : str
            Plot type
        markersize : int
            Marker size
        event_source : object
            Event source
        """

        # set markersize
        if plot_type in [
            "timeseries",
            "periodic",
            "scatter",
            "periodic-violin",
            "taylor",
            "fairmode-target",
            "fairmode-statsummary",
        ]:
            if isinstance(ax, dict):
                for sub_ax in ax.values():
                    for line in sub_ax.lines:
                        line.set_markersize(markersize)
            elif isinstance(ax, list):
                for sub_ax in ax:
                    for line in sub_ax.lines:
                        line.set_markersize(markersize)
            else:
                if plot_type == "taylor":
                    for line in self.plotting.taylor_polar_relevant_axis.lines:
                        line.set_markersize(markersize)
                else:
                    for line in ax.lines:
                        line.set_markersize(markersize)

            # update characteristics per plot type
            # this is made to keep the changes when selecting stations with lasso
            if plot_type in [
                "timeseries",
                "periodic",
                "scatter",
                "taylor",
                "fairmode-target",
                "fairmode-statsummary",
            ]:
                self.plot_characteristics[plot_type]["plot"]["markersize"] = markersize
            elif plot_type == "periodic-violin":
                self.plot_characteristics[plot_type]["plot"]["median"][
                    "markersize"
                ] = markersize

        elif plot_type == "map":
            # zero selected and unselected stations
            if event_source == self.interactive_elements[plot_type]["markersize_sl"][0]:
                # actually have zero selected stations currently?
                # if so, update active markersizes
                if len(self.absolute_selected_station_inds) == 0:
                    for collection in self.plot_axes["map"].collections:
                        if isinstance(
                            collection, matplotlib.collections.PathCollection
                        ):
                            markersizes = collection.get_sizes()
                            markersizes[:] = markersize
                            collection.set_sizes(markersizes)

                # actually have selected stations currently?
                # if so, update active opacities
                elif (len(self.absolute_non_selected_station_inds) > 0) & (
                    len(self.absolute_selected_station_inds) > 0
                ):
                    for collection in self.plot_axes["map"].collections:
                        if isinstance(
                            collection, matplotlib.collections.PathCollection
                        ):
                            markersizes = collection.get_sizes()
                            markersizes[
                                self.absolute_non_selected_station_inds
                            ] = markersize
                            collection.set_sizes(markersizes)

                # update characteristics per plot type for unselected stations
                self.plot_characteristics["map"]["marker_unselected"]["s"] = markersize

                # update characteristics per plot type for zero selected stations
                self.plot_characteristics["map"]["marker_zero_stations_selected"][
                    "s"
                ] = markersize

            # selected stations
            elif (
                event_source == self.interactive_elements[plot_type]["markersize_sl"][1]
            ):
                # actually have selected stations currently?
                # if so, update active opacities
                if len(self.absolute_selected_station_inds) > 0:
                    for collection in self.plot_axes["map"].collections:
                        if isinstance(
                            collection, matplotlib.collections.PathCollection
                        ):
                            markersizes = collection.get_sizes()
                            markersizes[
                                self.absolute_selected_station_inds
                            ] = markersize
                            collection.set_sizes(markersizes)

                # update characteristics per plot type
                self.plot_characteristics["map"]["marker_selected"]["s"] = markersize

        # draw changes
        self.figure.canvas.draw_idle()

        return None

    def update_opacity(self, ax, plot_type, opacity, event_source):
        """
        Update markers opacity for each plot type

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes
        plot_type : str
            Plot type
        opacity : float
            Opacity
        event_source : object
            Event source
        """

        # set opacity
        if plot_type == "map":
            # zero selected and unselected stations
            if event_source == self.interactive_elements[plot_type]["opacity_sl"][0]:
                # actually have zero selected stations currently?
                # if so, update active opacities
                if len(self.absolute_selected_station_inds) == 0:
                    for collection in self.plot_axes["map"].collections:
                        if isinstance(
                            collection, matplotlib.collections.PathCollection
                        ):
                            if Version(matplotlib.__version__) < Version("3.4"):
                                opacities = collection.get_facecolor()
                                opacities[:, -1] = opacity
                                collection.set_facecolor(opacities)
                            else:
                                opacities = collection.get_facecolor()[:, -1]
                                opacities[:] = opacity
                                collection.set_alpha(opacities)

                # actually have selected stations currently?
                # if so, update active opacities
                elif (len(self.absolute_non_selected_station_inds) > 0) & (
                    len(self.absolute_selected_station_inds) > 0
                ):
                    for collection in self.plot_axes["map"].collections:
                        if isinstance(
                            collection, matplotlib.collections.PathCollection
                        ):
                            if Version(matplotlib.__version__) < Version("3.4"):
                                opacities = collection.get_facecolor()
                                opacities[
                                    self.absolute_non_selected_station_inds, -1
                                ] = opacity
                                collection.set_facecolor(opacities)
                            else:
                                opacities = collection.get_facecolor()[:, -1]
                                opacities[
                                    self.absolute_non_selected_station_inds
                                ] = opacity
                                collection.set_alpha(opacities)

                # update characteristics per plot type for unselected stations
                self.plot_characteristics["map"]["marker_zero_stations_selected"][
                    "alpha"
                ] = opacity

                # update characteristics per plot type for zero selected stations
                self.plot_characteristics["map"]["marker_unselected"]["alpha"] = opacity

            # selected stations
            elif event_source == self.interactive_elements[plot_type]["opacity_sl"][1]:
                # actually have selected stations currently?
                # if so, update active opacities
                if len(self.absolute_selected_station_inds) > 0:
                    for collection in self.plot_axes["map"].collections:
                        if isinstance(
                            collection, matplotlib.collections.PathCollection
                        ):
                            if Version(matplotlib.__version__) < Version("3.4"):
                                opacities = collection.get_facecolor()
                                opacities[
                                    self.absolute_selected_station_inds, -1
                                ] = opacity
                                collection.set_facecolor(opacities)
                            else:
                                opacities = collection.get_facecolor()[:, -1]
                                opacities[self.absolute_selected_station_inds] = opacity
                                collection.set_alpha(opacities)

                # update characteristics per plot type
                self.plot_characteristics["map"]["marker_selected"]["alpha"] = opacity

        # redraw points
        self.figure.canvas.draw_idle()

        return None

    def update_linewidth(self, ax, plot_type, linewidth):
        """
        Update line widths for each plot type

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes
        plot_type : str
            Plot type
        linewidth : int
            Line width
        """

        # set linewidth
        if isinstance(ax, dict):
            for sub_ax in ax.values():
                for line in sub_ax.lines:
                    line.set_linewidth(linewidth)
        else:
            for line in ax.lines:
                if (plot_type == "scatter") and (
                    (list(line.get_xdata()) == [0, 0.5])
                    or (list(line.get_xdata()) == [0, 1])
                ):
                    continue
                else:
                    line.set_linewidth(linewidth)

        # update characteristics per plot type
        if plot_type == "periodic-violin":
            self.plot_characteristics[plot_type]["plot"]["median"][
                "linewidth"
            ] = linewidth
        elif plot_type == "timeseries":
            self.plot_characteristics[plot_type]["smooth"]["format"][
                "linewidth"
            ] = linewidth
        elif plot_type == "regression":
            self.plot_characteristics[plot_type]["regression"]["linewidth"] = linewidth
        else:
            self.plot_characteristics[plot_type]["plot"]["linewidth"] = linewidth

        # draw changes
        self.figure.canvas.draw_idle()

        return None

    def update_smooth_window(self, plot_type, smooth_window):
        """
        Update smooth window

        Parameters
        ----------
        plot_type : str
            Plot type
        smooth_window : int
            Smooth window
        """

        # update characteristics per plot type
        self.plot_characteristics[plot_type]["smooth"]["window"] = smooth_window

        # get index of smooth in plot options
        all_plot_options = self.plot_characteristics[plot_type]["plot_options"]
        index = all_plot_options.index("smooth")

        # remove smooth plot option
        self.timeseries_options.model().item(index).setCheckState(QtCore.Qt.Unchecked)

        # create smooth lines
        if smooth_window > 0:
            # add smooth plot option
            self.timeseries_options.model().item(index).setCheckState(QtCore.Qt.Checked)

        # draw changes
        self.figure.canvas.draw_idle()

        return None

    def update_smooth_min_points(self, plot_type, smooth_min_points):
        """
        Update smooth minimum points

        Parameters
        ----------
        plot_type : str
            Plot type
        smooth_min_points : int
            Smooth minimum points
        """

        # update characteristics per plot type
        self.plot_characteristics[plot_type]["smooth"][
            "min_points_percentage"
        ] = smooth_min_points

        # get window to check if we need to redraw
        window = self.timeseries_smooth_window_sl.value()

        # get index of smooth in plot options
        all_plot_options = self.plot_characteristics[plot_type]["plot_options"]
        index = all_plot_options.index("smooth")

        # remove smooth plot option
        self.timeseries_options.model().item(index).setCheckState(QtCore.Qt.Unchecked)

        if window > 0:
            # add smooth plot option
            self.timeseries_options.model().item(index).setCheckState(QtCore.Qt.Checked)

            # draw changes
            self.figure.canvas.draw_idle()

        return None

    def update_violin_widths(self, ax, plot_type, width):
        """
        Update violin widths for violin plots

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes
        plot_type : str
            Plot type
        width : int
            Width
        """

        # set violin widths
        if plot_type == "periodic-violin":
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
        self.plot_characteristics[plot_type]["plot"]["violin"]["widths"] = widths

        # draw changes
        self.figure.canvas.draw_idle()

        return None

    def save_axis_figure_dialog(self, plot_type, relevant_temporal_resolution=None):
        """
        Function to create the dialog box to save each plot figure

        Parameters
        ----------
        plot_type : str
            Plot type
        relevant_temporal_resolution : str
            Temporal resolution
        """

        default_filename = "{0}-{1}-{2}-{3}-{4}-{5}-{6}.png".format(
            self.read_instance.network[0],
            self.read_instance.species[0],
            self.read_instance.resolution,
            self.read_instance.start_date,
            self.read_instance.end_date,
            plot_type,
            str(relevant_temporal_resolution),
        )

        if relevant_temporal_resolution is None:
            default_filename = default_filename.split("-None")[0] + ".png"

        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        figure_path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self.read_instance,
            "Choose folder to save figure",
            default_filename,
            "All Files (*);;Figures (*.png)",
            options=options,
        )
        if figure_path:
            return figure_path

    def save_axis_figure_func(self):
        """
        Function to save each plot figure
        """

        # get option and plot names
        event_source = self.sender()
        plot_type = event_source.objectName().split("_save")[0]
        if plot_type in ["periodic_violin", "fairmode_target", "fairmode_statsummary"]:
            plot_type = plot_type.replace("_", "-")

        # set extent expansion
        for i, position in enumerate(
            [
                self.read_instance.position_1,
                self.read_instance.position_2,
                self.read_instance.position_3,
                self.read_instance.position_4,
                self.read_instance.position_5,
            ]
        ):
            if plot_type == position:
                if i + 1 == 1:
                    expand_x, expand_y = 1.4, 1.4
                elif i + 1 == 2:
                    if plot_type == "scatter":
                        expand_x, expand_y = 1.30, 1.25
                    else:
                        expand_x, expand_y = 1.2, 1.3
                else:
                    if plot_type == "scatter":
                        expand_x, expand_y = 1.25, 1.2
                    else:
                        expand_x, expand_y = 1.2, 1.2
                continue

        # remove titles
        for key in self.read_instance.active_dashboard_plots:
            if key != "None":
                if isinstance(self.plot_axes[key], dict):
                    for relevant_temporal_resolution, sub_ax in self.plot_axes[
                        key
                    ].items():
                        if relevant_temporal_resolution in ["hour"]:
                            sub_ax.set_title(
                                label="",
                                fontsize=self.plot_characteristics[key]["axis_title"][
                                    "fontsize"
                                ],
                                loc=self.plot_characteristics[key]["axis_title"]["loc"],
                            )
                else:
                    self.plot_axes[key].set_title(
                        label="",
                        fontsize=self.plot_characteristics[key]["axis_title"][
                            "fontsize"
                        ],
                        loc=self.plot_characteristics[key]["axis_title"]["loc"],
                    )

        # hide colourbar
        if plot_type != "map":
            self.plot_axes["cb"].set_visible(False)

        # draw changes
        self.figure.canvas.draw_idle()

        # make screenshot and save
        if isinstance(self.plot_axes[plot_type], dict):
            for relevant_temporal_resolution, sub_ax in self.plot_axes[
                plot_type
            ].items():
                extent = sub_ax.get_window_extent().transformed(
                    self.figure.dpi_scale_trans.inverted()
                )
                if relevant_temporal_resolution == "hour":
                    expand_x, expand_y = 1.15, 1.3
                elif relevant_temporal_resolution == "dayofweek":
                    expand_x, expand_y = 1.25, 1.3
                elif relevant_temporal_resolution == "month":
                    expand_x, expand_y = 1.15, 1.3

                # get folder where figure will be saved
                figure_path = self.save_axis_figure_dialog(
                    plot_type, relevant_temporal_resolution
                )

                # save figure
                if figure_path is not None:
                    self.figure.savefig(
                        figure_path, bbox_inches=extent.expanded(expand_x, expand_y)
                    )
        else:
            extent = (
                self.plot_axes[plot_type]
                .get_window_extent()
                .transformed(self.figure.dpi_scale_trans.inverted())
            )

            # get folder where figure will be saved
            figure_path = self.save_axis_figure_dialog(plot_type)

            # save figure
            if figure_path is not None:
                self.figure.savefig(
                    figure_path, bbox_inches=extent.expanded(expand_x, expand_y)
                )

        # add titles
        for key in self.read_instance.active_dashboard_plots:
            if key != "None":
                if isinstance(self.plot_axes[key], dict):
                    for relevant_temporal_resolution, sub_ax in self.plot_axes[
                        key
                    ].items():
                        if relevant_temporal_resolution in ["hour"]:
                            sub_ax.set_title(
                                **self.plot_characteristics[key]["axis_title"]
                            )
                else:
                    self.plot_axes[key].set_title(
                        **self.plot_characteristics[key]["axis_title"]
                    )

        # show colourbar
        self.plot_axes["cb"].set_visible(True)

        # draw changes
        self.figure.canvas.draw_idle()

        return None

    def save_axis_figure_data_dialog(self):
        """
        Function to create the dialog box to choose directory where each plot figure data will be saved
        """

        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self.read_instance, "Choose folder to save data", options=options
        )
        if path:
            return path

    def save_axis_figure_data_func(self):
        """
        Function to save each plot figure data
        """

        # get option and plot names
        event_source = self.sender()
        plot_type = event_source.objectName().split("_save")[0]
        if plot_type in ["periodic_violin", "fairmode_target", "fairmode_statsummary"]:
            plot_type = plot_type.replace("_", "-")
        plot_options = copy.deepcopy(self.current_plot_options[plot_type])

        tests_generate_output = False
        if plot_type == "map":
            labela = self.map_z1.currentText()
            labelb = self.map_z2.currentText()
        else:
            labela = ""
            labelb = ""

        # get folder where data will be saved
        path = self.save_axis_figure_data_dialog()

        # save figure
        if path is not None:
            download_plot_data_to_csv(
                self.read_instance,
                self,
                plot_type,
                plot_type,
                plot_options,
                path,
                self.read_instance.networkspeci,
                tests_generate_output,
                labela,
                labelb,
            )

    def update_aggregation_statistic(self):
        """
        Update general aggregation statistic
        """

        # get statistic
        self.read_instance.selected_statistic_aggregation = (
            self.read_instance.cb_statistic_aggregation.currentText()
        )

        # update statistic in memory
        self.read_instance.statistic_aggregation = (
            self.read_instance.selected_statistic_aggregation
        )

    def update_timeseries_aggregation_statistic(self):
        """
        Update timeseries aggregation statistic
        """

        # get statistic
        self.read_instance.selected_timeseries_statistic_aggregation = (
            self.timeseries_stat.currentText()
        )

        # update statistic in memory
        self.read_instance.timeseries_statistic_aggregation = (
            self.read_instance.selected_timeseries_statistic_aggregation
        )

    def update_option_on_combobox(self, event_source, index, uncheck=True):
        """
        Check or uncheck option in combobox dropdown

        Parameters
        ----------
        event_source : object
            Event source
        index : int
            Index
        uncheck : bool
            Check or uncheck combobox
        """

        self.read_instance.block_MPL_canvas_updates = True
        if uncheck:
            event_source.model().item(index).setCheckState(QtCore.Qt.Unchecked)
        else:
            event_source.model().item(index).setCheckState(QtCore.Qt.Checked)
        self.read_instance.block_MPL_canvas_updates = False

        return None
