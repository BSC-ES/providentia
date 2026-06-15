""" Class for the dashboard matplotlib navigation toolbar """

import copy
from enum import Enum
import os

import matplotlib
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from packaging.version import Version
from PyQt5 import QtGui, QtWidgets

from providentia.auxiliar import CURRENT_PATH, join
from .configuration import ProvConfiguration
from .configuration import load_conf
from .dashboard_elements import InputDialog
from .dashboard_interactivity import LassoSelector
from .fields_menus import multispecies_conf
from .plot_aux import get_map_extent
from .plot_formatting import harmonise_xy_lims_paradigm
from .read_aux import generate_file_trees
from .warnings_prv import show_message
from .writing import export_configuration, export_data_npz, export_netcdf


class _Mode(str, Enum):
    """
    Defines the interactive tool modes available for the user interface.
    """

    NONE = ""
    PAN = "pan/zoom"
    ZOOM = "zoom rect"
    LASSO = "lasso"

    def __str__(self):
        """
        Returns the string value of the current mode.

        Returns
        -------
        str
            The human-readable name of the mode.
        """
        return self.value

    @property
    def _navigate_mode(self):
        """
        Retrieves the navigation internal state name for the backend.

        Returns
        -------
        str or None
            The uppercase name of the mode if active, otherwise None.
        """
        return self.name if self is not _Mode.NONE else None


class NavigationToolbar(NavigationToolbar2QT):
    """Class that customises the Matplotlib navigation toolbar with specific buttons and bespoke functionality."""

    def __init__(self, read_instance=None, canvas_instance=None):
        """
        Initialises the toolbar, defines the toolitems, and sets custom icons.

        Parameters
        ----------
        read_instance : object, optional
            The source instance containing the data and configuration (default is None).
        canvas_instance : object, optional
            The Matplotlib canvas instance to which this toolbar is attached (default is None).
        """

        self.read_instance = read_instance
        self.canvas_instance = canvas_instance

        # only display wanted buttons
        NavigationToolbar2QT.toolitems = (
            (
                "Save data",
                "Save current instance of data and metadata",
                "",
                "save_data",
            ),
            (
                "Load data",
                "Load toolbar selections from configuration file",
                "",
                "conf_dialogs",
            ),
            (None, None, None, None),
            ("World", "Set world view", "world", "world"),
            ("Home", "Set to original view", "home", "home"),
            ("Back", "Back to previous view", "back", "back"),
            ("Forward", "Forward to next view", "forward", "forward"),
            (None, None, None, None),
            ("Zoom", "Zoom to rectangle", "zoom_to_rect", "zoom"),
            ("Pan", "Pan axes with left mouse, zoom with right", "move", "pan"),
            ("Lasso", "Select stations with lasso", "", "connect_lasso"),
            (None, None, None, None),
            ("Save figure", "Save the figure", "filesave", "save_figure"),
        )

        # allow access to methods of parent class NavigationToolbar2QT
        super(NavigationToolbar, self).__init__(canvas_instance, read_instance)

        # set toolbar icons
        actions = self.findChildren(QtWidgets.QAction)
        self._actions["save_data"].setIcon(
            QtGui.QIcon(join(CURRENT_PATH, "resources/save_icon.png"))
        )
        self._actions["conf_dialogs"].setIcon(
            QtGui.QIcon(join(CURRENT_PATH, "resources/conf_icon.png"))
        )
        self._actions["world"].setIcon(
            QtGui.QIcon(join(CURRENT_PATH, "resources/world_icon.png"))
        )
        self._actions["home"].setIcon(
            QtGui.QIcon(join(CURRENT_PATH, "resources/original_view_icon.png"))
        )
        self._actions["zoom"].setIcon(
            QtGui.QIcon(join(CURRENT_PATH, "resources/zoom_icon.png"))
        )
        self._actions["connect_lasso"].setIcon(
            QtGui.QIcon(join(CURRENT_PATH, "resources/lasso_icon.png"))
        )
        self._actions["save_figure"].setIcon(
            QtGui.QIcon(join(CURRENT_PATH, "resources/save_fig_icon.png"))
        )

        # allow lasso button to be checked
        self._actions["connect_lasso"].setCheckable(True)

    def _update_buttons_checked(self):
        """
        Synchronises the visual check states of toolbar buttons with the current interaction mode.
        """

        if "pan" in self._actions:
            self._actions["pan"].setChecked(self.mode == "pan/zoom")
        if "zoom" in self._actions:
            self._actions["zoom"].setChecked(self.mode == "zoom rect")
        if "connect_lasso" in self._actions:
            self._actions["connect_lasso"].setChecked(self.mode == "lasso")

    def save_data(self):
        """
        Triggers a file dialogue to export data or configuration settings to various formats.
        """

        filetypes = {"NetCDF": "nc", "Numpy file": "npz", "Configuration": "conf"}
        sorted_filetypes = sorted(filetypes.items())
        startpath = os.path.expanduser(matplotlib.rcParams["savefig.directory"])
        daterange = (
            self.read_instance.le_start_date.text()
            + "_"
            + self.read_instance.le_end_date.text()
        )
        default_name = "PRV_" + str(self.read_instance.species[0]) + "_" + daterange
        start = join(startpath, default_name)

        filter_ext = [
            "%s (%s)" % (name, "*.%s" % ext) for name, ext in sorted_filetypes
        ]
        filter_ext = ";;".join(filter_ext)

        # prompt with file name and extension
        fname, fext = QtWidgets.QFileDialog.getSaveFileName(
            None, "Choose a filename to save to", start, filter_ext
        )
        chose_npz = "npz" in fext
        chose_conf = "conf" in fext
        if fname:
            # save dir for next time, unless empty str (i.e., use cwd).
            if startpath != "":
                matplotlib.rcParams["savefig.directory"] = os.path.dirname(fname)
                try:
                    if chose_npz:
                        output = export_data_npz(
                            self.read_instance, fname, input_dialogue=True
                        )
                    elif chose_conf:
                        output = export_configuration(self.read_instance, fname)
                    else:
                        output = export_netcdf(
                            self.read_instance, fname, input_dialogue=True
                        )
                    if output is True:
                        msg = "The data was successfully saved in {}.".format(fname)
                        show_message(self.read_instance, msg)
                except Exception as e:
                    msg = "There was an error saving the file."
                    self.read_instance.logger.error(e)
                    show_message(self.read_instance, msg)

    def check_for_axis_limit_changes(self, previous_state, current_state):
        """
        Identifies which plot axes have had their limits modified and initiates visual synchronisation.

        Parameters
        ----------
        previous_state : dict
            A dictionary mapping axis objects to their prior (x_limit, y_limit) tuples.
        current_state : dict
            A dictionary mapping axis objects to their new (x_limit, y_limit) tuples.
        """

        # check which limit changed
        for axis, (prev_xlim, prev_ylim) in previous_state.items():
            new_xlim, new_ylim = current_state[axis]

            # find the different limit
            if prev_xlim != new_xlim or prev_ylim != new_ylim:
                self.harmonise_changed_axis(axis)

    def harmonise_changed_axis(self, axis):
        """
        Identifies the context of a modified axis and propagates limit changes to related plots.

        Parameters
        ----------
        axis : matplotlib.axes.Axes
            The specific axis object that has undergone a change in limits.
        """

        # iterate through each plot until the one related to the axis is found
        for plot_type, axes in self.canvas_instance.plot_axes.items():
            # for periodic plots, check the the axes inside the dictionary
            if plot_type == "periodic" or plot_type == "periodic-violin":
                for periodic_type in axes:
                    if axes[periodic_type] == axis:
                        break
            # compare the main axis for the rest of the plots
            elif axes == axis:
                break

        # apply harmonise to the plots with time
        if plot_type in ["periodic", "periodic-violin", "timeseries"]:
            plot_options = copy.deepcopy(
                self.canvas_instance.current_plot_options[plot_type]
            )
            harmonise_xy_lims_paradigm(
                self.read_instance,
                self.canvas_instance,
                self.canvas_instance.plot_axes[plot_type],
                plot_type,
                self.canvas_instance.plot_characteristics[plot_type],
                plot_options,
                relim=True,
                autoscale=False,
            )

    def world(self):
        """
        Resets the map view to global coordinates and updates the navigation history.
        """

        self.set_history_buttons()
        ax = self.canvas_instance.plot_axes["map"]
        ax.set_extent([-180, 180, -90, 90])
        self.push_current()
        self.canvas.draw_idle()

        # update map extent
        self.read_instance.map_extent = get_map_extent(self.canvas_instance)

    def home(self):
        """
        Restores the original view of all plots and synchronises limits across linked axes.
        """

        # get all the limits before doing clicking on home
        previous_state = {
            axis: (axis.get_xlim(), axis.get_ylim())
            for axis in self.canvas_instance.figure.axes
        }

        super().home()

        # get all the limits after doing clicking on home
        current_state = {
            axis: (axis.get_xlim(), axis.get_ylim())
            for axis in self.canvas_instance.figure.axes
        }

        # harmonise axis if needed
        self.check_for_axis_limit_changes(previous_state, current_state)

        # update map extent
        self.read_instance.map_extent = get_map_extent(self.canvas_instance)

    def back(self):
        """
        Reverts to the previous view in the navigation stack and synchronises all linked axes.
        """

        # get all the limits before doing clicking on back
        previous_state = {
            axis: (axis.get_xlim(), axis.get_ylim())
            for axis in self.canvas_instance.figure.axes
        }

        super().back()

        # get all the limits after doing clicking on back
        current_state = {
            axis: (axis.get_xlim(), axis.get_ylim())
            for axis in self.canvas_instance.figure.axes
        }

        # harmonise axis if needed
        self.check_for_axis_limit_changes(previous_state, current_state)

        # update map extent
        self.read_instance.map_extent = get_map_extent(self.canvas_instance)

    def forward(self):
        """
        Moves to the next view in the navigation stack and synchronises all linked axes.
        """

        # get all the limits before doing clicking on forward
        previous_state = {
            axis: (axis.get_xlim(), axis.get_ylim())
            for axis in self.canvas_instance.figure.axes
        }

        super().forward()

        # get all the limits after doing clicking on forward
        current_state = {
            axis: (axis.get_xlim(), axis.get_ylim())
            for axis in self.canvas_instance.figure.axes
        }

        # harmonise axis if needed
        self.check_for_axis_limit_changes(previous_state, current_state)

        # update map extent
        self.read_instance.map_extent = get_map_extent(self.canvas_instance)

    def drag_pan(self, event):
        """
        Handles the continuous panning of plot axes and synchronises related views in real-time.

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            The mouse event triggered by dragging the plot surface.
        """

        super().drag_pan(event)

        # harmonise axis if needed
        self.harmonise_changed_axis(event.inaxes)

        # update map extent
        self.read_instance.map_extent = get_map_extent(self.canvas_instance)

    def release_zoom(self, event):
        """
        Finalises the zoom-to-rectangle action and synchronises all linked axes to the new view.

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            The mouse release event that defines the final zoom boundary.
        """

        super().release_zoom(event)

        # harmonise axis if needed
        self.harmonise_changed_axis(event.inaxes)

        # update map extent
        self.read_instance.map_extent = get_map_extent(self.canvas_instance)

    def save_figure(self):
        """
        Validates data presence before triggering the standard dialogue to save the current plot as an image.
        """

        if (
            self.read_instance.le_minimum_value.text() == ""
            and self.read_instance.le_minimum_value.text() == ""
        ):
            # control that data was read
            msg = "The figure could not be created. Read the data first."
            show_message(self.read_instance, msg)
            return

        # inherit from qt
        super().save_figure(self)

    def conf_dialogs(self):
        """
        Manages the selection and loading of configuration files through sequential user dialogues.
        """

        conf_to_load = self.filename_dialog()

        # is user pressed cancel
        if conf_to_load is None:
            return

        try:
            load_conf(self.read_instance, fpath=conf_to_load)
            all_sections = self.read_instance.sub_opts.keys()

            if len(all_sections) == 1:
                okpressed = False
                selected_section = list(all_sections)[0]
            else:
                title = "Sections"
                msg = "Select section to load"
                dialog = InputDialog(self, title, msg, all_sections)
                selected_section, okpressed = dialog.selected_option, dialog.okpressed

            if okpressed or len(all_sections) == 1:
                self.reload_conf(selected_section, conf_to_load)

        except Exception as e:
            msg = "There was an error loading the configuration file."
            show_message(self.read_instance, msg)
            self.read_instance.logger.error(e)

    def connect_lasso(self):
        """
        Toggles the free-form lasso selection tool on the map for station selection.
        """

        if not self.canvas_instance.figure.canvas.widgetlock.available(self):
            self.set_message("lasso unavailable")
            self.mode = _Mode.NONE
            self.canvas_instance.lasso_active = False
            self._update_buttons_checked()
            return

        # if lasso button is pressed then activate lasso event
        if self._actions["connect_lasso"].isChecked():
            # update mode
            self.mode = _Mode.LASSO
            self.canvas_instance.lasso_active = True

            # connect lasso event
            if (Version(matplotlib.__version__) < Version("3.2")) or (
                self.read_instance.machine in ["mn5", "nord4"]
            ):
                blit = False
            else:
                blit = True
            self.canvas_instance.lasso_event = LassoSelector(
                self.canvas_instance.plot_axes["map"],
                onselect=self.canvas_instance.onlassoselect,
                useblit=blit,
                props=self.canvas_instance.lasso,
                button=[1],
            )

            # release canvas drawing (from self if owned by other toolbar buttons)
            if self.canvas_instance.figure.canvas.widgetlock.isowner(self):
                self.canvas_instance.figure.canvas.widgetlock.release(self)

            # ensure lasso is initialised correctly by forcing a draw
            self.canvas_instance.figure.canvas.draw_idle()

        # otherwise, deactivate lasso event
        else:
            # update mode
            self.mode = _Mode.NONE
            self.canvas_instance.lasso_active = False

            # disconnect lasso event
            self.canvas_instance.lasso_event.disconnect_events()

        # update checked buttons
        self._update_buttons_checked()

    def filename_dialog(self):
        """
        Opens a file selection dialogue for the user to locate and choose a configuration file.

        Returns
        -------
        str or None
            The absolute path to the selected file, or None if the dialogue was cancelled.
        """

        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(
            self.read_instance,
            "Choose configuration file to load",
            "",
            "All Files (*);;Python Files (*.py",
            options=options,
        )
        if filename:
            return filename

    def reload_conf(self, section, fpath):
        """
        Resets the current session and reinitialises the dashboard with parameters from a specific configuration section.

        Parameters
        ----------
        section : str
            The specific section name within the configuration file to be loaded.
        fpath : str
            The file path to the configuration file.
        """

        # save current active_dashboard_plots
        previous_active_dashboard_plots = copy.deepcopy(
            self.read_instance.active_dashboard_plots
        )

        # delete current active config attributes
        for k in self.read_instance.current_config:
            try:
                vars(self.read_instance).pop(k)
            except:
                pass

        # reinitialise default configuration variables
        # no commandline arguments considered as reloading config
        provconf = ProvConfiguration(self.read_instance, **{"dashboard": True})

        # get current GHOST version
        current_ghost_version = self.read_instance.ghost_version

        # update config and section attributes of instance
        self.read_instance.config = fpath
        self.read_instance.section = section
        self.read_instance.from_conf = True
        self.read_instance.current_config = self.read_instance.sub_opts[section]
        for k, val in self.read_instance.current_config.items():
            setattr(self.read_instance, k, provconf.parse_parameter(k, val))

        # now all variables have been parsed, check validity of those, throwing errors where necessary
        provconf.check_validity()

        # get new active_dashboard_plots and set it in memory to be the previous variable
        new_active_dashboard_plots = copy.deepcopy(
            self.read_instance.active_dashboard_plots
        )
        self.read_instance.active_dashboard_plots = previous_active_dashboard_plots

        # generate file trees if GHOST version has changed
        if current_ghost_version != self.read_instance.ghost_version:
            generate_file_trees(self.read_instance)

        # update species, models, qa & flags
        self.read_instance.config_bar_initialisation = True
        self.read_instance.update_configuration_bar_fields()
        self.read_instance.config_bar_initialisation = False

        # read and filter
        self.read_instance.handle_data_selection_update()

        # update multispecies menu from conf
        multispecies_conf(self.read_instance)

        # update active dashboard plots
        for position, (previous_plot_type, new_plot_type) in enumerate(
            zip(previous_active_dashboard_plots, new_active_dashboard_plots)
        ):
            if previous_plot_type != new_plot_type:
                self.read_instance.handle_layout_update(
                    new_plot_type, sender=position + 2
                )

        # reset from_conf variable to False after reading and filtering is complete
        self.read_instance.from_conf = False
