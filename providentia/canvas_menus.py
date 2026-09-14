""" Class for plot settings menus """

from functools import partial
import platform

from PyQt5 import QtCore, QtGui, QtWidgets
import yaml

from providentia.auxiliar import CURRENT_PATH, join
from .dashboard_elements import CheckableComboBox, ComboBox, MenuLineEdit
from .dashboard_elements import set_formatting

PROVIDENTIA_ROOT = "/".join(CURRENT_PATH.split("/")[:-1])
settings_dict = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings/internal/canvas_menus.yaml"))
)
# get operating system specific formatting
operating_system = platform.system()
if operating_system == "Darwin":
    formatting_dict = yaml.safe_load(
        open(join(PROVIDENTIA_ROOT, "settings/internal/stylesheet_mac.yaml"))
    )
elif operating_system == "Linux":
    formatting_dict = yaml.safe_load(
        open(join(PROVIDENTIA_ROOT, "settings/internal/stylesheet_linux.yaml"))
    )
elif operating_system in ["Windows", "MINGW32_NT", "MINGW64_NT"]:
    formatting_dict = yaml.safe_load(
        open(join(PROVIDENTIA_ROOT, "settings/internal/stylesheet_windows.yaml"))
    )


# how faded a slider is drawn while it is disabled. Faded rather than
# restyled, as a slider given a stylesheet stops being drawn by the platform
# altogether and so changes appearance when enabled as well - and the platform
# draws a disabled slider identically to an enabled one, leaving the sliders
# held by an automatic setting (map point sizing, histogram bins) looking as
# though they can still be dragged
DISABLED_SLIDER_OPACITY = 0.35

# gap between each settings menu slider and the value shown to its right - the
# same as the margin between a menu panel's edge and the controls inside it
# (see canvas_menus.yaml). The value's own width is whatever the widest value
# any slider in the menu can reach needs (see ValueSlider.fit_group())
SLIDER_READOUT_GAP = 10

# length and abbreviated unit of one step of each temporal resolution, for a
# slider counting steps to be shown as the duration it spans
TIMESTEP_DURATIONS = {
    "hourly": (1, "h"),
    "3hourly": (3, "h"),
    "6hourly": (6, "h"),
    "daily": (1, "d"),
    "monthly": (1, "mo"),
    "annual": (1, "yr"),
}


def get_smoothing_timestep(canvas_instance):
    """
    Get the temporal resolution of the points the timeseries smoothing window
    counts: the chunks when the timeseries is chunked, otherwise the data's
    own (see plot_options.smooth()).

    Parameters
    ----------
    canvas_instance : instance
        Canvas instance

    Returns
    -------
    str or None
        Temporal resolution, or None if it cannot be told yet
    """

    chunk_stat = getattr(canvas_instance, "timeseries_chunk_stat", None)
    chunk_resolution = getattr(canvas_instance, "timeseries_chunk_resolution", None)
    if (chunk_stat is not None) and (chunk_resolution is not None):
        chunk_texts = [chunk_stat.currentText(), chunk_resolution.currentText()]
        if not any(text in ["", "None"] for text in chunk_texts):
            return chunk_resolution.currentText()

    return getattr(canvas_instance.read_instance, "active_resolution", None)


def set_slider_enabled(slider, enabled):
    """
    Set whether a settings menu slider can be used, fading it while it cannot.

    Parameters
    ----------
    slider : QtWidgets.QSlider
        Slider to enable or disable
    enabled : bool
        Whether the slider can be used
    """

    slider.setEnabled(enabled)

    # the number showing the slider's value fades along with it
    widgets = [slider]
    if isinstance(slider, ValueSlider):
        widgets.append(slider.readout)

    for widget in widgets:
        if enabled:
            widget.setGraphicsEffect(None)
        else:
            faded = QtWidgets.QGraphicsOpacityEffect(widget)
            faded.setOpacity(DISABLED_SLIDER_OPACITY)
            widget.setGraphicsEffect(faded)

    return None


class ValueSlider(QtWidgets.QSlider):
    """
    Settings menu slider with its value, and its unit, shown to its right.

    The value follows the handle while it is dragged, although the plot is
    only redrawn once it is let go, so the value being moved to can be seen
    before it is applied. It also follows every other change of value,
    including those made in code with the slider's signals blocked (e.g. an
    automatic setting keeping the slider in step with what it has drawn), and
    is shown, hidden and raised along with the slider.

    The slider and its value share the width the slider is given in
    canvas_menus.yaml, the value ending where it does. Every slider in a menu
    gives its value the same width - that of the widest value any of them can
    reach - so that they all end at the same point, and are fitted again
    whenever what they can reach changes.
    """

    def __init__(self, parent, value_divisor=1, unit=None, get_timestep=None):
        """
        Initialise slider and the value shown beside it

        Parameters
        ----------
        parent : QtWidgets.QWidget
            Widget the slider and its value are drawn on
        value_divisor : int
            What the slider's integer position is divided by to give the value
            it sets (e.g. 10 for an opacity of 0-1 set in steps of 0.1), which
            is what is shown
        unit : str, optional
            Abbreviated unit shown after the value, or "timesteps" for a count
            of time steps, shown as the duration it spans (default is None,
            i.e. no unit)
        get_timestep : callable, optional
            Returns the temporal resolution of the steps a "timesteps" slider
            counts - see get_smoothing_timestep()
        """

        super().__init__(QtCore.Qt.Horizontal, parent)

        self.value_divisor = value_divisor
        self.unit = unit
        self.get_timestep = get_timestep

        # width shared by the slider and its value, the sliders fitted
        # together with it (a menu's) and the width given to the value - set
        # by SettingsMenu and fit_group()
        self.footprint_width = None
        self.group = [self]
        self.readout_width = 0

        self.readout = set_formatting(
            QtWidgets.QLabel(parent), formatting_dict["settings_sublabel"]
        )
        self.readout.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.readout.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        self.readout.hide()

        # while dragging the value is not changed until the handle is let go
        # (tracking is off), so the position being dragged to is shown instead
        self.sliderMoved.connect(self.show_value)
        self.show_value(self.value())

    def format_value(self, position):
        """
        Get the text for the value a slider position sets, with its unit.

        Parameters
        ----------
        position : int
            Slider position

        Returns
        -------
        str
            Value with its unit
        """

        if self.unit == "timesteps":
            # the duration the steps span, when their length is known
            timestep = self.get_timestep() if self.get_timestep else None
            step_length, abbreviation = TIMESTEP_DURATIONS.get(
                str(timestep).replace("_instantaneous", ""), (None, None)
            )
            if step_length is None:
                text = str(position)
            else:
                text = "{} {}".format(position * step_length, abbreviation)
        else:
            if self.value_divisor == 1:
                text = str(position)
            else:
                text = "{:.1f}".format(position / self.value_divisor)
            if self.unit == "%":
                text += "%"
            elif self.unit:
                text += " " + self.unit

        return text

    def show_value(self, position):
        """
        Show the value a slider position sets.

        Parameters
        ----------
        position : int
            Slider position
        """

        self.readout.setText(self.format_value(position))

        return None

    def refresh_readout(self):
        """
        Show the value again, and fit the menu's sliders to it, for a unit that
        has changed without the value doing so (e.g. the temporal resolution a
        count of steps is in).
        """

        self.show_value(self.sliderPosition())
        self.fit_group()

        return None

    def widest_value_width(self):
        """
        Get the width the widest value the slider can reach takes up.

        Returns
        -------
        int
            Width in pixels
        """

        # measured in the font the value is shown in, which a stylesheet only
        # sets once the label has been polished
        self.readout.ensurePolished()
        metrics = self.readout.fontMetrics()

        return max(
            metrics.horizontalAdvance(self.format_value(position))
            for position in (self.minimum(), self.maximum())
        )

    def fit_group(self):
        """
        Size the values of this slider's menu to the widest any of its sliders
        can reach, shortening the sliders to make room for them.
        """

        sliders = [slider for slider in self.group if slider.footprint_width]
        if not sliders:
            return None

        readout_width = max(slider.widest_value_width() for slider in sliders)
        for slider in sliders:
            slider.readout_width = readout_width
            slider.resize(
                slider.footprint_width - SLIDER_READOUT_GAP - readout_width,
                slider.height(),
            )
            slider.place_readout()

        return None

    def sliderChange(self, change):
        """
        Keep the value shown up to date with every change of value, whether or
        not the slider's signals are blocked, and the menu's sliders fitted to
        the values they can reach.
        """

        super().sliderChange(change)
        if change == QtWidgets.QAbstractSlider.SliderValueChange:
            self.show_value(self.value())
        elif change == QtWidgets.QAbstractSlider.SliderRangeChange:
            self.fit_group()

    def place_readout(self):
        """Put the value just to the right of the slider."""

        self.readout.setGeometry(
            self.x() + self.width() + SLIDER_READOUT_GAP,
            self.y(),
            self.readout_width,
            self.height(),
        )

        return None

    def moveEvent(self, event):
        super().moveEvent(event)
        self.place_readout()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.place_readout()

    # placed on these too, as a hidden slider is only sent the events for being
    # moved and resized once it is shown - and the menus are moved to their
    # plots while hidden (see Dashboard.update_qt_elements_geometry())
    def setGeometry(self, *args):
        super().setGeometry(*args)
        self.place_readout()

    def setVisible(self, visible):
        super().setVisible(visible)
        self.place_readout()
        if visible:
            self.refresh_readout()
        self.readout.setVisible(visible)

    def raise_(self):
        super().raise_()
        self.readout.raise_()


class SettingsMenu(object):
    def __init__(self, plot_type, canvas_instance):
        """
        Initialise object to create plot settings menu

        Parameters
        ----------
        plot_type : str
            Plot type
        canvas_instance : instance
            Canvas instance
        """

        self.canvas_instance = canvas_instance

        self.elements = list(settings_dict[plot_type].keys())
        self.buttons = {}
        self.labels = {}
        self.comboboxes = {}
        self.checkable_comboboxes = {}
        self.sliders = {}
        self.checkboxes = {}
        self.lineedits = {}
        # a plot type can define more than one container (e.g. a nested
        # sub-menu's own background panel, positioned separately from the
        # main one) - self.container keeps pointing at the last one built,
        # for the common single-container case
        self.containers = {}

        for element_name in self.elements:
            element_settings = settings_dict[plot_type][element_name]
            element_type = element_settings["element_type"]
            if element_type in [
                "button",
                "container",
                "label",
                "combobox",
                "checkable_combobox",
                "slider",
                "checkbox",
                "lineedit",
            ]:
                # Add element
                element = getattr(self, "add_" + element_type)(element_settings)

                # Add options as items to options combobox
                if element_name == "options":
                    if plot_type in [
                        "periodic_violin",
                        "fairmode_target",
                        "fairmode_statsummary",
                    ]:
                        plot_type_corr = plot_type.replace("_", "-")
                    else:
                        plot_type_corr = plot_type
                    element.addItems(
                        self.canvas_instance.plot_characteristics[plot_type_corr][
                            "plot_options"
                        ]
                    )

                # Apply common properties
                if "relative_position" in element_settings.keys():
                    element.move(
                        self.buttons["settings_button"].x()
                        + element_settings["relative_position"][0],
                        self.buttons["settings_button"].y()
                        + element_settings["relative_position"][1],
                    )
                if "size" in element_settings.keys():
                    element.resize(
                        element_settings["size"][0], element_settings["size"][1]
                    )
                    # a slider shares its width with the value shown beside it -
                    # see the fitting once every element has been added, below
                    if element_type == "slider":
                        element.footprint_width = element_settings["size"][0]
                if "fixed_width" in element_settings.keys():
                    element.setFixedWidth(element_settings["fixed_width"])
                if "style" in element_settings.keys():
                    element.setStyleSheet(element_settings["style"])
                if "object_name" in element_settings.keys():
                    element.setObjectName(element_settings["object_name"])

                # Hide element
                element.hide()

                # Save element in corresponding dictionary
                if element_type == "container":
                    self.container = element
                    self.containers[element_name] = element
                elif element_type in ["button", "label", "slider", "lineedit"]:
                    getattr(self, element_type + "s")[element_name] = element
                elif element_type in ["combobox", "checkable_combobox", "checkbox"]:
                    getattr(self, element_type + "es")[element_name] = element

            else:
                error = f"Error: Unknown element type: {element_type}"
                self.canvas_instance.read_instance.logger.error(error)

        # the menu's sliders are fitted together, so they all end at the same
        # point - see ValueSlider.fit_group()
        sliders = list(self.sliders.values())
        for slider in sliders:
            slider.group = sliders
        if sliders:
            sliders[0].fit_group()

    def add_button(self, element_settings):
        """
        Add button as settings menu

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QPushButton
            Button
        """

        button = set_formatting(
            QtWidgets.QPushButton(self.canvas_instance),
            formatting_dict[element_settings["formatting_dict"]],
        )
        # icon-only (the plot corner gear/save buttons) or text-only (e.g.
        # a nested sub-menu's nav button) - support either
        if "path" in element_settings.keys():
            button.setIcon(QtGui.QIcon(join(CURRENT_PATH, element_settings["path"])))
            button.setIconSize(
                QtCore.QSize(element_settings["size"][0], element_settings["size"][1])
            )
        if "text" in element_settings.keys():
            button.setText(element_settings["text"])
        button.clicked.connect(partial(self.connect, element_settings["function"]))

        # deliberately not click-focusable - taking focus on click makes a
        # field being edited commit during the press, so the redraw that
        # triggers runs before the release and the click never completes.
        # SettingsMenu.connect() applies any pending edit anyway
        button.setFocusPolicy(QtCore.Qt.NoFocus)

        return button

    def add_container(self, element_settings):
        """
        Add elements container

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QWidget
            Container
        """

        container = set_formatting(
            QtWidgets.QWidget(self.canvas_instance),
            formatting_dict[element_settings["formatting_dict"]],
        )
        # clicking the panel background takes focus, so a line edit being
        # edited inside it commits. Without this the background is
        # focus-transparent and clicking off a field changed nothing
        container.setFocusPolicy(QtCore.Qt.ClickFocus)
        container.raise_()

        return container

    def add_label(self, element_settings):
        """
        Add label

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QLabel
            Label
        """

        if "formatting_dict" in element_settings.keys():
            label = set_formatting(
                QtWidgets.QLabel(element_settings["text"], self.canvas_instance),
                formatting_dict[element_settings["formatting_dict"]],
            )
        else:
            label = QtWidgets.QLabel(element_settings["text"], self.canvas_instance)

        return label

    def add_combobox(self, element_settings):
        """
        Add combobox

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QComboBox
            Combobox
        """

        combobox = set_formatting(
            ComboBox(self.canvas_instance),
            formatting_dict[element_settings["formatting_dict"]],
        )
        combobox.currentTextChanged.connect(
            partial(self.connect, element_settings["function"])
        )

        return combobox

    def add_slider(self, element_settings):
        """
        Add slider

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        ValueSlider
            Slider
        """

        slider = ValueSlider(
            self.canvas_instance,
            value_divisor=element_settings.get("value_divisor", 1),
            unit=element_settings.get("value_unit"),
            get_timestep=partial(get_smoothing_timestep, self.canvas_instance),
        )
        slider.setTracking(element_settings["tracking"])
        slider.setTickInterval(element_settings["tick_interval"])
        if "minimum" in element_settings.keys():
            slider.setMinimum(int(element_settings["minimum"]))
        if "maximum" in element_settings.keys():
            slider.setMaximum(int(element_settings["maximum"]))
        if "value" in element_settings.keys():
            slider.setValue(int(element_settings["value"]))
        slider.valueChanged.connect(partial(self.connect, element_settings["function"]))

        return slider

    def add_checkable_combobox(self, element_settings):
        """
        Add checkable combobox

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QComboBox
            Combobox with options to check
        """

        checkable_combobox = set_formatting(
            CheckableComboBox(self.canvas_instance),
            formatting_dict[element_settings["formatting_dict"]],
        )
        checkable_combobox.currentTextChanged.connect(
            partial(self.connect, element_settings["function"])
        )

        return checkable_combobox

    def add_checkbox(self, element_settings):
        """
        Add checkbox

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QCheckbox
            Combobox
        """

        combobox = set_formatting(
            QtWidgets.QCheckBox(element_settings["text"], self.canvas_instance),
            formatting_dict[element_settings["formatting_dict"]],
        )
        combobox.stateChanged.connect(
            partial(self.connect, element_settings["function"])
        )

        # deliberately not click-focusable - taking focus on click makes a
        # field being edited commit during the press, so the redraw that
        # triggers runs before the release and the click never completes.
        # SettingsMenu.connect() applies any pending edit anyway
        combobox.setFocusPolicy(QtCore.Qt.NoFocus)

        return combobox

    def add_lineedit(self, element_settings):
        """
        Add line edit

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QLineEdit
            Line edit
        """

        lineedit = set_formatting(
            MenuLineEdit(self.canvas_instance),
            formatting_dict[element_settings["formatting_dict"]],
        )
        if "placeholder" in element_settings.keys():
            lineedit.setPlaceholderText(element_settings["placeholder"])
        # committed (Enter, or clicking away), not textChanged - typing
        # shouldn't redraw the map on every keystroke. MenuLineEdit's own
        # signal rather than editingFinished, which Qt suppresses on
        # focus-out when a validator considers the text intermediate
        lineedit.committed.connect(
            partial(self.connect, element_settings["function"])
        )

        return lineedit

    def get_elements(self):
        """
        Get elements inside menu settings

        Returns
        -------
        list
            Menu settings elements
        """

        sliders = list(self.sliders.values())
        comboboxes = list(self.comboboxes.values())
        labels = list(self.labels.values())
        checkable_comboboxes = list(self.checkable_comboboxes.values())
        checkboxes = list(self.checkboxes.values())
        lineedits = list(self.lineedits.values())

        return (
            list(self.containers.values())
            + sliders
            + comboboxes
            + labels
            + checkable_comboboxes
            + checkboxes
            + lineedits
        )

    def connect(self, function):
        """
        Connect element to functions in settings dictionary

        Parameters
        ----------
        function : str
            Function name
        """

        if hasattr(self.canvas_instance, "interactive_elements"):
            # Call function only after all elements have been added
            if self.canvas_instance.interactive_elements.keys() == settings_dict.keys():
                # apply any field holding an uncommitted value before this
                # control's handler runs, so a value typed and then abandoned
                # is never dropped. The commit calls each field's handler
                # directly, so this cannot re-enter here
                commit = getattr(
                    self.canvas_instance, "commit_map_pending_edits", None
                )
                if commit is not None:
                    commit()

                getattr(self.canvas_instance, function)()

        return None
