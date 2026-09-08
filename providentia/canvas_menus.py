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
        QtWidgets.QSlider
            Slider
        """

        slider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self.canvas_instance)
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
