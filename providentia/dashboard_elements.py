""" Functions and classes to create and format dashboard PyQt elements """

import copy
from functools import partial
import platform
from textwrap import wrap

import numpy as np
from PyQt5 import QtCore, QtWidgets, QtGui
import yaml

from providentia.auxiliar import CURRENT_PATH, join

PROVIDENTIA_ROOT = "/".join(CURRENT_PATH.split("/")[:-1])
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


def set_formatting(
    PyQt5_obj, format, valid_obj=None, disabled=False, extra_arguments={}
):
    """
    Function that takes a PyQt5 object and applies some defined formatting

    Parameters
    ----------
    PyQt5_obj : object
        PyQt5 element
    format : dict
        Format dictionary
    valid_obj : list
        PyQt5 element to format
    disabled : bool
        Whether we want to format the disabled version of the object
    extra_arguments : list
        Extra arguments

    Returns
    -------
    object
        PyQt5 element with new format
    """

    # initialise style
    full_defined_style = ""

    # iterate through formatting dictionary and apply defined font modifiers/object formatting values
    for obj_type in format:
        if valid_obj:
            if obj_type not in valid_obj:
                continue

        if len(extra_arguments) > 0:
            if obj_type in extra_arguments:
                cut_extra_arguments = extra_arguments[obj_type]
            else:
                cut_extra_arguments = copy.deepcopy(extra_arguments)
        else:
            cut_extra_arguments = {}

        defined_style = ""

        for format_name, format_val in format[obj_type].items():
            if format_name in cut_extra_arguments:
                format_val = cut_extra_arguments[format_name]
                del cut_extra_arguments[format_name]

            if format_name == "height":
                PyQt5_obj.setFixedHeight(int(format_val))
            elif format_name == "width":
                PyQt5_obj.setFixedWidth(int(format_val))
            elif format_name == "min-height":
                PyQt5_obj.setMinimumHeight(int(format_val))
            elif format_name == "min-width":
                PyQt5_obj.setMinimumWidth(int(format_val))
            elif format_name == "max-height":
                PyQt5_obj.setMaximumHeight(int(format_val))
            elif format_name == "max-width":
                PyQt5_obj.setMaximumWidth(int(format_val))
            else:
                defined_style += "{}: {};".format(format_name, format_val)

        # have remaining extra arguments to add?
        if len(cut_extra_arguments) > 0:
            for format_name, format_val in cut_extra_arguments.items():
                if format_name == "height":
                    PyQt5_obj.setFixedHeight(int(format_val))
                elif format_name == "width":
                    PyQt5_obj.setFixedWidth(int(format_val))
                elif format_name == "min-height":
                    PyQt5_obj.setMinimumHeight(int(format_val))
                elif format_name == "min-width":
                    PyQt5_obj.setMinimumWidth(int(format_val))
                elif format_name == "max-height":
                    PyQt5_obj.setMaximumHeight(int(format_val))
                elif format_name == "max-width":
                    PyQt5_obj.setMaximumWidth(int(format_val))
                else:
                    defined_style += "{}: {};".format(format_name, format_val)

        if disabled:
            defined_style = "{}:disabled {{ {} }} ".format(obj_type, defined_style)
        else:
            defined_style = "{} {{ {} }} ".format(obj_type, defined_style)
        full_defined_style += defined_style

    # apply style sheet
    PyQt5_obj.setStyleSheet(full_defined_style)

    return PyQt5_obj


def wrap_tooltip_text(tooltip_text, max_width, format_type):
    """
    Function which takes the text for a tooltip and wraps it by the screen pixel width.
    It does this by estimating the pixel width of the tooltip text (as formatted),
    and then gets the ratio exceedance over the screen pixel width.
    If there is an exceedance (i.e. > 1), the text is then broken into n max_char pieces
    based on the position of the first exceedance in the text
    (i.e. the part of the text which first exceeds the screen pixel width).

    Parameters
    ----------
    tooltip_text : str
        Tooltip text
    max_width : int
        Maximum width accepted
    format_type : str
        Element type

    Returns
    -------
    str
        Updated tooltip text
    """

    tooltip_label = set_formatting(
        QtWidgets.QLabel(text=tooltip_text),
        formatting_dict[format_type],
        valid_obj=["QToolTip"],
    )
    tooltip_width = (
        tooltip_label.fontMetrics().boundingRect(tooltip_label.text()).width()
    )
    if tooltip_width > max_width:
        ratio = tooltip_width / max_width
        max_char = int(np.floor((len(tooltip_text) / ratio) * 1.0))
        tooltip_text = "\n".join(wrap(tooltip_text, max_char))

    return tooltip_text


def center(window):
    """
    Center window
    Reference: https://wiki.qt.io/How_to_Center_a_Window_on_the_Screen

    Parameters
    ----------
    window : MessageBox
        Message box
    """

    window.setGeometry(
        QtWidgets.QStyle.alignedRect(
            QtCore.Qt.LeftToRight,
            QtCore.Qt.AlignCenter,
            window.size(),
            QtWidgets.qApp.desktop().availableGeometry(),
        )
    )


class ComboBox(QtWidgets.QComboBox):
    """Modify default class of PyQT5 combobox."""

    def __init__(self, parent=None):
        """
        Initialise class

        Parameters
        ----------
        parent : object
            Dashboard instance
        """

        super(ComboBox, self).__init__(parent)

        # setMaxVisibleItems only works if the box is editable
        # this creates a line edit that we need to overwrite
        self.setEditable(True)
        # self.AdjustToContents
        self.setSizeAdjustPolicy(self.AdjustToMinimumContentsLengthWithIcon)
        self.setMaxVisibleItems(20)
        self.AdjustToContents

        # overwrite default line edit by an invisible one
        self.lineEdit().setFrame(False)
        self.lineEdit().setReadOnly(True)
        self.currentTextChanged.connect(self.fixCursorPosition)

    def fixCursorPosition(self):
        """
        Move (invisible) cursor to first position to avoid cutting off the start
        """

        # apply only to comboboxes with text lengths of more than 8 chars
        if len(self.lineEdit().text()) >= 8:
            self.lineEdit().setCursorPosition(0)
            self.lineEdit().setFocus()

    def showPopup(self):
        """
        Show pop-up
        """

        # set index of selected choice to highlight it
        text = self.lineEdit().text()
        index = self.findText(text, QtCore.Qt.MatchFixedString)
        self.setCurrentIndex(index)

        # show pop-up
        super().showPopup()

        # increase the width of the elements on popup so they can be read
        self.view().setMinimumWidth(self.view().sizeHintForColumn(0) + 10)

        # add vertical scroll bar
        self.view().setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)


class CheckableComboBox(QtWidgets.QComboBox):
    def __init__(self, *args, **kwargs):
        """
        Initialise class
        """

        super().__init__(*args, **kwargs)

        # make the combo editable to set a custom text, but readonly
        self.setEditable(True)
        self.lineEdit().setReadOnly(True)
        self.lineEdit().setPlaceholderText("Select option/s:")
        self.setMaxVisibleItems(20)
        self.currentTextChanged.connect(self.fixCursorPosition)

        # make the lineedit the same color as QComboBox
        palette = QtWidgets.QApplication.palette()
        palette.setBrush(QtGui.QPalette.Base, palette.button())
        self.lineEdit().setPalette(palette)

        # update the text when an item is toggled
        self.model().dataChanged.connect(self.updateText)

        # hide and show popup when clicking the line edit
        self.lineEdit().installEventFilter(self)
        self.closeOnLineEditClick = False

        # prevent popup from closing when clicking on an item
        self.view().viewport().installEventFilter(self)

    def fixCursorPosition(self):
        """
        Move (invisible) cursor to first position to avoid cutting off the start
        """

        # apply only to comboboxes with text lengths of more than 8 chars
        if len(self.lineEdit().text()) >= 8:
            self.lineEdit().setCursorPosition(0)
            self.lineEdit().setFocus()

    def resizeEvent(self, event):
        """
        Resize event after updating text

        Parameters
        ----------
        event : QResizeEvent
            Resize event
        """

        # recompute text to elide as needed
        self.updateText()
        super().resizeEvent(event)

    def eventFilter(self, obj, event):
        """
        Custom multi-select dropdown where clicking items checks/unchecks them while keeping the
        popup visually smooth.

        Parameters
        ----------
        obj : object
            Element
        event : QtCore.QEvent
            Event
        """

        if obj == self.lineEdit():
            if event.type() == QtCore.QEvent.MouseButtonRelease:
                if self.closeOnLineEditClick:
                    self.hidePopup()
                else:
                    self.showPopup()
                return True
            return False

        if obj == self.view().viewport():
            if event.type() == QtCore.QEvent.MouseButtonRelease:
                index = self.view().indexAt(event.pos())
                if not index.isValid():
                    return False

                item = self.model().item(index.row())
                if not (item.flags() & QtCore.Qt.ItemIsEnabled):
                    return True

                # Toggle check state
                item.setCheckState(
                    QtCore.Qt.Unchecked
                    if item.checkState() == QtCore.Qt.Checked
                    else QtCore.Qt.Checked
                )

                # Save current visual state
                row_to_restore = index.row()
                scroll_value = self.view().verticalScrollBar().value()

                # Hide now; reopen after the layout/menu has settled
                self.hidePopup()

                def _reopen_and_restore():
                    # Block repaints during reopen to avoid flicker
                    self.setUpdatesEnabled(False)
                    self.view().setUpdatesEnabled(False)

                    self.showPopup()

                    # Restore scroll & selection
                    try:
                        self.view().verticalScrollBar().setValue(scroll_value)
                    except RuntimeError:
                        pass  # view may be recreated

                    if 0 <= row_to_restore < self.model().rowCount():
                        idx = self.model().index(row_to_restore, 0)
                        if idx.isValid():
                            self.view().setCurrentIndex(idx)
                            self.view().scrollTo(idx)

                    self.view().setUpdatesEnabled(True)
                    self.setUpdatesEnabled(True)

                # Defer to next cycle so geometry changes (from text update) are applied
                QtCore.QTimer.singleShot(0, _reopen_and_restore)
                return True

        return False

    def showPopup(self):
        """
        Custom show pop up
        """

        model = self.model()

        # Find the first enabled item
        first_enabled_index = -1
        for i in range(model.rowCount()):
            item = model.item(i)
            if item.flags() & QtCore.Qt.ItemIsEnabled:
                first_enabled_index = i
                break

        # Open the popup first
        super().showPopup()

        # Highlight the first enabled item in the popup view only
        if first_enabled_index != -1:
            idx = model.index(first_enabled_index, 0)
            self.view().selectionModel().clearSelection()  # clear any existing selection
            self.view().selectionModel().setCurrentIndex(
                idx, QtCore.QItemSelectionModel.SelectCurrent
            )
            self.view().scrollTo(idx)

        # Adjust popup width
        width = self.view().sizeHintForColumn(0) + 40
        self.view().setMinimumWidth(width)
        self.closeOnLineEditClick = True

    def hidePopup(self):
        """
        Custom hide pop up
        """

        super().hidePopup()
        self.startTimer(100)
        self.updateText()

    def timerEvent(self, event):
        """
        Stop timer and disable closing the popup

        Parameters
        ----------
        event : QtCore.QEvent
            Event
        """

        self.killTimer(event.timerId())
        self.closeOnLineEditClick = False

    def updateText(self):
        """
        Show checked elements as one string in line
        """

        texts = [
            self.model().item(i).text()
            for i in range(self.model().rowCount())
            if self.model().item(i).checkState() == QtCore.Qt.Checked
        ]
        self.lineEdit().setText(", ".join(texts))

    def addItem(self, text, data=None, enabled=True):
        """
        Add a single checkable item to the model

        Parameters
        ----------
        text : str
            Display text of the item
        data : int, None
            Associated data
        enabled : bool, default True
            If False, item is disabled and grayed out.
        """

        item = QtGui.QStandardItem()
        item.setText(text)
        item.setData(data if data is not None else text)
        flags = QtCore.Qt.ItemIsUserCheckable
        if enabled:
            flags |= QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable
        item.setFlags(flags)
        item.setData(QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)

        if not enabled:
            item.setData(
                QtGui.QBrush(QtGui.QColor(150, 150, 150)), QtCore.Qt.ForegroundRole
            )

        self.model().appendRow(item)

        # Highlight the first enabled item if nothing is selected
        if self.currentIndex() == -1:
            for i in range(self.model().rowCount()):
                first_item = self.model().item(i)
                if first_item.flags() & QtCore.Qt.ItemIsEnabled:
                    self.setCurrentIndex(i)
                    break

    def addItems(self, texts, datalist=None, enabled_list=None):
        """
        Add multiple checkable items to the model

        Parameters
        ----------
        texts : list of str
            Display texts for the items
        datalist : list
            Data for each item
        enabled_list : list
            Enabled state for each item
        """

        for i, text in enumerate(texts):
            data = datalist[i] if datalist and i < len(datalist) else None
            enabled = (
                enabled_list[i] if enabled_list and i < len(enabled_list) else True
            )
            self.addItem(text, data, enabled)

    def currentData(self, all=False):
        """
        Return item data from the model

        Parameters
        ----------
        all : bool
            If True, return data for all items
            If False, return data only for checked items

        Returns
        -------
        list
            List of item data values
        """

        return [
            self.model().item(i).data()
            for i in range(self.model().rowCount())
            if all or self.model().item(i).checkState() == QtCore.Qt.Checked
        ]


class QVLine(QtWidgets.QFrame):
    """
    Define class that generates vertical separator line
    """

    def __init__(self, parent=None):
        super(QVLine, self).__init__(parent)
        self.setFrameShape(QtWidgets.QFrame.VLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)


class Switch(QtWidgets.QPushButton):
    """Define class that generates switch buttons."""

    def __init__(self, parent=None):
        """
        Initialise class

        Parameters
        ----------
        parent : object
            Dashboard instance
        """

        super(Switch, self).__init__(parent)
        self.setCheckable(True)

    def paintEvent(self, event):
        """
        Define switch properties

        Parameters
        ----------
        event : QtCore.QEvent
            Event
        """

        # set switch properties
        radius = 9
        width = 20
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.translate(self.rect().center())

        # add grey border to switch main box
        painter.setPen(QtGui.QPen(QtCore.Qt.gray))

        # set white background
        painter.setBrush(QtCore.Qt.white)
        painter.drawRoundedRect(
            QtCore.QRect(-width, -radius, 2 * width, 2 * radius), radius, radius
        )

        # set colours and labels on switch
        label = "ON" if self.isChecked() else "OFF"
        bg_colour = QtCore.Qt.black if self.isChecked() else QtCore.Qt.gray
        text_colour = QtCore.Qt.white if self.isChecked() else QtCore.Qt.black

        # set switch background color
        painter.setBrush(QtGui.QBrush(bg_colour))

        # remove switch border color
        painter.setPen(QtGui.QPen(QtCore.Qt.NoPen))

        # change position depending on check
        sw_rect = QtCore.QRect(-radius, -radius, width + radius, 2 * radius)
        if not self.isChecked():
            sw_rect.moveLeft(-width)
        painter.drawRoundedRect(sw_rect, radius, radius)

        # add label (ON / OFF)
        painter.setPen(QtGui.QPen(text_colour))
        painter.drawText(sw_rect, QtCore.Qt.AlignCenter, label)


class MessageBox(QtWidgets.QWidget):
    def __init__(self, msg, parent=None, confirmation=False):
        """
        Initialise class

        Parameters
        ----------
        msg : str
            Text on message box
        parent : object
            Dashboard instance
        confirmation : bool
            Indicates whether we want to ask a question with Yes or No as an answer
        """

        super().__init__(parent)

        self.result = False
        msg_box = self.create_msg_box(msg, confirmation)
        if msg_box is not None:
            layout = QtWidgets.QVBoxLayout(self)
            layout.addWidget(msg_box)
            center(self)

    def create_msg_box(self, msg, confirmation=False):
        """
        Create message box

        Parameters
        ----------
        msg : str
            Text on message box
        confirmation : bool
            Indicates whether we want to ask a question with Yes or No as an answer
        """

        # add warning box
        msg_box = QtWidgets.QMessageBox()
        msg_box.setWindowTitle("Warning")
        msg_box.setText(msg)

        if confirmation:
            # add yes button
            yes_button = set_formatting(
                QtWidgets.QPushButton("Yes"), formatting_dict["popup_button"]
            )

            # add no button
            no_button = set_formatting(
                QtWidgets.QPushButton("No"), formatting_dict["popup_button"]
            )

            msg_box.addButton(yes_button, QtWidgets.QMessageBox.AcceptRole)
            msg_box.addButton(no_button, QtWidgets.QMessageBox.RejectRole)

        else:
            # add ok button
            ok_button = set_formatting(
                QtWidgets.QPushButton("OK"), formatting_dict["popup_button"]
            )

            msg_box.addButton(ok_button, QtWidgets.QMessageBox.AcceptRole)

        # create wrapper to center
        wrapper = partial(center, msg_box)
        QtCore.QTimer.singleShot(0, wrapper)

        result = msg_box.exec_()

        if confirmation:
            self.result = result == QtWidgets.QMessageBox.AcceptRole


class InputDialog(QtWidgets.QWidget):
    def __init__(self, read_instance, title, msg, options, parent=None):
        """
        Initialise class

        Parameters
        ----------
        read_instance : object
            Instance of class Dashboard or Report
        title : str
            Title of input dialog
        msg : str
            Text of input dialog
        options : list
            Dialog options
        parent : object
            Dashboard instance
        """

        super().__init__(parent)

        dialog = self.create_dialog_box(read_instance, title, msg, options)
        if dialog is not None:
            layout = QtWidgets.QVBoxLayout(self)
            layout.addWidget(dialog)

    def create_dialog_box(self, read_instance, title, msg, options):
        dialog = QtWidgets.QInputDialog(self)
        self.selected_option, self.okpressed = dialog.getItem(
            read_instance, title, msg, options, 0, False
        )
        if not self.okpressed:
            return


class CheckDialog(QtWidgets.QDialog):
    def __init__(self, items):
        super().__init__()

        self.setWindowTitle("Select Items")

        layout = QtWidgets.QVBoxLayout(self)

        self.list_widget = QtWidgets.QListWidget()

        # remove focus
        self.list_widget.setFocusPolicy(QtCore.Qt.NoFocus)

        for item_text in items:
            item = QtWidgets.QListWidgetItem(item_text)

            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.Unchecked)

            self.list_widget.addItem(item)

        self.list_widget.itemClicked.connect(self.toggle_item)

        layout.addWidget(self.list_widget)

        button = QtWidgets.QPushButton("OK")
        button.clicked.connect(self.accept)
        layout.addWidget(button)

    def toggle_item(self, item):
        item.setCheckState(
            QtCore.Qt.Unchecked
            if item.checkState() == QtCore.Qt.Checked
            else QtCore.Qt.Checked
        )

    def get_checked_items(self):
        checked = []

        for i in range(self.list_widget.count()):
            item = self.list_widget.item(i)

            if item.checkState() == QtCore.Qt.Checked:
                checked.append(item.text())

        return checked


def create_custom_cursor(size=24):
    """
    Create custom loading cursor (Providentia logo).

    Parameters
    ----------
    size : int
        Size of the cursor in pixels
    """

    pix = QtGui.QPixmap(join(PROVIDENTIA_ROOT, "assets/logo.png"))
    pix = pix.scaled(
        size, size, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation
    )
    return QtGui.QCursor(pix)


def set_cursor(cursor_function, function, size=24):
    """
    Set custom cursor when performing a function that takes time to execute,
    and restore it to normal when finished.

    Parameters
    ----------
    cursor_function : str
        Name of the function which currently owns the cursor, to avoid resetting the cursor if it is already set by another function
    function : str
        Name of the function which is trying to set the cursor
    size : int
        Size of the cursor in pixels

    Returns
    -------
    str
        Function which owns thr cursor, to be used for restoring the cursor to normal when function is finished
    """

    cursor = QtWidgets.QApplication.overrideCursor()
    if cursor is None:
        pix = create_custom_cursor(size)
        QtWidgets.QApplication.setOverrideCursor(pix)
        return function
    else:
        return cursor_function


def unset_cursor(cursor_function, function):
    """
    Unset custom cursor if the function which owns the cursor is the one trying to restore it,
    to avoid restoring the cursor to normal if it is still being used by another function.

    Parameters
    ----------
    cursor_function : str
        Name of the function which currently owns the cursor, to avoid resetting the cursor if it is already set by another function
    function : str
        Name of the function which is trying to set the cursor
    """
    if cursor_function == function:
        QtWidgets.QApplication.restoreOverrideCursor()
