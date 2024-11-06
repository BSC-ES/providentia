""" Functions to create and format PyQT elements"""


from functools import partial
import json
import os
import platform
import yaml

import numpy as np
from PyQt5 import QtCore, QtWidgets, QtGui
from textwrap import wrap

from providentia.auxiliar import CURRENT_PATH

PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
# get operating system specific formatting
operating_system = platform.system()
if operating_system == 'Darwin':
    formatting_dict = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/internal/stylesheet_mac.yaml')))
elif operating_system == 'Linux':
    formatting_dict = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/internal/stylesheet_linux.yaml')))
elif operating_system in ['Windows','MINGW32_NT','MINGW64_NT']:
    formatting_dict = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/internal/stylesheet_windows.yaml')))


def set_formatting(PyQt5_obj, format):
    """ Function that takes a PyQt5 object and applies some defined formatting. """

    # initialise style
    defined_style = ""

    # iterate through formatting dictionary and apply defined font modifiers/object formatting values
    for format_name, format_val in format.items():
        if format_name == 'font':
            
            defined_font = QtGui.QFont(format['font']['style'], int(format['font']['size']))

            if 'bold' in format[format_name].keys():
                defined_font.setBold(format['font']['bold'])
            
            if 'italic' in format[format_name].keys():
                defined_font.setItalic(format['font']['italic'])

            if 'underline' in format[format_name].keys():
                defined_font.setUnderline(format['font']['underline'])

            PyQt5_obj.setFont(defined_font)

            if hasattr(PyQt5_obj, 'lineEdit'):
                PyQt5_obj.lineEdit().setFont(defined_font)

        elif format_name == 'height':
            PyQt5_obj.setFixedHeight(format_val)

        elif format_name == 'width':
            PyQt5_obj.setFixedWidth(format_val)

        elif format_name == 'min_height':
            PyQt5_obj.setMinimumHeight(format_val)

        elif format_name == 'min_width':
            PyQt5_obj.setMinimumWidth(format_val)

        elif format_name == 'max_height':
            PyQt5_obj.setMaximumHeight(format_val)

        elif format_name == 'max_width':
            PyQt5_obj.setMaximumWidth(format_val)

        elif format_name == 'color':
            defined_style += "color: {};".format(format_val)

        elif format_name == 'background-color':
            defined_style += "background-color: {};".format(format_val)

        elif format_name == 'selection-color':
            defined_style += "selection-color: {};".format(format_val)

        elif format_name == 'selection-background-color':
            defined_style += "selection-background-color: {};".format(format_val)

        elif format_name == 'border':
            if format_val is not None:
                defined_style += "border: {}px solid {};".format(format['border']['size'], 
                                                                 format['border']['colour'])
            else:
                defined_style += "border: None;"

        elif format_name == 'border-radius':
            defined_style += "border-radius: {}px;".format(format_val)

        elif format_name == 'padding':
            defined_style += "padding: {}px;".format(format_val)

    # apply style sheet
    PyQt5_obj.setStyleSheet(defined_style)

    return PyQt5_obj


def wrap_tooltip_text(tooltip_text, max_width):
    """ Function which takes the text for a tooltip and wraps it by the screen pixel width.
        It does this by estimating the pixel width of the tooltip text (as formatted),
        and then gets the ratio exceedance over the screen pixel width.
        If there is an exceedance (i.e. > 1), the text is then broken into n max_char pieces
        based on the position of the first exceedance in the text
        (i.e. the part of the text which first exceeds the screen pixel width).
    """

    tooltip_label = set_formatting(QtWidgets.QLabel(text=tooltip_text), formatting_dict['tooltip'])
    tooltip_width = tooltip_label.fontMetrics().boundingRect(tooltip_label.text()).width()
    if tooltip_width > max_width:
        ratio = tooltip_width/max_width
        max_char = int(np.floor((len(tooltip_text)/ratio)*1.0))
        tooltip_text = '\n'.join(wrap(tooltip_text, max_char))

    return tooltip_text


def center(window):
    # Reference: https://wiki.qt.io/How_to_Center_a_Window_on_the_Screen

    window.setGeometry(
        QtWidgets.QStyle.alignedRect(
            QtCore.Qt.LeftToRight,
            QtCore.Qt.AlignCenter,
            window.size(),
            QtWidgets.qApp.desktop().availableGeometry(),
        )
    )

            
class ComboBox(QtWidgets.QComboBox):
    """ Modify default class of PyQT5 combobox. """

    def __init__(self, parent=None):

        super(ComboBox, self).__init__(parent)

        # setMaxVisibleItems only works if the box is editable
        # this creates a line edit that we need to overwrite
        self.setEditable(True)
        #self.AdjustToContents
        self.setSizeAdjustPolicy(self.AdjustToMinimumContentsLengthWithIcon)
        self.setMaxVisibleItems(20)
        self.AdjustToContents

        # overwrite default line edit by an invisible one
        self.lineEdit().setFrame(False)
        self.lineEdit().setReadOnly(True)
        self.currentTextChanged.connect(self.fixCursorPosition)

    def fixCursorPosition(self):
        """ Move (invisible) cursor to first position to avoid cutting off the start. """

        # apply only to comboboxes with text lengths of more than 8 chars
        if len(self.lineEdit().text()) >= 8:
            self.lineEdit().setCursorPosition(0)
            self.lineEdit().setFocus()
    
    def showPopup(self):
        """ Show pop-up. """

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
    """ Create combobox with multiple selection options.
        Reference: https://gis.stackexchange.com/questions/350148/qcombobox-multiple-selection-pyqt5. 
    """
    
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        # make the combo editable to set a custom text, but readonly
        self.setEditable(True)
        self.setMaxVisibleItems(20)
        self.lineEdit().setPlaceholderText('Select option/s:')
        self.lineEdit().setReadOnly(True)
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
        """ Move (invisible) cursor to first position to avoid cutting off the start. """

        # apply only to comboboxes with text lengths of more than 8 chars
        if len(self.lineEdit().text()) >= 8:
            self.lineEdit().setCursorPosition(0)
            self.lineEdit().setFocus()
    
    def resizeEvent(self, event):

        # recompute text to elide as needed
        self.updateText()
        super().resizeEvent(event)

    def eventFilter(self, object, event):

        if object == self.lineEdit():
            if event.type() == QtCore.QEvent.MouseButtonRelease:
                if self.closeOnLineEditClick:
                    self.hidePopup()
                else:
                    self.showPopup()
                return True
            return False

        if object == self.view().viewport():
            if event.type() == QtCore.QEvent.MouseButtonRelease:
                index = self.view().indexAt(event.pos())
                item = self.model().item(index.row())
                if item.checkState() == QtCore.Qt.Checked:
                    item.setCheckState(QtCore.Qt.Unchecked)
                else:
                    item.setCheckState(QtCore.Qt.Checked)
                return True

        return False

    def showPopup(self):

        super().showPopup()

        # increase the width of the elements on popup so they can be read
        self.view().setMinimumWidth(self.view().sizeHintForColumn(0) + 10)

        # when the popup is displayed, a click on the lineedit should close it
        self.closeOnLineEditClick = True

    def hidePopup(self):
        
        super().hidePopup()
        
        # used to prevent immediate reopening when clicking on the lineEdit
        self.startTimer(100)
        
        # refresh the display text when closing
        self.updateText()

    def timerEvent(self, event):
        
        # after timeout, kill timer, and reenable click on line edit
        self.killTimer(event.timerId())
        self.closeOnLineEditClick = False

    def updateText(self):
        
        texts = []
        for i in range(self.model().rowCount()):
            if self.model().item(i).checkState() == QtCore.Qt.Checked:
                texts.append(self.model().item(i).text())
        text = ", ".join(texts)

        # compute elided text (with "...")
        # metrics = QtGui.QFontMetrics(self.lineEdit().font())
        # elidedText = metrics.elidedText(text, QtCore.Qt.ElideRight, self.lineEdit().width())

        self.lineEdit().setText(text)

    def addItem(self, text, data=None):

        item = QtGui.QStandardItem()
        item.setText(text)
        if data is None:
            item.setData(text)
        else:
            item.setData(data)
        item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsUserCheckable)
        item.setData(QtCore.Qt.Unchecked, QtCore.Qt.CheckStateRole)
        self.model().appendRow(item)

    def addItems(self, texts, datalist=None):

        for i, text in enumerate(texts):
            try:
                data = datalist[i]
            except (TypeError, IndexError):
                data = None
            self.addItem(text, data)

    def currentData(self, all=False):
        
        # return the list of selected items data
        res = []
        for i in range(self.model().rowCount()):
            if not all:
                if self.model().item(i).checkState() == QtCore.Qt.Checked:
                    res.append(self.model().item(i).data())
            else:
                res.append(self.model().item(i).data())

        return res


class QVLine(QtWidgets.QFrame):
    """ Define class that generates vertical separator line. """

    def __init__(self, parent=None):
        super(QVLine, self).__init__(parent)
        self.setFrameShape(QtWidgets.QFrame.VLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)


class Switch(QtWidgets.QPushButton):
    """ Define class that generates switch buttons. """

    def __init__(self, parent=None):
        super(Switch, self).__init__(parent)
        self.setCheckable(True)

    def paintEvent(self, event):

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
        painter.drawRoundedRect(QtCore.QRect(-width, -radius, 2*width, 2*radius), radius, radius)

        # set colours and labels on switch
        label = "ON" if self.isChecked() else "OFF"
        bg_colour = QtCore.Qt.black if self.isChecked() else QtCore.Qt.gray
        text_colour = QtCore.Qt.white if self.isChecked() else QtCore.Qt.black

        # set switch background color
        painter.setBrush(QtGui.QBrush(bg_colour))

        # remove switch border color
        painter.setPen(QtGui.QPen(QtCore.Qt.NoPen))

        # change position depending on check
        sw_rect = QtCore.QRect(-radius, -radius, width + radius, 2*radius)
        if not self.isChecked():
            sw_rect.moveLeft(-width)
        painter.drawRoundedRect(sw_rect, radius, radius)

        # add label (ON / OFF)
        painter.setPen(QtGui.QPen(text_colour))
        painter.drawText(sw_rect, QtCore.Qt.AlignCenter, label)


class MessageBox(QtWidgets.QWidget):

    def __init__(self, msg, parent=None):

        super().__init__(parent)
        msg_box = self.create_msg_box(msg)
        if msg_box is not None:
            layout = QtWidgets.QVBoxLayout(self)
            layout.addWidget(msg_box)
            center(self)

    def create_msg_box(self, msg):

        # add warning box
        msg_box = QtWidgets.QMessageBox()
        msg_box.setWindowTitle("Warning")
        msg_box.setText(msg)

        # add ok button
        ok_button = set_formatting(QtWidgets.QPushButton("OK"), formatting_dict['button_popup'])
        msg_box.addButton(ok_button, QtWidgets.QMessageBox.AcceptRole)

        # create wrapper to center
        wrapper = partial(center, msg_box)
        QtCore.QTimer.singleShot(0, wrapper)
        msg_box.exec_()


class InputDialog(QtWidgets.QWidget):

    def __init__(self, read_instance, title, msg, options, parent=None):

        super().__init__(parent)
        
        dialog = self.create_dialog_box(read_instance, title, msg, options)
        if dialog is not None:
            layout = QtWidgets.QVBoxLayout(self)
            layout.addWidget(dialog)

    def create_dialog_box(self, read_instance, title, msg, options):

        dialog = QtWidgets.QInputDialog(self)
        self.selected_option, self.okpressed = dialog.getItem(read_instance, title, msg, options, 0, False)
        if not self.okpressed:
            return
