""" Auxiliary classes for Dashboard """

import os
import copy
import time
import json
from textwrap import wrap
import numpy as np
import matplotlib
from matplotlib.widgets import _SelectorWidget
from matplotlib.lines import Line2D
from PyQt5 import QtCore, QtWidgets, QtGui
from functools import partial
from .read_aux import get_default_qa
from .aux import show_message

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
formatting_dict = json.load(open(os.path.join(CURRENT_PATH, 'conf/stylesheet.json')))

def set_formatting(PyQt5_obj, format):
    """ Function that takes a PyQt5 object and applies some defined formatting. """

    # initialise style
    defined_style = ""

    # iterate through formatting dictionary and apply defined font modifiers/object formatting values
    for format_name, format_val in format.items():
        if format_name == 'font':
            
            defined_font = QtGui.QFont(format['font']['style'], format['font']['size'])

            if 'bold' in format[format_name].keys():
                defined_font.setBold(format['font']['bold'])
            
            if 'italic' in format[format_name].keys():
                defined_font.setItalic(format['font']['italic'])

            if 'underline' in format[format_name].keys():
                defined_font.setUnderline(format['font']['underline'])

            PyQt5_obj.setFont(defined_font)

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

class LassoSelector(_SelectorWidget):
    """
    Selection curve of an arbitrary shape.
    For the selector to remain responsive you must keep a reference to it.
    The selected path can be used in conjunction with `~.Path.contains_point`
    to select data points from an image.
    In contrast to `Lasso`, `LassoSelector` is written with an interface
    similar to `RectangleSelector` and `SpanSelector`, and will continue to
    interact with the Axes until disconnected.
    Example usage::
        ax = plt.subplot()
        ax.plot(x, y)
        def onselect(verts):
            print(verts)
        lasso = LassoSelector(ax, onselect)
    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The parent Axes for the widget.
    onselect : function
        Whenever the lasso is released, the *onselect* function is called and
        passed the vertices of the selected path.
    useblit : bool, default: True
        Whether to use blitting for faster drawing (if supported by the
        backend). See the tutorial :doc:`/tutorials/advanced/blitting`
        for details.
    props : dict, optional
        Properties with which the line is drawn, see `matplotlib.lines.Line2D`
        for valid properties. Default values are defined in ``mpl.rcParams``.
    button : `.MouseButton` or list of `.MouseButton`, optional
        The mouse buttons used for rectangle selection.  Default is ``None``,
        which corresponds to all buttons.
    """

    def __init__(self, canvas_instance, ax, onselect, useblit=True, props=None, button=None):
        super().__init__(ax, onselect, useblit=useblit, button=button)
        self.verts = None
        props = {
            **(props if props is not None else {}),
            # Note that self.useblit may be != useblit, if the canvas doesn't
            # support blitting.
            'animated': self.useblit, 'visible': False,
        }
        line = Line2D([], [], **props)
        self.ax.add_line(line)
        self._selection_artist = line
        self.canvas_instance = canvas_instance

    def _press(self, event):
        self.verts = [self._get_data(event)]
        self._selection_artist.set_visible(True)

    def _release(self, event):
        if self.verts is not None:
            self.verts.append(self._get_data(event))
            self.onselect(self.verts)
        self._selection_artist.set_data([[], []])
        self._selection_artist.set_visible(False)
        self.verts = None
        if self.canvas_instance.map_annotation_disconnect:
            self.canvas_instance.map_annotation_event = self.canvas_instance.figure.canvas.mpl_connect('motion_notify_event', self.canvas_instance.hover_map_annotation)
            self.canvas_instance.map_annotation_disconnect = False

    def _onmove(self, event):
        if not self.canvas_instance.map_annotation_disconnect:
            self.canvas_instance.figure.canvas.mpl_disconnect(self.canvas_instance.map_annotation_event)
            self.canvas_instance.map_annotation_disconnect = True
        if self.verts is None:
            return
        self.verts.append(self._get_data(event))
        self._selection_artist.set_data(list(zip(*self.verts)))

        self.update()

class ComboBox(QtWidgets.QComboBox):
    """ Modify default class of PyQT5 combobox to always dropdown from fixed
        position box position, stopping truncation of data. """

    def __init__(self, parent=None):

        super(ComboBox, self).__init__(parent)

        # setMaxVisibleItems only works if the box is editable
        # this creates a line edit that we need to overwrite
        self.setEditable(True)
        self.setMaxVisibleItems(15)

        # overwrite default line edit by an invisible one
        self.lineEdit().setFrame(False)
        self.lineEdit().setReadOnly(True)

        # simulate enter event to change the current text
        self.keyPressEvent(QtGui.QKeyEvent(QtCore.QEvent.KeyPress, QtCore.Qt.Key_Enter, QtCore.Qt.NoModifier))
        self.activated.connect(self.highlightCurrent)

    def highlightCurrent(self):
        """ Set index of selected choice to highlight it. """

        text = self.lineEdit().text()
        index = self.findText(text, QtCore.Qt.MatchFixedString)
        self.setCurrentIndex(index)

    def showPopup(self):
        """ Show pop-up. """

        super().showPopup()

        # add vertical scroll bar
        self.view().setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)

        # increase the width of the elements on popup so they can be read
        self.view().setMinimumWidth(self.view().sizeHintForColumn(0))

class CheckableComboBox(QtWidgets.QComboBox):
    """ Create combobox with multiple selection options.
        Reference: https://gis.stackexchange.com/questions/350148/qcombobox-multiple-selection-pyqt5. 
    """
    
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        # make the combo editable to set a custom text, but readonly
        self.setEditable(True)
        self.lineEdit().setPlaceholderText('Select option/s:')
        self.lineEdit().setReadOnly(True)

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
        metrics = QtGui.QFontMetrics(self.lineEdit().font())
        elidedText = metrics.elidedText(text, QtCore.Qt.ElideRight, self.lineEdit().width())
        self.lineEdit().setText(elidedText)

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

    def currentData(self):
        
        # return the list of selected items data
        res = []
        for i in range(self.model().rowCount()):
            if self.model().item(i).checkState() == QtCore.Qt.Checked:
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

class PopUpWindow(QtWidgets.QWidget):
    """ Class that generates generalised pop-up window. """

    def __init__(self, read_instance, menu_root, menu_levels, main_window_geometry):
        super(PopUpWindow, self).__init__()

        # add input arguments to self
        self.read_instance = read_instance
        self.menu_root = menu_root
        self.menu_levels = menu_levels
        self.main_window_geometry = main_window_geometry
        self.menu_current = menu_root
        for _, menu_level in enumerate(menu_levels):
            self.menu_current = self.menu_current[menu_level]

        # generate GUI window for root page in menu
        self.generate_window()

        # define stylesheet for tooltips
        self.setStyleSheet("QToolTip { font: %spt %s}" % (formatting_dict['tooltip']['font']['size'],
                                                          formatting_dict['tooltip']['font']['style']))

    def generate_window(self):
        """ Generate GUI window for current menu level. """

        # get current menu level keys
        menu_current_keys = list(self.menu_current.keys())
        
        # set window title
        self.setWindowTitle(self.menu_current['window_title'])

        # create parent layout
        parent_layout = QtWidgets.QVBoxLayout()
        parent_layout.setAlignment(QtCore.Qt.AlignTop)

        # define spacing/margin variables
        self.layout_spacing = 10
        parent_layout.setSpacing(self.layout_spacing)
        self.page_margin = 5
        parent_layout.setContentsMargins(self.page_margin,self.page_margin,self.page_margin,self.page_margin)

        # set page title
        title_label = set_formatting(QtWidgets.QLabel(self, text=self.menu_current['page_title']), formatting_dict['title_popup'])
        title_label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignTop)

        # add title to parent frame
        parent_layout.addWidget(title_label)

        # create layout for placing buttons horizontally (aligned left)
        button_row = QtWidgets.QHBoxLayout()
        button_row.setAlignment(QtCore.Qt.AlignHCenter)
        
        # initialise variable to update if have any buttons in button row
        self.have_buttons = False

        # check if need to create home/previous button
        if len(self.menu_levels) >= 1:
            self.have_buttons = True
            
            # create home button
            root_button = set_formatting(QtWidgets.QPushButton("HOME"), formatting_dict['button_popup'])
            root_button.setStyleSheet("color: blue;")
            root_button.setToolTip('Return to home menu page')
            root_button.setFixedWidth(80)
            button_row.addWidget(root_button)
            root_button.clicked.connect(self.root_page)
            
            # create previous button
            previous_button = set_formatting(QtWidgets.QPushButton("PREVIOUS"), formatting_dict['button_popup'])
            previous_button.setStyleSheet("color: green;")
            previous_button.setToolTip('Return to previous menu page')
            previous_button.setFixedWidth(80)
            button_row.addWidget(previous_button)
            previous_button.clicked.connect(self.previous_page)

        # check if need to create checkbox selection buttons
        if 'select_buttons' in menu_current_keys:
            self.have_buttons = True
            
            # need to create "Select All" button?
            if 'all' in self.menu_current['select_buttons']:
                select_all_button = set_formatting(QtWidgets.QPushButton("Select All"), formatting_dict['button_popup'])
                select_all_button.setFixedWidth(100)
                button_row.addWidget(select_all_button)
                select_all_button.clicked.connect(self.select_all)
            
            # need to create "Clear All" button?
            if 'clear' in self.menu_current['select_buttons']:
                clear_all_button = set_formatting(QtWidgets.QPushButton("Clear All"), formatting_dict['button_popup'])
                clear_all_button.setFixedWidth(100)
                button_row.addWidget(clear_all_button)
                clear_all_button.clicked.connect(self.clear_all)
            
            # need to create "Select Default" button?
            if 'default' in self.menu_current['select_buttons']:
                select_default_button = set_formatting(QtWidgets.QPushButton("Select Default"), formatting_dict['button_popup'])
                select_default_button.setFixedWidth(100)
                button_row.addWidget(select_default_button)
                select_default_button.clicked.connect(self.select_all_default)

        # add button row to parent layout (if have some buttons)
        if self.have_buttons:
            parent_layout.addLayout(button_row)

        # create dictionary to store current page status
        self.page_memory = {}

        # gather list of all different menu types that need to be plotted
        menu_types = []
        if 'navigation_buttons' in menu_current_keys:
            menu_types.append('navigation_buttons')
        if 'rangeboxes' in menu_current_keys:
            menu_types.append('rangeboxes')
        if 'checkboxes' in menu_current_keys:
            menu_types.append('checkboxes')
        if 'multispecies' in menu_current_keys:
            menu_types.append('multispecies')

        # have at least 1 menu type?
        if len(menu_types) > 0:

            # create grid of menu types (horizontally concatenating different menu type grids)
            horizontal_parent = self.create_grid(menu_types)

            # add grid to parent layout
            parent_layout.addLayout(horizontal_parent)
 
        # set finalised layout
        self.setLayout(parent_layout)

        # set geometry to match that of main window
        self.setGeometry(self.main_window_geometry)

        # show pop-up window
        self.show()

        # setup event to get selected checkbox indices when closing window
        quit_event = QtWidgets.QAction("Quit", self)
        quit_event.triggered.connect(self.closeEvent)
        
    def create_grid(self, menu_types):
        """ Create grid for each needed checkbox/rangebox/navigation button menu types, that wrap vertically
            and concatenate them horizontally together.
        """

        # create horizontal layout to place all menu types within
        horizontal_parent = QtWidgets.QHBoxLayout()
        
        # set spacing between different menu type grids
        horizontal_parent.setSpacing(25)

        # align grids to centre and top
        horizontal_parent.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignTop)
        
        # order appearance of menu types grids in menu (from left to right)
        menu_type_order_dict = {'navigation_buttons': 1, 'multispecies': 1, 'rangeboxes': 2, 'checkboxes': 3}
        menu_types = sorted(menu_types, key=menu_type_order_dict.__getitem__)
        
        # iterate through passed menu types
        for menu_type in menu_types:

            # create empty grid
            scroll_area_content = QtWidgets.QWidget()
            grid = QtWidgets.QGridLayout(scroll_area_content)

            # align grid to center and top
            grid.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignCenter)

            # initialise indices used for indexing objects in grid
            start_row_n = 0
            row_n = 0
            column_n = 0

            # get dictionary nested inside current menu for menu type
            menu_current_type = self.menu_current[menu_type]
            
            # add button to add row in multispecies
            if menu_type == 'multispecies':
                element_label = 'ADD ROW'
                add_row_button = set_formatting(QtWidgets.QPushButton(element_label), 
                                                formatting_dict['multispecies_popup'])
                add_row_button.clicked.connect(lambda: self.add_multispecies_widgets(grid))
                grid.addWidget(add_row_button, 0, 3, QtCore.Qt.AlignCenter)
                
            # if have no labels for current menu type, continue to next menu type
            if len(menu_current_type['labels']) == 0 and menu_type != 'multispecies':
                continue

            # get keys associated with current menu type in current menu level
            current_menu_keys = list(menu_current_type.keys())
            
            # create dictionary to store variables that save page status per menu type (also set formatting dict for each row in grid, and vertical spacing)
            if menu_type == 'checkboxes':
                row_format_dict = formatting_dict['checkbox_popup']
                grid_vertical_spacing = 0
                if ('keep_selected' in current_menu_keys) & ('remove_selected' in current_menu_keys):
                    self.page_memory['checkboxes'] = {'keep_selected':[], 'remove_selected':[], 'n_column_consumed':3, 
                                                      'ordered_elements':['keep_selected', 'remove_selected'], 
                                                      'widget':[QtWidgets.QCheckBox,
                                                                QtWidgets.QCheckBox]}
                elif 'keep_selected' in current_menu_keys:
                    self.page_memory['checkboxes'] = {'keep_selected':[], 'n_column_consumed':2, 
                                                      'ordered_elements':['keep_selected'], 
                                                      'widget':[QtWidgets.QCheckBox]}
                elif 'remove_selected' in current_menu_keys:
                    self.page_memory['checkboxes'] = {'remove_selected':[], 'n_column_consumed':2, 
                                                      'ordered_elements':['remove_selected'], 
                                                      'widget':[QtWidgets.QCheckBox]}
            elif menu_type == 'rangeboxes':
                row_format_dict = formatting_dict['rangebox_popup']
                grid_vertical_spacing = 3
                if ('current_lower' in current_menu_keys) & ('current_upper' in current_menu_keys):
                    self.page_memory['rangeboxes'] = {'current_lower': [], 'current_upper': [], 'apply_selected': [],
                                                      'n_column_consumed': 4, 
                                                      'ordered_elements': ['current_lower', 
                                                                           'current_upper', 
                                                                           'apply_selected'], 
                                                      'widget':[QtWidgets.QLineEdit, 
                                                                QtWidgets.QLineEdit, 
                                                                QtWidgets.QCheckBox]}
                elif 'current_lower' in current_menu_keys:
                    self.page_memory['rangeboxes'] = {'current_lower': [], 'n_column_consumed': 2, 
                                                      'ordered_elements': ['current_lower'], 
                                                      'widget': [QtWidgets.QLineEdit,
                                                                 QtWidgets.QCheckBox]}
                elif 'current_upper' in current_menu_keys:
                    self.page_memory['rangeboxes'] = {'current_upper': [], 'n_column_consumed': 2, 
                                                      'ordered_elements': ['current_upper'], 
                                                      'widget': [QtWidgets.QLineEdit,
                                                                 QtWidgets.QCheckBox]}
            elif menu_type == 'navigation_buttons':
                row_format_dict = formatting_dict['navigation_button_popup']
                grid_vertical_spacing = 3
                self.page_memory['navigation_buttons'] = {'buttons': [], 'n_column_consumed': 1, 
                                                          'ordered_elements': ['buttons'], 
                                                          'widget': [QtWidgets.QPushButton]}

            elif menu_type == 'multispecies':
                row_format_dict = formatting_dict['multispecies_popup']
                grid_vertical_spacing = 3
                self.page_memory['multispecies'] = {'network': [], 'species': [], 'matrix': [],
                                                    'current_lower': [], 'current_upper': [], 
                                                    'current_filter_species_fill_value': [],
                                                    'apply_selected': [],
                                                    'n_column_consumed': 7, 
                                                    'ordered_elements':['network', 
                                                                        'matrix',
                                                                        'species',
                                                                        'current_lower', 
                                                                        'current_upper',
                                                                        'current_filter_species_fill_value', 
                                                                        'apply_selected'], 
                                                    'widget':[QtWidgets.QComboBox,
                                                              QtWidgets.QComboBox,
                                                              QtWidgets.QComboBox,
                                                              QtWidgets.QLineEdit, 
                                                              QtWidgets.QLineEdit, 
                                                              QtWidgets.QLineEdit, 
                                                              QtWidgets.QCheckBox
                                                             ]}

            # if have more than 1 column per label, need column headers
            if len(self.page_memory[menu_type]['ordered_elements']) > 1:
                have_column_headers = True
                if menu_type == 'multispecies':
                    start_row_n = 2
                else:
                    start_row_n = 1
            else:
                have_column_headers = False

            # set vertical/horizontal grid spacing
            grid.setHorizontalSpacing(3)
            grid.setVerticalSpacing(grid_vertical_spacing)

            # calculate currently occupied vertical space
            occupied_vertical_space_before_grid = self.page_margin + formatting_dict['title_popup']['height'] + self.layout_spacing + self.page_margin
            if self.have_buttons:
                occupied_vertical_space_before_grid += (formatting_dict['button_popup']['height'] + self.layout_spacing)
            if have_column_headers:
                occupied_vertical_space_before_grid += (formatting_dict['column_header_label_popup']['height'] + grid_vertical_spacing)
            
            # add horizontal scrollbar spacing
            occupied_vertical_space_before_grid += row_format_dict['height']*1.5

            # initialise variable for tracking available vertical space when appending rows of grid
            currently_occupied_vertical_space = copy.deepcopy(occupied_vertical_space_before_grid)

            # iterate through all grid labels
            for label_ii, label in enumerate(menu_current_type['labels']):
                
                # evaluate if all available vertical space has been consumed
                row_available_space = self.main_window_geometry.height() - currently_occupied_vertical_space
                
                # if available space <= than row height, force a new column to be started
                if row_available_space <= (row_format_dict['height']):
                    column_n+=self.page_memory[menu_type]['n_column_consumed']
                    row_n = 0
                    currently_occupied_vertical_space = copy.deepcopy(occupied_vertical_space_before_grid)

                # set page subtitle (if exists)
                if menu_type in ['checkboxes', 'rangeboxes']:
                    if 'subtitles' in current_menu_keys:
                        if label_ii in menu_current_type['subtitle_inds']:

                            # add subtitle to frame
                            subtitle_label = set_formatting(QtWidgets.QLabel(self, text=menu_current_type['subtitles'][menu_current_type['subtitle_inds'].index(label_ii)]), formatting_dict['subtitle_popup'])
                            grid.addWidget(subtitle_label, start_row_n+row_n, column_n, QtCore.Qt.AlignCenter)
                            
                            # update occupied vertical space 
                            currently_occupied_vertical_space += formatting_dict['subtitle_popup']['height']
                            
                            # iterate row_n
                            row_n +=1

                # add label to left of checkboxes / rangeboxes, and also add a tooltip (if exists)
                if menu_type in ['checkboxes', 'rangeboxes']:
                    rangebox_label = set_formatting(QtWidgets.QLabel(self, text = str(label)), 
                                                    formatting_dict['rangebox_label_popup'])
                    if menu_type == 'rangeboxes':
                        if len(menu_current_type['tooltips']) > 0:
                            rangebox_label.setToolTip(wrap_tooltip_text(menu_current_type['tooltips'][label_ii], 
                                                                        self.main_window_geometry.width()))
                    grid.addWidget(rangebox_label, start_row_n+row_n, column_n, QtCore.Qt.AlignLeft)

                # create all elements in column, per row
                for (element_ii, element), widget in zip(enumerate(self.page_memory[menu_type]['ordered_elements']),
                                                         self.page_memory[menu_type]['widget']):

                    # if menu type ==  'rangeboxes' then add 1 to element ii, because placed a label in first column
                    if menu_type in ['checkboxes', 'rangeboxes']:
                        element_label = ''
                        element_ii += 1
                    elif menu_type in ['multispecies']:
                        element_ii -= 1
                    else:
                        element_label = label

                    # append widget to page memory dictionary
                    if menu_type == 'multispecies':
                        self.page_memory[menu_type][element].append(set_formatting(widget(), row_format_dict))
                    else:
                        self.page_memory[menu_type][element].append(set_formatting(widget(element_label), 
                                                                                   row_format_dict))

                    # set checkboxes / rangeboxes to previous values
                    if menu_type in ['checkboxes', 'rangeboxes']:
                        
                        # check if checkbox is currently selected, if so select it again
                        # map checkbox value first if necessary
                        if menu_type == 'checkboxes':
                            if 'map_vars' in current_menu_keys:
                                var_to_check = menu_current_type['map_vars'][label_ii]
                            else:
                                var_to_check = copy.deepcopy(label)
                            if var_to_check in menu_current_type[element]:
                                self.page_memory[menu_type][element][label_ii].setCheckState(QtCore.Qt.Checked)
                        
                        # set rangeboxes to previous set value (if any)
                        elif menu_type == 'rangeboxes':
                            if element != 'apply_selected':
                                self.page_memory[menu_type][element][label_ii].setText(menu_current_type[element][label_ii])
                                self.page_memory[menu_type][element][label_ii].setFixedWidth(75)
                            else:
                                if 'map_vars' in current_menu_keys:
                                    var_to_check = menu_current_type['map_vars'][label_ii]
                                else:
                                    var_to_check = copy.deepcopy(label)
                                if var_to_check in menu_current_type[element]:
                                    self.page_memory[menu_type][element][label_ii].setCheckState(QtCore.Qt.Checked)
                        
                    # if menu type == navigation_buttons
                    elif menu_type == 'navigation_buttons':
                        
                        # add tooltip
                        if len(menu_current_type['tooltips']) > 0:
                            self.page_memory[menu_type][element][label_ii].setToolTip(wrap_tooltip_text(menu_current_type['tooltips'][label_ii], 
                                                                                                        self.main_window_geometry.width()))
                        
                        # add connectivity to buttons
                        self.page_memory[menu_type][element][label_ii].clicked.connect(self.open_new_page)

                    # if menu type == multispecies
                    elif menu_type == 'multispecies':
                        
                        # add connectivity to apply button
                        if element == 'apply_selected':
                            self.page_memory[menu_type][element][label_ii].setObjectName('apply_' + str(label_ii))
                            self.page_memory[menu_type][element][label_ii].clicked.connect(self.handle_filter_species_change)

                        # format and connect for changes
                        else:
                            if element in ['network', 'species', 'matrix']:
                                # format
                                self.page_memory[menu_type][element][label_ii].setFixedWidth(100)
                                self.page_memory[menu_type][element][label_ii].AdjustToContents

                                # add connectivity to options
                                self.page_memory[menu_type][element][label_ii].setObjectName('comboboxes_' + str(label_ii))
                                self.page_memory[menu_type][element][label_ii].currentTextChanged.connect(self.handle_multispecies_params_change)
                            else:
                                # add connectivity to options
                                self.page_memory[menu_type][element][label_ii].setObjectName('texts_' + str(label_ii))
                                self.page_memory[menu_type][element][label_ii].textChanged.connect(self.handle_multispecies_params_change)

                    # add element to grid (aligned left)
                    if menu_type == 'multispecies':
                        grid.addWidget(self.page_memory[menu_type][element][label_ii], 
                                       start_row_n+row_n, element_ii+1, QtCore.Qt.AlignLeft)
                    else:
                        grid.addWidget(self.page_memory[menu_type][element][label_ii], 
                                       start_row_n+row_n, column_n+element_ii, QtCore.Qt.AlignLeft)
                
                # update multispecies filtering fields for each row
                if menu_type == 'multispecies':
                    self.update_multispecies_fields(label_ii)
                    self.read_instance.multispecies_initialisation = False

                # iterate row_n
                row_n +=1

                # add row vertical space to total occupied space
                currently_occupied_vertical_space += (row_format_dict['height'] + grid_vertical_spacing)

            # get values before closing the tab
            if menu_type == 'multispecies':
                self.set_multispecies_previous_fields()

            # add column headers to menu type grid if needed
            if have_column_headers:
                for column_number in np.arange(0, column_n+1, self.page_memory[menu_type]['n_column_consumed']):
                    if menu_type == 'checkboxes':
                        texts = ['K', 'R']
                        for i, text in enumerate(texts):
                            grid.addWidget(set_formatting(QtWidgets.QLabel(self, text=text), 
                                                          formatting_dict['column_header_label_popup']), 
                                           0, column_number+i+1, QtCore.Qt.AlignCenter)
                    elif menu_type == 'rangeboxes':
                        texts = ['Min', 'Max', 'A']
                        for i, text in enumerate(texts):
                            grid.addWidget(set_formatting(QtWidgets.QLabel(self, text=text), 
                                                          formatting_dict['column_header_label_popup']), 
                                           0, column_number+i+1, QtCore.Qt.AlignCenter)
                    elif menu_type == 'multispecies':
                        if len(self.menu_current['multispecies']['labels']) > 0:
                            texts = ['Network', 'Matrix', 'Species', 'Min', 'Max', 'Fill value', 'A']
                            for i, text in enumerate(texts):
                                grid.addWidget(set_formatting(QtWidgets.QLabel(self, text=text), 
                                                              formatting_dict['column_header_label_popup']), 
                                               1, i, QtCore.Qt.AlignCenter)

            # set horizontal scroll properties
            scroll_area = QtWidgets.QScrollArea() 
            scroll_area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            scroll_area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
            scroll_area.setWidgetResizable(True)
            scroll_area.setFrameShape(0)

            # add horizontal scroll
            scroll_area.setWidget(scroll_area_content)
            horizontal_parent.addWidget(scroll_area)

        # return horizontally concatenated menu type grids
        return horizontal_parent

    def open_new_page(self):
        """ Function to open new page in pop-up window. """

        # get selected navigation button text
        selected_navigation_button = self.sender().text()
        
        # add selected navigation button text to menu levels list
        self.menu_levels.append(selected_navigation_button)
        
        # create new pop-up page for selected navigation button
        self.new_window = PopUpWindow(self.read_instance, self.menu_root, self.menu_levels, self.main_window_geometry)
        
        # sleep briefly to allow new page to be generated
        time.sleep(0.1)
        
        # close current pop-up page
        self.close()

    def root_page(self):
        """ Function that returns pop-up window to root menu level page. """

        # create new pop-up page for root menu level
        self.new_window = PopUpWindow(self.read_instance, self.menu_root, [], self.main_window_geometry)
        
        # sleep briefly to allow new page to be generated
        time.sleep(0.1)
        
        # close current pop-up page
        self.close()

    def previous_page(self):
        """ Function that returns pop-up window to previous menu level page. """

        # create new pop-up page for previous menu level
        self.new_window = PopUpWindow(self.read_instance, self.menu_root, self.menu_levels[:-1], 
                                      self.main_window_geometry)
        
        # sleep briefly to allow new page to be generated
        time.sleep(0.1)
        
        # close current pop-up page
        self.close()

    def select_all(self):
        """ Function to select all checkboxes. """

        for element in self.page_memory['checkboxes']['ordered_elements']:
            for checkbox_ii, checkbox in enumerate(self.page_memory['checkboxes'][element]):
                self.page_memory['checkboxes'][element][checkbox_ii].setCheckState(QtCore.Qt.Checked)

    def clear_all(self):
        """ Function to clear all checkboxes. """

        for element in self.page_memory['checkboxes']['ordered_elements']:
            for checkbox_ii, checkbox in enumerate(self.page_memory['checkboxes'][element]):
                self.page_memory['checkboxes'][element][checkbox_ii].setCheckState(QtCore.Qt.Unchecked)

    def select_all_default(self):
        """ Function to select all default selected checkboxes. """

        # unselect all checkboxes first
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for checkbox_ii, checkbox in enumerate(self.page_memory['checkboxes'][element]):
                self.page_memory['checkboxes'][element][checkbox_ii].setCheckState(QtCore.Qt.Unchecked)

        # map default variables to positional indices
        if 'map_vars' in list(self.menu_current['checkboxes'].keys()):
            default_inds = [np.where(self.menu_current['checkboxes']['map_vars'] == default_var)[0][0] for default_var in self.menu_current['checkboxes']['remove_default']]
        else:
            default_inds = [np.where(self.menu_current['checkboxes']['labels'] == default_var)[0][0] for default_var in self.menu_current['checkboxes']['remove_default']]

        # now select only desired default checkboxes
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for default_ind in default_inds:
                self.page_memory['checkboxes'][element][default_ind].setCheckState(QtCore.Qt.Checked)

    def add_multispecies_widgets(self, grid):
        """ Function to add new line to pop-up window. """
        
        # get widget line
        label_ii = len(self.menu_current['multispecies']['labels'])

        # initalise labels
        if len(self.menu_current['multispecies']['labels']) == 0:
            texts = ['Network', 'Matrix', 'Species', 'Min', 'Max', 'Fill value', 'A']
            for i, text in enumerate(texts):
                grid.addWidget(set_formatting(QtWidgets.QLabel(self, text=text), 
                                              formatting_dict['column_header_label_popup']), 
                                1, i, QtCore.Qt.AlignCenter)

        # update menu_current labels
        # check if they exist since they might have been added with the function multispecies_conf
        if 'networkspeci_' + str(label_ii) not in self.menu_current['multispecies']['labels']:
            self.menu_current['multispecies']['labels'].append('networkspeci_' + str(label_ii))

        # add new line only when add_row button is pressed
        for (element_ii, element), widget in zip(enumerate(self.page_memory['multispecies']['ordered_elements']),
                                                           self.page_memory['multispecies']['widget']):
            
            # add element to memory
            self.page_memory['multispecies'][element].append(set_formatting(widget(), 
                                                                            formatting_dict['multispecies_popup']))
            grid.addWidget(self.page_memory['multispecies'][element][label_ii], label_ii + 2, element_ii, 
                           QtCore.Qt.AlignLeft)

            if element == 'apply_selected':
                self.page_memory['multispecies'][element][label_ii].setObjectName('apply_' + str(label_ii))
                self.page_memory['multispecies'][element][label_ii].clicked.connect(self.handle_filter_species_change)

            else:
                if element in ['network', 'species', 'matrix']:
                    # format
                    self.page_memory['multispecies'][element][label_ii].setFixedWidth(100)
                    self.page_memory['multispecies'][element][label_ii].AdjustToContents

                    # add connectivity to options
                    self.page_memory['multispecies'][element][label_ii].setObjectName('comboboxes_' + str(label_ii))
                    self.page_memory['multispecies'][element][label_ii].currentTextChanged.connect(self.handle_multispecies_params_change)
                else:
                    # add connectivity to options
                    self.page_memory['multispecies'][element][label_ii].setObjectName('texts_' + str(label_ii))
                    self.page_memory['multispecies'][element][label_ii].textChanged.connect(self.handle_multispecies_params_change)

            # add element to grid
            grid.addWidget(self.page_memory['multispecies'][element][label_ii], label_ii + 2, element_ii, 
                           QtCore.Qt.AlignLeft)
                                    
        # update multispecies filtering fields for each row
        self.read_instance.multispecies_initialisation = True
        self.update_multispecies_fields(label_ii)
        self.read_instance.multispecies_initialisation = False

    def handle_filter_species_change(self):
        """ Function to add or remote filter species after clicking on apply button. """
        
        # get widget line
        event_source = self.sender()
        label_ii = int(event_source.objectName().split('_')[1])
        
        if event_source.isChecked():
            self.update_filter_species(label_ii)
        else:
            self.update_filter_species(label_ii, add_filter_species=False)

    def update_multispecies_fields(self, label_ii):
        """ Update multispecies fields in tab.

            :param label_ii: Corresponding widget line in dashboard
            :type label_ii: int
        """

        # set variable to block interactive handling while updating config bar parameters
        self.read_instance.block_config_bar_handling_updates = True

        # set some default configuration values when initialising config bar
        if self.read_instance.multispecies_initialisation:
            
            # set initial selected config variables as set .conf files or defaults
            self.read_instance.selected_widget_network.update({label_ii: copy.deepcopy(self.read_instance.selected_network)})
            self.read_instance.selected_widget_matrix.update({label_ii: self.read_instance.parameter_dictionary[self.read_instance.selected_species]['matrix']})
            self.read_instance.selected_widget_species.update({label_ii: copy.deepcopy(self.read_instance.selected_species)})
            self.read_instance.selected_widget_lower.update({label_ii: str(np.nan)})
            self.read_instance.selected_widget_upper.update({label_ii: str(np.nan)})
            self.read_instance.selected_widget_filter_species_fill_value.update({label_ii: str(np.nan)})
            self.read_instance.selected_widget_apply.update({label_ii: False})

        # initialise/update fields - maintain previously selected values wherever possible
        # clear fields
        if len(self.menu_current['multispecies']['labels']) > 0:
            self.page_memory['multispecies']['network'][label_ii].clear()
            self.page_memory['multispecies']['matrix'][label_ii].clear()
            self.page_memory['multispecies']['species'][label_ii].clear()
            self.page_memory['multispecies']['current_lower'][label_ii].setText(str(np.nan))
            self.page_memory['multispecies']['current_upper'][label_ii].setText(str(np.nan))
            self.page_memory['multispecies']['current_filter_species_fill_value'][label_ii].setText(str(np.nan))
            self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Unchecked)

        # if have no available observational data, return from function, updating variable informing that have no data
        if len(self.read_instance.available_observation_data) == 0:
            self.read_instance.no_data_to_read = True
            # unset variable to allow interactive handling from now
            self.read_instance.block_config_bar_handling_updates = False
            return
        else:
            self.read_instance.no_data_to_read = False
        
        # get available networks for selected resolution
        available_networks = copy.deepcopy(list(self.read_instance.available_observation_data.keys()))
        for network in list(self.read_instance.available_observation_data.keys()):
            if self.read_instance.selected_resolution not in self.read_instance.available_observation_data[network].keys():
                available_networks.remove(network)
        
        # update network and set current text in network widget
        self.page_memory['multispecies']['network'][label_ii].addItems(available_networks)
        if self.read_instance.selected_widget_network[label_ii] in available_networks:
            self.page_memory['multispecies']['network'][label_ii].setCurrentText(self.read_instance.selected_widget_network[label_ii])
        else:
            self.read_instance.selected_widget_network.update({label_ii: self.page_memory['multispecies']['network'][label_ii].currentText()})

        # add available matrices
        available_matrices = sorted(self.read_instance.available_observation_data[self.read_instance.selected_widget_network[label_ii]][self.read_instance.selected_resolution])
        self.page_memory['multispecies']['matrix'][label_ii].addItems(available_matrices)
        
        # update matrix and set current text in matrix widget
        if self.read_instance.selected_widget_matrix[label_ii] in available_matrices:
            self.page_memory['multispecies']['matrix'][label_ii].setCurrentText(self.read_instance.selected_widget_matrix[label_ii])
        else:
            self.read_instance.selected_widget_matrix.update({label_ii: self.page_memory['multispecies']['matrix'][label_ii].currentText()})
                            
        # add available species
        available_species = sorted(self.read_instance.available_observation_data[self.read_instance.selected_widget_network[label_ii]][self.read_instance.selected_resolution][self.read_instance.selected_widget_matrix[label_ii]])
        self.page_memory['multispecies']['species'][label_ii].addItems(available_species)
        
        # update species and set current text in species widget
        if self.read_instance.selected_widget_species[label_ii] in available_species:
            self.page_memory['multispecies']['species'][label_ii].setCurrentText(self.read_instance.selected_widget_species[label_ii])
        else:
            self.read_instance.selected_widget_species.update({label_ii: self.page_memory['multispecies']['species'][label_ii].currentText()})
        
        # update current lower
        if self.read_instance.selected_widget_lower[label_ii] != str(np.nan):
            self.page_memory['multispecies']['current_lower'][label_ii].setText(str(self.read_instance.selected_widget_lower[label_ii]))

        # update current upper
        if self.read_instance.selected_widget_upper[label_ii] != str(np.nan):
            self.page_memory['multispecies']['current_upper'][label_ii].setText(str(self.read_instance.selected_widget_upper[label_ii]))

        # update current fill value
        if self.read_instance.selected_widget_filter_species_fill_value[label_ii] != str(np.nan):
            self.page_memory['multispecies']['current_filter_species_fill_value'][label_ii].setText(str(self.read_instance.selected_widget_filter_species_fill_value[label_ii]))

        # update apply
        if self.read_instance.selected_widget_apply[label_ii] == True:
            self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Checked)

        self.read_instance.block_config_bar_handling_updates = False

    def set_multispecies_previous_fields(self):
        """ Function to set the same multispecies filtering fields that we previously had upon reopening the tab.
        """

        # set variable to block interactive handling while updating config bar parameters
        self.read_instance.block_config_bar_handling_updates = True
        
        # set previous bounds and check apply_selected
        # do this only when previous bounds exist (after tab has been closed)
        if len(self.menu_current['multispecies']['previous_lower'].keys()) > 0:
            for label_ii, label in enumerate(self.menu_current['multispecies']['labels']):
                if (label_ii+1) <= len(self.menu_current['multispecies']['previous_lower'].keys()):
                    self.page_memory['multispecies']['current_lower'][label_ii].setText(
                        str(self.menu_current['multispecies']['previous_lower'][label_ii]))
                    self.page_memory['multispecies']['current_upper'][label_ii].setText(
                        str(self.menu_current['multispecies']['previous_upper'][label_ii]))
                    self.page_memory['multispecies']['current_filter_species_fill_value'][label_ii].setText(
                        str(self.menu_current['multispecies']['previous_filter_species_fill_value'][label_ii]))
                    if self.menu_current['multispecies']['previous_apply'][label_ii]:
                        self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Checked)
                    else:
                        self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Unchecked)

        self.read_instance.block_config_bar_handling_updates = False

    def handle_multispecies_params_change(self, changed_param):
        """ Function which handles interactive updates of multispecies filtering fields. """
        
        if (changed_param != '') and (not self.read_instance.block_config_bar_handling_updates):
            
            # get event origin source
            event_source = self.sender()

            # get line
            label_ii = int(event_source.objectName().split('_')[1])

            # reset bounds and apply values on widgets
            self.read_instance.selected_widget_apply.update({label_ii: False})

            # remove previous networkspeci from lists
            self.update_filter_species(label_ii, add_filter_species=False)

            # if network, matrix or species have changed then respective
            # current selection for the changed param
            if event_source == self.page_memory['multispecies']['network'][label_ii]:
                self.read_instance.selected_widget_network.update({label_ii: changed_param})
                self.read_instance.selected_widget_matrix.update({label_ii: sorted(list(
                    self.read_instance.available_observation_data[self.read_instance.selected_widget_network[label_ii]][self.read_instance.selected_resolution].keys()))[0]})
            elif event_source == self.page_memory['multispecies']['matrix'][label_ii]:
                self.read_instance.selected_widget_matrix.update({label_ii: changed_param})
                self.read_instance.selected_widget_species.update({label_ii: sorted(list(
                    self.read_instance.available_observation_data[self.read_instance.selected_widget_network[label_ii]][self.read_instance.selected_resolution][self.read_instance.selected_widget_matrix[label_ii]].keys()))[0]})
            elif event_source == self.page_memory['multispecies']['species'][label_ii]:
                self.read_instance.selected_widget_species.update({label_ii: changed_param})
            elif event_source == self.page_memory['multispecies']['current_lower'][label_ii]:
                self.read_instance.selected_widget_lower.update({label_ii: changed_param})
            elif event_source == self.page_memory['multispecies']['current_upper'][label_ii]:
                self.read_instance.selected_widget_upper.update({label_ii: changed_param})
            elif event_source == self.page_memory['multispecies']['current_filter_species_fill_value'][label_ii]:
                self.read_instance.selected_widget_filter_species_fill_value.update({label_ii: changed_param})


            # update multispecies filtering fields
            self.update_multispecies_fields(label_ii)

    def closeEvent(self, event):
        """ Function to get status of current page upon closing of pop-up window. """

        # take everything from page memory dictionary and put it back into menu level object

        # iterate through menu types
        for menu_type in list(self.page_memory.keys()):
            for element in self.page_memory[menu_type]['ordered_elements']:
                if menu_type == 'checkboxes':
                    selected_vars = []
                    for checkbox_ii, checkbox in enumerate(self.page_memory[menu_type][element]):
                        if checkbox.checkState() == QtCore.Qt.Checked:
                            # map selected position index to variable name
                            if 'map_vars' in list(self.menu_current[menu_type].keys()):
                                selected_vars.append(self.menu_current[menu_type]['map_vars'][checkbox_ii])
                            else:
                                selected_vars.append(self.menu_current[menu_type]['labels'][checkbox_ii])
                    # update previous selected variable
                    if element == 'keep_selected':
                        self.menu_current[menu_type]['previous_keep_selected'] = copy.deepcopy(self.menu_current[menu_type]['keep_selected'])
                    elif element == 'remove_selected':
                        self.menu_current[menu_type]['previous_remove_selected'] = copy.deepcopy(self.menu_current[menu_type]['remove_selected'])
                    # update selected variable
                    self.menu_current[menu_type][element] = selected_vars
                
                if menu_type == 'rangeboxes':
                    set_vals = []
                    selected_vars = []
                    for rangebox_ii, rangebox in enumerate(self.page_memory[menu_type][element]):
                        if element != 'apply_selected':
                            set_vals.append(rangebox.text())
                        else:
                            if rangebox.checkState() == QtCore.Qt.Checked:
                                # map selected position index to variable name
                                if 'map_vars' in list(self.menu_current[menu_type].keys()):
                                    selected_vars.append(self.menu_current[menu_type]['map_vars'][rangebox_ii])
                                else:
                                    selected_vars.append(self.menu_current[menu_type]['labels'][rangebox_ii])
                               
                    # update previous set variable
                    if element == 'current_lower':
                        self.menu_current[menu_type]['previous_lower'] = copy.deepcopy(self.menu_current[menu_type]['current_lower'])
                    elif element == 'current_upper':
                        self.menu_current[menu_type]['previous_upper'] = copy.deepcopy(self.menu_current[menu_type]['current_upper'])
                    elif element == 'apply_selected':
                        self.menu_current[menu_type]['previous_apply'] = copy.deepcopy(self.menu_current[menu_type]['apply_selected'])
                    
                    # update set value
                    if element != 'apply_selected':
                        self.menu_current[menu_type][element] = set_vals
                    else:
                        self.menu_current[menu_type][element] = selected_vars       

                elif menu_type == 'multispecies':
                    for label_ii, label in enumerate(self.menu_current[menu_type]['labels']):
                        # update previous set variable
                        if element == 'current_lower':
                            self.menu_current[menu_type]['previous_lower'].update({label_ii: self.page_memory['multispecies']['current_lower'][label_ii].text()})
                        elif element == 'current_upper':
                            self.menu_current[menu_type]['previous_upper'].update({label_ii: self.page_memory['multispecies']['current_upper'][label_ii].text()})
                        elif element == 'current_filter_species_fill_value':
                            self.menu_current[menu_type]['previous_filter_species_fill_value'].update({label_ii: self.page_memory['multispecies']['current_filter_species_fill_value'][label_ii].text()})
                        elif element == 'apply_selected':
                            if self.page_memory['multispecies']['apply_selected'][label_ii].isChecked():
                                self.menu_current[menu_type]['previous_apply'].update({label_ii: True})
                            else:
                                self.menu_current[menu_type]['previous_apply'].update({label_ii: False})

    def update_filter_species(self, label_ii, add_filter_species=True):
        """ Function to update filter species after launching the dashboard with a configuration file or 
            by editing the fields in the multispecies filtering tab in the dashboard. 

            :param instance: Instance of class ProvidentiaMainWindow
            :type instance: object
            :param label_ii: Corresponding widget line in dashboard
            :type label_ii: int
            :param add_filter_species: boolean to indicate if networkspeci has to be added or removed
            :type add_filter_species: boolean
        """

        # update previous filter species
        self.read_instance.previous_filter_species = copy.deepcopy(self.read_instance.filter_species)

        # get selected network, species and bounds
        network = self.read_instance.selected_widget_network[label_ii]
        speci = self.read_instance.selected_widget_species[label_ii]
        networkspeci = network + '|' + speci
        current_lower = self.read_instance.selected_widget_lower[label_ii]
        current_upper = self.read_instance.selected_widget_upper[label_ii]
        current_filter_species_fill_value = self.read_instance.selected_widget_filter_species_fill_value[label_ii]

        # get filter species after changes
        current_filter_species = [current_lower, current_upper, current_filter_species_fill_value]

        # if apply button is checked or filter_species in configuration file, add networkspecies in filter_species
        if add_filter_species:
            
            # do not add to filter_species if lower and upper bounds are nan
            if current_lower == str(np.nan) or current_upper == str(np.nan):
                msg = 'Data bounds cannot be empty.'
                show_message(self.read_instance, msg)
                self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Unchecked)
                return

            # only add to filter_species when lower bound if it contains :, > or >=
            if ('<' in current_lower):
                msg = 'Lower bound ({}) for {} cannot contain < or <=. '.format(current_lower, networkspeci)
                show_message(self.read_instance, msg)
                self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Unchecked)
                return
            elif (':' not in current_lower) and ('>' not in current_lower):
                msg = 'Lower bound ({}) for {} should contain > or >=. '.format(current_lower, networkspeci)
                show_message(self.read_instance, msg)
                self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Unchecked)
                return

            # only add to filter_species when upper bound if it contains :, < or <=
            if ('>' in current_upper):
                msg = 'Upper bound ({}) for {} cannot contain > or >=. '.format(current_upper, networkspeci)
                show_message(self.read_instance, msg)
                self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Unchecked)
                return
            elif (':' not in current_upper) and ('<' not in current_upper):
                msg = 'Upper bound ({}) for {} should contain < or <=. '.format(current_upper, networkspeci)
                show_message(self.read_instance, msg)
                self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Unchecked)
                return

            # add or update networkspeci
            # check selected lower and upper bounds and fill value are numbers or nan
            try:
                if networkspeci in self.read_instance.filter_species.keys():
                    if current_filter_species not in self.read_instance.filter_species[networkspeci]:
                        self.read_instance.filter_species[networkspeci].append(current_filter_species)
                else:
                    self.read_instance.filter_species[networkspeci] = [current_filter_species]

            # if any of the fields are not numbers, return from function
            except ValueError:
                msg = 'Warning: Data limit fields must be numeric.'
                show_message(self.read_instance, msg)
                self.page_memory['multispecies']['apply_selected'][label_ii].setCheckState(QtCore.Qt.Unchecked)
                return

            # get quality flags for species if the information is not available in qa_per_species
            if speci not in self.read_instance.qa_per_species:
                # get species in memory 
                species = copy.deepcopy(self.read_instance.species)
                filter_species = [val.split('|')[1] 
                                  for val in list(copy.deepcopy(self.read_instance.filter_species).keys())]
                qa_species = species + filter_species
                
                # add
                qa_species.append(speci)
                self.read_instance.qa_per_species = {speci:get_default_qa(self.read_instance, speci) 
                                                     for speci in qa_species}

        # if apply button is unchecked, remove networkspecies from filter_species
        else:
            # remove from filter_species
            filter_species_aux = copy.deepcopy(self.read_instance.filter_species)
            if networkspeci in filter_species_aux.keys():
                for networkspeci in filter_species_aux:
                    if current_filter_species in filter_species_aux[networkspeci]:
                        sub_networkspeci_ii = self.read_instance.filter_species[networkspeci].index(current_filter_species)
                        del self.read_instance.filter_species[networkspeci][sub_networkspeci_ii]
                        if len(self.read_instance.filter_species[networkspeci]) == 0:
                            del self.read_instance.filter_species[networkspeci]

            # remove from qa_per_species
            if ((speci in self.read_instance.qa_per_species) and 
                (networkspeci not in self.read_instance.filter_species.keys())):
                del self.read_instance.qa_per_species[speci]
