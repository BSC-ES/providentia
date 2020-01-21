""" Auxiliary classes for Dashboard """

from functools import partial
import numpy as np
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui

class ComboBox(QtWidgets.QComboBox):
    """Modify default class of PyQT5 combobox to always dropdown from fixed
    position box postion, stopping truncation of data"""

    def show_popup(self):
        """Shows popups"""
        QtWidgets.QComboBox.show_popup(self)
        self.view().parent().move(self.mapToGlobal(QtCore.QPoint()))


class QVLine(QtWidgets.QFrame):
    """Define class that generates vertical separator line"""

    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QtWidgets.QFrame.VLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)


class QHLine(QtWidgets.QFrame):
    """Define class that generates horizontal separator line"""

    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QtWidgets.QFrame.HLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)


class PopUpWindow(QtWidgets.QWidget):
    """Define class that generates generalised pop-up window"""

    def __init__(self, window_type='', window_titles=[], checkbox_labels=[],
                 default_checkbox_selection=[], selected_indices={}):
        super(PopUpWindow, self).__init__()

        # add input arguments to self
        self.window_type = window_type
        self.window_titles = window_titles
        self.checkbox_labels = checkbox_labels
        self.default_checkbox_selection = default_checkbox_selection
        self.selected_indices = selected_indices

        # create UI
        self.init_ui()

    def init_ui(self):
        """Initialise user interface"""

        # set window title
        self.setWindowTitle(self.window_type)

        # get pop up window dimensions
        # window_width = self.width()
        window_height = self.height()

        # create parent layout to hold N horizontally laid out windows
        parent_layout = QtWidgets.QHBoxLayout()
        parent_layout.setAlignment(QtCore.Qt.AlignTop)

        # get N of nested windows from length of window_titles list
        n_nested_windows = len(self.window_titles)

        # create frame to hold checkboxes
        self.checkboxes = [[] for _ in range(len(self.checkbox_labels))]

        # iterate through N nested windows, creating each one and placing in
        # parent frame accordingly
        for nested_window_n in range(n_nested_windows):

            # create nested parent layout to pull together a title, button
            # row, and grid of checkboxes
            nested_parent_layout = QtWidgets.QVBoxLayout()
            nested_parent_layout.setAlignment(QtCore.Qt.AlignTop)

            # define spacing/margin variables
            nested_parent_layout.setSpacing(10)
            nested_parent_layout.setContentsMargins(0, 0, 0, 0)

            # create title label
            title_label = \
                QtWidgets.QLabel(self,
                                 text=self.window_titles[nested_window_n])
            title_label.setAlignment(QtCore.Qt.AlignCenter)
            my_font = QtGui.QFont()
            my_font.setPointSize(16)
            title_label.setFont(my_font)

            # create row of buttons
            button_row = QtWidgets.QHBoxLayout()
            button_row.setAlignment(QtCore.Qt.AlignLeft)
            select_all_button = QtWidgets.QPushButton("Select All")
            select_all_button.setMinimumWidth(100)
            select_all_button.setMaximumWidth(100)
            clear_all_button = QtWidgets.QPushButton("Clear All")
            clear_all_button.setMinimumWidth(100)
            clear_all_button.setMaximumWidth(100)
            select_default_button = QtWidgets.QPushButton("Select Default")
            select_default_button.setMinimumWidth(120)
            select_default_button.setMaximumWidth(120)

            # order buttons in grid layout
            button_row.addWidget(select_all_button)
            button_row.addWidget(clear_all_button)
            button_row.addWidget(select_default_button)

            # add connectivity to buttons
            select_all_button.clicked.connect(partial(self.select_all, nested_window_n))
            clear_all_button.clicked.connect(partial(self.clear_all, nested_window_n))
            select_default_button.clicked.connect(partial(self.select_all_default, nested_window_n))

            # create grid of checkboxes
            checkbox_grid = QtWidgets.QGridLayout()

            # define spacing/margin variables
            checkbox_grid.setHorizontalSpacing(1)
            checkbox_grid.setVerticalSpacing(1)
            # checkbox_grid.setContentsMargins(0,0,0,0)

            # create checkboxes
            # force a new column to be started if the available vertical space
            # for each row in grid goes below a critical value (18.2) if
            # checkbox has been previously selected (without updating
            # network/resolution/species), then reselect it again
            row_n = 0
            column_n = 0
            current_selected_indices = \
                self.selected_indices[self.window_type][nested_window_n]

            for checkbox_index, val in enumerate(self.checkbox_labels[nested_window_n]):
                # checkbox_grid.setMaximumHeight(30)
                row_available_space = window_height/(row_n+1)
                if row_available_space < 15:
                    column_n += 2
                    row_n = 0
                self.checkboxes[nested_window_n].append(
                    QtWidgets.QCheckBox(val))
                if checkbox_index in current_selected_indices:
                    self.checkboxes[nested_window_n][checkbox_index].setCheckState(
                        QtCore.Qt.Checked)
                checkbox_grid.addWidget(self.checkboxes[nested_window_n][checkbox_index],
                                        row_n, column_n)
                # checkbox_grid.addWidget(self.checkboxes[nested_window_n][checkbox_index])
                row_n += 1

            # position title, button row and checkbox grid in nested parent layout

            # add title to nested parent frame
            nested_parent_layout.addWidget(title_label)

            # add button row to nested parent frame
            nested_parent_layout.addLayout(button_row)

            # add checkbox grid to nested parent frame
            nested_parent_layout.addLayout(checkbox_grid)

            # add nested parent layout to parent frame
            parent_layout.addLayout(nested_parent_layout)

            # add vertical separation line (if not last nested window that is
            # being iterated through)
            if (nested_window_n+1) != n_nested_windows:
                parent_layout.addWidget(QVLine())

        # set finalised layout
        self.setLayout(parent_layout)

        # maximise window to fit screen
        self.showMaximized()

        # setup event to get selected checkbox indices when closing window
        quit_event = QtWidgets.QAction("Quit", self)
        quit_event.triggered.connect(self.closeEvent)

    def select_all(self, nested_window_n):
        """function to select all checkboxes"""
        for checkbox_index, _ in enumerate(self.checkboxes[nested_window_n]):
            self.checkboxes[nested_window_n][checkbox_index].setCheckState(
                QtCore.Qt.Checked)

    def clear_all(self, nested_window_n):
        """function to clear all checkboxes"""
        for checkbox_index, _ in enumerate(self.checkboxes[nested_window_n]):
            self.checkboxes[nested_window_n][checkbox_index].setCheckState(QtCore.Qt.Unchecked)

    def select_all_default(self, nested_window_n):
        """function to select all default selected checkboxes"""
        # unselect all checkboxes first
        for checkbox_index, _ in enumerate(self.checkboxes[nested_window_n]):
            self.checkboxes[nested_window_n][checkbox_index].setCheckState(QtCore.Qt.Unchecked)

        # now select only desired default checkboxes
        for checkbox_index in self.default_checkbox_selection[nested_window_n]:
            self.checkboxes[nested_window_n][checkbox_index].setCheckState(QtCore.Qt.Checked)


    def closeEvent(self, event):
        """Function to get indices of selected checkboxes upon closing of pop-up window"""

        selected_indices = []

        for checkbox_index in range(len(self.checkboxes)):
            checked_indices = np.array([], dtype=np.uint8)
            for checkbox_index2 in range(len(self.checkboxes[checkbox_index])):
                if self.checkboxes[checkbox_index][checkbox_index2].checkState() == QtCore.Qt.Checked:
                    checked_indices = np.append(checked_indices, checkbox_index2)
            selected_indices.append(checked_indices)

        self.selected_indices[self.window_type] = selected_indices
