""" Auxiliary classes for Dashboard """

import copy
import time
from textwrap import wrap
import numpy as np
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui


# setup dictionary characterising formats for all GUI window objects (i.e. buttons, titles etc.)
formatting_dict = {'title_menu':               {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20, 'underline':True},
                   'label_menu':               {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'button_menu':              {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'checkbox_menu':            {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'combobox_menu':            {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'lineedit_menu':            {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'title_popup':              {'font':QtGui.QFont("SansSerif", 13),  'colour':'black', 'height':23, 'underline':True},
                   'button_popup':             {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'navigation_button_popup':  {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':20},
                   'checkbox_popup':           {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':18},
                   'rangebox_label_popup':     {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':18},
                   'rangebox_popup':           {'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':18, 'max_width':90},
                   'column_header_label_popup':{'font':QtGui.QFont("SansSerif", 9.5), 'colour':'black', 'height':18, 'italic':True},
                   'tooltip':                  {'font':QtGui.QFont("SansSerif", 8)}
                   }


def set_formatting(PyQt5_obj, formats):
    """function that takes a PyQt5 object and applies some defined formatting"""

    #first get defined font for object
    defined_font = formats['font']

    #iterate through formats dictionary and apply defined font modifiers/object formatting values
    for format_name, format_val in formats.items():
        if format_name == 'bold':
            defined_font.setBold(format_val)

        if format_name == 'italic':
            defined_font.setItalic(format_val)

        elif format_name == 'underline':
            defined_font.setUnderline(format_val)

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

        elif format_name == 'colour':
            PyQt5_obj.setStyleSheet("color: {};".format(format_val))

    # now apply font to object
    PyQt5_obj.setFont(defined_font)

    # return modified PyQt5 object
    return PyQt5_obj


def wrap_tooltip_text(tooltip_text, max_width):
    """function which takes the text for a tooltip and wraps it by the screen pixel width.
       It does this by estimating the pixel width of the tooltip text (as formatted),
       and then gets the ratio exceedance over the screen pixel width.
       If there is an exceedance (i.e. > 1), the text is then broken into n max_char pieces
       based on the position of the first exceedance in the text
       (i.e. the part of the text which first exceeds the screen pixel width)
    """

    tooltip_label = set_formatting(QtWidgets.QLabel(text=tooltip_text), formatting_dict['tooltip'])
    tooltip_width = tooltip_label.fontMetrics().boundingRect(tooltip_label.text()).width()
    if tooltip_width > max_width:
        ratio = tooltip_width/max_width
        max_char = int(np.floor((len(tooltip_text)/ratio)*1.0))
        tooltip_text = '\n'.join(wrap(tooltip_text, max_char))

    return tooltip_text


class ComboBox(QtWidgets.QComboBox):
    """Modify default class of PyQT5 combobox to always dropdown from fixed
    position box postion, stopping truncation of data"""

    def showPopup(self):
        """Shows popups"""
        QtWidgets.QComboBox.showPopup(self)
        self.view().parent().move(self.mapToGlobal(QtCore.QPoint()))
        #self.view().setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)


class QVLine(QtWidgets.QFrame):
    """Define class that generates vertical separator line"""

    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QtWidgets.QFrame.VLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)


class PopUpWindow(QtWidgets.QWidget):
    """Define class that generates generalised pop-up window"""

    def __init__(self, menu_root, menu_levels, main_window_geometry):
        super(PopUpWindow, self).__init__()

        #add input arguments to self
        self.menu_root = menu_root
        self.menu_levels = menu_levels
        self.main_window_geometry = main_window_geometry
        self.menu_current = menu_root
        for _, menu_level in enumerate(menu_levels):
            self.menu_current = self.menu_current[menu_level]

        #generate GUI window for root page in menu
        self.generate_window()

        #define stylesheet for tooltips
        self.setStyleSheet("QToolTip { font: %spt %s}"%(formatting_dict['tooltip']['font'].pointSizeF(), formatting_dict['tooltip']['font'].family()))

    def generate_window(self):

        """generate GUI window for current menu level"""

        #get current menu level keys
        menu_current_keys = list(self.menu_current.keys())

        #set window title
        self.setWindowTitle(self.menu_current['window_title'])

        #create parent layout
        parent_layout = QtWidgets.QVBoxLayout()
        parent_layout.setAlignment(QtCore.Qt.AlignTop)

        #define spacing/margin variables
        self.layout_spacing = 10
        parent_layout.setSpacing(self.layout_spacing)
        self.page_margin = 5
        parent_layout.setContentsMargins(self.page_margin,self.page_margin,self.page_margin,self.page_margin)

        #set page title
        title_label = set_formatting(QtWidgets.QLabel(self, text=self.menu_current['page_title']), formatting_dict['title_popup'])
        title_label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignTop)

        #add title to parent frame
        parent_layout.addWidget(title_label)

        #create layout for placing buttons horizontally (aligned left)
        button_row = QtWidgets.QHBoxLayout()
        button_row.setAlignment(QtCore.Qt.AlignLeft)
        #initialise variable to update if have any buttons in button row
        self.have_buttons = False

        #check if need to create home/previous button
        if len(self.menu_levels) >= 1:
            self.have_buttons = True
            #create home button
            root_button = set_formatting(QtWidgets.QPushButton("HOME"), formatting_dict['button_popup'])
            root_button.setStyleSheet("color: blue;")
            root_button.setToolTip('Return to home menu page')
            root_button.setFixedWidth(80)
            button_row.addWidget(root_button)
            root_button.clicked.connect(self.root_page)
            #create previous button
            previous_button = set_formatting(QtWidgets.QPushButton("PREVIOUS"), formatting_dict['button_popup'])
            previous_button.setStyleSheet("color: green;")
            previous_button.setToolTip('Return to previous menu page')
            previous_button.setFixedWidth(80)
            button_row.addWidget(previous_button)
            previous_button.clicked.connect(self.previous_page)

        #check if need to create checkbox selection buttons
        if 'select_buttons' in menu_current_keys:
            self.have_buttons = True
            #need to create "Select All" button?
            if 'all' in self.menu_current['select_buttons']:
                select_all_button = set_formatting(QtWidgets.QPushButton("Select All"), formatting_dict['button_popup'])
                select_all_button.setFixedWidth(100)
                button_row.addWidget(select_all_button)
                select_all_button.clicked.connect(self.select_all)
            #need to create "Clear All" button?
            if 'clear' in self.menu_current['select_buttons']:
                clear_all_button = set_formatting(QtWidgets.QPushButton("Clear All"), formatting_dict['button_popup'])
                clear_all_button.setFixedWidth(100)
                button_row.addWidget(clear_all_button)
                clear_all_button.clicked.connect(self.clear_all)
            #need to create "Select Default" button?
            if 'default' in self.menu_current['select_buttons']:
                select_default_button = set_formatting(QtWidgets.QPushButton("Select Default"), formatting_dict['button_popup'])
                select_default_button.setFixedWidth(100)
                button_row.addWidget(select_default_button)
                select_default_button.clicked.connect(self.select_all_default)

        #add button row to parent layout (if have some buttons)
        if self.have_buttons:
            parent_layout.addLayout(button_row)

        #create dictionary to store current page status
        self.page_memory = {}

        #gather list of all different menu types that need to be plotted
        menu_types = []
        if 'navigation_buttons' in menu_current_keys:
            menu_types.append('navigation_buttons')
        if 'rangeboxes' in menu_current_keys:
            menu_types.append('rangeboxes')
        if 'checkboxes' in menu_current_keys:
            menu_types.append('checkboxes')
        #have at least1 menu type?
        if len(menu_types) > 0:
            #create grid of menu types (horizontally concatenating different menu type grids)
            grid = self.create_grid(menu_types)
            #add grid to parent layout
            parent_layout.addLayout(grid)

        #set finalised layout
        self.setLayout(parent_layout)

        #set geometry to match that of main window
        self.setGeometry(self.main_window_geometry)

        #show pop-up window
        self.show()

        #setup event to get selected checkbox indices when closing window
        quit_event = QtWidgets.QAction("Quit", self)
        quit_event.triggered.connect(self.closeEvent)

    def create_grid(self, menu_types):
        """create grid for each needed checkbox/rangebox/navigation button menu types, that wrap vertically
           and concatenate them horizontally together
        """

        #create horizontal layout to place all menu types within
        horizontal_parent = QtWidgets.QHBoxLayout()
        #set spacing between different menu type grids
        horizontal_parent.setSpacing(25)
        #align grids to centre and top
        horizontal_parent.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignTop)

        #order appearance of menu types grids in menu (from left to right)
        menu_type_order_dict = {'navigation_buttons':1, 'rangeboxes':2, 'checkboxes':3}
        menu_types = sorted(menu_types, key=menu_type_order_dict.__getitem__)

        #iterate through passed menu types
        for menu_type in menu_types:

            #create empty grid
            grid = QtWidgets.QGridLayout()
            #align grid to top
            grid.setAlignment(QtCore.Qt.AlignTop)

            #initialise indices used for indexing objects in grid
            start_row_n = 0
            row_n = 0
            column_n = 0

            #get dictionary nested inside current menu for menu type
            menu_current_type = self.menu_current[menu_type]

            #if have no labels for current menu type, continue to next menu type
            if len(menu_current_type['labels']) == 0:
                continue

            #get keys associated with current menu type in current menu level
            current_menu_keys = list(menu_current_type.keys())

            #create dictionary to store variables that save page status per menu type (also set formatting dict for each row in grid, and vertical spacing)
            if menu_type == 'checkboxes':
                row_format_dict = formatting_dict['checkbox_popup']
                grid_vertical_spacing = 0
                if ('keep_selected' in current_menu_keys) & ('remove_selected' in current_menu_keys):
                    self.page_memory['checkboxes'] = {'keep_selected':[], 'remove_selected':[], 'n_column_consumed':3, 'ordered_elements':['keep_selected','remove_selected'], 'widget':QtWidgets.QCheckBox}
                elif 'keep_selected' in current_menu_keys:
                    self.page_memory['checkboxes'] = {'keep_selected':[], 'n_column_consumed':2, 'ordered_elements':['keep_selected'], 'widget':QtWidgets.QCheckBox}
                elif 'remove_selected' in current_menu_keys:
                    self.page_memory['checkboxes'] = {'remove_selected':[], 'n_column_consumed':2, 'ordered_elements':['remove_selected'], 'widget':QtWidgets.QCheckBox}
            elif menu_type == 'rangeboxes':
                row_format_dict = formatting_dict['rangebox_popup']
                grid_vertical_spacing = 3
                if ('current_lower' in current_menu_keys) & ('current_upper' in current_menu_keys):
                    self.page_memory['rangeboxes'] = {'current_lower':[], 'current_upper':[], 'n_column_consumed':3, 'ordered_elements':['current_lower','current_upper'], 'widget':QtWidgets.QLineEdit}
                elif 'current_lower' in current_menu_keys:
                    self.page_memory['rangeboxes'] = {'current_lower':[], 'n_column_consumed':2, 'ordered_elements':['current_lower'], 'widget':QtWidgets.QLineEdit}
                elif 'current_upper' in current_menu_keys:
                    self.page_memory['rangeboxes'] = {'current_upper':[], 'n_column_consumed':2, 'ordered_elements':['current_upper'], 'widget':QtWidgets.QLineEdit}
            elif menu_type == 'navigation_buttons':
                row_format_dict = formatting_dict['navigation_button_popup']
                grid_vertical_spacing = 3
                self.page_memory['navigation_buttons'] = {'buttons':[], 'n_column_consumed':1, 'ordered_elements':['buttons'], 'widget':QtWidgets.QPushButton}

            #if have more than 1 column per label, need column headers
            if len(self.page_memory[menu_type]['ordered_elements']) > 1:
                have_column_headers = True
                start_row_n = 1
            else:
                have_column_headers = False

            #set vertical/horizontal grid spacing
            grid.setHorizontalSpacing(3)
            grid.setVerticalSpacing(grid_vertical_spacing)

            #calculate currently occupied vertical space
            occupied_vertical_space_before_grid = self.page_margin + formatting_dict['title_popup']['height'] + self.layout_spacing + self.page_margin
            if self.have_buttons:
                occupied_vertical_space_before_grid += (formatting_dict['button_popup']['height'] + self.layout_spacing)
            if have_column_headers:
                occupied_vertical_space_before_grid += (formatting_dict['column_header_label_popup']['height'] + grid_vertical_spacing)

            #initialise variable for tracking available vertical space when appending rows of grid
            currently_occupied_vertical_space = copy.deepcopy(occupied_vertical_space_before_grid)

            #iterate through all grid labels
            for label_ii, label in enumerate(menu_current_type['labels']):

                #evaluate if all available vertical space has been consumed
                row_available_space = self.main_window_geometry.height() - currently_occupied_vertical_space
                #if available space <= than row height, force a new column to be started
                if row_available_space <= (row_format_dict['height']):
                    column_n+=self.page_memory[menu_type]['n_column_consumed']
                    row_n = 0
                    currently_occupied_vertical_space = copy.deepcopy(occupied_vertical_space_before_grid)

                #if menu type == 'rangeboxes', then need to add label to left of rangeboxes, also add a tooltip (if exist)
                if (menu_type == 'checkboxes') or (menu_type == 'rangeboxes'):
                    rangebox_label = set_formatting(QtWidgets.QLabel(self, text = label), formatting_dict['rangebox_label_popup'])
                    if menu_type == 'rangeboxes':
                        if len(menu_current_type['tooltips']) > 0:
                            rangebox_label.setToolTip(wrap_tooltip_text(menu_current_type['tooltips'][label_ii], self.main_window_geometry.width()))
                    grid.addWidget(rangebox_label, start_row_n+row_n, column_n, QtCore.Qt.AlignLeft)

                #create all elements in column, per row
                for element_ii, element in enumerate(self.page_memory[menu_type]['ordered_elements']):
                    #if menu type == 'rangeboxes' then add 1 to element ii, because placed a label in first column
                    if (menu_type == 'checkboxes') or (menu_type == 'rangeboxes'):
                        element_label = ''
                        element_ii+=1
                    else:
                        element_label = label

                    #append widget to page memory dictionary
                    self.page_memory[menu_type][element].append(set_formatting(self.page_memory[menu_type]['widget'](element_label), row_format_dict))

                    #put text label to left of keep checkbox/rangebox (rather than to right)
                    if menu_type in ['checkboxes','rangeboxes']:
                        #check if checkbox is currently selected, if so select it again
                        #map checkbox value first if necessary
                        if menu_type == 'checkboxes':
                            if 'map_vars' in current_menu_keys:
                                var_to_check = menu_current_type['map_vars'][label_ii]
                            else:
                                var_to_check = copy.deepcopy(label)
                            if var_to_check in menu_current_type[element]:
                                self.page_memory[menu_type][element][label_ii].setCheckState(QtCore.Qt.Checked)
                        #set rangeboxes to previous set value (if any)
                        elif menu_type == 'rangeboxes':
                            self.page_memory[menu_type][element][label_ii].setText(menu_current_type[element][label_ii])

                    #if menu type == navigation_buttons, add connectivity to buttons, and also add tooltip
                    elif menu_type == 'navigation_buttons':
                        if len(menu_current_type['tooltips']) > 0:
                            self.page_memory[menu_type][element][label_ii].setToolTip(wrap_tooltip_text(menu_current_type['tooltips'][label_ii], self.main_window_geometry.width()))
                        self.page_memory[menu_type][element][label_ii].clicked.connect(self.open_new_page)

                    #add element to grid (aligned left)
                    grid.addWidget(self.page_memory[menu_type][element][label_ii], start_row_n+row_n, column_n+element_ii, QtCore.Qt.AlignLeft)

                #iterate row_n
                row_n +=1

                #add row vertical space to total occupied space
                currently_occupied_vertical_space += (row_format_dict['height'] + grid_vertical_spacing)

            #add column headers to menu type grid if needed
            if have_column_headers:
                for column_number in np.arange(0, column_n+1, self.page_memory[menu_type]['n_column_consumed']):
                    if menu_type == 'checkboxes':
                        grid.addWidget(set_formatting(QtWidgets.QLabel(self, text='K'), formatting_dict['column_header_label_popup']), 0, column_number+1, QtCore.Qt.AlignCenter)
                        grid.addWidget(set_formatting(QtWidgets.QLabel(self, text='R'), formatting_dict['column_header_label_popup']), 0, column_number+2, QtCore.Qt.AlignCenter)
                    elif menu_type == 'rangeboxes':
                        grid.addWidget(set_formatting(QtWidgets.QLabel(self, text='Lower Bound'), formatting_dict['column_header_label_popup']), 0, column_number+1, QtCore.Qt.AlignCenter)
                        grid.addWidget(set_formatting(QtWidgets.QLabel(self, text='Upper Bound'), formatting_dict['column_header_label_popup']), 0, column_number+2, QtCore.Qt.AlignCenter)

            #add menu type grid to horizontal layout
            horizontal_parent.addLayout(grid)

        #return horizontally concatenated menu type grids
        return horizontal_parent

    #------------------------------------------------------------------------#
    #------------------------------------------------------------------------#
    #functions that handle callbacks upon clicking on buttons

    def open_new_page(self):
        """function to open new page in pop-up window"""

        #get selected navigation button text
        selected_navigation_button = self.sender().text()
        #add selected navigation button text to menu levels list
        self.menu_levels.append(selected_navigation_button)
        #create new pop-up page for selected navigation button
        self.new_window = PopUpWindow(self.menu_root, self.menu_levels, self.main_window_geometry)
        #sleep briefly to allow new page to be generated
        time.sleep(0.1)
        #close current pop-up page
        self.close()

    def root_page(self):
        """function that returns pop-up window to root menu level page"""

        #create new pop-up page for root menu level
        self.new_window = PopUpWindow(self.menu_root, [], self.main_window_geometry)
        #sleep briefly to allow new page to be generated
        time.sleep(0.1)
        #close current pop-up page
        self.close()

    def previous_page(self):
        """function that returns pop-up window to previous menu level page"""

        #create new pop-up page for previous menu level
        self.new_window = PopUpWindow(self.menu_root, self.menu_levels[:-1], self.main_window_geometry)
        #sleep briefly to allow new page to be generated
        time.sleep(0.1)
        #close current pop-up page
        self.close()

    def select_all(self):
        """function to select all checkboxes"""
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for checkbox_ii, checkbox in enumerate(self.page_memory['checkboxes'][element]):
                self.page_memory['checkboxes'][element][checkbox_ii].setCheckState(QtCore.Qt.Checked)

    def clear_all(self):
        """function to clear all checkboxes"""
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for checkbox_ii, checkbox in enumerate(self.page_memory['checkboxes'][element]):
                self.page_memory['checkboxes'][element][checkbox_ii].setCheckState(QtCore.Qt.Unchecked)

    def select_all_default(self):
        """function to select all default selected checkboxes"""
        #unselect all checkboxes first
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for checkbox_ii, checkbox in enumerate(self.page_memory['checkboxes'][element]):
                self.page_memory['checkboxes'][element][checkbox_ii].setCheckState(QtCore.Qt.Unchecked)

        #map default variables to positional indices
        if 'map_vars' in list(self.menu_current['checkboxes'].keys()):
            default_inds = [np.where(self.menu_current['checkboxes']['map_vars'] == default_var)[0][0] for default_var in self.menu_current['checkboxes']['remove_default']]
        else:
            default_inds = [np.where(self.menu_current['checkboxes']['labels'] == default_var)[0][0] for default_var in self.menu_current['checkboxes']['remove_default']]

        #now select only desired default checkboxes
        for element in self.page_memory['checkboxes']['ordered_elements']:
            for default_ind in default_inds:
                self.page_memory['checkboxes'][element][default_ind].setCheckState(QtCore.Qt.Checked)

    #------------------------------------------------------------------------#
    #------------------------------------------------------------------------#

    def closeEvent(self, event):

        """function to get status of current page upon closing of pop-up window"""

        #take everything from page memory dictionary and put it back into menu level object

        #iterate through menu types
        for menu_type in list(self.page_memory.keys()):
            for element in self.page_memory[menu_type]['ordered_elements']:
                #if menu type == 'checkboxes', get variable names of all checkboxes ticked
                if menu_type == 'checkboxes':
                    selected_vars = []
                    for checkbox_ii, checkbox in enumerate(self.page_memory[menu_type][element]):
                        if checkbox.checkState() == QtCore.Qt.Checked:
                            #map selected position index to variable name
                            if 'map_vars' in list(self.menu_current[menu_type].keys()):
                                selected_vars.append(self.menu_current[menu_type]['map_vars'][checkbox_ii])
                            else:
                                selected_vars.append(self.menu_current[menu_type]['labels'][checkbox_ii])
                    #update previous selected variable
                    if element == 'keep_selected':
                        self.menu_current[menu_type]['previous_keep_selected'] = copy.deepcopy(self.menu_current[menu_type]['keep_selected'])
                    elif element == 'remove_selected':
                        self.menu_current[menu_type]['previous_keep_selected'] = copy.deepcopy(self.menu_current[menu_type]['remove_selected'])
                    #update selected variable
                    self.menu_current[menu_type][element] = selected_vars
                #if menu type == 'rangeboxes', get current values of all rangeboxes
                if menu_type == 'rangeboxes':
                    set_vals = []
                    for rangebox in self.page_memory[menu_type][element]:
                        set_vals.append(rangebox.text())
                    #update previous set variable
                    if element == 'current_lower':
                        self.menu_current[menu_type]['previous_lower'] = copy.deepcopy(self.menu_current[menu_type]['current_lower'])
                    elif element == 'current_upper':
                        self.menu_current[menu_type]['previous_upper'] = copy.deepcopy(self.menu_current[menu_type]['current_upper'])
                    #update set value
                    self.menu_current[menu_type][element] = set_vals

