""" Settings menus """

import json
import os
import platform
import yaml

from functools import partial
from PyQt5 import QtCore, QtGui, QtWidgets 

from .dashboard_elements import CheckableComboBox, ComboBox
from .dashboard_elements import set_formatting

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
settings_dict = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/canvas_menus.yaml')))
# get operating system specific formatting
operating_system = platform.system()
if operating_system == 'Darwin':
    formatting_dict = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/internal/stylesheet_mac.yaml')))
elif operating_system == 'Linux':
    formatting_dict = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/internal/stylesheet_linux.yaml')))
elif operating_system in ['Windows','MINGW32_NT','MINGW64_NT']:
    formatting_dict = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/internal/stylesheet_windows.yaml')))

class SettingsMenu(object):

    def __init__(self, plot_type, canvas_instance):
        
        self.canvas_instance = canvas_instance

        self.elements = list(settings_dict[plot_type].keys())
        self.buttons = {}
        self.labels = {}
        self.comboboxes = {}
        self.checkable_comboboxes = {}
        self.sliders = {}

        for element_name in self.elements:
            element_settings = settings_dict[plot_type][element_name]
            element_type = element_settings['element_type']
            if element_type in ['button', 'container', 'label', 'combobox', 'checkable_combobox', 'slider']:

                # Add element
                element = getattr(self, 'add_' + element_type)(element_settings)

                # Add options as items to options combobox
                if element_name == 'options':
                    if plot_type in ['periodic_violin','fairmode_target','fairmode_statsummary']:
                        plot_type_corr = plot_type.replace('_','-')
                    else:
                        plot_type_corr = plot_type
                    element.addItems(self.canvas_instance.plot_characteristics[plot_type_corr]['plot_options']) 
                
                # Apply common properties
                if 'relative_position' in element_settings.keys():
                    element.move(self.buttons['settings_button'].x()+element_settings['relative_position'][0],
                                 self.buttons['settings_button'].y()+element_settings['relative_position'][1]) 
                if 'size' in element_settings.keys():
                    element.resize(element_settings['size'][0], element_settings['size'][1])
                if 'fixed_width' in element_settings.keys():
                    element.setFixedWidth(element_settings['fixed_width'])
                if 'style' in element_settings.keys():
                    element.setStyleSheet(element_settings['style'])
                if 'object_name' in element_settings.keys():
                    element.setObjectName(element_settings['object_name'])

                # Hide element
                element.hide()

                # Save element in corresponding dictionary
                if element_type == 'container':
                    self.container = element
                elif element_type in ['button', 'label', 'slider']:
                    getattr(self, element_type + 's')[element_name] = element
                elif element_type in ['combobox', 'checkable_combobox']:
                    getattr(self, element_type + 'es')[element_name] = element

            else:
                print(f'ERROR: Unknown element type: {element_type}')


    def add_button(self, element_settings):
        """ Add button as settings menu

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QPushButton
            Button
        """

        button = set_formatting(QtWidgets.QPushButton(self.canvas_instance), 
                                                      formatting_dict[element_settings['formatting_dict']])
        button.setIcon(QtGui.QIcon(os.path.join(CURRENT_PATH, element_settings['path'])))
        button.setIconSize(QtCore.QSize(element_settings['size'][0], element_settings['size'][1]))
        button.clicked.connect(partial(self.connect, element_settings['function']))
        
        return button 

    def add_container(self, element_settings):
        """ Add elements container

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QWidget
            Container
        """

        container = set_formatting(QtWidgets.QWidget(self.canvas_instance), 
                                                     formatting_dict[element_settings['formatting_dict']])
        container.raise_()

        return container

    def add_label(self, element_settings):
        """ Add label

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QLabel
            Label
        """
                
        if 'formatting_dict' in element_settings.keys():
            label = set_formatting(QtWidgets.QLabel(element_settings['text'], self.canvas_instance), 
                                   formatting_dict[element_settings['formatting_dict']])
        else:
            label = QtWidgets.QLabel(element_settings['text'], self.canvas_instance)

        return label

    def add_combobox(self, element_settings):
        """ Add combobox

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QComboBox
            Combobox
        """

        combobox = set_formatting(ComboBox(self.canvas_instance), formatting_dict[element_settings['formatting_dict']])
        combobox.currentTextChanged.connect(partial(self.connect, element_settings['function']))

        return combobox

    def add_slider(self, element_settings):
        """ Add slider

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
        slider.setTracking(element_settings['tracking'])
        slider.setTickInterval(element_settings['tick_interval'])
        if 'minimum' in element_settings.keys():
            slider.setMinimum(int(element_settings['minimum']))
        if 'maximum' in element_settings.keys():
            slider.setMaximum(int(element_settings['maximum']))
        if 'value' in element_settings.keys():
            slider.setValue(int(element_settings['value']))
        slider.valueChanged.connect(partial(self.connect, element_settings['function']))

        return slider

    def add_checkable_combobox(self, element_settings):
        """ Add checkable combobox

        Parameters
        ----------
        element_settings : dict
            Settings

        Returns
        -------
        QtWidgets.QComboBox
            Combobox with options to check
        """

        checkable_combobox = set_formatting(CheckableComboBox(self.canvas_instance), 
                                            formatting_dict[element_settings['formatting_dict']])
        checkable_combobox.currentTextChanged.connect(partial(self.connect, element_settings['function']))

        return checkable_combobox

    def get_elements(self):
        """ Get elements inside menu settings
        
        Returns
        -------
        list
            Menu settings elements
        """

        sliders = list(self.sliders.values())
        comboboxes = list(self.comboboxes.values())
        labels = list(self.labels.values())
        checkable_comboboxes = list(self.checkable_comboboxes.values())

        return [self.container] + sliders + comboboxes + labels + checkable_comboboxes

    def connect(self, function):
        """ Connect element to functions in settings dictionary """

        if hasattr(self.canvas_instance, 'interactive_elements'):
            # Call function only after all elements have been added
            if self.canvas_instance.interactive_elements.keys() == settings_dict.keys():
                getattr(self.canvas_instance, function)()

        return None