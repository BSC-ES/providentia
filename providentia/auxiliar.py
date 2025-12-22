""" Auxiliary functions """

import os
import socket
import sys
import time

import numpy as np

from .unit_converter import UnitConverter, get_molecular_mass

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


def join(*args):
    """
    Join paths making sure they are all in the right direction
    """
    return os.path.join(*args).replace("\\", "/")


def deep_merge(dict1, dict2):
    """
    Merge dictionaries recursively to avoid values getting replaced

    Parameters
    ----------
    dict1 : dict
        First dictionary
    dict2 : dict
        Second dictionary

    Returns
    -------
        Merged dictionary
    """

    for key, value in dict2.items():
        if key in dict1 and isinstance(dict1[key], dict) and isinstance(value, dict):
            # if both values are dictionaries merge them recursively
            deep_merge(dict1[key], value)
        else:
            # replace or add the value
            dict1[key] = value

    return dict1


def expand_plot_characteristics(plot_characteristics, mode):
    """
    Get values from active mode and expand generic plot characteristics

    Parameters
    ----------
    plot_characteristics : dict
        Plot characteristics dictionary
    mode : str
        Active mode

    Returns
    -------
    dict
        Plot characteristics dictionary
    """

    keys_to_remove = ["dashboard", "report", "library", "tests"]
    for plot_type in plot_characteristics:
        # get all keys
        plot_type_characteristics = plot_characteristics[plot_type]
        if mode in plot_type_characteristics:
            for key, value in plot_type_characteristics[mode].items():
                # if it already exists (under general settings), update values
                if key in plot_type_characteristics.keys():
                    # if it is a dict, merge general dictionary with and mode dictionary
                    if isinstance(plot_type_characteristics[key], dict):
                        plot_type_characteristics[key] = deep_merge(
                            plot_type_characteristics[key].copy(), value)
                    # if it is a list, extend list removing duplicates
                    elif isinstance(plot_type_characteristics[key], list):
                        plot_type_characteristics[key] = list(
                            set(plot_type_characteristics[key] + value))
                # if it doesn't exist, create new key
                else:
                    plot_type_characteristics[key] = value
                # remove options from plot options if they do not apply to active mode
                if key == 'plot_options':
                    if mode != 'report':
                        if 'obs' in value:
                            value.remove('obs')
                        if 'individual' in value:
                            value.remove('individual')
                        # remove bias option in dashboard only for map plot type
                        # since we select the statistics from other dropdowns
                        if ('bias' in value) and (plot_type[:4] == 'map'):
                            value.remove('bias')
                    if (mode == 'dashboard') and ('multispecies' in value):
                        value.remove('multispecies')
                    plot_type_characteristics[key] = value

        # remove mode keys
        for key in keys_to_remove:
            plot_type_characteristics.pop(key, None)

    return plot_characteristics


def pad_array(arr, length, pad_value=np.nan):
    """
    Pad array with pad value if its length is less than input length

    Parameters
    ----------
    arr : numpy.array
        Array to pad
    length : int
        Length the array will have after padding
    pad_value : int, float, optional
        Value to append if length is lesss than input length, by default np.nan

    Returns
    ------
    numpy.array
        Padded array
    """

    pad_size = max(0, length - len(arr))

    return np.pad(arr, (0, pad_size), constant_values=pad_value)


def get_machine():
    """
    Get machine where code is running

    Returns
    -------
    str
        Machine
    """

    # return github machine if tests are running in actions
    if os.getenv("GITHUB_ACTIONS") == "true":
        return "github"
    
    # get BSC machine name (if have one)
    machine = os.environ.get('BSC_MACHINE', None)

    # set current machine
    if machine is None:
        hostname = os.environ.get('HOSTNAME', '')
        
        # setup retrial system for getting ip address as occasionaly breaks
        retry = 0
        ip = None
        while True:
            try:
                ip = socket.gethostbyname(socket.gethostname())
                break
            except:
                if retry == 3:
                    break
                else:
                    retry+=1
                    time.sleep(1)
        if "bscearth" in hostname:
            machine = "workstation"
        elif "transfer" in hostname:
            machine = "storage5"
        elif "bscesdust02.bsc.es" in hostname:
            machine = "dust"
        elif ip == "84.88.185.205":
            machine = "oper"
        elif ip == "84.88.185.48":
            machine = "hub"
        else:
            machine = "local"

    return machine

class Tee:
    """
    Class to show interpolation output in terminal from notebook
    """

    def __init__(self, *streams):
        """
        Initialise the object

        Parameters
        ----------
        *streams : file-like objects
            One or more streams (e.g., sys.stdout, open file objects) where
            output should be duplicated
        """

        self.streams = streams

    def write(self, data):
        """
        Write data to all output streams

        Parameters
        ----------
        data : str
            The string to write to each stream
        """

        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self):
        """
        Flush all output streams ensuring that any buffered output is written immediately
        """

        for s in self.streams:
            s.flush()


def get_standard_parameters_by_speci(speci, ghost_version):
    """
    Get GHOST standard parameters dictionary

    Parameters
    ----------
    speci : str
        Speci to plot
    ghost_version : str
        GHOST version

    Returns
    -------
    dict
        GHOST standard parameters dictionary
    """

    sys.path.insert(1, join(CURRENT_PATH, 'dependencies/GHOST_standards/{}'.format(ghost_version)))
    from GHOST_standards import standard_parameters
    standard_parameters = standard_parameters
    
    # get cut of standard parameters for original speci to process
    for standard_parameter in standard_parameters.keys():
        if standard_parameters[standard_parameter]['bsc_parameter_name'] == speci:
            return standard_parameters[standard_parameter]


def unit_conversion(initial_units, final_units, standard_parameter_speci):
    """
    Use unit_converter class to return conversion object from initial to final units.

    Parameters
    ----------
    initial_units : str
        Input units
    final_units : str
        Output units
    standard_parameter_speci : dict
        GHOST standard parameters dictionary

    Returns
    -------
    obj, str
        Conversion object or error
    """

    # if units are unitless, then no need for conversion (i.e. conversion factor = 1.0)   
    if (final_units == 'unitless') or (final_units == '-') or (final_units == '1'):
        return type('UnitConverter', (object,), {'conversion_factor':1.0, 'output_standard_units':'unitless'})
    
    # determine chemical formula of species 
    if 'chemical_formula_charge' in list(standard_parameter_speci.keys()):
        speci_chemical_formula = standard_parameter_speci['chemical_formula_charge']
    else:
        speci_chemical_formula = standard_parameter_speci['chemical_formula']

    # get input (model) quantity for conversion
    if ('units_quantity' in list(standard_parameter_speci.keys())) and (initial_units == standard_parameter_speci['standard_units']):
        initial_quantity = standard_parameter_speci['units_quantity']
    else:
        conv_obj = UnitConverter(initial_units, initial_units, 1, species=speci_chemical_formula)
        initial_quantity = conv_obj.output_quantity

    # get output (observational) quantity for conversion
    if ('units_quantity' in list(standard_parameter_speci.keys())) and (final_units == standard_parameter_speci['standard_units']):
        final_quantity = standard_parameter_speci['units_quantity']
    else:
        conv_obj = UnitConverter(final_units, final_units, 1, species=speci_chemical_formula)
        final_quantity = conv_obj.output_quantity

    # unit converter module does not produce conversion factor for temperature, 
    # but both input and output units should be Kelvin (i.e. conversion factor = 1.0) 
    # if input not in K then return error
    if final_quantity == 'temperature': 
        if initial_units == 'K':
            return type('UnitConverter', (object,), {'conversion_factor':1.0, 'output_standard_units':'K'})
        else:
            error = "Error: Experiment units should be 'K', but are set as '{}'".format(initial_units)
            return error

    # initial and final quantities not equal (convert to observational units, 
    # standard_temperature=293.15, standard_pressure=1013.25)
    if initial_quantity != final_quantity:     
        
        # convert units
        input_units = {'temperature': 'K', 
                       'pressure': 'hPa', 
                       'molar_mass': 'kg mol-1', 
                       initial_quantity: initial_units}
        input_values = {'temperature': 293.15, 
                        'pressure': 1013.25, 
                        'molar_mass': get_molecular_mass(speci_chemical_formula), 
                        initial_quantity: 1.0}
        conv_obj = UnitConverter(input_units, final_units, input_values, 
                                 species=speci_chemical_formula, 
                                 input_quantity=initial_quantity, 
                                 output_quantity=final_quantity)
    
    # same quantity conversion
    else:
        conv_obj = UnitConverter(initial_units, final_units, 1.0, 
                                 species=speci_chemical_formula, 
                                 input_quantity=initial_quantity, 
                                 output_quantity=final_quantity) 
    
    # return conversion object
    return conv_obj


def get_conversion_factor(initial_units, final_units, standard_parameter_speci):
    """
    Get conversion factor to convert from initial to final units.
    Convenience wrapper for unit_conversion function.

    Parameters
    ----------
    initial_units : str
        Input units
    final_units : str
        Output units
    standard_parameter_speci : dict
        GHOST standard parameters dictionary

    Returns
    -------
    float, str
        Conversion factor or error
    """

    conv_obj = unit_conversion(initial_units, final_units, standard_parameter_speci)
    return conv_obj.conversion_factor


def get_standard_units(initial_units, standard_parameter_speci):
    """
    Get standardised units for given initial units.
    Convenience wrapper for unit_conversion function.

    Parameters
    ----------
    initial_units : str
        Input units
    standard_parameter_speci : dict
        GHOST standard parameters dictionary

    Returns
    -------
    str
        Output standard units or error
    """
    
    conv_obj = unit_conversion(initial_units, initial_units, standard_parameter_speci)
    return conv_obj.output_standard_units