import os
import numpy as np
import socket
import time

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


def join(*args):
    """Join paths making sure they are all in the right direction
    """
    return os.path.join(*args).replace("\\", "/")


def deep_merge(dict1, dict2):
    """Merge dictionaries recursively to avoid values getting replaced

    Parameters
    ----------
    dict1 : dict
        First dictionary
    dict2 : dict
        Second dictionary
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
    """Get values from active mode and expand generic plot characteristics

    Parameters
    ----------
    plot_characteristics : dict
        Plot characteristics dictionary
    mode : str
        Active mode
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
    """Pad array with pad value if its length is less than input length

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

    # get BSC machine name (if have one)
    machine = os.environ.get('BSC_MACHINE', None)

    # set current machine
    if machine is None:
        hostname = os.environ.get('HOSTNAME', '')
        
        #setup retrial system for getting ip address as occasionaly breaks
        retry = 0
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
        elif "bscesoper01.bsc.es" in hostname:
            machine = "oper"
        elif ip == "84.88.185.48":
            machine = "hub"
        else:
            machine = "local"

    return machine