import os

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
        Second dictionary_
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

    keys_to_remove = ["dashboard", "offline", "interactive", "tests"]
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
                    if mode != 'offline':
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
