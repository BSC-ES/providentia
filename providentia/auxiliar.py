import os

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

def join(*args):
    """Join paths making sure they are all in the right direction
    """
    return os.path.join(*args).replace("\\","/")


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
    for section in plot_characteristics:
        # get all keys
        section_data = plot_characteristics[section]
        if mode in section_data:
            for key, value in section_data[mode].items():
                # if it already exists (under general settings), update values
                if key in section_data.keys():
                    section_data[key].update(value)
                # if it doesn't exist, create new key
                else:
                    section_data[key] = value

        # remove mode keys
        for key in keys_to_remove:
            section_data.pop(key, None)

    return plot_characteristics