""" Auxiliary interpolation functions """

import copy
import os
import numpy as np


def get_aeronet_bin_radius_from_bin_variable(binvar):
    """
    Returns the specific AERONET particle radius associated with a given bin variable name.
    AERONET bins are not bins in the classical sense in that they have a radius min/max,
    rather they represent an instance on the distribution curve across all sizes.

    Parameters
    ----------
    binvar : str
        The AERONET bin variable identifier (e.g. 'vconcaerobin1').

    Returns
    -------
    float
        The radius value in micrometres corresponding to the size distribution instance.
    """

    aeronet_bin_variable_to_bin_radius = {
        "vconcaerobin1": 0.05,
        "vconcaerobin2": 0.065604,
        "vconcaerobin3": 0.086077,
        "vconcaerobin4": 0.112939,
        "vconcaerobin5": 0.148184,
        "vconcaerobin6": 0.194429,
        "vconcaerobin7": 0.255105,
        "vconcaerobin8": 0.334716,
        "vconcaerobin9": 0.439173,
        "vconcaerobin10": 0.576227,
        "vconcaerobin11": 0.756052,
        "vconcaerobin12": 0.991996,
        "vconcaerobin13": 1.301571,
        "vconcaerobin14": 1.707757,
        "vconcaerobin15": 2.240702,
        "vconcaerobin16": 2.939966,
        "vconcaerobin17": 3.857452,
        "vconcaerobin18": 5.061260,
        "vconcaerobin19": 6.640745,
        "vconcaerobin20": 8.71314,
        "vconcaerobin21": 11.432287,
        "vconcaerobin22": 15.00,
    }

    return aeronet_bin_variable_to_bin_radius[binvar]


def get_model_bin_radii(model_name):
    """
    Returns the particle bin edge radii and bin particle densities for a specified model.

    Parameters
    ----------
    model_name : str
        The name of the model to retrieve bin information for.

    Returns
    -------
    r_edges : list or None
        A list of bin edge radii if the model is recognised; otherwise None.
    rho_bins : list or None
        A list of particle densities (kg/m^3) if the model is recognised; otherwise None.
    """

    if "monarch" in model_name:
        r_edges = [0.1, 0.18, 0.3, 0.6, 1.0, 1.8, 3.0, 6.0, 10.0]
        rho_bins = [2500.0, 2500.0, 2500.0, 2500.0, 2650.0, 2650.0, 2650.0, 2650.0]
        return r_edges, rho_bins
    else:
        return None, None


def get_aeronet_model_bin(model_name, aeronet_bin_radius):
    """
    Identifies the specific model bin that encompasses a given AERONET bin radius.

    Parameters
    ----------
    model_name : str
        The name of the model to retrieve bin information for.
    aeronet_bin_radius : float
        The radius value from the AERONET distribution.

    Returns
    -------
    bin_index : int
        The index of the matching model bin.
    bin_min : float
        The minimum radius edge of the identified model bin.
    bin_max : float
        The maximum radius edge of the identified model bin.
    rho_bin : float
        The particle density (kg/m^3) associated with the identified model bin.
    """

    # get model bin raddi and bin rh
    r_edges, rho_bins = get_model_bin_radii(model_name)

    # get model bin index which contains AERONET bin radius instance
    if aeronet_bin_radius == r_edges[-1]:
        bin_index = len(r_edges) - 1
    else:
        bin_index = np.searchsorted(r_edges, aeronet_bin_radius, side="right") - 1

    return bin_index, r_edges[bin_index], r_edges[bin_index + 1], rho_bins[bin_index]


def get_model_to_aeronet_bin_transform_factor(model_name, rmin, rmax):
    """
    Calculates the factor to transform aerosol size distributions from model bins to AERONET's
    22 bin format, assuming a constant function.

    Parameters
    ----------
    model_name : str
        The name of the model being processed.
    rmin : float
        The minimum radius of the specific model bin.
    rmax : float
        The maximum radius of the specific model bin.

    Returns
    -------
    bin_transform_factor : float
    """

    # get bin integral (per model)
    if model_name == "monarch":
        bin_transform_factor = 1.0 / (np.log(rmax) - np.log(rmin))

    return bin_transform_factor


def check_for_ghost(network_name):
    """
    Determines whether a selected network is sourced from GHOST or not.

    Parameters
    ----------
    network_name : str
        The name of the network to be checked.

    Returns
    -------
    bool
        True if the network is a GHOST network, False otherwise.
    """

    if "/" in network_name:
        return False
    else:
        return True


def check_directory_existence(directory_tree_str, directories_not_to_test=None):
    """
    Validates, creates, and sets standardised permissions and ownership for a directory hierarchy.

    Parameters
    ----------
    directory_tree_str : str
        The full path of the directory tree to be processed.
    directories_not_to_test : str, optional
        A prefix of the directory tree that should be excluded from existence and permission checks.
    """

    # define special characters that will need to be escaped for execution of commands
    special_characters = ["(", ")"]

    # first check if instructed to not to check existence of part of directory tree str
    if directories_not_to_test is not None:
        # modify directory_tree_str to not include directories_not_to_test
        directory_tree_str = directory_tree_str.replace(directories_not_to_test, "")
        # set directory_str_to_test to be directories_not_to_test
        directory_str_to_test = copy.deepcopy(directories_not_to_test)
    else:
        # else, set directory_str_to_test as initially empty string
        directory_str_to_test = ""

    # split directory tree str
    directory_tree_str_split = directory_tree_str.split("/")

    # iterate through directory_tree_str_split
    for current_directory in directory_tree_str_split:
        # if current directory is empty string then do not check directory existence
        if current_directory != "":
            # add current_dictionary to directory_str_to_test
            directory_str_to_test = directory_str_to_test + "/" + current_directory

            # does directory_str_to_test exist?
            if not os.path.isdir(directory_str_to_test):
                # escape certain special characters in directory_str_to_test
                alt_directory_str_to_test = copy.deepcopy(directory_str_to_test)
                for ch in special_characters:
                    alt_directory_str_to_test = alt_directory_str_to_test.replace(
                        ch, "\{}".format(ch)
                    )

                # if not, create it
                os.system("mkdir {}".format(alt_directory_str_to_test))

                # give 770 permissions to directory
                os.system("chmod 770 {}".format(alt_directory_str_to_test))

                # make group owner bsc32
                os.system("chgrp bsc32 {}".format(alt_directory_str_to_test))


def set_file_permissions_ownership(file_str):
    """
    Sets standardised file permissions and group ownership for a specified file path.

    Parameters
    ----------
    file_str : str
        The path of the file to be updated with new permissions and ownership.
    """

    # define special characters that will need to be escaped for execution of commands
    special_characters = ["(", ")"]

    # escape certain special characters in file_str
    for ch in special_characters:
        file_str = file_str.replace(ch, "\{}".format(ch))

    # give 770 permissions to file
    os.system("chmod 770 {}".format(file_str))

    # make group owner bsc32
    os.system("chgrp bsc32 {}".format(file_str))


def findMiddle(input_len):
    """
    Determines the middle index or indices for a given length.

    Parameters
    ----------
    input_len : int
        The total length of the sequence to calculate the middle for.

    Returns
    -------
    middle : int or list
        The middle index as an integer for odd lengths, or a list of two middle indices for even lengths.
    """

    middle = float(input_len) / 2
    if middle % 2 != 0:
        return int(middle - 0.5)
    else:
        return [int(middle - 1), int(middle)]
