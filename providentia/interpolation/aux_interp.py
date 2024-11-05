""" Auxiliary functions. """

def get_aeronet_bin_radius_from_bin_variable(binvar):
    """ Return AERONET bin radius instance from bin variable. AERONET bins are not bins in the classical sense in that 
        they have a radius min/max, rather they represent an instance on the distribution curve across all sizes.      
        
        :param binvar: bin variable
        :type binvar: str
        :return: AERONET bin radius
        :rtype: float
    """
    
    aeronet_bin_variable_to_bin_radius = {'vconcaerobin1': 0.05,       'vconcaerobin2': 0.065604, 
                                          'vconcaerobin3': 0.086077,   'vconcaerobin4': 0.112939, 
                                          'vconcaerobin5': 0.148184,   'vconcaerobin6': 0.194429, 
                                          'vconcaerobin7': 0.255105,   'vconcaerobin8': 0.334716, 
                                          'vconcaerobin9': 0.439173,   'vconcaerobin10': 0.576227,
                                          'vconcaerobin11': 0.756052,  'vconcaerobin12': 0.991996, 
                                          'vconcaerobin13': 1.301571,  'vconcaerobin14': 1.707757, 
                                          'vconcaerobin15': 2.240702,  'vconcaerobin16': 2.939966, 
                                          'vconcaerobin17': 3.857452,  'vconcaerobin18': 5.061260, 
                                          'vconcaerobin19': 6.640745,  'vconcaerobin20': 8.71314,
                                          'vconcaerobin21': 11.432287, 'vconcaerobin22': 15.00}

    return aeronet_bin_variable_to_bin_radius[binvar]

def get_model_bin_radii(model_name):
    """ Return model bin edges radii, and bin relative humidity
        
        :param model_name: name of model
        :type model_name: str
        :return: model bin radii, model bin rh 
        :rtype: list, list
    """ 

    if 'monarch' in model_name:
        r_edges =[0.1, 0.18, 0.3, 0.6, 1.0, 1.8, 3.0, 6.0, 10.0]
        rho_bins = [2500.0, 2500.0, 2500.0, 2500.0, 2650.0, 2650.0, 2650.0, 2650.0]
        return r_edges, rho_bins
    else:
        return None, None

def get_aeronet_model_bin(model_name, aeronet_bin_radius):
    """ Return details of model bin which contains AERONET bin radius instance, and 
        
        :param model_name: name of model
        :type model_name: str
        :param aeronet_bin_radius: bin radius
        :type aeronet_bin_radius: float
        :return: model bin index, bin radius minimum, bin radius maximum, bin relative humidity
        :rtype: int, float, float, float
    """ 

    import numpy as np

    # get model bin raddi and bin rh
    r_edges, rho_bins = get_model_bin_radii(model_name)

    # get model bin index which contains AERONET bin radius instance
    if aeronet_bin_radius == r_edges[-1]:
        bin_index = len(r_edges)-1
    else:
        bin_index = np.searchsorted(r_edges, aeronet_bin_radius, side='right') - 1

    return bin_index, r_edges[bin_index], r_edges[bin_index+1], rho_bins[bin_index]

def get_model_to_aeronet_bin_transform_factor(model_name, rmin, rmax):
    """ Return factor which transforms aerosol size distribution data from model bins 
        to AERONET's 22 bins format, assuming a constant function.

        :param model_name: name of model
        :type model_name: str
        :param rmin: minimum bin radius
        :type rmin: float
        :param rmax: maximum bin radius
        :type rmax: float
        :return: transform factor
        :rtype: np.float32
    """

    import numpy as np

    # get bin integral (per model)
    if model_name == 'monarch':
        bin_transform_factor = 1.0/(np.log(rmax) - np.log(rmin))

    return bin_transform_factor

def check_for_ghost(network_name):
    """ Check whether the selected network comes from GHOST or not.
        All non-GHOST networks start with an asterisk at their name.

        :param network_name: network name
        :type network_name: str
        :return: True or False if network is from GHOST 
        :rtype: boolean
    """

    if '/' in network_name:
        return False
    else:
        return True
        
def check_directory_existence(directory_tree_str, directories_not_to_test=None):
    """ Iterate through provided directory tree string.
        First, check if directory exists, if not create it.
        Second, check if each directory has 770 permissions.
        Finally, also ensure group owner is bsc32.

        :param directory_tree_str: directory tree
        :type directory_tree_str: str
        :param directories_not_to_test: directories that will not be checked
        :type directories_not_to_test: str
    """

    import copy
    import os

    # define special characters that will need to be escaped for execution of commands      
    special_characters = ['(',')']

    # first check if instructed to not to check existence of part of directory tree str 
    if directories_not_to_test is not None:
        # modify directory_tree_str to not include directories_not_to_test
        directory_tree_str = directory_tree_str.replace(directories_not_to_test,'')
        # set directory_str_to_test to be directories_not_to_test
        directory_str_to_test = copy.deepcopy(directories_not_to_test)
    else:
        # else, set directory_str_to_test as initially empty string
        directory_str_to_test = ''

    # split directory tree str
    directory_tree_str_split = directory_tree_str.split('/')

    # iterate through directory_tree_str_split
    for current_directory in directory_tree_str_split:

        # if current directory is empty string then do not check directory existence
        if current_directory != '':

            # add current_dictionary to directory_str_to_test
            directory_str_to_test = directory_str_to_test + '/' + current_directory

            # does directory_str_to_test exist? 
            if os.path.isdir(directory_str_to_test) == False:
                
                # escape certain special characters in directory_str_to_test
                alt_directory_str_to_test = copy.deepcopy(directory_str_to_test)
                for ch in special_characters:
                    alt_directory_str_to_test = alt_directory_str_to_test.replace(ch,'\{}'.format(ch))

                # if not, create it
                os.system("mkdir {}".format(alt_directory_str_to_test))

                
                
                # give 770 permissions to directory
                os.system("chmod 770 {}".format(alt_directory_str_to_test))
                
                # make group owner bsc32
                os.system("chgrp bsc32 {}".format(alt_directory_str_to_test))
            os.system("mkdir /gpfs/scratch/bsc32/bsc032388/Providentia/logs/interpolation/arguments/pip")
            os.system("touch /gpfs/scratch/bsc32/bsc032388/Providentia/logs/interpolation/arguments/pip")

def set_file_permissions_ownership(file_str):
    """ Set 770 permissions for a newly written file and also ensure group owner is bsc32.

        :param file_str: name of the file to set permissions and owner
        :type file_str: str
    """

    import os
    
    # define special characters that will need to be escaped for execution of commands      
    special_characters = ['(',')']

    # escape certain special characters in file_str
    for ch in special_characters:
        file_str = file_str.replace(ch,'\{}'.format(ch))

    # give 770 permissions to file
    os.system("chmod 770 {}".format(file_str))

    # make group owner bsc32
    os.system("chgrp bsc32 {}".format(file_str))

def findMiddle(input_len):
    """ Find middle index/indices of a list.

        :param input_len: list of coordinates
        :type input_len: list
        :return: middle index of list 
        :rtype: int
    """

    middle = float(input_len)/2
    if middle % 2 != 0:
        return int(middle - .5)
    else:
        return [int(middle-1), int(middle)]