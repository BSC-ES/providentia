""" Auxiliar functions. """


def get_aeronet_bin_radius_from_bin_variable(binvar):
    """ Return AERONET bin radius instance from bin variable AERONET bins are not bins in the classical sense in that 
        they have a radius min/max, rather they represent an instance on the distribution curve across all sizes.      
        
        :param binvar: bin variable
        :type binvar: str
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


def check_for_ghost(network_name):
    """ Check whether the selected network comes from GHOST or not.
        All non-GHOST networks start with an asterisk at their name.

        :param network_name: network name
        :type network_name: str
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
                
                # make group owner Earth
                os.system("chgrp Earth {}".format(alt_directory_str_to_test))

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

    # make group owner Earth
    os.system("chgrp Earth {}".format(file_str))

def findMiddle(input_len):
    """ Find middle index/indices of a list.

        :param input_len: list of coordinates
        :type input_len: list
    """

    middle = float(input_len)/2
    if middle % 2 != 0:
        return int(middle - .5)
    else:
        return [int(middle-1), int(middle)]
        