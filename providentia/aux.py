import json
import numpy as np


def which_bounds(instance, species):
    """If there are bounds defined in a config file, fill that value,
    if it is withing the feasible bounds of the species"""

    lower = np.float32(instance.parameter_dictionary[species]['extreme_lower_limit'])
    upper = np.float32(instance.parameter_dictionary[species]['extreme_upper_limit'])

    if hasattr(instance, 'lower_bound'):
        if instance.lower_bound >= lower:
            lower = instance.lower_bound

    if hasattr(instance, 'upper_bound'):
        if instance.upper_bound <= upper:
            upper = instance.upper_bound

    return np.float32(lower), np.float32(upper)


def which_qa(instance, return_defaults=False):
    """Checks if the species we currently have selected belongs to the ones
    that have specific qa flags selected as default"""

    if return_defaults or (not hasattr(instance, 'qa')):
        if instance.selected_species in instance.qa_exceptions:
            return instance.specific_qa
        else:
            return instance.general_qa

    if hasattr(instance, 'qa'):
        # if conf has only 1 QA
        if isinstance(instance.qa, int):
            return [instance.qa]
        # if the QAs are written with their names
        elif instance.qa == "":
            return []
        elif isinstance(instance.qa, str):
            return [instance.standard_QA_name_to_QA_code[q.strip()] for q in instance.qa.split(",")]
        # return subset the user has selected in conf
        else:
            return list(instance.qa)


def which_flags(instance):
    """if there are flags coming from a config file, select those"""

    if hasattr(instance, 'flags'):
        # if conf has only one flag
        if isinstance(instance.flags, int):
            return [instance.flags]
        # if flags are writtern as strings
        elif instance.flags == "":
            return []
        elif isinstance(instance.flags, str):
            return [instance.standard_data_flag_name_to_data_flag_code[f.strip()]
                    for f in instance.flags.split(",")]
        else:
            return list(instance.flags)
    else:
        return []


def get_qa_codes(instance):
    """Retrieve QA codes from GHOST_standards using the qa flags' names.

    Specific flags are defined for the following species:
    ['WND_DIR_10M','WND_SPD_10M','RH_2M','PREC_ACCUM','SNOW_ACCUM',
    'SNOW_DEPTH','CEILING_HEIGHT','VIS_DIST','CLOUD_CVG','CLOUD_CVG_FRAC']"""

    # get names from json files
    specific_qa_names = json.load(open(
        "providentia/conf/default_flags.json"))['specific_qa']
    general_qa_names = json.load(open(
        "providentia/conf/default_flags.json"))['general_qa']
    # get codes
    specific_qa = [instance.standard_QA_name_to_QA_code[qa_name]
                   for qa_name in specific_qa_names]
    general_qa = [instance.standard_QA_name_to_QA_code[qa_name]
                  for qa_name in general_qa_names]
    # get difference of flags, needed later for updating default selection
    qa_diff = list(set(general_qa) - set(specific_qa))

    return specific_qa, general_qa, qa_diff
