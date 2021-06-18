import sys
import pandas as pd
from providentia import aux


class InitStandards:

    def __init__(self, obs_root, ghost_version):
        """ Read from ghost standards """
        sys.path.insert(1, '{}/GHOST_standards/{}'.format(obs_root, ghost_version))
        from GHOST_standards import standard_parameters, \
            get_standard_metadata, standard_data_flag_name_to_data_flag_code, \
            standard_QA_name_to_QA_code
        # modify standard parameter dictionary to have BSC standard parameter names as
        # keys (rather than GHOST)
        self.parameter_dictionary = dict()
        for _, param_dict in standard_parameters.items():
            self.parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict
        # get standard metadata dictionary
        self.standard_metadata = get_standard_metadata({'standard_units': ''})
        # create list of metadata variables to read (make global)
        self.metadata_vars_to_read = [key for key in self.standard_metadata.keys() if
                                      pd.isnull(self.standard_metadata[key]['metadata_type']) == False]
        self.metadata_dtype = [(key, self.standard_metadata[key]['data_type']) for key in
                               self.metadata_vars_to_read]
        self.standard_data_flag_name_to_data_flag_code = \
            standard_data_flag_name_to_data_flag_code
        self.standard_QA_name_to_QA_code = standard_QA_name_to_QA_code
        self.qa_exceptions = ['dir10', 'spd10', 'rho2', 'acprec', 'acsnow', 'si',
                              'cldbot', 'vdist', 'ccovmean', 'cfracmean']
        self.specific_qa, self.general_qa, self.qa_diff = aux.get_qa_codes(self)
