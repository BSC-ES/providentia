""" Class for the conversion of units """

from collections import Counter
import copy
import math
import re
from string import digits
import sys

import numpy as np
import pandas as pd
    
class UnitConverter:
    """
    Python class that converts units (specifically tailored for conversions of components in the field of 
    Atmospheric Composition). 

    Accepts mixture of SI base/derived/equivalent units and non-SI units as input/output arguments, and
    can handle complex formatted input/output units (i.e. without spaces, spaced, including '/' or 'per' terms).

    Converts units between same quantities (i.e. length --> length, cm --> m), 
    and converts between different quantities (i.e. mass concentration --> mole fraction, mg/m3 --> ppbv), 
    when given the correct input parameters, and have a set conversion formula.

    Takes the following input arguments:

    1. input_units -- string of input units for same quantity conversion, or a dictionary of input quantities and units 
    for converting between quantities. Quantity names should be given with underscore joins when more than 1 word. 
    e.g. 'kg m-2 s-1'
    e.g. {'molar_mass':'g mol-1', 'mole_fraction':'nmol mol-1'}

    2. output_units -- string of output units
    e.g. 'kg kg-1'

    3. input_value -- input value for same quantity conversion, or a dictionary of input quantities and values 
    for converting between quantities. Quantity names should be given with underscore joins when more than 1 word.
    e.g. 28.0
    e.g. {'temperature':273.15, 'pressure':1013.25} 

    4. precision **OPTIONAL -- precision of converted output / conversion factor (fixed number of decimal places).
    e.g. precision=2

    5. species **OPTIONAL -- species (in chemical notation), needed when converting units associated with an 
    element or molecule (e.g. ppmv N). Can pass pass molecular charges for ions when needed (e.g. NO3- or SO4--).  
    e.g. species='O3'
    e.g. species='SO4--'
    
    6. valence **OPTIONAL -- valence of species, needed when converting some units associated with ions (e.g. eq l-1).
    e.g. valence=1

    7. input_quantity **OPTIONAL -- input quantity, needed when wanting to return conversion_factor for conversions
    with more than 1 input variable, and for cases where multiple quantities can be represented by the same units.
    The quantity refers to the input variable for which a conversion factor to the output units is wanted.
    e.g. input_quantity='mass_concentration'

    8. output_quantity **OPTIONAL -- output quantity, needed for cases where multiple quantities can be represented 
    by the same units.
    e.g. output_quantity='momentum'

    ----------------------------------------------------------------------------------------
    
    How to use class

    -----------------------------------------------------
    Same quantity unit conversion (e.g. m3 --> dm3)

    1. Define input variables 
    input_units = 'm3'
    output units = 'dm3'
    input_value = 1.0

    2. Do conversion
    conv_obj = unit_converter.UnitConverter(input_units, output_units, input_value)

    -----------------------------------------------------
    Different unit quantity conversion (e.g. ug m-3 --> ppbv)

    1. Define input variables 
    input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mass_concentration':'ug m-3'}
    output_units = 'ppbv'
    input_values = {'temperature':273.0, 'pressure':1013.0, 'molar_mass':48.0, 'mass_concentration':10.0}

    2. Do conversion
    conv_obj = unit_converter.UnitConverter(input_units, output_units, input_values, 
                                            input_quantity='mass_concentration')

    ----------------------------------------------------------------------------------------
    Output Class Attributes

    conv_obj.converted_value                   --> Converted output value 
    conv_obj.conversion_factor                 --> Conversion factor to get from input value to output value (Note this is not valid for the temperature quantity conversions)
    conv_obj.input_cleaned_units               --> Cleaned version of input units 
    conv_obj.output_cleaned_units              --> Cleaned version of output units 
    conv_obj.input_reference_units             --> Input units as reference units (equivalent SI base units for given quantity)
    conv_obj.input_reference_simplified_units  --> Input units as simplified reference units (equivalent SI base units for given quantity)   
    conv_obj.output_reference_units            --> Output units as reference units (equivalent SI base units for given quantity)
    conv_obj.output_reference_simplified_units --> Output units as simplified reference units (equivalent SI base units for given quantity)
    conv_obj.input_standard_units              --> Input units as standard units (standardised form of units)
    conv_obj.input_standard_simplified_units   --> Input units as simplified standard units (standardised form of units)
    conv_obj.output_standard_units             --> Output units as standard units (standardised form of units)
    conv_obj.output_standard_simplified_units  --> Output units as simplified standard units (standardised form of units)
    conv_obj.input_quantity                    --> Quantity of input units 
    conv_obj.input_quantities                  --> Quantities of input units 
    conv_obj.output_quantity                   --> Quantity of output units
    """

    def __init__(self, input_units, output_units, input_value, precision=-999, species='', valence=np.nan,
                 input_quantity='', output_quantity=''):
        """
        Initialises the UnitConverter and executes the conversion logic based on provided units and values.
        
        Parameters
        ----------
        input_units : str or dict
            Input units string or dictionary of quantities and units for inter-quantity conversion.
        output_units : str
            Desired output units string.
        input_value : float or dict
            Input value or dictionary of values corresponding to input quantities.
        precision : int, optional
            Fixed number of decimal places for the result.
        species : str, optional
            Chemical species (e.g., 'O3') for mass-to-mole conversions.
        valence : int, optional
            Ion valence for equivalent unit conversions.
        input_quantity : str, optional
            Specific quantity being converted when multiple interpretations exist.
        output_quantity : str, optional
            Specific target quantity when multiple interpretations exist.
        """
        
        # set input units to self 
        if isinstance(input_units, dict):
            self.input_units = {k: str(v) for k, v in input_units.items()}
        else:
            self.input_units = str(input_units)
        
        # set output units to self 
        self.output_units = str(output_units)
        
        # set input value to self 
        if isinstance(input_value, dict):
            self.input_value = {k: float(v) for k, v in input_value.items()}
        else:
            self.input_value = float(input_value)
        
        # set precision to self 
        self.precision = int(precision)
        # set species to self
        self.species = str(species)
        # set valence to self
        if np.isnan(valence):
            self.valence = get_valence(self.species)
        else:
            self.valence = int(valence)
        
        # set input/output quantity to self
        self.user_input_quantity = str(input_quantity)
        self.user_output_quantity = str(output_quantity)

        # load key dictionaries into self
        self.load_dictionaries()

        # make conversion
        self.do_conversion()
    
    def do_conversion(self):
        """
        Executes the unit conversion logic, branching between scaling within the same quantity or calculating between different quantities.
        """

        # get prefix_standard names
        standard_prefix_names = []
        for prefix in self.standard_prefixes:
            standard_prefix_names.append(self.standard_prefixes[prefix]['standard_units'])
        
        # determine what type of conversion is being done
        # if input units are a string then want to convert between same quantity
        if isinstance(self.input_units, str):
            self.conversion_type = 'same_quantities'
        # elif input units are a dict then want to convert between quantities
        elif isinstance(self.input_units, dict):
            self.conversion_type = 'between_quantities'

        # to convert between same quantity, get conversion factor to convert between given input and output units to standard SI
        # the ratio between these conversion factors * the input value then gives the desired output value in the desired units   
        if self.conversion_type == 'same_quantities':

            # get input/output units into standardised format and get unified conversion factor to base SI
            self.standardise_unit_form(self.input_units, 'input')
            self.standardise_unit_form(self.output_units, 'output')

            # checks if units are of same quantity and can be converted between
            self.check_quantity_agreement()

            # get scaling ratio between input and output units
            self.scaling_factor = self.input_scaling_factor_to_reference/self.output_scaling_factor_to_reference
            # print('input/output scaling to reference :', self.input_scaling_factor_to_reference,self.output_scaling_factor_to_reference) 
            # print('scaling factor :', self.scaling_factor)

            # correct for temperature scaling, only scale for absolute temps. 
            if self.output_quantity == 'temperature': 
                if (self.input_standard_units[0] not in standard_prefix_names) & (self.output_standard_units[0] not in standard_prefix_names):
                    self.correct_temperature(self.input_value,self.input_standard_units,self.output_standard_units)
        
            # multiply input value by final calculated scaling factor 
            else:                  
                self.converted_value = self.input_value*self.scaling_factor
                # save out conversion factor 
                self.conversion_factor = self.scaling_factor
    
        #--------------------------------------------------------------
        # to convert between quantities, process each of the input units into standard SI units, 
        # and then calculate output value using conversion formula
        elif self.conversion_type == 'between_quantities':

            # put input quantity names into list
            self.input_quantities = list(self.input_units.keys())
            # add mathematical constants to input quantities
            for m_quantity in self.constants:
                self.input_quantities.append(m_quantity)
    
            # get output units into standardised format and get unified conversion factor to base SI
            self.standardise_unit_form(self.output_units, 'output')
 
            # check that have everything correctly configured to convert between input variables --> output variable
            self.check_conversion_validity()

            # list for appending input values, converted to SI
            self.SI_input_quantities = {}
    
            # add mathematical constants to SI_input_quantities
            for m_quantity in self.constants:
                self.SI_input_quantities[m_quantity] = self.constants[m_quantity]['value']
    
            # iterate through input units, getting into standardised format and get unified conversion factor and 
            # get unified conversion factor to base SI
            for quantity in self.input_quantities:
                if quantity not in self.constants:

                    # standardise unit form
                    self.standardise_unit_form(self.input_units[quantity], 'input')
    
                    # scaling factor
                    self.scaling_factor = self.input_scaling_factor_to_reference

                    # save scaling factor if want to return conversion factor
                    if self.user_input_quantity == quantity:
                        conversion_scaling_factor_1 = self.scaling_factor

                    # correct for temperature scaling to SI Kelvin
                    if quantity == 'temperature':
                        if self.input_standard_units[0] not in standard_prefix_names:
                            self.correct_temperature(self.input_value[quantity], self.input_standard_units,'K')
            
                    # convert to base SI
                    else:
                        self.converted_value = self.input_value[quantity]*self.scaling_factor

                    # add SI input parameter name:converted_value to dict
                    self.SI_input_quantities[quantity] = self.converted_value
    
            # if have a conversion input quantity, return conversion factor to get output quantity 
            # by: input_quantity*conversion_factor
            if self.user_input_quantity != '':
                formula = self.get_conversion_formula(formula_type='conversion_factor')
                conversion_scaling_factor_2 = eval(formula,self.SI_input_quantities)
                self.conversion_factor = conversion_scaling_factor_1 * conversion_scaling_factor_2
            
            # take converted SI quantities and calculate SI output quantity using appropriate formula
            self.converted_value = eval(self.formula,self.SI_input_quantities)
    
            # finally, convert output SI quantity into desired output units (if necessary)
            if self.output_standard_units != self.output_reference_units:
        
                self.standardise_unit_form(self.output_reference_units, 'input')

                # get scaling ratio between input (SI units) and output units
                self.scaling_factor = self.input_scaling_factor_to_reference/self.output_scaling_factor_to_reference
        
                # save scaling factor if want to return conversion factor
                if self.user_input_quantity != '':
                    conversion_scaling_factor_3 = self.scaling_factor
                    self.conversion_factor = self.conversion_factor * conversion_scaling_factor_3
            
                # correct for temperature scaling from K
                if self.output_quantity == 'temperature': 
                    if self.output_standard_units[0] not in standard_prefix_names:
                        self.correct_temperature(self.converted_value, 'K', self.output_standard_units)
                # convert to desired output units
                else:
                    self.converted_value = self.converted_value*self.scaling_factor
    
        # finish by rounding converted value/conversion factor to fixed number of decimal places (if desired)
        if self.precision != -999:
            self.converted_value = np.around(self.converted_value, self.precision)
            try:
                self.conversion_factor = np.around(self.conversion_factor, self.precision)
            except:
                pass
            
    def check_quantity_agreement(self):
        """
        Verifies that the input and output units represent the same physical quantity before proceeding with scaling.
        """

        #print('cleaned units are: {} --> {}'.format(self.input_cleaned_units,self.output_cleaned_units))
        #print('reference units are: {} --> {}'.format(self.input_reference_units,self.output_reference_units))
        #print('standard units are: {} --> {}'.format(self.input_standard_units,self.output_standard_units))
        
        # get SI quantities that are represented by input/output units
        self.input_quantity = self.get_quantity_representation(self.input_reference_units, self.input_reference_simplified_units, 'input', self.user_input_quantity)
        self.output_quantity = self.get_quantity_representation(self.output_reference_units, self.output_reference_simplified_units, 'output', self.user_output_quantity)
        
        # are quantities represented same between input/output?
        if self.input_quantity != self.output_quantity:
            sys.exit("Input and output quantities do not agree: '{}', '{}'. Check input and output units.".format(self.input_quantity, self.output_quantity))
            
    def check_conversion_validity(self):
        """
        Validates that the input dictionaries are the same length and retrieves the appropriate formula for inter-quantity conversion.
        """
    
        # check if input unit/value dicts have the same length
        if len(list(self.input_units.keys())) != len(list(self.input_value.keys())):
            sys.exit("Input unit and value dictionaries do not have same length.")
    
        # check that input unit/value dicts have same quantity names
        if set(list(self.input_units.keys())) != set(list(self.input_value.keys())):
            sys.exit("Input unit and value dictionaries do not have same defined quantity names.")

        # get SI quantities represented by output_units
        self.output_quantity = self.get_quantity_representation(self.output_reference_units, self.output_reference_simplified_units, 'output', self.user_output_quantity)

        # get conversion formula between input and output quantities
        self.formula = self.get_conversion_formula()
             
    def get_quantity_representation(self, units, units_simplified, unit_type, unit_quantity):
        """
        Identifies the physical quantity name associated with a set of reference units by checking SI base and derived definitions.
        """

        # set list for storing valid quantities
        valid_quantities = []

        # is highest level quantity represented a base SI unit?
        for quantity in self.SI_base_quantities:
            if (Counter(units) == Counter(self.SI_base_quantities[quantity]['base_units'])) or (Counter(units_simplified) == Counter(self.SI_base_quantities[quantity]['base_units'])):
                valid_quantities.append(quantity)
                
        # is highest level quantity represented a derived SI unit?
        for quantity in self.SI_derived_quantities:
            if (Counter(units) == Counter(self.SI_derived_quantities[quantity]['base_units'])) or (Counter(units_simplified) == Counter(self.SI_derived_quantities[quantity]['base_units'])):
                valid_quantities.append(quantity)
            
        # is highest level quantity represented a scientific constant?
        for quantity in self.constants:
            if (Counter(units) == Counter(self.constants[quantity]['base_units'])) or (Counter(units_simplified) == Counter(self.constants[quantity]['base_units'])):
                valid_quantities.append(quantity)

        # if have just 1 valid quantity, return that
        if len(valid_quantities) == 1:
            return valid_quantities[0] 
        
        # elif have more than 1 valid quantity, if user has given a input/output quantity, use that to select correct matching quantity
        elif len(valid_quantities) > 1:
            # if user has provided a quantity, and that is a valid quantity, then siomply take that
            if unit_quantity in valid_quantities:
                return unit_quantity
            # elif user has not provided a quantity the operation depends if are converting between the same quantity or not
            elif unit_quantity == '':
                # if have just unitless and another quantity, just take the non-unitless option
                if (len(valid_quantities) == 2) & ('unitless' in valid_quantities):
                    valid_quantities.remove('unitless')
                    return valid_quantities[0]
                # if are converting between same quantities, as cannot give a specific quantity, simply return simplified units as quantity
                elif self.conversion_type == 'same_quantities':
                    return units_simplified
                # if are converting between different quantities then throw error as the specific quantity is needed to be known for conversion 
                else:
                    sys.exit("Multiple valid quantities found for {} units: '{}'. Valid quantites are: {}. You must specify {}_quantity.".format(unit_type, units, valid_quantities, unit_type))
            else:
                sys.exit("Specified {}_quantity: '{}' not valid for {}_units: '{}'. Valid quantities are: {}.".format(unit_type, unit_quantity, unit_type, units, valid_quantities))

        # elif unable to find quantity for reference units, simply return simplified units as quantity
        elif len(valid_quantities) == 0:
            return units_simplified

    def get_conversion_formula(self, formula_type='standard'):
        """ 
        Retrieves the appropriate mathematical formula for either standard output calculation or conversion factor derivation.

        Parameters
        ----------
        formula_type : str
            Type of formula to retrieve; either 'standard' to get the formula to obtain the output quantity, 
            or 'conversion_factor' to get the formula to obtain the conversion factor

        Returns
        -------
        formula : str
            The string representation of the mathematical formula to be evaluated.
        """

        # is the conversion formula standard (i.e. formula to get just the output quantity)
        if formula_type == 'standard':                   
            # do you have a formula for converting to desired output quantity?
            available_conversion_formulae_quantities = list(self.conversion_formulae.keys())
    
            # if don't have formula, return error message
            if self.output_quantity not in available_conversion_formulae_quantities:
                sys.exit("Do not have formula to calculate desired output quantity: '{}'.".format(self.output_quantity))
            
            else:
                # check if have given required input variables for calculating output quantity using one of formulae
                # defined in conversion_formulae dictionary
                quantity_formulae = self.conversion_formulae[self.output_quantity]
                remove_digits = str.maketrans('', '', digits)
                for formula in quantity_formulae:
                    formula_quantities = formula.translate(remove_digits).replace('*','$').replace('/','$').replace('(','').replace(')','').replace('.','').strip().split('$')
                    formula_quantities = [formula_quantity for formula_quantity in formula_quantities if formula_quantity]
                    # test to see if have all of formula quantities required in input quantities
                    have_quantities = set(formula_quantities).issubset(self.input_quantities)
                    if have_quantities:
                        return formula

                # if do not valid input quantities for one of the formulae, return error message
                sys.exit("Do not have correct input quantities to calculate desired output quantity: '{}'.".format(self.output_quantity))

        # is the conversion formula a factor (i.e. formula to get factor between an input and output quantity)
        elif formula_type == 'conversion_factor':                   
    
            # do you have a formula for getting conversion_factor
            available_conversion_formulae_quantities = list(self.conversion_formulae.keys())
    
            # if don't have formula, return error message
            if ('{}>{}--conversion_factor'.format(self.user_input_quantity, self.output_quantity) 
                not in available_conversion_formulae_quantities):
                sys.exit("Do not have defined formula to calculate conversion factor between input quantity: '{}' and output quantity: '{}'.".format(self.user_input_quantity, self.output_quantity))
            
            else:
                # check if have given required input variables for calculating conversion_factor using one of formulae 
                # defined in conversion_formulae dictionary
                quantity_formulae = self.conversion_formulae['{}>{}--conversion_factor'.format(self.user_input_quantity,self.output_quantity)]
                remove_digits = str.maketrans('', '', digits)
                for formula in quantity_formulae:
                    formula_quantities = formula.translate(remove_digits).replace('*','$').replace('/','$').replace('(','').replace(')','').replace('.','').strip().split('$')
                    formula_quantities = [formula_quantity for formula_quantity in formula_quantities if formula_quantity]
                    # test to see if have all of formula quantities required in input quantities
                    have_quantities = set(formula_quantities).issubset(self.input_quantities)
                    if have_quantities:
                        return formula

                # if do not valid input quantities for one of the formulae, return error message
                sys.exit("Do not have correct input quantities to calculate conversion factor between input quantity: '{}' and output quantity: '{}'.".format(self.user_input_quantity, self.output_quantity))
             
    def correct_temperature(self, input_value, input_standard_unit, output_standard_unit):
        """
        Adjusts temperature values by applying both multiplicative scaling and additive offsets for different thermometric scales.
        """

        # kelvin --> kelvin
        if ('K' in input_standard_unit) & ('K' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
            
        # celsius --> celsius
        if ('°C' in input_standard_unit) & ('°C' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
        
        # fahrenheit --> fahrenheit
        if ('°F' in input_standard_unit) & ('°F' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
        
        # rankine --> rankine
        if ('°R' in input_standard_unit) & ('°R' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
        
        # kelvin --> fahrenheit
        if ('K' in input_standard_unit) & ('°F' in output_standard_unit):
            self.converted_value = (input_value * self.scaling_factor) - 459.67
        
        # kelvin --> rankine
        if ('K' in input_standard_unit) & ('°R' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
        
        # kelvin --> celsius
        if ('K' in input_standard_unit) & ('°C' in output_standard_unit):
            self.converted_value = input_value - 273.15
        
        # celsius --> kelvin
        if ('°C' in input_standard_unit) & ('K' in output_standard_unit):
            self.converted_value = input_value + 273.15
            
        # celsius --> fahrenheit
        if ('°C' in input_standard_unit) & ('°F' in output_standard_unit):
            self.converted_value = (input_value * self.scaling_factor) + 32.0
            
        # celsius --> rankine
        if ('°C' in input_standard_unit) & ('°R' in output_standard_unit):
            self.converted_value = (input_value + 273.15) * self.scaling_factor
        
        # fahrenheit --> kelvin
        if ('°F' in input_standard_unit) & ('K' in output_standard_unit):
            self.converted_value = (input_value + 459.67) * self.scaling_factor
            
        # fahrenheit --> celsius
        if ('°F' in input_standard_unit) & ('°C' in output_standard_unit):
            self.converted_value = (input_value - 32.0) * self.scaling_factor
    
        # fahrenheit --> rankine
        if ('°F' in input_standard_unit) & ('°R' in output_standard_unit):
            self.converted_value = input_value + 459.67
            
        # rankine --> kelvin
        if ('°R' in input_standard_unit) & ('K' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
    
        # rankine --> celsius
        if ('°R' in input_standard_unit) & ('°C' in output_standard_unit):
            self.converted_value = (input_value - 491.67) * self.scaling_factor
            
        # rankine --> fahrenheit
        if ('°R' in input_standard_unit) & ('°F' in output_standard_unit):
           self.converted_value = input_value - 459.67
           
        return
    
    def standardise_unit_form(self, units, unit_type):
        """
        Decomposes and standardises raw unit strings into cleaned, reference SI, and formatted standard forms while calculating scaling factors.

        Parameters
        ----------
        units : str
            The raw unit string to be processed (e.g., 'ug/m3' or 'kg m-2 s-1').
        unit_type : str
            Indicator of whether the units are 'input' or 'output'.
        """

        # strip all leading/trailing white spaces from string
        units = units.strip()

        # throw error if units string is empty
        if units == '':
            sys.exit("The {} units are empty. If you want to set empty units, set it as: 'unitless'".format(unit_type))

        # ensure mass densities and mole fractions per element are correctly formatted
        for element in ['C','CaCO3','Cl','N','S']:
            units = units.replace('g {}/m3'.format(element), 'g{}/m3'.format(element))  
            units = units.replace('g {} m-3'.format(element), 'g{} m-3'.format(element))
            units = units.replace('mol {}/mol'.format(element), 'mol{}/mol'.format(element))  
            units = units.replace('mol {} mol-1'.format(element), 'mol{} mol-1'.format(element))
            units = units.replace('g {}/l'.format(element), 'g{}/l'.format(element))  
            units = units.replace('g {} l-1'.format(element), 'g{} l-1'.format(element))         

        # list for appending reference units (to which all units are scaled around)
        reference_units = []
        
        # list for appending standardised units
        standard_units = []
        
        # list for appending scaling factors for each unit to get to equivalent SI base unit  
        scaling_factors_to_reference = []

        # list for storing dependencies for calculation of non-SI units scaling factors
        non_SI_unit_dependencies = []
         
        #-----------------------------------------------------
        # remove '/' or 'per' strings from input unit string
        # track all units which occur after '/' or 'per' strings as being in the denominator
        if '/' in units:
            split_per =  units.split('/')
        else:
            split_per =  units.split('per')
        
        # join split_per back together to 1 string (ensuring there is 1 character whitespace between separate units)
        split_per = [u.strip() for u in split_per]
        units = ' '.join(split_per)
        
        # get denominator indices
        denominator_inds = list(range(len(split_per[0])+1,len(units))) 

        #-----------------------------------------------------
        # separate units string into separate distinct units (i.e. separate quantities)
        
        # pass through all combinations of SI/non-SI quantities/+prefix units and find all matches with unit string
        # take longest string matches out iteratively of string until original string is empty
        
        match_list = []

        # SI base units
        for quantity in self.SI_base_quantities:
        
            for prefix in self.standard_prefixes:
                for prefix_variant in self.standard_prefixes[prefix]['units']:
                    # does unit string contain one of the SI base units + prefix?
                    if prefix_variant+self.SI_base_quantities[quantity]['base_units'] in units:
                        match_list.append(prefix_variant+self.SI_base_quantities[quantity]['base_units'])
                        scaling_factors_to_reference.append(self.standard_prefixes[prefix]['factor']) 
                        reference_units.append(self.SI_base_quantities[quantity]['base_units'])
                        standard_units.append(self.standard_prefixes[prefix]['standard_units']+self.SI_base_quantities[quantity]['base_units'])
                        non_SI_unit_dependencies.append([])
                
                    for equiv_unit in self.SI_base_quantities[quantity]['equiv_units']:
                        # does unit string contain one of the SI base equivalent units + prefix?
                        if prefix_variant+equiv_unit in units:
                            match_list.append(prefix_variant+equiv_unit)
                            scaling_factors_to_reference.append(self.standard_prefixes[prefix]['factor']) 
                            reference_units.append(self.SI_base_quantities[quantity]['base_units'])
                            standard_units.append(self.standard_prefixes[prefix]['standard_units']+self.SI_base_quantities[quantity]['base_units'])
                            non_SI_unit_dependencies.append([])
            
        # SI derived equivalent units
        for quantity in self.SI_derived_quantities:

            for equiv_unit in self.SI_derived_quantities[quantity]['equiv_units']:
                for prefix in self.standard_prefixes:
                    for prefix_variant in self.standard_prefixes[prefix]['units']:
                        # does unit string contain one of the SI derived equivalent units + prefix?
                        if prefix_variant+equiv_unit in units:
                            match_list.append(prefix_variant+equiv_unit)
                            scaling_factors_to_reference.append(self.standard_prefixes[prefix]['factor'])
                            reference_units.append(self.SI_derived_quantities[quantity]['base_units'])
                            standard_units.append(self.standard_prefixes[prefix]['standard_units']+self.SI_derived_quantities[quantity]['standard_units'])
                            non_SI_unit_dependencies.append([])

        # non-SI units
        for quantity in self.non_SI_quantities:

            for quantity_variant in self.non_SI_quantities[quantity]:
                for unit in self.non_SI_quantities[quantity][quantity_variant]['units']:                        
                    for prefix in self.standard_prefixes:
                        # do not test 'grams' with prefix 'kilo', as already have that as base unit (do not want double matches)
                        if (quantity_variant == 'gram') & (prefix == 'kilo'):
                            continue
                        for prefix_variant in self.standard_prefixes[prefix]['units']:
                            # does unit string contain one of the non-SI units+prefix?
                            if prefix_variant+unit in units:
                                match_list.append(prefix_variant+unit)
                                scaling_factors_to_reference.append(self.non_SI_quantities[quantity][quantity_variant]['factor'] * self.standard_prefixes[prefix]['factor'])
                                standard_units.append(self.standard_prefixes[prefix]['standard_units']+self.non_SI_quantities[quantity][quantity_variant]['standard_units'])
                                try:
                                    reference_units.append(self.SI_base_quantities[quantity]['base_units'])
                                except:
                                    reference_units.append(self.SI_derived_quantities[quantity]['base_units'])
                                if 'dependencies' in self.non_SI_quantities[quantity][quantity_variant]:
                                    non_SI_unit_dependencies.append(self.non_SI_quantities[quantity][quantity_variant]['dependencies'])
                                else:
                                    non_SI_unit_dependencies.append([])

        # if have no matches then exit
        if len(match_list) == 0:
            sys.exit("The {} units: '{}' are not recognised. This may be due to incorrect capitalisations or spacings.".format(unit_type, units))

        # sort match list by length of string (biggest first), then alphabetically
        # also sort scaling factors/unit lists with match list 
        match_list, reference_units, standard_units, scaling_factors_to_reference, non_SI_unit_dependencies = map(list, zip(*sorted(zip(match_list, 
            reference_units, standard_units, scaling_factors_to_reference, non_SI_unit_dependencies), key=lambda item: (-len(item[0]), item[0]))))

        # lists for storing accepted match units
        accepted_cleaned_units = []
        accepted_reference_units = []
        accepted_standard_units = []
 
        # list to store unit multiplication scaling factors per matched unit 
        mult_factors = []

        # list for storing lists of inds for each separate accepted matched unit in original input string   
        accepted_unit_inds_in_str = []      
 
        # create separate copy of original units for cutting down iteratively
        cut_units = copy.deepcopy(units)
        
        # create array for all indices which remain unmatached in cut_units 
        unit_unmatched_inds = list(range(len(cut_units))) 

        # iteratively get inds of matches in original unit string
        for mm, match_unit in enumerate(match_list):

            # only proceed if match still exists in cut original input units
            if match_unit in cut_units:
                # get all inds of match occurances in unit_string
                match_inds = [i for i in list(range(len(units))) if units.startswith(match_unit, i)]

                # iterate through match inds, if ind not already part of another distinct unit string then proceed
                for match_ind in match_inds:

                    # set str for match unit with exponent 
                    match_unit_with_exponent = copy.deepcopy(match_unit)

                    if not any(match_ind in x for x in accepted_unit_inds_in_str):
                        relevant_inds = list(range(match_ind, match_ind+len(match_unit)))
                        # evaluate if next characters are numeric after relevant inds (i.e exponent)
                        match_unit_exponent = ''
                        while True:
                            next_relevant_ind = relevant_inds[-1]+1
                            if next_relevant_ind >= len(units):
                                break    
                            if (units[next_relevant_ind] == '-') or (units[next_relevant_ind].isdigit()):
                                relevant_inds.append(next_relevant_ind)
                                match_unit_with_exponent += units[next_relevant_ind]
                            else:
                                break 

                        # remove matched inds from unit_unmatched_inds array
                        unit_unmatched_inds = sorted([x for x in unit_unmatched_inds if x not in relevant_inds])                        
                        
                        # remove matched units from original units string            
                        cut_units = ''
                        for x in unit_unmatched_inds:
                            cut_units += units[x]
                        
                        # append indices of matched unit indices in original input string
                        accepted_unit_inds_in_str.append(relevant_inds)
 
                        # create dict of relevant cleaned/reference/standard units
                        relevant_units = {'clean':match_unit_with_exponent, 'reference':reference_units[mm], 'standard':standard_units[mm]}

                        # throw error if have non-SI unit dependencies not met
                        dependencies = non_SI_unit_dependencies[mm]
                        if ('species' in dependencies) & ('valence' in dependencies):
                            if (self.species == '') & (np.isnan(self.valence)):
                                sys.exit(" The {} units: '{}' requires both 'species' and 'valence' to be passed as input arguments.".format(unit_type, units))
                        elif 'species' in dependencies:
                            if self.species == '':
                                sys.exit("The {} units: '{}' requires 'species' to be passed as an input argument.".format(unit_type, units)) 
                        elif 'valence' in dependencies:
                            if np.isnan(self.valence):
                                sys.exit("The {} units: '{}' requires 'valence' to be passed as an input argument.".format(unit_type, units))  

                        # evaluate if all of relevant indices are in denominator indices list
                        if all(elem in denominator_inds for elem in relevant_inds): 
                            # if so, then invert all matched units exponents (as in denominator) in cleaned units, reference and standard units 
                            invert_exponents = True
                        # if only some of relevant indices are in denominator indices, then throw error as all of them should be
                        elif any(elem in denominator_inds for elem in relevant_inds):
                            sys.exit("The {} units: '{}' syntax is invalid.".format(unit_type, units))
                        # else, do not need to invert exponents
                        else:
                            invert_exponents = False            

                        # iterate through all cleaned/reference/standard units character by character to set all exponents correctly
                        # this involves raising the reference/standard units to the power of the cleaned unit exponent
                        # also this involves inverting exponents if units are in the denominator, adding -1 exponents if no exponents present
                        for relevant_unit_key in ['clean','reference','standard']:
                            relevant_unit = relevant_units[relevant_unit_key]

                            # if relevant unit is a compound with a number then do not check for exponent
                            if 'CaCO3' in relevant_unit:
                                if relevant_unit_key == 'clean':
                                    match_unit_exponent = 1
                                continue

                            # create empty strings for appending new altered units
                            new_unit = ''
                            current_symbol = ''
                            current_exponent = ''

                            # iterate through characters in unit string
                            for char in relevant_unit:
                                # is current character a space?
                                if char == ' ':
                                    # if current symbol is not empty string, append it to new_unit str
                                    if current_symbol != '':
                                        new_unit+=current_symbol
                                        current_symbol = ''
                                        # determine whether to add exponent to previous symbol (did not have one)
                                        # handling inverting exponents, and raising to power of cleaned unit exponent
                                        if relevant_unit_key == 'clean': 
                                            if invert_exponents:
                                                new_exponent = -1
                                            else:
                                                new_exponent = 1
                                        else:
                                            new_exponent = copy.deepcopy(match_unit_exponent)
                                        # if new exponent != 1, then add it
                                        if new_exponent != 1:
                                            new_unit+=str(new_exponent)
                                    # if current exponent is not empty string, determine whether to invert it, 
                                    # whether to raise it to cleaned unit exponent, and then append it to new_unit str (if != 1)
                                    if current_exponent != '': 
                                        if relevant_unit_key == 'clean':
                                            if invert_exponents:
                                                new_exponent = -int(current_exponent)
                                            else:
                                                new_exponent = int(current_exponent)
                                        else:
                                            new_exponent = int(current_exponent)*match_unit_exponent
                                        # if new exponent != 1, then add it
                                        if new_exponent != 1:
                                            new_unit+=str(new_exponent)
                                        current_exponent = ''
                                    # append space to new_unit str
                                    new_unit+=' '

                                # determine if current character is part of exponent
                                elif (char == '-') or (char.isdigit()):
                                    # if current symbol is not empty string, append it to new_unit str  
                                    if current_symbol != '':
                                        new_unit+=current_symbol
                                        current_symbol = ''
                                    # append exponent character
                                    current_exponent+=char

                                # otherwise current character is symbol
                                else:
                                    # if current exponent is not empty string, determine whether to invert it,
                                    # whether to raise it to cleaned unit exponent, and then append it to new_unit str (if != 1)
                                    if current_exponent != '':
                                        if relevant_unit_key == 'clean':
                                            if invert_exponents:
                                                new_exponent = -int(current_exponent)
                                            else:
                                                new_exponent = int(current_exponent)
                                        else:
                                            new_exponent = int(current_exponent)*match_unit_exponent
                                        # if new exponent != 1, then add it
                                        if new_exponent != 1:
                                            new_unit+=str(new_exponent)
                                        current_exponent = ''
                                    # append symbol character
                                    current_symbol+=char

                            # if have current symbol remaining that have not been handled, append it to new_unit str
                            if current_symbol != '':
                                new_unit+=current_symbol
                                # determine whether to add exponent to previous symbol (did not have one)
                                # handling inverting exponents, and raising to power of cleaned unit exponent
                                if relevant_unit_key == 'clean':
                                    if invert_exponents:
                                        new_exponent = -1
                                    else:
                                        new_exponent = 1
                                else:
                                    new_exponent = copy.deepcopy(match_unit_exponent)
                                # if new exponent != 1, then add it
                                if new_exponent != 1:
                                    new_unit+=str(new_exponent)

                            # if have current exponent remaining that have not been handled, determine whether to invert it,
                            # whether to raise it to cleaned unit exponent, and then append it to new_unit str (if != 1)
                            if current_exponent != '':
                                if relevant_unit_key == 'clean':
                                    if invert_exponents:
                                        new_exponent = -int(current_exponent)
                                    else:
                                        new_exponent = int(current_exponent)
                                else:
                                    new_exponent = int(current_exponent)*match_unit_exponent
                                # if new exponent != 1, then add it
                                if new_exponent != 1:
                                    new_unit+=str(new_exponent)  

                            # if relevant_unit_key is == 'clean', set match unit exponent
                            if relevant_unit_key == 'clean':
                                match_unit_exponent = copy.deepcopy(new_exponent)

                            # write modified new_unit back into relevant_units dict
                            relevant_units[relevant_unit_key] = copy.deepcopy(new_unit)

                        # append now fully corrected relevant units to accepted units lists
                        accepted_cleaned_units.append(relevant_units['clean'])
                        accepted_reference_units.append(relevant_units['reference'])
                        accepted_standard_units.append(relevant_units['standard'])

                        # generate new scaling factor to reference, accounting for exponent 
                        # (relevant_scaling_factor**exponent)
                        relevant_scaling_factor_to_reference = scaling_factors_to_reference[mm]
                        # if exponent == 1, do no modification to scaling factor, append to mult_factors
                        if match_unit_exponent == 1:
                            mult_factors.append(relevant_scaling_factor_to_reference) 
                        # otherwise generate new scaling factor and append to mult_factors
                        else:  
                            mult_factors.append(relevant_scaling_factor_to_reference**match_unit_exponent)

        # check all remaining characters from cut units
        # if string is not empty or all whitespace, units are not valid, returning invalid warning
        if (cut_units != '') & (not cut_units.isspace()):
            sys.exit("Parts of the {} units: '{}' are not recognised. This may be due to incorrect capitalisations or spacings.".format(unit_type, units))
   
        # sort each of list indices by first integer in each list
        # also sort scaling factors/unit lists equivalently
        accepted_unit_inds_in_str, accepted_cleaned_units, accepted_reference_units, accepted_standard_units, mult_factors = map(list, 
            zip(*sorted(zip(accepted_unit_inds_in_str, accepted_cleaned_units, accepted_reference_units, accepted_standard_units, mult_factors), 
            key=lambda x: x[0][0])))      

        # put together cleaned unit string (adding in spaces)
        cleaned_unit = ''
        for ii, accepted_cleaned_unit in enumerate(accepted_cleaned_units):
            cleaned_unit += accepted_cleaned_unit
            if ii != (len(accepted_cleaned_units)-1):
                cleaned_unit += ' '

        # put together reference unit string
        reference_unit = ' '.join(accepted_reference_units)

        # put together standard unit string
        standard_unit = ' '.join(accepted_standard_units)

        # get reference simplified units
        reference_simplified_unit = self.simplify_units(reference_unit)

        # get standard simplified units
        standard_simplified_unit = self.simplify_units(standard_unit)

        # multiply all scaling factors together in multiplication factor list
        scaling_factor_to_reference = np.prod(mult_factors)

        # if scaling_factor is NaN then it is because the species formula is needed to be input for conversion
        if pd.isnull(scaling_factor_to_reference):
            sys.exit("The {} units: '{}' requires 'species' to be passed as an input argument.".format(unit_type, units)) 

        # save output to self
        if unit_type == 'input':
            self.input_cleaned_units                 = cleaned_unit
            self.input_reference_units               = reference_unit
            self.input_reference_simplified_units    = reference_simplified_unit
            self.input_standard_units                = standard_unit
            self.input_standard_simplified_units     = standard_simplified_unit
            self.input_scaling_factor_to_reference   = scaling_factor_to_reference
        
        elif unit_type == 'output':
            self.output_cleaned_units                = cleaned_unit
            self.output_reference_units              = reference_unit
            self.output_reference_simplified_units   = reference_simplified_unit
            self.output_standard_units               = standard_unit
            self.output_standard_simplified_units    = standard_simplified_unit
            self.output_scaling_factor_to_reference  = scaling_factor_to_reference
                
        return 
        
    def simplify_units(self, units):
        """
        Reduces complex unit strings by summing exponents for identical base units and removing terms with zero exponents.

        Parameters
        ----------
        units : str
            The unit string to be simplified, containing space-separated terms with or without exponents.

        Returns
        -------
        simplified : str
            The simplified unit string maintaining the original order of first appearance.
        """

        exponents = {}
        order = []  # first occurrence order
        pattern = re.compile(r'(°?[a-zA-Z]+)([-+]?\d+)?')  # allow ° at start

        for part in units.split():
            match = pattern.fullmatch(part)
            if match:
                base = match.group(1)
                exp = int(match.group(2)) if match.group(2) else 1
                if base not in exponents:
                    order.append(base)
                exponents[base] = exponents.get(base, 0) + exp

        # Remove zero exponents
        exponents = {u: e for u, e in exponents.items() if e != 0}

        # Rebuild simplified string in original order
        simplified_parts = []
        for u in order:
            if u in exponents:
                part = f"{u}" if exponents[u] == 1 else f"{u}{exponents[u]}"
                simplified_parts.append(part)

        simplified = ' '.join(simplified_parts)
        return simplified

    def load_dictionaries(self):
        """
        Defines dictionaries for SI prefixes, SI base quantities, SI derived quantities, non-SI quantities, constants, 
        and conversion formulae to support unit parsing and conversion.
        """


        # dictionary of standard prefixes with associated multiplying factors
        self.standard_prefixes = {
                                  'yotta':{'factor':1e24,  'standard_units':'Y',  'units':['yotta','Yotta','Y']},
                                  'zetta':{'factor':1e21,  'standard_units':'Z',  'units':['zetta','Zetta','Z']},
                                  'exa':  {'factor':1e18,  'standard_units':'E',  'units':['exa','Exa','E']},
                                  'peta': {'factor':1e15,  'standard_units':'P',  'units':['peta','Peta','P']},
                                  'tera': {'factor':1e12,  'standard_units':'T',  'units':['tera','Tera','T']},
                                  'giga': {'factor':1e9,   'standard_units':'G',  'units':['giga','Giga','G']},
                                  'mega': {'factor':1e6,   'standard_units':'M',  'units':['mega','Mega','M']},
                                  'kilo': {'factor':1e3,   'standard_units':'k',  'units':['kilo','Kilo','k']},
                                  'hecto':{'factor':1e2,   'standard_units':'h',  'units':['hecto','Hecto','h']},
                                  'deca': {'factor':1e1,   'standard_units':'da', 'units':['deca','Deca','deka','Deka',
                                                                                           'da']},
                                  '':     {'factor':1.0,   'standard_units':'',   'units':['']},
                                  'deci': {'factor':1e-1,  'standard_units':'d',  'units':['deci','Deci','d']},
                                  'centi':{'factor':1e-2,  'standard_units':'c',  'units':['centi','Centi','c']},
                                  'milli':{'factor':1e-3,  'standard_units':'m',  'units':['milli','Milli','m']},
                                  'micro':{'factor':1e-6,  'standard_units':'u',  'units':['micro','Micro','µ','u']},
                                  'nano': {'factor':1e-9,  'standard_units':'n',  'units':['nano','Nano','n']},
                                  'pico': {'factor':1e-12, 'standard_units':'p',  'units':['pico','Pico','p']},
                                  'femto':{'factor':1e-15, 'standard_units':'f',  'units':['femto','Femto','f']},
                                  'atto': {'factor':1e-18, 'standard_units':'a',  'units':['atto','Atto','a']},
                                  'zepto':{'factor':1e-21, 'standard_units':'z',  'units':['zepto','Zepto','z']},
                                  'yocto':{'factor':1e-24, 'standard_units':'y',  'units':['yocto','Yocto','y']}}


        # define SI base quantities with units 
        self.SI_base_quantities = {
                                   'amount of substance':{'unit_name':'mole',
                                                          'base_units':'mol',
                                                          'equiv_units':['moles','mols','mole']},
                                   'electric current'   :{'unit_name':'ampere',
                                                          'base_units':'A',
                                                          'equiv_units':['amp','amps','amperes']},
                                   'length'             :{'unit_name':'metre',
                                                          'base_units':'m',
                                                          'equiv_units':['metres','meters','meter','metre']},
                                   'luminous intensity' :{'unit_name':'candela',
                                                          'base_units':'cd',
                                                          'equiv_units':['candela','candelas']},
                                   'mass'               :{'unit_name':'kilogram',
                                                          'base_units':'kg',
                                                          'equiv_units':['kilograms','kilogram','kg']},
                                   'temperature'        :{'unit_name':'kelvin',
                                                          'base_units':'K',
                                                          'equiv_units':['kelvin','Kelvin','kelvins','Kelvins']},
                                   'time'               :{'unit_name':'second',
                                                          'base_units':'s',
                                                          'equiv_units':['secs','seconds','sec','second']}}                     
 
        # define common SI derived quantities with units (formed by powers, products or quotients of the base units and 
        # are unlimited in number)
        self.SI_derived_quantities = {'unitless':                     {'unit_name':'unitless', 
                                                                       'base_units':'',
                                                                       'standard_units':'',
                                                                       'equiv_units':['unitless']},
                                      'absorbed_dose':                {'unit_name':'gray',                           
                                                                       'base_units':'m2 s-2',
                                                                       'standard_units':'Gy',
                                                                       'equiv_units':['Gy','gray','grays']},
                                      'absorption_coefficient':       {'unit_name':'reciprocal metre',
                                                                       'base_units':'m-1',
                                                                       'standard_units':'m-1',
                                                                       'equiv_units':['none']},
                                      'acceleration':                 {'unit_name':'metre per second squared',
                                                                       'base_units':'m s-2',
                                                                       'standard_units':'m s-2',
                                                                       'equiv_units':['none']},
                                      'action':                       {'unit_name':'joule second',
                                                                       'base_units':'m2 kg s-1', 
                                                                       'standard_units':'m2 kg s-1', 
                                                                       'equiv_units':['none']},
                                      'angular_acceleration':         {'unit_name':'radian per second squared', 
                                                                       'base_units':'s-2',
                                                                       'standard_units':'s-2', 
                                                                       'equiv_units':['none']},
                                      'angular_momentum':             {'unit_name':'newton metre second',
                                                                       'base_units':'m2 kg s-1',
                                                                       'standard_units':'m2 kg s-1', 
                                                                       'equiv_units':['none']},
                                      'angular_velocity':             {'unit_name':'radian per second',
                                                                       'base_units':'s-1', 
                                                                       'standard_units':'s-1', 
                                                                       'equiv_units':['none']},
                                      'angle':                        {'unit_name':'radian',
                                                                       'base_units':'m m-1',
                                                                       'standard_units':'rad',
                                                                       'equiv_units':['rad','rads','radian','radians','unitless']},
                                      'area':                         {'unit_name':'square metre', 
                                                                       'base_units':'m2',
                                                                       'standard_units':'m2',
                                                                       'equiv_units':['none']},
                                      'backscattering_coefficient':   {'unit_name':'reciprocal metre',
                                                                       'base_units':'m-1',
                                                                       'standard_units':'m-1',
                                                                       'equiv_units':['none']},
                                      'catalytic_activity':           {'unit_name':'katal',
                                                                       'base_units':'mol s-1',
                                                                       'standard_units':'kat', 
                                                                       'equiv_units':['kat','katal','katals']},
                                      'column_number_density':        {'unit_name':'number per square metre',
                                                                       'base_units':'m-2',
                                                                       'standard_units':'m-2',
                                                                       'equiv_units':['none']},
                                      'electrical_capacitance':       {'unit_name':'farad', 
                                                                       'base_units':'s4 A2 kg-1 m-2',
                                                                       'standard_units':'F',
                                                                       'equiv_units':['F','farad','farads']},
                                      'electric_charge':              {'unit_name':'coulomb',
                                                                       'base_units':'s A',
                                                                       'standard_units':'C',
                                                                       'equiv_units':['C','coulomb','coulombs']},
                                      'electrical_conductance':       {'unit_name':'siemens',
                                                                       'base_units':'s3 A2 kg-1 m-2',
                                                                       'standard_units':'S',
                                                                       'equiv_units':['S','siemens']},
                                      'electrical_inductance':        {'unit_name':'henry',
                                                                       'base_units':'kg m2 s-2 A-2',
                                                                       'standard_units':'H',
                                                                       'equiv_units':['H','henry','henries']},
                                      'electrical_resistance':        {'unit_name':'ohm', 
                                                                       'base_units':'kg m2 s-3 A-2',
                                                                       'standard_units':'Ω',
                                                                       'equiv_units':['Ω','ohm','ohms']},
                                      'energy':                       {'unit_name':'joule', 
                                                                       'base_units':'kg m2 s-2',
                                                                       'standard_units':'J', 
                                                                       'equiv_units':['J','joule','joules']},
                                      'equivalent_dose':              {'unit_name':'sievert', 
                                                                       'base_units':'m2 s-2', 
                                                                       'standard_units':'Sv', 
                                                                       'equiv_units':['Sv','sievert','sieverts']},
                                      'frequency':                    {'unit_name':'hertz', 
                                                                       'base_units':'s-1', 
                                                                       'standard_units':'Hz', 
                                                                       'equiv_units':['Hz','hertz']},
                                      'force':                        {'unit_name':'newton', 
                                                                       'base_units':'kg m s-2', 
                                                                       'standard_units':'N', 
                                                                       'equiv_units':['N','newton','newtons']},
                                      'heat_capacity':                {'unit_name':'joule per kelvin', 
                                                                       'base_units':'m2 kg s-2 K-1', 
                                                                       'standard_units':'m2 kg s-2 K-1',
                                                                       'equiv_units':['none']},
                                      'illuminance':                  {'unit_name':'lux',
                                                                       'base_units':'cd m-2', 
                                                                       'standard_units':'lx', 
                                                                       'equiv_units':['lx','lux','luxes']},
                                      'jerk':                         {'unit_name':'metre per second cubed',
                                                                       'base_units':'m s-3',
                                                                       'standard_units':'m s-3', 
                                                                       'equiv_units':['none']},
                                      'luminance':                    {'unit_name':'candela per square metre',
                                                                       'base_units':'cd m-2', 
                                                                       'standard_units':'cd m-2', 
                                                                       'equiv_units':['none']},
                                      'luminous_flux':                {'unit_name':'lumen', 
                                                                       'base_units':'cd', 
                                                                       'standard_units':'lm',
                                                                       'equiv_units':['lm','lumen','lumens']},
                                      'magnetic_flux':                {'unit_name':'weber',
                                                                       'base_units':'kg m2 s-2 A-1',
                                                                       'standard_units':'Wb',
                                                                       'equiv_units':['Wb','weber','webers']},
                                      'magnetic_field_strength':      {'unit_name':'tesla',
                                                                       'base_units':'kg s-2 A-1', 
                                                                       'standard_units':'T', 
                                                                       'equiv_units':['T','tesla']},
                                      'mass_absorption_cross_section':{'unit_name':'area per mass',
                                                                       'base_units':'m2 kg-1',
                                                                       'standard_units':'m2 kg-1',
                                                                       'equiv_units':['none']},
                                      'mass_area_density':            {'unit_name':'kilogram per square metre',
                                                                       'base_units':'kg m-2',
                                                                       'standard_units':'kg m-2',
                                                                       'equiv_units':['none']},
                                      'mass_concentration':           {'unit_name':'kilogram per cubic metre',
                                                                       'base_units':'kg m-3', 
                                                                       'standard_units':'kg m-3',
                                                                       'equiv_units':['none']},
                                      'mass_flux':                    {'unit_name':'kilogram per square metre per second',
                                                                       'base_units':'kg m-2 s-1',
                                                                       'standard_units':'kg m-2 s-1',
                                                                       'equiv_units':['none']},                                  
                                      'mass_fraction':                {'unit_name':'kilogram per kilogram',
                                                                       'base_units':'kg kg-1', 
                                                                       'standard_units':'kg kg-1',
                                                                       'equiv_units':['unitless']},
                                      'molarity':                     {'unit_name':'mole per cubic metre', 
                                                                       'base_units':'mol m-3', 
                                                                       'standard_units':'M', 
                                                                       'equiv_units':['Molar','molar','M']},
                                      'molar_heat_capacity':          {'unit_name':'joule per kelvin mole',
                                                                       'base_units':'m2 kg s-2 K-1 mol-1', 
                                                                       'standard_units':'m2 kg s-2 K-1 mol-1', 
                                                                       'equiv_units':['none']},
                                      'molar_mass':                   {'unit_name':'kilogram per mole',
                                                                       'base_units':'kg mol-1', 
                                                                       'standard_units':'kg mol-1', 
                                                                       'equiv_units':['none']},
                                      'molar_volume':                 {'unit_name':'cubic metre per mole', 
                                                                       'base_units':'m3 mol-1', 
                                                                       'standard_units':'m3 mol-1', 
                                                                       'equiv_units':['none']},
                                      'momentum':                     {'unit_name':'newton second', 
                                                                       'base_units':'m kg s-1', 
                                                                       'standard_units':'m kg s-1',  
                                                                       'equiv_units':['none']},
                                      'mole_fraction':                {'unit_name':'mole per mole', 
                                                                       'base_units':'mol mol-1', 
                                                                       'standard_units':'mol mol-1', 
                                                                       'equiv_units':['unitless']},
                                      'number_density':               {'unit_name':'number per cubic metre',
                                                                       'base_units':'m-3',
                                                                       'standard_units':'m-3',
                                                                       'equiv_units':['none']},
                                      'optical_depth':                {'unit_name':'extinction coefficient per path length', 
                                                                       'base_units':'m m-1', 
                                                                       'standard_units':'m m-1', 
                                                                       'equiv_units':['unitless']},
                                      'pH':                           {'unit_name':'mole per cubic metre',
                                                                       'base_units':'mol m-3 mol-1 m3',
                                                                       'standard_units':'mol m-3 mol-1 m3',
                                                                       'equiv_units':['unitless']},
                                      'population_density':           {'unit_name':'population per square metre',
                                                                       'base_units':'m-2',
                                                                       'standard_units':'m-2',
                                                                       'equiv_units':['none']},
                                      'power':                        {'unit_name':'watt', 
                                                                       'base_units':'kg m2 s-3',
                                                                       'standard_units':'W',
                                                                       'equiv_units':['W','watt','watts']},
                                      'pressure':                     {'unit_name':'pascal',
                                                                       'base_units':'kg m-1 s−2',
                                                                       'standard_units':'Pa', 
                                                                       'equiv_units':['Pa','pascal','pascals']},
                                      'radioactivity':                {'unit_name':'becquerel',
                                                                       'base_units':'s-1', 
                                                                       'standard_units':'Bq',
                                                                       'equiv_units':['Bq','becquerel','becquerels']},
                                      'refractive_index':             {'unit_name':'metre per second per metre per second',
                                                                       'base_units':'m s-1 m-1 s',
                                                                       'standard_units':'m s-1 m-1 s',
                                                                       'equiv_units':['unitless']},
                                      'relative_humidity':            {'unit_name':'pascal per pascal',
                                                                       'base_units':'kg m-1 s-2 kg-1 m s2',
                                                                       'standard_units':'kg m-1 s-2 kg-1 m s2',
                                                                       'equiv_units':['unitless']},
                                      'snap':                         {'unit_name':'metre per second to the fourth', 
                                                                       'base_units':'m s-4',
                                                                       'standard_units':'m s-4', 
                                                                       'equiv_units':['none']},
                                      'scattering_coefficient':       {'unit_name':'reciprocal metre',
                                                                       'base_units':'m-1',
                                                                       'standard_units':'m-1',
                                                                       'equiv_units':['none']},
                                      'single_scattering_albedo':     {'unit_name':'scattering coefficient per extinction coefficient', 
                                                                       'base_units':'m m-1',
                                                                       'standard_units':'m m-1',
                                                                       'equiv_units':['unitless']},
                                      'solid_angle':                  {'unit_name':'steradian', 
                                                                       'base_units':'m2 m-2',
                                                                       'standard_units':'sr',
                                                                       'equiv_units':['sr','steradian','steradians','unitless']},
                                      'specific_heat_capacity':       {'unit_name':'joule per kilogram kelvin', 
                                                                       'base_units':'m2 s-2 K-1',
                                                                       'standard_units':'m2 s-2 K-1',
                                                                       'equiv_units':['none']},
                                      'specific_volume':              {'unit_name':'cubic metre per kilogram', 
                                                                       'base_units':'m3 kg-1', 
                                                                       'standard_units':'m3 kg-1', 
                                                                       'equiv_units':['none']},
                                      'spectral_irradiance':          {'unit_name':'watt per cubic metre', 
                                                                       'base_units':'kg m-1 s-3', 
                                                                       'standard_units':'kg m-1 s-3', 
                                                                       'equiv_units':['none']},
                                      'temperature_gradient':         {'unit_name':'kelvin per metre', 
                                                                       'base_units':'K m-1', 
                                                                       'standard_units':'K m-1', 
                                                                       'equiv_units':['none']},
                                      'thermal_resistance':           {'unit_name':'kelvin per watt', 
                                                                       'base_units':'K s3 m-2 kg-1',
                                                                       'standard_units':'K s3 m-2 kg-1', 
                                                                       'equiv_units':['none']},
                                      'thermal_expansion_coefficient':{'unit_name':'reciprocal kelvin',
                                                                       'base_units':'K-1', 
                                                                       'standard_units':'K-1',
                                                                       'equiv_units':['none']},
                                      'torque':                       {'unit_name':'newton metre',
                                                                       'base_units':'m2 kg s-2', 
                                                                       'standard_units':'m2 kg s-2',
                                                                       'equiv_units':['none']},
                                      'velocity':                     {'unit_name':'metre per second',
                                                                       'base_units':'m s-1', 
                                                                       'standard_units':'m s-1', 
                                                                       'equiv_units':['none']},
                                      'voltage':                      {'unit_name':'volt',
                                                                       'base_units':'kg m2 s-3 A-1', 
                                                                       'standard_units':'V', 
                                                                       'equiv_units':['V','volt','volts']},
                                      'volume':                       {'unit_name':'cubic metre',
                                                                       'base_units':'m3',
                                                                       'standard_units':'m3', 
                                                                       'equiv_units':['none']},
                                      'volume_concentration':         {'unit_name':'cubic metre per square metre',
                                                                       'base_units':'m3 m-2',
                                                                       'standard_units':'m3 m-2',
                                                                       'equiv_units':['none']},
                                      'volume_scattering_function':   {'unit_name':'reciprocal metre per steradian',
                                                                       'base_units':'m-1 m2 m-2',
                                                                       'standard_units':'m-1 m2 m-2', 
                                                                       'equiv_units':['none']},
                                      'volumetric_flow':              {'unit_name':'cubic metre per second', 
                                                                       'base_units':'m3 s-1', 
                                                                       'standard_units':'m3 s-1', 
                                                                       'equiv_units':['none']},
                                      'wavenumber':                   {'unit_name':'reciprocal metre', 
                                                                       'base_units':'m-1',
                                                                       'standard_units':'m-1',
                                                                       'equiv_units':['none']},
                                      'yank':                         {'unit_name':'newton per second',
                                                                       'base_units':'m kg s-3',
                                                                       'standard_units':'m kg s-3', 
                                                                       'equiv_units':['none']}}

        # important constants
        self.constants = {
                          'avogadro_constant' : {'value':6.022140857e23, 
                                                 'base_units':'mol-1',               
                                                 'equiv_units':'L'},
                          'boltzmann_constant': {'value':1.38064852e-23, 
                                                 'base_units':'m2 kg s-2 K-1',       
                                                 'equiv_units':'k'},
                          'molar_gas_constant': {'value':8.3144598,      
                                                 'base_units':'kg m2 s-2 K-1 mol-1', 
                                                 'equiv_units':'R'},
                          'molar_mass_dry_air': {'value':0.02896470,
                                                 'base_units':'kg mol-1',
                                                 'equiv_units':'none'},
                          'pi':                 {'value':3.141592653589793,
                                                 'base_units':'none',
                                                 'equiv_units':'unitless'},
                          'planck_constant'   : {'value':6.62606979e-34, 
                                                 'base_units':'m2 kg s-1',           
                                                 'equiv_units':'h'}}

        # map non-SI units to SI quantities with conversion factors
        # can be more 1 type of non-SI units for each SI quantity
        # factor is the number neccessary to multiply by to get to base units
        self.non_SI_quantities = {
                                  'amount of substance':{'equivalent':              {'factor':1./self.valence,
                                                                                     'standard_units':'Eq',
                                                                                     'units':['eq','Eq','equiv','Equiv',
                                                                                              'equivalent','Equivalent',
                                                                                              'equivalents','Equivalents'],
                                                                                     'dependencies':['valence']},
                                                         'mole of carbon':          {'factor':1.0/get_N_atoms(self.species,'C'), 
                                                                                     'standard_units':'molC',
                                                                                     'units':['molC','moleC'],
                                                                                     'dependencies':['species']},                                                            
                                                         'mole of chlorine':        {'factor':1.0/get_N_atoms(self.species,'Cl'), 
                                                                                     'standard_units':'molCl',          
                                                                                     'units':['molCl','moleCl'],
                                                                                     'dependencies':['species']},
                                                         'mole of hydrogen':        {'factor':1.0/get_N_atoms(self.species,'H'), 
                                                                                     'standard_units':'molH',          
                                                                                     'units':['molH','moleH'],
                                                                                     'dependencies':['species']},                                                            
                                                         'mole of nitrogen':        {'factor':1.0/get_N_atoms(self.species,'N'), 
                                                                                     'standard_units':'molN',       
                                                                                     'units':['molN','moleN'],
                                                                                     'dependencies':['species']},                                                                              
                                                         'mole of sulphur':         {'factor':1.0/get_N_atoms(self.species,'S'),    
                                                                                     'standard_units':'molS',        
                                                                                     'units':['molS','moleS'],
                                                                                     'dependencies':['species']}},
                                  'angle':              {'degree':                  {'factor':math.pi/180.0, 
                                                                                     'standard_units':'°',
                                                                                     'units':['deg','degrees','°','angular degrees']}},
                                  'area':               {'acre':                    {'factor':4046.8564224,
                                                                                     'standard_units':'acre',
                                                                                     'units':['acre','acres']}},
                                  'energy':             {'electronvolt':            {'factor':1.6021766208e-19,
                                                                                     'standard_units':'eV',
                                                                                     'units':['eV','electronvolt','electronvolts']}},
                                  'length':             {'mile':                    {'factor':1609.344,
                                                                                     'standard_units':'miles',
                                                                                     'units':['mi', 'miles', 'mile']},                                                       
                                                         'feet':                    {'factor':0.3048,                                     
                                                                                     'standard_units':'ft',             
                                                                                     'units':['feet','ft']},                                                               
                                                         'inch':                    {'factor':0.0254,                                     
                                                                                     'standard_units':'inch',       
                                                                                     'units':['inch','inches']},                                                                             
                                                         'yard':                    {'factor':0.9144,                                       
                                                                                     'standard_units':'yd',          
                                                                                     'units':['yd','yds','yards','yard']},                                              
                                                         'astronomical unit':       {'factor':149597870700,                                   
                                                                                     'standard_units':'au',           
                                                                                     'units':['au']}},
                                  'mass':               {'gram':                    {'factor':1e-3,
                                                                                     'standard_units':'g', 
                                                                                     'units':['g', 'gram', 'grams']},                                                        
                                                         'metric ton':              {'factor':1e3,                                         
                                                                                     'standard_units':'ton',            
                                                                                     'units':['ton','tons','tonne','tonnes', 't']},                                             
                                                         'pound':                   {'factor':0.45359237,                                 
                                                                                     'standard_units':'lb',         
                                                                                     'units':['pounds','lb','lbs']},                                                                         
                                                         'ounce':                   {'factor':0.02834952,                                    
                                                                                     'standard_units':'oz',          
                                                                                     'units':['oz','ozs','ounce','ounces']},                                            
                                                         'unified atomic mass unit':{'factor':1.660539040e-27,                                
                                                                                     'standard_units':'u',            
                                                                                     'units':['u','Da']},                                                                      
                                                         'gram of carbon':          {'factor':1e-3*(get_molecular_mass(self.species)/(get_molecular_mass('C')*get_N_atoms(self.species,'C'))), 
                                                                                     'standard_units':'gC',          
                                                                                     'units':['gC','gramC','gramsC'],
                                                                                     'dependencies':['species']},                                                      
                                                         'gram of chlorine':        {'factor':1e-3*(get_molecular_mass(self.species)/(get_molecular_mass('Cl')*get_N_atoms(self.species,'Cl'))), 
                                                                                     'standard_units':'gCl',         
                                                                                     'units':['gCl','gramCl','gramsCl'],
                                                                                     'dependencies':['species']},                                                                                                 
                                                         'gram of nitrogen':        {'factor':1e-3*(get_molecular_mass(self.species)/(get_molecular_mass('N')*get_N_atoms(self.species,'N'))), 
                                                                                     'standard_units':'gN',          
                                                                                     'units':['gN','gramN','gramsN'],
                                                                                     'dependencies':['species']},
                                                         'gram of calcium carbonate':{'factor':1e-3*((get_molecular_mass(self.species)/self.valence)/(get_molecular_mass('CaCO3')/2)), # valence of calcium carbonate is 2
                                                                                     'standard_units':'gCaCO3',          
                                                                                     'units':['gCaCO3','gramCaCO3','gramsCaCO3'],
                                                                                     'dependencies':['species','valence']},                                                     
                                                         'gram of sulphur':         {'factor':1e-3*(get_molecular_mass(self.species)/(get_molecular_mass('S')*get_N_atoms(self.species,'S'))), 
                                                                                     'standard_units':'gS',           
                                                                                     'units':['gS','gramS','gramsS'],
                                                                                     'dependencies':['species']}},       
                                  # note - volumetric fractions and the mole fractions of the components of an ideal gas mixture are interchangeable --> does not hold for liquids
                                  'mole_fraction':      {'ppmv':                    {'factor':1e-6, 
                                                                                     'standard_units':'umol mol-1', 
                                                                                     'units':['ppmv','ppm','Partspermillion',
                                                                                              'partspermillion', 'ppmV']},                   
                                                         'ppbv':                    {'factor':1e-9,
                                                                                     'standard_units':'nmol mol-1',
                                                                                     'units':['ppbv','ppb','Partsperbillion',
                                                                                              'partsperbillion','ppbV']}, 
                                                         'pptv':                    {'factor':1e-12, 
                                                                                     'standard_units':'pmol mol-1', 
                                                                                     'units':['pptv','ppt','Partspertrillion',
                                                                                              'partspertrillion','pptV']}, 
                                                         'ppmv of carbon':          {'factor':1e-6/(get_N_atoms(self.species,'C')), 
                                                                                     'standard_units':'umolC mol-1', 
                                                                                     'units':['ppmvC','ppmC','partspermillioncarbon',
                                                                                              'PartspermillionCarbon','ppmVC'],
                                                                                     'dependencies':['species']}, 
                                                         'ppmv of chlorine':        {'factor':1e-6/(get_N_atoms(self.species,'Cl')), 
                                                                                     'standard_units':'umolCl mol-1', 
                                                                                     'units':['ppmvCl','ppmCl',
                                                                                              'partspermillionchlorine',
                                                                                              'PartspermillionChlorine',
                                                                                              'ppmVCl'],
                                                                                     'dependencies':['species']}, 
                                                         'ppmv of nitrogen':        {'factor':1e-6/(get_N_atoms(self.species,'N')), 
                                                                                     'standard_units':'umolN mol-1', 
                                                                                     'units':['ppmvN','ppmN',
                                                                                              'partspermillionnitrogen',
                                                                                              'PartspermillionNitrogen',
                                                                                              'ppmVN'],
                                                                                     'dependencies':['species']}, 
                                                         'ppmv of sulphur':         {'factor':1e-6/(get_N_atoms(self.species,'S')), 
                                                                                     'standard_units':'umolS mol-1', 
                                                                                     'units':['ppmvS','ppmS',
                                                                                              'partspermillionsulphur',
                                                                                              'PartspermillionSulphur',
                                                                                              'partspermillionsulfur',
                                                                                              'PartspermillionSulfur',
                                                                                              'ppmVS'],
                                                                                     'dependencies':['species']}, 
                                                         'ppbv of carbon':          {'factor':1e-9/(get_N_atoms(self.species,'C')),       
                                                                                     'standard_units':'nmolC mol-1', 
                                                                                     'units':['ppbvC','ppbC',
                                                                                              'partsperbillioncarbon',
                                                                                              'PartsperbillionCarbon',
                                                                                              'ppbVC'],
                                                                                     'dependencies':['species']}, 
                                                         'ppbv of chlorine':        {'factor':1e-9/(get_N_atoms(self.species,'Cl')),  
                                                                                     'standard_units':'nmolCl mol-1', 
                                                                                     'units':['ppbvCl','ppbCl',
                                                                                              'partsperbillionchlorine',
                                                                                              'PartsperbillionChlorine',
                                                                                              'ppbVCl'],
                                                                                     'dependencies':['species']}, 
                                                         'ppbv of nitrogen':        {'factor':1e-9/(get_N_atoms(self.species,'N')), 
                                                                                     'standard_units':'nmolN mol-1', 
                                                                                     'units':['ppbvN','ppbN',
                                                                                              'partsperbillionnitrogen',
                                                                                              'PartsperbillionNitrogen',
                                                                                              'ppbVN'], 
                                                                                     'dependencies':['species']}, 
                                                         'ppbv of sulphur':         {'factor':1e-9/(get_N_atoms(self.species,'S')), 
                                                                                     'standard_units':'nmolS mol-1', 
                                                                                     'units':['ppbvS','ppbS',
                                                                                              'partsperbillionsulphur',
                                                                                              'PartsperbillionSulphur',
                                                                                              'partsperbillionsulfur',
                                                                                              'PartsperbillionSulfur',
                                                                                              'ppbVS'],
                                                                                     'dependencies':['species']}, 
                                                         'pptv of carbon':          {'factor':1e-12/(get_N_atoms(self.species,'C')), 
                                                                                     'standard_units':'pmolC mol-1', 
                                                                                     'units':['pptvC','pptC',
                                                                                              'partspertrillioncarbon',
                                                                                              'PartspertrillionCarbon',
                                                                                              'pptVC'],
                                                                                     'dependencies':['species']}, 
                                                         'pptv of chlorine':        {'factor':1e-12/(get_N_atoms(self.species,'Cl')), 
                                                                                     'standard_units':'pmolCl mol-1', 
                                                                                     'units':['pptvCl','pptCl',
                                                                                              'partspertrillionchlorine',
                                                                                              'PartspertrillionChlorine',
                                                                                              'pptVCl'],
                                                                                     'dependencies':['species']}, 
                                                         'pptv of nitrogen':        {'factor':1e-12/(get_N_atoms(self.species,'N')), 
                                                                                     'standard_units':'pmolN mol-1', 
                                                                                     'units':['pptvN','pptN',
                                                                                              'partspertrillionnitrogen',
                                                                                              'PartspertrillionNitrogen',
                                                                                              'pptVN'],
                                                                                     'dependencies':['species']}, 
                                                         'pptv of sulphur':         {'factor':1e-12/(get_N_atoms(self.species,'S')), 
                                                                                     'standard_units':'pmolS mol-1', 
                                                                                     'units':['pptvS','pptS',
                                                                                              'partspertrillionsulphur',
                                                                                              'PartspertrillionSulphur',
                                                                                              'partspertrillionsulfur',
                                                                                              'PartspertrillionSulfur',
                                                                                              'pptVS'],
                                                                                     'dependencies':['species']}},
                                  'pressure':           {'atmosphere':              {'factor':101325, 
                                                                                     'standard_units':'atm', 
                                                                                     'units':['atm', 'atmospheres',
                                                                                              'atmosphere']},                                         
                                                         'bar':                     {'factor':1e5, 
                                                                                     'standard_units':'bar',
                                                                                     'units':['bar']},                                                                     
                                                         'torr':                    {'factor':133.3224, 
                                                                                     'standard_units':'Torr',
                                                                                     'units':['Torr','torr']}, 
                                                         'psi':                     {'factor':6894.76, 
                                                                                     'standard_units':'psi',
                                                                                     'units':['psi']}},  
                                  'temperature':        {'celsius':                 {'factor':1.0, 
                                                                                     'standard_units':'°C',
                                                                                     'units':['degC','celsius','°C',
                                                                                              '°Celsius','degreesC',
                                                                                              'degCelsius',
                                                                                              'degreesCelsius']}, 
                                                         'rankine':                 {'factor':5./9., 
                                                                                     'standard_units':'°R', 
                                                                                     'units':['degR','rankine','°R',
                                                                                              '°Rankine','degreesR',
                                                                                              'degRankine',
                                                                                              'degreesRankine']}, 
                                                         'fahrenheit':              {'factor':5./9.,
                                                                                     'standard_units':'°F',
                                                                                     'units':['degF','fahrenheit', '°F', 
                                                                                              '°Fahrenheit', 'degreesF', 
                                                                                              'degFahrenheit', 
                                                                                              'degreesFahrenheit']}},
                                  'time':               {'minute':                  {'factor':60.0,
                                                                                     'standard_units':'mins', 
                                                                                     'units':['min','mins','minutes',
                                                                                              'minute']}, 
                                                         'hour':                    {'factor':3600.0,
                                                                                     'standard_units':'hours',
                                                                                     'units':['hr','hrs','hours',
                                                                                              'hour']},
                                                         'day':                     {'factor':86400.0, 
                                                                                     'standard_units':'days', 
                                                                                     'units':['day','days']},
                                                         'year':                    {'factor':31536000.0, 
                                                                                     'standard_units':'years', 
                                                                                     'units':['year','years', 'yr']}},
                                  'volume':             {'litre':                   {'factor':1e-3,
                                                                                     'standard_units':'l', 
                                                                                     'units':['l','litre','litres',
                                                                                              'liters','L']}, 
                                                         'imperial pint':           {'factor':0.000568261, 
                                                                                     'standard_units':'imperial_pints', 
                                                                                     'units':['pt','pint','pints']},
                                                         'imperial gallon':         {'factor':0.00454609, 
                                                                                     'standard_units':'gal',
                                                                                     'units':['gal','gallon','gals',
                                                                                              'gallons']}}}

        # quantity conversion formulae
        self.conversion_formulae = {
                                    'mole_fraction':                                                             ["mass_concentration*((molar_gas_constant*temperature)/(molar_mass*pressure))",        
                                                                                                                  "number_density*((molar_gas_constant*temperature)/(avogadro_constant*pressure))",
                                                                                                                  "mass_fraction*(molar_mass_dry_air/molar_mass)"],
                                    'mass_fraction':                                                             ["mole_fraction*(molar_mass/molar_mass_dry_air)"],
                                    'mass_concentration':                                                        ["mole_fraction*((molar_mass*pressure)/(molar_gas_constant*temperature))",
                                                                                                                  "(number_density*molar_mass)/avogadro_constant",
                                                                                                                  "molarity*molar_mass"],
                                    'molarity':                                                                  ["mass_concentration/molar_mass"],
                                    'number_density':                                                            ["mole_fraction*((avogadro_constant*pressure)/(molar_gas_constant*temperature))",
                                                                                                                  "(mass_concentration*avogadro_constant)/molar_mass"],
                                    'backscattering_coefficient':                                                ["volume_scattering_function*(2.*pi)"],
                                    'scattering_coefficient':                                                    ["volume_scattering_function*(4.*pi)"],
                                    'volume_scattering_function':                                                ["backscattering_coefficient/(2.*pi)",
                                                                                                                  "scattering_coefficient/(4.*pi)"],
                                    'mass_concentration>mole_fraction--conversion_factor':                       ["((molar_gas_constant*temperature)/(molar_mass*pressure))"],
                                    'mole_fraction>mass_concentration--conversion_factor':                       ["((molar_mass*pressure)/(molar_gas_constant*temperature))"],
                                    'molarity>mass_concentration--conversion_factor':                            ["molar_mass"],
                                    'mass_concentration>molarity--conversion_factor':                            ["1./molar_mass"],
                                    'mass_fraction>mole_fraction--conversion_factor':                            ["molar_mass_dry_air/molar_mass"],
                                    'mole_fraction>mass_fraction--conversion_factor':                            ["molar_mass/molar_mass_dry_air"],
                                    'volume_scattering_function>backscattering_coefficient--conversion_factor':  ["2.*pi"],
                                    'volume_scattering_function>scattering_coefficient--conversion_factor':      ["4.*pi"], 
                                    'backscattering_coefficient>volume_scattering_function--conversion_factor':  ["1./(2.*pi)"],
                                    'scattering_coefficient>volume_scattering_function--conversion_factor':      ["1./(4.*pi)"]
        }

def get_molecular_mass(molecule):
    """
    Calculate the molecular mass of a chemical compound in kg mol-1.

    Parameters
    ----------
    molecule : str
        The chemical formula of the molecule to be processed.

    Returns
    -------
    float
        The total molecular mass in kilograms per mole, or NaN if the input is empty.
    """

    # if molecule is empty string then return molecular mass as np.NaN
    if molecule == '':
        return np.nan

    # define atomic masses per element in g mol-1
    atomic_masses = {
                     "H":1.00794,"He":4.002602,"Li":6.941,"Be":9.012182,"B":10.811,"C":12.011,"N":14.00674,"O":15.9994,
                     "F":18.9984032,"Ne":20.1797,"Na":22.989768,"Mg":24.3050,"Al":26.981539,"Si":28.0855,"P":30.973762,
                     "S":32.066,"Cl":35.4527,"Ar":39.948,"K":39.0983,"Ca":40.078,"Sc":44.955910,"Ti":47.88,"V":50.9415,
                     "Cr":51.9961,"Mn":54.93805,"Fe":55.847,"Co":58.93320,"Ni":58.6934,"Cu":63.546,"Zn":65.39,
                     "Ga":69.723,"Ge":72.61,"As":74.92159,"Se":78.96,"Br":79.904,"Kr":83.80,"Rb":85.4678,"Sr":87.62,
                     "Y":88.90585,"Zr":91.224,"Nb":92.90638,"Mo":95.94,"Tc":98,"Ru":101.07,"Rh":102.90550,"Pd":106.42,
                     "Ag":107.8682,"Cd":112.411,"In":114.82,"Sn":118.710,"Sb":121.757,"Te":127.60,"I":126.90447,
                     "Xe":131.29,"Cs":132.90543,"Ba":137.327,"La":138.9055,"Ce":140.115,"Pr":140.90765,"Nd":144.24,
                     "Pm":145,"Sm":150.36,"Eu":151.965,"Gd":157.25,"Tb":158.92534,"Dy":162.50,"Ho":164.93032,
                     "Er":167.26,"Tm":168.93421,"Yb":173.04,"Lu":174.967,"Hf":178.49,"Ta":180.9479,"W":183.85,
                     "Re":186.207,"Os":190.2,"Ir":192.22,"Pt":195.08,"Au":196.96654,"Hg":200.59,"Tl":204.3833,
                     "Pb":207.2,"Bi":208.98037,"Po":209,"At":210,"Rn":222,"Fr":223,"Ra":226.0254,"Ac":227,"Th":232.0381,
                     "Pa":213.0359,"U":238.0289,"Np":237.0482,"Pu":244,"Am":243,"Cm":247,"Bk":247,"Cf":251,"Es":252,
                     "Fm":257,"Md":258,"No":259,"Lr":260,"Rf":261,"Db":262,"Sg":263,"Bh":262,"Hs":265,"Mt":266}


    s = re.findall('([A-Z][a-z]?)([0-9]*)', molecule)
    compoundweight = 0

    # take sum of all element atomic masses (in g mol-1) 
    for element, count in s:
        count = int(count or '1')
        compoundweight += atomic_masses[element] * count

    # convert to kg mol-1
    compoundweight = compoundweight/1e3

    return compoundweight

def get_N_atoms(molecule, element):
    """
    Determine the integer count of a specific element within a chemical formula.

    Parameters
    ----------
    molecule : str
        The chemical formula of the molecule to be analysed.
    element : str
        The chemical symbol of the element to count.

    Returns
    -------
    int or float
        The number of atoms of the specified element, or NaN if the input is invalid or missing.
    """

    # if molecule or element is empty string then return number of atoms as np.NaN
    if (molecule == '') or (element == ''):
        return np.nan

    parse_list = re.findall(r'[A-Z][a-z]*|\d+', re.sub(r'[A-Z][a-z]*(?![\da-z])', r'\g<0>1', molecule))
    # if element not in molecule, return np.nan 
    try:
        n_atoms_element = parse_list[parse_list.index(element)+1]
    except:
        return np.nan   

    return int(n_atoms_element)

def get_valence(molecule):
    """
    Retrieve the chemical valence or absolute charge for a specific molecule.

    Parameters
    ----------
    molecule : str
        The chemical formula or ion symbol to be looked up.

    Returns
    -------
    int or float
        The valence value associated with the molecule, or NaN if not found.
    """

    # if molecule is empty string then return valence as np.NaN
    if molecule == '':
        return np.nan

    valence = {
        'Al':        3,  # total dissolved Al, mainly Al³⁺ and hydrolyzed forms
        'As':        3,  # total As, dominated by As³⁺ (arsenite) and As⁵⁺ (arsenate)
        'C20H12':    1,  # neutral organic compound (benzo[a]pyrene)
        'Ca++':      2,  # explicitly Ca²⁺ ion
        'CaCO3':     2,  # calcium carbonate, CaCO3
        'Cd':        2,  # Cd²⁺ is the stable aqueous form
        'Cl-':       1,  # chloride anion Cl⁻
        'Co':        2,  # Co²⁺ is the common dissolved oxidation state
        'Cr':        3,  # total dissolved Cr, mostly Cr³⁺ (some Cr⁶⁺ possible)
        'Cu':        2,  # Cu²⁺ dominant in water; Cu⁺ rarely stable
        'F-':        1,  # fluoride anion F⁻
        'Fe':        3,  # total Fe, mostly Fe³⁺ (Fe²⁺ minor)
        'H+':        1,  # hydrogen ion H⁺
        'Hg':        2,  # Hg²⁺ is dominant soluble form
        'K+':        1,  # potassium cation K⁺
        'Mg++':      2,  # magnesium cation Mg²⁺
        'Mn':        2,  # total Mn, dominated by Mn²⁺
        'CH4O3S':    1,  # methanesulfonic acid (monoprotic acid)
        'Na+':       1,  # sodium cation Na⁺
        'NH4+':      1,  # ammonium cation NH₄⁺
        'NH4NO3':    1,  # neutral salt, effectively 1 for mass balance
        'Ni':        2,  # total Ni, mainly Ni²⁺ in water
        'NO3-':      1,  # nitrate anion NO₃⁻
        'Pb':        2,  # total dissolved Pb, predominantly Pb²⁺
        'PO4---':    3,  # phosphate anion PO₄³⁻
        'Se':        2,  # total Se, dominated by SeO₄²⁻ (Se(VI))
        'SO4--':     2,  # sulfate anion SO₄²⁻
        'V':         5,  # total V, mainly as vanadate V(V)
        'Zn':        2,  # Zn²⁺ is stable aqueous form
    }

    if molecule not in valence:
        return np.nan
    else:
        return valence[molecule]