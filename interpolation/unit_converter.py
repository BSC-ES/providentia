#WRITTEN BY DENE BOWDALO

###------------------------------------------------------------------------------------###
###IMPORT MODULES
###------------------------------------------------------------------------------------###

from collections import Counter
import numpy as np
import re
import sys

###------------------------------------------------------------------------------------###

class convert_units(object):

    '''class that converts units 

       accepts mixture of SI base/derived/equivalent units and non-SI units as input/output arguments
       can handle complex formatted input/output units (i.e. without spaces, spaced, including '/' or 'per' terms)

       converts units between same quantities (i.e. length --> length, cm --> m)
       and converts between different quantities (i.e. mass density --> mole fraction, mg/m3 --> ppbv) **when given the correct parameters needed/ and have formula

       **units which are solely per something e.g. frequency (s-1) are prefixed with 'xx' to prevent confusion between units
       i.e. ms -1 (per miliseconds) could be taken to be m s-1 (metres per second) or vice versa in processing  

       requires input of: 
       1. input units -- either single string for same quantity conversion, or a dictionary of input quantity name and units for converting between quantities {'quantity':unit1, 'quantity2':unit2}  
       **quantity names should be given with underscore joins when more 1 word 

       2. output units

       3. input value/s -- either a single value for same quantity conversion, or a dictionary of input quantity name and values for converting between quantities {'quantity':value, 'quantity2':value2} 
       **quantity names should be given with underscore joins when more 1 word 

       4. Level of precision of converted output/conversion factor (fixed number of decimal places) **OPTIONAL 

       5. Measured species (in chemical notation). **OPTIONAL -- used when needing to convert units associated with an element (e.g. ppmv N) 

       6. Conversion input quantity **OPTIONAL -- used when wanting to return conversion factor for conversions with more than 1 input variables, i.e. this quantity refers to the input variable for which a conversion factor to the output units is wanted

       ----------------------------------------------------------------------------------------
       
       How to use class

       -----------------------------
       Same quantity unit conversion (m3 --> dm3)

       1. Define input variables 
       input_units = 'm3'
       output units = 'dm3'
       input_value = 1.0

       2. Create object instance
       conv_obj = converter.convert_units(input_units, output_units, input_value)

       3. Do conversion
       conv_obj.do_conversion()

       4. Processed output now stored in class 
       conv_obj.converted_value   --> Converted output value 
       conv_obj.conversion_factor --> Conversion factor to get from input value to output value (*Note this is not valid for temperature parameters)

       -----------------------------
       #Different unit quantity conversion (ug m-3 --> ppbv)

       1. Define input variables 
       input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'kg mol-1', 'mass_density':'ug m-3'}
       output_units = 'ppbv'
       input_values = {'temperature':273.0, 'pressure':1013.0, 'molar_mass':48.0, 'mass_density':10.0}

       2. Create object instance
       conv_obj = converter.convert_units(input_units, output_units, input_values, conversion_input_quantity='mass_density')

       3. Do conversion
       conv_obj.do_conversion()

       4. Processed output now stored in class 
       conv_obj.converted_value   --> Converted output value 
       conv_obj.conversion_factor --> Conversion factor to get from input value to output value (*Note this is not valid for temperature parameters)
    '''

    def __init__(self,input_units,output_units,input_value,precision=-999,measured_species='',conversion_input_quantity=''):
        if isinstance(input_units, dict):
            self.input_units = {k: str(v) for k, v in input_units.items()}
        else:
            self.input_units = str(input_units)
        
        self.output_units = str(output_units)
        
        if isinstance(input_value, dict):
            self.input_value = {k: float(v) for k, v in input_value.items()}
        else:
            self.input_value = float(input_value)
        
        self.precision                 = int(precision)
        self.measured_species          = str(measured_species)
        self.conversion_input_quantity = str(conversion_input_quantity)
        
        #load key dictionaries into self
        self.load_dictionaries()
        
    
    def do_conversion(self):

        '''main conversion module'''
        
        #get prefix_standard names
        standard_prefix_names = []
        for prefix in self.standard_prefixes:
            standard_prefix_names.append(self.standard_prefixes[prefix]['standard_units'])
        
        #handle errors
        try:
        
            #if input unit is just a string, get conversion factor to convert between given input and output units to standard SI
            #the ratio between these conversion factors * the input value then gives the desired output value in the desired units   
    
            if isinstance(self.input_units, str):

                #get input/output units into standardised format and get unified conversion factor to base SI
                self.standardise_unit_form(self.input_units,'input')
                self.standardise_unit_form(self.output_units,'output')

                #checks if units are of same quantity and can be converted between
                self.check_quantity_agreement()

                #get scaling ratio between input and output units
                self.scaling_factor = self.input_scaling_factor_to_reference/self.output_scaling_factor_to_reference
    
                #correct for temperature scaling, only scale for absolute temps. 
                if (self.output_represented_quantity == 'temperature') & (self.input_standard_units[0] not in standard_prefix_names) & (self.output_standard_units[0] not in standard_prefix_names):
                    self.correct_temperature(self.input_value,self.input_standard_units,self.output_standard_units)
        
                #multiply input value by final calculated scaling factor 
                else:                  
                    self.converted_value = self.input_value*self.scaling_factor
                    #save out conversion factor 
                    self.conversion_factor = self.scaling_factor
    
            #--------------------------------------------------------------
            #if input units are in dict, want to convert between quantities.
            #process each of the input units into standard SI units, and then calculate output value using conversion formula

            elif isinstance(self.input_units, dict):
    
                #put input quantity names into list
                self.input_quantities = self.input_units.keys()
                #add mathematical constants to input quantities
                for m_quantity in self.constants:
                    self.input_quantities.append(m_quantity)
    
                #get output units into standardised format and get unified conversion factor to base SI
                self.standardise_unit_form(self.output_units,'output')
    
                #check that have formula that converts between input variables --> output variable
                self.check_conversion_validity()
    
                #list for appending input values, converted to SI
                self.SI_input_quantities = {}
    
                #add mathematical constants to SI_input_quantities
                for m_quantity in self.constants:
                    self.SI_input_quantities[m_quantity] = self.constants[m_quantity]['value']
    
                #iterate through input units, getting into standardised format and get unified conversion factor and get unified conversion factor to base SI
                for quantity in self.input_quantities:
                    if quantity not in self.constants:
                    
                        self.standardise_unit_form(self.input_units[quantity],'input')
    
                        #scaling factor
                        self.scaling_factor = self.input_scaling_factor_to_reference

                        #save scaling factor if want to return conversion factor
                        if self.conversion_input_quantity == quantity:
                            conversion_scaling_factor_1 = self.scaling_factor

                        #correct for temperature scaling to SI Kelvin
                        if (quantity == 'temperature') & (self.input_standard_units[0] not in standard_prefix_names):
                            self.correct_temperature(self.input_value[quantity],self.input_standard_units,'K')
            
                        #convert to base SI
                        else:
                            self.converted_value = self.input_value[quantity]*self.scaling_factor

                        #add SI input parameter name:converted_value to dict
                        self.SI_input_quantities[quantity] = self.converted_value
    
                #if have a conversion input quantity, return conversion factor to get output quantity by: input_quantity*conversion_factor
                if self.conversion_input_quantity != '':
                    formula = self.conversion_formulae['{}>{}--conversion_factor'.format(self.conversion_input_quantity,self.output_represented_quantity)]
                    conversion_scaling_factor_2 = eval(formula,self.SI_input_quantities)
                    self.conversion_factor = conversion_scaling_factor_1 * conversion_scaling_factor_2
            
                #take converted SI quantities and calculate SI output quantity using appropriate formula
                self.converted_value = eval(self.formula,self.SI_input_quantities)
    
                #finally, convert output SI quantity into desired output units (if necessary)
                if self.output_standard_units != self.output_reference_units:
        
                    self.standardise_unit_form(self.output_reference_units,'input')

                    #get scaling ratio between input (SI units) and output units
                    self.scaling_factor = self.input_scaling_factor_to_reference/self.output_scaling_factor_to_reference
        
                    #save scaling factor if want to return conversion factor
                    if self.conversion_input_quantity != '':
                        conversion_scaling_factor_3 = self.scaling_factor
                        self.conversion_factor = self.conversion_factor * conversion_scaling_factor_3
            
                    #correct for temperature scaling from K
                    if (self.output_represented_quantity == 'temperature') & (self.output_standard_units[0] not in standard_prefix_names):
                        self.correct_temperature(self.converted_value,'K',self.output_standard_units)
                    #convert to desired output units
                    else:
                        self.converted_value = self.converted_value*self.scaling_factor
    
            #finish by rounding converted value/conversion factor to fixed number of decimal places (if desired)
            if self.precision != -999:
                self.converted_value = np.around(self.converted_value, self.precision)
                try:
                    self.conversion_factor = np.around(self.conversion_factor, self.precision)
                except:
                    pass
            return
            
        except:
            sys.exit(self.error_message)
        
            return
        
    def check_quantity_agreement(self):

       '''check if input quantity is same as output quantity'''

        #print 'cleaned units are: {} --> {}'.format(self.input_cleaned_units,self.output_cleaned_units)
        #print 'reference units are: {} --> {}'.format(self.input_reference_units,self.output_reference_units)
        #print 'standard units are: {} --> {}'.format(self.input_standard_units,self.output_standard_units)
        
        #get SI quantities that are represented by input/output units
        self.input_represented_quantity = self.get_quantity_representation(self.input_reference_units)
        self.output_represented_quantity = self.get_quantity_representation(self.output_reference_units)
        
        #are quantities represented same between input/output?
        if self.input_represented_quantity != self.output_represented_quantity:
            self.error_message = 'Input and Output Quantities do not agree: {}, {}'.format(self.input_represented_quantity, self.output_represented_quantity)
            return 1+'a'
            
    def check_conversion_validity(self):

        '''check if have formula for converting to output quantities, and correct input quantities for doing so'''
    
        #check if input unit/value dicts have the same length
        if len(self.input_units.keys()) != len(self.input_value.keys()):
            self.error_message = 'Input unit and value dictionaries do not have same length'
            return 1+'a'
        
        #check that input unit/value dicts have same quantity names
        if set(self.input_units.keys()) != set(self.input_value.keys()):
            self.error_message = 'Input unit and value dictionaries do not have same quantity names'
            return 1+'a'  
        
        #get SI quantities represented by output_units
        self.output_represented_quantity = self.get_quantity_representation(self.output_reference_units)
        
        #Firstly, do you have a formula for converting to desired output quantity?
        available_conversion_formulae_quantities = self.conversion_formulae.keys()
    
        #if don't have formula, return error message
        if self.output_represented_quantity not in available_conversion_formulae_quantities:
            self.error_message = 'Do not have formula to calculate desired output quantity'
            return 1+'a'
        else:
            #next, check if have given required input variables for calculating output quantity using one of formulae defined in conversion_formulae dictionary
            quantity_formulae = self.conversion_formulae[self.output_represented_quantity]
            for formula in quantity_formulae:
                formula_quantities = formula.replace('*','$').replace('/','$').replace('(','').replace(')','').split('$')
                #test to see if have all of formula quantities required in input quantities
                have_quantities = set(formula_quantities).issubset(self.input_quantities)
                if have_quantities == True:
                    self.formula = formula
                    return    
            
            #if do not valid input quantities for one of the formula return error message
            print formula_quantities
            self.error_message = 'Do not have correct input quantities to calculate desired output quantity'
            return 1+'a'
             
    def get_quantity_representation(self, units):

        '''get quantity represented by given units'''

        #is highest level quantity represented a base SI unit?
        for quantity in self.SI_base_quantities:
            if Counter(units) == Counter(self.SI_base_quantities[quantity]['base_units']):
                return quantity
                
        #is highest level quantity represented a derived SI unit?
        for quantity in self.SI_derived_quantities:
            if Counter(units) == Counter(self.SI_derived_quantities[quantity]['base_units']):
                return quantity
            
        #is highest level quantity represented a scientific constant?
        for quantity in self.constants:
            if Counter(units) == Counter(self.constants[quantity]['base_units']):
                return quantity
                
        #otherwise undefined derived units, rearrange to put in alphabetical order, preserving spaces
        reference_unit_split = units.split(' ')
        reference_unit_split_sorted = sorted(reference_unit_split)
        return ' '.join(reference_unit_split_sorted)
                                   
    def correct_temperature(self,input_value,input_standard_unit,output_standard_unit):
        
       '''if output quantity is solely a temperature, shift scales if required'''

        #kelvin --> kelvin
        if ('K' in input_standard_unit) & ('K' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
            
        #celsius --> celsius
        if ('°C' in input_standard_unit) & ('°C' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
        
        #fahrenheit --> fahrenheit
        if ('°F' in input_standard_unit) & ('°F' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
        
        #rankine --> rankine
        if ('°R' in input_standard_unit) & ('°R' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
        
        #kelvin --> fahrenheit
        if ('K' in input_standard_unit) & ('°F' in output_standard_unit):
            self.converted_value = (input_value * self.scaling_factor) - 459.67
        
        #kelvin --> rankine
        if ('K' in input_standard_unit) & ('°R' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
        
        #kelvin --> celsius
        if ('K' in input_standard_unit) & ('°C' in output_standard_unit):
            self.converted_value = input_value - 273.15
        
        #celsius --> kelvin
        if ('°C' in input_standard_unit) & ('K' in output_standard_unit):
            self.converted_value = input_value + 273.15
            
        #celsius --> fahrenheit
        if ('°C' in input_standard_unit) & ('°F' in output_standard_unit):
            self.converted_value = (input_value * self.scaling_factor) + 32.0
            
        #celsius --> rankine
        if ('°C' in input_standard_unit) & ('°R' in output_standard_unit):
            self.converted_value = (input_value + 273.15) * self.scaling_factor
        
        #fahrenheit --> kelvin
        if ('°F' in input_standard_unit) & ('K' in output_standard_unit):
            self.converted_value = (input_value + 459.67) * self.scaling_factor
            
        #fahrenheit --> celsius
        if ('°F' in input_standard_unit) & ('°C' in output_standard_unit):
            self.converted_value = (input_value - 32.0) * self.scaling_factor
    
        #fahrenheit --> rankine
        if ('°F' in input_standard_unit) & ('°R' in output_standard_unit):
            self.converted_value = input_value + 459.67
            
        #rankine --> kelvin
        if ('°R' in input_standard_unit) & ('K' in output_standard_unit):
            self.converted_value = input_value * self.scaling_factor
    
        #rankine --> celsius
        if ('°R' in input_standard_unit) & ('°C' in output_standard_unit):
            self.converted_value = (input_value - 491.67) * self.scaling_factor
            
        #rankine --> fahrenheit
        if ('°R' in input_standard_unit) & ('°F' in output_standard_unit):
           self.converted_value = input_value - 459.67
           
        return

    def get_molecular_mass(self, species):

        '''get molecular mass of a specific molecule in kg mol-1'''

        s = re.findall('([A-Z][a-z]?)([0-9]*)', species)
        compoundweight = 0

        #take sum of all element atomic masses (in g mol-1) 
        for element, count in s:
            count = int(count or '1')
            compoundweight += self.atomic_masses[element] * count
    
        #convert to kg mol-1
        compoundweight = compoundweight/1e3
    
        return compoundweight

    def get_N_atoms(self, molecule, element):

        '''get integer number of atoms of a desired element in a molecule.'''

        parse_list = re.findall(r'[A-Z][a-z]*|\d+', re.sub('[A-Z][a-z]*(?![\da-z])', r'\g<0>1', molecule))
        try:
            n_atoms_element = parse_list[parse_list.index(element)+1]
        except:
            n_atoms_element = 1
    
        return int(n_atoms_element)
    
    def standardise_unit_form(self, units, unit_type):

        '''return standardised units 

        cleaned units (cleaned up version of text input)
        reference units (equivalent SI base units for given quantity)
        standard units (input units in standardised form) 
        '''

        #first strip all white spaces from string
        units = units.replace(' ', '')
        
        #list for appending reference units (to which all units are scaled around)
        reference_units = []
        
        #list for appending standardised units
        standard_units = []
        
        #list for appending scaling factors for each unit to get to equivalent SI base unit  
        scaling_factors_to_reference = []
         
        #2 separate lists to separate out operation required with scaling factors, mult_factors = factors to multiply together, div_factors = factors to iteratively divide through
        mult_factors = []
        div_factors =  []
        
        #-----------------------------------------------------
        #remove '/' or 'per' and add negative exponents
        #split units into separate list items on '/' or 'per' (1 list item if none exist)
        if '/' in units:
            split_per =  units.split('/')
        else:
            split_per =  units.split('per')
        #modify (add exponents) to only units directly after / (list indices != 0)
        #add exponent in standard notation to unit (replacing '/' or 'per')
        for ii, per in enumerate(split_per):
            if ii != 0:                
                #get int of exponent (if exists) 
                exponent = re.sub("\D", "", per)
                unit_str = re.sub('[^a-zA-Z]+', "", per)
                #if no existing exponent, then exponent is -1
                if exponent == '':
                    split_per[ii] = unit_str+'-1'
                #if have exponent, make it negative
                else:
                    split_per[ii] = unit_str+'-'+exponent
        #join split_per back together to 1 string
        units = ''.join(split_per)
        
        #-----------------------------------------------------
        #separate units string into separate distinct units (i.e. separate quantities)
        
        #pass through all combinations of SI/non-SI quantities/+prefix units and find all matches with unit string
        #take longest string matches out iteratively of string until original string is empty
        
        match_list = []
        
        #SI base units
        for quantity in self.SI_base_quantities:
            
            #does unit string contain one of the SI base units?
            if self.SI_base_quantities[quantity]['base_units'] in units:
                match_list.append(self.SI_base_quantities[quantity]['base_units'])
                scaling_factors_to_reference.append(1.0)
                reference_units.append(self.SI_base_quantities[quantity]['base_units'])
                standard_units.append(self.SI_base_quantities[quantity]['base_units'])
                #print 1,  self.SI_base_quantities[quantity]['base_units']
            
            for equiv_unit in self.SI_base_quantities[quantity]['equiv_units']:
                #does unit string contain one of the SI base equivalent units?
                if equiv_unit in units:
                    match_list.append(equiv_unit)
                    scaling_factors_to_reference.append(1.0)
                    reference_units.append(self.SI_base_quantities[quantity]['base_units'])
                    standard_units.append(self.SI_base_quantities[quantity]['base_units'])
                    #print 2, equiv_unit
        
            for prefix in self.standard_prefixes:
                for prefix_variant in self.standard_prefixes[prefix]['units']:
                    #does unit string contain one of the SI base units + prefix?
                    if prefix_variant+self.SI_base_quantities[quantity]['base_units'] in units:
                        match_list.append(prefix_variant+self.SI_base_quantities[quantity]['base_units'])
                        scaling_factors_to_reference.append(self.standard_prefixes[prefix]['factor']) 
                        reference_units.append(self.SI_base_quantities[quantity]['base_units'])
                        standard_units.append(self.standard_prefixes[prefix]['standard_units']+self.SI_base_quantities[quantity]['base_units'])
                        #print 3, prefix_variant+self.SI_base_quantities[quantity]['base_units']
                
                    for equiv_unit in self.SI_base_quantities[quantity]['equiv_units']:
                        #does unit string contain one of the SI base equivalent units + prefix?
                        if prefix_variant+equiv_unit in units:
                            match_list.append(prefix_variant+equiv_unit)
                            scaling_factors_to_reference.append(self.standard_prefixes[prefix]['factor']) 
                            reference_units.append(self.SI_base_quantities[quantity]['base_units'])
                            standard_units.append(self.standard_prefixes[prefix]['standard_units']+self.SI_base_quantities[quantity]['base_units'])
                            #print 4, prefix_variant+equiv_unit
            
        #SI derived equivalent units
        for quantity in self.SI_derived_quantities:
            for equiv_unit in self.SI_derived_quantities[quantity]['equiv_units']:
                #does unit string contain one of the SI derived equivalent units?
                if equiv_unit in units:
                    match_list.append(equiv_unit)
                    scaling_factors_to_reference.append(1.0)
                    reference_units.append(self.SI_derived_quantities[quantity]['base_units'])
                    standard_units.append(self.SI_derived_quantities[quantity]['standard_units'])
                    #print 5, equiv_unit
            
                for prefix in self.standard_prefixes:
                    for prefix_variant in self.standard_prefixes[prefix]['units']:
                        #does unit string contain one of the SI derived equivalent units + prefix?
                        if prefix_variant+equiv_unit in units:
                            match_list.append(prefix_variant+equiv_unit)
                            scaling_factors_to_reference.append(self.standard_prefixes[prefix]['factor'])
                            reference_units.append(self.SI_derived_quantities[quantity]['base_units'])
                            standard_units.append(self.standard_prefixes[prefix]['standard_units']+self.SI_derived_quantities[quantity]['standard_units'])
                            #print 6, prefix_variant+equiv_unit
            
        #non-SI units
        for quantity in self.non_SI_quantities:
            for quantity_variant in self.non_SI_quantities[quantity]:
                for unit in self.non_SI_quantities[quantity][quantity_variant]['units']:
                    #does unit string contain one of the non-SI units?
                    if unit in units:
                        match_list.append(unit)
                        scaling_factors_to_reference.append(self.non_SI_quantities[quantity][quantity_variant]['factor'])
                        standard_units.append(self.non_SI_quantities[quantity][quantity_variant]['standard_units'])
                        try:
                            reference_units.append(self.SI_base_quantities[quantity]['base_units'])
                        except:
                            reference_units.append(self.SI_derived_quantities[quantity]['base_units'])
                        #print 7, unit
                        
                        
                    for prefix in self.standard_prefixes:
                        #do not test 'grams' with prefix 'kilo', as already have that as base unit (do not want double matches)
                        if (quantity_variant == 'gram') & (prefix == 'kilo'):
                            continue
                        for prefix_variant in self.standard_prefixes[prefix]['units']:
                            #does unit string contain one of the non-SI units+prefix?
                            if prefix_variant+unit in units:
                                match_list.append(prefix_variant+unit)
                                scaling_factors_to_reference.append(self.non_SI_quantities[quantity][quantity_variant]['factor'] * self.standard_prefixes[prefix]['factor'])
                                standard_units.append(self.standard_prefixes[prefix]['standard_units']+self.non_SI_quantities[quantity][quantity_variant]['standard_units'])
                                try:
                                    reference_units.append(self.SI_base_quantities[quantity]['base_units'])
                                except:
                                    reference_units.append(self.SI_derived_quantities[quantity]['base_units'])
                                #print 8, prefix_variant+unit
         
         
        #sort match list by length of string (biggest first), then alphabetically
        #sort scaling factors/unit lists with match list 
        xyzs= zip(match_list, scaling_factors_to_reference, reference_units, standard_units)
        xyzs.sort(key=lambda item: (-len(item[0]), item[0]))
        match_list = [x for x, y, z, s in xyzs]
        scaling_factors_to_reference = [y for x, y, z, s in xyzs]
        reference_units = [z for x, y, z, s in xyzs]
        standard_units = [s for x, y, z, s in xyzs]
        
        #list for storing lists of inds for each separate unit  
        separate_unit_inds = []
        scaling_factor_inds = []
        
        #create separate copy of original units for cutting down iteratively
        cut_units = units
      
        #iteratively get inds of matches in original unit string
        for mm,match in enumerate(match_list):

            #check match is not all derived unit with negative exponent (invalid)
            #if so, skip to next match
            #stops things like --> ms-1, where should be m s-1
            try:
                if int(units.replace(match,'',1)) < 0:
                    continue
            except:
                pass

            #only proceed if match still exists in cut original units
            if match in cut_units:
                #get all inds of match occurances in unit_string
                match_inds = [i for i in range(len(units)) if units.startswith(match, i)]
                #iterate through match inds, if ind not already part of another distinct unit string then proceed
                for match_ind in match_inds:
                    if any(match_ind in x for x in separate_unit_inds) == False:
                        #get range of inds for each separate unit match, related to the original string
                        valid_inds = range(match_ind, match_ind+len(match))
                        separate_unit_inds.append(valid_inds)
                        scaling_factor_inds.append(mm)
                        cut_units = cut_units.replace(match,'')
                        #break
    
        #strip alphabet characters from cut units, if string is not empty, units are not valid, returning invalid warning
        remaining_units =  " ".join(re.findall("[a-zA-Z]+", cut_units)) 
        if remaining_units != '':
            self.error_message = 'Some or all of {} units are not recognised'.format(unit_type)
            return 1+'a'
    
        #cut scaling factors/unit lists to only be valid matches
        scaling_factors_to_reference = np.array(scaling_factors_to_reference)[scaling_factor_inds] 
        reference_units = np.array(reference_units)[scaling_factor_inds] 
        standard_units = np.array(standard_units)[scaling_factor_inds] 
    
        #sort each of list indices by first integer in each list
        #sort scaling factors/unit lists equivalently
        xyzs = zip(separate_unit_inds, scaling_factors_to_reference, reference_units, standard_units)
        xyzs.sort(key=lambda x: x[0][0])
        separate_unit_inds = [x for x, y, z, s in xyzs]
        scaling_factors_to_reference = [y for x, y, z, s in xyzs]  
        reference_units = [z for x, y, z, s in xyzs] 
        standard_units = [s for x, y, z, s in xyzs] 
        
        #append gaps between list of ind lists to anterior list (exponents)
        #also separate scaling factors into multiplication and division lists based on sign of exponent, and modify scaling factors for power of exponent
    
        for ii, current_list in enumerate(separate_unit_inds): 
            current_list_start = current_list[0]
            
            #test if difference between current list start ind & prev list end ind is > 1
            #if so then add integers missing to prev list inds
            #do not test on first iteration (don't have prev list)
            if ii != 0:
                if (current_list_start - prev_list_end) > 1:
                    missing_inds = range(prev_list_end+1,current_list_start)
                    #get exponent
                    exponent = ''
                    for m_i in missing_inds:
                        separate_unit_inds[ii-1].append(m_i)
                        exponent+=units[m_i]
                        
                    #append exponent to reference/standard units
                    reference_units[ii-1] = reference_units[ii-1]+exponent 
                    standard_units[ii-1]  = standard_units[ii-1]+exponent
                        
                    #get just numbers from exponent
                    exponent_num = int(re.sub("\D", "", exponent))
                    #check if exponent is negative
                    #if so, append scaling factor to division list
                    if '-' in exponent:
                        div_factors.append(scaling_factors_to_reference[ii-1]**exponent_num)
                    #otherwise then exponent is positive
                    #append scaling factor to multiplication list 
                    else:
                        mult_factors.append(scaling_factors_to_reference[ii-1]**exponent_num)
                
                #no exponent, do no modification to scaling factor, append to mult_factors
                else:
                    mult_factors.append(scaling_factors_to_reference[ii-1])
                    
            prev_list_end = current_list[-1]    
 
            #if last element, get difference between len of original unit string and ind+1, append ind differences
            if ii == (len(separate_unit_inds)-1):
                if (prev_list_end+1) != len(units):
                    missing_inds = range(prev_list_end+1,len(units))
                    #get exponent
                    exponent = ''
                    for m_i in missing_inds:
                        separate_unit_inds[ii].append(m_i)
                        exponent+=units[m_i]
                        
                    #append exponent to reference/standard units
                    reference_units[ii] = reference_units[ii]+exponent  
                    standard_units[ii]  = standard_units[ii]+exponent  
                    
                    #get just numbers from exponent
                    exponent_num = int(re.sub("\D", "", exponent))
                    #check if exponent is negative
                    #if so, append scaling factor to division list
                    if '-' in exponent:
                        div_factors.append(scaling_factors_to_reference[ii]**exponent_num)
                    #otherwise then exponent is positive
                    #append scaling factor to multiplication list 
                    else:
                        mult_factors.append(scaling_factors_to_reference[ii]**exponent_num)
                                
                #no exponent, do no modification to scaling factor, append to mult_factors
                else:
                    mult_factors.append(scaling_factors_to_reference[ii])            
            
                        
        #put together cleaned unit string
        cleaned_unit = ''
        for ii, current_list in enumerate(separate_unit_inds):
            cleaned_unit = cleaned_unit + units[current_list[0]:(current_list[-1]+1)]
            if ii != (len(separate_unit_inds)-1):
                cleaned_unit += ' '
                
        #put together reference unit string
        reference_unit = ' '.join(reference_units)
                
        #put together standard unit string
        standard_unit = ' '.join(standard_units)
        
        
        #multiply all scaling factors together in multiplication list
        scaling_factor_to_reference = np.prod(mult_factors)
        
        #divide calculated factor iteratively by each of scaling factors in division list to get final unified scaling factor to reference
        for div_fact in div_factors:
            scaling_factor_to_reference = scaling_factor_to_reference / div_fact
    
        #save output to self
        if unit_type == 'input':
            self.input_cleaned_units                 = cleaned_unit
            self.input_reference_units               = reference_unit
            self.input_standard_units                = standard_unit
            self.input_scaling_factor_to_reference   = scaling_factor_to_reference
        
        elif unit_type == 'output':
            self.output_cleaned_units                = cleaned_unit
            self.output_reference_units              = reference_unit
            self.output_standard_units               = standard_unit
            self.output_scaling_factor_to_reference  = scaling_factor_to_reference
                
        return 
        
    def load_dictionaries(self):
        
        '''define all key dictionaries needed for unit conversion'''

        #dictionary of standard prefixes with associated multiplying factors
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
        'deca': {'factor':1e1,   'standard_units':'da', 'units':['deca','Deca','deka','Deka','da']},
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


        #define SI base quantities with units 
        self.SI_base_quantities = {
        'amount of substance':{'unit_name':'mole',     'base_units':'mol', 'equiv_units':['moles','mols','mole']},
        'electric current'   :{'unit_name':'ampere',   'base_units':'A',   'equiv_units':['amp','amps','amperes']},
        'length'             :{'unit_name':'metre',    'base_units':'m',   'equiv_units':['metres','meters','meter','metre']},
        'luminous intensity' :{'unit_name':'candela',  'base_units':'cd',  'equiv_units':['candela','candelas']},
        'mass'               :{'unit_name':'kilogram', 'base_units':'kg',  'equiv_units':['kilograms','kilogram','kg']},
        'temperature'        :{'unit_name':'kelvin',   'base_units':'K',   'equiv_units':['kelvin','Kelvin','kelvins','Kelvins']},
        'time'               :{'unit_name':'second',   'base_units':'s',   'equiv_units':['secs','seconds','sec','second']},
        #bnak unit only used to prevent confusion between units for rates (i.e. frequency - s-1)
        'blank'              :{'unit_name':'blank',    'base_units':'xx',  'equiv_units':['none']}}                     
 
        #define common SI derived quantities with units (formed by powers, products or quotients of the base units and are unlimited in number)
        self.SI_derived_quantities ={
        'absorbed_dose':                {'unit_name':'gray',                           'base_units':'m2 s-2',              'standard_units':'Gy',                  'equiv_units':['Gy','gray','grays']},
        'acceleration':                 {'unit_name':'metre per second squared',       'base_units':'m s-2',               'standard_units':'m s-2',               'equiv_units':['none']},
        'action':                       {'unit_name':'joule second',                   'base_units':'m2 kg s-1',           'standard_units':'m2 kg s-1',           'equiv_units':['none']},
        'angular_acceleration':         {'unit_name':'radian per second squared',      'base_units':'xx s-2',              'standard_units':'s-2',                 'equiv_units':['none']},
        'angular_momentum':             {'unit_name':'newton metre second',            'base_units':'m2 kg s-1',           'standard_units':'m2 kg s-1',           'equiv_units':['none']},
        'angular_velocity':             {'unit_name':'radian per second',              'base_units':'xx s-1',              'standard_units':'s-1',                 'equiv_units':['none']},
        'angle':                        {'unit_name':'radian',                         'base_units':'m m-1',               'standard_units':'rad',                 'equiv_units':['rad','rads','radian','radians']},
        'area':                         {'unit_name':'square metre',                   'base_units':'m2',                  'standard_units':'m2',                  'equiv_units':['none']},
        'area_density':                 {'unit_name':'kilogram per square metre',      'base_units':'m-2 kg',              'standard_units':'m-2 kg',              'equiv_units':['none']},
        'catalytic_activity':           {'unit_name':'katal',                          'base_units':'s-1 mol',             'standard_units':'kat',                 'equiv_units':['kat','katal','katals']},
        'electrical_capacitance':       {'unit_name':'farad',                          'base_units':'kg-1 m-2 s4 A2',      'standard_units':'F',                   'equiv_units':['F','farad','farads']},
        'electric_charge':              {'unit_name':'coulomb',                        'base_units':'s A',                 'standard_units':'C',                   'equiv_units':['C','coulomb','coulombs']},
        'electrical_conductance':       {'unit_name':'siemens',                        'base_units':'kg-1 m-2 s3 A2',      'standard_units':'S',                   'equiv_units':['S','siemens']},
        'electrical_inductance':        {'unit_name':'henry',                          'base_units':'kg m2 s-2 A-2',       'standard_units':'H',                   'equiv_units':['H','henry','henries']},
        'electrical_resistance':        {'unit_name':'ohm',                            'base_units':'kg m2 s-3 A-2',       'standard_units':'Ω',                   'equiv_units':['Ω','ohm','ohms']},
        'energy':                       {'unit_name':'joule',                          'base_units':'kg m2 s-2',           'standard_units':'J',                   'equiv_units':['J','joule','joules']},
        'equivalent_dose':              {'unit_name':'sievert',                        'base_units':'m2 s-2',              'standard_units':'Sv',                  'equiv_units':['Sv','sievert','sieverts']},
        'frequency':                    {'unit_name':'hertz',                          'base_units':'xx s-1',              'standard_units':'Hz',                  'equiv_units':['Hz','hertz']},
        'force':                        {'unit_name':'newton',                         'base_units':'kg m s-2',            'standard_units':'N',                   'equiv_units':['N','newton','newtons']},
        'heat_capacity':                {'unit_name':'joule per kelvin',               'base_units':'m2 kg s-2 K-1',       'standard_units':'m2 kg s-2 K-1',       'equiv_units':['none']},
        'jerk':                         {'unit_name':'metre per second cubed',         'base_units':'m s-3',               'standard_units':'m s-3',               'equiv_units':['none']},
        'illuminance':                  {'unit_name':'lux',                            'base_units':'m-2 cd',              'standard_units':'lx',                  'equiv_units':['lx','lux','luxes']},
        'luminous_flux':                {'unit_name':'lumen',                          'base_units':'cd',                  'standard_units':'lm',                  'equiv_units':['lm','lumen','lumens']},
        'magnetic_flux':                {'unit_name':'weber',                          'base_units':'kg m2 s-2 A-1',       'standard_units':'Wb',                  'equiv_units':['Wb','weber','webers']},
        'magnetic_field_strength':      {'unit_name':'tesla',                          'base_units':'kg s-2 A-1',          'standard_units':'T',                   'equiv_units':['T','tesla']},
        'mass_density':                 {'unit_name':'kilogram per cubic metre',       'base_units':'m-3 kg',              'standard_units':'m-3 kg',              'equiv_units':['none']},
        'molarity':                     {'unit_name':'mole per cubic metre',           'base_units':'mol cm-3',            'standard_units':'M' ,                  'equiv_units':['Molar','molar','M']},
        'molar_heat_capacity':          {'unit_name':'joule per kelvin mole',          'base_units':'m2 kg s-2 K-1 mol-1', 'standard_units':'m2 kg s-2 K-1 mol-1', 'equiv_units':['none']},
        'molar_mass':                   {'unit_name':'kilogram per mole',              'base_units':'kg mol-1',            'standard_units':'kg mol-1',            'equiv_units':['none']},
        'molar_volume':                 {'unit_name':'cubic metre per mole',           'base_units':'m3 mol-1',            'standard_units':'m3 mol-1',            'equiv_units':['none']},
        'momentum':                     {'unit_name':'newton second',                  'base_units':'m kg s-1',            'standard_units':'m kg s-1' ,           'equiv_units':['none']},
        'mole_fraction':                {'unit_name':'mole per mole',                  'base_units':'mol mol-1',           'standard_units':'mol mol-1',           'equiv_units':['none']},
        'number_density':               {'unit_name':'number per cubic metre',         'base_units':'xx m-3',              'standard_units':'m-3',                 'equiv_units':['none']},
        'power':                        {'unit_name':'watt',                           'base_units':'kg m2 s-3',           'standard_units':'W',                   'equiv_units':['W','watt','watts']},
        'pressure':                     {'unit_name':'pascal',                         'base_units':'kg m-1 s−2',          'standard_units':'Pa',                  'equiv_units':['Pa','pascal','pascals']},
        'radioactivity':                {'unit_name':'becquerel',                      'base_units':'xx s-1',              'standard_units':'Bq',                  'equiv_units':['Bq','becquerel','becquerels']},
        'snap':                         {'unit_name':'metre per second to the fourth', 'base_units':'m s-4',               'standard_units':'m s-4',               'equiv_units':['none']},
        'solid_angle':                  {'unit_name':'steradian',                      'base_units':'m2 m-2',              'standard_units':'sr',                  'equiv_units':['sr','steradian','steradians']},
        'specific_heat_capacity':       {'unit_name':'joule per kilogram kelvin',      'base_units':'m2 s-2 K-1',          'standard_units':'m2 s-2 K-1',          'equiv_units':['none']},
        'specific_volume':              {'unit_name':'cubic metre per kilogram',       'base_units':'m3 kg-1',             'standard_units':'m3 kg-1',             'equiv_units':['none']},
        'spectral_irradiance':          {'unit_name':'watt per cubic metre',           'base_units':'m-1 kg s-3',          'standard_units':'m-1 kg s-3',          'equiv_units':['none']},
        'temperature_gradient':         {'unit_name':'kelvin per metre',               'base_units':'m-1 K',               'standard_units':'m-1 K',               'equiv_units':['none']},
        'thermal_resistance':           {'unit_name':'kelvin per watt',                'base_units':'m-2 kg-1 s3 K',       'standard_units':'m-2 kg-1 s3 K',       'equiv_units':['none']},
        'thermal_expansion_coefficient':{'unit_name':'reciprocal kelvin',              'base_units':'xx K-1',              'standard_units':'K-1',                 'equiv_units':['none']},
        'torque':                       {'unit_name':'newton metre',                   'base_units':'m2 kg s-2',           'standard_units':'m2 kg s-2',           'equiv_units':['none']},
        'velocity':                     {'unit_name':'metre per second',               'base_units':'m s-1',               'standard_units':'m s-1',               'equiv_units':['none']},
        'voltage':                      {'unit_name':'volt',                           'base_units':'kg m2 s-3 A-1',       'standard_units':'V',                   'equiv_units':['V','volt','volts']},
        'volume':                       {'unit_name':'cubic metre',                    'base_units':'m3',                  'standard_units':'m3',                  'equiv_units':['none']},
        'volumetric_flow':              {'unit_name':'cubic metre per second',         'base_units':'m3 s-1',              'standard_units':'m3 s-1',              'equiv_units':['none']},
        'wavenumber':                   {'unit_name':'reciprocal metre',               'base_units':'xx m-1',              'standard_units':'m-1',                 'equiv_units':['none']},
        'yank':                         {'unit_name':'newton per second',              'base_units':'m kg s-3',            'standard_units':'m kg s-3',            'equiv_units':['none']}}


        #map non-SI units to SI quantities with conversion factors
        #can be more 1 type of non-SI units for each SI quantity
        self.non_SI_quantities = {
        'area':         {'acre':        {'factor':4046.8564224,     'standard_units':'acre',       'units':['acre','acres']}},
        'energy':       {'electronvolt':{'factor':1.6021766208e-19, 'standard_units':'eV',         'units':['eV','electronvolt','electronvolts']}},
        'length':       {'mile':        {'factor':1609.344,         'standard_units':'miles',      'units':['mi','miles']},                                                                   'feet':          {'factor':0.3048,      'standard_units':'ft',             'units':['feet','ft']},                                                                    'inch':            {'factor':0.0254,     'standard_units':'in',         'units':['in','inch','inches']},                                                                               'yard':          {'factor':0.9144,                                                                        'standard_units':'yd',          'units':['yd','yds','yards','yard']},                                              'astronomical unit':       {'factor':149597870700,                                                                  'standard_units':'au',          'units':['au']}},
        'mass':         {'gram':        {'factor':1e-3,             'standard_units':'g',          'units':['g','gram','grams']},                                                             'metric ton':    {'factor':1e3,         'standard_units':'ton',            'units':['ton','tons','tonne','tonnes']},                                                  'pound':           {'factor':0.45359237, 'standard_units':'lb',         'units':['pounds','lb','lbs']},                                                                                'ounce':         {'factor':0.02834952,                                                                    'standard_units':'oz',          'units':['oz','ozs','ounce','ounces']},                                            'unified atomic mass unit':{'factor':1.660539040e-27,                                                               'standard_units':'u',           'units':['u','Da']},                                                                   'gram of nitrogen':{'factor':1e-3*(self.get_molecular_mass(self.measured_species)/(get_molecular_mass('N')*self.get_N_atoms(self.measured_species,'N'))), 'standard_units':'gN',          'units':['gN','gramN','gramsN']},                                                                                                    'gram of carbon':{'factor':1e-3*(self.get_molecular_mass(self.measured_species)/(self.get_molecular_mass('C')*self.get_N_atoms(self.measured_species,'C'))), 'standard_units':'gC',          'units':['gC','gramC','gramsC']},                                                  'gram of sulphur': {'factor':1e-3*(self.get_molecular_mass(self.measured_species)/(self.get_molecular_mass('S')*self.get_N_atoms(self.measured_species,'S'))), 'standard_units':'gS',          'units':['gS','gramS','gramsS']}},       
        #note - volumetric fractions and the mole fractions of the components of an ideal gas mixture are interchangeable --> does not hold for liquids
        'mole_fraction':{'ppmv':        {'factor':1e-6,             'standard_units':'umol mol-1', 'units':['ppmv','ppm','Partspermillion','partspermillion','ppmV']},                        'ppbv':          {'factor':1e-9,        'standard_units':'nmol mol-1',     'units':['ppbv','ppb','Partsperbillion','partsperbillion','ppbV']},                        'pptv':            {'factor':1e-12,      'standard_units':'pmol mol-1', 'units':['pptv','ppt','Partspertrillion','partspertrillion','pptV']},                                          'ppmv of carbon':{'factor':1e-6/(self.get_N_atoms(self.measured_species,'C')), 'standard_units':'umolC mol-1', 'units':['ppmvC','ppmC','partspermillioncarbon','PartspermillionCarbon','ppmVC']}, 'ppmv of nitrogen':        {'factor':1e-6/(self.get_N_atoms(self.measured_species,'N')), 'standard_units':'umolN mol-1', 'units':['ppmvN','ppmN','partspermillionnitrogen','PartspermillionNitrogen','ppmVN']}, 'ppmv of sulphur': {'factor':1e-6/(self.get_N_atoms(self.measured_species,'S')),                                                                                                                                     'standard_units':'umolS mol-1', 'units':['ppmvS','ppmS','partspermillionsulphur','PartspermillionSulphur','partspermillionsulfur','PartspermillionSulfur','ppmVS']}, 'ppbv of carbon':{'factor':1e-9/(self.get_N_atoms(self.measured_species,'C')),                                                                                                                                     'standard_units':'nmolC mol-1', 'units':['ppbvC','ppbC','partsperbillioncarbon','PartsperbillionCarbon','ppbVC']}, 'ppbv of nitrogen':{'factor':1e-9/(self.get_N_atoms(self.measured_species,'N')),                                                                                                                                     'standard_units':'nmolN mol-1', 'units':['ppbvN','ppbN','partsperbillionnitrogen','PartsperbillionNitrogen','ppbVN']}, 'ppbv of sulphur': {'factor':1e-9/(self.get_N_atoms(self.measured_species,'S')), 'standard_units':'nmolS mol-1', 'units':['ppbvS','ppbS','partsperbillionsulphur','PartsperbillionSulphur','partsperbillionsulfur','PartsperbillionSulfur','ppbVS']}, 'pptv of carbon':{'factor':1e-12/(self.get_N_atoms(self.measured_species,'C')), 'standard_units':'pmolC mol-1', 'units':['pptvC','pptC','partspertrillioncarbon','PartspertrillionCarbon','pptVC']}, 'pptv of nitrogen':{'factor':1e-12/(self.get_N_atoms(self.measured_species,'N')), 'standard_units':'pmolN mol-1', 'units':['pptvN','pptN','partspertrillionnitrogen','PartspertrillionNitrogen','pptVN']}, 'pptv of sulphur': {'factor':1e-12/(self.get_N_atoms(self.measured_species,'S')), 'standard_units':'pmolS mol-1', 'units':['pptvS','pptS','partspertrillionsulphur','PartspertrillionSulphur','partspertrillionsulfur','PartspertrillionSulfur','pptVS']}},
        'pressure':     {'atmosphere':  {'factor':101325,           'standard_units':'atm',        'units':['atm', 'atmospheres']},                                                           'bar':           {'factor':1e5,         'standard_units':'bar',            'units':['bar']},                                                                          'torr':            {'factor':133.3224,   'standard_units':'Torr',       'units':['Torr','torr']},                                                                                      'psi':           {'factor':6894.76,                                                                       'standard_units':'psi',         'units':['psi']}},  
        'temperature':  {'celsius':     {'factor':1.0,              'standard_units':'°C',         'units':['degC','dC','celsius','°C','°Celsius','degreesC','degCelsius','degreesCelsius']}, 'rankine':       {'factor':5./9.,       'standard_units':'°R',             'units':['degR','dR','rankine','°R','°Rankine','degreesR','degRankine','degreesRankine']}, 'fahrenheit':      {'factor':5./9.,      'standard_units':'°F',         'units':['degF', 'dF', 'fahrenheit', '°F', '°Fahrenheit', 'degreesF', 'degFahrenheit', 'degreesFahrenheit']}},
        'time':         {'minute':      {'factor':60.0,             'standard_units':'mins',       'units':['min','mins','minutes']},                                                         'hour':          {'factor':3600.0,      'standard_units':'hours',          'units':['h','hr','hrs','hours']},                                                         'day':             {'factor':86400.0,    'standard_units':'days',       'units':['day','days']}},
        'volume':       {'litre':       {'factor':1e-3,             'standard_units':'l',          'units':['l','litres','liters','L']},                                                      'imperial pint': {'factor':0.000568261, 'standard_units':'imperial_pints', 'units':['pt','pint','pints']},                                                            'imperial gallon': {'factor':0.00454609, 'standard_units':'gal',        'units':['gal','gallon','gals','gallons']}}}


        #important constants
        self.constants = {
        'avogadro_constant' : {'value':6.022140857e23, 'base_units':'mol-1',               'equiv_units':'L'},
        'boltzmann_constant': {'value':1.38064852e-23, 'base_units':'m2 kg s-2 K-1',       'equiv_units':'k'},
        'molar_gas_constant': {'value':8.3144598,      'base_units':'kg m2 s-2 K-1 mol-1', 'equiv_units':'R'},
        'planck_constant'   : {'value':6.62606979e-34, 'base_units':'m2 kg s-1',           'equiv_units':'h'}}

        #quantity conversion formulae
        self.conversion_formulae = {
        'mole_fraction':                                ["mass_density*((molar_gas_constant*temperature)/(molar_mass*pressure))",        "number_density/((avogadro_constant*pressure)/(molar_gas_constant*temperature))"],
        'mass_density':                                 ["mole_fraction/((molar_gas_constant*temperature)/(molar_mass*pressure))",       "(number_density*molar_mass)/avogadro_constant"],
        'number_density':                               ["mole_fraction*((avogadro_constant*pressure)/(molar_gas_constant*temperature))","(mass_density*avogadro_constant)/molar_mass"],
        'mass_density>mole_fraction--conversion_factor':"((molar_gas_constant*temperature)/(molar_mass*pressure))",
        'mole_fraction>mass_density--conversion_factor':"mole_fraction/((molar_gas_constant*temperature)/(molar_mass*pressure))"}

        #define atomic masses per element in g mol-1
        self.atomic_masses = {
        "H":1.00794,"He":4.002602,"Li":6.941,"Be":9.012182,"B":10.811,"C":12.011,"N":14.00674,"O":15.9994,"F":18.9984032,"Ne":20.1797,"Na":22.989768,"Mg":24.3050,"Al":26.981539,"Si":28.0855,"P":30.973762,"S":32.066,"Cl":35.4527,"Ar":39.948,
        "K":39.0983,"Ca":40.078,"Sc":44.955910,"Ti":47.88,"V":50.9415,"Cr":51.9961,"Mn":54.93805,"Fe":55.847,"Co":58.93320,"Ni":58.6934,"Cu":63.546,"Zn":65.39,"Ga":69.723,"Ge":72.61,"As":74.92159,"Se":78.96,"Br":79.904,"Kr":83.80,
        "Rb":85.4678,"Sr":87.62,"Y":88.90585,"Zr":91.224,"Nb":92.90638,"Mo":95.94,"Tc":98,"Ru":101.07,"Rh":102.90550,"Pd":106.42,"Ag":107.8682,"Cd":112.411,"In":114.82,"Sn":118.710,"Sb":121.757,"Te":127.60,"I":126.90447,"Xe":131.29,
        "Cs":132.90543,"Ba":137.327,"La":138.9055,"Ce":140.115,"Pr":140.90765,"Nd":144.24,"Pm":145,"Sm":150.36,"Eu":151.965,"Gd":157.25,"Tb":158.92534,"Dy":162.50,"Ho":164.93032,"Er":167.26,"Tm":168.93421,"Yb":173.04,"Lu":174.967,
        "Hf":178.49,"Ta":180.9479,"W":183.85,"Re":186.207,"Os":190.2,"Ir":192.22,"Pt":195.08,"Au":196.96654,"Hg":200.59,"Tl":204.3833,"Pb":207.2,"Bi":208.98037,"Po":209,"At":210,"Rn":222,"Fr":223,"Ra":226.0254,"Ac":227,"Th":232.0381,
        "Pa":213.0359,"U":238.0289,"Np":237.0482,"Pu":244,"Am":243,"Cm":247,"Bk":247,"Cf":251,"Es":252,"Fm":257,"Md":258,"No":259,"Lr":260,"Rf":261,"Db":262,"Sg":263,"Bh":262,"Hs":265,"Mt":266}