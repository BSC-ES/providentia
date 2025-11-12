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
conv_obj = unit_converter.convert_units(input_units, output_units, input_value)

-----------------------------------------------------
Different unit quantity conversion (e.g. ug m-3 --> ppbv)

1. Define input variables 
input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mass_concentration':'ug m-3'}
output_units = 'ppbv'
input_values = {'temperature':273.0, 'pressure':1013.0, 'molar_mass':48.0, 'mass_concentration':10.0}

2. Do conversion
conv_obj = unit_converter.convert_units(input_units, output_units, input_values, 
                                        input_quantity='mass_concentration')

----------------------------------------------------------------------------------------
Output Class Attributes

conv_obj.converted_value                   --> Converted output value 
conv_obj.conversion_factor                 --> Conversion factor to get from input value to output value (*Note this is not valid for temperature parameters)
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