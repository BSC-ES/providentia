Python program that converts units (specifically tailored for conversions of components in the field of Atmospheric Composition) 

accepts mixture of SI base/derived/equivalent units and non-SI units as input/output arguments
can handle complex formatted input/output units (i.e. without spaces, spaced, including '/' or 'per' terms)

converts units between same quantities (i.e. length --> length, cm --> m)
and converts between different quantities (i.e. mass density --> mole fraction, mg/m3 --> ppbv) **when given the correct parameters needed/ and have formula

**units which are solely per something e.g. frequency (s-1) are prefixed with 'xx' to prevent confusion between units
i.e. ms-1 (per miliseconds) could be taken to be m s-1 (metres per second) or vice versa in processing  

requires input of: 
1. input units -- either single string for same quantity conversion, or a dictionary of input quantity name and units for converting between quantities {'quantity':unit1, 'quantity2':unit2}  
**quantity names should be given with underscore joins when more 1 word 

2. output units

3. input value/s -- either a single value for same quantity conversion, or a dictionary of input quantity name and values for converting between quantities {'quantity':value, 'quantity2':value2} 
**quantity names should be given with underscore joins when more 1 word 

4. Level of precision of converted output/conversion factor (fixed number of decimal places) **OPTIONAL 

5. Measured species (in chemical notation). **OPTIONAL -- needed when converting units associated with an element (e.g. ppmv N) 

6. Conversion input quantity **OPTIONAL -- used when wanting to return conversion factor for conversions with more than 1 input variables, i.e. this quantity refers to the input variable for which a conversion factor to the output units is wanted

**Examples**

**Same quantity unit conversion (m3 --> dm3)**

1. Define input variables
```
input_units = 'm3'
output units = 'dm3'
input_value = 1.0
```

2. Do conversion
```
conv_obj = unit_converter.convert_units(input_units, output_units, input_value)
conv_obj.do_conversion()
```

**Different unit quantity conversion (ug m-3 --> ppbv)**

1. Define input variables
```
input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mass_density':'ug m-3'}
output_units = 'ppbv'
input_values = {'temperature':273.0, 'pressure':1013.0, 'molar_mass':48.0, 'mass_density':10.0}
```

2. Do conversion
```
conv_obj = unit_converter.convert_units(input_units, output_units, input_values, conversion_input_quantity='mass_density')
conv_obj.do_conversion()
```

**Output Class Attributes**

```
conv_obj.converted_value             --> Converted output value 
conv_obj.conversion_factor           --> Conversion factor to get from input value to output value (*Note this is not valid for temperature parameters)
conv_obj.input_cleaned_units         --> Cleaned version of input units 
conv_obj.output_cleaned_units        --> Cleaned version of output units 
conv_obj.input_reference_units       --> Input units as reference units (equivalent SI base units for given quantity)
conv_obj.output_reference_units      --> Output units as reference units (equivalent SI base units for given quantity)
conv_obj.input_standard_units        --> Input units as standard units (standardised form of units)
conv_obj.output_standard_units       --> Output units as standard units (standardised form of units)
conv_obj.input_represented_quantity  --> Quantity of input units
conv_obj.output_represented_quantity --> Quantity of output units
```
