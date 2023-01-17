import numpy as np
import sys
import unit_converter

test_dict = {'1':{'orig_unit':'g', 'orig_val':500, 'new_unit':'kg', 'new_val':0.5},
             '2':{'orig_unit':'K', 'orig_val':3, 'new_unit':'degC', 'new_val':-270.15},
             '3':{'orig_unit':'degF', 'orig_val':50, 'new_unit':'K', 'new_val':283.15},
             '4':{'orig_unit':'m/s2', 'orig_val':3000, 'new_unit':'km/s2', 'new_val':3},
             '5':{'orig_unit':'kg/m3', 'orig_val':18, 'new_unit':'pg/um3', 'new_val':0.018},
             '6':{'orig_unit':'W/J', 'orig_val':25, 'new_unit':'mW/kJ', 'new_val':2.5e+7},
             '7':{'orig_unit':'inches', 'orig_val':3, 'new_unit':'cm', 'new_val':7.62},
             '8':{'orig_unit':'feet', 'orig_val':5, 'new_unit':'m', 'new_val':1.524}
            }
             
for idx in sorted(list(test_dict.keys())):
    conv_obj = unit_converter.convert_units(test_dict[idx]['orig_unit'],test_dict[idx]['new_unit'],test_dict[idx]['orig_val'])
    conv_obj.do_conversion()
    if np.isclose(conv_obj.converted_value,test_dict[idx]['new_val']) == False:
        sys.exit('{} {} to {} should be {}, is {}'.format(test_dict[idx]['orig_val'],test_dict[idx]['orig_unit'],test_dict[idx]['new_unit'],test_dict[idx]['new_val'],conv_obj.converted_value))
    else:
        print('{} {} to {} {} correctly converted'.format(test_dict[idx]['orig_val'],test_dict[idx]['orig_unit'],test_dict[idx]['new_val'],test_dict[idx]['new_unit']))

orig_unit = 'ug/m3'
orig_val = 1.9957
new_unit = 'ppbv'
new_val = 1.0
input_units ={'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mass_density':orig_unit}
input_values = {'temperature':293.15, 'pressure':1013.25, 'molar_mass':47.998, 'mass_density':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, conversion_input_quantity='mass_density')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppbv'
orig_val = 1.0
new_unit = 'ug/m3'
new_val = 1.9957
input_units ={'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':293.15, 'pressure':1013.25, 'molar_mass':47.998, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, conversion_input_quantity='mole_fraction')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 150.0
new_unit = 'mg/m3'
new_val = 294.51
input_units ={'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.0, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, conversion_input_quantity='mole_fraction')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 150.0
new_unit = 'mgC/m3'
new_val = 240.0
input_units ={'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.0, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, conversion_input_quantity='mole_fraction', measured_species='C3H8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'mg/m3'
orig_val = 294.51
new_unit = 'mgC/m3'
new_val = 240.0
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='C3H8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'mgC/m3'
orig_val = 240.0
new_unit = 'mg/m3'
new_val = 294.51
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='C3H8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppbv'
orig_val = 10.0
new_unit = 'umol/mol'
new_val = 0.01
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='C3H8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ug/m3'
orig_val = 13.0
new_unit = 'mg m-3'
new_val = 0.013
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='C3H8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmvC'
orig_val = 3
new_unit = 'ppmv'
new_val = 1
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='C3H8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmvC'
orig_val = 1
new_unit = 'umol mol-1'
new_val = 1
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='CH8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 1
new_unit = 'ppmvC'
new_val = 3
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='C3H8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'umolC mol-1'
orig_val = 3
new_unit = 'ppmv'
new_val = 1
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='C3H8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'umolC mol-1'
orig_val = 1
new_unit = 'umol mol-1'
new_val = 1
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='CH8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 1
new_unit = 'umolC mol-1'
new_val = 3
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,measured_species='C3H8')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 450.0
new_unit = 'mg/m3'
new_val = 884.0
input_units ={'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.01, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, conversion_input_quantity='mole_fraction')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 450.0
new_unit = 'mg/L'
new_val = 0.884
input_units ={'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.01, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, conversion_input_quantity='mole_fraction')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'mg/m3'
orig_val = 884.0
new_unit = 'ppmv'
new_val = 450.0
input_units ={'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mass_density':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.01, 'mass_density':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, conversion_input_quantity='mass_density')
conv_obj.do_conversion()
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))
