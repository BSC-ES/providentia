import numpy as np
import os
import sys

sys.path.insert(1, os.path.join(sys.path[0], '..'))
import unit_converter

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Same quantity unit conversion (e.g. m3 --> dm3)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

print('\nSAME QUANTITY UNIT CONVERSION TESTS')

test_dict = {'1': {'orig_unit':'g', 'orig_val':500, 'new_unit':'kg', 'new_val':0.5, 'input_quantity':'mass', 'output_quantity':'mass'},
             '2': {'orig_unit':'K', 'orig_val':3, 'new_unit':'degC', 'new_val':-270.15, 'input_quantity':'temperature', 'output_quantity':'temperature'},
             '3': {'orig_unit':'degF', 'orig_val':50, 'new_unit':'K', 'new_val':283.15, 'input_quantity':'temperature', 'output_quantity':'temperature'},
             '4': {'orig_unit':'m/s2', 'orig_val':3000, 'new_unit':'km/s2', 'new_val':3, 'input_quantity':'acceleration', 'output_quantity':'acceleration'},
             '5': {'orig_unit':'kg/m3', 'orig_val':18, 'new_unit':'pg/um3', 'new_val':0.018, 'input_quantity':'mass_concentration', 'output_quantity':'mass_concentration'},
             '6': {'orig_unit':'W/J', 'orig_val':25, 'new_unit':'mW/kJ', 'new_val':2.5e+7, 'input_quantity':'frequency', 'output_quantity':'frequency'},
             '7': {'orig_unit':'inches', 'orig_val':3, 'new_unit':'cm', 'new_val':7.62, 'input_quantity':'length', 'output_quantity':'length'},
             '8': {'orig_unit':'feet', 'orig_val':5, 'new_unit':'m', 'new_val':1.524, 'input_quantity':'length', 'output_quantity':'length'},
             '9': {'orig_unit':'s-1', 'orig_val':1000, 'new_unit':'ms-1', 'new_val':1, 'input_quantity':'radioactivity', 'output_quantity':'radioactivity'},
             '10':{'orig_unit':'Hz', 'orig_val':3, 'new_unit':'MHz', 'new_val':0.000003, 'input_quantity':'frequency', 'output_quantity':'frequency'},
             '11':{'orig_unit':'m-1', 'orig_val':25, 'new_unit':'um-1', 'new_val':0.000025, 'input_quantity':'backscattering_coefficient', 'output_quantity':'backscattering_coefficient'},
             '12':{'orig_unit':'umol mol-1', 'orig_val':100.0, 'new_unit':'mmol mol-1', 'new_val':0.1, 'input_quantity':'mole_fraction', 'output_quantity':'mole_fraction'},
             '13':{'orig_unit':'kg kg-1', 'orig_val':1.0, 'new_unit':'g kg-1', 'new_val':1000.0, 'input_quantity':'mass_fraction', 'output_quantity':'mass_fraction'},
             '14':{'orig_unit':'m-1', 'orig_val':6.0, 'new_unit':'um-1', 'new_val':6e-6, 'input_quantity':'absorption_coefficient', 'output_quantity':'absorption_coefficient'},
             '15':{'orig_unit':'um-1', 'orig_val':6e-6, 'new_unit':'m-1', 'new_val':6.0, 'input_quantity':'scattering_coefficient', 'output_quantity':'scattering_coefficient'},
             '16':{'orig_unit':'m2 kg-1', 'orig_val':15.0, 'new_unit':'km2 g-1', 'new_val':1.5e-8, 'input_quantity':'mass_absorption_cross_section', 'output_quantity':'mass_absorption_cross_section'},
             '17':{'orig_unit':'m m-1', 'orig_val':22.0, 'new_unit':'mm mm-1', 'new_val':22.0, 'input_quantity':'single_scattering_albedo', 'output_quantity':'single_scattering_albedo'},
             '18':{'orig_unit':'unitless', 'orig_val':22.0, 'new_unit':'unitless', 'new_val':22.0, 'input_quantity':'single_scattering_albedo', 'output_quantity':'single_scattering_albedo'},
             '19':{'orig_unit':'m m-1', 'orig_val':22.0, 'new_unit':'m km-1', 'new_val':2.2e4, 'input_quantity':'single_scattering_albedo', 'output_quantity':'single_scattering_albedo'}, 
             '20':{'orig_unit':'mol m-3 mol-1 m3', 'orig_val':7.0, 'new_unit':'nmol mm-3 nmol-1 mm3', 'new_val':7.0, 'input_quantity':'pH', 'output_quantity':'pH'},
             '21':{'orig_unit':'unitless', 'orig_val':7.0, 'new_unit':'unitless', 'new_val':7.0, 'input_quantity':'pH', 'output_quantity':'pH'},
             '22':{'orig_unit':'kg m-2 s-1', 'orig_val':13.0, 'new_unit':'g mm-2 hour-1', 'new_val':46.8, 'input_quantity':'mass_flux', 'output_quantity':'mass_flux'},
             '23':{'orig_unit':'kg m-1 s-2 kg-1 m s2', 'orig_val':84.0, 'new_unit':'g mm-1 min-2 g-1 mm min2', 'new_val':84.0, 'input_quantity':'relative_humidity', 'output_quantity':'relative_humidity'},
             '24':{'orig_unit':'unitless', 'orig_val':84.0, 'new_unit':'unitless', 'new_val':84.0, 'input_quantity':'relative_humidity', 'output_quantity':'relative_humidity'},
             '25':{'orig_unit':'m m-1', 'orig_val':64.0, 'new_unit':'um um-1', 'new_val':64.0, 'input_quantity':'optical_depth', 'output_quantity':'optical_depth'},
             '26':{'orig_unit':'unitless', 'orig_val':64.0, 'new_unit':'unitless', 'new_val':64.0, 'input_quantity':'optical_depth', 'output_quantity':'optical_depth'},
             '27':{'orig_unit':'m m-1', 'orig_val':64.0, 'new_unit':'um m-1', 'new_val':6.4e7, 'input_quantity':'optical_depth', 'output_quantity':'optical_depth'},
             '28':{'orig_unit':'um hours-1 um-1 hours', 'orig_val':28.0, 'new_unit':'m s-1 m-1 s', 'new_val':28.0, 'input_quantity':'refractive_index', 'output_quantity':'refractive_index'},
             '29':{'orig_unit':'unitless', 'orig_val':28.0, 'new_unit':'unitless', 'new_val':28.0, 'input_quantity':'refractive_index', 'output_quantity':'refractive_index'},
             '30':{'orig_unit':'m3 m-2', 'orig_val':70.0, 'new_unit':'mm3 km-2', 'new_val':7e16, 'input_quantity':'volume_concentration', 'output_quantity':'volume_concentration'},
             '31':{'orig_unit':'hours', 'orig_val':24.0, 'new_unit':'minutes', 'new_val':1440.0, 'input_quantity':'time', 'output_quantity':'time'},
             '32':{'orig_unit':'angular degrees', 'orig_val':270.0, 'new_unit':'radians', 'new_val':4.7124, 'input_quantity':'angle', 'output_quantity':'angle'},   
             '33':{'orig_unit':'m s-1', 'orig_val':3.0, 'new_unit':'um min-1', 'new_val':1.8e8, 'input_quantity':'velocity', 'output_quantity':'velocity'},
             '34':{'orig_unit':'m-2', 'orig_val':25.0, 'new_unit':'km-2', 'new_val':2.5e7, 'input_quantity':'population_density', 'output_quantity':'population_density'},
             '35':{'orig_unit':'km-2', 'orig_val':2.0, 'new_unit':'nm-2', 'new_val':2e-24, 'input_quantity':'column_number_density', 'output_quantity':'column_number_density'},
             '36':{'orig_unit':'m3 s-1', 'orig_val':21.0, 'new_unit':'nm3 min-1', 'new_val':1.26e30, 'input_quantity':'volumetric_flow', 'output_quantity':'volumetric_flow'},
             '37':{'orig_unit':'m2', 'orig_val':17.0, 'new_unit':'mm2', 'new_val':1.7e7, 'input_quantity':'area', 'output_quantity':'area'},
             '38':{'orig_unit':'unitless', 'orig_val':1.0, 'new_unit':'unitless', 'new_val':1.0, 'input_quantity':'unitless', 'output_quantity':'unitless'}
            }

for idx in list(test_dict.keys()):
    conv_obj = unit_converter.convert_units(test_dict[idx]['orig_unit'],test_dict[idx]['new_unit'],test_dict[idx]['orig_val'],input_quantity=test_dict[idx]['input_quantity'],output_quantity=test_dict[idx]['output_quantity'])
    if np.isclose(conv_obj.converted_value,test_dict[idx]['new_val']) == False:
        sys.exit('{} {} to {} should be {}, is {} (given quantities)'.format(test_dict[idx]['orig_val'],test_dict[idx]['orig_unit'],test_dict[idx]['new_unit'],test_dict[idx]['new_val'],conv_obj.converted_value))
    else:
        print('{} {} to {} {} correctly converted (given quantities)'.format(test_dict[idx]['orig_val'],test_dict[idx]['orig_unit'],test_dict[idx]['new_val'],test_dict[idx]['new_unit']))

    conv_obj = unit_converter.convert_units(test_dict[idx]['orig_unit'],test_dict[idx]['new_unit'],test_dict[idx]['orig_val'])
    if np.isclose(conv_obj.converted_value,test_dict[idx]['new_val']) == False:
        sys.exit('{} {} to {} should be {}, is {} (no given quantities)'.format(test_dict[idx]['orig_val'],test_dict[idx]['orig_unit'],test_dict[idx]['new_unit'],test_dict[idx]['new_val'],conv_obj.converted_value))
    else:
        print('{} {} to {} {} correctly converted (no given quantities)'.format(test_dict[idx]['orig_val'],test_dict[idx]['orig_unit'],test_dict[idx]['new_val'],test_dict[idx]['new_unit']))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Valence dependent unit conversion (e.g. eq l-1 --> mol l-1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

print('\nVALENCE DEPENDENT UNIT CONVERSION TESTS')

orig_unit = 'eq l-1'
orig_val = 1.0
new_unit = 'mol l-1'
new_val = 1.0
conv_obj = unit_converter.convert_units(orig_unit, new_unit, orig_val, species='NO3-', valence=1)
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'eq l-1'
orig_val = 1.0
new_unit = 'mol l-1'
new_val = 0.5
conv_obj = unit_converter.convert_units(orig_unit, new_unit, orig_val, species='SO4--')
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Species dependent unit conversion (e.g. mgC/m3 --> mg/m3)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

print('\nSPECIES DEPENDENT UNIT CONVERSION TESTS')

orig_unit = 'mg/m3'
orig_val = 294.51
new_unit = 'mgC/m3'
new_val = 240.0
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='C3H8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'mgC/m3'
orig_val = 240.0
new_unit = 'mg/m3'
new_val = 294.51
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='C3H8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppbv'
orig_val = 10.0
new_unit = 'umol/mol'
new_val = 0.01
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='C3H8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ug/m3'
orig_val = 13.0
new_unit = 'mg m-3'
new_val = 0.013
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='C3H8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmvC'
orig_val = 3
new_unit = 'ppmv'
new_val = 1
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='C3H8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmvC'
orig_val = 1
new_unit = 'umol mol-1'
new_val = 1
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='CH8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 1
new_unit = 'ppmvC'
new_val = 3
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='C3H8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'umolC mol-1'
orig_val = 3
new_unit = 'ppmv'
new_val = 1
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='C3H8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'umolC mol-1'
orig_val = 1
new_unit = 'umol mol-1'
new_val = 1
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='CH8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 1
new_unit = 'umolC mol-1'
new_val = 3
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='C3H8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'eq l-1'
orig_val = 1.0
new_unit = 'mol l-1'
new_val = 0.5
conv_obj = unit_converter.convert_units(orig_unit, new_unit, orig_val, species='SO4', valence=2)
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'mgCaCO3 l-1'
orig_val = 15.0
new_unit = 'mg l-1'
new_val = 15.0
conv_obj = unit_converter.convert_units(orig_unit,new_unit,orig_val,species='CaCO3', valence=2) 
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Different unit quantity conversion (e.g. ug m-3 --> ppbv)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

print('\nDIFFERENT QUANTITY UNIT CONVERSION TESTS')

orig_unit = 'ug/m3'
orig_val = 1.9957
new_unit = 'ppbv'
new_val = 1.0
input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mass_concentration':orig_unit}
input_values = {'temperature':293.15, 'pressure':1013.25, 'molar_mass':47.998, 'mass_concentration':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='mass_concentration', output_quantity='mole_fraction')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppbv'
orig_val = 1.0
new_unit = 'ug/m3'
new_val = 1.9957
input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':293.15, 'pressure':1013.25, 'molar_mass':47.998, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='mole_fraction', output_quantity='mass_concentration')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 150.0
new_unit = 'mg/m3'
new_val = 294.51
input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.0, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='mole_fraction', output_quantity='mass_concentration')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 150.0
new_unit = 'mgC/m3'
new_val = 240.0
input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.0, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='mole_fraction', output_quantity='mass_concentration', species='C3H8')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 450.0
new_unit = 'mg/m3'
new_val = 884.0
input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.01, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='mole_fraction', output_quantity='mass_concentration')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'ppmv'
orig_val = 450.0
new_unit = 'mg/L'
new_val = 0.884
input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.01, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='mole_fraction', output_quantity='mass_concentration')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'mg/m3'
orig_val = 884.0
new_unit = 'ppmv'
new_val = 450.0
input_units = {'temperature':'K', 'pressure':'hPa', 'molar_mass':'g mol-1', 'mass_concentration':orig_unit}
input_values = {'temperature':273.15, 'pressure':1013.25, 'molar_mass':44.01, 'mass_concentration':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='mass_concentration', output_quantity='mole_fraction')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'mol/l'
orig_val = 15.0
new_unit = 'mg/l'
new_val = 660150.0
input_units = {'molar_mass':'g mol-1', 'molarity':orig_unit}
input_values = {'molar_mass':44.01, 'molarity':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='molarity', output_quantity='mass_concentration')
if np.isclose(conv_obj.converted_value,new_val,atol=1.0) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'eq l-1'
orig_val = 1.0
new_unit = 'mg l-1'
new_val = 48035.0
input_units = {'molar_mass':'g mol-1', 'molarity':orig_unit}
input_values = {'molar_mass':96.07, 'molarity':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='molarity', output_quantity='mass_concentration', species='SO4--')
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'kg kg-1'
orig_val = 5.7981605e-08
new_unit = 'nmol mol-1'
new_val = 34.98961
input_units = {'molar_mass':'g mol-1', 'mass_fraction':orig_unit}
input_values = {'molar_mass':47.998, 'mass_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='mass_fraction', output_quantity='mole_fraction')
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'nmol mol-1'
orig_val = 34.98961
new_unit = 'kg kg-1'
new_val = 5.7981605e-08
input_units = {'molar_mass':'g mol-1', 'mole_fraction':orig_unit}
input_values = {'molar_mass':47.998, 'mole_fraction':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='mole_fraction', output_quantity='mass_fraction')
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'm-1 sr-1'
orig_val = 1.0
new_unit = 'Mm-1'
new_val = 6.283185e6
input_units = {'volume_scattering_function':orig_unit}
input_values = {'volume_scattering_function':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='volume_scattering_function', output_quantity='backscattering_coefficient')
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'm-1 sr-1'
orig_val = 1.0
new_unit = 'Mm-1'
new_val = 1.2566370e7
input_units = {'volume_scattering_function':orig_unit}
input_values = {'volume_scattering_function':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='volume_scattering_function', output_quantity='scattering_coefficient')
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'Mm-1'
orig_val = 6.283185e6
new_unit = 'm-1 sr-1'
new_val = 1.0
input_units = {'backscattering_coefficient':orig_unit}
input_values = {'backscattering_coefficient':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='backscattering_coefficient', output_quantity='volume_scattering_function')
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))

orig_unit = 'Mm-1'
orig_val = 1.2566370e7
new_unit = 'm-1 sr-1'
new_val = 1.0
input_units = {'scattering_coefficient':orig_unit}
input_values = {'scattering_coefficient':orig_val}
conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, input_quantity='scattering_coefficient', output_quantity='volume_scattering_function')
if np.isclose(conv_obj.converted_value,new_val) == False:
    sys.exit('{} {} to {} should be {}, is {}'.format(orig_val,orig_unit,new_unit,new_val,conv_obj.converted_value))
else:
    print('{} {} to {} {} correctly converted'.format(orig_val,orig_unit,new_val,new_unit))