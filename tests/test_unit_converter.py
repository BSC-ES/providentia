import numpy as np
import pytest

from .aux_functions import check_unit_conversion
import providentia.unit_converter as unit_converter

possibilities = [
    ('g', 500, 'kg', 0.5, 'mass', 'mass'),
    ('K', 3, 'degC', -270.15, 'temperature', 'temperature'),
    ('degF', 50, 'K', 283.15, 'temperature', 'temperature'),
    ('m/s2', 3000, 'km/s2', 3, 'acceleration', 'acceleration'),
    ('kg/m3', 18, 'pg/um3', 0.018, 'mass_concentration', 'mass_concentration'),
    ('W/J', 25, 'mW/kJ', 2.5e7, 'frequency', 'frequency'),
    ('inches', 3, 'cm', 7.62, 'length', 'length'),
    ('feet', 5, 'm', 1.524, 'length', 'length'),
    ('s-1', 1000, 'ms-1', 1, 'radioactivity', 'radioactivity'),
    ('Hz', 3, 'MHz', 0.000003, 'frequency', 'frequency'),
    ('m-1', 25, 'um-1', 0.000025, 'backscattering_coefficient', 'backscattering_coefficient'),
    ('umol mol-1', 100.0, 'mmol mol-1', 0.1, 'mole_fraction', 'mole_fraction'),
    ('kg kg-1', 1.0, 'g kg-1', 1000.0, 'mass_fraction', 'mass_fraction'),
    ('m-1', 6.0, 'um-1', 6e-6, 'absorption_coefficient', 'absorption_coefficient'),
    ('um-1', 6e-6, 'm-1', 6.0, 'scattering_coefficient', 'scattering_coefficient'),
    ('m2 kg-1', 15.0, 'km2 g-1', 1.5e-8, 'mass_absorption_cross_section', 'mass_absorption_cross_section'),
    ('m m-1', 22.0, 'mm mm-1', 22.0, 'single_scattering_albedo', 'single_scattering_albedo'),
    ('unitless', 22.0, 'unitless', 22.0, 'single_scattering_albedo', 'single_scattering_albedo'),
    ('m m-1', 22.0, 'm km-1', 2.2e4, 'single_scattering_albedo', 'single_scattering_albedo'),
    ('mol m-3 mol-1 m3', 7.0, 'nmol mm-3 nmol-1 mm3', 7.0, 'pH', 'pH'),
    ('unitless', 7.0, 'unitless', 7.0, 'pH', 'pH'),
    ('kg m-2 s-1', 13.0, 'g mm-2 hour-1', 46.8, 'mass_flux', 'mass_flux'),
    ('kg m-1 s-2 kg-1 m s2', 84.0, 'g mm-1 min-2 g-1 mm min2', 84.0, 'relative_humidity', 'relative_humidity'),
    ('unitless', 84.0, 'unitless', 84.0, 'relative_humidity', 'relative_humidity'),
    ('m m-1', 64.0, 'um um-1', 64.0, 'optical_depth', 'optical_depth'),
    ('unitless', 64.0, 'unitless', 64.0, 'optical_depth', 'optical_depth'),
    ('m m-1', 64.0, 'um m-1', 6.4e7, 'optical_depth', 'optical_depth'),
    ('um hours-1 um-1 hours', 28.0, 'm s-1 m-1 s', 28.0, 'refractive_index', 'refractive_index'),
    ('unitless', 28.0, 'unitless', 28.0, 'refractive_index', 'refractive_index'),
    ('m3 m-2', 70.0, 'mm3 km-2', 7e16, 'volume_concentration', 'volume_concentration'),
    ('hours', 24.0, 'minutes', 1440.0, 'time', 'time'),
    ('angular degrees', 270.0, 'radians', 4.7124, 'angle', 'angle'),
    ('m s-1', 3.0, 'um min-1', 1.8e8, 'velocity', 'velocity'),
    ('m-2', 25.0, 'km-2', 2.5e7, 'population_density', 'population_density'),
    ('km-2', 2.0, 'nm-2', 2e-24, 'column_number_density', 'column_number_density'),
    ('m3 s-1', 21.0, 'nm3 min-1', 1.26e30, 'volumetric_flow', 'volumetric_flow'),
    ('m2', 17.0, 'mm2', 1.7e7, 'area', 'area'),
    ('unitless', 1.0, 'unitless', 1.0, 'unitless', 'unitless')
]
@pytest.mark.parametrize("orig_unit, orig_val, new_unit, new_val, input_quantity, output_quantity", possibilities)
def test_same_quantity_unit_conversion(orig_unit, orig_val, new_unit, new_val, input_quantity, output_quantity):
    """Same quantity unit conversion."""
    
    conv_obj = unit_converter.convert_units(orig_unit, new_unit, orig_val, input_quantity=input_quantity, output_quantity=output_quantity)
    check_unit_conversion(conv_obj, new_val, orig_val, orig_unit, new_unit)
    
    conv_obj = unit_converter.convert_units(orig_unit, new_unit, orig_val)
    check_unit_conversion(conv_obj, new_val, orig_val, orig_unit, new_unit)


possibilities = [
    # Valence dependent unit conversion
    ('eq l-1', 1.0, 'mol l-1', 1.0, 'NO3-', 1), 
    ('eq l-1', 1.0, 'mol l-1', 0.5, 'SO4--', np.nan), 
    # Species dependent unit conversion
    ('mg/m3', 294.51, 'mgC/m3', 240.0, 'C3H8', np.nan), 
    ('mgC/m3', 240.0, 'mg/m3', 294.51, 'C3H8', np.nan), 
    ('ppbv', 10.0, 'umol/mol', 0.01, 'C3H8', np.nan), 
    ('ug/m3', 13.0, 'mg m-3', 0.013, 'C3H8', np.nan), 
    ('ppmvC', 3, 'ppmv', 1, 'C3H8', np.nan), 
    ('ppmvC', 1, 'umol mol-1', 1, 'CH8', np.nan), 
    ('ppmv', 1, 'ppmvC', 3, 'C3H8', np.nan), 
    ('umolC mol-1', 3, 'ppmv', 1, 'C3H8', np.nan), 
    ('umolC mol-1', 1, 'umol mol-1', 1, 'CH8', np.nan),
    ('ppmv', 1, 'umolC mol-1', 3, 'C3H8', np.nan),
    ('eq l-1', 1.0, 'mol l-1', 0.5, 'SO4', 2),
    ('mgCaCO3 l-1', 15.0, 'mg l-1', 15.0, 'CaCO3', 2)
]
@pytest.mark.parametrize("orig_unit, orig_val, new_unit, new_val, species, valence", possibilities)
def test_valence_species_unit_conversion(orig_unit, orig_val, new_unit, new_val, species, valence):
    """Valence and species dependent unit conversions."""

    conv_obj = unit_converter.convert_units(orig_unit, new_unit, orig_val, species=species, valence=valence)
    check_unit_conversion(conv_obj, new_val, orig_val, orig_unit, new_unit, atol=1.0)


possibilities = [
    ('ug/m3', 1.9957, 'ppbv', 1.0, '', 
     {'temperature': 'K', 'pressure': 'hPa', 'molar_mass': 'g mol-1', 'mass_concentration': 'ug/m3'},
     {'temperature': 293.15, 'pressure': 1013.25, 'molar_mass': 47.998, 'mass_concentration': 1.9957},
     'mass_concentration', 'mole_fraction'),
    ('ppbv', 1.0, 'ug/m3', 1.9957, '',
     {'temperature': 'K', 'pressure': 'hPa', 'molar_mass': 'g mol-1', 'mole_fraction': 'ppbv'},
     {'temperature': 293.15, 'pressure': 1013.25, 'molar_mass': 47.998, 'mole_fraction': 1.0},
     'mole_fraction', 'mass_concentration'),
    ('ppmv', 150.0, 'mg/m3', 294.51, '',
     {'temperature': 'K', 'pressure': 'hPa', 'molar_mass': 'g mol-1', 'mole_fraction': 'ppmv'},
     {'temperature': 273.15, 'pressure': 1013.25, 'molar_mass': 44.0, 'mole_fraction': 150.0},
     'mole_fraction', 'mass_concentration'),
    ('ppmv', 150.0, 'mgC/m3', 240.0, 'C3H8',
     {'temperature': 'K', 'pressure': 'hPa', 'molar_mass': 'g mol-1', 'mole_fraction': 'ppmv'},
     {'temperature': 273.15, 'pressure': 1013.25, 'molar_mass': 44.0, 'mole_fraction': 150.0},
     'mole_fraction', 'mass_concentration'),
    ('ppmv', 450.0, 'mg/m3', 884.0, '',
     {'temperature': 'K', 'pressure': 'hPa', 'molar_mass': 'g mol-1', 'mole_fraction': 'ppmv'},
     {'temperature': 273.15, 'pressure': 1013.25, 'molar_mass': 44.01, 'mole_fraction': 450.0},
     'mole_fraction', 'mass_concentration'),
    ('ppmv', 450.0, 'mg/L', 0.884, '',
     {'temperature': 'K', 'pressure': 'hPa', 'molar_mass': 'g mol-1', 'mole_fraction': 'ppmv'},
     {'temperature': 273.15, 'pressure': 1013.25, 'molar_mass': 44.01, 'mole_fraction': 450.0},
     'mole_fraction', 'mass_concentration'),
    ('mg/m3', 884.0, 'ppmv', 450.0, '',
     {'temperature': 'K', 'pressure': 'hPa', 'molar_mass': 'g mol-1', 'mass_concentration': 'mg/m3'}, 
     {'temperature': 273.15, 'pressure': 1013.25, 'molar_mass': 44.01, 'mass_concentration': 884.0},
     'mass_concentration', 'mole_fraction'),
    ('mol/l', 15.0, 'mg/l', 660150.0, '', 
     {'molar_mass': 'g mol-1', 'molarity': 'mol/l'},
     {'molar_mass': 44.01, 'molarity': 15.0},
     'molarity', 'mass_concentration'),
    ('eq l-1', 1.0, 'mg l-1', 48035.0, 'SO4--',
     {'molar_mass': 'g mol-1', 'molarity': 'eq l-1'}, 
     {'molar_mass': 96.07, 'molarity': 1.0},
     'molarity', 'mass_concentration'),
    ('kg kg-1', 5.7981605e-08, 'nmol mol-1', 34.98961, '',
     {'molar_mass': 'g mol-1', 'mass_fraction': 'kg kg-1'}, 
     {'molar_mass': 47.998, 'mass_fraction': 5.7981605e-08},
     'mass_fraction', 'mole_fraction'),
    ('nmol mol-1', 34.98961, 'kg kg-1', 5.7981605e-08, '',
     {'molar_mass': 'g mol-1', 'mole_fraction': 'nmol mol-1'},
     {'molar_mass': 47.998, 'mole_fraction': 34.98961},
     'mole_fraction', 'mass_fraction'),
    ('m-1 sr-1', 1.0, 'Mm-1', 6.283185e6, '',
     {'volume_scattering_function': 'm-1 sr-1'}, 
     {'volume_scattering_function': 1.0}, 
     'volume_scattering_function', 'backscattering_coefficient'),
    ('m-1 sr-1', 1.0, 'Mm-1', 1.2566370e7, '',
     {'volume_scattering_function': 'm-1 sr-1'},
     {'volume_scattering_function': 1.0},
     'volume_scattering_function', 'scattering_coefficient'),
    ('Mm-1',  6.283185e6, 'm-1 sr-1', 1.0, '',
     {'backscattering_coefficient': 'Mm-1'},
     {'backscattering_coefficient':  6.283185e6},
     'backscattering_coefficient', 'volume_scattering_function'),
    ('Mm-1', 1.2566370e7, 'm-1 sr-1', 1.0, '',
     {'scattering_coefficient': 'Mm-1'}, 
     {'scattering_coefficient':  1.2566370e7},
     'scattering_coefficient', 'volume_scattering_function')
]
@pytest.mark.parametrize("orig_unit, orig_val, new_unit, new_val, species, input_units, input_values, input_quantity, output_quantity", 
                         possibilities)
def test_different_quantity_unit_conversion(orig_unit, orig_val, new_unit, new_val, species, input_units,
                                            input_values, input_quantity, output_quantity):
    """Different quantity unit conversion"""

    conv_obj = unit_converter.convert_units(input_units, new_unit, input_values, 
                                            input_quantity=input_quantity, 
                                            output_quantity=output_quantity,
                                            species=species)
    check_unit_conversion(conv_obj, new_val, orig_val, orig_unit, new_unit, atol=1.0)