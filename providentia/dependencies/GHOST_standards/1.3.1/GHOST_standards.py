import numpy as np

###--------------------------------------------------------------------------------------------------###
###MEASURED PARAMETER STANDARDISATIONS
###--------------------------------------------------------------------------------------------------###

#Define dictionary with details pertaining to standardised measured parameters
 
standard_parameters ={  
'O3':                  {'long_parameter_name':'ozone',                                                                                                                            'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'44201', 'airbase_name':'Ozone (air)',                          'airbase_code':'00007', 'ebas_parameter_name':'ozone',                        'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'007', 'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'O3',      'bsc_parameter_name':'sconco3'         },
'NO':                  {'long_parameter_name':'nitrogen monoxide',                                                                                                                'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':5000.0,   'extreme_upper_monthly_median':250.0,  'aqs_code':'42601', 'airbase_name':'Nitrogen monoxide (air)',              'airbase_code':'00038', 'ebas_parameter_name':'nitrogen_monoxide',            'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'008', 'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NO',      'bsc_parameter_name':'sconcno'         },
'NO2':                 {'long_parameter_name':'nitrogen dioxide',                                                                                                                 'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':5000.0,   'extreme_upper_monthly_median':250.0,  'aqs_code':'42602', 'airbase_name':'Nitrogen dioxide (air)',               'airbase_code':'00008', 'ebas_parameter_name':'nitrogen_dioxide',             'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'006', 'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NO2',     'bsc_parameter_name':'sconcno2'        },
'CO':                  {'long_parameter_name':'carbon monoxide',                                                                                                                  'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':10000.0,  'extreme_upper_monthly_median':3000.0, 'aqs_code':'42101', 'airbase_name':'Carbon monoxide (air)',                'airbase_code':'00010', 'ebas_parameter_name':'carbon_monoxide',              'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'005', 'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'CO',      'bsc_parameter_name':'sconcco'         },
'CH4':                 {'long_parameter_name':'methane',                                                                                                                          'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':50000.0,  'extreme_upper_monthly_median':5000.0, 'aqs_code':'',      'airbase_name':'Methane (air)',                        'airbase_code':'00041', 'ebas_parameter_name':'methane',                      'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'CH4',     'bsc_parameter_name':'sconcch4'        },
'C6H6':                {'long_parameter_name':'benzene',                                                                                                                          'matrix':'gas',     'standard_units':'pmol mol-1',      'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':10000.0,  'extreme_upper_monthly_median':1000.0, 'aqs_code':'',      'airbase_name':'Benzene (air)',                        'airbase_code':'00020', 'ebas_parameter_name':'benzene',                      'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C6H6',    'bsc_parameter_name':'sconcc6h6'       },
'ISOPRENE':            {'long_parameter_name':'isoprene',                                                                                                                         'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':0.1,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'43243', 'airbase_name':'Isoprene (air)',                       'airbase_code':'00451', 'ebas_parameter_name':'isoprene',                     'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C5H8',    'bsc_parameter_name':'sconcisop'       },
'SO2':                 {'long_parameter_name':'sulphur dioxide',                                                                                                                  'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':5000.0,   'extreme_upper_monthly_median':300.0,  'aqs_code':'',      'airbase_name':'Sulphur dioxide (air)',                'airbase_code':'00001', 'ebas_parameter_name':'sulphur_dioxide',              'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'004', 'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO2',     'bsc_parameter_name':'sconcso2'        },
'NH3':                 {'long_parameter_name':'ammonia',                                                                                                                          'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':3000.0,   'extreme_upper_monthly_median':300.0,  'aqs_code':'',      'airbase_name':'Ammonia (air)',                        'airbase_code':'00035', 'ebas_parameter_name':'ammonia',                      'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH3',     'bsc_parameter_name':'sconcnh3'        },
'HNO3':                {'long_parameter_name':'nitric acid',                                                                                                                      'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Nitric acid (air)',                    'airbase_code':'00050', 'ebas_parameter_name':'nitric_acid',                  'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'HNO3',    'bsc_parameter_name':'sconchno3'       },
'PAN':                 {'long_parameter_name':'peroxyacetyl nitrate',                                                                                                             'matrix':'gas',     'standard_units':'nmol mol-1',      'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':200.0,    'extreme_upper_monthly_median':20.0,   'aqs_code':'',      'airbase_name':'Peroxyacetyl nitrate (air)',           'airbase_code':'00034', 'ebas_parameter_name':'',                             'ebas_matrix':'air',     'ebas_priority_units':['mol/mol','/m3'], 'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C2H3NO5', 'bsc_parameter_name':'sconcpan'        },
'PM10':                {'long_parameter_name':'PM10 mass',                                                                                                                        'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':10000.0,  'extreme_upper_monthly_median':1000.0, 'aqs_code':'',      'airbase_name':'Particulate matter < 10 µm (aerosol)', 'airbase_code':'00005', 'ebas_parameter_name':'pm10_mass',                    'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'pm10'            },
'PM2.5':               {'long_parameter_name':'PM2.5 mass',                                                                                                                       'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':10000.0,  'extreme_upper_monthly_median':1000.0, 'aqs_code':'',      'airbase_name':'Particulate matter < 2.5 m (aerosol)', 'airbase_code':'06001', 'ebas_parameter_name':'pm25_mass',                    'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'pm2p5'           },
'PM1':                 {'long_parameter_name':'PM1 mass',                                                                                                                         'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':10000.0,  'extreme_upper_monthly_median':1000.0, 'aqs_code':'',      'airbase_name':'Particulate matter < 1 µm (aerosol)',  'airbase_code':'06002', 'ebas_parameter_name':'pm1_mass',                     'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'pm1'             },
'Al':                  {'long_parameter_name':'aluminium',                                                                                                                        'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'aluminium (aerosol)',                  'airbase_code':'00604', 'ebas_parameter_name':'aluminium',                    'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Al',      'bsc_parameter_name':'sconcal'         },
'As':                  {'long_parameter_name':'arsenic',                                                                                                                          'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Arsenic (aerosol)',                    'airbase_code':'00018', 'ebas_parameter_name':'arsenic',                      'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'As',      'bsc_parameter_name':'sconcas'         },
'BC':                  {'long_parameter_name':'black carbon',                                                                                                                     'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'black_carbon (aerosol)',               'airbase_code':'00391', 'ebas_parameter_name':'',                             'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'sconcbc'         },
'eBC':                 {'long_parameter_name':'equivalent black carbon',                                                                                                          'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'equivalent_black_carbon',      'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'sconcebc'        },
'C':                   {'long_parameter_name':'carbon',                                                                                                                           'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'total_carbon',                 'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'sconcc'          },
'Ca++':                {'long_parameter_name':'calcium',                                                                                                                          'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'calcium (aerosol)',                    'airbase_code':'00629', 'ebas_parameter_name':'calcium',                      'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ca',      'bsc_parameter_name':'sconcca'         },
'Cd':                  {'long_parameter_name':'cadmium',                                                                                                                          'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Cadmium (aerosol)',                    'airbase_code':'00014', 'ebas_parameter_name':'cadmium',                      'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cd',      'bsc_parameter_name':'sconccd'         },
'Cl-':                 {'long_parameter_name':'chloride',                                                                                                                         'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'chloride (aerosol)',                   'airbase_code':'00631', 'ebas_parameter_name':'chloride',                     'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cl',      'bsc_parameter_name':'sconccl'         },
'Co':                  {'long_parameter_name':'cobalt',                                                                                                                           'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Cobalt (aerosol)',                     'airbase_code':'00064', 'ebas_parameter_name':'cobalt',                       'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Co',      'bsc_parameter_name':'sconccobalt'     },
'Cr':                  {'long_parameter_name':'chromium',                                                                                                                         'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Chromium (aerosol)',                   'airbase_code':'00016', 'ebas_parameter_name':'chromium',                     'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cr',      'bsc_parameter_name':'sconccr'         },
'Cu':                  {'long_parameter_name':'copper',                                                                                                                           'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Copper (aerosol)',                     'airbase_code':'00073', 'ebas_parameter_name':'copper',                       'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cu',      'bsc_parameter_name':'sconccu'         },
'DUST':                {'long_parameter_name':'dust',                                                                                                                             'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'sconcdu'         },
'EC':                  {'long_parameter_name':'elemental carbon',                                                                                                                 'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Elemental carbon (aerosol)',           'airbase_code':'00771', 'ebas_parameter_name':'elemental_carbon',             'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'sconcec'         },
'Fe':                  {'long_parameter_name':'iron',                                                                                                                             'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'Iron (aerosol)',                       'airbase_code':'00065', 'ebas_parameter_name':'iron',                         'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Fe',      'bsc_parameter_name':'sconcfe'         },
'Hg':                  {'long_parameter_name':'mercury',                                                                                                                          'matrix':'aerosol', 'standard_units':'pg m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Mercury (aerosol)',                    'airbase_code':'00013', 'ebas_parameter_name':'mercury',                      'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Hg',      'bsc_parameter_name':'sconchg'         },
'K+':                  {'long_parameter_name':'potassium',                                                                                                                        'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'potassium (aerosol)',                  'airbase_code':'00657', 'ebas_parameter_name':'potassium',                    'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'K',       'bsc_parameter_name':'sconck'          },
'Mg++':                {'long_parameter_name':'magnesium',                                                                                                                        'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'magnesium (aerosol)',                  'airbase_code':'00659', 'ebas_parameter_name':'magnesium',                    'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mg',      'bsc_parameter_name':'sconcmg'         },
'Mn':                  {'long_parameter_name':'manganese',                                                                                                                        'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Manganese (aerosol)',                  'airbase_code':'00017', 'ebas_parameter_name':'manganese',                    'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mn',      'bsc_parameter_name':'sconcmn'         },
'MSA':                 {'long_parameter_name':'methanesulfonic acid',                                                                                                             'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'methanesulfonic_acid',         'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'CH4O3S',  'bsc_parameter_name':'sconcmsa'        },
'Na+':                 {'long_parameter_name':'sodium',                                                                                                                           'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'sodium (aerosol)',                     'airbase_code':'00668', 'ebas_parameter_name':'sodium',                       'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Na',      'bsc_parameter_name':'sconcna'         },
'NH4+':                {'long_parameter_name':'ammonium',                                                                                                                         'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Particulate ammonium (aerosol)',       'airbase_code':'00045', 'ebas_parameter_name':'ammonium',                     'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4',     'bsc_parameter_name':'sconcnh4'        },
'NH4NO3':              {'long_parameter_name':'ammonium nitrate',                                                                                                                 'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'ammonium_nitrate',             'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4NO3',  'bsc_parameter_name':'sconcnh4no3'     },
'Ni':                  {'long_parameter_name':'nickel',                                                                                                                           'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Nickel (aerosol)',                     'airbase_code':'00015', 'ebas_parameter_name':'nickel',                       'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ni',      'bsc_parameter_name':'sconcni'         },
'NO3-':                {'long_parameter_name':'nitrate',                                                                                                                          'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Particulate nitrate (aerosol)',        'airbase_code':'00046', 'ebas_parameter_name':'nitrate',                      'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NO3',     'bsc_parameter_name':'sconcno3'        },
'OC':                  {'long_parameter_name':'organic carbon',                                                                                                                   'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Organic carbon (aerosol)',             'airbase_code':'00772', 'ebas_parameter_name':'organic_carbon',               'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'sconcoc'         },
'Pb':                  {'long_parameter_name':'lead',                                                                                                                             'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Lead (aerosol)',                       'airbase_code':'00012', 'ebas_parameter_name':'lead',                         'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Pb',      'bsc_parameter_name':'sconcpb'         },
'Se':                  {'long_parameter_name':'selenium',                                                                                                                         'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Selenium (aerosol)',                   'airbase_code':'00048', 'ebas_parameter_name':'selenium',                     'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Se',      'bsc_parameter_name':'sconcse'         },
'SO4--':               {'long_parameter_name':'sulphate',                                                                                                                         'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Particulate sulphate (aerosol)',       'airbase_code':'00047', 'ebas_parameter_name':'sulphate_total',               'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'sconcso4'        },
'SO4--_NSS':           {'long_parameter_name':'sulphate: non-sea salt',                                                                                                           'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_corrected',           'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'sconcso4nss'     },
'SO4--_SS':            {'long_parameter_name':'sulphate: sea salt',                                                                                                               'matrix':'aerosol', 'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_seasalt',             'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'sconcso4ss'      },
'V':                   {'long_parameter_name':'vanadium',                                                                                                                         'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Vanadium (aerosol)',                   'airbase_code':'00049', 'ebas_parameter_name':'vanadium',                     'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'V',       'bsc_parameter_name':'sconcv'          },
'Zn':                  {'long_parameter_name':'zinc',                                                                                                                             'matrix':'aerosol', 'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':5000.0,   'extreme_upper_monthly_median':500.0,  'aqs_code':'',      'airbase_name':'Zinc (aerosol)',                       'airbase_code':'00063', 'ebas_parameter_name':'zinc',                         'ebas_matrix':'aerosol', 'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Zn',      'bsc_parameter_name':'sconczn'         },
'PM10_Al':             {'long_parameter_name':'PM10 aluminium',                                                                                                                   'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'aluminium in PM10 (aerosol)',          'airbase_code':'05604', 'ebas_parameter_name':'aluminium',                    'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Al',      'bsc_parameter_name':'pm10al'          },
'PM10_As':             {'long_parameter_name':'PM10 arsenic',                                                                                                                     'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Arsenic in PM10 (aerosol)',            'airbase_code':'05018', 'ebas_parameter_name':'arsenic',                      'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'As',      'bsc_parameter_name':'pm10as'          },
'PM10_BC':             {'long_parameter_name':'PM10 black carbon',                                                                                                                'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'black_carbon in PM10 (aerosol)',       'airbase_code':'05391', 'ebas_parameter_name':'',                             'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm10bc'          },
'PM10_eBC':            {'long_parameter_name':'PM10 equivalent black carbon',                                                                                                     'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'equivalent_black_carbon',      'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm10ebc'         },
'PM10_C':              {'long_parameter_name':'PM10 carbon',                                                                                                                      'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'total_carbon',                 'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm10c'           },
'PM10_Ca++':           {'long_parameter_name':'PM10 calcium',                                                                                                                     'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'calcium in PM10 (aerosol)',            'airbase_code':'05629', 'ebas_parameter_name':'calcium',                      'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ca',      'bsc_parameter_name':'pm10ca'          },
'PM10_Cd':             {'long_parameter_name':'PM10 cadmium',                                                                                                                     'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Cadmium in PM10 (aerosol)',            'airbase_code':'05014', 'ebas_parameter_name':'cadmium',                      'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cd',      'bsc_parameter_name':'pm10cd'          },
'PM10_Cl-':            {'long_parameter_name':'PM10 chloride',                                                                                                                    'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'chloride in PM10 (aerosol)',           'airbase_code':'05631', 'ebas_parameter_name':'chloride',                     'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cl',      'bsc_parameter_name':'pm10cl'          },
'PM10_Co':             {'long_parameter_name':'PM10 cobalt',                                                                                                                      'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Cobalt in PM10 (aerosol)',             'airbase_code':'05064', 'ebas_parameter_name':'cobalt',                       'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Co',      'bsc_parameter_name':'pm10cobalt'      },
'PM10_Cr':             {'long_parameter_name':'PM10 chromium',                                                                                                                    'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Chromium in PM10 (aerosol)',           'airbase_code':'05016', 'ebas_parameter_name':'chromium',                     'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cr',      'bsc_parameter_name':'pm10cr'          },
'PM10_Cu':             {'long_parameter_name':'PM10 copper',                                                                                                                      'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Copper in PM10 (aerosol)',             'airbase_code':'05073', 'ebas_parameter_name':'copper',                       'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cu',      'bsc_parameter_name':'pm10cu'          },
'PM10_DUST':           {'long_parameter_name':'PM10 dust',                                                                                                                        'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'pm10du'          },
'PM10_EC':             {'long_parameter_name':'PM10 elemental carbon',                                                                                                            'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Elemental carbon in PM10 (aerosol)',   'airbase_code':'05771', 'ebas_parameter_name':'elemental_carbon',             'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm10ec'          },
'PM10_Fe':             {'long_parameter_name':'PM10 iron',                                                                                                                        'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'Iron in PM10 (aerosol)',               'airbase_code':'05065', 'ebas_parameter_name':'iron',                         'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Fe',      'bsc_parameter_name':'pm10fe'          },
'PM10_Hg':             {'long_parameter_name':'PM10 mercury',                                                                                                                     'matrix':'pm10',    'standard_units':'pg m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Mercury in PM10 (aerosol)',            'airbase_code':'05013', 'ebas_parameter_name':'mercury',                      'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Hg',      'bsc_parameter_name':'pm10hg'          },
'PM10_K+':             {'long_parameter_name':'PM10 potassium',                                                                                                                   'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'potassium in PM10 (aerosol)',          'airbase_code':'05657', 'ebas_parameter_name':'potassium',                    'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'K',       'bsc_parameter_name':'pm10k'           },
'PM10_Mg++':           {'long_parameter_name':'PM10 magnesium',                                                                                                                   'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'magnesium in PM10 (aerosol)',          'airbase_code':'05659', 'ebas_parameter_name':'magnesium',                    'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mg',      'bsc_parameter_name':'pm10mg'          },
'PM10_Mn':             {'long_parameter_name':'PM10 manganese',                                                                                                                   'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Manganese in PM10 (aerosol)',          'airbase_code':'05017', 'ebas_parameter_name':'manganese',                    'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mn',      'bsc_parameter_name':'pm10mn'          },
'PM10_MSA':            {'long_parameter_name':'PM10 methanesulfonic acid',                                                                                                        'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'methanesulfonic_acid',         'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'CH4O3S',  'bsc_parameter_name':'pm10msa'         },
'PM10_Na+':            {'long_parameter_name':'PM10 sodium',                                                                                                                      'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'sodium in PM10 (aerosol)',             'airbase_code':'05668', 'ebas_parameter_name':'sodium',                       'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Na',      'bsc_parameter_name':'pm10na'          },
'PM10_NH4+':           {'long_parameter_name':'PM10 ammonium',                                                                                                                    'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Ammonium in PM10 (aerosol)',           'airbase_code':'05045', 'ebas_parameter_name':'ammonium',                     'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4',     'bsc_parameter_name':'pm10nh4'         },
'PM10_NH4NO3':         {'long_parameter_name':'PM10 ammonium nitrate',                                                                                                            'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'ammonium_nitrate',             'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4NO3',  'bsc_parameter_name':'pm10nh4no3'      },
'PM10_Ni':             {'long_parameter_name':'PM10 nickel',                                                                                                                      'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Nickel in PM10 (aerosol)',             'airbase_code':'05015', 'ebas_parameter_name':'nickel',                       'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ni',      'bsc_parameter_name':'pm10ni'          },
'PM10_NO3-':           {'long_parameter_name':'PM10 nitrate',                                                                                                                     'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Nitrate in PM10 (aerosol)',            'airbase_code':'05046', 'ebas_parameter_name':'nitrate',                      'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NO3',     'bsc_parameter_name':'pm10no3'         },
'PM10_OC':             {'long_parameter_name':'PM10 organic carbon',                                                                                                              'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Organic carbon in PM10 (aerosol)',     'airbase_code':'05772', 'ebas_parameter_name':'organic_carbon',               'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm10oc'          },
'PM10_Pb':             {'long_parameter_name':'PM10 lead',                                                                                                                        'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Lead in PM10 (aerosol)',               'airbase_code':'05012', 'ebas_parameter_name':'lead',                         'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Pb',      'bsc_parameter_name':'pm10pb'          },
'PM10_Se':             {'long_parameter_name':'PM10 selenium',                                                                                                                    'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Selenium in PM10 (aerosol)',           'airbase_code':'05048', 'ebas_parameter_name':'selenium',                     'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Se',      'bsc_parameter_name':'pm10se'          },
'PM10_SO4--':          {'long_parameter_name':'PM10 sulphate',                                                                                                                    'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'sulphate in PM10 (aerosol)',           'airbase_code':'05047', 'ebas_parameter_name':'sulphate_total',               'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'pm10so4'         },
'PM10_SO4--_NSS':      {'long_parameter_name':'PM10 sulphate : non-sea salt',                                                                                                     'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_corrected',           'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'pm10so4nss'      },
'PM10_SO4--_SS':       {'long_parameter_name':'PM10 sulphate : sea salt',                                                                                                         'matrix':'pm10',    'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_seasalt',             'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'pm10so4ss'       },
'PM10_V':              {'long_parameter_name':'PM10 vanadium',                                                                                                                    'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Vanedium in PM10 (aerosol)',           'airbase_code':'05049', 'ebas_parameter_name':'vanadium',                     'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'V',       'bsc_parameter_name':'pm10v'           },
'PM10_Zn':             {'long_parameter_name':'PM10 zinc',                                                                                                                        'matrix':'pm10',    'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':5000.0,   'extreme_upper_monthly_median':500.0,  'aqs_code':'',      'airbase_name':'Zinc in PM10 (aerosol)',               'airbase_code':'05063', 'ebas_parameter_name':'zinc',                         'ebas_matrix':'pm10',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Zn',      'bsc_parameter_name':'pm10zn'          },
'PM2.5_Al':            {'long_parameter_name':'PM2.5 aluminium',                                                                                                                  'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'aluminium in PM2.5 (aerosol)',         'airbase_code':'01604', 'ebas_parameter_name':'aluminium',                    'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Al',      'bsc_parameter_name':'pm2p5al'         },
'PM2.5_As':            {'long_parameter_name':'PM2.5 arsenic',                                                                                                                    'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Arsenic in PM2.5 (aerosol)',           'airbase_code':'01018', 'ebas_parameter_name':'arsenic',                      'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'As',      'bsc_parameter_name':'pm2p5as'         },
'PM2.5_BC':            {'long_parameter_name':'PM2.5 black carbon',                                                                                                               'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'black_carbon in PM2.5 (aerosol)',      'airbase_code':'01391', 'ebas_parameter_name':'',                             'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm2p5bc'         },
'PM2.5_eBC':           {'long_parameter_name':'PM2.5 equivalent black carbon',                                                                                                    'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'equivalent_black_carbon',      'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm2p5ebc'        },
'PM2.5_C':             {'long_parameter_name':'PM2.5 carbon',                                                                                                                     'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'total_carbon',                 'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm2p5c'          },
'PM2.5_Ca++':          {'long_parameter_name':'PM2.5 calcium',                                                                                                                    'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'calcium in PM2.5 (aerosol)',           'airbase_code':'01629', 'ebas_parameter_name':'calcium',                      'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ca',      'bsc_parameter_name':'pm2p5ca'         },
'PM2.5_Cd':            {'long_parameter_name':'PM2.5 cadmium',                                                                                                                    'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Cadmium in PM2.5 (aerosol)',           'airbase_code':'01014', 'ebas_parameter_name':'cadmium',                      'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cd',      'bsc_parameter_name':'pm2p5cd'         },
'PM2.5_Cl-':           {'long_parameter_name':'PM2.5 chloride',                                                                                                                   'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'chloride in PM2.5 (aerosol)',          'airbase_code':'01631', 'ebas_parameter_name':'chloride',                     'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cl',      'bsc_parameter_name':'pm2p5cl'         },
'PM2.5_Co':            {'long_parameter_name':'PM2.5 cobalt',                                                                                                                     'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Cobalt in PM2.5 (aerosol)',            'airbase_code':'01064', 'ebas_parameter_name':'cobalt',                       'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Co',      'bsc_parameter_name':'pm2p5cobalt'     },
'PM2.5_Cr':            {'long_parameter_name':'PM2.5 chromium',                                                                                                                   'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Chromium in PM2.5 (aerosol)',          'airbase_code':'01016', 'ebas_parameter_name':'chromium',                     'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cr',      'bsc_parameter_name':'pm2p5cr'         },
'PM2.5_Cu':            {'long_parameter_name':'PM2.5 copper',                                                                                                                     'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Copper in PM2.5 (aerosol)',            'airbase_code':'01073', 'ebas_parameter_name':'copper',                       'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cu',      'bsc_parameter_name':'pm2p5cu'         },
'PM2.5_DUST':          {'long_parameter_name':'PM2.5 dust',                                                                                                                       'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'pm2p5du'         },
'PM2.5_EC':            {'long_parameter_name':'PM2.5 elemental carbon',                                                                                                           'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Elemental carbon in PM2.5 (aerosol)',  'airbase_code':'01771', 'ebas_parameter_name':'elemental_carbon',             'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm2p5ec'         },
'PM2.5_Fe':            {'long_parameter_name':'PM2.5 iron',                                                                                                                       'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'Iron in PM2.5 (aerosol)',              'airbase_code':'01065', 'ebas_parameter_name':'iron',                         'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Fe',      'bsc_parameter_name':'pm2p5fe'         },
'PM2.5_Hg':            {'long_parameter_name':'PM2.5 mercury',                                                                                                                    'matrix':'pm2.5',   'standard_units':'pg m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Mercury in PM2.5 (aerosol)',           'airbase_code':'01013', 'ebas_parameter_name':'mercury',                      'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Hg',      'bsc_parameter_name':'pm2p5hg'         },
'PM2.5_K+':            {'long_parameter_name':'PM2.5 potassium',                                                                                                                  'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'potassium in PM2.5 (aerosol)',         'airbase_code':'01657', 'ebas_parameter_name':'potassium',                    'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'K',       'bsc_parameter_name':'pm2p5k'          },
'PM2.5_Mg++':          {'long_parameter_name':'PM2.5 magnesium',                                                                                                                  'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'magnesium in PM2.5 (aerosol)',         'airbase_code':'01659', 'ebas_parameter_name':'magnesium',                    'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mg',      'bsc_parameter_name':'pm2p5mg'         },
'PM2.5_Mn':            {'long_parameter_name':'PM2.5 manganese',                                                                                                                  'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Manganese in PM2.5 (aerosol)',         'airbase_code':'01017', 'ebas_parameter_name':'manganese',                    'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mn',      'bsc_parameter_name':'pm2p5mn'         },
'PM2.5_MSA':           {'long_parameter_name':'PM2.5 methanesulfonic acid',                                                                                                       'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'methanesulfonic_acid',         'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'CH4O3S',  'bsc_parameter_name':'pm2p5msa'        },
'PM2.5_Na+':           {'long_parameter_name':'PM2.5 sodium',                                                                                                                     'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'sodium in PM2.5 (aerosol)',            'airbase_code':'01668', 'ebas_parameter_name':'sodium',                       'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Na',      'bsc_parameter_name':'pm2p5na'         },
'PM2.5_NH4+':          {'long_parameter_name':'PM2.5 ammonium',                                                                                                                   'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Ammonium in PM2.5 (aerosol)',          'airbase_code':'01045', 'ebas_parameter_name':'ammonium',                     'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4',     'bsc_parameter_name':'pm2p5nh4'        },
'PM2.5_NH4NO3':        {'long_parameter_name':'PM2.5 ammonium nitrate',                                                                                                           'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'ammonium_nitrate',             'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4NO3',  'bsc_parameter_name':'pm2p5nh4no3'     },
'PM2.5_Ni':            {'long_parameter_name':'PM2.5 nickel',                                                                                                                     'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Nickel in PM2.5 (aerosol)',            'airbase_code':'01015', 'ebas_parameter_name':'nickel',                       'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ni',      'bsc_parameter_name':'pm2p5ni'         },
'PM2.5_NO3-':          {'long_parameter_name':'PM2.5 nitrate',                                                                                                                    'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Nitrate in PM2.5 (aerosol)',           'airbase_code':'01046', 'ebas_parameter_name':'nitrate',                      'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NO3',     'bsc_parameter_name':'pm2p5no3'        },
'PM2.5_OC':            {'long_parameter_name':'PM2.5 organic carbon',                                                                                                             'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Organic carbon in PM2.5 (aerosol)',    'airbase_code':'01772', 'ebas_parameter_name':'organic_carbon',               'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm2p5oc'         },
'PM2.5_Pb':            {'long_parameter_name':'PM2.5 lead',                                                                                                                       'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Lead in PM2.5 (aerosol)',              'airbase_code':'01012', 'ebas_parameter_name':'lead',                         'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Pb',      'bsc_parameter_name':'pm2p5pb'         },
'PM2.5_Se':            {'long_parameter_name':'PM2.5 selenium',                                                                                                                   'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Selenium in PM2.5 (aerosol)',          'airbase_code':'01048', 'ebas_parameter_name':'selenium',                     'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Se',      'bsc_parameter_name':'pm2p5se'         },
'PM2.5_SO4--':         {'long_parameter_name':'PM2.5 sulphate',                                                                                                                   'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'sulphate in PM2.5 (aerosol)',          'airbase_code':'01047', 'ebas_parameter_name':'sulphate_total',               'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'pm2p5so4'        },
'PM2.5_SO4--_NSS':     {'long_parameter_name':'PM2.5 sulphate : non-sea salt',                                                                                                    'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_corrected',           'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'pm2p5so4nss'     },
'PM2.5_SO4--_SS':      {'long_parameter_name':'PM2.5 sulphate : sea salt',                                                                                                        'matrix':'pm2.5',   'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_seasalt',             'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'pm2p5so4ss'      },
'PM2.5_V':             {'long_parameter_name':'PM2.5 vanadium',                                                                                                                   'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Vanedium in PM2.5 (aerosol)',          'airbase_code':'01049', 'ebas_parameter_name':'vanadium',                     'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'V',       'bsc_parameter_name':'pm2p5v'          },
'PM2.5_Zn':            {'long_parameter_name':'PM2.5 zinc',                                                                                                                       'matrix':'pm2.5',   'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':5000.0,   'extreme_upper_monthly_median':500.0,  'aqs_code':'',      'airbase_name':'Zinc in PM2.5 (aerosol)',              'airbase_code':'01063', 'ebas_parameter_name':'zinc',                         'ebas_matrix':'pm25',    'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Zn',      'bsc_parameter_name':'pm2p5zn'         },
'PM1_Al':              {'long_parameter_name':'PM1 aluminium',                                                                                                                    'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'aluminium',                    'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Al',      'bsc_parameter_name':'pm1al'           },
'PM1_As':              {'long_parameter_name':'PM1 arsenic',                                                                                                                      'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'arsenic',                      'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'As',      'bsc_parameter_name':'pm1as'           },
'PM1_BC':              {'long_parameter_name':'PM1 black carbon',                                                                                                                 'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm1bc'           },
'PM1_eBC':             {'long_parameter_name':'PM1 equivalent black carbon',                                                                                                      'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'equivalent_black_carbon',      'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm1ebc'          },
'PM1_C':               {'long_parameter_name':'PM1 carbon',                                                                                                                       'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'total_carbon',                 'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm1c'            },
'PM1_Ca++':            {'long_parameter_name':'PM1 calcium',                                                                                                                      'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'calcium',                      'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ca',      'bsc_parameter_name':'pm1ca'           },
'PM1_Cd':              {'long_parameter_name':'PM1 cadmium',                                                                                                                      'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'cadmium',                      'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cd',      'bsc_parameter_name':'pm1cd'           },
'PM1_Cl-':             {'long_parameter_name':'PM1 chloride',                                                                                                                     'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'chloride',                     'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cl',      'bsc_parameter_name':'pm1cl'           },
'PM1_Co':              {'long_parameter_name':'PM1 cobalt',                                                                                                                       'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'cobalt',                       'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Co',      'bsc_parameter_name':'pm1cobalt'       },
'PM1_Cr':              {'long_parameter_name':'PM1 chromium',                                                                                                                     'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'chromium',                     'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cr',      'bsc_parameter_name':'pm1cr'           },
'PM1_Cu':              {'long_parameter_name':'PM1 copper',                                                                                                                       'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'copper',                       'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cu',      'bsc_parameter_name':'pm1cu'           },
'PM1_DUST':            {'long_parameter_name':'PM1 dust',                                                                                                                         'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'pm1du'           },
'PM1_EC':              {'long_parameter_name':'PM1 elemental carbon',                                                                                                             'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'elemental_carbon',             'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm1ec'           },
'PM1_Fe':              {'long_parameter_name':'PM1 iron',                                                                                                                         'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'iron',                         'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Fe',      'bsc_parameter_name':'pm1fe'           },
'PM1_Hg':              {'long_parameter_name':'PM1 mercury',                                                                                                                      'matrix':'pm1',     'standard_units':'pg m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'mercury',                      'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Hg',      'bsc_parameter_name':'pm1hg'           },
'PM1_K+':              {'long_parameter_name':'PM1 potassium',                                                                                                                    'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'potassium',                    'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'K',       'bsc_parameter_name':'pm1k'            },
'PM1_Mg++':            {'long_parameter_name':'PM1 magnesium',                                                                                                                    'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'magnesium',                    'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mg',      'bsc_parameter_name':'pm1mg'           },
'PM1_Mn':              {'long_parameter_name':'PM1 manganese',                                                                                                                    'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'manganese',                    'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mn',      'bsc_parameter_name':'pm1mn'           },
'PM1_MSA':             {'long_parameter_name':'PM1 methanesulfonic acid',                                                                                                         'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'methanesulfonic_acid',         'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'CH4O3S',  'bsc_parameter_name':'pm1msa'          },
'PM1_Na+':             {'long_parameter_name':'PM1 sodium',                                                                                                                       'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sodium',                       'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Na',      'bsc_parameter_name':'pm1na'           },
'PM1_NH4+':            {'long_parameter_name':'PM1 ammonium',                                                                                                                     'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'ammonium',                     'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4',     'bsc_parameter_name':'pm1nh4'          },
'PM1_NH4NO3':          {'long_parameter_name':'PM1 ammonium nitrate',                                                                                                             'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'ammonium_nitrate',             'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4NO3',  'bsc_parameter_name':'pm1nh4no3'       },
'PM1_Ni':              {'long_parameter_name':'PM1 nickel',                                                                                                                       'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'nickel',                       'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ni',      'bsc_parameter_name':'pm1ni'           },
'PM1_NO3-':            {'long_parameter_name':'PM1 nitrate',                                                                                                                      'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'nitrate',                      'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NO3',     'bsc_parameter_name':'pm1no3'          },
'PM1_OC':              {'long_parameter_name':'PM1 organic carbon',                                                                                                               'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'organic_carbon',               'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'C',       'bsc_parameter_name':'pm1oc'           },
'PM1_Pb':              {'long_parameter_name':'PM1 lead',                                                                                                                         'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'lead',                         'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Pb',      'bsc_parameter_name':'pm1pb'           },
'PM1_Se':              {'long_parameter_name':'PM1 selenium',                                                                                                                     'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'selenium',                     'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Se',      'bsc_parameter_name':'pm1se'           },
'PM1_SO4--':           {'long_parameter_name':'PM1 sulphate',                                                                                                                     'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_total',               'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'pm1so4'          },
'PM1_SO4--_NSS':       {'long_parameter_name':'PM1 sulphate : non-sea salt',                                                                                                      'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_corrected',           'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'pm1so4nss'       },
'PM1_SO4--_SS':        {'long_parameter_name':'PM1 sulphate : sea salt',                                                                                                          'matrix':'pm1',     'standard_units':'ug m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_seasalt',             'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'pm1so4ss'        },
'PM1_V':               {'long_parameter_name':'PM1 vanadium',                                                                                                                     'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'vanadium',                     'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'V',       'bsc_parameter_name':'pm1v'            },
'PM1_Zn':              {'long_parameter_name':'PM10 zinc',                                                                                                                        'matrix':'pm1',     'standard_units':'ng m-3',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':5000.0,   'extreme_upper_monthly_median':500.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'zinc',                         'ebas_matrix':'pm1',     'ebas_priority_units':['/m3'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Zn',      'bsc_parameter_name':'pm1zn'           },
'Al_PREC':             {'long_parameter_name':'aluminium in precipitaton',                                                                                                        'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'aluminium (precip)',                   'airbase_code':'00605', 'ebas_parameter_name':'aluminium',                    'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Al',      'bsc_parameter_name':'precal'          },
'As_PREC':             {'long_parameter_name':'arsenic in precipitaton',                                                                                                          'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Arsenic (precip)',                     'airbase_code':'02018', 'ebas_parameter_name':'arsenic',                      'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'As',      'bsc_parameter_name':'precas'          },
'Ca++_PREC':           {'long_parameter_name':'calcium in precipitaton',                                                                                                          'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'calcium (precip)',                     'airbase_code':'00630', 'ebas_parameter_name':'calcium',                      'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ca',      'bsc_parameter_name':'precca'          },
'Cd_PREC':             {'long_parameter_name':'cadmium in precipitaton',                                                                                                          'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Cadmium (precip)',                     'airbase_code':'02014', 'ebas_parameter_name':'cadmium',                      'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cd',      'bsc_parameter_name':'preccd'          },
'Cl-_PREC':            {'long_parameter_name':'chloride in precipitaton',                                                                                                         'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'chloride (precip)',                    'airbase_code':'00632', 'ebas_parameter_name':'chloride',                     'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cl',      'bsc_parameter_name':'preccl'          },
'Co_PREC':             {'long_parameter_name':'cobalt in precipitaton',                                                                                                           'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Cobalt (precip)',                      'airbase_code':'02064', 'ebas_parameter_name':'cobalt',                       'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Co',      'bsc_parameter_name':'preccobalt'      },
'Cr_PREC':             {'long_parameter_name':'chromium in precipitaton',                                                                                                         'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'Chromium (precip)',                    'airbase_code':'02016', 'ebas_parameter_name':'chromium',                     'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cr',      'bsc_parameter_name':'preccr'          },
'Cu_PREC':             {'long_parameter_name':'copper in precipitaton',                                                                                                           'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Copper (precip)',                      'airbase_code':'02073', 'ebas_parameter_name':'copper',                       'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cu',      'bsc_parameter_name':'preccu'          },
'Fe_PREC':             {'long_parameter_name':'iron in precipitaton',                                                                                                             'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':10.0,   'extreme_lower_limit':0.0,     'extreme_upper_limit':20000.0,  'extreme_upper_monthly_median':2000.0, 'aqs_code':'',      'airbase_name':'Iron (precip)',                        'airbase_code':'02065', 'ebas_parameter_name':'iron',                         'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Fe',      'bsc_parameter_name':'precfe'          },
'Hg_PREC':             {'long_parameter_name':'mercury in precipitaton',                                                                                                          'matrix':'precip',  'standard_units':'ng l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Mercury (precip)',                     'airbase_code':'02013', 'ebas_parameter_name':'mercury',                      'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Hg',      'bsc_parameter_name':'prechg'          },
'K+_PREC':             {'long_parameter_name':'potassium in precipitaton',                                                                                                        'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'potassium (precip)',                   'airbase_code':'00658', 'ebas_parameter_name':'potassium',                    'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'K',       'bsc_parameter_name':'preck'           },
'Mg++_PREC':           {'long_parameter_name':'magnesium in precipitaton',                                                                                                        'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'magnesium (precip)',                   'airbase_code':'00660', 'ebas_parameter_name':'magnesium',                    'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mg',      'bsc_parameter_name':'precmg'          },
'Mn_PREC':             {'long_parameter_name':'manganese in precipitaton',                                                                                                        'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Manganese (precip)',                   'airbase_code':'02017', 'ebas_parameter_name':'manganese',                    'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mn',      'bsc_parameter_name':'precmn'          },
'MSA_PREC':            {'long_parameter_name':'methanesulfonic acid in precipitaton',                                                                                             'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'methanesulfonic_acid',         'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'CH4O3S',  'bsc_parameter_name':'precmsa'         },
'Na+_PREC':            {'long_parameter_name':'sodium in precipitaton',                                                                                                           'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'sodium (precip)',                      'airbase_code':'00669', 'ebas_parameter_name':'sodium',                       'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Na',      'bsc_parameter_name':'precna'          },
'NH4+_PREC':           {'long_parameter_name':'ammonium in precipitaton',                                                                                                         'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'ammonium (precip)',                    'airbase_code':'00664', 'ebas_parameter_name':'ammonium',                     'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4',     'bsc_parameter_name':'precnh4'         },
'NH4NO3_PREC':         {'long_parameter_name':'ammonium nitrate in precipitaton',                                                                                                 'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'ammonium_nitrate',             'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4NO3',  'bsc_parameter_name':'precnh4no3'      },
'Ni_PREC':             {'long_parameter_name':'nickel in precipitaton',                                                                                                           'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Nickel (precip)',                      'airbase_code':'02015', 'ebas_parameter_name':'nickel',                       'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ni',      'bsc_parameter_name':'precni'          },
'NO3-_PREC':           {'long_parameter_name':'nitrate in precipitaton',                                                                                                          'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'nitrate (precip)',                     'airbase_code':'00666', 'ebas_parameter_name':'nitrate',                      'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NO3',     'bsc_parameter_name':'precno3'         },
'Pb_PREC':             {'long_parameter_name':'lead in precipitaton',                                                                                                             'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Lead (precip)',                        'airbase_code':'02012', 'ebas_parameter_name':'lead',                         'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Pb',      'bsc_parameter_name':'precpb'          },
'Se_PREC':             {'long_parameter_name':'selenium in precipitaton',                                                                                                         'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':0.2,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':100.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'selenium',                     'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Se',      'bsc_parameter_name':'precse'          },
'SO4--_PREC':          {'long_parameter_name':'sulphate in precipitaton',                                                                                                         'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'sulphate (precip)',                    'airbase_code':'00719', 'ebas_parameter_name':'sulphate_total',               'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'precso4'         },
'SO4--_NSS_PREC':      {'long_parameter_name':'sulphate: non-sea salt in precipitaton',                                                                                           'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'sulphate_corrected (precip)',          'airbase_code':'00720', 'ebas_parameter_name':'sulphate_corrected',           'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'precso4nss'      },
'SO4--_SS_PREC':       {'long_parameter_name':'sulphate: sea salt in precipitaton',                                                                                               'matrix':'precip',  'standard_units':'mg l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_seasalt',             'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'precso4ss'       },
'V_PREC':              {'long_parameter_name':'vanadium in precipitaton',                                                                                                         'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'vanadium (precip)',                    'airbase_code':'00728', 'ebas_parameter_name':'vanadium',                     'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'V',       'bsc_parameter_name':'precv'           },
'Zn_PREC':             {'long_parameter_name':'zinc in precipitaton',                                                                                                             'matrix':'precip',  'standard_units':'ug l-1',          'axis_label':'concentration', 'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Zinc (precip)',                        'airbase_code':'02063', 'ebas_parameter_name':'zinc',                         'ebas_matrix':'precip',  'ebas_priority_units':['/m3','/l'],      'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Zn',      'bsc_parameter_name':'preczn'          },
'Al_WETDEP':           {'long_parameter_name':'aluminium wet deposition flux',                                                                                                    'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'aluminium (precip)',                   'airbase_code':'00605', 'ebas_parameter_name':'aluminium',                    'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Al',      'bsc_parameter_name':'wetal'           },
'As_WETDEP':           {'long_parameter_name':'arsenic wet deposition flux',                                                                                                      'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Arsenic (precip)',                     'airbase_code':'02018', 'ebas_parameter_name':'arsenic',                      'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'As',      'bsc_parameter_name':'wetas'           },
'Ca++_WETDEP':         {'long_parameter_name':'calcium wet deposition flux',                                                                                                      'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'calcium (precip)',                     'airbase_code':'00630', 'ebas_parameter_name':'calcium',                      'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ca',      'bsc_parameter_name':'wetca'           },
'Cd_WETDEP':           {'long_parameter_name':'cadmium wet deposition flux',                                                                                                      'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Cadmium (precip)',                     'airbase_code':'02014', 'ebas_parameter_name':'cadmium',                      'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cd',      'bsc_parameter_name':'wetcd'           },
'Cl-_WETDEP':          {'long_parameter_name':'chloride wet deposition flux',                                                                                                     'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'chloride (precip)',                    'airbase_code':'00632', 'ebas_parameter_name':'chloride',                     'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cl',      'bsc_parameter_name':'wetcl'           },
'Co_WETDEP':           {'long_parameter_name':'cobalt wet deposition flux',                                                                                                       'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Cobalt (precip)',                      'airbase_code':'02064', 'ebas_parameter_name':'cobalt',                       'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Co',      'bsc_parameter_name':'wetcobalt'       },
'Cr_WETDEP':           {'long_parameter_name':'chromium wet deposition flux',                                                                                                     'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Chromium (precip)',                    'airbase_code':'02016', 'ebas_parameter_name':'chromium',                     'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cr',      'bsc_parameter_name':'wetcr'           },
'Cu_WETDEP':           {'long_parameter_name':'copper wet deposition flux',                                                                                                       'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Copper (precip)',                      'airbase_code':'02073', 'ebas_parameter_name':'copper',                       'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Cu',      'bsc_parameter_name':'wetcu'           },
'Fe_WETDEP':           {'long_parameter_name':'iron wet deposition flux',                                                                                                         'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Iron (precip)',                        'airbase_code':'02065', 'ebas_parameter_name':'iron',                         'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Fe',      'bsc_parameter_name':'wetfe'           },
'Hg_WETDEP':           {'long_parameter_name':'mercury wet deposition flux',                                                                                                      'matrix':'wetdep',  'standard_units':'pg m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Mercury (precip)',                     'airbase_code':'02013', 'ebas_parameter_name':'mercury',                      'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Hg',      'bsc_parameter_name':'wethg'           },
'K+_WETDEP':           {'long_parameter_name':'potassium wet deposition flux',                                                                                                    'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'potassium (precip)',                   'airbase_code':'00658', 'ebas_parameter_name':'potassium',                    'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'K',       'bsc_parameter_name':'wetk'            },
'Mg++_WETDEP':         {'long_parameter_name':'magnesium wet deposition flux',                                                                                                    'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'magnesium (precip)',                   'airbase_code':'00660', 'ebas_parameter_name':'magnesium',                    'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mg',      'bsc_parameter_name':'wetmg'           },
'Mn_WETDEP':           {'long_parameter_name':'manganese wet deposition flux',                                                                                                    'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Manganese (precip)',                   'airbase_code':'02017', 'ebas_parameter_name':'manganese',                    'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Mn',      'bsc_parameter_name':'wetmn'           },
'MSA_WETDEP':          {'long_parameter_name':'methanesulfonic acid wet deposition flux',                                                                                         'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'methanesulfonic_acid',         'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'CH4O3S',  'bsc_parameter_name':'wetmsa'          },
'Na+_WETDEP':          {'long_parameter_name':'sodium wet deposition flux',                                                                                                       'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'sodium (precip)',                      'airbase_code':'00669', 'ebas_parameter_name':'sodium',                       'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Na',      'bsc_parameter_name':'wetna'           },
'NH4+_WETDEP':         {'long_parameter_name':'ammonium wet deposition flux',                                                                                                     'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'ammonium (precip)',                    'airbase_code':'00664', 'ebas_parameter_name':'ammonium',                     'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4',     'bsc_parameter_name':'wetnh4'          },
'NH4NO3_WETDEP':       {'long_parameter_name':'ammonium nitrate wet deposition flux',                                                                                             'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'ammonium_nitrate',             'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NH4NO3',  'bsc_parameter_name':'wetnh4no3'       },
'Ni_WETDEP':           {'long_parameter_name':'nickel wet deposition flux',                                                                                                       'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Nickel (precip)',                      'airbase_code':'02015', 'ebas_parameter_name':'nickel',                       'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Ni',      'bsc_parameter_name':'wetni'           },
'NO3-_WETDEP':         {'long_parameter_name':'nitrate wet deposition flux',                                                                                                      'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'nitrate (precip)',                     'airbase_code':'00666', 'ebas_parameter_name':'nitrate',                      'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'NO3',     'bsc_parameter_name':'wetno3'          },
'Pb_WETDEP':           {'long_parameter_name':'lead wet deposition flux',                                                                                                         'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Lead (precip)',                        'airbase_code':'02012', 'ebas_parameter_name':'lead',                         'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Pb',      'bsc_parameter_name':'wetpb'           },
'Se_WETDEP':           {'long_parameter_name':'selenium wet deposition flux',                                                                                                     'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'selenium',                     'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Se',      'bsc_parameter_name':'wetse'           },
'SO4--_WETDEP':        {'long_parameter_name':'sulphate wet deposition flux',                                                                                                     'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'sulphate (precip)',                    'airbase_code':'00719', 'ebas_parameter_name':'sulphate_total',               'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'wetso4'          },
'SO4--_NSS_WETDEP':    {'long_parameter_name':'sulphate: non-sea salt wet deposition flux',                                                                                       'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'sulphate_corrected (precip)',          'airbase_code':'00720', 'ebas_parameter_name':'sulphate_corrected',           'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'wetso4nss'       },
'SO4--_SS_WETDEP':     {'long_parameter_name':'sulphate: sea salt wet deposition flux',                                                                                           'matrix':'wetdep',  'standard_units':'ug m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'sulphate_seasalt',             'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'SO4',     'bsc_parameter_name':'wetso4ss'        },
'V_WETDEP':            {'long_parameter_name':'vanadium wet deposition flux',                                                                                                     'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'vanadium (precip)',                    'airbase_code':'00728', 'ebas_parameter_name':'vanadium',                     'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'V',       'bsc_parameter_name':'wetv'            },
'Zn_WETDEP':           {'long_parameter_name':'zinc wet deposition flux',                                                                                                         'matrix':'wetdep',  'standard_units':'ng m-2 min-1',    'axis_label':'flux',          'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':2000.0,   'extreme_upper_monthly_median':200.0,  'aqs_code':'',      'airbase_name':'Zinc (precip)',                        'airbase_code':'02063', 'ebas_parameter_name':'zinc',                         'ebas_matrix':'wetdep',  'ebas_priority_units':['/m2'],           'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'Zn',      'bsc_parameter_name':'wetzn'           },
'WND_DIR_10M':         {'long_parameter_name':'10m wind direction',                                                                                                               'matrix':'meteo',   'standard_units':'angular degrees', 'axis_label':'direction',     'minimum_resolution':15.0,   'extreme_lower_limit':1.0,     'extreme_upper_limit':360.0,    'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'direction_angle',    'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'dir10'           }, 
'WND_SPD_10M':         {'long_parameter_name':'10m wind speed',                                                                                                                   'matrix':'meteo',   'standard_units':'m s-1',           'axis_label':'speed',         'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':900.0,    'extreme_upper_monthly_median':50.0,   'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'speed',              'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'spd10'           },  
'TEMP_2M':             {'long_parameter_name':'2m air temperature',                                                                                                               'matrix':'meteo',   'standard_units':'K',               'axis_label':'temperature',   'minimum_resolution':1.0,    'extreme_lower_limit':-658.85, 'extreme_upper_limit':891.15,   'extreme_upper_monthly_median':323.15, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'air_temp',           'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'t2'              },  
'RH_2M':               {'long_parameter_name':'average relative humidity',                                                                                                        'matrix':'meteo',   'standard_units':'unitless',        'axis_label':'humidity (%)',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':1000.0,   'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'avg_rh',             'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'rho2'            },  
'SST':                 {'long_parameter_name':'sea surface temperature',                                                                                                          'matrix':'meteo',   'standard_units':'K',               'axis_label':'temperature',   'minimum_resolution':1.0,    'extreme_lower_limit':223.15,  'extreme_upper_limit':723.15,   'extreme_upper_monthly_median':50.0,   'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'sst',                'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'sst'             },  
'DEW_PT_2M':           {'long_parameter_name':'2m dew point',                                                                                                                     'matrix':'meteo',   'standard_units':'K',               'axis_label':'temperature',   'minimum_resolution':1.0,    'extreme_lower_limit':-708.85, 'extreme_upper_limit':641.15,   'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'dpoint_temp',        'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'td2'             },  
'PRES_2M':             {'long_parameter_name':'2m atmospheric pressure',                                                                                                          'matrix':'meteo',   'standard_units':'hPa',             'axis_label':'pressure',      'minimum_resolution':10.0,   'extreme_lower_limit':450.0,   'extreme_upper_limit':10900.0,  'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'station_p',          'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'pshltr'          },  
'SLP':                 {'long_parameter_name':'sea level pressure',                                                                                                               'matrix':'meteo',   'standard_units':'hPa',             'axis_label':'pressure',      'minimum_resolution':10.0,   'extreme_lower_limit':8600.0,  'extreme_upper_limit':10900.0,  'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'sea_level_pressure', 'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'slp'             },  
'PREC_ACCUM':          {'long_parameter_name':'liquid precipitation accumulation',                                                                                                'matrix':'meteo',   'standard_units':'mm',              'axis_label':'accumulation',  'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':9998.0,   'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'prec_depth',         'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'acprec'          },  
'SNOW_ACCUM':          {'long_parameter_name':'snow accumulation',                                                                                                                'matrix':'meteo',   'standard_units':'cm',              'axis_label':'accumulation',  'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1200.0,   'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'snow_acc',           'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'acsnow'          },  
'SNOW_DEPTH':          {'long_parameter_name':'snow depth',                                                                                                                       'matrix':'meteo',   'standard_units':'cm',              'axis_label':'depth',         'minimum_resolution':1.0,    'extreme_lower_limit':0.0,     'extreme_upper_limit':1200.0,   'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'si',                 'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'si'              },  
'CEILING_HEIGHT':      {'long_parameter_name':'the height above ground level of the lowest cloud or obscuring phenomena layer aloft with 5/8 or more summation total sky cover.', 'matrix':'meteo',   'standard_units':'m',               'axis_label':'height',        'minimum_resolution':200.0,  'extreme_lower_limit':0.0,     'extreme_upper_limit':22000.0,  'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'ceiling',            'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'cldbot'          },  
'VIS_DIST':            {'long_parameter_name':'visibility distance',                                                                                                              'matrix':'meteo',   'standard_units':'m',               'axis_label':'distance',      'minimum_resolution':100.0,  'extreme_lower_limit':0.0,     'extreme_upper_limit':160000.0, 'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'distance',           'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'vdist'           },  
'CLOUD_CVG':           {'long_parameter_name':'cloud coverage, taken as maximum cloud coverage across cloud layers',                                                              'matrix':'meteo',   'standard_units':'oktas',           'axis_label':'coverage',      'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':8.0,      'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'cloud',              'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'ccovmean'        },  
'CLOUD_CVG_FRAC':      {'long_parameter_name':'cloud coverage fraction, taken as maximum cloud coverage across cloud layers',                                                     'matrix':'meteo',   'standard_units':'unitless',        'axis_label':'coverage (%)',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':100.0,    'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'cloud_frac',         'aeronet_parameter_name':'',                             'aeronet_matrix':'',          'chemical_formula':'',        'bsc_parameter_name':'cfracmean'       },  
'AOD_500nm':           {'long_parameter_name':'aerosol optical depth at 500nm',                                                                                                   'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 500nm',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':20.0,     'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'Total_AOD_500nm[tau_a]',       'aeronet_matrix':'oneill',    'chemical_formula':'',        'bsc_parameter_name':'od500aero'       },  
'AOD_500nm_COARSE':    {'long_parameter_name':'coarse mode aerosol optical depth at 500nm',                                                                                       'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 500nm',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':20.0,     'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'Coarse_Mode_AOD_500nm[tau_c]', 'aeronet_matrix':'oneill',    'chemical_formula':'',        'bsc_parameter_name':'od500aerocoarse' },  
'AOD_500nm_FINE':      {'long_parameter_name':'fine mode aerosol optical depth at 500nm',                                                                                         'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 500nm',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':20.0,     'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'Fine_Mode_AOD_500nm[tau_f]',   'aeronet_matrix':'oneill',    'chemical_formula':'',        'bsc_parameter_name':'od500aerofine'   },  
'FINE_MODE_FRAC_500nm':{'long_parameter_name':'fine mode fraction at 500nm',                                                                                                      'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'',              'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':1.0,      'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'FineModeFraction_500nm[eta]',  'aeronet_matrix':'oneill',    'chemical_formula':'',        'bsc_parameter_name':'fm500frac'       },  
'AOD_380nm':           {'long_parameter_name':'aerosol optical depth at 380nm',                                                                                                   'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 380nm',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':20.0,     'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'AOD_380nm',                    'aeronet_matrix':'directsun', 'chemical_formula':'',        'bsc_parameter_name':'od380aero'       },  
'AOD_440nm':           {'long_parameter_name':'aerosol optical depth at 440nm',                                                                                                   'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 440nm',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':20.0,     'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'AOD_440nm',                    'aeronet_matrix':'directsun', 'chemical_formula':'',        'bsc_parameter_name':'od440aero'       },  
'AOD_550nm':           {'long_parameter_name':'aerosol optical depth at 550nm',                                                                                                   'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 550nm',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':20.0,     'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'AOD_550nm',                    'aeronet_matrix':'directsun', 'chemical_formula':'',        'bsc_parameter_name':'od550aero'       },  
'AOD_675nm':           {'long_parameter_name':'aerosol optical depth at 675nm',                                                                                                   'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 675nm',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':20.0,     'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'AOD_675nm',                    'aeronet_matrix':'directsun', 'chemical_formula':'',        'bsc_parameter_name':'od675aero'       },  
'AOD_870nm':           {'long_parameter_name':'aerosol optical depth at 870nm',                                                                                                   'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 870nm',  'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':20.0,     'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'AOD_870nm',                    'aeronet_matrix':'directsun', 'chemical_formula':'',        'bsc_parameter_name':'od870aero'       },  
'AOD_1020nm':          {'long_parameter_name':'aerosol optical depth at 1020nm',                                                                                                  'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 1020nm', 'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':20.0,     'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'AOD_1020nm',                   'aeronet_matrix':'directsun', 'chemical_formula':'',        'bsc_parameter_name':'od1020aero'      },  
'AE_440-870nm':        {'long_parameter_name':'angstrom exponent between 440 and 870 nm',                                                                                         'matrix':'aerosol', 'standard_units':'unitless',        'axis_label':'AOD at 1020nm', 'minimum_resolution':np.NaN, 'extreme_lower_limit':0.0,     'extreme_upper_limit':5.0,      'extreme_upper_monthly_median':np.NaN, 'aqs_code':'',      'airbase_name':'',                                     'airbase_code':'',      'ebas_parameter_name':'',                             'ebas_matrix':'',        'ebas_priority_units':[''],              'naps_code':'',    'ncdc_isd_parameter_name':'',                   'aeronet_parameter_name':'440-870_Angstrom_Exponent',    'aeronet_matrix':'directsun', 'chemical_formula':'',        'bsc_parameter_name':'ae440-870aero'   }  
}

###--------------------------------------------------------------------------------------------------###
###NETWORK STANDARDISATIONS
###--------------------------------------------------------------------------------------------------###

#Define dictionary with details pertaining to standardised networks.
#information regarding the hierarchy of preferred data levels is detailed, by data network, and returns the ranking of a specific data level for data of a specific network.
#Data with a data level of a high defined ranking, is preferentially selected over data with a low data level ranking, when there is data duplication.       

standard_networks = {
'AERONET_v3':       {'file_storage_type':'all_per_day',        'external_metadata_files':['AERONET_META.csv'],                                          'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'2.0':3, '1.5':2, '1.0':1}},    
'CANADA_NAPS':      {'file_storage_type':'all_per_year',       'external_metadata_files':[],                                                            'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'1.0':1}},
'CAPMoN':           {'file_storage_type':'all_per_year',       'external_metadata_files':[],                                                            'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'2.0':1}},
'EANET':            {'file_storage_type':'site_per_year',      'external_metadata_files':[],                                                            'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'1.0':1}},
'EBAS':             {'file_storage_type':'site_per_multiyear', 'external_metadata_files':['EBAS_META_MANUAL.csv'],                                      'metadata_data_join_variable_name':'reference', 'data_level_hierarchy':{'5.0':8, '4.5':7, '4.0':6, '3.0':5, '2.0':4, '1.5':3, '1.0':2, '0.0':1}},
'EEA_AIRBASE':      {'file_storage_type':'site_per_year',      'external_metadata_files':[],                                                            'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'1.0':1}},
'EEA_AQ_eReporting':{'file_storage_type':'site_per_year',      'external_metadata_files':['EEA_AQ_eReporting_META.csv','EEA_AQ_eReporting_AIDE_D.csv'], 'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'E1a':2,'E2a':1}},
'NCDC_ISD':         {'file_storage_type':'site_per_year',      'external_metadata_files':['isd-history.csv'],                                           'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'1.0':1}}, 
'SEARCH':           {'file_storage_type':'all_per_year',       'external_metadata_files':[],                                                            'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'1.0':1}},
'US_EPA_AQS':       {'file_storage_type':'all_per_year',       'external_metadata_files':[],                                                            'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'1.0':1}},
'US_EPA_CASTNET':   {'file_storage_type':'all_per_year',       'external_metadata_files':[],                                                            'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'1.0':1}},
'WMO_GAW':          {'file_storage_type':'site_per_year',      'external_metadata_files':[],                                                            'metadata_data_join_variable_name':'',          'data_level_hierarchy':{'1.0':1}},
}

###--------------------------------------------------------------------------------------------------###
###TEMPORAL RESOLUTION STANDARDISATIONS
###--------------------------------------------------------------------------------------------------###

#Define dictionary multiple variables associated with standard temporal resolution output codes 

standard_temporal_resolutions = {
'H':  {'temporal_resolution_datetime_frequency':'H',  'temporal_resolution_path':'hourly',                'descriptive_temporal_resolution':'hours'},
'3H': {'temporal_resolution_datetime_frequency':'3H', 'temporal_resolution_path':'3hourly',               'descriptive_temporal_resolution':'hours'},
'6H': {'temporal_resolution_datetime_frequency':'6H', 'temporal_resolution_path':'6hourly',               'descriptive_temporal_resolution':'hours'},
'D':  {'temporal_resolution_datetime_frequency':'D',  'temporal_resolution_path':'daily',                 'descriptive_temporal_resolution':'days'},
'M':  {'temporal_resolution_datetime_frequency':'MS', 'temporal_resolution_path':'monthly',               'descriptive_temporal_resolution':'months'},
'HI': {'temporal_resolution_datetime_frequency':'H',  'temporal_resolution_path':'hourly_instantaneous',  'descriptive_temporal_resolution':'hours'},
'3HI':{'temporal_resolution_datetime_frequency':'3H', 'temporal_resolution_path':'3hourly_instantaneous', 'descriptive_temporal_resolution':'hours'},
'6HI':{'temporal_resolution_datetime_frequency':'6H', 'temporal_resolution_path':'6hourly_instantaneous', 'descriptive_temporal_resolution':'hours'}
}

###--------------------------------------------------------------------------------------------------###
###DATA VARIABLE STANDARDISATIONS
###--------------------------------------------------------------------------------------------------###

#Define dictionary with details pertaining to standardised data variables.
#These data variables include time, the measured parameter variable, as well as variables which inform of the quality of each separate measurement.

def get_standard_data(parameter_details):

    standard_data = {
    #TIME
    'time':                                             {'value':[], 'standard_name':'time',                                              'long_name':'time',                                              'units':'minutes since 0001-01-01 00:00:00', 'data_type':np.uint32,  'string_format':np.NaN, 'description':'Integer time in minutes since 0001-01-01 00:00 UTC. Time given refers to the start of the time window the measurement is representative of (temporal resolution).'},        

    #MEASUREMENT VALUES
    parameter_details['bsc_parameter_name']:            {'value':[], 'standard_name':parameter_details['long_parameter_name'],            'long_name':parameter_details['long_parameter_name'],            'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN, 'description':'Measured value of surface {} for the stated temporal resolution.'.format(parameter_details['long_parameter_name'])}, 

    #PARAMETERS INFORMING QUALITY OF MEASUREMENT 
    'temporal_resolution':                              {'value':[], 'standard_name':'temporal resolution',                               'long_name':'temporal measurement resolution',                   'units':'minutes',                           'data_type':np.uint16,  'string_format':np.NaN, 'description':'Temporal resolution that a specific measurement is representative of. Reported as an integer string (in minutes). If measured resolution is instantaneous, the given integer is 0.'}, 
    'reported_lower_limit_of_detection_per_measurement':{'value':[], 'standard_name':'reported lower limit of detection per measurement', 'long_name':'reported lower limit of detection per measurement', 'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN, 'description':'Reported lower limit of detection of measurement methodology, for a specific measurement, in {}.'.format(parameter_details['standard_units'])}, 
    'reported_upper_limit_of_detection_per_measurement':{'value':[], 'standard_name':'reported upper limit of detection per measurement', 'long_name':'reported upper limit of detection per measurement', 'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN, 'description':'Reported upper limit of detection of measurement methodology, for a specific measurement, in {}.'.format(parameter_details['standard_units'])}, 
    'reported_uncertainty_per_measurement':             {'value':[], 'standard_name':'reported measurement uncertainty per measurement',  'long_name':'reported measurement uncertainty per measurement',  'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN, 'description':'Reported measurement uncertainty (±) of methodology, for a specific measurement. In principal this refers to the inherent uncertainty on every measurement as a function of the quadratic addition of the accuracy and precision metrics (at the same confidence interval), but is often reported incosistently e.g. being solely determined from random errors (i.e. just the measuremental precision). This is given in absolute terms in {}.'.format(parameter_details['standard_units'])},
    'derived_uncertainty_per_measurement':              {'value':[], 'standard_name':'derived measurement uncertainty per measurement',   'long_name':'derived measurement uncertainty per measurement',   'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN, 'description':'Derived measurement uncertainty (±) of methodology, for a specific measurement. This is calculated through the quadratic addition of reported (or if not available, documented) accuracy and precision metrics. This is given in absolute terms in {}.'.format(parameter_details['standard_units'])},
    'day_night_code':                                   {'value':[], 'standard_name':'day/night code',                                    'long_name':'day/night code per measurement',                    'units':'unitless',                          'data_type':np.uint8,   'string_format':np.NaN, 'description':'Binary indication if measurement was made during the day or night. Day=0, Night=1. The classification is made by calculating the solar elevation angle for a latitude/longitude/measurement height at a mid-measurement window timestamp. If the solar elevation angle is > 0, it is classed as daytime, otherwise it is nightime. Classification is 255 if cannot be made.'}, 
    'weekday_weekend_code':                             {'value':[], 'standard_name':'weekday/weekend code',                              'long_name':'weekday/weekend code per measurement',              'units':'unitless',                          'data_type':np.uint8,   'string_format':np.NaN, 'description':'Binary indication if measurement was made during the weekday or weekend. Weekday=0, Weekend=1. The classification is made by evaluating if the local-time of a mid-measurement window timestamp falls on a weekday or on the weekend. Classification is 255 if cannot be made.'}, 
    'season_code':                                      {'value':[], 'standard_name':'season code',                                       'long_name':'season code per measurement',                       'units':'unitless',                          'data_type':np.uint8,   'string_format':np.NaN, 'description':'Code decreeing if measurement was made during the Spring, Summer, Autumn or Winter Seasons. Spring=0, Summer=1, Autumn=2, Winter=3. The classification is made by evaluating which season the local-time of a mid-measurement window timestamp falls in. Classification is 255 if cannot be made.'}, 
    'flag':                                             {'value':[], 'standard_name':'flags',                                             'long_name':'data reporter provided standardised flags',         'units':'unitless',                          'data_type':np.object,  'string_format':np.NaN, 'description':'List of associated data flags per measurement, indicating the data quality of a specific measurement, provided by the native data reporter.'},
    'qa':                                               {'value':[], 'standard_name':'qa',                                                'long_name':'qa standardised data flags',                        'units':'unitless',                          'data_type':np.object,  'string_format':np.NaN, 'description':'List of derived quality assurance flags per measurement, informing on the quality associated with each observation, determined by multiple scientific quality control/assurance checks.'},
    }

    return standard_data

###--------------------------------------------------------------------------------------------------###
###METADATA VARIABLE STANDARDISATIONS
###--------------------------------------------------------------------------------------------------###

#Define dictionary with details pertaining to standardised metadata variables.
#These metadata variables refers to information which provides qualitative information about measurements, typically across large swathes of time.

def get_standard_metadata(parameter_details):

    standard_metadata = {    
    #VERSION OF GLOBALLY HARMONISED OBSERVATIONAL SURFACE TREATMENT (GHOST)
    'GHOST_version':{'value':[], 'standard_name':'GHOST version', 'long_name':'Globally Harmonised Observational Surface Treatment (GHOST) version', 'units':'unitless', 'data_type':np.object, 'string_format':'upper_short', 'metadata_type':'STATION MISCELLANEOUS', 'description':'Version of the Globally Harmonised Observational Surface Treatment (GHOST).'}, 
 
    #NETWORK STATION REFERENCE 
    'station_reference':{'value':[], 'standard_name':'station reference', 'long_name':'station reference identifier', 'units':'unitless', 'data_type':np.object, 'string_format':'short', 'metadata_type':'STATION MISCELLANEOUS', 'description':'reference ID for station.'}, 

    #NETWORK PROVIDED STATION INFORMATION
    'station_timezone':                 {'value':[], 'standard_name':'station timezone',                  'long_name':'station timezone',                                'units':'unitless',                   'data_type':np.object,  'string_format':'short',        'metadata_type':'STATION MISCELLANEOUS',   'description':'Name of the local timezone that the measuring station is located in. This is automatically generated using Timezone Finder Python package (taking longitude and latitude as inputs).'},
    'latitude':                         {'value':[], 'standard_name':'latitude',                          'long_name':'latitude',                                        'units':'decimal degrees North',      'data_type':np.float64, 'string_format':np.NaN,         'metadata_type':'STATION POSITION',        'description':'Geodetic latitude of measuring instrument, in decimal degrees North.'},
    'longitude':                        {'value':[], 'standard_name':'longitude',                         'long_name':'longitude',                                       'units':'decimal degrees East',       'data_type':np.float64, 'string_format':np.NaN,         'metadata_type':'STATION POSITION',        'description':'Geodetic longitude of measuring instrument, in decimal degrees East.'},
    'altitude':                         {'value':[], 'standard_name':'altitude',                          'long_name':'altitude relative to mean sea level',             'units':'m',                          'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION POSITION',        'description':'Altitude of the ground level at the station, relative to the mean sea level datum, in metres.'},
    'sampling_height':                  {'value':[], 'standard_name':'sampling height',                   'long_name':'sampling height relative to ground',              'units':'m',                          'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION POSITION',        'description':'Height above the ground of the inlet/instrument/sampler, in metres.'},
    'measurement_altitude':             {'value':[], 'standard_name':'measurement altitude',              'long_name':'measurement altitude relative to mean sea level', 'units':'m',                          'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION POSITION',        'description':'Altitude of the inlet/instrument/sampler, relative to the mean sea level datum, in metres.'},
    'distance_to_building':             {'value':[], 'standard_name':'distance to building',              'long_name':'distance to the nearest building',                'units':'m',                          'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION MISCELLANEOUS',   'description':'Distance to the nearest building of the inlet/instrument/sampler, in metres.'},
    'distance_to_kerb':                 {'value':[], 'standard_name':'distance to kerb',                  'long_name':'distance to the street kerb',                     'units':'m',                          'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION MISCELLANEOUS',   'description':'Distance to the street kerb of the inlet/instrument/sampler, in metres.'},
    'distance_to_junction':             {'value':[], 'standard_name':'distance to junction',              'long_name':'distance to the nearest road junction',           'units':'m',                          'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION MISCELLANEOUS',   'description':'Distance to the nearest road junction of the inlet/instrument/sampler, in metres.'},
    'distance_to_source':               {'value':[], 'standard_name':'distance to source',                'long_name':'distance to the main emission source',            'units':'km',                         'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION MISCELLANEOUS',   'description':'Distance to the main emission source (see variable: "main_emission_source") of the inlet/instrument/sampler, in kilometres.'},
    'street_width':                     {'value':[], 'standard_name':'street width',                      'long_name':'width of the street',                             'units':'m',                          'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION MISCELLANEOUS',   'description':'Width of the street where measurements are being made (if applicable), in metres.'},
    'street_type':                      {'value':[], 'standard_name':'street type',                       'long_name':'type of street',                                  'units':'unitless',                   'data_type':np.object,  'string_format':'title_short',  'metadata_type':'STATION MISCELLANEOUS',   'description':'Type of street where measurements are being made (if applicable).'},
    'daytime_traffic_speed':            {'value':[], 'standard_name':'daytime traffic speed',             'long_name':'average daytime speed of passing traffic',        'units':'km hr-1',                    'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION MISCELLANEOUS',   'description':'Average daytime speed of the passing traffic where measurements are being made (if applicable), in kilometres per hour.'},
    'daily_passing_vehicles':           {'value':[], 'standard_name':'daily passing vehicles',            'long_name':'average daily number of passing vehicles',        'units':'unitless',                   'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION MISCELLANEOUS',   'description':'Average number of vehicles passing daily.'},
    'horizontal_geodetic_datum':        {'value':[], 'standard_name':'horizontal geodetic datum',         'long_name':'horizontal geodetic datum',                       'units':'unitless',                   'data_type':np.object,  'string_format':'upper_short',  'metadata_type':'STATION MISCELLANEOUS',   'description':"Name of the horizontal geodetic datum used in defining geodetic latitudes and longitudes on the Earth's surface. A geodetic datum defines the shape and size of the earth, as well as the reference point for the various coordinate systems used in mapping the earth. Datums may be global, meaning that they represent the whole Earth, or they may be local, meaning that they represent an ellipsoid best-fit to only a portion of the Earth. An example of a global datum is the World Geodetic System (WGS84), the default datum used for the Global Positioning System."}, 
    'data_level':                       {'value':[], 'standard_name':'data level',                        'long_name':'data level',                                      'units':'unitless',                   'data_type':np.object,  'string_format':'short',        'metadata_type':'STATION MISCELLANEOUS',   'description':'Data level of data reported. This varies with network.'},
    'climatology':                      {'value':[], 'standard_name':'climatology',                       'long_name':'climatology',                                     'units':'unitless',                   'data_type':np.object,  'string_format':'upper_short',  'metadata_type':'STATION MISCELLANEOUS',   'description':'Name of the climatology of which the observations pertain to.'},
    'station_name':                     {'value':[], 'standard_name':'station name',                      'long_name':'station name',                                    'units':'unitless',                   'data_type':np.object,  'string_format':'title_short',  'metadata_type':'STATION MISCELLANEOUS',   'description':'Name of station where the measurement was conducted.'},
    'city':                             {'value':[], 'standard_name':'city',                              'long_name':'city',                                            'units':'unitless',                   'data_type':np.object,  'string_format':'title_short',  'metadata_type':'STATION MISCELLANEOUS',   'description':'Name of the city the station is located in.'},
    'country':                          {'value':[], 'standard_name':'country',                           'long_name':'country',                                         'units':'unitless',                   'data_type':np.object,  'string_format':'title_short',  'metadata_type':'STATION MISCELLANEOUS',   'description':'Name of the country the station is located in.'},
    'administrative_country_division_1':{'value':[], 'standard_name':'administrative country division 1', 'long_name':'administrative country division 1',               'units':'unitless',                   'data_type':np.object,  'string_format':'title_short',  'metadata_type':'STATION MISCELLANEOUS',   'description':'Name of the first (i.e. largest) country administrative division in which the station lies, e.g. countries within soverign state, state, province, county etc. These are defined for the purposes of managing of land and the affairs of people. This is automatically generated using Reverse Geocoder Python package (taking longitude and latitude as inputs).'},
    'administrative_country_division_2':{'value':[], 'standard_name':'administrative country division 2', 'long_name':'administrative country division 2',               'units':'unitless',                   'data_type':np.object,  'string_format':'title_short',  'metadata_type':'STATION MISCELLANEOUS',   'description':'Name of the second (i.e. second largest) country administrative division in which the station lies, e.g. countries within soverign state, state, province, county etc. These are defined for the purposes of managing of land and the affairs of people. This is automatically generated using Reverse Geocoder Python package (taking longitude and latitude as inputs).'},
    'population':                       {'value':[], 'standard_name':'population',                        'long_name':'population',                                      'units':'unitless',                   'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION MISCELLANEOUS',   'description':'Population size of the nearest urban settlement.'},
    'representative_radius':            {'value':[], 'standard_name':'representative_radius',             'long_name':'representative_radius',                           'units':'km',                         'data_type':np.float32, 'string_format':np.NaN,         'metadata_type':'STATION CLASSIFICATIONS', 'description':'Radius of representativity of the measurements made (i.e. for what distance scale around the sampling point would the measurements be very similar?), given in kilometres.  A quantitative version of the "measurement_scale" classification.'},
    'network':                          {'value':[], 'standard_name':'network',                           'long_name':'data network',                                    'units':'unitless',                   'data_type':np.object,  'string_format':'short',        'metadata_type':'STATION MISCELLANEOUS',   'description':'The name of the network which reports data for the specific station in question.'},
    'associated_networks':              {'value':[], 'standard_name':'other associated networks',         'long_name':'other associated networks',                       'units':'unitless',                   'data_type':np.object,  'string_format':'long',         'metadata_type':'STATION MISCELLANEOUS',   'description':'String pair of associated network name and station reference. Format: network1:station_reference1;network2:station_reference2'},

    #NETWORK PROVIDED STANDARDISED CLASSIFICATIONS
    'standardised_network_provided_area_classification':   {'value':[], 'standard_name':'standardised network provided area classification',    'long_name':'standardised network provided area classification',    'units':'unitless', 'data_type':np.object, 'string_format':'lower_short', 'metadata_type':'STATION CLASSIFICATIONS', 'description':'Standardised network provided classification, describing type of area a measurement station is situated in.'},                
    'standardised_network_provided_station_classification':{'value':[], 'standard_name':'standardised network provided station classification', 'long_name':'standardised network provided station classification', 'units':'unitless', 'data_type':np.object, 'string_format':'lower_short', 'metadata_type':'STATION CLASSIFICATIONS', 'description':'Standardised network provided classification, categorising the type of air measured by a station.'},                
    'standardised_network_provided_main_emission_source':  {'value':[], 'standard_name':'standardised network provided main emission source',   'long_name':'standardised network provided main emission source',   'units':'unitless', 'data_type':np.object, 'string_format':'lower_short', 'metadata_type':'STATION CLASSIFICATIONS', 'description':'Standardised network provided classification, describing the main emission source influencing air measured at a station.'},                
    'standardised_network_provided_land_use':              {'value':[], 'standard_name':'standardised network provided land use',               'long_name':'standardised network provided land use type',          'units':'unitless', 'data_type':np.object, 'string_format':'lower_short', 'metadata_type':'STATION CLASSIFICATIONS', 'description':'Standardised network provided classification, describing the dominant land use in the area of the reporting station.'},                
    'standardised_network_provided_terrain':               {'value':[], 'standard_name':'standardised network provided terrain',                'long_name':'standardised network provided terrain type',           'units':'unitless', 'data_type':np.object, 'string_format':'lower_short', 'metadata_type':'STATION CLASSIFICATIONS', 'description':'Standardised network provided classification, describing the dominant terrain in the area of the reporting station.'},                 
    'standardised_network_provided_measurement_scale':     {'value':[], 'standard_name':'standardised network provided measurement scale',      'long_name':'standardised network provided measurement scale',      'units':'unitless', 'data_type':np.object, 'string_format':'lower_short', 'metadata_type':'STATION CLASSIFICATIONS', 'description':'Standardised network provided classification, a denotation of the geographic scope of the air quality measurements made.'},                         

    #GLOBALLY GRIDDED CLASSIFICATIONS
    'ESDAC_Iwahashi_landform_classification':           {'value':[], 'standard_name':'ESDAC Iwahashi landform classification',            'long_name':'ESDAC Iwahashi landform classification',                      'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'European Soil Data Centre (ESDAC) Iwahashi landform classification. The classification presents relief classes which are classified using an unsupervised nested-means algorithms and a three part geometric signature. Slope gradient, surface texture and local convexity are calculated based on the SRTM30 digital elevation model, within a given window size and classified according to the inherent data set properties. This is a dynamic landform classification method. Native resolution of 0.0083 x 0.0083 degrees. A correction for costal sites is made: if the native class is "water", then the modal classification of the neighbouring grid boxes is used instead (lowest code kept preferentially in case of a tie). If the site is truly an "ocean" site, all the surrounding gridcells will be water also, and therefore the class will be maintained as "water".'},
    'ESDAC_modal_Iwahashi_landform_classification_5km': {'value':[], 'standard_name':'ESDAC modal Iwahashi landform classification 5km',  'long_name':'ESDAC modal Iwahashi landform classification in 5km radius',  'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal European Soil Data Centre (ESDAC) Iwahashi landform classification in radius of 5km around station location.'},
    'ESDAC_modal_Iwahashi_landform_classification_25km':{'value':[], 'standard_name':'ESDAC modal Iwahashi landform classification 25km', 'long_name':'ESDAC modal Iwahashi landform classification in 25km radius', 'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal European Soil Data Centre (ESDAC) Iwahashi landform classification in radius of 25km around station location.'},
    'ESDAC_Meybeck_landform_classification':            {'value':[], 'standard_name':'ESDAC Meybeck landform classification',             'long_name':'ESDAC Meybeck landform classification',                       'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'European Soil Data Centre (ESDAC) Meybeck landform classification. The classification presents relief classes which are calculated based on the relief roughness. Roughness and elevation are classified based on a digital elevation model according to static thresholds, with a given window size. This is a static landform classification method. Native resolution of 0.0083 x 0.0083 degrees. A correction for costal sites is made: if the native class is "water", then the modal classification of the neighbouring grid boxes is used instead (lowest code kept preferentially in case of a tie). If the site is truly an "ocean" site, all the surrounding gridcells will be water also, and therefore the class will be maintained as "water".'},
    'ESDAC_modal_Meybeck_landform_classification_5km':  {'value':[], 'standard_name':'ESDAC modal Meybeck landform classification 5km',   'long_name':'ESDAC modal Meybeck landform classification in 5km radius',   'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal European Soil Data Centre (ESDAC) Meybeck landform classification in radius of 5km around station location.'},
    'ESDAC_modal_Meybeck_landform_classification_25km': {'value':[], 'standard_name':'ESDAC modal Meybeck landform classification 25km',  'long_name':'ESDAC modal Meybeck landform classification in 25km radius',  'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal European Soil Data Centre (ESDAC) Meybeck landform classification in radius of 25km around station location.'},
    'GHSL_settlement_model_classification':             {'value':[], 'standard_name':'GHSL settlement model classification',              'long_name':'GHSL settlement model classification',                        'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) settlement model classification (technical label: GHS_SMOD_POPMT_GLOBE_R2019A). The classification delineates and classify settlement typologies via a logic of population size, population and built-up area densities as a refinement of the ‘degree of urbanization’ method described by EUROSTAT. The classification is derived by using the GHS_POP_MT_GLOBE_R2019A and GHS_BUILT_LDSMT_GLOBE_R2018A products. The GHS Settlement Model grid is an improvement of the GHS Settlement Grid (R2016A) introducing a more detailed classification of settlements in two levels, also called ‘refined degree of urbanization’. The Settlement Model is provided at detailed level (Second Level - L2). The First Level, as a porting of the Degree of Urbanization adopted by EUROSTAT can be obtained aggregating L2. Native resolution of 1.0 x 1.0 kilometres.'},   
    'GHSL_modal_settlement_model_classification_5km':   {'value':[], 'standard_name':'GHSL modal settlement model classification 5km',    'long_name':'GHSL modal settlement model classification in 5km radius',    'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal Global Human Settlement Layer (GHSL) settlement model classification in radius of 5km around station location.'},   
    'GHSL_modal_settlement_model_classification_25km':  {'value':[], 'standard_name':'GHSL modal settlement model classification 25km',   'long_name':'GHSL modal settlement model classification in 25km radius',   'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal Global Human Settlement Layer (GHSL) settlement model classification in radius of 25km around station location.'},   
    'Joly-Peuch_classification_code':                   {'value':[], 'standard_name':'Joly-Peuch classification code',                    'long_name':'Joly-Peuch classification code',                              'units':'unitless', 'data_type':np.float32, 'string_format':np.NaN,        'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Joly-Peuch European classification code (range of 1-10) designed to objectively stratify stations between those diplaying rural and urban signatures (most rural == 1, most urban == 10). This classification is objectively made per species. The species that this is done for are: O3, NO2, SO2, CO, PM10, PM2.5. See reference here: https://www.sciencedirect.com/science/article/abs/pii/S1352231011012088'},        
    'Koppen-Geiger_classification':                     {'value':[], 'standard_name':'Koppen-Geiger classification',                      'long_name':'Koppen-Geiger classification',                                'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Koppen-Geiger classification, classifying the global climates into 5 main groups (30 total groups with subcategories). Native resolution of 0.0083 x 0.0083 degrees. A correction for costal sites is made: if the native class is "water", then the modal classification of the neighbouring grid boxes is used instead (lowest code kept preferentially in case of a tie). If the site is truly an "ocean" site, all the surrounding gridcells will be water also, and therefore the class will be maintained as "water". See citation: Beck, H.E., N.E. Zimmermann, T.R. McVicar, N. Vergopolan, A. Berg, E.F. Wood: Present and future Köppen-Geiger climate classification maps at 1-km resolution, Nature Scientific Data, 2018.'},
    'Koppen-Geiger_modal_classification_5km':           {'value':[], 'standard_name':'Koppen-Geiger modal classification 5km',            'long_name':'Koppen-Geiger classification',                                'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal Koppen-Geiger classification in radius of 5km around station location.'},
    'Koppen-Geiger_modal_classification_25km':          {'value':[], 'standard_name':'Koppen-Geiger modal classification 25km',           'long_name':'Koppen-Geiger classification',                                'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal Koppen-Geiger classification in radius of 25km around station location.'},
    'MODIS_MCD12C1_v6_IGBP_land_use':                   {'value':[], 'standard_name':'MODIS MCD12C1 v6 IGBP land use',                    'long_name':'MODIS MCD12C1 v6 IGBP land use',                              'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Majority land use class from the Moderate Resolution Imaging Spectroradiometer (MODIS) Land Cover Climate Modeling Grid (CMG) MCD12C1 version 6, using the International Geosphere-Biosphere Programme (IGBP) classification. Native resolution of 0.05 x 0.05 degrees. See dataset user guide here: https://lpdaac.usgs.gov/documents/101/MCD12_User_Guide_V6.pd. A correction for costal sites is made: if the native class is "water bodies", then the modal classification of the neighbouring grid boxes is used instead (lowest code kept preferentially in case of a tie). If the site is truly an "ocean" site, all the surrounding gridcells will be water also, and therefore the class will be maintained as "water".'},
    'MODIS_MCD12C1_v6_modal_IGBP_land_use_5km':         {'value':[], 'standard_name':'MODIS MCD12C1 v6 IGBP modal land use 5km',          'long_name':'MODIS MCD12C1 v6 IGBP modal land use in 5km radius',          'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal land use in radius of 5km around the station location from the Moderate Resolution Imaging Spectroradiometer (MODIS) Land Cover Climate Modeling Grid (CMG) MCD12C1 version 6, using the International Geosphere-Biosphere Programme (IGBP) classification.'},
    'MODIS_MCD12C1_v6_modal_IGBP_land_use_25km':        {'value':[], 'standard_name':'MODIS MCD12C1 v6 IGBP modal land use 25km',         'long_name':'MODIS MCD12C1 v6 IGBP modal land use in 25km radius',         'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal land use in radius of 25km around the station location from the Moderate Resolution Imaging Spectroradiometer (MODIS) Land Cover Climate Modeling Grid (CMG) MCD12C1 version 6, using the International Geosphere-Biosphere Programme (IGBP) classification.'},
    'MODIS_MCD12C1_v6_UMD_land_use':                    {'value':[], 'standard_name':'MODIS MCD12C1 v6 UMD land use',                     'long_name':'MODIS MCD12C1 v6 UMD land use',                               'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Majority land use class from the Moderate Resolution Imaging Spectroradiometer (MODIS) Land Cover Climate Modeling Grid (CMG) MCD12C1 version 6, using the University of Maryland (UMD) classification. Native resolution of 0.05 x 0.05 degrees. See dataset user guide here: https://lpdaac.usgs.gov/documents/101/MCD12_User_Guide_V6.pd. A correction for costal sites is made: if the native class is "water bodies", then the modal classification of the neighbouring grid boxes is used instead (lowest code kept preferentially in case of a tie). If the site is truly an "ocean" site, all the surrounding gridcells will be water also, and therefore the class will be maintained as "water".'},
    'MODIS_MCD12C1_v6_modal_UMD_land_use_5km':          {'value':[], 'standard_name':'MODIS MCD12C1 v6 UMD modal land use 5km',           'long_name':'MODIS MCD12C1 v6 UMD modal land use in 5km radius',           'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal land use in radius of 5km around the station location from the Moderate Resolution Imaging Spectroradiometer (MODIS) Land Cover Climate Modeling Grid (CMG) MCD12C1 version 6, using the University of Maryland (UMD) classification.'},
    'MODIS_MCD12C1_v6_modal_UMD_land_use_25km':         {'value':[], 'standard_name':'MODIS MCD12C1 v6 UMD modal land use 25km',          'long_name':'MODIS MCD12C1 v6 UMD modal land use in 25km radius',          'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal land use in radius of 25km around the station location from the Moderate Resolution Imaging Spectroradiometer (MODIS) Land Cover Climate Modeling Grid (CMG) MCD12C1 version 6, using the University of Maryland (UMD) classification.'},
    'MODIS_MCD12C1_v6_LAI':                             {'value':[], 'standard_name':'MODIS MCD12C1 v6 LAI',                              'long_name':'MODIS MCD12C1 v6 Leaf Area Index',                            'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Majority Leaf Area Index class from the Moderate Resolution Imaging Spectroradiometer (MODIS) Land Cover Climate Modeling Grid (CMG) MCD12C1 version 6. Native resolution of 0.05 x 0.05 degrees. See dataset user guide here: https://lpdaac.usgs.gov/documents/101/MCD12_User_Guide_V6.pd. A correction for costal sites is made: if the native class is "water bodies", then the modal classification of the neighbouring grid boxes is used instead (lowest code kept preferentially in case of a tie). If the site is truly an "ocean" site, all the surrounding gridcells will be water also, and therefore the class will be maintained as "water".'},
    'MODIS_MCD12C1_v6_modal_LAI_5km':                   {'value':[], 'standard_name':'MODIS MCD12C1 v6 LAI modal 5km',                    'long_name':'MODIS MCD12C1 v6 modal Leaf Area Index in 5km radius',        'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal Leaf Area Index in radius of 5km around the station location from the Moderate Resolution Imaging Spectroradiometer (MODIS) Land Cover Climate Modeling Grid (CMG) MCD12C1 version 6.'},
    'MODIS_MCD12C1_v6_modal_LAI_25km':                  {'value':[], 'standard_name':'MODIS MCD12C1 v6 LAI modal 25km',                   'long_name':'MODIS MCD12C1 v6 modal Leaf Area Index in 25km radius',       'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Modal Leaf Area Index in radius of 25km around the station location from the Moderate Resolution Imaging Spectroradiometer (MODIS) Land Cover Climate Modeling Grid (CMG) MCD12C1 version 6.'},
    'WMO_region':                                       {'value':[], 'standard_name':'WMO region',                                        'long_name':'WMO region code',                                             'units':'unitless', 'data_type':np.object,  'string_format':'title_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'World Meteorological Organization (WMO) region of station. The available regions are: Africa, Asia, South America, "Northern America, Central America and the Caribbean", South-West Pacific, Europe and Antarctica.'},
    'WWF_TEOW_terrestrial_ecoregion':                   {'value':[], 'standard_name':'WWF TEOW terrestrial ecoregion',                    'long_name':'WWF TEOW terrestrial ecoregion',                              'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Terrestrial Ecoregions of the World (TEOW) World Wildlife Foundation (WWF) classification. There are 825 terrestrial ecoregions. Ecoregions are relatively large units of land containing distinct assemblages of natural communities and species, with boundaries that approximate the original extent of natural communities prior to major land-use change. See citation: Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N., Underwood, E. C., DAmico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth. Bioscience 51(11):933-938.'},
    'WWF_TEOW_biogeographical_realm':                   {'value':[], 'standard_name':'WWF TEOW biogeographical realm',                    'long_name':'WWF TEOW biogeographical realm',                              'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Terrestrial Ecoregions of the World (TEOW) World Wildlife Foundation (WWF) classification. There are 8 biogeographical realms. See citation: Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N., Underwood, E. C., DAmico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth. Bioscience 51(11):933-938.'},
    'WWF_TEOW_biome':                                   {'value':[], 'standard_name':'WWF TEOW biome',                                    'long_name':'WWF TEOW biome',                                              'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Terrestrial Ecoregions of the World (TEOW) World Wildlife Foundation (WWF) classification. There are 14 biomes. See citation: Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N., Underwood, E. C., DAmico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth. Bioscience 51(11):933-938.'},    
    'UMBC_anthrome_classification':                     {'value':[], 'standard_name':'UMBC anthrome classification',                      'long_name':'UMBC anthrome classification',                                'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'University of Maryland Baltimore County (UMBC) anthrome classification, describing the anthropogenic land use (for the year 2000). There are 20 distinct classifications. Native resolution of 0.0833 x 0.0833 degrees. A correction for costal sites is made: if the native anthrome class is "water", then the modal classification of the neighbouring grid boxes is used instead (lowest code kept preferentially in case of a tie). If the site is truly an "ocean" site, all the surrounding gridcells will be water also, and therefore the class will be maintained as "water".'},
    'UMBC_modal_anthrome_classification_5km':           {'value':[], 'standard_name':'UMBC modal anthrome classification 5km',            'long_name':'UMBC modal anthrome classification in 5km radius',            'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'University of Maryland Baltimore County (UMBC) modal anthrome classification in radius of 5km around station location.'},
    'UMBC_modal_anthrome_classification_25km':          {'value':[], 'standard_name':'UMBC modal anthrome classification 25km',           'long_name':'UMBC modal anthrome classification in 25km radius',           'units':'unitless', 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'University of Maryland Baltimore County (UMBC) modal anthrome classification in radius of 25km around station location.'},
          
    #GLOBALLY GRIDDED PRODUCTS
    'EDGAR_v4.3.2_annual_average_BC_emissions':                        {'value':[], 'standard_name':'EDGAR v4.3.2 annual average BC emissions',                         'long_name':'EDGAR v4.3.2 annual average black carbon emissions',                                                          'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average BC emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'EDGAR_v4.3.2_annual_average_CO_emissions':                        {'value':[], 'standard_name':'EDGAR v4.3.2 annual average CO emissions',                         'long_name':'EDGAR v4.3.2 annual average carbon monoxide emissions',                                                       'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average CO emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'EDGAR_v4.3.2_annual_average_NH3_emissions':                       {'value':[], 'standard_name':'EDGAR v4.3.2 annual average NH3 emissions',                        'long_name':'EDGAR v4.3.2 annual average ammonia emissions',                                                               'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average NH3 emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'EDGAR_v4.3.2_annual_average_NMVOC_emissions':                     {'value':[], 'standard_name':'EDGAR v4.3.2 annual average NMVOC emissions',                      'long_name':'EDGAR v4.3.2 annual average non-methane volatile organic compound emissions',                                 'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average NMVOC emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'EDGAR_v4.3.2_annual_average_NOx_emissions':                       {'value':[], 'standard_name':'EDGAR v4.3.2 annual average NOx emissions',                        'long_name':'EDGAR v4.3.2 annual average emissions of nitrogen oxides',                                                    'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average NOx emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'EDGAR_v4.3.2_annual_average_OC_emissions':                        {'value':[], 'standard_name':'EDGAR v4.3.2 annual average OC emissions',                         'long_name':'EDGAR v4.3.2 annual average organic carbon emissions',                                                        'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average OC emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'EDGAR_v4.3.2_annual_average_PM10_emissions':                      {'value':[], 'standard_name':'EDGAR v4.3.2 annual average PM10 emissions',                       'long_name':'EDGAR v4.3.2 annual average PM10 emissions',                                                                  'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average PM10 emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'EDGAR_v4.3.2_annual_average_biogenic_PM2.5_emissions':            {'value':[], 'standard_name':'EDGAR v4.3.2 annual average biogenic PM2.5 emissions',             'long_name':'EDGAR v4.3.2 annual average biogenic PM2.5 emissions',                                                        'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average biogenic PM2.5 emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'EDGAR_v4.3.2_annual_average_fossilfuel_PM2.5_emissions':          {'value':[], 'standard_name':'EDGAR v4.3.2 annual average fossil fuel PM2.5 emissions',          'long_name':'EDGAR v4.3.2 annual average fossil fuel PM2.5 emissions',                                                     'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average fossil fuel PM2.5 emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'EDGAR_v4.3.2_annual_average_SO2_emissions':                       {'value':[], 'standard_name':'EDGAR v4.3.2 annual average SO2 emissions',                        'long_name':'EDGAR v4.3.2 annual average sulphur dioxide emissions',                                                       'units':'kg m-2 s-1', 'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'EDGAR v4.3.2 annual average SO2 emissions, in kilograms per squared metre per second. Native resolution of 0.1 x 0.1 degrees.'},         
    'ASTER_v3_altitude':                                               {'value':[], 'standard_name':'ASTER v3 altitude',                                                'long_name':'ASTER v3 altitude, relative to EGM96 geoid datum',                                                            'units':'m',          'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Altitude from ASTER v3 digital elevation model, relative to EGM96 geoid datum, in metres. The dataset was generated using 1,880,306 Level-1A scenes (taken from the NASA TERRA spacecraft) acquired between March 1, 2000 and November 30, 2013. The ASTER GDEM was created by stacking all individual cloud-masked scene DEMs and non-cloud-masked scene DEMs, then applying various algorithms to remove abnormal data. A statistical approach is not always effective for anomaly removal in areas with a limited number of images. Several existing reference DEMs were used to replace residual anomalies caused by the insufficient number of stacked scenes. In addition to ASTER GDEM, the ASTER Global Water Body Database (ASTWBD) was generated as a by-product to correct elevation values of water body surfaces like sea, rivers, and lakes. The ASTWBD was applied to GDEM to provide proper elevation values for water body surfaces. The sea and lake have a flattened elevation value. The river has a stepped-down elevation value from the upper stream to the lower stream. Native resolution of 30m x 30m.'},         
    'ETOPO1_altitude':                                                 {'value':[], 'standard_name':'ETOPO1 altitude',                                                  'long_name':'ETOPO1 altitude, relative to sea level datum',                                                                'units':'m',          'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Altitude from ETOPO1 digital elevation model, relative to sea level datum, in metres. Native resolution of 0.0166 x 0.0166 degrees. A correction for coastal sites is made: if the derived altitude is <= -5 m, the maximum altitude of the neighbouring grid boxes will be used instead. If all neighbouring grid boxes have altitudes <= -5 m, the original value will be retained.'},         
    'ETOPO1_max_altitude_difference_5km':                              {'value':[], 'standard_name':'ETOPO1 max altitude difference 5km',                               'long_name':'ETOPO1 maximum altitude difference between the ETOPO1_altitude and all ETOPO1 altitudes in 5km radius',       'units':'m',          'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Altitude difference between the ETOPO1_altitude, and the minimum ETOP1 altitude in a radius of 5 km around the station location, in metres.'},         
    'GHSL_built_up_area_density':                                      {'value':[], 'standard_name':'GHSL built up area density',                                       'long_name':'GHSL built up area density',                                                                                  'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) built up area density (technical label: GHS_BUILT_LDSMT_GLOBE_R2018A), in units of built-up area fraction per gridcell (0-100). The product is a multitemporal information layer on built-up presence as derived from Landsat image collections (GLS1975, GLS1990, GLS2000, and ad-hoc Landsat 8 collection 2013/2014). Native resolution of 0.25 x 0.25 kilometres.'},   
    'GHSL_average_built_up_area_density_5km':                          {'value':[], 'standard_name':'GHSL average built up area density 5km',                           'long_name':'GHSL average built up area density in 5km radius',                                                            'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) average built up area density in a radius of 5km around the station location.'},   
    'GHSL_average_built_up_area_density_25km':                         {'value':[], 'standard_name':'GHSL average built up area density 25km',                          'long_name':'GHSL average built up area density in 25km radius',                                                           'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) average built up area density in a radius of 25km around the station location.'},   
    'GHSL_max_built_up_area_density_5km':                              {'value':[], 'standard_name':'GHSL max built up area density 5km',                               'long_name':'GHSL max built up area density in 5km radius',                                                                'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) max built up area density in a radius of 5km around the station location.'},   
    'GHSL_max_built_up_area_density_25km':                             {'value':[], 'standard_name':'GHSL max built up area density 25km',                              'long_name':'GHSL max built up area density in 25km radius',                                                               'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) max built up area density in a radius of 25km around the station location.'},   
    'GHSL_population_density':                                         {'value':[], 'standard_name':'GHSL population density',                                          'long_name':'GHSL population density',                                                                                     'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) population density (technical label: GHS_POP_MT_GLOBE_R2019A), in populus per squared kilometre. It depicts the distribution of population, expressed as the number of people per cell. Residential population estimates for target years 1975, 1990, 2000 and 2015 provided by CIESIN GPWv4.10 were disaggregated from census or administrative units to grid cells, informed by the distribution and density of built-up as mapped in the GHSL global layer per corresponding epoch. Native resolution of 0.25 x 0.25 kilometres.'},   
    'GHSL_average_population_density_5km':                             {'value':[], 'standard_name':'GHSL average population density 5km',                              'long_name':'GHSL average population density in 5km radius',                                                               'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) average population density in a radius of 5km around the station location.'},   
    'GHSL_average_population_density_25km':                            {'value':[], 'standard_name':'GHSL average population density 25km',                             'long_name':'GHSL average population density in 25km radius',                                                              'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) average population density in a radius of 25km around the station location.'},   
    'GHSL_max_population_density_5km':                                 {'value':[], 'standard_name':'GHSL max population density 5km',                                  'long_name':'GHSL max population density in 5km radius',                                                                   'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) max population density in a radius of 5km around the station location.'},   
    'GHSL_max_population_density_25km':                                {'value':[], 'standard_name':'GHSL max population density 25km',                                 'long_name':'GHSL max population density in 25km radius',                                                                  'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Global Human Settlement Layer (GHSL) max population density in a radius of 25km around the station location.'},   
    'GPW_population_density':                                          {'value':[], 'standard_name':'GPW population density',                                           'long_name':'GPW population density',                                                                                      'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Gridded Population of the World (GPW), population density, in populus per squared kilometre, from either version 3 and 4 of the provided gridded datasets, dependent on the data year: v3 (1990-2000), v4 (2000-2015). Native resolution of 0.04166 x 0.04166 for v3 data; native resolution of 0.0083 x 0.0083 degrees for v4 data.'},         
    'GPW_average_population_density_5km':                              {'value':[], 'standard_name':'GPW average population density 5km',                               'long_name':'GPW average population density in 5km radius',                                                                'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Gridded Population of the World (GPW), average population density in a radius of 5 km around the station location.'},         
    'GPW_average_population_density_25km':                             {'value':[], 'standard_name':'GPW average population density 25km',                              'long_name':'GPW average population density in 25km radius',                                                               'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Gridded Population of the World (GPW), average population density in a radius of 25 km around the station location.'},         
    'GPW_max_population_density_5km':                                  {'value':[], 'standard_name':'GPW max population density 5km',                                   'long_name':'GPW average population density in 5km radius',                                                                'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Gridded Population of the World (GPW), maximum population density in a radius of 5 km around the station location.'},         
    'GPW_max_population_density_25km':                                 {'value':[], 'standard_name':'GPW max population density 25km',                                  'long_name':'GPW average population density in 25km radius',                                                               'units':'xx km–2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Gridded Population of the World (GPW), maximum population density in a radius of 25 km around the station location.'},         
    'NOAA-DMSP-OLS_v4_nighttime_stable_lights':                        {'value':[], 'standard_name':'NOAA-DMSP-OLS v4 nighttime stable lights',                         'long_name':'NOAA DMSP-OLS version 4 nighttime stable lights',                                                             'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'National Oceanic and Atmospheric Administration (NOAA), Defense Meteorological Satellite Program - Operational Linescane System (DMSP-OLS) version 4 nighttime stable lights. Native resolution of 0.0083 x 0.0083 degrees. The values represent a brightness index ranging from 0 to 63. The sensor saturates at a value of 63.'},         
    'NOAA-DMSP-OLS_v4_average_nighttime_stable_lights_5km':            {'value':[], 'standard_name':'NOAA-DMSP-OLS v4 average nighttime stable lights 5km',             'long_name':'NOAA DMSP-OLS version 4 average nighttime stable lights in 5km radius',                                       'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'National Oceanic and Atmospheric Administration (NOAA), Defense Meteorological Satellite Program - Operational Linescane System (DMSP-OLS) version 4 average nighttime stable lights in 5km radius of measurement station. The values represent a brightness index ranging from 0 to 63. The sensor saturates at a value of 63.'},         
    'NOAA-DMSP-OLS_v4_average_nighttime_stable_lights_25km':           {'value':[], 'standard_name':'NOAA-DMSP-OLS v4 average nighttime stable lights 25km',            'long_name':'NOAA DMSP-OLS version 4 average nighttime stable lights in 25km radius',                                      'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'National Oceanic and Atmospheric Administration (NOAA), Defense Meteorological Satellite Program - Operational Linescane System (DMSP-OLS) version 4 average nighttime stable lights in 25km radius of measurement station. The values represent a brightness index ranging from 0 to 63. The sensor saturates at a value of 63.'},         
    'NOAA-DMSP-OLS_v4_max_nighttime_stable_lights_5km':                {'value':[], 'standard_name':'NOAA-DMSP-OLS v4 max nighttime stable lights 5km',                 'long_name':'NOAA DMSP-OLS version 4 maximum nighttime stable lights in 5km radius',                                       'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'National Oceanic and Atmospheric Administration (NOAA), Defense Meteorological Satellite Program - Operational Linescane System (DMSP-OLS) version 4 maximum nighttime stable lights in 5km radius of measurement station. The values represent a brightness index ranging from 0 to 63. The sensor saturates at a value of 63.'},         
    'NOAA-DMSP-OLS_v4_max_nighttime_stable_lights_25km':               {'value':[], 'standard_name':'NOAA-DMSP-OLS v4 max nighttime stable lights 25km',                'long_name':'NOAA DMSP-OLS version 4 maximum nighttime stable lights in 25km radius',                                      'units':'unitless',   'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'National Oceanic and Atmospheric Administration (NOAA), Defense Meteorological Satellite Program - Operational Linescane System (DMSP-OLS) version 4 maximum nighttime stable lights in 25km radius of measurement station. The values represent a brightness index ranging from 0 to 63. The sensor saturates at a value of 63.'},         
    'OMI_level3_column_annual_average_NO2':                            {'value':[], 'standard_name':'OMI level3 column annual average NO2',                             'long_name':'OMI level3 column annual average nitrogen dioxide',                                                           'units':'xx cm-2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'AURA Ozone monitoring instrument (OMI) level3 column annual average NO2, in molecules per squared centimetres. Native resolution of 0.25 x 0.25 degrees.'},         
    'OMI_level3_column_cloud_screened_annual_average_NO2':             {'value':[], 'standard_name':'OMI level3 column cloud screened annual average NO2',              'long_name':'OMI level3 column cloud screened annual average nitrogen dioxide',                                            'units':'xx cm-2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'AURA Ozone monitoring instrument (OMI) level3 column cloud screened (where cloud fraction is less than 30 percent) annual average NO2, in molecules per squared centimetres. Native resolution of 0.25 x 0.25 degrees.'},         
    'OMI_level3_tropospheric_column_annual_average_NO2':               {'value':[], 'standard_name':'OMI level3 tropospheric column annual average NO2',                'long_name':'OMI level3 tropospheric column annual average nitrogen dioxide',                                              'units':'xx cm-2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'AURA Ozone monitoring instrument (OMI) level3 tropospheric column annual average NO2, in molecules per squared centimetres. Native resolution of 0.25 x 0.25 degrees.'},         
    'OMI_level3_tropospheric_column_cloud_screened_annual_average_NO2':{'value':[], 'standard_name':'OMI level3 tropospheric column cloud screened annual average NO2', 'long_name':'OMI level3 tropospheric column cloud screened annual average nitrogen dioxide',                               'units':'xx cm-2',    'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'AURA Ozone monitoring instrument (OMI) level3 tropospheric column cloud screened (where cloud fraction is less than 30 percent) annual average NO2, in molecules per squared centimetres. Native resolution of 0.25 x 0.25 degrees.'},         
    'GSFC_coastline_proximity':                                        {'value':[], 'standard_name':'GSFC coastline proximity',                                         'long_name':'GSFC proximity to the coastline',                                                                             'units':'km',         'data_type':np.float32, 'string_format':np.NaN, 'metadata_type':'GLOBALLY GRIDDED CLASSIFICATIONS', 'description':'Proximity to the coastline provided by the NASA Goddard Space Flight Center (GSFC) Ocean Color Group, in kilometres, produced using the Generic Mapping Tools package. Native resolution of 0.01 x 0.01 degrees. Negative distances represent locations over land (including land-locked bodies of water), while positive distances represent locations over the ocean. There is an uncertainty of up to 1 km in the computed distance at any given point.'},         

    #MEASUREMENT INFORMATION
    'primary_sampling_type':                                   {'value':[], 'standard_name':'primary sampling type',                                    'long_name':'standardised primary sampling type',                                    'units':'unitless',                          'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Standardised primary sampling type.'},
    'primary_sampling_instrument_name':                        {'value':[], 'standard_name':'primary sampling instrument name',                         'long_name':'standardised primary sampling instrument name',                         'units':'unitless',                          'data_type':np.object,  'string_format':'short',       'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Standardised name of the primary sampling instrument (if no specific instrument is used, or known, this is the standardised primary sampling type).'},  
    'primary_sampling_instrument_documented_flow_rate':        {'value':[], 'standard_name':'primary sampling instrument documented flow rate',         'long_name':'primary sampling instrument documented sampling flow rate',             'units':'l min-1',                           'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Volume (litres) of fluid which passes to the primary sampling instrument, per unit time (minutes), as given in instrumental manual/documentation. Can be a range: e.g. 1.0-3.0.'},  
    'primary_sampling_instrument_reported_flow_rate':          {'value':[], 'standard_name':'primary sampling instrument reported flow rate',           'long_name':'primary sampling instrument reported sampling flow rate',               'units':'l min-1',                           'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Volume (litres) of fluid which passes to the primary sampling instrument, per unit time (minutes), as given in metadata. Can be a range: e.g. 1.0-3.0.'},
    'primary_sampling_process_details':                        {'value':[], 'standard_name':'primary sampling process details',                         'long_name':'primary sampling process details',                                      'units':'unitless',                          'data_type':np.object,  'string_format':'lower_long',  'metadata_type':np.NaN,                            'description':'Miscellaneous details regarding assumptions made in the standardisation of the primary sampling type/instrument.'},
    'primary_sampling_instrument_manual_name':                 {'value':[], 'standard_name':'primary sampling instrument manual name',                  'long_name':'primary sampling instrument manual name',                               'units':'unitless',                          'data_type':np.object,  'string_format':'short',       'metadata_type':np.NaN,                            'description':'Path to the location in the esarchive of the manual for the specific primary sampling instrument.'},
    'primary_sampling_further_details':                        {'value':[], 'standard_name':'primary sampling further details',                         'long_name':'primary sampling further details',                                      'units':'unitless',                          'data_type':np.object,  'string_format':'lower_long',  'metadata_type':np.NaN,                            'description':'Further associated details regarding the specifics of the primary sampling instrument/type.'},
    'sample_preparation_types':                                {'value':[], 'standard_name':'sample preparation types',                                 'long_name':'standardised sample preparation types',                                 'units':'unitless',                          'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Standardised sample preparation types utilised in the measurement process. Mutiple types are separated by ";".'},               
    'sample_preparation_techniques':                           {'value':[], 'standard_name':'sample preparation techniques',                            'long_name':'standardised specific preparation techniques',                          'units':'unitless',                          'data_type':np.object,  'string_format':'short',       'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Standardised sample preparation techniques utilised in the measurement process. Mutiple names are separated by ";".'},               
    'sample_preparation_process_details':                      {'value':[], 'standard_name':'sample preparation process details',                       'long_name':'sample preparation process details',                                    'units':'unitless',                          'data_type':np.object,  'string_format':'lower_long',  'metadata_type':np.NaN,                            'description':'Miscellaneous details regarding assumptions made in the standardisation of the sample preparation types/techniques. Multiple details specific to different types are separated by ";".'},
    'sample_preparation_further_details':                      {'value':[], 'standard_name':'sample preparation further details',                       'long_name':'sample preparation further details',                                    'units':'unitless',                          'data_type':np.object,  'string_format':'lower_long',  'metadata_type':np.NaN,                            'description':'Further associated details regarding the specifics of the sample preparation types/instruments. Multiple details specific to different types are separated by ";".'},
    'measurement_methodology':                                 {'value':[], 'standard_name':'measurement methodology',                                  'long_name':'standardised measurement methodology',                                  'units':'unitless',                          'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Standardised name of the measurement methodology.'},               
    'measuring_instrument_name':                               {'value':[], 'standard_name':'measuring instrument name',                                'long_name':'standardised measuring instrument name',                                'units':'unitless',                          'data_type':np.object,  'string_format':'short',       'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Standardised name of the measuring instrument.'},  
    'measuring_instrument_sampling_type':                      {'value':[], 'standard_name':'measuring instrument sampling type',                       'long_name':'standardised sampling type of the measuring instrument',                'units':'unitless',                          'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Standardised name of the measuring instrument sampling type.'},     
    'measuring_instrument_documented_flow_rate':               {'value':[], 'standard_name':'measuring instrument documented flow rate',                'long_name':'measuring instrument documented flow rate',                             'units':'l min-1',                           'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Volume (litres) of fluid which passes to the measuring instrument, per unit time (minutes), as given in instrumental manual/documentation. Can be a range: e.g. 1.0-3.0.'}, 
    'measuring_instrument_reported_flow_rate':                 {'value':[], 'standard_name':'measuring instrument reported flow rate',                  'long_name':'measuring instrument reported flow rate',                               'units':'l min-1',                           'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Volume (litres) of fluid which passes to the measuring instrument, per unit time (minutes), as given in metadata. Can be a range: e.g. 1.0-3.0.'},  
    'measuring_instrument_process_details':                    {'value':[], 'standard_name':'measuring instrument process details',                     'long_name':'measuring instrument process details',                                  'units':'unitless',                          'data_type':np.object,  'string_format':'lower_long',  'metadata_type':np.NaN,                            'description':'Miscellaneous details regarding assumptions made in the standardisation of the measurement methodology/instrument.'},
    'measuring_instrument_manual_name':                        {'value':[], 'standard_name':'measuring instrument manual name',                         'long_name':'measuring instrument manual name',                                      'units':'unitless',                          'data_type':np.object,  'string_format':'short',       'metadata_type':np.NaN,                            'description':'Path to the location in the esarchive of the manual for the specific measuring instrument.'},
    'measuring_instrument_further_details':                    {'value':[], 'standard_name':'measuring instrument further details',                     'long_name':'measuring instrument further details',                                  'units':'unitless',                          'data_type':np.object,  'string_format':'lower_long',  'metadata_type':np.NaN,                            'description':'Further associated details regarding the specifics of the measurement methodology/instrument.'},
    'measuring_instrument_reported_units':                     {'value':[], 'standard_name':'measuring instrument reported units',                      'long_name':'measuring instrument reported measurement units',                       'units':'unitless',                          'data_type':np.object,  'string_format':'short',       'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Units that the measured parameter are natively reported in.'}, 
    'measuring_instrument_reported_lower_limit_of_detection':  {'value':[], 'standard_name':'measuring instrument reported lower limit of detection',   'long_name':'measuring instrument reported lower limit of detection',                'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN,        'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Lower limit of detection of measurement methodology, as given in metadata.'},
    'measuring_instrument_documented_lower_limit_of_detection':{'value':[], 'standard_name':'measuring instrument documented lower limit of detection', 'long_name':'measuring instrument documented lower limit of detection',              'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN,        'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Lower limit of detection of measurement methodology, as given in the instrumental manual/documentation.'}, 
    'measuring_instrument_reported_upper_limit_of_detection':  {'value':[], 'standard_name':'measuring instrument reported upper limit of detection',   'long_name':'measuring instrument reported upper limit of detection',                'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN,        'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Upper limit of detection of measurement methodology, as given in metadata.'},
    'measuring_instrument_documented_upper_limit_of_detection':{'value':[], 'standard_name':'measuring instrument documented upper limit of detection', 'long_name':'measuring instrument documented upper limit of detection',              'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN,        'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Upper limit of detection of measurement methodology, as given in the instrumental manual/documentation.'}, 
    'measuring_instrument_reported_uncertainty':               {'value':[], 'standard_name':'measuring instrument reported measurement uncertainty',    'long_name':'measuring instrument reported measurement uncertainty',                 'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Measurement uncertainty (±), as given in metadata. In principal this refers to the inherent uncertainty on every measurement as a function of the quadratic addition of the accuracy and precision metrics (at the same confidence interval), but is often reported incosistently e.g. being solely determined from random errors (i.e. just the measuremental precision). It can be given in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%); or a percentage quantity after a fixed limit (i.e. 0.5%>=50).'},
    'measuring_instrument_documented_uncertainty':             {'value':[], 'standard_name':'measuring instrument documented measurement uncertainty',  'long_name':'measuring instrument documented measurement uncertainty',               'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Measurement uncertainty (±), as given in the instrumental manual/documentation. In principal this refers to the inherent uncertainty on every measurement as a function of the quadratic addition of the accuracy and precision metrics (at the same confidence interval), but is often reported incosistently e.g. being solely determined from random errors (i.e. just the measuremental precision). This can be given in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%); or a percentage quantity after a fixed limit (i.e. 0.5%>=50).'},
    'measuring_instrument_reported_accuracy':                  {'value':[], 'standard_name':'measuring instrument reported measurement accuracy',       'long_name':'measuring instrument reported measurement accuracy',                    'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Measurement accuracy (±), as given in metadata. Accuracy describes the difference between the measurement and the actual value of the part that is measured. It includes: Bias (a measure of the difference between the true value and the observed value of a part -- If the “true” value is unknown, it can be calculated by averaging several measurements with the most accurate measuring equipment available) and Linearity (a measure of how the size of the part affects the bias of a measurement system -- It is the difference in the observed bias values through the expected range of measurement). This can be given as in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%); or a percentage quantity after a fixed limit (i.e. 0.5%>=50).'}, 
    'measuring_instrument_documented_accuracy':                {'value':[], 'standard_name':'measuring instrument documented measurement accuracy',     'long_name':'measuring instrument documented measurement accuracy',                  'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Measurement accuracy (±), as given in the instrumental manual/documentation. Accuracy describes the difference between the measurement and the actual value of the part that is measured. It includes: Bias (a measure of the difference between the true value and the observed value of a part -- If the “true” value is unknown, it can be calculated by averaging several measurements with the most accurate measuring equipment available) and Linearity (a measure of how the size of the part affects the bias of a measurement system -- It is the difference in the observed bias values through the expected range of measurement). This can be given as in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%); or a percentage quantity after a fixed limit (i.e. 0.5%>=50).'}, 
    'measuring_instrument_reported_precision':                 {'value':[], 'standard_name':'measuring instrument reported measurement precision',      'long_name':'measuring instrument reported measurement precision',                   'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Measurement precision (±), as given in metadata. Precision describes the variation you see when you measure the same part repeatedly with the same device. It includes the following two types of variation: Repeatability (variation due to the measuring device -- it is the variation observed when the same operator measures the same part repeatedly with the same device) and Reproducibility (variation due to the operators and the interaction between operator and part -- It is the variation of the bias observed when different operators measure the same parts using the same device). This can be given as in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%); or a percentage quantity after a fixed limit (i.e. 0.5%>=50).'}, 
    'measuring_instrument_documented_precision':               {'value':[], 'standard_name':'measuring instrument documented measurement precision',    'long_name':'measuring instrument documented measurement precision',                 'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Measurement precision (±), as given in instrumental manual/documentation. Precision describes the variation you see when you measure the same part repeatedly with the same device. It includes the following two types of variation: Repeatability (variation due to the measuring device -- it is the variation observed when the same operator measures the same part repeatedly with the same device) and Reproducibility (variation due to the operators and the interaction between operator and part -- It is the variation of the bias observed when different operators measure the same parts using the same device). This can be given as in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%); or a percentage quantity after a fixed limit (i.e. 0.5%>=50).'}, 
    'measuring_instrument_reported_zero_drift':                {'value':[], 'standard_name':'measuring instrument reported zero drift',                 'long_name':'measuring instrument reported zero drift',                              'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Zero drift of measuring instrument per unit of time, as given in metadata. Zero drift (or baseline drift) refers to the shifting of the whole calibration by the same amount caused by slippage or due to undue warming up of the electronic circuits. It is reported as the maximum possible drift per unit of time in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%/day); or a percentage quantity after a fixed limit (i.e. 0.5%>=50/day).'}, 
    'measuring_instrument_documented_zero_drift':              {'value':[], 'standard_name':'measuring instrument documented zero drift',               'long_name':'measuring instrument documented zero drift',                            'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Zero drift of measuring instrument per unit of time, as given in instrumental manual/documentation. Zero drift (or baseline drift) refers to the shifting of the whole calibration by the same amount caused by slippage or due to undue warming up of the electronic circuits. It is reported as the maximum possible drift per unit of time in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%/day); or a percentage quantity after a fixed limit (i.e. 0.5%>=50/day).'}, 
    'measuring_instrument_reported_span_drift':                {'value':[], 'standard_name':'measuring instrument reported span drift',                 'long_name':'measuring instrument reported span drift',                              'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Span drift of measuring instrument per unit of time, as given in metadata. Span drift (or sensitivity drift) refers to when there is proportional change in the indication of an instrument all along the upward scale, hence higher calibrations end up being shifted more than lower calibrations. It is reported as the maximum possible drift per unit of time in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%/day); or a percentage quantity after a fixed limit (i.e. 0.5%>=50/day).'}, 
    'measuring_instrument_documented_span_drift':              {'value':[], 'standard_name':'measuring instrument documented span drift',               'long_name':'measuring instrument documented span drift',                            'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Span drift of measuring instrument per unit of time, as given in instrumental manual/documentation. Span drift (or sensitivity drift) refers to when there is proportional change in the indication of an instrument all along the upward scale, hence higher calibrations end up being shifted more than lower calibrations. It is reported as the maximum possible drift per unit of time in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%/day); or a percentage quantity after a fixed limit (i.e. 0.5%>=50/day).'}, 
    'measuring_instrument_reported_zonal_drift':               {'value':[], 'standard_name':'measuring instrument reported zonal drift',                'long_name':'measuring instrument reported zonal drift',                             'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Zonal drift of measuring instrument per unit of time, as given in metadata. Zonal drift refers to when drift occurs only over a portion of the full scale or span of an instrument, while the remaining portion of the scale remains unaffected. It is reported as the maximum possible drift per unit of time in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%/day); or a percentage quantity after a fixed limit (i.e. 0.5%>=50/day).'}, 
    'measuring_instrument_documented_zonal_drift':             {'value':[], 'standard_name':'measuring instrument documented zonal drift',              'long_name':'measuring instrument documented zonal drift',                           'units':parameter_details['standard_units'], 'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Zonal drift of measuring instrument per unit of time, as given in instrumental manual/documentation. Zonal drift refers to when drift occurs only over a portion of the full scale or span of an instrument, while the remaining portion of the scale remains unaffected. It is reported as the maximum possible drift per unit of time in absolute terms; as a percentage; the greater of either an absolute value or percentage (i.e. 25.0/0.5%/day); or a percentage quantity after a fixed limit (i.e. 0.5%>=50/day).'}, 
    'measuring_instrument_reported_measurement_resolution':    {'value':[], 'standard_name':'measuring instrument reported measurement resolution',     'long_name':'measuring instrument reported measurement resolution',                  'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN,        'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Measurement resolution, as given in metadata. The measurement resolution is defined as the smallest change or increment in the measured quantity that the instrument can detect. However it is often reported inconsistently, often being simply the number of digits an instrument can display, which does not relate to the actual physical resolution of the instrument.'}, 
    'measuring_instrument_documented_measurement_resolution':  {'value':[], 'standard_name':'measuring instrument documented measurement resolution',   'long_name':'measuring instrument documented measurement resolution',                'units':parameter_details['standard_units'], 'data_type':np.float32, 'string_format':np.NaN,        'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Measurement resolution, as given in instrumental manual/documentation. The measurement resolution is defined as the smallest change or increment in the measured quantity that the instrument can detect. However it is often reported inconsistently, often being simply the number of digits an instrument can display, which does not relate to the actual physical resolution of the instrument.'}, 
    'measuring_instrument_reported_absorption_cross_section':  {'value':[], 'standard_name':'measuring instrument reported absorption cross section',   'long_name':'measuring instrument reported absorption cross section',                'units':'cm2',                               'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Assumed molecule cross-section for parameter being measured (in cm2/molecule), as given in metadata. This field is only used for parameters being measured using optical methods, where a molecule cross section is assumed for processing the measurement values. Physically it is the effective area of the molecule that photon needs to traverse in order to be absorbed. The larger the absorption cross section, the easier it is to photoexcite the molecule. Can be a range: e.g. 1e-15-1.5e-15.'}, 
    'measuring_instrument_documented_absorption_cross_section':{'value':[], 'standard_name':'measuring instrument documented absorption cross section', 'long_name':'measuring instrument documented absorption cross section',              'units':'cm2',                               'data_type':np.object,  'string_format':'lower_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Assumed molecule cross-section for parameter being measured (in cm2/molecule), as given in instrumental manual/documentation. This field is only used for parameters being measured using optical methods, where a molecule cross section is assumed for processing the measurement values. Physically it is the effective area of the molecule that photon needs to traverse in order to be absorbed. The larger the absorption cross section, the easier it is to photoexcite the molecule. Can be a range: e.g. 1e-15-1.5e-15.'}, 
    'measuring_instrument_inlet_information':                  {'value':[], 'standard_name':'measuring instrument inlet information',                   'long_name':'measuring instrument measurement inlet information',                    'units':'unitless',                          'data_type':np.object,  'string_format':'lower_long',  'metadata_type':np.NaN,                            'description':'Description of sampling inlet of the measuring instrument.'}, 
    'measuring_instrument_calibration_scale':                  {'value':[], 'standard_name':'measuring instrument calibration scale',                   'long_name':'measuring instrument calibration scale name',                           'units':'unitless',                          'data_type':np.object,  'string_format':'upper_short', 'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'Name of calibration scale used for the calibration of the measuring instrument'},
    'network_provided_volume_standard_temperature':            {'value':[], 'standard_name':'network provided volume standard temperature',             'long_name':'network provided volume standard temperature',                          'units':'K',                                 'data_type':np.float64, 'string_format':np.NaN,        'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'The temperature (in Kelvin) associated with the volume of the sampled gas (which varies with temperature and pressure). This volume is typially normalised in-instrument to a standard temperature and pressure. These standard values typically follow network/continental/global standards (e.g. European Union) for the measured component. If no in-instrument normalisation is done then the reported temperature should be reported as the internal temperature of the instrument (i.e. the measurement conditions). If no numbers are reported explicitly per measurement, then the sample gas temperature is assumed to be the known network standard temperature for the measured component.'},  
    'network_provided_volume_standard_pressure':               {'value':[], 'standard_name':'network provided volume standard pressure',                'long_name':'network provided volume standard pressure ',                            'units':'hPa',                               'data_type':np.float64, 'string_format':np.NaN,        'metadata_type':'MEASUREMENT PROCESS INFORMATION', 'description':'The pressure (in hPa) associated with the volume of the sampled gas (which varies with temperature and pressure). This volume is typially normalised in-instrument to a standard temperature and pressure. These standard values typically follow network/continental/global standards (e.g. European Union) for the measured component. If no in-instrument normalisation is done then the reported pressure should be reported as the internal pressure of the instrument (i.e. the measurement conditions). If no numbers are reported explicitly per measurement, then the sample gas pressure is assumed to be the known network standard pressure for the measured component.'},  

    #CONTACT INFORMATION
    'principal_investigator_name':         {'value':[], 'standard_name':'principal investigator name',          'long_name':'principal investigator name',          'units':'unitless', 'data_type':np.object, 'string_format':'title_short', 'metadata_type':'STATION MISCELLANEOUS', 'description':'Full name of the principal scientific investigator for the specific reported data.'}, 
    'principal_investigator_institution':  {'value':[], 'standard_name':'principal investigator institution',   'long_name':'principal investigator institution',   'units':'unitless', 'data_type':np.object, 'string_format':'title_short', 'metadata_type':'STATION MISCELLANEOUS', 'description':'Institution of the principal scientific investigator for the specific reported data.'}, 
    'principal_investigator_email_address':{'value':[], 'standard_name':'principal investigator email address', 'long_name':'principal investigator email address', 'units':'unitless', 'data_type':np.object, 'string_format':'lower_short', 'metadata_type':'STATION MISCELLANEOUS', 'description':'Email address of the principal scientific investigator for the specific reported data.'}, 
    'contact_name':                        {'value':[], 'standard_name':'contact name',                         'long_name':'contact name',                         'units':'unitless', 'data_type':np.object, 'string_format':'title_short', 'metadata_type':'STATION MISCELLANEOUS', 'description':'Full name of the principal data contact for the specific reported data.'}, 
    'contact_institution':                 {'value':[], 'standard_name':'contact institution',                  'long_name':'contact institution',                  'units':'unitless', 'data_type':np.object, 'string_format':'title_short', 'metadata_type':'STATION MISCELLANEOUS', 'description':'Institution of the principal data contact for the specific reported data.'}, 
    'contact_email_address':               {'value':[], 'standard_name':'contact email address',                'long_name':'contact email address',                'units':'unitless', 'data_type':np.object, 'string_format':'lower_short', 'metadata_type':'STATION MISCELLANEOUS', 'description':'Email address of the principal data contact for the specific reported data.'}, 

    #TIMESTAMPS
    'meta_update_stamp':  {'value':[], 'standard_name':'metdata update timestamp', 'long_name':'metdata update timestamp', 'units':'unitless', 'data_type':np.uint32, 'string_format':np.NaN, 'metadata_type':np.NaN, 'description':'Time stamp of metadata updates in integer minutes from 0001-01-01 00:00 UTC.'}, 
    'data_download_stamp':{'value':[], 'standard_name':'data download timestamp',  'long_name':'data download timestamp',  'units':'unitless', 'data_type':np.uint32, 'string_format':np.NaN, 'metadata_type':np.NaN, 'description':'Time stamp of date/time of data download in integer minutes from 0001-01-01 00:00 UTC.'}, 

    #FURTHER DETAIL
    'network_sampling_details':     {'value':[], 'standard_name':'network sampling details',      'long_name':'network sampling details',      'units':'unitless', 'data_type':np.object, 'string_format':'lower_long', 'metadata_type':np.NaN, 'description':'Extra details provided by the reporting network about the sampling methods employed.'},
    'network_uncertainty_details':  {'value':[], 'standard_name':'network uncertainty details',   'long_name':'network uncertainty details',   'units':'unitless', 'data_type':np.object, 'string_format':'lower_long', 'metadata_type':np.NaN, 'description':'Extra details provided by the reporting network about the uncertainties involved with the measurement methods employed.'},
    'network_maintenance_details':  {'value':[], 'standard_name':'network maintenance details',   'long_name':'network maintenance details',   'units':'unitless', 'data_type':np.object, 'string_format':'lower_long', 'metadata_type':np.NaN, 'description':'Extra details provided by the reporting network about the operational maintenance done at the station.'},
    'network_qa_details':           {'value':[], 'standard_name':'network qa details',            'long_name':'network qa details',            'units':'unitless', 'data_type':np.object, 'string_format':'lower_long', 'metadata_type':np.NaN, 'description':'Extra details provided by the reporting network about the in-network quality assurance of measurements.'},
    'network_miscellaneous_details':{'value':[], 'standard_name':'network miscellaneous details', 'long_name':'network miscellaneous details', 'units':'unitless', 'data_type':np.object, 'string_format':'lower_long', 'metadata_type':np.NaN, 'description':'Extra miscellanous details provided by the reporting network.'},
    
    #WARNINGS
    'process_warnings':             {'value':[], 'standard_name':'process warnings',              'long_name':'process warnings',              'units':'unitless', 'data_type':np.object, 'string_format':'lower_long', 'metadata_type':np.NaN, 'description':'Warnings accumulated through GHOST processing regarding the data that should be considered.'}    
    }

    return standard_metadata

###--------------------------------------------------------------------------------------------------###
###DATA REPORTER PROVIDED DATA FLAG STANDARDISATIONS
###--------------------------------------------------------------------------------------------------###

#Define dictionary associating standardised data reporter provided data flag names with unique data flag codes

#These standardised data flags refer specifically to information/or lack of information provided by the data provider to fill the 'flag' field,
#and not to post-processing/quality control of observations (given separately by the ‘qa’ field).
    
standard_data_flag_name_to_data_flag_code = {

#Basic Data Flags
#-----------------------------------------------------
'Valid Data': 0,

'Preliminary Data': 1,

'Missing Data': 2, 

'Invalid Data - Unspecified': 3,

'Un-Flagged Data': 4,
#-----------------------------------------------------


#Estimated Data Flags 
#-----------------------------------------------------
'Estimated Data - Unspecified': 10,

'Estimated Data - Measured Negative Value': 11,

'Estimated Data - No Value Detected': 12,

'Estimated Data - Value Below Detection Limit': 13,

'Estimated Data - Value Above Detection Limit': 14,

'Estimated Data - Value Substituted from Secondary Monitor': 15,

'Estimated Data - Multiple Parameters Aggregated': 16,
#-----------------------------------------------------


#Extreme/Irregular Data Flags
#-----------------------------------------------------
'Extreme/Irregular Data - Unspecified': 20,

'Data Does Not Meet Internal Network Quality Control Criteria': 21,

'High Variability of Data': 22,

'Irregular Data Manually Screened and Accepted': 23,

'Irregular Data Manually Screened and Rejected': 24,

'Negative Value': 25,

'No Value Detected': 26,
    
'Reconstructed/Recalculated Data': 27,

'Value Close to Detection Limit': 28,

'Value Below Acceptable Range': 29,

'Value Above Acceptable Range': 30,

'Value Below Detection Limit': 31,

'Value Above Detection Limit': 32,

#-----------------------------------------------------


#Measurement Issue Data Flags
#-----------------------------------------------------
'Measurement Issue - Unspecified': 40,

'Chemical Issue': 41,

'Erroneous Sampling Operation': 42,

'Extreme Internal Instrument Meteorological Conditions': 43,

'Extreme Ambient Laboratory Meteorological Conditions': 44,

'Extreme External Meteorological Conditions': 45,

'Extreme Sample Transport Conditions': 46,

'Invalid Flow Rate': 47,

'Human Error': 48,

'Matrix Effect': 49,

'Mechanical Issue/Non-Operational Equipment': 50,

'No Technician': 51,

'Operational Maintenance Check Issue': 52,

'Physical Issue With Filter': 53,

'Power Failure': 54,

'Sample Diluted for Analysis': 55,

'Unmeasured Key Meteorological Parameter': 56,

#-----------------------------------------------------


#Operational Maintenance Data Flags
#-----------------------------------------------------
'Operational Maintenance - Unspecified': 60,

'Calibration': 61,

'Accuracy Check': 62,

'Blank Check': 63,

'Detection Limits Check': 64,

'Precision Check': 65,

'Retention Time Check': 66,

'Span Check': 67,

'Zero Check': 68,

'Instrumental Inspection': 69,

'Instrumental Repair': 70,

'Quality Control Audit': 71,
#-----------------------------------------------------


#Data Formatting/Processing Issue Data Flags
#-----------------------------------------------------
'Data Formatting/Processing Issue': 80,

'Corrected Data Formatting/Processing Issue': 81,
#-----------------------------------------------------


#Aggregation/Representation Flags
#-----------------------------------------------------
'Aggregation/Representation Issue - Unspecified': 90,

'Data Window Completeness < 90%': 91,

'Data Window Completeness < 75%': 92,

'Data Window Completeness < 66%': 93,

'Data Window Completeness < 50%': 94,

'Data Window Completeness < 25%': 95,

'>= 75% of Measurements in Window Below Detection Limit': 96,

'>= 50% of Measurements in Window Below Detection Limit': 97,
#-----------------------------------------------------


# Weather Condition Flags
#-----------------------------------------------------
'No Significant Weather': 100,

'Precipitation - Unspecified Intensity': 101,

'Precipitation - Light': 102,

'Precipitation - Moderate': 103,

'Precipitation - Heavy': 104,

'Drizzle - Unspecified Intensity': 105,

'Drizzle - Light': 106,

'Drizzle - Moderate': 107,

'Drizzle - Heavy': 108,

'Freezing Drizzle - Unspecified Intensity': 109,

'Freezing Drizzle - Light': 110,

'Freezing Drizzle - Moderate': 111,

'Freezing Drizzle - Heavy': 112,

'Rain - Unspecified Intensity': 113,

'Rain - Light': 114,

'Rain - Moderate': 115,

'Rain - Heavy': 116,

'Rain Shower/s - Unspecified Intensity': 117,

'Rain Shower/s - Light': 118,

'Rain Shower/s - Moderate': 119,

'Rain Shower/s - Heavy': 120,

'Freezing Rain - Unspecified Intensity': 121,

'Freezing Rain - Light': 122,

'Freezing Rain - Moderate': 123,

'Freezing Rain - Heavy': 124,

'Freezing Rain Shower/s - Unspecified Intensity': 125,

'Freezing Rain Shower/s - Light': 126,

'Freezing Rain Shower/s - Moderate': 127,

'Freezing Rain Shower/s - Heavy': 128,

'Snow - Unspecified Intensity': 129,

'Snow - Light': 130,

'Snow - Moderate': 131,

'Snow - Heavy': 132,

'Snow Shower/s - Unspecified Intensity': 133,

'Snow Shower/s - Light': 134,

'Snow Shower/s - Moderate': 135,

'Snow Shower/s - Heavy': 136,

'Hail - Unspecified Intensity': 137,

'Hail - Light': 138,

'Hail - Moderate': 139,

'Hail - Heavy': 140,

'Hail Shower/s - Unspecified Intensity': 141,

'Hail Shower/s - Light': 142,

'Hail Shower/s - Moderate': 143,

'Hail Shower/s - Heavy': 144,

'Ice Pellets - Unspecified Intensity': 145,

'Ice Pellets - Light': 146,

'Ice Pellets - Moderate': 147,

'Ice Pellets - Heavy': 148,

'Ice Pellets Shower/s - Unspecified Intensity': 149,

'Ice Pellets Shower/s - Light': 150,

'Ice Pellets Shower/s - Moderate': 151,

'Ice Pellets Shower/s - Heavy': 152,

'Snow Pellets - Unspecified Intensity': 153,

'Snow Pellets - Light': 154,

'Snow Pellets - Moderate': 155,

'Snow Pellets - Heavy': 156,

'Snow Pellets Shower/s - Unspecified Intensity': 157,

'Snow Pellets Shower/s - Light': 158,

'Snow Pellets Shower/s - Moderate': 159,

'Snow Pellets Shower/s - Heavy': 160,

'Snow Grains - Unspecified Intensity': 161,

'Snow Grains - Light': 162,

'Snow Grains - Moderate': 163,

'Snow Grains - Heavy': 164,

'Diamond Dust - Unspecified Intensity': 165,

'Diamond Dust - Light': 166,

'Diamond Dust - Moderate': 167,

'Diamond Dust - Heavy': 168,

'Glaze': 169,

'Rime': 170,

'Thunderstorm': 171,

'Funnel Cloud/s': 172,

'Squalls': 173,

'Tropical Cyclone (Cyclone/Hurricane/Typhoon)': 174,

'Duststorm': 175,

'Sandstorm': 176,

'Dust/Sand Whirls': 177,

'High Winds': 178,
#-----------------------------------------------------


# Atmospheric Obscuration/Local Area Contamination Flags
#-----------------------------------------------------

'No Atmospheric Obscuration': 180,

'Atmospheric Obscuration - Unknown': 181,

'Dust': 182,

'Blowing Dust': 183,

'Drifting Dust': 184,

'Sand': 185,

'Blowing Sand': 186,

'Drifting Sand': 187,

'Blowing Snow': 188,

'Drifting Snow': 189,

'Fog': 190,

'Freezing Fog': 191,

'Ground Fog': 192,

'Ice Fog': 193,

'Haze': 194,

'Mist': 195,

'Sea Spray': 196,

'Smoke': 197,

'Volcanic Ash': 198,

'No Local Contamination': 199,

'Local Contamination - Unspecified': 200,

'Agricultural Contamination': 201,

'Bird-Dropping Contamination': 202,

'Construction Contamination': 203,

'Industrial Contamination': 204,

'Insect Contamination': 205,

'Internal Laboratory/Instrument Contamination': 206,

'Pollen/Leaf Contamination': 207,

'Traffic Contamination': 208,
#-----------------------------------------------------


#Exceptional Event Data Flags
#-----------------------------------------------------
'Exceptional Event - Unspecified': 210,

#Natural Events
#---------------------------

'Seismic Activity': 211,

'Stratospheric Ozone Intrusion': 212,

'Volcanic Eruptions': 213, 

'Wildfire': 214,

#Anthropogenically Induced Events
#---------------------------
'Chemical Spill/Industrial Accident': 220,

'Cleanup After a Major Disaster': 221,

'Demolition': 222,

'Fireworks': 223,
    
'Infrequent Large Gathering': 224,

'Terrorist Act': 225,
#-----------------------------------------------------


# Meteorological Measurement Flags
#-----------------------------------------------------

'Visibility Distance Unlimited': 230,

'Ceiling Height Unlimited': 231
#-----------------------------------------------------

}

###--------------------------------------------------------------------------------------------------###
###QUALITY ASSURANCE FLAGS STANDARDISATIONS
###--------------------------------------------------------------------------------------------------###

#Define dictionary associating standardised quality assurance (QA) flag names with unique QA flag codes

standard_QA_name_to_QA_code = {

#Most Common QA Flags
#-----------------------------------------------------
#Missing Measurement
#Measurement is missing (i.e. NaN).
'Missing Measurement': 0,

#Infinite Value
#Value is infinite -- happens when data values are outside of the range that the float32 data type can handle (-3.4E+38 to +3.4E+38).
'Infinite Value': 1,

#Negative Measurement
#Measurement is negative in absolute terms.
'Negative Measurement': 2,

#Zero Measurement
#Have measurement equal to zero.
'Zero Measurement': 3,

#Not Maximum Data Quality Level
#Measurement is not of highest data quality level available from data provider
'Not Maximum Data Quality Level': 4,

#Preliminary Data
#Measurement has been flagged by data provider to be preliminary in nature (i.e. pending review)
'Preliminary Data': 5,

#Invalid Data Provider Flags - GHOST Decreed
#Measurements are associated with data quality flags given by the data provider which have been decreed by the GHOST project architects to suggest the measurements are associated with substantial uncertainty/bias
'Invalid Data Provider Flags - GHOST Decreed': 6,

#Invalid Data Provider Flags - Network Decreed
#Measurements are associated with data quality flags given by the data provider which have been decreed by the reporting network to suggest the measurements are associated with substantial uncertainty/bias 
'Invalid Data Provider Flags - Network Decreed': 7,

#No Valid Data to Average
#After screening by key QA flags, no valid data remains to average in temporal window
'No Valid Data to Average': 8,

#Duplicate Station
#Station has been decreed to be a duplicate (i.e. reporting the same data as from another network, but the data from another network has been preferred)   
'Duplicate Station': 9,
#-----------------------------------------------------


#Measurement Methodology Flags
#-----------------------------------------------------
#Methodology Not Mapped
#The measurement methodology used has not yet been mapped to standardised dictionaries of measurement methodologies.
'Methodology Not Mapped': 10,

#Assumed Primary Sampling
#A level of assumption has been made in determining the primary sampling type.
'Assumed Primary Sampling': 11,

#Assumed Sample Preparation
#A level of assumption has been made in determining the sample preparation.
'Assumed Sample Preparation': 12,

#Assumed Measurement Methodology
#A level of assumption has been made in determining the measurement methodology. 
'Assumed Measurement Methodology': 13,

#Unknown Primary Sampling Type
#The specific name of the primary sampling type is unknown.
'Unknown Primary Sampling Type': 14,

#Unknown Primary Sampling Instrument
#The specific name of the primary sampling instrument is unknown.
'Unknown Primary Sampling Instrument': 15,

#Unknown Sample Preparation Type
#The specific name of the sample preparation type is unknown.
'Unknown Sample Preparation Type': 16,

#Unknown Sample Preparation Technique
#The specific name of the sample preparation technique is unknown.
'Unknown Sample Preparation Technique': 17, 

#Unknown Measurement Method
#The specific name of the measurement method is unknown.
'Unknown Measurement Method': 18, 

#Unknown Measuring Instrument
#The specific name of measuring instrument is unknown. 
'Unknown Measuring Instrument': 19,

#Erroneous Primary Sampling
#The primary sampling is not appropriate to prepare the specific parameter for subsequent measurement.
'Erroneous Primary Sampling': 20,

#Erroneous Sample Preparation
#The sample preparation is not appropriate to prepare the specific parameter for subsequent measurement.
'Erroneous Sample Preparation': 21,

#Erroneous Measurement Method
#The measurement methodology used is not known to be able to measure the specific parameter. Only do check when known (or have assumed method).
'Erroneous Measurement Methodology': 22,

#Invalid QA Measurement Method 
#The specific measurement methodology has been decreed not to conform to QA standards as the method is not sufficiently proven/ subject to substantial biases/uncertainty. Only do check when known (or have assumed method).
'Invalid QA Measurement Methodology': 23,
#-----------------------------------------------------


#Sample Gas Volume Flags
#-----------------------------------------------------
#Sample Gas Volume - Network Standard
#The volume of the sampled gas (which varies with temperature and pressure) is assumed to be for a known standard network defined temperature and pressure for a component. This happens when no sample gas temperature/pressure are reported explicitly per measurement.
'Sample Gas Volume - Network Standard': 30,

#Sample Gas Volume - Unknown
#The volume of the sampled gas (which varies with temperature and pressure) is unknown. This happens when no sample gas temperature/pressure are reported explicitly per measurement AND there is no known standard network defined temperature and pressure for a component.
'Sample Gas Volume - Unknown': 31,

#Unit Conversion - Network Standard Sample Gas Volume Assumption
#Unit conversion has been done assuming the volume of the sampled gas (which varies with temperature and pressure) is for a known standard network defined temperature and pressure for a component.
'Unit Conversion - Network Standard Sample Gas Volume Assumption': 32,

#Unit Conversion - Educated Guess Sample Gas Volume Assumption
#Unit conversion has been done making an educated guess at the sampled gas temperature/pressure measurements are associated with.
'Unit Conversion - Educated Guess Sample Gas Volume Assumption': 33,
#-----------------------------------------------------


#Positional Metadata Doubt
#-----------------------------------------------------
#Station Position Doubt - DEM Decreed
#There exists significant doubt about the accuracy of the station position, as when cross-referencing the reported station altitude with the altitude from the ASTER v3 DEM (30 metre resolution) the absolute difference returned is >= 50 metres.
'Station Position Doubt - DEM Decreed': 40,

#Station Position Doubt - Manually Decreed
#There exists significant doubt about the accuracy of the station position, determined from empirical/word of mouth evidence
'Station Position Doubt - Manually Decreed': 41,
#-----------------------------------------------------


#Weather / Atmospheric Obscuration / Local Contamination / Exceptional Event Flags
#-----------------------------------------------------
#Local Precipitation
#There is some form of precipitation present at the measuring location at the time of measurement, reported by the data providers
'Local Precipitation': 50,

#Local Extreme Weather
#There is some form of extreme weather (heavy precipitation, funnel clouds, squalls, tropical cyclones, duststorm, sandstorm, dust/sand whirls, high winds) present at the measuring location at the time of measurement, reported by the data providers
'Local Extreme Weather': 51,

#Local Atmospheric Obscuration
#There is some form of atmospheric obscuration present at the measuring location at the time of measurement, reported by the data providers
'Local Atmospheric Obscuration': 52,

#Local Contamination
#There is some form of local contamination present at the measuring location at the time of measurement, reported by the data providers
'Local Contamination': 53,

#Local Exceptional Event
#There is some kind of exceptional event (both natural and anthropogenically induced) active at the measuring location at the time of measurement, reported by the data providers
'Local Exceptional Event': 54,


#Quality Control Checks
#-----------------------------------------------------
#Non-Integer Local Timezone (relative to UTC)
#Local timezone has been determined to be non-integer for time of measurement, relative to UTC.
'Non-Integer Local Timezone (relative to UTC)': 60,

#Below Documented Lower Limit of Detection
#Measurement is below or equal to the instrumental documented lower limit of detection.
'Below Documented Lower Limit of Detection': 61,

#Below Reported Lower Limit of Detection
#Measurement is below or equal to the network reported lower limit of detection. 
'Below Reported Lower Limit of Detection': 62,

#Below Preferential Lower Limit of Detection
#Measurement is below or equal to the preferential lower limit of detection (reported, and then documented). 
'Below Preferential Lower Limit of Detection': 63,

#Above Documented Upper Limit of Detection
#Measurement is above or equal to the instrumental documented upper limit of detection. 
'Above Documented Upper Limit of Detection': 64,

#Above Reported Upper Limit of Detection
#Measurement is above or equal to the network reported upper limit of detection. 
'Above Reported Upper Limit of Detection': 65,

#Above Preferential Upper Limit of Detection
#Measurement is above or equal to the preferential upper limit of detection (reported, and then documented). 
'Above Preferential Upper Limit of Detection': 66,

#Insufficient Measurement Resolution - Documented
#The documented resolution of measurement is coarser than a set limit (variable by measured parameter).
'Insufficient Measurement Resolution - Documented': 67,

#Insufficient Measurement Resolution - Reported
#The reported resolution of measurement is coarser than a set limit (variable by measured parameter).
'Insufficient Measurement Resolution - Reported': 68,

#Insufficient Measurement Resolution - Preferential
#The preferential resolution of measurement (reported, and then documented) is coarser than a set limit (variable by measured parameter).
'Insufficient Measurement Resolution - Preferential': 69,

#Insufficient Measurement Resolution - Empirical
#The resolution of measurement is analysed month by month. If the minimum difference between observations is coarser than a set limit (variable by measured parameter), measurements are flagged.
'Insufficient Measurement Resolution - Empirical': 70,

#Persistent Recurring Values - 5/6
#Check for persistently recurring values. Check is done by using a moving window of 9 measurements. If 5/6 (i.e. 83.33%) of values in the window are the same then the entire window is flagged.
'Persistent Recurring Values - 5/6': 71,

#Persistent Recurring Values - 9/12
#Check for persistently recurring values. Check is done by using a moving window of 12 measurements. If 9/12 (i.e. 75%) of values in the window are the same, then the entire window is flagged.
'Persistent Recurring Values - 9/12': 72,

#Persistent Recurring Values - 16/24
#Check for persistently recurring values. Check is done by using a moving window of 24 measurements. If 16/24 (i.e. 66.66%) of values in the window are the same, then the entire window is flagged.
'Persistent Recurring Values - 16/24': 73,

#Monthly Fractional Unique Values <= 10% 
#Fraction of unique values in a UTC month is <= 10% of total valid values in period 
'Monthly Fractional Unique Values <= 10%': 74,

#Monthly Fractional Unique Values <= 30% 
#Fraction of unique values in a UTC month is <= 30% of total valid values in period 
'Monthly Fractional Unique Values <= 30%': 75,

#Monthly Fractional Unique Values <= 50% 
#Fraction of unique values in a UTC month is <= 50% of total valid values in period 
'Monthly Fractional Unique Values <= 50%': 76,

#Monthly Fractional Unique Values <= 70% 
#Fraction of unique values in a UTC month is <= 70% of total valid values in period 
'Monthly Fractional Unique Values <= 70%': 77,

#Monthly Fractional Unique Values <= 90% 
#Fraction of unique values in a UTC month is <= 90% of total valid values in period 
'Monthly Fractional Unique Values <= 90%': 78,

#Data Outlier - Exceeds Scientifically Decreed Lower/Upper Limit
#The measured value is below or greater than scientifically feasible lower/upper limits (variable by parameter).
'Data Outlier - Exceeds Scientifically Decreed Lower/Upper Limit': 79,

#Data Outlier - Monthly Median Exceeds Scientifically Decreed Upper Limit 
#The median of the measurements in a month is greater than a scientifically feasible limit (variable by parameter).
'Data Outlier - Monthly Median Exceeds Scientifically Decreed Upper Limit': 80,

#Data Outlier - Network Decreed
#Data has been reported to be an outlier through data flags by the network data reporters (and not manually checked and verified as valid).
'Data Outlier - Network Decreed': 81,

#Data Outlier - Manually Decreed
#Data has been found and decreed manually to be an outlier.
'Data Outlier - Manually Decreed': 82,

#Possible Data Outlier - Monthly Adjusted Boxplot
#Measured value exceeds adjusted boxplot inner fence (lower or upper) of monthly data, therefore is a possible data outlier.
'Possible Data Outlier - Monthly Adjusted Boxplot': 83,

#Probable Data Outlier - Monthly Adjusted Boxplot
#Measured value exceeds adjusted boxplot outer fence (lower or upper) of monthly data, therefore is a probable data outlier.
'Probable Data Outlier - Monthly Adjusted Boxplot': 84,

#Monthly Distribution Consistency 
#QA flags which describe how consistent a monthly distribution is with other distributions for the same month, across the years.  
#This is done by calculating how consistent a monthly distribution of measurements is with the distributions for the same month of the neighbouring 5 years (backwards and forwards in time), as well as calculating how this metric of consistency differs from all other calculated metrics for the same month, across the time series.
#Kernel density estimation (kde) is used to estimate the probability density function (PDF) for the month in question, as well as the neighbouring year months.  
#The intersection is then taken between the month in question's PDF and the surrounding months PDF's, for which 1.0 - the median intersection value is taken as the metric of consistency (perfect intersection = 0.0, no intersection = 1.0). 
#Only the neighbouring 5 years are used, to limit the impact of any long-term trends biasing the comparison. 
#Additionally determined is to what extent does the median PDF intersection metric per month differ from the typical (median) intersection score for the same month across the entire timerseries (0.0 = perfectly consistent or better than the median intersection score, 1.0 = intersection score is as far below as possible the median intersection score). 
#This metric provides a weighting for the multi-annual variance of the monthly distributions. If there is a low intersection score, but there are typically always low intersection scores due to high variance across the years, then the difference score will be close to 0.0. Conversely, if there is a low intersection score, but there are typically always high intersection scores, then the difference score will be closer to 1.0, as it is outlying.           
#These metrics are then both summed to give a total classification score of the normality of the monthly distribution (2.0 = Extremely different from typical case. 0.0 = Extremely consistent with typical case).
#The classification score is split into 10 zones (in range increments of 0.2) from the most consistent monthly distributions in Zone 1 (score of 0.0 to 0.2), to the least consistent monthly distributions in Zone 10 (score of 1.8 to 2.0).
#All months for which a total classification score can be determined are flagged with a flag for the respective classification score zone.
#If 2/3, 4/6, 8/12 consecutive months are classed as zone 6 or higher, then it is suspected there is a systematic reason for the non-typical distributions, and the entire periods are flagged. 
#Monthly Distribution Consistency - Zone 1
'Monthly Distribution Consistency - Zone 1': 85,

#Monthly Distribution Consistency - Zone 2
'Monthly Distribution Consistency - Zone 2': 86,

#Monthly Distribution Consistency - Zone 3
'Monthly Distribution Consistency - Zone 3': 87,

#Monthly Distribution Consistency - Zone 4
'Monthly Distribution Consistency - Zone 4': 88,

#Monthly Distribution Consistency - Zone 5
'Monthly Distribution Consistency - Zone 5': 89,

#Monthly Distribution Consistency - Zone 6
'Monthly Distribution Consistency - Zone 6': 90,

#Monthly Distribution Consistency - Zone 7
'Monthly Distribution Consistency - Zone 7': 91,

#Monthly Distribution Consistency - Zone 8
'Monthly Distribution Consistency - Zone 8': 92,

#Monthly Distribution Consistency - Zone 9
'Monthly Distribution Consistency - Zone 9': 93,

#Monthly Distribution Consistency - Zone 10
'Monthly Distribution Consistency - Zone 10': 94,

#Monthly Distribution Consistency - Unclassified
'Monthly Distribution Consistency - Unclassified': 95,

#Systematic Inconsistent Monthly Distributions - 2/3 Months >= Zone 6
#2 out of 3 months' distributions are classed as Zone 6 or higher, suggesting there is potentially systematic reason for the inconsistent distributions across the 3 months.
'Systematic Inconsistent Monthly Distributions - 2/3 Months >= Zone 6': 96,

#Systematic Inconsistent Monthly Distributions - 4/6 Months >= Zone 6
#4 out of 6 months' distributions are classed as Zone 6 or higher, suggesting there is potentially systematic reason for the inconsistent distributions across the 6 months.
'Systematic Inconsistent Monthly Distributions - 4/6 Months >= Zone 6': 97,

#Systematic Inconsistent Monthly Distributions - 8/12 Months >= Zone 6
#8 out of 12 months' distributions are classed as Zone 6 or higher, suggesting there is potentially systematic reason for the inconsistent distributions across the 12 months.
'Systematic Inconsistent Monthly Distributions - 8/12 Months >= Zone 6': 98,

#-----------------------------------------------------


}

###--------------------------------------------------------------------------------------------------###
###DEFINE METADATA WARNING MODIFICATION FIELDS 
###--------------------------------------------------------------------------------------------------###

#define fields for which a warning will be wrote upon modification, associated with the names of the fields as they should be written in the text warning.

metadata_modification_field_to_warning_name = {'latitude':                                            'Latitude',
                                               'longitude':                                           'Longitude',
                                               'altitude':                                            'Altitude',
                                               'sampling_height':                                     'Sampling Height',
                                               'measurement_altitude':                                'Measurement Altitude',
                                               'standardised_network_provided_area_classification':   'Standardised Network Provided Area Classification',
                                               'standardised_network_provided_station_classification':'Standardised Network Provided Station Classification',
                                               'standardised_network_provided_main_emission_source':  'Standardised Network Provided Main Emission Source',
                                               'standardised_network_provided_land_use':              'Standardised Network Provided Land Use',
                                               'standardised_network_provided_terrain':               'Standardised Network Provided Terrain',
                                               'standardised_network_provided_measurement_scale':     'Standardised Network Provided Measurement Scale',
                                               'representative_radius':                               'Representative Radius',
                                               'primary_sampling_type':                               'Primary Sampling Type',
                                               'primary_sampling_instrument_name':                    'Primary Sampling Instrument Name',
                                               'sample_preparation_types':                            'Sample Preparation Types',
                                               'sample_preparation_techniques':                       'Sample Preparation Techniques',
                                               'measurement_methodology':                             'Measurement Methodology',
                                               'measuring_instrument_name':                           'Measuring Instrument Name'
                                              }

###--------------------------------------------------------------------------------------------------###
###MEASUREMENT PROCESS STANDARDISATIONS
###--------------------------------------------------------------------------------------------------###

#----------------------------------------------------------------------------------------#
#STANDARD SAMPLING INFORMATION
#----------------------------------------------------------------------------------------#
#define dictionary that contains information pertaining to standardised sampling types and specific sampling instrumentation

#The "sampling type" refers to type of sampling used by either the measurement instrument, or for 2 (or more) step-processes, in the sample preparation phase 

standard_sampling_types = {
#------------------------------------
#low volume continuous/high volume continuous

#Ambient air is continuously drawn in by a pump of some description. 

#The rate that air is sampled can be particularly important for the measurements of particulate matter.
#Many particulate matter measurements utilise a 2-step process, where standalone samplers first sample and filter air, and then the filtered products are measured by separate instrument.
#The standalone samplers can be separated into 2 discrete groups: low volume samplers & high volume samplers.

#Many instruments also have built-in continuous samplers. These are in-almost all cases low-volume, when not otherwise stated.

#--------
#low volume continuous

#Ambient air is continuously drawn in using a low volume sampler (typically sampling < ~24,000 litres of air over a 24-hour period)
#Some of these samplers have in-built filters, designed to specifically retain certain compounds.

'low volume continuous':{
    'instruments':{
        'APM PMS-104':            {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5'],                                                             'instrument_further_details':'See more detail here: http://apm.co.kr/sub2/view.php?id=20&ca_id=&page='},
        'Atmoservice PNS3D15'  :  {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],                                                       'instrument_further_details':'See more detail here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6096885/'},
        'Comde Derenda LVS 3.1':  {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],             'documented_flow_rate':'16.6667-58.3333', 'instrument_manual_name':'Comde_Derenda_MVS6.1_Specs'},
        'Comde Derenda MVS 6.1':  {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],             'documented_flow_rate':'16.6667-91.6667', 'instrument_manual_name':'Comde_Derenda_MVS6.1_Specs'},
        'Comde Derenda PNS 16T':  {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],             'documented_flow_rate':'16.6667-91.6667', 'instrument_manual_name':'Comde_Derenda_PNS_16T_Specs.pdf'},
        'Elecos APM-1':           {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],                                                       'instrument_further_details':'See more detail here: https://uk-air.defra.gov.uk/assets/documents/reports/empire/quarg/quarg_11.pdf', 'primary_sampling_process_details':'Make assumption sampling instrument is low volume continuous.', 'primary_sampling_assumption':True},
        'IVL P-Model':            {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],                                                       'instrument_further_details':'Not much detail provided, but believe this sampler uses different sampling inlets which are then named as P-Model S10 (for PM10) P-Model S2.5 (for PM2.5). Group different referenced "S" names all under "IVL P-Model". See detail here for S10 version: https://www.aces.su.se/reflab/wp-content/uploads/2016/11/ACES_Report_4.pdf; See detail here for S2.5 version: https://pubs.rsc.org/en/content/articlelanding/2017/em/c7em00122c#!divAbstract'},
        'KNF IP20-T':             {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],             'documented_flow_rate':'11.5',            'instrument_further_details':'See more detail here: https://www.overstocklabequipment.com/product-p/un811ktp.htm'},
        'Leckel SEQ 47-50':       {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],             'documented_flow_rate':'16.66-50',        'instrument_manual_name':'Leckel_SEQ_47-50_Specs.pdf'},
        'Leckel LVS3':            {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],             'documented_flow_rate':'16.66-38.33',     'instrument_manual_name':'Leckel_LVS3_MVS6_Specs.pdf'},
        'Leckel MVS6':            {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],             'documented_flow_rate':'38.33-58.33',     'instrument_manual_name':'Leckel_LVS3_MVS6_Specs.pdf'},
        'MCZ LVS1':               {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5'],                   'documented_flow_rate':'16.66-38.33',     'instrument_manual_name':'MCZ_LVS1_Specs.pdf'},
        'MCZ LVS16':              {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5'],                   'documented_flow_rate':'16.66-38.33',     'instrument_manual_name':'MCZ_LVS16_Specs.pdf'},
        'NILU EK':                {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','SO2'],                     'documented_flow_rate':'10.0',            'instrument_further_details':'See detail here: http://products.nilu.no/ProductsDivision/AirSamplers/SequentialAirSamplerTypeEK.aspx'},
        'NILU FK':                {'sample_preparation_type':'filter',           'isolated_parameters':['NO2','SO2','PM10','PM2.5','PM1'], 'documented_flow_rate':'0.5-5.0',         'instrument_further_details':'See more detail here:https://www.nilu.no/wp-content/uploads/dnn/38-2000-bs.pdf'},
        'NILU SS2000':            {'sample_preparation_type':'filter',           'isolated_parameters':['NO2','SO2'],                      'documented_flow_rate':'0.5-1.0',         'instrument_manual_name':'NILU_SS2000_Specs.pdf'},
        'Riemer Autocan 2/4':     {'sample_preparation_type':'filter',                                                                                                               'instrument_further_details':'See more details here: http://fp7.actris.eu/Portals/97/deliverables/PU/WP4_D4.7_M32.pdf'},
        'Riemer SPF4':            {'sample_preparation_type':'filter',           'isolated_parameters':['SO2','NH3','HNO3'],                                                         'instrument_further_details':'See more detail here: http://www.riemer-mt.com/low.volume.php'},
        'Rivm Universal Sampler': {'sample_preparation_type':'sorbent trapping',                                                           'documented_flow_rate':'0.1-0.65',        'instrument_further_details':'See more details here: https://www.rivm.nl/bibliotheek/rapporten/723101055.pdf'},
        'TCR Terora Skypost Gas': {'sample_preparation_type':'sorbent trapping',                                                           'documented_flow_rate':'0.5-10.0',        'instrument_manual_name':'TCR_Tecora_Skypost_Gas_Specs.pdf'},
        'TCR Tecora Sentinel':    {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],                                                       'instrument_manual_name':'TCR_Tecora_Sentinel_Specs.pdf', 'primary_sampling_assumption':True, 'primary_sampling_process_details':'Make assumption sampling instrument is low volume continuous.'},
        'Thermo Andersen FH95-KF':{'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],             'documented_flow_rate':'38.33',           'instrument_further_details':'Not many details but what there is, is here: https://core.ac.uk/download/pdf/38621632.pdf'},
        'Thermo Partisol 2000':   {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5'],                   'documented_flow_rate':'16.6667',         'instrument_manual_name':'Thermo_Partisol_2000_Manual.pdf'},
        'Thermo Partisol 2025':   {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5','PM1'],             'documented_flow_rate':'16.6667',         'instrument_manual_name':'Thermo_Partisol_2025_Specs.pdf'},
        'Thermo Partisol 2000i':  {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5'],                   'documented_flow_rate':'5-18',            'instrument_manual_name':'Thermo_Partisol_2000i_Manual.pdf'},
        'Thermo Partisol 2025i':  {'sample_preparation_type':'filter',           'isolated_parameters':['PM10','PM2.5'],                   'documented_flow_rate':'10-19',           'instrument_manual_name':'Thermo_Partisol_2025i_Specs.pdf.pdf'},
        'Thermo Partisol 2300':   {'sample_preparation_type':'filter',           'isolated_parameters':['SO2','PM10','PM2.5'],             'documented_flow_rate':'5-18',            'instrument_manual_name':'Thermo_Partisol_2300_Manual.pdf.pdf'}
    }
},

#--------
#high volume continuous

#Ambient air is continuously drawn in using a high volume sampler instrumentation (typically sampling > ~100,000 litres of air over a 24-hour period)
#Some of these samplers have in-built filters, designed to specifically retain certain compounds.

'high volume continuous':{
    'instruments':{
        'Digitel DH-77':             {'sample_preparation_type':'filter', 'isolated_parameters':['PM10','PM2.5','PM1'], 'documented_flow_rate':'100-1000',         'instrument_manual_name':'Digitel_DH-77_Specs.pdf'},
        'Digitel DHA-80':            {'sample_preparation_type':'filter', 'isolated_parameters':['PM10','PM2.5','PM1'], 'documented_flow_rate':'100-1000',         'instrument_manual_name':'Digitel_DHA-80_Manual.pdf'},
        'MCV CAV-A/mb':              {'sample_preparation_type':'filter', 'isolated_parameters':['PM10','PM2.5'],       'documented_flow_rate':'333.333-1333.333', 'instrument_further_details':'See more details here: http://www.mcvsa.com/Productes/Atmosfera/CaptadordealtovolumenMCVCAVAmb/tabid/114/Default.aspx'},
        'MCV CAV-A/msb':             {'sample_preparation_type':'filter', 'isolated_parameters':['PM10','PM2.5'],       'documented_flow_rate':'250-666.667',      'instrument_further_details':'See more details here: http://www.mcvsa.com/Productes/Atmosfera/CaptadordealtovolumensecuencialMCVCAVAMSb/tabid/139/language/es-ES/Default.aspx'},
        'Sierra 305-2000':           {'sample_preparation_type':'filter', 'isolated_parameters':['PM10','PM2.5','PM1'],                                            'instrument_further_details':'See more details here: https://pubs.acs.org/doi/pdf/10.1021/es00095a726'},
        'Thermo MFC-PM10':           {'sample_preparation_type':'filter', 'isolated_parameters':['PM10'],               'documented_flow_rate':'566,34-1699.01',   'instrument_manual_name':'Thermo MFC-PM10_Specs.pdf'},
        'Thermo VFC-PM10':           {'sample_preparation_type':'filter', 'isolated_parameters':['PM10'],               'documented_flow_rate':'1132.67',          'instrument_manual_name':'Thermo VFC-PM10_Specs.pdf'},
        'Tisch Environmental 6070D': {'sample_preparation_type':'filter', 'isolated_parameters':['PM10'],               'documented_flow_rate':'1020-1240',        'instrument_manual_name':'Tisch_Environmental_6070D_Specs.pdf'},     
        'Tisch Environmental 6070DX':{'sample_preparation_type':'filter', 'isolated_parameters':['PM10'],               'documented_flow_rate':'1020-1240',        'instrument_manual_name':'Tisch_Environmental_6070D_Specs.pdf'}     
        }
},

#------------------------------------
#injection

#The measuring instrument is injected with a limited quantity of air
#The limited injected quantity is typically pre-processed to aid the detection of a specific species.

'injection':{
    'instruments':{}
},

#------------------------------------
#continuous injection

#The measuring instrument is periodically injected with limited quantities of air.
#The injected quantities can either be from continuous automated collection, or from pre-processed loaded samples, which are periodically injected.

'continuous injection':{
    'instruments':{}
},

#------------------------------------
#passive              

#Air is not drawn in, rather the air sampled is the ambient air which interacts which measurement apparatus

'passive':{
    'instruments':{}
},

#------------------------------------
#remote            

#The measuring instrument does not actively sample air, but uses advanced optical techniques to measure parameters in air over long distances

'remote':{
    'instruments':{}
},

#------------------------------------
#manual              
#No instrument is used to determine measured values, they are determined manually 
#e.g. this applies for some colorimetric methods where measured values are derived manually via the colour of reagent after reaction with a compound of interest

'manual':{
    'instruments':{}
},

#------------------------------------
#unknown

'unknown':{
    'instruments':{}
}
}
    

#----------------------------------------------------------------------------------------#
#STANDARD SAMPLING PREPARATION INFORMATION
#----------------------------------------------------------------------------------------#

#define dictionary that contains information pertaining to standardised sampling preparation types
#Sample Preparation refers to all preparations made to the sample before the ultimate measurement 

standard_sampling_preparation_types = {
#------------------------------------
#flask

#Sample is collected/stored in measurement flasks/canisters from ambient air, or filled by pump
#The canisters can be opened for a short instantaneous window, or can be opened periodically over a long window to get a longer representative sample of air. 

'flask':{
    'specific_techniques':{}
},

#------------------------------------
#bag

#Sample is collected/stored in gas sampling bags (typically teflon) from ambient air, or filled by pump.
#These bags are a cheap alternative to canisters, with much reduced stability times.  

'bag':{
    'specific_techniques':{}
},


#------------------------------------
#preconcentration       

#This refers to the process of concentrating a sample before analysis, so that trace components won't be overlooked. 
#This is done typically through absorption of the sample onto a cooled, sorbent-packed trap before thermal desorption to transfer very quickly the sample to the analytical system. 
#Can be more than 1 trap, with different adsorbents.

'preconcentration':{
    'specific_techniques':{}
},

#------------------------------------
#filter                 

#Air is passed through a filtering system, selectively retaining compound(s) of interest

'filter':{
    'specific_techniques':{}
},

#------------------------------------
#filter pack

#Air is passed through a filter pack, selectively retaining compound(s) of interest
#The filter packs can contain different filters to target the retention of different compounds. 

'filter pack':{
    'specific_techniques':{
        '1-Stage Filter Pack':{},
        '2-Stage Filter Pack':{},
        '3-Stage Filter Pack':{}
    }
},

#------------------------------------
#denuder

#Air is passed through a denuder before analysis to selectively retain compound(s) of interest 
#a denuder is cylindrical or annular conduit or tube internally coated with a reagent that selectively reacts with a stable flow of gas

'denuder':{
    'specific_techniques':{
        'CEH DELTA':  {'isolated_parameters':['NH3','HNO3','SO2','HCl','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'],'preparation_further_details':'See specifications here: https://www.ceh.ac.uk/services/delta-active-sampler-system'},
        'Riemer DEN2':{                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    'preparation_further_details':'See here for confirmation of technique: https://ebas-submit.nilu.no/Submit-Data/Data-Reporting/Comments/Modules/all-instrument-types-manufacturers-models'},
        'UBA Olaf':   {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    'preparation_further_details':'See here for confirmation of technique: https://ebas-submit.nilu.no/Submit-Data/Data-Reporting/Comments/Category/Aerosol/Aerosol-chemical-speciation/inorganic-aerosol-chemistry/instrument-manufacturer'}
        }
},

#------------------------------------
#sorbent trapping           

#Sample is passed through a sorbent material to trap and retain the compound(s) of interest.
#Diffusive samplers are a specific example of the utilisation of sorbent trapping.

'sorbent trapping':{
    'specific_techniques':{

#--------
#diffusive sampler

#Diffussive samplers are a specific utilisation of sorbent trapping to passively trap species over long time periods, typically in urban areas.
#Diffusive samplers are an inexpensive method for measuring atmospheric composition over periods from one day to several weeks.
#These samplers operate on the principle of a sorbent trapping and retaining a compound of interest. 
#Sorption often takes place through reaction with a reagent (i.e. typically TEA for NO2).  
#After a suitable sample time, the trapped compound can then be analysed using colorimetry/spectrophotometry/ion chromatography/GC methods 

#Diffusive samplers come in various forms, e.g: Palmes, Passam, Ogawa Badge, and Radiello radial sampler tubes.
#The way these tubes differ is associated with the way they sample air and deal with variations in meteorological conditions. 

#This method is typically associated with very large uncertainties. 
#These uncertainties mainly come from:
#1. The non-specificity of trapped gases (i.e. HONO and PAN are readily uptaken when attempting to trap NO2).
#2. Variability in temperature, relative humidity, and air velocity directly impacting upon uptake rate of the desired gas, and effectiveness of the sorbent.
#Uncertainties associated with this methodology are typically reported between 10-35%: http://publications.jrc.ec.europa.eu/repository/bitstream/JRC51106/reqno_jrc51106_eur_23793.pdf[1].pdf

#Unless additional information is given, make assumption all diffusive samplers trap gases through sorbent trapping, followed by subsequent analysis through ion chromatography (by injection).

        'Diffusive Sampler':{'preparation_further_details':'See Reference for NO2: https://ec.europa.eu/jrc/en/publication/eur-scientific-and-technical-research-reports/review-application-diffusive-samplers-measurement-nitrogen-dioxide-ambient-air-european; reference for NH3:https://www.sciencedirect.com/science/article/pii/S1352231012001136; see reference for VOCs: http://www.siremlab.com/pages/wms/pdf/Engineering-Issue-Paper-2015-Passive-Samplers-for-Investigations-of-Air-Quality.pdf'}
    }
},


#------------------------------------
#reagent reaction

#Air is reacted with a liquid/solid chemical reagent to allow subsequent measurement of a specific compound

#The variety of types of passive sampling phases are numerous. These are described in detail below. 
#The name of the passive sampling method is given as the instrument name.
#However, yypically only passive sampling phase information is given, with no detail about the analysis phase. 
#Unless additional information is given, make assumption all passive sampling methods involve chemical reaction with a reagent, followed by subsequent analysis through spectrophotometry (by injection).

'reagent reaction':{
    'specific_techniques':{

#--------
#Griess-Saltzman

#NO2 reacts with water producing nitrous acid.
#This is reacted with an aromatic amine, to form a diazonium salt.
#Next, the solution is reacted with organic coupling agent (modified Saltzman reagent - 1960): 0.5% sulfanilic acid, 5% acetic acid, 0.005% NEDA(N-(1-Naphthyl)ethylenediamine dihydrochloride), to form a deeply coloured azo-dye.
#If an extended sampling time is required, there may be some bleaching by SO2. Saltzman recommended the addition of acetone to the reagent to prevent this.
#The conversion of NO2 to azo-dye is not quantitative, therefore a factor is introduced to represent the conversion efficiency under a given set of experimental conditions.
#The absorbance of the solution is then measured spectrophotometrically.
#Bennett (1930) gave stated a slighted altered reagent for the specific measurement of NO.

        'Griess-Saltzman':{'isolated_parameters':['NO','NO2'], 'preparation_further_details':'See References: Colorimetric Azo Dye Methods for the Atmospheric Analysis of Nitrogen Dioxide; Historical Development, Szonntagh, E.L, 1979; https://www3.epa.gov/ttn/naaqs/standards/nox/previous/1982-NOx-AQCD.pdf'},

#--------
#Lyshkow 

#A modification of the Griess-Saltzman method, using a modified reagent to produce an azo-dye: 0.15% Sulfanilamide, 1.5% tartaric acid, 0.005% NEDA, 0.005%  2-Naphthol-3, 6-disulphonic acid disodium salt. 
#The absorbance of the solution is measured spectrophotometrically, as with Griess-Saltzman.
#The reagent as with the Griess-Saltzman method can be slightly modified for the specific measurement of NO.

        'Lyshkow':{'isolated_parameters':['NO','NO2'], 'preparation_further_details':'See References: Colorimetric Azo Dye Methods for the Atmospheric Analysis of Nitrogen Dioxide; Historical Development, Szonntagh, E.L, 1979; https://www3.epa.gov/ttn/naaqs/standards/nox/previous/1982-NOx-AQCD.pdf'},

#--------
#Jacobs-Hochheiser 

#Designed by Jacobs and Hochheiser (1958) to avoid the azo-dye bleaching by SO2 that occurs in the Griess-Saltzman method.
#NO2 is absorbed in a solution of sodium hydroxide.
#Butanol is added as a surfactant to improve gas transfer.
#After sampling any sulfite from absorbed SO2 is oxidised with hydrogen peroxide.
#Next, the solution is acidified with phosphoric acid.
#A diazotizer and coupling agent are then added to produce an azo-dye.
#The absorbance of the solution is then measured spectrophotometrically.
    
        'Jacobs-Hochheiser':{'isolated_parameters':['NO2'], 'preparation_further_details':'See References: Colorimetric Azo Dye Methods for the Atmospheric Analysis of Nitrogen Dioxide; Historical Development, Szonntagh, E.L, 1979; https://www3.epa.gov/ttn/naaqs/standards/nox/previous/1982-NOx-AQCD.pdf'},

#--------
#Sodium Arsenite 

#NO2 is absorbed in an alkaline solution of sodium arsenite.
#Any sulfite in the solution is oxidised with peroxide.
#The solution is then the solution is acidified with phosphoric acid.
#An azo-dye is then formed through addition of sulfanilamide and NEDA.

        'Sodium Arsenite':{'isolated_parameters':['NO2'], 'preparation_further_details':'See details: https://www3.epa.gov/ttn/naaqs/standards/nox/previous/1982-NOx-AQCD.pdf'},

#--------
#TEA 

#NO2 is absorbed in a solution of triethanolamine and n-butanol surfactant.
#The solution is then reacted with the modified Saltzman reagent (1960) to form a deeply coloured azo-dye.
#The absorbance of the solution is then measured spectrophotometrically.

        'TEA':{'isolated_parameters':['NO2'], 'preparation_further_details':'See References: Comparison of Several Methods for Nitrogen Dioxide and Sulfur Dioxide in Metro Manila Air, Quirit, L.L, Hernandez, K.N, Lee, B.J, 2008; https://www3.epa.gov/ttn/naaqs/standards/nox/previous/1982-NOx-AQCD.pdf'},

#--------
#TGS-ANSA 

#Ambient air is bubbled through a solution containing triethanolamine, guaiacol and sodium metabilsulfite.
#Through this process, NO2 is converted to nitrite ion.
#8-Anilinonaphthalene-1-sulfonic acid ammonium salt (ANSA) and sulfanilamide are then added to the solution.
#The absorbance of the solution is then measured spectrophotometrically (at 550nm).

        'TGS-ANSA':{'isolated_parameters':['NO2'], 'preparation_further_details':'See details: https://www3.epa.gov/ttn/naaqs/standards/nox/previous/1982-NOx-AQCD.pdf'},

#--------
#Sodium Phenolate 

#A solution of alkaline phenol and hypochlorite react with ammonia to form indophenol blue.
#The solution is buffered at a pH of 9.5 with a borate buffer in order to decrease hydrolysis of cyanates and organic nitrogen compounds.
#The sample is then is then distilled into a solution of boric acid.
#The blue color formed is intensified with sodium nitroprusside and measured colorimetrically. 

        'Sodium Phenolate':{'isolated_parameters':['NH3'], 'preparation_further_details':'See further details:https://www.epa.gov/sites/production/files/2015-06/documents/epa-350.1.pdf'},

#--------
#Nessler 

#A dechlorinating agent is added to distilled water until the distillate shows no trace of ammonia. 
#NaOH is added to the distilled water until the pH is 9.5.     
#NH3 is absorbed from ambient air into the solution.
#The solution is buffered at a pH of 9.5 with a borate buffer in order to decrease hydrolysis of cyanates and organic nitrogen compounds.
#The sample is then is then distilled into a solution of boric acid.
#The distillate is then reacted with Nessler reagent.
#The absorbance of the solution is then measured spectrophotometrically (at 425nm).

        'Nessler':{'isolated_parameters':['NH3'], 'preparation_further_details':'See further details:https://www.umass.edu/mwwp/pdf/epa350_2NH3titration.pdf'},

#--------
#pararosaniline (also referred to as the West and Gaeke Method)

#SO2 is absorbed from air in a solution of potasium tetrachloromercuate. 
#A dichlorosulfitomercuate complex, which resists oxidation in air, is formed. This complex is stable to strong oxidants (i.e. ozone, NOy).
#This complex is then reacted with pararosaniline and formaldehyde to form intensely coloured C19H19N3 methyl sulphonic acid. 
#The absorbance of the solution is then measured spectrophotometrically.

        'Pararosaniline':{'isolated_parameters':['SO2'], 'preparation_further_details':'Also called West and Gaeke Method. See further details:https://nepis.epa.gov/Exe/ZyNET.exe/9101GLZ9.txt?ZyActionD=ZyDocument&Client=EPA&Index=Prior%20to%201976&Docs=&Query=&Time=&EndTime=&SearchMethod=1&TocRestrict=n&Toc=&TocEntry=&QField=&QFieldYear=&QFieldMonth=&QFieldDay=&UseQField=&IntQFieldOp=0&ExtQFieldOp=0&XmlQuery=&File=D%3A%5CZYFILES%5CINDEX%20DATA%5C70THRU75%5CTXT%5C00000022%5C9101GLZ9.txt&User=anonymous&Password=anonymous&SortMethod=h%7C-&MaximumDocuments=1&FuzzyDegree=0&ImageQuality=r75g8/r75g8/x150y150g16/i425&Display=hpfr&DefSeekPage=x&SearchBack=ZyActionL&Back=ZyActionS&BackDesc=Results%20page&MaximumPages=1&ZyEntry=12'},

#--------
#hydrogen peroxide 
#SO2 is absorbed from air in a solution of 3% hydrogen peroxide.
#This produces stable, non-volatile sulphuric acid. 
#SO3 is filtered in an 80% solution of isopropyl alcohol.
#Aliquots from the impinger are analysed by titration with barium perchlorate and 2% thorin indicator to the pink-orange end point.

        'Hydrogen Peroxide':{'isolated_parameters':['SO2'], 'preparation_further_details':'See further details:https://nepis.epa.gov/Exe/ZyNET.exe/9101GLZ9.txt?ZyActionD=ZyDocument&Client=EPA&Index=Prior%20to%201976&Docs=&Query=&Time=&EndTime=&SearchMethod=1&TocRestrict=n&Toc=&TocEntry=&QField=&QFieldYear=&QFieldMonth=&QFieldDay=&UseQField=&IntQFieldOp=0&ExtQFieldOp=0&XmlQuery=&File=D%3A%5CZYFILES%5CINDEX%20DATA%5C70THRU75%5CTXT%5C00000022%5C9101GLZ9.txt&User=ANONYMOUS&Password=anonymous&SortMethod=h%7C-&MaximumDocuments=1&FuzzyDegree=0&ImageQuality=r75g8/r75g8/x150y150g16/i425&Display=hpfr&DefSeekPage=x&SearchBack=ZyActionL&Back=ZyActionS&BackDesc=Results%20page&MaximumPages=1&ZyEntry=17 ; https://www.tandfonline.com/doi/pdf/10.1080/00022470.1972.10469715'},

#------
#detection tube
    
#Gas detection tubes operate on a chemical reaction between the vapour-phase compound and a liquid or solid detecting reagent, which is supported on an inert matrix. 
#The most common types of reactions are:
#   -Acid-base reactions These include reactions of acidic gases like HCl and HF with bases, and reaction of alkaline vapors such as ammonia with an acid in the tube. A dye present in the tube changes color as the pH changes on exposure to the vapors.
#   -Reduction-oxidation (Red-ox) reactions: These generate an oxidized or reduced compound, which has a different color. The chlorine tube uses oxidative coupling of colorless o-toluidine to form an orange azo-dye. White di-iodine pentoxide is reduced by CO and many organic vapors to form deep brown-colored iodine. Orange chromium (VI) is reduced by many organic compounds to form brown or green-colored Cr(III) compounds.
#   -Ligand-exchange reactions: These generate new complexes that are more colored than the starting reagents. The most notable is the conversion of white lead acetate to brown-black lead sulfide in the detection of H2S. In the case of phosphine, the exchange of PH3 for the chlorine ligand of HgCl2 releases HCl, which then causes a pH- dependent dye-color change.
#   -Pre-layers or Pre-tubes: These are used to condition the sample by controlling humidity, removing interferences, or transforming the analyte to another detectable compound. Examples include drying agents in NH3 and HCl tubes, organic removal by charcoal or oxidation in selective CO tubes, and oxidation of NO to NO2 in the nitrogen oxides tube.

#Most detection tubes are length-of-stain types. In these tubes, the reaction of the gas with the supported reagent is fast, compared to the transport of the bulk air sample through the tube. Therefore, all of the detected vapors are reacted within the tube. As a result, there is not a strong dependence of the readings on the rate at which the gas is sampled. However, a very high flow rate can cause some smearing to a high reading. 
#Conversely, low flow rates are less likely to affect the stain length, but can give low readings by concentrating the colored products in a shorter section of the tube. In cases of flow extremes, errors outside the standard 25% accuracy can be produced.
#Detection tubes are usually calibrated using piston hand pumps. The flow during a single pump stroke initially rises sharply and then decays exponentially. The best accuracy is therefore obtained when the flow through the tube mimics this profile.
#The advantages of detection tubes over other analytical methods are simplicity of use, rapid response, low cost, and very low maintenance. Because tubes are typically pre-calibrated, no calibration equipment is necessary. 
#Errors are prevented by directly marking the calibration information on each tube, and accuracy is further ensured by controlling the volume of air sampled. 
#Air sampling can also be performed using piston pumps, which latch into a precisely defined position to fix the volume. These pumps pull a strong vacuum initially and thus create substantially higher flowrate than the bellows pumps. Piston pumps generate a high flow initially followed by an approximately exponential decay, whereas bellows pumps provide a more steady flow initially followed by the slow decay. 
#The difference in flow patterns means that the pumps cannot be interchanged between types. For example, piston pumps sometime cause a smearing of the color stain when used on tubes originally developed for bellows pumps. This occurs because the higher flow rates do not allow enough contact time to give sharp endpoints when a piston pump is used.   

#Uncertainty on measurements is typically high with this method, typically reported between 10-25 %: https://www.honeywellanalytics.com/~/media/honeywell-analytics/products/colorimetric-gas-detection-tubes/documents/hagastubeshandbook667802enlr.pdf?la=en

#Unless additional information is given, make assumption all detection tubes trap and produce characteristic colour through reagent reaction, followed by subsequent analysis through colorimetry (manual).
        'Detection Tube':{'preparation_further_details':'See details: https://www.honeywellanalytics.com/~/media/honeywell-analytics/products/colorimetric-gas-detection-tubes/documents/hagastubeshandbook667802enlr.pdf?la=en'}
    }
},

#------------------------------------
#intermediate measurement

#A measurement is made by instrument/manually which feeds later into the ultimate measurement  

#---No specific sampling preparation name
'intermediate measurement':{
    'specific_techniques':{}
},

#------------------------------------
#unknown

'unknown':{
    'specific_techniques':{}
}
}


#----------------------------------------------------------------------------------------#
#STANDARD MEASUREMENT METHOD INFORMATION
#----------------------------------------------------------------------------------------#

#define dictionary that contains information pertaining to standardised measurement methodologies and specific instrumentation
#the "measurement methodology" refers to the name of the method used by the measuring instrument to determine the measured values for a parameter)

standard_measurement_methodologies = {
#------------------------------------
#ultraviolet photometry

#Operates on the principle that a specific species efficiently absorbs light at a known wavelength in the UV range. This is the case for ozone, at 253.65nm.
#The degree to which the UV light is absorbed by a specific species is directly related to the species concentration as described by the Beer-Lambert Law (I/Io = e−KLC; K = molecular absorption coefficient at STP (308 cm-1 atm-1 for O3), L = optical path length of cell, C = species concentration , I = light intensity of sample gas, Io = light intensity of sample without measured species (reference gas) ).

#Typical measurement operation for ozone:
#1. Ambient air is continuously drawn through the analyser by a vacuum pump. 
#2. Sample is split into two flow paths. One path incorporates a scrubber which selectively removes ozone, while the other path does not. 
#3. Both flow paths are connected to a solenoid valve, which switches at a fixed interval to allow either the sample gas stream (I), or the scrubbed, reference gas stream (Io), to flow through a quartz tube (cell) of accurately known length. 
#4. Typically a mercury vapour lamp at one end of the quartz cell produces a monochromatic beam of ultraviolet light at 254 nm. A vacuum diode at the opposite end of the cell measures the intensity of transmitted light (I/Io). 
#5. The intensities of 254 nm UV light transmitted through the sample (I) and ozone-free reference gas streams (Io) are related to the concentration (C) of ozone in the sample gas stream according to the Beer-Lambert Law.
#6. The analyzer’s microprocessor calculates the absolute concentration in molecules cm-3, and is then converted to a mixing ratio using in-instrument measured sample temperature and pressure.

#Known interferences: Gaseous hydrocarbons with strong absorption at 254 nm, such as aromatic hydrocarbons (i.e., benzene and substituted benzene rings).

#The greatest source of uncertaity regarding this measurement type is the cross sectional value used, for which the conventional value has been called into question: see further details: https://www.bipm.org/en/bipm/chemistry/gas-metrology/ozone/ozone-cross-section.html

'ultraviolet photometry':{'abbreviation':'UVP', 'sampling_type':'low volume continuous', 'measured_parameters':['O3'], 'qa_accepted_parameters':['O3'], 
    'instruments':{
        '2B 202':                   {'documented_lower_limit_of_detection':'3.0',  'documented_upper_limit_of_detection':'250000.0', 'documented_precision':'1.5/2.0%',  'documented_accuracy':'1.5/2.0%', 'documented_zero_drift':'2.0/day',    'documented_span_drift':'1.0%/day',       'documented_resolution':'0.1',   'documented_absorption_cross_section':'1.15e-17',   'documented_flow_rate':'1.0',     'instrument_manual_name':'2B_202_Manual.pdf'},
        '2B 205':                   {'documented_lower_limit_of_detection':'2.0',  'documented_upper_limit_of_detection':'250000.0', 'documented_precision':'1.0/2.0%',  'documented_accuracy':'1.0/2.0%', 'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/day',       'documented_resolution':'0.1',   'documented_absorption_cross_section':'1.15e-17',   'documented_flow_rate':'1.8',     'instrument_manual_name':'2B_205_Manual.pdf'},
        '2B 211':                   {'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':'2000.0',   'documented_precision':'0.5/1.0%',  'documented_accuracy':'1.0/2.0%', 'documented_zero_drift':'1.0/day',    'documented_span_drift':'0.5%/day',       'documented_resolution':'0.1',   'documented_absorption_cross_section':'1.15e-17',   'documented_flow_rate':'2.0',     'instrument_manual_name':'2B_211_Manual.pdf'},
        '2B 106-L':                 {'documented_lower_limit_of_detection':'3.0',  'documented_upper_limit_of_detection':'100000.0', 'documented_precision':'1.5/2.0%',  'documented_accuracy':'1.5/2.0%', 'documented_zero_drift':'3.0/day',    'documented_span_drift':'1.0%/day',       'documented_resolution':'0.1',   'documented_absorption_cross_section':'1.15e-17',   'documented_flow_rate':'1.0',     'instrument_manual_name':'2B_106-L_Manual.pdf'},
        '2B 106-OEM-L':             {'documented_lower_limit_of_detection':'3.0',  'documented_upper_limit_of_detection':'100000.0', 'documented_precision':'1.5/2.0%',  'documented_accuracy':'1.5/2.0%', 'documented_zero_drift':'3.0/day',    'documented_span_drift':'1.0%/day',       'documented_resolution':'0.1',   'documented_absorption_cross_section':'1.15e-17',   'documented_flow_rate':'1.0',     'instrument_manual_name':'2B_106-L_Manual.pdf'},
        '2B POM':                   {'documented_lower_limit_of_detection':'3.0',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'1.5/2.0%',  'documented_accuracy':'1.5/2.0%', 'documented_zero_drift':'2.0/day',    'documented_span_drift':'1.0%/day',       'documented_resolution':'0.1',   'documented_absorption_cross_section':'1.15e-17',   'documented_flow_rate':'0.8',     'instrument_manual_name':'2B_POM_Manual.pdf'},
        'CSI 3100':                 {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:https://www.healtheffects.org/system/files/ResearchReport-70.pdf'},
        'Dasibi 1003-AH':           {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Dasibi 1003-PC':           {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Dasibi 1003-RS':           {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Dasibi 1008-AH':           {'documented_lower_limit_of_detection':'2.0',  'documented_upper_limit_of_detection':'1000.0',   'documented_precision':'1.0',       'documented_accuracy':'1.0',                                                                                                                       'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'2.0',     'instrument_manual_name':'Dasibi_1008_Manual.pdf'},
        'Dasibi 1008-HC':           {'documented_lower_limit_of_detection':'2.0',  'documented_upper_limit_of_detection':'1000.0',   'documented_precision':'1.0',       'documented_accuracy':'1.0',                                                                                                                       'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'2.0',     'instrument_manual_name':'Dasibi_1008_Manual.pdf'},
        'Dasibi 1008-PC':           {'documented_lower_limit_of_detection':'2.0',  'documented_upper_limit_of_detection':'1000.0',   'documented_precision':'1.0',       'documented_accuracy':'1.0',                                                                                                                       'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'2.0',     'instrument_manual_name':'Dasibi_1008_Manual.pdf'},
        'Dasibi 1008-PS':           {'documented_lower_limit_of_detection':'2.0',  'documented_upper_limit_of_detection':'1000.0',   'documented_precision':'1.0',       'documented_accuracy':'1.0',                                                                                                                       'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'2.0',     'instrument_manual_name':'Dasibi_1008_Manual.pdf'},
        'Dasibi 1008-RS':           {'documented_lower_limit_of_detection':'2.0',  'documented_upper_limit_of_detection':'1000.0',   'documented_precision':'1.0',       'documented_accuracy':'1.0',                                                                                                                       'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'2.0',     'instrument_manual_name':'Dasibi_1008_Manual.pdf'},
        'Dasibi 1108 W/GEN':        {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:https://books.google.es/books?id=SuIsDwAAQBAJ&pg=PA48&lpg=PA48&dq=dasibi+1108&source=bl&ots=tnXwozqtBT&sig=ACfU3U0bA3efWduqeNFD2g2zk-RqP7cQFg&hl=en&sa=X&ved=2ahUKEwjRk8bX2f7fAhUD4OAKHSr1A5wQ6AEwDHoECAAQAQ#v=onepage&q=dasibi%201108&f=false'},
        'DKK-TOA GUX-113E':         {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:https://books.google.es/books?id=KpQhwCLgVgAC&pg=PA209&lpg=PA209&dq=GUX-113E&source=bl&ots=e1rn8FsfBV&sig=Sd0wYb5DAjG4GZBrpsFPAuZ6Uqo&hl=en&sa=X&ved=0ahUKEwjbivvTq67bAhVDuRQKHY50ASYQ6AEIMjAB#v=onepage&q=GUX-113E&f=false'},
        'DKK-TOA GUX-313E':         {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'1000.0',                                                                         'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/day',                                                                                                                              'instrument_manual_name':'DKK-TOA_GUX-313E_Specs.pdf'} ,
        'Dylec 1100':               {'documented_lower_limit_of_detection':'0.1',  'documented_upper_limit_of_detection':'500.0',    'documented_precision':'0.5%',      'documented_accuracy':'0.5%',     'documented_zero_drift':'0.5%/month', 'documented_span_drift':'0.5%/month',                                                                                          'documented_flow_rate':'1.5-2.0', 'instrument_further_details':'Specifications online:https://translate.google.com/translate?hl=en&sl=ja&u=http://www.dylec.co.jp/01ozone/m1100.html&prev=search'},
        'Dylec 1200':               {                                              'documented_upper_limit_of_detection':'20000.0',  'documented_precision':'0.5%',      'documented_accuracy':'0.5%',     'documented_zero_drift':'0.5%/month', 'documented_span_drift':'0.5%/month',                                                                                          'documented_flow_rate':'1.5-2.0', 'instrument_further_details':'Specifications online:https://translate.google.com/translate?hl=en&sl=ja&u=http://www.dylec.co.jp/01ozone/m1200.html&prev=search'},
        'Ebara EG-2001F':           {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here: http://www.data.jma.go.jp/gmd/env/ghg_obs/en/station/station_minamitorishima.html'},
        'Ebara EG-2001FTP':         {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here: http://www.data.jma.go.jp/gmd/env/ghg_obs/en/station/station_minamitorishima.html'},
        'Ebara EG-3000F':           {                                              'documented_upper_limit_of_detection':'1000.0',   'documented_precision':'5.0',       'documented_accuracy':'5.0',      'documented_zero_drift':'5.0/month',  'documented_span_drift':'5.0/month',                                                                                           'documented_flow_rate':'1.5',     'instrument_manual_name':'Ebara_EG-3000F_Manual.pdf'},
        'Ecotech EC9810':           {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'20000.0',  'documented_precision':'1.0/1.0%',                                    'documented_zero_drift':'1.0/day',    'documented_span_drift':'0.5%/day',       'documented_resolution':'0.001', 'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'0.5',     'instrument_manual_name':'Ecotech_EC9810_Manual.pdf'},
        'Ecotech Serinus10':        {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'20000.0',  'documented_precision':'0.5/0.2%',  'documented_accuracy':'1.0%',     'documented_zero_drift':'0.3/day',    'documented_span_drift':'0.3/0.5%/7days',                                                                                      'documented_flow_rate':'0.5',     'instrument_manual_name':'Ecotech_Serinus10_Specs.pdf'},
        'Environics 300':           {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Environnement SA O331M':   {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:https://www.geo.fu-berlin.de/met/ag/trumf/Lehre/Lehrveranstaltungen_alt/Dok_2006/Vorlesung_SchwarzesDreieck_20060503.pdf'},
        'Environnement SA O341M':   {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Environnement SA O342e':   {'documented_lower_limit_of_detection':'0.2',  'documented_upper_limit_of_detection':'10000.0',                                      'documented_accuracy':'1.0%',     'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5%/day',                                                                                            'documented_flow_rate':'1.0',     'instrument_manual_name':'EnvironnementSA_O342e_Specs.pdf'},
        'Environnement SA O342M':   {'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'10000.0',                                      'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0/week',   'documented_span_drift':'1.0%/week',                                                                                           'documented_flow_rate':'1.0',     'instrument_manual_name':'EnvironnementSA_O342M_Specs.pdf'},
        'Horiba APOA300':           {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:https://books.google.es/books?id=qmWXBwAAQBAJ&pg=PA88&lpg=PA88&dq=horiba+apna+300&source=bl&ots=dBQJ0TEai1&sig=plLm6SAUIZwuIhy37JzI8nZXccQ&hl=en&sa=X&ved=2ahUKEwihxOS9uvDfAhWGnhQKHYzoAPcQ6AEwCHoECAIQAQ#v=onepage&q=horiba%20apna%20300&f=false'},
        'Horiba APOA350':           {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:https://books.google.es/books?id=WqvzBQAAQBAJ&pg=SL13-PA35&lpg=SL13-PA35&dq=Horiba+APOA350&source=bl&ots=Ltdp9qBqrc&sig=qKEZbhp9Vd0UPEbMIL6Al8CKITM&hl=en&sa=X&ved=0ahUKEwjcjJKai9bbAhXH7BQKHa_gDvQ4ChDoAQgtMAE#v=onepage&q=Horiba%20APOA350&f=false'},
        'Horiba APOA360':           {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'1.0%',      'documented_accuracy':'1.0%',     'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5/day',                                                                                             'documented_flow_rate':'0.8',     'instrument_manual_name':'Horiba_APMA360_Specs.pdf'},
        'Horiba APOA370':           {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'100000.0', 'documented_precision':'1.0%',      'documented_accuracy':'1.0%',     'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5/day',                                                                                             'documented_flow_rate':'0.7',     'instrument_manual_name':'Horiba_APOA370_Specs.pdf'},
        'MCV 48AUV':                {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://www.academia.edu/18561546/Temporal_patterns_of_surface_ozone_levels_in_different_habitats_of_the_North_Western_Mediterranean_basin'},
        'PCI LC-12':                {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Philips K50094':           {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/philipsK50094/view?facet=HTML+Representation'},   
        'Sabio 6030':               {'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'0.5%',      'documented_accuracy':'0.5%',     'documented_zero_drift':'1.0/month',  'documented_span_drift':'1.0%/month',                                                                                          'documented_flow_rate':'0.5-1.0', 'instrument_manual_name':'Sabio_6030_Specs.pdf'},
        'Seres OZ 2000G':           {'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':'10000.0',                                      'documented_accuracy':'1.0%',     'documented_zero_drift':'1.5/week',   'documented_span_drift':'1.0%/month',                                                                                          'documented_flow_rate':'1.0',     'instrument_manual_name':'Seres_OZ2000G_Specs.pdf'},
        'SIR S-5014':               {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here: http://www.ecomission.co.za/products/model_s_5014.htm'},
        'Tanabyte 722':             {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'1.0',       'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0%/week',  'documented_span_drift':'1.0%/week',      'documented_resolution':'0.1',                                                       'documented_flow_rate':'0.5-2.0', 'instrument_manual_name':'Tanabyte_72x_Specs.pdf'},
        'Tanabyte 723':             {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'1.0',       'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0%/week',  'documented_span_drift':'1.0%/week',      'documented_resolution':'0.1',                                                       'documented_flow_rate':'0.5-2.0', 'instrument_manual_name':'Tanabyte_72x_Specs.pdf'},
        'Tanabyte 724':             {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'1.0',       'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0%/week',  'documented_span_drift':'1.0%/week',      'documented_resolution':'0.1',                                                       'documented_flow_rate':'0.5-2.0', 'instrument_manual_name':'Tanabyte_72x_Specs.pdf'},
        'Tanabyte 725':             {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'1.0',       'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0%/week',  'documented_span_drift':'1.0%/week',      'documented_resolution':'0.1',                                                       'documented_flow_rate':'0.5-2.0', 'instrument_manual_name':'Tanabyte_72x_Specs.pdf'},
        'Tanabyte 726':             {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'1.0',       'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0%/week',  'documented_span_drift':'1.0%/week',      'documented_resolution':'0.1',                                                       'documented_flow_rate':'0.5-2.0', 'instrument_manual_name':'Tanabyte_72x_Specs.pdf'},
        'Tanabyte 734':             {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'1.0',       'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0%/week',  'documented_span_drift':'1.0%/week',      'documented_resolution':'0.1',                                                       'documented_flow_rate':'0.5-2.0', 'instrument_manual_name':'Tanabyte_72x_Specs.pdf'},
        'Teledyne 400':             {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},   
        'Teledyne 400A':            {'documented_lower_limit_of_detection':'0.6',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'0.5%',      'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/day',                                        'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'0.8',     'instrument_manual_name':'Teledyne_400A_Manual.pdf'},   
        'Teledyne 400E':            {'documented_lower_limit_of_detection':'0.6',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'0.5%',      'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/day',                                        'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'0.8',     'instrument_manual_name':'Teledyne_400E_Manual.pdf'},   
        'Teledyne ML8810':          {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Teledyne ML9810':          {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Teledyne ML9810A':         {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/ML9810A/view'},
        'Teledyne ML9810B':         {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Teledyne ML9811':          {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Teledyne ML9812':          {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Teledyne T204-O3':         {'documented_lower_limit_of_detection':'2.0',  'documented_upper_limit_of_detection':'1000.0',   'documented_precision':'0.5%>=100', 'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0%/day',   'documented_span_drift':'1.0%/day',                                        'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'0.8',     'instrument_manual_name':'Teledyne_T204_Manual.pdf'},
        'Teledyne T400':            {'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'10000.0',  'documented_precision':'0.5%>=100', 'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/day',                                                                                            'documented_flow_rate':'0.8',     'instrument_manual_name':'Teledyne_T400_Manual.pdf'},
        'Thermo 49':                {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Thermo 49C':               {'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':'200000.0', 'documented_precision':'1.0',       'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/month',                                      'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'1.0-3.0', 'instrument_manual_name':'Thermo_49C_Manual.pdf'},   
        'Thermo 49CPS':             {'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':'5000.0',                                       'documented_accuracy':'1.0%',     'documented_zero_drift':'100.0/day',  'documented_span_drift':'1.0%/day',                                                                                            'documented_flow_rate':'1.0-3.0', 'instrument_manual_name':'Thermo_49CPS_Specs.pdf'},
        'Thermo 49i':               {'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'200000.0', 'documented_precision':'1.0',       'documented_accuracy':'1.0%',     'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/month',                                      'documented_absorption_cross_section':'1.1476e-17', 'documented_flow_rate':'1.0-3.0', 'instrument_manual_name':'Thermo_49i_Manual.pdf,Thermo_49i_Specs.pdf'},
        'Thermo 49w':               {                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:https://www.rivm.nl/en/Documents_and_publications/Scientific/Reports/1993/juni/Acceptance_procedure_and_testresults_of_the_ozone_analyzers_model_49W_made_by_Thermo_Environmental_Instruments_on_behalf_of_the_Dutch_National_Air_Quality_Monitoring_Network'},
        'Wedding & Associates 1010':{                                                                                                                                                                                                                                                                                                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet photometry here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'}
    }
},

#------------------------------------
#visible photometry

#Operates on the principle that a specific species efficiently absorbs light at a known wavelength in the visible range. This is the case for NO2, at 405nm
#The degree to which the visible light is absorbed by a specific species is directly related to the species concentration as described by the Beer-Lambert Law (C = 1/Lσ * ln(Io/I) ; σ = absorption cross section (6.06×10-19 cm2 molec-1 for NO2 at 405nm), L = optical path length of cell, C = species concentration , I = light intensity of sample gas, Io = light intensity of sample without measured species (reference gas) ).

#Typical measurement operation for NO2 (and indirectly NO):
#1. Ambient air is continuously drawn through the instrument by a vacuum pump. 
#2. An NO2 scrubber valve alternately bypasses and sends the sample air through a heated NO2 scrubber to remove all NO2 in the sample. 
#3. The NO2-scrubbed/unscrubbed air alternately passes through the optical cell (of known length) and the cell flow meter. Alternate switching of the NO2 scrubber valve once every 5 seconds allows the measurement of a light intensity in the absence of NO2 and presence of NO2. 
#4. A light-emitting diode (LED) emits a monochromatic beam of light at 405nm. A photodiode at the other end of the cell measures the intensity of transmitted light (I/Io).
#5. The intensities of 254 nm light transmitted through the sample (I) and NO2-free reference gas streams (Io) are related to the concentration (C) of ozone in the sample gas stream according to the Beer-Lambert Law.
#6. The analyser's microprocessor calculates the absolute concentration in molecules cm-3, and is then converted to a mixing ratio using in-instrument measured sample temperature and pressure.
#7. NO is measured indirectly by bypassing the NO2 scrubber and measuring the light intensity while adding (I) or not adding (Io) ozone to convert NO to NO2 according to the well-known reaction: NO + O3 → NO2 + O2

#Known interferences: Water vapour, small particles (< 5 um)

'visible photometry':{'abbreviation':'VP', 'sampling_type':'low volume continuous', 'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO','NO2','NOx'], 
    'instruments':{
        '2B 405nm':{'documented_lower_limit_of_detection':'1.0', 'documented_upper_limit_of_detection':{'NO':'2000.0','NO2':'10000.0'}, 'documented_precision':'0.5/0.5%', 'documented_resolution':'0.1', 'documented_absorption_cross_section':'6.06e-19', 'instrument_manual_name':'2B_405nm_Manual.pdf'}
    }
},
#--------------------------------------------------------------------------------------###
#ethylene chemiluminescence 

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of ozone with ethylene, leading to an excited molecule.
#The return to a fundamental electronic state of the excited molecules is made by luminous radiation in a specific spectrum, which can be measured.
#The concentration of sample ozone is directly proportional to the intensity of light emitted.
#The broadband emission is detected using a photomultiplier tube (at 440 nm for ethylene + ozone).

#Known interferences: water vapour

'ethylene chemiluminescence':{'abbreviation':'ECL', 'sampling_type':'low volume continuous', 'measured_parameters':['O3'], 'qa_accepted_parameters':['O3'], 
    'instruments':{
        'Beckman 950A':    {                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ethylene chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Bendix 8002':     {'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':'200.0', 'instrument_further_details':'No manual or specifications available but confirmation it is uses ethylene chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf ; details about limit of detection:https://books.google.es/books?id=G0EYUql3_98C&pg=PA283&lpg=PA283&dq=thermo+43a&source=bl&ots=k2RUJXlAVJ&sig=yIxG08AoSyuqedtS6dp-O7O04X4&hl=en&sa=X&ved=0ahUKEwj17qqB3OfbAhVCVRQKHanFDlo4ChDoAQgoMAA#v=onepage&q=thermo%2043a&f=false'},
        'CSI 2000':        {                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ethylene chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'McMillan 1100-1': {                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ethylene chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'McMillan 1100-2': {                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ethylene chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'McMillan 1100-3': {                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ethylene chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Meloy 0A325-2R':  {                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ethylene chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Meloy 0A350-2R':  {                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ethylene chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'},
        'Teledyne ML8410E':{                                                                                             'instrument_further_details':'No manual or specifications available but confirmation it is uses ethylene chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'}
    }
},

#--------------------------------------------------------------------------------------###
#eosin-Y chemiluminescence 

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of ozone with eosin-Y, leading to an excited molecule.
#The return to a fundamental electronic state of the excited molecules is made by luminous radiation in a specific spectrum, which can be measured.
#The concentration of sample ozone is directly proportional to the intensity of light emitted.
#The broadband emission is detected using a photomultiplier tube (at ~560 nm for eosin-Y + ozone).

#Known interferences: water vapour

'eosin-Y chemiluminescence':{'abbreviation':'EYCL', 'sampling_type':'low volume continuous', 'measured_parameters':['O3'], 'qa_accepted_parameters':['O3'], 
    'instruments':{
        'Luminox LOZ-3':{'instrument_further_details':'No manual or specifications available but confirmation it is uses eosin-Y chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'}
    }
},

#--------------------------------------------------------------------------------------###
#rhodamine B chemiluminescence 

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of ozone with rhodamine B, leading to an excited molecule.
#The return to a fundamental electronic state of the excited molecules is made by luminous radiation in a specific spectrum, which can be measured.
#The concentration of sample ozone is directly proportional to the intensity of light emitted.
#The broadband emission is detected using a photomultiplier tube (at ~580 nm for rhodamine B + ozone).

#Known interferences: water vapour

'rhodamine B chemiluminescence':{'abbreviation':'RBC', 'sampling_type':'low volume continuous', 'measured_parameters':['O3'], 'qa_accepted_parameters':['O3'], 
    'instruments':{
        'Philips PW9771': {'instrument_further_details':'No manual or specifications available but confirmation it is uses rhodamine B chemiluminescence here:http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf'}
    }
},

#------------------------------------
#chemiluminescence (internal molybdenum converter)

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of NO with ozone (NO+O3 --> NO2*+O2).
#The return to a fundamental electronic state of the excited NO2* molecules is made by luminous radiation in a 600-3000 nm spectrum (NO2* --> NO2 + hv), which can be measured.
#This method is designed specifically to directly measure NO directly (and NOx/NO2 indirectly). O3 can also be measured by inverting the measurement technique (using canister of NO for reaction).

#Typical measurement operation for NO (and NOx/NO2 indirectly):
#1. Sample air is drawn into the reaction cell via two separate (alternating) paths; the NO and NOx channels.
#2. NO in the first path (NO channel) reacts with ozone (from ozone generator) according to the following reaction: NO+O3 --> NO2*+O2.
#3. The excited NO2 molecules quickly returns to the ground state, releasing the excess energy. This release takes the form of a quantum of light (hv). The distribution of wavelengths for these quanta range between 600 and 3000 nm, with a peak at about 1200 nm (NO2* --> NO2 +hv (1200nm))
#4. Luminescence from the reaction is detected using a photomultiplier tube (PMT). Photons enter the PMT and strike a negatively charged photo cathode causing it to emit electrons. These electrons are accelerated by an applied high voltage and multiplied through a sequence of similar acceleration steps (dynodes) until a useable current signal is generated. The more light present (in this case photons given off by the chemiluminescent reaction described above), the more current is produced. Therefore the more NO present in the reaction cell the more current is produced by the PMT.
#5. Before entering the PMT, light is passed through a high-pass optical filter only transparent to wavelengths of light above 645nm. The narrowness of this band of sensitivity allows avoids the potential of extraneous light and radiation that might interfere with the measurement (e.g. some oxides of sulfur can also be chemiluminescent emitters when in contact with O3 but give off light at much shorter wavelengths (usually around 260nm to 480nm)).
#6. All things being constant (temperature, pressure, amount of ozone present, etc.), the relationship between the amount of NO present in the reaction cell and the amount of light emitted from the reaction is very linear. If more NO is present, more IR light is produced. 
#7. In addition, sometimes the excited NO2 collides with other gaseous molecules in the reaction cell chamber or even the molecules of the reaction cell walls and transfers its excess energy to this collision partner (represented by M) without emitting any light at all. In fact, by far the largest portion of the excited NO2 returns to the ground state this way, leaving only a few percent yield of usable chemiluminescence (NO2* + M --> NO2 + M). The probability of a collision between the NO2* molecule and a collision partner M increases proportionally with the reaction cell pressure. This non-radiating collision with the NO2* molecules is usually referred to as third body quenching. Even under the best conditions only about 20% of the NO2 that is formed by the key chemiluminescence reaction is in the excited state. In order to maximize chemiluminescence, the reaction cell is maintained at reduced pressure (thereby reducing the amount of available collision partners) and is supplied with a large, constant excess of ozone (about 3000-5000 ppm) from the internal ozone generator.
#8. The second path (NOx channel) travels through a delay loop and a heated molybdenum converter (at 315degC). The heated molybdenum reacts with NO2 in the sample gas and produces a NO gas and a variety of molybdenum. (xNO2 + yMo --> xNO + MyOz (at 315degC)).
#9. Once the NO2 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#10. By converting the NO2 in the sample gas into NO, the analyser indirectly measures total NOx content of the sample gas (i.e. the NO present + the converted NO2 present). 
#11. By switching the sample gas stream in and out of the molybdenum converter every 6 - 10 seconds, the analyser is able to quasi-continuously measure both the NO and the total NOx content.
#12. The NO2 concentration is not directly measured but calculated by subtracting the known NO content of the sample gas from the known NOx content.

#Known interferences: Water vapour (above 20 ppmv), 3rd body quenching (CO2, SOx), other NOy species converted to NO by molybdenum converter (PAN, ethyl nitrate, ethyl nitrite, HONO, HNO3, methyl nitrate, n-propyl nitrate, n-butyl nitrate, nitrocresol, NH3), other species undergoing chemiluminescence with O3 (SOx).

#Not QA accepted measurement technique for NO2/NOx due to bias in measuring NO2, due to molybdenum converters.

'chemiluminescence (internal molybdenum converter)':{'abbreviation':'CL(IMC)', 'sampling_type':'low volume continuous', 'measured_parameters':['NO','NO2','NOx','O3'], 'qa_accepted_parameters':['NO','O3'], 
    'instruments':{
        'Airpointer-NOx':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'20000.0',                                       'documented_precision':'1.0/1.0%', 'documented_accuracy':'1.0%>=100.0', 'documented_zero_drift':'0.4/day',   'documented_span_drift':'1.0%>=100/day',   'documented_flow_rate':'1.0',     'instrument_further_details':'cannot find converter documentation (assume internal molybdenum converter), see other details here: https://www.airpointer.com/#sc-tabs-1547659231751|sc-tabs-1547659231301', 'measuring_assumption':True},
        'Beckman 952A':          {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://www.tandfonline.com/doi/pdf/10.1080/00022470.1979.10470816?needAccess=true', 'measuring_assumption':True},
        'Bendix 8101-B':         {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://www.tandfonline.com/doi/pdf/10.1080/00022470.1980.10465149', 'measuring_assumption':True},
        'Bendix 8101-C':         {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://www.ce.berkeley.edu/sites/default/files/assets/users/kirchstetter/Kirchstetter%20et%20al%2C%20Nitrous%20Acid%20in%20Motor%20Vehicle%20Exhaust%2C%20EST%201996.pdf', 'measuring_assumption':True},
        'CSI 1600':              {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://nepis.epa.gov/Exe/ZyNET.exe/30003UBG.txt?ZyActionD=ZyDocument&Client=EPA&Index=1986%20Thru%201990&Docs=&Query=&Time=&EndTime=&SearchMethod=1&TocRestrict=n&Toc=&TocEntry=&QField=&QFieldYear=&QFieldMonth=&QFieldDay=&UseQField=&IntQFieldOp=0&ExtQFieldOp=0&XmlQuery=&File=D%3A%5CZYFILES%5CINDEX%20DATA%5C86THRU90%5CTXT%5C00000006%5C30003UBG.txt&User=ANONYMOUS&Password=anonymous&SortMethod=h%7C-&MaximumDocuments=1&FuzzyDegree=0&ImageQuality=r75g8/r75g8/x150y150g16/i425&Display=hpfr&DefSeekPage=x&SearchBack=ZyActionL&Back=ZyActionS&BackDesc=Results%20page&MaximumPages=1&ZyEntry=1', 'measuring_assumption':True},
        'Dasibi 2108':           {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter:https://books.google.es/books?id=-I7pCAAAQBAJ&pg=PA372&lpg=PA372&dq=DASIBI+2108+manual&source=bl&ots=ZwNR_ShKjO&sig=yUbYYfnJXqD8eAJovJ6pNjaa468&hl=en&sa=X&ved=0ahUKEwj_jL_SxqzbAhXrB8AKHXvKDkUQ6AEIWTAI#v=onepage&q=DASIBI%202108%20manual&f=false', 'measuring_assumption':True},
        'DKK-TOA GLN-114E':      {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'cannot find converter documentation but states instrument is previous version of DKK-TOA GLN-314E here: https://www.toadkk.co.jp/english/product/env/saleend/ap.html. This uses internal molybdenum converter, so assume same here', 'measuring_assumption':True},
        'DKK-TOA GLN-314E':      {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'2000.0',                                                                                                                'documented_zero_drift':'0.6/day',   'documented_span_drift':'1.0%/day',                                          'instrument_manual_name':'DKK-TOA_GLN-314E_Specs.pdf'},
        'Ecophysics CLD700AL':   {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':{'NO':'100000.0','NO2':'1000.0','NOx':'1000.0'},                                    'documented_accuracy':'1.0%',        'documented_zero_drift':'0.0',                                                  'documented_flow_rate':'0.7',     'instrument_manual_name':'Ecophysics_CLDAL700_specs.pdf'},
        'Ecophysics CLD770AL':   {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but information is found about converter here:https://www.atmos-chem-phys-discuss.net/acp-2014-187/acpd-14-10135-2014.pdf'},
        'Ecotech EC9841T':       {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'2000.0',                                        'documented_precision':'0.1/0.5%',                                      'documented_zero_drift':'0.1/day',   'documented_span_drift':'1.0%/day',        'documented_flow_rate':'0.64',    'instrument_manual_name':'Ecotech_EC9841T_Specs.pdf', 'measuring_instrument_further_details':'Converter information not given in Ecotech_EC9841T_Specs.pdf, but found online: https://books.google.es/books?id=PknpYLJjsgIC&pg=PT101&lpg=PT101&dq=ecotech+EC9841T+molybdenum&source=bl&ots=p34J8Jais1&sig=smuKt3IiSRlcdvCjufYUF1Wh2iQ&hl=en&sa=X&ved=0ahUKEwiSwNqUoefaAhVMORQKHX5kBrcQ6AEIKTAA#v=onepage&q=ecotech%20EC9841T%20molybdenum&f=false'},
        'Ecotech Serinus40':     {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'20000.0',                                       'documented_precision':'0.4/0.5%', 'documented_accuracy':'1.0%',        'documented_zero_drift':'0.4/day',   'documented_span_drift':'1.0%/7days',      'documented_flow_rate':'0.3',     'instrument_manual_name':'Ecotech_Serinus40_Manual.pdf'},
        'Environnement SA AC30M':{'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://helda.helsinki.fi/bitstream/handle/10138/1157/cont71.pdf?sequence=3&isAllowed=y', 'measuring_assumption':True},
        'Environnement SA AC31M':{'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.35', 'documented_upper_limit_of_detection':'10000.0',                                                                          'documented_accuracy':'1.0%',        'documented_zero_drift':'1.0/week',  'documented_span_drift':'1.0%/day',        'documented_flow_rate':'0.6',     'instrument_manual_name':'EnvironnementSA_AC31M_Specs.pdf'},
        'Environnement SA AC32M':{'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'50000.0',                                       'documented_precision':'1.0%',     'documented_accuracy':'1.0%',        'documented_zero_drift':'1.0/day',   'documented_span_drift':'1.0%/7days',      'documented_flow_rate':'0.66',    'instrument_manual_name':'EnvironnementSA_AC32M_Manual.pdf'},
        'Environnement SA AC32e':{'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.2',  'documented_upper_limit_of_detection':'10000.0',                                       'documented_precision':'1.0%',     'documented_accuracy':'1.0%',        'documented_zero_drift':'1.0/day',   'documented_span_drift':'1.0/day',         'documented_flow_rate':'0.66',    'instrument_manual_name':'EnvironnementSA_AC32e_Specs.pdf'},
        'Horiba APNA300':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://books.google.es/books?id=qmWXBwAAQBAJ&pg=PA88&lpg=PA88&dq=horiba+apna+300&source=bl&ots=dBQJ0TEai1&sig=plLm6SAUIZwuIhy37JzI8nZXccQ&hl=en&sa=X&ved=2ahUKEwihxOS9uvDfAhWGnhQKHYzoAPcQ6AEwCHoECAIQAQ#v=onepage&q=horiba%20apna%20300&f=false', 'measuring_assumption':True},
        'Horiba APNA350':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: http://www.wseas.us/e-library/conferences/2009/timisoara/SSE2/SSE2-11.pdf', 'measuring_assumption':True},
        'Horiba APNA360':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'10000.0',                                       'documented_precision':'1.0%',     'documented_accuracy':'1.0%',        'documented_zero_drift':'0.5/day',   'documented_span_drift':'0.5/day',         'documented_flow_rate':'0.8',     'instrument_manual_name':'Horiba_APNA360_Manual.pdf'},
        'Horiba APNA370':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'10000.0',                                       'documented_precision':'1.0%',     'documented_accuracy':'1.0%',        'documented_zero_drift':'0.5/day',   'documented_span_drift':'0.5/day',         'documented_flow_rate':'0.8',     'instrument_manual_name':'Horiba_APNA370_Manual.pdf'},
        'MCV 30QL':              {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://studylib.es/doc/4716520/expedient-de-contractació-núm.-tipus-de-contracte--servei...', 'measuring_assumption':True},
        'Meloy NA530R':          {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=MELOY+NA530R&source=bl&ots=-uI9Lv61Zl&sig=Av_D54Gtv1qh8qVrAa4qLde7_cw&hl=en&sa=X&ved=0ahUKEwjg1q2-8uvaAhUGhiwKHcCxAFYQ6AEIKzAB#v=onepage&q=MELOY%20NA530R&f=false', 'measuring_assumption':True},
        'Philips K50034':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://www.google.co.uk/url?sa=t&rct=j&q=&esrc=s&source=web&cd=11&ved=2ahUKEwipieOHk73nAhUd8uAKHdEdA1oQFjAKegQIAxAB&url=http%3A%2F%2Ffbd.beun.edu.tr%2Findex.php%2Fzkufbd%2Farticle%2Fdownload%2F53%2F100&usg=AOvVaw0l7fQAqr1Pealzx0MaZ6d0 ', 'measuring_assumption':True},
        'Philips K50102':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/philipsK50102/view', 'measuring_assumption':True},
        'Philips K50235':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://www.arpae.it/cms3/documenti/_cerca_doc/aria/rimini/rn_valloni_tesi_laurea_mar2010.pdf', 'measuring_assumption':True},
        'Philips PW9762/02':     {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://books.google.es/books?id=T5B9CAAAQBAJ&pg=PA415&lpg=PA415&dq=Philips+PW9762/02&source=bl&ots=cPN-nS_ZqB&sig=hPPtzoZLDIjgkvnmuhNb8mB0AC8&hl=en&sa=X&ved=0ahUKEwj0it6zuqzbAhUMCMAKHdI-AikQ6AEIKzAB#v=onepage&q=Philips%20PW9762%2F02&f=false', 'measuring_assumption':True},
        'Seres 2000G':           {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':'20000.0',                                                                          'documented_accuracy':'1.0%',        'documented_zero_drift':'1.0%/week', 'documented_span_drift':'1.0%/week',       'documented_flow_rate':'0.5',     'instrument_manual_name':'Seres_2000G_Specs.pdf'},
        'SIR S-5012':            {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'Converter information not given in a manual or specs, but found online: http://www.ecomission.co.za/products/model_s_5012.htm'},
        'Teledyne 200':          {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/API200/view', 'measuring_assumption':True},
        'Teledyne 200A':         {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'20000.0',                                       'documented_precision':'0.5%',     'documented_accuracy':'1.0%',        'documented_zero_drift':'0.5/day',   'documented_span_drift':'0.5%/7days',      'documented_flow_rate':'0.5',     'instrument_manual_name':'Teledyne_200A_Manual.pdf'},
        'Teledyne 200AU':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'2000.0',                                        'documented_precision':'0.5%>=50', 'documented_accuracy':'1.0%/0.1',    'documented_zero_drift':'0.1/day',   'documented_span_drift':'0.05/0.5%/7days', 'documented_flow_rate':'1.0',     'instrument_manual_name':'Teledyne_200AU_Manual.pdf'},
        'Teledyne 200E':         {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'20000.0',                                       'documented_precision':'0.5%',     'documented_accuracy':'1.0%',        'documented_zero_drift':'0.5/day',   'documented_span_drift':'0.5%/7days',      'documented_flow_rate':'0.5',     'instrument_manual_name':'Teledyne_200E_Manual.pdf'},
        'Teledyne 200EU':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'2000.0',                                        'documented_precision':'0.5%>=50', 'documented_accuracy':'1.0%/0.1',    'documented_zero_drift':'0.1/day',   'documented_span_drift':'0.05/0.5%/day',   'documented_flow_rate':'1.0',     'instrument_manual_name':'Teledyne_200EU_addendum.pdf'}, 
        'Teledyne 265E':         {'measured_parameters':['O3'],             'qa_accepted_parameters':['O3'], 'documented_lower_limit_of_detection':'0.3',  'documented_upper_limit_of_detection':'2000.0',                                        'documented_precision':'0.5%',     'documented_accuracy':'1.0%',        'documented_zero_drift':'0.5',       'documented_span_drift':'0.5%',            'documented_flow_rate':'0.5',     'instrument_manual_name':'Teledyne_265E_Addendum.pdf'},
        'Teledyne ML2041':       {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://cambridgeshireinsight.org.uk/wp-content/uploads/2017/12/Cambs-City-Council-Air-Quality-Annual-Status-Report-2015.pdf', 'measuring_assumption':True},
        'Teledyne ML8440E':      {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but information is found about converter here: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=21&ved=0ahUKEwjLg9-In-raAhXFshQKHdzgBfE4FBAWCCgwAA&url=https%3A%2F%2Fnagoya.repo.nii.ac.jp%2F%3Faction%3Drepository_uri%26item_id%3D20273%26file_id%3D17%26file_no%3D1&usg=AOvVaw1jJ6Phm7rzGeqXxCppKjS7'},
        'Teledyne ML8840':       {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':'200.0',                                                                                                                                                                                                                                   'instrument_further_details':'No manual or specifications available but information is found about converter here: https://books.google.es/books?id=OFfsCAAAQBAJ&pg=PA40&lpg=PA40&dq=Monitor+Labs+8840&source=bl&ots=tm9D_6gJ-t&sig=bwhco8GlW_p1yft9aESRbijGUN0&hl=en&sa=X&ved=0ahUKEwjU7KravKzbAhUjB8AKHZRJC-YQ6AEIYDAJ#v=onepage&q=Monitor%20Labs%208840&f=false ; details about limit of detection:https://books.google.es/books?id=G0EYUql3_98C&pg=PA283&lpg=PA283&dq=thermo+43a&source=bl&ots=k2RUJXlAVJ&sig=yIxG08AoSyuqedtS6dp-O7O04X4&hl=en&sa=X&ved=0ahUKEwj17qqB3OfbAhVCVRQKHanFDlo4ChDoAQgoMAA#v=onepage&q=thermo%2043a&f=false'},
        'Teledyne ML8841':       {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter:https://books.google.es/books?id=qmWXBwAAQBAJ&pg=PA88&lpg=PA88&dq=Monitor+Labs+8841&source=bl&ots=dBOJ1Rybm3&sig=VOhLyY-EPom-jDAcT7ST9qsQN7M&hl=en&sa=X&ved=0ahUKEwjB077wxKzbAhUjCMAKHXN_BE4Q6AEIVTAJ#v=onepage&q=Monitor%20Labs%208841&f=false', 'measuring_assumption':True},
        'Teledyne ML9841':       {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: https://www.arb.ca.gov/aaqm/qa/qa-manual/vol4/measurement_of_nox.pdf',              'measuring_assumption':True},
        'Teledyne ML9841A':      {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/ML8941A/view', 'measuring_assumption':True},
        'Teledyne ML9841B':      {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/ML9841B/view', 'measuring_assumption':True},
        'Teledyne ML9841T':      {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter: http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/ML9841T/view', 'measuring_assumption':True},
        'Teledyne T200':         {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'20000.0',                                       'documented_precision':'0.5%>=50', 'documented_accuracy':'1.0%',        'documented_zero_drift':'0.5/day',   'documented_span_drift':'0.5%/day',        'documented_flow_rate':'0.5',     'instrument_manual_name':'Teledyne_T200_Manual.pdf'},
        'Teledyne T200U':        {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'2000.0',                                        'documented_precision':'0.5%>=5' , 'documented_accuracy':'1.0%/0.1',    'documented_zero_drift':'0.1/day',   'documented_span_drift':'0.05/0.5%/day',   'documented_flow_rate':'1.1',     'instrument_manual_name':'Teledyne_T200U_Addendum.pdf'},
        'Teledyne T204-NOx':     {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'20000.0',                                       'documented_precision':'0.5%>=50', 'documented_accuracy':'1.0%',        'documented_zero_drift':'0.5/day',   'documented_span_drift':'0.5%/day',        'documented_flow_rate':'0.5',     'instrument_manual_name':'Teledyne_T204_Manual.pdf'},
        'Teledyne T265':         {'measured_parameters':['O3'],             'qa_accepted_parameters':['O3'], 'documented_lower_limit_of_detection':'0.3',  'documented_upper_limit_of_detection':'2000.0',                                        'documented_precision':'0.5%',     'documented_accuracy':'1.0%',        'documented_zero_drift':'0.5/day',   'documented_span_drift':'0.5%/day',        'documented_flow_rate':'0.5',     'instrument_manual_name':'Teledyne_T265_Addendum.pdf'},
        'Thermo 14B':            {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_manual_name':'Thermo_14B_Operation.pdf'},
        'Thermo 14B/E':          {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_manual_name':'Thermo_14B/E_Operation.pdf'},
        'Thermo 14D/E':          {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses chemiluminescence here, assume internal molybdenum converter:https://books.google.es/books?id=fXysxVKVzyQC&pg=PA137&lpg=PA137&dq=Thermo+14D/E&source=bl&ots=8ks5Vf7mrS&sig=r0-I08UWw6bVqUPPbIavs2GlIHk&hl=en&sa=X&ved=0ahUKEwiiot2juKzbAhXKKMAKHbE3D-AQ6AEINzAB#v=onepage&q=Thermo%2014D%2FE&f=false', 'measuring_assumption':True},
        'Thermo 42':             {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'cannot find any information. AQS states it as continuous for NO2 which is almost always chemiluminescence with internal molybdenum converter. Assume this.', 'measuring_assumption':True},
        'Thermo 42C':            {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'50.0', 'documented_upper_limit_of_detection':'5000000.0',                                                                        'documented_accuracy':'1.0%',        'documented_zero_drift':'50.0/day',  'documented_span_drift':'1.0%/day',        'documented_flow_rate':'0.025',   'instrument_manual_name':'Thermo_42C_Manual.pdf'},
        'Thermo 42C-TL':         {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'200.0',                                                                            'documented_accuracy':'1.0%',        'documented_zero_drift':'0.0',       'documented_span_drift':'1.0%/day',        'documented_flow_rate':'1.0',     'instrument_manual_name':'Thermo_42C-TL_Manual.pdf'},
        'Thermo 42i':            {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.4',  'documented_upper_limit_of_detection':'100000.0',                                                                         'documented_accuracy':'1.0%',        'documented_zero_drift':'0.4/day',   'documented_span_drift':'1.0%/day',        'documented_flow_rate':'0.6-0.8', 'instrument_manual_name':'Thermo_42i_Manual.pdf'},
        'Thermo 42i-LS':         {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'10.0', 'documented_upper_limit_of_detection':'500000.0',                                                                         'documented_accuracy':'1.0%',        'documented_zero_drift':'5.0/day',   'documented_span_drift':'1.0%/day',        'documented_flow_rate':'0.025',   'instrument_manual_name':'Thermo_42i-LS_Manual.pdf'},
        'Thermo 42i-TL':         {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'], 'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'200.0',                                                                            'documented_accuracy':'1.0%',        'documented_zero_drift':'0.0',       'documented_span_drift':'1.0%/day',        'documented_flow_rate':'1.0',     'instrument_manual_name':'Thermo_42i-TL_Manual.pdf'},
        'Thermo 42S':            {'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO'],                                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but information is found about converter here: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/98JD01202'}
    }
},

#------------------------------------
#chemiluminescence (external molybdenum converter)

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of NO with ozone (NO+O3 --> NO2*+O2).
#The return to a fundamental electronic state of the excited NO2* molecules is made by luminous radiation in a 600-3000 nm spectrum (NO2* --> NO2 + hv), which can be measured.
#This method is designed specifically to directly measure NO (and NOy indirectly).

#Typical measurement operation for NO (and NOy indirectly)
#1. Sample air is drawn into the reaction cell via two separate (alternating) paths; the NO and NOy channels.
#2. NO in the first path (NO channel) reacts with ozone (from ozone generator) according to the following reaction: NO+O3 --> NO2*+O2.
#3. The excited NO2 molecules quickly returns to the ground state, releasing the excess energy. This release takes the form of a quantum of light (hv). The distribution of wavelengths for these quanta range between 600 and 3000 nm, with a peak at about 1200 nm (NO2* --> NO2 +hv (1200nm))
#4. Luminescence from the reaction is detected using a photomultiplier tube (PMT). Photons enter the PMT and strike a negatively charged photo cathode causing it to emit electrons. These electrons are accelerated by an applied high voltage and multiplied through a sequence of similar acceleration steps (dynodes) until a useable current signal is generated. The more light present (in this case photons given off by the chemiluminescent reaction described above), the more current is produced. Therefore the more NO present in the reaction cell the more current is produced by the PMT.
#5. Before entering the PMT, light is passed through a high-pass optical filter only transparent to wavelengths of light above 645nm. The narrowness of this band of sensitivity allows avoids the potential of extraneous light and radiation that might interfere with the measurement (e.g. some oxides of sulfur can also be chemiluminescent emitters when in contact with O3 but give off light at much shorter wavelengths (usually around 260nm to 480nm)).
#6. All things being constant (temperature, pressure, amount of ozone present, etc.), the relationship between the amount of NO present in the reaction cell and the amount of light emitted from the reaction is very linear. If more NO is present, more IR light is produced. 
#7. In addition, sometimes the excited NO2 collides with other gaseous molecules in the reaction cell chamber or even the molecules of the reaction cell walls and transfers its excess energy to this collision partner (represented by M) without emitting any light at all. In fact, by far the largest portion of the excited NO2 returns to the ground state this way, leaving only a few percent yield of usable chemiluminescence (NO2* + M --> NO2 + M). The probability of a collision between the NO2* molecule and a collision partner M increases proportionally with the reaction cell pressure. This non-radiating collision with the NO2* molecules is usually referred to as third body quenching. Even under the best conditions only about 20% of the NO2 that is formed by the key chemiluminescence reaction is in the excited state. In order to maximize chemiluminescence, the reaction cell is maintained at reduced pressure (thereby reducing the amount of available collision partners) and is supplied with a large, constant excess of ozone (about 3000-5000 ppm) from the internal ozone generator.
#8. The second path (NOy channel) before entering the instrument enters an externally mounted heated molybdenum converter (at 315degC), with less than 10cm tubing between the converter and sample entry. The heated molybdenum reacts with NOy in the sample gas and produces a NO gas and a variety of molybdenum (3NOy + Mo --> 3NO + MoO3 (at 315degC)). Minimising the transit time and surface contact area between the sample inlet and converter enables the conversion of labile components of NOy.
#9. Once the NOy in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#10. By converting the NOy in the sample gas into NO, the analyser indirectly measures total NOy content of the sample gas (i.e. the NO present + the converted NO2 present). 
#11. By switching the sample gas stream in and out of the molybdenum converter every 6 - 10 seconds, the analyser is able to quasi-continuously measure both the NO and the total NOy content.

#Known interferences: Water vapour (above 20 ppmv), 3rd body quenching (CO2, SOx), other species undergoing chemiluminescence with O3 (SOx).

'chemiluminescence (external molybdenum converter)':{'abbreviation':'CL(EMC)', 'sampling_type':'low volume continuous', 'measured_parameters':['NO','NOy'], 'qa_accepted_parameters':['NO','NOy'], 
    'instruments':{
        'Ecotech EC9843':    {'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'2000.0', 'documented_precision':'0.05/0.5%', 'documented_accuracy':'1.0%',     'documented_zero_drift':'0.1/day', 'documented_span_drift':'1.0%/day',      'documented_flow_rate':'0.64', 'instrument_manual_name':'Ecotech_EC9843_Specs.pdf, Ecotech_EC9843_Manual.pdf'},
        'Teledyne 200EU-NOy':{'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'2000.0', 'documented_precision':'0.5%>=50',  'documented_accuracy':'1.0%/0.1', 'documented_zero_drift':'0.1/day', 'documented_span_drift':'0.05/0.5%/day', 'documented_flow_rate':'1.0',  'instrument_manual_name':'Teledyne_200EU-NOy_addendum.pdf'},
        'Thermo 42C-NOy':    {'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'200.0',                                      'documented_accuracy':'1.0%',     'documented_zero_drift':'0.0',     'documented_span_drift':'1.0%/day',      'documented_flow_rate':'1.0',  'instrument_manual_name':'Thermo_42C-NOy_Manual.pdf'},
        'Thermo 42i-NOy':    {'documented_lower_limit_of_detection':'0.05', 'documented_upper_limit_of_detection':'1000.0',                                     'documented_accuracy':'1.0%',     'documented_zero_drift':'0.0',     'documented_span_drift':'1.0%/day',      'documented_flow_rate':'1.0',  'instrument_manual_name':'Thermo_42i-NOy_Manual.pdf'}
    }
},

#------------------------------------
#chemiluminescence (internal photolytic converter) (CL(IPC))

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of NO with ozone (NO+O3 --> NO2*+O2).
#The return to a fundamental electronic state of the excited NO2* molecules is made by luminous radiation in a 600-3000 nm spectrum (NO2* --> NO2 + hv), which can be measured.
#This method is designed specifically to directly measure NO (and NOx/NO2 indirectly). 

#Typical measurement operation for NO (and NOx/NO2 indirectly):
#1. Sample air is drawn into the reaction cell via two separate (alternating) paths; the NO and NOx channels.
#2. NO in the first path (NO channel) reacts with ozone (from ozone generator) according to the following reaction: NO+O3 --> NO2*+O2.
#3. The excited NO2 molecules quickly returns to the ground state, releasing the excess energy. This release takes the form of a quantum of light (hv). The distribution of wavelengths for these quanta range between 600 and 3000 nm, with a peak at about 1200 nm (NO2* --> NO2 +hv (1200nm))
#4. Luminescence from the reaction is detected using a photomultiplier tube (PMT). Photons enter the PMT and strike a negatively charged photo cathode causing it to emit electrons. These electrons are accelerated by an applied high voltage and multiplied through a sequence of similar acceleration steps (dynodes) until a useable current signal is generated. The more light present (in this case photons given off by the chemiluminescent reaction described above), the more current is produced. Therefore the more NO present in the reaction cell the more current is produced by the PMT.
#5. Before entering the PMT, light is passed through a high-pass optical filter only transparent to wavelengths of light above 645nm. The narrowness of this band of sensitivity allows avoids the potential of extraneous light and radiation that might interfere with the measurement (e.g. some oxides of sulfur can also be chemiluminescent emitters when in contact with O3 but give off light at much shorter wavelengths (usually around 260nm to 480nm)).
#6. All things being constant (temperature, pressure, amount of ozone present, etc.), the relationship between the amount of NO present in the reaction cell and the amount of light emitted from the reaction is very linear. If more NO is present, more IR light is produced. 
#7. In addition, sometimes the excited NO2 collides with other gaseous molecules in the reaction cell chamber or even the molecules of the reaction cell walls and transfers its excess energy to this collision partner (represented by M) without emitting any light at all. In fact, by far the largest portion of the excited NO2 returns to the ground state this way, leaving only a few percent yield of usable chemiluminescence (NO2* + M --> NO2 + M). The probability of a collision between the NO2* molecule and a collision partner M increases proportionally with the reaction cell pressure. This non-radiating collision with the NO2* molecules is usually referred to as third body quenching. Even under the best conditions only about 20% of the NO2 that is formed by the key chemiluminescence reaction is in the excited state. In order to maximize chemiluminescence, the reaction cell is maintained at reduced pressure (thereby reducing the amount of available collision partners) and is supplied with a large, constant excess of ozone (about 3000-5000 ppm) from the internal ozone generator.
#8. The second path (NOx channel) travels through a delay loop and a photolytic converter. As the sample gas passes through the converter chamber it is exposed to blue light at specific wavelengths (350-420 nm) from a photolytic light source (blue LEDs, metal halide). This selectively converts NO2 to NO with negligible radiant heating or interference from other gases (NO2 +hv (350-420nm) --> NO + O).
#9. Once the NO2 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#10. By converting the NO2 in the sample gas into NO, the analyser indirectly measures total NOx content of the sample gas (i.e. the NO present + the converted NO2 present). 
#11. By switching the sample gas stream in and out of the photolytic converter every 6 - 10 seconds, the analyser is able to quasi-continuously measure both the NO and the total NOx content.
#12. The NO2 concentration is not directly measured but calculated by subtracting the known NO content of the sample gas from the known NOx content.

#Known interferences: Water vapour (above 20 ppmv), 3rd body quenching (CO2, SOx), thermal decomposition of PAN to NO2 within the photolysis cell

'chemiluminescence (internal photolytic converter)':{'abbreviation':'CL(IPC)', 'sampling_type':'low volume continuous', 'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO','NO2','NOx'], 
    'instruments':{
        'Ecophysics CraNOx':{'documented_lower_limit_of_detection':'0.025',                               'documented_upper_limit_of_detection':'1000.0',                                    'documented_accuracy':'1.0%',                                                                            'documented_flow_rate':'2.7', 'instrument_manual_name':'Ecophysics_CraNOx_Specs.pdf'},
        'Teledyne T200P':   {'documented_lower_limit_of_detection':'0.4',                                 'documented_upper_limit_of_detection':'4000.0', 'documented_precision':'0.5%>=50', 'documented_accuracy':'1.0%',     'documented_zero_drift':'0.5/day', 'documented_span_drift':'0.5%/day', 'documented_flow_rate':'0.5', 'instrument_manual_name':'Teledyne_T200P_Manual.pdf'},
        'Teledyne T200UP':  {'documented_lower_limit_of_detection':{'NO':'0.05','NO2':'0.1','NOx':'0.1'}, 'documented_upper_limit_of_detection':'2000.0', 'documented_precision':'0.5%>=5',  'documented_accuracy':'1.0%/0.1', 'documented_zero_drift':'0.1/day', 'documented_span_drift':'0.5%/day', 'documented_flow_rate':'1.1', 'instrument_manual_name':'Teledyne_T200UP_addendum.pdf'}
    }
},

#------------------------------------
#chemiluminescence (internal molybdenum and quartz converters) (CL(IMQC))

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of NO with ozone (NO+O3 --> NO2*+O2).
#The return to a fundamental electronic state of the excited NO2* molecules is made by luminous radiation in a 600-3000 nm spectrum (NO2* --> NO2 + hv), which can be measured.
#This method is designed specifically to directly measure NO (and NH3/NOx/NO2 indirectly). 

#Typical measurement operation for NO (and NH3/NOx/NO2 indirectly)
#1. Sample air is drawn into the reaction cell via three separate (alternating) paths; the NO, TNx and NOx channels.
#2. NO in the first path (NO channel) reacts with ozone (from ozone generator) according to the following reaction: NO+O3 --> NO2*+O2.
#3. The excited NO2 molecules quickly returns to the ground state, releasing the excess energy. This release takes the form of a quantum of light (hv). The distribution of wavelengths for these quanta range between 600 and 3000 nm, with a peak at about 1200 nm (NO2* --> NO2 +hv (1200nm))
#4. Luminescence from the reaction is detected using a photomultiplier tube (PMT). Photons enter the PMT and strike a negatively charged photo cathode causing it to emit electrons. These electrons are accelerated by an applied high voltage and multiplied through a sequence of similar acceleration steps (dynodes) until a useable current signal is generated. The more light present (in this case photons given off by the chemiluminescent reaction described above), the more current is produced. Therefore the more NO present in the reaction cell the more current is produced by the PMT.
#5. Before entering the PMT, light is passed through a high-pass optical filter only transparent to wavelengths of light above 645nm. The narrowness of this band of sensitivity allows avoids the potential of extraneous light and radiation that might interfere with the measurement (e.g. some oxides of sulfur can also be chemiluminescent emitters when in contact with O3 but give off light at much shorter wavelengths (usually around 260nm to 480nm)).
#6. All things being constant (temperature, pressure, amount of ozone present, etc.), the relationship between the amount of NO present in the reaction cell and the amount of light emitted from the reaction is very linear. If more NO is present, more IR light is produced. 
#7. In addition, sometimes the excited NO2 collides with other gaseous molecules in the reaction cell chamber or even the molecules of the reaction cell walls and transfers its excess energy to this collision partner (represented by M) without emitting any light at all. In fact, by far the largest portion of the excited NO2 returns to the ground state this way, leaving only a few percent yield of usable chemiluminescence (NO2* + M --> NO2 + M). The probability of a collision between the NO2* molecule and a collision partner M increases proportionally with the reaction cell pressure. This non-radiating collision with the NO2* molecules is usually referred to as third body quenching. Even under the best conditions only about 20% of the NO2 that is formed by the key chemiluminescence reaction is in the excited state. In order to maximize chemiluminescence, the reaction cell is maintained at reduced pressure (thereby reducing the amount of available collision partners) and is supplied with a large, constant excess of ozone (about 3000-5000 ppm) from the internal ozone generator.
#8. In the second path (TNx channel), sample air is passed through an internal heated quartz converter (at 980 degC). The heated quartz reacts with NH3 and NO2 in the sample gas and produces a NO gas (4NH3 + 5O2 --> 4NO + 6H2O).
#9. Once the NO2 and NH3 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#10. By converting the NO2 and NH3 in the sample gas into NO, the analyser indirectly measures total TNx content of the sample gas. 
#11. The third path (NOx channel) travels through a delay loop and a heated molybdenum converter (at 315degC). The heated molybdenum reacts with NO2 in the sample gas and produces a NO gas and a variety of molybdenum. (xNO2 + yMo --> xNO + MyOz (at 315degC)).
#12. Once the NO2 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#13. By converting the NO2 in the sample gas into NO, the analyser measures total NOx content of the sample gas (i.e. the NO present + the converted NO2 present). 
#14. By switching the sample gas stream between the 3 streams frequently, the analyser is able to quasi-continuously measure both the NO, the total NOx content and TNx.
#15. NO2 is indirectly calculated by subtracting the known NO content of the sample gas from the known NOx content.
#16. NH3 is indirectly calculated by subtracting the known NOx content of the sample gas from the known TNx content.

#Known interferences: Water vapour (above 20 ppmv), 3rd body quenching (CO2, SOx), other NOy species converted to NO by molybdenum converter (PAN, ethyl nitrate, ethyl nitrite, HONO, HNO3, methyl nitrate, n-propyl nitrate, n-butyl nitrate, nitrocresol, NH3), other species undergoing chemiluminescence with O3 (SOx).

#Not QA accepted measurement technique for NO2/NOx due to bias in measuring NO2, due to molybdenum converters.

'chemiluminescence (internal molybdenum and quartz converters)':{'abbreviation':'CL(IMQC)', 'sampling_type':'low volume continuous', 'measured_parameters':['NH3','NO','NO2','NOx'], 'qa_accepted_parameters':['NH3','NO'], 
    'instruments':{
        'Ecotech EC9842': {'documented_lower_limit_of_detection':'0.5', 'documented_upper_limit_of_detection':'2000.0', 'documented_precision':'0.5%', 'documented_accuracy':'1.0%', 'documented_zero_drift':'1.0/day', 'documented_span_drift':'0.5%/day', 'documented_resolution':'0.1',  'documented_flow_rate':'0.000355', 'instrument_manual_name':'Ecotech_EC9842_Specs.pdf'},
        'Ecotech EC9842T':{'documented_lower_limit_of_detection':'0.4',                                                 'documented_precision':'2.0%',                               'documented_zero_drift':'0.1/day', 'documented_span_drift':'1.0%/day',                                                                    'instrument_further_details':'Details found online:http://aep.alberta.ca/air/legislation-and-policy/air-monitoring-directive/documents/AmbientAirMonitoringSpecifications-Sep2014.pdf'}
    }
},


#------------------------------------
#chemiluminescence (internal molybdenum converter and external quartz converter) (CL(IMC-EQC))

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of NO with ozone (NO+O3 --> NO2*+O2).
#The return to a fundamental electronic state of the excited NO2* molecules is made by luminous radiation in a 600-3000 nm spectrum (NO2* --> NO2 + hv), which can be measured.
#This method is designed specifically to directly measure NO (and NH3/NOx/NO2 indirectly). 

#Typical measurement operation for NO (and NH3/NOx/NO2 indirectly)
#1. Sample air is drawn into the reaction cell via three separate (alternating) paths; the NO, TNx and NOx channels.
#2. NO in the first path (NO channel) reacts with ozone (from ozone generator) according to the following reaction: NO+O3 --> NO2*+O2.
#3. The excited NO2 molecules quickly returns to the ground state, releasing the excess energy. This release takes the form of a quantum of light (hv). The distribution of wavelengths for these quanta range between 600 and 3000 nm, with a peak at about 1200 nm (NO2* --> NO2 +hv (1200nm))
#4. Luminescence from the reaction is detected using a photomultiplier tube (PMT). Photons enter the PMT and strike a negatively charged photo cathode causing it to emit electrons. These electrons are accelerated by an applied high voltage and multiplied through a sequence of similar acceleration steps (dynodes) until a useable current signal is generated. The more light present (in this case photons given off by the chemiluminescent reaction described above), the more current is produced. Therefore the more NO present in the reaction cell the more current is produced by the PMT.
#5. Before entering the PMT, light is passed through a high-pass optical filter only transparent to wavelengths of light above 645nm. The narrowness of this band of sensitivity allows avoids the potential of extraneous light and radiation that might interfere with the measurement (e.g. some oxides of sulfur can also be chemiluminescent emitters when in contact with O3 but give off light at much shorter wavelengths (usually around 260nm to 480nm)).
#6. All things being constant (temperature, pressure, amount of ozone present, etc.), the relationship between the amount of NO present in the reaction cell and the amount of light emitted from the reaction is very linear. If more NO is present, more IR light is produced. 
#7. In addition, sometimes the excited NO2 collides with other gaseous molecules in the reaction cell chamber or even the molecules of the reaction cell walls and transfers its excess energy to this collision partner (represented by M) without emitting any light at all. In fact, by far the largest portion of the excited NO2 returns to the ground state this way, leaving only a few percent yield of usable chemiluminescence (NO2* + M --> NO2 + M). The probability of a collision between the NO2* molecule and a collision partner M increases proportionally with the reaction cell pressure. This non-radiating collision with the NO2* molecules is usually referred to as third body quenching. Even under the best conditions only about 20% of the NO2 that is formed by the key chemiluminescence reaction is in the excited state. In order to maximize chemiluminescence, the reaction cell is maintained at reduced pressure (thereby reducing the amount of available collision partners) and is supplied with a large, constant excess of ozone (about 3000-5000 ppm) from the internal ozone generator.
#8. In the second path (TNx channel), sample air is passed through an externally mounted heated quartz converter (at 980 degC). The heated quartz reacts with NH3 and NO2 in the sample gas and produces a NO gas (4NH3 + 5O2 --> 4NO + 6H2O).
#9. Once the NO2 and NH3 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#10. By converting the NO2 and NH3 in the sample gas into NO, the analyser indirectly measures total TNx content of the sample gas. 
#11. The third path (NOx channel) travels through a delay loop and a heated molybdenum converter (at 315degC). The heated molybdenum reacts with NO2 in the sample gas and produces a NO gas and a variety of molybdenum. (xNO2 + yMo --> xNO + MyOz (at 315degC)).
#12. Once the NO2 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#13. By converting the NO2 in the sample gas into NO, the analyser measures total NOx content of the sample gas (i.e. the NO present + the converted NO2 present). 
#14. By switching the sample gas stream between the 3 streams frequently, the analyser is able to quasi-continuously measure both the NO, the total NOx content and TNx.
#15. NO2 is indirectly calculated by subtracting the known NO content of the sample gas from the known NOx content.
#16. NH3 is indirectly calculated by subtracting the known NOx content of the sample gas from the known TNx content.

#Known interferences: Water vapour (above 20 ppmv), 3rd body quenching (CO2, SOx), other NOy species converted to NO by molybdenum converter (PAN, ethyl nitrate, ethyl nitrite, HONO, HNO3, methyl nitrate, n-propyl nitrate, n-butyl nitrate, nitrocresol, NH3), other species undergoing chemiluminescence with O3 (SOx).

#Not QA accepted measurement technique for NO2/NOx due to bias in measuring NO2, due to molybdenum converters.

'chemiluminescence (internal molybdenum converter and external quartz converter)':{'abbreviation':'CL(IMC-EQC)', 'sampling_type':'low volume continuous', 'measured_parameters':['NH3','NO','NO2','NOx'], 'qa_accepted_parameters':['NH3','NO'], 
    'instruments':{
        'Environnement SA AC32e-CNH3':{'documented_lower_limit_of_detection':{'NO':'0.2','NO2':'0.2','NOx':'0.2','NH3':'1.0'}, 'documented_upper_limit_of_detection':{'NO':'10000.0','NO2':'10000.0','NOx':'10000.0','NH3':'1000.0'}, 'documented_precision':'1.0%', 'documented_accuracy':'1.0%', 'documented_zero_drift':'1.0/day', 'documented_span_drift':'1.0/7days', 'documented_flow_rate':'0.66', 'instrument_manual_name':'EnvironnementSA_AC32e-CNH3_Specs.pdf'}
    }
},  

   
#------------------------------------
#chemiluminescence (internal molybdenum and stainless steel converters) (CL(IMSC))

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of NO with ozone (NO+O3 --> NO2*+O2).
#The return to a fundamental electronic state of the excited NO2* molecules is made by luminous radiation in a 600-3000 nm spectrum (NO2* --> NO2 + hv), which can be measured.
#This method is designed specifically to directly measure NO (and NH3/NOx/NO2 indirectly). 

#Typical measurement operation for NO (and NH3/NOx/NO2 indirectly)
#1. Sample air is drawn into the reaction cell via three separate (alternating) paths; the NO, TNx and NOx channels.
#2. NO in the first path (NO channel) reacts with ozone (from ozone generator) according to the following reaction: NO+O3 --> NO2*+O2.
#3. The excited NO2 molecules quickly returns to the ground state, releasing the excess energy. This release takes the form of a quantum of light (hv). The distribution of wavelengths for these quanta range between 600 and 3000 nm, with a peak at about 1200 nm (NO2* --> NO2 +hv (1200nm))
#4. Luminescence from the reaction is detected using a photomultiplier tube (PMT). Photons enter the PMT and strike a negatively charged photo cathode causing it to emit electrons. These electrons are accelerated by an applied high voltage and multiplied through a sequence of similar acceleration steps (dynodes) until a useable current signal is generated. The more light present (in this case photons given off by the chemiluminescent reaction described above), the more current is produced. Therefore the more NO present in the reaction cell the more current is produced by the PMT.
#5. Before entering the PMT, light is passed through a high-pass optical filter only transparent to wavelengths of light above 645nm. The narrowness of this band of sensitivity allows avoids the potential of extraneous light and radiation that might interfere with the measurement (e.g. some oxides of sulfur can also be chemiluminescent emitters when in contact with O3 but give off light at much shorter wavelengths (usually around 260nm to 480nm)).
#6. All things being constant (temperature, pressure, amount of ozone present, etc.), the relationship between the amount of NO present in the reaction cell and the amount of light emitted from the reaction is very linear. If more NO is present, more IR light is produced. 
#7. In addition, sometimes the excited NO2 collides with other gaseous molecules in the reaction cell chamber or even the molecules of the reaction cell walls and transfers its excess energy to this collision partner (represented by M) without emitting any light at all. In fact, by far the largest portion of the excited NO2 returns to the ground state this way, leaving only a few percent yield of usable chemiluminescence (NO2* + M --> NO2 + M). The probability of a collision between the NO2* molecule and a collision partner M increases proportionally with the reaction cell pressure. This non-radiating collision with the NO2* molecules is usually referred to as third body quenching. Even under the best conditions only about 20% of the NO2 that is formed by the key chemiluminescence reaction is in the excited state. In order to maximize chemiluminescence, the reaction cell is maintained at reduced pressure (thereby reducing the amount of available collision partners) and is supplied with a large, constant excess of ozone (about 3000-5000 ppm) from the internal ozone generator.
#8. In the second path (TNx channel), sample air is first drawn through an internal heated stainless steel converter (at 825 degC). The heated stainless steel reacts with NH3 and NO2 in the sample gas and produces a NO gas (4NH3 + 5O2 --> 4NO + 6H2O).
#9. Once the NO2 and NH3 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#10. By converting the NO2 and NH3 in the sample gas into NO, the analyser indirectly measures total TNx content of the sample gas. 
#11. The third path (NOx channel) travels through a delay loop and a heated molybdenum converter (at 315degC). The heated molybdenum reacts with NO2 in the sample gas and produces a NO gas and a variety of molybdenum. (xNO2 + yMo --> xNO + MyOz (at 315degC)).
#12. Once the NO2 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#13. By converting the NO2 in the sample gas into NO, the analyser measures total NOx content of the sample gas (i.e. the NO present + the converted NO2 present). 
#14. By switching the sample gas stream between the 3 streams frequently, the analyser is able to quasi-continuously measure both the NO, the total NOx content and TNx.
#15. NO2 is indirectly calculated by subtracting the known NO content of the sample gas from the known NOx content.
#16. NH3 is indirectly calculated by subtracting the known NOx content of the sample gas from the known TNx content.

#Known interferences: Water vapour (above 20 ppmv), 3rd body quenching (CO2, SOx), other NOy species converted to NO by molybdenum converter (PAN, ethyl nitrate, ethyl nitrite, HONO, HNO3, methyl nitrate, n-propyl nitrate, n-butyl nitrate, nitrocresol, NH3), other species undergoing chemiluminescence with O3 (SOx).

#Not QA accepted measurement technique for NO2/NOx due to bias in measuring NO2, due to molybdenum converters.

'chemiluminescence (internal molybdenum and stainless steel converters)':{'abbreviation':'CL(IMSC)', 'sampling_type':'low volume continuous', 'measured_parameters':['NH3','NO','NO2','NOx'], 'qa_accepted_parameters':['NH3','NO'], 
    'instruments':{
        'Teledyne 201A':{'documented_lower_limit_of_detection':'1.0', 'documented_upper_limit_of_detection':'2000.0', 'documented_precision':'1.0%', 'documented_accuracy':'2.0%', 'documented_zero_drift':'2.0/day', 'documented_span_drift':'1.0%/week', 'documented_flow_rate':'1.0', 'instrument_manual_name':'Teledyne_201A_Manual.pdf'}
    }
},
#------------------------------------
#chemiluminescence (internal molybdenum converter and external stainless steel converter) (CL(IMC-ESC))

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of NO with ozone (NO+O3 --> NO2*+O2).
#The return to a fundamental electronic state of the excited NO2* molecules is made by luminous radiation in a 600-3000 nm spectrum (NO2* --> NO2 + hv), which can be measured.
#This method is designed specifically to directly measure NO (and NH3/NOx/NO2 indirectly). 

#Typical measurement operation for NO (and NH3/NOx/NO2 indirectly)
#1. Sample air is drawn into the reaction cell via three separate (alternating) paths; the NO, TNx and NOx channels.
#2. NO in the first path (NO channel) reacts with ozone (from ozone generator) according to the following reaction: NO+O3 --> NO2*+O2.
#3. The excited NO2 molecules quickly returns to the ground state, releasing the excess energy. This release takes the form of a quantum of light (hv). The distribution of wavelengths for these quanta range between 600 and 3000 nm, with a peak at about 1200 nm (NO2* --> NO2 +hv (1200nm))
#4. Luminescence from the reaction is detected using a photomultiplier tube (PMT). Photons enter the PMT and strike a negatively charged photo cathode causing it to emit electrons. These electrons are accelerated by an applied high voltage and multiplied through a sequence of similar acceleration steps (dynodes) until a useable current signal is generated. The more light present (in this case photons given off by the chemiluminescent reaction described above), the more current is produced. Therefore the more NO present in the reaction cell the more current is produced by the PMT.
#5. Before entering the PMT, light is passed through a high-pass optical filter only transparent to wavelengths of light above 645nm. The narrowness of this band of sensitivity allows avoids the potential of extraneous light and radiation that might interfere with the measurement (e.g. some oxides of sulfur can also be chemiluminescent emitters when in contact with O3 but give off light at much shorter wavelengths (usually around 260nm to 480nm)).
#6. All things being constant (temperature, pressure, amount of ozone present, etc.), the relationship between the amount of NO present in the reaction cell and the amount of light emitted from the reaction is very linear. If more NO is present, more IR light is produced. 
#7. In addition, sometimes the excited NO2 collides with other gaseous molecules in the reaction cell chamber or even the molecules of the reaction cell walls and transfers its excess energy to this collision partner (represented by M) without emitting any light at all. In fact, by far the largest portion of the excited NO2 returns to the ground state this way, leaving only a few percent yield of usable chemiluminescence (NO2* + M --> NO2 + M). The probability of a collision between the NO2* molecule and a collision partner M increases proportionally with the reaction cell pressure. This non-radiating collision with the NO2* molecules is usually referred to as third body quenching. Even under the best conditions only about 20% of the NO2 that is formed by the key chemiluminescence reaction is in the excited state. In order to maximize chemiluminescence, the reaction cell is maintained at reduced pressure (thereby reducing the amount of available collision partners) and is supplied with a large, constant excess of ozone (about 3000-5000 ppm) from the internal ozone generator.
#8. In the second path (TNx channel), sample air is first drawn through an externally mounted heated stainless steel converter (at 750 degC). The heated stainless steel reacts with NH3 and NO2 in the sample gas and produces a NO gas (4NH3 + 5O2 --> 4NO + 6H2O).
#9. Once the NO2 and NH3 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#10. By converting the NO2 and NH3 in the sample gas into NO, the analyser indirectly measures total TNx content of the sample gas. 
#11. The third path (NOx channel) travels through a delay loop and a heated molybdenum converter (at 315degC). The heated molybdenum reacts with NO2 in the sample gas and produces a NO gas and a variety of molybdenum. (xNO2 + yMo --> xNO + MyOz (at 315degC)).
#12. Once the NO2 in the sample gas has been converted to NO, it is routed to the reaction cell where it undergoes the chemiluminescence reaction and is measured by the PMT.
#13. By converting the NO2 in the sample gas into NO, the analyser measures total NOx content of the sample gas (i.e. the NO present + the converted NO2 present). 
#14. By switching the sample gas stream between the 3 streams frequently, the analyser is able to quasi-continuously measure both the NO, the total NOx content and TNx.
#15. NO2 is indirectly calculated by subtracting the known NO content of the sample gas from the known NOx content.
#16. NH3 is indirectly calculated by subtracting the known NOx content of the sample gas from the known TNx content.

#Known interferences: Water vapour (above 20 ppmv), 3rd body quenching (CO2, SOx), other NOy species converted to NO by molybdenum converter (PAN, ethyl nitrate, ethyl nitrite, HONO, HNO3, methyl nitrate, n-propyl nitrate, n-butyl nitrate, nitrocresol, NH3), other species undergoing chemiluminescence with O3 (SOx).

#Not QA accepted measurement technique for NO2/NOx due to bias in measuring NO2, due to molybdenum converters.

'chemiluminescence (internal molybdenum converter and external stainless steel converter)':{'abbreviation':'CL(IMC-ESC)', 'sampling_type':'low volume continuous', 'measured_parameters':['NH3','NO','NO2','NOx'], 'qa_accepted_parameters':['NH3','NO'], 
    'instruments':{
        'Teledyne 201E':  {'documented_lower_limit_of_detection':'1.0', 'documented_upper_limit_of_detection':'2000.0',   'documented_accuracy':{'NH3':'2.0%','NO':'1.0%','NO2':'1.0%','NOx':'1.0%'}, 'documented_zero_drift':'2.0/day', 'documented_span_drift':'1.0%/day', 'documented_flow_rate':'0.5-1.0', 'instrument_manual_name':'Teledyne_201E_Manual.pdf'},
        'Thermo 17C':     {                                                                                                                                                                                                                                                                                    'instrument_further_details':'Converter information found here:https://archive.epa.gov/nrmrl/archive-etv/web/pdf/01_vs_thermo.pdf'},
        'Thermo 17i':     {'documented_lower_limit_of_detection':'1.0', 'documented_upper_limit_of_detection':'100000.0', 'documented_accuracy':'1.0%',                                               'documented_zero_drift':'1.0/day', 'documented_span_drift':'1.0%/day', 'documented_flow_rate':'0.6',     'instrument_manual_name':'Thermo_17i_Manual.pdf'}
    }
},
#------------------------------------
#sulphur chemiluminescence - gas chromatography (SC-GC)

#Chemiluminescence is the emission of light (luminescence), as the result of a chemical reaction. Chemiluminescence occurs as a result of the reaction of SO and ozone (SO+O3 --> SO2*+O2).
#The return to a fundamental electronic state of the excited SO2* molecules is made by luminous radiation in a specific spectrum (SO2* --> SO2 + hv), which can be measured (by a photomultiplier tube).
#The concentration of sample total sulphur is directly proportional to the intensity of light emitted.
#This method can be mixed with gas chromatography (GC) to allow to determination of specific sulphur compounds (i.e. SO2).
#The gas sample is passed through a GC column before being ultimately measured by the photomultiplier detector.

#Gas chromatography (GC) is a method used for separating and analysing compounds that can be vaporized without decomposition.
#A sample solution is injected into a instrument, entering a gas stream which transports the sample (mobile phase) into a separation tube known as the "column".
#The mobile phase is a carrier gas, usually an inert gas such as helium or an unreactive gas such as nitrogen. Helium remains the most commonly used carrier gas in about 90% of instruments although hydrogen is preferred for improved separations.
#The column consists of a microscopic layer of liquid or polymer on an inert solid support a microscopic layer of liquid or polymer on an inert solid support (stationary phase), inside a piece of glass or metal tubing, placed inside a piece of glass or metal tubing.
#Once inside the column, the gaseous compounds being analysed interact with the walls of the column coated with a stationary phase. This causes each compound to elute at a different time, known as the retention time of the compound. 
#The comparison of retention times is what gives GC its analytical usefulness.
#If greater separation of compounds is required, multiple distinct columns can be used for this purpose.

#See details of method for SO2: http://www.nrcresearchpress.com/doi/pdf/10.1139/v84-073, http://www.interline.nl/media/1000035/specifiekzwavelchromatograph.pdf'

'sulphur chemiluminescence - gas chromatography':{'abbreviation':'SC-GC', 'sampling_type':'low volume continuous', 'measured_parameters':['SO2'], 'qa_accepted_parameters':['SO2'], 
    'instruments':{}
},

#------------------------------------
#flame photometric detection (FPD)

#Many elements give characteristic emission when burned in flame. Absorption of energy from the flame allows a ground state atom or molecule to reach an excited state. 
#The atom/molecule may return to the ground state through emission of light (luminescence), which can be subsequently measured (by a photomultiplier tube). This process is a chemiluminescent process (where luminescence occurs as the result of a chemical reaction).
#The concentration of the sample gas is directly proportional to the intensity of light emitted.
#This method has been applied for the measurement of sulphur containing species (i.e. SO2).
#For specific measurement of solely SO2, the sample gas must be scrubbed of other sulphur species prior to measurement, and the photomultiplier detector measures emission centred near 394nm.

#Known interferences: Other sulphur compounds

#See details of method for SO2: https://nepis.epa.gov/Exe/ZyNET.exe/9101A5EI.txt?ZyActionD=ZyDocument&Client=EPA&Index=1976%20Thru%201980&Docs=&Query=&Time=&EndTime=&SearchMethod=1&TocRestrict=n&Toc=&TocEntry=&QField=&QFieldYear=&QFieldMonth=&QFieldDay=&UseQField=&IntQFieldOp=0&ExtQFieldOp=0&XmlQuery=&File=D%3A%5CZYFILES%5CINDEX%20DATA%5C76THRU80%5CTXT%5C00000026%5C9101A5EI.txt&User=ANONYMOUS&Password=anonymous&SortMethod=h%7C-&MaximumDocuments=1&FuzzyDegree=0&ImageQuality=r75g8/r75g8/x150y150g16/i425&Display=hpfr&DefSeekPage=x&SearchBack=ZyActionL&Back=ZyActionS&BackDesc=Results%20page&MaximumPages=1&ZyEntry=22

'flame photometric detection':{'abbreviation':'FPD', 'sampling_type':'low volume continuous', 'measured_parameters':['SO2'], 'qa_accepted_parameters':['SO2'],  
    'instruments':{
        'Bendix 8303':    {'instrument_further_details':'No manual or specifications available but confirmation it is uses flame photometric detection here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'},
        'Meloy SA185-2A': {'instrument_further_details':'No manual or specifications available but confirmation it is uses flame photometric detection here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'},
        'Meloy SA285E':   {'instrument_further_details':'No manual or specifications available but confirmation it is uses flame photometric detection here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'},
        'Teledyne ML8450':{'instrument_further_details':'No manual or specifications available but confirmation it is uses flame photometric detection here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'}
    }
},

#------------------------------------
#flame ionisation detection (FID)

#This method is based on the principle of the generation of an electrical current that is proportional to the rate of ion formation, dependent on the concentrations of species in the sample gas.
#The method is typically the standard detection method for column separated hydrocarbons, however the method is also sensitive to almost all compounds, mostly combustible ones. 
#There are, however, a few compounds to which the method has very little, if any, sensitivity, including: O2, N2, SO2, NO, N2O, NO2, NH3, CO, CO2, and H2O.
#When used as a standalone method (without prior chromatography) it is non-specific for mixed gases in the sample.

#Typical measurement operation for hydrocarbons:
#1. The sample gas is introduced into a hydrogen flame inside the flame ionisation detector.
#2. Any hydrocarbons in the sample will produce ions when they are burnt. 
#3. Ions are detected using a metal collector which is biased with a high DC voltage. 
#4. The current across this collector is thus proportional to the rate of ionisation which in turn depends upon the concentration of hydrocarbons in the sample gas.

'flame ionisation detection':{'abbreviation':'FID', 'sampling_type':'low volume vontinuous', 'measured_parameters':['C2H6','CH4','ISOPRENE'], 'qa_accepted_parameters':[], 
    'instruments':{}
},

#------------------------------------
#selective combustion - flame ionisation detection (SC-FID)

#This method combines flame ionisation detection with selctive combustion to enable to measurement of specific hydrocarbons in a mixed gaseous sample.
#This is most commonly employed for the measurement of CH4, but can be done for othwer hydrocarbons also. 

#Typical measurement operation for hydrocarbons:
#1. The sample gas is split into 2 lines
#2. On the first line, the sample gas is introduced into a hydrogen flame inside the flame ionisation detector (reference measurement).
#2. Any hydrocarbons in the sample will produce ions when they are burnt. 
#3. Ions are detected using a metal collector which is biased with a high DC voltage. 
#4. The current across this collector is thus proportional to the rate of ionisation which in turn depends upon the concentration of hydrocarbons in the sample gas.
#5. On the second line, the sample gas is first passed though a selective catalytic combustion unit, leaving only a specific hydrocarbon in the sample (oxidising the others), for example CH4
#6. This is then measured using FID, and thus gives the concentration of the isolated hydrocarbon
#7. If the isolated hydrocarbon measured on the second line is CH4, then line1 - line2 will give the NMHC concentration

'selective combustion - flame ionisation detection':{'abbreviation':'SC-FID', 'sampling_type':'low volume vontinuous', 'measured_parameters':['C2H6','CH4','ISOPRENE'], 'qa_accepted_parameters':['C2H6','CH4','ISOPRENE'], 
    'instruments':{
        'Horiba APHA350':{'measured_parameters':['CH4'], 'qa_accepted_parameters':['CH4'],                                                                                                                                                                     'instrument_further_details':'Details about instrument found here: https://www.umweltbundesamt.de/sites/default/files/medien/publikation/long/2702.pdf'}, 
        'Horiba APHA360':{'measured_parameters':['CH4'], 'qa_accepted_parameters':['CH4'], 'documented_precision':'1.0%', 'documented_accuracy':'1.0%',                                     'documented_span_drift':'0.5%/week', 'documented_flow_rate':'0.5', 'instrument_manual_name':'Horiba_APMA360_Specs.pdf'}, 
        'Horiba APHA370':{'measured_parameters':['CH4'], 'qa_accepted_parameters':['CH4'], 'documented_precision':'2.0%', 'documented_accuracy':'1.0%', 'documented_zero_drift':'2.0%/day', 'documented_span_drift':'2.0%/day',  'documented_flow_rate':'0.9', 'instrument_manual_name':'Horiba_APHA370_Manual.pdf'} 
    }
},

#------------------------------------
#conductimetry (CD)

#Involves the absorption of a specific gas species in deionized water, to produce an acid, which can be detected by a conductivity cell.
#Through the use of the deionized water reagent, measurements are susceptible to interferences from CO2, salt aerosols, acid mists and basic gases. 
#Addition of hydrogen peroxide to the deionized water minimises interference (through reduced solubility of CO2 within the water).
#This method been used extensively to measure SO2.

#Known interferences: Any gas that can form or remove ions (NO2, NH3, HCl, Cl2 are the worst known interferants).

#Not QA accepted measurement technique for SO2 due to dominance of other more accurate techniques in measuring this species

#See details of method for SO2:https://books.google.es/books?id=DiRevS__HUkC&pg=PA6&lpg=PA6&dq=SO2+conductimetric&source=bl&ots=wTME59BaGo&sig=D8wxHHP99GihBGlBJv1DXCUxcis&hl=en&sa=X&ved=0ahUKEwjmloGg7-fbAhXDzRQKHdtIAwAQ6AEIZjAH#v=onepage&q=SO2%20conductimetric&f=false'

'conductimetry':{'abbreviation':'CD', 'sampling_type':'low volume continuous', 'measured_parameters':['NH3','SO2'], 'qa_accepted_parameters':['NH3'], 
    'instruments':{
        'Asarco 500':{                                'measured_parameters':['SO2'], 'qa_accepted_parameters':[], 'instrument_further_details':'cannot find any information. AQS states it as conductance. Assume conductimetry as detailed in general method reference below', 'measuring_assumption':True},
        'ECN AMOR':  {'sample_preparation':'Denuder', 'measured_parameters':['NH3'],                              'instrument_further_details':'See reference here:https://www.sciencedirect.com/science/article/pii/096016869390280C'}
    }
},

#------------------------------------
#second derivative spectroscopy (SDS)

#Derivative spectroscopy involves plotting the first, second or higher order derivatives of a spectrum with respect to wavelength. 
#Usually this is obtained by a microprocessor connected in series with a spectrophotometer, which computes the derivative with respect to time as the spectrum is scanned at constant speed.
#The "true" wavelength derivative is linearly related to the time derivative recorded, the magnitude of which is directly affected by scanning speed and spectral band width.
#The derivative process provides two general advantages: first, an effective enhancement of resolution, which can be useful to separate two or more components with overlapping spectra; second, a discrimination in favour of the sharpest features of a spectrum, used to eliminate interferences by broad band constituents. However, this process results in a decrease in the signal to noise ratio.
#Both advantages and disadvantages increase with the derivative order. Generally, the second derivative is more useful than the first ones.

#Not QA accepted measurement technique for NO/NO2/SO2 due to dominance of other more accurate techniques in measuring these species

#See details about method for NO,NO2,NH3,SO2:https://www.osapublishing.org/as/abstract.cfm?uri=as-53-11-1352

'second derivative spectroscopy':{'abbreviation':'SDS', 'sampling_type':'low volume continuous', 'measured_parameters':['NO','NO2','NH3','SO2'], 'qa_accepted_parameters':['NH3'], 
    'instruments':{
        'Ecotech SM1000':{'measured_parameters':['SO2'], 'qa_accepted_parameters':[], 'instrument_further_details':'No manual or specifications available but confirmation it is uses second derivative spectroscopy here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'}
    }
},

#------------------------------------
#coulometry (CM)

#Measures the current necessary to maintain a low concentration of halide and halogen at equilibrium. 
#Typical coulometric methods use solutions of potassium bromide (halide) and bromine (halogen) in dilute sulphuric acid.  
#The concentration of bromine is measured as a voltage between an indicator electrode and a bromide reference electrode.
#When the measured species is present, the bromine concentration is decreased through reaction. A feedback system then compensates by generator bromine at a generator electrode.
#The current required to regenerate the bromine is directly proportional to the concentration of the measured species.
#The method is extremely sensitive to interference through reaction with other undesired species (typically NO, NO2, O3, sulphur species, Cl2).
#Other species are typically removed through pre-filters.
#This method been used to measure NO2 and SO2.

#Known interferences: Reaction with other gases, typically: NO, NO2, O3, sulphur species, Cl2.

#Not QA accepted measurement technique for NO2 or NOx due to dominance of other more accurate techniques in measuring these species

#See details about method for NO2: https://www.tandfonline.com/doi/pdf/10.1080/00022470.1962.10468113; and for SO2: https://www.tandfonline.com/doi/pdf/10.1080/00022470.1968.10469104

'coulometry':{'abbreviation':'CM', 'sampling_type':'low volume continuous', 'measured_parameters':['NO2','NOx','SO2'], 'qa_accepted_parameters':[], 
    'instruments':{
        'Philips PW9700':{'measured_parameters':['SO2'], 'documented_lower_limit_of_detection':'4.0', 'instrument_further_details':'No manual or specifications available but confirmation it is uses coulometry here, and gives detection limit:https://books.google.es/books?id=KGQyBwAAQBAJ&pg=PT423&lpg=PT423&dq=PHILIPS+PW9700&source=bl&ots=lbloftODx4&sig=uBDyugqVzVHAhH2K21mDH7q2wYc&hl=en&sa=X&ved=0ahUKEwjxgYGZwOfbAhWHORQKHWEiBCAQ6AEIMTAB#v=onepage&q=PHILIPS%20PW9700&f=false'},
        'Philips PW9755':{'measured_parameters':['SO2'],                                              'instrument_further_details':'No manual or specifications available but confirmation it is uses coulometry here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'}
    }
},
#------------------------------------
#polarography (PO)

#An electrochemical technique where the cell current is measured as a function of time and as a function of the potential between the indicator and reference electrodes. 
#The working electrode is a dropping mercury (or other liquid conductor) electrode and unstirred solutions are used. 
#The potential is varied using pulses of linearly increasing amplitude (one pulse during each drop lifetime) and the current is sampled before and after each voltage pulse.

#Not QA accepted measurement technique for NO/NO2/SO2 due to dominance of other more accurate techniques in measuring these species

#See details about method for NO, NO2 and SO2: https://www.sciencedirect.com/science/article/pii/0039914079800041

'polarography':{'abbreviation':'PO', 'sampling_type':'injection', 'measured_parameters':['NO','NO2','SO2'], 'qa_accepted_parameters':[], 
    'instruments':{}
},

#------------------------------------
#capillary electrophoresis (CE)

#Capillary electrophoresis (CE) is a family of electrokinetic separation methods performed in submillimeter diameter capillaries and in micro and nanofluidic channels.
#In capillary electrophoresis methods, analytes migrate through electrolyte solutions under the influence of an electric field. 
#Analytes can be separated according to ionic mobility and/or partitioning into an alternate phase via non-covalent interactions. 
#Additionally, analytes may be concentrated or "focused" by means of gradients in conductivity and pH.

#Not QA accepted measurement technique for SO2 due to dominance of other more accurate techniques in measuring this species
#See details about method for SO2: https://journals.sagepub.com/doi/abs/10.1177/1082013213494900

#See details for ions here: https://www.waters.com/waters/library.htm?locale=en_AU&lid=1822645&cid=511436 and https://www.isterre.fr/english/platforms-services-1235/platform-synthesis-and-characterization/synthesis-and-solid-liquid-interactions/article/capillary-electrophoresis-waters.html

'capillary electrophoresis':{'abbreviation':'CE', 'sampling_type':'injection', 'measured_parameters':['SO2','HNO3','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'qa_accepted_parameters':['HNO3','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 
    'instruments':{}
},

#------------------------------------
#ultraviolet fluorescence (UVF)

#Fluorescence is the emission of light (luminescence) by a substance that has absorbed light or other electromagnetic radiation (excitation). In most cases, the emitted light has a longer wavelength, and therefore lower energy, than the absorbed radiation.
#Ultraviolet fluorescence specifically refers to the process of a species being excited by ultraviolet light (10nm to 400nm) (species + hv (UV) --> species*).
#The return to a fundamental electronic state of the excited species is made by luminous radiation on a longer wavelength spectrum (species* --> species + hv).
#SO2 is efficiently excited to an excited state (SO2*) by UV light between 190 nm-230 nm, subsequently fluorescing light at approximately 330nm (SO2* --> SO2 +hv(330nm)), which can be measured.

#Typical measurement operation for SO2:
#1. The sample is drawn into the sample bulkhead (typically the sample flows through a hydrocarbon “kicker,” which removes hydrocarbons from the sample by forcing the hydrocarbon molecules to permeate through the tube wall).
#2. The sample then flows into the measurement cell, where pulsating UV light from a light source (typically mercury or zinc vapour lamp) between 190nm-230nm excites the SO2 molecules by the reaction: SO2 + hv (190nm-230nm) --> SO2* . A vacuum diode, UV detector that converts UV light to a DC current is used to measure the intensity of the excitation UV source lamp.
#3. Typically a band pass filter limits the wavelength of the UV source light to approximately 214 nm, and therefore the excitation reaction becomes: SO2 + hv (214nm) --> SO2* 
#4. The amount SO2 converted to SO2* in the sample chamber is dependent on the average intensity of the UV light (Ia) and not its peak intensity because the intensity of UV light is not constant in every part of the sample chamber. Some of the photons are absorbed by the SO2 as the light travels through the sample gas. The equation for defining the average intensity of the UV light (Ia) is: Ia = I0 [1 − exp(− ax(SO2))], Where: I0 = Intensity of the excitation UV light, a = The absorption coefficient of SO2 (a constant), SO2 = Concentration of SO2 in the sample chamber, x = The distance between the UV source and the SO2 molecule(s) being affected (path length). 
#5. SO2* then fluoresces light at a longer (lower energy) wavelength centred at 330nm: SO2* --> SO2 + hv (330nm).
#6. The amount of detectable UV given off by the decay of the SO2* is affected by the rate at which this reaction occurs (k). F = k(SO2*), Where: F = the amount of fluorescent light given off, k = The rate at which the SO2* decays into SO2, SO2* = Amount of excited-state SO2 in the sample chamber. Therefore: (SO2*) -kf --> SO2 + hv (330nm)
#7. The function (k) is affected by the temperature of the gas. The warmer the gas, the faster the individual molecules decay back into their ground state and the more photons of UV light are given off per unit of time. Given that the absorption rate of SO2 (a) is constant, the amount of fluorescence (F) is a result of: a) The amount of SO2* created which is affected by the variable factors: concentration of SO2; intensity of UV light (I0); path length of the UV light(x) and b) The amount of fluorescent light created which is affected by the variable factors: the amount of SO2* present and the rate of decay (k) which changes based on the temperature of the gas.
#8. When and the intensity of the light (I0) is known; path length of excited light is short (x); the temperature of the gas is known and compensated for so that the rate of SO2*decay is constant (k). and; no interfering conditions are present (such as interfering gases or stray light); the amount of fluorescent light emitted (F) is directly related to the concentration of the SO2 in the Sample Chamber.
#9. A Photo Multiplier Tube (PMT) detects the UV given off by the SO2* decay (330 nm) and outputs an analog signal. Several focusing lenses and optical filters ensure that both the PMT and source lamp vacuum diode are exposed to an optimum amount of only the right wavelengths of UV. To further assure that the PMT only detects light given off by decaying SO2* the pathway of the excitation UV and field of view of the PMT are perpendicular to each other and the inside surfaces of the sample chamber are coated with a layer of black teflon that absorbs stray light.
#10. The net result of the careful instrumental precision is any variation in UV fluorescence can be directly attributed to changes in the concentration of SO2 in the sample gas.

#Known interferences: 3rd body quenching (NO, CO2, O2, H2O), light pollution, UV absorption by ozone, other species undergoing ultraviolet fluorescence (poly-nuclear aromatics, NO)

'ultraviolet fluorescence':{'abbreviation':'UVF', 'sampling_type':'low volume continuous', 'measured_parameters':['SO2'], 'qa_accepted_parameters':['SO2'], 
    'instruments':{
        'Airpointer-SO2':        {'documented_lower_limit_of_detection':'0.5',   'documented_upper_limit_of_detection':'10000.0',   'documented_precision':'1.0/1.0%',    'documented_accuracy':'1.0%>=100.0',                  'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%>=100/day',                                  'documented_flow_rate':'0.5',       'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here: https://www.airpointer.com/#sc-tabs-1547659231751|sc-tabs-1547659231301'},
        'Beckman 953':           {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'},  
        'Dasibi 4108':           {'documented_lower_limit_of_detection':'1.0',   'documented_upper_limit_of_detection':'2000.0',    'documented_precision':'1.0',                                                                                                                                                                                                                   'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here, and details about limit of detection:https://books.google.es/books?id=ZNvqCAAAQBAJ&pg=PA155&lpg=PA155&dq=Dasibi+4108&source=bl&ots=rxhgzzOLrY&sig=W2zmIMXsFjAcNUXGLgTCYViACR0&hl=en&sa=X&ved=0ahUKEwj8xOXN2ufbAhXDchQKHRVWBksQ6AEIVjAI#v=onepage&q=Dasibi%204108&f=false'},
        'DKK-TOA GFS-32':        {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:http://itepsrv1.itep.nau.edu/itep_course_downloads/DAI%20resources/old%20resource%20files/Reference%20and%20Equivalent%20Methods%201-2007.pdf'},
        'DKK-TOA GFS-112E':      {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://www.federalregister.gov/documents/2000/01/18/00-1083/ambient-air-monitoring-reference-and-equivalent-methods-designation-of-a-new-equivalent-method-for'},
        'DKK-TOA GFS-312E':      {'documented_lower_limit_of_detection':'0.2',   'documented_upper_limit_of_detection':'1000.0',                                                                                                'documented_zero_drift':'0.5/day',    'documented_span_drift':'1.0%/day',                                                                           'instrument_manual_name':'DKK-TOA_GFS-312E_Specs.pdf'},
        'Ecotech AM2020':        {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'},  
        'Ecotech EC9850':        {'documented_lower_limit_of_detection':'0.5',   'documented_upper_limit_of_detection':'20000.0',   'documented_precision':'0.5/1.0%',    'documented_accuracy':'1.0%',                         'documented_zero_drift':'1.0/day',    'documented_span_drift':'0.5%/day',      'documented_resolution':'0.001', 'documented_flow_rate':'0.5',       'instrument_manual_name':'Ecotech_EC9850_Manual.pdf'},
        'Ecotech EC9850T':       {'documented_lower_limit_of_detection':'0.2',   'documented_upper_limit_of_detection':'2000.0',    'documented_precision':'0.5/0.5%',                                                          'documented_zero_drift':'0.2/day',    'documented_span_drift':'0.5%/day',                                       'documented_flow_rate':'0.5',       'instrument_manual_name':'Ecotech_EC9850T_Specs.pdf'},
        'Ecotech Serinus50':     {'documented_lower_limit_of_detection':'0.3',   'documented_upper_limit_of_detection':'20000.0',   'documented_precision':'0.5/0.15%',   'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.5/day',    'documented_span_drift':'1.0%/day',                                       'documented_flow_rate':'0.75',      'instrument_manual_name':'Ecotech_Serinus50_Manual.pdf'},
        'Environnement SA AF20M':{                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://books.google.es/books?id=sZLqCAAAQBAJ&pg=PA73&lpg=PA73&dq=environnement+af20m&source=bl&ots=OOOfZQbtkW&sig=ACfU3U0bdsNJH5j2XhFzrkF2tnkb6e2HKQ&hl=en&sa=X&ved=2ahUKEwjAidHQtf_fAhXD8uAKHcYOBRYQ6AEwAnoECAgQAQ#v=onepage&q=environnement%20af20m&f=false'},
        'Environnement SA AF21M':{                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2586140/'},
        'Environnement SA AF22M':{'documented_lower_limit_of_detection':'0.4',   'documented_upper_limit_of_detection':'10000.0',                                         'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5%/day',                                       'documented_flow_rate':'0.3',       'instrument_manual_name':'Environnement_SA_AF22M_Specs.pdf'},
        'Horiba APSA350':        {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://books.google.es/books?id=qmWXBwAAQBAJ&pg=PA88&lpg=PA88&dq=apsa350&source=bl&ots=dBQKXVwgi2&sig=ACfU3U033XGNb0kOfeN06k0HDr7Oot-Imw&hl=en&sa=X&ved=2ahUKEwi0jbGJsP_fAhWC2OAKHTYyCU44ChDoATAAegQIABAB#v=onepage&q=apsa350&f=false'},
        'Horiba APSA360':        {'documented_lower_limit_of_detection':'0.5',   'documented_upper_limit_of_detection':'10000.0',   'documented_precision':'1.0%',        'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5/day',                                        'documented_flow_rate':'0.8',       'instrument_manual_name':'Horiba_APMA360_Specs.pdf'},
        'Horiba APSA370':        {'documented_lower_limit_of_detection':'0.5',   'documented_upper_limit_of_detection':'10000.0',   'documented_precision':'1.0%',        'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5/day',                                        'documented_flow_rate':'0.7',       'instrument_manual_name':'Horiba_APSA370_Specs.pdf'},
        'MCV 64FUV':             {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'cannot find any information. EIONET states it as automatic for SO2 which is almost always ultraviolet fluorescencer. Assume this.', 'measuring_assumption':True},  
        'Meloy SA700':           {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'},  
        'Philips K50033':        {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&ved=2ahUKEwiaspe-19jnAhWI4IUKHVswA1oQFjADegQICBAB&url=http%3A%2F%2Ffbd.beun.edu.tr%2Findex.php%2Fzkufbd%2Farticle%2Fdownload%2F53%2F100&usg=AOvVaw0l7fQAqr1Pealzx0MaZ6d0'},
        'Seres SF2000G':         {'documented_lower_limit_of_detection':'1.0',   'documented_upper_limit_of_detection':'20000.0',                                         'documented_accuracy':'1.0%',                         'documented_zero_drift':'2.0/week',   'documented_span_drift':'1.0%/month',                                     'documented_flow_rate':'0.5-0.833', 'instrument_manual_name':'Seres_SF2000G_Specs.pdf'},
        'SIR S-5001':            {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:http://www.eufar.net/instruments/423'},
        'Teledyne 100':          {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'cannot find any information. AQS states it as ultraviolet fluorescence. Assume this', 'measuring_assumption':True},
        'Teledyne 100A':         {'documented_lower_limit_of_detection':'0.4',   'documented_upper_limit_of_detection':'20000.0',   'documented_precision':'0.5%',        'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5%/7days',                                     'documented_flow_rate':'0.65',      'instrument_manual_name':'Teledyne_100A_Manual.pdf'},
        'Teledyne 100AH':        {'documented_lower_limit_of_detection':'100.0', 'documented_upper_limit_of_detection':'5000000.0', 'documented_precision':'0.5%',        'documented_accuracy':'1.0%',                         'documented_zero_drift':'1000.0/day', 'documented_span_drift':'0.5%/7days',                                     'documented_flow_rate':'0.65',      'instrument_manual_name':'Teledyne_100AH_Manual.pdf'},
        'Teledyne 100E':         {'documented_lower_limit_of_detection':'0.4',   'documented_upper_limit_of_detection':'20000.0',   'documented_precision':'0.5%>=50.0',  'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5%/day',                                       'documented_flow_rate':'0.65',      'instrument_manual_name':'Teledyne_100E_Manual.pdf'},
        'Teledyne 100EH':        {'documented_lower_limit_of_detection':'100.0', 'documented_upper_limit_of_detection':'5000000.0', 'documented_precision':'0.5%',        'documented_accuracy':'1.0%',                         'documented_zero_drift':'1000.0/day', 'documented_span_drift':'0.5%/7days',                                     'documented_flow_rate':'0.7',       'instrument_manual_name':'Teledyne_100EH_Addendum.pdf'},
        'Teledyne 100EU':        {'documented_lower_limit_of_detection':'0.05',  'documented_upper_limit_of_detection':'20000.0',   'documented_precision':'0.5%',        'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.2/day',    'documented_span_drift':'0.5%/7days',                                     'documented_flow_rate':'0.65',      'instrument_manual_name':'Teledyne_100EU_Addendum.pdf'},
        'Teledyne 101A':         {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:http://pagina.jccm.es/medioambiente/rvca/pdfs/INFORME_UNIDAD-MOVIL-MORA_2_15%20.pdf'},
        'Teledyne T100':         {'documented_lower_limit_of_detection':'0.4',   'documented_upper_limit_of_detection':'20000.0',   'documented_precision':'0.5%>=50.0',  'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5%/day',                                       'documented_flow_rate':'0.65',      'instrument_manual_name':'Teledyne_T100_Manual.pdf'},
        'Teledyne T100H':        {'documented_lower_limit_of_detection':'200.0', 'documented_upper_limit_of_detection':'5000000.0', 'documented_precision':'0.5%>=10000', 'documented_accuracy':'1.0%',                         'documented_zero_drift':'1000.0/day', 'documented_span_drift':'0.5%/day',                                       'documented_flow_rate':'0.7',       'instrument_manual_name':'Teledyne_T100H_Addendum.pdf'},
        'Teledyne T100U':        {'documented_lower_limit_of_detection':'0.05',  'documented_upper_limit_of_detection':'20000.0',   'documented_precision':'0.5%',        'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.2/day',    'documented_span_drift':'0.5%/day',                                       'documented_flow_rate':'0.65',      'instrument_manual_name':'Teledyne_T100U_Addendum.pdf'},
        'Teledyne T101':         {'documented_lower_limit_of_detection':'0.4',   'documented_upper_limit_of_detection':'20000.0',   'documented_precision':'0.5%>=50.0',  'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.5/day',    'documented_span_drift':'0.5%/day',                                       'documented_flow_rate':'0.65',      'instrument_manual_name':'Teledyne_T101_Manual.pdf'},
        'Teledyne ML8850':       {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://books.google.es/books?id=BODqCAAAQBAJ&pg=PA137&lpg=PA137&dq=MONITOR+LABS+8850&source=bl&ots=TUXOyv5p5D&sig=Ult9_8r_m8xwlzmZIHNu3a3JaMs&hl=en&sa=X&ved=0ahUKEwjsu52E2ufbAhXBvhQKHVNNCrkQ6AEIdzAM#v=onepage&q=MONITOR%20LABS%208850&f=false'},
        'Teledyne ML8850S':      {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://www.der.wa.gov.au/images/documents/your-environment/air/publications/perth-photochemical-smog-study-1996.pdf'},
        'Teledyne ML9850':       {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:http://www.epa.ie/licences/lic_eDMS/090151b28005fa5c.pdf'},
        'Teledyne ML9850B':      {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:http://www.epa.ie/licences/lic_eDMS/090151b28005fa5c.pdf'},
        'Thermo 43':             {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA249&lpg=PA249&dq=monitor+labs+8450&source=bl&ots=-uIdNzc1_l&sig=LqRK0OiodCMg5Sow1wlXQgIVLx0&hl=en&sa=X&ved=0ahUKEwjCy8eNvOfbAhVBuhQKHZTJBUoQ6AEIWDAG#v=onepage&q=monitor%20labs%208450&f=false'},  
        'Thermo 43A':            {'documented_lower_limit_of_detection':'1.0',   'documented_upper_limit_of_detection':'100.0',                                                                                                                                                                                                                                                     'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here, and details about limit of detection:https://books.google.es/books?id=G0EYUql3_98C&pg=PA283&lpg=PA283&dq=thermo+43a&source=bl&ots=k2RUJXlAVJ&sig=yIxG08AoSyuqedtS6dp-O7O04X4&hl=en&sa=X&ved=0ahUKEwj17qqB3OfbAhVCVRQKHanFDlo4ChDoAQgoMAA#v=onepage&q=thermo%2043a&f=false'},
        'Thermo 43B':            {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here:http://www.southwestohioair.org/UserFiles/Servers/Server_3788196/File/EnvironmentalServices/AirQuality/Monitoring/AirMonitoringNetwork16-17.pdf'},
        'Thermo 43C':            {'documented_lower_limit_of_detection':'0.5',   'documented_upper_limit_of_detection':'100000.0',  'documented_precision':'1.0/1.0%',    'documented_accuracy':'1.0%',                         'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/day',                                       'documented_flow_rate':'0.5',       'instrument_manual_name':'Thermo_43C_Specs.pdf'},
        'Thermo 43C-TLE':        {'documented_lower_limit_of_detection':'0.05',  'documented_upper_limit_of_detection':'1000.0',    'documented_precision':'0.2/1.0%',                                                          'documented_zero_drift':'0.2/day',    'documented_span_drift':'1.0%/day',                                       'documented_flow_rate':'0.5',       'instrument_manual_name':'Thermo_43C-TLE_Specs.pdf'},
        'Thermo 43H':            {'documented_lower_limit_of_detection':'10.0',  'documented_upper_limit_of_detection':'20000.0',                                                                                                                                                                                                                                                   'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here, and details about detection limits: https://www.yumpu.com/fr/document/read/7309580/catalogue-cleanair-europe/36'},
        'Thermo 43i':            {'documented_lower_limit_of_detection':'0.5',   'documented_upper_limit_of_detection':'100000.0',                                        'documented_accuracy':'1.0%',                         'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/day',                                       'documented_flow_rate':'0.5',       'instrument_manual_name':'Thermo_43i_Manual.pdf'},
        'Thermo 43i-TLE':        {'documented_lower_limit_of_detection':'0.05',  'documented_upper_limit_of_detection':'1000.0',                                          'documented_accuracy':'1.0%',                         'documented_zero_drift':'0.2/day',    'documented_span_drift':'1.0%/day',                                       'documented_flow_rate':'0.5',       'instrument_manual_name':'Thermo_43i-TLE_Manual.pdf'},
        'Thermo 43S':            {                                                                                                                                                                                                                                                                                                                                                  'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here: https://books.google.es/books?id=JqpwDwAAQBAJ&pg=PA232&lpg=PA232&dq=thermo+43s&source=bl&ots=h9MazeIX98&sig=ACfU3U0F3I9zCn-qI-pM7Np-yN5w7AxRGg&hl=en&sa=X&ved=2ahUKEwiwh9ri1f7fAhXh8eAKHfd9CxYQ6AEwCXoECAkQAQ#v=onepage&q=thermo%2043s&f=false'},
        'Thermo 43W':            {'documented_lower_limit_of_detection':'1.0',   'documented_upper_limit_of_detection':'750.0',                                                                                                                                                                                                                                                     'instrument_further_details':'No manual or specifications available but confirmation it is uses ultraviolet fluorescence here, and details about limit of detection:https://books.google.es/books?id=G0EYUql3_98C&pg=PA283&lpg=PA283&dq=thermo+43a&source=bl&ots=k2RUJXlAVJ&sig=yIxG08AoSyuqedtS6dp-O7O04X4&hl=en&sa=X&ved=0ahUKEwj17qqB3OfbAhVCVRQKHanFDlo4ChDoAQgoMAA#v=onepage&q=thermo%2043a&f=false'},
        'Thermo 450i':           {'documented_lower_limit_of_detection':'1.5',   'documented_upper_limit_of_detection':'100000.0',                                        'documented_accuracy':'1.0%<100000.0 5.0%>=100000.0', 'documented_zero_drift':'1.0/day',    'documented_span_drift':'1.0%/day',                                       'documented_flow_rate':'1.0',       'instrument_manual_name':'Thermo_450i_Manual.pdf'}
    }
},

#------------------------------------
#thermal reduction - ultraviolet fluorescence (TR-UVF)

#Indirect method for the measurement of sulphate (SO4--) through the conversion of SO4-- to SO2 by thermal reduction, which is then measured by ultraviolet fluorescence.
#The sulphate is converted by drawing a continuous stream of sample across a hot reactive surface that reduces sulphate particles in the sample stream to sulfur dioxide gas. 
#The concentrations of sulphate in the ambient air are then quantified by comparing the signal produced when aerosol-laden sample is drawn directly into the converter to a background signal that is produced when the sample stream is run through a high-efficiency particulate aerosol filter that removes the sulfate before conversion. 
#The difference in signal between the filtered and unfiltered sample can be attributed to sulfur dioxide that is formed from sulfate particles in the unfiltered sample stream. 
#By routinely switching between the filtered and unfiltered sample streams, the instrument readings can be continuously adjusted or corrected for changes in background signal that might be produced by traces of SO2 or other interfering gases. 
#This frequent adjustment for changes in background signal improves the system stability, which in turn improves the limit of detection.


'thermal reduction - ultraviolet fluorescence':{'abbreviation':'TR-UVF', 'sampling_type':'low volume continuous', 'measured_parameters':['SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'qa_accepted_parameters':['SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 
    'instruments':{
        'Thermo 5020i':{'documented_lower_limit_of_detection':'0.5', 'documented_upper_limit_of_detection':'4000.0', 'documented_accuracy':'1.0%', 'documented_zero_drift':'1.0/day', 'documented_span_drift':'1.0%/day', 'documented_flow_rate':'0.4-0.5', 'instrument_manual_name':'Thermo_5020i_Manual.pdf'}
    }
},

#------------------------------------
#laser-induced fluorescence (LIF)

#Fluorescence is the emission of light (luminescence) by a substance that has absorbed light or other electromagnetic radiation (excitation). In most cases, the emitted light has a longer wavelength, and therefore lower energy, than the absorbed radiation.
#Laser-induced fluorescence (LIF) is a method where fluorescence is induced by an atom or molecule being excited  by the absorption of laser light.
#The emission light is detected using a photomutiplier tube. The amount of a species present in the sample is proportional to the intensity of the emission light of a specific wavelength and quantification is possible by using calibration data.
#The LIF method has been used mostly in the measurement of hydroxyl and hydroperoxyl radicals, but is also used for measurement of NO/NO2.

#Known Interferences: Laser power modulation (laser excitation could generate OH radicals biasing the measurement). 3rd body quenching.

#See details of method for NO/NO2 here:https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/1999GL900015 and here:http://etheses.whiterose.ac.uk/11668/1/Amy%20Foulds_Final%20MSc%20Thesis.pdf

'laser-induced fluorescence':{'abbreviation':'LIF', 'sampling_type':'low volume continuous', 'measured_parameters':['NO','NO2','NOx'], 'qa_accepted_parameters':['NO','NO2','NOx'], 
    'instruments':{}
},

#------------------------------------
#vacuum ultraviolet resonance fluorescence (VURF)

#Fluorescence is the emission of light (luminescence) by a substance that has absorbed light or other electromagnetic radiation (excitation). In most cases, the emitted light has a longer wavelength, and therefore lower energy, than the absorbed radiation.
#Vacuum Ultraviolet Resonance Fluorescence is a method where a reasonance lamp excited by reasonance fluorescence discharge (in combiantion with an optical filter) produces photons in the ultraviolet which react with a sample gas, inducing fluorescence in the vacuum ultraviolet (UV radiaion 10nm-200nm), subsequently detected by a photomultiplier tube.
#This method can be employed for CO, with the reasonance lamp emitting UV light between 145nm-151nm, with fluorescence occurring between 160nm-190nm.

#Known interferences: Water vapour, drifts in lamp intensity, continuum raman scattering by O2. 

'vacuum ultraviolet resonance fluorescence':{'abbreviation':'VURF', 'sampling_type':'low volume continuous', 'measured_parameters':['CO'], 'qa_accepted_parameters':['CO'], 
    'instruments':{
        'Aerolaser AL5001':{'documented_lower_limit_of_detection':'0.8', 'documented_upper_limit_of_detection':'100000.0', 'measuring_instrument_manual_name':'Aerolaser_AL5001_Specs.pdf'},
        'Aerolaser AL5002':{'documented_lower_limit_of_detection':'0.8', 'documented_upper_limit_of_detection':'100000.0', 'measuring_instrument_manual_name':'Aerolaser_AL5002_Specs.pdf'}
    }
},

#------------------------------------
#cavity ringdown spectroscopy (CRDS)

#Based on absorption spectroscopy, Cavity Ringdown Spectroscopy works by attuning light rays to the unique molecular fingerprint of the sample species. 
#By measuring the time it takes the light to fade or "ring-down", you receive an accurate molecular count in milliseconds. 
#The time of light decay, in essence, provides an exact, non-invasive, and rapid means to detect contaminants in the air, in gases, and even in the breath.
#The method is typically employed to measure CO and other greenhouse gases (e.g CO2, CH4, H2O).

#Typical measurement operation:
#1. A Continuous Wave (CW) diode laser emits a directed beam of light energy through an ultra-high reflective mirror into the absorption cell (cavity).
#2. The light reflects back and forth between two ultra-high reflective mirrors multiple times, up to a total path length of 100 kilometers.
#3. Once the photodiode detector “sees” a preset level of light energy, the light source is shuttered or diverted from the cavity.
#4. On each successive pass, a small amount of light or ring-down signal emits through the second mirror and is sensed by the light detector.
#5. Once the light "rings down", the detector achieves a point of zero light energy in milliseconds, and the measurement is complete.    
#6. The computer-controlled system tunes the laser off the absorption peak for the sample species to determine the τ empty value, equivalent to a zero baseline correction. It tunes back to the absorption peak to determine the τ value, dependent on the sample species concentration.
#7. The concentration of the sample species is directly calculated using Beer’s Law. The measured value constitutes an absolute measurement and is unaffected by losses outside the ring-down cavity.

#Known interferences: Water vapour, CO2 and particulates.

'cavity ringdown spectroscopy':{'abbreviation':'CRDS', 'sampling_type':'low volume continuous', 'measured_parameters':['CO','CO2','CH4','H20'], 'qa_accepted_parameters':['CO','CO2','CH4','H20'], 
    'instruments':{
        'LGR Off-Axis ICOS':{'measured_parameters':['CO'],              'qa_accepted_parameters':['CO'],              'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':'4000.0',                                                                                                                                    'documented_precision':'0.05',                                             'documented_span_drift':'1.0/day',                                                                                               'instrument_manual_name':'LGR_Off-Axis_ICOS_Specs.pdf'},
        'Picarro G1301':    {'measured_parameters':['CO2','CH4','H2O'], 'qa_accepted_parameters':['CO2','CH4','H2O'],                                                                                                                                                                                                                                                                                                                                                                                                               'documented_flow_rate':'0.4',  'instrument_manual_name':'Picarro_G1301_Manual.pdf'},
        'Picarro G2301':    {'measured_parameters':['CO2','CH4','H2O'], 'qa_accepted_parameters':['CO2','CH4','H2O'],                                               'documented_upper_limit_of_detection':{'CO2':'1000000.0','CH4':'20000.0','H2O':np.NaN},                                                                                            'documented_precision':{'CO2':'25.0','CH4':'0.22','H2O':'30000.0'},        'documented_span_drift':{'CO2':'120.0/day','CH4':'1.0/day','H2O':'100000.0'},                     'documented_flow_rate':'0.4',  'instrument_manual_name':'Picarro_G2301_Specs.pdf'},
        'Picarro G2302':    {'measured_parameters':['CO','CO2','H2O'],  'qa_accepted_parameters':['CO','CO2','H2O'],                                                'documented_upper_limit_of_detection':{'CO':'5000.0','CO2':'1000000.0','H2O':np.NaN},                 'documented_uncertainty':{'CO':'2.0','CO2':np.NaN,'H2O':np.NaN},             'documented_precision':{'CO':'2.0','CO2':'50.0','H2O':'50000.0'},          'documented_span_drift':{'CO':'15.0/day','CO2':'150.0/day','H2O':'100000.0'},                     'documented_flow_rate':'0.4',  'instrument_manual_name':'Picarro_G2302_Specs.pdf'},
        'Picarro G2303':    {'measured_parameters':['CO'],              'qa_accepted_parameters':['CO'],                                                                                                                                                                                                                                                                                                                                                                                                                                                           'instrument_further_details':'EIONET reports instrumnent as a CRDS instrument but cannot find other documentation, assume it is a CRDS instrument as other Picarro instruments are: http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/picarroG2303/view', 'measuring_assumption':True},
        'Picarro G2401':    {                                                                                                                                       'documented_upper_limit_of_detection':{'CO':'5000.0','CO2':'1000000.0','CH4':'20000.0','H2O':np.NaN}, 'documented_uncertainty':{'CO':'2.0','CO2':'50.0','CH4':'1.0','H2O':np.NaN}, 'documented_precision':{'CO':'1.0','CO2':'10.0','CH4':'0.3','H2O':np.NaN}, 'documented_span_drift':{'CO':'10.0/day','CO2':'100.0/day','CH4':'1.0/day','H2O':'100000.0/day'}, 'documented_flow_rate':'0.4',  'instrument_manual_name':'Picarro_G2401_Specs.pdf'}
    }
},
#------------------------------------
#cavity attenuated phase shift spectroscopy (CAPS)

#Operates on the principle that a specific species efficiently absorbs light at a known wavelength. This is the case for NO2, at 450nm.
#The degree to which the light is absorbed by a specific species is directly related to the species concentration as described by the Beer-Lambert Law (A = εLC; A = Absorbance (mol litres-1), ε = Molar absorptivity (litres mol-1 cm-1), L = mean optical path length of cell, C = species concentration).
#Emitted light in an optical cell is reflected back and forth between two mirrors, building intensity and running a very long path length. The long path length extends the “time” or “life” of the photon, thus providing ample time to measure absorbance when a species is present.
#The method is typically employed for direct measurement of NO2.

#Typical measurement operation for NO2:
#1. The sample is drawn into the sample bulkhead, before being passed into the measurement cell. 
#2. Light is emitted from a blue ultraviolet (UV) light emitting diode (LED) centered at 450 nm (a prominent absorption band for NO2). 
#3. The measurement cell contains high reflectivity mirrors located at either end to provide an extensive optical path length. The optical cell resides in a temperature controlled oven. The oven raises the ambient temperature of the sample gas to 45 degrees Celsius. This mitigates the formation of moisture on the surfaces of the mirrors while also minimiSing changes in the absorption coefficient due to temperature fluctuations.
#4. The CAPS method is unique in that it applies the fundamental optical absorption law in the frequency domain, rather than using relative changes in light intensity as the primary signal. 
#5. UV light from the modulating high intensity LED enters a near confocal optical cell through the rear of mirror A. The intensity of the light, as observed by a vacuum photodiode detector, which is also modulating at a slightly different frequency, located behind Mirror B, builds exponentially in the cell while the LED is ON. The opposite is true when the LED is OFF. Because both mirrors are highly reflective at 450 nm, a prominent absorption band for NO2, the light takes a considerable amount of time to plateau in the absence of the absorbing gas. 
#6. When NO2 is present, the mean path length traveled by the light is significantly reduced. This has two effects on the observed intensity as measured by the detector: a) The light plateau intensity level is lower, b) The light intensity plateaus sooner. Thus, an observed phase shift from the modulating LED is detected. The phase shift is largest when measuring zero air and decreases when NO2 is present.
#7. Both the LED and the Detector are modulated ON and OFF such that the observed signal has a much lower frequency, equal to the difference between the modulated frequencies and is referred to as a beat frequency. The system hardware and software take advantage of this, as it makes it easier to post process the signal using a micro controller. The technique is known as heterodyning.
#8. The instrument translates the phase shift from the presence of absorbing gas into a concentration measurement. Typical absorption techniques of other analysers take a reference and measure value of the light intensity “level” in order to derive concentration and compensate for source drift. Using the CAPS technique the amount of phase shift remains constant for a given concentration, even if the LED drifts over time. The measurement approach offers many advantages over traditional chemiluminescence analysers, such as faster response (single gas stream), lower noise at span and more importantly greater specificity.

#Known interferences: direct spectral interference with photochemically produced 1,2-dicarbonyl species (e.g., glyoxal, methylglyoxal)

'cavity attenuated phase shift spectroscopy':{'abbreviation':'CAPS', 'sampling_type':'low volume continuous', 'measured_parameters':['NO2'], 'qa_accepted_parameters':['NO2'], 
    'instruments':{
        'Environnement SA AS32M':{'documented_lower_limit_of_detection':'0.1',  'documented_upper_limit_of_detection':'1000.0', 'documented_precision':'0.05',    'documented_accuracy':'4.0%', 'documented_zero_drift':'0.75/day', 'documented_span_drift':'1.5/day',                                'instrument_manual_name':'EnvironnementSA_AS32M_Specs.pdf'},
        'Teledyne T500U':        {'documented_lower_limit_of_detection':'0.04', 'documented_upper_limit_of_detection':'1000.0', 'documented_precision':'0.5%>=5', 'documented_accuracy':'1.0%', 'documented_zero_drift':'0.1/day',  'documented_span_drift':'0.5%/day', 'documented_flow_rate':'0.9', 'instrument_manual_name':'Thermo_T500U_Manual.pdf'}
    }
},

#------------------------------------
#differential optical absorption spectroscopy (DOAS)

#The basic principle used in Differential optical absorption spectroscopy (DOAS) is absorption spectroscopy. DOAS allows the quantitative determination of multiple atmospheric trace gas concentrations by recording and evaluating the characteristic absorption structures (lines or bands) of the trace gas molecules along an absorption path of known length in the open atmosphere, following the Beer-Lambert law (I = Io e−KLC; K = molecular absorption coefficient at STP, L = optical path length of cell, C = species concentration , I = light intensity of sample gas, Io = light intensity of sample without measured species (reference gas) ).
#The wavelength of light where a distinct absorption peak occurs is determined for analyte. A wavelength on either side of the absorption peak is next determined. The intensity of a light source at wavelength is measured and then the intensity is measured again after the light passes through the analyte. The difference of the intensities is proportional to the concentration to the analyte.
#DOAS is a long path measuring technique. Measurements can be made in an optical pathway from 1 to 10 kilometers.
#The method measures the average concentration of a species along the path length, not for any single molecule.
#In the real atmosphere, multiple effects contribute to the overall attenuation of the light. In particular, aerosols and clouds scatter light and thereby reduce the intensity of the direct beam while increasing intensities measured in other directions. Also, there rarely is only one absorber relevant at a given wavelength and this also needs to be accounted for. 
#The solution to this problem lies in the use of measurements at several wavelengths. Each molecule has a characteristic absorption spectrum (its spectral fingerprint) and therefore, simultaneous measurements at different wavelengths enable the separation of the contributions of the different absorbers. This is what DOAS does.
#Scattering by aerosols also needs accounting for. Their extinction cross-sections can be approximated by power laws (λ**-4 for Rayleigh scattering and λ**-1..0 for Mie scattering), the coefficients of which can be determined in the fit. The basic principle behind the separation of aerosol extinction and trace gas absorption is that the latter are identified using those parts of their absorption cross-sections that vary rapidly with wavelength. 
#The more slowly varying parts of the absorption can not be separated from extinction by aerosols which is why the absorption cross-sections are often high pass filtered before use in a DOAS retrieval.
#In addition to their effect on the spectral distribution of the intensity measured, aerosols can also have a significant impact on the light path of scattered light. This has to be modelled explicitly when computing the airmass factors for the light path correction.
#The DOAS technique is characterised by the following: 
#(a) measuring the transmitted light intensity over a relatively (compared to the width of an absorption band) broad spectral interval; 
#(b) high‐pass filtering of the spectra to obtain a differential absorption signal and eliminating broad‐band extinction processes such as Rayleigh and Mie scattering (RS and MS); 
#(c) quantitative determination of trace column densities by matching the observed spectral signatures to prerecorded (reference) spectra by, for instance, least‐squares methods.
#DOAS instruments are often divided into two main groups: passive and active ones. The active DOAS system such as longpath(LP)-systems and cavity-enhanced(CE) DOAS systems have their own light-source, whereas passive ones use the sun as their light source, e.g. MAX(Multi-axial)-DOAS. Also the moon can be used for night-time DOAS measurements, but here usually direct light measurements need to be done instead of scattered light measurements as it is the case for passive DOAS systems such as the MAX-DOAS.

#Typical measurement operation:
#1. The positioning of the emitter and the receiver define the monitoring path. The light source in the emitter is a high-pressure xenon lamp. This type of light source radiates an almost smooth spectrum ranging from approximately 200 nm up to 500 nm, and from about 1200 to 3000 nm. Within these ranges a number of gaseous substances show specific absorption spectra. The lamp spectrum is however not perfectly smooth, but the remaining "hills" in the spectral output are being taken care of in the evaluation.
#2. The emitted light beam is directed towards the receiver, and on its way the intensity is affected by scattering and absorption in molecules and particles.
#3. To ascertain the initial intensity (Io) of the light source across the path, there are typically two approaches: (a) Use of measurements at different light paths. If it is possible to take two measurements with the same initial intensity but at different light path lengths L, the dependence on the initial intensity cancels if one looks at the ratio of the two measurements. This approach is e.g. used in zenith-sky DOAS, where a measurement at low sun (long light path through the atmosphere) is analysed with a measurement taken at high sun (short light path). (b) Use of measurements at different wavelengths. In this approach, one exploits the fact that many molecules have structured spectra and light at different wavelengths experiences different absorption strengths. If the initial intensity does not vary with wavelength, one can look at the ratio of two measurements at different wavelengths to remove the dependence on Io. This approach is often used in long-path measurements which employ a lamp. It also is one reason to call this method differential as changes in absorption are used.
#4. From the receiver the captured light is led via an opto-fibre to the analyser. The function of the fibre is only to avoid exposing the opto-analyser to dust, high humidity, temperature variations, etc.
#5. When the light reaches the analyser, it enters a spectrometer. Inside the spectrometer, a grating refracts the light into its wavelength components. The refracted light is then projected onto a rapid scanning slit in front of a photo-multiplier detector or an infrared sensitive diode, where a selected part of the spectrum is detected. The scanning slit device makes it possible to record all wavelengths separately.
#6. As the grating is moveable, any chosen part of the spectrum can be detected. The wavelength window can thus be optimised for a certain component, with respect to parameters such as sensitivity and interfering pollutants. Approximately 100 scans per second are recorded.
#7. The current from the detector is converted into digital signals by a 12 bit analogue-to-digital converter, and the signal is stored and accumulated in a multi-channel register. The detected spectrum is typically 40 nm wide in the UV range and approximately twice as wide in the IR range. Each scan is digitised into 1000 points.
#8. Each pollutant is monitored during a time period entered by the operator. When the data accumulation is finished, the evaluation process is started. At the same time the next data accumulation period starts.

#Known interferences: Interference by miscellaneous atmospheric constituents. Heavy rain and fog, and even high humidity. Atmospheric turbulence, such as that from thermal-induced effects, can distort reflections. Anything that interrupts the path of the laser will cause some interference (i.e., animals, cars, planes, etc.). 

'differential optical absorption spectroscopy':{'abbreviation':'DOAS', 'sampling_type':'remote', 'measured_parameters':['CO','NO','NO2','NOx','O3','SO2','NH3','C6H6'], 'qa_accepted_parameters':['CO','NO','NO2','NOx','O3','SO2','NH3','C6H6'], 
    'instruments':{
        'Environnement SA SANOA':{'measured_parameters':['NO','NO2','NOx','O3','SO2','NH3'], 'qa_accepted_parameters':['NO','NO2','NOx','O3','SO2','NH3'], 'documented_lower_limit_of_detection':{'NO':'1.5','NO2':'0.6','NOx':'1.5','O3':'0.6','SO2':'0.2','NH3':'3.5'}, 'documented_accuracy':'1.0%', 'instrument_manual_name':'EnvironnementSA_SANOA_Specs.pdf'},
        'OPSIS AR500':           {                                                                                                                                                                                                                                                                      'instrument_manual_name':'OPSIS_AR500_Manual.pdf'},
        'OPSIS AR550':           {                                                                                                                                                                                                                                                                      'instrument_manual_name':'OPSIS_AR500_Manual.pdf'},
        'OPSIS AR600':           {                                                                                                                                                                                                                                                                      'instrument_manual_name':'OPSIS_AR500_Manual.pdf'},
        'OPSIS AR650':           {                                                                                                                                                                                                                                                                      'instrument_manual_name':'OPSIS_AR500_Manual.pdf'}
    }
},
#------------------------------------
#electrochemical membrane diffusion (EMD)

#Methodology developed specifically for measurement of NH3: https://www.sciencedirect.com/science/article/pii/S0003267003012650?via%3Dihub
#Gas is sampled in a sampler comprising two opposite channels separated by a gas permeable, water repellent polypropylene membrane. 
#Subsequently, the acid sample solution is pumped into a selector where an alkaline solution is added to ionize all sampled ambient acid gasses, resulting in an enhanced selectivity. 
#In the selector, the ammonia can diffuse through a second membrane into a purified water stream where an electrolyte conductivity sensor quantifies the resulting ammonium concentration. 
#The realized system is shown to be selective enough not to be influenced by normal ambient carbon dioxide concentrations. 
#Experiments with a gas flow of 3 ml/min, containing ammonia concentrations ranging from 9.8 to 0.3 ppm in a nitrogen carrier flow, into a 15 μl/min sample solution flow and finally into a 5 μl/min purified water stream have been carried out and show that the system is sensitive to ammonia concentration below 1 ppmv.

#Known interferences: CO2

#See details of method here:https://www.sciencedirect.com/science/article/pii/S0003267003012650?via%3Dihub

'electrochemical membrane diffusion': {'abbreviation':'EMD', 'sampling_type':'low volume continuous', 'measured_parameters':['NH3'], 'qa_accepted_parameters':['NH3'], 
    'instruments':{}
},

#------------------------------------
#photoacoustic spectroscopy (PS)

#Gaseous samples are continuously drawn through a measurement cell, where they are interrogated spectroscopically with the output radiation of a carbon dioxide/tunable quantum cascade laser. 
#This is used in conjunction with a highly sensitive acoustic detector to measure trace gas concentrations. Laser tuning is accomplished using proprietary algorithms that do not require the use of a spectrometer or other wavelength measurement device. 
#The laser is tuned to a known gaseous absorption line and passed through a cell containing the gas sample. If the trace gas species is present, the gas sample will be slightly heated. This heating can be measured very accurately and linearly by microphones in the cell and the amplitude of the electrical signal from the microphones correlate with the trace gas concentration. 
#If there is no trace gas present, there will be no signal from the microphone. This is, therefore, the highly desirable "zero-background measurement". 
#Traditional spectroscopic methods measure the total power of the laser beam with and without a sample in the path. For very weak attenuations that result from trace gas concentrations, this process requires taking the difference of two large numbers, each with a finite uncertainty, to compute a small real number. This complication is avoided using this method.
#This method has been used to measure NH3.

'photoacoustic spectroscopy':{'abbreviation':'PS', 'sampling_type':'low volume continuous', 'measured_parameters':['NH3'], 'qa_accepted_parameters':['NH3'], 
    'instruments':{
        'Pranalytica Nitrolux100': {'documented_lower_limit_of_detection':'0.1', 'documented_upper_limit_of_detection':'300.0',  'documented_precision':'0.1/10%', 'documented_accuracy':'0.1/10%', 'documented_zero_drift':'0.1/week', 'documented_span_drift':'10.0%/month', 'documented_flow_rate':'1.6', 'instrument_manual_name':'Pranalytica_Nitrolux_Specs.pdf,Pranalytica_Nitrolux100_Specs.pdf'},
        'Pranalytica Nitrolux200': {'documented_lower_limit_of_detection':'0.2', 'documented_upper_limit_of_detection':'500.0',  'documented_precision':'0.2/10%', 'documented_accuracy':'0.2/10%', 'documented_zero_drift':'0.2/week', 'documented_span_drift':'10.0%/month', 'documented_flow_rate':'1.6', 'instrument_manual_name':'Pranalytica_Nitrolux_Specs.pdf',},
        'Pranalytica Nitrolux1000':{'documented_lower_limit_of_detection':'1.0', 'documented_upper_limit_of_detection':'2000.0', 'documented_precision':'1.0/10%', 'documented_accuracy':'1.0/10%', 'documented_zero_drift':'1.0/week', 'documented_span_drift':'10.0%/month', 'documented_flow_rate':'1.6', 'instrument_manual_name':'Pranalytica_Nitrolux_Specs.pdf'}
    }
},

#------------------------------------
#non-dispersive infrared absorption (luft) (NDIR-L)

#Operates on the principle that a specific species efficiently absorbs IR light at a known wavelength in the IR range. This is the case for CO, at 4.7um
#The degree to which the IR light is absorbed by a specific species is directly related to the species concentration as described by the Beer-Lambert Law (I/Io = e−KLC; K = molecular absorption coefficient at STP, L = optical path length of cell, C = species concentration , I = light intensity of sample gas, Io = light intensity of sample without measured species (reference gas) ).

#Typical measurement operation for CO:
#1. The sample is drawn into the through the sample bulkhead
#2. IR radiation from a IR source (at 4.7 um) is passed through a rotating shutter creating a series of pulses.
#3. These pulses are split simultaneously into 2 measurement cells. One cell contains the sample gas, termed the sample stream (I), and cell contains a non-absorbing gas at the specific wavelength (i.e N2), termed the reference stream (Io).
#4. Both the pulses then exit the cells and fall on a IR photo-detector. 
#5. The intensities of 4.7 um light transmitted through the sample (I) and CO-free reference beams (Io) are related to the concentration (C) of CO in the sample gas stream according to the Beer-Lambert Law.
#6. The analyser’s microprocessor calculates the absolute concentration in molecules cm-3, and is then converted to a mixing ratio using in-instrument measured sample temperature and pressure.

#Known Interferences:
#Water vapour and CO2

'non-dispersive infrared absorption (luft)':{'abbreviation':'NDIR-L', 'sampling_type':'low volume continuous', 'measured_parameters':['CO','CH4','CO2','NH3','NO','SO2','H2O'], 'qa_accepted_parameters':['CO'], 
    'instruments':{
        'Beckman 866':      {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (luft) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Bendix 8501-5CA':  {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (luft) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Horiba 300E':      {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (luft) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Horiba 300SE':     {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (luft) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Horiba AQM10':     {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (luft) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Horiba AQM11':     {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (luft) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Horiba AQM12':     {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (luft) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Massachusetts CO1':{'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'cannot find any information. AQS states it as non-dispersive infrared absorption. Assume this is correct, and assume luft type', 'measuring_assumption':True},
        'MSA 202S':         {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (luft) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'}, 
        'Signal 418':       {                              'documented_lower_limit_of_detection':{'CO':'1000.0','CH4':'2000.0','CO2':'300.0','NH3':'1000.0','NO':'5000.0','SO2':'2000.0','H2O':'5000.0'}, 'documented_precision':'0.1%', 'documented_accuracy':'0.5%', 'documented_zero_drift':'0.5%/day', 'documented_span_drift':'0.5%/day', 'documented_flow_rate':'0.2-2.0', 'instrument_manual_name':'Signal_418_Specs.pdf'},
        'Teledyne ML8310':  {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (luft) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Teledyne ML8830':  {'measured_parameters':['CO'],                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption here (assume luft):https://books.google.es/books?id=7eL1Rd8hB-0C&pg=PA418&lpg=PA418&dq=monitor+labs+8830&source=bl&ots=GlI29t_UB9&sig=uqlixhqiWdhaQf8BvZ5GQmu_1gw&hl=en&sa=X&ved=0ahUKEwi2rNOCweLbAhVDrRQKHVf7AsoQ6AEIVDAF#v=onepage&q=monitor%20labs%208830&f=false', 'measuring_assumption':True}
    }
},
#------------------------------------
#non-dispersive infrared absorption (gas-filter correlation) (NDIR-GFC)

#Operates on the principle that a specific species efficiently absorbs IR light at a known wavelength in the IR range. This is the case for CO, at 4.7um
#The degree to which the IR light is absorbed by a specific species is directly related to the species concentration as described by the Beer-Lambert Law (I/Io = e−KLC; K = molecular absorption coefficient at STP, L = optical path length of cell, C = species concentration , I = light intensity of sample gas, Io = light intensity of sample without measured species (reference gas) ).
#The methodology also applies the gas-filter correlation technique. This is applied as water vapour and CO2 absorb light at 4.7 um as well as CO, and ensures measurements are specific for CO.

#Typical measurement operation for CO:
#1. The sample is drawn into the through the sample bulkhead
#2. IR radiation from a IR source (at 4.7 um) is passed through a multi-pass cell length with sample gas. The sample cell uses mirrors at each end to reflect the IR beam back and forth through the sample gas a number of times. The total length that the reflected light travels is directly related to the intended sensitivity of the instrument. The lower the concentrations the instrument is designed to detect, the longer the light path must be in order to create detectable levels of attenuation.
#A gas-filter correlation wheel is combined with this system. This wheel contains three parts to increase measurement accuracy: CO, N2 and the mask:
#    The CO window contains saturation (40 %) of CO which acts as a reference beam – absorbing a known amount of light.
#    The N2 window, containing 100 % N2, does not absorb IR at 4.7 um at all and is used during normal CO measurement.
#    The mask totally blocks the light source and is used to determine background signals, the strength of other signals relative to each other and the background.
#3. The IR source is passed through the rotating filter wheel before entering the multi-pass cell: 
#    When the light beam is intercepted by the CO portion of the wheel, the carbon monoxide, which is at relatively high concentration, absorbs all wavelengths that are co-specific, creating and emanating light beam that is “CO blind”. This “optically scrubbed” portion of the beam is designated the Reference beam (Io).
#    The nitrogen-intercepted portion of the beam, which is “CO sensitive”, is designated the sample beam (I). 
#    In the absence of CO no attenuation of the Io and I portion of the beam will occur.
#    Species other than CO will cause an equal attenuation of both I and Io portions of the beam. 
#    If CO is present in the air being sampled, then the beam portion generated by the CO side of the wheel will experience no attenuation, but the beam portion generated by the N2 portion of the wheel will be attenuated to the degree dictated by the level of CO concentration.
#    When the light beam is intercepted by the mask, this is labelled the "dark portion". This provides a zero light reference point to compensate for the "dark current" of the detector.
#4. A band pass filter is fitted on the exit of the sample cell, allowing only 4.7 um wavelength light to pass.
#5. The infrared radiation then exits the optical cell and falls on an IR photo-detector. 
#6. The intensities of 4.7 um light transmitted through the sample (I) and CO-free reference beams (Io) are related to the concentration (C) of CO in the sample gas stream according to the Beer-Lambert Law.
#7. The analyser's microprocessor records the imbalance between the I and Io beams portions and performs a data linearisation, calculating the absolute concentration in molecules cm-3, and is then converted to a mixing ratio using in-instrument measured sample temperature and pressure.

#Known Interferences:
#Water vapour and CO2

'non-dispersive infrared absorption (gas-filter correlation)':{'abbreviation':'NDIR-GFC', 'sampling_type':'low volume continuous', 'measured_parameters':['CO'], 'qa_accepted_parameters':['CO'], 
    'instruments':{
        'Dasibi 3003':           {                                                                                                                                                                                                                                                                                                                 'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (gas-filter correlation) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Dasibi 3008':           {'documented_lower_limit_of_detection':'100.0', 'documented_upper_limit_of_detection':'1000000.0',  'documented_precision':'100.0',     'documented_accuracy':'1.0%',                           'documented_zero_drift':'200.0/day', 'documented_span_drift':'1.0%/day',        'documented_flow_rate':'1.0',     'instrument_manual_name':'Dasibi_3008_Manual.pdf'},
        'DKK-TOA GFC-311E':      {'documented_lower_limit_of_detection':'400.0', 'documented_upper_limit_of_detection':'100000.0',                                                                                               'documented_zero_drift':'100.0/day', 'documented_span_drift':'2.0%/day',                                          'instrument_manual_name':'DKK-TOA_GFC-311E_Specs.pdf'},
        'Ecotech EC9830T':       {'documented_lower_limit_of_detection':'20.0',  'documented_upper_limit_of_detection':'20000.0',    'documented_precision':'20/0.1%',   'documented_accuracy':'1.0%',                           'documented_zero_drift':'50.0/day',  'documented_span_drift':'0.5%/day',        'documented_flow_rate':'1.0',     'instrument_manual_name':'Ecotech_EC9830T_Specs.pdf'},
        'Ecotech Serinus30':     {'documented_lower_limit_of_detection':'40.0',  'documented_upper_limit_of_detection':'200000.0',   'documented_precision':'20/0.1%',   'documented_accuracy':'1.0%',                           'documented_zero_drift':'100.0/day', 'documented_span_drift':'0.5%/day',        'documented_flow_rate':'1.0',     'instrument_manual_name':'Ecotech_Serinus30_Manual.pdf'},
        'Environnement SA CO10M':{                                                                                                                                                                                                                                                                                                                 'instrument_further_details':':No manual or specifications available but confirmation it is uses non-dispersive infrared absorption here (assume gas-filter correlation, as latest version of instrument uses it):http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/enviroCO10M/view', 'measuring_assumption':True},
        'Environnement SA CO11M':{                                                                                                                                                                                                                                                                                                                 'instrument_further_details':':No manual or specifications available but confirmation it is uses non-dispersive infrared absorption here (assume gas-filter correlation, as next version of instrument uses it):http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf', 'measuring_assumption':True},
        'Environnement SA CO12M':{'documented_lower_limit_of_detection':'40.0',  'documented_upper_limit_of_detection':'200000.0',                                       'documented_accuracy':'1.0%',                           'documented_zero_drift':'500.0/day', 'documented_span_drift':'0.5%/day',        'documented_flow_rate':'1.0',     'instrument_manual_name':'EnvironnementSA_CO12M_Specs.pdf'},
        'Maihak Unor 6N':        {                                                                                                                                                                                                                                                                                                                 'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption here (assume gas-filter correlation):https://books.google.es/books?id=ZDb3koqwApMC&pg=PA468&lpg=PA468&dq=maihak+unor+6n&source=bl&ots=EP38JON-AR&sig=H6eppHv6VGHfbrPaqjOBREST0mE&hl=en&sa=X&ved=2ahUKEwiLmbzS1O_fAhXiShUIHQBqArsQ6AEwC3oECAEQAQ#v=onepage&q=maihak%20unor%206n&f=false', 'measuring_assumption':True}, 
        'Philips K50109':        {                                                                                                                                                                                                                                                                                                                 'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (gas-filter correlation):https://www.isiaq.org/docs/papers/472.pdf'}, 
        'SIR S-5006':            {'documented_lower_limit_of_detection':'40.0',  'documented_upper_limit_of_detection':'200000.0',   'documented_precision':'100.0',     'documented_accuracy':'1.0%',                           'documented_zero_drift':'10.0/day',  'documented_span_drift':'0.5%/day',        'documented_flow_rate':'1.0',     'instrument_manual_name':'SIR_S-5006_Specs.pdf'},
        'Teledyne 300':          {                                                                                                                                                                                                                                                                                                                 'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption here (assume gas-filter correlation, as other instruments with very similar codes use it):http://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf', 'measuring_assumption':True}, 
        'Teledyne 300E':         {'documented_lower_limit_of_detection':'40.0',  'documented_upper_limit_of_detection':'1000000.0',  'documented_precision':'0.5%/200',  'documented_accuracy':'1.0%',                           'documented_zero_drift':'100.0/day', 'documented_span_drift':'100.0/0.5%/day',  'documented_flow_rate':'0.8',     'instrument_manual_name':'Teledyne_300E_Manual.pdf'},
        'Teledyne 300EM':        {'documented_lower_limit_of_detection':'200.0', 'documented_upper_limit_of_detection':'5000000.0',  'documented_precision':'1.0%/1000', 'documented_accuracy':'1.0%',                           'documented_zero_drift':'500.0/day', 'documented_span_drift':'500.0',           'documented_flow_rate':'0.8',     'instrument_manual_name':'Teledyne_300E_Manual.pdf'},
        'Teledyne 300EU':        {'documented_lower_limit_of_detection':'20.0',  'documented_upper_limit_of_detection':'100000.0',   'documented_precision':'0.5%',      'documented_accuracy':'1.0%',                           'documented_zero_drift':'20.0/day',  'documented_span_drift':'0.5%/day',        'documented_flow_rate':'1.8',     'instrument_manual_name':'Teledyne_300EU_Addendum.pdf'},
        'Teledyne ML8930':       {                                                                                                                                                                                                                                                                                                                 'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption here (assume gas-filter correlation): http://dd.eionet.europa.eu/vocabularyconcept/aq/measurementequipment/ML8930/view?facet=HTML+Representation', 'measuring_assumption':True},        
        'Teledyne ML9830':       {                                                                                                                                                                                                                                                                                                                 'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (gas-filter correlation) here:https://twobtech.com/pdf/211%20FEM%20Modification_June2018.pdf'},    
        'Teledyne ML9830B':      {                                                                                                                                                                                                                                                                                                                 'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (gas-filter correlation) here:https://twobtech.com/pdf/211%20FEM%20Modification_June2018.pdf'},    
        'Teledyne T300':         {'documented_lower_limit_of_detection':'40.0',  'documented_upper_limit_of_detection':'1000000.0',  'documented_precision':'0.5%/200',  'documented_accuracy':'1.0%',                           'documented_zero_drift':'100.0/day', 'documented_span_drift':'0.5%/day',        'documented_flow_rate':'0.8',     'instrument_manual_name':'Teledyne_T300_Manual.pdf'},
        'Teledyne T300M':        {'documented_lower_limit_of_detection':'200.0', 'documented_upper_limit_of_detection':'1000000.0',  'documented_precision':'1.0%/1000', 'documented_accuracy':'1.0%',                           'documented_zero_drift':'500.0/day', 'documented_span_drift':'500.0/day',       'documented_flow_rate':'0.8',     'instrument_manual_name':'Teledyne_T300_Manual.pdf'},
        'Teledyne T300U':        {'documented_lower_limit_of_detection':'20.0',  'documented_upper_limit_of_detection':'100000.0',   'documented_precision':'0.5%',      'documented_accuracy':'1.0%',                           'documented_zero_drift':'20.0/day',  'documented_span_drift':'0.5%>=500.0/day', 'documented_flow_rate':'1.8',     'instrument_manual_name':'Teledyne_T300U_Addendum.pdf'},
        'Thermo 48':             {                                                                                                                                                                                                                                                                                                                 'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (gas-filter correlation) here:https://books.google.es/books?id=yCDGFlP3toUC&pg=PA248&lpg=PA248&dq=HORIBA+AQM-10&source=bl&ots=-uIdMsd0Yg&sig=U34lJHqy1fXno_AwV36rZIKCwuE&hl=en&sa=X&ved=0ahUKEwjvlL2IxuLbAhUJthQKHZCQB90Q6AEIQjAC#v=onepage&q=HORIBA%20AQM-10&f=false'},
        'Thermo 48C':            {'documented_lower_limit_of_detection':'40.0',  'documented_upper_limit_of_detection':'10000000.0', 'documented_precision':'100.0',     'documented_accuracy':'1.0%',                           'documented_zero_drift':'100.0/day', 'documented_span_drift':'1.0%/day',        'documented_flow_rate':'0.5-2.0', 'instrument_manual_name':'Thermo_48C_Specs.pdf'},
        'Thermo 48C-TL':         {                                                                                                                                                                                                                                                                                                                 'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (gas-filter correlation) here:http://www.broward.org/Environment/AirQuality/AirMonitoring/Pages/AmbientAirMonitoringInstruments.aspx'},
        'Thermo 48i':            {'documented_lower_limit_of_detection':'40.0',  'documented_upper_limit_of_detection':'10000000.0',                                     'documented_accuracy':'1.0%<1000000.0 2.5%>=1000000.0', 'documented_zero_drift':'100.0/day', 'documented_span_drift':'1.0%/day',        'documented_flow_rate':'1.0',     'instrument_manual_name':'Thermo_48i_Manual.pdf'},
        'Thermo 48i-TLE':        {'documented_lower_limit_of_detection':'40.0',  'documented_upper_limit_of_detection':'1000000.0',                                      'documented_accuracy':'1.0%/40.0',                      'documented_zero_drift':'100.0/day', 'documented_span_drift':'1.0%/day',        'documented_flow_rate':'0.5',     'instrument_manual_name':'Thermo_48i-TLE_Manual.pdf'}
    }
},

#------------------------------------
#non-dispersive infrared absorption (cross-flow modulation) (NDIR-CFM)

#Operates on the principle that a specific species efficiently absorbs IR light at a known wavelength in the IR range. This is the case for CO, at 4.7um.
#The degree to which the IR light is absorbed by a specific species is directly related to the species concentration as described by the Beer-Lambert Law (I/Io = e−KLC; K = molecular absorption coefficient at STP, L = optical path length of cell, C = species concentration , I = light intensity of sample gas, Io = light intensity of sample without measured species (reference gas) ).
#The methodology also applies the cross-flow modulation technique. This is applied as water vapour and CO2 absorb light at 4.7 um as well as CO, and ensures measurements are specific for CO.

#Typical measurement operation for CO:
#1. The sample is drawn into the through the sample bulkhead
#2. IR radiation from a IR source (at 4.7 um) is passed through a measurement cell.
#3. Fixed amounts of the sample gas (I) and the reference gas (Io) are injected alternately into the measurement cell (cross-flow modulation). This has the advantage over the gas-filter correlation technique that, in principle, when analyzing small amounts of gas there is no generation of zero-drift. An additional advantage is that the elimination of rotary sectors precludes the need for optical adjustment. These features assure greatly improved stability over long periods of measurement. The reference gas is generated by purging the sample through an oxidation process, where an oxidizing catalyst burns the CO to CO2. These features eliminate interference from other elements, resulting in highly accurate measurements.
#4. The infrared radiation then exits the optical cell and falls on an IR photo-detector. 
#5. The intensities of 4.7 um light transmitted through the sample (I) and CO-free reference beams (Io) are related to the concentration (C) of CO in the sample gas stream according to the Beer-Lambert Law.
#6. The analyzer’s microprocessor calculates the absolute concentration in molecules cm-3, and is then converted to a mixing ratio using in-instrument measured sample temperature and pressure.

#Known Interferences:
#Water vapour and CO2

'non-dispersive infrared absorption (cross-flow modulation)':{'abbreviation':'NDIR-CFM', 'sampling_type':'low volume continuous', 'measured_parameters':['CO'], 'qa_accepted_parameters':['CO'], 
    'instruments':{
        'Horiba APMA300':{                                                                                                                                                                                                                                                                   'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption here (assume cross-flow modulation): https://books.google.es/books?id=qmWXBwAAQBAJ&pg=PA88&lpg=PA88&dq=horiba+apna+300&source=bl&ots=dBQJ0TEai1&sig=plLm6SAUIZwuIhy37JzI8nZXccQ&hl=en&sa=X&ved=2ahUKEwihxOS9uvDfAhWGnhQKHYzoAPcQ6AEwCHoECAIQAQ#v=onepage&q=horiba%20apna%20300&f=false', 'measuring_assumption':True},
        'Horiba APMA350':{                                                                                                                                                                                                                                                                   'instrument_further_details':'No manual or specifications available but confirmation it is uses non-dispersive infrared absorption (cross-flow modulation) here: http://www.wseas.us/e-library/conferences/2009/timisoara/SSE2/SSE2-11.pdf'},
        'Horiba APMA360':{'documented_lower_limit_of_detection':'20.0', 'documented_upper_limit_of_detection':'100000.0', 'documented_precision':'1.0%', 'documented_accuracy':'1.0%', 'documented_zero_drift':'20.0/day', 'documented_span_drift':'20.0/day', 'documented_flow_rate':'1.5', 'instrument_manual_name':'Horiba_APMA360_Specs.pdf'}, 
        'Horiba APMA370':{'documented_lower_limit_of_detection':'20.0', 'documented_upper_limit_of_detection':'100000.0', 'documented_precision':'1.0%', 'documented_accuracy':'1.0%', 'documented_zero_drift':'20.0/day', 'documented_span_drift':'20.0/day', 'documented_flow_rate':'1.5', 'instrument_manual_name':'Horiba_APMA370_Manual.pdf'}
    }
},

#------------------------------------
#dual isotope fluorescence (DIF)  

#This method is an adaption of non-dispersive infrared absorption for the measurement of CO.
#A broadband of infrared source is used to stimulate fluorescence of C16O and C18O, contained in a sealed cell.
#The stimulated radiation is passed first through a filter that alternately passes radiation from only one of the isotopic sources and is then passed through the sample cell to a solid state detector.
#Because almost all atmospheric CO consists of the C16O isotope, only radiation produced by the C16O source will be absorbed by the CO in the sample.
#The radiation produced by the C18O source is used as a reference beam to account for interference with H2O, or other interferants.

#Known Interferences:
#Water vapour and CO2

#See further details: https://books.google.es/books?id=V0oEzMp7axAC&pg=PA385&lpg=PA385&dq=dual+isotope+fluorescence+carbon+monoxide&source=bl&ots=azm5b16AqC&sig=pDZOFxw9LqZxQMsW8wOSxKv1YSE&hl=en&sa=X&ved=0ahUKEwiR7o2iyOLbAhVEtRQKHSFECGEQ6AEINjAB#v=onepage&q=dual%20isotope%20fluorescence%20carbon%20monoxide&f=false

'dual isotope fluorescence':{'abbreviation':'DIF', 'sampling_type':'low volume continuous', 'measured_parameters':['CO'], 'qa_accepted_parameters':['CO'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography (GC)

#Gas chromatography (GC) is a method used for separating and analysing compounds that can be vaporized without decomposition.
#A sample solution is injected into a instrument, entering a gas stream which transports the sample (mobile phase) into a separation tube known as the "column".
#The mobile phase is a carrier gas, usually an inert gas such as helium or an unreactive gas such as nitrogen. Helium remains the most commonly used carrier gas in about 90% of instruments although hydrogen is preferred for improved separations.
#The column consists of a microscopic layer of liquid or polymer on an inert solid support a microscopic layer of liquid or polymer on an inert solid support (stationary phase), inside a piece of glass or metal tubing, placed inside a piece of glass or metal tubing.
#Once inside the column, the gaseous compounds being analysed interact with the walls of the column coated with a stationary phase. This causes each compound to elute at a different time, known as the retention time of the compound. 
#The comparison of retention times is what gives GC its analytical usefulness.
#If greater separation of compounds is required, multiple distinct columns can be used for this purpose.

#Almost all methods involving GC are typically paired with a detector of some description to measure the GC separated components.
#These detectors are various in number, and are outlined in detail below.
#Typically GC-detector analysis is used for the separation of short lived hydrocarbon species. 

#Typically sample processing is undertaken before the passing of samples into the GC. 
#Typically this emcompasses the preconcentration, or passing samples through sorbent traps.

#If no other information is given, it is assumed all GC-detector analysis is done via injection and do not operate online (i.e. not continuously).

#------------------------------------
#gas chromatography - electron capture detection (GC-ECD)

#The Electron Capture Detector (ECD) is a selective detector for electro-negative compounds, especially halogens.
#A beta-ray radio active source which can ionize the carrier gas is located in the detector. 
#A current is produced between two electrodes in the detector supplied with a potential difference and this is monitored as a continuous background current. 
#When there are electro-negative components present in the carrier gas, the background current is reduced because these components capture electrons.
#This change in the background current is registered as 'holes' in the background. 
#Since it concerns very small holes, the background current should not be too high. 
#Therefore the (linear) dynamic range of the ECD is often limited.

'gas chromatography - electron capture detection': {'abbreviation':'GC-ECD', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE','PAN'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE','PAN'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - flame ionisation detection (GC-FID)

#Flame ionisation detection operates in the exact same way as it does as a standalone method, only in this instance coupled to a GC.
#It is based on the principle of the generation of an electrical current that is proportional to the rate of ion formation, dependent on the concentrations of species in the sample gas.

#A limitation of GC-FID is that the retention time is merely characteristic of a compound and certainly not adjudicative or confirmatory of the specificity of that compound.
#This means a peak at a given retention time is not a unique qualitative measure (to the exclusion of every other compound in the universe).

'gas chromatography - flame ionisation detection':{'abbreviation':'GC-FID', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - photoionisation detection (GC-PID)

#Gas Chromatography - Photoionisation Detection (GC-PID) is a technique used to analyse a wide range of aromatic hydrocarbons and other organic compounds. 
#A typical application is the analysis of hydrocarbon pollution of water. The PID uses ultraviolet light to ionize the components exiting the column. 
#The ions are collected by electrodes and the current generated measures the concentration. 
#As common with other GC techniques, a carrier gas is required with low Water and Oxygen impurities since Water and Oxygen can interact with the stationary phase and cause significant problems such as high baseline noise and column bleed in the output gas chromatogram which both reduces the analyser sensitivity and decreases column lifetime.  
#Helium and Nitrogen are normally used. 
#A make-up gas is used in GC-PID at the exit of the column to increase the flow rate into the detector. 
#The routine calibration of the analyser using a calibration mixture is common.

'gas chromatography - photoionisation detection':{'abbreviation':'GC-PID', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - mercuric oxide reduction detection (GC-HgO) 

#Mercuric oxide (HgO) reduction is a detection method available for hydrogen, carbon monoxide, and other reducing gases.
#It uses mercuric oxide reaction tube and a mercury lamp in a heated UV detector cell. 
#The reaction tube is heated to 260-300oC. 
#Located immediately downstream of the reaction tube, the UV detector cell is heated to 170degC. 
#The UV detector cell is equipped with a mercury lamp and a UV photodiode. When a reducing gas such as carbon monoxide elutes from the GC column, it reacts with the HgO to form gaseous mercury vapor, which is then swept into the UV cell. 
#The gaseous mercury absorbs the UV light from the mercury lamp as it flows through the cell. The change in transmittance is converted by the data system into an absorbance output, which is linearly proportional to the amount of reducing gas.
#More information here: https://srigc.com/cn/downloads/41/RGD.pdf

'gas chromatography - mercuric oxide reduction detection':{'abbreviation':'GC-HgO', 'sampling_type':'injection', 'measured_parameters':['CO'], 'qa_accepted_parameters':['CO'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - dual flame ionisation detection (GC-DFID)

#To limit the main limitation of GC-FID, i.e. lack of specificity, it is possible to add another GC column and analyse the sample concurrently on both columns.
#This can be done by adding a y-splitter that will take the single injection made into the injector port and divide the sample into two different pathways with one part of the sample going to one column for analysis and the second part going to another column for analysis.
#The two different columns are given two different stationary phases, which provides different separation of the analytes (both in terms of retention time and possibly even eluting order). 
#This different eluting order and different retention times only minimises, but does not entirely eliminate the possibility of co-elution as again, the resolving (separating) power of the method is only determined by one non-unique measure which is the two retention times. 
#Even though there is a change in the eluting order potentially and the retention times are different based upon the stationary phase composition, again, remains the basic limitation of GC-FID remains in that a retention time is merely characteristic of the analyte, but is certainly not confirmatory.
#With dual detectors, the concentration of the analyte(s) of interest is determined from two calibration curves, and can serve as quantitative confirmation.

'gas chromatography - dual flame ionisation detection':{'abbreviation':'GC-DFID', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - fourier transform infrared spectroscopy (GC-FTIR)

#Fourier Transform Infrared Spectroscopy (FTIR) is an analytical technique used to identify organic (and in some cases inorganic) materials. 
#This technique measures the absorption of infrared radiation by the sample material versus wavelength. 
#The infrared absorption bands identify molecular components and structures.
#When a material is irradiated with infrared radiation, absorbed IR radiation usually excites molecules into a higher vibrational state. 
#The wavelength of light absorbed by a particular molecule is a function of the energy difference between the at-rest and excited vibrational states. 
#The wavelengths that are absorbed by the sample are characteristic of its molecular structure.
#The FTIR spectrometer uses an interferometer to modulate the wavelength from a broadband infrared source. 
#A detector measures the intensity of transmitted or reflected light as a function of its wavelength. 
#The signal obtained from the detector is an interferogram, which must be analyzed with a computer using Fourier transforms to obtain a single-beam infrared spectrum. 
#The FTIR spectra are usually presented as plots of intensity versus wavenumber (in cm-1). Wavenumber is the reciprocal of the wavelength. 
#The intensity can be plotted as the percentage of light transmittance or absorbance at each wavenumber.

'gas chromatography - fourier transform infrared spectroscopy':{'abbreviation':'GC-FTIR', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - mass spectrometry (GC-MS)

#One of the most common detector types is a mass spectrometer.
#Mass spectrometry is the study of matter through the formation of gas-phase ions (with or without fragmentation) that are detected and characterised by their mass to charge ratios (m/z) and relative abundances. 
#The technique basically studies the effect of ionizing energy on molecules. It depends upon chemical reactions in the gas phase in which sample molecules are consumed during the formation of ionic and neutral species.
#A mass spectrometer generates multiple ions from the sample under investigation, it then separates them according to their specific mass-to-charge ratio (m/z), and then records the relative abundance of each ion type.
#The first step in the mass spectrometric analysis of compounds is the production of gas phase ions of the compound, basically by electron ionization. 
#This molecular ion undergoes fragmentation. Each primary product ion derived from the molecular ion, in turn, undergoes fragmentation, and so on. 
#The ions are separated in the mass spectrometer according to their mass-to-charge ratio, and are detected in proportion to their abundance. 
#A mass spectrum of the molecule is thus produced. It displays the result in the form of a plot of ion abundance versus mass-to-charge ratio. 
#Ions provide information concerning the nature and the structure of their precursor molecule. 
#In the spectrum of a pure compound, the molecular ion, if present, appears at the highest value of m/z (followed by ions containing heavier isotopes) and gives the molecular mass of the compound.

#There are a variety of different mass spectrometer types (almost never reported in observational data). 
#The most important difference among various mass spectrometers is the type of the analyser. 
#It is a fundamental characteristic of an instrument, and can not be modified. 
#The main analyzer types are: quadrupole, ion trap, time-of-flight (TOF), sector and Fourier-transform (FT) type instruments. 
#Some analysers are inherently capable for MS–MS (like ion trap or FT-MS), some others can be connected in tandem, offering MS–MS capabilities. 
#The most common analyzers (quadrupoles, sectors and in most respect ion traps) pass one ion (a given m / z ratio) at any one time, and filter out all other ions (which are lost for analysis).

#To obtain a mass spectrum, the analyser scans through the mass range of interest. 
#The consequence of this scanning is, that the sensitivity of detecting a given ion is much less compared to the case, when only a given mass is monitored (when the analyser stands still to pass only one ion). 
#These analyzers can be used in two ways: Either monitoring a given ion (with full sensitivity), or scanning the mass range (with loss of sensitivity) to obtain a mass spectrum. 
#Other analyzer types (TOF, FT-MS, in some respect ion traps and special sector instruments) can simultaneously detect all ions in a mass range, so there is no need for scanning. 
#In this case the full sensitivity is obtained irrespectively if one ion or a full mass spectrum is studied. 
#At present quadrupole-type instruments are most widespread. These are simple to use, cheap and rugged, good for basic GC–MS studies and for quantitation.

#In many instances the phrase 'mass selective detector' or 'MSD' is be used in place of mass spectrometer.
#However, there is no fundamental difference between mass selective detection and mass spectrometry; the name 'mass selective detection' is part of a marketing strategy, and implies a simple, relatively cheap instrument.
#Therefore all reported instances of MSD are grouped under mass spectrometry.

#Determination of Hg: https://www.atmos-meas-tech.net/9/2195/2016/amt-9-2195-2016.pdf
#Determination of Ni: https://www.longdom.org/open-access/gcms-determination-of-the-isotope-distribution-of-nickel-2157-7064-1000397.pdf
#Determination of OC: https://assets.thermofisher.com/TFS-Assets/CMD/Application-Notes/ANCCSGCTOCAN_0311.pdf

'gas chromatography - mass spectrometry':{'abbreviation':'GC-MS', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 
    'instruments':{}
},

#------------------------------------
#pyrolysis - gas chromatography - mass spectrometry (Py-GC-MS)

#For the measurement of black carbon, GCMS can be used following pyrolysis.
#See more information here: https://www.sciencedirect.com/science/article/abs/pii/S014663800800185X and https://digital.csic.es/handle/10261/55429

'pyrolysis - gas chromatography - mass spectrometry':{'abbreviation':'Py-GC-MS', 'sampling_type':'injection', 'measured_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 'qa_accepted_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - direct temperature resolved mass spectrometry (GC-DTMS) 

#Direct Temperature Resolved Mass Spectrometry comprises of thermal depolymerization of compounds into smaller fragments that are either immediately analyzed by mass spectrometry.
#This technique has the advantage that compounds with high boiling points are prevented from condensation on 'cold' surfaces.
#The method is capable of detecting both extractable and non-extractable more complex organic fractions.
#It is typically used to give detailed information about a broad range of organic compounds such as lipids, waxes, terpenoids, polynuclear aromatic compounds, oligosaccharides, small peptides, protein fragments, and larger more condensed polymeric structures resulting from charred foodstuffs.

'gas chromatography - direct temperature resolved mass spectrometry':{'abbreviation':'GC-DTMS', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - mass spectrometry - flame ionisation detection (GC-MS-FID)

#In some instances, the simultaneous application of different detectors is extremely advantageous.
#For example, the combination of an a flame ionisation detector (which has a wide quantitation response range) with a mass spectrometer (which obtains qualitative information of unknown compounds) enables identification of volatile and semi-volatile compounds at trace levels.
#A splitter is applied after the GC column, splitting output to an FID and MS detector. 
#Using such a system enables users to obtain both a TIC (Total Ion Chromatogram) and FID chromatogram at the same time, giving the advantages of both detection types simultaneously.
#By acquiring abundant information at the same time, the detector splitting system saves both cost for analytical instruments and columns, and operators' time.

'gas chromatography - mass spectrometry - flame ionisation detection':{'abbreviation':'GC-MS-FID', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - mass spectrometry - photoionisation detection (GC-MS-PID)

#In the same way as GC-MS-FID, simultaneous detection using MS and PID detectors through splitting of GC output gives simultaneous advantage of both detection types. 

'gas chromatography - mass spectrometry - photoionisation detection':{'abbreviation':'GC-MS-PID', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - flame ionisation detection - electron capture detection (GC-FID-ECD)

#In the same way as GC-MS-FID, simultaneous detection using FID and ECD detectors through splitting of GC output gives simultaneous advantage of both detection types. 

'gas chromatography - flame ionisation detection - electron capture detection':{'abbreviation':'GC-FID-ECD', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE','PAN'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE','PAN'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - flame ionisation detection - photoionisation detection (GC-FID-PID)

#In the same way as GC-MS-FID, simultaneous detection using FID and PID detectors through splitting of GC output gives simultaneous advantage of both detection types. 

'gas chromatography - flame ionisation detection - photoionisation detection':{'abbreviation':'GC-FID-PID', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - fourier transform infrared spectroscopy - mass spectrometry (GC-FTIR-MS) 

#In the same way as GC-MS-FID, simultaneous detection using FTIR and MS detectors through splitting of GC output gives simultaneous advantage of both detection types. 

'gas chromatography - fourier transform infrared spectroscopy - mass spectrometry':{'abbreviation':'GC-FTIR-MS', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},

#------------------------------------
#gas chromatography - unknown detection (GC-?)

'gas chromatography - unknown detection':{'abbreviation':'GC-?', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'],  
    'instruments':{}
},

#------------------------------------
#high performance liquid chromatography - mass spectrometry (HPLC-MS)

#Liquid chromatography (LC) is a separation process used to isolate the individual components of a mixture. 
#This process involves mass transfer of a sample through a polar mobile phase and non-polar stationary phase.
#The device is a column packed with the porous medium made of a granular solid material (i.e., stationary phase), such as polymers and silica, where the sample is injected and the solvent (i.e., mobile phase) passes to transport the sample.
#When a sample is injected, it is adsorbed on the stationary phase, and the solvent passes through the column to separate the compounds one by one, based on their relative affinity to the packing materials and the solvent. 
#The component with the most affinity to the stationary phase is the last to separate. This is because high affinity corresponds to more time to travel to the end of the column.

#High-performance liquid chromatography (HPLC), also known as high-pressure liquid chromatography, is an advanced type of LC. 
#HPLC is amenable to a wide range of applications, such as pharmaceuticals and food analysis. It is especially useful for low or non-volatile organic compounds, which cannot be handled with gas chromatography.
#The difference between traditional LC and HPLC is that the solvent in LC travels by the force of gravity. 
#In HPLC, the solvent travels under high pressure obtained by means of a pump to overcome the pressure drop in the packed column, which reduces the time of separation. 

#As with gas chromatography, liquid chromatography can be combined with mass spectrometry detection.
#Combining the two analytical methods reduces experimental error and improves accuracy. This technique is very useful in applications that involve a huge number of compounds, such as environmental effluents.

'high performance liquid chromatography - mass spectrometry':{'abbreviation':'HPLC-MS', 'sampling_type':'injection', 'measured_parameters':['Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'qa_accepted_parameters':['Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'],
    'instruments':{}
},

#------------------------------------
#high performance liquid chromatography - ultraviolet detection (HPLC-UV)

#The most commonly used detector in liquid chromatography is ultra violet absorption
#A UV detector uses a molecules ability to absorb ultraviolet light – detectors can use a single wavelength or variable wavelengths which can be more sensitive. 
#By passing UV light through the individual components of an eluting sample mix and measuring the amount of UV light absorbed by each component, you can determine the individual quantities of each component.

'high performance liquid chromatography - ultraviolet detection':{'abbreviation':'HPLC-UV', 'sampling_type':'injection', 'measured_parameters':['CH4','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','NH4+','PM10_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'qa_accepted_parameters':['CH4','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'],
    'instruments':{}
},

#------------------------------------
#high performance liquid chromatography - fluorescence detection (HPLC-FLD)

#Fluorescence detection can be a strong alternative to ultraviolet or other detectors for some compounds. 

'high performance liquid chromatography - fluorescence detection':{'abbreviation':'HPLC-FLD', 'sampling_type':'injection', 'measured_parameters':['Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'qa_accepted_parameters':['Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'],
    'instruments':{}
},

#------------------------------------
#proton transfer reaction

#An alternative to using gas chromatography for the separation of compounds is proton transfer reaction. 
#The main advantage of proton transfer reaction is that in most of cases it does not induce any significant fragmentation, and when it does only few fragments are produced. 
#Soft-ionisation results in only one or two characteristic ions, and the matrix of signals is much less complex than in other mass spectrometric techniques.
#Upon entering a PTR device, gas-phase molecules are ionised (by an ionisation source), transferring a proton (H3O+, O2+, NO+, or others) to the molecules of interest.
#For the proton transfer reaction to occur, the ionised molecule must have a higher proton affinity than the transferred proton (i.e. H3O+).
#Fast switching between reagent ions can also be used to separate isomeric compounds (molecules containing same atoms, but having different conformation), which have different reactivity with the various reagent ions.
#The ions produced in this reaction are subsequently detected by a mass analyser.
#Using PTR, comes with the advantage that ambient air can be directly sampled, without any prior sample processing (i.e. preconcentration). 

#------------------------------------
#proton transfer reaction - mass spectrometry (PTR-MS)

#Almost always PTR is coupled with a mass spectrometer for detection (PTR-MS).
#Proton Transfer Reaction - Mass Spectrometry characterises compounds according to their mass-to-charge ratio. 
#A great advantage of PTR-MS is that common air constituents (N2, O2) have low proton affinities and do not react, thus no diluting buffer gas is required.
#PTR-MS is the most widely used method for real-time monitoring of volatile organic compounds (VOCs) at low concentrations.

#Unless another specific detector is stated, where Proton Transfer Reaction is given as a methodology, it will be assumed that it uses Mass Spectrometry detection.

'proton transfer reaction - mass spectrometry':{'abbreviation':'PTR-MS', 'sampling_type':'injection', 'measured_parameters':['CH4','C6H6','ISOPRENE'], 'qa_accepted_parameters':['CH4','C6H6','ISOPRENE'], 
    'instruments':{}
},
#------------------------------------
#Passive methods (colorimetry, spectrophotometry, ion chromatography, titration)

#It is not always possible to measure atmospheric composition using a automatic, continuous methodology. 
#Typically automated instrumentation is expensive, needs dedicated laboratory space, and requires regular measuremental maintenance.
#There is a growing trend in the use of passive methodologies, which are much cheaper, require less maintenance and can be used in high volume. 

#These methods typically involve 2 or more sample steps. Typically these are:
#1. Air is sampled passively, and a desired compound(s) is first trapped (through either physical or chemical means) (passive phase).
#2. The trapped species concentration is determined in a laboratory, using a suitable methodology (typically done in batch) (analysis phase). 

#The in-lab analysis phase methodologies are almost always one of the following: 
#1. colorimetry
#2. spectrophotometry
#3. ion chromatography 
#4. titration

#---------------------
#colorimetry (CO)

#Colorimetry fundamentally refers to the determination of species concentrations based on the colour of a solution.
#Intrinsically, this is based on Beer-Lambert's law, according to which the absorption of light transmitted through the medium is directly proportional to the medium concentration.
#However, the way this colour is extrapolated to give a concentration, adds complexity to the definition of this method.
#In the most basic of senses the colour can be analysed manually, with the concentration directly extrapolated through a colour chart, or given through the length of colour stain on a glass tube from reagent reaction.
#The determination can be done in a more specific, automated fashion through the use of colorimeters.
#In a colorimeter, a beam of light with a specific wavelength is passed through a solution via a series of lenses, which navigate the colored light to the measuring device. 
#This analyses the color compared to an existing standard. A microprocessor then calculates the absorbance or percent transmittance.
#If the concentration of the solution is greater, more light will be absorbed, which can be identified by measuring the difference between the amount of light at its origin and that after passing the solution.
#Colorimeters are usually portable and use LED light sources and color filters. As a result, they operate at fixed wavelengths and can only accommodate tests that incorporate those wavelengths.   

#Unless additional information is given, make assumption colorimetry method uses instrumental injection of sample species, following chemical reaction with a reagent.

#Not QA accepted measurement technique for multiple species due to dominance of other more accurate techniques in measuring these species

'colorimetry':{'abbreviation':'CO', 'sampling_type':'injection', 'measured_parameters':['O3','NO','NO2','NOx','NH3','PAN','CO','SO2','HNO3','Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['PAN','NH3','HNO3','Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{}
},
#---------------------
#spectrophotometry (SP)

#Spectrophotometry is also based on Beer-Lambert's law, where the absorption of light transmitted through the medium is directly proportional to the medium concentration.
#A spectrophotometer is a photometer (a device for measuring light intensity) that can measure intensity as a function of the color, or more specifically, the wavelength of light.  
#Spectrophotometers are usually bench top instruments and use light sources that can produce a range of wavelengths, use monochromators to select for a desired wavelength. As a result, spectrophotometers can be used for a broad range of tests. 

#The distinction between spectrophotometry and colorimetry is not always clear. 
#Colorimeters and spectrophotometers both measure sample absorbance to determine analyte concentrations. 
#The main distinction is that a spectrophotometer can be set to measure % transmittance or absorbance over a wide range of wavelengths, whereas a colorimeter determines absorbance at specific, visible, wavelengths.

#Unless additional information is given, make assumption spectrophotometry method uses instrumental injection of sample species, following chemical reaction with a reagent.

#Not QA accepted measurement technique for multiple species due to dominance of other more accurate techniques in measuring these species

'spectrophotometry':{'abbreviation':'SP', 'sampling_type':'injection', 'measured_parameters':['O3','NO','NO2','NOx','NH3','PAN','CO','SO2','HNO3','Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['PAN','NH3','HNO3','Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{ 
        'Perkin Elmer 18 UV/VIS':{'instrument_further_details':'No manual or specifications available but confirmation it is uses spectrophotometry here: https://www.recycledgoods.com/perkin-elmer-lambda-18-uv-vis-spectrophotometer/'},
        'Shimadzu UV-1201':      {'instrument_further_details':'No manual or specifications available but confirmation it is uses spectrophotometry here: http://www.medwow.com/med/spectrophotometer/shimadzu/uv-1201/12483.model-spec'},
        'Thermo Aquakem 200':    {'instrument_manual_name':'Thermo_Aquekem_200_250_600_Specs.pdf'},
        'Thermo Aquakem 250':    {'instrument_manual_name':'Thermo_Aquekem_200_250_600_Specs.pdf'},
        'Thermo Aquakem 600':    {'instrument_manual_name':'Thermo_Aquekem_200_250_600_Specs.pdf'},
    }
},
#---------------------
#ion chromatography (IC)

#Ion chromatography (or ion-exchange chromatography) is a chromatography process that separates ions and polar molecules based on their affinity to the ion exchanger. 
#It works on almost any kind of charged molecule—including large proteins, small nucleotides, and amino acids. However, ion chromatography must be done in conditions that are one unit away from the isoelectric point of a protein. 
#The two types of ion chromatography are anion-exchange and cation-exchange. 
#Cation-exchange chromatography is used when the molecule of interest is positively charged. The molecule is positively charged because the pH for chromatography is less than the pI. In this type of chromatography, the stationary phase is negatively charged and positively charged molecules are loaded to be attracted to it. 
#Anion-exchange chromatography is when the stationary phase is positively charged and negatively charged molecules (meaning that pH for chromatography is greater than the pI) are loaded to be attracted to it.

#Typical Operation:
#1. A sample of the mixture to be analyzed (analyte) is injected into a carrier fluid (the eluent).
#2. The combination is passed through a column containing a stationary fixed material (adsorbent).
#3. Compounds contained in the analyte are then partitioned between the stationary adsorbent and the moving eluent/analyte mixture.
#4. Different dissolved materials adhere to the adsorbent with different forces. The ones that adhere strongly are moved through the adsorbent more slowly as the eluent flows over them. As the eluent flows through the column the components of the analyte will move down the column at different speeds and therefore separate from one another.
#5. A detector is used to analyze the output at the end of the column. There are various detectors, typically they use an electrical conductivity detector.
#6. Each time analyte molecules/ions emerge from the chromatography column the detector generates a measurable signal which is usually printed out as a peak on the chromatogram.
#7. A suppressor is used to reduce the background conductance of the eluent and at the same time enhance the conductance of the sample ions.

#Not QA accepted measurement technique for NO/NO2/NOx/SO2 due to dominance of other more accurate techniques in measuring these species

#See details of method for NH3:https://www.atmos-meas-tech.net/7/2733/2014/amt-7-2733-2014.html'

'ion chromatography':{'abbreviation':'IC', 'sampling_type':'injection', 'measured_parameters':['NO','NO2','NOx','PAN','NH3','SO2','HNO3','Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['PAN','NH3','HNO3','Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{
        'Dionex DX-120':             {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          'instrument_manual_name':'Dionex_DX-120_Manual.pdf'},
        'Dionex DX-600':             {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          'instrument_further_details':'No manual or specifications available but confirmation it is uses ion chromatography here: https://books.google.es/books?id=9Yq9BgAAQBAJ&pg=PA97&lpg=PA97&dq=ic+dionex+600&source=bl&ots=lUbV9Iwed7&sig=ACfU3U3gUbLCaUNYglSBZQeHrWmVITMVoA&hl=en&sa=X&ved=2ahUKEwigsonXt9roAhUVAWMBHafJAXwQ6AEwEHoECAkQMQ#v=onepage&q=ic%20dionex%20600&f=false'},
        'Dionex 2010i':              {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          'instrument_further_details':'No manual or specifications available but confirmation it is uses ion chromatography here: https://books.google.es/books?id=fortCAAAQBAJ&pg=PA395&lpg=PA395&dq=dionex+2010i&source=bl&ots=qObKPXBKMt&sig=ACfU3U13WHkWJOW7eD9NSNHuggZFjQnGVA&hl=en&sa=X&ved=2ahUKEwiZ2obZiKLgAhW4BGMBHYu7DHMQ6AEwDXoECAEQAQ#v=onepage&q=dionex%202010i&f=false'},
        'Metrohm 881 compact IC pro':{                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          'instrument_manual_name':'Metrohm_881_compact_IC_pro.pdf'},
        'Metrohm MARGA':             {'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['HCl','HNO2','HNO3','NH3','SO2','Ca++','PM10_Ca++','PM2.5_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','K+','PM10_K+','PM2.5_K+','Mg++','PM10_Mg++','PM2.5_Mg++','Na+','PM10_Na+','PM2.5_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','NO3-','PM10_NO3-','PM2.5_NO3-','SO4--','PM10_SO4--','PM2.5_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'qa_accepted_parameters':['HCl','HNO2','HNO3','NH3','SO2','Ca++','PM10_Ca++','PM2.5_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','K+','PM10_K+','PM2.5_K+','Mg++','PM10_Mg++','PM2.5_Mg++','Na+','PM10_Na+','PM2.5_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','NO3-','PM10_NO3-','PM2.5_NO3-','SO4--','PM10_SO4--','PM2.5_SO4--'], 'documented_lower_limit_of_detection':{'HCl':np.NaN,'HNO2':np.NaN,'HNO3':np.NaN,'NH3':np.NaN,'SO2':np.NaN,'Ca++':'0.09','PM10_Ca++':'0.09','PM2.5_Ca++':'0.09','Cl-':'0.01','PM10_Cl-':'0.01','PM2.5_Cl-':'0.01','K+':'0.09','PM10_K+':'0.09','PM2.5_K+':'0.09','Mg++':'0.06','PM10_Mg++':'0.06','PM2.5_Mg++':'0.06','Na+':'0.05','PM10_Na+':'0.05','PM2.5_Na+':'0.05','NH4+':'0.05','PM10_NH4+':'0.05','PM2.5_NH4+':'0.05','NO3-':'0.05','PM10_NO3-':'0.05','PM2.5_NO3-':'0.05','SO4--':'0.04','PM10_SO4--':'0.04','PM2.5_SO4--':'0.04','SO4--_NSS':'0.04','PM10_SO4--_NSS':'0.04','PM2.5_SO4--_NSS':'0.04','SO4--_SS':'0.04','PM10_SO4--_SS':'0.04','PM2.5_SO4--_SS':'0.04'}, 'documented_flow_rate':'16.67', 'instrument_manual_name':'Metrohm_MARGA_Specs.pdf'},
        'Thermo ICS-1100':           {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          'instrument_manual_name':'Thermo_ICS-1100_Manual.pdf'},
        'Thermo ICS-2000':           {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          'instrument_manual_name':'Thermo_ICS-2000_Manual.pdf'},
        'Thermo ICS-2100':           {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          'instrument_manual_name':'Thermo_ICS-2100_Manual.pdf'},
        'Thermo ICS-3000':           {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          'instrument_manual_name':'Thermo_ICS-3000_Manual.pdf'},
        'Thermo ICS-5000':           {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          'instrument_manual_name':'Thermo_ICS-5000_Manual.pdf'}
    }
},

#---------------------
#titration (TI)

#Titration is a common laboratory method of quantitative chemical analysis that is used to determine the concentration of an identified analyte. 
#Since volume measurements play a key role in titration, it is also known as volumetric analysis. 
#A reagent, called the titrant or titrator is prepared as a standard solution. 
#A known concentration and volume of titrant reacts with a solution of analyte or titrand to determine concentration. 
#The volume of titrant reacted is called titration volume.
#Titration can be also done automatically, in a continuous process.

#Typical measurement operation:
#1. A typical titration begins with a beaker or Erlenmeyer flask containing a very precise volume of the analyte and a small amount of indicator (such as phenolphthalein) placed underneath a calibrated burette or chemistry pipetting syringe containing the titrant. 
#2. Small volumes of the titrant are then added to the analyte and indicator until the indicator changes color in reaction to the titrant saturation threshold, reflecting arrival at the endpoint of the titration. 
#3. Depending on the endpoint desired, single drops or less than a single drop of the titrant can make the difference between a permanent and temporary change in the indicator.    
#4. When the endpoint of the reaction is reached, the volume of reactant consumed is measured and used to calculate the concentration of analyte by: Ca = Ct*Vt*M / Va  ;where Ca is the concentration of the analyte, typically in molarity; Ct is the concentration of the titrant, typically in molarity; Vt is the volume of the titrant used, typically in liters; M is the mole ratio of the analyte and reactant from the balanced chemical equation; and Va is the volume of the analyte used, typically in liters.

#Unless additional information is given, make assumption titration uses instrumental injection of sample species, following chemical reaction with a reagent.

#Not QA accepted measurement technique for SO2 due to dominance of other more accurate techniques in measuring this species

#See more details about method for SO2:https://www3.epa.gov/ttnemc01/ctm/ctm013B.pdf

'titration':{'abbreviation':'TI', 'sampling_type':'injection', 'sample_preparation_type':'reagent reaction', 'measured_parameters':['SO2'], 'qa_accepted_parameters':[], 
    'instruments':{}
},

##------------------------------------
#aerosol mass spectrometry (AMS)

#This method is a variation of the Mass Spectrometry method to analyse gas-phase species, specific for aerosols.
#In Aerosol Mass Spectrometry, aerosol particles are – together with a small fraction (10−7) of a carrier gas – directed onto the tungsten vaporizer, typically heated to a temperature of 550–600 ◦C.
#Non-refractory material flash-vaporises, and the emerging vapor is electron-impact-ionised (70 eV) for subsequent analysis with a quadrupole or a time-of-flight mass spectrometer.
#In addition to the measurement of the aerosol beam (“beam-open” mass spectrum), also the instrument background is measured with the aerosol beam blocked (“beam-closed” mass spectrum). 
#The particle contribution is calculated from the difference of both measurements (“difference” spectrum), assuming that all particle components vaporize quickly compared to the beam-open–beam-closed cycle length. 
#For calculation of mass concentrations of species from the difference spectrum, each individual m/z is associated with one or several of these species.
#For the standard analysis of the unit mass resolution spectra, these associations are listed in the “frag table”, which needs to be adapted for special measurement situations by the user. 
#Measurements with the high-resolution time-of-flight AMS (HR-ToF-AMS) allow the quantification of individual signals in the spectra with certain elemental composition.
#For measurements of the continental or urban background aerosol the assumptions behind this procedure are typically well met and the AMS provides robust quantitative information on the sub-micron aerosol composition. 
#However, under certain conditions, e.g., when measuring close to anthropogenic sources or in the marine environment, the limitations of these assumptions are sometimes reached and the standard analysis could result in misinterpretation of the mass spectra.
#The method can be performed online (sampling ambient air) or offline (through injection). 
#It is assumed unless otherwise stated that instrumental sampling is done via injection.

'aerosol mass spectrometry':{'abbreviation':'AMS', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'qa_accepted_parameters':['Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 
    'instruments':{
        'Aerodyne CTOF AMS':       {'measured_parameters':['PM2.5_Cl-','PM1_Cl-','PM2.5_NO3-','PM1_NO3-','PM2.5_NH4+','PM1_NH4+','PM2.5_SO4--','PM1_SO4--'], 'qa_accepted_parameters':['PM2.5_Cl-','PM1_Cl-','PM2.5_NO3-','PM1_NO3-','PM2.5_NH4+','PM1_NH4+','PM2.5_SO4--','PM1_SO4--','PM2.5_SO4--_NSS','PM1_SO4--_NSS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'documented_lower_limit_of_detection':{'PM2.5_Cl-':np.NaN,'PM1_Cl-':np.NaN,'PM2.5_NO3-':'0.0012','PM1_NO3-':'0.0012','PM2.5_NH4+':'0.024','PM1_NH4+':'0.024','PM2.5_SO4--':'0.0024','PM1_SO4--':'0.0024','PM2.5_SO4--_NSS':'0.0024','PM1_SO4--_NSS':'0.0024','PM2.5_SO4--_SS':'0.0024','PM1_SO4--_SS':'0.0024'}, 'documented_flow_rate':'0.85',  'instrument_manual_name':'Aerodyne_AMS_Specs.pdf'},
        'Aerodyne HRTOF AMS':      {'measured_parameters':['PM2.5_Cl-','PM1_Cl-','PM2.5_NO3-','PM1_NO3-','PM2.5_NH4+','PM1_NH4+','PM2.5_SO4--','PM1_SO4--'], 'qa_accepted_parameters':['PM2.5_Cl-','PM1_Cl-','PM2.5_NO3-','PM1_NO3-','PM2.5_NH4+','PM1_NH4+','PM2.5_SO4--','PM1_SO4--','PM2.5_SO4--_NSS','PM1_SO4--_NSS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'documented_lower_limit_of_detection':{'PM2.5_Cl-':np.NaN,'PM1_Cl-':np.NaN,'PM2.5_NO3-':'0.0029','PM1_NO3-':'0.0029','PM2.5_NH4+':'0.058','PM1_NH4+':'0.058','PM2.5_SO4--':'0.0058','PM1_SO4--':'0.0058','PM2.5_SO4--_NSS':'0.0058','PM1_SO4--_NSS':'0.0058','PM2.5_SO4--_SS':'0.0058','PM1_SO4--_SS':'0.0058'}, 'documented_flow_rate':'0.85',  'instrument_manual_name':'Aerodyne_AMS_Specs.pdf'},
        'Aerodyne quadrupole ACSM':{'measured_parameters':['PM2.5_Cl-','PM1_Cl-','PM2.5_NO3-','PM1_NO3-','PM2.5_NH4+','PM1_NH4+','PM2.5_SO4--','PM1_SO4--'], 'qa_accepted_parameters':['PM2.5_Cl-','PM1_Cl-','PM2.5_NO3-','PM1_NO3-','PM2.5_NH4+','PM1_NH4+','PM2.5_SO4--','PM1_SO4--','PM2.5_SO4--_NSS','PM1_SO4--_NSS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'documented_lower_limit_of_detection':{'PM2.5_Cl-':'0.02','PM1_Cl-':'0.02','PM2.5_NO3-':'0.02','PM1_NO3-':'0.02','PM2.5_NH4+':'0.5','PM1_NH4+':'0.5','PM2.5_SO4--':'0.04','PM1_SO4--':'0.04','PM2.5_SO4--_NSS':'0.04','PM1_SO4--_NSS':'0.04','PM2.5_SO4--_SS':'0.04','PM1_SO4--_SS':'0.04'},                     'documented_flow_rate':'0.085', 'instrument_manual_name':'Aerodyne_quadrupole_ACSM_Specs.pdf'},
        'Aerodyne eTOF ACSM':      {'measured_parameters':['PM2.5_Cl-','PM1_Cl-','PM2.5_NO3-','PM1_NO3-','PM2.5_NH4+','PM1_NH4+','PM2.5_SO4--','PM1_SO4--'], 'qa_accepted_parameters':['PM2.5_Cl-','PM1_Cl-','PM2.5_NO3-','PM1_NO3-','PM2.5_NH4+','PM1_NH4+','PM2.5_SO4--','PM1_SO4--','PM2.5_SO4--_NSS','PM1_SO4--_NSS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'documented_lower_limit_of_detection':{'PM2.5_Cl-':'0.003','PM1_Cl-':'0.003','PM2.5_NO3-':'0.007','PM1_NO3-':'0.007','PM2.5_NH4+':'0.06','PM1_NH4+':'0.06','PM2.5_SO4--':'0.006','PM1_SO4--':'0.006','PM2.5_SO4--_NSS':'0.006','PM1_SO4--_NSS':'0.006','PM2.5_SO4--_SS':'0.006','PM1_SO4--_SS':'0.006'},         'documented_flow_rate':'0.085', 'instrument_manual_name':'Aerodyne_eTOF_ACSM_Specs.pdf'},    
    }
},
#------------------------------------
#gravimetry (GR)

#An air pump draws ambient air at a constant flow rate into a specially shaped inlet where particulate matter is separated into size fractions. 
#Particulate matter is then collected on a filter. Each filter is weighed before and after use, to determine the net mass gain due to collected matter. 
#The total volume of air filtered is known from the constant air flow, and the difference in filter weights is used to calculate the particulate matter concentration in micrograms per cubic meter (μg/m3) of air.

#Known Interference:
#Particulate matter may be lost during filter handling and weighing procedures, especially if filter is exposed to warming.
#Gaseous species may contaminate filters.
#Humidity and absorbed water may be difficult to control both during operations and when handling filters.
#Removing filters and transporting to a lab for analysis may affect results. 
#Meteorological conditions may affect flow rate.

'gravimetry':{'abbreviation':'GR', 'sampling_type':'manual', 'measured_parameters':['PM10','PM2.5','PM1'], 'qa_accepted_parameters':['PM10','PM2.5','PM1'], 
    'instruments':{}
},

#------------------------------------
#tapered element oscillating microbalance - gravimetry (TEOM-GR)

#The TEOM - Gravimetric method essentially involves a  true “gravimetric” instrument that draws (then heats) ambient air through a filter at constant flow rate, continuously weighing the filter and calculating near real-time mass concentrations of particulate matter
#The weighing principle used in the tapered element oscillating microbalance (TEOM) TEOM mass transducer is similar to that of a laboratory microbalance in that the mass change detected by the sensor is the result of the measurement of a change in a parameter (in this case, frequency) that is directly coupled via a physical law (or from first principles).

#Typical operation:
#1. Air is drawn through a tapered glass element with a filter attached. 
#2. The element oscillates according to a characteristic frequency, that decreases as mass accumulates on the attached filter. 
#3. Measurement of the change in frequency converts to measurement of the accumulated mass.

'tapered element oscillating microbalance - gravimetry':{'abbreviation':'TEOM-GR', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['PM10','PM2.5','PM1'], 'qa_accepted_parameters':['PM10','PM2.5','PM1'], 
    'instruments':{
        'Thermo 1400':  {                                                                                                                                                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses tapered element oscillating microbalance - gravimetry here: https://uk-air.defra.gov.uk/assets/documents/reports/empire/lsoman/LSO_manual_2012_18_R&P_TEOM1400.pdf'},
        'Thermo 1400a': {                                                                                                                                                                                                                                                                              'documented_flow_rate':'3.0',     'instrument_manual_name':'Thermo_1400a_Manual.pdf'},
        'Thermo 1400ab':{                                                                                   'documented_lower_limit_of_detection':'0.06', 'documented_upper_limit_of_detection':'5000000', 'documented_precision':'1.5', 'documented_accuracy':'0.75%', 'documented_resolution':'0.1', 'documented_flow_rate':'0.5-4.0', 'instrument_manual_name':'Thermo_1400ab_Specs.pdf'},
        'Thermo 1405':  {                                                                                                                                 'documented_upper_limit_of_detection':'1000000', 'documented_precision':'2.0', 'documented_accuracy':'0.75%', 'documented_resolution':'0.1', 'documented_flow_rate':'3.0',     'instrument_manual_name':'Thermo_1405_Specs.pdf'},
        'Thermo 1405D': {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'],                                                                                                                                                                                            'documented_flow_rate':'3.0',     'instrument_manual_name':'Thermo_1405D_Manual.pdf'}
    }
},
#------------------------------------
#tapered element oscillating microbalance - filter dynamics measurement system - gravimetry (TEOM-FDMS-GR)

#The Filter Dynamics Measurement System (FDMS) is a hybrid of the TEOM system which is able to measure the semi-volatile fraction of airborne particulate matter. 
#This can be important in the study of primary and secondary PM and to ensure there are no losses of these semi volatile fractions due to sampling conditions. 

'tapered element oscillating microbalance - filter dynamics measurement system - gravimetry':{'abbreviation':'TEOM-FDMS-GR', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['PM10','PM2.5','PM1'], 'qa_accepted_parameters':['PM10','PM2.5','PM1'], 
    'instruments':{
        'Thermo 1405-DF':{'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'], 'documented_upper_limit_of_detection':'1000000', 'documented_precision':'2.0', 'documented_resolution':'0.1', 'documented_flow_rate':'1.67-12.0', 'instrument_manual_name':'Thermo_1405-DF_Specs.pdf'},
        'Thermo 1405F':  {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'],                                                                                                               'documented_flow_rate':'3.0',       'instrument_manual_name':'Thermo_1405F_Manual.pdf'}
    }
},

#------------------------------------
#beta-attenuation (BA)

#Beta particles (electrons with energies in the 0.01 to 0.1 MeV range) are attenuated according to an approximate exponential function when they pass through particulate deposits on a filter tape. 
#Automated samplers (analyzers) use a continuous filter tape, first measuring the attenuation by the unexposed tape, and then measuring the attenuation after the tape has passed through the ambient air flow. 
#The attenuation measurement converts to a measure of the mass on the filter, so that the filters do not require later laboratory analysis for the mass variable. For some devices, the beta particle source is 14C.

#Known Interference:
#Particulate matter may be lost due to filter tape advance and vibration, especially if filter is exposed to warming.
#Gaseous species may contaminate filters.
#Humidity and absorbed water may be difficult to control during operations. 
#Meteorological conditions may affect flow rate.
#Although on-site real-time mass measurement offers significant improvements over the filter removal and laboratory analysis process, the beta emission and detection process present additional on-site maintenance requirements.

'beta-attenuation':{'abbreviation':'BA', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['PM10','PM2.5','PM1'], 'qa_accepted_parameters':['PM10','PM2.5','PM1'], 
    'instruments':{
        'Andersen FH62 I-N':          {                                                                                                                                                                                                                             'documented_accuracy':'2.0',                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://www3.epa.gov/ttn/amtic/files/ambient/inorganic/mthd-1-1.pdf'},
        'Andersen FH621-N':           {                                                                                   'documented_lower_limit_of_detection':'2.0', 'documented_upper_limit_of_detection':'10000',                                               'documented_accuracy':'2.0',                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://cdiac.ess-dive.lbl.gov/programs/NARSTO/Compendium_combined.pdf; see details about detection limits here: https://books.google.es/books?id=gDJEBAAAQBAJ&pg=PT323&lpg=PT323&dq=andersen+fh621-n&source=bl&ots=fU7hu5zKUU&sig=ACfU3U1TC8JhHedeoUMYXqfHdPIfrrN0rw&hl=en&sa=X&ved=2ahUKEwi--YrC_oPgAhVi5OAKHaf2CEcQ6AEwAnoECAcQAQ#v=onepage&q=andersen%20fh621-n&f=false'},
        'Dasibi 7001':                {                                                                                                                                                                                                                                                                                                                                                                                                                                         'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: http://www.sistemapiemonte.it/ambiente/srqa/stazioni/pdf/151.pdf'},
        'Eberline FH62-1':            {                                                                                                                                                                                                                             'documented_accuracy':'2.0',                                                                                                                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: http://adsabs.harvard.edu/abs/2012AtmEn..54...18Gf'},
        'Environnement SA MP101M':    {                                                                                   'documented_lower_limit_of_detection':'0.5', 'documented_upper_limit_of_detection':'10000',                                                                                                                                                                                                                   'documented_flow_rate':'16.666667',     'instrument_manual_name':'Environnement_SA_MP101M_Specs.pdf'},
        'Environnement SA MPSI100':   {                                                                                                                                                                                                                                                                                                                                                                                                                                         'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://eurekamag.com/pdf/008/008897460.pdf'},
        'FAI Swam 5a':                {                                                                                   'documented_lower_limit_of_detection':'1.0', 'documented_upper_limit_of_detection':'10000',                                                                             'documented_uncertainty':'1.0',                                                                                                       'documented_flow_rate':'13.33-41.6667', 'instrument_manual_name':'FAI_Swam_5a_specs.pdf', 'measuring_instrument_further_details':'See more details here: http://fai-instruments.net/swam-dual-channel-2/?lang=en'}, 
        'Horiba APDA350E':            {                                                                                                                                                                                                                                                                                                                                                                                                                                         'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://books.google.es/books?id=wYUUzLCIMzsC&pg=PA112&lpg=PA112&dq=horiba+apda+350e&source=bl&ots=uvOvFFezdD&sig=ACfU3U041o7NZz92BUIgpBh1x0CLKwY5lQ&hl=en&sa=X&ved=2ahUKEwjfo6_x4ojgAhUFQhoKHeuAAlgQ6AEwCHoECAEQAQ#v=onepage&q=horiba%20apda%20350e&f=false'},
        'Horiba APDA351E':            {                                                                                                                                                                                                                                                                                                                                                                                                                                         'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://pdfs.semanticscholar.org/fc4c/3b522b04166381468138e590ac4b866759d3.pdf'},
        'Horiba APDA371':             {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'], 'documented_lower_limit_of_detection':'4.8', 'documented_upper_limit_of_detection':'10000',                                                                                                             'documented_resolution':'0.1',                                                                        'documented_flow_rate':'16.666667',     'instrument_manual_name':'Horiba_APDA371_Specs.pdf'},
        'Horiba FH621-R':             {                                                                                                                                                                                                                                                                                                                                                                                                                                         'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://scindeks-clanci.ceon.rs/data/pdf/1451-9372/2008/1451-93720801005J.pdf'},
        'Kimoto 168S':                {                                                                                                                                                                                                                                                                                                                                                                                                                                         'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://www.dri.edu/images/stories/editors/eafeditor/Watsonetal2010SCAQMDFugDustReport.pdf'},
        'Kimoto SPM-613D':            {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'],                                                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://www3.epa.gov/ttnamti1/files/2006conference/tisch.pdf'},
        'Kimoto SPM-613F':            {                                                                                   'documented_lower_limit_of_detection':'2.0',  'documented_upper_limit_of_detection':'1000',  'documented_precision':'2.0%',               'documented_accuracy':'3.0%',                                                                'documented_zero_drift':'3.0/day', 'documented_span_drift':'2.0%/day', 'documented_flow_rate':'16.666667',     'instrument_manual_name':'Kimoto_SPM-613-F_Manual.pdf'},
        'Met One BAM1020':            {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'], 'documented_lower_limit_of_detection':'4.8',  'documented_upper_limit_of_detection':'10000',                                                                            'documented_uncertainty':'1.0', 'documented_resolution':'0.24',                                                                       'documented_flow_rate':'16.7',          'instrument_manual_name':'Met_One_BAM1020_Manual.pdf'}, 
        'OPSIS SM200':                {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'], 'documented_lower_limit_of_detection':'0.5',  'documented_upper_limit_of_detection':'1000',                                                                                                                                                                                                                   'documented_flow_rate':'16.666667',     'instrument_manual_name':'OPSIS_SM200_Specs'},
        'Thermo 5014i':               {                                                                                   'documented_lower_limit_of_detection':'4.0',  'documented_upper_limit_of_detection':'10000', 'documented_precision':'2.0<80.0 5.0>=80.0', 'documented_accuracy':'5.0%',                                 'documented_resolution':'0.1',                                                                        'documented_flow_rate':'16.666667',     'instrument_manual_name':'Thermo_5014i_Manual.pdf'},
        'Thermo FH62 C14':            {                                                                                   'documented_lower_limit_of_detection':'4.0',  'documented_upper_limit_of_detection':'10000', 'documented_precision':'2.0',                                                                              'documented_resolution':'1.0',                                                                        'documented_flow_rate':'0.6',           'instrument_manual_name':'Thermo_FH62_C14_Specs.pdf'},
        'Thermo Andersen FH62 I-R':   {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'], 'documented_lower_limit_of_detection':'3.0',  'documented_upper_limit_of_detection':'10000',                                              'documented_accuracy':'1.0',                                  'documented_resolution':'1.0',                                                                        'documented_flow_rate':'16.666667',     'instrument_further_details':'Instrument was called Eberline FH62 I-R previously, specifications for that instrument are here: https://www.breckland.gov.uk/media/1176/The-Detailed-Assessment-for-PM10-2004/pdf/The_Detailed_Assessment_for_PM10_2004.pdf'},
        'Verewa F701':                {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'], 'documented_lower_limit_of_detection':'5.0',  'documented_upper_limit_of_detection':'10000',                                                                                                                                                                                                                  'documented_flow_rate':'16.666667',     'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: http://www.ecomonitoring.pl/pylomierz%20-%20imisja.html'},
        'Verewa F701-20':             {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'], 'documented_lower_limit_of_detection':'1.0',  'documented_upper_limit_of_detection':'10000',                                              'documented_accuracy':'2.0%',                                                                                                                                       'documented_flow_rate':'16.666667',     'instrument_manual_name':'Verewa_F701-20_Manual.pdf'},
        'Verewa F703':                {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'],                                                                                                                                                                                                                                                                                                                                                       'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://books.google.es/books?id=zx6YBwAAQBAJ&pg=SL13-PA29&lpg=SL13-PA29&dq=Verewa+F703&source=bl&ots=zS_cDXxd_x&sig=ACfU3U2QDSpJnu-ns_aORl3mN9x_SrGYeQ&hl=en&sa=X&ved=2ahUKEwjuhNf38IbgAhXxoXEKHeP6BDoQ6AEwC3oECAQQAQ#v=onepage&q=Verewa%20F703&f=false'},
        'Wedding and Associates Beta':{                                                                                                                                                                                                                                                                                                                                                                                                                                         'instrument_further_details':'No manual or specifications available but confirmation it is uses beta-attenuation here: https://uk-air.defra.gov.uk/assets/documents/reports/empire/quarg/quarg_11.pdf'}
    }
},

#------------------------------------
#impaction (IM)

#Impaction is a methodology employed for the measurement of particulate matter where some information on the particle size distribution is desired
#An impactor has two co-linear plates; one acts as a collection surface for the particles, the other one has small nozzle or nozzles in it to control the flow velocity. 
#The sample flow is first led through the nozzles to achieve a certain, exact flow velocity. 
#After the sample passes through the nozzles it is turned sharply in front of the collection plate in which case particles larger than stage cut diameter cannot follow the flow stream lines but are impacted onto the collection plate. 
#Particles smaller that the stage cut diameter continue to the following impactor stages where they are further size classified and collected. 
#By changing the dimensions in each impactor stage, different sized particles can be collected on different impactor stages. Typically there are stage cut points of 10, 2.5 and 1.0 µm.

#Before the measurement a collection substrate is placed on each of the impactor stages and this substrate is weighed before and after the measurement to determine particle size distribution. 
#Depending on the used substrate material, a chemical composition of the collected particles may also be determined after the sample collection.
#Typically PM10 impactor particles >1.0 µm are collected on 25 mm substrates and particles <1.0 µm are collected on a 47 mm filter. 
#Different material substrates and filters can be used in the impactor depending on the used particle analysis method.

'impaction':{'abbreviation':'IM', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['PM10','PM2.5','PM1'], 'qa_accepted_parameters':['PM10','PM2.5','PM1'], 
    'instruments':{
        'Dekati PM10 Impactor':{'documented_flow_rate':'10.0-30.0', 'instrument_manual_name':'Dekati_PM10_Impactor_Specs.pdf'}
    }
},
#------------------------------------
#nephelometry (NP)

#Nephelometry is based on the bulk scattering of visible light by an air sample. It is a mature technology and provides reliable measurements of the spectral (aerosol) scattering coefficient of the air (unit m‐1).
#An integrating nephelometer combines the light from all scattering angles to produce a spectral scattering coefficient for the wavelength of light used.
#The signal obtained by integrating nephelometer measurements is closely proportional to the cosine weighted scattering function, which is equal to the scattering part of the Lambert‐Beer extinction coefficient. 
#Thus the primary metric measured is the particle scattering coefficient that is usually given in reciprocal megameters (Mm‐1). 
#This ability to provide one important basic physical characteristic of the aerosol makes integrating nephelometers particularly valuable in climate change research.

#This method uses the light scattered by tiny particles (nephelometry) to determine the particulate matter concentration in the ambient air directly and continuously.
#A highly sensitive scattered light sensor lies at the heart of the applied measuring method. The light emitted by an intensity stabilized laser diode illuminates a measuring space defined by the optical path. 
#The light scattered by all the aerosol particles inside this measuring space is captured by a semiconductor photodetector positioned at an angle of 90°. 
#After amplification, the outcome of this measurement is made available as a voltage signal.
#The signal is directly proportional to the mass concentration of the aerosol in the measuring space.

'nephelometry':{'abbreviation':'NP', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['PM10','PM2.5','PM1'], 'qa_accepted_parameters':['PM10','PM2.5','PM1'], 
    'instruments':{}
},
#------------------------------------
#light scattering photometry (LSP)

#Light scattering photometers do not integrate the light from all scattered angles, but only view a portion of the available scattering angles. 
#Like nephelometry, it produces a bulk scattering coefficient of an air sample with similar time resolution and range. In general they are calibrated with a standard test aerosol such as Arizona Test Dust, but typically these test aerosols have significantly different optical properties from real ambient air aerosols.
#Again, therefore, assumptions about the type of aerosol present, or frequent calibration by other methods, is needed to convert the scattering coefficient to a mass concentration, and uncertainties tend to be higher than for nephelometers.
#See here for distinction between nephelometry and light scattering photometry: https://db-airmontech.jrc.ec.europa.eu/download/PM_Mass_MMTI_OpticalDevices.pdf

'light scattering photometry':{'abbreviation':'LSP', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['PM10','PM2.5','PM1'], 'qa_accepted_parameters':['PM10','PM2.5','PM1'], 
    'instruments':{
        'Comde Derenda APM-2':{'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'],                                                                                                                                                                                                                                                                         'instrument_manual_name':'Comde_Derenda_APM-2_Specs.pdf'},
        'Turnkey Osiris':     {                                                                                   'documented_lower_limit_of_detection':'0.01', 'documented_upper_limit_of_detection':'6000.0',                                                                                                        'documented_resolution':'0.1', 'documented_flow_rate':'0.6',       'instrument_manual_name':'Turnkey_Osiris_Specs.pdf'},
        'Unitec LSPM10':      {'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'],                                               'documented_upper_limit_of_detection':'20000.0', 'documented_accuracy':'1.0%', 'documented_zero_drift':'0.1%/day', 'documented_span_drift':'0.1%/day', 'documented_resolution':'0.1', 'documented_flow_rate':'15.0-40.0', 'instrument_manual_name':'Unitec_LSPM10_Manual.pdf'}
    }
},

#------------------------------------
#optical particle counter (OPC)

#Optical particle counters differ from nephelometry and light scattering photometry in being based on measuring the light scattered by individual particles. 
#The light scattered by individual particles is picked up as (typically) a discrete signal, so that the size of each particle can be gauged from the intensity of the signal, while the number concentration of particles can be gauged from the number of signals and the air flow rate through the light beam volume. 
#In general, therefore, a size spectrum of particles will be obtained, so that information about specific size fractions, such as PM10 and PM2.5, can be obtained simultaneously. 
#These can then be converted to mass concentrations by making assumptions about the particle density.
#Uncertainties arise from the determination of the size of the volume of air being examined by the laser, scattering from more than one particle at the same time (coincidence), and the conversions from intensity to size and then to mass. 
#Nevertheless, performance can be comparable with other automatic methods for PM10 and PM2.5.

'optical particle counter':{'abbreviation':'OPC', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['PM10','PM2.5','PM1'], 'qa_accepted_parameters':['PM10','PM2.5','PM1'], 
    'instruments':{
        'Grimm 107':       {                                                 'documented_flow_rate':'1.2',  'instrument_manual_name':'Grimm_107_Specs.pdf', 'measuring_instrument_further_details':'Portable version of Grimm 180.'},
        'Grimm 180':       {                                                 'documented_flow_rate':'1.2',  'instrument_manual_name':'Grimm_180_Manual.pdf'},
        'Grimm 190':       {                                                                                'instrument_further_details':'No manual or specifications available but confirmation it is uses optical particle counter here: https://core.ac.uk/download/pdf/52132153.pdf'},
        'Grimm 365':       {                                                 'documented_flow_rate':'1.2',  'instrument_manual_name':'Grimm_365_Specs.pdf'},
        'Palas Fidas 200': {'documented_upper_limit_of_detection':'10000.0', 'documented_flow_rate':'5.0',  'instrument_manual_name':'Palas_Fidas_200_Specs.pdf'},
        'Palas Fidas 200E':{'documented_upper_limit_of_detection':'10000.0', 'documented_flow_rate':'5.0',  'instrument_manual_name':'Palas_Fidas_200_Specs.pdf'},
        'Palas Fidas 200S':{'documented_upper_limit_of_detection':'10000.0', 'documented_flow_rate':'5.0',  'instrument_manual_name':'Palas_Fidas_200_Specs.pdf'}
    }
},

#------------------------------------
#beta-attenuation - nephelometry (BA-NP)

#Hybrid of beta-attenuation and nephelometry methods enhancing the overall performance.
#Has accuracy of the beta-attenuation method with high time resolution of a nephelometer.
#see https://www3.epa.gov/ttnamti1/files/2006conference/goohssharp.pdf for more information

'beta-attenuation - nephelometry':{'abbreviation':'BA-NP', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['PM10','PM2.5','PM1'], 'qa_accepted_parameters':['PM10','PM2.5','PM1'], 
    'instruments':{
        'Thermo 5030 SHARP':{'measured_parameters':['PM10','PM2.5'], 'qa_accepted_parameters':['PM10','PM2.5'], 'documented_lower_limit_of_detection':'0.5', 'documented_upper_limit_of_detection':'10000', 'documented_precision':'2.0<80.0 5.0>=80.0', 'documented_accuracy':'5.0%', 'documented_span_drift':'0.02%/day', 'documented_flow_rate':'16.666667', 'instrument_manual_name':'Thermo_5030_SHARP_Specs.pdf'}
    }
},

#------------------------------------
#thermal-optical analysis (TOA)

#Thermal-optical analysis (TOA) has been widely used to separate carbonaceous aerosols from ambient and source samples into two components, organic and elemental carbon. 
#This method uses volatility to separate groups of carbon, and laser monitoring to correct for the transformation of non-absorbing carbon into pyrolytic carbon that absorbs light.
#Thermal analysis has been used to classify carbonaceous aerosols in both ambient and source samples into light-absorbing carbon (LAC) and organic carbon (OC).
#Because thermal analysis does not yield a true value of LAC, the refractory component it measures takes a different name: elemental carbon (EC). (This has no relationship to a chemist’s definition of “elemental.”) 
#Thermal methods rely on the assumption that organic carbon is volatile, and the strongly light absorbing carbon is stable or refractory at elevated temperatures. 
#It assumes that LAC in a sample is unchanged until exposed to both elevated temperatures and an oxidant, so that carbon driven off at low temperatures or in the absence of oxygen must be OC. 
#Ideally, TOA would be a separation-and-detection technique similar to chromatography, but failures of separation have hampered its interpretation.
#Some OC undergoes pyrolysis or charring in the absence of oxygen, forming pyrolyzed carbon (“charring”) instead of volatilizing (Huntzicker et al., 1982). 
#For that reason, the sample is also monitored optically. The sample usually becomes blacker as charring occurs in an inert atmosphere, and the blackness decreases as the pyrolyzed carbon burns later in the analysis. 
#This change is detected by monitoring either transmittance or reflectance (sometimes both). When the filter’s transmittance or reflectance returns to the original value, it is assumed that the remaining material is EC. 
#This point is colloquially termed the OC/EC split and the combined analysis is known as thermal-optical analysis (TOA).
#Investigations of the thermal-optical method have examined “round-robin” tests of different protocols, or sequences of temperature and analysis environments (e.g., Cadle et al., 1990; Countess et al., 1990; Schmid et al., 2001). 
#Studies have provided the differences between OC/EC split produced by various temperature programs (e.g., Chow et al., 2001; Subrama- nian et al., 2006), and quantification of positive and negative artifacts (e.g., Cadle et al., 1983; Huebert and Charlson, 1996;Kirchstetter et al., 2001). 
#Some examples of successful results have been agreement on measurement of total carbon (Schmid et al., 2001) and reproducible measurements of EC and OC from careful and consistent application of thermal programs (Schauer et al., 2003). 
#A wide variety of temperature programs has been used in TOA, often with differing results (Schmid et al., 2001; Watson et al., 2005). The combination of temperature program and optical analysis best suited to analyzing different types of samples is still a subject of debate. 
#Many papers have pointed out challenges in interpreting TOA results (Cadle et al., 1983; Reid et al., 1998; Chow et al., 2001; Conny et al., 2003; Subramanian et al., 2006), and the method has both detractors and supporters.

#See more details on method here: https://www.tandfonline.com/doi/pdf/10.1080/02786820802360690?needAccess=true

'thermal-optical analysis':{'abbreviation':'TOA', 'sampling_type':'injection', 'measured_parameters':['C','PM10_C','PM2.5_C','PM1_C','EC','PM10_EC','PM2.5_EC','PM1_EC','OC','PM10_OC','PM2.5_OC','PM1_OC'], 'qa_accepted_parameters':['C','PM10_C','PM2.5_C','PM1_C','EC','PM10_EC','PM2.5_EC','PM1_EC','OC','PM10_OC','PM2.5_OC','PM1_OC'], 
    'instruments':{
        'DRI Model 2001':     {'instrument_manual_name':'DRI_Model_2001_Manual.pdf'},
        'Sunset OCEC Model 4':{'instrument_manual_name':'Sunset_OCEC_4_Specs.pdf'},
        'Sunset OCEC Model 5':{'instrument_manual_name':'Sunset_OCEC_5l_Manual.pdf', 'instrument_further_details':'See more details here:http://www.sunlab.com/wp-content/uploads/Lab-Instrument-brochure.pdf'}
    }
},

#------------------------------------
#aethalometer (ATH)

#The aethalometer is an instrument designed for the measurement of black carbon deposit collected continuously deposited on quartz fiber filter (designed by Magee Scientific).
#The 'aethalometer' is the dominant market instrument (designed by Magee Scientific), of which there are multiple versions. 
#It measures the concentration of optically absorbing (‘black’) suspended particulates in a gas colloid stream; commonly visualised as smoke or haze, often seen in ambient air under polluted conditions. 
#The word aethalometer is derived from the Classical Greek verb ‘aethaloun’, meaning ‘to blacken with soot’.
#The gas stream (frequently ambient air) passes through a filter material which traps the suspended particulates, creating a deposit of increasing density. 
#A light beam projected through the deposit is attenuated by those particles which are absorbing (‘black’) rather than scattering (‘white’).
#Measurements are made at successive regular time intervals. The increase in attenuation from one measurement to the next is proportional to the increase in the density of optically absorbing material on the filter: which, in turn, is proportional to the concentration of the material in the sampled air stream. 
#The sample is collected as a spot on a roll of filter tape. 
#When the density of the deposit spot reaches a pre-set limit, the tape advances to a fresh spot and the measurements continue. 
#Measurement of the sample gas flow rate and knowledge of the instrument’s optical and mechanical characteristics permit a calculation of the average concentration of absorbing particles in the gas stream during the sampling period. 
#Aethalometers may operate on time-base periods as rapid as 1 second, providing quasi-real-time data. Comparison of aethalometer data with other physical and chemical analyses allows the output to be expressed as a concentration of black carbon.
#More information here: https://www.tandfonline.com/doi/full/10.1080/027868290901972

'aethalometer':{'abbreviation':'ATH', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 'qa_accepted_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 
    'instruments':{
        'Magee AE16':{                                                                                                                                                                                                                                                                                                                                                                         'instrument_further_details':'No manual or specifications available but confirmation instrument is aethalometer here: https://www.arm.gov/publications/tech_reports/handbooks/aeth_handbook.pdf'},
        'Magee AE22':{'measured_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 'qa_accepted_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'],                                                                                               'documented_resolution':'0.1',  'documented_flow_rate':'2.0-6.0', 'instrument_manual_name':'Magee_AE31_Specs.pdf'},
        'Magee AE31':{'measured_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 'qa_accepted_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'],                                                                                               'documented_resolution':'0.1',  'documented_flow_rate':'2.0-6.0', 'instrument_manual_name':'Magee_AE31_Specs.pdf'},
        'Magee AE33':{'measured_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 'qa_accepted_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 'documented_lower_limit_of_detection':'0.005', 'documented_upper_limit_of_detection':'100.0', 'documented_resolution':'0.03', 'documented_flow_rate':'2.0-5.0', 'instrument_manual_name':'Magee_AE33_Specs.pdf'}
    }
},

#------------------------------------
#multi angle absorption photometer (MAAP)

#The MAAP instrument (introduced with the Thermo-Scientific model 5012) was specifically developed to reduce the uncertainties in black carbon measurements caused by aerosol scattering.
#The MAAP measures black carbon mass loadings based on aerosol optical absorption. In common with the Aethalometer this measurement is made on particles collected on a filter substrate, however there are some important differences in how the instruments operate, as well as with the algorithms used to extract the absorption coefficient and mass loading.
#With the MAAP the optical absorption coefficient of aerosol collected on a filter is determined by radiative transfer considerations which include multiple scattering effects and absorption enhancement due to reflections from the filter. 
#This calculation is based on the transmitted and reflected phase functions which are defined by directly measured values of transmission, direct and diffuse back scattering. 
#Full details of the algorithms used are given by Petzold et al (2002 and 2004). 
#Lab work showed that the diffuse back scatter component was a strong function of the fraction of scattering aerosol, being decreased (relative to the filter) for low scattering fractions and increased for high scattering fractions. 
#This lab work also informed the choice of angles for the detectors used.

#See here for more details: http://www.cas.manchester.ac.uk/restools/instruments/aerosol/maap/

#See references:
#Petzold A, Kramer H, and Schönlinner M, Continuous Measurement of Atmospheric Black Carbon Using a Multi-angle Absorption Photometer. Environmental Science and Pollution Research, Special Issue 4, 78-82, 2002.
#Petzold A and Schönlinner M, The Multi-angle absorption photometer – A new method for the measurement of aerosol light absorption and atmospheric black carbon. Journal of Aerosol Science, 35, 421-441, 2004.

'multi angle absorption photometer':{'abbreviation':'MAAP', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 'qa_accepted_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'],
    'instruments':{
        'Thermo 5012':{'documented_lower_limit_of_detection':'0.02', 'documented_upper_limit_of_detection':'180.0', 'documented_flow_rate':'16.7', 'instrument_manual_name':'Thermo_5012_Manual.pdf'}
    }
},

#------------------------------------
#particulate soot absorption photometer (PSAP)

#The PSAP instrument (manufactured by Radiance Research Inc.) measures the change in optical transmission of a filter as aerosol is deposited on the filter over time.
#The PSAP measures black carbon mass loadings based on aerosol optical absorption. In common with the Aethalometer this measurement is made on particles collected on a filter substrate, however there are some important differences in how the instruments operate, as well as with the algorithms used to extract the absorption coefficient and mass loading.
#This measurement is made at 3 wavelengths (467, 530 and 660nm), the instrument has a 10s time resolution and typically operates with a flow rate of around 1 litre/min. 
#The change in transmission from one measurement to the next is related to the optical absorption coefficient (bap m-1) of the aerosol deposited in this time. 
#In turn this can be related to black carbon mass by dividing by the specific absorption for black carbon at the wavelength of the instrument. 
#In the case of the PSAP (530nm channel) this is given as 10m2/g.
#Comparisons between the PSAP aerosol optical absorption measurement and reference techniques presented by Bond et al (1999), show that the PSAP also responds to pure scattering particles and external mixtures of scattering and absorbing particles. 
#This is because back scatter from scattering particles on the filter reduces the transmitted intensity. Bond et al show that the apparent absorption coefficient of purely scattering particles as measured by the PSAP is about 2% of their scattering coefficient. 
#Additionally even when the PSAP measurement is corrected for errors in flow rate, collection spot area and the presence of scattering particles, the same series of experiments showed that PSAP optical absorption coefficients were about 20% higher than the reference measurement for the type of aerosol used. 
#Thus PSAP data must be carefully calibrated and corrected based on measurements of aerosol scattering coefficient.

#See reference: Bond TC, Anderson TL, Campbell D, Calibration and intercomparison of filter-based measurements of visible light absorption by aerosols. Aerosol Science and Technology 30, 582–600, 1999

'particulate soot absorption photometer':{'abbreviation':'PSAP', 'sampling_type':'low volume continuous', 'sample_preparation_type':'filter', 'measured_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'], 'qa_accepted_parameters':['BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC'],
    'instruments':{
        'Radiance PSAP':{'documented_flow_rate':'1.0', 'instrument_manual_name':'Radiance_PSAP_Manual.pdf'}
    }
},

#------------------------------------
#flame atomic absorption spectroscopy (F-AAS)

#Atomic absorption spectroscopy (AAS) detects elements in either liquid or solid samples through the application of characteristic wavelengths of electromagnetic radiation from a light source. 
#Individual elements will absorb wavelengths differently, and these absorbances are measured against standards. In effect, AAS takes advantage of the different radiation wavelengths that are absorbed by different atoms.
#In AAS, analytes are first atomized so that their characteristic wavelengths are emitted and recorded. Then, during excitation, electrons move up one energy level in their respective atoms when those atoms absorb a specific energy.
#As electrons return to their original energy state, they emit energy in the form of light. 
#This light has a wavelength that is characteristic of the element. Depending on the light wavelenth and its intensity, specific elements can be detected and their concentrations measured.
#Applying the Beer-Lambert law directly in AA spectroscopy is difficult due to variations in the atomization efficiency from the sample matrix, and nonuniformity of concentration and path length of analyte atoms.
#Flame atomic absorption spectrometry (F-AAS) is where the sample is first atomised using flames, principally a air-acetylene flame with a temperature of about 2300 °C and a nitrous oxide system (N2O)-acetylene flame with a temperature of about 2700 °C. 
#The latter flame, in addition, offers a more reducing environment, being ideally suited for analytes with high affinity to oxygen.
#Liquid or dissolved samples are typically used with flame atomizers. 
#The sample solution is aspirated by a pneumatic analytical nebulizer, transformed into an aerosol, which is introduced into a spray chamber, where it is mixed with the flame gases and conditioned in a way that only the finest aerosol droplets (< 10 μm) enter the flame. 
#This conditioning process is responsible that only about 5% of the aspirated sample solution reaches the flame, but it also guarantees a relatively high freedom from interference.
#See here for more information: https://sites.chem.colostate.edu/diverdi/C433/experiments/iron%20analysis/references/aas%20methods%20-%20agilent.pdf

'flame atomic absorption spectroscopy':{'abbreviation':'F-AAS', 'sampling_type':'injection', 'measured_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{
        'GBC 932':                  {'instrument_further_details':'No manual or specifications available but confirmation it is uses flame atomic absorption spectroscopy here: http://www.speciation.net/Database/Instruments/GBC-Scientific-Equipment-Ltd/932-Plus--Atomic-Absorption-Spectrometer-;i3193'},
        'Perkin Elmer AAnalyst 600':{'instrument_manual_name':'Perkin_Elmer_AAnalyst_600_Specs.pdf'},
        'Perkin Elmer AAnalyst 700':{'instrument_manual_name':'Perkin_Elmer_AAnalyst_700_Specs.pdf'},
        'Perkin Elmer AAnalyst 800':{'instrument_manual_name':'Perkin_Elmer_AAnalyst_800_Specs.pdf'},
        'Varian 240 AA':            {'instrument_manual_name':'Varian_240_AA_Specs.pdf'}
    }
},

#------------------------------------
#graphite furnace atomic absorption spectroscopy (GF-AAS)

#Graphite furnace atomic absorption spectroscopy (GF-AAS) (also known as Electrothermal Atomic Absorption Spectroscopy (ET-AAS)) is a type of spectrometry that uses a graphite-coated furnace to vaporize the sample. 
#In GF-AAS, samples are deposited in a small graphite or pyrolytic carbon coated graphite tube, which can then be heated to vaporize and atomize the analyte. 
#The main advantages of the graphite furnace comparing to flame atomic absorption are the following:
#-The detection limits for the graphite furnace fall in the ppb range for most elements
#-Interference problems are minimized with the development of improved instrumentation
#-The graphite furnace can determine most elements measurable by flame atomic absorption in a wide variety of matrices.

'graphite furnace atomic absorption spectroscopy':{'abbreviation':'GF-AAS', 'sampling_type':'injection', 'measured_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{
        'Perkin Elmer AAnalyst 600': {'instrument_manual_name':'Perkin_Elmer_AAnalyst_600_Specs.pdf'},
        'Perkin Elmer AAnalyst 700': {'instrument_manual_name':'Perkin_Elmer_AAnalyst_700_Specs.pdf'},
        'Perkin Elmer AAnalyst 800': {'instrument_manual_name':'Perkin_Elmer_AAnalyst_800_Specs.pdf'},
        'Unicam Solar 989 QZ':       {'instrument_further_details':'No manual or specifications available but confirmation it is uses graphite furnace atomic absorption spectroscopy here: https://www.sedgeochem.uni-bremen.de/aas.html'},
        'Varian SpectrAA 400 GTA-96':{'instrument_further_details':'No manual or specifications available but confirmation it is uses graphite furnace atomic absorption spectroscopy here: http://www.labwrench.com/?equipment.view/equipmentNo/20158/Varian/SpectrAA-400/'},
        'Varian SpectrAA 880z':      {'instrument_further_details':'No manual or specifications available but confirmation it is uses graphite furnace atomic absorption spectroscopy here: http://www.medwow.com/med/atomic-absorption/varian/spectraa-880z/2927.model-spec'}, 
        'Varian 240 AA':{'instrument_manual_name':'Varian_240_AA_Specs.pdf'}
        }
},

#------------------------------------
#cold vapour atomic absorption spectroscopy (CV-AAS)

#Cold vapour atomic absorption spectroscopy (CV-AAS) is a method limited to the determination of mercury, due to it being the only metallic element to have a large enough vapor pressure at ambient temperature. 
#Because of this, it has an important use in determining organic mercury compounds in samples and their distribution in the environment. 
#The method initiates by converting mercury into Hg2+ by oxidation from nitric and sulfuric acids, followed by a reduction of Hg2+ with tin(II) chloride. 
#The mercury, is then swept into a long-pass absorption tube by bubbling a stream of inert gas through the reaction mixture. 
#The concentration is determined by measuring the absorbance of this gas at 253.7 nm. Detection limits for this technique are in the parts-per-billion range making it an excellent mercury detection atomization method.
#See here for determination of Cd using this method also: https://www.ncbi.nlm.nih.gov/pubmed/19071344

'cold vapour atomic absorption spectroscopy':{'abbreviation':'CV-AAS', 'sampling_type':'injection', 'measured_parameters':['Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg'], 'qa_accepted_parameters':['Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg'], 
    'instruments':{
        'Varian 240 AA':{'instrument_manual_name':'Varian_240_AA_Specs.pdf'}
    }
},

#------------------------------------
#hydride generation atomic absorption spectroscopy (HG-AAS)

#The hydride generation atomic absorption spectroscopy (HG-AAS) method is specialised for the measurement of specific elements. 
#The technique provides a means of introducing samples containing arsenic, antimony, selenium, bismuth, and lead into an atomizer in the gas phase. 
#With these elements, hydride atomization enhances detection limits by a factor of 10 to 100 compared to alternative methods. 
#Hydride generation occurs by adding an acidified aqueous solution of the sample to a 1% aqueous solution of sodium borohydride, all of which is contained in a glass vessel. 
#The volatile hydride generated by the reaction that occurs is swept into the atomization chamber by an inert gas, where it undergoes decomposition. 
#This process forms an atomized form of the analyte, which can then be measured by absorption or emission spectrometry.
#See more information here: https://www.sciencedirect.com/topics/chemistry/hydride-generation

'hydride generation atomic absorption spectroscopy':{'abbreviation':'HG-AAS', 'sampling_type':'injection', 'measured_parameters':['As','PM10_As','PM2.5_As','PM1_As','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se'], 'qa_accepted_parameters':['As','PM10_As','PM2.5_As','PM1_As','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se'], 
    'instruments':{
        'Varian 240 AA':{'instrument_manual_name':'Varian_240_AA_Specs.pdf'}
    }
},

#------------------------------------
#flame atomic emission spectroscopy (F-AES)

#Atomic emission spectroscopy (AES) is a method of chemical analysis that uses the intensity of light emitted from a flame, plasma, arc, or spark at a particular wavelength to determine the quantity of an element in a sample. 
#The wavelength of the atomic spectral line gives the identity of the element while the intensity of the emitted light is proportional to the number of atoms of the element.

#AAS and AES are similar methods but the differences are as follows:
#-Light Source: In AAS, a monochromatic light source is used to provide energy for the excitation of electrons. In the case of AES, it is a flame that is often used.
#-Atomization: In AAS, there is a separate chamber for atomization of the sample. However, in AES, atomization takes place step by step upon the introduction of the sample to the flame.
#-Principle of operation: In AAS, when monochromatic light is bombarded through the sample the atoms absorb energy, and the extent of absorption is recorded. In AES, the sample which gets atomised in the flame then absorbs the energy through the electrons which get excited. Later this energy is released upon the relaxation of the atoms and is measured by the instrument as the emitted energy.

#In flame atomic emission spectroscopy (F-AES) a sample of a material (analyte) is brought into the flame as a gas, sprayed solution, or directly inserted into the flame by use of a small loop of wire, usually platinum. 
#The heat from the flame evaporates the solvent and breaks chemical bonds to create free atoms. The thermal energy also excites the atoms into excited electronic states that subsequently emit light when they return to the ground electronic state. 
#Each element emits light at a characteristic wavelength, which is dispersed by a grating or prism and detected in the spectrometer.

'flame atomic emission spectroscopy':{'abbreviation':'F-AES', 'sampling_type':'injection', 'measured_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{}
},

#------------------------------------
#inductively coupled plasma atomic emission spectroscopy (ICP-AES)

#Inductively coupled plasma atomic emission spectroscopy (ICP-AES), also referred to as inductively coupled plasma optical emission spectrometry (ICP-OES), is an analytical technique used for the detection of chemical elements. 
#It is a type of emission spectroscopy that uses the inductively coupled plasma to produce excited atoms and ions that emit electromagnetic radiation at wavelengths characteristic of a particular element. 
#It is a flame technique with a flame temperature in a range from 6000 to 10,000 K. 
#The intensity of this emission is indicative of the concentration of the element within the sample.

'inductively coupled plasma atomic emission spectroscopy':{'abbreviation':'ICP-AES', 'sampling_type':'injection', 'measured_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{
        'Thermo IRIS Intrepid II':{'instrument_further_details':'No manual or specifications available but confirmation it is uses inductively coupled plasma atomic emission spectroscopy here: http://www.speciation.net/Database/Instruments/Thermo-Scientific/IRIS-Intrepid-II-;i123'},
        'Varian Vista-Pro':       {'instrument_manual_name':'Varian_Vista-Pro_Specs.pdf'}
    }
},

#------------------------------------
#cold vapour atomic fluorescence spectroscopy (CV-AFS)

#Cold vapour atomic fluorescence spectroscopy (CV-AFS) is a subset of the analytical technique known as atomic fluorescence spectroscopy (AFS).
#Used in the measurement of trace amounts of volatile heavy metals such as mercury, cold vapour AFS makes use of the unique characteristic of mercury that allows vapour measurement at room temperature. 
#Free mercury atoms in a carrier gas are excited by a collimated ultraviolet light source at a wavelength of 253.7 nanometres. 
#The excited atoms re-radiate their absorbed energy (fluoresce) at this same wavelength. Unlike the directional excitation source, the fluorescence is omnidirectional and may thus be detected using a photomultiplier tube or UV photodiode.
#Gold coated traps may be used to collect mercury in ambient air or other media. 
#The traps are then heated, releasing the mercury from the gold while passing argon through the cartridge. 
#This preconcentrates the mercury, increasing sensitivity, and also transfers the mercury into an inert gas.
#See more information for measurement of Hg here: https://ec.europa.eu/environment/air/pdf/pp_mercury4.pdf

'cold vapour atomic fluorescence spectroscopy':{'abbreviation':'CV-AFS', 'sampling_type':'injection', 'measured_parameters':['Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg'], 'qa_accepted_parameters':['Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg'], 
    'instruments':{}
},

#------------------------------------
#inductively coupled plasma mass spectrometry (ICP-MS)

#Inductively coupled plasma mass spectrometry (ICP-MS) is a type of mass spectrometry which is capable of detecting metals and several non-metals at concentrations as low as one part per quadrillion on non-interfered low-background isotopes. 
#This is achieved by ionizing the sample with inductively coupled plasma and then using a mass spectrometer to separate and quantify those ions.
#Compared to atomic absorption spectroscopy, ICP-MS has greater speed, precision, and sensitivity. 
#However, compared with other types of mass spectrometry, such as thermal ionization mass spectrometry (TIMS) and glow discharge mass spectrometry (GD-MS), ICP-MS introduces many interfering species: argon from the plasma, component gases of air that leak through the cone orifices, and contamination from glassware and the cones.

'inductively coupled plasma mass spectrometry':{'abbreviation':'ICP-MS', 'sampling_type':'injection', 'measured_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{
        'Agilent 7500':          {'instrument_manual_name':'Agilent_7500_Specs.pdf'},
        'Agilent 7800':          {'instrument_manual_name':'Agilent_7800_Specs.pdf'},
        'Agilent 7900':          {'instrument_manual_name':'Agilent_7900_Specs.pdf'},
        'Agilent 8900':          {'instrument_manual_name':'Agilent_8900_Specs.pdf'},
        'Perkin Elmer ELAN 6000':{'instrument_further_details':'No manual or specifications available but confirmation it is uses inductively coupled plasma mass spectrometry here: http://www.speciation.net/Database/Instruments/PerkinElmer-SCIEX/ELAN-6000-;i533'},
        'Thermo iCAP RQ':        {'instrument_manual_name':'Thermo_iCAP_RQ_Specs.pdf'}
    }
},

#------------------------------------
#X-ray fluorescence spectroscopy (XRFS)

#X-Ray Fluorescence Spectroscopy (XRFS) is an analytical technique that uses the interaction of x-rays with a material to determine its elemental composition. XRF is suitable for solids, liquids and powders, and in most circumstances is non-destructive.
#There are two main XRF methodologies - energy dispersive (EDXRF) and wavelength dispersive (WDXRF). Each method has its own advantages and disadvantages.
#The range of detectable elements varies according to instrument configuration and set up, but typically EDXRF covers all elements from sodium (Na) to uranium (U), whilst WDXRF can extend this down to beryllium (Be). Concentrations can range from 100% down to ppm and in some cases sub-ppm levels. 
#Limits of detection depend upon the specific element and the sample matrix, but as a general rule, heavier elements will have better detection limits.
#XRF spectroscopy is the technique of analyzing the fluorescent X-Rays in order to gain information on the elemental composition of a particular material.
#The key components of a typical XRF spectrometer are: 1. Source of X-Rays used to irradiate the sample. 2. Sample. 3. Detection of the emitted fluorescent X-Rays. 4. The resulting XRF spectrum shows intensity of X-Rays (usually in counts per second) as a function of energy (usually in eV).

'X-ray fluorescence spectroscopy':{'abbreviation':'XRFS', 'sampling_type':'injection', 'measured_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{}
},

#------------------------------------
#continuous flow analysis (CFA)

#In continuous flow analysis (CFA), a sample is injected into a flowing carrier solution passing rapidly through small-bore tubing. 
#The sample is mixed with a reagent, which reacts with the sample to develop a color and determine the sample concentration (typically by spectrophotometry).
#The use of carefully controlled flow conditions ensures that the color development reaction is reproducible, so that the color measurement need not wait until the reaction has gone to completion.
#There are several technologies of CFA such as segmented flow analysis (SFA), which uses turbulent flow conditions that allow for complete sample dispersion. 
#Other technologies include flow injected analysis (FIA) and sequential injection analysis, which uses laminar flow existing in the narrow-bore tubing to mix with the reagent to eliminate the need for air bubble partitioning.
#Most continuous flow analyzers depend on color reactions using a flow through photometer, however, also methods have been developed that use ISE, flame photometry, ICAP, fluorometry, and so forth.

'continuous flow analysis':{'abbreviation':'CFA', 'sampling_type':'injection', 'sample_preparation_type':'reagent reaction', 'measured_parameters':['NH3','Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 'qa_accepted_parameters':['NH3','Al','PM10_Al','PM2.5_Al','PM1_Al','As','PM10_As','PM2.5_As','PM1_As','BC','PM10_BC','PM2.5_BC','PM1_BC','eBC','PM10_eBC','PM2.5_eBC','PM1_eBC','C','PM10_C','PM2.5_C','PM1_C','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cd','PM10_Cd','PM2.5_Cd','PM1_Cd','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Co','PM10_Co','PM2.5_Co','PM1_Co','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','Hg','PM10_Hg','PM2.5_Hg','PM1_Hg','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','Mn','PM10_Mn','PM2.5_Mn','PM1_Mn','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','MSA','PM10_MSA','PM2.5_MSA','PM1_MSA','Na+','PM10_Na+','PM2.5_Na+','PM1_Na+','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NH4NO3','PM10_NH4NO3','PM2.5_NH4NO3','PM1_NH4NO3','Ni','PM10_Ni','PM2.5_Ni','PM1_Ni','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','OC','PM10_OC','PM2.5_OC','PM1_OC','Pb','PM10_Pb','PM2.5_Pb','PM1_Pb','Se','PM10_Se','PM2.5_Se','PM1_Se','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS','V','PM10_V','PM2.5_V','PM1_V','Zn','PM10_Zn','PM2.5_Zn','PM1_Zn'], 
    'instruments':{
       'Bran+Luebbe Autoanalyzer 3':{'measured_parameters':['NH3','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','K+','PM10_K+','PM2.5_K+','PM1_K+','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-'],                                                                                                                                                                                                                                                                                                                                           'qa_accepted_parameters':['NH3','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','K+','PM10_K+','PM2.5_K+','PM1_K+','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-'],                                                                                                                                                                                                                                                                                                                                      'instrument_manual_name':'Bran+Luebbe_Autoanalyzer_3_Specs.pdf'},
       'Foss FIAstar 5000':         {'measured_parameters':['NH3','Al','PM10_Al','PM2.5_Al','PM1_Al','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-'],                                                                                                                                                                                                                                                     'qa_accepted_parameters':['NH3','Al','PM10_Al','PM2.5_Al','PM1_Al','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','NH4+','PM10_NH4+','PM2.5_NH4+','PM1_NH4+','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-'],                                                                                                                                                                                                                                                'instrument_manual_name':'Foss_FIAstar_5000_Specs.pdf'},
       'Thermo Gallery':            {'measured_parameters':['SO2','NH3','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'],'qa_accepted_parameters':['NH3','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'instrument_manual_name':'Thermo_Gallery_Specs.pdf'},  
       'Thermo Gallery Plus':       {'measured_parameters':['SO2','NH3','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'],'qa_accepted_parameters':['NH3','Ca++','PM10_Ca++','PM2.5_Ca++','PM1_Ca++','Cl-','PM10_Cl-','PM2.5_Cl-','PM1_Cl-','Cr','PM10_Cr','PM2.5_Cr','PM1_Cr','Cu','PM10_Cu','PM2.5_Cu','PM1_Cu','Fe','PM10_Fe','PM2.5_Fe','PM1_Fe','K+','PM10_K+','PM2.5_K+','PM1_K+','Mg++','PM10_Mg++','PM2.5_Mg++','PM1_Mg++','NO3-','PM10_NO3-','PM2.5_NO3-','PM1_NO3-','SO4--','PM10_SO4--','PM2.5_SO4--','PM1_SO4--','SO4--_NSS','PM10_SO4--_NSS','PM2.5_SO4--_NSS','PM1_SO4--_NSS','SO4--_SS','PM10_SO4--_SS','PM2.5_SO4--_SS','PM1_SO4--_SS'], 'instrument_manual_name':'Thermo_Gallery_Specs.pdf'}  
    }
},

#------------------------------------
#photometry - direct (P-D)

#remote photometry involves measurements of spectral sun/moon irradiance (direct) and sky radiances (sky) 
#A sensor head fitted with 25 cm collimators is attached to a 40 cm robot base which systematically points the sensor head at the Sun, the sky and the moon according to a preprogrammed routine.
#The direct sun measurements are made in eight spectral bands requiring approximately 10 seconds. 
#Eight interference filters at wavelengths of 340, 380, 440, 500, 670, 870, 940 and 1020 nm are located in a filter wheel which is rotated by a direct drive stepping motor. 
#The 940 nm channel is used for column water abundance determination. A preprogrammed sequence of measurements is taken by these instruments starting at an air mass of 7 in the morning and ending at an air mass of 7 in the evening.
#Optical depth is calculated from spectral extinction of direct beam radiation at each wavelength based on the Beer-Bouguer Law. Attenuation due to Rayleigh scatter, and absorption by ozone, and gaseous pollutants is estimated and removed to isolate the aerosol optical depth (AOD). 
#A sequence of three such measurements are taken 30 seconds apart creating a triplet observation per wavelength.
#During the large air mass periods direct sun measurements are made at 0.25 air mass intervals, while at smaller air masses the interval between measurements is typically 15 minutes. 
#The time variation of clouds is usually greater than that of aerosols causing an observable variation in the triplets that can be used to screen clouds in many cases. 
#Additionally the 15-minute interval allows a longer temporal frequency check for cloud contamination.

#The latest models can also perform nighttime measurements of the spectral lunar irradiance.
#See paper here for detailed description of the methodology verified for a specific instrument:
#https://www.atmos-meas-tech.net/9/631/2016/amt-9-631-2016.pdf

#See paper here for more details on the methodology applied at a specific station:
#https://www.atmos-meas-tech.net/10/3007/2017/amt-10-3007-2017.pdf

'photometry - direct':{'abbreviation':'P-D', 'sampling_type':'remote', 'measured_parameters':['AOD_500nm', 'AOD_500nm_COARSE', 'AOD_500nm_FINE', 'FINE_MODE_FRAC_500nm', 'AOD_380nm', 'AOD_440nm', 'AOD_550nm', 'AOD_675nm', 'AOD_870nm', 'AOD_1020nm', 'AE_440-870nm'], 'qa_accepted_parameters':['AOD_500nm', 'AOD_500nm_COARSE', 'AOD_500nm_FINE', 'FINE_MODE_FRAC_500nm', 'AOD_380nm', 'AOD_440nm', 'AOD_550nm', 'AOD_675nm', 'AOD_870nm', 'AOD_1020nm', 'AE_440-870nm'],
    'instruments':{
        'Cimel CE318-1':{                               'instrument_further_details':'No manual or specifications available but confirmation it is uses photometry - direct here: https://aeronet.gsfc.nasa.gov/new_web/system_descriptions_instrument.html'}, 
        'Cimel CE318-2':{                               'instrument_further_details':'No manual or specifications available but confirmation it is uses photometry - direct here: https://aeronet.gsfc.nasa.gov/new_web/system_descriptions_instrument.html'},
        'Cimel CE318-N':{'documented_precision':'0.1%', 'instrument_manual_name':'Cimel_CE318N_Manual.pdf'},
        'Cimel CE318-T':{'documented_precision':'0.1%', 'instrument_manual_name':'Cimel_CE318T_Manual.pdf'}
    }
},

#------------------------------------
#photometry - sky (P-S)

#In addition to the direct solar irradiance measurements remote photometric instruments makes with a field of view of 1.2 degrees, these instruments measure the sky radiance in four spectral bands (440, 670, 870 and 1020 nm) along the solar principal plane (i.e., at constant azimuth angle, with varied scattering angles) up to nine times a day and along the solar almucantar (i.e., at constant elevation angle, with varied azimuth angles) up to six times a day. 
#The approach is to acquire aureole and sky radiances observations through a large range of scattering angles from the sun through a constant aerosol profile to retrieve size distribution, phase function and aerosol optical depth. 
#More than eight almucantar sequences are made daily at an optical air mass of 4, 3, 2 and 1.7 both morning and afternoon. Sky radiance measurements are inverted with the Dubovik and Nakajima inversions to provide aerosol properties of size distribution and phase function over the particle size range of 0.1 to 5 um.
#The robot mounted sensor head is parked pointed near-nadir when idle to prevent contamination of the optical windows from rain and foreign particles.

'photometry - sky':{'abbreviation':'P-S', 'sampling_type':'remote', 'measured_parameters':[], 'qa_accepted_parameters':[],
    'instruments':{
        'Cimel CE318-1':{                               'instrument_further_details':'No manual or specifications available but confirmation it is uses photometry - direct here: https://aeronet.gsfc.nasa.gov/new_web/system_descriptions_instrument.html'}, 
        'Cimel CE318-2':{                               'instrument_further_details':'No manual or specifications available but confirmation it is uses photometry - direct here: https://aeronet.gsfc.nasa.gov/new_web/system_descriptions_instrument.html'},
        'Cimel CE318-N':{'documented_precision':'0.1%', 'instrument_manual_name':'Cimel_CE318N_Manual.pdf'},
        'Cimel CE318-T':{'documented_precision':'0.1%', 'instrument_manual_name':'Cimel_CE318T_Manual.pdf'}
    }
},

#------------------------------------
#unknown methodology

'unknown':{'abbreviation':'?', 'measured_parameters':[], 'qa_accepted_parameters':[], 
    'instruments':{}
}
}

###--------------------------------------------------------------------------------------------------###
#DEFINE STANDARD COUNTRY DEFINITION INFORMATION
###--------------------------------------------------------------------------------------------------###

#dictionary that stores country ISO3166-1 codes with associated country names

country_ISO3166_to_name = {
'AF': 'Afghanistan',
'AX': 'Åland Islands',
'AL': 'Albania',
'DZ': 'Algeria',
'AS': 'American Samoa',
'AD': 'Andorra',
'AO': 'Angola',
'AI': 'Anguilla',
'AQ': 'Antarctica',
'AG': 'Antigua and Barbuda',
'AR': 'Argentina',
'AM': 'Armenia',
'AW': 'Aruba',
'AU': 'Australia',
'AT': 'Austria',
'AZ': 'Azerbaijan',
'BS': 'Bahamas',
'BH': 'Bahrain',
'BD': 'Bangladesh',
'BB': 'Barbados',
'BY': 'Belarus',
'BE': 'Belgium',
'BZ': 'Belize',
'BJ': 'Benin',
'BM': 'Bermuda',
'BT': 'Bhutan',
'BO': 'Bolivia',
'BQ': 'Bonaire, Sint Eustatius and Saba',
'BA': 'Bosnia and Herzegovina',
'BW': 'Botswana',
'BV': 'Bouvet Island',
'BR': 'Brazil',
'IO': 'British Indian Ocean Territory',
'BN': 'Brunei',
'BG': 'Bulgaria',
'BF': 'Burkina Faso',
'BI': 'Burundi',
'KH': 'Cambodia',
'CM': 'Cameroon',
'CA': 'Canada',
'CV': 'Cape Verde',
'KY': 'Cayman Islands',
'CF': 'Central African Republic',
'TD': 'Chad',
'CL': 'Chile',
'CN': 'China',
'CX': 'Christmas Island',
'CC': 'Cocos Islands',
'CO': 'Colombia',
'KM': 'Comoros',
'CG': 'Republic of the Congo',
'CD': 'Democratic Republic of the Congo',
'CK': 'Cook Islands',
'CR': 'Costa Rica',
'CI': "Côte d'Ivoire",
'HR': 'Croatia',
'CU': 'Cuba',
'CW': 'Curaçao',
'CY': 'Cyprus',
'CZ': 'Czech Republic',
'DK': 'Denmark',
'DJ': 'Djibouti',
'DM': 'Dominica',
'DO': 'Dominican Republic',
'EC': 'Ecuador',
'EG': 'Egypt',
'SV': 'El Salvador',
'GQ': 'Equatorial Guinea',
'ER': 'Eritrea',
'EE': 'Estonia',
'ET': 'Ethiopia',
'FK': 'Falkland Islands',
'FO': 'Faroe Islands',
'FJ': 'Fiji',
'FI': 'Finland',
'FR': 'France',
'GF': 'French Guiana',
'PF': 'French Polynesia',
'TF': 'French Southern Territories',
'GA': 'Gabon',
'GM': 'Gambia',
'GE': 'Georgia',
'DE': 'Germany',
'GH': 'Ghana',
'GI': 'Gibraltar',
'GR': 'Greece',
'GL': 'Greenland',
'GD': 'Grenada',
'GP': 'Guadeloupe',
'GU': 'Guam',
'GT': 'Guatemala',
'GG': 'Guernsey',
'GN': 'Guinea',
'GW': 'Guinea-Bissau',
'GY': 'Guyana',
'HT': 'Haiti',
'HM': 'Heard Island and McDonald Islands',
'VA': 'Vatican',
'HN': 'Honduras',
'HK': 'Hong Kong',
'HU': 'Hungary',
'IS': 'Iceland',
'IN': 'India',
'ID': 'Indonesia',
'IR': 'Iran',
'IQ': 'Iraq',
'IE': 'Ireland',
'IM': 'Isle of Man',
'IL': 'Israel',
'IT': 'Italy',
'JM': 'Jamaica',
'JP': 'Japan',
'JE': 'Jersey',
'JO': 'Jordan',
'KZ': 'Kazakhstan',
'KE': 'Kenya',
'KI': 'Kiribati',
'KP': 'North Korea',
'KR': 'South Korea',
'KW': 'Kuwait',
'KG': 'Kyrgyzstan',
'LA': "Lao People's Democratic Republic",
'LV': 'Latvia',
'LB': 'Lebanon',
'LS': 'Lesotho',
'LR': 'Liberia',
'LY': 'Libya',
'LI': 'Liechtenstein',
'LT': 'Lithuania',
'LU': 'Luxembourg',
'MO': 'Macao',
'MK': 'North Macedonia',
'MG': 'Madagascar',
'MW': 'Malawi',
'MY': 'Malaysia',
'MV': 'Maldives',
'ML': 'Mali',
'MT': 'Malta',
'MH': 'Marshall Islands',
'MQ': 'Martinique',
'MR': 'Mauritania',
'MU': 'Mauritius',
'YT': 'Mayotte',
'MX': 'Mexico',
'FM': 'Micronesia',
'MD': 'Moldova',
'MC': 'Monaco',
'MN': 'Mongolia',
'ME': 'Montenegro',
'MS': 'Montserrat',
'MA': 'Morocco',
'MZ': 'Mozambique',
'MM': 'Myanmar',
'NA': 'Namibia',
'NR': 'Nauru',
'NP': 'Nepal',
'NL': 'Netherlands',
'NC': 'New Caledonia',
'NZ': 'New Zealand',
'NI': 'Nicaragua',
'NE': 'Niger',
'NG': 'Nigeria',
'NU': 'Niue',
'NF': 'Norfolk Island',
'MP': 'Northern Mariana Islands',
'NO': 'Norway',
'OM': 'Oman',
'PK': 'Pakistan',
'PW': 'Palau',
'PS': 'Palestinian Territory',
'PA': 'Panama',
'PG': 'Papua New Guinea',
'PY': 'Paraguay',
'PE': 'Peru',
'PH': 'Philippines',
'PN': 'Pitcairn',
'PL': 'Poland',
'PT': 'Portugal',
'PR': 'Puerto Rico',
'QA': 'Qatar',
'RE': 'Réunion',
'RO': 'Romania',
'RU': 'Russia',
'RW': 'Rwanda',
'BL': 'Saint Barthélemy',
'SH': 'Saint Helena, Ascension and Tristan da Cunha',
'KN': 'Saint Kitts and Nevis',
'LC': 'Saint Lucia',
'MF': 'Saint Martin',
'PM': 'Saint Pierre and Miquelon',
'VC': 'Saint Vincent and the Grenadines',
'WS': 'Samoa',
'SM': 'San Marino',
'ST': 'Sao Tome and Principe',
'SA': 'Saudi Arabia',
'SN': 'Senegal',
'RS': 'Serbia',
'SC': 'Seychelles',
'SL': 'Sierra Leone',
'SG': 'Singapore',
'SX': 'Sint Maarten',
'SK': 'Slovakia',
'SI': 'Slovenia',
'SB': 'Solomon Islands',
'SO': 'Somalia',
'ZA': 'South Africa',
'GS': 'South Georgia and the South Sandwich Islands',
'SS': 'South Sudan',
'ES': 'Spain',
'LK': 'Sri Lanka',
'SD': 'Sudan',
'SR': 'Suriname',
'SJ': 'Svalbard and Jan Mayen',
'SZ': 'Swaziland',
'SE': 'Sweden',
'CH': 'Switzerland',
'SY': 'Syria',
'TW': 'Taiwan',
'TJ': 'Tajikistan',
'TZ': 'Tanzania',
'TH': 'Thailand',
'TL': 'Timor-Leste',
'TG': 'Togo',
'TK': 'Tokelau',
'TO': 'Tonga',
'TT': 'Trinidad and Tobago',
'TN': 'Tunisia',
'TR': 'Turkey',
'TM': 'Turkmenistan',
'TC': 'Turks and Caicos Islands',
'TV': 'Tuvalu',
'UG': 'Uganda',
'UA': 'Ukraine',
'AE': 'United Arab Emirates',
'GB': 'United Kingdom',
'US': 'United States',
'UM': 'United States Minor Outlying Islands',
'UY': 'Uruguay',
'UZ': 'Uzbekistan',
'VU': 'Vanuatu',
'VE': 'Venezuela',
'VN': 'Vietnam',
'VG': 'British Virgin Islands',
'VI': 'U.S. Virgin Islands',
'WF': 'Wallis and Futuna',
'EH': 'Western Sahara',
'XK': 'Kosovo',
'YE': 'Yemen',
'ZM': 'Zambia',
'ZW': 'Zimbabwe',
'CS': 'Serbia and Montenegro'
}

#dictionary that stores country FIPS 10-4 codes with associated country names

country_FIPS104_to_name = {
'AN': 'Andorra',
'AE': 'United Arab Emirates',
'AF': 'Afghanistan',
'AC': 'Antigua and Barbuda',
'AV': 'Anguilla',
'AL': 'Albania',
'AM': 'Armenia',
'AO': 'Angola',
'AY': 'Antarctica',
'AR': 'Argentina',
'AQ': 'American Samoa',
'AU': 'Austria',
'AS': 'Australia',
'AA': 'Aruba',
'AJ': 'Azerbaijan',
'BK': 'Bosnia and Herzegovina',
'BB': 'Barbados',
'BG': 'Bangladesh',
'BE': 'Belgium',
'UV': 'Burkina Faso',
'BU': 'Bulgaria',
'BA': 'Bahrain',
'BY': 'Burundi',
'BN': 'Benin',
'TB': 'Saint Barthélemy',
'BD': 'Bermuda',
'BX': 'Brunei',
'BL': 'Bolivia',
'BR': 'Brazil',
'BF': 'Bahamas',
'BT': 'Bhutan',
'BV': 'Bouvet Island',
'BC': 'Botswana',
'BO': 'Belarus',
'BH': 'Belize',
'CA': 'Canada',
'CH': 'China',
'CK': 'Cocos Islands',
'CG': 'Democratic',
'CT': 'Central African Republic',
'CF': 'Republic of the Congo',
'CG': 'Democratic Republic of the Congo',
'SZ': 'Switzerland',
'IV': "Côte d'Ivoire",
'CW': 'Cook Islands',
'CI': 'Chile',
'CM': 'Cameroon',
'CO': 'Colombia',
'CS': 'Costa Rica',
'CU': 'Cuba',
'CV': 'Cape Verde',
'UC': 'Curaçao',
'KT': 'Christmas Island',
'CY': 'Cyprus',
'EZ': 'Czech Republic',
'GM': 'Germany',
'DJ': 'Djibouti',
'DA': 'Denmark',
'DO': 'Dominica',
'DR': 'Dominican Republic',
'AG': 'Algeria',
'EC': 'Ecuador',
'EN': 'Estonia',
'EG': 'Egypt',
'WI': 'Western Sahara',
'ER': 'Eritrea',
'SP': 'Spain',
'ET': 'Ethiopia',
'FI': 'Finland',
'FG': 'French Guiana',
'FJ': 'Fiji',
'FK': 'Falkland Islands',
'FM': 'Micronesia',
'FO': 'Faroe Islands',
'FR': 'France',
'FS': 'French Southern Territories',
'GB': 'Gabon',
'UK': 'United Kingdom',
'GJ': 'Grenada',
'GG': 'Georgia',
'FG': 'French',
'GK': 'Guernsey',
'GH': 'Ghana',
'GI': 'Gibraltar',
'GL': 'Greenland',
'GA': 'Gambia',
'GV': 'Guinea',
'GP': 'Guadeloupe',
'EK': 'Equatorial Guinea',
'GR': 'Greece',
'SX': 'South Georgia and the South Sandwich Islands',
'GT': 'Guatemala',
'GQ': 'Guam',
'PU': 'Guinea-Bissau',
'GY': 'Guyana',
'HK': 'Hong Kong',
'HM': 'Heard Island and McDonald Islands',
'HO': 'Honduras',
'HR': 'Croatia',
'HA': 'Haiti',
'HU': 'Hungary',
'ID': 'Indonesia',
'EI': 'Ireland',
'IS': 'Israel',
'IM': 'Isle of Man',
'IN': 'India',
'IO': 'British Indian Ocean Territory',
'IZ': 'Iraq',
'IR': 'Iran',
'IC': 'Iceland',
'IT': 'Italy',
'JE': 'Jersey',
'JM': 'Jamaica',
'JO': 'Jordan',
'JA': 'Japan',
'KE': 'Kenya',
'KG': 'Kyrgyzstan',
'CB': 'Cambodia',
'KR': 'Kiribati',
'CN': 'Comoros',
'SC': 'Saint Kitts and Nevis',
'KN': 'North Korea',
'KS': 'South Korea',
'KV': 'Kosovo',
'KU': 'Kuwait',
'CJ': 'Cayman Islands',
'KZ': 'Kazakhstan',
'LA': 'Laos',
'LE': 'Lebanon',
'ST': 'Saint Lucia',
'LS': 'Liechtenstein',
'CE': 'Sri Lanka',
'LI': 'Liberia',
'LT': 'Lesotho',
'LH': 'Lithuania',
'LU': 'Luxembourg',
'LG': 'Latvia',
'LY': 'Libya',
'MO': 'Morocco',
'MN': 'Monaco',
'MD': 'Moldova',
'MJ': 'Montenegro',
'RN': 'Saint Martin',
'MA': 'Madagascar',
'RM': 'Marshall Islands',
'MK': 'North Macedonia',
'ML': 'Mali',
'BM': 'Myanmar',
'MG': 'Mongolia',
'MC': 'Macao',
'CQ': 'Northern Mariana Islands',
'MB': 'Martinique',
'MR': 'Mauritania',
'MH': 'Montserrat',
'MT': 'Malta',
'MP': 'Mauritius',
'MV': 'Maldives',
'MI': 'Malawi',
'MX': 'Mexico',
'MY': 'Malaysia',
'MZ': 'Mozambique',
'WA': 'Namibia',
'NC': 'New Caledonia',
'NG': 'Niger',
'NF': 'Norfolk Island',
'NI': 'Nigeria',
'NU': 'Nicaragua',
'NL': 'Netherlands',
'NO': 'Norway',
'NP': 'Nepal',
'NR': 'Nauru',
'NE': 'Niue',
'NZ': 'New Zealand',
'MU': 'Oman',
'PM': 'Panama',
'PE': 'Peru',
'FP': 'French Polynesia',
'PP': 'Papua New Guinea',
'RP': 'Philippines',
'PK': 'Pakistan',
'PL': 'Poland',
'SB': 'Saint Pierre and Miquelon',
'PC': 'Pitcairn',
'RQ': 'Puerto Rico',
'WE': 'Palestinian Territory',
'PO': 'Portugal',
'PS': 'Palau',
'PA': 'Paraguay',
'QA': 'Qatar',
'RE': 'Réunion',
'RO': 'Romania',
'RI': 'Serbia',
'RS': 'Russia',
'RW': 'Rwanda',
'SA': 'Saudi Arabia',
'BP': 'Solomon Islands',
'SE': 'Seychelles',
'SU': 'Sudan',
'OD': 'South Sudan',
'SW': 'Sweden',
'SN': 'Singapore',
'SH': 'Saint Helena, Ascension and Tristan da Cunha',
'SI': 'Slovenia',
'SV': 'Svalbard and Jan Mayen',
'LO': 'Slovakia',
'SL': 'Sierra Leone',
'SM': 'San Marino',
'SG': 'Senegal',
'SO': 'Somalia',
'NS': 'Suriname',
'TP': 'Sao Tome and Principe',
'ES': 'El Salvador',
'NN': 'Sint Maarten',
'SY': 'Syria',
'WZ': 'Swaziland',
'TK': 'Turks and Caicos Islands',
'CD': 'Chad',
'FS': 'French',
'TO': 'Togo',
'TH': 'Thailand',
'TI': 'Tajikistan',
'TL': 'Tokelau',
'TT': 'Timor-Leste',
'TX': 'Turkmenistan',
'TS': 'Tunisia',
'TN': 'Tonga',
'TU': 'Turkey',
'TD': 'Trinidad and Tobago',
'TV': 'Tuvalu',
'TW': 'Taiwan',
'TZ': 'Tanzania',
'UP': 'Ukraine',
'UG': 'Uganda',
'UM': 'United States Minor Outlying Islands',
'US': 'United States',
'UY': 'Uruguay',
'UZ': 'Uzbekistan',
'VT': 'Vatican',
'VC': 'Saint Vincent and the Grenadines',
'VE': 'Venezuela',
'VI': 'British Virgin Islands',
'VQ': 'U.S. Virgin Islands',
'VM': 'Vietnam',
'NH': 'Vanuatu',
'WF': 'Wallis and Futuna',
'WS': 'Samoa',
'YM': 'Yemen',
'MF': 'Mayotte',
'SF': 'South Africa',
'ZA': 'Zambia',
'ZI': 'Zimbabwe',
'YI': 'Serbia and Montenegro'}
###--------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------###