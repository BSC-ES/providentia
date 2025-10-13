# CAMS

Providentia's download mode supports downloading experiment data from **Atmosphere Data Store** datasets provided by CAMS:

[https://ads.atmosphere.copernicus.eu/datasets](https://ads.atmosphere.copernicus.eu/datasets)

## Available datasets

Currently, four CAMS datasets are available for download:

1. [**CAMS European Air Quality Forecasts**](#1-cams-european-air-quality-forecasts)
2. [**CAMS Global Atmospheric Composition Forecasts**](#2-cams-global-atmospheric-composition-forecasts)
3. [**CAMS European Air Quality Reanalyses**](#3-cams-european-air-quality-reanalyses)
4. [**CAMS Global Reanalysis (EAC4)**](#4-cams-global-reanalysis-eac4)

## 1. CAMS European Air Quality Forecasts

[Dataset Link](https://ads.atmosphere.copernicus.eu/datasets/cams-europe-air-quality-forecasts?tab=overview)

Using this dataset you can either download:

* **Regional Analysis data**, or
* **Regional Forecast data**

This dataset contains multi-level data.

This dataset provides hourly data only. You must set `resolution = hourly` in your configuration file; otherwise, the download will not work.

### Mandatory fields for regional analysis data

```ini
experiment = cams_analysis_<model>
domain = regional
resolution = hourly
```

### Mandatory fields for regional forecast data

```ini
experiment = cams_forecast_<model>
domain = regional
resolution = hourly
```

> Replace `<model>` with one of the available model names listed below.

### Available models

Use the following model names (same as in the API):

* `ensemble`
* `chimere`
* `dehm`
* `emep`
* `euradim`
* `gemaq`
* `lotos`
* `match`
* `minni`
* `mocage`
* `monarch`
* `silam`

### Available variables (Species)

These are the available species:

* `ammonia`
* `carbon_monoxide`
* `dust`
* `formaldehyde`
* `glyoxal`
* `nitrogen_dioxide`
* `nitrogen_monoxide`
* `non_methane_vocs`
* `ozone`
* `particulate_matter_10um`
* `particulate_matter_2.5um`
* `pm10_sea_salt_dry`
* `pm2.5_ammonium`
* `pm2.5_nitrate`
* `pm2.5_sulphate`
* `sulphur_dioxide`
* `total_elementary_carbon`

> Providentia can only read species in GHOST format. The mapping is provided in the [GHOST–CAMS species mapping](#ghostcams-species-mapping) section.

### Assumptions

Providentia assumes the following fixed values when downloading data:

* `time = 00:00-23:00`
* `leadtime_hour = 0-96`
* `level = 0`

> All available data in this dataset is multi-level, but Providentia retrieves only the ground level.

### Example configuration files and corresponding API requests

#### Example 1: CAMS Regional Analysis (DEHM Model, sconcno2)

```ini
[cams_forecast_dehm-regional] 
start_date = 20250701 
end_date = 20250702
species = sconcno2
experiments = cams_forecast_dehm
domain = regional
resolution = hourly
interpolated = False
```

```python
dataset = 'cams-europe-air-quality-forecasts'
request = {
'variable' : 'nitrogen_dioxide',
'leadtime_hour': [
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
        '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23',
        '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34',
        '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45',
        '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56',
        '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67',
        '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78',
        '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89',
        '90', '91', '92', '93', '94', '95', '96'
    ],
'type' : 'forecast',
'date' : '2025-07-01/2025-07-01',
'model' : 'dehm',
'level' : '0',
'time' : ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
        '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'],
'data_format' : 'netcdf_zip',
}
```

#### Example 2: CAMS Regional Forecast (Ensemble Model, sconcno2)

```ini
[cams_forecast_ensemble-regional]
start_date = 20250701
end_date = 20250702
species = sconcno2
experiments = cams_forecast_ensemble
domain = regional
resolution = hourly
interpolated = False
```

```python
dataset = 'cams-europe-air-quality-forecasts'
request = {
    'variable': 'nitrogen_dioxide',
    'leadtime_hour': [
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
        '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23',
        '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34',
        '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45',
        '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56',
        '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67',
        '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78',
        '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89',
        '90', '91', '92', '93', '94', '95', '96'
    ],
    'type': 'forecast',
    'date': '2025-07-01/2025-07-01',
    'model': 'ensemble',
    'level': '0',
    'time': ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
        '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'],
    'data_format': 'netcdf_zip',
}
```

## 2. CAMS Global Atmospheric Composition Forecasts

[Dataset Link](https://ads.atmosphere.copernicus.eu/datasets/cams-global-atmospheric-composition-forecasts?tab=overview)

Using this dataset you can download **Global Forecast data**.

This dataset contains both single-level and multi-level data.

This dataset provides hourly and 3hourly data. In your configuration file, you must set `resolution = hourly` or `resolution = 3hourly`, depending on whether the data is single-level or multi-level. Otherwise, the download will fail.

### Mandatory fields for global forecast single-level data

```ini
experiment = cams_forecast
domain = global
resolution = hourly
```

### Mandatory fields for global forecast multi-level data

```ini
experiment = cams_forecast
domain = global
resolution = 3hourly
```

### Available variables (Species)

Whether the data is single or multi depends on the selected species. These are the available species:

#### Single-level:

* `2m_dewpoint_temperature`
* `2m_temperature`
* `asymmetry_factor_1020nm`
* `asymmetry_factor_440nm`
* `cloud_base_height`
* `dust_aerosol_optical_depth_550nm`
* `mean_sea_level_pressure`
* `particulate_matter_10um`
* `particulate_matter_1um`
* `particulate_matter_2.5um`
* `relative_humidity`
* `sea_surface_temperature`
* `single_scattering_albedo_1020nm`
* `single_scattering_albedo_440nm`
* `single_scattering_albedo_550nm`
* `snow_depth`
* `surface_pressure`
* `total_absorption_aerosol_optical_depth_1020nm`
* `total_absorption_aerosol_optical_depth_440nm`
* `total_absorption_aerosol_optical_depth_550nm`
* `total_aerosol_optical_depth_1020nm`
* `total_aerosol_optical_depth_380nm`
* `total_aerosol_optical_depth_440nm`
* `total_aerosol_optical_depth_500nm`
* `total_aerosol_optical_depth_550nm`
* `total_fine_mode_aerosol_optical_depth_500nm`
* `total_precipitation`
* `visibility`

#### Multi-level:

* `ammonia`
* `ammonium`
* `attenuated_backscatter_due_to_aerosol_532nm_from_ground`
* `carbon_monoxide`
* `ethane`
* `ethanol`
* `ethene`
* `formaldehyde`
* `fraction_of_cloud_cover`
* `glyoxal`
* `hydrogen_chloride`
* `hydrogen_fluoride`
* `isoprene`
* `lead`
* `methane`
* `methane_sulfonic_acid`
* `nitrate`
* `nitric_acid`
* `nitrogen_dioxide`
* `nitrogen_monoxide`
* `ozone`
* `peroxyacetyl_nitrate`
* `propane`
* `propene`
* `sulphate_aerosol_mixing_ratio`
* `sulphur_dioxide`

> Providentia can only read species in GHOST format. The mapping is provided in the [GHOST–CAMS species mapping](#ghostcams-species-mapping) section.

### Assumptions

Providentia assumes the following fixed values when downloading data:

#### For single-level data:

* `time = 00:00`
* `leadtime_hour = 0-120` (for every hour)
* `type = forecast`

#### For multi-level data:

* `time = 00:00`
* `leadtime_hour = 0-120` (for every three hours)
* `type = forecast`
* `model_level = 60` (up to and including 2019-07-09) `137` (for data after 2019-07-09, when the levels increased)

> For the multi-level data in this dataset, Providentia retrieves only the ground level.

### Example configuration files and corresponding API requests

#### Example 1: CAMS Global Forecast (sst, Single-level, hourly)

```ini
[cams_forecast-global-hourly]
start_date = 20171201 
end_date = 20171202
species = sst
experiments = cams_forecast
domain = global
resolution = hourly
interpolated = False
```

```python
dataset = 'cams-global-atmospheric-composition-forecasts'
request = {
'variable' : 'sea_surface_temperature',
'leadtime_hour' : ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120'],
'type' : 'forecast',
'date' : '2017-12-01/2017-12-01',
'time' : '00:00',
'data_format' : 'netcdf_zip',
}
```

#### Example 2: CAMS Global Forecast (sconcno2, Multi-level, 3hourly)

```ini
[cams_forecast-global-3hourly]
start_date = 20200601  
end_date = 20200602
species = sconcno2
experiments = cams_forecast
domain = global
resolution = 3hourly
interpolated = False
```

```python
dataset = 'cams-global-atmospheric-composition-forecasts'
request = {
'variable' : 'nitrogen_dioxide',
'leadtime_hour' : ['0', '3', '6', '9', '12', '15', '18', '21', '24', '27', '30', '33', '36', '39', '42', '45', '48', '51', '54', '57', '60', '63', '66', '69', '72', '75', '78', '81', '84', '87', '90', '93', '96', '99', '102', '105', '108', '111', '114', '117', '120'],
'type' : 'forecast',
'date' : '2020-06-01/2020-06-01',
'model_level' : '137',
'time' : '00:00',
'data_format' : 'netcdf_zip',
}
```

## 3. CAMS European Air Quality Reanalyses

[Dataset Link](https://ads.atmosphere.copernicus.eu/datasets/cams-europe-air-quality-reanalyses?tab=overview)

Using this dataset you can download **Regional Reanalysis data**.

This dataset contains multi-level data.

This dataset provides hourly data only. You must set `resolution = hourly` in your configuration file; otherwise, the download will not work.

### Mandatory fields for regional reanalysis data

```ini
experiment = cams_reanalysis_<model>_<stream>
domain = regional
resolution = hourly
```

> Replace `<model>` and `<stream>` with one of the available model names and streams listed below.

### Available models

Use the following model names (same as in the API):

* `ensemble`
* `chimere`
* `emep`
* `lotos`
* `match`
* `minni`
* `mocage`
* `monarch`
* `silam`
* `euradim`
* `dehm`
* `gemaq`

### Available streams

There is two available streams for the regional reanalysis dataset:

* `validated` available for the period **2014–2021 (inclusive)**
* `interim` available for the period **2021–2024 (inclusive)**

### Available variables (Species)

These are the available species:

* `ammonia`
* `carbon_monoxide`
* `dust`
* `formaldehyde`
* `glyoxal`
* `nitrogen_dioxide`
* `nitrogen_monoxide`
* `non_methane_vocs`
* `ozone`
* `particulate_matter_10um`
* `particulate_matter_2.5um`
* `pm10_sea_salt_dry`
* `sulphur_dioxide`
* `total_elementary_carbon`

> Providentia can only read species in GHOST format. The mapping is provided in the [GHOST–CAMS species mapping](#ghostcams-species-mapping) section.

### Assumptions

Providentia assumes the following fixed values when downloading data:

* `level = 0`

> All available data in this dataset is multi-level, but Providentia retrieves only the ground level.

### Example configuration files and corresponding API requests

#### Example 1: CAMS Regional Reanalysis (Ensemble Model, sconcno2, Interim Reanalysis)

```ini
[cams_reanalysis_ensemble_interim-regional]
start_date = 20240101
end_date = 20240201
species = sconcno2
experiments = cams_reanalysis_ensemble_interim
domain = regional
resolution = hourly
interpolated = False
```

```python
dataset = 'cams-europe-air-quality-reanalyses'
request = {
'variable' : 'nitrogen_dioxide',
'year' : '2024',
'month' : '01',
'type' : 'interim_reanalysis',
'model' : 'ensemble',
'level' : '0',
}
```

#### Example 2: CAMS Regional Reanalysis (Ensemble Model, sconcno2, Validated Reanalysis)

```ini
[cams_reanalysis_ensemble_validated-regional] 
start_date = 20200101 
end_date = 20200301
species = sconcno2
experiments = cams_reanalysis_ensemble_validated
domain = regional
resolution = hourly
interpolated = False
```

```python
dataset = 'cams-europe-air-quality-reanalyses'
request = {
'variable' : 'nitrogen_dioxide',
'year' : '2020',
'month' : '01',
'type' : 'validated_reanalysis',
'model' : 'ensemble',
'level' : '0',
}
```

## 4. CAMS Global Reanalysis (EAC4)

[Dataset Link](https://ads.atmosphere.copernicus.eu/datasets/cams-global-reanalysis-eac4?tab=overview)

Using this dataset you can download **Global Reanalysis data**.

This dataset contains both single-level and multi-level data.

This dataset provides 3hourly data only. You must set `resolution = 3hourly` in your configuration file; otherwise, the download will not work.

### Mandatory fields for regional reanalysis data

```ini
experiment = cams_reanalysis
domain = global
resolution = 3hourly
```

### Available variables (Species)

Whether the data is single or multi depends on the selected species. These are the available species:

#### Single-level:

* `2m_dewpoint_temperature`
* `2m_temperature`
* `dust_aerosol_optical_depth_550nm`
* `mean_sea_level_pressure`
* `particulate_matter_10um`
* `particulate_matter_1um`
* `particulate_matter_2.5um`
* `relative_humidity`
* `sea_surface_temperature`
* `snow_depth`
* `surface_pressure`
* `total_aerosol_optical_depth_550nm`

#### Multi-level:

* `ammonia`
* `ammonium`
* `carbon_monoxide`
* `ethane`
* `ethanol`
* `ethene`
* `formaldehyde`
* `fraction_of_cloud_cover`
* `isoprene`
* `lead`
* `methane_sulfonic_acid`
* `nitrate`
* `nitric_acid`
* `nitrogen_dioxide`
* `nitrogen_monoxide`
* `ozone`
* `peroxyacetyl_nitrate`
* `propane`
* `propene`
* `sulphate_aerosol_mixing_ratio`
* `sulphur_dioxide`

> Providentia can only read species in GHOST format. The mapping is provided in the [GHOST–CAMS species mapping](#ghostcams-species-mapping) section.

### Assumptions

Providentia assumes the following fixed values when downloading data:

#### For single-level data:

* `time = 00:00-21:00` (for every three hours)

#### For multi-level data:

* `time = 00:00-21:00` (for every three hours)
* `model_level = 60` 

> For the multi-level data in this dataset, Providentia retrieves only the ground level.

### Example configuration files and corresponding API requests

#### Example 1: CAMS Global Reanalysis (sst, Single-level)

```ini
[cams_reanalysis-global-single]
start_date = 20181101
end_date = 20181130
species = sst
experiments = cams_reanalysis
domain = global
resolution = 3hourly
interpolated = False
```

```python
dataset = 'cams-global-reanalysis-eac4'
request = {
'variable' : 'sea_surface_temperature',
'date' : '2018-11-01/2018-11-29',
'time' : ['00:00', '03:00', '06:00', '09:00', '12:00', '15:00', '18:00', '21:00'],
'data_format' : 'netcdf_zip',
}
```

#### Example 2: CAMS Global Reanalysis (sconcno2, Multi-level)

```ini
[cams_reanalysis-global-multi]
start_date = 20181101
end_date = 20181130
species = sconcno2
experiments = cams_reanalysis
domain = global
resolution = 3hourly
interpolated = False
```

```python
dataset = 'cams-global-reanalysis-eac4'
request = {
'variable' : 'nitrogen_dioxide',
'date' : '2018-11-01/2018-11-29',
'model_level' : '60',
'time' : ['00:00', '03:00', '06:00', '09:00', '12:00', '15:00', '18:00', '21:00'],
'data_format' : 'netcdf_zip',
}
```

## GHOST–CAMS species mapping

Here is the mapping between GHOST variable names and CAMS species:

```
absod1020aero: total_absorption_aerosol_optical_depth_1020nm
absod440aero: total_absorption_aerosol_optical_depth_440nm
absod550aero: total_absorption_aerosol_optical_depth_550nm
acprec: total_precipitation
asy1020aero: asymmetry_factor_1020nm
asy440aero: asymmetry_factor_440nm
cfracmax: fraction_of_cloud_cover
cldbot: cloud_base_height
lbsco532: attenuated_backscatter_due_to_aerosol_532nm_from_ground
od1020aero: total_aerosol_optical_depth_1020nm
od380aero: total_aerosol_optical_depth_380nm
od440aero: total_aerosol_optical_depth_440nm
od500aero: total_aerosol_optical_depth_500nm
od500aerofine: total_fine_mode_aerosol_optical_depth_500nm
od550aero: total_aerosol_optical_depth_550nm
od550du: dust_aerosol_optical_depth_550nm
pm1: particulate_matter_1um
pm10: particulate_matter_10um
pm10du: dust
pm10so4ss: pm10_sea_salt_dry
pm2p5: particulate_matter_2.5um
pm2p5ec: total_elementary_carbon
pm2p5nh4: pm2.5_ammonium
pm2p5no3: pm2.5_nitrate
pm2p5so4: pm2.5_sulphate
pshltr: surface_pressure
rh2: relative_humidity
sca1020aero: single_scattering_albedo_1020nm
sca440aero: single_scattering_albedo_440nm
sca550aero: single_scattering_albedo_550nm
sconcc2h4: ethene
sconcc2h6: ethane
sconcc3h6: propene
sconcc3h8: propane
sconcch4: methane
sconcco: carbon_monoxide
sconcetoh: ethanol
sconcglyox: glyoxal
sconchcho: formaldehyde
sconchcl: hydrogen_chloride
sconchf: hydrogen_fluoride
sconchno3: nitric_acid
sconcisop: isoprene
sconcmsa: methane_sulfonic_acid
sconcnh3: ammonia
sconcnh4: ammonium
sconcnmvoc: non_methane_vocs
sconcno: nitrogen_monoxide
sconcno2: nitrogen_dioxide
sconcno3: nitrate
sconco3: ozone
sconcpan: peroxyacetyl_nitrate
sconcpb: lead
sconcso2: sulphur_dioxide
sconcso4: sulphate_aerosol_mixing_ratio
si: snow_depth
slp: mean_sea_level_pressure
sst: sea_surface_temperature
t2: 2m_temperature
td2: 2m_dewpoint_temperature
vdist: visibility
```
