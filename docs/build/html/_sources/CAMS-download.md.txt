# CAMS

Providentia's download mode supports downloading model output data provided by **CAMS** from the **Atmosphere Data Store**:

[https://ads.atmosphere.copernicus.eu/datasets](https://ads.atmosphere.copernicus.eu/datasets)

An ECMWF account is required to download the files. For a tutorial on how to create and set up your account, see the [Account setup in Atmosphere Data Store](#ADS) page.

Data requests may take some time to complete. After starting a request from Providentia, you can monitor its progress by following the steps described in the [Check CAMS download](Check-CAMS-download) page.

## Available CAMS datasets

Currently, four CAMS datasets are available for download:

1. [**CAMS European Air Quality Forecasts**](#1-cams-european-air-quality-forecasts)
2. [**CAMS Global Atmospheric Composition Forecasts**](#2-cams-global-atmospheric-composition-forecasts)
3. [**CAMS European Air Quality Reanalyses**](#3-cams-european-air-quality-reanalyses)
4. [**CAMS Global Reanalysis (EAC4)**](#4-cams-global-reanalysis-eac4)

(how-to-enable-ads-download)=
## How to enable Atmosphere Data Store download 

To download CAMS data from the Atmosphere Data Store, include one of the following model names along with the domain. For details about available models, streams, and specific resolutions, see the corresponding sections below:

* `cams_analysis_<model>-regional` — [CAMS European Air Quality Forecasts](#1-cams-european-air-quality-forecasts)
* `cams_forecast_<model>-regional` — [CAMS European Air Quality Forecasts](#1-cams-european-air-quality-forecasts)
* `cams_forecast-global` — [CAMS Global Atmospheric Composition Forecasts](#2-cams-global-atmospheric-composition-forecasts)
* `cams_reanalysis_<model>_<stream>-regional` — [CAMS European Air Quality Reanalyses](#3-cams-european-air-quality-reanalyses)
* `cams_reanalysis-global` — [CAMS Global Reanalysis (EAC4)](#4-cams-global-reanalysis-eac4)

Once you have selected the dataset and specified the model and domain, make sure to set:

```ini
dl_interpolated = False
```

or

Answer `n` to the prompt:

_"Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded)"_

## 1. CAMS European Air Quality Forecasts

[Dataset Link](https://ads.atmosphere.copernicus.eu/datasets/cams-europe-air-quality-forecasts?tab=overview)

Using this dataset you can either download:

* **Regional Analysis data**, or
* **Regional Forecast data**

This dataset contains multi-level data.

Only hourly data is available in this dataset. You must set `resolution = hourly` in your configuration file; otherwise, the download will not work.

### Mandatory fields for regional analysis data

```ini
model = cams_analysis_<model>
domain = regional
resolution = hourly
```

### Mandatory fields for regional forecast data

```ini
model = cams_forecast_<model>
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

* `pm10`
* `pm10du`
* `pm10so4ss`
* `pm2p5`
* `pm2p5ec`
* `pm2p5nh4`
* `pm2p5no3`
* `pm2p5so4`
* `sconcco`
* `sconcglyox`
* `sconchcho`
* `sconcnh3`
* `sconcnmvoc`
* `sconcno`
* `sconcno2`
* `sconco3`
* `sconcso2`

> Providentia can only read species in GHOST format. If you want to know the mapping from CAMS variables to GHOST species, please refer to the [CAMS-GHOST species mapping](#cams-ghost-species-mapping) section.

### Assumptions

Providentia assumes the following fixed values when downloading data:

* `time = 00:00-23:00`
* `leadtime_hour = 0-96`
* `level = 0`

> All available data in this dataset is multi-level, but Providentia retrieves only the ground level.

### Example configuration files and corresponding API requests

#### Example 1: CAMS Regional Analysis (DEHM Model, sconcno2)

```ini
[cams_analysis_dehm-regional] 
start_date = 20250701 
end_date = 20250702
species = sconcno2
model = cams_analysis_dehm
domain = regional
resolution = hourly
dl_interpolated = False
```

```python
dataset = 'cams-europe-air-quality-forecasts'
request = {
'variable' : 'nitrogen_dioxide',
'leadtime_hour' : '0',
'type' : 'analysis',
'date' : '2025-07-01/2025-07-01',
'model' : 'dehm',
'level' : '0',
'time' : ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'],
'data_format' : 'netcdf_zip',
}
```

#### Example 2: CAMS Regional Forecast (Ensemble Model, sconcno2)

```ini
[cams_forecast_ensemble-regional]
start_date = 20250701
end_date = 20250702
species = sconcno2
model = cams_forecast_ensemble
domain = regional
resolution = hourly
dl_interpolated = False
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
    'time': '00:00',
    'data_format': 'netcdf_zip',
}
```

## 2. CAMS Global Atmospheric Composition Forecasts

[Dataset Link](https://ads.atmosphere.copernicus.eu/datasets/cams-global-atmospheric-composition-forecasts?tab=overview)

Using this dataset you can download **Global Forecast data**.

This dataset contains both single-level and multi-level data.

hourly and 3hourly data is available in this dataset. In your configuration file, you must set `resolution = hourly` or `resolution = 3hourly`, depending on whether the data is single-level or multi-level. Otherwise, the download will fail.

### Mandatory fields for global forecast single-level data

```ini
model = cams_forecast
domain = global
resolution = hourly
```

### Mandatory fields for global forecast multi-level data

```ini
model = cams_forecast
domain = global
resolution = 3hourly
```

### Available variables (Species)

Whether the data is single or multi depends on the selected species. These are the available species:

#### Single-level:

* `absod1020aero`
* `absod440aero`
* `absod550aero`
* `acprec`
* `asy1020aero`
* `asy440aero`
* `cldbot`
* `od1020aero`
* `od380aero`
* `od440aero`
* `od500aero`
* `od500aerofine`
* `od550aero`
* `od550du`
* `pm1`
* `pm10`
* `pm2p5`
* `pshltr`
* `rh2`
* `sca1020aero`
* `sca440aero`
* `sca550aero`
* `si`
* `slp`
* `sst`
* `t2`
* `td2`
* `vdist`

#### Multi-level:

* `cfracmax`
* `lbsco532`
* `sconcc2h4`
* `sconcc2h6`
* `sconcc3h6`
* `sconcc3h8`
* `sconcch4`
* `sconcco`
* `sconcetoh`
* `sconcglyox`
* `sconchcho`
* `sconchcl`
* `sconchf`
* `sconchno3`
* `sconcisop`
* `sconcmsa`
* `sconcnh3`
* `sconcnh4`
* `sconcno`
* `sconcno2`
* `sconcno3`
* `sconco3`
* `sconcpan`
* `sconcpb`
* `sconcso2`
* `sconcso4`

> Providentia can only read species in GHOST format. If you want to know the mapping from CAMS variables to GHOST species, please refer to the [CAMS-GHOST species mapping](#cams-ghost-species-mapping) section.

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
model = cams_forecast
domain = global
resolution = hourly
dl_interpolated = False
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
model = cams_forecast
domain = global
resolution = 3hourly
dl_interpolated = False
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

Only hourly data is available in this dataset. You must set `resolution = hourly` in your configuration file; otherwise, the download will not work.

### Mandatory fields for regional reanalysis data

```ini
model = cams_reanalysis_<model>_<stream>
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

* `pm10`
* `pm10du`
* `pm10so4ss`
* `pm2p5`
* `pm2p5ec`
* `sconcco`
* `sconcglyox`
* `sconchcho`
* `sconcnh3`
* `sconcnmvoc`
* `sconcno`
* `sconcno2`
* `sconco3`
* `sconcso2`

> Providentia can only read species in GHOST format. If you want to know the mapping from CAMS variables to GHOST species, please refer to the [CAMS-GHOST species mapping](#cams-ghost-species-mapping) section.

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
model = cams_reanalysis_ensemble_interim
domain = regional
resolution = hourly
dl_interpolated = False
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
model = cams_reanalysis_ensemble_validated
domain = regional
resolution = hourly
dl_interpolated = False
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

Only 3hourly data is available in this dataset. You must set `resolution = 3hourly` in your configuration file; otherwise, the download will not work.

### Mandatory fields for regional reanalysis data

```ini
model = cams_reanalysis
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

* `cfracmax`
* `sconcc2h4`
* `sconcc2h6`
* `sconcc3h6`
* `sconcc3h8`
* `sconcco`
* `sconcetoh`
* `sconchcho`
* `sconchno3`
* `sconcisop`
* `sconcmsa`
* `sconcnh3`
* `sconcnh4`
* `sconcno`
* `sconcno2`
* `sconcno3`
* `sconco3`
* `sconcpan`
* `sconcpb`
* `sconcso2`
* `sconcso4`

> Providentia can only read species in GHOST format. If you want to know the mapping from CAMS variables to GHOST species, please refer to the [CAMS-GHOST species mapping](#cams-ghost-species-mapping) section.

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
model = cams_reanalysis
domain = global
resolution = 3hourly
dl_interpolated = False
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
model = cams_reanalysis
domain = global
resolution = 3hourly
dl_interpolated = False
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

## CAMS-GHOST species mapping

Here is the mapping between CAMS variables and GHOST species names:

```
2m_dewpoint_temperature: td2
2m_temperature: t2
ammonia: sconcnh3
ammonium: sconcnh4
asymmetry_factor_1020nm: asy1020aero
asymmetry_factor_440nm: asy440aero
attenuated_backscatter_due_to_aerosol_532nm_from_ground: lbsco532
carbon_monoxide: sconcco
cloud_base_height: cldbot
dust: pm10du
dust_aerosol_optical_depth_550nm: od550du
ethane: sconcc2h6
ethanol: sconcetoh
ethene: sconcc2h4
formaldehyde: sconchcho
fraction_of_cloud_cover: cfracmax
glyoxal: sconcglyox
hydrogen_chloride: sconchcl
hydrogen_fluoride: sconchf
isoprene: sconcisop
lead: sconcpb
mean_sea_level_pressure: slp
methane: sconcch4
methane_sulfonic_acid: sconcmsa
nitrate: sconcno3
nitric_acid: sconchno3
nitrogen_dioxide: sconcno2
nitrogen_monoxide: sconcno
non_methane_vocs: sconcnmvoc
ozone: sconco3
particulate_matter_10um: pm10
particulate_matter_1um: pm1
particulate_matter_2.5um: pm2p5
peroxyacetyl_nitrate: sconcpan
pm10_sea_salt_dry: pm10so4ss
pm2.5_ammonium: pm2p5nh4
pm2.5_nitrate: pm2p5no3
pm2.5_sulphate: pm2p5so4
propane: sconcc3h8
propene: sconcc3h6
relative_humidity: rh2
sea_surface_temperature: sst
single_scattering_albedo_1020nm: sca1020aero
single_scattering_albedo_440nm: sca440aero
single_scattering_albedo_550nm: sca550aero
snow_depth: si
sulphate_aerosol_mixing_ratio: sconcso4
sulphur_dioxide: sconcso2
surface_pressure: pshltr
total_absorption_aerosol_optical_depth_1020nm: absod1020aero
total_absorption_aerosol_optical_depth_440nm: absod440aero
total_absorption_aerosol_optical_depth_550nm: absod550aero
total_aerosol_optical_depth_1020nm: od1020aero
total_aerosol_optical_depth_380nm: od380aero
total_aerosol_optical_depth_440nm: od440aero
total_aerosol_optical_depth_500nm: od500aero
total_aerosol_optical_depth_550nm: od550aero
total_elementary_carbon: pm2p5ec
total_fine_mode_aerosol_optical_depth_500nm: od500aerofine
total_precipitation: acprec
visibility: vdist
```
