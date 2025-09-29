# CAMS Download

Providentia's download mode supports downloading experiment data from **Atmosphere Data Store** datasets provided by CAMS:

[https://ads.atmosphere.copernicus.eu/datasets](https://ads.atmosphere.copernicus.eu/datasets)

## Available Datasets

Currently, four CAMS datasets are available for download:

1. **CAMS European Air Quality Forecasts**
2. **CAMS Global Atmospheric Composition Forecasts**
3. **CAMS European Air Quality Reanalyses**
4. **CAMS Global Reanalysis (EAC4)**

## 1. CAMS European Air Quality Forecasts

[Dataset Link](https://ads.atmosphere.copernicus.eu/datasets/cams-europe-air-quality-forecasts?tab=download)

You can download either:

* **Analysis data**, or
* **Forecast data**

This dataset provides hourly data only. You must set `resolution = hourly` in your configuration file; otherwise, the download will not work.

### Mandatory Fields for Analysis Data

```ini
experiment = cams_analysis_<model>
domain = regional
resolution = hourly
```

### Mandatory Fields for Forecast Data

```ini
experiment = cams_forecast_<model>
domain = regional
resolution = hourly
```

> Replace `<model>` with one of the available model names listed below.

### Available Models

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

### Available Variables (Species)

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

> Providentia can only read species in GHOST format. The mapping is provided in the [GHOST–CAMS Species Mapping](#ghostcams-species-mapping) section.

### Assumptions

Providentia assumes the following fixed values when downloading data:

* `level = 0`
* `time = 00:00`
* `leadtime_hour = 0-96`

### Example Configuration Files and Corresponding API Requests

#### Example 1: CAMS Analysis (Ensemble)

```ini
[cams_analysis_ensemble-regional]
start_date = 20250701
end_date = 20250801
species = sconcno2
experiments = cams_analysis_ensemble
domain = regional
resolution = hourly
interpolated = false
```

```python
dataset = 'cams-europe-air-quality-forecasts'
request = {
    'variable': 'nitrogen_dioxide',
    'leadtime_hour': '0',
    'type': 'analysis',
    'date': '2025-07-01/2025-07-31',
    'model': 'ensemble',
    'level': '0',
    'time': '00:00',
    'data_format': 'netcdf_zip',
}
```

#### Example 2: CAMS Forecast (Ensemble)

```ini
[cams_forecast_ensemble-regional]
start_date = 20250701
end_date = 20250702
species = sconcno2
experiments = cams_forecast_ensemble
domain = regional
resolution = hourly
interpolated = false
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

## GHOST–CAMS Species Mapping

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
