# ERA5

Providentia's download mode supports downloading model output data provided by **ERA5** from the **Climate Data Store** (CDL):

[https://ads.atmosphere.copernicus.eu/datasets](https://ads.atmosphere.copernicus.eu/datasets)

and the **Simulation and Data Laboratory 'Climate Science'** (SDL) at the Jülich Supercomputing Centre:

[https://datapub.fz-juelich.de/slcs/](https://datapub.fz-juelich.de/slcs/)

An ECMWF account is required to download the files. For a tutorial on how to create and set up your account, see the [ECMWF Account setup](#ADS) page.

Data requests may take some time to complete. After starting a request from Providentia, you can monitor its progress by following the steps described in the [Check ECMWF download](Check-CAMS-download) page.

## Available ERA5 datasets

Currently, two ERA5 datasets are available for download:

1. [**ERA5 hourly data on single levels from 1940 to present**](#1-era5-hourly-data-on-single-levels-from-1940-to-present) (CDL)
2. [**Reanalysis Tropopause Data Repository**](#2-reanalysis-tropopause-data-repository) (SDL)

(how-to-enable-era5-download)=
## How to enable ERA5 download 

To download ERA5 data, include one of the following model names along with the domain:

* `era5_reanalysis-global` - [**ERA5 hourly data on single levels from 1940 to present**](#1-era5-hourly-data-on-single-levels-from-1940-to-present)
* `era5_tropopause-global` - [**Reanalysis Tropopause Data Repository**](#2-reanalysis-tropopause-data-repository)

Once you have selected the dataset and specified the model and domain, make sure to set:

```ini
dl_interpolated = False
```

or

Answer `n` to the prompt:

_"Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded)"_


## 1. ERA5 hourly data on single levels from 1940 to present

[Dataset Link](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview)

Using this dataset you can download **Global Reanalysis data**.

This dataset contains single-level data.

Only hourly data is available in this dataset. You must set `resolution = hourly` in your configuration file; otherwise, the download will not work.

### Mandatory fields for regional reanalysis data

```ini
model = era5_reanalysis
domain = global
resolution = hourly
```

### Available variables (Species)

These are the available species:

* `acprec`
* `acsnow`
* `aswin`
* `aswincs`
* `aswindir`
* `aswindircs`
* `aswtoa`
* `blh`
* `cld`
* `cldbot`
* `clddf`
* `cldtop`
* `dir10`
* `gust10`
* `photi`
* `pshltr`
* `si`
* `slp`
* `spd10`
* `sst`
* `t2`
* `td2`
* `u10`
* `v10`

> Providentia can only read species in GHOST format. If you want to know the mapping from CDS ERA5 variables to GHOST species, please refer to the [CDS ERA5-GHOST species mapping](#cds-era5-ghost-species-mapping) section.

### Product Species

There are five derived variables computed from two ERA5 variables. When downloading one of these, the two original variables are also downloaded, formatted and saved together with the final result.

#### cld

### Product Species

There are five derived variables computed from two ERA5 variables. When downloading one of these, the two original variables are also downloaded, formatted, and saved together with the final result.

#### cld

**cld** is computed from **aswin** and **aswincs**, which in CDS are:

- surface_solar_radiation_downwards  
- surface_solar_radiation_downward_clear_sky  

When creating **cld**, both variables are also downloaded and saved.

**Formula:**

\[
\text{cld} = \frac{\text{aswin}}{\text{aswincs}}
\]

#### clddf

**clddf** is computed from **aswindir** and **aswindircs**, which in CDS are:

- total_sky_direct_solar_radiation_at_surface  
- clear_sky_direct_solar_radiation_at_surface  

When creating **clddf**, both variables are also downloaded and saved.

**Formula:**

\[
\text{clddf} = \frac{\text{aswindir}}{\text{aswindircs}}
\]

#### dir10

#### dir10

**dir10** is computed from **u10** and **v10**, which in CDS are:

- 10m_u_component_of_wind  
- 10m_v_component_of_wind  

When downloading **dir10**, both variables are also downloaded and saved.

**Formula:**

\[
\text{dir10} = \left( \frac{180}{\pi} \cdot \arctan2(-u10, -v10) \right)
\]

\[
\text{if } \text{dir10} < 0 \;\Rightarrow\; \text{dir10} = \text{dir10} + 360
\]

#### photi

**photi** is computed from **aswin** and **aswtoa**, which in CDS are:

- surface_solar_radiation_downwards  
- toa_incident_solar_radiation  

When creating **photi**, both variables are also downloaded and saved.

**Formula:**

\[
\text{photi} = \frac{\text{aswin}}{\text{aswtoa}}
\]

#### spd10

**spd10** is computed from **u10** and **v10**, which in CDS are:

- 10m_u_component_of_wind  
- 10m_v_component_of_wind  

When downloading **spd10**, both variables are also downloaded and saved.

**Formula:**

\[
\text{spd10} = \sqrt{u10^2 + v10^2}
\]

## 2. Reanalysis Tropopause Data Repository

[Dataset Link](https://datapub.fz-juelich.de/slcs/tropopause/index.html)

Using this dataset you can download **Global Tropopause Reanalysis data**.

This dataset contains single-level data.

Only hourly data is available in this dataset. You must set `resolution = hourly` in your configuration file; otherwise, the download will not work.

### Mandatory fields for regional reanalysis data

```ini
model = era5_tropopause
domain = global
resolution = hourly
```

### Available variables (Species)

These are the available species:

* `tphclp`
* `tphdyn`
* `tphwmo1`
* `tphwmo2`

> Providentia can only read species in GHOST format. If you want to know the mapping from SDL ERA5 variables to GHOST species, please refer to the [SDL ERA5-GHOST species mapping](#sdl-era5-ghost-species-mapping) section.

## CDS ERA5-GHOST species mapping

Here is the mapping between CDS ERA5 variables and GHOST species names:

```
acprec : total_precipitation
acsnow : snowfall
aswin : surface_solar_radiation_downwards
aswincs : surface_solar_radiation_downward_clear_sky
aswindir : total_sky_direct_solar_radiation_at_surface
aswindircs : clear_sky_direct_solar_radiation_at_surface
aswtoa : toa_incident_solar_radiation
blh : boundary_layer_height
cld : ['aswin', 'aswincs']
cldbot : cloud_base_height
clddf : ['aswindir', 'aswindircs']
cldtop : high_cloud_cover
dir10 : ['u10', 'v10']
gust10 : instantaneous_10m_wind_gust
photi : ['aswin', 'aswtoa']
pshltr : surface_pressure
si : snow_depth
slp : mean_sea_level_pressure
spd10 : ['u10', 'v10']
sst : sea_surface_temperature
t2 : 2m_temperature
td2 : 2m_dewpoint_temperature
u10 : 10m_u_component_of_wind
v10 : 10m_v_component_of_wind
```

## SDL ERA5-GHOST species mapping

Here is the mapping between SDL ERA5 variables and GHOST species names:

```
tphclp: clp_z # cold point height
tphdyn: dyn_z # dynamical tropopause height
tphwmo1: wmo_1st_z # WMO 1st tropopause height
tphwmo2: wmo_2nd_z # WMO 2nd tropopause height
```