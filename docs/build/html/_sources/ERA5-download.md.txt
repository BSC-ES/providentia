# ERA5

Providentia's download mode supports downloading model output data provided by **ERA5** from the **Climate Data Store** (CDS):

[https://ads.atmosphere.copernicus.eu/datasets](https://ads.atmosphere.copernicus.eu/datasets)

and the **Simulation and Data Laboratory 'Climate Science'** (SDL) at the Jülich Supercomputing Centre:

[https://datapub.fz-juelich.de/slcs/](https://datapub.fz-juelich.de/slcs/)

An ECMWF account is required when downloading ERA5 data from the CDS. For a tutorial on how to create and set up your account, see the [Start downloading ECMWF data](#Start-downloading-ECMWF-data) page.

Data requests may take some time to complete. After starting a request from Providentia, you can monitor its progress by following the steps described in the [Check ECMWF download](#Check-ECMWF-download) page.

## Available ERA5 datasets

Currently, two ERA5 datasets are available for download:

1. [**ERA5 hourly data on single levels from 1940 to present**](#1-era5-hourly-data-on-single-levels-from-1940-to-present) (CDS)
2. [**Reanalysis Tropopause Data Repository**](#2-reanalysis-tropopause-data-repository) (SDL)

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

(download-a-specific-area-era5)=
## Download a specific area
If you do not want to download the entire NetCDF file, you can download data for a specific geographic area by including `longitude` and `latitude` as configuration fields.

For example, to download data covering Europe in Providentia:

```ini
longitude = -28, 53
latitude = 35, 72
```

The equivalent area for the request would be the following:

```
{
    ...
    "area": [72, -28, 35, 53],
    ...
}
```

## 1. ERA5 hourly data on single levels from 1940 to present

[Dataset Link](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview)

Using this dataset you can download **Global Reanalysis data**.

This dataset contains single-level data.

Only hourly global data is available in this dataset. You must set `resolution = hourly` and `domain = global` in your configuration file; otherwise, the download will not work.

### Mandatory fields

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
* `auvin`
* `blh`
* `cldaf`
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

### Dataset temporal coverage

The dataset covers the period from January 1940 to the present.

### Fixed Download Settings

Providentia assumes the following fixed values when downloading data:

* `product_type = reanalysis`

### Product Species

There are five derived variables computed from two ERA5 variables. When downloading one of these, the two original variables are also downloaded, formatted and saved together with the final result.

#### cldaf

**cldaf** is computed from **aswin** and **aswincs**, which in CDS are:

- surface_solar_radiation_downwards  
- surface_solar_radiation_downward_clear_sky  

When creating **cldaf**, both variables are also downloaded and saved.

**Formula:**

$$
\begin{aligned}
\text{cld} &= \frac{\text{aswin}}{\text{aswincs}}
\end{aligned}
$$

#### clddf

**clddf** is computed from **aswindir** and **aswindircs**, which in CDS are:

- total_sky_direct_solar_radiation_at_surface  
- clear_sky_direct_solar_radiation_at_surface  

When creating **clddf**, both variables are also downloaded and saved.

**Formula:**

$$
\begin{aligned}
\text{clddf} &= \frac{\text{aswindir}}{\text{aswindircs}}
\end{aligned}
$$

#### dir10

**dir10** is computed from **u10** and **v10**, which in CDS are:

- 10m_u_component_of_wind  
- 10m_v_component_of_wind  

When downloading **dir10**, both variables are also downloaded and saved.

**Formula:**

$$
\begin{aligned}
\text{dir10} &= \frac{180}{\pi}\arctan2(-u10,\,-v10)
\end{aligned}
$$

$$
\begin{aligned}
\text{if } \text{dir10} &< 0 \\
\text{dir10} &= \text{dir10} + 360
\end{aligned}
$$

#### photi

**photi** is computed from **aswin** and **aswtoa**, which in CDS are:

- surface_solar_radiation_downwards  
- toa_incident_solar_radiation  

When creating **photi**, both variables are also downloaded and saved.

**Formula:**

$$
\begin{aligned}
\text{photi} &= \frac{\text{aswin}}{\text{aswtoa}}
\end{aligned}
$$

#### spd10

**spd10** is computed from **u10** and **v10**, which in CDS are:

- 10m_u_component_of_wind  
- 10m_v_component_of_wind  

When downloading **spd10**, both variables are also downloaded and saved.

**Formula:**

$$
\begin{aligned}
\text{spd10} &= \sqrt{u10^2 + v10^2} 
\end{aligned}
$$

## Example configuration file

```ini
[ERA5-CDS] 
start_date = 20220101
end_date = 20220201
species = t2
model = era5_reanalysis-global
resolution = hourly
dl_interpolated = False
```

## 2. Reanalysis Tropopause Data Repository

[Dataset Link](https://datapub.fz-juelich.de/slcs/tropopause/index.html)

Using this dataset you can download **Global Tropopause Reanalysis data**.

This dataset contains single-level data.

Only hourly global data is available in this dataset. You must set `resolution = hourly` and `domain = global` in your configuration file; otherwise, the download will not work.

### Mandatory fields

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

### Dataset temporal coverage

The dataset covers the period from January 2000 to the present.

### Fixed Download Settings

Providentia downloads tropopause data from the `v2/era5` directory.

## Example configuration file

```ini
[ERA5-SDL] 
start_date = 20220101
end_date = 20220201
species = tphclp
model = era5_tropopause-global
resolution = hourly
dl_interpolated = False
```

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
auvin: downward_uv_radiation_at_the_surface
blh : boundary_layer_height
cldaf : [aswin, aswincs]
cldbot : cloud_base_height
clddf : [aswindir, aswindircs]
cldtop : high_cloud_cover
dir10 : [u10, v10]
gust10 : instantaneous_10m_wind_gust
photi : [aswin, aswtoa]
pshltr : surface_pressure
si : snow_depth
slp : mean_sea_level_pressure
spd10 : [u10, v10]
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