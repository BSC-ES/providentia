# BDRC

Providentia's download mode supports downloading model output data from the **Barcelona Dust Regional Center** (BDRC) through its THREDDS Data Server.

[https://dust.aemet.es/thredds/catalog/dataRoot/catalog.html](https://dust.aemet.es/thredds/catalog/dataRoot/catalog.html)

BDRC credentials are required to download the files. To request credentials for data download please fill in the [Contact form](https://dust.aemet.es/contact-info) with a short bio and briefly explain your motivation. 

## Available BDRC datasets

Currently, only one BDRC dataset is available for download:

1. [**BDRC THREDDS Public Data**](#bdrc-thredds-public-data)

## How to enable BDRC download 

To download data from the Barcelona Dust Regional Center, include the following model name along with the domain. For details about available models,  see the section below:

* `bdrc_<model>-regional` - [**BDRC THREDDS Public Data**](#bdrc-thredds-public-data)

Once you have specified the model and domain, make sure to set:

```ini
dl_interpolated = False
```

or

Answer `n` to the prompt:

_"Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded)"_

## BDRC THREDDS Public Data

[Dataset Link](https://dust.aemet.es/thredds/catalog/dataRoot/catalog.html)

Using this dataset you can download **Regional data**.

Only 3hourly regional data is available in this dataset. You must set `resolution = 3hourly` and `domain = regional` in your configuration file; otherwise, the download will not work.

### Mandatory fields

```ini
model = bdrc_<model>
domain = regional
resolution = 3hourly
```

### Available models

Use the following model names:

* `aladin`
* `bsc_dream8b`
* `cams_ifs`
* `dream8_cams`
* `ema_regcm4`
* `icon_art`
* `lotos_euros`
* `metoffice_um`
* `mocage`
* `monarch`
* `multi_model`
* `nasa_geos`
* `ncep_gefs`
* `noa`
* `silam`
* `wrf_nemo`
* `zamg_wrf_chem`

### Available variables (Species)

The only available species is:

* `od550du`
* `sconcdu`

You may also specify related species in the [mapping_species.yaml](interpolation-mapping-species) file (e.g. `od550aero`, `od500aerocoarse` or other mapped aliases). Regardless of the requested alias, the downloaded data will always correspond to `od550du`.

## Example configuration files

#### Example 1: ZAMG WRF CHEM MODEL
```ini
[BDRC-zamg_wrf_chem] 
model = bdrc_zamg_wrf_chem
species = od550aero
resolution = 3hourly
start_date = 20230201
end_date = 20230202
domain = regional
dl_interpolated = False
```

#### Example 2: WRF NEMO MODEL
```ini
[BDRC-wrf_nemo] 
model = bdrc_wrf_nemo
species = od550aero
resolution = 3hourly
start_date = 20230211
end_date = 20230212
domain = regional
dl_interpolated = False
```

## .env file

An `.env` file will appear in the Providentia root directory when using the download mode. It is designed to store specific user preferences.

   - **BDRC_USER:** This setting specifies the username used to connect to the BDRC. It can be any valid username.
   - **BDRC_PWD:** This setting allows you to save the password needed for connecting to BDRC 

These values can be changed directly on the `.env` file and also be updated by Providentia during the next run.
