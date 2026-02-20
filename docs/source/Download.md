# Overview

Providentia's download mode retrieves modelled and observational data from BSC systems and external sources (CAMS, Zenodo, ACTRIS) for local use.

## Getting started

To start downloading data, simply add `--download` or `--dl` as a launch option along with the **mandatory** configuration file on the command line:

```
./bin/providentia --config='/path/to/file/example.conf' --download
```

```
./bin/providentia --config='/path/to/file/example.conf' --dl
```

This will get the data that needs to be downloaded from your configuration file and save it into the directories specified in `settings/data_paths.yaml` for `local`.

The download mode fetches all the content specified in your configuration file across all sections. To only run one specific section, add the `--section` parameter to the command.

## Types of downloads

Providentia supports four types of downloads. For detailed instructions, please visit the respective pages:  

1. **Download from BSC HPC Machines** 
   - Downloads GHOST and non-GHOST data and model outputs from BSC HPC machines. You must have a BSC account to access this feature.
   - How to get this type of download:  
     - For GHOST networks, answer `y` to the prompt:  
       _Do you want to download observational data from the BSC remote machine? (Otherwise, GHOST observational data will be retrieved from Zenodo)_  
       or set `dl_ghost_source` to `bsc`.  
     - For non-GHOST networks and interpolated/non-interpolated model data, no special action is required.
   - To see more information, check the [BSC download page](BSC-download).

2. **Download of GHOST network data from Zenodo**  
   - Downloads GHOST networks from the [GHOST Zenodo webpage](https://zenodo.org/records/10637450).  
   - How to get this type of download: answer `n` to the HPC prompt: _Do you want to download observational data from the BSC remote machine? (Otherwise, GHOST observational data will be retrieved from Zenodo)_ or set `dl_ghost_source` to `zenodo`.
   - To see more information, check the [Zenodo download page](Zenodo-download).

3. **Download of ACTRIS network data from Thredds**  
   - Downloads observational networks from [ACTRIS Thredds](https://thredds.nilu.no/thredds/catalog.html).  
   - How to get this type of download: write `actris/actris` on the `network` field in your configuration file.
   - To see more information, check the [ACTRIS download page](ACTRIS-download).

4. **Download of non-interpolated model data from the Atmosphere Data Store (ADS)**  
   - Downloads model outputs from the [Atmosphere Data Store](https://ads.atmosphere.copernicus.eu/datasets). You must have an ECMWF account to access this feature.
   - How to get this type of download: specify the model as `cams_analysis`, `cams_forecast` or `cams_reanalysis` in your configuration, and set `dl_interpolated` to `False`.
   - To see more information, check the [CAMS download page](CAMS-download).

## Download configuration fields

All parameters that can be used in the download configuration files can be found in the [Shared Parameters](shared-parameters) or [Download Parameters](download-parameters) sections of the Configuration Fields page.

### Automation of the download

When running downloads, the questions presented during a download can be skipped by setting the appropriate variables. This allows downloads to be fully automated without any user input.

Each of these variables corresponds directly to one of the questions asked during a manual download. 

| Variable | Original Question | Expected Values |
|-------------------- | -------------------- | -------------------- |
| `dl_overwrite`     | _There are some files that were already downloaded in a previous download, do you want to overwrite them ([y]/n)?_ | `True` (overwrite existing files) or `False` (keep existing files)   |
| `dl_ghost_source`  | _Do you want to download observational data from the BSC remote machine? (Otherwise, GHOST observational data will be retrieved from Zenodo) ([y]/n)_ | `bsc` (download from BSC remote machine) or `zenodo` (retrieve from Zenodo)  |
| `dl_interpolated`  | _Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded) ([y]/n)_ | `True` (download interpolated) or `False` (download non-interpolated)                                |
| `dl_mode` | _Which type of data do you want to download? Observational, modelled or both? ([both]/obs/mod)_ | `obs` (download observations), `mod` (download models) or `both` (download both) |
| `dl_thredds_update` | _File containing information of the files available in Thredds for {actris_parameter} ({info_path}) already exists. Do you want to update it (y/[n])?_ | `True` or `False`     |                 |
| `network_type` | _Do you want to download all the GHOST networks? (Otherwise all the non-GHOST networks will be downloaded) ([y]/n)_   | `ghost` (use all GHOST networks) or `non-ghost` (use all non-GHOST networks)                |
## Using wildcards

You can use the `*` wildcard in the following fields to automatically select all available values:

- `network`, `observation`, `framework` 
- `model`, `models`, `experiments`, `experiment`  
- `species`  
- `resolution`  
- `start_date`
- `end_date`  

**Note:** Using wildcards may result in large downloads, so use with caution.