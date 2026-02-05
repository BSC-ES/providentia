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

2. **Download of network from Zenodo**  
   - Downloads GHOST networks from the [GHOST Zenodo webpage](https://zenodo.org/records/10637450).  
   - How to get this type of download: answer `n` to the HPC prompt: _Do you want to download observational data from the BSC remote machine? (Otherwise, GHOST observational data will be retrieved from Zenodo)_ or set `dl_ghost_source` to `zenodo`.
   - To see more information, check the [Zenodo download page](Zenodo-download).

3. **Download of network from ACTRIS (Thredds)**  
   - Downloads observational networks from [ACTRIS Thredds](https://thredds.nilu.no/thredds/catalog.html).  
   - How to get this type of download: put `actris/actris` in the `network` field in your configuration.
   - To see more information, check the [ACTRIS download page](ACTRIS-download).

4. **Download of non-interpolated model data from the Atmosphere Data Store (ADS)**  
   - Downloads model outputs from the [Atmosphere Data Store](https://ads.atmosphere.copernicus.eu/datasets).  
   - How to get this type of download: specify the model as `cams_analysis`, `cams_forecast` or `cams_reanalysis` in your configuration, and set `dl_interpolated` to `False`.
   - To see more information, check the [CAMS download page](CAMS-download).

### Download configuration fields

Only the following configuration fields are used during download. All required fields must be provided.

| Variable        | Description                              | Required | Default |
|-----------------|------------------------------------------|----------|---------|
| `ghost_version`         | GHOST version used when a GHOST network is selected | No | 1.5 |
| `network`, `observation`, `framework` | Observation network to use               | Yes      | —       |
| `model`, `models`, `experiments`, `experiment` |  Model ID(s) to be interpolated           | No       | —       |
| `domain` | Domain of the model, can be indicated in the model field (e.g. `regional`, `global`) | No | — |
| `ensemble` | Ensemble of the model, can be indicated in the model field (e.g. `000`, `001`) | No | — |
| `species`         | Species to load (e.g. `sconco3`, `pm10`)| Yes      | —       |
| `resolution`      | Observation data resolution (e.g. `hourly`, `daily`) | Yes | —       |
| `model_resolution` | Model resolution if different from observations | No | Same as `resolution` |
| `start_date`      | Start date of download (`YYYYMM`)   | Yes      | —       |
| `end_date`        | End date of download (`YYYYMM`)     | Yes      | —       |
| `filter_species`  | Optional filter to select specific species | No     | —       |


## Automation of the download

In order to add the download to your scripts or if you just want to make it without the user input, here are all the variables you need to have

| Variable  | Description    | Original Question | Expected Values |
|--------------------|----------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| `dl_overwrite`     | Indicates whether previously downloaded files should be overwritten. | _There are some files that were already downloaded in a previous download, do you want to overwrite them ([y]/n)?_ | `True` (overwrite existing files) or `False` (keep existing files)   |
| `dl_ghost_source`  | Determines where GHOST observations are downloaded from. | _Do you want to download observational data from the BSC remote machine? (Otherwise, GHOST observational data will be retrieved from Zenodo) ([y]/n)_ | `bsc` (download from BSC remote machine) or `zenodo` (retrieve from Zenodo)  |
| `dl_interpolated`  | Specifies whether the interpolated versions of the model output should be downloaded. | _Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded) ([y]/n)_ | `True` (download interpolated) or `False` (download non-interpolated)                                |
| `dl_mode` | Selects what to download when both observations and model output are present in the configuration file. | _Which type of data do you want to download? Observational, modelled or both? ([both]/obs/mod)_      | `obs` (download observations), `mod` (download models) or `both` (download both) |
| `network_type` | Determines whether to use all GHOST or all non-GHOST networks when the observation field uses the `*` wildcard.  | _Do you want to download all the GHOST networks? (Otherwise all the non-GHOST networks will be downloaded) ([y]/n)_   | `ghost` (use all GHOST networks) or `non-ghost` (use all non-GHOST networks)                |

## Using wildcards

You can use the `*` wildcard in the following fields to automatically select all available values:

- `network`, `observation`, `framework` 
- `model`, `models`, `experiments`, `experiment`  
- `species`  
- `resolution`  
- `start_date`
- `end_date`  

**Note:** Using wildcards may result in large downloads, so use with caution.