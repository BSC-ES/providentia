# Configuration fields

Configuration fields determine how Providentia behaves during a run using a configuration file. Below is a full list of all available parameters organised by mode.

These fields can also be set via command line arguments, for more details, see the [Command line configuration page](Command-line-configuration).

## Common parameters

Some of these fields are required depending on the mode. If a parameter required by a given mode is missing, Providentia will fail.

| Parameter | Required in | Description |
| ------ | ------ | ------ |
| network | Dashboard, Report, Interpolation | Network you want to load observations from. Can be multiple (e.g. `CAPMoN, EBAS`). Adding a wild card (\*) is going to expand to certain variables (vconcaerobin* → vconcaerobin1, vconcaerobin2, etc.). |
| species | All modes | Species to load. Can be multiple (e.g. `sconco3, sconcno2`). |
| resolution | All modes | Temporal resolution of the observations you want to load (e.g. `3hourly`). |
| start_date | All modes | Comparison start date in YYYYMMDD format (e.g. `20170101`). |
| end_date | All modes | Comparison end date in YYYYMMDD format (e.g. `20180601`). |
| model | Dashboard, Report, Interpolation | ID of model. The model IDs can be mapped to different names by adding a list of alternative names after the model IDs (e.g. `mod1, mod2 (altmod1, altmod2)`). After interpolation model IDs will consist of 3 separate strings [ID-domain-ensemble]. The domain and ensemble  |
| model_resolution | Optional | Model resolution if different from observations | 
| domain | Optional | Domain of the model. Can be multiple. |
| ensemble | Optional | Ensemble member number or ensemble statistic of the model. defaults to all members available. Can be multiple. |  
| forecast | Optional | Controls how forecast data is handled. Can be multiple (e.g. `day`, `daily`, `combined`). If wanting to limit to specific forecast days, then add the forecast day to the option, (e.g. `day1`, `daily2`, `combined3`). Different options cannot be mixed however (i.e. `day` options cannot be set with `daily` or `combined` options). This variable must be set to a valid value when performing interpolation for forecast data to be interpolated.|
| filter_species | Optional | Filter read species by other species data within a data range. The first value set is the lower bound to filter by, and the second value the upper bound. Place a sign before each bound value to inform if the filter should be inclusive or exclusive of the bound, e.g. `<` or `<=`. If not wishing to set either the lower or upper bounds, a `:` can be used. Optionally, a fill value can also be given as a third value to impose what the filtered data is set to, by default this is `NaN`. Multiple filters can be  set together separated by a comma (e.g. `network1:species1 (>lowerlim, <=upperlim, fillvalue), network2:species2 (:, <upperlim)`). |

## Analysis and visualisation modes (Dashboard, Report, Library)

Apart from the common parameters, these are the fields used by all analysis and visualisation modes ([Dashboard](Dashboard), [Report](Report), [Library](Library)). All parameters in this section are optional.

| Parameter | Description |
| ------ | ------ |
| statistic_mode | Statistic mode: Temporal&#124;Spatial (default), Spatial&#124;Temporal or Flattened. |
| statistic_aggregation | Aggregation statistic, e.g. Median. |
| timeseries_statistic_aggregation | Timeseries aggregation statistic, e.g. Median. |
| periodic_statistic_mode | Periodic statistic mode: Independent (default), Cycle. |
| periodic_statistic_aggregation | Periodic aggregation statistic, e.g. Mean (default). |
| temporal_colocation | Boolean variable to set if you want to temporally colocate the observation and model data. |
| spatial_colocation | Boolean variable to set if you want to spatially colocate the observation and model data across multiple species. |
| plot_characteristics_filename | The path to the file containing the plot characteristics. |
| observations_data_label | Alias for observational data |
| lower_bound | Filter out data lower than this set limit. If multiple species are being read then this can either be one value, setting the same limit across species, or multiple values per species (e.g. `3, 4, 5`). |
| upper_bound | Filter out data above this set limit. If multiple species are being read then this can either be one value, setting the same limit across species, or multiple values per species (e.g. `3, 4, 5`). |
| map_extent | Set the map plot extents with the syntax: minimum longitude, maximum longitude, minimum latitude, maximum latitude. |
| remove_extreme_stations | Type of extreme stations removal, from the options given in `remove_extreme_stations.yaml`. |
| resampling_resolution | Resolution you want to resample your data to. Options: `hourly`, `3hourly`, `6hourly`, `daily`, `monthly`, `annual`.|

## Dashboard 

This parameter is used only in the [Dashboard mode](Dashboard). It is optional.

| Parameter | Description |
| ------ | ------ |
| active_dashboard_plots | Plots that will be active in the dashboard once it is launched (e.g. `timeseries, periodic-violin, scatter, distribution`). |

## Report 

These parameters are used only in the [Report mode](Report). All of them are optional.

| Parameter | Description |
| ------ | ------ |
| report_type | Type of report to generate that defines which plots the report will contain, from the options given in `report_plots.yaml`. |
| report_summary | Boolean variable to set if you wish to make specific plots for each station in subsection. |
| report_stations | Boolean variable to set if you wish to make summary plots across station subsection. |
| report_title | The header in the first page of the report (as in the PDF). |
| report_filename | The filename of the report or the path to create the report (as in the PDF). |
| harmonise_stations | Boolean variable to set if you wish to harmonise axes limits across stations for stations report. |
| harmonise_summary | Boolean variable to set if you wish to harmonise axes limits across subsections for summary report. |

If the number of networks and species are **both** multiple but not equal, Providentia will throw the error `Error: The number of "network" and "species" fields is not the same.` and the user will be required to clearly specify which networks and species they want. For example, this would not be accepted:

```
network = EBAS, EEA_AQ_eReporting
species = sconco3, sconcno2, sconcso2
```

But this would:
```
network = EBAS, EBAS, EBAS, EEA_AQ_eReporting, EEA_AQ_eReporting, EEA_AQ_eReporting
species = sconco3, sconcno2, sconcso2, sconco3, sconcno2, sconcso2
```

## Interpolation 

These parameters are used only in the [Interpolation mode](Interpolation). All of them are optional.

| Parameter | Description |
| ------ | ------ |
| interp_n_neighbours | The number of nearest neighbours to use in the interpolation of model output to observational stations. If not set, this defaults to `4`.|
| interp_spinup_timesteps | Number of initial timesteps skipped for model spin-up. If not set, this defaults to `0`.|
| interp_model_downsampling | Sets the statistic for the downsampling of the model resolution to the observational resolution. Options: `mean`, `median`.|
| interp_model_upsampling | Sets the method for the upsampling of the model resolution to the observational resolution. `fill` linearly fills between measurements, and `gaps` sets NaN values for times that the model does not have.|
| interp_multiprocessing | Boolean variable to set if you wish to use multiprocessing instead of greasy to interpolate on HPC machines. |
| network_type | Determines whether to use all GHOST or all non-GHOST networks when the observation field uses the `*` wildcard.  | 

## Download

These parameters are used only in the [Download mode](Download). All of them are optional.

| Parameter | Description |
| ------ | ------ |
| dl_overwrite    | Indicates whether previously downloaded files should be overwritten. | _There are some files that were already downloaded in a previous download, do you want to overwrite them ([y]/n)?_ | 
| dl_ghost_source  | Determines where GHOST observations are downloaded from. |
| dl_interpolated  | Specifies whether the interpolated versions of the model output should be downloaded. |  
| dl_mode | Selects what to download when both observations and model output are present in the configuration file. | 
| network_type | Determines whether to use all GHOST or all non-GHOST networks when the observation field uses the `*` wildcard.  | 