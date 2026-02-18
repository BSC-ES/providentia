# Configuration fields

Configuration fields determine how Providentia behaves during a run using a configuration file. Below is a full list of all available parameters organised by mode.

These fields can also be set via command line arguments, for more details, see the [Command line configuration page](Command-line-configuration).

(shared-parameters)=
## Shared parameters

Some of these fields are required depending on the mode. If a parameter required by a given mode is missing, Providentia will fail.

| Parameter | Required in | Description | Default
| ------ | ------ | ------ | ------ |
| `network`, `observation`, `framework` | Dashboard, Report, Interpolation | Network(s) to load observations from. Multiple values are allowed (e.g. `CAPMoN, EBAS`). Wildcards (`*`) expand to multiple variables (e.g. `vconcaerobin*` → `vconcaerobin1`, `vconcaerobin2`, etc.). For GHOST networks, the selection is dictated by GHOST. For non-GHOST networks, available options are defined under the `nonghost_available_networks` key in `available_inputs.yaml` and can be modified by the user. | — |
| `model`, `models`, `experiments`, `experiment` | Dashboard, Report, Interpolation | Model(s). Model IDs can optionally include domain and/or ensemble information in the following formats: `modelID`, `modelID-ensemble`, `modelID-domain` or `modelID-domain-ensemble`. Also, `model`, `domain` and/or `ensemble` can be specified separately. Model IDs can also be mapped to alternative names (aliases) by appending them in parentheses after the ID (e.g. `mod1-dom-ens, mod2-dom-ens (altmod1, altmod2)`). | — |
| `species` | All modes | Species to load. Can be multiple (e.g. `sconco3, sconcno2`). Dictated by GHOST. See the [Available Species](Available-species) page for options. | — |
| `resolution` | All modes | Temporal resolution of the observations to load (e.g. `hourly`, `daily`). For GHOST networks, the resolution is dictated by GHOST. For non-GHOST networks, available options are defined under the `nonghost_available_resolutions` key in `available_inputs.yaml` and can be modified by the user. || — |
| `start_date` | All modes | Comparison start date in `YYYYMMDD` format or `YYYYMM` when interpolation is enabled (e.g. `20170101`). | — |
| `end_date` | All modes | Comparison end date in `YYYYMMDD` format or `YYYYMM` when interpolation is enabled (e.g. `20180601`). | — |
| `ghost_version` | Optional | GHOST version used when a GHOST network is selected. | `1.5` |
| `domain` | Optional | Domain of the model (e.g. `regional`, `global`). When multiple model IDs and multiple ensembles/domains are provided, all possible combinations of model, domain and ensemble will be used. Options are defined under the `available_domains` key in `available_inputs.yaml` and can be modified by the user. | All available |
| `ensemble` | Optional | Ensemble of the model (e.g. `000`, `001`). When multiple model IDs and multiple ensembles/domains are provided, all possible combinations of model, domain and ensemble will be used. | All available, except in interpolation mode where the default is `000` |
| `forecast` | Optional | Controls how forecast data is handled: `day`, `daily`, `combined`, `dayN`, `dailyN`, `combinedN`. Multiple values of the same type can be provided (e.g. `day`, `day2`, `day3`), but different types cannot be mixed (i.e. `day` options cannot be combined with `daily` or `combined`). To limit to specific forecast days, append the day number to the option (e.g. `day1`, `daily2`, `combined3`). This variable must be set to a valid value when performing interpolation for forecast data. | — |
| `filter_species` | Optional | Filter read species by other species data within a data range. The first value set is the lower bound to filter by, and the second value the upper bound. Place a sign before each bound value to inform if the filter should be inclusive or exclusive of the bound, (e.g. `<` or `<=`). If not wishing to set either the lower or upper bounds, a `:` can be used. Optionally, a fill value can also be given as a third value to impose what the filtered data is set to, by default this is `NaN`. Multiple filters can be  set together separated by a comma (e.g. `network1:species1 (>lowerlim, <=upperlim, fillvalue), network2:species2 (:, <upperlim)`). | — |
| `ghost_root` | Optional | Root directory for GHOST observations, overwrites `data_paths.yaml` | From `data_paths.yaml` |
| `nonghost_root` | Optional | Root directory for non-GHOST observations, overwrites `data_paths.yaml` | From `data_paths.yaml` |
| `mod_root` | Optional | Root directory for interpolated model data, overwrites `data_paths.yaml` | From `data_paths.yaml` |
| `mod_to_interp_root` | Optional | Root directory for non-interpolated model data, overwrites `data_paths.yaml` | From `data_paths.yaml` |
| `config_dir` | Optional | Path to all configuration files. | `configurations/` |
| `cartopy_data_dir` | Optional | Cartopy data directory. | In HPC: `/gpfs/projects/bsc32/software/rhel/9.2/software/Cartopy/0.23.0-foss-2023b-Python-3.11.5/lib/python3.11/site-packages/cartopy/data`. In local: Downloaded from the internet on the fly. |

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

(visualization-parameters)=
## Parameters for analysis and visualization modes (Dashboard, Report, Library)

Apart from the common parameters, these are the fields used by all analysis and visualisation modes ([Dashboard](Dashboard), [Report](Report), [Library](Library)). All parameters in this section are **optional**.

| Parameter | Description | Default |
| ------ | ------ | ------ |
| `statistic_mode` | Statistic mode: `Temporal\|Spatial`, `Spatial\|Temporal`, `Flattened`. | `Temporal\|Spatial` |
| `statistic_aggregation` | Aggregation statistic: `Median`, `Mean`, `p1`, `p5`, `p10`, `p25`, `p75`, `p90`, `p95`, `p99`. | Depends on `statistic_mode`: `Median` if `Temporal\|Spatial` or `Spatial\|Temporal`; no aggregation if `Flattened` |
| `timeseries_statistic_aggregation` | Timeseries aggregation statistic: `Median`, `Mean`, `p1`, `p5`, `p10`, `p25`, `p75`, `p90`, `p95`, `p99`. | `Median` |
| `periodic_statistic_mode` | Periodic statistic mode: `Independent`, `Cycle`. | `Independent` |
| `periodic_statistic_aggregation` | Periodic aggregation statistic: `Median`, `Mean`, `p1`, `p5`, `p10`, `p25`, `p75`, `p90`, `p95`, `p99`. | `Median` |
| `temporal_colocation` | Boolean variable to set if you want to temporally colocate the observation and model data. | `True` |
| `temporal_colocation_active` | Boolean variable | `False` |
| `spatial_colocation` | Boolean variable to set if you want to spatially colocate the observation and model data across multiple species. | `True` |
| `spatial_colocation_tolerance` | Spatial colocation tolerance to match stations by `longitudes`/`latitudes` and/or `measurement_altitudes` (in metres) | `19.053` |
| `spatial_colocation_validation` | Boolean variable to validate spatial colocation intersections via position using `spatial_colocation_tolerance` | `True` |
| `spatial_colocation_validation_tolerance` | Spatial colocation validation tolerance to validate station reference/station name match of stations by longitude/latitude position (in metres) | `10000.0` |
| `spatial_colocation_station_reference` | Boolean variable to indicate the use of `station_reference` variable for spatial colocation |  `True` |
| `spatial_colocation_station_name` | Boolean variable to indicate usage of `station_name` variable for spatial colocation |  `True` |
| `spatial_colocation_longitude_latitude` | Boolean variable to indicate the use of `longitude` and `latitude` variables for spatial colocation |  `True` |
| `spatial_colocation_measurement_altitude` | Boolean variable to indicate the use of `measurement_altitude` variable for spatial colocation |  `True` |
| `plot_characteristics_filename` | The path to the file containing the plot characteristics. | — |
| `observations_data_label` | Alias for observational data | `observations` |
| `lower_bound` | Filter out data lower than this set limit. If multiple species are being read then this can either be one value, setting the same limit across species or multiple values per species (e.g. `3, 4, 5`). | — |
| `upper_bound` | Filter out data above this set limit. If multiple species are being read then this can either be one value, setting the same limit across species or multiple values per species (e.g. `3, 4, 5`). | — |
| `map_extent` | Set the map plot extents with the syntax: minimum longitude, maximum longitude, minimum latitude, maximum latitude (e.g. `-30, 50, 20, 90`). | `[-180, 180, -90, 90]` in Dashboard |
| `remove_extreme_stations` | Type of extreme stations removal, from the options given in `remove_extreme_stations.yaml`. | — |
| `resampling_resolution` | Resolution you want to resample your data to: `hourly`, `3hourly`, `6hourly`, `daily`, `monthly`, `annual`. | — |

(dashboard-parameters)=
## Dashboard parameters

This parameter is used only in the [Dashboard mode](Dashboard). It is **optional**.

| Parameter | Description | Default |
| ------ | ------ | ------ |
| `active_dashboard_plots` | Plots that will be active in the dashboard once it is launched (e.g. `timeseries, periodic-violin, scatter, distribution`). | `timeseries, statsummary, distribution, periodic` |

(report-parameters)=
## Report parameters

These parameters are used only in the [Report mode](Report). All of them are **optional**.

| Parameter | Description | Default |
| ------ | ------ | ------ |
| `report_type` | Type of report to generate that defines which plots the report will contain, from the options given in `report_plots.yaml`. | `standard` |
| `report_summary` | Boolean variable to set if you wish to make specific plots for each station in subsection. | `True` |
| `report_stations` | Boolean variable to set if you wish to make summary plots across station subsection. | `False` |
| `report_title` | The header in the first page of the report (as in the PDF). | `Providentia Report` |
| `report_filename` | The filename of the report or the path to create the report (as in the PDF). | `PROVIDENTIA_Report` |
| `harmonise_stations` | Boolean variable to set if you wish to harmonise axes limits across stations for stations report. | `True` |
| `harmonise_summary` | Boolean variable to set if you wish to harmonise axes limits across subsections for summary report. | `True` |

(interpolation-parameters)=
## Interpolation parameters

These parameters are used only in the [Interpolation mode](Interpolation). All of them are **optional**.

| Parameter | Description | Default |
| ------ | ------ | ------ |
| `interp_n_neighbours` | Number of nearest neighbours used for interpolation | `4` |
| `interp_reverse_vertical_orientation` | Reverse vertical order of model levels | `False` |
| `interp_chunk_size` | Minimum number of jobs per interpolation chunk | `16` |
| `interp_job_array_limit` | Maximum number of chunks in the job array | `100` |
| `interp_multiprocessing` | Use multiprocessing instead of Greasy on HPC systems | `False` |
| `interp_spinup_timesteps` | Number of initial timesteps skipped for model spin-up | `0` |
| `interp_model_downsampling` | Statistic for the downsampling of the model resolution to the observational resolution: `mean`, `median`. | `mean` |
| `interp_model_upsampling` | Method for the upsampling of the model resolution to the observational resolution: `fill`, `gaps`. `fill` linearly fills between measurements, and `gaps` sets NaN values for times that the model does not have. | `fill` |
| `network_type` | Determines whether to use all GHOST or all non-GHOST networks when the `observation` field uses the `*` wildcard.  | — |
| `model_resolution` | Optional | Model resolution if different from observations. | Same as `resolution` |

(download-parameters)=
## Download parameters

These parameters are used only in the [Download mode](Download). All of them are **optional**.

| Parameter | Description | Default |
| ------ | ------ | ------ |
| `dl_overwrite` | Indicates whether previously downloaded files should be overwritten: `True`, `False`. | — |
| `dl_ghost_source`  | Determines where GHOST observations are downloaded from: `bsc`, `zenodo`. | — |
| `dl_interpolated`  | Specifies whether the interpolated versions of the model output should be downloaded: `True`, `False`. | — |
| `dl_mode` | Selects what to download when both observations and model output are present in the configuration file: `obs`, `mod`, `both`. | — |
| `network_type` | Determines whether to use all GHOST or all non-GHOST networks when the observation field uses the `*` wildcard: `ghost`, `non-ghost`.  | — |
| `dl_timeout` | Sets the timeout (in seconds) for downloads from HPC systems, covering interpolated and non-interpolated model data as well as GHOST and non-GHOST observations. | `180` |
| `model_resolution` | Optional | Model resolution if different from observations. | Same as `resolution` |