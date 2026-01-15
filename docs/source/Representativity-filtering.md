# Representativity filtering

One major limitation associated with observations is the amount of gaps often found between measurements. If these observations are directly compared with typically complete model data, this would impose a significant bias upon the comparison. Filtering by representativity filters provides a way to control the temporal robustness of the observations for evaluation. 

Providentia has multiple such filters. By default these are **NOT** activated (i.e. set at 0). 

## Native vs averaged filters 

Providentia has two types of representativity filters, those which calculate the representativity percent across the native resolution (which can be variable), and those which calculate it at the average resolution, e.g. the resolution being used in Providentia such as hourly, 3hourly etc. The native filters are only available when using GJOST data, and are differentiated between, with the native representativity variables all containing the `native` string.

 The first is the `representativity_percent` which is used to set the minimum % of representativity needed for a station's measurements to be valid in a specific time period (i.e. daily, monthly, etc.). If for instance we wanted to check the representativity of hourly station data each day, we would use the `daily_representativity_percent` filter. If we set this at 50%, this would mean any daily periods where there are <50% of observations missing, would be set to be NaN. By default, this is always set at 0%.

If wanting to apply filters to the entire time period in memory, rather than shorter windows, this can be done using the `all_representativity_percent` and `all_max_gap_percent` filters.
Any stations which are completely NaN after filtering are removed from the map.

Whatever the current temporal `resolution` is of the current data in memory, the standard representativity filters will be available for all coarser resolutions, e.g. for hourly resolution data, the standard filters will be available for `daily`, `monthly` and `all` resolutions. 

When using GHOST data there are some extra filters available which can be used to assess the representativity of the native resolution data, rather than solely the AC standard averaged temporal resolutions. For example if measurements are natively measured every 10 minutes and these are then averaged to hourly resolution, then the representativity checks will assess the representativity of this native data within the desired time period. When GHOST data is used there is therefore native versions of the standard filters. Additionally because of the nature of using native data, there are additional filters at the current resolution in use, e.g. for hourly data: `hourly_native_representativity_percent`. These filters can be distinguished from the standard forms by the occurrence of the native word within the filter name.    

## Available filters

The available representativity filters in Providentia (when using GHOST data), per temporal resolution are as follows:

### hourly, hourly_instantaneous

```
hourly_native_representativity_percent = 0
daily_native_representativity_percent = 0
daily_representativity_percent = 0
monthly_native_representativity_percent = 0
monthly_representativity_percent = 0
all_representativity_percent = 0
```

### 3hourly, 6hourly, 3hourly_instantaneous, 6hourly_instantaneous, daily

```
daily_native_representativity_percent = 0
monthly_native_representativity_percent = 0
monthly_representativity_percent = 0
all_representativity_percent = 0
```

### monthly

```
monthly_native_representativity_percent = 0
all_representativity_percent = 0
```