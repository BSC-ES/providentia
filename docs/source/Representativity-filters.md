# Representativity filters

One major limitation associated with observations is the amount of gaps often found between measurements.

If these observations are directly compared with model data, this would impose a significant bias upon the comparison.   

The representativity filters provides one way to control the temporal robustness of the observations for evaluation. 

Providentia has two types of representativity filters. The first is the `representativity_percent` which is used to set the minimum % of representativity needed for a station's measurements to be valid in a specific time period (i.e. daily, monthly, etc.). If for instance we wanted to check the representativity of hourly station data each day, we would use the `daily_representativity_percent` filter. If we set this at 50%, this would mean any daily periods where there are <50% of observations missing, would be set to be NaN. By default, this is always set at 0%.

The second type of representativity filter is the `max_gap_percent`. This is used to set the **maximum** permitted gap allowed for a station's measurements in a specific time period. For example if we wanted to check the max gap of hourly station data each day, we would use the `daily_max_gap_percent` filter. If we set this at 50%, this would mean any daily periods where there is a continuous gap of >50%, would be set to be NaN.

If wanting to apply filters to the entire time period in memory, rather than shorter windows, this can be done using the `all_representativity_percent` and `all_max_gap_percent` filters.
Any stations which are completely NaN after filtering are removed from the map.

Whatever the current temporal `resolution` is of the current data in memory, the standard representativity filters will be available for all coarser resolutions, e.g. for hourly resolution data, the standard filters will be available for `daily`, `monthly` and `all` resolutions. 

When using GHOST data there are some extra filters available which can be used to assess the representativity of the native resolution data, rather than solely the AC standard averaged temporal resolutions. For example if measurements are natively measured every 10 minutes and these are then averaged to hourly resolution, then the representativity checks will assess the representativity of this native data within the desired time period. When GHOST data is used there is therefore native versions of the standard filters. Additionally because of the nature of using native data, there are additional filters at the current resolution in use, e.g. for hourly data: `hourly_native_representativity_percent`. These filters can be distinguished from the standard forms by the occurrence of the native word within the filter name.    

## Defaults

The default representativity filters set in Providentia (when using GHOST data), per temporal `resolution` are as follows:

### hourly, hourly_instantaneous

```
hourly_native_representativity_percent = 0
hourly_native_max_gap_percent = 100 
daily_native_representativity_percent = 0
daily_representativity_percent = 0
daily_native_max_gap_percent = 100
daily_max_gap_percent = 100
monthly_native_representativity_percent = 0
monthly_representativity_percent = 0
monthly_native_max_gap_percent = 100
monthly_max_gap_percent = 100
all_representativity_percent = 0
all_max_gap_percent = 100
```

### 3hourly, 6hourly, 3hourly_instantaneous, 6hourly_instantaneous, daily

```
daily_native_representativity_percent = 0
daily_native_max_gap_percent = 100
monthly_native_representativity_percent = 0
monthly_representativity_percent = 0
monthly_native_max_gap_percent = 100 
monthly_max_gap_percent = 100
all_representativity_percent = 0
all_max_gap_percent = 100
```

### monthly

```
monthly_native_representativity_percent = 0
monthly_native_max_gap_percent = 100
all_representativity_percent = 0
all_max_gap_percent = 100
```