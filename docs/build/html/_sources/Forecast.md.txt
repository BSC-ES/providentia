# Forecast

Providentia is designed to handle multiple types of model data, including free-running, diagnostic, analysis, reanalysis, and forecast output.

For the forecast case, there are a few important details that you need to know in order to fully exploit the data, which this page will go into. 

## Interpolation

If the model data is stored as a forecast, that is, multiple forecast days are provided for each calendar day, then to ensure that all forecast days are interpolated, the `forecast` variable in the configuration file must be set to a valid value (defined in the next section).

If the `forecast` variable is set to a valid value, no matter what it is, then the interpolation will be done for all forecast days. If the `forecast` variable is not set or set to a value that is not valid, then the interpolation will happen for only the first forecast day (i.e. treating the data as a non-forecast type).

## Forecast options

```
forecast = daily
```

The forecast variable can be set to the following options: `day`, `daily`, `combined`. Setting the variable like this will utilise all of the available forecast days. It is also possible to combine a number with the forecast option to limit data to a specific forecast day, e.g. `day1`, `daily2`, `combined3`, and multiple options can set at once, e.g. `day1, day2, day3`. Different options cannot be mixed however, i.e. `day` options cannot be set with `daily` or `combined` options.

The effects of each forecast option are now explained:  

* `day` is used when wanting to compare different forecast days. The effect this has in Providentia is essentially to treat each forecast day as an independent model. Thus statistics for each forecast day can be directly compared.  

* `combined` is used when wanting to integrate data across all forecast days, providing statistics which represent the average state across all forecast days. The timeseries plot shows the average timeseries across all forecast days, for the entire time range.

* `daily` is also used when wanting to integrate data across all forecast days, however the timeseries plot is transformed to display the averaged data per forecast day, i.e. the plotted timeseries data for forecast day 1 in this case would be the average day 1 forecast across the entire time range.  

If no forecast option is set for a model with forecast data, data from the first forecast day will simply be taken to represent the model (i.e. treating the data as a non-forecast type).

It is possible to compare different models where one has forecast data and the other does not, in these cases the non-forecast data is treated as day 1 of a forecast for comparison purposes. Be aware however, if you are loading multiple species, and a model for one species has forecast data, but the same model does not for another, then this could cause Providentia to crash as it is not designed to handle this.

When using the dashboard, models that have have interpolated forecast data are automatically detected. When clicking on the **MODELS** button on the menu bar, if the model has forecast data, then when you select the model, an extra drop-down menu will appear to the right. This will allow you to select the forecast option, i.e. `combined`, `daily` or `day`. If no option is set, data from the first forecast day will be taken to represent the model (i.e. treating the data as a non-forecast type). If a forecast option is selected, another drop-down menu will then appear to the right of that, allowing you to select specific forecast days. If no specific forecast days are selected but a forecast option is, then all forecast days will be utilised. 

 


