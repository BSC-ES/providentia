# Plot types and options

Users can choose the plots from the options below.

## Plot types

The standard plot types consist of:

- Map (`map`)
- Metadata summary (`metadata`)
- Timeseries (`timeseries`)
- Periodic plot (`periodic`)
- Periodic violin plot (`periodic-violin`)
- Box plot (`boxplot`)
- Distribution plot (`distribution`)
- Scatter plot (`scatter`)
- Heat map (`heatmap`) - Not available in dashboard
- Table that gives one statistic per subsection per model (`table`) - Not available in dashboard
- Table that gives multiple statistics per model (`statsummary`)
- Taylor Diagram (`taylor`)
- FAIRMODE target plot (`fairmode-target`)
- FAIRMODE statistics summary plot (`fairmode-statsummary`)

Some of these plots are created for specific statistics, namely: `map`, `periodic`, `heatmap`, `taylor` and `table`. In the report and library modes, this statistic is defined by aggregating the `-[stat]` field to the plot type. `[stat]` should be replaced with any of the base statistic names (e.g. p5, Mean), or model bias names (e.g. r2, RMSE). For example to show the median values spatially, `map-p50` would be set as the plot name, or `map-r2` to show the coefficient of determination. The available statistic names are documented in `settings/basic_stats.yaml` and `settings/model_bias_stats.yaml`. For the Taylor diagram, only `r`and `r2`can be used.

The timeseries can also be used to show how statistics vary in time. In order to do this, we need to add `-[stat]` and the temporal resolution after the plot type name (e.g. `timeseries-Mean-daily`, `timeseries-r2-monthly`, `timeseries-r-annual`).

For the `metadata` plot the metadata displayed is set to a default list of metadata fields. For the `statsummary` plot the statistics displayed are set to a default list of absolute and bias statistics. These default options can be changed in `settings/plot_characteristics.yaml`.

### Map (`map`)

![Map](uploads/Map.png)

### Metadata (`metadata`)

![Metadata](uploads/Metadata.png)

### Timeseries (`timeseries`)

![Timeseries](uploads/Timeseries.png)

### Periodic (`periodic`)

![Periodic](uploads/Periodic.png)

### Periodic violin (`periodic-violin`)

![Periodic-violin](uploads/Periodic-violin.png)

### Boxplot (`boxplot`)

![Boxplot](uploads/Boxplot.png)

### Distribution (`distribution`)

![Distribution](uploads/Distribution.png)

### Scatter plot (`scatter`)

![Scatter](uploads/Scatter.png)

### Heatmap (`heatmap`)

![Heatmap](uploads/Heatmap.png)

### Table (`table`)

![Table](uploads/Table.png)

### Statistics summary (`statsummary`)

![Statsummary](uploads/Statsummary.png)

### Taylor diagram (`taylor`)

![taylor](uploads/taylor.png)

### FAIRMODE target (`fairmode-target`)

![fairmode-target](uploads/fairmode-target.png)

### FAIRMODE statistics summary (`fairmode-statsummary`)

![fairmode-statsummary](uploads/fairmode-statsummary.png)

## Plot options

It is possible to create advanced plots by adding one or more of the following words to each basic plot type or choosing the options in the dashboard:

### Only show observations (`_obs`)

The extension `_obs` allows users to only show observations in their plots.

![obs](uploads/obs.jpg)

Plot types: 
- Dashboard: Not applicable 
- Report: `map`, `timeseries`, `periodic`, `periodic-violin`, `distribution`, `boxplot`
- Library: Not applicable 

### Split the plots by label (`_individual`)

The extension `_individual` allows users to disaggregate the plots and see the plots by models, individually. This can help to visualise the results in a clear way when multiple models have been selected.

![individual](uploads/individual.jpg)

Plot types: 
- Dashboard: Not applicable
- Report: `timeseries`, `periodic`, `periodic-violin`, `distribution`, `scatter`, `boxplot`, `taylor`, `fairmode-target`, `fairmode-statsummary`
- Library: Not applicable

### Add annotations (`_annotate`)

If the configuration option `_annotate` is added, a box will be created on the plots to show several statistical data. The style and position of this box, as well as the statistics, can be defined by the user in `plot_characteristics.yaml` under `settings` by changing the parameter `annotate_stats` per plot type.

![annotate](uploads/annotate.jpg)

Plot types: 
- Dashboard: `map`, `timeseries`, `periodic`, `periodic-violin`, `distribution`, `scatter`, `boxplot`, `taylor`, `fairmode-target`
- Report: `map`, `timeseries`, `periodic`, `distribution`, `scatter`, `boxplot`, `heatmap`, `taylor`, `fairmode-target`
- Library: `map`, `timeseries`, `periodic`, `distribution`, `scatter`, `boxplot`, `heatmap`, `taylor`, `fairmode-target`

### Get the bias of the data (`_bias`)

Alternatively the plots can be modified to show, rather than the absolute observational vs model values, the bias between these pairings. This is done by adding `_bias` to the base plot names, for example: `distribution_bias` or `periodic-Max_bias`.

![bias](uploads/bias.jpg)

Plot types: 
- Dashboard: `timeseries`, `periodic`, `distribution`, `statsummary`
- Report: `map`, `timeseries`, `periodic`, `distribution`, `heatmap`, `table`, `statsummary`
- Library: `timeseries`, `periodic`, `distribution`, `heatmap`, `table`, `statsummary`

### Add a smooth line to the timeseries (`_smooth`)

Adding the option `_smooth` to the `timeseries` plot will plot a smooth line over the timeseries.

![smooth](uploads/smooth.jpg)

Plot types: 
- Dashboard: `timeseries`
- Report: `timeseries`
- Library: `timeseries`

### Add a regression line to the scatter plot (`_regression`)

Adding the option `_regression` will plot the linear regression between observations and model.

![regression](uploads/regression.jpg)

Plot types: 
- Dashboard: `scatter`
- Report: `scatter`
- Library: `scatter`

### Hide points and only show regression / smooth lines (`_hidedata`)

The option `_hidedata` needs to be accompanied by `_smooth` in the `timeseries` plot and by `_regression` in the `scatter` plot.

![hidedata](uploads/hidedata.png)

Plot types: 
- Dashboard: `timeseries`, `scatter`
- Report: `timeseries`, `scatter`
- Library: `timeseries`, `scatter`

### Make the scale logarithmic in x axis (`_logx`)

Adding the options `_logx` will set the x axis to be logarithmically scaled. 

![logs](uploads/logs.jpg)

Plot types: 
- Dashboard: `distribution`, `scatter`
- Report: `distribution`, `scatter`
- Library:  `distribution`, `scatter`

### Make the scale logarithmic in y axis (`_logy`)

Adding the options `_logy` will set the y axis to be logarithmically scaled.

Plot types: 
- Dashboard: `timeseries`, `periodic`, `periodic-violin`, `distribution`, `scatter`, `boxplot`
- Report: `timeseries`, `periodic`, `periodic-violin`, `distribution`, `scatter`, `boxplot`
- Library: `timeseries`, `periodic`, `periodic-violin`, `distribution`, `scatter`, `boxplot`

### Get plot by more than one network species (`_multispecies`)

Incorporate all read species in the plot type. 

![multispecies](uploads/multispecies.jpg)

Plot types: 
- Dashboard: Not applicable
- Report: `boxplot`, `heatmap`, `table`, `statsummary`
- Library: `boxplot`, `heatmap`, `table`, `statsummary`

### Show the model grid in the maps (`_domain`)

Adding `_domain` will add the model grid on top of the map.

![domain](uploads/domain.png)

Plot types: 
- Dashboard: `map`
- Report: `map`
- Library: `map`

### Add threshold line (`_threshold`)

Adding `_threshold` will add a line indicating the exceedances. These exceedances are set in `settings/exceedances.yaml`.

![threshold](uploads/threshold.png)

Plot types: 
- Dashboard: `timeseries`, `periodic`, `periodic-violin`, `distribution`, `scatter`, `boxplot`
- Report: `timeseries`, `periodic`, `periodic-violin`, `distribution`, `scatter`, `boxplot`
- Library: `timeseries`, `periodic`, `periodic-violin`, `distribution`, `scatter`, `boxplot`

### Normalise boxplot (`_normalise`)

Adding `_normalise` will normalise the boxplot. 

![normalise](uploads/normalise.png)

Plot types: 
- Dashboard: Not applicable
- Report: `boxplot`
- Library: `boxplot`