# Plot types and options

This page is only useful to those who want to create their own reports with the offline version. Users can choose which plots these reports will have, as you will see below.

## Plot types

The standard plot types consist of: maps (`map`), metadata summary (`metadata`), timeseries (`timeseries`), periodic plots (`periodic`), periodic violin plots (`periodic-violin`), box plots (`boxplot`), distribution plots (`distribution`), scatter plots (`scatter`), heat maps (`heatmap`), tables that give one stat per subsection per experiment (`table`) and tables that give multiple stats per experiment (`statsummary`).

Some of these plots involve calculating a specific statistic, namely: `map`, `periodic`, `heatmap` and `table`. This statistic is defined by aggregating the `-[stat]` field to the plot type will make a plot of the required type for the specific type of statistic set. `[stat]` should be replaced with any of the base statistic names (e.g. p5, Mean), or experiment bias names. For example to show the median values spatially, `map-p50` would be set as the plot name, or `map-r2` to show the coefficient of determination. The available statistic names are documented in `providentia/settings/basic_stats_dict.json` and `providentia/settings/experiment_bias_stats_dict.json` for reference.

For the `metadata` plot the metadata displayed is set to a default list of metadata fields. For the `statsummary` plot the statistics displayed are set to a default list of absolute and bias statistics. These default options can be changed in either the
`providentia/settings/plot_characteristics_dashboard.json` and `providentia/settings/plot_characteristics_offline.json` files, depending on which mode Providentia is being ran in.

### Map (`map`)

![Map](uploads/c09bd5457c62747147ceed0eae681c66/Map.png)

### Metadata (`metadata`)

![Metadata](uploads/4b70f94292b6bb60c03942b00dbe4aa0/Metadata.png)

### Timeseries (`timeseries`)

![Timeseries](uploads/616c4d754acb3ea71db03d7a2c2c03fe/Timeseries.png)

### Periodic (`periodic`)

![Periodic](uploads/14e16fa930ff86811775d2c6ca3b3810/Periodic.png)

### Periodic violin (`periodic-violin`)

![Periodic_violin](uploads/bcd66dd9858dc03c3bde612cbb5e3ba6/Periodic_violin.png)

### Boxplot (`boxplot`)

![Boxplot](uploads/d84ef502fe31b9fa3a0770e058b9158e/Boxplot.png)

### Distribution (`distribution`)

![Distribution](uploads/2ad285ec3f4e142208c4d087a4df6c06/Distribution.png)

### Scatter plot (`scatter`)

![Scatter](uploads/aa8bd0c96f7fc8e1dc8f285210c98cbd/Scatter.png)

### Heatmap (`heatmap`)

![Heatmap](uploads/a2fe2f27d53c2e5e1a32eeed9b07ec56/Heatmap.png)

### Table (`table`)

![Table](uploads/58bdc6481b7bbbdac66f53fa5c43b576/Table.png)

### Statistics summary (`statsummary`)

![Statsummary](uploads/88448c8e7dd9ef8515c03d5363df00da/Statsummary.png)

## Plot options

It is possible to create advanced plots by adding one or more of the following words to each basic plot type or choosing the options in the dashboard:

### Only show observations (`_obs`)

The extension `_obs` allows users to only show observations in their plots.

![Screenshot_from_2022-10-04_16-13-50](uploads/316de7338fe068ef9913e6474f20258f/Screenshot_from_2022-10-04_16-13-50.jpg)

### Split the plots by label (`_individual`)

The extension `_individual` allows users to disaggregate the plots and see the plots by experiments, individually. This can help to visualise the results in a clear way when multiple experiments have been selected.

![Screenshot_from_2022-10-04_16-13-53](uploads/42a70a0027d4397555345d9a73202d29/Screenshot_from_2022-10-04_16-13-53.jpg)

### Add annotations (`_annotate`)

If the configuration option `_annotate` is added, a box will be created on the plots to show several statistical data. The style and position of this box, as well as the statistics, can be defined by the user in plot characteristics under ``settings`` by changing the parameter `annotate_stats`.

![Screenshot_from_2022-10-04_16-17-04](uploads/906536e687e29f3e8a4ccf760485d884/Screenshot_from_2022-10-04_16-17-04.jpg)

### Get the bias of the data (`_bias`)

Alternatively the plots can be modified to show, rather than the absolute observational vs experiment values, the bias between these pairings. This is done by adding `_bias` to the base plot names, for example: `distribution_bias` or `periodic-Max_bias`.

![Screenshot_from_2022-10-04_16-10-41](uploads/f54cd53c29707eb80c5ec05e853d187d/Screenshot_from_2022-10-04_16-10-41.jpg)

### Add a smooth line to the timeseries (`_smooth`)

Adding the option `_smooth` to the `timeseries` plot will plot a smooth line over the timeseries.

![Screenshot_from_2022-10-04_16-17-13](uploads/a4bb25564231b55545a4878fb2da8fba/Screenshot_from_2022-10-04_16-17-13.jpg)

### Add a regression line to the scatter plot (`_regression`)

Adding the option `_regression` will plot the linear regression between observations and experiment.

![Screenshot_from_2022-10-04_16-17-08](uploads/d12382d6634b7ba63b7a2c04ca49cabe/Screenshot_from_2022-10-04_16-17-08.jpg)

### Make the scale logarithmic (`_logx` / `_logy`)

Adding the options `_logx` or `_logy` will set the desired axis to be logarithmically scaled. 

![Screenshot_from_2022-10-04_16-17-10](uploads/23248ef58007c22ab1021878cae710e5/Screenshot_from_2022-10-04_16-17-10.jpg)

### Get plot by more than one network species (`_multispecies`)

Incorporate all read species in the plot type. 

![Screenshot_from_2022-10-04_16-14-01](uploads/2e4902e7d1be1d79db47b5818c7288a2/Screenshot_from_2022-10-04_16-14-01.jpg)

### Hide points and only show regression / smooth lines (`_hidedata`)

The option `_hidedata` needs to be accompanied by `_smooth` in the `timeseries` plot and by `_regression` in the `scatter` plot.

![hidedata](uploads/341d06f77cebecffc9802d0190c0de91/hidedata.png)

### Show the model grid in the maps (`_domain`)

Adding `_domain` will add the model grid on top of the map.

![domain](uploads/4dd8180f18b54e2a3c6d4e0e39827d0c/domain.png)