# Overview

The library was designed for users to exploit the powerful backend of Providentia, allowing the use of the Providentia functionalities directly in your own scripts or within [Jupyter notebooks](Notebooks).

## Getting started

If the script you are running is not inside the Providentia home directory, in order to import from Providentia, then it is necessary to add these lines to your code at the start:      

```
import sys
sys.path.append(path)
```

where `path` is the path where your Providentia code exists.     

You will also need to load the necessary modules for it to function, for this you can activate the local environment that is created the first time you run the tool:

```
conda activate providentia-env_v[version]
```

where `version` is the relevant version of Providentia, e.g. `3.0.2`.

Then the Providentia library can be used in any script, through simply importing it as a module:

```
import providentia as prv
```

After the import you will have full access to the backend of Providentia and all of its functions as described in the following sections.

NOTE: In the future we will allow for Providentia to be imported as a module on the BSC machines, making this step redundant. 

## Library configuration fields

All parameters that can be used in the library configuration files can be found in the [Shared Parameters](shared-parameters) or [Analysis and Visualization Parameters](visualization-parameters) sections of the Configuration Fields page.

## Features

Importing Providentia gives access to the entirety of the Providentia backend. This includes use of all of the different Providentia modes directly in your own scripts, such as interpolation, report etc, as well as all the smaller functions that Providentia performs, such as reading, filtering, plotting, writing, etc. through the library mode. 

Tutorials on how to use the library mode can be found under the folder [tutorials](https://github.com/BSC-ES/providentia/tree/master/tutorials).

### Accessing the Providentia modes

Each of the different Providentia modes can be used through the Providentia object, as well as would normally execute them via the command line. 

After importing Providentia, you need to create the object using a configuration file:
```
provi = prv.Providentia('interactive_template.conf')
```

Providentia looks for configuration files in the `configurations` folder by default, if your configuration file is saved elsewhere, pass the complete path.

The configuration will become part of the object, therefere you can call the modes as explained below. You can pass any arguments you wish when calling each function to override what is set in the loaded configuration file. 

When using Providentia inside a Python script, instead of through the terminal or in a Jupyter notebook, we need to call its functions with the main guard:

```
if __name__ == "__main__":
    provi = prv.Providentia('interactive_template.conf')
    # Other functions (load, plot, etc.)
```

**Download**

To download data:

```
provi.download()
```

**Interpolation**

To run the interpolation:

```
provi.interpolation()
```

**Dashboard**

To open the dashboard:

```
provi.dashboard()
```

**Report**

To generate a report:

```
provi.report()
```

### Library: Reading and filtering data from a configuration file

The Library mode was specifically designed to access to all Providentia sub-functions such as reading, filtering, plotting, etc.

As a first step, data needs to be read and filtered appropriately, for this we use the funcion `load`:

```
provi.load()
```

where `provi` is the object instance previously defined, which can access all of the class variables and methods.  

Only one section-subsection pair in the configuration file can be read and filtered by the library at one time, with the exception in the creation of plots with the `multispecies` option (see plotting section). Where there is more than one pair defined in the file, the specific pair wished to be read can be set by passing the section and/or subsection arguments, e.g.:

```
provi = prv.Providentia('interactive_template.conf', section='OTHERSECTION', subsection='OTHERSUBSECTION')
provi.load()
```

If there are multiple section and subsection pairs in the configuration file, and a specific pair is not explicitly set to be read, then simply the first pair is taken to be read.

When wanting to change any data or filters applied to the data, simply update the configuration file and reinitiate the class instance. 

If wanting to overwrite any arguments in the configuration file, directly in the script, each argument can be simply passed when initiating the class instance, e.g.:

```
provi = prv.Providentia('interactive_template.conf', network='EANET')
provi.load()
```

### Library: Reset filter

After data has been read and filtered, it will stay filtered unless reset. This can be done by the following method: 

```
provi.reset()
```  

If at any point wanting to return to the state when the Interactive class was initialised, this can be done by adding the following argument to the method:

```
provi.reset(initialise=True)
```

### Library: Applying filters

If wanting to apply a filter not set in the configuration file, this can be done using the following method:

```
provi.filter(field, ...)
```

where `field` is the field to filter by. The fields to filter by can be representativity fields, period fields or metadata fields.

If the field is numeric, lower and upper limits to retain data between can be set as follows:

```
provi.filter(field, lower=28, upper=31)
```

NOTE: It is not mandatory to pass `lower` and `upper` arguments together.

If the field is textual, values to keep and remove associated data with can be set as follows:

```
provi.filter(field, keep='Spain', remove='')
```

NOTE: It is not mandatory to pass `keep` and `remove` arguments together.

If multiple values are wanted to be removed concurrently, the arguments passed should be lists. 

```
provi.filter(field, remove=['Spain','France'])
```           

For the specific case of representativity fields, the argument passed should be `limit`:

```
provi.filter(rep_field, limit=20)
```           

### Library: Filtering for specific station/s

If wanting to filter data for a specific station or stations, the following convenience method can also be used (rather than using `apply_filter`):

```
provi.filter_station(station)
```

where `station` corresponds to the `station_reference` of one station of interest or can be a list of multiple stations given their `station_reference` values. This will then subsequently mean the data in memory is filtered for the relevant station or stations.

### Library: Print active configuration file

If wanting to visually print the configuration file that is currently active this can be done by:
```
provi.print_config()
```

If wanting to visually print any other existing configuration, this can be done by passing the `conf` or `config` argument: 
```
provi.print_config(conf='another-configuration.conf')
```

### Library: Calculate statistics

In order to calculate statistics for a dataset, this can done via the following method:

```
stat_calc = provi.statistic(stat, labela='OBS')
```

where `stat` is the statistic wished to be calculated, and `labela` is the name of the observations/model data (this can be an alias set in the configuration file or the original dataset name). The statistic returned will be one summary value, which can be formulated in differing ways, e.g. taking a median time series across all stations and then calculating the statistic or calculating the statistic across all stations and then taking the median. This formulation can be set by changing the statistical mode in the configuration file, for which more information can be found [here](Statistics).

If wanting to calculate a bias statistic, then an additional dataset is needed to be compared against, which can be set via the following argument:

```
labelb='EMEP'
```

where `EMEP` in this case is the name of the dataset that is wished to be compared against.

For bias statistics for which a subtraction is involved, it is always done as `datasetb - dataseta`.

If wanting to calculate statistics at each individual station then this can be set by adding the following argument:

```
per_station=True
```

### Library: Saving data

The data which has been read and filtered can be saved out using the following method: 

```
provi.save(format='nc')
```

where `format` is the format of the saved data, and can be `nc` (netCDF), `np` (numpy) or `conf` (Providentia configuration file). In the case of `conf`, rather than the data being saved, a Providentia configuration file is generated, which when loaded would return the exact same data in memory.

The filename of the saved data can be set by passing the following argument:

```
fname='/mypath/myfilename'
```

If `fname` is not provided then one is generated automatically, and saved in the directory `saved_data`.

More details can be found in the [Saved file formats section](Saved-file-formats).
 
### Library: Accessing data in memory

As well as being able to save data out, data can be returned directly in memory in specific formats, via the following method: 

```
data = provi.data(format='xr')
```

where `format` is the format of the returned data, and can be `nc` (netCDF), `np` (numpy) or `xr` (xarray).

The returned variable will contain all read data and metadata variables available. Therefore, if reading a GHOST network you can expect the metadata to be more substantial.

### Library: Accessing specific variable in memory

If wanting to extract a specific variable in memory, rather than the entire read dataset, this can be done using the following method:

```
var_data = provi.variable(var='myvar')
```

where `var` is the name of the variable wished to be read.

### Library: Plotting 

All plot types that are available in other modes of Providentia are available through the library, with the addition of the option to plot the legend as a standalone plot. All plots can be made using the method:

```
provi.plot(plot_type)
```

where `plot_type` is the plot type wished to be made (all types stated later in this section).

#### Plot customisation

All plots have been tailored to appear nice in Jupyter notebooks, but any plot settings can be modified by the user, per plot type, in `settings/plot_characteristics.yaml`.

On top of that, there are numerous arguments to the plot method that can be used to customise the plots, each of which are subsequently detailed.

##### Species

In non-multispecies plots you can choose the species that you want to plot using the `networkspeci` argument, e.g.:

```
networkspeci='EEA|sconco3'
```

##### Legend

For user convenience a legend has been integrated into each appropriate plot type. If this is not wanted, then this can be deactivated using the following argument to the plot method:

```
legend=False
```

If wanting to remove the observations data label from the legend, then this can be done by passing the following argument:

```
set_obs_legend=False
```

##### Colourbar

For some plots a colourbar is integrated into the plot, e.g. map. This can be deactivated by passing the following argument:

```
cb=False
```

##### Limiting plotted data

By default, all read data will be plotted. If wanting to limit which data is plotted, this can be done by passing a list of the relevant data labels to be plotted, e.g.:

```
data_labels=['EMEP', 'MONARCH']
```

NOTE: These labels should be the dataset aliases, if set.

##### Setting title, x or y labels

The title, xlabel and ylabel will be set automatically (if required) for all plot types, however if wanting to overwrite them this can be done by passing the following arguments:

```
title='My custom title'        
xlabel='My custom xlabel'        
ylabel='My custom ylabel'         
```

##### Custom formatting

To have more granular control of the plot formatting, and overwrite specific variables set in `settings/plot_characteristics.yaml`, the `format` argument can be used to pass a dictionary which sets specific plot type format variables to overwrite, e.g.:

```
format={"figsize": [14,7], "xtick_params": {"labelsize": 22}}
```

##### Plot options

As in reports, plot options can be set for each plot type. There are 2 ways to set the plot options. One way is passing them all through a list, e.g.:

```
plot_options=['annotate', 'bias']
```

The other way is to pass each individual plot option as an argument e.g.:

```
annotate=True, bias=True
```

The available plot options per plot type are given in the sections below.

##### Saving / returning plot object in memory

Rather than viewing the plot on the screen, it can be returned in memory for more editing or saved to a file.

To return the plot in memory, this can be done as follows:

```
plot_obj = provi.plot(plot_type, return_plot=True)
```

where `plot_obj` is the object of the plot returned.

To save the plot, the following argument can be passed:

```
save='fname'
```

where `save` is the name of the file to save to. If save is set as: `save=True`, then the filename will be automatically generated and saved in the `plots` directory.

#### Saving plot data

Use `save_data` to store the plot elements as CSV files.

```
plot_obj = provi.plot(plot_type, save_data=True)
```

By default, the files will be saved in the `saved_data` folder, users can change the path by passing `save_data_path`.

```
plot_obj = provi.plot(plot_type, save_data=True, save_data_path='path')
```

#### Plot types

In this section the available plot types will de detailed, with any specific subtleties involved in calling the plotting method detailed. Available plot options are also stated per plot type. 

##### legend

The legend plot can be made as follows:  

```
provi.plot('legend')
```

It has no available plot options.

##### metadata

The metadata plot can be made as follows:    

```
provi.plot('metadata')
```

It has no available plot options.

##### map

The map plot can be made as follows:   

```
provi.plot('map-stat', labela='OBS')
```

where stat is the statistic to be plotted per station on the map, and `labela` is the data label of the dataset to calculate the statistics with.

By adding `labelb` to the plot arguments, you can turn the plot type into a bias plot:

```
provi.plot('map-stat', labela='OBS', labelb='MONARCH')
```

For bias statistics for which a subtraction is involved, it is always done as `datasetb - dataseta`.

The map extent is by default automatically set based on what is plotted on the map. This can be also manually controlled by passing the argument:

```
map_extent = [lonmin, lonmax, latmin, latmax]
```

where `lonmin`, `lonmax`, `latmin`, and `latmax` are the rectangular bounds of the map extent, e.g.: `map_extent = [-180, 180, -90, 90]`

The available plot options are: `annotate` and `domain`.
 
##### timeseries

The timeseries plot can be made as follows:   

```
provi.plot('timeseries')
```

By adding a statistic and resolution, we can get statistical timeseries to see the variation of a statistic with time:

```
provi.plot('timeseries-stat-resolution')
```

The available plot options are: `annotate`, `bias`, `hidedata`, `logy` and `smooth`.

##### periodic

The periodic plot can be made as follows:   

```
provi.plot('periodic-stat')
```

where `stat` is the statistic that is wanted to be plotted. 

The available plot options are: `annotate`, `bias` and `logy`.

##### periodic-violin

The periodic-violin plot can be made as follows: 

```
provi.plot('periodic-violin')
```

The available plot options are: `annotate`and `logy`.

##### distribution

The distribution plot can be made as follows: 

```
provi.plot('distribution')
```

The available plot options are: `annotate`, `bias`, `logx` and `logy`.

##### scatter

The scatter plot can be made as follows: 

```
provi.plot('scatter')
```

The available plot options are: `annotate`, `hidedata`, `logx`, `logy` and `regression`.

##### boxplot

The boxplot can be made as follows: 

```
provi.plot('boxplot')
```

The available plot options are: `annotate`, `logy` and `multispecies`.

##### heatmap

The heatmap plot can be made as follows: 

```
provi.plot('heatmap-stat')
```

where `stat` is the statistic that is wanted to be plotted. 

The available plot options are: `annotate`, `bias` and `multispecies`.

##### table

The table plot can be made as follows: 

```
provi.plot('table')
```

The available plot options are: `bias` and `multispecies`.

##### statsummary

The statsummary plot can be made as follows: 

```
provi.plot('statsummary')
```

The available plot options are: `bias` and `multispecies`.

##### taylor

The Taylor diagram can be made as follows: 

```
provi.plot('taylor-stat')
```

where `stat` is the statistic that is wanted to be plotted, from `r` or `r2`. 

The available plot option is: `annotate`.

##### fairmode-target

The FAIRMODE target plot can be made as follows: 

```
provi.plot('fairmode-target')
```

The available plot option is: `annotate`.

##### fairmode-statsummary

The FAIRMODE statistic summary plot can be made as follows: 

```
provi.plot('fairmode-statsummary')
```

It has no available plot options.