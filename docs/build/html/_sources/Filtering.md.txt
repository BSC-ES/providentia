# Filtering

In Providentia there are substantial number of options when it comes to filtering data. Here we go through each of these one by one, and explain how they can be applied through the dashboard, or via configuration file.

## Dashboard filtering

ALL of the filters explained in the following sections can be selected on the dashboard menu bar, either in the **Data Selection** section, or in the **Filters** section. 

A general important rule to remember is when things are modified in the **Data Selection** section, then data must be re-read using the **READ** button for the filters to be applied.

When things are modified in the **Filters** section, then the **FILTER** button must be clicked to apply the filters. These filters can then be removed (back to defaults) by clicking the **RESET** button. 

## Metadata

The most common type of data filtering is by metadata, e.g. for country, longitude, altitude etc. 

Metadata can either by numerical or text and the method for filtering for these slightly varies, as described in the following sections. There will be substantially more metadata fields available when using GHOST data compared to not.  

One common question is what metadata available variables exist for the loaded data that can be filtered by. These can be seen empirically by looking through the sub-menus under the **META** button on the dashboard, with all available variables are organised into 5 categories: `station position`, `station classifications`, `station miscellaneous`, `globally gridded classifications`, and `measurement process information`. These are also defined for GHOST data in this [publication](https://doi.org/10.5194/essd-16-4417-2024).

The subsequent question is what fields are available for each metadata variable, specifically for text metadata. The best solution for this is to empirically check by opening the dashboard and navigating to the variable you are interested in and viewing the available fields.

### Numerical metadata

Numerical metadata is filtered by setting a lower and upper bound of values to retain, for example only keeping stations with latitudes greater or equal to 30°N, but less or equal to 72°N: 

```
latitude = 30, 72
```

By default on the dashboard, the variables will display the full available range for the variable, e.g. minimum and maximum latitude.

### Text metadata

Text metadata is filtered by either setting to **keep** or **remove** specific fields for a variable, for example keeping Spanish and French stations, or removing just UK stations: 

```
country = keep: Spain, France
country = remove: United Kingdom
```

On the dashboard all available fields for a variable will be displayed with associated keep and remove checkboxes per field. Simply check the fields you wish you keep or remove.

## QA flags

Observations of atmospheric species are taken by scientists in the real-world where issues with instruments, meteorology, or even human error can mean observations are subject to significant biases, which left alone could impose significant biases for evaluations with model data. Fortunately when this occurs observations are typically flagged so they can be screened out.

There are two sets of QA flags available in Providentia for filtering observations when using GHOST or ACTRIS data.

* **flags** - flags that relate to standardised data flags taken from the data provider.
* **qa** - flags that relate to GHOST performed quality control checks.

These variables can be set explicitely (using codes or names), for example:
```
flags = 0, 1, 2
qa = 0, 1, 2
flags = Valid Data, Preliminary Data, Missing Data
qa = Missing Measurement, Infinite Value, Negative Measurement
flags = 
qa = 
```

Leaving the variables blank means they will be not applied at all. It can also be set to add or subtract from the default fields, for example:
```
add_qa = 3, 10
subtract_qa = 3, 10
add_flags = 3, 10
subtract_flags = 3, 10
``` 

When using the dashboard by clicking on the **FLAGS** or **QA** buttons on the menu bar, pop-up menus will appear giving the option to interactively select fields.

See [here](QA-filtering) for more in-depth information on filtering by QA or data flags, as well as full definitions of all available flags that can be filtered by. 

## Representativity

One major limitation associated with observations is the amount of gaps often found between measurements. If these observations are directly compared with typically complete model data, this would impose a significant bias upon the comparison. Filtering by representativity filters provides a way to control the temporal robustness of the observations for evaluation. 

Providentia has multiple such filters. By default these are **NOT** activated (i.e. set at 0). 

```
hourly_native_representativity_percent = 20
```

On the dashboard the representativity variables can be modified via the **% REP** button on the menu bar.

See [here](Representativity-filtering) for more in-depth information on filtering by representativity.

## Periods

It is often desired to select certain periods over a time range, for example just daytime data, or data in the summer. Providentia gives an easy way to filter in such a way using the `period` variable when using GHOST data. 

The available period fields that that can be selected are: `Daytime`, `Nighttime`, `Weekday`, `Weekend`, `Spring`, `Summer`, `Autumn`, `Winter`.The availability of some of these fields is contingent on the temporal resolution that is active, for example **Daytime** data cannot be selected for monthly data.

When wanting to apply these via configuration file the syntax is the same as the text metadata filtering with a small caveat, that you can both keep and remove fields at the same time by using a **double pipe "||"** between the keep and remove definitions to distinguish between them, for example:

```
period = keep: Winter
period = remove: Daytime
period = keep: Spring, Summer || remove: Weekday
```

On the dashboard the period fields to keep or remove can be accessed via the **PERIOD** button on the menu bar.

## Bounds

Often it is desired to remove values which exceed certain extreme bounds, as it is known that data should appear at such extremes. These bounds will be by default active in Providentia, with extreme bounds associated with a given species taken from definitions in GHOST.

These bounds can also be modified in the configuration file as follows:

```
lower_bound = 10
upper_bound = 1000
```

Additionally in the case of multiple species, then multiple bounds can be given for the number of species that exist, for example:

```
species = sconco3, pm2p5, pm10
lower_bound = 10, 20, 30
upper_bound = 1000, 2000, 3000
```
On the dashboard these bounds can be set through the lower and upper bounds rangeboxes on the menu bar.

## Calibration factor

Which not explicitely filtering per se, a commonly used feature to correct for model biases is by way of a calibration factor. This is used when a model has a known bias and you want to correct for it.

It is defined simply as a factor that is applied to model data. It can be set via `calibration_factor` in the configuration file, in the following ways:

- To add:
```
calibration_factor = +10
```
- To subtract:
```
calibration_factor = -10
```
- To multiply:
```
calibration_factor = *10
```
- To divide:
```
calibration_factor = /10
```

The calibration factor can also be defined independently for diferent models:
```
calibration_factor = a54s-regional-000 (*0.62), a4xf-regional-000 (*0.51)
```

It can also be defined independently for different species:
```
network = EEA_AQ_eReporting
species = pm2p5, pm10
calibration_factor = a54s-regional-000 (*0.62, *0.4), a4xf-regional-000 (*0.52, *0.9)
```

On the dashboard there is no current way to interactively set the calibration factor, so data must be loaded by configuration file first if wanting to apply this.