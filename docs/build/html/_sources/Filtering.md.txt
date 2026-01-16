# Filtering

In Providentia there are substantial number of options when it comes to filtering data. Here we go through each of these one by one, and explain how they can be applied through the dashboard, or via configuration file.

## Dashboard filtering

ALL of the filters explained in the following sections can be selected on the dashboard menu bar, either in the **Data Selection** section, or in the **Filters** section. 

A general important rule to remember is when things are modified in the **Data Selection** section, then data must be re-read using the **READ** button for the filters to be applied.

When things are modified in the **Filters** section, then the **FILTER** button must be clicked to apply the filters. These filters can then be removed (back to defaults) by clicking the **RESET** button. 

(filtering_metadata)= 
## Metadata

The most common type of data filtering is by metadata, e.g. for country, longitude, altitude etc. 

Metadata can either by numerical or text and the method for filtering for these slightly varies, as will be described. There will be substantially more metadata fields available when using GHOST data compared to when not.  

One common question is what metadata available variables exist for the loaded data that can be filtered by. These can be seen empirically by looking through the sub-menus under the **META** button on the dashboard, with all available variables are organised into 5 categories: `station position`, `station classifications`, `station miscellaneous`, `globally gridded classifications`, and `measurement process information`. These are also defined for GHOST data in this [publication](https://doi.org/10.5194/essd-16-4417-2024).

The subsequent question is what fields are then available for each metadata variable, specifically for text metadata. The best solution for this is to empirically check by opening the dashboard and navigating to the variable you are interested in and viewing the available fields.

### Numerical metadata

Numerical metadata is filtered by setting a lower and upper bound of values to retain, for example only keeping stations with latitudes greater or equal to 30°N, but less or equal to 72°N: 

```
latitude = 30, 72
```

By default on the dashboard metadata menu, the variables will display the full available range for the variable, e.g. minimum and maximum latitude.

### Text metadata

Text metadata is filtered by either setting to **keep** or **remove** specific fields for a variable, for example keeping Spanish and French stations, or removing just UK stations: 

```
country = keep: Spain, France
country = remove: United Kingdom
```

On the dashboard metadata menu all available fields for a variable will be displayed, with associated keep (K) and remove (R) checkboxes per field. Simply check the fields you wish you keep or remove for that variable.

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

Leaving the variables blank means they will be not applied at all. 

It can also be set to add or subtract from the default fields, for example:
```
add_qa = 3, 10
subtract_qa = 3, 10
add_flags = 3, 10
subtract_flags = 3, 10
``` 

When using the dashboard by clicking on the **FLAGS** or **QA** buttons on the menu bar, pop-up menus will appear giving the option to interactively select fields.

See [here](QA-filtering) for more in-depth information on filtering by QA or data flags, as well as full definitions of all available flags that can be filtered by. 

## Multispecies filtering

Multispecies filtering refers to the ability to filter loaded species data by the values of another species. For example, when performing investigations of dust in the atmosphere it is a common practice to filter the AOD by the Angstrom exponent, to isolate values associated with dust.

Multispecies filtering can be set in the configuration files using the `filter_species` variable, where we define the network and species to filter by, the lower range to filter associated data by, the upper range to filter associated data by, and the value to set associated data to (default is nan). If you want to leave one of the filter ranges open-ended, then use a colon (:). Spatial colocation must be active in order to apply multispecies filtering (it is active by default).

```
[All]
network = AERONET_v3_lev1.5
species = od550aero
filter_species = AERONET_v3_lev1.5:ae440-870aero (>0.6, :, nan)
spatial_colocation = True
```

On the dashboard this can be replicated under the **MULTI** button on the menu bar. The apply (A) checkbox must be checked to apply the filter, and then data re-read by hitting the **READ** button on the main menu bar.

See [here](Multispecies-filtering) for more in-depth information on multspecies filtering.

(representativity)= 
## Representativity

One major limitation often associated with observations is the amount of gaps between measurements. If these observations are directly compared with typically complete model data, this would impose a significant bias upon the comparison. Filtering by representativity filters provides a way to control the temporal robustness of the observations for evaluation. 

Providentia has multiple such filters available.

### Native vs averaged filters 

Providentia has two types of representativity filters, those which calculate the representativity percent of observations at the **native** measured resolution (which can be variable), and those which calculate it at an **averaged** resolution, i.e. the data resolution being used in Providentia such as hourly, 3hourly etc. 

The **native** filters are only available when using GHOST data, and will always return a lower value than the **averaged** ones, as an hourly period with just 30 minutes being represented with measurements would be classed as being 100% represented in the average sense, whereas it would be 50% in the native sense.

### Available filters

For both **native** and **averaged** filters, it is possible to filter by the representativity across different periods. For example, data for each day with representativity less than a certain amount can be screened out (i.e. set as nan).

The periods that representativity can be filtered by are: **hour**, **day**, **month**, and **year**. The hourly case only applies for **native** data (when the data resolution is hourly).

There is an additional filter for the **averaged** case, where the representativity of data across across the entire time range can be filtered by, thus removing entire stations when the representativity is less than a certain percent. This is set by the `all_representativity_percent` variable in the configuration file, and the `All` line under the **% REP** button on the dashboard.

The **native** and **averaged** representativity filters can be set in the configuration file as follows, with the difference in the naming convention being that the native fields contain the `native` string. By default all fields are set to be 0% (i.e. the filters are not active).

```
# Native representativity filters
hourly_native_representativity_percent = 0
daily_native_representativity_percent = 0
monthly_native_representativity_percent = 0
annual_native_representativity_percent = 0

# Averaged representativity filters
daily_representativity_percent = 0
monthly_representativity_percent = 0
annual_representativity_percent = 0
all_representativity_percent = 0
```

As you change the temporal resolution of the data to be coarser, some of these fields will become unavailable. For example when using monthly data, you can not filter by daily representativity.

On the dashboard the filters can be set under the **% REP** button on the menu bar.

(periods)= 
## Periods

It is often desired to select or remove certain periods over a time range, for example just keep daytime data, or remove summertime data. Providentia gives an easy way to filter in such a way using the `period` variable when using GHOST data. 

The available period fields that that can be selected are: `Daytime`, `Nighttime`, `Weekday`, `Weekend`, `Spring`, `Summer`, `Autumn`, `Winter`. The availability of some of these fields is contingent on the active temporal resolution of data, for example **Daytime** values cannot be selected when monthly data is loaded.

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

Additionally in the case of multiple species, then multiple bounds can be given for the number of species that are loaded, for example:

```
species = sconco3, pm2p5, pm10
lower_bound = 10, 20, 30
upper_bound = 1000, 2000, 3000
```
On the dashboard these bounds can be set through the lower and upper bounds rangeboxes on the menu bar.

## Calibration factor

While not explicitely filtering per se, a commonly used feature to correct for model biases is by way of a calibration factor. This is used when a model has a known bias and you want to correct for it.

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