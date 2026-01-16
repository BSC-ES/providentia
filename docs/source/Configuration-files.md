# Configuration files

Configuration files form the backbone of operation of Providentia. All modes use these files for configuring Providentia exactly how you want it to run. We go through here exactly how they work.

## Overview

A basic configuration file looks like the following:

```
[PRV_sconco3_a365]
network = EBAS
species = sconco3
resolution = hourly
start_date = 20180101
end_date = 20180601
model = cams61_chimere_ph2-eu-000, cams61_monarch_ph2-eu-000
temporal_colocation = True
spatial_colocation = True
report_type = standard
report_summary = True
report_stations = False
report_title = Report
report_filename = PROVIDENTIA_Report
```

This configuration can also have subsections, as in:

```
[PRV_sconco3_a365]
network = EBAS
species = sconco3
resolution = hourly
start_date = 20180101
end_date = 20180601
model = cams61_chimere_ph2-eu-000, cams61_monarch_ph2-eu-000
temporal_colocation = True
spatial_colocation = True
report_type = standard
report_summary = True
report_stations = False
report_filename = PROVIDENTIA_Report
report_title = Report

    [[Barcelona]]
    latitude = 39.8, 41.8
    longitude = 1.5, 2.5

    [[Madrid]]
    latitude = 39.57, 42.2
    longitude = -4.57, -2.42
```
## Sections

It is **mandatory** to define sections to launch Providentia. Their names must be wrapped in brackets (`[ ]`) and cannot include interpuncts (`·`).

A section can be thought of containing all the general information needed for an analysis, e.g. `species`, `network`, `start_date`, `end_date` etc.

For a list of available fields and their descriptions that can be set per section, please refer to the [configuration fields](Configuration-fields) page. The field names in the configuration files are case-insensitive, meaning that if you use capitals or lower-case (providing that the spelling is correct), then the field will be read properly.

On top of this, numerous [filter fields](Filtering) can be also be set to filter data in a variety ways.

## Subsections

After defining the sections, the user can optionally create subsections and set specific information for each one. In order to do this, it is important to wrap the subsection names in double brackets (`[[ ]]`).

These subsections **MUST** be located under the sections, and can be thought of as subsets that we want to compare in the analsysis, e.g. different regions, different QA methods etc. 

Above we defined two subsections as examples of different regions that we want to compare, **Barcelona** and **Madrid**. We select data just in those regions by keeping data within a certain longitude and latitude range, as follows:

```
[[Barcelona]]
latitude = 39.8, 41.8
longitude = 1.5, 2.5

[[Madrid]]
latitude = 39.57, 42.2
longitude = -4.57, -2.42
```

There are a wide number of [filter fields](Filtering) that can be used to filter data. If any variables are repeated between sections and subsections, the subsection variable will overwrite the information set in the section.  

The specific functionality of the sections and subsections differs between the different modes. For the reports, a report will be created per section, and then individual plots for each plot type will be made for each subsection, in each report. For the dashboard and when using Providentia as a library, you will only be able to load one section or subsection at a time.

## Colocation

When performing evaluations of model data with observations, it is of high importance to ensure you are comparing apples with apples, rather than apples with oranges.

One way that evaluations can often be biased is due to gapped observations being compared with non-gapped model data. This can be resolved by ensuring both observational and model data is equally temporally gapped, called **temporal colocation**.

When loading multiple species, the number of available stations per species will most likely be different, therefore unless this is controlled for, this will lead to biases when comparing statistics across species. This can be resolved by ensuring only stations which are available for all species are retained, called **spatial colocation**.

Please see the [colocation](Colocation) page for more information on how to apply these.

## Filtering

In both subsections and subsections, there are a vast number of options that can be applied to filter data from the configuration file, based on metadata, lower/upper bounds, representativity, QA etc. Please see the [filtering](Filtering) page for the full guide on how to do so.

