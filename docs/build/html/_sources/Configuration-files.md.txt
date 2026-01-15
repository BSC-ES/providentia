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
    period = keep: Spring, Daytime

    [[Madrid]]
    latitude = 39.57, 42.2
    longitude = -4.57, -2.42
    area_classification = keep: rural || remove: urban-suburban
```
## Sections

It is **mandatory** to define sections to launch Providentia. Their names must be wrapped in brackets (`[ ]`) and cannot include interpuncts (`·`).

A section can be though of containing all the general information needed for an analysis, e.g. `species`, `network`, `start_date`, `end_date` etc.

For a full list of available fields and their descriptions per section, please refer to the [Configuration Fields](Configuration-fields) page.

In the reports, one report is generated per section.In the dashboard and when using Providentia as a library, you will only be able to load one section or subsection at a time.

## Subsections

After defining the sections, the user can optionally create subsections and set specific information for each one. In order to do this, it is important to wrap the subsection names in double brackets (`[[ ]]`).

These subsections **MUST** be located under the sections, and can be thought of as subsets that we want to compare in the analsysis, e.g. different regions, different QA methods etc. 

Above we define two subsections as examples of different areas that we want to compare, Barcelona and Madrid. For these separate subsections, it is necessary to define the latitude and longitude which describe the area in order to delimit the analysis to the desired area. In the example, we used the following parameters for Madrid:
```
[[Madrid]]
latitude = 39.57, 42.2
longitude = -4.57, -2.42
```

We can also get the data for all the stations in a country, in this case Spain, by writing:
```
[[Spain]]
country = keep: Spain
```

Additionally, we can apply metadata filters. In the subsection for Madrid, we constrained the stations according to the `area_classification` criteria:
```
area_classification = keep: rural || remove: urban-suburban
```

If you do not know which names can be used to filter the metadata, you should use the dashboard to select the fields, export the configuration {ref}`exporting feature <exporting-configuration-files-and-data>` and copy the fields that you want from there.

`NOTE: The field names in the configuration files are case-insensitive, meaning that if you use capitals or small-cases (providing that the spelling is correct), then the field will be read properly.`

## Colocation

When performing evaluations of model data with observations, it is of high importance to ensure you are comparing apples with apples, rather than apples with oranges.

One way that evaluations can often be biased is due to gapped observations being compared with non-gapped model data. This can be resolved by ensuring both observational and model data is equally temporally gapped, called **temporal colocation**.

When loading multiple species, the number of available stations per species will most likely be different, therefore unless this is controlled for, this will lead to biases when comparing statistics across species. This can be resolved by ensuring only stations which are available for all species are retained, called **spatial colocation**.

Please see the [Colocation](Colocation) page for more information on how to apply these.

## Filtering

In both subsections and subsections, there are a vast number of options that can be applied to filter data from the configuration file, based on metadata, lower/upper bounds, representativity, QA etc. Please see the [Filtering](Filtering) page for the full guide on how to do so.

