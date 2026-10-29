# Gridded data

Besides being interpolated to the station locations, a model can be loaded as gridded data: its original output, read on the model grid.

Gridded data can be used in the dashboard, report and library modes. It is only drawn on the map: it is not used in any other plot type or in the statistics.

## Loading gridded data

### Dashboard

In the **MODELS** menu, each model has an interpolated and a gridded column. A model can be checked in both columns at the same time, to compare its interpolated values at the stations with its grid. Only one gridded model can be loaded at a time: if a second one is checked, it is unchecked again and a warning names the gridded model already checked.

If only gridded data is loaded (no observations or interpolated models), the map statistic menus and the station selection options are disabled, as there are no stations to show. Gridded models are never offered in the map **Dataset 1** and **Dataset 2** menus, as they have no data at the stations: the lead days drawn on the grid are chosen with **Grid lead days** in the map settings menu (see [Forecast data](#forecast-data)).

### Configuration files (dashboard, report and library)

Add `::gridded` at the end of the model name:

```ini
model = cams61_emep_ph2-eu-000::gridded
```

To load the same model both interpolated and gridded, list it twice. Models without a tag are read as interpolated, which can also be set explicitly with `::interpolated`:

```ini
model = cams61_emep_ph2-eu-000, cams61_emep_ph2-eu-000::gridded
```

Only one gridded model can be loaded at a time. If more are set, the first one is used. 

### Where gridded files are read from

Gridded files are searched for in the paths defined in `settings/interp_models.yaml`, and in `mod_to_interp_root`, defined in `settings/data_paths.yaml`.

Gridded data is only read when clicking **READ** or when loading a configuration file. If some days of the selected period have no valid gridded data (e.g. a forecast run or a file is missing), a warning is shown, as they are not drawn on the grid.

## Forecast data

Gridded data is read for the same forecast days selected for its interpolated version (e.g. days 2 and 4 for `day2, day4`, or the selected days for `daily` and `combined`), never for the ones of other models. If no forecast option is set, only the first forecast day is read, as for interpolated data. When forecast options are set in the model name, a gridded model loaded with its interpolated version must have the same one as it (e.g. `a59g-regional-000-day3, a59g-regional-000-day3::gridded`). The forecast option and lead days used are shown next to the gridded model in the terminal when reading, e.g. `a59g-regional-000::gridded (day option with data from lead days 2 and 4)`, or `a59g-regional-000::gridded (with data from lead day 1)` when no forecast option is applied to it.

Lead day N of a forecast run covers the 24 hours that start N-1 days after the run starts, e.g. lead day 1 of the run of `2023010100` covers 2023-01-01 00:00 to 2023-01-01 23:00. Forecast runs are expected to be 24 hours apart: a missing run leaves days without data for its lead days (see above).

See [Forecast](Forecast.md) for the forecast options.

The table below uses the same model `fc`, loaded both interpolated and gridded, with forecast data for lead days 1 to N. For each forecast mode, it shows what happens when no lead days are set (all are loaded), when only lead day 1 is set, and when lead days 1 and 3 are set:

| Forecast mode | Selected lead days | Data selected for the map              | Lead days drawn on the grid | Legend                                                   |
|---------------|--------------------|----------------------------------------|-----------------------------|----------------------------------------------------------|
| None          | None               | `fc`                                   | Day 1                       | `fc`: square and dot                                     |
| `day`         | None               | `fc-dayX` selected in `z1` or `labela` | Day X only                  | `fc-dayX`: square and dot<br>Other `fc-day`: dot only    |
| `day`         | 1                  | `fc-day1`                              | Day 1                       | `fc-day1`: square and dot                                |
| `day`         | 1 and 3            | `fc-day1` or `fc-day3` selected in `z1` or `labela` | Day of the selected entry (1 or 3) | Selected entry: square and dot<br>Other `fc-day`: dot only    |
| `daily`       | None               | `fc-daily`                             | Average of days 1 to N      | `fc-daily`: square and dot                               |
| `daily`       | 1                  | `fc-daily`                             | Day 1                       | `fc-daily`: square and dot                               |
| `daily`       | 1 and 3            | `fc-daily`                             | Average of days 1 and 3     | `fc-daily`: square and dot                               |
| `combined`    | None               | `fc-combined`                          | Average of days 1 to N      | `fc-combined`: square and dot                            |
| `combined`    | 1                  | `fc-combined`                          | Day 1                       | `fc-combined`: square and dot                            |
| `combined`    | 1 and 3            | `fc-combined`                          | Average of days 1 and 3     | `fc-combined`: square and dot                            |

The table describes report and library modes, where the data selected for the map (`z1` in reports, `labela` in library) chooses the lead day drawn. In the dashboard, the grid shows the lead days checked in **Grid lead days** instead, whatever is selected for the map (see "Dashboard: grid lead days" below).

`X` is any number from 1 to N. In the legend, a filled dot marks the interpolated data and an empty square the gridded data. Gridded models are only shown in the legend next to a map drawn over their grid: in the dashboard, and in library mode on maps with the `domain` option (the only maps with a legend). In other library plots, in the library `legend` plot, and in reports, where the legend is only shown on plots of station data, gridded models are left out of the legend. Each map also describes its grid in its title (see [Map title](#map-title)).

With `day`, each lead day is a separate entry, so the grid always shows a single lead day: the one of the selected entry. With `daily` and `combined`, all lead days are merged into one entry, so the grid shows the average of the loaded lead days, which is that same day when only one is loaded.

When the observations or another model are selected for the map, they have no lead day of their own, so the grid shows the average of the loaded lead days (lead day 1 if no forecast option is set). For example, with `day2, day4`, the grid under the `fc-day2` map shows lead day 2, while the grid under the observations map shows the average of lead days 2 and 4.

`daily` and `combined` only differ in the timeseries plot: on the map, both the stations and the grid show the mean over all loaded lead days.

If a model is loaded gridded only, it has its own forecast options, and its grid uses their lead days. With `day`, it gets one label per lead day, e.g. `fc-day2 (gridded)` and `fc-day4 (gridded)` for `day2, day4`, and each lead day gets its own map in reports, showing the observations over the grid of that lead day, as the gridded model has no data at the stations. With `daily` and `combined`, it keeps its label (e.g. `fc (gridded)`) and the grid shows the average of the selected lead days. As it has no data at the stations, it does not switch the plots of station data to the forecast layout (e.g. the `daily` timeseries per forecast day): only models with forecast data at the stations do. Where gridded models are shown in the legend, its entry is shown with a square only. In `day` mode, only the lead days drawn on the grid are shown.

### Dashboard: grid lead days

In the dashboard, the **Grid lead days** drop-down in the map settings menu lists the lead days read for the gridded model (`1` if it is not read with a forecast option), and the grid shows the average of the checked ones. It does not depend on the data selected in **Dataset 1**: for example, with `day2, day4`, the grid shows the average of lead days 2 and 4 under every map until one of them is unchecked. The map title and the legend follow the checked lead days. If none is checked, no grid is drawn.

After reading, all lead days are checked if the models have changed (e.g. a model added or its forecast days changed in the **MODELS** menu), otherwise the previous selection is kept.

### Period

In the dashboard, the grid shows the period selected for the map. In report and library modes, it shows all the loaded period.

## Map title

When a grid is drawn under a map, the map title names the gridded model and the lead days drawn, so it is clear what the grid shows even when it belongs to another model. The gridded model is shown by its alias if one is set in the configuration file. Maps without a grid (e.g. bias maps) keep their usual title.

* **Report**: on its own line, under the data label:
  ```
  fc-day2
  fc (gridded, lead day 2)
  SECTION (1751 stations)
  ```
* **Library**: at the end of the title, after a `+`, e.g. `observations Mean at 1751 stations + fc (gridded, mean of lead days 2 and 4)`. It is not added when a title is passed with `title`.
* **Dashboard**: on its own line, under the number of selected stations, e.g. `fc (gridded, mean of lead days 2 and 4)`.

## Stations over the grid

When a grid is drawn under the stations, the station markers get an edge to tell them apart from the grid. On maps without a grid, the stations are drawn without it. The edge is set with the `stations_over_grid` variable, under `plot` in the `map` section of the plot characteristics file, for each mode:

```yaml
"plot":
  {
    "stations": { "s": 20, "picker": 1, "zorder": 5, "linewidth": 0 },
    "stations_over_grid": { "edgecolor": "black", "linewidth": 0.4 },
    "grid": { "zorder": 2, "alpha": null },
  },
```

`stations_over_grid` is applied on top of `stations` when a grid is drawn, so any marker property can be changed over the grid. If it is not set (e.g. in a custom plot characteristics file made before it was added), the stations look the same with and without a grid.

## Limitations

* Only the `Mean` statistic can be drawn on the grid. Other statistics show a message and no grid is drawn.
* The grid is not drawn on bias maps (when two datasets are selected for the map), as it is an absolute field.
* Resampling and temporal colocation are not applied to gridded data. A warning is shown in the dashboard if they are active.
* Filters only apply to stations, so they do not change the grid.
