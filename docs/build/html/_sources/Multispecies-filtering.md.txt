# Multispecies filtering

Multispecies filtering refers to the ability to filter loaded species data by the values of another species.

For example, when performing investigations of dust in the atmosphere it is a common practice to filter the AOD by the Angstrom exponent, to isolate values associated with dust. For instance, we can set the AOD to be nan when the Angstrom exponent is above 0.6, removing all non-dust AOD.

If we take a look at the timeseries for one station, we can see what this actually means:

![multispecies-filtering](uploads/multispecies-filtering.jpg)

It can be observed that there are less data points for AOD after the filtering is applied, with the removal of the data points happening when the Angstrom exponent is above 0.6.

Multispecies filtering can be set in the configuration files using the `filter_species` variable, where we define the network and species to filter by, data bounds to filter by, and an optional fill value. The first value set is the lower bound to filter by, and the second value the upper bound. Place a sign before each bound value to inform if the filter should be inclusive or exclusive of the bound, (e.g. `>` or `>=`). If no sign is set then it is assumed the bound is inclusive, i.e. `>=`. If not wishing to set either the lower or upper bounds, a `:` can be used. Optionally, a fill value can also be given as a third value to impose what the filtered data is set to, by default this is `NaN`. Spatial colocation must be active in order to apply multispecies filtering (it is active by default).

```
[All]
network = AERONET_v3_lev1.5
species = od550aero
filter_species = AERONET_v3_lev1.5:ae440-870aero (>0.6, :, nan)
spatial_colocation = True
```

On the dashboard this can be replicated under the **SPECIES** button on the menu bar. The apply (A) checkbox must be checked to apply the filter, and then data re-read by hitting the **READ** button on the main menu bar.

![multi](uploads/multi.png)

It is also possible to filter data by more than one species at a time, for example:

```
network = nasa-aeronet/directsun_v3-lev15
species = od550aero
filter_species = nasa-aeronet/directsun_v3-lev15:ae440-870aero (>0.75, <=1.2, nan), nasa-aeronet/directsun_v3-lev15:ae440-870aero (>1.2, :, 0)
spatial_colocation = True
```