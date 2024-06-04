# Colocation

There are two types of colocations that can be used in Providentia:

- Temporal colocation
- Spatial colocations

## Temporal colocation

The temporal colocation can be used to remove the gaps where the data of the observations or experiments are missing. Once the colocation is turned on, the user has access to more plot types (i.e. scatter plot and Taylor diagram).

**Without temporal colocation:**

![Screenshot_from_2024-06-03_15-16-50](uploads/7802fcb048934f876df0400a4638c2c5/Screenshot_from_2024-06-03_15-16-50.png)

**With temporal colocation:**

![Screenshot_from_2024-06-03_15-16-55](uploads/b74674db95135dab801e5d2b5f1f1b67/Screenshot_from_2024-06-03_15-16-55.png)

## Spatial colocation

When loading more than one species we may want to ensure that the stations that we have in our network/s measure all the species that we have selected. To do this, we need to activate the spatial colocation.

**Without spatial colocation:**

![no_spatial_colocation](uploads/5b547da5c0fa6da8473c611af5f30fbe/no_spatial_colocation.png)

**With spatial colocation:**

![spatial_colocation](uploads/ad31643908ef4be7c482de3ea678ebb7/spatial_colocation.png)

If you see that the number of stations is not equal, but has been reduced, it might be that the missing stations have NaN values for the current period.

In the dashboard, the spatial colocation is always on by default since we do not load multiple species at the same time, except when we use the multispecies filtering.