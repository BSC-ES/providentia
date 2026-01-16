# Saved file formats

Currently, it is possible to export three different types of files: configuration files, Numpy and NetCDF files.

## .conf file format (Providentia)

Users can download their configurations (.conf) and use the corresponding files to launch Providentia again. You can use your configuration files in the dashboard, report, interpolation and download modes through the command line or make use of the load button.

## Numeric file formats

The following table summarizes the variables that are available after exporting data into numpy (.npz) files and NetCDF (.nc) files.

As Providentia is capable of loading multiple network and species simultaneously, some of these variables are provided per [network]-[species].

### Numpy (.npz file format)

#### Variables

| Variable | Description |
| ------ | ------ |
| [network]-[species]_data | Values of the desired species for both observations and models |
| [network]-[species]_ghost_data | GHOST data variables used for additional filtering | 
| [network]-[species]_metadata | Metadata of the observations which varies per month gives as a multidimensional array |
| time | Time in given resolution from the start date |
| data_labels | Labels associated with each data array, e.g. observations, model_1, etc. |
| ghost_data_variables | The names of the GHOST data variables used for additional filtering |
| resolution | Temporal resolution of data |
| start_date | Start date of data |
| end_date | End date of data |
| temporal_colocation | Boolean stating if observations and models have been temporally colocated |
| spatial_colocation | Boolean stating if data has been spatially colocated across [network]-[species] |
| filter_species | Data ranges per species used filter read data |
| ghost_version | Version of GHOST |

#### Loading the data

Loading a .npz file in python is done simply by:

```python
In [1]: import numpy as np                                                                           
In [2]: obs = np.load("/home/bsc32/bsc32099/PRV_sconco3_20160101_20160601.npz", allow_pickle=True)
```
Note it is necessary that the `allow_pickle` argument is set as True.

To investigate the variables that the loaded .npz has inside it, we can use the "files" method:
                 
```python
In [3]: obs.files                                                                                    
Out[3]: ['EBAS-sconco3_ghost_data', 'EBAS-sconco3_data', 'EBAS-sconco3_metadata'...]
```

Values for a data variable are returned by: 

```python
In [4]: data = obs['EBAS-sconco3_data']                                                                    
```

Metadata access is special in the .npz files. The metadata variable names can be returned by:

```python
In [5]: metadata_vars = obs['metadata'].dtype.names
```

Any specific metadata field can be accessed by using one of the metadata variable names: 

```python
In [6]: latitude = obs['metadata']['latitude']
```

### NetCDF (.nc file format)

#### Variables

| Variable | Description |
| ------ | ------ |
| [network]-[species]_data | Values of the desired species for both observations and models |
| [network]-[species]_ghost_data | GHOST data variables used for additional filtering | 
| [network]-[species]_[metadata_var] | Metadata of the observations which varies per month given per variable |
| [network]-[species]_qa | Quality assurance flags, GHOST performed quality control checks |
| [network]-[species]_flags | Data flags, standardised flags taken from the data provider |
| time | Time in given resolution from the start date |
| data_labels | Labels associated with each data array, e.g. observations, model_1, etc. |
| ghost_data_variables | The names of the GHOST data variables used for additional filtering |

resolution, start_date, end_date, temporal_colocation, spatial_colocation, filter_species and ghost_version are stored as attributes of [network]-[species]_data.

#### Loading the data

You can read these files as you would usually do, typically using the library netCDF4:

``` python
from netCDF4 import Dataset 
dataset = Dataset("PRV_sconco3_20160101_20170101.nc")
```

Or xarray:
``` python
import xarray as xr 
dataset = xr.open_dataset("PRV_sconco3_20160101_20170101.nc")
```