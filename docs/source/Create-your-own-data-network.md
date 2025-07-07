# Create your own data network

If you have data from observations in an Excel or CSV file and want to create netCDF files that Providentia can visualise, you can use [NES](https://earth.bsc.es/gitlab/es/NES/-/wikis/home), a library created by the Earth Sciences Department to easily read and write netCDF files containing geospatial data in different projections (regular, rotated, Mercator, etc.).

You can install it by cloning the repository:

```
git clone https://earth.bsc.es/gitlab/es/NES.git
```

And, in order to use the latest features, go to the master branch:

```
git checkout master
```

You can learn more about how to create your own datasets with the following tutorials:

- Create hourly datasets from Port Barcelona: https://earth.bsc.es/gitlab/es/NES/-/blob/master/tutorials/2.Creation/2.4.Create_Points_Port_Barcelona.ipynb

- Create monthly datasets from CSIC:
https://earth.bsc.es/gitlab/es/NES/-/blob/master/tutorials/2.Creation/2.5.Create_Points_CSIC.ipynb

Providentia will be able to find these datasets if they are stored in the path specified in the `nonghost_root` key under the corresponding machine in the `settings/data_paths.yaml` file. The name of the network must be included in the `nonghost_available_networks` list in `settings/init_prov.yaml`. For example, in the nord3v2 machine, to get the data from CSIC stored at this location: `/esarchive/obs/csic/csic/monthly/sconcnh3/`, you would need to modify the configuration files as follows.

Update `data_paths.yaml` so it stores the dataset directory:

```
"nord3v2": {
        "nonghost_root": "/esarchive/obs",
        ...
    },

```

Update `init_prov.yaml` by adding the corresponding data network to the `nonghost_available_networks` list:

```
'nonghost_available_networks': [..., 'csic/csic']

```
If you then launch the dashboard and select any period in 2019, you will see that you can select CSIC as a network. If you want to use a configuration file to launch it or create reports, then you will need to write:

```
network = csic/csic

```