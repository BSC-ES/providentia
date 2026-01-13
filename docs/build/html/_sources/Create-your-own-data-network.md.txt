# Create your own data network

If you have data from observations in an Excel or CSV file and want to create netCDF files that Providentia can visualise, you can use [NES](https://pypi.org/project/nes/), a library created by the Earth Sciences Department to easily read and write netCDF files containing geospatial data in different projections (regular, rotated, Mercator, etc.).

You can install it by cloning the repository:

```
pip install nes
```

You can learn how to create your own datasets with NES following the tutorials `5. Create hourly datasets from Port Barcelona` and `6. Create monthly datasets from CSIC` in the [tutorials folder](https://github.com/BSC-ES/providentia/tree/master/tutorials)

Providentia will be able to find these datasets if they are stored in the path specified in the `nonghost_root` key under the corresponding machine in the `settings/data_paths.yaml` file. The name of the network must be included in the `nonghost_available_networks` list in `settings/available_inputs.yaml`.