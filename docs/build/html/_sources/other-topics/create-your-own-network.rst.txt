Create your own data network
============================

If you have data from observations in an Excel or CSV file and want to create netCDF files that Providentia can visualise, you can use NES, a library created by the Earth Sciences Department to easily read and write netCDF files containing geospatial data in different projections (regular, rotated, Mercator, etc.).

You can install it by cloning the repository:

::

    git clone https://earth.bsc.es/gitlab/es/NES.git

And, in order to use the latest features, go to the master branch:

::

    git checkout master

You can learn more about how to create your own datasets with the following tutorials:

- Create hourly datasets from Port Barcelona: https://earth.bsc.es/gitlab/es/NES/-/blob/master/tutorials/2.Creation/2.4.Create_Points_Port_Barcelona.ipynb
- Create monthly datasets from CSIC: https://earth.bsc.es/gitlab/es/NES/-/blob/master/tutorials/2.Creation/2.5.Create_Points_CSIC.ipynb

Providentia will be able to find these datasets if they are saved in the directory /esarchive/obs/. The data path must be included in providentia/settings/nonghost_files.json. As an example, the path to the data from CSIC (/esarchive/obs/csic/csic/monthly/sconcnh3/) correspond in this dictionary to:

::

    "csic/csic": {
        "monthly": {
            "sconcnh3": []
            }
        }

If you launch the dashboard and select any period in 2019, you will see that you can select CSIC as a network. If you want to use a configuration file to launch it or create offline reports, then you will need to write:

::
    
    network = csic/csic