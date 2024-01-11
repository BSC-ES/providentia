Saved file formats
==================

Currently, it is possible to export three different types of files: configuration files and Numpy and NetCDF files.

.conf file format (Providentia)
-------------------------------

Users can download their configurations (.conf) and use the corresponding files to create offline reports or launch Providentia again. In Run the tool you can learn how to load your configurations through the command line or from the dashboard.

Saved variables (Numpy and NetCDF)
----------------------------------

The following table details the available variables that are available after exporting data in numpy (.npz) files and NetCDF (.nc) files.

As Providentia is capable of loading multiple network and species simultaneously, some of these variables are provided per [network]-[species].

.. list-table:: 
   :widths: 20 20 20 40
   :header-rows: 1

   * - Variable
     - .npz
     - .nc
     - Description
   * - [network]-[species]_data
     - ✓
     - ✓ 
     - Values of the desired species for both observations and experiments
   * - [network]-[species]_ghost_data
     - ✓
     - ✓
     - GHOST data variables used for additional filtering
   * - [network]-[species]_metadata
     - ✓
     - ✗
     - Metadata of the observations which varies per month (given per variable in .nc files)
   * - time
     - ✓
     - ✓
     - Time in given resolution from the start date
   * - data_labels
     - ✓
     - ✓
     - Labels associated with each data array, e.g. observations, experiment_1, etc.
   * - ghost_data_variables
     - ✓
     - ✓
     - The names of the GHOST data variables used for additional filtering
   * - resolution
     - ✓
     - ✗
     - Temporal resolution of data (attribute of [network]-[species]_ghost_data in .nc files)
   * - start_date
     - ✓
     - ✗
     - Start date of data (attribute of [network]-[species]_ghost_data in .nc files)
   * - end_date
     - ✓
     - ✗
     - End date of data (attribute of [network]-[species]_ghost_data in .nc files)
   * - temporal_colocation
     - ✓
     - ✗
     - Boolean stating if observations and experiments have been temporally colocated (attribute of [network]-[species]_ghost_data in .nc files)
   * - spatial_colocation
     - ✓
     - ✗
     - Boolean stating if data has been spatially colocated across [network]-[species] (attribute of [network]-[species]_ghost_data in .nc files)
   * - filter_species
     - ✓
     - ✗
     - Data ranges per species used filter read data (attribute of [network]-[species]_ghost_data in .nc files)
   * - ghost_version
     - ✓
     - ✗
     - Version of GHOST (attribute of [network]-[species]_ghost_data in .nc files)

Note for metadata, in the .nc files the information is saved per variable (e.g. ), whereas in the .npz files it is saved in multidimensional array.


.npz file format (Numpy)
------------------------

Loading the data
^^^^^^^^^^^^^^^^

Loading a .npz file in python is done simply by:

::

    In [1]: import numpy as np                                                                           
    In [2]: obs = np.load("/home/bsc32/bsc32099/PRV_sconco3_20160101_20160601.npz", allow_pickle=True)

Note It is necessary that the allow_pickle argument is set as True.

To investigate the variables that the loaded .npz has inside it, we can use the "files" method:

::

    In [3]: obs.files                                                                                    
    Out[3]: ['EBAS-sconco3_ghost_data', 'EBAS-sconco3_data', 'EBAS-sconco3_metadata'...]

Accessing the data
^^^^^^^^^^^^^^^^^^
Values for a data variable are returned by:

::

    In [4]: data = obs['EBAS-sconco3_data']                                                                    

Accessing the metadata
^^^^^^^^^^^^^^^^^^^^^^

Metadata access is a special in the .npz files. The metadata variable names can be returned by:

::

    In [5]: metadata_vars = obs['metadata'].dtype.names

Any specific metadata field can be accessed by using one of the metadata variable names:

::

    In [6]: latitude = obs['metadata']['latitude']

.nc file format (NetCDF)
------------------------

You can read these files as you would usually do:

::

    ncdf = Dataset("PRV_sconco3_20160101_20170101.nc")