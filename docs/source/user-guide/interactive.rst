================
Interactive mode
================

The interactive mode was designed for users to exploit the powerful backend of Providentia, allowing the use of the Providentia functions, as well as allowing for the running in Jupyter notebooks. 

Providentia interactive can be used in any script, through simply importing the Interactive class as follows:   

::

   from providentia import Interactive

Using Providentia in your own scripts
-------------------------------------

If your script is not inside the Providentia directory, in order to import from Providentia, then it is necessary to add these lines to your code first:            

::
   
   import sys
   sys.path.append(provdir)

      
where ``provdir`` is the path where your Providentia code exists.

You will also need to load the necessary modules for it to function, for this you can run:

::
   
   source provdir/load_modules.sh

After the import you will have full access to the backend of Providentia and all of its functions as described in the following sections.

NOTE: In the future we will allow for Providentia to be imported as a module on the BSC machines, making this step redundant. 

Starting a Jupyter notebook
---------------------------

A Jupyter notebook that can run Providentia interactive can also be now easily launched with the following command:

::

   ./bin/providentia --interactive         

This will either start a job if you are using a machine with SLURM, or directly open Jupyter notebooks. In the first case, a file named ``interactive.out`` will be created. Information in the file will then allow for a connection to a Jupyter notebook to be set up. Firstly, an SSH tunnel from the local machine needs to be set up, by pasting a given command into a Linux/MacOS terminal or equivalent (i.e. PuTTY on Windows), e.g.:

::

   Create an SSH tunnel via terminal on your local machine:
   ssh -N -L 8825:s04r1b14:8825 bsc32025@nord4.bsc.es

Secondly, the notebook can than be accessed via a web browser following information from the bottom of the file, e.g.:     

::

   To access the server, open this file in a browser:
         file:///.statelite/tmpfs/gpfs/home/bsc32/bsc32025/.local/share/jupyter/runtime/jpserver-800644-open.html
      Or copy and paste one of these URLs:
         http://s04r1b14:8825/lab?token=57a994dfaecf92e5aefd1d5dbc5a9831fb8060f873fd8181
      or http://127.0.0.1:8825/lab?token=57a994dfaecf92e5aefd1d5dbc5a9831fb8060f873fd8181

One additional step is also needed inside the notebook if wishing to embed plots:     

::

   %matplotlib inline
              
This is a magic function that allows for the rendering of figures directly in the notebook.        

NOTE: This must set after the importing of the Interactive class.

A template Jupyter notebook, demonstrating some of the Providentia interactive features can be found in ``notebooks/interactive_template.ipynb``.

Features
--------

The interactive version gives access to the entirety of the Providentia backend. This therefore includes reading, filtering, plotting, writing, etc. After importing the Providentia Interactive class, the class methods are then accessible. In the following sections, each of these methods will be described. 

Reading and filtering data from a .conf file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As a first step, data needs to be read and filtered appropriately.

In order to let Providentia know what data is wanted to be read, and what filters to apply, this is all defined in a .conf file (the same as used in the dashboard and offline modes). A class instance is created by the following line: 

::
   
   provi = Interactive(conf='interactive_template.conf')

where ``provi`` is the class instance, which can access all of the class variables and methods.  

Only one section-subsection pair in the .conf can be read and filtered by Providentia interactive at one time, with the exception in the creation of plots with the ``multispecies`` option (see plotting section). Where there is more than one pair defined in the file, the specific pair wished to be read can be set by passing the section and/or subsection arguments, e.g.:

::

   provi = Interactive(conf='interactive_template.conf', section='OTHERSECTION', subsection='OTHERSUBSECTION')

If there are multiple section and subsection pairs in the .conf file, and a specific pair is not explicitly set to be read, then simply the first pair is taken to be read.

When wanting to change any data or filters applied to the data, simply update the appropriate the .conf and reinitiate the class instance. 

If wanting to overwrite any arguments in the .conf file, directly in the script, each argument can be simply passed when initiating the class instance, e.g.:

::

   provi = Interactive(conf='interactive_template.conf', network='EANET')


Reset filter
^^^^^^^^^^^^

After data has been read and filtered, it will stay filtered until it is reset. This can be done by the following method: 

::

   provi.reset_filter()  

If at any point wanting to return to the state when the Interactive class was initialised, this can be done by adding the following argument to the method:

::

   initialise=True 

Applying specific filter
^^^^^^^^^^^^^^^^^^^^^^^^

If wanting to apply a filter not set in the .conf, this can be done using the following method:

::

   provi.apply_filter(field, ...)

where ``field`` is the field to filter by. The fields to filter by can be representativity fields, period fields or metadata fields.

If the field is numeric, lower and upper limits to retain data between can be set as follows:

::
   
   provi.apply_filter(field, lower=28, upper=31)


NOTE: it is not mandatory to pass ``lower`` and ``upper`` arguments together.

If the field is textual, values to keep and remove associated data with can be set as follows:

::

   provi.apply_filter(field, keep='Spain', remove='')


If multiple values are wanted to be removed concurrently, the arguments passed should be lists.

::
   
   provi.apply_filter(field, remove=['Spain','France'])


NOTE: it is not mandatory to pass ``keep`` and ``remove`` arguments together.

For the specific case of representativity fields, the argument passed should be ``limit``:

::
   
   provi.apply_filter(rep_field, limit=20)


Filtering for specific station/s
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If wanting to filter data for a specific station or stations, the following convenience method can also be used (rather than using ``apply_filter``):

::

   provi.select_station(station)


where ``station`` can be a str ``station_reference`` of one station of interest, or can be a list of multiple stations of interest. This will then subsequently mean the data in memory is filtered for the relevant station or stations.

Print active .conf file
^^^^^^^^^^^^^^^^^^^^^^^

If wanting to visually print the .conf file that is currently active this can be done by:

::

   provi.print_config()


If wanting to visually print any other existing .conf, this can be done by passing the ``conf`` or ``config`` argument: 

::

   conf='important.conf'

Calculate statistics
^^^^^^^^^^^^^^^^^^^^

In order to calculate statistics for a dataset, this can done via the following method:

::
   
   stat_calc = provi.calculate_stat(stat, labela='OBS')

where ``stat`` is the statistic wished to be calculated, and ``labela`` is the name of the observations/experiment data (this can be an alias set in the .conf or the original dataset name). The statistic returned will be one summary value, which can be formulated in differing ways, e.g. taking a median time series across all stations and then calculating the statistic, or calculating the statistic across all stations and then taking the median. This formulation can be set by changing the statistical mode in the .conf, for which more information can be found [here](Statistics).

If wanting to calculate a bias statistic, then an additional dataset is needed to be compared against, which can be set via the following argument:

::

   labelb='EMEP'

where ``EMEP`` in this case is the name of the dataset that is wished to be compared against.

For bias statistics for which a subtraction is involved, it is always done as ``datasetb - dataseta``.

If wanting to calculate statistics at each individual station then this can be set by adding the following argument:

::

   per_station=True

Saving data
^^^^^^^^^^^

The data which has been read and filtered can be saved out using the following method: 

::

   provi.save(format='nc')

where ``format`` is the format of the saved data, and can be ``nc`` (netCDF), ``np`` (numpy) or ``conf`` (Providentia .conf file). In the case of ``conf``, rather than the data being saved, a Providentia .conf file is generated, which when loaded would return the exact same data in memory.

The filename of the saved data can be set by passing the following argument:

::

   fname='/mypath/myfilename'

If ``fname`` is not provided then one is generated automatically, and saved in the directory ``saved_data``.
 
Accessing data in memory
^^^^^^^^^^^^^^^^^^^^^^^^

As well as being able to save data out, data can be returned directly in memory in specific formats, via the following method: 

::

   data = provi.get_data(format='xr')

where ``format`` is the format of the returned data, and can be ``nc`` (netCDF), ``np`` (numpy) or ``xr`` (xarray).

The returned variable will contain all read data and metadata variables available. Therefore, if reading a GHOST network you can expect the metadata to be more substantial.

Accessing specific variable in memory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If wanting to extract a specific variable in memory, rather than the entire read dataset, this can be done using the following method:

::

   var_data = provi.get_var(var='myvar')

where ``var`` is the name of the variable wished to be read.

Plotting 
^^^^^^^^

All plot types that are available in other modes of Providentia are available in the interactive mode, with the addition of the option to plot the legend as a standalone plot. All plots can be made using the method:

::

   provi.make_plot(plot_type)

where ``plot_type`` is the plot type wished to be made (all types stated later in this section).

All plots have been tailored to appear nice in Jupyter notebooks, but any plot settings can be modified by the user, per plot type, in:      
``settings/plot_characteristics_interactive.yaml``.

Plot customisation
^^^^^^^^^^^^^^^^^^

There are numerous arguments to the plot method that can be used to customise the plots, each of which are subsequently detailed.

**Legend**

For user convenience a legend has been integrated into each appropriate plot type. If this is not wanted, then this can be deactivated using following argument to the plot method: 

::

   legend=False

If wanting to remove the observations data label from the legend, then this can be done by passing the following argument:

::

   set_obs_legend=False

**Colourbar**

For some plots a colourbar is integrated into the plot, e.g. map. This can be deactivated by passing the following argument:

::

   cb=False

**Limiting plotted data**

By default, all read data will be plotted. If wanting to limit which data is plotted, this can be done by passing a list of the relevant data labels to be plotted, e.g.:

::

   data_labels=['EMEP', 'MONARCH']

NOTE: these labels should be the dataset aliases, if set.

**Setting title, x or y labels**

The title, xlabel and ylabel will be set automatically (if required) for all plot types, however if wanting to overwrite them this can be done by passing the following arguments:

::

   title='My custom title'        
   xlabel='My custom xlabel'        
   ylabel='My custom ylabel'         

**Custom formatting**

To have more granular control of the plot formatting, and overwrite specific variables set in ``settings/plot_characteristics_interactive.yaml``, the ``format`` argument can be used to pass a dictionary which sets specific plot type format variables to overwrite, e.g.:

::

   format={"figsize": [14,7], "xtick_params": {"labelsize": 22}}

**Plot options**

As in the offline version, plot options can be set for each plot type. There are 2 ways to set the plot options. One way is passing them all through a list, e.g.:          

::

   plot_options=['annotate', 'bias']

The other way is to pass each individual plot option as an argument e.g.:

::

   annotate=True

The available plot options per plot type are given in the sections below.

**Saving / returning plot object in memory**

Rather than viewing the plot on the screen, it can be returned in memory, or saved to a file.

To return the plot in memory, this can be done as follows:

::

   plot_obj = provi.make_plot(plot_type, return_plot=True)

where ``plot_obj`` is the object of the plot returned.

To save the plot, the following argument can be passed:

::

   save='fname'

where ``save`` is the name of the file to save to. If save is set as: ``save=True``, then the filename will be automatically generated and saved in the ``plots`` directory.

Plot types
^^^^^^^^^^

In this section the available plot types will de detailed, with any specific subtleties involved in calling the plotting method detailed. Available plot options are also stated per plot type. 

**legend**
The legend plot can be made as follows:   

::
   
   provi.make_plot('legend')

It has no available plot options.

**metadata**

The metadata plot can be made as follows:       

::

   provi.make_plot('metadata')

It has no available plot options.

**map**

The map plot can be made as follows:   

::

   provi.make_plot('map-stat', labela='OBS')

where stat is the statistic to be plotted per station on the map, and ``labela`` is the data label of the dataset to calculate the statistics with.

By adding `labelb` to the plot arguments, you can turn the plot type into a bias plot:

::

   provi.make_plot('map-stat', labela='OBS', labelb='MONARCH')

For bias statistics for which a subtraction is involved, it is always done as ``datasetb - dataseta``.

The map extent is by default automatically set based on what is plotted on the map. This can be also manually controlled by passing the argument:

::

   map_extent = [lonmin, lonmax, latmin, latmax]

where ``lonmin``, ``lonmax``, ``latmin``, and ``latmax`` are the rectangular bounds of the map extent, e.g.: ``map_extent = [-180, 180, -90, 90]``.

The available plot options are:      
``annotate``, ``domain``.
 
**timeseries**

The timeseries plot can be made as follows:   

::

   provi.make_plot('timeseries')

The available plot options are:      
``annotate``, ``bias``, ``hidedata``, ``logy``, ``smooth``

**periodic**

The periodic plot can be made as follows:   

::

   provi.make_plot('periodic-stat')

where ``stat`` is the statistic that is wanted to be plotted. 

The available plot options are:      
``annotate``, ``bias``, ``logy``

**periodic-violin**

The periodic-violin plot can be made as follows: 

::

   provi.make_plot('periodic-violin')


The available plot options are:      
``annotate``, ``logy``

**distribution**

The distribution plot can be made as follows: 

::

   provi.make_plot('distribution')


The available plot options are:   
``annotate``, ``bias``, ``logx``, ``logy``

**scatter**

The scatter plot can be made as follows: 

::

   provi.make_plot('scatter')

The available plot options are:   
``annotate``, ``hidedata``, ``logx``, ``logy``, ``regression``

**boxplot**

The boxplot can be made as follows: 

::

   provi.make_plot('boxplot')

The available plot options are:   
``annotate``, ``logy``, ``multispecies``

**heatmap**

The heatmap plot can be made as follows: 

::

   provi.make_plot('heatmap-stat')

where ``stat`` is the statistic that is wanted to be plotted. 

The available plot options are:   
``annotate``, ``bias``, ``multispecies``

**table**

The table plot can be made as follows: 

::

   provi.make_plot('table')

The available plot options are:   
``bias``, ``multispecies``

**statsummary**

The statsummary plot can be made as follows: 

::

   provi.make_plot('statsummary')

The available plot options are:   
``bias``, ``multispecies``

Functions
---------

.. automodule:: providentia.interactive
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: read, filter, make_colourbar, make_legend, set_config