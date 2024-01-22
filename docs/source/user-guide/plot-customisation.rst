Plot customisation
==================

Editing the plot style from the configuration files
---------------------------------------------------

If you want to edit the style of your plots, you will need to edit it in one of the following files and launch the tool again:

- Plots in the dashboard: ``/providentia/settings/plot_characteristics_dashboard.json``
- Plots in the offline reports: ``/providentia/settings/plot_characteristics_offline.json``

Setting custom bounds and cmap per species
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users can define the color and bounds of the colorbar (cmap, vmin and vmax) per species using a dictionary, with the keys being the names of the species inside ``basic_stats.json`` and ``experiment_bias_stats.json``. An example can be seen in the image below:

::

  "Mean":        {"function": "calculate_mean", 
                  "order": 0, 
                  "label": "Mean", 
                  "arguments": {}, 
                  "units": "[measurement_units]", 
                  "minimum_bias": [0.0],
                  "vmin_absolute": {"sconco3": 0, "sconcno2": 0},
                  "vmax_absolute": {"sconco3": 20, "sconcno2": 5}, 
                  "vmin_bias": {}, 
                  "vmax_bias": {},
                  "cmap_absolute": "viridis",
                  "cmap_bias": "RdYlBu_r"},

If they define the cmap, they will need to give a complete list of cmap options for each of the species that they load or otherwise a warning will appear. For vmin and vmax, they can define the bounds for some species and the rest will take the data minimum and maximum values.

Remove extreme stations by their statistical values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to automatically remove stations that have certain statistical values, you will need to add your criteria in the file `remove_extreme_stations.json`. An example of this exists for CAMS:

::

  "CAMS": {"r": ["<0.3"],
           "NMB": ["<-100.0", ">100.0"],
           "NRMSE": [">100.0"]}

The statistics can be general, across all components, or they can be specific per component, for example:

::

  "CAMS": {"sconco3": {"r": ["<0.3"],
                       "NMB": ["<-100.0", ">100.0"],
                       "NRMSE": [">100.0"]},
           "sconcno2": {"r": ["<0.55"],
                        "NMB": ["<-20.0", ">20.0"],
                        "NRMSE": [">200.0"]}}

You will also need to add the variable `remove_extreme_stations` in your configuration file:

::

  remove_extreme_stations=CAMS

Calculating exceedances
^^^^^^^^^^^^^^^^^^^^^^^

In Providentia ``exceedances`` is available in the list of available statistics. How it is currently implemented is simplistic, but users can simply state a threshold/limit value per component, and each instance where values exceed this limit will be counted. Therefore the exceeedances statistic simply gives the number of instances. The threshold values can be set in the file ``settings/exceedances.json``, as so:

::
  
  {"sconco3": 90.21, 
   "sconcno2": 106.38}

Editing the plot style in the dashboard
---------------------------------------

Changing the plot style
^^^^^^^^^^^^^^^^^^^^^^^

The style of the plots can be edited by clicking on the burger menus and changing the settings.

.. image:: ../images/plot-customization/burger-menu.png
  :alt: Burger menu

Legend picking
^^^^^^^^^^^^^^

Clicking on the legend labels will remove or add data to each of the plots. If the label appears in bold, the data will be visible. If not, it will disappear.

.. image:: ../images/plot-customization/legend-picking.png
  :alt: Legend picking

Changing the statistics
^^^^^^^^^^^^^^^^^^^^^^^

The statistics in the ``statsummary`` can be updated from the burger menu.

.. image:: ../images/plot-customization/statistics-change.png
  :alt: Changing the statistics

Information on hover
^^^^^^^^^^^^^^^^^^^^
Most plots show information when hovering over them. Take a look for instance at the distribution plot:

.. image:: ../images/plot-customization/info-hover.png
  :alt: Information on hover

Smoothing
^^^^^^^^^
It is possible to add a smoothing line to the timeseries plot and make the points disappear. In order to achieve this, you will need to increase the smooth window, which by default is 0 and bring the marker size down to 0. You can also use the plot option ``hidedata`` to hide the points.

.. image:: ../images/plot-customization/smoothing.png
  :alt: Smooth line
