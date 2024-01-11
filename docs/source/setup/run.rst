Run the tool
============

Running the tool on a BSC machine
---------------------------------

If you are on Nord3v2 or Marenostrum4, you should request an interactive session:

::

    ./bin/providentia

After some seconds you will have entered onto an allocated node and the dashboard will be launched.
By default, a wall time of 2 hours is requested, with 12 CPUs and 30Gb of total memory. This can be modified as desired using the bash options. You can check the available options with:

::

    ./bin/providentia -u

In order to run the tool, it is necessary to load the modules. These have already been loaded using the first command.

**Launch the dashboard**

As explained, you can launch the dashboard by simply running:

::

    ./bin/providentia


If you already have a configuration file, you can specify its path in the command line with the argument --config:

::
    
    ./bin/providentia --config='/path/to/file/example.conf'


A pop-up window will immediately appear where you can choose the section or subsection of interest. After that, the graphical window of Providentia will appear and you can begin using the tool.

**Generate an offline report**

If you have your configuration file ready, you can generate reports by running Providentia with the following command:

::
    
    ./bin/providentia --config='/path/to/file/example.conf' --offline


The mandatory commands are:

* --config: specify the path of your configuration file
* --offline: for running offline, without launching the dashboard

You can also launch the dashboard or get a report for only one section by using the option  --section. In order to indicate subsections, you will need to write the section name, followed by a hyphen (-) and the subsection name.

The reports will be saved under the folder reports. You can add a path in the report_filename of the configuration file to change the default directory.

Running the tool on a personal machine
--------------------------------------

If you do not have access to the BSC machines, you will need to define the directories where your data is stored. You can do this by editing the file configuration.py and defining ghost_root, nonghost_root and exp_root. For instance:

::
    
    'local_ghost_root': '/data/providentia/obs/ghost/'
    'local_nonghost_root': '/data/providentia/obs/nonghost/'
    'local_exp_root': '/data/providentia/exp/'

You should download the data that you need from the BSC systems to the local machine. The tool has set locations where it grabs the observational/experiment data on the **esarchive**:

* GHOST observational data: `/esarchive/obs/ghost`
* Non-GHOST observational data: `/esarchive/obs/`
* Interpolated experiment data: `/esarchive/recon/exp_interp`

The entire data collection comprises **several terabytes**, therefore it is not recommended to download the entire archive.

If you are running the tool on a personal machine, you will also need to install all the modules by yourself. To do this you can create a virtual environment with Conda and install the necessary packages by:

::
    
    conda create -n providentia-env python=3.9
    conda activate providentia-env
    conda install -c conda-forge cartopy
    conda install -c conda-forge jupyterlab
    pip install -r requirements.txt

Once the environment is activated and all modules are installed, you can open the dashboard and generate reports as explained above.

Enjoy!