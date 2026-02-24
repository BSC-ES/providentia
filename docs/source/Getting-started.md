# Getting started

If you have access to the HPC machines at Barcelona Supercomputing Center (BSC), the first thing you need to decide is whether you want to use Providentia on a supercomputer (MN5 or Nord4) or on your local computer.

We recommend working on local machines to everyone, including the users at BSC, because the interactive features of the dashboard are faster and you do not need to wait in queue to get resources and use the software. The only disadvantage is that the data (models and observations) stored on HPC cannot be accessed directly and need to be downloaded onto your local machine using the [download](Download) mode in advance. If you do not want to download the data and instead you prefer to use an HPC machine for your analysis, we recommend reading the Wiki section [Connection setup](Connection-setup).

If you do not have access to the machines, you will only be able to use the download mode to get CAMS model data and observations from limited sources, i.e. Zenodo for GHOST and NILU Thredds for ACTRIS. If you want to use your own data, consider checking the tutorial [on how to format model data](https://github.com/BSC-ES/providentia/blob/master/tutorials/2.%20Formatting%20model%20data.ipynb) and reading the section [Create your own data network](Create-your-own-data-network) to process and create observational netCDF files that Providentia can read.

## Prerequisites

Providentia works best on Linux and macOS. For these machines, these are the prerequisites:

### 1. Git
Install Git by following the instructions here: [https://git-scm.com/install/linux](https://git-scm.com/install/linux)

### 2. Conda
Install conda. We recommend Miniconda because it is lightweight. You can download the .sh file **for Linux** from [the official website](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2) after creating an account. If you get an error and your VPN is active, make sure to deactivate it before downloading it.

From the terminal on the WSL run:

```bash
cd Downloads
bash Miniconda3-latest-Linux-x86_64.sh
```

The software has not been designed to work on Windows. If you are a Windows user you have three options:
- Run it from a Windows Subsystem for Linux (WSL) if you have Windows Pro. Read the instructions [here](WSL). (**Recommended**)
- Use a virtual machine software like Oracle VirtualBox to access Linux. Read the instructions [here](VM).
- Run it in Git Bash. Read the instructions [here](Git-Bash).

## Cloning the project

Use the following command to get a copy of the repository in your machine:

```
git clone https://github.com/bsc-es/providentia
```

When you have finished cloning the repository from Github, you are automatically in the branch `master`. It is recommended to use that branch as it contains the latest features and bug fixes.

## Running the tool the first time

Once cloned, you should be able to open the dashboard by running this command from your terminal:

```
cd providentia
./bin/providentia
```

The first time the software runs in a local machine it will create a conda environment called `providentia-env_v[version]` with all the modules needed. If you encountered any other problem, feel free to [contact us](Home).

In HPC, the environment is not created by the user as it is stored in a shared folder. Every time we run Providentia on HPC, a wall time of 2 hours is requested, with 12 CPUs and 30Gb of total memory. This can be modified as desired using the bash options. You can check the available options with:

```
./bin/providentia --usage
``` 

## Accessing the data

When you open the dashboard on a local machine the first time, you don't see anything on the dropdowns and you need to place the data into a local directory. By default, the data is read from `/home/{user}/data/providentia`. If for some reason you want to store it elsewhere you can edit the paths in `settings/data_paths.yaml`. 

### Data directory tree and filename conventions

The datasets need to be saved following a very specific directory tree. The download mode takes care of that when saving the files, more details can be found in the [download](Download) section. However, if you are using your own data you will need to take that into account.

By default, in the folder `/home/{user}/data/providentia` (or your preferred) there should be three folders: 
- `mod`: Interpolated model data as in: {GHOST version} -> {model}_{domain}_{ensemble} -> {resolution} -> {species} -> {network} -> {species}_{year}{month}.nc.
- `mod_to_interp`: Model data to interpolate as in: {model} -> {domain} -> {resolution} -> {species} -> {species}_{year}{month}.nc.
- `obs`: Observation datasets. For GHOST as in: ghost -> {network} -> {GHOST version} -> {resolution} -> {species} -> {species}_{year}{month}.nc. For non-GHOST as in: nonghost -> {provider} -> {network} -> {resolution} -> {species} -> {species}_{year}{month}.nc.

As observed, datasets must be saved per month, independently of their temporal resolution. An example of a working directory tree is the following:

```
├── mod
│   └── 1.5
│       └── cams61_emep_ph2-eu-000
│           └── hourly
│               └── sconcno2
│                   └── eea-eionet
│                       ├── sconcno2_201801.nc
│                       ├── sconcno2_201802.nc
│                       ├── sconcno2_201803.nc
│                       ├── sconcno2_201804.nc
│                       ├── sconcno2_201805.nc
│                       ├── sconcno2_201806.nc
│                       ├── sconcno2_201807.nc
│                       ├── sconcno2_201808.nc
│                       ├── sconcno2_201809.nc
│                       ├── sconcno2_201810.nc
│                       ├── sconcno2_201811.nc
│                       └── sconcno2_201812.nc
├── mod_to_interp
│   └── cams61_emep_ph2
│       └── eu
│           └── hourly
│               └── sconcno2
│                   ├── sconcno2_201801.nc
│                   ├── sconcno2_201802.nc
│                   ├── sconcno2_201803.nc
│                   ├── sconcno2_201804.nc
│                   ├── sconcno2_201805.nc
│                   ├── sconcno2_201806.nc
│                   ├── sconcno2_201807.nc
│                   ├── sconcno2_201808.nc
│                   ├── sconcno2_201809.nc
│                   ├── sconcno2_201810.nc
│                   ├── sconcno2_201811.nc
│                   └── sconcno2_201812.nc
└── obs
    ├── ghost
    │   └── EBAS
    │       └── 1.5
    │           └── hourly
    │               └── sconcno2
    │                   ├── sconcno2_201801.nc
    │                   ├── sconcno2_201802.nc
    │                   ├── sconcno2_201803.nc
    │                   ├── sconcno2_201804.nc
    │                   ├── sconcno2_201805.nc
    │                   ├── sconcno2_201806.nc
    │                   ├── sconcno2_201807.nc
    │                   ├── sconcno2_201808.nc
    │                   ├── sconcno2_201809.nc
    │                   ├── sconcno2_201810.nc
    │                   ├── sconcno2_201811.nc
    │                   └── sconcno2_201812.nc
    └── nonghost
        └── eea
            └── eionet
                └── hourly
                    └── sconcno2
                        ├── sconcno2_201801.nc
                        ├── sconcno2_201802.nc
                        ├── sconcno2_201803.nc
                        ├── sconcno2_201804.nc
                        ├── sconcno2_201805.nc
                        ├── sconcno2_201806.nc
                        ├── sconcno2_201807.nc
                        ├── sconcno2_201808.nc
                        ├── sconcno2_201809.nc
                        ├── sconcno2_201810.nc
                        ├── sconcno2_201811.nc
                        └── sconcno2_201812.nc
```

If you are running Providentia on HPC, you will already see that there are options to choose from in the menu on the top. The data is being read from the paths specified in `settings/data_paths.yaml`.

## Providentia internal directories

When cloning the Providentia repository from GitHub, the project is automatically created with a predefined directory structure. These internal directories are required for the correct execution of Providentia and are used to store configuration files, intermediate outputs, results and visualizations.

Providentia uses the following directories to store essential files during execution:

```
├── configurations/     Configuration files required to run all Providentia modes.
├── logs/               Output files generated during interpolation.
├── notebooks/          Template notebooks.
├── plots/              Saved plots when Providentia is used as a module.
├── reports/            Generated reports.
├── saved_data/         Configuration, NetCDF and NumPy files when used as a module.
└── settings/           Files that configure various aspects of Providentia.
```

The `settings/` directory contains the main configuration files that control Providentia’s behavior. These files can be modified by the user to adapt and configure Providentia’s behavior according to their needs.

```
├── basic_stats.yaml                 Defines the stats properties
├── color_palettes.yaml              Defines the possible color palettes
├── data_paths.yaml                  Defines dataset paths categorized by machine.              	
├── experiment_bias_stats.yaml       Defines the experiment stats properties
├── init_prov.yaml                   Stores initialization settings, including non-ghost available networks and resolutions.             
├── interp_experiments.yaml          Specifies locations for non-interpolated experiments. 
├── remove_extreme_stations.yaml     Defines criteria for filtering stations if you want to automatically remove them. 
├── exceedances.yaml                 Stores threshold values for exceedance statistics             	 
├── fairmode.yaml                    Stores configurations for the Fairmode plots.                
├── plot_characteristics.yaml        Configures plot appearance settings.	
└── report_plots.yaml                Defines plot types per report type.        
```

## Statistics

Before explaining how to use each mode, it is important to note that statistics are computed in numerous ways, depending on the user’s needs. A thorough explanation can be found in the [Statistics](Statistics) section.

## Launching the dashboard

As explained, you can launch the dashboard by simply running:

```
./bin/providentia
```

If you want to define which data is loaded in advance, you can use a configuration file. Some examples can be found under the folder `configurations`, for more details read the section [Configuration files](Configuration-files). Once you have it, you can specify its path in the command line with the argument `--config`:

```
./bin/providentia --config='/path/to/file/example.conf'
```

If you have multiple sections or subsections, a pop-up window will immediately appear where you can choose the section or subsection of interest. After that, the graphical window of Providentia will appear and you can begin using the tool.

An initial set of plots will be displayed, including the timeseries, distribution, statistics summary, and periodic plots. To take full advantage of Providentia, you can explore the wide range of plotting options described in [Plot types and options](Plot-types-and-options). We also recommend reading the [Plot customisation](Plot-customisation) section.

More details can be found in the [dashboard section](Dashboard).

## Generating a report

With the configuration file you can also generate PDF reports. In order to do this, you should use the argument `report`:

```
./bin/providentia --config=/path/to/file/example.conf --report
```

If your configuration file is inside the folder `configurations`, you don't need to specify the full path:

```
./bin/providentia --config=example.conf --report
```

You can launch the dashboard or get a report for only one section by using the option  `--section`. In order to indicate subsections, you will need to write the section name, followed by an interpunct (·) and the subsection name.

```
./bin/providentia --config=example.conf --report --section=All·France
```

The reports will be saved under the folder `reports`. You can add a path in the `report_filename` of the configuration file to change the default directory.

More details can be found in the [report section](Report).

## Using Providentia backend functions

Providentia can be imported and used in your own Python scripts. Some examples on how to use Providentia's backend functions can be found in the [tutorials folder](https://github.com/BSC-ES/providentia/tree/master/tutorials).

Also, a Jupyter notebook with an active conda environment can be launched with the following command:

```
./bin/providentia --notebook
```

More details can be found in the [library section](Library).

## Interpolating your model data to observations

If you want to visualise data from your model, you will need to interpolate it to the network. Using a configuration file, you can start interpolating your model data to your desired observational network.
```
./bin/providentia --config=example.conf --interpolate
```

More details can be found in the [interpolation section](Interpolation).

Enjoy!