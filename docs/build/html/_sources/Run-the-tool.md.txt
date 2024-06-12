# Run the tool

## Running the tool on a BSC machine

If you are on Nord3v2 or Marenostrum 5, you should request an interactive session:

```
./bin/providentia
``` 

After some seconds you will have entered onto an allocated node and the dashboard will be launched. If you want to launch it multiple times and avoid waiting in the queue you can use the debug mode. To do this, you will need to run the tool using the `--debug` flag and then rerun it:

```
./bin/providentia --debug
./bin/providentia
```

By default, a wall time of 2 hours is requested, with 12 CPUs and 30Gb of total memory. This can be modified as desired using the bash options. You can check the available options with:

```
./bin/providentia --usage
``` 

### Launch the dashboard

As explained, you can launch the dashboard by simply running:

```
./bin/providentia
```

If you already have a [configuration file](Configuration-files), you can specify its path in the command line with the argument `--config`:

```
./bin/providentia --config='/path/to/file/example.conf'
```

A pop-up window will immediately appear where you can choose the section or subsection of interest. After that, the graphical window of Providentia will appear and you can begin using the tool. 

### Generate an offline report

If you have your [configuration file](Configuration-files) ready, you can generate reports by running Providentia with the following command:

```
./bin/providentia --config='/path/to/file/example.conf' --offline
```

The mandatory commands are:

* `--config`: specify the path of your [configuration file](Configuration-files)
* `--offline`: for running offline, without launching the dashboard

You can also launch the dashboard or get a report for only one section by using the option  `--section`. In order to indicate subsections, you will need to write the section name, followed by a hyphen (-) and the subsection name.

The reports will be saved under the folder `reports`. You can add a path in the `report_filename` of the configuration file to change the default directory.

### Using interactive mode

A Jupyter notebook can be launched with the following command:
```
./bin/providentia --interactive
```        

## Running the tool on a personal machine

If you do not have access to the BSC machines, you will need to define the directories where your data is stored. You can do this by editing the file `settings/data_paths.yaml` and defining `ghost_root`, `nonghost_root` and `exp_root` under `local`. For instance:

```
'ghost_root': '/data/providentia/obs/ghost/'
'nonghost_root': '/data/providentia/obs/nonghost/'
'exp_root': '/data/providentia/exp/'
```

You should download the data that you need from the BSC systems to the local machine. To do this, you can use the [download mode](Download).

If it is the first time that you launch Providentia, an environment called `providentia-env` will be created with all the necessary dependencies. If for some reason you want to create the environment from scratch, you can use:

```
conda env create -f environment.yaml
```

You might get a warning like:

```
WARNING conda.models.version:get_matcher(556): Using .* with relational operator is superfluous and deprecated and will be removed in a future version of conda. Your spec was 1.6.0.*, but conda is ignoring the .* and treating it as 1.6.0
```

This can be removed by updating conda:

```
conda update conda
conda install -n base conda=24.4.0 conda-build=24.3.0
```

Check what the latest versions of [conda](https://github.com/conda/conda/releases) and [conda-build](https://github.com/conda/conda-build/releases) are.

Enjoy!