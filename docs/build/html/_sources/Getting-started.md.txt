# Getting started

The first thing you need to decide is whether you want to use Providentia on HPC or on your local computer. 

We always recommend working on local machines because the interactive features of the dashboard are faster and you do not need to wait in queue to use the software. The only disadvantage is that the data (experiments and observations) stored on HPC cannot be accessed. To solve this we developed the [download mode](Download), which allows you to download data from HPC directly onto your local machine.

If you do not want to download the data and instead you prefer to use an HPC machine for your analysis, we recommend reading the Wiki section [Connection setup](Connection-setup)

## Cloning the project

Independently of the machine, use the following command to get a copy of the repository in your local or HPC machine:

```
git clone https://earth.bsc.es/gitlab/ac/Providentia.git
```

When you have finished cloning the repository from Gitlab, you are automatically in the branch `master`. It is recommended to use that branch as it contains the latest features and bug fixes.

## Running the tool the first time

Before running the tool, we suggest checking if conda is installed on your machine. If it is installed, Providentia will create a conda environment on your machine called `providentia-env_v[version]` and install everything that's needed the first time you run it. You can test this by doing:

```
./bin/providentia
```

If the dashboard opens, it worked!

If it didn't and it is because conda is missing, you can follow the steps on [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to install it. If you encountered any other problem, feel free to [contact us](https://earth.bsc.es/gitlab/ac/providentia/-/wikis/home#contact-persons).

By default, a wall time of 2 hours is requested, with 12 CPUs and 30Gb of total memory. This can be modified as desired using the bash options. You can check the available options with:

```
./bin/providentia --usage
``` 

## Accessing the data

If you are running the tool on HPC, you will already see that there are options to choose from in the menu on the top. However, if you opened the dashboard on a local machine, you won't see anything on the dropdowns and you will need to download the data into a directory. By default, the data will be stored in `/home/{user}/data/providentia`, if for some reason you want to save it elsewhere you can edit the path in `settings/data_paths.yaml`.

## Launching the dashboard

As explained, you can launch the dashboard by simply running:

```
./bin/providentia
```

If you already have a [configuration file](Configuration-files), you can specify its path in the command line with the argument `--config`:

```
./bin/providentia --config='/path/to/file/example.conf'
```

If you have defined multiple sections or subsections, a pop-up window will immediately appear where you can choose the section or subsection of interest. After that, the graphical window of Providentia will appear and you can begin using the tool. 

## Generating a report

If you have your [configuration file](Configuration-files) ready, you can generate reports by running Providentia with the following command:

```
./bin/providentia --config='/path/to/file/example.conf' --report
```

The mandatory commands are:

* `--config`: specify the path of your [configuration file](Configuration-files)
* `--report`: for creating a report, without launching the dashboard

You can also launch the dashboard or get a report for only one section by using the option  `--section`. In order to indicate subsections, you will need to write the section name, followed by a hyphen (-) and the subsection name.

The reports will be saved under the folder `reports`. You can add a path in the `report_filename` of the configuration file to change the default directory.

## Using Providentia backend functions

A Jupyter notebook can be launched with the following command:

```
./bin/providentia --notebook
```        

Some examples on how to use Providentia's backend functions can be found in [the notebooks folder](https://earth.bsc.es/gitlab/ac/providentia/-/tree/master/notebooks).

Providentia can also be imported and used in your own Python scripts.

## Interpolating your model data to observations

Using a [configuration file](Configuration-files), you can start interpolating your model data to your desired observational network.
```
./bin/providentia --config='/path/to/file/example.conf' --interpolate
```

More details can be found in the [interpolation section](Interpolation).

## Redirecting log output to a file  

Providentia allows saving its output to a log file using the `--logfile` option:  

```bash
./bin/providentia --logfile <filename>
```

More details [here](https://earth.bsc.es/gitlab/ac/Providentia/-/wikis/Redirecting-output-to-a-file).

Enjoy!