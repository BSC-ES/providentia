# Command line configuration

Command-line configuration allows users to quickly override default or configuration file values at runtime without modifying the configuration files.

All parameters that are allowed in the configuration file can also be passed via the command line. For a full list, please visit the [Configuration fields page](Configuration-fields).

However, there are some parameters that can only be set through the command line. These are:

| Parameter | Description | Default |
| ------ | ------ | ------ |
| report, reports, offline | Start [reports](Report). | — |
| download, dl | Start [downloads](Download). | — |
| interpolation, interp, interpolate | Start [interpolations](Interpolation). | — |
| notebook, nb, jupyter | Open a [Jupyter notebook](Notebooks) to use Providentia as a [library](Library). | — |
| clean | Removes log files. | — |
| debug | Start [debug mode](Running-the-tool-on-debug). | — |
| cores | Number of cores. | 12 |
| time | Requested time. | 02:00:00 |
| jobname | Job name. | PRV |
| queue | Job queue. | debug |
| mem | Memory. | 20Gb |
| constraint | Memory constraint. | medmem (to use nodes with 64GB of memory) |
| version, V | Print version of Providentia. | — |
| logfile | Redirects output to a file. More info in the [Redirecting output to a file page](Redirecting-output-to-a-file). | — |
| conf, config | Configuration file path. | — |
| config_dir | Path to all configuration files. | — |
| section | Section within configuration file. | — |
| cartopy_data_dir | Cartopy data directory. | In MN5: `/gpfs/projects/bsc32/software/rhel/9.2/software/Cartopy/0.23.0-foss-2023b-Python-3.11.5/lib/python3.11/site-packages/cartopy/data`. In other HPC: `/gpfs/projects/bsc32/software/rhel/7.5/ppc64le/POWER9/software/Cartopy/0.17.0-foss-2018b-Python-3.7.0/lib/python3.7/site-packages/Cartopy-0.17.0-py3.7-linux-ppc64le.egg/cartopy/data`. Locally: Downloaded from the internet on the fly. |
| generate_file_tree, gft | Generate file tree to update data directories | — |
| disable_file_tree, dft | Disable file tree to update data directories | — |

You can also see all the available options by typing `./bin/providentia --usage` and `./bin/providentia --help`.

## Use cases

If a parameter is set both in the configuration file and on the command line, the command-line value takes precedence and will override the configuration file. For example:

```
./bin/providentia --network=EBAS
```

The final value used would be `EBAS`.

To provide a parameter with multiple values, do not include spaces between them. For example:

```
./bin/providentia --active_dashboard_plots=timeseries,metadata,periodic-violin,boxplot
```

To specify a subsection, add the name of the parent section followed by an interpunct (·) before the subsection name, like in:

```
./bin/providentia --conf=configurations/test.conf --section=SECTIONA·Spain
```