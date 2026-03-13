# Command line configuration

Command-line configuration allows users to quickly override default or configuration file values at runtime without modifying the configuration files.

All parameters that are allowed in the configuration file can also be passed via the command line. For a full list, please visit the [Configuration fields page](Configuration-fields).

However, there are some parameters that can only be set through the command line. These are:

| Parameter | Description | Default |
| ------ | ------ | ------ |
| `report`, `reports`, `offline` | Start [reports](Report). | — |
| `download`, `dl` | Start [downloads](Download). | — |
| `interpolation`, `interpolate`, `interp` | Start [interpolations](Interpolation). | — |
| `notebook`, `nb`, `jupyter` | Open a [Jupyter notebook](Notebooks) to use Providentia as a [library](Library). | — |
| `debug` | Start [debug mode](Running-the-tool-on-debug). | — |
| `clean` | Removes log files. | — |
| `conf`, `config` | Configuration file path or configuration file name if the file is stored in `providentia/configurations`. | — |
| `section` | Section within configuration file. | — |
| `logfile` | Redirects output to a file. More info in the [Redirecting output to a file](Redirecting-output-to-a-file)  page. | — |
| `cores` | Number of cores. | `12` |
| `time` | Requested time. | `02:00:00` |
| `jobname` | Job name. | `PRV` |
| `queue` | Job queue. | `debug` |
| `mem` | Memory. | `20Gb` |
| `constraint` | Memory constraint. | `medmem` (to use nodes with 64GB of memory) |
| `version`, `V` | Print version of Providentia. | — |
| `generate_file_tree`, `gft` | Generate file tree to update data directories. | `False` |
| `disable_file_tree`, `dft` | Disable file tree to update data directories. | `False` |
| `cpus_per_task` | Number of CPUs per task. | `12` |

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