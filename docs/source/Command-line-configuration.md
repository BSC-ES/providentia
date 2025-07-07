# Command line configuration

There are some parameters that may be configured by command line when executing Providentia, although none are mandatory. 

These parameters can be modified when running the bash script to launch Providentia. You can see all the available options by typing `./bin/providentia --usage` and `./bin/providentia --help`.

Namely, these are:

| Parameter | Description | Default |
| ------ | ------ | ------ |
| cores | Number of cores. | 12 |
| time | Requested time. | 02:00:00 |
| jobname | Job name. | PRV |
| queue | Job queue. | debug |
| mem | Memory. | 20Gb |
| constraint | Memory constraint. | medmem (to use nodes with 64GB of memory) |
| version | Print version of Providentia. | |
| clean | Removes log files. | |
| logfile | Redirects output to a file. More info in the [wiki page](https://earth.bsc.es/gitlab/ac/Providentia/-/wikis/Redirecting-output-to-a-file). | |
| debug | Start [debug mode](https://earth.bsc.es/gitlab/ac/Providentia/-/wikis/Interactive-mode#starting-a-jupyter-notebook). | |
| interactive | Open a Jupyter notebook to [use Providentia as a library](https://earth.bsc.es/gitlab/ac/Providentia/-/wikis/Run-the-tool#running-the-tool-on-a-bsc-machine). | |
| conf | Configuration file path. | |
| config | Configuration file path. | |
| config_dir | Path to all configuration files. | |
| section | Section within configuration file. | |
| ghost_version | GHOST version. | 1.5 |
| cartopy_data_dir | Cartopy data directory. | In MN5: `/gpfs/projects/bsc32/software/rhel/9.2/software/Cartopy/0.23.0-foss-2023b-Python-3.11.5/lib/python3.11/site-packages/cartopy/data`. In other HPC: `/gpfs/projects/bsc32/software/rhel/7.5/ppc64le/POWER9/software/Cartopy/0.17.0-foss-2018b-Python-3.7.0/lib/python3.7/site-packages/Cartopy-0.17.0-py3.7-linux-ppc64le.egg/cartopy/data`. Locally: Downloaded from the internet on the fly. |
| generate_file_tree | Generate file tree to update data directories | |
| file_tree | Generate file tree to update data directories | |
| gft | Generate file tree to update data directories | |
| offline | Start [reports](https://earth.bsc.es/gitlab/ac/Providentia/-/wikis/Run-the-tool#generate-a-report). | |
| network | Network you want to load observations from. Can be multiple (e.g. `CAPMoN, EBAS`). Adding a wild card (\*) is going to expand to certain variables (vconc* → vconc1, vconc2, etc.). | EBAS |
| species | Species to load. Can be multiple (e.g. `sconco3, sconcno2`). | sconco3 |
| resolution | Temporal resolution of the observations you want to load (e.g. `3hourly`). | hourly |
| start_date | Comparison start date in YYYYMMDD format (e.g. `20170101`). | 20180101 |
| end_date | Comparison end date in YYYYMMDD format (e.g. `20180601`). | 20190101 |
| experiments |  ID of interpolated experiment using providentia-interpolation. The experiment IDs can be mapped to different names by adding a list of alternative names after the experiment IDs (e.g. `exp1, exp2 (altexp1, altexp2)`). | |
| temporal_colocation | Boolean variable to set if you want to temporally colocate the observation and experiment data. | False |
| spatial_colocation | Boolean variable to set if you want to spatially colocate the observation and experiment data across multiple species. | True |
| filter_species | Filter read species by other species data within a data range (can be multiple) (e.g. `network1:species1 (lowerlim, upperlim), network2:species2 (lowerlim, upperlim)`). | |
| lower_bound | Data lower limit | |
| upper_bound | Data upper limit | |
| report_type | Type of report to generate that defines which plots the report will contain, from the options given in `report_plots.yaml`. | standard |
| report_summary | Boolean variable to set if you wish to make specific plots for each station in subsection.  | True |
| report_stations | Boolean variable to set if you wish to make summary plots across station subsection. | False |
| report_title | The header in the first page of the report (as in the PDF). | Providentia Report |
| report_filename | The filename of the report or the path to create the report (as in the PDF). | PROVIDENTIA_Report |
| map_extent | Set the map plot extents with the syntax: minimum longitude, maximum longitude, minimum latitude, maximum latitude. | |
| active_dashboard_plots | Plots that will be active in the dashboard once it is launched (e.g. `timeseries, periodic-violin, scatter, distribution`). | timeseries, statsummary, distribution, periodic |
| plot_characteristics_filename | The path to the file containing the plot characteristics. | |

If you want to use the parameter `active_dashboard_plots`, it is important that you leave no spaces between the plot types as in:

```
./bin/providentia --active_dashboard_plots=timeseries,metadata,periodic-violin,boxplot
```

To specify a subsection, add the name of the parent section followed by an interpunct (·) before the subsection name, like in:

```
./bin/providentia --conf=configurations/test.conf --section=SECTIONA·Spain
```