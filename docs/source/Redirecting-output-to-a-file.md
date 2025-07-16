# Redirecting Output to a File

Providentia provides the option to save its output in a log file. 

## Supported Modes

The logging feature is available in the following modes:
- **Dashboard**
- **Report**
- **Notebook**
- **Download**

It is **not available** in the **Interpolation** mode.

## Usage

To enable logging, use the `--logfile` argument when running Providentia from the command line in **Dashboard**, **Report**, and **Download** modes. 

To enable logging in a notebook, set the `logfile` argument when calling `Interactive`, as explained below. In this case, each Providentia object will be linked to a single log file.

### Default Logging
```bash
./bin/providentia --logfile
```
In the notebook mode, the command looks like this:

```python
Interactive(conf='debug.conf',logfile=True)
```

This command will create a log file inside the `logs` folder within the directory corresponding to the active mode. For example, logs generated in dashboard mode will be saved in `logs/dashboard`. The log files are named with a timestamp in the format `%Y%m%d%H%M%S.log`, such as `20250313123045.log`.

### Custom Log Filename
If you want to specify a custom filename, you can provide it as an argument:
```bash
./bin/providentia --logfile=custom_filename
```
In the notebook mode, it is like this:

```python
Interactive(conf='debug.conf',logfile='custom_filename')
```
This will save the log file with the name `custom_filename` inside the default `logs` folder.

### Custom Log File Path
You can also define a custom file path, absolute or relative:
```bash
./bin/providentia --logfile=/custom/path/custom_filename
```
In the notebook mode, it is like this:

```python
Interactive(conf='debug.conf',logfile='/custom/path/custom_filename')
```
This will save the log file in `/custom/path/` with the name `custom_filename`. Ensure that the specified path exists, Providentia will not create non-existent directories.
