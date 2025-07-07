# Download
Now it is possible to download data directly into your machine outside HPC via the download mode. With this, you won't need to spend time waiting for a job to be allocated for you and the tool interactive features will run smoothly as there will be no x11 forwarding.

To download data you just need to add `--download` or `--dl` to your command:

```
./bin/providentia --config='/path/to/file/example.conf' --download
```
OR

```
./bin/providentia --config='/path/to/file/example.conf' --dl
```

This will get the data that needs to be downloaded from your configuration file and save it into the directories specified in `settings/data_paths.yaml` for `local`.

The download mode fetches all the content specified in your configuration file across all sections, similar to how it functions in the reports.


## Types of downloads
There are currently four types of downloads available:

1. **Download of Network from HPC Machines:**
   - Downloads ghost and non-ghost networks from the storage5 machine (or mn5 if storage5 is down).
   - Saves them within the `ghost_root` and `nonghost_root` paths specified in `settings/data_paths.yaml`.
   - **How to get this type of download:** Answer `y` to the _Do you want to download from BSC remote machine?_ question or change the `BSC_DL_CHOICE` value to `y`.

2. **Download of Network from Zenodo:**
   - Downloads networks from the [GHOST Zenodo webpage](https://zenodo.org/records/10637450).
   - No need for a BSC HPC machine account, as the download is from the internet.
   - To download a network, use the network name exactly as it appears in the zip file. For example, use `EBAS-COLOSSAL_tursk` instead of `EBAS`.
   - **How to get this type of download:** Answer `n` to the _Do you want to download from BSC remote machine?_ question or change the `BSC_DL_CHOICE` value to `n`.

3. **Download of Interpolated Experiments:**
   - Downloads interpolated experiments from the storage5 machine (or mn5 if storage5 is down).
   - Saves them in the `exp_root` paths specified in `settings/data_paths.yaml`.
   - **How to get this type of download:** Set the `interpolated` field in the configuration file to `True` or don't even set it since that is the default.

4. **Download of Non-interpolated Experiments:**
   - For local downloads:
       - Downloads non-interpolated experiments from the storage5 machine (or mn5 if storage5 is down).
       - The preferred path is specified in `interp_experiments.yaml`. If the experiment is not listed or no path or experiment domain is found, the system will look in the `exp_to_interp_root` path specified in `settings/data_paths.yaml`.
   - For HPC downloads:
       - In this case, the download (more accurately, the copy) will transfer non-interpolated experiments from esarchive to the gpfs folder found in `exp_to_interp_root` in the `settings/data_paths.yaml` file.
       - It only copies from the paths that are specified in `interp_experiments.yaml`. 
   - **How to get this type of download:** Set the `interpolated` field in the configuration file to `False`.

*To know how to define an experiment in `interp_experiments.yaml`, please visit [this page](https://earth.bsc.es/gitlab/ac/Providentia/-/wikis/Interpolation#define-experiments).

## Where to download from

As for today download mode is available:
- **In local** for all 4 modes:
    1. Download of Network from HPC Machines.
    2. Download of Network from Zenodo.
    3. Download of Interpolated Experiments.
    4. Download of Non-interpolated Experiment (download from esarchive or gpfs to local).
- **In storage5/nord4** for:
    1. Download of Non-interpolated Experiment (copy from esarchive to gpfs).

## Download Requirements Summary

The table below summarizes the basic requirements to download data from BSC:

| Download Type                  | Available from Local | Available from HPC Machines | Need to have a BSC account to an HPC machine |
|--------------------------------|----------------------|-----------------------------|---------------------------------------------|
| Network from HPC Machines      | ✓                    | ✗                           | ✓                                           |
| Network from Zenodo            | ✓                    | ✗                           | ✗                                           |
| Interpolated Experiment        | ✓                    | ✗                           | ✓                                           |
| Non-interpolated Experiment    | ✓                    | ✓                           | ✓                                           |

## Fields

Fields work a bit different in this mode, as only the parameters below will be used, some of them without being mandatory:

| Parameter       | Mandatory in Download Mode | Mandatory in Other Modes |
| --------------- | -------------------------- | ------------------------ |
| network         | ✓                           | ✓                        |
| species         | ✗                           | ✓                        |
| resolution      | ✗                           | ✓                       |
| start_date      | ✓                           | ✓                        |
| end_date        | ✓                           | ✓                        |
| experiments     | ✗                           | ✗                        |
| filter_species  | ✗                           | ✗                        |

Other important changes from the other modes are:

1. **Wildcard (*) keyword:** If you want to download all available networks or experiments, simply use the wildcard `*` in the `networks` or `experiments` field. Be careful with this.
2. **Default keyword:** The keyword `default` functions the same as leaving the field empty, meaning it will retrieve all possible values for that parameter.

## .env file

An `.env` file will appear in the Providentia root directory when using the download mode. It is designed to store specific user preferences.

   - **BSC_DL_CHOICE:** (values: y/n) This setting controls where GHOST network data will be downloaded from. If set to `n`, the data will be downloaded from the [Zenodo](https://zenodo.org/records/10637450) webpage. If set to `y`, it will download from the remote machine `transfer1.bsc.es` by default.
   - **PRV_USER:** This setting specifies the username used to connect to the remote machines. It can be any valid username, e.g.:bsc000000.
   - **PRV_PWD:** This setting allows you to save the password needed for connecting to remote machines. If you have properly generated a SSH key, there will be no need to input a password.
   - **OVERWRITE:** (values: y/n) This setting determines whether to automatically overwrite (`y`) or not overwrite (`n`) any existing files in your local directories.

These values can be changed directly on the `.env` file and also be updated by Providentia during the next run.