# Download
Now it is possible to download data directly into your machine outside HPC via the download mode. With this, you won't need to spend time waiting for a job to be allocated for you and the tool interactive features will run smoothly as there will be no x11 forwarding.

To download data you just need to add `--download` to your command:

```
./bin/providentia --config='/path/to/file/example.conf' --download
```

This will get the data that needs to be downloaded from your configuration file and save it into the directories specified in `settings/data_paths` for `local`.

The download mode fetches all the content specified in your configuration file across all sections, similar to how it functions in offline mode.

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