# Download

Now it is possible to download data directly into your machine outside HPC via the download mode. With this, you won't need to spend time waiting for a job to be allocated for you and the tool interactive features will run smoothly as there will be no x11 forwarding.

To download data you just need to add `--download` to your command:

```
./bin/providentia --config='/path/to/file/example.conf' --download
```

This will get the data that needs to be downloaded from your configuration file and save it into the directories specified in `settings/data_paths` for `local`.

If the data is GHOST, it will be grabbed from [Zenodo](https://zenodo.org/records/10637450) and for non-GHOST data, this will come from the `esarchive`.