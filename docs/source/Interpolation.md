# Interpolation

Since the release of version 2.4.0 in June 2024, our old tool Providentia Interpolation has been merged into Providentia.

These are some details that you need to be aware of:
- The old configuration files can no longer be used, instead you will be able to use the ones from Providentia to interpolate and visualise your data.
- The file `defined_experiments.py` has been removed and there is no need to define the full paths to the experiments anymore. However, the experiment ids need to be passed into the `/settings/experiment_names.yaml`.
- The code that runs the interpolation is found under the folder `interpolation`. The logs will be created there as usual.

Now that you know these points, you can run the interpolation by doing:
```
./bin/providentia --config='/path/to/file/example.conf' --interpolate
```

This will generate logs in different folders:

- Submission logs: To give information on the slurm and Greasy submissions to the HPC machines.
- Management logs: To inform about the interpolation overall. Most errors will appear here.
- Interpolation logs: To get information about individual interpolations and know how long it took to do them.