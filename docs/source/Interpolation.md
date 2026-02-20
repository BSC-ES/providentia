# Interpolation

This mode allows users to spatially interpolate model output against available observational stations, allowing it to be subsequently evaluated in **Providentia**. 

Interpolation consists of spatially interpolating gridded model outputs to observational station locations using a nearest-neighbour approach.

The Interpolation mode allows to interpolate models against **GHOST** and **non-GHOST** observations.

## Getting started

To start an interpolation, you need to add either `--interp`, `--interpolate` or `--interpolation` as a launch option along with the **mandatory** configuration file on the command line:

```
./bin/providentia --config='/path/to/file/example.conf' --interp 
```   
```
./bin/providentia --config='/path/to/file/example.conf' --interpolate 
```         
```
./bin/providentia--config='/path/to/file/example.conf' --interpolation 
```     

The interpolation mode fetches all the content specified in your configuration file across all sections. To only run one specific section, add the `--section` parameter to the command.

In terms of performance, we recommend running Providentia Interpolation in MareNostrum5.

## Interpolation configuration fields

All parameters that can be used in the interpolation configuration files can be found in the [Shared Parameters](shared-parameters) or [Interpolation Parameters](interpolation-parameters) sections of the Configuration Fields page.

## Execution details

### Local users

For local execution, the interpolation runs in the background using multiprocessing.
Be aware that this can be demanding on the computer, so ensure your machine has sufficient resources before running.


### BSC HPC users

Upon submission, a first job named `PRV` will start the submission process which will make a job called `PRVI` appear in the SLURM queue, shortly afterwards a job array named `PRVI_$SLURMJOBID` (containing the jobs for all the defined variable combinations) will be submitted.

When all jobs have been completed (or there is a failure) the `PRVI` job will exit the queue.


### Interpolation considerations

Providentia is highly flexible when handling interpolation between model and observational data, for example in mapping species, adjusting for different temporal resolutions and using wildcards to select multiple values at once.

#### Mapping species

When checking if an model is stored in a location with the corresponding domain, resolution, and species, consider that the species might not always be listed under the same name.

The file `settings/internal/mapping_species.yaml` contains a dictionary mapping original species names to their alternative names. 

Note that the mapping species file is only used when the species name from the configuration file is not found in the expected location, meaning Povidentia first looks for the species written in the configuration file. If it is not found, it then searches for the corresponding mapped species in `mapping_species.yaml`.

#### Different temporal resolutions between observations and model

When you have observational and model data with different temporal resolutions, Providentia is very adaptable to try and ensure that an interpolation takes place.

For each temporal resolution you are wishing to interpolate to, Providentia will go through a series of steps:
1. It will first check to see if you both observations and model data at that resolution. If there are no observations at the resolution, the interpolation will not be performed.
2. If there are observations but no model data, Providentia will next check if there is model data at a finer resolution available. If there is, it will then downsample the model data to the coarser resolution of the observations. 
3. If there is no finer model data available, it will next check if there is model data at a coarser resolution available. If there is, it will then upsample the model data to the finer resolution of the observations.
4. If there is no finer or coarser model data available, the interpolation will not be performed.

The downsampling or upsampling of the model data that Providentia performs can be controlled via a few variables.

The statistic for the downsampling of model data to a coarser observational resolution can be set via the `interp_model_downsampling` variable. The valid options are: `mean` and `median`, with the default being `mean`.

```
interp_model_downsampling = mean
```

The method for the upsampling of model data to a finer observational resolution can be set via the `interp_model_upsampling` variable. The valid options are: `fill` and `gaps`, with the default being `fill`. `fill` linearly fills between measurements, and `gaps` sets NaN values for times that the model does not have.

```
interp_model_upsampling = fill
```

#### Using wildcards

You can use the `*` wildcard in the following fields to automatically select all available values:

- `network`, `observation`, `framework` 
- `model`, `models`, `experiments`, `experiment`  
- `species`  
- `resolution`  
- `start_date`
- `end_date`  

**Note:** Using wildcards may result in large numbers of interpolations, so use with caution.

## Logs

Every time an interpolation is done, logs are saved in the `logs/interpolation` folder.

To check the status/output of an interpolation job, the following log files are created on submission in different directories:

* ### Management logs

   These logs provide an overview of the interpolation process. Most errors will appear here.
   
   Located in the `logs/interpolation/management_logs` folder, look for an `$SLURMJOBID.out` file.

* ### Submission logs

   These logs contain information about the Slurm and Greasy submissions to the HPC machines.
   
   Found in the `logs/interpolation/submission_logs` folder, search for an `$SLURMJOBID.out` file.

* ### Interpolation logs

   These logs give information about individual interpolations and how long it took to do them.
   
   Found in the `logs/interpolation/interpolation_logs` folder, for each individual interpolation, new directories are created with the structure `{model}/{species}/{network}/{resolution}`. Inside these directories, logs for each month are stored as `{YYYYMM}_{exit_code}.out`. If successful, the exit code will be 0.

## Input data

### Observational data

Observational data is read from the directories defined in `settings/data_paths.yaml`, with `ghost_root` for GHOST observations and `nonghost_root` for non-GHOST observations.

If no network can be located under `ghost_root` or `nonghost_root`, the interpolation will fail during submission.

### Model data

Providentia locates model data differently depending on whether it is run locally or on BSC HPC systems.

For local and HPC executions, model data is firstly located using the paths defined in `settings/data_paths.yaml`.  
By default, models are expected to be found under the `mod_to_interp_root` directory defined in that file.

If no model can be located under `mod_to_interp_root`, the interpolation will fail during submission.

### Model data (HPC-specific)

On BSC HPC systems, some models are stored in fixed locations that cannot be easily moved.  
For this reason, in addition to `settings/data_paths.yaml`, HPC users can define model locations in `settings/interp_models.yaml`.

When running an interpolation on HPC, Providentia searches for model data in the following order:

1. `mod_to_interp_root` defined in `data_paths.yaml`
2. Paths defined in `interp_models.yaml`
3. If the model is not found, the interpolation fails during submission

(define-models)=
#### Defining models in `interp_models.yaml`

The `settings/interp_models.yaml` file contains a dictionary of default relevant models grouped by type, which contains the list of model names and their possible storage paths.

If a model is located in one of the predefined paths (for example `/esarchive/exp/monarch/`), it only needs to be added to the corresponding model list.

If the model is stored elsewhere, you must define:
- The **model type**
- The **model name**
- The **model storage directory**, excluding the model name

```
"example_model_type": {
        "models": ["example_model_name"],
        "paths": [ 
            "/example/model/path"
        ]
}
```

You can find this exact template at the end of the `interp_models.yaml` file.

##### Model directory structure

When adding a new model to a directory, if you want it to be read from Providentia, the subdirectories inside the model storage directory must follow this structure: `{model_name}/{domain}/{resolution}/{species}`. For example: `cams61_monarch_ph3/eu/hourly/sconco3`.

There can be multiple paths to the same model, and you can add them to the list of paths. The order is important: the first path that works on the machine will be used.

There's normally two location types of model data:

* **gpfs**: Accessible by the ***MareNostrum5/Nord4*** machines.
* **esarchive**: Accessible by the ***Nord4*** machine.

If you are using a machine that allows both types of paths, it is recommended to list your `gpfs` paths first. This is because when reading data from the `esarchive`, a major limitation on the read time is the transfer speed between the 2 machines, reading directly from the `gpfs`  directory circumvents this therefore.

## Output data

Interpolated model data is written to the directory defined by `mod_root` in `settings/data_paths.yaml`.  
The default value of this path depends on the execution environment.

- **BSC HPC users**:  
  `/gpfs/projects/bsc32/AC_cache/recon/exp_interp/`

- **Local users**:  
  `~/data/providentia/mod`

This can be changed by updating `mod_root` or editing `settings/data_paths.yaml`.