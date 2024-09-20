# Interpolation

Starting from the release of version 2.4.0 in August 2024, our old tool Providentia Interpolation will be merged into Providentia and it will be now refered as the Interpolation mode.

This mode allows users to spatially interpolate experiment output against available observational stations to be viewable in **Providentia**. 

The Interpolation mode allows to interpolate experiments against **GHOST** and **non-GHOST** observations.

Contained in: `/gpfs/projects/bsc32/AC_cache/recon/exp_interp/` exists an archive of previously interpolated experiment output, depending on the version. If the experiment you want to analyse is not there, or if you have re-run the experiment, add new variables etc. then this guide will aid you in interpolating the experiment output.

## Starting an Interpolation

The interpolation has been designed to only run on the **MareNostrum5/Nord3v2** machines, so first the user must be logged onto either machine.

To start an interpolation,  you need to add either --interp, --interpolate, or --interpolation as a launch option along with the **mandatory** configuration file in the command line. This will initiate the interpolation process.

```
./bin/providentia --interp --config='/path/to/file/example.conf'
```   
```
./bin/providentia --interpolate  --config='/path/to/file/example.conf'
```         
```
./bin/providentia --interpolation --config='/path/to/file/example.conf'
```         
Upon submission, a first job named `PRV` will start the submission process which will make a job called `PRVI` appear in the SLURM queue, shortly afterwards a job array named `PRVI_$SLURMJOBID` (containing the jobs for all the defined variable combinations) will be submitted.

When all jobs have been completed (or there is a failure) the `PRVI` job will exit the queue.

In terms of performance, we recommend running Providentia Interpolation in MareNostrum5.

## Checking output  

The code that runs the interpolation is found under the `providentia/interpolation` folder. All the logs will be created there too.

To check the status/output of an interpolation job, the following log files are created on submission in different directories:

* ### Management logs

   These logs provide an overview of the interpolation process. Most errors will appear here.
   
   Located in the `providentia/interpolation/management_logs` folder, look for an `$SLURMJOBID.out` file.

* ### Submission logs

   These logs contain information about the Slurm and Greasy submissions to the HPC machines.
   
   Found in the `providentia/interpolation/submission_logs` folder, search for an `$SLURMJOBID.out` file.

* ### Interpolation logs

   These logs give information about individual interpolations and how long it took to do them.
   
    Found in the `providentia/interpolation/interpolation_logs` folder, for each individual interpolation, new directories are created with the structure `{experiment}/{species}/{network}/{resolution}`. Inside these directories, logs for each month are stored as `{YYYYMM}_{exit_code}.out`. If successful, the exit code will be 0.

 
If the interpolation is succesful, the resulted interpolated files are stored under the directory : `/gpfs/projects/bsc32/AC_cache/recon/exp_interp/` followed by the latest version of GHOST used by the interpolation.

## Define experiments

All the experiments that are runned in the interpolation need to be defined in `/settings/interp_experiments.yaml`.

This file contains a dictionary of default relevant experiments grouped by type, in where you find the list of names of the experiments and the list of possible paths for them. 

If you have a new experiment to interpolate and it is located on one of the examples paths such as `/esarchive/exp/monarch/`, you only need to add it to its list of experiments.

If the experiment is not located in one of the predefined paths, you will need to add the **experiment type**, the **experiment name**, and the **experiment storage directory** (excluding the experiment name). For example:

```
"example_experiment_type": {
        "experiments": ["example_experiment_name"],
        "paths": [ 
            "/example/experiment/path"
        ]
```

You can find this exact template at the end of the interp_experiments.yaml file.

When adding a new experiment, the subdirectories inside the experiment storage directory must follow this structure: `{experiment_name}/{domain}/{resolution}/{species}`. For example: `cams61_monarch_ph3/eu/hourly/sconco3`.

There can be multiple paths to the same experiments, and you can add them to the list of paths. The order is important: the first path that works on the machine will be used.

There's normally two types the location of experiment data:

* **gpfs**: Accessible by the ***MareNostrum5/Nord3v2*** machines.
* **esarchive**: Accessible by the ***Nord3v2*** machine.

If you are using a machine that allows both types of paths, it is recommended to list your `gpfs` paths first. This is because when reading data from the `esarchive`, a major limitation on the read time is the transfer speed between the 2 machines, reading directly from the `gpfs`  directory circumvents this therefore.