#!/bin/bash

#SBATCH --job-name=SubmitG                                         
#SBATCH --ntasks=1                                                
#SBATCH −−cpus−per−task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00                                            

srun --cpu-bind=core --output=/esarchive/recon/ghost_interp/process_job_arrays/management_logs/$SLURM_JOB_ID.out --error=/esarchive/recon/ghost_interp/process_job_arrays/management_logs/$SLURM_JOB_ID.err python -u /esarchive/recon/ghost_interp/process_job_arrays/GHOST_experiment_interpolation_submission.py $SLURM_JOB_ID $PWD
