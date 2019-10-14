#!/bin/bash

#SBATCH --job-name=PRVI_S                                         
#SBATCH --ntasks=1                                                
#SBATCH −−cpus−per−task=1
#SBATCH --time=48:00:00                                            
#SBATCH --qos=debug
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

#load modules
source load_modules.sh

srun --cpu-bind=core --output=management_logs/$SLURM_JOB_ID.out --error=management_logs/$SLURM_JOB_ID.err python -u experiment_interpolation_submission.py $SLURM_JOB_ID
