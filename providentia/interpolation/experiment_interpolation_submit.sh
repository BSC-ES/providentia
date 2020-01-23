#!/bin/bash

#SBATCH --job-name=PRVI                                         
#SBATCH --ntasks=1                                                
#SBATCH --time=48:00:00                                            
#SBATCH --qos=bsc_es
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

#load modules
source load_modules.sh

srun --output=management_logs/$SLURM_JOB_ID.out python -u experiment_interpolation_submission.py $SLURM_JOB_ID
