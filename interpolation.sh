#!/bin/bash

#SBATCH --job-name=PRVI                                         
#SBATCH --ntasks=1                                                
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# load modules
source load_modules.sh

srun --output=providentia/interpolation/management_logs/$SLURM_JOB_ID.out python -u providentia/interpolation/experiment_interpolation_submission.py $SLURM_JOB_ID $*
