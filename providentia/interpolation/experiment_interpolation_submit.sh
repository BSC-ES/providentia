#!/bin/bash

#SBATCH --job-name=PRVI                                         
#SBATCH --ntasks=1                                                
#SBATCH --output=providentia/interpolation/submission_logs/%J.out
#SBATCH --error=/dev/null

echo "Loading modules..."

# load modules
source load_modules.sh

srun --output=providentia/interpolation/management_logs/$SLURM_JOB_ID.out python -u providentia/interpolation/experiment_interpolation_submission.py $SLURM_JOB_ID $*

