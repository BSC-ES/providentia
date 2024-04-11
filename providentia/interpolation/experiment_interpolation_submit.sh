#!/bin/bash

#SBATCH --job-name=PRVI                                         
#SBATCH --ntasks=1                                                
#SBATCH --output=submission_logs/%J.out
#SBATCH --error=/dev/null

echo "Loading modules..."

# load modules
source load_modules.sh

srun --output=management_logs/$SLURM_JOB_ID.out python -u experiment_interpolation_submission.py $SLURM_JOB_ID
