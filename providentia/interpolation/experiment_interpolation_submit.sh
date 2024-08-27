#!/bin/bash

#SBATCH --job-name=PRVI                                         
#SBATCH --ntasks=1                                                
#SBATCH --output=providentia/interpolation/submission_logs/%J.out
#SBATCH --error=/dev/null

JOB_ID="--slurm_job_id=$SLURM_JOB_ID"

srun --output=providentia/interpolation/management_logs/$SLURM_JOB_ID.out python -u providentia/interpolation/experiment_interpolation_submission.py $@ $JOB_ID
