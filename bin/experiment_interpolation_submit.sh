#!/bin/bash

#SBATCH --job-name=PRVI                                         
#SBATCH --ntasks=1                                                
#SBATCH --output=logs/interpolation/submission_logs/%J.out
#SBATCH --error=/dev/null

JOB_ID="--slurm_job_id=$SLURM_JOB_ID"

srun --output=logs/interpolation/management_logs/$SLURM_JOB_ID.out python -u -c "from providentia.main import main; main()" $@ $JOB_ID