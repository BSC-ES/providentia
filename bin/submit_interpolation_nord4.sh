#!/bin/bash

#SBATCH --job-name=PRVI                                         
#SBATCH --ntasks=1                                                
#SBATCH --output=logs/interpolation/submission_logs/%J.out
#SBATCH --error=logs/interpolation/submission_logs/%J.err
#SBATCH --account=bsc32
#SBATCH --qos=bsc_es

# load modules from conda environment
module load Miniconda3/24.7.1-0
eval "$(conda shell.bash hook)"
conda config --set env_prompt '($(basename {default_env})) '
if ! { conda config --show-sources | grep '/gpfs/projects/bsc32/repository/apps/conda_envs/'; } >/dev/null 2>&1; then
    conda config --append envs_dirs /gpfs/projects/bsc32/repository/apps/conda_envs/
fi
if { conda env list | grep 'providentia-env_v3.1.0-nord4'; } >/dev/null 2>&1; then 
    echo "INFO: Activating conda environment in /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_v3.1.0-nord4..."
    conda activate providentia-env_v3.1.0-nord4
else 
    echo "ERROR: Environment not found in /gpfs/projects/bsc32/repository/apps/conda_envs/, contact the developers."
fi

JOB_ID="--slurm_job_id=$SLURM_JOB_ID"

# Disable CPU binding to avoid affinity conflicts with the CPUs allocated to the job step when launching interpolation from dashboard
srun --cpu-bind=none --output=logs/interpolation/management_logs/$SLURM_JOB_ID.out python3 -u -c "from providentia.main import main; main()" $@ $JOB_ID