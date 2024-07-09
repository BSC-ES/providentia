#!/bin/bash

#SBATCH --job-name=PRVI                                         
#SBATCH --ntasks=1                                                
#SBATCH --output=providentia/interpolation/submission_logs/%J.out
#SBATCH --error=/dev/null
#SBATCH --account=bsc32
#SBATCH --qos=gp_bsces

# load modules from conda environment
module load anaconda
eval "$(conda shell.bash hook)"
conda config --set env_prompt '($(basename {default_env})) '
if ! { conda config --show-sources | grep '/gpfs/projects/bsc32/repository/apps/conda_envs/'; } >/dev/null 2>&1; then
    conda config --append envs_dirs /gpfs/projects/bsc32/repository/apps/conda_envs/
fi
if { conda env list | grep 'providentia-env_2.4.0'; } >/dev/null 2>&1; then 
    echo "Activating conda environment in /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_2.4.0..."
    conda activate providentia-env_2.4.0
else 
    echo "Environment not found in /gpfs/projects/bsc32/repository/apps/conda_envs/"
fi

# load greasy
echo "Loading dependencies..."
source load_modules.sh

srun --output=providentia/interpolation/management_logs/$SLURM_JOB_ID.out python -u providentia/interpolation/experiment_interpolation_submission.py $SLURM_JOB_ID $*