#!/bin/bash

#SBATCH --job-name=PRVI                                         
#SBATCH --ntasks=1                                                
#SBATCH --output=providentia/interpolation/submission_logs/%J.out
#SBATCH --error=/dev/null
#SBATCH --account=bsc32
#SBATCH --qos=gp_bsces

# load modules from conda environment
active_env=false
eval "$(conda shell.bash hook)"
if { conda env list | grep 'providentia-env'; } >/dev/null 2>&1; then 
    echo "Activating environment..."
    conda activate providentia-env
    active_env=true
else
    # create environment if it does not exist if we are in glogin4 (with internet)
    echo "Environment providentia-env does not exist."; 
    node=$(hostname -s)
    if [ "${node}" == "glogin4" ]; then
        echo "Creating environment..."
        conda create -n providentia-env -y python=3.9.16
        conda activate providentia-env
        conda install -c conda-forge cartopy -y
        conda install -c conda-forge jupyterlab -y
        pip install -r requirements.txt
        active_env=true
    else
        echo "Please connect to MN5 using glogin4 to create the environment automatically."
    fi
fi

# load greasy
echo "Loading greasy..."
module load greasy/2.2.4.1

srun --output=providentia/interpolation/management_logs/$SLURM_JOB_ID.out python -u providentia/interpolation/experiment_interpolation_submission.py $SLURM_JOB_ID $*