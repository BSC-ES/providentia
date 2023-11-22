#!/bin/bash
#SBATCH --ntasks 12
#SBATCH --qos debug
#SBATCH --time 02:00:00
#SBATCH --job-name PROVI
#SBATCH --output log_jupyter-notebook-%J.out
#SBATCH --error log_jupyter-notebook-%J.err
# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
# print tunneling instructions jupyter-log
echo -e "
MacOS or linux terminal command to create your ssh tunnel
ssh -N -L ${port}:${node}:${port} ${user}@nord4.bsc.es
Use a Browser on your local machine to go to:
localhost:${port}  (prefix w/ https:// if using password)
"
# load modules or conda environments here
source load_modules.sh
module load jupyterlab/3.0.9-foss-2019b-Python-3.7.4
export PYTHONPATH=/esarchive/obs/ghost/scripts/Providentia-interactive:${PYTHONPATH}
# DON'T USE ADDRESS BELOW.
# DO USE TOKEN BELOW
jupyter-lab --no-browser --port=${port} --ip=${node}
