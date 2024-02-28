#!/bin/bash
#SBATCH --output interactive.out

# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
# print tunneling instructions jupyter-log
echo -e "
Create an SSH tunnel via terminal on your local machine:
ssh -N -L ${port}:${node}:${port} ${user}@nord4.bsc.es
"
export PYTHONPATH=$(pwd):${PYTHONPATH}
# DON'T USE ADDRESS BELOW.
# DO USE TOKEN BELOW
jupyter-lab --no-browser --port=${port} --ip=${node} --ServerApp.allow_origin="*" --ServerApp.ip="0.0.0.0"