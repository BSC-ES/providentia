#!/bin/bash
#SBATCH --output notebook.out

# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
login=$(printf $SLURM_SUBMIT_HOST | awk -F 'login' '{print $2}')
if [ "${BSC_MACHINE}" == "nord3v2" ]; then
    bsc_hostname="nord${login}"
elif [ "${BSC_MACHINE}" == "nord4" ]; then
    bsc_hostname="n4login${login}"
elif [ "${BSC_MACHINE}" = "amd" ]; then
    bsc_hostname="amdlogin1"
elif [ "${BSC_MACHINE}" = "mn5" ]; then
    bsc_hostname="glogin${login}"
fi

# print tunneling instructions jupyter-log
echo -e "
Create an SSH tunnel via terminal on your local machine:
ssh -N -L ${port}:${node}:${port} ${user}@${bsc_hostname}.bsc.es
"
export PYTHONPATH=$(pwd):${PYTHONPATH}
jupyter-lab --no-browser --port=${port} --ip=${node} --ServerApp.allow_origin="*"