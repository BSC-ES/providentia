#!/bin/bash
#SBATCH --output interactive.out

# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
login=$(printf $HOSTNAME | awk -F 'login' '{print $2}')
if [ "${BSC_MACHINE}" == "nord3v2" ]; then
    bsc_hostname="nord$login"
elif [ "${BSC_MACHINE}" == "mn4" ]; then
    bsc_hostname="mn$login"
elif [ "${BSC_MACHINE}" = "amd" ]; then
    bsc_hostname="amdlogin1"
elif [ "${BSC_MACHINE}" = "power" ]; then
    bsc_hostname="plogin$login"
elif [ "${ip}" == "84.88.185.48" ]; then
    bsc_hostname="$HOSTNAME"
fi

# print tunneling instructions jupyter-log
echo -e "
Create an SSH tunnel via terminal on your local machine:
ssh -N -L ${port}:${node}:${port} ${user}@${bsc_hostname}.bsc.es
"
export PYTHONPATH=$(pwd):${PYTHONPATH}
# DON'T USE ADDRESS BELOW.
# DO USE TOKEN BELOW
jupyter-lab --no-browser --port=${port} --ip=${node}
