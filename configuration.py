#GHOST Interactive Configuration File 
###------------------------------------------------------------------------------------###

###------------------------------------------------------------------------------------###
###IMPORT BUILT-IN MODULES
###------------------------------------------------------------------------------------###

import os
import re
import socket
import subprocess
import sys

###------------------------------------------------------------------------------------###
###SETTING SOME NECESSARY DETAILS --> DO NOT MODIFY!
###------------------------------------------------------------------------------------###

#get hostname
hostname = socket.getfqdn()
#get available N CPUs
if 'power.cte' in hostname:
    bash_command = 'squeue -h -o "%C"'
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    available_CPUs = int(re.findall(r'\d+', str(output))[0])
else:
    available_CPUs = int(os.cpu_count())

#set cartopy data directory (needed on CTE-POWER as has no external internet connection)
if 'power.cte' in hostname:
    cartopy_data_dir = '/esarchive/scratch/Earth/fbeninca/python/data'
#on all machines except CTE-POWER, pull from internet
else:
    cartopy_data_dir = ''

###------------------------------------------------------------------------------------###
###DEFINE NUMBER OF CPUs TO UTILISE
###------------------------------------------------------------------------------------###

#Define number of CPUs to process on (leave empty to automatically utilise all available CPUs)
#NOTE: if this value is set higher than the actual number of CPUs available, then the max number of CPUs is used. 
n_CPUs = ''
if (n_CPUs == '') or (int(n_CPUs) > available_CPUs):
    n_CPUs = available_CPUs

###------------------------------------------------------------------------------------###
###DEFINE OBSERVATIONAL/EXPERIMENT ROOT DATA DIRECTORIES
###------------------------------------------------------------------------------------###

#Define observational root data directory (if undefined it is automatically taken from the BSC machine the tool is ran on)
obs_root = ''
#set observational root data directory if left undefined
if obs_root == '':
    #running on workstation?
    if 'bscearth' in hostname:
        obs_root = '/esarchive/obs/ghost'
    #running on CTE-POWER?
    elif 'power.cte' in hostname:
        obs_root = '/gpfs/projects/bsc32/AC_cache/obs/ghost'
    #can not recognise machine? --> exit with message
    else:
        sys.exit('GHOST Interactive cannot be executed as observational root directory is empty, and the machine the tool is being run on is not recognised.')
    
#Define experiment root data directory
exp_root = ''
#set experiment root data directory if left undefined
if exp_root == '':
    #running on workstation?
    if 'bscearth' in hostname:
        exp_root = '/esarchive/recon/ghost_interp'
    #running on CTE-POWER?
    elif 'power.cte' in hostname:
        exp_root = '/gpfs/projects/bsc32/AC_cache/recon/ghost_interp'
    #can not recognise machine? --> exit with message
    else:
        sys.exit('GHOST Interactive cannot be executed as experiment root directory is empty, and the machine the tool is being run on is not recognised.')
    
###------------------------------------------------------------------------------------###
###DEFINE COLOURMAPS (see all options here: https://matplotlib.org/examples/color/colormaps_reference.html)
###------------------------------------------------------------------------------------###

#set the sequential colourmap (used to view absolute values)
#one of the perceptually uniform sequential colourmaps is recommended here: viridis, plasma, inferno, magma
sequential_colourmap = 'viridis'

#set the warm sequential colourmap (used to view absolute biases going from zero bias to a maximum bias)
sequential_colourmap_warm = 'Reds'

#set the diverging colourmap (used to evaluate differences)
diverging_colourmap = 'bwr'
