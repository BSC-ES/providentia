# Providentia interpolation configuration file

# define the QOS (Quality of Service) used to manage jobs on the SLURM system
# OPTIONS: bsc_es (max walltime of 2 days), class_a (max walltime of 3 days), class_b (max walltime of 3 days),
#          prace (max walltime of 3 days), debug (max walltime of 2 hours)
#  'default' = 'bsc_es'
qos = 'bsc_es'

# For MN5
# qos="gp_bsces"


# set GHOST version data to work with
#  'default' = latest GHOST version
GHOST_version = 'default'

# define N nearest neighbours to use for interpolation
#  'default' = 4
n_neighbours_to_find = 'default'

# define date range to process
start_date = '202309'  # YYYYMM   START FROM THIS POINT
end_date   = '202311'  # YYYYMM   GO UP TO THIS POINT

# define all experiments (by ID) to process
experiments_to_process = ['a6gn']

# define complete list of desired species to process (filtered based on availability of species in experiment 
# directories), e.g. ['sconco3','sconcno2']
# for binned parameters, all bins can be processed by replacing the bin number with a '*', e.g. 'vconcaerobin*' 
#  'default' = all species defined in relevant version of GHOST_standards.py:
species_to_process = ['sconco3']

# define complete list of desired grid_types to process (this is later filtered
# based on availability of grid types in experiment directories)
#  'default': global, regional, eu, ip, cat, mad, can, bcn
grid_types_to_process = ['regional']

# define options for interpolation of ensemble simulations, this can be a specific
# ensemble number, or a specific ensemble statistic
#  **IGNORED FOR NON-ENSEMBLE SIMULATIONS --> leave as 'default'
#  --> for ensemble members, simply give the member number, or for all members input 'allmembers'
#  --> for ensemble statistics simply define them as 'stat_' + the stat name provided in the filename,
#  e.g. 'stat_av' for ensemble average, 'stat_av_an' for ensemble analysis average, or 'stat_av_inc' 
#       for ensemble analysis - first guess average
#  'default' = all ensemble members 
ensemble_options = ['default']

# define complete list of observational networks to interpolate against
#  --> in case of GHOST, write providers name (e.g. EBAS, EEA_AQ_eReporting, etc.).
#      Available GHOST networks are given here: https://earth.bsc.es/gitlab/ac/GHOST 
#  --> in case of esarchive non-GHOST, write the two following folders after /esarchive/obs/ that correspond to the 
#      provider you need. (e.g. nilu/ebas, eea/eionet, nasa-aeronet/oneill_v3-lev15)
#  'default' = all processed GHOST networks in relevant version  
networks_to_interpolate_against = ['eea/eionet']

# define complete list of interpolated model temporal resolutions to output (interpolating against equivalent temporal 
# resolution observational files)
# in the case of instantaneous resolutions, if unavailable from the model, a temporal average resolution is attempted 
# to be taken instead, e.g. for 'hourly_instantaneous', 'hourly' is taken instead
# in the case of non-instantaneous resolutions, if unavailable from the model, an instantaneous resolution is attempted 
# to be taken instead, e.g. for 'hourly', 'hourly_instantaneous' is taken instead
#  'default' = hourly, hourly_instantaneous, daily, monthly
temporal_resolutions_to_output = ['hourly']
