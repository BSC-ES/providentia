#Providentia interpolation configuration file

###--------------------------------------------------------------------------------------------------###
###SETTING SOME KEY VARIABLES --> DO NOT MODIFY UNLESS YOU KNOW WHAT YOU ARE DOING!
###--------------------------------------------------------------------------------------------------###

#set GHOST version data to work with
GHOST_version = '1.2'

#define the QOS (Quality of Service) used to manage jobs on the SLURM system
#OPTIONS: bsc_es (max walltime of 2 days), prace (max walltime of 3 days), debug (max walltime of 2 hours)
qos = 'debug'

#define the chunk size to process tasks (i.e. the size of CPU chunks used to process tasks)
#this should not exceed the number of CPUs per node on machine (160 virtual cores/40 physical cores on power9, 48 physical cores on marenostrum4)
chunk_size = 20

#define job array limit to process tasks (i.e. the number of CPU chunks allowed to be submitted simultaneously)
job_array_limit = 100

#define if multithreading should be used when processing tasks (i.e. use all threads per physical CPU)
multithreading = False

###--------------------------------------------------------------------------------------------------###
#MODIFY INTERPOLATION VARIABLES HERE
###--------------------------------------------------------------------------------------------------###

#define N nearest neighbours to use for interpolation
n_neighbours_to_find = 4

#define date range to process
start_date = '197001' #YYYYMM   START FROM THIS POINT 
end_date =   '201901' #YYYYMM   GO UP TO THIS POINT    

#define all experiments (by ID) to process
experiments_to_process = ['a1vv','a2jj']

#define complete list of desired species to process (this is later filtered based on availability of species in experiment directories)
species_to_process = ['sconcco','sconchno3','sconcnh3','sconcisop','sconcno','sconcno2','sconco3','sconcpan','sconcso2',
                      'pm10','pm2p5','pm1',
                      'sconcal','sconcas','sconcbc','sconcc','sconcca','sconccd','sconccl','sconccobalt','sconccr','sconccu','sconcec','sconcfe','sconchg','sconck','sconcmg','sconcmn','sconcna','sconcnh4','sconcni','sconcno3','sconcoc','sconcpb','sconcse','sconcso4','sconcv','sconczn',
                      'pm10al','pm10as','pm10bc','pm10c','pm10ca','pm10cd','pm10cl','pm10cobalt','pm10cr','pm10cu','pm10ec','pm10fe','pm10hg','pm10k','pm10mg','pm10mn','pm10na','pm10nh4','pm10ni','pm10no3','pm10oc','pm10pb','pm10se','pm10so4','pm10v','pm10zn',
                      'pm2p5al','pm2p5as','pm2p5bc','pm2p5c','pm2p5ca','pm2p5cd','pm2p5cl','pm2p5cobalt','pm2p5cr','pm2p5cu','pm2p5ec','pm2p5fe','pm2p5hg','pm2p5k','pm2p5mg','pm2p5mn','pm2p5na','pm2p5nh4','pm2p5ni','pm2p5no3','pm2p5oc','pm2p5pb','pm2p5se','pm2p5so4','pm2p5v','pm2p5zn',
                      'pm1al','pm1as','pm1bc','pm1c','pm1ca','pm1cd','pm1cl','pm1cobalt','pm1cr','pm1cu','pm1ec','pm1fe','pm1hg','pm1k','pm1mg','pm1mn','pm1na','pm1nh4','pm1ni','pm1no3','pm1oc','pm1pb','pm1se','pm1so4','pm1v','pm1zn',                            
                      'acprec','dir','dir10','spd','spd10','cldbot','vdist','t','t2','td','td2','slp','acsnow','si','p','pshltr','p10','sst','stc','stc10','stc40','stc100','stc200','ccovmean','cfracmean','rh','rho2',
                      'od500aero','od500aerocorse','od500aerofine','fm500frac','od380aero','od440aero','od550aero','od675aero','od870aero','od1020aero','ae440_870'] 

#define complete list of desired grid_types to process (this is later filtered based on availability of grid types in experiment directories)
grid_types_to_process = ['regional','eu','ip','cat','mad','can']

#define complete list of desired model temporal resolutions to process (this is currently only allowed to be only 'hourly')
model_temporal_resolutions_to_process = ['hourly']

#define ensemble member numbers to process (leave empty to process all members) --> this field is ignored for non-ensemble simulations
model_ensemble_members_to_process = []

#define complete list of GHOST observational networks to interpolate against
#GHOST_networks_to_interpolate_against  = ['EBAS','EEA_AQ_eReporting','NCDC_ISD']
GHOST_networks_to_interpolate_against  = ['EEA_AQ_eReporting']

#define complete lost of interpolated model temporal resolutions to output (interpolating against equivalent temporal resolution observational files)
temporal_resolutions_to_output = ['hourly', 'hourly_instantaneous', 'daily', 'monthly']