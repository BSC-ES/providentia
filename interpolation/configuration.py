#define date range to process
start_date = '197001' #YYYYMM   START FROM THIS POINT 
end_date =   '201901' #YYYYMM   GO UP TO THIS POINT    

#define all experiments (by ID) to process
experiments_to_process = ['a1vv']

#define complete list of desired species to process (this is later filtered based on availability of species in experiment directories)
species_to_process = ['sconcco','sconchno3','sconcnh3','sconcisop','sconcno','sconcno2','sconco3','sconcpan','sconcso2',
                      'pm10','pm2p5','pm1',
                      'sconcal','sconcas','sconcbc','sconcc','sconcca','sconccd','sconccl','sconccobalt','sconccr','sconccu','sconcec','sconcfe','sconchg','sconck','sconcmg','sconcmn','sconcna','sconcnh4','sconcni','sconcno3','sconcoc','sconcpb','sconcse','sconcso4','sconcv','sconczn',
                      'pm10al','pm10as','pm10bc','pm10c','pm10ca','pm10cd','pm10cl','pm10cobalt','pm10cr','pm10cu','pm10ec','pm10fe','pm10hg','pm10k','pm10mg','pm10mn','pm10na','pm10nh4','pm10ni','pm10no3','pm10oc','pm10pb','pm10se','pm10so4','pm10v','pm10zn',
                      'pm2p5al','pm2p5as','pm2p5bc','pm2p5c','pm2p5ca','pm2p5cd','pm2p5cl','pm2p5cobalt','pm2p5cr','pm2p5cu','pm2p5ec','pm2p5fe','pm2p5hg','pm2p5k','pm2p5mg','pm2p5mn','pm2p5na','pm2p5nh4','pm2p5ni','pm2p5no3','pm2p5oc','pm2p5pb','pm2p5se','pm2p5so4','pm2p5v','pm2p5zn',
                      'pm1al','pm1as','pm1bc','pm1c','pm1ca','pm1cd','pm1cl','pm1cobalt','pm1cr','pm1cu','pm1ec','pm1fe','pm1hg','pm1k','pm1mg','pm1mn','pm1na','pm1nh4','pm1ni','pm1no3','pm1oc','pm1pb','pm1se','pm1so4','pm1v','pm1zn']                            

#define complete list of desired grid_types to process (this is later filtered based on availability of grid types in experiment directories)
grid_types_to_process = ['regional','eu','ip','cat','mad','can']

#define complete list of desired model temporal resolutions to process (this is currently only allowed to be only 'hourly')
model_temporal_resolutions_to_process = ['hourly']

#define complete list of GHOST observational networks to interpolate against
#GHOST_networks_to_interpolate_against  = ['EBAS','EIONET']
GHOST_networks_to_interpolate_against  = ['EIONET']

#define complete lost of interpolated model temporal resolutions to output (interpolating against equivalent temporal resolution observational files)
temporal_resolutions_to_output = ['hourly', 'daily', 'monthly']

#define N nearest neighbours to use for interpolation
n_neighbours_to_find = 4

#define qos (quality of service): bsc_es, debug
qos = 'debug'