#define all experiment IDs, with associated file directories on 'esarchive' or optionally on 'gpfs'
defined_experiments_dictionary = {#----------------------------------------------------------------------------#
                                   #MONARCH experiments
                                   #b007 - previous Europe region benchmark
                                   'b007':{'esarchive':'/esarchive/exp/nmmb-bsc-ctm/b007/'},
                                   #a1wd - current Europe region benchmark 
                                   'a1wd':{'esarchive':'/esarchive/exp/monarch/a1wd/'},
                                   #a1vv - radiation/photolysis solved every 5 minutes
                                   'a1vv':{'esarchive':'/esarchive/exp/monarch/a1vv/'},
                                   #a1vw - horizontal diffusion smag = 0.4 
                                   'a1vw':{'esarchive':'/esarchive/exp/monarch/a1vw/'},
                                   #a1x8 - MEGAN NO soil emissions
                                   'a1x8':{'esarchive':'/esarchive/exp/monarch/a1x8/'},
                                   #a1xa - 0.1x0.1 degree horizontal resolution
                                   'a1xa':{'esarchive':'/esarchive/exp/monarch/a1xa/'},
				                   #a1xf - CAMS-81 gridded temporal profiles
                                   'a1xf':{'esarchive':'/esarchive/exp/monarch/a1xf/'},
                                   #a21m - CAMS-81 v2.1 gridded temporal profiles	
                                   'a21m':{'esarchive':'/esarchive/exp/monarch/a21m/'},
                                   #a1yq - new SOA Simple scheme, speciation primary emissions HERMES  
				                   'a1yq':{'esarchive':'/esarchive/exp/monarch/a1yq/'},	
                                   #a1yz - effect of dry deposition on aerosols VDRY=1
                                   'a1yz':{'esarchive':'/esarchive/exp/monarch/a1yz/'},
                                   #a1zc - CO boundary conditions in volume mixing ratio
                                   'a1zc':{'esarchive':'/esarchive/exp/monarch/a1zc/'},		
                                   #a28j
                                   'a28j':{'esarchive':'/esarchive/exp/monarch/a28j/'},	
                                   #----------------------------------------------------------------------------#
                                   #CALIOPE OPERATIONAL experiments
                                   'b012':{'esarchive':'/esarchive/exp/wrf-hermes-cmaq/b012/'},
                                   #----------------------------------------------------------------------------#                         
                                   #AUTO-WRF-HERMES-CMAQ experiments
                                   'a148':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a148/'},
                                   #----------------------------------------------------------------------------#
                                   #CAMS 50 experiments
                                   #CAMS 50 ensemble median
                                   'cams50':{        'esarchive':'/esarchive/recon/ecmwf/cams50/'},
                                   #CAMS 50 CHIMERE 
                                   'cams50_chimere':{'esarchive':'/esarchive/recon/ecmwf/cams50_chimere/'},
                                   #CAMS 50 EMEP 
                                   'cams50_emep':{   'esarchive':'/esarchive/recon/ecmwf/cams50_emep/'},
                                   #CAMS 50 MATCH 
                                   'cams50_match':{  'esarchive':'/esarchive/recon/ecmwf/cams50_match/'},
                                   #CAMS 50 MOCAGE 
                                   'cams50_mocage':{ 'esarchive':'/esarchive/recon/ecmwf/cams50_mocage/'},
                                   #CAMS 50 SILAM 
                                   'cams50_silam':{  'esarchive':'/esarchive/recon/ecmwf/cams50_silam/'}
                                  }
