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
                                   #other
                                   'a1wdwodu':{'esarchive':'/esarchive/exp/monarch/a1wdwodu/'},
                                   'a21m':{'esarchive':'/esarchive/exp/monarch/a21m/'},
                                   'a22r':{'esarchive':'/esarchive/exp/monarch/a22r/'},
                                   'a22s':{'esarchive':'/esarchive/exp/monarch/a22s/'},
                                   'a28j':{'esarchive':'/esarchive/exp/monarch/a28j/'},
                                   'a28j_v2':{'esarchive':'/esarchive/exp/monarch/a28j_v2/'},
                                   'a28jwodu':{'esarchive':'/esarchive/exp/monarch/a28jwodu/'},
                                   'a29d':{'esarchive':'/esarchive/exp/monarch/a29d/'},
                                   'a29d_v2':{'esarchive':'/esarchive/exp/monarch/a29d_v2/'},
                                   'a29dwodu':{'esarchive':'/esarchive/exp/monarch/a29dwodu/'},
                                   'a2ba':{'esarchive':'/esarchive/exp/monarch/a2ba/'},
                                   'a2bg':{'esarchive':'/esarchive/exp/monarch/a2bg/'},
                                   'a2bk':{'esarchive':'/esarchive/exp/monarch/a2bk/'},
                                   'a2bz':{'esarchive':'/esarchive/exp/monarch/a2bz/'},
                                   'a22r':{'esarchive':'/esarchive/exp/monarch/a22r/'},
                                   #----------------------------------------------------------------------------#
                                   #CALIOPE OPERATIONAL experiments
                                   'b012':{'esarchive':'/esarchive/exp/wrf-hermes-cmaq/b012/'},
                                   #----------------------------------------------------------------------------#                         
                                   #AUTO-WRF-HERMES-CMAQ experiments
                                   'a1ay':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1ay/'}, 
                                   'a1uo':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1uo/'},
                                   'a148':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a148/'},
                                   'a1wz':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1wz/'},
                                   'a1y6':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1y6/'},
                                   'a1y7':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1y7/'},
                                   'a1y8':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1y8/'},
                                   'a1y9':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1y9/'},
                                   'a1ya':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1ya/'},
                                   'a1zk':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1zk/'},
                                   'a1zm':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a1zm/'},
                                   'a26h':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a26h/'},
                                   'a26j':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a26j/'},
                                   'a26o':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a26o/'},
                                   'a272':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a272/'},  
                                   'a29h':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a29h/'},
                                   'a29i':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a29i/'},
                                   'a2a1':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a2a1/'},
                                   'a2a6':{'esarchive':'/esarchive/exp/auto-wrf-hermes-cmaq/a2a6/'},
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