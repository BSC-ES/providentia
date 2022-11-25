import numpy as np
from netCDF4 import Dataset
import glob
import pandas as pd
import copy

# script to generate Providentia .conf files for CAMS Interim Reanalysis (IRA),
# and CAMS Validated Reanalysis (VRA)

# please indicate the type of station files to use (IRA or VRA)  
process_type = 'IRA'
# please indicate the year of station files to process (YYYY)
station_year = '2022'
# please indicate the GHOST version
ghost_version = '1.4'
# please indicate the species to process (as a list)
species = ['O3', 'NO2', 'CO', 'SO2', 'PM10', 'PM25']
# please indicate the GHOST network to use
network = 'EEA_AQ_eReporting'
# please indicate the temporal resolution of data to use (hourly, daily or monthly)
resolution = 'hourly'
# please indicate the start_date of data to read (YYYYMMDD)
start_date = '20200101'
# please indicate the end_date of data to read (YYYYMMDD)
end_date = '20210101'

def cross_check(stations_unfound, ghost_stations, associated_networks, mapped_stations):
   
    found_stations = []

    for stn in stations_unfound:

        #modified ghost references
        for ghost_stn in ghost_stations:
            if stn in ghost_stn:
                found_stations.append(stn)
                mapped_stations[stn] = ghost_stn
                break

        if stn not in found_stations:

            #associated networks
            for associated_network_ii, associated_network in enumerate(associated_networks):
                if stn in associated_network:
                    found_stations.append(stn) 
                    mapped_stations[stn] = ghost_stations[associated_network_ii]
                    break

    for stn in found_stations:
        stations_unfound.remove(stn)

    return stations_unfound, mapped_stations

def create_config_file(fname, section, network, resolution, species, start_date, end_date, mapped_stations):
  
    mapped_stations_str = ', '.join([str(elem) for elem in sorted(list(mapped_stations.values()))])

    config_file = open('generated_confs/{}'.format(fname), "w")
    config_file.write('[{}]\n'.format(section))
    config_file.write('network = {}\n'.format(network))
    config_file.write('species = {}\n'.format(species))
    config_file.write('resolution = {}\n'.format(resolution))
    config_file.write('start_date = {}\n'.format(start_date))
    config_file.write('end_date = {}\n'.format(end_date))
    config_file.write('station_reference = keep:{} ||'.format(mapped_stations_str))
    config_file.close()

species_map = {'CO':'sconcco','NO2':'sconcno2','O3':'sconco3','SO2':'sconcso2','PM10':'pm10','PM25':'pm2p5'}

if process_type == 'IRA':
    station_files = sorted(glob.glob('/esarchive/obs/ineris/eionet-cams2_40-ira/original_files/stations_{}/stations*'.format(station_year)))
elif process_type == 'VRA':
    station_files = sorted(glob.glob('/esarchive/obs/ineris/eionet-cams2_40-vra/original_files/stations_{}/stations*'.format(station_year)))

# iterate through species
for speci in species:
    station_files_cut = [station_file for station_file in station_files if '.{}.'.format(speci) in station_file]

    speci_mapped = species_map[speci] 

    for station_file in station_files_cut:

        if 'assimilation' in station_file:
            eval_type = 'assimilation'
        elif 'validation' in station_file:
            eval_type = 'validation'

        station_file_read = pd.read_csv(station_file, sep=" ", header=None)
        stations = station_file_read[0].tolist()
        stations_unfound = copy.deepcopy(stations)
        mapped_stations = {}

        # get set number of station file
        setno = station_file.split('_')[-1]

        # create section name
        section = '{}-{}'.format(process_type,eval_type)

        # set fname of conf file
        conf_fname = 'ecmwf_cams2_40_{}_{}_{}_background_{}_{}_{}.conf'.format(speci,process_type,eval_type,setno, start_date, end_date)

        # get relevant ghost files
        ghost_files = sorted(glob.glob('/gpfs/projects/bsc32/AC_cache/obs/ghost/{}/{}/{}/{}/{}*'.format(network,ghost_version,resolution,speci_mapped,speci_mapped)))

        for ghost_file in ghost_files:
            root = Dataset(ghost_file)
            ghost_station_references = root['station_reference'][:]
            ghost_associated_networks =  root['associated_networks'][:]

            stations_unfound, mapped_stations = cross_check(stations_unfound, ghost_station_references, ghost_associated_networks, mapped_stations)

            #replace turkey codes
            stations_unfound = ['TR{}'.format(stn[1:]) if (stn[0] == 'T') & (stn[1].isdigit()) else stn for stn in stations_unfound]
            stations_unfound, mapped_stations = cross_check(stations_unfound, ghost_station_references, ghost_associated_networks, mapped_stations)

            root.close()

        create_config_file(conf_fname, section, network, resolution, speci_mapped, start_date, end_date, mapped_stations)
