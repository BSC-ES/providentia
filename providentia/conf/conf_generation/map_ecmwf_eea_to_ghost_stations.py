import numpy as np
from netCDF4 import Dataset
import glob
import pandas as pd
import copy

species_map = {'CO':'sconcco','NO2':'sconcno2','O3':'sconco3','SO2':'sconcso2','PM10':'pm10','PM25':'pm2p5'}

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

def create_config_file(fname, eval_type, network, resolution, matrix, species, start_date, end_date, mapped_stations):
  
    mapped_stations_str = ', '.join([str(elem) for elem in sorted(list(mapped_stations.values()))])

    config_file = open('generated_confs/{}'.format(fname), "w")
    config_file.write('[{}]\n'.format(eval_type))
    config_file.write('selected_network = {}\n'.format(network))
    config_file.write('selected_resolution = {}\n'.format(resolution))
    config_file.write('selected_matrix = {}\n'.format(matrix))
    config_file.write('selected_species = {}\n'.format(species))
    config_file.write('start_date = {}\n'.format(start_date))
    config_file.write('end_date = {}\n'.format(end_date))
    config_file.write('station_reference = keep:{};'.format(mapped_stations_str))
    config_file.close()

stations_files = sorted(glob.glob('/esarchive/obs/ecmwf/eea-cams50/original_files/stations/*'))

for stations_file in stations_files:

    print(stations_file)
    if 'assimilation' in stations_file:
        eval_type = 'assimilation'
    elif 'validation' in stations_file:
        eval_type = 'validation'

    stations_file_read = pd.read_csv(stations_file, sep=" ", header=None)
    stations = stations_file_read[0].tolist()
    stations_unfound = copy.deepcopy(stations)
    mapped_stations = {}

    species = species_map[stations_file.split('.')[1]]
    if species == 'pm10':
        matrix = 'pm10'
    elif species == 'pm2p5':
        matrix = 'pm2.5'
    else:
        matrix = 'gas'

    conf_fname = 'ecmwf_cams2_40_eea_{}_background_{}_set10.conf'.format(species,eval_type)

    ghost_files = sorted(glob.glob('/gpfs/projects/bsc32/AC_cache/obs/ghost/EEA_AQ_eReporting/1.4/hourly/{}/{}*'.format(species,species)))

    for ghost_file in ghost_files:
        root = Dataset(ghost_file)
        ghost_station_references = root['station_reference'][:]
        ghost_associated_networks =  root['associated_networks'][:]

        stations_unfound, mapped_stations = cross_check(stations_unfound, ghost_station_references, ghost_associated_networks, mapped_stations)

        #replace turkey codes
        stations_unfound = ['TR{}'.format(stn[1:]) if (stn[0] == 'T') & (stn[1].isdigit()) else stn for stn in stations_unfound]
        stations_unfound, mapped_stations = cross_check(stations_unfound, ghost_station_references, ghost_associated_networks, mapped_stations)

        root.close()

    print(len(stations_unfound), stations_unfound)
    print()

    create_config_file(conf_fname, eval_type, 'EEA_AQ_eReporting', 'hourly', matrix, species, '20180101', '20190101', mapped_stations)

