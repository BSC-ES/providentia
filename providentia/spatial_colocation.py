""" Functions for the spatial colocation """

import itertools
import numpy as np
import pyproj
from scipy.spatial import cKDTree


def spatial_colocation_nonghost(station_references, longitudes, latitudes):
    """ Given multiple species, return intersecting indices for matching stations across species (per network/species)
        for non-GHOST data.

        This is done by 
            1. Cross-checking the station references between species to get matching station_references
            2. Cross-checking matching longitude / latitude coordinates to a tolerance of 19m difference

        The tolerance is calculated by allowing for a tolerance of 19.053m in the 3 independent x,y,z dimensions, 
        as is done in GHOST to distinguish unique stations.
        Using Pythagoras in 3D √(11**2 +11**2 + 11**2) = 19.053.

        :param station_references: dictionary of station references per network/species
        :type station_references: dict
        :param longitudes: dictionary of longitudes per network/species
        :type longitudes: dict
        :param latitudes: dictionary of latitudes per network/species
        :type latitudes: dict
        :return: intersecting indices per network/species
        :rtype: dict
    """

    # get indices of intersection of station references across species  
    intersecting_station_references = list(set.intersection(*map(set,list(station_references.values()))))
    intersecting_indices = {}
    for networkspecies in station_references:
        intersecting_indices[networkspecies] = np.array([list(station_references[networkspecies]).index(ref) 
                                                        for ref in intersecting_station_references], dtype=np.int64)

    # set variable for first networkspecies
    firstnetworkspecies = list(intersecting_indices.keys())[0]

    # if have zero intersecting indices across species, then return with warning message
    if len(intersecting_indices[firstnetworkspecies]) == 0:
        print('Warning: No intersecting stations across networks/species')
        return intersecting_indices

    # if non-intersecting indices unaccounted for across species, 
    # then attempt to resolve them by matching longitudes / latitudes
    if len(intersecting_indices[firstnetworkspecies]) != len(station_references[firstnetworkspecies]):

        # set tolerance for matching longitudes and latitudes in metres
        tolerance = 19.053

        # get non-intersecting indices, longitudes and latitudes across speci
        non_intersecting_indices = {networkspecies: np.setdiff1d(np.arange(len(station_references[networkspecies])), intersecting_indices[networkspecies]) for networkspecies in station_references}
        non_intersecting_longitudes = {networkspecies: longitudes[networkspecies][non_intersecting_indices[networkspecies]] for networkspecies in longitudes}
        non_intersecting_latitudes = {networkspecies: latitudes[networkspecies][non_intersecting_indices[networkspecies]] for networkspecies in latitudes}

        # get non-intersecting station longitudes and latitudes for first speci
        firstnetworkspecies_longitudes = non_intersecting_longitudes[firstnetworkspecies]
        firstnetworkspecies_latitudes = non_intersecting_latitudes[firstnetworkspecies]

        # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
        # (Earth Centred, Earth Fixed) coordinates assuming WGS84 datum and ellipsoid, and that all heights equal zero
        # ECEF coordinates represent positions (in metres) as X, Y, Z coordinates, approximating the earth surface as an ellipsoid of revolution
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        firstnetworkspecies_x, firstnetworkspecies_y, firstnetworkspecies_z = pyproj.transform(lla, ecef, 
            firstnetworkspecies_longitudes, firstnetworkspecies_latitudes, 
            np.zeros(len(firstnetworkspecies_longitudes)), radians=False)
        
        # merge coordinates to 3D array
        firstnetworkspecies_xyz = np.column_stack((firstnetworkspecies_x, firstnetworkspecies_y, firstnetworkspecies_z))

        # iterate through all other speci, and get intersections (within tolerance) of longitudes and latitudes 
        # with first speci longitudes and latitudes
        pairwise_intersect_inds = {firstnetworkspecies:[]}
        for networkspecies in non_intersecting_longitudes:

            if networkspecies == firstnetworkspecies:
                continue

            nextnetworkspecies_longitudes = non_intersecting_longitudes[networkspecies]
            nextnetworkspecies_latitudes = non_intersecting_latitudes[networkspecies]

            # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
            nextnetworkspecies_x, nextnetworkspecies_y, nextnetworkspecies_z = pyproj.transform(lla, ecef, 
                nextnetworkspecies_longitudes, nextnetworkspecies_latitudes, 
                np.zeros(len(nextnetworkspecies_longitudes)), radians=False)
            
            # merge coordinates to 3D array
            nextnetworkspecies_xyz = np.column_stack((nextnetworkspecies_x, nextnetworkspecies_y, nextnetworkspecies_z))            

            # get all indices of next speci xyz coords, within tolerance of each first speci xyz coords
            idx = cKDTree(nextnetworkspecies_xyz).query_ball_point(firstnetworkspecies_xyz, tolerance)

            # get all indices where have non-duplicated and duplicated matched indices
            unique_idx, unique_idx_counts = np.unique(list(itertools.chain(*idx)), return_counts=True)
            nondup_idx = unique_idx[np.where(unique_idx_counts == 1)[0]]
            dup_idx = unique_idx[np.where(unique_idx_counts > 1)[0]]

            # gather list of indices for first speci and next speci of already resolved within tolerance indices
            fs_wtol_inds = []
            ns_wtol_inds = []

            # resolve all duplicated matched indices
            idx, unresolved_dup_idx, fs_wtol_inds, ns_wtol_inds = resolve_duplicate_spatial_colocation_matches(idx, 
                                                                    nondup_idx, dup_idx, 
                                                                    firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                                    fs_wtol_inds, ns_wtol_inds)

            # pass through again to resolve all unresovered duplicated indices after first pass
            if len(unresolved_dup_idx) > 0:
                idx, _, fs_wtol_inds, ns_wtol_inds = resolve_duplicate_spatial_colocation_matches(idx, 
                                                    nondup_idx, unresolved_dup_idx, 
                                                    firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                    fs_wtol_inds, ns_wtol_inds)

            # iterate though next speci within tolerance indices, per each first speci coord 
            for idx_ii, idx_l in enumerate(idx):
                # if position has already been resolved (through resolving duplicates) then continue
                if idx_ii in fs_wtol_inds:
                    continue
                
                # no matches, then append nothing
                elif len(idx_l) == 0:
                    continue

                # just 1 match, then append
                elif len(idx_l) == 1:
                    fs_wtol_inds.append(idx_ii) 
                    ns_wtol_inds.append(idx_l[0])

                # more than 1 match, then find which match has closest distance, then append
                else:
                    # find the dists between all relevant first speci xyz coords and next speci duplicate xyz coord
                    idx_l = np.array(idx_l)
                    dists = cKDTree([firstnetworkspecies_xyz[idx_ii]]).query(nextnetworkspecies_xyz[idx_l], k=1)[0]
                    ordered_idx_l = idx_l[np.argsort(dists)]
                    fs_wtol_inds.append(idx_ii) 
                    ns_wtol_inds.append(ordered_idx_l[0])

            # order matched indices for first species in ascending order, and order next speci indices in smae way
            if len(fs_wtol_inds) > 0:
                fs_wtol_inds, ns_wtol_inds = list(zip(*sorted(zip(fs_wtol_inds, ns_wtol_inds))))

                # set indices where first species differences are within tolerance, i.e. intersecting 
                pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)] = non_intersecting_indices[firstnetworkspecies][np.array(fs_wtol_inds)]
                pairwise_intersect_inds[firstnetworkspecies].extend(non_intersecting_indices[firstnetworkspecies][np.array(fs_wtol_inds)])
            
                # get indices where next species differences are within tolerance, i.e. intersecting 
                pairwise_intersect_inds[networkspecies] = non_intersecting_indices[networkspecies][np.array(ns_wtol_inds)]
            else:
                pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)] = np.array([], dtype=np.int64)
                pairwise_intersect_inds[networkspecies] = np.array([], dtype=np.int64)

        # get indices (for first networkspecies) where longitude and latitudes intersect across all species
        pairwise_intersect_inds_unique, counts = np.unique(pairwise_intersect_inds[firstnetworkspecies], return_counts=True)
        pairwise_intersect_inds[firstnetworkspecies] = pairwise_intersect_inds_unique[counts == (len(longitudes)-1)]

        if len(pairwise_intersect_inds[firstnetworkspecies]) > 0:
            # get specific intersect indices across all species, for rest of species
            for networkspecies in non_intersecting_longitudes:
                if networkspecies == firstnetworkspecies:
                    continue
                _, species_intersect_inds, _ = np.intersect1d(pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)], pairwise_intersect_inds[firstnetworkspecies], return_indices=True)
                pairwise_intersect_inds[networkspecies] = pairwise_intersect_inds[networkspecies][species_intersect_inds]

            # append newly found intersecting indices to previously found intersect inds
            for networkspecies in non_intersecting_longitudes:
                intersecting_indices[networkspecies] = np.array(np.append(intersecting_indices[networkspecies], pairwise_intersect_inds[networkspecies]), dtype=np.int64)

    return intersecting_indices


def spatial_colocation_ghost(longitudes, latitudes, measurement_altitudes):
    """ Given multiple species, return intersecting indices for matching stations across species (per network/species)
        for GHOST data.

        This is done by cross-checking matching longitude / latitude / measurement altitudes coordinates 
        to a tolerance of 19.053m difference.
        This tolerance is calculated by allowing for a tolerance of 11m in the 3 independent x,y,z dimensions, 
        as is done in GHOST to distinguish unique stations.
        Using Pythagoras in 3D √(11**2 +11**2 + 11**2) = 19.053.

        A current limitation is that at one station there can be several measurement methods, which
        are represented as unique stations in GHOST. Currently, if these stations have the same measurement 
        position, simply the first of these stations will be preferentially chosen as a match.
        This could be better done by prioritising first by method, when have multiple matches.

        :param longitudes: dictionary of longitudes per network/species
        :type longitudes: dict
        :param latitudes: dictionary of latitudes per network/species
        :type latitudes: dict
        :param measurement_altitudes: dictionary of measurement altitudes per network/species
        :type measurement_altitudes: dict
        :return: intersecting indices per network/species
        :rtype: dict
    """

    # set tolerance for matching longitudes / latitudes / measurement_altitudes in metres
    tolerance = 19.053

    # set variable for first networkspecies
    firstnetworkspecies = list(longitudes.keys())[0]
    
    # get station coordinates for firstnetworkspecies
    firstnetworkspecies_longitudes = longitudes[firstnetworkspecies]
    firstnetworkspecies_latitudes = latitudes[firstnetworkspecies] 
    firstnetworkspecies_measurement_altitudes = measurement_altitudes[firstnetworkspecies]

    # convert longitudes / latitudes / measurement_altitudes in geographic coordinates to cartesian ECEF 
    # (Earth Centred, Earth Fixed) coordinates assuming WGS84 datum and ellipsoid, and that all heights equal zero
    # ECEF coordinates represent positions (in metres) as X, Y, Z coordinates, approximating the earth surface as an ellipsoid of revolution
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    firstnetworkspecies_x, firstnetworkspecies_y, firstnetworkspecies_z = pyproj.transform(lla, ecef, 
        firstnetworkspecies_longitudes, firstnetworkspecies_latitudes, firstnetworkspecies_measurement_altitudes,
        radians=False)

    # merge coordinates to 3D array
    firstnetworkspecies_xyz = np.column_stack((firstnetworkspecies_x, firstnetworkspecies_y, firstnetworkspecies_z))

    # iterate through all other speci, and get intersections (within tolerance) of 
    # longitudes / latitudes / measurement_altitudes, with first speci longitudes and latitudes
    pairwise_intersect_inds = {firstnetworkspecies:[]}
    for networkspecies in longitudes:

        if networkspecies == firstnetworkspecies:
            continue

        nextnetworkspecies_longitudes = longitudes[networkspecies]
        nextnetworkspecies_latitudes = latitudes[networkspecies]
        nextnetworkspecies_measurement_altitudes = measurement_altitudes[networkspecies]

        # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
        nextnetworkspecies_x, nextnetworkspecies_y, nextnetworkspecies_z = pyproj.transform(lla, 
            ecef, nextnetworkspecies_longitudes, nextnetworkspecies_latitudes, nextnetworkspecies_measurement_altitudes,
            radians=False)
        
        # merge coordinates to 3D array
        nextnetworkspecies_xyz = np.column_stack((nextnetworkspecies_x, nextnetworkspecies_y, nextnetworkspecies_z))            

        # get all indices of next speci xyz coords, within tolerance of each first speci xyz coords
        idx = cKDTree(nextnetworkspecies_xyz).query_ball_point(firstnetworkspecies_xyz, tolerance)

        # get all indices where have non-duplicated and duplicated matched indices
        unique_idx, unique_idx_counts = np.unique(list(itertools.chain(*idx)), return_counts=True)
        nondup_idx = unique_idx[np.where(unique_idx_counts == 1)[0]]
        dup_idx = unique_idx[np.where(unique_idx_counts > 1)[0]]

        # gather list of indices for first speci and next speci of already resolved within tolerance indices
        fs_wtol_inds = []
        ns_wtol_inds = []

        # resolve all duplicated matched indices
        idx, unresolved_dup_idx, fs_wtol_inds, ns_wtol_inds = resolve_duplicate_spatial_colocation_matches(idx, 
                                                                nondup_idx, dup_idx, 
                                                                firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                                fs_wtol_inds, ns_wtol_inds)

        # pass through again to resolve all unresovered duplicated indices after first pass
        if len(unresolved_dup_idx) > 0:
            idx, _, fs_wtol_inds, ns_wtol_inds = resolve_duplicate_spatial_colocation_matches(idx, 
                                                   nondup_idx, unresolved_dup_idx, 
                                                   firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                   fs_wtol_inds, ns_wtol_inds)

        # iterate though next speci within tolerance indices, per each first speci coord 
        for idx_ii, idx_l in enumerate(idx):
            # if position has already been resolved (through resolving duplicates) then continue
            if idx_ii in fs_wtol_inds:
                continue
            
            # no matches, then append nothing
            elif len(idx_l) == 0:
                continue

            # just 1 match, then append
            elif len(idx_l) == 1:
                fs_wtol_inds.append(idx_ii) 
                ns_wtol_inds.append(idx_l[0])

            # more than 1 match, then find which match has closest distance, then append
            else:
                # find the dists between all relevant first speci xyz coords and next speci duplicate xyz coord
                idx_l = np.array(idx_l)
                dists = cKDTree([firstnetworkspecies_xyz[idx_ii]]).query(nextnetworkspecies_xyz[idx_l], k=1)[0]
                ordered_idx_l = idx_l[np.argsort(dists)]
                fs_wtol_inds.append(idx_ii) 
                ns_wtol_inds.append(ordered_idx_l[0])

        # order matched indices for first species in ascending order, and order next speci indices in smae way
        if len(fs_wtol_inds) > 0:
            fs_wtol_inds, ns_wtol_inds = list(zip(*sorted(zip(fs_wtol_inds, ns_wtol_inds))))

            # set indices where first species differences are within tolerance, i.e. intersecting 
            pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)] = np.array(fs_wtol_inds)
            pairwise_intersect_inds[firstnetworkspecies].extend(np.array(fs_wtol_inds))
        
            # get indices where next species differences are within tolerance, i.e. intersecting 
            pairwise_intersect_inds[networkspecies] = np.array(ns_wtol_inds)
        else:
            pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)] = np.array([], dtype=np.int64)
            pairwise_intersect_inds[networkspecies] = np.array([], dtype=np.int64)

    # get indices (for first networkspecies) where longitude, latitudes and measurement_altitudes intersect across all species
    pairwise_intersect_inds_unique, counts = np.unique(pairwise_intersect_inds[firstnetworkspecies], return_counts=True)
    pairwise_intersect_inds[firstnetworkspecies] = pairwise_intersect_inds_unique[counts == (len(longitudes)-1)]

    # get specific intersect indices across all species, for rest of species
    intersecting_indices = {}
    for networkspecies in longitudes:
        if networkspecies == firstnetworkspecies:
            intersecting_indices[networkspecies] = np.array(pairwise_intersect_inds[networkspecies], dtype=np.int64)
        else:
            _, species_intersect_inds, _ = np.intersect1d(pairwise_intersect_inds['{}_{}'.format(firstnetworkspecies, networkspecies)], pairwise_intersect_inds[firstnetworkspecies], return_indices=True)
            intersecting_indices[networkspecies] = np.array(pairwise_intersect_inds[networkspecies][species_intersect_inds], dtype=np.int)

    return intersecting_indices


def resolve_duplicate_spatial_colocation_matches(idx, nondup_idx, dup_idx, 
                                                 firstnetworkspecies_xyz, nextnetworkspecies_xyz,
                                                 fs_wtol_inds, ns_wtol_inds):

    """ Function that resolves duplicate indices found during spatial colocation of 2 species.

        In spatial colocation it is neccessary to match each stations geographically within a specific
        tolerance, for 2 different species. 

        In some cases the stations within the tolerance can match for multiple stations. 
        In order to resolve this, for each duplicated station, it is set to match with the closest station to it
        in terms of 3D distance. This is done iteratively until there are no more duplicates.

        :param idx: per first networkspeci xyz coords, a list of indices of next networkspeci xyz coords within tolerance
        :type idx: array 
        :param nondup_idx: next networkspeci idx coords that are not duplicated across idx array 
        :type nondup_idx: array
        :param dup_idx: next networkspeci idx coords that are duplicated across idx array 
        :type dup_idx: array
        :param firstnetworkspecies_xyz: ECEF coordinates for first networkspeci stations
        :type firstnetworkspecies_xyz: array
        :param nextnetworkspecies_xyz: ECEF coordinates for next networkspeci stations
        :type nextnetworkspecies_xyz: array
        :param fs_wtol_inds: first networkspeci station indices within tolerance (i.e. have paired match)
        :type fs_wtol_inds: list
        :param ns_wtol_inds: next networkspeci station indices within tolerance (i.e. have paired match)
        :type ns_wtol_inds: list
        :return: idx, unresolved_dup_idx, fs_wtol_inds, ns_wtol_inds 
        :rtype: array, list, list, list
    """

    # resolve all duplicated matched indices by finding for which index the distance is closest
    unresolved_dup_idx = []
    for dup_index in dup_idx:

        # get all relevant first speci indices for which contain the duplicate next speci index
        relevant_idx = np.array([idx_ii for idx_ii, idx_l in enumerate(idx) if dup_index in idx_l])
        # find the dists between all relevant first speci xyz coords and next speci duplicate xyz coord
        dists = cKDTree([nextnetworkspecies_xyz[dup_index]]).query(firstnetworkspecies_xyz[relevant_idx], k=1)[0]
        # order the relevant first speci indices by dists (closest first)
        ordered_relevant_idx = relevant_idx[np.argsort(dists)]
        # iterate through ordered relevant first speci indices until find an ind for which can claim duplicate ind
        # keep going in order iteratively until have exhausted all options for duplicate ind
        for ordered_relevant_index_ii, ordered_relevant_index in enumerate(ordered_relevant_idx):
        
            # if current relevant first speci index has no competing indices, then append matched index,
            # remove duplicate index from other match lists, and then break out of iteration
            if len(idx[ordered_relevant_index]) == 1:
                fs_wtol_inds.append(ordered_relevant_index)
                ns_wtol_inds.append(dup_index)
                for next_ordered_relevant_index in ordered_relevant_idx[ordered_relevant_index_ii+1:]:
                    idx[next_ordered_relevant_index].remove(dup_index)
                break

            # otherwise, 
            # find the dists between the relevant first speci xyz coord and next speci duplicate xyz coords
            else:
                idx_l = np.array(idx[ordered_relevant_index])
                dists = cKDTree([firstnetworkspecies_xyz[ordered_relevant_index]]).query(nextnetworkspecies_xyz[idx_l], k=1)[0]
                ordered_idx_l = idx_l[np.argsort(dists)]

                # if the closest index is the dup index, or is another non-duplicate, then append it
                # if the closest index is the dup index, then remove duplicate index from other match lists 
                # and then break out of iteration
                if (ordered_idx_l[0] == dup_index) or (ordered_idx_l[0] in nondup_idx):
                    fs_wtol_inds.append(ordered_relevant_index)
                    ns_wtol_inds.append(ordered_idx_l[0])
                    if ordered_idx_l[0] == dup_index:
                        for next_ordered_relevant_index in ordered_relevant_idx[ordered_relevant_index_ii+1:]:
                            idx[next_ordered_relevant_index].remove(dup_index)
                        break

                # else, if the closest index is another duplicated index,
                # then cannot currently resolve position, so come back to this later
                else:
                    unresolved_dup_idx.append(dup_index)
                    break

    return idx, unresolved_dup_idx, fs_wtol_inds, ns_wtol_inds 
