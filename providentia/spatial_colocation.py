""" Class for spatial colocation """

import copy
import itertools
from packaging.version import Version
import sys

import numpy as np
import pyproj
import scipy
from scipy.spatial import cKDTree

class SpatialColocation:

    """ Given multiple species, return intersecting indices for matching stations across species (per network/species)
        for GHOST/non-GHOST data.

        This is done by: 
            1. Cross-checking the station references between species to get matching station references
            2. Cross-checking the station names between species to get matching station names
            3. Cross-checking the position coordinates (longitude / latitude / measurement altitude (if available)) 
               to a given tolerance

        The user has the option to configure by which variables they perform the spatial colocation, but by default all 
        ways are active. Each check is performed in order while non-intersecting stations exist across species.

        If any station references / names per speci are duplicated, then do not count them as intersecting.
        This can happen for example when at one station, several measurement methods are used for measurement.
        These stations would eventually be matched in the position check, with simply the first available station 
        being taken as the match. This is a current limitation, and could be better done by prioritising first by 
        method, when have multiple matches.

        Because station names are not unique and in theory matching names can be located at different places in the 
        world, validation is used to ensure the stations are genuine matches. This is done by checking the 
        longitude / latitude of the stations are within a certain validation tolerance.

        The default position check tolerance is calculated by allowing for a tolerance of 11m in the 3 independent 
        x,y,z dimensions, as is done in GHOST to distinguish unique stations. 
        Using Pythagoras in 3D √(11**2 +11**2 + 11**2) = 19.053m
    """

    def __init__(self, read_instance):

        self.read_instance = read_instance

        # get list of all networkspecies 
        self.networkspecies = list(self.read_instance.station_references.keys())

        # set variable for first networkspeci
        self.firstnetworkspeci = self.networkspecies[0]

        # do spatial colocation
        self.intersecting_indices = {networkspeci:[] for networkspeci in self.networkspecies}

        # intialise non-intersection indices
        self.get_non_intersections()

        # check if have measurement_altitudes or not
        if not self.read_instance.station_measurement_altitudes:
            if self.read_instance.spatial_colocation_measurement_altitude:
                print("Warning: spatial_colocation is not using measurement_altitude as no valid values were found for this variable.")
            self.read_instance.spatial_colocation_measurement_altitude = False 

        # check if have station_names or not
        if not self.read_instance.station_names:
            if self.read_instance.spatial_colocation_station_name:
                print("Warning: spatial_colocation is not using station_name as no valid values were found for this variable.")
            self.read_instance.spatial_colocation_station_name = False 

        # if wanting to use measurement altitudes for spatial colocation, but not lons/lats, then return as it is not possible
        if (self.read_instance.spatial_colocation_measurement_altitude) & (not self.read_instance.spatial_colocation_longitude_latitude):
            print("Warning: spatial_colocation is set to False, as spatial_colocation_longitude_latitude must be set to True if spatial_colocation_measurement_altitude is True.")
            return

        # if not wanting to use any variables for spatial colcoation, then return as there is nothing to be done
        elif (not self.read_instance.spatial_colocation_station_reference) & (not self.read_instance.spatial_colocation_station_name) & (not self.read_instance.spatial_colocation_longitude_latitude) & (not self.read_instance.spatial_colocation_measurement_altitude):
            print("Warning: spatial_colocation is set to False, as have no active variables to perform colocation.")
            return

        # define coordinate systems for transformation
        if Version(pyproj.__version__) >= Version("2.0.0"):
            self.lla = {"proj": "latlong", "ellps": "WGS84", "datum": "WGS84"}
            self.ecef = {"proj": "geocent", "ellps": "WGS84", "datum": "WGS84"}
            self.transformer = pyproj.Transformer.from_crs(self.lla, self.ecef)
        else:
            self.lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
            self.ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')

        # by station reference?
        if self.read_instance.spatial_colocation_station_reference:
            self.by_station_reference()

        # by station name (if still have indices unaccounted for)?
        if self.read_instance.spatial_colocation_station_name:
            self.by_station_name()

        # by measurement position: longitude, latitude and option for measurement_altitude (if still have indices unaccounted for)?
        if self.read_instance.spatial_colocation_longitude_latitude:
            self.by_position()

    def by_station_reference(self):

        # for GHOST data, remove measurement method str and duplicate station number from station references
        if self.read_instance.reading_ghost:
            station_references_cut = {k: np.array(['_'.join(ref.split('_')[:-1]) for ref in v]) for k,v in self.read_instance.station_references.items()}
        else:
            station_references_cut = copy.deepcopy(self.read_instance.station_references)

        # get intersecting station references
        intersecting_station_references = list(set.intersection(*map(set,list(station_references_cut.values()))))

        # only continue wih check if have intersecting station references
        if len(intersecting_station_references) > 0:

            # if any station references per speci are duplicated, then do not count them as intersecting
            refs_to_remove = []
            for networkspeci in station_references_cut:
                unique_refs, unique_counts = np.unique(station_references_cut[networkspeci], return_counts=True)
                refs_to_remove.extend(unique_refs[np.where(unique_counts > 1)[0]])
            refs_to_remove = np.unique(refs_to_remove)
            intersecting_station_references = np.array([ref for ref in intersecting_station_references if ref not in refs_to_remove])

            # get indices of intersecting station references per speci
            for networkspeci in station_references_cut:

                station_references_sorted = np.argsort(station_references_cut[networkspeci])
                intersect_refs_pos = np.searchsorted(station_references_cut[networkspeci][station_references_sorted], intersecting_station_references)
                self.intersecting_indices[networkspeci] = np.array(np.append(self.intersecting_indices[networkspeci], 
                                                                             station_references_sorted[intersect_refs_pos]),
                                                                             dtype=np.int64)

            # validate intersections by longitude/latitude
            self.validate_intersections()

            # sort intersecting indices
            self.sort_intersecting_indices()

            # get non-intersection indices
            self.get_non_intersections()

    def by_station_name(self):

        # only do check if have remaining non-intersecting indices
        if (len(self.intersecting_indices[self.firstnetworkspeci]) != len(self.read_instance.station_references[self.firstnetworkspeci])):

            # get intersecting station names
            intersecting_station_names = list(set.intersection(*map(set,list(self.non_intersecting_station_names.values()))))

            # only continue wih check if have intersecting station names
            if len(intersecting_station_names) > 0:

                # if any station names per speci are duplicated, then do not count them as intersecting
                names_to_remove = []
                for networkspeci in self.read_instance.station_names:
                    unique_names, unique_counts = np.unique(self.read_instance.station_names[networkspeci], return_counts=True)
                    names_to_remove.extend(unique_names[np.where(unique_counts > 1)[0]])

                names_to_remove = np.unique(names_to_remove)
                intersecting_station_names = np.array([name for name in intersecting_station_names if name not in names_to_remove])

                # get indices of intersecting station names per speci
                for networkspeci in self.read_instance.station_names:

                    station_names_sorted = np.argsort(self.read_instance.station_names[networkspeci])
                    intersect_names_pos = np.searchsorted(self.read_instance.station_names[networkspeci][station_names_sorted], intersecting_station_names)
                    self.intersecting_indices[networkspeci] = np.array(np.append(self.intersecting_indices[networkspeci], 
                                                                                 station_names_sorted[intersect_names_pos]),
                                                                                 dtype=np.int64)

                # validate intersections by longitude/latitude
                self.validate_intersections()

                # sort intersecting indices
                self.sort_intersecting_indices()

                # get non-intersection indices
                self.get_non_intersections()
            
    def by_position(self):
        
        # only do check if have remaining non-intersecting indices
        if (len(self.intersecting_indices[self.firstnetworkspeci]) != len(self.read_instance.station_references[self.firstnetworkspeci])):
            
            # get non-intersecting station coordinates for first speci
            firstnetworkspeci_longitudes = self.non_intersecting_longitudes[self.firstnetworkspeci]
            firstnetworkspeci_latitudes = self.non_intersecting_latitudes[self.firstnetworkspeci]
            if self.read_instance.spatial_colocation_measurement_altitude:
                firstnetworkspeci_measurement_altitudes = self.non_intersecting_measurement_altitudes[self.firstnetworkspeci]
            else:
                firstnetworkspeci_measurement_altitudes = np.zeros(len(firstnetworkspeci_longitudes))

            # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
            # (Earth Centred, Earth Fixed) coordinates assuming WGS84 datum and ellipsoid, and that all heights equal zero
            # ECEF coordinates represent positions (in metres) as X, Y, Z coordinates, approximating the earth surface as an ellipsoid of revolution
            if Version(pyproj.__version__) >= Version("2.0.0"):
                firstnetworkspeci_x, firstnetworkspeci_y, firstnetworkspeci_z = self.transformer.transform(
                    firstnetworkspeci_longitudes, firstnetworkspeci_latitudes, 
                    firstnetworkspeci_measurement_altitudes, radians=False)
            else:
                firstnetworkspeci_x, firstnetworkspeci_y, firstnetworkspeci_z = pyproj.transform(
                    self.lla, self.ecef, firstnetworkspeci_longitudes, firstnetworkspeci_latitudes, 
                    firstnetworkspeci_measurement_altitudes, radians=False)

            # merge coordinates to 3D array
            self.firstnetworkspeci_xyz = np.column_stack((firstnetworkspeci_x, firstnetworkspeci_y, firstnetworkspeci_z))

            # iterate through all other speci, and get intersections (within tolerance) of coordinates with first speci coordinates
            pairwise_intersect_inds = {self.firstnetworkspeci:[]}
            for networkspeci in self.non_intersecting_longitudes:

                if networkspeci == self.firstnetworkspeci:
                    continue

                nextnetworkspeci_longitudes = self.non_intersecting_longitudes[networkspeci]
                nextnetworkspeci_latitudes = self.non_intersecting_latitudes[networkspeci]
                if self.read_instance.spatial_colocation_measurement_altitude:
                    nextnetworkspeci_measurement_altitudes = self.non_intersecting_measurement_altitudes[networkspeci]
                else:
                    nextnetworkspeci_measurement_altitudes = np.zeros(len(nextnetworkspeci_longitudes))

                # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
                if Version(pyproj.__version__) >= Version("2.0.0"):
                    nextnetworkspeci_x, nextnetworkspeci_y, nextnetworkspeci_z = self.transformer.transform( 
                        nextnetworkspeci_longitudes, nextnetworkspeci_latitudes, 
                        nextnetworkspeci_measurement_altitudes, radians=False)
                else:
                    nextnetworkspeci_x, nextnetworkspeci_y, nextnetworkspeci_z = pyproj.transform(self.lla, self.ecef, 
                        nextnetworkspeci_longitudes, nextnetworkspeci_latitudes, 
                        nextnetworkspeci_measurement_altitudes, radians=False)
                
                # merge coordinates to 3D array
                self.nextnetworkspeci_xyz = np.column_stack((nextnetworkspeci_x, nextnetworkspeci_y, nextnetworkspeci_z))            

                # get all indices of next speci xyz coords, within tolerance of each first speci xyz coords
                if Version(scipy.__version__) < Version("1.6.0"):
                    self.idx = cKDTree(self.nextnetworkspeci_xyz).query_ball_point(self.firstnetworkspeci_xyz, self.read_instance.spatial_colocation_tolerance)
                else:
                    self.idx = cKDTree(self.nextnetworkspeci_xyz).query_ball_point(self.firstnetworkspeci_xyz, self.read_instance.spatial_colocation_tolerance, workers=-1)

                # resolve all duplicated matched indices
                self.resolve_duplicate_spatial_colocation_matches()

                # iterate though next speci within tolerance indices, per each first speci coord 
                for idx_ii, idx_l in enumerate(self.idx):
                    # if position has already been resolved (through resolving duplicates) then continue
                    if idx_ii in self.fs_wtol_inds:
                        continue
                    
                    # no matches, then append nothing
                    elif len(idx_l) == 0:
                        continue

                    # just 1 match, then append
                    elif len(idx_l) == 1:
                        self.fs_wtol_inds.append(idx_ii) 
                        self.ns_wtol_inds.append(idx_l[0])

                    # more than 1 match, this shouldn't happen, so throw error
                    else:
                        sys.exit('Error: Spatial colocation could not resolve duplicate matches between station coordinates.')

                # have some interesects between species?
                if len(self.fs_wtol_inds) > 0:

                    # set indices where first species differences are within tolerance, i.e. intersecting 
                    pairwise_intersect_inds['{}_{}'.format(self.firstnetworkspeci, networkspeci)] = self.non_intersecting_indices[self.firstnetworkspeci][np.array(self.fs_wtol_inds)]
                    pairwise_intersect_inds[self.firstnetworkspeci].extend(self.non_intersecting_indices[self.firstnetworkspeci][np.array(self.fs_wtol_inds)])
                
                    # get indices where next species differences are within tolerance, i.e. intersecting 
                    pairwise_intersect_inds[networkspeci] = self.non_intersecting_indices[networkspeci][np.array(self.ns_wtol_inds)]
                # if have no intersects, then append empty array
                else:
                    pairwise_intersect_inds['{}_{}'.format(self.firstnetworkspeci, networkspeci)] = np.array([], dtype=np.int64)
                    pairwise_intersect_inds[networkspeci] = np.array([], dtype=np.int64)

            # get indices (for first networkspeci) where coordinates intersect across all species
            pairwise_intersect_inds_unique, counts = np.unique(pairwise_intersect_inds[self.firstnetworkspeci], return_counts=True)
            pairwise_intersect_inds[self.firstnetworkspeci] = pairwise_intersect_inds_unique[counts == (len(self.read_instance.station_longitudes)-1)]

            # get specific intersect indices across all species, for rest of species
            if len(pairwise_intersect_inds[self.firstnetworkspeci]) > 0:
                for networkspeci in self.non_intersecting_longitudes:
                    if networkspeci == self.firstnetworkspeci:
                        continue
                    _, species_intersect_inds, _ = np.intersect1d(pairwise_intersect_inds['{}_{}'.format(self.firstnetworkspeci, networkspeci)], pairwise_intersect_inds[self.firstnetworkspeci], return_indices=True)
                    pairwise_intersect_inds[networkspeci] = pairwise_intersect_inds[networkspeci][species_intersect_inds]

                # append newly found intersecting indices to previously found intersect inds
                for networkspeci in self.non_intersecting_longitudes:
                    self.intersecting_indices[networkspeci] = np.array(np.append(self.intersecting_indices[networkspeci], pairwise_intersect_inds[networkspeci]), dtype=np.int64)

            # sort intersecting indices
            self.sort_intersecting_indices()

            # get non-intersection indices
            self.get_non_intersections()

    def validate_intersections(self):

        """ Double check each station colocation by longitude / latitude position.
            If any of the the differences in longitude / latitude for a station across the relevant 
            longitudes / latitudes exceed a certain tolerance, then drop station.
        """

        # only do validation if appropriate variable is True
        if self.read_instance.spatial_colocation_validation:

            # iterate through stations in first networkspeci, and get longitudes / latitudes of each station across
            # networkspecies
            stn_inds_to_drop = np.array([], dtype=np.int64)
            for stn_ii in range(len(self.intersecting_indices[self.firstnetworkspeci])):
                station_lons = np.array([])
                station_lats = np.array([])
                for networkspeci in self.networkspecies:
                    station_lons = np.append(station_lons, self.read_instance.station_longitudes[networkspeci][self.intersecting_indices[networkspeci][stn_ii]])
                    station_lats = np.append(station_lats, self.read_instance.station_latitudes[networkspeci][self.intersecting_indices[networkspeci][stn_ii]])
                
                # convert speci longitude and latitudes in geographic coordinates to cartesian ECEF 
                station_z = np.zeros(len(station_lons))
                if Version(pyproj.__version__) >= Version("2.0.0"):
                    station_lons_ecef, station_lats_ecef, _ = self.transformer.transform(
                        station_lons, station_lats, station_z, radians=False)
                else:
                    station_lons_ecef, station_lats_ecef, _ = pyproj.transform(
                        self.lla, self.ecef, station_lons, station_lats, station_z, radians=False)

                # if the maximum difference in station longitudes or latitudes exceeds a certain tolerance, then drop
                # the station form intersections
                lon_diff_max = np.max(np.abs(station_lons_ecef - station_lons_ecef[:, np.newaxis]))
                lat_diff_max = np.max(np.abs(station_lats_ecef - station_lats_ecef[:, np.newaxis]))
                if (lon_diff_max > self.read_instance.spatial_colocation_validation_tolerance) or (lat_diff_max > self.read_instance.spatial_colocation_validation_tolerance):
                    stn_inds_to_drop = np.append(stn_inds_to_drop, stn_ii)

            # update intersecting indices
            for networkspeci in self.networkspecies:
                self.intersecting_indices[networkspeci] = np.array([inter_ii for stn_ii,inter_ii in enumerate(self.intersecting_indices[networkspeci]) if stn_ii not in stn_inds_to_drop], dtype=np.int64)

    def sort_intersecting_indices(self):
        """ get the order of the first networkspecies in ascending order, and order other networkspecies in the 
            same order. 
        """
        
        # sort first networkspecies in ascending order, and order other networkspecies in the same order
        sorted_inds = np.argsort(self.intersecting_indices[self.networkspecies[0]])
        for networkspeci in self.networkspecies:
            self.intersecting_indices[networkspeci] = self.intersecting_indices[networkspeci][sorted_inds]

    def get_non_intersections(self):

        """ get the indices of stations where no-intersections have been found per networkspecies, and
            associated metadata variables.
        """

        # get the indices of stations where no-intersections have been found per networkspecies
        self.non_intersecting_indices = {networkspeci: np.setdiff1d(np.arange(len(self.read_instance.station_references[networkspeci])), self.intersecting_indices[networkspeci]) for networkspeci in self.networkspecies}

        # get associated metadata variables for stations with no-intersections
        if self.read_instance.spatial_colocation_station_name:
            self.non_intersecting_station_names = {networkspeci: self.read_instance.station_names[networkspeci][self.non_intersecting_indices[networkspeci]] for networkspeci in self.read_instance.station_names}
        if self.read_instance.spatial_colocation_longitude_latitude:
            self.non_intersecting_longitudes = {networkspeci: self.read_instance.station_longitudes[networkspeci][self.non_intersecting_indices[networkspeci]] for networkspeci in self.read_instance.station_longitudes}
            self.non_intersecting_latitudes = {networkspeci: self.read_instance.station_latitudes[networkspeci][self.non_intersecting_indices[networkspeci]] for networkspeci in self.read_instance.station_latitudes}
        if self.read_instance.spatial_colocation_measurement_altitude:
            self.non_intersecting_measurement_altitudes = {networkspeci: self.read_instance.station_measurement_altitudes[networkspeci][self.non_intersecting_indices[networkspeci]] for networkspeci in self.read_instance.station_measurement_altitudes}

    def resolve_duplicate_spatial_colocation_matches(self):

        """ Method that resolves duplicate indices found during spatial colocation of 2 species.

            In spatial colocation it is neccessary to match each stations geographically within a specific
            tolerance, for 2 different species. 

            In some cases the stations within the tolerance can match for multiple stations. 
            In order to resolve this, for each duplicated station, it is set to match with the closest station to it
            in terms of 3D distance. This is done iteratively until there are no more duplicates.
        """

        # determine which indices of next species are found to match with mutiple first speci coords, and which not
        unique_idx, unique_idx_counts = np.unique(list(itertools.chain(*self.idx)), return_counts=True)
        nondup_idx = unique_idx[np.where(unique_idx_counts == 1)[0]]
        dup_idx = unique_idx[np.where(unique_idx_counts > 1)[0]]

        # resolve all duplicated matched indices by finding for which index the distance is closest
        # keep on iterating until you have no more duplicates, or it is impossible to resolve them
        previous_n_dups = len(dup_idx) 

        # gather list of indices for first speci and next speci of already resolved within tolerance indices
        self.fs_wtol_inds = []
        self.ns_wtol_inds = []

        while True:

            # before iterating over duplicate next speci indices, iterate over each of first speci coord with multiple matches
            # resolve simple cases where the closest distance is a non-duplicate index
            for idx_ii, idx_l in enumerate(self.idx):
                if len(idx_l) > 1:
                    idx_l = np.array(idx_l) 
                    # find the dists between duplicate matched next speci xyz coords and first speci coord
                    if Version(scipy.__version__) < Version("1.6.0"):
                        dists = cKDTree([self.firstnetworkspeci_xyz[idx_ii]]).query(self.nextnetworkspeci_xyz[idx_l], k=1)[0]
                    else:
                        dists = cKDTree([self.firstnetworkspeci_xyz[idx_ii]]).query(self.nextnetworkspeci_xyz[idx_l], k=1, workers=-1)[0]
                    # order indices by closest dists
                    ordered_idx_l = idx_l[np.argsort(dists)]
                    if ordered_idx_l[0] in nondup_idx:
                        self.fs_wtol_inds.append(idx_ii)
                        self.ns_wtol_inds.append(ordered_idx_l[0])
                        self.idx[idx_ii] = [ordered_idx_l[0]]

            # update nondup_idx and dup_idx
            unique_idx, unique_idx_counts = np.unique(list(itertools.chain(*self.idx)), return_counts=True)
            nondup_idx = unique_idx[np.where(unique_idx_counts == 1)[0]]
            dup_idx = unique_idx[np.where(unique_idx_counts > 1)[0]]

            # create a list of storing matched indices that cannot resolve on this pass
            unresolved_dup_idx = []

            # iterate over duplicate next speci indices
            for dup_index in dup_idx:

                # get all relevant first speci indices for which contain the duplicate next speci index
                relevant_idx = np.array([idx_ii for idx_ii, idx_l in enumerate(self.idx) if dup_index in idx_l])

                # if have no relevant_idx values, it is because in previous iterations the dup_index was resolved
                # continue to the next dup_index
                if len(relevant_idx) == 0:
                    continue

                # find the dists between all relevant first speci xyz coords and next speci duplicate xyz coord
                if Version(scipy.__version__) < Version("1.6.0"):
                    dists = cKDTree([self.nextnetworkspeci_xyz[dup_index]]).query(self.firstnetworkspeci_xyz[relevant_idx], k=1)[0]
                else:
                    dists = cKDTree([self.nextnetworkspeci_xyz[dup_index]]).query(self.firstnetworkspeci_xyz[relevant_idx], k=1, workers=-1)[0]
                # order the relevant first speci indices by dists (closest first)
                ordered_relevant_idx = relevant_idx[np.argsort(dists)]
                # iterate through ordered relevant first speci indices until find an ind for which can claim duplicate ind
                # keep going in order iteratively until have exhausted all options for duplicate ind
                for ordered_relevant_index_ii, ordered_relevant_index in enumerate(ordered_relevant_idx):
                
                    # if current relevant first speci index has no competing indices, then append matched index,
                    # remove duplicate index from other match lists, and then break out of iteration
                    if len(self.idx[ordered_relevant_index]) == 1:
                        self.fs_wtol_inds.append(ordered_relevant_index)
                        self.ns_wtol_inds.append(dup_index)
                        for next_ordered_relevant_index in ordered_relevant_idx[ordered_relevant_index_ii+1:]:
                            self.idx[next_ordered_relevant_index].remove(dup_index)
                        break

                    # otherwise, find the dists between the relevant first speci xyz coord and next speci duplicate xyz coords
                    else:
                        idx_l = np.array(self.idx[ordered_relevant_index])
                        if Version(scipy.__version__) < Version("1.6.0"):
                            dists = cKDTree([self.firstnetworkspeci_xyz[ordered_relevant_index]]).query(self.nextnetworkspeci_xyz[idx_l], k=1)[0]
                        else:
                            dists = cKDTree([self.firstnetworkspeci_xyz[ordered_relevant_index]]).query(self.nextnetworkspeci_xyz[idx_l], k=1, workers=-1)[0]
                        ordered_idx_l = idx_l[np.argsort(dists)]

                        # if the closest index is the dup index, or is another non-duplicate, then append it
                        # if the closest index is the dup index, then remove duplicate index from other match lists 
                        # and then break out of iteration
                        if (ordered_idx_l[0] == dup_index) or (ordered_idx_l[0] in nondup_idx):
                            self.fs_wtol_inds.append(ordered_relevant_index)
                            self.ns_wtol_inds.append(ordered_idx_l[0])
                            if ordered_idx_l[0] == dup_index:
                                self.idx[ordered_relevant_index] = [ordered_idx_l[0]]
                                for next_ordered_relevant_index in ordered_relevant_idx[ordered_relevant_index_ii+1:]:
                                    self.idx[next_ordered_relevant_index].remove(dup_index)
                                break
                            elif ordered_idx_l[0] in nondup_idx:
                                self.idx[ordered_relevant_index] = [ordered_idx_l[0]]
                            
                        # else, if the closest index is another duplicated index,
                        # then cannot currently resolve position, so come back to this later
                        else:
                            unresolved_dup_idx.append(dup_index)
                            break

            # if have no more duplicated matched indices to resolve, then break out of loop and function
            if len(unresolved_dup_idx) == 0:
                return
            # otherwise keep on iterating until there are no more duplicates, or it is impossible to resolve duplicates
            else:
                # break out of function if have same number of duplicates that there was on the previous loop (i.e. impossible to resolve duplicates)
                if previous_n_dups == len(unresolved_dup_idx):
                    return

                # update nondup_idx and dup_idx
                unique_idx, unique_idx_counts = np.unique(list(itertools.chain(*self.idx)), return_counts=True)
                nondup_idx = unique_idx[np.where(unique_idx_counts == 1)[0]]
                dup_idx = unique_idx[np.where(unique_idx_counts > 1)[0]]
                previous_n_dups = len(dup_idx)