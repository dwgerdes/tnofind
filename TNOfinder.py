import numpy as np
import csv
from scipy.spatial import KDTree
import networkx as nx
from KBO import *
from linktools import *

class TNOfinder(object):
    def __init__(self, objcat, band = None, exclude_objids=None, date_start=None, look_ahead_nights=30,
                nominal_distance=60):

# objcat is a Catalog (for now) which must contain at least the following columns:
#     date (as a DateTime object)
#     nite
#     ra, dec in radians or ephem.Angle objects
#     band
#     expnum
#     exptime (in seconds)
#     objid (a unique identifier, such as snobjid)
#     mag
#     ccd
        self.objects = objcat
        if band is not None: self.objects = Catalog(obj for obj in self.objects if obj.band == band)
        if exclude_objids is not None: self.objects = Catalog(obj for obj in self.objects if obj.objid not in exclude_objids)
#
        if date_start is not None: self.objects = Catalog(obj for obj in self.objects if obj.date>=date_start)
        self.objects.add_constant('obscode', 807)
        self.objects.add_constant('err', astrometric_err)
        self.objects.orderby('date')
        self.band = band
        self.look_ahead_nights = look_ahead_nights   # how many visits to look ahead for linking
        self.nominal_distance = nominal_distance  # AU
        self.vmax = 150 # arcsec/day
        self.nites = sorted(set(point.nite for point in self.objects))
        self.para = exposure_parallax()        # get dlon, dlat, vlon, vlat for each exposure
        self.good_triplets = []
        self.graph = nx.Graph()
        self.candidates = []
        
    
    def nites_between(self, nite1, nite2):
        '''
        Returns the time in days between two nites, where nite2>nite1.
        '''
        return ephem.date(pretty_nite(nite2))-ephem.date(pretty_nite(nite1))
        
    def cosine_cut(self, displacement_arcsec):
        '''
        The requirement on directional alignment between two points and that expected from earth reflex motion.
        It's more generous for lower values of the displacement (in arcseconds).
        '''
        if displacement_arcsec<400:
            return 0.7
        elif 400<displacement_arcsec<=800:
            return 0.7 + (0.95-0.7)/(800-400)*(displacement_arcsec-400)
        else:
            return 0.95
    
    def tno_like(self, point1, point2, lon1, lat1, dlon, dlat):
        '''
        returns True if point1 and point2 are in a TNO-like relationship to each other, i.e. displacement
        magnitude and direction consistent with earth reflex motion. 
        point1 = point in earlier visit
        point2 = point in subsequent visit that we are trying to link to point1.
        lon1, lat1 = ecliptic coords of point1
        dlon, dlat = expected change in direction from one visit to next: dlon = nextvisit.dlon - thisvisit.dlon, etc.
        '''
        lon2, lat2 = Ecliptic(Equatorial(point2.ra, point2.dec)).get()
        displacement = ephem.separation((lon2, lat2), (lon1, lat1))
        norm = np.sqrt(np.cos(lat1)**2*dlon**2 + dlat**2)
        dot = np.cos(lat1)**2*(lon2 - lon1)*dlon + (lat2 - lat1)*dlat
        if point2.date != point1.date:
            velocity = displacement/(point2.date - point1.date)
        else: 
            velocity=9999
        cosine = dot/(norm*displacement)
        displacement_asec = displacement*180/np.pi*3600
        return True if velocity<self.vmax and cosine>self.cosine_cut(displacement_asec) else False
       
    def link_obj(self, point):
        nites = self.nites
        thisnite = point.nite
        look_ahead_nites = sorted([n for n in nites if 0<self.nites_between(thisnite,n)<self.look_ahead_nights])
        lon, lat = Ecliptic(Equatorial(point.ra, point.dec)).get()
        this_dlon, this_dlat = self.para[point.expnum]['dlon'], self.para[point.expnum]['dlat']
        next_obj = []
        for next_nite in look_ahead_nites:
            deltaT = self.nites_between(thisnite,next_nite)
            next_dlon, next_dlat, next_vlon, next_vlat = parallax(point.ra, point.dec, point.date+deltaT)
            dlon, dlat = next_dlon-this_dlon, next_dlat-this_dlat
            search_center = Equatorial(Ecliptic(lon+dlon/self.nominal_distance, lat+dlat/self.nominal_distance))
            deltaR = np.sqrt(dlon**2*np.cos(lat)**2 + dlat**2)/self.nominal_distance
#            print ephem.separation((search_center.ra, search_center.dec), (point.ra, point.dec))*3600*180/np.pi/deltaT
            sep_max = self.vmax*np.pi/(180*3600)*(deltaT)/2
            current_objects = [obj for obj in self.objects if obj.nite == next_nite \
                               and (self.band is None or obj.band == self.band) \
                                and ephem.separation((search_center.ra, search_center.dec), (obj.ra, obj.dec))<deltaR]
            # now the real work: for each object, test to see
            # if it's consistent with being the next point in
            # a KBO trajectory, i.e. consistent in direction and displacement with earth parallax.
            for point2 in current_objects:
                next_dlon, next_dlat = self.para[point2.expnum]['dlon'], self.para[point2.expnum]['dlat']
                dlon, dlat = next_dlon-this_dlon, next_dlat-this_dlat
                if self.tno_like(point, point2, lon, lat, dlon, dlat): next_obj.append(point2)
        return next_obj
    
    def find_triplets(self):
        cat1 = self.objects
        current_nite = 0
        good_triplets = []
        for obj1 in cat1:
            if obj1.nite>current_nite:
                self.report_state(current_nite, good_triplets)
                print 'Linking points from nite: ', obj1.nite
                current_nite = obj1.nite
            next_points = self.link_obj(obj1)
            for obj2 in next_points:
#                print '.',
                next_next_points = self.link_obj(obj2)
                for obj3 in next_next_points:
#                    print '-',
                    triple=Catalog([obj1, obj2, obj3])
                    orbit = Orbit(triple)
                    if orbit.chisq<2 and orbit.ndof==1 and orbit.elements['a']>20:    # good triplet
                        good_triplets.append(triple)
        self.good_triplets = good_triplets
        self.candidates = self.tnocands(good_triplets)
        return good_triplets

    def make_graph(self, triplets):
        G = nx.Graph()
        for triplet in triplets:
            observations = sorted([obs.objid for obs in triplet])
            G.add_nodes_from(observations)
            G.add_edge(observations[0], observations[1])
            G.add_edge(observations[1], observations[2])
        return G

    def subgraphs(self, graph):
        return [s for s in nx.connected_component_subgraphs(graph)]

    def connected_observations(self, subgraph):
        objids = [node for node in nx.nodes_iter(subgraph)]
        points = [p for p in self.objects if p.objid in objids]
        cat = Catalog(points)
        cat.add_constant('obscode', 807)
        cat.add_constant('err', astrometric_err)
        return cat
    
    def tnocands(self, triplets):
        subs = self.subgraphs(self.make_graph(triplets))
        observations = [self.connected_observations(sub) for sub in subs]
        orbits = [Orbit(obs) for obs in observations]
        return zip(orbits, observations)
    
    def report_state(self, nite, triplets):
        cands = self.tnocands(triplets)
        print 'Status after nite ', nite
        print 'Good triplets found: ', len(cands)
        for cand in cands:
            orbit, observations = cand[0], cand[1]
            print 'Chi^2, ndof: ', orbit.chisq, orbit.ndof
            print 'Elements: ', orbit.elements
            print 'Element errors:', orbit.elements_errs
            print 'Observations: '
            print observations.writes()
            print
        print '--------------------------------------------------------------------------------------------'