import numpy as np
import csv
from scipy.spatial import KDTree
import networkx as nx
from KBO import *
import linkutils
import itertools
import uuid
import os.path
import cPickle as pickle
from TNOcandidate import TNOcandidate
from Triplet import Triplet
from pympler import tracker


class TNOfinder(object):
    def __init__(self, objcat, bands = None, exclude_objids=None, date_start=None, look_ahead_nights=30,
                nominal_distance=60, astrometric_err=0.15, runid=None, classifier=None):

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
        self._version = '1.0.0'
        self.astrometric_err = astrometric_err
        self.set_catalog(objcat, bands, exclude_objids, date_start)        
        self.look_ahead_nights = look_ahead_nights   # how many visits to look ahead for linking
        self.nominal_distance = nominal_distance  # AU
        self.vmax = 150 # arcsec/day
        self.para = linkutils.exposure_parallax()        # get dlon, dlat, vlon, vlat for each exposure
        self.good_triplets = []
        self.graph = nx.Graph()
        self.candidates = []
        self.linkpoints = []
        self.runid=runid        # identifier for this linking run
        self.clf = classifier
#        self.memory_tracker = tracker.SummaryTracker()
 #       self.memory_tracker.print_diff()

        
    
    def set_catalog(self, cat, bands=None, exclude_objids=None, date_start=None):
        self.objects = cat
        self.bands = bands
        if exclude_objids is not None: self.objects = Catalog(obj for obj in self.objects if obj.objid not in exclude_objids)
        if bands is not None: self.objects = Catalog(obj for obj in self.objects if obj.band in bands)
        self.nites = sorted(set(point.nite for point in self.objects))
        if date_start is not None: self.objects = Catalog(obj for obj in self.objects if obj.date>=date_start)
        self.objects.add_constant('obscode', 807)
        self.objects.add_constant('err', self.astrometric_err)
        self.objects.orderby('expnum')


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

    def getDatePos(self, objid):
        '''
        Return date, ra, dec for a given objid
        '''
        point = [p for p in self.objects if p.objid==objid][0]
        ra = point.ra
        dec = point.dec
        date = point.date
        return date, ra, dec

    def vEcliptic(self, id1, id2):
        '''
        Compute the components of velocity between two points in ecliptic coords. Units: arcsec/day.
        '''
        date1, ra1, dec1 = self.getDatePos(id1)
        date2, ra2, dec2 = self.getDatePos(id2)
        ecl1 = ephem.Ecliptic(ephem.Equatorial(ra1, dec1))
        ecl2 = ephem.Ecliptic(ephem.Equatorial(ra2, dec2))
        lat1, lon1 = ecl1.lat, ecl1.lon
        lat2, lon2 = ecl2.lat, ecl2.lon
        dt = date2-date1
        dlon = np.mod(np.abs(lon2-lon1), 2*np.pi)
        if dlon > np.pi: dlon = 2*np.pi - dlon
        vlon = dlon/dt*np.cos(lat1)*180/np.pi*3600
        vlat = (lat2-lat1)/dt*180/np.pi*3600
        return vlon, vlat

    def make_vmap(self, linkmap, report_interval=100):
        '''
        Stuffs all the velocities for the points in the link table into a dictionary that can be exported.
        '''
        vlist = []
        npts=len(linkmap.keys())
        i=0
        for k in linkmap.keys():
            if report_interval>0:
                if i % report_interval == 0: print 'Building velocity table for point ', i, ' of ', npts
                i+=1
            v = [self.vEcliptic(k, ind2) for ind2 in linkmap[k]]
            vlist.append(dict(zip(linkmap[k], v)))
        return dict(zip(linkmap.keys(), vlist))
    

    def find_linkable(self, report_interval=100):
        linkable_ids = {}
        i=0
        npts = len(self.objects)
        for point in self.objects:
            if report_interval>0:
                if i % report_interval == 0: print 'Building link table for point ', i, ' of ', npts
                i+=1
            linkable_points = self.link_obj(point)
            linkable_ids[point.objid]=[p.objid for p in linkable_points]
        return linkable_ids

    def tno_like(self, point1, point2, debug=False):
        lon1, lat1 = Ecliptic(Equatorial(point1.ra, point1.dec)).get()
        lon2, lat2 = Ecliptic(Equatorial(point2.ra, point2.dec)).get()
        displacement = ephem.separation((lon2, lat2), (lon1, lat1))
        displacement_asec = displacement*180/np.pi*3600
        if point2.date != point1.date:
            velocity = displacement_asec/(point2.date - point1.date)
        else: 
            velocity=9999
        this_dlon, this_dlat, this_vlon, this_vlat = linkutils.parallax(point1.ra, point1.dec, point1.date) 
        next_dlon, next_dlat, next_vlon, next_vlat = linkutils.parallax(point2.ra, point2.dec, point2.date)   
        dlon, dlat = next_dlon-this_dlon, next_dlat-this_dlat
        dot = np.cos(lat1)**2*(lon2 - lon1)*dlon + (lat2 - lat1)*dlat
        norm = np.sqrt(np.cos(lat1)**2*dlon**2 + dlat**2)
        cosine = dot/(norm*displacement)
        TNOlike = True if velocity<self.vmax and cosine>self.cosine_cut(displacement_asec) else False
        if debug:
            linkinfo = {'v':velocity, 'cos':cosine, 'dot':dot, 'displacement':displacement_asec, 'cut':self.cosine_cut(displacement_asec),
                     'point1':point1, 'point2':point2, 'lon1':lon1, 'lon2':lon2, 'lat1':lat1, 'lat2':lat2, 'dlon':dlon, 'dlat':dlat, 'norm':norm,
                     'TNOlike':TNOlike}
            self.linkpoints.append(linkinfo)
        return TNOlike
       
    def link_obj(self, point, verbose=False):
        nites = self.nites
        thisnite = point.nite
        if verbose: print 'Linking point ', point.objid, point.date, thisnite, point.ra, point.dec, point.band
        look_ahead_nites = sorted([n for n in nites if 0<self.nites_between(thisnite,n)<self.look_ahead_nights])
        lon, lat = Ecliptic(Equatorial(point.ra, point.dec)).get()
#        this_dlon, this_dlat = self.para[point.expnum]['dlon'], self.para[point.expnum]['dlat']
        this_dlon, this_dlat, this_vlon, this_vlat = linkutils.parallax(point.ra, point.dec, point.date) 
        next_obj = []
        for next_nite in look_ahead_nites:
            if verbose: print 'Linking target nite: ', next_nite
            deltaT = self.nites_between(thisnite,next_nite)
            next_dlon, next_dlat, next_vlon, next_vlat = linkutils.parallax(point.ra, point.dec, point.date+deltaT)
            dlon, dlat = next_dlon-this_dlon, next_dlat-this_dlat
            search_center = Equatorial(Ecliptic(lon+dlon/self.nominal_distance, lat+dlat/self.nominal_distance))
            deltaR = np.sqrt(dlon**2*np.cos(lat)**2 + dlat**2)/self.nominal_distance
#            print ephem.separation((search_center.ra, search_center.dec), (point.ra, point.dec))*3600*180/np.pi/deltaT
            sep_max = self.vmax*np.pi/(180*3600)*(deltaT)/2
            current_objects = [obj for obj in self.objects if obj.nite == next_nite \
                               and (self.bands is None or obj.band in self.bands) \
                                and ephem.separation((search_center.ra, search_center.dec), (obj.ra, obj.dec))<deltaR]
            if verbose: print '   Found ', len(current_objects), ' points in search window'
            # now the real work: for each object, test to see
            # if it's consistent with being the next point in
            # a KBO trajectory, i.e. consistent in direction and displacement with earth parallax.
            if True:
                for point2 in current_objects:
                    if verbose: print 'Examining point ', point2.objid, point2.date, thisnite, point2.ra, point2.dec, point.band, ' ... ', 
    #                next_dlon, next_dlat = self.para[point2.expnum]['dlon'], self.para[point2.expnum]['dlat']
    #                dlon, dlat = next_dlon-this_dlon, next_dlat-this_dlat
    #                if self.tno_like(point, point2, lon, lat, dlon, dlat, debug=True): 
                    if self.tno_like(point, point2, debug=False):
                        if verbose: print 'Point is tno_like...'
                        next_obj.append(point2)
                    else:
                        if verbose: print 'Point NOT tno_like...'
        return next_obj

    def getTriplets(self, linkmap, everyN=1):
        alltriplets = []
        ntriplets = 0
        npts = 0
        for ind1 in linkmap.keys():
            if npts % everyN == 0:
                next_inds = linkmap[ind1]
                for ind2 in next_inds:
                    next_next_inds = linkmap[ind2]
                    for ind3 in next_next_inds:
                        triplet = [ind1, ind2, ind3]
                        alltriplets.append(triplet)
                        ntriplets +=1
            npts +=1
        return alltriplets

    def find_triplets(self, dataframe, allow_unbound=False, verbose=False, link_table=None, vmap_table=None):
        current_nite = 0
#        self.memory_tracker.print_diff()
        good_triplets = []
        good_ids = []
        if link_table is None:
            if verbose:
                print '*** TNOfinder stage 1: Building link table ***'
            linkable_ids = self.find_linkable()
            linkmap_out = os.path.join('linker_output',str(self.runid)+'_linkmap.pickle')
            with open(linkmap_out,'wb') as f:
                if verbose:
                    print '   ----> Writing link table to file ', linkmap_out
                pickle.dump(linkable_ids, f)
        else:
            if verbose:
                print '*** TNOfinder stage 1: Reading link table from file '+link_table+' ***'
            with open(link_table,'rb') as f:
                linkable_ids = pickle.load(f)
        if verbose:
            print '*** TNOfinder: Building link table complete ***'
            print

        if vmap_table is None:
            if verbose:
                print '*** TNOfinder stage 2: Building velocity map ***'
            vmap = self.make_vmap(linkable_ids)
            vmap_out = os.path.join('linker_output',str(self.runid)+'_vmap.pickle')
            with open(vmap_out,'wb') as fout:
                if verbose:
                    print '   ----> Writing velocity map to file ', vmap_out
                pickle.dump(vmap, fout)
        else:
            if verbose:
                print '*** TNOfinder stage 2: Reading velocity map from file ' + vmap_table + ' ***'
            with open(vmap_table,'rb') as fv:
                vmap = pickle.load(fv)


        if verbose:
            print '*** TNOfinder: Building velocity table complete ***'
            print
            print '*** TNOfinder stage 3: Finding triplets ***'
        ipts = 0
        if True:   # this is here in case we only want to build the link and velocity tables. Could handle this scenario through a config flag...
            for obj1 in self.objects:
                if obj1.nite>current_nite:
                    if verbose:
                        if current_nite>0: 
                            self.report_state(current_nite, good_triplets)
#                            self.memory_tracker.print_diff()
                        print 'Linking points from nite: ', obj1.nite
                    current_nite = obj1.nite 
                next_points = [p for p in self.objects if p.objid in linkable_ids[obj1.objid]]
                for obj2 in next_points:
    #                if verbose: print '.',
                    next_next_points = [p for p in self.objects if p.objid in linkable_ids[obj2.objid]]
                    for obj3 in next_next_points:
    #                    if verbose: print '-',
                        ipts +=1
                        #
                        #  Insert triplet-checking here (reject the really bad ones before orbit-fitting):
                        idlist = [obj1.objid, obj2.objid, obj3.objid]
                        t = Triplet(idlist, vmap, dataframe, classifier=self.clf)
                        if t.checkTriplet2() == 1 or self.clf is None:  # good triplet according to machine-learning classifier if we have one
                            triple=Catalog([obj1, obj2, obj3])
                            orbit = Orbit(triple)
                            if (orbit.chisq<5 and orbit.ndof==1 and orbit.elements['a']>self.nominal_distance/2):
                                good_triplets.append(triple)
                                self.good_triplets.append(triple)
#                                if verbose: 
#                                    print 'New good triplet:'
#                                    self.report_state(current_nite, [good_triplets[-1]])
                            if allow_unbound and (orbit.ndof==0 and orbit.elements['a']>self.nominal_distance/2):
                                good_triplets.append(triple)
                                self.good_triplets.append(triple)
        if verbose:
            print '*** TNOfinder: Triplet search complete. '
            print
            print '*** TNOfinder stage 4: merging triplets, building final candidate list ***'
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
        cat.add_constant('err', self.astrometric_err)
        return cat

    def distinct_chains(self, cat):
        '''
        A catalog built from connected observations may contain more than one point per exposure, for example if 
        a point links to two closely-spaced points in a subsequent exposure. This routine returns a list of all catalogs 
        that can be built from the input catalog containing only one point per exposure. Note that if the catalog
        does not have any duplicate exposure numbers, the result is a list containing only the input catalog
        '''
        exps = set([p.expnum for p in cat])
        objids = [[p.objid for p in cat if p.expnum==e] for e in exps]
        combos = [oid for oid in itertools.product(*objids)]   # itertools wizardry
        return [Catalog([p for p in cat if p.objid in combo]) for combo in combos]
 
    def tnocands(self, triplets):
        subs = self.subgraphs(self.make_graph(triplets))
        cats = [self.connected_observations(sub) for sub in subs]
        observations = []
        for cat in cats:
            unique_cats = self.distinct_chains(cat)
            for u in unique_cats:
                observations.append(u)
        orbits = [Orbit(obs) for obs in observations]
        return zip(orbits, observations)
    
    def report_state(self, nite, triplets):
        cands = self.tnocands(triplets)
        print 'Status after nite ', nite
        print 'Good triplets found: ', len(cands)
        for cand in cands:
            orbit, observations = cand[0], cand[1]
            if orbit.chisq/orbit.ndof<2:
                tc = TNOcandidate(orbit, observations, runid=self.runid, csvname=os.path.join('linker_output',str(self.runid)+'_rolling_'+str(nite)))
                print 'Chi^2, ndof: ', round(orbit.chisq,3), orbit.ndof
                print 'Elements: ', orbit.elements
                print 'Element errors:', orbit.elements_errs
                print 'Observations: '
                print observations.writes()
                print
        print '--------------------------------------------------------------------------------------------'
