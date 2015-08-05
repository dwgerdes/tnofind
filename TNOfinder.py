class ObjectLinker(object):
    def __init__(self, target, objcat, band = None, exclude_snobjids=None, date_start=None, look_ahead_visits=2):

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
        if exclude_snobjids is not None: self.objects = Catalog(obj for obj in self.objects if obj.objid not in exclude_snobjids)
#
        if date_start is not None: self.objects = Catalog(obj for obj in self.objects if obj.date>=date_start)
        self.objects.add_constant('obscode', 807)
        self.objects.add_constant('err', astrometric_err)
        self.objects.orderby('date')
        self.field = fields[target]
        self.band = band
        self.look_ahead_visits = look_ahead_visits   # how many visits to look ahead for linking
        self.vmax = 150    # arcsec/day
        self.nominal_distance = 40  # AU
#        self.trees_by_nite, self.points_by_nite = trees_points_by_nite(objcat)
        
    def points_near(self, center, r, nite):
        x = [center.ra, center.dec]
        near_ids = self.trees_by_nite[nite].query_ball_point(x, r)
        near_points = [self.points_by_nite[nite][i] for i in near_ids]
        return near_points
        
    
    def cosine_cut(self, displacement):
        '''
        The requirement on directional alignment between two points and that expected from earth reflex motion.
        It's more generous for lower values of the displacement (in arcseconds).
        '''
        if displacement<400:
            return 0.7
        elif 400<displacement<=800:
            return 0.7 + (0.95-0.7)/(800-400)*(displacement-400)
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
#        norm = np.sqrt(np.cos(self.centerlat)**2*dlon**2 + dlat**2)
#        dot = np.cos(self.centerlat)**2*(lon2 - lon1)*dlon + (lat2 - lat1)*dlat
        norm = np.sqrt(np.cos(lat1)**2*dlon**2 + dlat**2)
        dot = np.cos(lat1)**2*(lon2 - lon1)*dlon + (lat2 - lat1)*dlat
        if point2.date != point1.date:
            velocity = displacement/(point2.date - point1.date)
        else: 
            velocity=9999
        cosine = dot/(norm*displacement)
        # The old way:
#        return True if cosine > np.cos(20*np.pi/180) and velocity < self.vmax*np.pi/(180*3600) else False 
#       The new way accounts for the noise in cosine for smaller values of displacement.
        displacement_asec = displacement*180/np.pi*3600
        return True if velocity<self.vmax and cosine>self.cosine_cut(displacement_asec) else False
       
    def link_obj(self, point):
        visits = self.field.visits
        visits.orderby('nite')
        thisvisit = visits[point.nite]
        this_dlon, this_dlat = para[point.expnum]['dlon'], para[point.expnum]['dlat']
        this_vlon, this_vlat = para[point.expnum]['vlon'], para[point.expnum]['vlat']
        nextvisit = thisvisit
        next_obj = []
        for i in range(self.look_ahead_visits):
            
            # get the next visit
            try:
                nextvisit = next(visit for visit in visits if visit.nite > nextvisit.nite)
            except StopIteration: pass
            if nextvisit is None: 
                print "nextvisit is None!"
                return []
        
            days_between_visits = ephem.date(pretty_nite(nextvisit.nite))-ephem.date(pretty_nite(thisvisit.nite))
            if days_between_visits>0:
                deltalon, deltalat = this_vlon/self.nominal_distance*days_between_visits*np.pi/180, this_vlat/self.nominal_distance*days_between_visits*np.pi/180
                lon, lat = Ecliptic(Equatorial(point.ra, point.dec)).get()
                search_center = Equatorial(Ecliptic(lon+deltalon, lat+deltalat))
#                print ephem.separation((search_center.ra, search_center.dec), (point.ra, point.dec))*3600*180/np.pi/days_between_visits
                sep_max = self.vmax*np.pi/(180*3600)*(days_between_visits)/2
                current_objects = [obj for obj in self.objects if obj.nite == nextvisit.nite \
                                   and (self.band is None or obj.band == self.band) \
                                    and ephem.separation((search_center.ra, search_center.dec), (obj.ra, obj.dec))<sep_max]
# This method uses KDtrees. Our current implementation is slower than the list-comprehension method above.        
#                current_objects = self.points_near(search_center, sep_max, nextvisit.nite)
                # now the real work: for each object, test to see
                # if it's consistent with being the next point in
                # a KBO trajectory.
                dlon, dlat = nextvisit.dlon - thisvisit.dlon, nextvisit.dlat - thisvisit.dlat
                for point2 in current_objects:
                    next_dlon, next_dlat = para[point2.expnum]['dlon'], para[point2.expnum]['dlat']
                    dlon, dlat = next_dlon-this_dlon, next_dlat-this_dlat
                    if self.tno_like(point, point2, lon, lat, dlon, dlat): next_obj.append(point2)
        return next_obj
