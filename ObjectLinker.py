from KBO import *
import pickle
import ephem
import numpy as np

class ObjectLinker(object):
    def __init__(self, field, band = None, filename = None, exclude_snobjids=None, date_start=None):
        if filename is None: filename = 'objects/kbo_data_Y2_' + field[1] + '.obj'
        with open(filename) as objfile:
            self.objects = Catalog(pickle.load(objfile))
#            self.objects = Catalog(filename, date=DateTime,ra=hours,dec=degrees,expnum=int,exptime=float,band=str,ccdnum=int,mag=float,\
#                                   ml_score=float,snobjid=int,fwhm=float,t_eff=float)

        self.objects.rename(date='DATE_OBS', ra='RA', dec='DEC', band='BAND', mag='MAG', ccd ='CCDNUM',
                            snobjid='SNOBJID', expnum = 'EXPNUM', exptime='EXPTIME')
        self.objects.refactor('date', toDateTime)
        self.objects.refactor('ra', lambda ra: hours(ra * np.pi/180))
        self.objects.refactor('dec', lambda dec: degrees(dec * np.pi/180))
        if band is not None: self.objects = Catalog(obj for obj in self.objects if obj.band == band)
        if exclude_snobjids is not None: self.objects = Catalog(obj for obj in self.objects if obj.snobjid not in exclude_snobjids)
# Exclude g-band to speed it up; they are lower S/N anyway.
        self.objects = Catalog(obj for obj in self.objects if obj.band in ['r','i'])
#
        if date_start is not None: self.objects = Catalog(obj for obj in self.objects if obj.date>=date_start)
        self.objects.add_constant('obscode', 807)
        self.objects.add_constant('err', 0.2)
        #self.objects.refactor('ra', hours)
        #self.objects.refactor('dec', degrees)

        self.objects.orderby('date')
        self.field = fields[field]
        self.band = band
        
    def sep(self, point1, point2):
        # returns separation between two points in arcsec
        return ephem.separation((point1.ra, point1.dec), (point2.ra, point2.dec))*180/np.pi*3600
        
    def link(self, snobjid):

        # get the object associated with this snobjid
        point = self.objects[snobjid]
        
        # link object and return
        next_snobjid = link_obj(point)
        return next_snobjid
       
    def link_obj(self, point):
        visits = self.field.visits
        visits.orderby('nite')
        thisvisit = visits[get_nite(point.date)]
        if thisvisit is None:
            print 'This point\'s visit not found!'
            return []
        nextvisit = thisvisit
        next_obj = []
        for i in range(3):
            
            # get the next visit
            try:
                nextvisit = next(visit for visit in visits if visit.nite > nextvisit.nite)
            except StopIteration: pass
            if nextvisit is None: 
                print "nextvisit is None!"
                return []
            
            current_objects = [obj for obj in self.objects if get_nite(obj.date) == get_nite(nextvisit.date) \
                               and (self.band is None or obj.band == self.band) ]
            # now the real work: for each object, test to see
            # if it's consistent with being the next point in
            # a KBO trajectory.
            lon, lat = Ecliptic(Equatorial(point.ra, point.dec)).get()
            centerlat = Ecliptic(field.center).lat
            for obj in current_objects:
                objlon, objlat = Ecliptic(Equatorial(obj.ra, obj.dec)).get()
                dlon, dlat = nextvisit.dlon - thisvisit.dlon, nextvisit.dlat - thisvisit.dlat
                displacement = ephem.separation((objlon, objlat), (lon, lat))
                norm = np.sqrt(np.cos(centerlat)**2*dlon**2 + dlat**2)
                dot = np.cos(centerlat)**2*(objlon - lon)*dlon + (objlat - lat)*dlat
                if obj.date != point.date:
                 velocity = displacement/(obj.date - point.date)
                else: 
                 break
                cosine = dot/(norm*displacement)
                if cosine > np.cos(20 * np.pi/180) and velocity < 150 * np.pi/648000:
                    next_obj.append(obj)
        return next_obj
