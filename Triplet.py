
import pandas as pd 
import ephem
import numpy as np
from linkutils import *

def getRow(ind, df, objid_col='snobjid'):
    '''
    ind is the objid or equivalent identifier
    df is the pandas dataframe
    This is a utility routine to return the row of the dataframe containing this object id. 
    '''
    return df.loc[df[objid_col]==ind]

class Triplet(object):
    def __init__(self, triple, vlist, df, classifier=None, objid_col='snobjid'):
        '''
        Container class for triplet info.
        triple is a list of three objids
        '''
        self.ids = triple
        self.vlon12, self.vlat12 = vlist[triple[0]][triple[1]]
        self.vlon23, self.vlat23 = vlist[triple[1]][triple[2]]
        self.dfRows = dict(zip(self.ids,[getRow(t, df, objid_col=objid_col) for t in self.ids]))
        self.ra = dict(zip(self.ids, [self.dfRows[t]['ra'].values[0]*np.pi/180 for t in self.ids]))
        self.dec = dict(zip(self.ids, [self.dfRows[t]['dec'].values[0]*np.pi/180 for t in self.ids]))
        self.date = dict(zip(self.ids, [self.dfRows[t]['date_obs'].values[0] for t in self.ids]))
        self.cos12 = self.reflex(self.ids[0], self.ids[1])
        self.cos23 = self.reflex(self.ids[1], self.ids[2])
        self.chisq_cut = 5
#        self.chisq, self.ndof = self.checkTriplet()
        self.clf = classifier   # ML classifier (from scikit-learn)
        
    def dvlon(self):
        return 2*(self.vlon12-self.vlon23)/(self.vlon12+self.vlon23)

    def dvlat(self):
        return 2*(self.vlat12-self.vlat23)/(self.vlat12+self.vlat23)

    def dcos(self):
        return self.cos12 - self.cos23
    
    def ML_data(self):
        '''
        variables for triplet classification
        '''
        return np.array(np.array([self.dvlon(), self.dcos()]))
        
    def checkTriplet(self):
        '''
        Full orbit fit
        '''
        ralist = [ephem.hours(self.ra[self.ids[i]]) for i in range(3)]
        declist = [ephem.degrees(self.dec[self.ids[i]]) for i in range(3)]
        datelist = [ephem.date(self.date[self.ids[i]]) for i in range(3)]
        orbit = Orbit(dates=datelist, ra=ralist, dec=declist, obscode=[807, 807, 807], err=0.15)
        return orbit.chisq, orbit.ndof
    
    def checkTriplet2(self):
        if self.clf is not None:
            return self.clf.predict(self.ML_data())[0]
        else:
            return -1
        
    def ML_prob(self):
        if self.clf is not None:
            return self.clf.predict_proba(self.ML_data())[0]
    
    def checkTriplet3(self):
        '''
        Dummy routine to compare timing with other triplet validation methods
        '''
        return 1, 1           # placeholder to compare timing with orbit fitting
    
    def reflex(self, ind1, ind2, debug=False):
        date1, ra1, dec1 = ephem.date(self.date[ind1]), ephem.hours(self.ra[ind1]), ephem.degrees(self.dec[ind1])
        date2, ra2, dec2 = ephem.date(self.date[ind2]), ephem.hours(self.ra[ind2]), ephem.degrees(self.dec[ind2])
        lon1, lat1 = ephem.Ecliptic(ephem.Equatorial(ra1, dec1)).get()
        lon2, lat2 = ephem.Ecliptic(ephem.Equatorial(ra2, dec2)).get()
        displacement = ephem.separation((lon2, lat2), (lon1, lat1))
        displacement_asec = displacement*180/np.pi*3600
        if date2 != date1:
            velocity = displacement_asec/(date2 - date1)
        else: 
            velocity=9999
        this_dlon, this_dlat, this_vlon, this_vlat = parallax(ra1, dec1, date1) 
        next_dlon, next_dlat, next_vlon, next_vlat = parallax(ra2, dec2, date2)   
        dlon, dlat = next_dlon-this_dlon, next_dlat-this_dlat
        dot = np.cos(lat1)**2*(lon2 - lon1)*dlon + (lat2 - lat1)*dlat
        norm = np.sqrt(np.cos(lat1)**2*dlon**2 + dlat**2)
        cosine = dot/(norm*displacement)
        return cosine
    
    
    def isGood(self):
        return True if self.ndof>0 and self.chisq/self.ndof<self.chisq_cut else False
    
    def isPure(self, fakeid_col='snfake_id'):
        # Does triplet all come from the same fakeid?
        fakeids = [self.dfRows[t][fakeid_col].values[0] for t in self.ids]
        return True if fakeids[0]==fakeids[1]==fakeids[2] and fakeids[0]>0 else False