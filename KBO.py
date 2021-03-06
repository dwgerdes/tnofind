from __future__ import division
import numpy as np
import ephem
from itertools import count
from ephem import hours, degrees, Equatorial, Ecliptic
from Orbit import Orbit
from Catalog import Catalog, DateTime, TimeDelta, Point
from DECamField import DECamField
from DECamExposure import DECamExposure
from ccdBounds import *
from MPCRecord import MPCRecord
import easyaccess as ea

sidereal_year = 365.256363
sidereal_rate = ephem.degrees(2*np.pi/sidereal_year)


fields = Catalog('fields.csv', name=str, centerra=float, centerdec=float, opposition=float, visitspath=str, orderedby='name')
fields.add_property('center', lambda pt: Equatorial(hours(pt.centerra), degrees(pt.centerdec)))
#fields.add_property('visits', lambda pt: Catalog(pt.visitspath, nite=int, date=float, dlon=float, dlat=float, vlon=float, vlat=float))
fields.refactor('opposition', DateTime)
#for field in fields: field.visits.refactor('date', DateTime)

def _exp_contains(self, ra1, dec1): return DECamField(self.ra, self.dec).contains(ra1, dec1)
def _exp_ellipse(self): return DECamField(self.ra, self.dec).ellipse()

exposures = Catalog('exposures.csv', expnum=int, date=float, ra=float, dec=float, exptime=float, nite=int, band=str, tag=str, object=str,
    fwhm_asec=float, t_eff=float, ellipticity=float, skybrightness=float, accepted=str, analyst=str, analyst_comment=str, lastchanged_time=str)
exposures.refactor('date', lambda date: ephem.date(date))
exposures.refactor('ra', lambda ra: hours(ra))
exposures.refactor('dec', lambda dec: degrees(dec))
exposures.add_function('contains', _exp_contains)
exposures.add_function('ellipse', _exp_contains)
exposures_by_nite = exposures.groupby(lambda pt: pt.nite)
        

def toDateTime(ISOdate):
    '''
    Transform an ISO date string into a DateTime object.
    '''
    datestring = ' '.join(ISOdate.split('T'))
    return DateTime(datestring)

def pretty_nite(nite):
    '''
    Return a string in yyyy/mm/dd format from a
    "nite number" of the form yyyymmdd.
    '''
    nitestr = str(nite)
    nitestr = "/".join([nitestr[:4], nitestr[4:6], nitestr[6:]])
    return nitestr

def sexagesimal(ra, dec):
    '''
    Return the sexagesimal representation, as pyEphem
    objects, of the coordinates (ra, dec), which are
    taken to be in decimal degree format.
    '''
    return hours(float(ra) * np.pi/180), degrees(float(dec) * np.pi/180)

def decimal(ra, dec):
    '''
    Return the decimal degree representation of the
    coordinates (ra, dec), which are taken to be in
    radians (or pyEphem objects stored that way).
    '''
    return hours(ra) * 180/np.pi, degrees(dec) * 180/np.pi

def days_between(start, stop):
    '''
    Return a generator object that iterates over
    the days between start and stop as DateTime objects.
    '''
    startd, stopd = DateTime(start), DateTime(stop)
    d = startd
    while d < stopd:
        yield d
        d = DateTime(d + 1)

def compute_chip(rockra, rockdec, expra, expdec):
    '''
    Given the ra and dec of a point and of the center
    of an exposure, find the CCD containing that point.
    
    Returns a pair of the CCD name and number.
    '''
    deltara = 180/np.pi*ephem.degrees(rockra-expra).znorm  # compute difference in degrees (normalized between -180, +180)
    deltadec = 180/np.pi*ephem.degrees(rockdec-expdec).znorm  # the 180/pi is because ephem.Angle objects are natively in radians
    ccdname = 'None'
    for k in ccdBounds:
        if deltara > ccdBounds[k][0] and deltara < ccdBounds[k][1] and deltadec > ccdBounds[k][2] and deltadec < ccdBounds[k][3]:
            ccdname = k
    return ccdname, ccdNum[ccdname]

def find_exposures(target_ra, target_dec):
    '''
    Find exposures containing a given (ra, dec).
    
    Returns a list of dictionaries with keys 'expnum', 'ccd', 'date', 'band'
    expnum is the exposure number
    ccd is the (integer) chip id
    date is the date as a DateTime object
    band is the band (g,r,i,z).
    '''
    match = []
    for exp in exposures:
        ccdname, ccdnum = compute_chip(target_ra, target_dec, exp.ra, exp.dec)
        this_match = {'expnum': exp.expnum, 'ccd': ccdnum, 'date': exp.date, 'band': exp.band}
        if ccdnum != -99: match.append(this_match)
    return match

def find_exposures_by_nite(nite, target_ra, target_dec, snfields = True):
    '''
    Find exposures on a given night containing a given (ra, dec).
    
    Returns a list of dictionaries with keys 'expnum', 'ccd', 'date', 'band'
    expnum is the exposure number
    ccd is the (integer) chip id
    date is the date as a DateTime object
    band is the band (g,r,i,z).
    '''
    match = []
    for exp in exposures_by_nite[nite]:
        if (snfields or not exp.object.startswith('DES supernova hex')):
            ccdname, ccdnum = compute_chip(target_ra, target_dec, exp.ra, exp.dec)
            this_match = {'expnum': exp.expnum, 'ccd': ccdnum, 'date': exp.date, 'band': exp.band}
            if ccdnum != -99: match.append(this_match)
    return match

def count(start = 0):
    i = start
    while True:
        i += 1
        yield i
       
def MPCobservation(point, temp_designation='       ', packed_designation='     '):
    newobject=False
    if packed_designation=='       ': newobject=True
    if newobject==False and packed_designation=='      ':
        print 'MPCobservation Error, must supply packed designation'
    rec = MPCRecord(obsnum=packed_designation, MPprovisional=temp_designation, obsdate=ephem.date(point.date+point.exptime*ephem.second/2), ra_obs_J2000=point.ra, dec_obs_J2000=point.dec,
                    mag=point.mag, band=point.band, newobject=newobject)
    return rec.record

def absolute_magnitude(orbit, apparent_mag, date_obs):
    ''' computes the absolute magnitude of an object, given its orbit and apparent magnitude on date_obs
    '''
    body = orbit.ellipticalBody()
    body.compute(date_obs)
    d_BS = body.sun_distance
    d_BE = body.earth_distance
    d_ES = 1.0
    cos_chi = (d_BE**2 + d_BS**2 - d_ES**2)/(2*d_BE*d_BS)
    chi = np.arccos(cos_chi)
#    P = (2/3)*((1-chi/np.pi)*np.cos(chi) + (1/np.pi)*np.sin(chi))
    P=1
    H = apparent_mag - 2.5*np.log10(d_BS**2*d_BE**2/(P*d_ES**4))
    return H

def db_connect(section='desoper'):
    db = ea.connect(section=section)
    cursor = db.cursor()
    return cursor

def convertDB(row):
# converts a query result to our standard format for observation files
    date_obs = toDateTime(row['DATE_OBS'])
    ra = ephem.hours(ephem.degrees(str(row['RA'])))
    dec = ephem.degrees(str(row['DEC']))
    myrow = {'date':date_obs, 'ra':ra, 'dec':dec, 'expnum':int(row['EXPNUM']), 'exptime':float(row['EXPTIME']), 'band':row['BAND'],
             'ccdnum':int(row['CCDNUM']), 'mag':round(float(row['MAG']),3), 
            'ml_score':float(row['ML_SCORE']), 'snobjid':int(row['SNOBJID'])}
    return myrow


def exposure_midpoint(obs, field):
    # Computes the midpoint of the exposure, accounting for the fact that some bands/fields are stacked. 
    nstack = 1
    if field in ['X3', 'C3']:
        if obs.band in ['g','r']:
            nstack = 3
        elif obs.band == 'i':
            nstack = 5
        elif obs.band == 'z':
            nstack = 11
    elif field in ['E1', 'E2', 'C1', 'C2', 'S1', 'S2', 'X1', 'X2'] and obs.band=='z':
        nstack = 2
    elif field=='wide':
        nstack=1
    return ephem.date(obs.date + nstack*ephem.second*obs.exptime/2)


def wrap_degrees(theta):
    '''
    Express angle as -180<theta<180 rather than 0 < theta < 360
    
    '''
    return theta if 0<theta<180 else theta-360

