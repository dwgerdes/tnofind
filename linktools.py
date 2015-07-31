from __future__ import division
import numpy as np
from KBO import *

def genNoise(exposure, nfakes=100):
    '''
    Returns a Catalog of fakes in the specified DECamExposure. 
    '''
    cat = Catalog(date=float,expnum=int,exptime=float, 
                   ra=float,dec=float,band=str,ccd=int,mag=float,nite=int,objid=int, fakeid=int)
    
    delta = ephem.degrees('1.5')
    points = np.random.uniform(-delta, delta, size=(nfakes*10,2))
    Nfiducial=0
    for i in range(len(points)):
        fake_ra = exposure.ra + points[i][0]
        fake_dec = exposure.dec + points[i][1]
        ccdName, ccdNum = compute_chip(fake_ra, fake_dec, exposure.ra, exposure.dec)
        if ccdNum>-99 and Nfiducial<nfakes:
            cat.append(ra=fake_ra, dec=fake_dec, date=exposure.date, expnum=exposure.expnum,
                      exptime=exposure.exptime, band=exposure.band, ccd=ccdNum, mag=20, objid=exposure.expnum*10000+i,
                       nite=exposure.nite, fakeid=-1)
            Nfiducial+=1
#    cat.refactor('date',toDateTime)
#    cat.refactor('ra',  lambda ra: hours(ra))
#    cat.refactor('dec', lambda dec: degrees(dec))
    return cat

    # returns a list of the fields that a given point is in.
def get_field(point):
    in_fields = []
    for field in fields:
#        print field.name, point.ra, point.dec, field.center.ra, field.center.dec
        ccdName, ccdNum = compute_chip(point.ra, point.dec, field.center.ra, field.center.dec) 
        if ccdNum>-99: in_fields.append(field.name)
    return in_fields

# Utility to return the list of nites a given fake_id was observed in a given field
def get_nites(fake_cat, fake_id, field):
    this_fake = [p for p in fake_cat if p.fakeid==fake_id and field in get_field(p)]    # observations matching this fake
    nites = sorted(set([p.nite for p in this_fake]))
    return nites

# defines a "detectable fake"
def isDetectable(fake_cat, fake_id, field):
    return True if len(get_nites(fake_cat, fake_id, field))>2 else False

# Gets the nites associated with each exposure if not present already in our input catalog.
# Requires a database connection.
def query_nites(cat_in):
    nites={}
    expnites={}
    exps = sorted(set([p.expnum for p in cat_in]))
    cursor=db_connect()
    for exp in exps:
        query = 'select e.nite from exposure e where e.expnum = '+str(exp)
        cursor.execute(query)
        rows = cursor.fetchall()
        expnites[exp]=int(rows[0][0])
    for p in cat_in:
        nites[p]=expnites[p.expnum]
    return nites

def t_opp(ra, dec):
    '''
    Computes the approximate time of opposition for the given (ra, dec)
    '''
    t0 = ephem.date('2015/01/01')
    lon, lat = Ecliptic(Equatorial(ra, dec)).get()
    sun = ephem.Sun()
    sun.compute(t0)
    lon_s, lat_s = Ecliptic(Equatorial(sun.ra, sun.dec)).get()
    delta_lon = lon-lon_s
    # first approximation
    t180 = t0 + sidereal_year/2 + delta_lon/sidereal_rate
    # correction
    sun.compute(t180)
    lon_s, lat_s = Ecliptic(Equatorial(sun.ra, sun.dec)).get()
    t180 = t180 + (lon-lon_s + np.pi)/sidereal_rate
    return ephem.date(t180)

def parallax(ra, dec, date):
    '''
    Modeled on anomalies() function in KBO.py.
    
    Calculate the shifts and rate of change in ecliptic coordinates expected for a stationary object 
    at the specified location on the specified date.
    
    date should be a DateTime object, pyEphem
    date, floating-point JD, or compatible string.
    Results are returned as a pair (dlon, dlat),
    in arbitrary units. To get appropriately scaled
    results, multiply by the reciprocal of the object's
    distance from the sun in AU.
    '''
    date = DateTime(date)
    omegat = sidereal_rate * (date - t_opp(ra, dec))
    lat = Ecliptic(Equatorial(ra, dec)).lat
    dlon = -np.sin(omegat)/np.cos(lat)
    dlat = np.cos(omegat)*np.sin(lat)
    vlon = -np.cos(omegat)/np.cos(lat)
    vlat = -np.sin(omegat)*np.sin(lat)
    return dlon, dlat, vlon, vlat

def exposure_parallax():
    p={}
    for exp in exposures:
        dlon, dlat, vlon, vlat = parallax(exp.ra, exp.dec, exp.date)
        p[exp.expnum]={'dlon':dlon, 'dlat':dlat, 'vlon':vlon, 'vlat':vlat}
    return p