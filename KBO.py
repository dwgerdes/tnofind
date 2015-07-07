from __future__ import division
import numpy as np
import ephem
from itertools import count
#from desdb import Connection
from ephem import hours, degrees, Equatorial, Ecliptic
from Orbit import Orbit
from Catalog import Catalog, DateTime, TimeDelta, Point
from DECamField import DECamField
from DECamExposure import DECamExposure
from ccdBounds import *
from MPCRecord import MPCRecord
import easyaccess as ea

sidereal_rate = ephem.degrees(2*np.pi/365.256363)

# adding 4.944 hours puts the sun at exactly 180 degrees ecliptic longitude
autumnal_equinox = DateTime(ephem.next_equinox('2014-07-24') + 4.944*ephem.hour)

fields = Catalog('fields.csv', name=str, centerra=float, centerdec=float, opposition=float, visitspath=str, orderedby='name')
fields.add_property('center', lambda pt: Equatorial(hours(pt.centerra), degrees(pt.centerdec)))
fields.add_property('visits', lambda pt: Catalog(pt.visitspath, nite=int, date=float, dlon=float, dlat=float, vlon=float, vlat=float))
fields.refactor('opposition', DateTime)
for field in fields: field.visits.refactor('date', DateTime)

def _exp_contains(self, ra1, dec1): return DECamField(self.ra, self.dec).contains(ra1, dec1)
def _exp_ellipse(self): return DECamField(self.ra, self.dec).ellipse()

exposures = Catalog('exposures.csv', expnum=int, date=DateTime, ra=hours, dec=degrees, exptime=float, band=str, tag=str, object=str)
exposures.add_function('contains', _exp_contains)
exposures.add_function('ellipse', _exp_contains)
expquality = Catalog('exposure_quality.csv',expnum=int,date=DateTime,band=str,object=str,accepted=str,t_eff=float,fwhm_asec=float,ellipticity=float,skybrightness=float)

def get_nite(date):
    '''
    Get a "nite number" of the form yyyymmdd
    from a DateTime or pyEphem date object.
    '''
    stdate = str(date)
    year, month, daytime = stdate.split('/')
    day, time = daytime.split()
    hour = time.split(':')[0]
    month = month.zfill(2)
    day = day.zfill(2)
    nite = int(year + month + day)
    return nite if int(hour) >= 12 else nite - 1

exposures_by_nite = exposures.groupby(lambda pt: get_nite(pt.date))

def good_visits(field):
    goodvisits = Catalog(nite=int,date=float,band=str,t_eff=float,fwhm_asec=float,ellipticity=float,skybrightness=float)
    for visit in field.visits:
        exps = [e for e in expquality if get_nite(e.date)==visit.nite and field.name in e.object]
        good = True
        for e in exps:
            if e.accepted == 'False': good = False

        

def snob_query(rock, date, rng):
    '''
    Return an SQL query string that looks for SN difference imaging objects near the predicted position of rock on the specified date.
    
    rock should be a pyEphem Orbit object.
    date should be a DateTime, pyEphem date, or
    floating-point Julian Day.
    rng specifies the search range in both ra and dec,
    in (decimal) degrees.
    Fields returned are date_obs, ra, dec, expnum,
    exptime, band, ccdnum, mag, pixelx, pixely,
    snobjid, ml_real_bogus_score.
    '''
    pos = rock.predict_pos(date)
    ra, dec = pos['ra'] * 180/np.pi, pos['dec'] * 180/np.pi
    nite = get_nite(date)
    query = "select e.date_obs, o.ra, o.dec, e.expnum, e.exptime, o.band, o.ccdnum, o.mag, o.ml_real_bogus_score, o.snobjid " \
    "from snobs_legacy o join exposure e on o.exposureid = e.expnum where e.nite = " + \
    str(nite) + " and o.ra between " + str(ra - rng) + " and " + str(ra + rng) + " and o.dec between " + \
    str(dec - rng) + " and " + str(dec + rng) + " order by e.date_obs"
    return query

def object_query(rock, date, rng):
    '''
    Return an SQL query string that looks for wide survey objects near the predicted position of rock on the specified date.
    
    rock should be a pyEphem Orbit object.
    date should be a DateTime, pyEphem date, or
    floating-point Julian Day.
    rng specifies the search range in both ra and dec,
    in (decimal) degrees.
    Fields returned are date_obs, ra, dec, expnum,
    exptime, band, ccd, mag_psf, xwin_image, ywin_image,
    tag.
    '''
    pos = rock.predict_pos(date)
    ra, dec = pos['ra'] * 180/np.pi, pos['dec'] * 180/np.pi
    nite = get_nite(date)
    query = "select e.date_obs, o.ra, o.dec, e.expnum, e.exptime, i.band, i.ccd, o.mag_psf, o.xwin_image, o.ywin_image, t.tag from " \
    "SVA1_FINALCUT o, image i, exposure e, runtag t where i.exposureid = e.id and o.imageid = i.id and i.run = t.run and " \
    "e.nite = " + str(nite) + " and o.ra between " + str(ra - rng) + " and " + str(ra + rng) + " and o.dec between " + \
    str(dec - rng) + " and " + str(dec + rng) + " order by e.date_obs;"
    return query

def field_query(field, band='i'):
    '''
    Return an SQL query string that looks for the nights when a given field was visited.
    
    This looks for i-band exposures by default.
    Fields returned are nite, date_obs, expnum, object.
    '''
    query = "select distinct e.nite, e.date_obs, e.expnum, e.object from exposure e where e.object like 'DES supernova hex SN-" + \
    field + "%' and e.band='" + band + "' order by e.date_obs"
    return query

def anomalies(field, date):
    '''
    Calculate the shifts in ecliptic coordinates expected for a stationary object in the specified field on the specified date.
    
    date should be a DateTime object, pyEphem
    date, floating-point JD, or compatible string.
    Results are returned as a pair (dlon, dlat),
    in arbitrary units. To get appropriately scaled
    results, multiply by the reciprocal of the object's
    distance from the sun in AU.
    '''
    field = fields[field]
    date = DateTime(date)
    omegat = sidereal_rate * (date - field.opposition)
    lat = Ecliptic(field.center).lat
    dlon = -np.sin(omegat)/np.cos(lat)
    dlat = np.cos(omegat)*np.sin(lat)
    return dlon, dlat

def velocities(field, date):
    '''
    Calculate the expected rates of change of ecliptic coordinates expected for a stationary object in the specified field on the specified date.
    
    date should be a DateTime object, pyEphem
    date, floating-point JD, or compatible string.
    Results are returned as a pair (vlon, vlat),
    in arbitrary units. To get appropriately scaled
    results, multiply by the sidereal rate and then
    by the reciprocal of the object's distance from
    the sun in AU.
    '''
    field = fields[field]
    date = DateTime(date)
    omegat = sidereal_rate * (date - field.opposition)
    lat = Ecliptic(field.center).lat
    vlon = -np.cos(omegat)/np.cos(lat)
    vlat = -np.sin(omegat)*np.sin(lat)
    return vlon, vlat

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
        if exp.contains(target_ra, target_dec):
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
        if exp.contains(target_ra, target_dec) and (snfields or not exp.object.startswith('DES supernova hex')):
            ccdname, ccdnum = compute_chip(target_ra, target_dec, exp.ra, exp.dec)
            this_match = {'expnum': exp.expnum, 'ccd': ccdnum, 'date': exp.date, 'band': exp.band}
            if ccdnum != -99: match.append(this_match)
    return match

def count(start = 0):
    i = start
    while True:
        i += 1
        yield i
 
    
    #def __init__(self, obsnum='     ', MPprovisional='       ', discovery=' ', note1=' ', 
    #             note2='C', obsdate=ephem.date('2000/01/01'), ra_obs_J2000=ephem.hours(0), dec_obs_J2000=ephem.degrees(0), 
    #             mag=99, band='r', observatoryCode='W84', newobject=True):       
def MPCobservation(point, temp_designation='       ', packed_designation='     '):
    newobject=False
    if packed_designation=='       ': newobject=True
    if newobject==False and packed_designation=='      ':
        print 'MPCobservation Error, must supply packed designation'
    rec = MPCRecord(obsnum=packed_designation, MPprovisional=temp_designation, obsdate=ephem.date(point.date+point.exptime*ephem.second/2), ra_obs_J2000=point.ra, dec_obs_J2000=point.dec,
                    mag=point.mag, band=point.band, newobject=newobject)
    return rec.record

def absolute_magnitude(orbit, apparent_mag, date_obs):
    # computes the absolute magnitude of an object, given its orbit and apparent magnitude on date_obs
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
    if field.name in ['X3', 'C3']:
        if obs.band in ['g','r']:
            nstack = 3
        elif obs.band == 'i':
            nstack = 5
        elif obs.band == 'z':
            nstack = 11
    elif field.name in ['E1', 'E2', 'C1', 'C2', 'S1', 'S2', 'X1','X2'] and obs.band=='z':
        nstack = 2
    return ephem.date(obs.date + nstack*ephem.second*obs.exptime/2)


