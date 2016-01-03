# encoding: utf-8
# cython: profile=False
# filename: linkutils.pyx

from libc.math cimport sin, cos
from KBO import DateTime, exposures
import numpy as np
from ephem import Ecliptic, Equatorial, hours, degrees, date, Sun

cimport cython

cdef double sidereal_year = 365.256363
cdef double sidereal_rate = degrees(2*np.pi/sidereal_year)

@cython.profile(False)
cdef inline double t_opp(double ra, double dec):
    '''
    Computes the approximate time of opposition for the given (ra, dec)
    '''
    cdef double t0, lon, lat, lon_s, lat_s, delta_lon, t180

    t0 = date('2015/01/01')
    lon, lat = Ecliptic(Equatorial(ra, dec)).get()
    sun = Sun()
    sun.compute(t0)
    lon_s, lat_s = Ecliptic(Equatorial(sun.ra, sun.dec)).get()
    delta_lon = lon-lon_s
    # first approximation
    t180 = t0 + sidereal_year/2 + delta_lon/sidereal_rate
    # correction
    sun.compute(t180)
    lon_s, lat_s = Ecliptic(Equatorial(sun.ra, sun.dec)).get()
    t180 = t180 + (lon-lon_s + np.pi)/sidereal_rate
    return date(t180)

def parallax(double ra, double dec, double date):
    '''
    
    Calculate the shifts and rate of change in ecliptic coordinates expected for a stationary object 
    at the specified location on the specified date.
    
    date should be a DateTime object, pyEphem
    date, floating-point JD, or compatible string.
    Results are returned as a pair (dlon, dlat),
    in arbitrary units. To get appropriately scaled
    results, multiply by the reciprocal of the object's
    distance from the sun in AU.
    '''
    cdef double omegat, dlon, dlat, vlon, vlat

    date = DateTime(date)
    omegat = sidereal_rate * (date - t_opp(ra, dec))
    lat = Ecliptic(Equatorial(ra, dec)).lat
    dlon = -sin(omegat)/cos(lat)
    dlat = cos(omegat)*sin(lat)
    vlon = -cos(omegat)/cos(lat)
    vlat = -sin(omegat)*sin(lat)
    return dlon, dlat, vlon, vlat


def exposure_parallax():
    p={}
    for exp in exposures:
        dlon, dlat, vlon, vlat = parallax(exp.ra, exp.dec, exp.date)
        p[exp.expnum]={'dlon':dlon, 'dlat':dlat, 'vlon':vlon, 'vlat':vlat}
    return p

def test():
    a,b,c,d = parallax(hours('3:30:00'),degrees('0'),date('2014/01/01'))
    para = exposure_parallax()
    return a,b,c,d
