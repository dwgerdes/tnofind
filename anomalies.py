from __future__ import division
import numpy as np

'''
Utility routines to convert among mean, true, and eccentric anomalies. 
'''

def MeanToEccenAnomaly(e, M):
    '''
    Convert mean anomaly to eccentric anomaly. All anomalies are in RADIANS.
    Adapted from Octave routine written by Luis Baars.
    '''
    if e>=1:
        print 'MeanToEccenAnomaly does not support parabolic or hyperbolic orbits!'
        return -999
    else:
## Make sure M lies between -pi and pi
        M = np.mod(M, 2*np.pi)
        if M>np.pi:
            M -= 2*np.pi

        if (-np.pi < M < 0) or (M > np.pi):
            E = M - e
        else:
            E = M + e

        Enew = E
        zero = 1e-06
        first = True
        while (first or np.abs(Enew-E)>zero):
            first = False
            E = Enew
            Enew = E + (M - E + e*np.sin(E))/(1 - e*np.cos(E))
        E = Enew
        return E

def EccenToTrueAnomaly(e, E):
    '''
    Convert eccentric anomaly to true anomaly. All anomalies are in RADIANS.
    Adapted from Octave routine written by Luis Baars.
    '''     
    if e>=1:
        print 'EccenToTrueAnomaly does not support parabolic or hyperbolic orbits!'
        f = -999
    else:
# Calc the true anomaly
        sinf = np.sin(E)*np.sqrt(1 - e**2)/(1 - e * np.cos(E));
        cosf = (np.cos(E) - e)/(1 - e * np.cos(E));
        f = np.arctan2(sinf, cosf);
    return f

def MeanToTrueAnomaly(e, M):
    '''
    Convert mean anomaly to true anomaly. All anomalies are in RADIANS.
    Adapted from Octave routine written by Luis Baars.
    '''   
    if e>=1:
        print 'MeanToTrueAnomaly does not support parabolic or hyperbolic orbits!'
        f = -999
    else:
        E = MeanToEccenAnomaly(e, M)
        f = EccenToTrueAnomaly(e, E)
    return f

def TrueToEccenAnomaly(e, f):
    '''
    Convert true anomaly to eccentric anomaly. All anomalies are in RADIANS.
    Adapted from Octave routine written by Luis Baars.
    '''   
    if e>=1:
        print 'TrueToEccenAnomaly does not support parabolic or hyperbolic orbits!'
        E = -999
    else:
## Calc the eccentric anomaly
        sinE = np.sin(f)*np.sqrt(1 - e**2)/(1 + e * np.cos(f))
        cosE = (e + np.cos(f))/(1 + e * np.cos(f))
        E = np.arctan2(sinE, cosE)
    return E

def EccenToMeanAnomaly(e, E):
    '''
    Convert eccentric anomaly to mean anomaly. All anomalies are in RADIANS.
    Adapted from Octave routine written by Luis Baars.
    '''   
    if e>=1:
        print 'EccenToMeanAnomaly does not support parabolic or hyperbolic orbits!'
        M = -999
    else:
## Calc the eccentric anomaly
        M = E - e * np.sin(E)
    return M


def TrueToMeanAnomaly(e, f):
    '''
    Convert true anomaly to mean anomaly. All anomalies are in RADIANS.
    Adapted from Octave routine written by Luis Baars.
    '''   
    if e>=1:
        print 'TrueToMeanAnomaly does not support parabolic or hyperbolic orbits!'
        M = -999
    else:
        E = TrueToEccenAnomaly(e, f)
        M = EccenToMeanAnomaly(e, E)
    return M
