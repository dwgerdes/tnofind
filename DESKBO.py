#!/usr/bin/env python

from __future__ import division

import ephem
import numpy as np
from Orbit import Orbit
from KBO import Catalog, DateTime, hours, degrees
import webcolors

GM = 4.*np.pi*np.pi/1.0000378 # solar gravitation, from orbfit
combinedMass = GM * 1.00134  #Alter GM to account for total SS mass, from orbfit
DAY = 1/365.25

class DESKBO(object):
    def __init__(self, obsfile, name=None, field=None, pltcolor='k'):
#        print 'Reading file ', obsfile
        self.name=name
        self.field=field
        self.obsfile=obsfile
        self.observations = Catalog(self.obsfile, date=DateTime, ra=hours, dec=degrees, expnum=int, exptime=float, band=str, ccdnum=int, mag=float, \
                                    ml_score=float, snobjid=int, orderedby='date')
        self.observations.add_constant('obscode', 807)
        self.observations.add_constant('err', 0.15)
        self.orbit = Orbit(self.observations)   # Warning: somehow orbit.predict_pos() does not work from here.
        self.body = self.orbit.ellipticalBody(name=self.name)
        self.pltcolor=pltcolor
        self.elements, self.errs = self.orbit.get_elements()
        obsdates = sorted([o.date for o in self.observations])
        self.body.compute(obsdates[0])
        self.discoveryDistance = self.body.sun_distance
        self.size = KBOsize(self)
        self.epoch = self.orbit.jd0
        self.mu = (self.epoch-self.elements['top'])*np.sqrt(combinedMass/self.elements['a']**3)*DAY*180/np.pi
        if self.mu<0: self.mu += 360.0
        
    def predict(self, d, obscode=807):
        # d is an ephem.Date object
        self.body.compute(d)
        pos = {'ra':self.body.a_ra, 'dec':self.body.a_dec}
        return pos
    
def DEScands(fields=None):

    prefix = '/Users/gerdes/TNO/'
    cands = [
        {'obsfile':prefix+'cands/VR113/VR113.csv', 'name':'2012 VR113', 'field':'X2', 'pltcol':np.array(webcolors.name_to_rgb('indianred'))/255},    # 1
        {'obsfile':prefix+'cands/VS113/VS113.csv', 'name':'2012 VS113', 'field':'X3', 'pltcol':np.array(webcolors.name_to_rgb('coral'))/255},        # 2
        {'obsfile':prefix+'cands/VU113/VU113.csv', 'name':'2012 VU113', 'field':'S2', 'pltcol':np.array(webcolors.name_to_rgb('indigo'))/255},       # 3
        {'obsfile':prefix+'cands/VV113/VV113.csv', 'name':'2012 VV113', 'field':'S2', 'pltcol':np.array(webcolors.name_to_rgb('green'))/255},        # 4
        {'obsfile':prefix+'cands/WD36/WD36.csv', 'name':'(451657) 2012 WD36', 'field':'S2', 'pltcol':np.array(webcolors.name_to_rgb('crimson'))/255},         # 5
        {'obsfile':prefix+'cands/YO9/YO9.csv', 'name':'2012 YO9','field':'S1','pltcol':np.array(webcolors.name_to_rgb('limegreen'))/255},        # 23
        {'obsfile':prefix+'cands/QO95/QO95.csv', 'name':'2013 QO95','field':'X1','pltcol':np.array(webcolors.name_to_rgb('royalblue'))/255},         # 6
        {'obsfile':prefix+'cands/QP95/QP95_matches.csv', 'name':'2013 QP95','field':'X1','pltcol':np.array(webcolors.name_to_rgb('magenta'))/255},           # 7
        {'obsfile':prefix+'cands/RB98/RB98.csv', 'name':'2013 RB98','field':'X2','pltcol':np.array(webcolors.name_to_rgb('blue'))/255},              # 8
        {'obsfile':prefix+'cands/RD98/RD98.csv', 'name':'2013 RD98','field':'X1','pltcol':np.array(webcolors.name_to_rgb('deeppink'))/255},          # 9
        {'obsfile':prefix+'cands/RF98/RF98.csv', 'name':'2013 RF98','field':'X3','pltcol':np.array(webcolors.name_to_rgb('khaki'))/255},             # 20
        {'obsfile':prefix+'cands/RG98/RG98.csv', 'name':'2013 RG98','field':'C3','pltcol':np.array(webcolors.name_to_rgb('cornflowerblue'))/255},    # 22
        {'obsfile':prefix+'cands/RM98/RM98.csv', 'name':'2013 RM98','field':'Stripe82','pltcol':np.array(webcolors.name_to_rgb('darkcyan'))/255},         #38
        {'obsfile':prefix+'cands/SE99/SE99.csv', 'name':'2013 SE99','field':'X3','pltcol':np.array(webcolors.name_to_rgb('crimson'))/255},          # 11
        {'obsfile':prefix+'cands/TH159/TH159.csv', 'name':'2013 TH159','field':'X3','pltcol':np.array(webcolors.name_to_rgb('lavender'))/255},         #24
        {'obsfile':prefix+'cands/TJ159/TJ159.csv', 'name':'TJ159','field':'Stripe82','pltcol':np.array(webcolors.name_to_rgb('cornflowerblue'))/255},         #39
        {'obsfile':prefix+'cands/TV158/TV158_matches.csv', 'name':'(437360) 2013 TV158','field':'X1','pltcol':np.array(webcolors.name_to_rgb('green'))/255},          # 10
        {'obsfile':prefix+'cands/VD24/VD24.csv', 'name':'2013 VD24','field':'X3','pltcol':np.array(webcolors.name_to_rgb('mediumpurple'))/255},      #27
        {'obsfile':prefix+'cands/QL441/QL441.csv', 'name':'2014 QL441','field':'X2','pltcol':np.array(webcolors.name_to_rgb('lightslategray'))/255}, # 12
        {'obsfile':prefix+'cands/QM441/QM441.csv', 'name':'2014 QM441','field':'S2','pltcol':np.array(webcolors.name_to_rgb('chocolate'))/255},      # 13
        {'obsfile':prefix+'cands/QN441/QN441.csv', 'name':'2014 QN441','field':'X2','pltcol':np.array(webcolors.name_to_rgb('limegreen'))/255},      # 14
        {'obsfile':prefix+'cands/QO441/QO441.csv', 'name':'2014 QO441','field':'X3','pltcol':np.array(webcolors.name_to_rgb('aquamarine'))/255},     # 15
        {'obsfile':prefix+'cands/QP441/QP441.csv', 'name':'2014 QP441','field':'X3','pltcol':np.array(webcolors.name_to_rgb('darkcyan'))/255},     # 16
        {'obsfile':prefix+'cands/QR441/QR441.csv', 'name':'2014 QR441','field':'E1','pltcol':np.array(webcolors.name_to_rgb('indianred'))/255},         #26
        {'obsfile':prefix+'cands/QS441/QS441.csv', 'name':'2014 QS441','field':'X1','pltcol':np.array(webcolors.name_to_rgb('indianred'))/255},         #30
        {'obsfile':prefix+'cands/QU441/QU441.csv', 'name':'2014 QU441','field':'S1','pltcol':np.array(webcolors.name_to_rgb('indianred'))/255},         #36
        {'obsfile':prefix+'cands/SB349/SB349.csv', 'name':'2014 SB349','field':'S2','pltcol':np.array(webcolors.name_to_rgb('indianred'))/255},         #36
        {'obsfile':prefix+'cands/SZ348/SZ348.csv', 'name':'2014 SZ348','field':'C3','pltcol':np.array(webcolors.name_to_rgb('blue'))/255},            # 18
        {'obsfile':prefix+'cands/TT85/TT85.csv', 'name':'2014 TT85','field':'S2','pltcol':np.array(webcolors.name_to_rgb('indianred'))/255},         # 17
        {'obsfile':prefix+'cands/TU85/TU85.csv', 'name':'2014 TU85','field':'X1','pltcol':np.array(webcolors.name_to_rgb('cornflowerblue'))/255},    # 21
        {'obsfile':prefix+'cands/UF224/UF224.csv', 'name':'2014 UF224','field':'X1','pltcol':np.array(webcolors.name_to_rgb('chocolate'))/255},       # 19
        {'obsfile':prefix+'cands/VT37/VT37.csv', 'name':'2014 VT37','field':'X3','pltcol':np.array(webcolors.name_to_rgb('mediumorchid'))/255},    #25
        {'obsfile':prefix+'cands/PD312/PD312.csv', 'name':'2015 PD312','field':'X2','pltcol':np.array(webcolors.name_to_rgb('indianred'))/255},         #30
        {'obsfile':prefix+'cands/PF312/PF312.csv', 'name':'2015 PF312','field':'S2','pltcol':np.array(webcolors.name_to_rgb('indianred'))/255},         #29
        {'obsfile':prefix+'cands/PK312/PK312.csv', 'name':'2015 PK312','field':'X3','pltcol':np.array(webcolors.name_to_rgb('indianred'))/255},          #34
        {'obsfile':prefix+'cands/QT11/QT11.csv', 'name':'2015 QT11','field':'X3','pltcol':np.array(webcolors.name_to_rgb('indianred'))/255}         #39
#        {'obsfile':prefix+'cands/Y3/test.csv', 'name':'2015 test','field':'X2','pltcol':np.array(webcolors.name_to_rgb('darkcyan'))/255},         #38
        ]
    rocks = [DESKBO(c['obsfile'], name=c['name'], field=c['field'],pltcolor=c['pltcol']) for c in cands]
    if fields is not None:
        rocks = [r for r in rocks if r.field in fields]
    return rocks

def known_snobjids():
    # Returns a list of snobjids associated with known objects. Useful if you want to exclude already-known objects from a linking pass, e.g.
    known = []
    for rock in DEScands():
        for obs in rock.observations:
            if obs.snobjid is not None and obs.snobjid>0:
                known.append(obs.snobjid)
    return known

def KBOmag(rock, band='r', exclude_nites=None):
    all_nites = set([get_nite(o.date) for o in rock.observations if o.exptime>90])
    if exclude_nites is not None:
        nites = sorted([n for n in all_nites if n not in exclude_nites])
    else:
        nites = sorted([n for n in all_nites])
    mags = [o.mag for o in rock.observations if o.band==band and o.ml_score>0.4 and o.mag>0]
    return np.average(mags)

def KBOcolor(rock, band1='r', band2='i', exclude_nites=None):
    # Returns the average color for the rock.
    # If exclude_nites is none, all available nites are used.
    # NB: solar colors are:
    #   g-r = 0.44 +/- 0.02
    #   r-i = 0.11 +/- 0.02
    #   i-z = 0.03 +/- 0.02
    all_nites = set([get_nite(o.date) for o in rock.observations if o.exptime>90])
    if exclude_nites is not None:
        nites = sorted([n for n in all_nites if n not in exclude_nites])
    else:
        nites = sorted([n for n in all_nites])
    cols = []
    for n in nites:
        exps = [o for o in rock.observations if get_nite(o.date)==n]
        mags_1 = np.array([e.mag for e in exps if e.band==band1])
        mags_2 = np.array([e.mag for e in exps if e.band==band2])
        if len(mags_1)>0:
            mag_1 = np.average(mags_1)
        else: mag_1 = 0
        if len(mags_2)>0:
            mag_2 = np.average(mags_2)
        else: mag_2 = 0
        col=-99
        if mag_1>0 and mag_2>0: col = mag_1-mag_2
        if col != -99:
            cols.append(col)
#            print n, col
    return cols, np.average(cols), np.std(cols)/np.sqrt(len(cols))

def KBOsize(rock, albedo=0.05):
    mag_r_sun = -27.1
    phase = 1
    size=[]
    for obs in rock.observations:
        if obs.band == 'r' and obs.mag>0 and obs.ml_score>0.4:
            mag_r = obs.mag
            rock.body.compute(obs.date)
            d_earth = rock.body.earth_distance
            d_sun = rock.body.sun_distance
            radius = np.sqrt(2.25e16*d_sun**2*d_earth**2/(albedo*phase))*10**(0.2*(mag_r_sun-23.9))
            size.append(2*radius)
    return np.average(size)
    
    
    
def main():
    rock = DESKBO('../cands/QO441/QO441.csv', name='DES KBO')
    cur_pos = rock.predict(ephem.now())
    print cur_pos['ra'], cur_pos['dec']
    print 'Mag g: ', KBOmag(rock, band='g', exclude_nites=[])
    print 'Mag r: ', KBOmag(rock, band='r', exclude_nites=[])
    print 'Mag i: ', KBOmag(rock, band='i', exclude_nites=[])
    print 'Mag z: ', KBOmag(rock, band='z', exclude_nites=[])
    print 'Estimated size: ', rock.size

    
if __name__ == '__main__':
    main()
