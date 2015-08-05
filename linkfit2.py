#!/usr/bin/env python

from __future__ import division
import numpy as np
import ephem
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
from KBO import *
from Orbit import Orbit
from ObjectLinker import ObjectLinker
from DESKBO import DEScands

force_link=False

def main():
    cand_file = sys.argv[-1]
    if force_link is True:
        force_point_cat = Catalog('../force_link_point.csv', date=DateTime, ra=hours, dec=degrees, expnum=int, exptime=float, \
                                  band=str, ccdnum=int,mag=float,ml_score=float,snobjid=int,fwhm=float,t_eff=float)
        print len(force_point_cat)
        print 'Forcing link with point: '
        for p in force_point_cat:
            print p.date, p.ra, p.dec, p.expnum, p.ccdnum, p.band, p.mag
    cur = 0
    band = None
    while cur < len(sys.argv):
        arg = sys.argv[cur]
        if arg in ['S1','S2','X1','X2','X3','E1','E2','C1','C2','C3','F1','F2','F3','F4','F5','F6','F7']:
            field = fields[arg]
            print "Found field", field.name, "with", len(field.visits), "visits"
        if arg in ['g','r','i','z']:
            band = arg
        cur += 1
    excluded_snobjids = []
    for rock in DEScands():
        for obs in rock.observations:
            if obs.snobjid is not None and obs.snobjid>0:
                excluded_snobjids.append(obs.snobjid)
    linker = ObjectLinker(field.name, band, filename=cand_file, exclude_snobjids=excluded_snobjids, date_start=ephem.date('2013/8/19'))
    print "Read", len(linker.objects), "candidate points from file", cand_file
    outfilename = field.name + '.cands'
    print "Writing to output file", outfilename + '...'
    with open(outfilename, 'wb') as outfile:
        outfile.write("") # blank the file
    outfile2name=field.name+'_good_candidates.txt'
    with open(outfile2name,'wb') as outfile2:
        outfile2.write("")
    i = 0
    pos1=[]
    pos2=[]
    pos3=[]
    plt.figure(1,facecolor='w', edgecolor='w', figsize=(8,8))
    if force_link is True:
        cat1 = force_point_cat
    else:
        cat1 = linker.objects
    for obj1 in cat1:
#    for obj1 in force_point_cat:
        print "Linking object", i, "of", len(linker.objects), "from", pretty_nite(get_nite(obj1.date)), "(snobjid", str(obj1.snobjid) + ")",
        next_points = linker.link_obj(obj1)
#        print 'Object 1 has possible links to ', len(next_points), ' points.'
        for obj2 in next_points:
            print ".",
            next_next_points = linker.link_obj(obj2)
#            print 'Object 2 has possible links to ', len(next_next_points), ' points.'
            for obj3 in next_next_points:
                print '-',
                triple = Catalog([obj1, obj2, obj3])
                ra = [obj1.ra, obj2.ra, obj3.ra]
                dec=[obj1.dec, obj2.dec, obj3.dec]
                dates = [exposure_midpoint(obj1, field), exposure_midpoint(obj2, field), exposure_midpoint(obj3, field)]
                print ra, dec, dates
#                print zip(ra, dec, dates, [obj1.mag, obj2.mag, obj3.mag])
                if np.std([obj1.mag, obj2.mag, obj3.mag])<0.3:
                    orbit = Orbit(ra=ra, dec=dec, dates=dates, obscode=[807,807,807],obserr=0.2)
    #                if orbit.elements['a']>5: 
                    if orbit.elements['a']>18 and orbit.elements['i']<80 and orbit.chisq<2 and np.std([obj1.mag,obj2.mag,obj3.mag])<0.3:
                        if obj1.snobjid not in pos1: pos1.append(obj1.snobjid)
                        if obj2.snobjid not in pos2: pos2.append(obj2.snobjid)
                        if obj3.snobjid not in pos3: pos3.append(obj3.snobjid)
                        if obj1.snobjid in pos2 or obj1.snobjid in pos3 or obj2.snobjid in pos3:
                            with open(outfile2name,'ab') as of2:
                                of.write(str(ephem.now()))
                                for name in ['a','e','i']:
                                    element, error = orbit.elements[name], orbit.elements_errs[name]
                                    of2.write((name + ': ' + str(element) + u' \u00b1 ' + str(error) + '\n').encode('utf8'))
                                of2.write('-'*64 + '\n')
                                of2.write(('Chi^2, ndof = '+str(orbit.chisq)+' '+str(orbit.ndof)+'\n').encode('utf8'))
                                of2.write(triple.writes() + '\n')
                                of2.write('='*64 + '\n')
                        with open(outfilename, 'ab') as outfile:
                            perihelion = orbit.perihelion()
                            peri, perierr = perihelion
                            a, aerr = orbit.elements['a'], orbit.elements_errs['a']
                            inc, incerr = orbit.elements['i'], orbit.elements_errs['i']
                            plt.scatter(np.log10(a), inc, alpha=0.7,color='b',s=6)
                            outfile.write(str(ephem.now()))
                            outfile.write((": " + str(a) + u' \u00b1 ' + str(aerr) + '\n').encode('utf8'))
                            outfile.write((": " + str(inc) + u' \u00b1 ' + str(incerr) + '\n').encode('utf8'))
                            for name in ['a','e','i']:
                                element, error = orbit.elements[name], orbit.elements_errs[name]
                                outfile.write((name + ': ' + str(element) + u' \u00b1 ' + str(error) + '\n').encode('utf8'))
                            outfile.write('-'*64 + '\n')
                            outfile.write(('Chi^2, ndof = '+str(orbit.chisq)+' '+str(orbit.ndof)+'\n').encode('utf8'))
                            outfile.write(triple.writes() + '\n')
                            outfile.write('='*64 + '\n')
                        
        print
        i += 1
#    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#    rc('text', usetex=True)
    plt.xlim([1,2.8])
    plt.ylim([0,180])
    plt.xlabel('log$_{10}$(semi-major axis) (AU)', fontsize=16)
    plt.ylabel('Inclination (deg.)', fontsize=16)
    plt.title('Y2 data: '+field.name+' supernova field', fontsize=20)
    plt.show()
        
if __name__ == '__main__':
    main()
