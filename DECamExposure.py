#!/usr/bin/env python

import os
import ephem
from DECamField import DECamField

class DECamExposure(object):
#
    def __init__(self, expnum=0, UTobs='2013-01-01 00:00:00', exptime=0, band='r', ra=ephem.degrees(0), dec=ephem.degrees(0), nite=20130101, tag='None', obj='None'):
        self.expnum = expnum
        self.UTobs = ephem.date(UTobs)
        self.exptime = exptime
        self.band = band
        self.ra = ephem.degrees(ra)
        self.dec = ephem.degrees(dec)
        self.tag = tag
        self.obj = obj
	self.nite = nite
    
    def contains(self, ra1, dec1):
        # returns True if the point (ra1, dec1) lies inside the field
        return DECamField(self.ra, self.dec).contains(ra1, dec1)
    
    def ellipse(self):
        return DECamField(self.ra, self.dec).ellipse()
    
    def dump(self):
        print 'ExpID: \t', self.expnum
        print 'UTobs: \t', self.UTobs
        print 'Exptime: \t', self.exptime
        print 'Band: \t', self.band
        print 'RA: \t', self.ra
        print 'DEC: \t', self.dec
        print 'Tag: \t', self.tag
        print 'Tile: \t', self.obj
        
    def local_files(self, rootdir):
        # Searches rootdir and its subdirectories for files (not directories) of the form DECam_nnnnnnnn_cc.* where nnnnnnnn is the expnum
        a = os.walk(rootdir)
        flist = []
        for root, dirs, files in a:
            for f in files:
                if str(self.expnum) in f and 'DECam_' in f:
                    flist.append(os.path.join(root, f))
        return flist
    
    def local_nulls(self, rootdir):
        # Searches rootdir and its subdirectories for a directory containing 'null_nnnnnnnn' where nnnnnnn is the expnum,
        # and makes a list of the files it contains
        a = os.walk(rootdir)
        nlist = []
        for root, dirs, files in a:
            for d in dirs:
                if 'null_'+str(self.expnum) in d:
                    for r, dirs2, files2 in os.walk(os.path.join(root, d)):
                        for f in files2:
                            nlist.append(os.path.join(r,f))
        return nlist
    
def main():
    pass

if __name__=="__main__":
    main()
