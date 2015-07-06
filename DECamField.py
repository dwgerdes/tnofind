#!/usr/bin/env python

import numpy as np
import ephem
from matplotlib.patches import Ellipse

class DECamField(object):
#
    def __init__(self, ra=ephem.degrees('0:0.0'), dec=ephem.degrees('0:0.0')):
        self.ra = ra
        self.dec = dec

    def ellipse(self):
        # An approximation to the DECam field of view, suitable e.g. for plotting
        semimajor_deg = ephem.degrees('1.08')
        semiminor_deg = ephem.degrees('0.98')
        center = (self.ra, self.dec)
        rotation = 0
        return Ellipse(center, 2*semimajor_deg, 2*semiminor_deg, rotation)
    
    def contains(self, ra1, dec1):
        # returns True if the point (ra1, dec1) lies inside the field
        radiff = ra1-self.ra
        if radiff > 2*np.pi: radiff -= 2*np.pi
        return radiff**2/ephem.degrees('1.08')**2 + (dec1-self.dec)**2/ephem.degrees('0.98')**2 <= 1.0
    
    

def main():
 # Try some SN fields:
    C1 = ephem.readdb("C1,f,03:37:05.83,-27:06:41.8,23.0,2000")
    C1.compute()
    C1field = DECamField(C1.ra, C1.dec)
    print C1field.ra, C1field.dec
    print C1field.contains(C1.ra, C1.dec+ephem.degrees('0.979'))
    print C1field.contains(ephem.degrees('04:00:00'),ephem.degrees('-30:00:00'))

if __name__=="__main__":
    main()