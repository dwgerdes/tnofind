#!/usr/bin/env python

import ephem


class MPCRecord(object):
    
    #
    # Defines a minor planet optical observation in the standard 80-column MPC format.
    # See: http://www.minorplanetcenter.net/iau/info/OpticalObs.html
    # For format of packed provisional designations see http://www.minorplanetcenter.net/iau/info/PackedDes.html
    # Observation codes for note1: see http://www.minorplanetcenter.net/iau/info/ObsNote.html
    # note2 defaults to 'C' for a CCD observation.
    # observatoryCode defaults to W84 for Tololo/DECam, others at
    # http://www.minorplanetcenter.net/iau/lists/ObsCodes.html
    
    
    def __init__(self, obsnum='     ', MPprovisional='       ', discovery=' ', note1=' ', 
                 note2='C', obsdate=ephem.date('2000/01/01'), ra_obs_J2000=ephem.hours(0), dec_obs_J2000=ephem.degrees(0), 
                 mag=99, band='r', observatoryCode='W84', newobject=True):
            self.spaces  = []
            a = ''
            one_space = ' '
            for i in range(20):
                self.spaces.append(a)
                a += one_space
    
            self.newobject = newobject 
            self.obsnum = self.setobsnum(obsnum)
            self.MPprovisional = self.setMPprovisional(MPprovisional)
            if discovery==' ' or discovery=='*':
                self.discovery = discovery
            else:
                print 'Error, unrecognized discovery flag, should be * or one blank space'
            if len(note1)==1:
                self.note1 = note1
            else:
                print 'Error, note1 must have length 1'
            if len(note2)==1:
                self.note2 = note2
            else:
                print 'Error, note2 must have length 1'
            self.obsdate = self.setObsdate(obsdate)                  # YYYY MM DD.dddddd, with precision normally given to 0.00001d (about 0.8s). For DES data should use 0.000001d
            self.ra_obs_J2000 = self.setRa(ra_obs_J2000)
            self.dec_obs_J2000 = self.setDec(dec_obs_J2000)
            self.mag_band = self.setMagBand(mag, band)
            self.observatoryCode = observatoryCode
            self.record = self.obsnum+self.MPprovisional+self.discovery+self.note1+self.note2+self.obsdate+self.ra_obs_J2000 \
                          +self.dec_obs_J2000+self.spaces[9]+self.mag_band+self.spaces[6]+self.observatoryCode
    
    def setobsnum(self, num):
        if self.newobject is True:
            thisobs=self.spaces[5]
        else:
            if len(num)==5:
                thisobs = num
            elif len(num)>5:
# Assume we are working with the unpacked designation, and pack it.
                thisobs=self.packedDesignation(num)
            else:
                while len(num)<5:
                    num = '0'+num
                thisobs=num
        return thisobs
                    
            
    def setMPprovisional(self, prov):
        thisprov=self.spaces[7]
        if self.newobject is True:     # up first six characters can be whatever we want
            if len(prov)<7:
                thisprov = prov+self.spaces[7-len(prov)]
        else:
            if ' ' in prov:                # assume we've got the unpacked designation, and pack it
                thisprov = self.packedDesignation(prov)
            elif ' ' not in prov and len(prov)==7: # assume it's already packed, leave it alone
                thisprov = prov
            else:                          # don't recognize this
                print 'Error, unrecognized MP designation'
                thisprov = self.spaces[7]
        assert len(thisprov)==7
        return thisprov
            
    def packedDesignation(self, prov):
        """
        Returns the packed form of the provisional designation, see http://www.minorplanetcenter.net/iau/info/PackedDes.html
        """
        year, code = prov.split(' ')
        assert (len(year)==4 and len(code)<6 and len(code)>1)
        century = year[:2]
        if century=='18':
            packed = 'I'
        elif century=='19':
            packed='J'
        elif century=='20':
            packed='K'
        else:
            print 'Error, incorrect year in MP provisional designation'
        packed += year[-2:]
#       Now the code
        half_month_1 = code[0]
        half_month_2 = code[1]
        if len(code)==2:
            cycle = '00'
        elif len(code)==3:
            cycle = '0'+code[-1]
        elif len(code)==4:
            cycle=code[-2:]
        else:
            cycle_count = int(code[2:])
            if cycle_count <= 99:
                cycle = '0'+str(cycle_count)
            else:
                head = int(str(cycle_count)[:2])
                assert head>=10
                alpha='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
                if head<62:
                    cycle=alpha[head-10]
                else:
                    print 'Error, invalid cycle count'
                cycle = cycle+code[-1]     # last digit of cycle count goes at end
        packed += half_month_1 + cycle + half_month_2
        return packed
    
    def check_packed(self):
        full = ['1995 XA', '1995 XL1', '1995 FB13', '1998 SQ108', '1998 SV127', '1998 SS162', '2099 AZ193', '2008 AA360', '2007 TA418', '2012 VP113']
        packed =['J95X00A','J95X01L', 'J95F13B','J98SA8Q','J98SC7V','J98SG2S', 'K99AJ3Z', 'K08Aa0A', 'K07Tf8A', 'K12VB3P' ]
        pairs = zip(full, packed)
        for p in pairs:
            pk = packedDesignation(p[0])
            if pk != p[1]:
                print 'Error, expected ', p[1], ', got ', pk
                
    def setObsdate(self, d):    # assume d is an ephem.date object
        ddd=d.triple()
        day = "%.6f" % ddd[2]
        if len(str(ddd[1]))==1:
            mm = '0'+str(ddd[1])
        else:
            mm = str(ddd[1])
        if len(day)==8: day='0'+day
        yyyymmdd = str(ddd[0])+' '+mm+' '+day
        assert len(yyyymmdd) == 17
        return yyyymmdd
    
    def setRa(self, ra):      # assume ra is an ephem.Angle
        hms = str(ra).split(':')
        hh = hms[0]
        if len(hh)==1: hh='0'+hh
        mm = hms[1]
        if len(mm)==1: mm='0'+mm
        ss = hms[2]
        if len(ss) == 4: ss='0'+ss
        ss=ss+' '
        ra_mpc = hh+' '+mm+' '+ss
        assert len(ra_mpc)==12
        return ra_mpc
        
    def setDec(self, dec):
        dms = str(dec).split(':')
        d = dms[0]
        if d[0] != '-':
            if len(d)==1:
                dd='+0'+d
            else:
                dd='+'+d
        else:
            if len(d)==2:
                dd='-0'+d[1:]
            else:
                dd=d
        assert len(dd)==3
        mm = dms[1]
        if len(mm)==1: mm='0'+mm
        ss = dms[2]
        if len(ss) == 3: ss='0'+ss
        ss=ss+' '
        dec_mpc = dd+' '+mm+' '+ss
        assert len(dec_mpc)==12
        return dec_mpc
    
    def setMagBand(self, mag, band):
        if (mag<10 or mag>30):
            mpcmag=self.spaces[5]
        else:
            mpcmag = "%4.1f" % mag + ' '
        assert len(mpcmag)==5
        assert len(band)==1 and band in ['g','r','i','z','Y']
        return mpcmag+band
        
        
        
    
    
def main():
    pass

if  __name__ == '__main__':
    main()
            
