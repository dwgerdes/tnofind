
#!/usr/bin/env python

from __future__ import division
import pandas as pd
import numpy as np
from ephem import Ecliptic, Equatorial, hours, degrees, date, separation
from multiprocessing import Pool, Manager, Queue, cpu_count
import time
import json
import argparse
import cPickle as pickle
import os
from linkmap import make_cat

def getDatePos(objid):
    '''
    Return date, ra, dec for a given objid
    '''
    point = df.loc[df['objid']==objid]
    ra = point['ra'].values[0]
    dec = point['dec'].values[0]
    date = point['date'].values[0]
    return date, ra, dec

def vEcliptic(id1, id2):
    '''
    Compute the components of velocity between two points in ecliptic coords. Units: arcsec/day.
    '''
    date1, ra1, dec1 = getDatePos(id1)
    date2, ra2, dec2 = getDatePos(id2)
    ecl1 = Ecliptic(Equatorial(ra1, dec1))
    ecl2 = Ecliptic(Equatorial(ra2, dec2))
    lat1, lon1 = ecl1.lat, ecl1.lon
    lat2, lon2 = ecl2.lat, ecl2.lon
    dt = date2-date1
    dlon = np.mod(np.abs(lon2-lon1), 2*np.pi)
    if dlon > np.pi: dlon = 2*np.pi - dlon
    vlon = dlon/dt*np.cos(lat1)*180/np.pi*3600
    vlat = (lat2-lat1)/dt*180/np.pi*3600
    return vlon, vlat

def get_vmap(args):
    objid, q = args
    q.put(objid)
    return {objid:dict(zip(linkmap[objid], [vEcliptic(objid, linkid) for linkid in linkmap[objid]]))}


if __name__ == '__main__':

    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--conf", required=True, help="path to the JSON configuration file")
    args = vars(ap.parse_args())
    conf = json.load(open(args["conf"]))
    infile = conf["infile"]
    output_dir = conf["output_dir"]
    runid = conf["runid"]
    Ncpu = conf["Ncpu"]

    df = make_cat(infile)
    linkmap_in = os.path.join(output_dir,str(runid)+'_linkmap.pickle')
    vmap_out= os.path.join(output_dir,str(runid)+'_vmap.pickle')

    vmap = {}
    linkmap = pickle.load(open(linkmap_in, 'rb'))
    print 'Read linkmap file ', linkmap_in
    objids = linkmap.keys()
    print 'Begin vmap generation, Ncpu = ', Ncpu
    pool = Pool(Ncpu)
    manager = Manager()
    queue = manager.Queue()
    args = [(key, queue) for key in objids]
    result = pool.map_async(get_vmap, args)

        # monitor loop
    size_old = 0
    while True:
        if result.ready():
            break
        else:
            size = queue.qsize()
            if size > size_old:
                print 'vmap progress: ', size,' of ', len(objids), ' objids completed (', round(size/len(objids)*100,1),'%)'
                size_old = size
            time.sleep(5)

    for d in result.get():
        vmap.update(d)


    print 'Done with vmap generation!'

    with open(vmap_out,'wb') as f:
        pickle.dump(vmap, f)
    print 'Wrote vmap to file ', vmap_out

