#!/usr/bin/env python

from __future__ import division
import matplotlib.pyplot as plt
import sys
import os.path
from KBO import *
import numpy as np
import pandas as pd
import easyaccess as ea
from DESKBO import DEScands
from TNOfinder import TNOfinder
from TNOcandidate import TNOcandidate
import json
import argparse

def make_cat(infile): #, t_eff_min=0.3, link_bands=["r","i"], exptime_min=90):
	cat_in = Catalog(infile, date=str, ra=str, dec=str, nite=int, expnum=int, exptime=float, 
		band=str,ccdnum=int,mag=float,snobjid=int,ml_score=float, snfake_id=int)
	'''
	cat_in.refactor('ra',  lambda ra: hours(degrees(ra)))
	cat_in.refactor('dec', lambda dec: degrees(dec))
	cat_in.refactor('date', lambda date: toDateTime(date))  
	cat_in.rename(objid='snobjid')
	cat_in.rename(ccd='ccdnum')
	cat_in.rename(fakeid='snfake_id')
	cat_in.add_constant('obscode', 807)
	cat_in.add_constant('err', 0.15)
	cat_in.orderby('expnum')
	cat_in = Catalog(p for p in cat_in if p.band in link_bands and p.exptime>=exptime_min)
	expnums_in = set([p.expnum for p in cat_in])
	exposures_in = [e for e in exposures if e.expnum in expnums_in]
	good_exps = [e.expnum for e in exposures_in if e.t_eff>t_eff_min]
	cat_good = Catalog([p for p in cat_in if p.expnum in good_exps])
	return cat_good
	'''
	return cat_in

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("-c", "--conf", required=True, help="path to the JSON configuration file")
	args = vars(ap.parse_args())
	conf = json.load(open(args["conf"]))
	infile = conf["infile"]
	runid = conf["runid"]
	look_ahead_nights = conf["look_ahead_nights"]
	nominal_distance = conf["nominal_distance"]
	link_bands = conf["link_bands"]
	t_eff_min = conf["t_eff_min"]
	exptime_min = conf["exptime_min"]
	cat = make_cat(infile) # link_bands=link_bands, t_eff_min=t_eff_min, exptime_min=exptime_min)
	for p in cat._points[:10]:
		print p.date, p.ra, p.dec, p.expnum, p.exptime, p.ccd
	print 'Using input file ', infile
	print 'Linker run id: ', runid
	print 'Length of input catalog: ', len(cat) 
	finder = TNOfinder(cat, look_ahead_nights=look_ahead_nights, nominal_distance=nominal_distance, exclude_objids=None, runid=runid)
	triplets = finder.find_triplets(verbose=True)    # THIS CAN TAKE AWHILE
	cands = finder.tnocands(triplets)
	for cand in cands:
		orbit = cand[0]
		observations = cand[1]
		tnocand = TNOcandidate(orbit, observations, runid=runid, csvname=os.path.join('linker_output',runid))


if __name__=='__main__':
	main()