#!/usr/bin/env python

from __future__ import division
import matplotlib.pyplot as plt
import sys
import os.path
from KBO import *
import numpy as np
import pandas as pd
import collections
import easyaccess as ea
from DESKBO import DEScands
from TNOfinder import TNOfinder
from TNOcandidate import TNOcandidate

def make_cat(infile):
	cat_in = Catalog(infile, date=str, ra=str, dec=str, nite=int, expnum=int, exptime=float, \
	                    band=str,ccd=int,mag=float,snobjid=int,pixelx=float,pixely=float,ml_score=float)
	cat_in.refactor('ra',  lambda ra: hours(degrees(ra)))
	cat_in.refactor('dec', lambda dec: degrees(dec))
	cat_in.refactor('date', lambda date: toDateTime(date))  
	cat_in.rename(objid='snobjid')
	cat_in.add_constant('obscode', 807)
	cat_in.add_constant('err', 0.15)
	cat_in.add_constant('fakeid',0)
	cat_in.orderby('expnum')
	cat_in = Catalog(p for p in cat_in if p.band in ['i','r'] and p.exptime>=90)
	cat_in.orderby('expnum')
	expnums_in = set([p.expnum for p in cat_in])
	exposures_in = [e for e in exposures if e.expnum in expnums_in]
	good_exps = [e.expnum for e in exposures_in if e.t_eff>0.3]
	cat_good = Catalog([p for p in cat_in if p.expnum in good_exps])
	return cat_good

def main():
	infile = sys.argv[-2]
	runid = sys.argv[-1]
	cat = make_cat(infile)
	print 'Using input file ', infile
	print 'Linker run id: ', runid
	print 'Length of input catalog: ', len(cat) 
	finder = TNOfinder(cat, look_ahead_nights=30, nominal_distance=40, exclude_objids=None, runid=runid)
	triplets = finder.find_triplets(verbose=True)    # THIS CAN TAKE AWHILE
	cands = finder.tnocands(triplets)
	for cand in cands:
		orbit = cand[0]
		observations = cand[1]
		tnocand = TNOcandidate(orbit, observations, runid=runid, csvname=os.path.join('linker_output',runid))


if __name__=='__main__':
	main()