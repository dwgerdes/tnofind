#!/usr/bin/env python

from __future__ import division
import pandas as pd
import numpy as np
from KBO import pretty_nite, toDateTime
from ephem import Ecliptic, Equatorial, hours, degrees, date, separation
import linkutils
from scipy.spatial import KDTree
import time
from multiprocessing import Pool, Manager, Queue, cpu_count
import json
import argparse
import cPickle as pickle
import os

# Replace these lines with input from a config file

#infile = 'wsdiff_catalogs/scrambled_nofakes/wsdiff_S82_nofakes_riz_scrambled.csv'
#infile = 'wsdiff_catalogs/real_data/wsdiff_S82_riz_nofakes.csv'
infile = 'wsdiff_catalogs/fakes_only/wsdiff_S82_riz_fakes.csv'

linkmap_out = 'wsdiff_catalogs/fakes_only/linkmap_fast_test.pickle'

arcmin = np.pi/(180*60)  # one arcminute in radians
df_exps = pd.read_csv('exposures.csv')   # master list of all exposures

def make_cat(infile):
	df =  pd.read_csv(infile)
	# rename some columns:
	df.rename(columns={'date_obs':'date', 'snobjid':'objid', 'snfake_id':'fakeid'}, inplace=True)
	df['ra'] = df.apply(lambda row: hours(degrees(str(row['ra']))), axis=1)
	df['dec'] = df.apply(lambda row: degrees(str(row['dec'])), axis=1)
	df['date'] = df.apply(lambda row: toDateTime(row['date']), axis=1)
	return df

def wrap_ra(r):
    return r if 0<r<np.pi else r-2*np.pi

def build_kdtree(df):
    '''
    Builds a KDtree from the input dataframe
    '''
    x = np.array([[wrap_ra(r) for r in df['ra'].values], 
                  df['dec'].values]).transpose()
    return KDTree(x)

def nites_between(nite1, nite2):
    '''
    Returns the time in days between two nites, where nite2>nite1.
    '''
    return date(pretty_nite(nite2))-date(pretty_nite(nite1))


def apply_offset(df, offset_ra, offset_dec):
	'''
	Applies an ra, dec offset to the coordinates in a dataframe. The dataframe is assumed
	to have coords that are ephem.Angle objects, and offsets are in radians.
	'''
	df_offset = df[['ra','dec','objid']]
	df_offset['ra'] += offset_ra
	df_offset['dec'] += offset_dec
	return df_offset


def cosine_cut(displacement_arcsec):
    '''
    The requirement on directional alignment between two points and that expected from earth reflex motion.
    It's more generous for lower values of the displacement (in arcseconds).
    '''
    if displacement_arcsec<400:
        return 0.7
    elif 400<displacement_arcsec<=800:
        return 0.7 + (0.95-0.7)/(800-400)*(displacement_arcsec-400)
    else:
            return 0.95

def tno_like(point1, point2, vmax=150, debug=False):
	#
	#  vmax is max velocity in arsec/day
	#
    lon1, lat1 = Ecliptic(Equatorial(point1['ra'], point1['dec'])).get()
    lon2, lat2 = Ecliptic(Equatorial(point2['ra'], point2['dec'])).get()
    displacement = separation((lon2, lat2), (lon1, lat1))
    displacement_asec = displacement*180/np.pi*3600
    if point2['date'] != point1['date']:
        velocity = displacement_asec/(point2['date'] - point1['date'])
    else: 
        velocity=9999
    this_dlon, this_dlat, this_vlon, this_vlat = linkutils.parallax(point1['ra'], point1['dec'], point1['date']) 
    next_dlon, next_dlat, next_vlon, next_vlat = linkutils.parallax(point2['ra'], point2['dec'], point2['date'])   
    dlon, dlat = next_dlon-this_dlon, next_dlat-this_dlat
    dot = np.cos(lat1)**2*(lon2 - lon1)*dlon + (lat2 - lat1)*dlat
    norm = np.sqrt(np.cos(lat1)**2*dlon**2 + dlat**2)
    cosine = dot/(norm*displacement)
    TNOlike = True if velocity<vmax and cosine>cosine_cut(displacement_asec) else False
    return TNOlike


def get_links(args):
	this_nite, q = args
	'''
	The key routine for building the linkmap. Examines all points in this_nite and identifies points in 
	exposures within look_ahead_nites of this_nite with separation consistent with parallax motion of a distant
	solar system object.
	'''
	this_nites_expnums = np.unique(df_nites[this_nite]['expnum'].values)
	link_nites = [nite for nite in nites if 0<nites_between(this_nite, nite)<=look_ahead_nites]
	# Compute parallax offsets from each exposure in current nite to the target nites
	print 'Processing nite: ', this_nite
	print '    Nite contains', len(df_nites[this_nite]),' points in ', len(this_nites_expnums), ' exposures'
	print '    ', len(link_nites),' nites available for linking'
	linkmap_nite = {}
	for exp in this_nites_expnums: 
		exposure = df_exps.loc[df_exps['expnum']==exp]
 		exp_ra, exp_dec = hours(exposure['ra'].values[0]), degrees(exposure['dec'].values[0])
		exp_date = date(exposure['date'].values[0])
		lon, lat = Ecliptic(Equatorial(exp_ra, exp_dec)).get()
		this_dlon, this_dlat, this_vlon, this_vlat = linkutils.parallax(exp_ra, exp_dec, exp_date)
		this_exps_points = df_nites[this_nite].loc[df_nites[this_nite]['expnum']==exp] # points in this exposure
		this_exps_points.index = range(len(this_exps_points))   # reindex starting from zero
		this_exps_objids = this_exps_points['objid'].values     # objids in this exposure
		linkmap_exp = dict(zip(this_exps_objids, [[] for i in range(len(this_exps_objids))]))
		for next_nite in link_nites:
			deltaT = nites_between(this_nite, next_nite)
			next_dlon, next_dlat, next_vlon, next_vlat = linkutils.parallax(exp_ra, exp_dec, exp_date+deltaT)
			dlon, dlat = next_dlon-this_dlon, next_dlat-this_dlat
			search_center = Equatorial(Ecliptic(lon+dlon/nominal_distance, lat+dlat/nominal_distance))
			deltaR = np.sqrt(dlon**2*np.cos(lat)**2 + dlat**2)/nominal_distance
			offset_ra, offset_dec = hours(search_center.ra - exp_ra), degrees(search_center.dec - exp_dec)
			shifted_points = apply_offset(this_exps_points, offset_ra, offset_dec)  # apply parallax offset
			this_exps_tree = build_kdtree(shifted_points)  # make KDtree from the shifted points
			near_list = this_exps_tree.query_ball_tree(trees[next_nite], deltaR)
			for ipoint in range(len(near_list)):
			    if near_list[ipoint] != []:
				    point1 = this_exps_points.iloc[ipoint]
				    for jpoint in near_list[ipoint]:
						point2 = df_nites[next_nite].iloc[jpoint]
						if tno_like(point1, point2):
							objid1 = point1['objid']
							objid2 = point2['objid']
							linkmap_exp[objid1].append(objid2)
		linkmap_nite.update(linkmap_exp)
	q.put(this_nite)
	return linkmap_nite

if __name__ == '__main__':

	ap = argparse.ArgumentParser()
	ap.add_argument("-c", "--conf", required=True, help="path to the JSON configuration file")
	args = vars(ap.parse_args())
	conf = json.load(open(args["conf"]))
	infile = conf["infile"]
	output_dir = conf["output_dir"]
	runid = conf["runid"]
	look_ahead_nites = conf["look_ahead_nights"]
	nominal_distance = conf["nominal_distance"]
	Ncpu = conf["Ncpu"]
#	t_eff_min = conf["t_eff_min"]
#	exptime_min = conf["exptime_min"]

	df = make_cat(infile)
	nites = np.unique(df['nite'].values)
	df_nites = dict(zip(nites, [df[df['nite']==n] for n in nites]))  # dictionary of nightly dataframes. Key is nite.
	for n in nites:
   		df_nites[n].index = range(len(df_nites[n]))  # reindex each nightly df
	trees  = dict(zip(nites, [build_kdtree(df_nites[n]) for n in nites]))
   	print 'Done building KDtrees'

   	print 'Begin linkmap generation, Ncpu = ', Ncpu
	linkmap = {}
	pool = Pool(Ncpu)
	manager = Manager()
	queue = manager.Queue()
	args = [(n, queue) for n in nites]
	result = pool.map_async(get_links, args)

	# monitor loop
	size_old = 0
	while True:
		if result.ready():
			break
		else:
			size = queue.qsize()
			if size > size_old:
				print 'linkmap progress: ', size,' of ', len(nites), ' nites completed'
				size_old = size
			time.sleep(0.1)

	for d in result.get():
		linkmap.update(d)


	print 'Done with linkmap generation!'

	linkmap_out = os.path.join(output_dir,str(runid)+'_linkmap.pickle')
	with open(linkmap_out,'wb') as f:
		pickle.dump(linkmap, f)
	print 'Wrote linkmap to file ', linkmap_out
