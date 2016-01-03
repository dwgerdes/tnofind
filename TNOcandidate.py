import numpy as np
import pandas as pd
from KBO import *
import itertools
import uuid
import os.path
import easyaccess as ea
import sqlalchemy
from orbfit_tools.orbfit_tools import mean_anomaly

class TNOcandidate(object):
    def __init__(self, orbit, observations, runid, csvname=None, writedb=False):
    	self.orbit = orbit
    	self.observations = observations
    	self.runid = runid                    # linker run id 
    	self.id = uuid.uuid4()                # unique id for this candidate
    	self.engine = self.desdb_engine() if writedb else None
    	self.df_obs = self.obs2df()
    	self.df_link = self.link2df()
    	self.df_orbit = self.orbit2df()
    	if csvname is not None:
    		csv_obs = str(csvname+'_obs.csv')
    		csv_link = str(csvname+'_link.csv')
    		csv_orbit = str(csvname+'_orbit.csv')
    		if os.path.isfile(csv_obs):  # if file already exists, append to it and don't write the header.
    			self.df_obs.to_csv(csv_obs, header=False, index=False, mode='a')  
    		else:
    			self.df_obs.to_csv(csv_obs, index=False)
    		if os.path.isfile(csv_link):  # if file already exists, append to it and don't write the header.
    			self.df_link.to_csv(csv_link, header=False, index=False, mode='a')  
    		else:
    			self.df_link.to_csv(csv_link, index=False)
    		if os.path.isfile(csv_orbit):  # if file already exists, append to it and don't write the header.
    			self.df_orbit.to_csv(csv_orbit, header=False, index=False, mode='a')  
    		else:
    			self.df_orbit.to_csv(csv_orbit, index=False)


    def desdb_engine(self, section='desoper', verbose=True):
    	'''
    	Establishes the connetion to the DESDB. Eventually can get this info from the user's .desservices file
    	'''
    	engine = sqlalchemy.create_engine("oracle+cx_oracle://gerdes:thUmp2ski@(DESCRIPTION = \
    		(LOAD_BALANCE=on) (FAILOVER=ON) (ADDRESS = (PROTOCOL = TCP)(HOST = leovip148.ncsa.uiuc.edu)(PORT = 1521)) \
    		(CONNECT_DATA = (SERVER = DEDICATED) (SERVICE_NAME = "+section+")))")
    	return engine

    def obs2df(self):
    	'''
    	Insert the Catalog of observations into a pandas dataframe
    	'''
    	# These next two lines are a kludge. They should be part of the input catalog.
    	self.observations.add_constant('obscode', 807)
    	self.observations.add_constant('fakeid', 0)
    	obs_cols = ['date_obs', 'ra', 'dec', 'expnum', 'exptime', 'band', 'ccd', 'mag', 'nite', 'objid', 'fakeid', 'date_added', 
    	            'obscode']
    	index = np.arange(len(self.observations))  
    	df = pd.DataFrame(columns=obs_cols, index=index)
    	ind=0
    	for obs in self.observations:
    		obsdata = np.array([obs.date, obs.ra, obs.dec, obs.expnum, obs.exptime, obs.band, obs.ccd, obs.mag, obs.nite, obs.objid, obs.fakeid, 
    			ephem.now(), obs.obscode])
    		df.ix[ind] = obsdata
    		ind+=1
    	return df

    def obs2db(self, obs_df, section='desoper', table='tnobs'):
    	'''
    	Write observation catalog to the database
    	'''
    	engine = sqlalchemy.create_engine("oracle+cx_oracle://gerdes:thUmp2ski@(DESCRIPTION = \
    		(LOAD_BALANCE=on) (FAILOVER=ON) (ADDRESS = (PROTOCOL = TCP)(HOST = leovip148.ncsa.uiuc.edu)(PORT = 1521)) \
    		(CONNECT_DATA = (SERVER = DEDICATED) (SERVICE_NAME = "+section+")))")
    	obs_df.to_sql(table, engine, if_exists='replace',index=False)

    def link2df(self):
    	'''
    	Insert the linked observations into a pandas dataframe. Use id, objid to join to the orbits and observations dataframes or db tables.
    	'''
    	link_cols = ['id','linkrun','objid']
    	index = np.arange(len(self.observations))
    	df = pd.DataFrame(columns=link_cols, index=index)
    	ind=0
    	for obs in self.observations:
    		linkdata = np.array([self.id, self.runid, obs.objid])
    		df.ix[ind] = linkdata
    		ind+=1
    	return df

    def orbit2df(self):
    	'''
    	Insert the orbit details into a pandas dataframe. 
    	'''
    	orbit_cols = ['id','linkrun', 'chisq', 'ndof', 'a', 'e', 'i', 
    	'aop', 'node', 'peri_jd', 'peri_date', 'epoch_jd',
    	'mean_anomaly', 'period', 'period_err', 
    	'a_err', 'e_err', 'i_err', 'aop_err', 'node_err', 'peri_err', 'lat0', 'lon0',
    	'xBary','yBary','zBary',
    	'abg_a', 'abg_b', 'abg_g', 
    	'abg_adot', 'abg_bdot', 'abg_gdot', 
    	'abg_a_err', 'abg_b_err', 'abg_g_err',
    	'abg_adot_err', 'abg_bdot_err', 'abg_gdot_err']
    	index = np.arange(1)
    	df = pd.DataFrame(columns=orbit_cols, index=index)

    	elements, errs = self.orbit.get_elements()
    	elements_abg, errs_abg = self.orbit.get_elements_abg()
    	epoch = self.orbit.jd0

    	peri_jd = -999 if np.isnan(elements['top']) else elements['top']
    	peri_date = -999 if np.isnan(elements['top']) else ephem.date(elements['top']-2415020)
    	peri_jd_err = -999 if (np.isnan(elements['top']) or np.isnan(errs['top'])) else errs['top']

    	orbit_data = np.array([self.id, self.runid, round(self.orbit.chisq,2), self.orbit.ndof, round(elements['a'],3), round(elements['e'],4), round(elements['i'],4),
    		round(elements['aop'],2), round(elements['lan'],3), round(peri_jd,2), peri_date, round(epoch,2),
    		mean_anomaly(elements, epoch), elements['a']**1.5, 3*np.sqrt(elements['a'])/2*errs['a'],
    		round(errs['a'],2), round(errs['e'],4), round(errs['i'],4), round(errs['aop'],2), round(errs['lan'],3), round(peri_jd_err,2), round(self.orbit.lat0,4), round(self.orbit.lon0,4),
    		round(self.orbit.xBary,4), round(self.orbit.yBary,4), round(self.orbit.zBary,4), 
    		elements_abg['a'], elements_abg['b'], elements_abg['g'],
    		elements_abg['adot'], elements_abg['bdot'], elements_abg['gdot'],
    		errs_abg['a'], errs_abg['b'], errs_abg['g'],
    		errs_abg['adot'], errs_abg['bdot'], errs_abg['gdot']
    		])
    	df.ix[0] = orbit_data
    	return df




