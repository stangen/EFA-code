#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 11:48:09 2018

@author: stangen
"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
from nicks_files.operational_cfsv2 import get_cfsv2_ensemble
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import efa_files.cfs_utilities_st as ut
import efa_files.madis_example.madis_utilities as mt
import time
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from efa_xray.state.ensemble import EnsembleState
from efa_xray.observation.observation import Observation
from efa_xray.assimilation.ensrf import EnSRF
import efa_xray.postprocess.postprocess as pp
from do_efa import Run_efa



ensemble_type = 'ecmwf'
variables = ['T2M', 'ALT']
ob_type = 'T2M'
#change this later
start_date = '20130401_0000'
end_date = '20130401_0000'

var_dict = {
        'ALT' : 'alts',
        'T2M' : 'temp'
        }

m_dict = {
           '01' : 'jan',
           '02' : 'feb',
           '03' : 'mar',
           '04' : 'apr',
           '05' : 'may',
           '06' : 'jun',
           '07' : 'jul',
           '08' : 'aug',
           '09' : 'sep',
           '10' : 'oct',
           '11' : 'nov',
           '12' : 'dec'           
           }

dates = mt.make_datelist(start_date, end_date, torf=True, timestep=86400)
for datestr in dates:
    year = datestr[0:4]
    month = datestr[4:6]
    day = datestr[6:8]
    hour = datestr[9:11]
    
    
    
    
    efa = Run_efa(datestr,ensemble_type,variables,m_dict,var_dict[ob_type])
    statecls, lats, lons, elevs = efa.load_data()
    
    obs = efa.load_obs(6)
    
    observations = []
    ob_all = []
    hx = []
    variance = []
    hx_mean = []
    se = []
    i = 0
    while i < 100:
    #for ob in obs:
    #ob = obs[0]
    #print(ob)
        ob=obs[i]
        ob_split = ob.split(',')
        #get the lat/lon of the station
        ob_lat = float(ob_split[1])
        ob_lon = float(ob_split[2])
        #get longitude positive-definite- ie -130 lon is 230 E
    #    if ob_lon < 0:
    #        ob_lon = ob_lon + 360
        ob_elev = float(ob_split[3])
        ob_time = mt.timestamp2utc(float(ob_split[4]))
        ob_value = float(ob_split[5])
        
        #call function to check elevation to see whether or not to assimilate ob
        TorF = efa.check_elevation(lats,lons,elevs,ob,ob_lat,ob_lon,ob_elev)
        
        obser = Observation(value=ob_value, time=ob_time,
                        lat=ob_lat,lon=ob_lon, obtype=ob_type, localize_radius=1000,
                        assimilate_this=TorF, error=1)
        
        if TorF == True:
            ob_all.append(ob)
            observations.append(ob_value)
            hx_oneob = obser.estimate(statecls)
            hx.append(hx_oneob)
            hx_variance = np.var(hx_oneob)
            #check for outlier data, 4 standard deviations from mean
            std_dev = np.sqrt(hx_variance)
            #if abs(ob_value-np.mean(hx_oneob)) < 10*std_dev:
            variance.append(hx_variance)
            hx_mean.append(np.mean(hx_oneob))
            se.append((np.mean(hx_oneob)-ob_value)**2)
            
        print(len(variance))
        i = i +1
    
    hx = np.array(hx)
    variance = np.array(variance)
    hx_mean = np.array(hx_mean)
    se = np.array(se)
    
    ##calculate the mean squared error
    mse = np.mean(se)
    print(mse)
    ##calculate the mean variance
    mean_var = np.mean(variance)
    print(mean_var)
        
    
    
    
        
    #thing = pp.obs_assimilation_statistics(statecls, observations)
    
        
