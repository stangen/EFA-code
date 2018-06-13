#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 16:19:04 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
import numpy as np
from nicks_files.operational_cfsv2 import get_cfsv2_ensemble
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import efa_files.cfs_utilities_st as ut
import time
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from EFA.efa_xray.state.ensemble import EnsembleState
from EFA.efa_xray.observation.observation import Observation
from EFA.efa_xray.assimilation.ensrf import EnSRF

class Run_efa():
    
    def __init__(self, datestr,ens_type,vrbls,ob_type):
        self.datestr = datestr
        self.y = datestr[0:4]
        self.m = datestr[4:6]
        self.d = datestr[6:8]
        self.h = datestr[9:11]
        self.ens_type = ens_type
        self.vrbls = vrbls
        self.m_dict = {
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
        self.ob_type = ob_type
        


    
    
    
    def load_data(self):
        # directory where the ensemble of all times is - ADDED HARDCODED FILENAME SPECIFIC FOR 1 FILE AT THE END
        infile = '/home/disk/hot/stangen/Documents/atms544/ensembles/'+self.ens_type+'/'+self.m_dict[self.m]+self.y+'/'+self.y+'-'+self.m+'-'+self.d+'_'+self.h+'_'+self.ens_type+'4var.nc' 
        
        print('loading netcdf file: '+self.y+self.m+self.d+'_'+self.h+'00')
        # loading/accessing the netcdf data            
        ncdata = Dataset(infile,'r')
        #print(ncdata.variables.keys())
        times = ncdata.variables['time']
        ftimes = num2date(times[:],
                          times.units)
        lats = ncdata.variables['lat'][:]
        lons = ncdata.variables['lon'][:]
        mems = ncdata.variables['ens'][:]
        #print(ncdata.variables)
        
        # storing the variable data in a dict (state?)
        allvars = {}
        for var in self.vrbls:
            allvars[var] = (['validtime','y','x','mem'],
                            ncdata.variables[var][:])
        lonarr, latarr = np.meshgrid(lons, lats)
        
        # Package into an EnsembleState object knowing the state and metadata
        statecls = EnsembleState.from_vardict(allvars,
                                      {'validtime' : ftimes,
                                       'lat' : (['y','x'], latarr),
                                       'lon' : (['y','x'], lonarr),
                                       'mem' : mems,
                                       })
        
        
        # directory where the orography file is
        orography = '/home/disk/hot/stangen/Documents/ensembles/orography/2013-04-01_00_'+self.ens_type+'.nc'
        orog_data = Dataset(orography,"r")
        #print(of.variables)
        elevs = orog_data.variables['orog'][0,:] # Elevation of ecmwf or eccc         
        
        return statecls, lats, lons, elevs

        
    def load_obs(self, forecast_hour):
        
        #get observations from n hours after the model was initialized
        dt0 = datetime.strptime(self.datestr,'%Y%m%d_%H%M')
        dt = dt0.replace(minute = 0, second=0, microsecond=0)
        dt = dt + timedelta(0,3600*forecast_hour)
        
        #convert back to strings
        dtstr = dt.strftime('%Y%m%d_%H00')
        dty = dtstr[0:4]
        dtm = dtstr[4:6]
        dtd = dtstr[6:8]
        dth = dtstr[9:11]
        
        # directory where the observations are
        obs_file = '/home/disk/hot/stangen/Documents/EFA/surface_obs/MADIS/'+dty+dtm+'/combined_'+self.ob_type+'/'+self.ob_type+'_'+dty+dtm+dtd+'_'+dth+'00.txt'
        print('loading '+self.ob_type+ ' from '+dty+dtm+dtd+'_'+dth+'00')
        #print(obs_file)
        f1 = open(obs_file, 'r')
        obs = f1.readlines()
        
        return obs
        
    
    
    def check_elevation(self,lats,lons,elevs,ob,ob_lat,ob_lon,ob_elev):
        #make a lat/lon array for comparison with station lat/lon
        lonarr, latarr = np.radians(np.meshgrid(lons, lats))
        
        R = 6371. # Radius of earth in km
        #convert degrees to radians
        ob_lat = np.radians(ob_lat)
        ob_lon = np.radians(ob_lon)
        #radian distance between points
        dlat = ob_lat-latarr
        dlon = ob_lon-lonarr
        #Haversine formula to calculate the great circle distance between ob and each lat/lon point (in km)
        a = np.sin(dlat/2)**2 + np.cos(ob_lat) * np.cos(latarr) * np.sin(dlon/2)**2
        dist = 2*R*np.arcsin(np.sqrt(a))
        #flatten the distance array to 1D
        dist_flat = dist.flatten()
        
        #find nearest 4 points, in order from closest to farthest
        #argpartition just puts the 4 smallest indices in front, array (1,k) puts the 
        #first four in order, the rest are not necessarily sorted
        k = 4
        idx = np.argpartition(dist_flat, (1,k))
        
        #check if the 4 closest points have elevations which differ by more than 300 meters
        elev_flat = elevs.flatten()
        elev_four_closest = elev_flat[idx[:k]]
        elev_diff = abs(ob_elev-elev_four_closest)
        #print(elev_diff)
        if elev_diff.max() > 300:
            TorF = False
        else:
            TorF = True
        
        return TorF
    
