#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 13:30:52 2017

@author: stangen
"""

"""
Created on Fri Sep  1 15:21:40 2017

@author: stangen
"""

"""
Goal: Determine if EFA creates more accurate forecasts 6-10 days in the
future. 
To do this: 
    1. Find difference between GEFS mean (prior) and observed for various
    parameters, starting with 500 mb height.
    2. Find difference between EFA-adjusted GEFS (posterior) and observed for 
    various parameters, starting with 500 mb height.
    3. Determine if EFA reduces error 6-10 days into the future.
    
"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
from nicks_files.operational_cfsv2 import get_cfsv2_ensemble
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import nicks_files.cfs_utilities as ut
import time
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from efa_xray.state.ensemble import EnsembleState
from efa_xray.observation.observation import Observation
from efa_xray.assimilation.ensrf import EnSRF

# Filepath of the precip analysis
orog1_path = '/home/disk/hot/stangen/Documents/prior_ensembles/ecmwf/201304/2013-04-10_00_ecmwf_T2M_ALT.nc'

# Filepath of the observed (analysis)
orog2_path = '/home/disk/hot/stangen/Documents/posterior_ensembles/ob_update_self/loc_1000/ecmwf/201304/2013-04-10_00_ecmwf_T2M_ALT.nc'
  
orog3_path = '/home/disk/hot/stangen/Documents/prior_ensembles/ecmwf/201304/2013-04-10_12_ecmwf_T2M_ALT.nc'             
# only variable I am using is 500mb height            
vrbls=['T2M']

# Load the prior data            
with Dataset(orog1_path,'r') as ncdata:
    #print(ncdata.variables)
    times = ncdata.variables['time']
    ftimes = num2date(times[:],times.units)
    lats = ncdata.variables['lat'][:]
    lons = ncdata.variables['lon'][:]
    for vrbl in vrbls:
        orog1 = ncdata.variables[vrbl][:]
        orog1_units = ncdata.variables[vrbl].units


# Load the analysis data
with Dataset(orog2_path, 'r') as ncdata:
    orog2 = ncdata.variables[vrbls[0]][:] 
    times2 = ncdata.variables['time']

#    T500_post = ncdata.variables['T500'][:]
#    T2M_post = ncdata.variables['T2M'][:]
with Dataset(orog3_path, 'r') as ncdata:
    orog3 = ncdata.variables[vrbls[0]][:] 
    times3 = ncdata.variables['time']

prior_mean = np.mean(orog1,axis=3)    
post_mean = np.mean(orog2,axis=3)
analysis_mean = np.mean(orog3,axis=3)
post_prior_diff= post_mean[1,:,:]-prior_mean[1,:,:]
anal_prior_diff = prior_mean[1,:,:]-analysis_mean[0,:,:]
anal_post_diff = post_mean[1,:,:]-analysis_mean[0,:,:]


# Plot the difference between the prior and posterior
#map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
map = Basemap(projection='merc',llcrnrlat=25,urcrnrlat=35,\
            llcrnrlon=-105,urcrnrlon=-95,lat_ts=30,resolution='c')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
#map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
#map.fillcontinents(color='coral',lake_color='aqua')
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0,360,1))
map.drawparallels(np.arange(-90,90,1))

# compute native map projection coordinates of lat/lon grid.
lon, lat = np.meshgrid(lons, lats)
x, y = map(lon, lat)

time_ind = 2
#cs = map.contourf(x,y,Z500_diff[0,:,:,0])#,levels=np.arange(-90, 91, 30), cmap=plt.cm.RdBu_r,linewidths=1.5, extend='both')
cs = map.contourf(x,y,post_prior_diff[:,:])
# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(orog1_units)

plt.title('{}'.format(ftimes[time_ind]))
plt.show()
    