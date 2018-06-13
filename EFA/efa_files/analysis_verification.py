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
precip_analysis_path = '/home/disk/hot/stangen/Documents/GEFS/analysis/' + \
                        '2017090600_2017091600/ensembles/precip/2017090606_10days.nc'

# Filepath of the observed (analysis)
analysis_path = '/home/disk/hot/stangen/Documents/GEFS/analysis/' + \
                '2017090600_2017091600/ensembles/2017090600_10days.nc'
                
# only variable I am using is 500mb height            
vrbls=['Z500']

# Load the prior data            
with Dataset(precip_analysis_path,'r') as ncdata:
    times = ncdata.variables['time']
    ftimes = num2date(times[:],
                      times.units)
    lats = ncdata.variables['lat'][:]
    lons = ncdata.variables['lon'][:]
    for vrbl in vrbls:
        Z500_precip_anl = ncdata.variables[vrbl][:]
        Z500_units = ncdata.variables[vrbl].units

    T2M_units = ncdata.variables['T2M'].units
    T2M_precip_anl = ncdata.variables['T2M'][:]

# Load the analysis data
with Dataset(analysis_path, 'r') as ncdata:
    Z500_anl = ncdata.variables['Z500'][:] 
#    T500_post = ncdata.variables['T500'][:]
#    T2M_post = ncdata.variables['T2M'][:]
    T2M_anl = ncdata.variables['T2M'][:]
    
differenceT2M = T2M_precip_anl[0,:,:]-T2M_anl[1,:,:]
difference500mb = Z500_precip_anl[0,:,:]-Z500_anl[1,:,:]

# Plot the difference between the prior and posterior
map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
#map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))

# compute native map projection coordinates of lat/lon grid.
lon, lat = np.meshgrid(lons, lats)
x, y = map(lon, lat)

time_ind = 20
#cs = map.contourf(x,y,Z500_diff[0,:,:,0])#,levels=np.arange(-90, 91, 30), cmap=plt.cm.RdBu_r,linewidths=1.5, extend='both')
cs = map.contourf(x,y,difference500mb)
# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(Z500_units)

plt.title('{}'.format(ftimes[time_ind]))
plt.show()
    