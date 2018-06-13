#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

# Load in the EFA-produced forecast and the observed.

# Filepath of the prior
efa_prior_path = '/home/disk/hot/stangen/Documents/GEFS/ensembles' + \
                 '/2017090600/2017090600_21mem_10days.nc'
            
# Filepath of the posterior
efa_post_path = '/home/disk/hot/stangen/Documents/GEFS/posterior' + \
                '/2017090600/2017090600_21mem_10days_efa.nc'
                
# Filepath of the observed (analysis)
analysis_path = '/home/disk/hot/stangen/Documents/GEFS/analysis/' + \
                '2017090600_2017091600/ensembles/2017090600_10days.nc'
                
# only variable I am using is 500mb height            
vrbls=['Z500']

# Load the prior data            
with Dataset(efa_prior_path,'r') as ncdata:
    times = ncdata.variables['time']
    ftimes = num2date(times[:],
                      times.units)
    lats = ncdata.variables['lat'][:]
    lons = ncdata.variables['lon'][:]
    mems = ncdata.variables['ens'][:]
    for vrbl in vrbls:
        Z500_prior = ncdata.variables[vrbl][:]
        Z500_units = ncdata.variables[vrbl].units
#    T500_prior = ncdata.variables['T500'][:]
#    T500_units = ncdata.variables['T500'].units
#    T2M_prior = ncdata.variables['T2M'][:]
    MSLP_prior = ncdata.variables['MSLP'][:]
    MSLP_units = ncdata.variables['MSLP'].units
#    PWAT_prior = ncdata.variables['PWAT'][:]
#    PWAT_units = ncdata.variables['PWAT'].units

# Load the posterior data    
with Dataset(efa_post_path, 'r') as ncdata:
    Z500_post = ncdata.variables['Z500'][:] 
#    T500_post = ncdata.variables['T500'][:]
#    T2M_post = ncdata.variables['T2M'][:]
    MSLP_post = ncdata.variables['MSLP'][:]
#    PWAT_post = ncdata.variables['PWAT'][:]

# Load the analysis data
with Dataset(analysis_path, 'r') as ncdata:
    Z500_anl = ncdata.variables['Z500'][:] 
#    T500_post = ncdata.variables['T500'][:]
#    T2M_post = ncdata.variables['T2M'][:]
    MSLP_anl = ncdata.variables['MSLP'][:]
    
   
## Subtract the difference between the prior/post and analysis for each point.
#Z500_prior_diff = Z500_prior-Z500_anl
#Z500_post_diff = Z500_post-Z500_anl
##T500_diff = T500_post-T500_prior
##T2M_diff = T2M_post-T2M_prior
#MSLP_prior_diff = MSLP_prior-MSLP_anl
#MSLP_post_diff = MSLP_post-Z500_anl
##PWAT_diff = PWAT_post-PWAT_prior

# Get the prior and posterior ensemble means for Z500
Z500_prior_ensmean = np.mean(Z500_prior, axis=3)
Z500_post_ensmean = np.mean(Z500_post, axis=3)

# Calculate the difference between Z500 prior/post ensmean and analysis
Z500_prior_ensmean_diff = Z500_prior_ensmean-Z500_anl
Z500_post_ensmean_diff = Z500_post_ensmean-Z500_anl
#Z500_ensdiff = {}
#for i in range(0,21):
#    Z500_ensdiff[i] = Z500_post[40,:,:,i]-Z500_prior[40,:,:,i]

# Get the prior and posterior ensemble means for MSLP
MSLP_prior_ensmean = np.mean(MSLP_prior, axis=3)
MSLP_post_ensmean = np.mean(MSLP_post, axis=3)

# Calculate the difference between MSLP prior/post ensmean and analysis
MSLP_prior_ensmean_diff = MSLP_prior_ensmean-MSLP_anl
MSLP_post_ensmean_diff = MSLP_post_ensmean-MSLP_anl

# Get the prior and posterior ensemble means for PWAT
#PWAT_prior_ensmean = np.mean(PWAT_prior, axis=3)
#PWAT_post_ensmean = np.mean(PWAT_post, axis=3)
#PWAT_ensmean_diff = PWAT_post_ensmean-PWAT_prior_ensmean


# The average absolute error difference (with specific lat/lon areas)
error_prior_Z500 = np.mean(abs(Z500_prior_ensmean_diff[:,100:150,240:330]))
error_post_Z500 = np.mean(abs(Z500_post_ensmean_diff[:,100:150,240:330]))



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
cs = map.contourf(x,y,Z500_prior_ensmean_diff[time_ind,:,:])
# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(Z500_units)

plt.title('{}'.format(ftimes[time_ind]))
plt.show()