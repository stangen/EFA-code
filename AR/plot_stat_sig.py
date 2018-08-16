#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 15:31:26 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
import numpy as np
from datetime import datetime, timedelta
import EFA.efa_files.cfs_utilities_st as ut
import surface_obs.madis_example.madis_utilities as mt
import os
import glob
import sys
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from EFA.efa_xray.observation.observation import Observation
from EFA.efa_xray.assimilation.ensrf import EnSRF
from EFA.duplicate_madaus.load_data import Load_Data
import EFA.duplicate_madaus.efa_functions as ef

"""
Script to plot where correlation between ob estimate and state
is statistically significant.
"""

#date/time of ensemble initialization
forecast_time = datetime(2015,11,11,0)
#ensemble type
ens = 'ecmwf'
#variable we want to look at
vrbl= 'IVT'
#latitude of ob location
ob_lat = 46
#longitude of ob location
ob_lon = -165
#condidence threshold
confidence = 90
#prior variables
prior_var = ['IWV','IVT','D-IVT']

#plotting variables
s = 35
n = 50
w = -180
e = -115
lat_ts = 40
figsize1 = 18
figsize2 = 12


#end forecast hour string
efh = '54hrs'

#convert to grid indices, add 1 to the right endpoints so python includes the last index
l = int(abs(-180-w)*2)
r = int(abs(-180-e)*2)+1
t = int((90-n)*2)
b = int((90-s)*2)+1

#convert lists of variables into strings for loading files
prior_var_str = ef.var_string(prior_var)

#get the times for the forecast
fy = forecast_time.strftime('%Y')
fm = forecast_time.strftime('%m')
fd = forecast_time.strftime('%d')
fh = forecast_time.strftime('%H')

# Filepath of the prior forecast
prior_path = '/home/disk/hot/stangen/Documents/prior_ensembles/'+ens+'/'+fy+fm+'/'+fy+'-'+fm+'-'+fd+'_'+fh+'_'+efh+'_'+prior_var_str+'.nc'

# Load the prior data 
print('loading prior netCDF: '+fy+fm+fd+fh)           
with Dataset(prior_path,'r') as ncdata:
    #print(ncdata.variables)
    times = ncdata.variables['time']
    ftimes = num2date(times[:],times.units)
    lats = ncdata.variables['lat'][:]
    lons = ncdata.variables['lon'][:]
    nens = len(ncdata.variables['ens'][:])
    prior = ncdata.variables[vrbl][:]
    prior_units = ncdata.variables[vrbl].units

#dummy elevations so closest_points function works    
elevs = np.zeros([len(lats),len(lons)])
#get observation time (12 hours after forecast)
ob_time = forecast_time + timedelta(seconds=3600*12)
#get ensemble estimate of observation (interp)    
interp, TorF = ef.closest_points(ob_lat,ob_lon,lats,lons,0,
                           elevs,ob_time,ftimes,prior,need_interp=True)

r_crit = ef.get_r_crit()
#get ranks/rank perts of ob estimate and state
ye_rank = np.argsort(interp)
ye_rank = np.argsort(ye_rank)

#time indices
time_ind = [2,4,6,8]

for i in time_ind:
    prior_mean = np.mean(prior[i,:],axis=-1)
    prior_mean_region = prior_mean[t:b,l:r]
    
    #sort once to obtain the indices of a sorted array
    prior_1sort = np.argsort(prior[i,:],axis=-1)
    #sort again to obtain the rank - double argsort is what causes code to run slowly
    prior_rank = np.argsort(prior_1sort,axis=-1)
    
    ye_rank_pert = ye_rank - np.mean(ye_rank)
    prior_rank_pert = prior_rank - np.mean(prior_rank,axis=-1)[:,:,None]
    #calculate covariance of ranks, then spearman correlation
    cov_rank = np.dot(prior_rank_pert,ye_rank_pert)/len(interp)
    spearmanr = cov_rank/(np.std(ye_rank_pert)*np.std(prior_rank_pert,axis=-1))
    #find which values are below the critical r value (and therefore fail)
    r_pass = abs(spearmanr) >= r_crit[len(interp)][confidence]
#        print(r_crit[len(ye)][self.localize_radius])
    r_pass_region = r_pass[t:b,l:r]
    #plot AR variable and where ob estimate is statistically significant
    
    fig = plt.figure(figsize=(figsize1,figsize2))
    ax1 = fig.add_subplot(111)
    # Plot the difference between the prior and posterior
    #map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
    m1 = Basemap(projection='merc',llcrnrlat=s,urcrnrlat=n,\
                llcrnrlon=w,urcrnrlon=e,lat_ts=lat_ts,resolution='c')
    # draw coastlines, country boundaries, fill continents.
    m1.drawcoastlines(linewidth=1.25)
    m1.drawcountries(linewidth=1.25)
    m1.drawstates()
    # draw lat/lon grid lines every half degree.
    #labels[left,right,top,bottom]  1=True 0=False
    m1.drawmeridians(np.arange(0,360,5),labels=[1,0,0,1],fontsize=10,linewidth=0.5)
    m1.drawparallels(np.arange(-90,90,1),labels=[1,0,0,1],fontsize=10,linewidth=0.5)
    # compute native map projection coordinates of lat/lon grid.
    lon, lat = np.meshgrid(lons, lats)
    #x, y = m1(lon, lat)
    x, y = m1(lon[t:b,l:r],lat[t:b,l:r])
    
    cs1 = m1.contourf(x,y,prior_mean_region)
    m1.contour(x,y,r_pass_region)
#    cs1 = m1.contourf(x,y,anl_mean_region,cont_int)
#    cs3 = m1.contour(x,y,anl_mean_region,cont_int,colors='white',linewidths=0.5)
#    cs2 = m1.contour(x,y,prior_mean_region,cont_int,colors='black',linewidths=1)
#    plt.clabel(cs3, inline=0, fontsize=10, colors='white',fmt='%.0f')
#    plt.clabel(cs2, inline=0, fontsize=12,fmt='%.0f')
#    cbar1 = m1.colorbar(cs1, location='right',pad="3%")
#    cbar1.set_label(varstring+' '+varunit,fontsize=12)
#    plt.title(ens_dict[ens]+' analysis (colorfill) and prior forecast (black contours) '+am+'/'+ad+'/'+ay+' '+ah+'Z, forecast initialized '+fm+'/'+fd+'/'+fy+' '+fh+'Z')

    
    
   # kcov[r_fail] = 0








    