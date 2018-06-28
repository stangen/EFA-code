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
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time


f = open('/home/disk/hot/stangen/Documents/surface_obs/MADIS/201304/combined_T2M/T2M_20130410_0600.txt','r')
lines = f.readlines()
f.close()

name = []; alts = [];lts = []; lngs = []; tms=[];
for l in lines:
    temp = l.split(',')
    name.append(temp[0])
    lts.append(temp[1])
    lngs.append(temp[2])
    tms.append(temp[4])
    alts.append(temp[5])

lts = np.float64(lts)
lngs = np.float64(lngs)   
tms = np.float64(tms) 
alts = np.float64(alts)
#set min/max for colorbar
#minn = round(min(alts)-2,0)
#maxx = round(max(alts)+2,0)
minn = 265
maxx = 305

# Filepath of the prior forecast
prior_path = '/home/disk/hot/stangen/Documents/prior_ensembles/ecmwf/201304/2013-04-10_00_ecmwf_T2M_ALT.nc'

# Filepath of the posterior forecast
post_path = '/home/disk/hot/stangen/Documents/posterior_ensembles/ob_update_self/loc_1000/ecmwf/201304/2013-04-10_00_ecmwf_T2M_ALT.nc'
 
# Filepath of the analysis at the desired time 
analysis_path = '/home/disk/hot/stangen/Documents/prior_ensembles/ecmwf/201304/2013-04-10_12_ecmwf_T2M_ALT.nc'             
# only variable I am using is 500mb height            
vrbls=['T2M']

# Load the prior data            
with Dataset(prior_path,'r') as ncdata:
    #print(ncdata.variables)
    times = ncdata.variables['time']
    ftimes = num2date(times[:],times.units)
    lats = ncdata.variables['lat'][:]
    lons = ncdata.variables['lon'][:]
    nens = len(ncdata.variables['ens'][:])
    for vrbl in vrbls:
        prior = ncdata.variables[vrbl][:]
        prior_units = ncdata.variables[vrbl].units


# Load the posterior data
with Dataset(post_path, 'r') as ncdata:
    post = ncdata.variables[vrbls[0]][:] 
    times2 = ncdata.variables['time']

# Load the analysis data
with Dataset(analysis_path, 'r') as ncdata:
    analysis = ncdata.variables[vrbls[0]][:] 
    times3 = ncdata.variables['time']
    atimes = num2date(times3[:],times3.units)

#get the ensemble mean and variance of the prior, posterior, and analysis
prior_mean = np.mean(prior,axis=3)    
post_mean = np.mean(post,axis=3)
analysis_mean = np.mean(analysis,axis=3)
prior_variance = np.var(prior,axis=3,ddof=1)
post_variance = np.var(post,axis=3,ddof=1)
analysis_variance=np.var(analysis,axis=3,ddof=1)
#differences between combos of prior, post, and analysis ens means at a certain time
post_prior_diff= post_mean[1,:,:]-prior_mean[1,:,:]
anal_prior_diff = prior_mean[1,:,:]-analysis_mean[0,:,:]
anal_post_diff = post_mean[1,:,:]-analysis_mean[0,:,:]

#calculate the covariance of the ensemble with a certain point (KMKN, ~32N, 98.5W)
prior_perts = prior-prior_mean[:,:,:,None]
point_prior = prior[1,116,163,:]
point2_prior = prior[1,115,162,:]
point_prior_perts = point_prior-np.mean(point_prior)
point_prior_var = np.var(prior[1,116,163,:],ddof=1)
cov_point = np.dot(prior_perts,point_prior_perts)/(nens-1)
point_cov = np.dot(point_prior,point2_prior)/(nens-1)

fig = plt.figure(figsize=(18,12))
ax1 = plt.subplot(111)
# Plot the difference between the prior and posterior
#map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
m1 = Basemap(projection='merc',llcrnrlat=27,urcrnrlat=40,\
            llcrnrlon=-103,urcrnrlon=-90,lat_ts=30,resolution='c')
# draw coastlines, country boundaries, fill continents.
m1.drawcoastlines(linewidth=1.25)
m1.drawcountries(linewidth=1.25)
m1.drawstates()

#Convert lat/lon to mercator coordinates
xalti,yalti = m1(lngs,lts)

#map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
m1.drawmapboundary(fill_color='aqua')
#map.fillcontinents(color='coral',lake_color='aqua')
# draw lat/lon grid lines every 30 degrees.
m1.drawmeridians(np.arange(0,360,.5))
m1.drawparallels(np.arange(-90,90,.5))

# compute native map projection coordinates of lat/lon grid.
lon, lat = np.meshgrid(lons, lats)
x, y = m1(lon, lat)

time_ind = 2
#cs = map.contourf(x,y,Z500_diff[0,:,:,0])#,levels=np.arange(-90, 91, 30), cmap=plt.cm.RdBu_r,linewidths=1.5, extend='both')
cs = m1.contour(x,y,prior_mean[1,:,:],20,vmin=minn, vmax=maxx)
cs1 = m1.contourf(x,y,post_mean[1,:,:],20,vmin=minn, vmax=maxx)
plt.clabel(cs, inline=1, fontsize=10)
#cs1 = m1.contour(x,y,post_variance[2,:,:],[286],vmin=minn, vmax=maxx)
#cs3 = m1.contour(x,y,analysis_mean[0,:,:],[286],vmin=minn, vmax=maxx)
#plt.clabel(cs3, inline=1, fontsize=10)
#cs2 = m1.scatter(xalti,yalti)
for i, txt in enumerate(alts):
    ax1.annotate(txt, (xalti[i],yalti[i]),clip_on=True)
cs2= m1.scatter(xalti,yalti,c=alts,zorder=4,s=50,vmin=minn,vmax=maxx,marker="o",edgecolor='k',cmap=plt.cm.gist_rainbow)

# Add Colorbar
cbar = m1.colorbar(cs1, location='bottom', pad="10%")
cbar2 = m1.colorbar(cs2,fraction=0.023)
cbar.set_label(prior_units)

plt.title('{}'.format(ftimes[time_ind]))
plt.show()
    