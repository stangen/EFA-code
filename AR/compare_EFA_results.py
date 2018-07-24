#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:24:41 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
import numpy as np
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time
import EFA.duplicate_madaus.efa_functions as ef

ens = 'ncep'
loc_rad = '1000'
forecast_time = datetime(2015,11,12,0) # when the forecast was initialized
analysis_time = datetime(2015,11,13,0) # when to compare with analysis
oberrvar = [1,10,100,1000, 'ensvarpractice']
efh = '54hrs'
grid = [-180,180,90,0,3]
prior_var = ['QF850','D-QF850']
#plotting variables
s = 35
n = 50
w = -180
e = -115
lat_ts = 40
figsize1 = 18
figsize2 = 12

cont_int = range(0,300,20)

post_var = ['QF850']
vrbl= 'QF850'

point_lon = -99#163
point_lat = 32#116

#get the times for the analysis
ay = analysis_time.strftime('%Y')
am = analysis_time.strftime('%m')
ad = analysis_time.strftime('%d')
ah = analysis_time.strftime('%H')

#get the times for the forecast
fy = forecast_time.strftime('%Y')
fm = forecast_time.strftime('%m')
fd = forecast_time.strftime('%d')
fh = forecast_time.strftime('%H')

#convert lists of variables into strings for loading files
prior_var_str = ef.var_string(prior_var)
grid_str = ef.var_string(grid)

#get forecast index- i.e how many timesteps after initialization we are looking at
timediff = analysis_time-forecast_time
tdd = timediff.days
tds = timediff.seconds
tdh = tdd*24+tds/3600
#timesteps of forecast output are 6 hours apart
timeind = int(tdh/6)

#dictionary for ensemble types, for use in plotting
ens_dict = {'ncep': 'GEFS',
            'eccc': 'CMC',
            'ecmwf': 'ECMWF'}

# Filepath of the analysis at the desired time 
analysis_path = '/home/disk/hot/stangen/Documents/prior_ensembles/'+ens+'/'+ay+am+'/'+ay+'-'+am+'-'+ad+'_'+ah+'_'+efh+'_'+prior_var_str+'.nc'             
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
    
# Load the analysis data
print('loading analysis netCDF: '+ay+am+ad+ah)
with Dataset(analysis_path, 'r') as ncdata:
    analysis = ncdata.variables[vrbl][:]
    times3 = ncdata.variables['time']
    atimes = num2date(times3[:],times3.units)
 
#means of analysis and prior, at analysis time
anl_mean = np.mean(analysis[0,:,:],axis=-1)
prior_mean = np.mean(prior[timeind,:,:],axis=-1)

#error of prior
prior_error = prior_mean-anl_mean

#-----------------------------------------------------------------------------#
#plot the analysis, with the prior as contours
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
x, y = m1(lon, lat)

cs1 = m1.contourf(x,y,anl_mean,20)
cs2 = m1.contour(x,y,prior_mean,cont_int,colors='black',linewidths=1)
plt.clabel(cs2, inline=0, fontsize=12,fmt='%.0f')
cbar1 = m1.colorbar(cs1, location='right',pad="3%")
cbar1.set_label('Moisture Flux (g/kg*m/s)',fontsize=12)
plt.title(ens_dict[ens]+' analysis and prior forecast on '+am+'/'+ad+'/'+ay+' '+ah+'Z, forecast initialized on '+fm+'/'+fd+'/'+fy+' '+fh+'Z')

#-----------------------------------------------------------------------------#
#plot the prior forecast error, with contour as the forecast
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

cs1 = m1.contourf(x,y,prior_error,20)
cs2 = m1.contour(x,y,prior_mean,cont_int,colors='black',linewidths=1)
cs3 = m1.contour(x,y,prior_error,[0],colors='white',linewidths=1)
plt.clabel(cs2, inline=0, fontsize=12,fmt='%.0f')#,cmap=plt.cm.Reds)
cbar1 = m1.colorbar(cs1, location='right',pad="3%")
cbar1.set_label('Moisture Flux (g/kg*m/s)',fontsize=12)
plt.title(ens_dict[ens]+' error of '+am+'/'+ad+'/'+ay+' '+ah+'Z prior forecast, initialized on '+fm+'/'+fd+'/'+fy+' '+fh+'Z')

#-----------------------------------------------------------------------------#


#loop through each observation error variance
for oev in oberrvar:
    post_varstr = vrbl+str(oev)
    # Filepath of the posterior forecast
    post_path = '/home/disk/hot/stangen/Documents/posterior_ensembles/gridded/ob_update_all/inf_none/loc_'+loc_rad+'/'+ens+'/'+fy+fm+'/'+fy+'-'+fm+'-'+fd+'_'+fh+'_'+efh+'_'+grid_str+'_'+post_varstr+'.nc'
    
    print('loading posterior netCDF: '+fy+fm+fd+fh+' '+str(oev)+' ob err var')
    # Load the posterior data
    with Dataset(post_path, 'r') as ncdata:
        post = ncdata.variables[vrbl][:]
        times2 = ncdata.variables['time']
    
    post_mean = np.mean(post[timeind,:,:],axis=-1)    
    post_error = post_mean-anl_mean
    
    #change in absolute error
    error_diff = abs(post_error)-abs(prior_error)
    prior_post_diff = post_mean-prior_mean
      
#-----------------------------------------------------------------------------#
    #plot the posterior forecast error, with contour as the forecast
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
    
    #cs1 = m1.contourf(x,y,post_error)
    cs1 = m1.contourf(x,y,anl_mean,20)
    cs2 = m1.contour(x,y,post_mean,cont_int,colors='black',linewidths=1)
    plt.clabel(cs2, inline=0, fontsize=12,fmt='%.0f')#,cmap=plt.cm.Reds)
    cs3 = m1.contour(x,y,prior_mean,cont_int,colors='grey',linewidths=1)
    plt.clabel(cs3, inline=0, fontsize=12,fmt='%.0f')
    cbar1 = m1.colorbar(cs1, location='right',pad="3%")
    cbar1.set_label('Moisture Flux (g/kg*m/s)',fontsize=12)
    plt.title(ens_dict[ens]+' analysis (colorfill), post (black), prior (grey) of '+am+'/'+ad+'/'+ay+' '+ah+'Z ob err var '+str(oev)+' forecast, initialized on '+fm+'/'+fd+'/'+fy+' '+fh+'Z')

#-----------------------------------------------------------------------------#
    #plot the difference between prior and posterior forecasts
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
    
    cs1 = m1.contourf(x,y,prior_post_diff)
    cbar1 = m1.colorbar(cs1, location='right',pad="3%")
    cbar1.set_label('Moisture Flux (g/kg*m/s)',fontsize=12)    
    plt.title(ens_dict[ens]+' posterior minus prior, '+am+'/'+ad+'/'+ay+' '+ah+'Z forecast, '+str(oev)+' ob err var , initialized on '+fm+'/'+fd+'/'+fy+' '+fh+'Z')

#-----------------------------------------------------------------------------#
    #plot the change in forecast error (error posterior - error prior)
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

    cs1 = m1.contourf(x,y,error_diff,20) 
    cs4 = m1.contour(x,y,error_diff,[0],colors='white',linewidths=1)
    #cs1 = m1.contourf(x,y,prior_post_diff)     
    cbar1 = m1.colorbar(cs1, location='right',pad="3%")
    cbar1.set_label('Moisture Flux (g/kg*m/s)',fontsize=12)
    plt.title(ens_dict[ens]+' change in absolute error of '+am+'/'+ad+'/'+ay+' '+ah+'Z forecast, '+str(oev)+' ob err var , initialized on '+fm+'/'+fd+'/'+fy+' '+fh+'Z')
