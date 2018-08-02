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

ens = 'ecmwf' #ncep, eccc, ecmwf
loc_rad = '1000'
forecast_time = datetime(2015,11,10,12) # when the forecast was initialized
analysis_time = datetime(2015,11,11,12) # when to compare with analysis
oberrvar = ['0-1','1','10','100']
#oberrvar = [1,10,100,1000, 'ensvar', 250, 500, 750]
varlist = ['TCW','TCW','TCW','TCW'] #used only when loading ob update self
efh = '54hrs'
grid = [-180,180,90,0,3]
prior_var = ['TCW']#['QF850','D-QF850']#

obd = 'self' #ob update 'all' or 'self'? 

#variable we want to look at
vrbl= 'TCW'#'QF850'#
#what observation we used to update the ensemble
ob_type = 'TCW' #used if ob update is 'all'

#plotting variables
s = 35
n = 50
w = -180
e = -115
lat_ts = 40
figsize1 = 18
figsize2 = 12

if vrbl == 'QF850':
    cont_int = range(0,300,20)
    varstring = 'moisture flux'
    varunit = '(g/kg*m/s)'
elif vrbl == 'TCW':
    cont_int = range(0,60,3)
    varstring = 'total column water'
    varunit = '(mm)'
elif vrbl == 'IVT':
    cont_int = range(0,1000,50)
    varstring = 'integrated vapor flux'
    varunit = '(kg/m/s)'
elif vrbl == 'IWV':
    cont_int = range(0,50,2)
    varstring = 'integrated water vapor'
    varunit = '(mm)'



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

#convert to grid indices, add 1 to the right endpoints so python includes the last index
l = int(abs(-180-w)*2)
r = int(abs(-180-e)*2)+1
t = int((90-n)*2)
b = int((90-s)*2)+1

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
anl_mean_region = anl_mean[t:b,l:r]
prior_mean = np.mean(prior[timeind,:,:],axis=-1)
prior_mean_region = prior_mean[t:b,l:r]

#error of prior
prior_error = prior_mean-anl_mean
prior_error_region = prior_error[t:b,l:r]

#mean absolute error
weights = np.cos(np.radians(lats))
mae_prior = np.mean(np.average(abs(prior_error),axis=0,weights=weights))
#mean absolute error of region in plot
weights_region = np.cos(np.radians(lats[t:b]))
mae_prior_region = np.mean(np.average(abs(prior_error[t:b,l:r]),axis=0,weights=weights_region))
print('Mean absolute error of the prior forecast: ',mae_prior)
print('Mean absolute error of the prior forecast in the plotted region: ',mae_prior_region)

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
#x, y = m1(lon, lat)
x, y = m1(lon[t:b,l:r],lat[t:b,l:r])

cs1 = m1.contourf(x,y,anl_mean_region,cont_int)
cs3 = m1.contour(x,y,anl_mean_region,cont_int,colors='white',linewidths=0.5)
cs2 = m1.contour(x,y,prior_mean_region,cont_int,colors='black',linewidths=1)
plt.clabel(cs3, inline=0, fontsize=10, colors='white',fmt='%.0f')
plt.clabel(cs2, inline=0, fontsize=12,fmt='%.0f')
cbar1 = m1.colorbar(cs1, location='right',pad="3%")
cbar1.set_label(varstring+' '+varunit,fontsize=12)
plt.title(ens_dict[ens]+' analysis (colorfill) and prior forecast (black contours) '+am+'/'+ad+'/'+ay+' '+ah+'Z, forecast initialized '+fm+'/'+fd+'/'+fy+' '+fh+'Z')


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

cs1 = m1.contourf(x,y,prior_error_region,20)
cs2 = m1.contour(x,y,prior_mean_region,cont_int,colors='black',linewidths=1)
cs3 = m1.contour(x,y,prior_error_region,[0],colors='white',linewidths=1)
plt.clabel(cs2, inline=0, fontsize=12,fmt='%.0f')#,cmap=plt.cm.Reds)
cbar1 = m1.colorbar(cs1, location='right',pad="3%")
cbar1.set_label(varstring+' error '+varunit,fontsize=12)
plt.title(ens_dict[ens]+' prior (black) and error (colorfill) of '+am+'/'+ad+'/'+ay+' '+ah+'Z forecast, initialized '+fm+'/'+fd+'/'+fy+' '+fh+'Z')

#-----------------------------------------------------------------------------#


#loop through each observation error variance
for oev in oberrvar:    
    vrbl2 = vrbl
    if obd == 'all':
        #one observation error variance per file
        post_varstr = ob_type+str(oev)
    elif obd == 'self':
        post_varstr = ef.var_num_string(varlist, oberrvar)
        #if self-updating, add ob error variance to variable name for access in netCDF
        vrbl2 = vrbl+str(oev)
    # Filepath of the posterior forecast
    post_path = '/home/disk/hot/stangen/Documents/posterior_ensembles/gridded/ob_update_'+obd+'/inf_none/loc_'+loc_rad+'/'+ens+'/'+fy+fm+'/'+fy+'-'+fm+'-'+fd+'_'+fh+'_'+efh+'_'+grid_str+'_'+post_varstr+'.nc'
    
    print('loading posterior netCDF: '+fy+fm+fd+fh+' '+str(oev)+' ob err var')
    # Load the posterior data
    with Dataset(post_path, 'r') as ncdata:
        post = ncdata.variables[vrbl2][:]
        times2 = ncdata.variables['time']
    
    post_mean = np.mean(post[timeind,:,:],axis=-1)  
    post_mean_region = post_mean[t:b,l:r]
    post_error = post_mean-anl_mean
    post_error_region = post_error[t:b,l:r]
    
    #change in absolute error
    error_diff = abs(post_error)-abs(prior_error)
    error_diff_region = error_diff[t:b,l:r]
    #difference between posterior and prior (change in forecast)
    prior_post_diff = post_mean-prior_mean
    prior_post_diff_region = prior_post_diff[t:b,l:r]
    
    mae_post = np.mean(np.average(abs(post_error),axis=0,weights=weights))
    mae_post_region = np.mean(np.average(abs(post_error_region),axis=0,weights=weights_region))
    print('Mean absolute error of the posterior forecast: ',mae_post)
    print('Mean absolute error of the posterior forecast in the plotted region: ',mae_post_region)
    print('Change in absolute error: ',mae_post-mae_prior)
    print('Change in absolute error: ',mae_post_region-mae_prior_region)
    
    #set the vmin and vmax so that the divergent colorbar of change in error is centered around 0
    #get just the region we're interested in
    error_diff_region = error_diff[t:b,l:r]
    #find largest absolute value of difference in error
    vmax = np.max(abs(error_diff_region))
    vmin = 0-vmax
      
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
    cs1 = m1.contourf(x,y,anl_mean_region,cont_int,cmap=plt.cm.jet)
    cs2 = m1.contour(x,y,post_mean_region,cont_int,colors='black',linewidths=1)
    plt.clabel(cs2, inline=0, fontsize=12,fmt='%.0f')#,cmap=plt.cm.Reds)
    cs3 = m1.contour(x,y,prior_mean_region,cont_int,colors='black',linestyles='dashed',linewidths=1)
    plt.clabel(cs3, inline=0, fontsize=12,fmt='%.0f')
    cbar1 = m1.colorbar(cs1, location='right',pad="3%")
    cbar1.set_label('analysis '+varstring+' '+varunit,fontsize=12)
    plt.title(ens_dict[ens]+' analysis (colorfill), prior (dashed), posterior (solid) of '+am+'/'+ad+'/'+ay+' '+ah+'Z ob err var '+str(oev)+' forecast, initialized on '+fm+'/'+fd+'/'+fy+' '+fh+'Z')

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

    cs1 = m1.contourf(x,y,error_diff_region,20,vmin=vmin,vmax=vmax,cmap=plt.cm.seismic) 
    #cs2 = m1.contour(x,y,error_diff,colors='white',linewidths=1)
    cs4 = m1.contour(x,y,error_diff_region,[0],colors='white',linewidths=2)
    cs3 = m1.contour(x,y,prior_post_diff_region,10,colors='black',linewidths=1)

    #plt.clabel(cs2, inline=0, fontsize=12, fmt='%.0f')
    plt.clabel(cs3, inline=0, fontsize=12)#,fmt='%.0f')
    #cs1 = m1.contourf(x,y,prior_post_diff)     
    cbar1 = m1.colorbar(cs1, location='right',pad="3%")
    cbar1.set_label(varstring+' error change '+varunit,fontsize=12)
    plt.title(ens_dict[ens]+' change in '+varstring+' (black) and change in absolute error (colorfill) of '+am+'/'+ad+'/'+ay+' '+ah+'Z forecast, '+str(oev)+' ob err var , initialized on '+fm+'/'+fd+'/'+fy+' '+fh+'Z')
    
