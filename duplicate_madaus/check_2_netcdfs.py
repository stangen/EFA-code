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
import efa_functions as ef

ens = 'ecmwf'
loc_rad = '100'
base_time = datetime(2013,4,10,0)
time_ind = 1
prior_var = ['T2M','ALT']
post_var = ['T2M']
vrbl= 'T2M'
point_lon = -99#163
point_lat = 32#116

#number of Hours After Initialization- corresponds with the analysis time
hai = str(time_ind*6)

savedir = '/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/plots/'

assim_time = base_time + timedelta(hours=6)

analysis_time = base_time + timedelta(hours=6*time_ind)

bty = base_time.strftime('%Y')
btm = base_time.strftime('%m')
btd = base_time.strftime('%d')
bth = base_time.strftime('%H')

bstr = bty+btm+btd+bth

aty = analysis_time.strftime('%Y')
atm = analysis_time.strftime('%m')
atd = analysis_time.strftime('%d')
ath = analysis_time.strftime('%H')

aaty = assim_time.strftime('%Y')
aatm = assim_time.strftime('%m')
aatd = assim_time.strftime('%d')
aath = assim_time.strftime('%H')

prior_var_string = ef.var_string(prior_var)
post_var_string = ef.var_string(post_var)


# Filepath of the prior forecast
prior_path = '/home/disk/hot/stangen/Documents/prior_ensembles/'+ens+'/'+bty+btm+'/'+bty+'-'+btm+'-'+btd+'_'+bth+'_'+ens+'_'+prior_var_string+'.nc'

# Filepath of the posterior forecast
post_path = '/home/disk/hot/stangen/Documents/posterior_ensembles/ob_update_self/loc_'+loc_rad+'/'+ens+'/'+bty+btm+'/'+bty+'T2M-'+btm+'-'+btd+'_'+bth+'_'+ens+'_'+post_var_string+'.nc'
 
# Filepath of the analysis at the desired time 
analysis_path = '/home/disk/hot/stangen/Documents/prior_ensembles/'+ens+'/'+aty+atm+'/'+aty+'-'+atm+'-'+atd+'_'+ath+'_'+ens+'_'+prior_var_string+'.nc'             
# only variable I am using is 500mb height            

#load the observations 6 hours in to the forecast (the assimilated observations)
f = open('/home/disk/hot/stangen/Documents/surface_obs/MADIS/'+aaty+aatm+'/combined_'+vrbl+'/'+vrbl+'_'+aaty+aatm+aatd+'_'+aath+'00.txt','r')
obs_assim = f.readlines()
f.close()

name = []; ob_var = [];lts = []; lngs = []; tms=[];
for l in obs_assim:
    temp = l.split(',')
    name.append(temp[0])
    lts.append(temp[1])
    lngs.append(temp[2])
    tms.append(temp[4])
    ob_var.append(temp[5])
    
#load the observations for our analysis time (i.e. 6, 12, 18, etc hours into the forecast)
f = open('/home/disk/hot/stangen/Documents/surface_obs/MADIS/'+aty+atm+'/combined_'+vrbl+'/'+vrbl+'_'+aty+atm+atd+'_'+ath+'00.txt','r')
obs_analysis = f.readlines()
f.close()

anl_name = []; anl_ob_var = [];anl_lts = []; anl_lngs = []; anl_tms=[];
for l in obs_analysis:
    temp = l.split(',')
    anl_name.append(temp[0])
    anl_lts.append(temp[1])
    anl_lngs.append(temp[2])
    anl_tms.append(temp[4])
    anl_ob_var.append(temp[5])

#convert to numpy arrays
lts = np.float64(lts)
lngs = np.float64(lngs)   
tms = np.float64(tms) 
ob_var = np.float64(ob_var)
anl_lts = np.float64(anl_lts)
anl_lngs = np.float64(anl_lngs)   
anl_tms = np.float64(anl_tms) 
anl_ob_var = np.float64(anl_ob_var)

#set min/max for colorbar
#minn = round(min(alts)-2,0)
#maxx = round(max(alts)+2,0)
minn = 265
maxx = 305

# Load the prior data            
with Dataset(prior_path,'r') as ncdata:
    #print(ncdata.variables)
    times = ncdata.variables['time']
    ftimes = num2date(times[:],times.units)
    lats = ncdata.variables['lat'][:]
    lons = ncdata.variables['lon'][:]
    nens = len(ncdata.variables['ens'][:])
    prior = ncdata.variables[vrbl][:]
    prior_units = ncdata.variables[vrbl].units


# Load the posterior data
with Dataset(post_path, 'r') as ncdata:
    post = ncdata.variables[vrbl][:]
    times2 = ncdata.variables['time']
    
#get the ensemble mean and variance of the prior, posterior, and analysis
prior_mean = np.mean(prior,axis=3)    
post_mean = np.mean(post,axis=3)

prior_variance = np.var(prior,axis=3,ddof=1)
post_variance = np.var(post,axis=3,ddof=1)
#differences between combos of prior, post, and analysis ens means at a certain time
post_prior_diff= post_mean[time_ind,:,:]-prior_mean[time_ind,:,:]

if time_ind != 1:
    # Load the analysis data
    with Dataset(analysis_path, 'r') as ncdata:
        analysis = ncdata.variables[vrbl][:]
        times3 = ncdata.variables['time']
        atimes = num2date(times3[:],times3.units)

    analysis_mean = np.mean(analysis,axis=3)
    analysis_variance=np.var(analysis,axis=3,ddof=1)
    anal_prior_diff = prior_mean[time_ind,:,:]-analysis_mean[0,:,:]
    anal_post_diff = post_mean[time_ind,:,:]-analysis_mean[0,:,:]



#calculate the covariance of the ensemble with a certain point (KMKN, ~32N, 98.5W)

yy = int((90-point_lat)*2)
xx = int(abs(-180-point_lon)*2)
prior_perts = prior-prior_mean[:,:,:,None]
point_prior = prior[1,yy,xx,:]
point2_prior = prior[1,115,162,:]
point_prior_perts = point_prior-np.mean(point_prior)
point_prior_var = np.var(point_prior_perts,ddof=1)
cov_point = np.dot(prior_perts,point_prior_perts)/(nens-1)
point_cov = np.dot(point_prior,point2_prior)/(nens-1)

#kalman gain matrix
kal_gain = cov_point/(point_prior_var+1)

for j in range(0,5):   

    fig = plt.figure(figsize=(18,12))
    ax1 = fig.add_subplot(111)
    # Plot the difference between the prior and posterior
    #map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
    m1 = Basemap(projection='merc',llcrnrlat=27,urcrnrlat=38,\
                llcrnrlon=-103,urcrnrlon=-92,lat_ts=30,resolution='c')
    # draw coastlines, country boundaries, fill continents.
    m1.drawcoastlines(linewidth=1.25)
    m1.drawcountries(linewidth=1.25)
    m1.drawstates()
    
    #Convert lat/lon to mercator coordinates
    xassimob,yassimob = m1(lngs,lts)
    xanlob,yanlob = m1(anl_lngs,anl_lts)
    
    xpoint,ypoint=m1(point_lon,point_lat)
    #map.fillcontinents(color='coral',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    m1.drawmapboundary(fill_color='aqua')
    
    #map.fillcontinents(color='coral',lake_color='aqua')
    # draw lat/lon grid lines every half degree.
    #labels[left,right,top,bottom]  1=True 0=False
    m1.drawmeridians(np.arange(0,360,1),labels=[1,0,0,1],fontsize=10,linewidth=0.5)
    m1.drawparallels(np.arange(-90,90,1),labels=[1,0,0,1],fontsize=10,linewidth=0.5)
    
    # compute native map projection coordinates of lat/lon grid.
    lon, lat = np.meshgrid(lons, lats)
    x, y = m1(lon, lat)
    

    #prior variance, prior mean values, and observation values 6 hours into forecast
    if j == 0:
        #cs = m1.contourf(x,y,post_mean[time_ind,:,:],20,vmin=minn, vmax=maxx)
        cs = m1.contour(x,y,prior_mean[1,:,:],20,vmin=minn, vmax=maxx,cmap=plt.cm.gist_rainbow)
        plt.clabel(cs, inline=0, fontsize=12,fmt='%.1f')#,cmap=plt.cm.Reds)
        cs1 = m1.contourf(x,y,prior_variance[1,:,:],20,cmap=plt.cm.YlOrRd)
        #cbar = m1.colorbar(cs, location='bottom', pad="10%")   

        for i, txt in enumerate(ob_var):
            ax1.annotate(txt, (xassimob[i],yassimob[i]),clip_on=True)
        cs2= m1.scatter(xassimob,yassimob,c=ob_var,zorder=4,s=50,vmin=minn,vmax=maxx,marker="o",edgecolor='k',cmap=plt.cm.gist_rainbow)        
        cbar2 = m1.colorbar(cs2,fraction=0.023) 
        cbar2.set_label('ob temperature (K)',fontsize=12)
        cbar1 = m1.colorbar(cs1, location='bottom',pad="3%")
        cbar1.set_label('Variance (hPa$^{2}$)',fontsize=12)
        plt.title('Prior variance, prior mean, and obs at time of assim.',weight='bold',fontsize=16, y=1.01)       
        plt.show()
        
        fig.savefig(savedir+bstr+'_6hr_var_prior_mean_obs.png',frameon=False,bbox_inches='tight')
    
    #Kalman gain at analysis time from an "ob" at one ensemble point taken 6 hours into forecast, 
    #prior mean values, and ob values 6 hours into forecast
    if j ==1:
        #cs = map.contourf(x,y,Z500_diff[0,:,:,0])#,levels=np.arange(-90, 91, 30), cmap=plt.cm.RdBu_r,linewidths=1.5, extend='both')
        cs = m1.contour(x,y,prior_mean[1,:,:],20,vmin=minn, vmax=maxx,cmap=plt.cm.gist_rainbow)#,cmap=plt.cm.RdBu_r,) #contours of prior mean
        #cs1 = m1.contourf(x,y,post_mean[1,:,:],20,vmin=minn, vmax=maxx) #color fill of post mean
        plt.clabel(cs, inline=0, fontsize=12,fmt='%.1f')
        #cs1 = m1.contourf(x,y,prior_variance[1,:,:])
        #cs1 = m1.contourf(x,y,cov_point[1,:,:])#,vmin=minn, vmax=maxx)
        cs1 = m1.contourf(x,y,kal_gain[time_ind,:,:])
        #cs3 = m1.contour(x,y,analysis_mean[0,:,:],[286],vmin=minn, vmax=maxx)
        #plt.clabel(cs3, inline=1, fontsize=10)
        #cs2 = m1.scatter(xalti,yalti)
        m1.scatter(xpoint,ypoint,s=100, zorder=4,c='black',marker='o',edgecolor='k')
        for i, txt in enumerate(ob_var):
            ax1.annotate(txt, (xassimob[i],yassimob[i]),clip_on=True)
        cs2= m1.scatter(xassimob,yassimob,c=ob_var,zorder=4,s=50,vmin=minn,vmax=maxx,marker="o",edgecolor='k',cmap=plt.cm.gist_rainbow)
        # Add Colorbar
              
        #cbar = m1.colorbar(cs, location='bottom', pad="10%")
        cbar2 = m1.colorbar(cs2,fraction=0.023)
        cbar2.set_label('ob temperature (K)',fontsize=12)
        cbar = m1.colorbar(cs1, location='bottom',pad='3%') 
        cbar.set_label('Kalman Gain',fontsize=12)         
        plt.title(str(ftimes[time_ind])+' Kalman gain & innovation',weight='bold',fontsize=16, y=1.01)
        plt.show()
        
        fig.savefig(savedir+bstr+'_'+hai+'hr_Kal_gain_innov.png',frameon=False,bbox_inches='tight')

    #Posterior ensemble mean and observation values at analysis time  
    if j == 2:
        m1.fillcontinents(color='coral',lake_color='aqua')
        #cs = m1.contourf(x,y,post_mean[time_ind,:,:],20,vmin=minn, vmax=maxx)
        cs = m1.contour(x,y,post_mean[time_ind,:,:],20,vmin=minn, vmax=maxx,cmap=plt.cm.gist_rainbow)
        plt.clabel(cs, inline=0, fontsize=12,fmt='%.1f')
        for i, txt in enumerate(anl_ob_var):
            ax1.annotate(txt, (xanlob[i],yanlob[i]),clip_on=True)
        cs2= m1.scatter(xanlob,yanlob,c=anl_ob_var,zorder=4,s=50,vmin=minn,vmax=maxx,marker="o",edgecolor='k',cmap=plt.cm.gist_rainbow)
        #cbar = m1.colorbar(cs, location='bottom', pad="10%")   
        cbar2 = m1.colorbar(cs2,fraction=0.023)
        cbar2.set_label('ob temperature (K)',fontsize=12)
        plt.title(str(ftimes[time_ind])+' Posterior ensemble mean & obs',weight='bold',fontsize=16, y=1.01)
        plt.show()
        
        fig.savefig(savedir+bstr+'_'+hai+'hr_loc'+loc_rad+'_post_mean_obs.png',frameon=False,bbox_inches='tight')
    
    #Change in forecast after assimilation at analysis time: posterior - prior    
    if j ==3:
        cs = m1.contourf(x,y,post_prior_diff[:,:])
        cbar = m1.colorbar(cs, location='bottom', pad="3%")
        cbar.set_label('Difference (K)',fontsize=12)
        plt.title(str(ftimes[time_ind])+' Change in ensemble mean after assimilation',weight='bold',fontsize=16, y=1.01)
        plt.show()
    
        fig.savefig(savedir+bstr+'_'+hai+'hr_loc'+loc_rad+'_change_ens_mean.png',frameon=False,bbox_inches='tight')
    #286 K contour of prior, posterior, and analysis ensemble means    
    if j == 4: 
        m1.fillcontinents(color='coral',lake_color='aqua')
        m1.contour(x,y,prior_mean[time_ind,:,:],[286], colors='r')
        m1.contour(x,y,post_mean[time_ind,:,:],[286],colors='k')
        if time_ind !=1:
            m1.contour(x,y,analysis_mean[0,:,:],[286],colors='b')
        for i, txt in enumerate(anl_ob_var):
            ax1.annotate(txt, (xanlob[i],yanlob[i]),clip_on=True)
        cs2= m1.scatter(xanlob,yanlob,c=anl_ob_var,zorder=4,s=50,vmin=minn,vmax=maxx,marker="o",edgecolor='k',cmap=plt.cm.gist_rainbow)
        cbar2 = m1.colorbar(cs2,fraction=0.023)
        cbar2.set_label('ob temperature (K)',fontsize=12)
        plt.title(str(ftimes[time_ind])+' 286K isotherms: red=prior, black=posterior',weight='bold',fontsize=16, y=1.01)
        plt.show()
        
        fig.savefig(savedir+bstr+'_'+hai+'hr_loc'+loc_rad+'_286K_isotherms.png',frameon=False,bbox_inches='tight')

    