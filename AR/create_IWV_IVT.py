#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 10:23:02 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
from datetime import datetime
import numpy as np
import EFA.efa_files.cfs_utilities_st as ut
import surface_obs.madis_example.madis_utilities as mt
import EFA.duplicate_madaus.efa_functions as ef
import os

# ecmwf, eccc for euro/canadian ensembles, ncep
ensemble_type = ['eccc']
#start and end date to get ensembles. 
start_date = datetime(2015,11,10,12) #YYYY,m,d,h
end_date = datetime(2015,11,10,12)
hourstep = 12 #how often you want a new forecast initialization, usually 12 hr

#variables with names coming from the raw TIGGE- see get_tigge_data if unsure of names.
#the order matters to make filename match exactly. 
surf_variables = ['D2M','SP','U10','V10']
upper_variables = ['Q','U','V']
levels = ['1000','925','850','700','500','300']
end_variables = ['IWV','IVT','D-IVT']

#a list of dates to loop through to load each forecast initialized on these dates
dates = mt.make_datetimelist(start_date,end_date,hourstep)  

surf_str = ef.var_string(surf_variables)+'_sfc.nc'
upper_str = ef.var_string(upper_variables)+'_'+ef.var_string(levels)+'_pl.nc'
nvars = len(end_variables)

g = 9.80665

outvar_string = ef.var_string(end_variables)

"""
This function creates a netCDF from a raw TIGGE netCDF. The main purpose
of this is to change surface pressure to altimeter setting, calculate
6-hourly precipitation, and to rename/shorten variable names. Unfortunately
using the float32 format for the variables uses twice the memory of the 
raw TIGGE int16 format. 
"""
for ens in ensemble_type:
    for date in dates:

        y,m,d,h = ef.dt_str_timedelta(date)
        
        print('Working on '+d+'_'+h+' '+ens)
        
        #rename TIGGE variable names
        vardict = {
                'T2M' : 't2m',
                'D2M' : 'd2m',
                'SP' : 'sp',
                'U10' : 'u10',
                'V10' : 'v10',
                'U' : 'u',
                'V' : 'v',
                'Q' : 'q',
                'ALT' : 'sp',
                'P6HR' : 'tp',
                'TCW' : 'tcw',
                'elev' : 'orog',
                'lat' : 'latitude',
                'lon' : 'longitude',
                'time' : 'time',
                'mem' : 'number'
            
                }
        
        outdir = '/home/disk/hot/stangen/Documents/prior_ensembles/'+ens+'/'+y+m+'/'
    
        #Create output directories if they don't yet exit
        if (os.path.isdir(outdir)):
            pass
        else:
            os.makedirs(outdir)
        
        surf_file = '/home/disk/hot/stangen/Documents/tigge_ensembles/'+ens+'/'+y+m+'/'+y+'-'+m+'-'+d+'_'+h+'_'+ens+'_'+surf_str
        upper_file = '/home/disk/hot/stangen/Documents/tigge_ensembles/'+ens+'/'+y+m+'/'+y+'-'+m+'-'+d+'_'+h+'_'+ens+'_'+upper_str
        
        #Load data about the surface ensemble
        surf = Dataset(surf_file, 'r')     
        # Shape of surf[variable] is nvars, ntimes, nmems, nlats, nlons
        #print(ncdata.variables.keys())
        # Find the indices corresponding to the start and end times
        tunit = surf.variables[vardict['time']].units
        ftimes = num2date(surf.variables[vardict['time']][:],tunit)
        nmems = len(surf.dimensions[vardict['mem']])
        #nmems4name = len(memfiles)
        ntimes = len(ftimes)
        nlats = len(surf.dimensions[vardict['lat']])
        nlons = len(surf.dimensions[vardict['lon']])
        
        # For the metadata, need a list of locations
        lats = surf.variables[vardict['lat']][:][:,None]
        lons = surf.variables[vardict['lon']][:][None,:]
        #And an array of ensemble members
        memarr = np.arange(1,nmems+1)
        
        p_surf = surf.variables['sp'][:]
        u_surf = surf.variables['u10'][:]
        v_surf = surf.variables['v10'][:]
        wind_surf = np.sqrt(u_surf**2+v_surf**2)

        
        #should have same lat/lon, members, and times as the surface
        #ensemble, so no need to load them here.
        #Shape of aloft[variable] is ntimes x nmems x nlevs x nlats x nlons
        #levs is ascending the atmosphere as index increases
        aloft = Dataset(upper_file, 'r')
        
        #time range of ensemble
        ftime_diff = ftimes[-1]-ftimes[0]
        tr = int((ftime_diff.days)*24 + (ftime_diff.seconds)/3600)
        tr_str = str(tr)+'hrs'
        
        # Allocate the state array
        print('Allocating the state vector array...')
        #state = np.zeros((nvars,ntimes,nlats,nlons,nmems))
        state = np.zeros((nvars,ntimes,nmems,nlats,nlons))
        
        #get surface vapor pressure (hPa, or mb)
        e = 6.112 * np.exp((17.67*(surf.variables['d2m'][:]-273.15))/((surf.variables['d2m'][:]-273.15)+243.5))
        #specific humidity (q) (vapor pressure and surface pressures in hPa)
        q_surf = 0.622 * e /(p_surf/100-0.378*e)
        
        #loop through each pressure level aloft
        for i, lev in enumerate(levels):
            #convert pressure to integer in pascals
            p_aloft = int(lev)*100
            q_aloft = aloft.variables['q'][:,:,i,:,:]
            u_aloft = aloft.variables['u'][:,:,i,:,:]
            v_aloft = aloft.variables['v'][:,:,i,:,:]
            
   #--------first time through the loop: we need to get the stuff less than 1000 hPa
            if i == 0:
                print('Working on surface pressure < 1000 hPa')
                #make array of boolean where surface pressure is greater than 1000 hPa.
                #this will create an array of the same shape of the ensemble, with
                #True where the condition is met, and False otherwise. 
                p_surf_less_1000 = p_surf >= p_aloft
 
                #calculate IWV of this column
                #where the condition is not true, it will create values of 0
                dp = np.subtract(p_surf,p_aloft,where=p_surf_less_1000)
                q_mean = np.add(q_surf,q_aloft,where=p_surf_less_1000)/2
                IWV = np.multiply(q_mean,dp,where=p_surf_less_1000)/g
                
                state[0,:,:,:,:] = np.add(state[0,:,:,:,:],IWV)
                #calculate IVT of this column                
                u_wind_mean = np.add(u_surf,u_aloft,where=p_surf_less_1000)/2
                v_wind_mean = np.add(v_surf,v_aloft,where=p_surf_less_1000)/2                
                u_IVT = np.multiply(IWV,u_wind_mean,where=p_surf_less_1000)
                v_IVT = np.multiply(IWV,v_wind_mean,where=p_surf_less_1000)
               
            else:
                print('Working on',lev)
                p_aloft_lower = int(levels[i-1])*100
                #if the surface pressure is less than both layers, no need to add
                #for that gridbox, since the pressure levels are below ground there.
                
    #-----------if surface pressure is between the 2 pressure levels
                p_surf_between = (p_surf >= p_aloft) & (p_surf < p_aloft_lower)
                #average surface values with values of upper layer, and get 0s
                #where this isn't the case
                dp_between = np.subtract(p_surf,p_aloft,where=p_surf_between)
                q_mean = np.add(q_surf,q_aloft,where=p_surf_between)/2
                IWV_between = np.multiply(q_mean,dp_between,where=p_surf_between)/g               
                state[0,:,:,:,:] = np.add(state[0,:,:,:,:],IWV_between)
                
                u_wind_mean_between = np.add(u_surf,u_aloft,where=p_surf_between)/2
                v_wind_mean_between = np.add(v_surf,v_aloft,where=p_surf_between)/2
                u_IVT_between = np.multiply(IWV_between,u_wind_mean_between,where=p_surf_between)
                v_IVT_between = np.multiply(IWV_between,v_wind_mean_between,where=p_surf_between)
                
    #-----------if surface pressure is greater than both pressure levels-
                #the two layers are both above ground               
                p_surf_greater = (p_surf > p_aloft) & (p_surf > p_aloft_lower)
                
                #dp is a scalar in this case
                dp_aloft = p_aloft_lower-p_aloft
                q_aloft_lower = aloft.variables['q'][:,:,i-1,:,:]
                q_mean = np.add(q_aloft_lower,q_aloft,where=p_surf_greater)/2
                IWV_greater = np.multiply(q_mean,dp_aloft,where=p_surf_greater)/g
                state[0,:,:,:,:] = np.add(state[0,:,:,:,:],IWV_greater)
                
                u_aloft_lower = aloft.variables['u'][:,:,i-1,:,:]
                v_aloft_lower = aloft.variables['v'][:,:,i-1,:,:]
                
                u_wind_mean_greater = np.add(u_aloft_lower,u_aloft,where=p_surf_greater)/2
                v_wind_mean_greater = np.add(v_aloft_lower,v_aloft,where=p_surf_greater)/2
                u_IVT_greater = np.multiply(IWV_greater,u_wind_mean_greater,where=p_surf_greater)
                v_IVT_greater = np.multiply(IWV_greater,v_wind_mean_greater,where=p_surf_greater)     
                
                #add the IVT from each category of surface pressure to layer pressure
                u_IVT = u_IVT + u_IVT_between + u_IVT_greater
                v_IVT = v_IVT + v_IVT_between + v_IVT_between
                                
        #find magnitude of IVT 
        state[1,:,:,:,:] = np.sqrt(u_IVT**2+v_IVT**2)
        #find direction of IVT
        state[2,:,:,:,:] = np.arctan2(v_IVT,u_IVT)*180/np.pi
        
        #get state so that it is nvars x ntimes x nlats x nlons x nmems
        state = np.swapaxes(state, 2, 3)
        state = np.swapaxes(state, 3, 4)
        
        print('Writing to netcdf...')
        # Convert times back to integers
        valid_times = date2num(ftimes,tunit)
        
        # Write ensemble forecast to netcdf - change name here
        #dset = Dataset(outdir+y+'-'+m+'-'+d+'_'+h+'_'+ens_type+'_'+outvar_string+'.nc','w')
        dset = Dataset(outdir+y+'-'+m+'-'+d+'_'+h+'_'+tr_str+'_'+outvar_string+'.87nc','w')
        dset.createDimension('time',None)
        dset.createDimension('lat',nlats)
        dset.createDimension('lon',nlons)
        dset.createDimension('ens',nmems)
        dset.createVariable('time','i4',('time',))
        dset.createVariable('lat',np.float32,('lat',))
        dset.createVariable('lon',np.float32,('lon'))
        dset.createVariable('ens','i4',('ens',))
        dset.variables['time'].units = tunit
        dset.variables['lat'].units = 'degrees_north'
        dset.variables['lon'].units = 'degrees_east'
        dset.variables['ens'].units = 'member_number'
        dset.variables['time'][:] = np.array(valid_times)
        dset.variables['lat'][:] = lats
        dset.variables['lon'][:] = lons
        dset.variables['ens'][:] = memarr
        for v,var in enumerate(end_variables):
            #var = vardict[var]
            print('Writing variable {}'.format(var))
            dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
            dset.variables[var].units = ef.get_units(var)
            dset.variables[var][:] = state[v,:,:,:,:]
        #completes writing the file
        dset.close()

        