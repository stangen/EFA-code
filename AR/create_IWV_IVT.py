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
start_date = datetime(2015,11,10,0) #YYYY,m,d,h
end_date = datetime(2015,11,10,0)
hourstep = 12 #how often you want a new forecast initialization, usually 12 hr

#variables with names coming from the raw TIGGE- see get_tigge_data if unsure of names.
#the order matters to make filename match exactly. 
surf_variables = ['D2M','SP','U10','V10']
upper_variables = ['Q','U','V']
levels = ['1000','925','850','700','500','300']
end_variables = ['IVT','IWV']

#sfc for surface, pl for elevated
lev = 'sfc'

#variables with names I want to have after it is processed
variables = ['TCW']#['ALT','T2M']##, 'P6HR', 'TCW']

#a list of dates to loop through to load each forecast initialized on these dates
dates = mt.make_datetimelist(start_date,end_date,hourstep)  

surf_str = ef.var_string(surf_variables)+'_sfc.nc'
upper_str = ef.var_string(upper_variables)+'_'+ef.var_string(levels)+'_pl.nc'
nvars = len(end_variables)

"""
This function creates a netCDF from a raw TIGGE netCDF. The main purpose
of this is to change surface pressure to altimeter setting, calculate
6-hourly precipitation, and to rename/shorten variable names. Unfortunately
using the float32 format for the variables uses twice the memory of the 
raw TIGGE int16 format. 
"""
for ens in ensemble_type:
    for date in dates:

        y = date.strftime('%Y')
        m = date.strftime('%m')
        d = date.strftime('%d')
        h = date.strftime('%H')
        
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
        # Shape of surf is nvars, ntimes, nmems, nlats, nlons
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
        
        #should have same lat/lon, members, and times as the surface
        #ensemble, so no need to load them here.
        #Shape of aloft is ntimes x nmems x nlevs x nlats x nlons
        #levs is ascending the atmosphere as index increases
        aloft = Dataset(upper_file, 'r')
        
        #time range of ensemble
        ftime_diff = ftimes[-1]-ftimes[0]
        tr = int((ftime_diff.days)*24 + (ftime_diff.seconds)/3600)
        tr_str = str(tr)+'hrs'
        
        # Allocate the state array
        print('Allocating the state vector array...')
        state = np.zeros((nvars,ntimes,nlats,nlons,nmems))
        
        #get surface vapor pressure (hPa, or mb)
        e = 6.112 * np.exp((17.67*(surf.variables['d2m'][:]-273.15))/((surf.variables['d2m'][:]-273.15)+243.5))
        #specific humidity (q) 
        q = 0.622 * e /(surf.variables['sp'][:]-0.378*e)
        