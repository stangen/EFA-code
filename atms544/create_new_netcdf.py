#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 12:22:51 2018

@author: stangen
"""
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime, timedelta
from subprocess import check_output
import numpy as np
import efa_files.cfs_utilities_st as ut
from efa_xray.state.ensemble import EnsembleState
import os

def create_new_netcdf(y,m,d,h,ens_type, vrbls):

    m_dict = {
               '01' : 'jan',
               '02' : 'feb',
               '03' : 'mar',
               '04' : 'apr',
               '05' : 'may',
               '06' : 'jun',
               '07' : 'jul',
               '08' : 'aug',
               '09' : 'sep',
               '10' : 'oct',
               '11' : 'nov',
               '12' : 'dec'           
               }
    
    vardict = {
                'T2M' : 't2m',
                'ALT' : 'sp',
                'P6HR' : 'tp',
                'TCW' : 'tcw',
                'elev' : 'orog',
                'lat' : 'latitude',
                'lon' : 'longitude',
                'time' : 'time',
                'mem' : 'number'
            
            }
    
    # This is the input directory for the raw TIGGE netcdf
    indir  = '/home/disk/hot/stangen/Documents/ensembles/'+ens_type+'/'+m_dict[m]+y+'/'+y+'-'+m+'-'+d+'_'+h+'_'+ens_type+'_T_SP_Precip_TCW_240hr.nc'
    # This is the directory for the orography file
    orography = '/home/disk/hot/stangen/Documents/ensembles/orography/2013-04-01_00_'+ens_type+'.nc'
    # This is the output directory for the netcdf with altimeter setting 
    outdir = '/home/disk/hot/stangen/Documents/atms544/ensembles/'+ens_type+'/'+m_dict[m]+y+'/'
    
    #Create output directories if they don't yet exit
    if (os.path.isdir(outdir)):
        pass
    else:
        os.makedirs(outdir)
    
    #Read the ensemble netcdf file
    ncdata = Dataset(indir,'r')
    
    # Shape of ncdata is nvars, ntimes, nmems, nlats, nlons
    #print(ncdata.variables.keys())
    # Find the indices corresponding to the start and end times
    tunit = ncdata.variables[vardict['time']].units
    ftimes = num2date(ncdata.variables[vardict['time']][:],tunit)
    nmems = len(ncdata.dimensions[vardict['mem']])
    #nmems4name = len(memfiles)
    ntimes = len(ftimes)
    nvars = len(vrbls)
    nlats = len(ncdata.dimensions[vardict['lat']])
    nlons = len(ncdata.dimensions[vardict['lon']])
    
    # Allocate the state array
    print('Allocating the state vector array...')
    state = np.zeros((nvars,ntimes,nlats,nlons,nmems))
    
    # For the metadata, need a list of locations
    lats = ncdata.variables[vardict['lat']][:][:,None]
    lons = ncdata.variables[vardict['lon']][:][None,:]
    # Do a 2d mesh of lat and lon
    lonarr, latarr = np.meshgrid(lons, lats)
    
    #And an array of ensemble members
    memarr = np.arange(1,nmems+1)
    
    #Read the orography netcdf file- for calculating altimeter setting
    orogdata = Dataset(orography, 'r')
    #print(orogdata.variables)
    elev = orogdata.variables[vardict['elev']][0,:,:]    
    
    # Now to populate the state array
    for v, var in enumerate(vrbls):
        field = ncdata.variables[vardict[var]][:,:,:,:]#[tbeg:tend,:,:]
        #print(field)
        print('Adding variable {}'.format(var))
        #convert surface pressure to altimeter setting
        if var=='ALT': 
            #Convert surface pressure to altimeter setting in mb
            #pressure in netcdf file is in pascals
            presinHg = field*.00029528744
            field = (presinHg/((288-0.0065*elev)/288)**5.2561)/.029528744
        if var=='P6HR':
            #create dummy field to facilitate subtracting of total precipitation
            #t-1 from t without saving over t, so the next subtraction still works
            field2 = np.zeros((ntimes,nmems,nlats,nlons))
            for t in range(0,ntimes):
                if t == 0:
                    field2[t,:,:,:] = field[t,:,:,:]
                #subtract previous time's precipitation to get 6 hour precipitation
                elif t > 0:
                    field2[t,:,:,:] = field[t,:,:,:] - field[t-1,:,:,:]
            #reassign 6 hour precip to field
            field = field2
        # make the ensemble dimension at the end of state
        field = np.swapaxes(field, 1, 3)
        field = np.swapaxes(field, 1, 2)
         # Populate its component of the state array
        state[v,:,:,:,:] = field
        
    print('Writing to netcdf...')
        # Convert times back to integers
    valid_times = date2num(ftimes,tunit)
    
    # Write ensemble forecast to netcdf - change name here
    dset = Dataset(outdir+y+'-'+m+'-'+d+'_'+h+'_'+ens_type+'4var.nc','w')
    dset.createDimension('time',None)
    dset.createDimension('lat',nlats)
    dset.createDimension('lon',nlons)
    dset.createDimension('ens',nmems)
    dset.createVariable('time','i4',('time',))
    dset.createVariable('lat',np.float64,('lat',))
    dset.createVariable('lon',np.float64,('lon'))
    dset.createVariable('ens','i4',('ens',))
    dset.variables['time'].units = tunit
    dset.variables['lat'].units = 'degrees_north'
    dset.variables['lon'].units = 'degrees_east'
    dset.variables['ens'].units = 'member_number'
    dset.variables['time'][:] = np.array(valid_times)
    dset.variables['lat'][:] = lats
    dset.variables['lon'][:] = lons
    dset.variables['ens'][:] = memarr
    for v,var in enumerate(vrbls):
        #var = vardict[var]
        print('Writing variable {}'.format(var))
        dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
        dset.variables[var].units = ut.get_units(var)
        dset.variables[var][:] = state[v,:,:,:,:]
    # Free up memory held by the state array
    del state



# form 'YYYY'
year = '2013'
# form 'mmm'
month = '04'
# form 'DD'
d = range(1,2)
day_string = []
for day in d:
    day = format(day, '02d')
    day = str(day)
    day_string.append('%s' % (day))
# form 'HH'
hour = ['00']#,'12']
# ecmwf, eccc for euro/canadian ensembles, ncep
ensemble_type = ['ecmwf']
#variables with names I want to have after it is processed
variables = ['T2M', 'ALT', 'P6HR', 'TCW']

for ens in ensemble_type:
    for d in day_string:
        for h in hour:
            print('Working on '+d+'_'+h+' '+ens)
            create_new_netcdf(year,month,d,h,ens,variables)
        
        
