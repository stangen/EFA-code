#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 12:22:51 2018

@author: stangen
"""
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime
import numpy as np
import EFA.efa_files.cfs_utilities_st as ut
import surface_obs.madis_example.madis_utilities as mt
import EFA.duplicate_madaus.efa_functions as ef
import os

#--------------Change these----------------------------------------------------
#start and end date to get ensembles. 
start_date = datetime(2015,11,10,0) #YYYY,m,d,h
end_date = datetime(2015,11,15,12)
hourstep = 12 #how often you want a new forecast initialization, usually 12 hr
# ecmwf, eccc for euro/canadian ensembles, ncep
ensemble_type = ['ncep']
#variables with names coming from the raw TIGGE- see get_tigge_data, or tigge 
#output filename if unsure of names.
#the order matters to make filename match exactly. 
in_variables = ['TCW']#['T2M','SP']#
#sfc for surface, pl for aloft
lev = 'sfc'
#variables with names I want to have after it is processed
variables = ['TCW']#['ALT','T2M']##, 'P6HR', 'TCW']
#a list of dates to loop through to load each forecast initialized on these dates
dates = mt.make_datetimelist(start_date,end_date,hourstep)   
#------------------------------------------------------------------------------


def create_new_netcdf(date,ens_type,in_vrbls,vrbls):
    """
    This function creates a netCDF from a raw TIGGE netCDF. The main purpose
    of this is to change surface pressure to altimeter setting, calculate
    6-hourly precipitation, and to rename/shorten variable names. Unfortunately
    using the float32 format for the variables uses twice the memory of the 
    raw TIGGE int16 format. 
    """

    y = date.strftime('%Y')
    m = date.strftime('%m')
    d = date.strftime('%d')
    h = date.strftime('%H')
    
    print('Working on '+d+'_'+h+' '+ens_type)
    
    #rename TIGGE variable names
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
    
    #build a var string corresponding with naming convention of tigge files
    invar_string = ef.var_string(in_vrbls)
    
    #build a var string corresponding with naming convention of output files
    outvar_string = ef.var_string(vrbls)
    
    
    # This is the input directory for the raw TIGGE netcdf
    indir  = '/home/disk/hot/stangen/Documents/tigge_ensembles/'+ens_type+'/'+y+m+'/'+y+'-'+m+'-'+d+'_'+h+'_'+ens_type+'_'+invar_string+'_'+lev+'.nc'
    # This is the directory for the orography file
    orography = '/home/disk/hot/stangen/Documents/tigge_ensembles/orography/2013-04-01_00_'+ens_type+'.nc'
    # This is the output directory for the netcdf with altimeter setting 
    outdir = '/home/disk/hot/stangen/Documents/prior_ensembles/'+ens_type+'/'+y+m+'/'
    
    #Create output directories if they don't yet exit
    if (os.path.isdir(outdir)):
        pass
    else:
        os.makedirs(outdir)
    
    #Read the ensemble netcdf file
    ncdata = Dataset(indir,'r')
    
    # Shape of ncdata is nvars, ntimes, nmems, nlats, nlons
    #print(ncdata.variables.keys())
    tunit = ncdata.variables[vardict['time']].units
    ftimes = num2date(ncdata.variables[vardict['time']][:],tunit)
    nmems = len(ncdata.dimensions[vardict['mem']])
    ntimes = len(ftimes)
    nvars = len(vrbls)
    nlats = len(ncdata.dimensions[vardict['lat']])
    nlons = len(ncdata.dimensions[vardict['lon']])
    
    #time range of ensemble, for naming of file
    ftime_diff = ftimes[-1]-ftimes[0]
    tr = int((ftime_diff.days)*24 + (ftime_diff.seconds)/3600)
    tr_str = str(tr)+'hrs'
    
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
       
    # Now to populate the state array
    for va, var in enumerate(vrbls):
        #reason index 0:2 is checked is that QF and D-QF can take place at different
        #levels (QF850, D-QF850), so to avoid specifying each possible variable name, 
        #just check 1st 2 chars.
        if var[0:2] not in ['QF','D-']:
            field = ncdata.variables[vardict[var]][:,:,:,:]
            #print(field)
            print('Adding variable {}'.format(var))
            #convert surface pressure to altimeter setting
            if var=='ALT': 
                #Read the orography netcdf file- for calculating altimeter setting
                #has a time index, even when only one time is gotten from TIGGE- requires
                #indexing like [0,:,:] to get this first (and only) time.
                orogdata = Dataset(orography, 'r')
                elev = orogdata.variables[vardict['elev']][0,:,:]  
                #Convert surface pressure to altimeter setting in mb
                #pressure in netcdf file is in pascals
                presinmb = field/100
                field = presinmb/((288-0.0065*elev)/288)**5.2561
            #find 6-hourly precipitation
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
        #magnitude of moisture flux
        if var[0:2] =='QF':
            q = ncdata.variables['q'][:,:,:,:]
            u = ncdata.variables['u'][:,:,:,:]
            v = ncdata.variables['v'][:,:,:,:]
            #moisture flux is qV, to find magnitude find distance from origin
            #to point, multiply by q, multiply by 1000 to get in units of g/kg.
            field = q*np.sqrt(u**2+v**2)*1000
        #direction of moisture flux (-180 to 180, unit circle degrees)    
        if var[0:2] =='D-':
            u = ncdata.variables['u'][:,:,:,:]
            v = ncdata.variables['v'][:,:,:,:]
            field = np.arctan2(v,u)*180/np.pi
                #print(field.shape)
        # make the ensemble dimension at the end of state
        field = np.swapaxes(field, 1, 3)
        field = np.swapaxes(field, 1, 2)
        # Populate its component of the state array
        state[va,:,:,:,:] = field
        
    print('Writing to netcdf...')
    # Convert times back to integers
    valid_times = date2num(ftimes,tunit)
    
    # Write ensemble forecast to netcdf - change name here
    #dset = Dataset(outdir+y+'-'+m+'-'+d+'_'+h+'_'+ens_type+'_'+outvar_string+'.nc','w')
    dset = Dataset(outdir+y+'-'+m+'-'+d+'_'+h+'_'+tr_str+'_'+outvar_string+'44.nc','w')
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
    for v,var in enumerate(vrbls):
        #var = vardict[var]
        print('Writing variable {}'.format(var))
        dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
        dset.variables[var].units = ef.get_units(var)
        dset.variables[var][:] = state[v,:,:,:,:]
    #completes writing the file
    dset.close()


#--------Call the function to make a new netCDF--------------------------------
for ens in ensemble_type:
    for d in dates:
        create_new_netcdf(d,ens,in_variables,variables)
        
        
