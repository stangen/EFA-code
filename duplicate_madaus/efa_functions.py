#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 16:23:50 2018

@author: stangen
"""

import numpy as np
from netCDF4 import Dataset, num2date, date2num
import EFA.efa_files.cfs_utilities_st as ut

def check_elev(lats,lons,elevs,ob_lat,ob_lon,ob_elev,k=4):
    """
    Function that takes in several variables to perform a terrain check 
    on observations. 
    lats: 1-d (masked) array of latitudes from ens netCDF file
    lons: 1-d (masked) array of longitudes from ens netCDF file
    elevs: nlat x nlon grid of elevations of the ens netCDF file
    ob_lat: latitude of the observation
    ob_lon: longitude of the observation
    ob_elev: elevation of the observation
    k: number of gridpoints to compare with, default is 4
    
    Returns: True or False
    """
    
    #make a lat/lon array for comparison with station lat/lon
    lonarr, latarr = np.radians(np.meshgrid(lons, lats))
    
    R = 6371. # Radius of earth in km
    #convert degrees to radians
    ob_lat = np.radians(ob_lat)
    ob_lon = np.radians(ob_lon)
    #radian distance between points
    dlat = ob_lat-latarr
    dlon = ob_lon-lonarr
    #Haversine formula to calculate the great circle distance between ob and each lat/lon point (in km)
    a = np.sin(dlat/2)**2 + np.cos(ob_lat) * np.cos(latarr) * np.sin(dlon/2)**2
    #print(a.min()) #was to check if negative values were really occurring, don't think so.
    dist = 2*R*np.arcsin(np.sqrt(a))
    #flatten the distance array to 1D
    dist_flat = dist.flatten()
    

    
    
    #get the indices of the k closest gridpoints
    idx = dist_flat.argsort(axis=None)[:k]
    #print(dist_flat[idx]) #prints distances of 4 closest gridpoints

    #find nearest 4 points, in order from closest to farthest
    #argpartition just puts the 4 smallest indices in front, array (1,k) puts the 
    #first four in order, the rest are not necessarily sorted

#   idx = np.argpartition(dist_flat, (1,k))
#   print(idx[0:4])
    
    #check if the 4 closest points have elevations which differ by more than 300 meters
    elev_flat = elevs.flatten()
    elev_four_closest = elev_flat[idx]
    elev_diff = abs(ob_elev-elev_four_closest)
    #print(elev_diff)
    if elev_diff.max() > 300:
        TorF = False
    else:
        TorF = True
    
    return TorF

def var_string(vrbls):
    """
    Function which takes in a list of variables and returns a string, formatted
    like var1_var2_var3, for use in saving files. 
    """
    var_string = ''
    for v in vrbls[:-1]:  
        var_string = var_string+v+'_'
    var_string = var_string+vrbls[-1]
    return var_string

def make_netcdf(state,outfile):
    """
    Function which creates/saves a netCDF file
    """
    #tunit='seconds since 1970-01-01'
    tunit='hours since 1900-01-01'
    # Write ensemble forecast to netcdf
    with Dataset(outfile,'w') as dset:
            dset.createDimension('time',None)
            dset.createDimension('lat',state.ny())
            dset.createDimension('lon',state.nx())
            dset.createDimension('ens',state.nmems())
            dset.createVariable('time','i4',('time',))
            dset.createVariable('lat',np.float32,('lat',))
            dset.createVariable('lon',np.float32,('lon'))
            dset.createVariable('ens','i4',('ens',))
            dset.variables['time'].units = tunit
            dset.variables['lat'].units = 'degrees_north'
            dset.variables['lon'].units = 'degrees_east'
            dset.variables['ens'].units = 'member_number'
            dset.variables['time'][:] = date2num(state.ensemble_times(),tunit)
            dset.variables['lat'][:] = state['lat'].values[:,0]
            dset.variables['lon'][:] = state['lon'].values[0,:]
            dset.variables['ens'][:] = state['mem'].values
            for var in state.vars():
                print('Writing variable {}'.format(var))
                dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
                dset.variables[var].units = ut.get_units(var)
                dset.variables[var][:] = state[var].values