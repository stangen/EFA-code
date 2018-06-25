#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 16:23:50 2018

@author: stangen
"""

import numpy as np
from netCDF4 import Dataset, num2date, date2num
import EFA.efa_files.cfs_utilities_st as ut

def closest_points(ob_lat, ob_lon, lats, lons, ob_elev, elevs, ob_time=None, 
                   times=None, variable=None, k=4, need_interp=False):
    """
    Function which uses the Haversine formula to compute the distances from one 
    observation lat/lon to each gridpoint. It finds the indices of the n closest 
    gridpoints (default 4). It then calls functions to check elevation and/or
    interpolate the ensemble to the observation location.
    
    ob_lat: latitude of the observation
    ob_lon: longitude of the observation
    lats: 1-d (masked) array of latitudes from ens netCDF file
    lons: 1-d (masked) array of longitudes from ens netCDF file
    ob_elev: required for check elevation function
    elevs: required for check elevation function
    ob_time: required for interpolate function
    times: required for interpolate function
    variable: required for interpolate function   
    k = number of points to find
    need_interp: if False, call function to check the elevation of the 4 
    nearest gridpoints. If True, call function to interpolate the ensemble 
    to the ob location/time.
    returns: output from check elevation or interpolate functions.
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
    dist = 2*R*np.arcsin(np.sqrt(a))
    #flatten the distance array to 1D
    dist_flat = dist.flatten()

    #find nearest 4 points (indices), in order from closest to farthest
    closest_4 = dist_flat.argsort(axis=None)[:k]
    #print(dist_flat[closest_4])
    #do we need to do interpolation?
    if need_interp==True:
        #call function to check elevation
        TorF = check_elev(closest_4,elevs,ob_elev)
        #if it passes terrain check, we need to interpolate
        if TorF==True:
            #4 closest distances
            distances = dist_flat[closest_4]
            #return to lat/lon gridbox indices
            closey, closex = np.unravel_index(closest_4, lonarr.shape)
            interp = interpolate(closey,closex,distances,times,variable,ob_time)
            return interp
        #if it fails terrain check, no need to interpolate, so return empty array
        elif TorF==False:
            return np.array(())
    #if no need to do interpolation, just call/return terrain check function
    elif need_interp==False:
        return check_elev(closest_4,elevs,ob_elev)


def interpolate(closey, closex, distances, times, variable, ob_time):
    """
    Function which interpolates the state to the observation location for a given variable
    
    closey = 4 closest lat indices of the ensemble
    closex = 4 closest lon indices of the ensemble
    distances = distance from 4 closest gridpoints to the observation, for use in weights
    times = forecast times of the ensemble, in datetime format
    variable = ensemble values of variable (ALT or T2M)
    ob_time = observation time, in datetime format
    returns: value of ensemble interpolated to ob location.

    """
    
    spaceweights = np.zeros(distances.shape)
    if (distances < 1.0).sum() > 0:
        spaceweights[distances.argmin()] = 1
    else:
        # Here, inverse distance weighting (for simplicity)
        spaceweights = 1.0 / distances
        spaceweights /= spaceweights.sum()
    # Get weights in time
    #all the valid times in the ensemble
    valids = times
    timeweights = np.zeros(valids.shape)
    # Check if we are outside the valid time range
    if (ob_time < valids[0]) or (ob_time > valids[-1]):
        print("Interpolation is outside of time range in state!")
        return None
    # Find where we are in this list
    #index after the time of the observation
    lastdex = (valids >= ob_time).argmax()
    # If we match a particular time value, then
    # this is just an identity
    if valids[lastdex] == ob_time:
        # Just make a one at this time
        timeweights[lastdex] = 1
    else:
        # Linear interpolation
        diff  = (valids[lastdex] - valids[lastdex-1])
        #often going to be 21600 seconds
        totsec = diff.seconds
        #calculate time difference between time after and time of observation
        #the abs will make this positive definite, which is okay since
        #the difference will always be negative
        thisdiff = abs(ob_time - valids[lastdex])
        thissec = thisdiff.seconds
        # Put in appropriate weights
        timeweights[lastdex-1] = float(thissec) / totsec
        timeweights[lastdex] = 1.0 - (float(thissec)/totsec)
    # Now that we have the weights, do the interpolation
    #an ntimes x 4 x nens array
    interp = variable[:,closey,closex,:]
    # Do a dot product with the time weights
    # And with the space weights
    if len(interp.shape) == 3:
        interp = (timeweights[:,None,None] * interp).sum(axis=0)
    else:
        interp = (timeweights[:,None,None,None] * interp).sum(axis=0)

    if len(interp.shape) == 3:
        interp = (spaceweights[:,None,None] * interp).sum(axis=1)
    else:
        interp = (spaceweights[:,None] * interp).sum(axis=0)
    # Return estimate from all ensemble members
    return interp

def check_elev(idx,elevs,ob_elev):
    """
    Function that takes in several variables to perform a terrain check 
    on observations- is the ens elevation of the 4 nearest gridpoints within 
    300 m of the ob elevation? 
    idx: index of 4 closest gridpoints of flattened lat/lon array
    elevs: nlat x nlon grid of elevations of the ens netCDF file
    ob_elev: elevation of the observation    
    Returns: True or False- True if passes terrain check, false if it doesn't.
    """
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