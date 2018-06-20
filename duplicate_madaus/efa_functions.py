#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 16:23:50 2018

@author: stangen
"""

import numpy as np

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
    print(a.min())
    dist = 2*R*np.arcsin(np.sqrt(a))
    #flatten the distance array to 1D
    dist_flat = dist.flatten()
    

    
    
    #get the indices of the k closest gridpoints
    idx = dist_flat.argsort(axis=None)[:k]
    print(dist_flat[idx])

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