#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 16:00:52 2018

@author: stangen
"""

from __future__ import print_function
import numpy as np
import pandas as pd
import xarray
import xarray.ufuncs as xu
import netCDF4
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, time
import pytz
from copy import deepcopy
import madis_utilities as mt

def nearest_points(lat, lon, npt=1):
    """
    Use the lat-lon arrays to return a list of indices
    of the nearest npt points to the given lat-lon
    """
    # Use sin of lat lon to handle periodic
    # and not worry about if we are in negative
    # degrees
    dist = xu.hypot(xu.sin(xu.radians(lat)) -
             xu.sin(xu.radians(lat)),\
             xu.cos(xu.radians(lon)) - 
             xu.cos(xu.radians(lon)))
    # Get indices of the flattened array
    nearest_raw = dist.argsort(axis=None)[:npt]
    # Convert back to 2-d coords
    nearest = np.unravel_index(nearest_raw, ['lat'].shape)
    return nearest

def distance_to_point(self, lat, lon):
        """
        Use Haversine formula to estimate distances from all
        gridpoints to a given location (lat, lon)
        """
        R = 6371. # Radius of earth in km
        lat = np.radians(lat)
        lon = np.radians(lon)
        dlat = lat - xu.radians(self['lat'].values)
        dlon = lon - xu.radians(self['lon'].values)
        a = xu.sin(dlat/2)**2 + xu.cos(lat) * xu.cos(xu.radians(self['lat'].values)) * \
                xu.sin(dlon/2)**2
        c = 2 * xu.arctan2(xu.sqrt(a), xu.sqrt(1.0-a))
        return R*c

base_dir = '/home/disk/hot/stangen/Documents/ensembles/orography/2013-04-01_00_eccc.nc'

of = netCDF4.Dataset(base_dir,"r")
#print(of.variables)
lats = of.variables['latitude'][:] # latitude of MADIS observation
lngs = of.variables['longitude'][:] # longitude of MADIS observation
elevs = of.variables['orog'][0,:] # Elevation of MADIS observation 
epoch = of.variables['time'][:] # Time of MADIS Observation
#make a lat/lon array for comparison with station lat/lon
lonarr, latarr = np.radians(np.meshgrid(lngs, lats))

fname = '/home/disk/hot/stangen/Documents/EFA/surface_obs/MADIS/201304/combined_alts/alts_20130430_0000.txt'

f1 = open(fname, 'r')
lines = f1.readlines()

s1 = lines[1000]

s1_split = s1.split(',')
#get the lat/lon of the station
lat = float(s1_split[1])
lon = float(s1_split[2])

R = 6371. # Radius of earth in km
#convert degrees to radians
lat = np.radians(lat)
lon = np.radians(lon)
#radian distance between points
dlat = lat-latarr
dlon = lon-lonarr
#Haversine formula to calculate the great circle distance between ob and each lat/lon point
a = np.sin(dlat/2)**2 + np.cos(lat) * np.cos(latarr) * np.sin(dlon/2)**2
dist = 2*R*np.arcsin(np.sqrt(a))
#flatten the distance array to 1D
dist_flat = dist.flatten()

#find nearest 4 points, in order from closest to farthest
#argpartition just puts the 4 smallest indices in front, array (1,k) puts the 
#first four in order, the rest are not necessarily sorted
k = 4
idx = np.argpartition(dist_flat, (1,k))

#check if the 4 closest points have elevations which differ by more than 300 meters
elev_flat = elevs.flatten()
elev_four_closest = elev_flat[idx[:k]]
elev_ob = float(s1_split[3])
elev_diff = abs(elev_ob-elev_four_closest)
if elev_diff.max() > 300:
    print('bad')
else:
    print('good')
#for elev in elev_four_closest:
#    if abs(elev_ob-elev > 300):
#        print('bad')
#    else:
#        print('good')

#a = xu.sin(dlat/2)**2 + xu.cos(lat) * xu.cos(latarr) * xu.sin(dlon/2)**2
#c = 2 * xu.arctan2(xu.sqrt(a), xu.sqrt(1.0-a))
#dist = R*c
#shortxu = dist.min()

short_dist_value = dist.min()
shortindex = np.argmin(dist_flat)
shortindexvalue = shortindex.min()

#latlonidx = b.index(min(b))

print(short_dist_value)
print(dist_flat[shortindex])


four_closest_distances = dist_flat[idx[:k]]

#print(dist)
#return R*c


#distance from 

#thing = nearest_points(lat,lon)
#for line in f1:
#    print(line)

#Create a function which takes in ecmwf or eccc, observation, returns true or false? 
#to know whether or not to assimilate the ob, if it passes terrain check