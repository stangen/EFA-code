#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 13:53:47 2018

@author: stangen
"""
import numpy as np

def get_ob_points(left,right,top,bottom,spacing):
    """
    Function which selects gridpoints for observations, based on the lat/lon
    boundaries and the spacing of the gridpoints (in deg lat/lon) at the lowest
    longitude. Causes less observation gridpoints to be selected at higher 
    latitudes than lower latitudes, since there is less distance between
    gridpoints at higher latitudes.
    """
    
#    left = np.radians(left)
#    right = np.radians(right)
#    top = np.radians(top)
#    bottom = np.radians(bottom)
    
    left = -180
    right = 180
    bottom = 0
    top = 90
    spacing = 2
    
    #create list of latitudes from given box
    lats = range(bottom,top+spacing,spacing)
    #find number of longitudes, based on latitude and size of box
    nlons = (right-left)/spacing*np.cos(np.radians(lats))
    lons = []
    #find longitudes equally spaced within the box
    for n in nlons:
        lons.append(np.linspace(left,right,num=int(round(n)),endpoint=False))
        
    #generate lat/lon pairs
    latlon = []
    for i,lat in enumerate(lats):
        for lon in lons[i]:
            latlon.append([lat,lon])
            
    return latlon
        
    
    
    