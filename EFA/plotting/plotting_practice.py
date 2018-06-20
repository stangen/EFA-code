#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 15:38:12 2017

@author: stangen
"""

#Get the required packages
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#Import the ncfile, assign a file handle, indicate read-only
my_example_nc_file = '/home/disk/hot/stangen/Documents/GEFS/netcdf/pgbf_2017081400_00.nc'
fh = Dataset(my_example_nc_file, mode='r')
#Print the variables to see what we have available
print(fh.variables)

#Get lat,lon, 500mb height
lons = fh.variables['longitude'][:]
lats = fh.variables['latitude'][:]
Z500 = fh.variables['HGT_500mb'][:]
Z5002 = fh.variables['HGT_500mb'][0,:,:]

time = fh.variables['time'][:]

Z500_units = fh.variables['HGT_500mb'].units

fh.close()


#print(Z500_00Z)

Z500m = np.mean(Z500)

#Get some parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()


m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=0,urcrnrlon=360,lat_ts=20,resolution='c')

# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

# Plot Data
cs = m.pcolor(xi, yi, Z5002)
#cs = m.contour(xi, yi, Z5002,15,linewidths=1.5)

# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(0., 361., 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(Z500_units)

# Add Title
plt.title('DJF Maximum Temperature')

plt.show()