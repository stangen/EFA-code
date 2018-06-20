#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 14:47:07 2017

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

map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
#map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))

# compute native map projection coordinates of lat/lon grid.
lon, lat = np.meshgrid(lons, lats)
x, y = map(lon, lat)

cs = map.contour(x,y,Z5002,15,linewidths=1.5)

# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(Z500_units)

plt.title('contour lines over filled continent background')
plt.show()