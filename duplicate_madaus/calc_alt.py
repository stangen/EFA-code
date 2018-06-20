#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 11:18:25 2018

@author: stangen
"""

import netCDF4 

#variables required for calculating altimeter:
elev = 4226*.3048#2822*.3048#1288#4226*.3048#860.1456 # in meters
station_alt = 1019.98 #station altimeter in inches
sp_mb = 864.49#907.29#864.49# #surface pressure in mb (hPa)
sp_in = sp_mb*.029528744#.02953#

#Method 1 of calculating altimeter:
#upon looking at actual observations, it appears this is the formula
#the NWS uses to convert surface pressure to altimeter setting:
m1 = (sp_in/((288-0.0065*elev)/288)**5.2561)/.029528744

m1_2 = sp_mb/((288-0.0065*elev)/288)**5.2561

#method 2 of calculating altimeter: 
m2 = (sp_mb-.3)*(1+(.0065*elev/288)*(1013.25/(sp_mb-.3))**0.190284)**(1/0.190284)

m3 = (sp_mb)*(1+(.0065*elev/288)*(1013.25/(sp_mb))**0.190284)**(1/0.190284)

m4 = (sp_mb)*(1+(.0065*elev/288.15)*(1013.25/(sp_mb))**0.190284)**(1/0.190284)


base_dir = "/home/disk/hot/stangen/Documents/surface_obs/MADIS/201806/raw/metar.20180614_0000"
f1 = netCDF4.Dataset(base_dir,"r")
variableNames = f1.variables.keys()
alt = f1.variables['altimeter'][:]
stn_name = f1.variables['stationName'][:,0:4] # Name of MADIS observation station
epoch = f1.variables['timeObs'][:] # Time of MADIS Observation
elevs = f1.variables['elevation'][:] # Elevation of MADIS observation 




stns = []
for n in stn_name:
    	#decode list of bytes and join them to get the station identifier string
        stns.append(b''.join(n).decode('utf-8'))

sta_alts = []
sta_times = []
sta_elev = []
for i,s in enumerate(stns):
    if s =='KSLC':
        sta_alts.append(alt[i])
        sta_times.append(epoch[i])
        sta_elev.append(elevs[i])
    
#slp = 
