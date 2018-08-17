#!/usr/bin/env python

#Import Python Libraries
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import sys
import os
from datetime import datetime
import EFA.duplicate_madaus.efa_functions as ef

#------------------Change these------------------------------------------------
#date (year, month, day, hour)
date = datetime(2015,11,12,12)
#gridded or madis observations?
ob_type = 'gridded' 
#varable type of observation
variable = 'IVT'

#-------These are if gridded obs-------
#ensemble type
ens = 'eccc'
#grid- western lat, eastern lat, southern lon, northern lon, degree spacing
grid = '-180_180_90_0_3'

#plot values with colors, or just locations?
values = 'colors' #'colors' or 'locations'

#Set lat/lng bounding box for basemap plot
minLat = 0
maxLat = 85
minLng = -180
maxLng = 180

#save figure
save = False
#------------------------------------------------------------------------------

#dictionaries for plotting
obtype_dict = {
        'madis' : 'Metar & Maritime',
        'gridded' : 'Gridded'
        }

var_dict = {
        'T2M' : '2-Meter Temperature',
        'ALT' : 'Altimeter',
        'IWV' : 'Integrated Water Vapor',
        'IVT' : 'Integrated Vapor Transport'
        }

var_unit = {
        'T2M' : '(K)',
        'ALT' : 'hPa',
        'IWV' : 'mm',
        'IVT' : 'kg/m/s'
        
        }

save_dir = '/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/plots/'

#get date strings from datetime object for accessing file
y, m, d, h = ef.dt_str_timedelta(date)
dstr = y+m+d
time = h+'00'


#Open MADIS text data file
#f = open('/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/plots/2013040106obs_allmaritime.txt','r')
if ob_type == 'madis':
    f = open("/home/disk/hot/stangen/Documents/surface_obs/MADIS/"+dstr[0:6]+"/combined_"+variable+"/"+variable+"_"+dstr+"_"+time+".txt","r")
elif ob_type == 'gridded':
    f = open("/home/disk/hot/stangen/Documents/gridded_obs/"+ens+"/"+dstr[0:6]+"/"+variable+"/"+dstr+"_"+time+"_"+grid+".txt","r")

lines = f.readlines()
f.close()

#variable to count how many observations are in North America and Europe
#this was to show that EFA results using MADIS obs probably show how EFA
#affects these regions specifically, instead of a global average.
N_A=0
Euro=0
#Read CSV data from file and extract data within bounding box
lats = []; lngs = []; var = []
for l in lines:
    temp = l.split(',')
    if ((float(temp[1]) >= minLat) and (float(temp[1]) <= maxLat) and (float(temp[2]) >= minLng) and (float(temp[2]) <= maxLng)):
        lats.append(temp[1])
        lngs.append(temp[2])
        var.append(temp[5])
    if((float(temp[1])>= 15) and (float(temp[2]) <=-45) and (float(temp[2]) >=-170)):
        N_A+=1
    if((float(temp[1])>=30) and (float(temp[2]) >=-15) and (float(temp[2]) <=45)):
        Euro+=1

#Convert list to numpy array
lats = np.float64(lats)
lngs = np.float64(lngs)
var = np.float64(var)

#set min/max for colorbar
minn = round(min(var)-2,0)
maxx = round(max(var)+2,0)

fig = plt.figure(figsize=(18,12))
ax1 = plt.subplot(111)
#Generate basemap
m1=Basemap(projection='cyl',llcrnrlat=minLat,urcrnrlat=maxLat,
   llcrnrlon=minLng,urcrnrlon=maxLng,resolution='c')

#Draw Latitude Lines
#labels[left,right,top,bottom]  1=True 0=False
parallels = np.arange(-90,90,2)
m1.drawparallels(parallels,labels=[0,0,0,0],fontsize=10,linewidth=0.5)

meridians = np.arange(0,360,2)
m1.drawmeridians(meridians,labels=[0,0,0,0],fontsize=10,linewidth=0.5)

#Draw Map
m1.drawmapboundary(fill_color='aqua')
m1.drawcoastlines()
m1.drawcountries()
m1.drawstates()
m1.fillcontinents(color='coral',lake_color='aqua')

#Convert lat/lon to mercator coordinates
xalti,yalti = m1(lngs,lats)

if values == 'locations':
    m1.scatter(xalti,yalti,c='b',zorder=4,s=10,vmin=minn,vmax=maxx,marker="o",edgecolor='k')
elif values == 'colors':
    m1.scatter(xalti,yalti,c=var,zorder=4,s=20,vmin=minn,vmax=maxx,marker="o",edgecolor='k',cmap=plt.cm.gist_rainbow)
    #Plot colorbar
    cbar = plt.colorbar(fraction=0.023)
    cbar.ax.set_title(var_unit[variable],y=1.02)

#Set plot and colorbar titles
plt.title(obtype_dict[ob_type]+" "+var_dict[variable]+" Obs, "+time+" UTC "+dstr[4:6]+"/"+dstr[6:8]+"/"+dstr[0:4],weight='bold',fontsize=16)
#Save figure
if save == True:
    plt.savefig(save_dir+"_"+ob_type+"_"+variable+"_"+values+".png",frameon=False,bbox_inches='tight')

print('Observations in North America: ',N_A)
print('Observations in Europe: ',Euro)
