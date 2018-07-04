#!/usr/bin/env python

#Import Python Libraries
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import numpy as np
import sys
import os

#Retrieve data from command line args
#dstr = sys.argv[1]
dstr = '20130401'
time = '0600'
var = 'alts'

save_dir = '/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/plots/'

#Open MADIS text data file
f = open('/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/plots/2013040100obs_allmaritime.txt','r')
#f = open("/home/disk/hot/stangen/Documents/EFA/surface_obs/MADIS/"+dstr[0:6]+"/combined_"+var+"/"+var+"_"+dstr+"_"+time+".txt","r")
#f = open("/home/disk/hot/stangen/Documents/EFA/surface_obs/MADIS/"+dstr[0:6]+"/raw/metar_alts_"+dstr+"_"+time+".txt","r")

lines = f.readlines()
f.close()

#Set lat/lng bounding box for basemap plot
minLat = 0
maxLat = 85
minLng = -180
maxLng = 180

N_A=0
Euro=0
#Read CSV data from file and extract data within bounding box
lats = []; lngs = []; alts = []
for l in lines:
    temp = l.split(',')
    if ((float(temp[1]) >= minLat) and (float(temp[1]) <= maxLat) and (float(temp[2]) >= minLng) and (float(temp[2]) <= maxLng)):
        lats.append(temp[1])
        lngs.append(temp[2])
        alts.append(temp[5])
    if((float(temp[1])>= 15) and (float(temp[2]) <=-45) and (float(temp[2]) >=-170)):
        N_A+=1
    if((float(temp[1])>=30) and (float(temp[2]) >=-15) and (float(temp[2]) <=45)):
        Euro+=1

#Convert list to numpy array
lats = np.float64(lats)
lngs = np.float64(lngs)
alts = np.float64(alts)

#set min/max for colorbar
minn = round(min(alts)-2,0)
maxx = round(max(alts)+2,0)

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

#Convert lat lon to x,y coordinates for plotting
m1.scatter(xalti,yalti,c='b',zorder=4,s=10,vmin=minn,vmax=maxx,marker="o")#,edgecolor='w',cmap=plt.cm.gist_rainbow)
#Plot colorbar
#cbar = plt.colorbar(fraction=0.023)
#Set plot and colorbar titles
plt.title("Metar & Stationary Maritime Altimeter Obs, "+time+" UTC "+dstr[4:6]+"/"+dstr[6:8]+"/"+dstr[0:4],weight='bold',fontsize=16)
#cbar.ax.set_title('(K)',y=1.02)
#Save figure
plt.savefig(save_dir+"madis_sfc_allmaritime.png",frameon=False,bbox_inches='tight')

print(N_A)
print(Euro)
