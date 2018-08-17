#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 14:03:03 2018

@author: stangen

Function to create a txt file with a time-series of observations for one location

"""

from surface_obs.madis_example.madis_utilities import make_datelist

#----------------_Change these as necessary------------------------------------
#Station name
sta_name = 'KSEA'
#Observation type
ob_type = 'T2M'
#------------------------------------------------------------------------------

save_dir = '/home/disk/hot/stangen/Documents/surface_obs/MADIS/one_station/'

s_date = '20130401_0000'
e_date = '20130511_1800'

#time-series to load observations
d_list = make_datelist(start_date=s_date, end_date=e_date)

#loop through each observation
for t in d_list:
    #strings for loading the appropriate ob file based on time
    y = t[0:4]
    m = t[4:6]
    d = t[6:8]
    h = t[8:10]
    obs_file = '/home/disk/hot/stangen/Documents/surface_obs/MADIS/'+y+m+'/combined_'+ob_type+'/'+ob_type+'_'+y+m+d+'_'+h+'00.txt'
    #read in the observations
    f1 = open(obs_file, 'r')
    obs = f1.readlines()
    #loop through each observation location in the observation file for a specific date
    for o in obs:
        if sta_name in o:
            #Write to a new observation file
            f = open(save_dir+sta_name+'_'+ob_type+'.txt', 'a')
            f.write(o)
            f.close()