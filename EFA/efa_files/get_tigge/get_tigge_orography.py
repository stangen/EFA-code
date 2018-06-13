#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 13:19:47 2018

@author: stangen
"""

import os
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

vardict = {
           '01' : 'jan',
           '02' : 'feb',
           '03' : 'mar',
           '04' : 'apr',
           '05' : 'may',
           '06' : 'jun',
           '07' : 'jul',
           '08' : 'aug',
           '09' : 'sep',
           '10' : 'oct',
           '11' : 'nov',
           '12' : 'dec'           
           }

# Create the year/month/day string for use in code
y = '2013'
m = '04'
d = range(1,3)
date_string = []
for day in d:
    day = format(day, '02d')
    day = str(day)
    date_string.append('%s-%s-%s' % (y, m, day))
print(date_string)

save_dir = '/home/disk/hot/stangen/Documents/ensembles/orography/'

#Create directories if they don't yet exit
if (os.path.isdir(save_dir)):
    pass
else:
    os.mkdir(save_dir)

def retreive_tigge_data():
    dates = date_string
    times = ['00']
    for date in dates:
        for time in times:            
            target = save_dir+date+'_'+time+'_ncep.nc' #% (vardict[m], y, date, time)
            tigge_pf_sfc_request(date, time, target)
            
def tigge_pf_sfc_request(date, time, target):
    server.retrieve({
        "class": "ti",
        "dataset": "tigge",
        "date": date,
        "expver": "prod",
        "grid": "0.5/0.5",
        "area": "90/-180/0/179.5",
        "levtype": "sfc",
        "origin": "kwbc",
        "format": "netcdf",
        "param": "228002",
        "step": "0",
        "target": target,
        "time": time,
        "type": "cf",        
    })
    
retreive_tigge_data()