#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 17:05:53 2017

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
# Some missing data towards the end of October
y = '2013'
m = '04'
d = range(1,31)
date_string = []
for day in d:
    day = format(day, '02d')
    day = str(day)
    date_string.append('%s-%s-%s' % (y, m, day))
print(date_string)

save_dir = '/home/disk/hot/stangen/Documents/ensembles/ncep/%s%s/' % (vardict[m], y)

#Create directories if they don't yet exit
if (os.path.isdir(save_dir)):
    pass
else:
    os.mkdir(save_dir)

def retreive_tigge_data():
    dates = date_string
    times = ['00', '12']
    for date in dates:
         for time in times:
             target = save_dir+date+'_'+time+'_ncep_T_SP.nc' 
             tigge_pf_sfc_request(date, time, target)
             #target = save_dir+date+'_'+time+'_ncep_pl.nc' #%s_%s_ncep_pl.nc' % (vardict[m], y, date, time)
             #tigge_pf_pl_request(date, time, target)
             
#full ecmwf from 0-240 hr should be about 1 Gb for 2 variables (t2m and pwat)
#ECMWF sfc: 167 is t2m, 136 is total column water
def tigge_pf_sfc_request(date, time, target):
    server.retrieve({
        "class": "ti",
        "dataset": "tigge",
        "date": date,
        "expver": "prod",
        "grid": "0.5/0.5",
        "area": "90/-180/0/179.5",
        "levtype": "sfc",
        "number": "1/TO/20",
        "origin": "kwbc",
        "format": "netcdf",
        "param": "134/167",
        "step": "0/TO/240/BY/6",
        "target": target,
        "time": time,
        "type": "pf",        
    })
    
# 156 is 500 mb height
def tigge_pf_pl_request(date, time, target):
    server.retrieve({
        "class": "ti",
        "dataset": "tigge",
        "date": date,
        "expver": "prod",
        "grid": "0.5/0.5",
        "area": "90/-180/0/179.5",
        "levelist": "500",
        "levtype": "pl",
        "number": "1/TO/20",
        "origin": "kwbc",
        "format": "netcdf",
        "param": "156",
        "step": "0/TO/240/BY/6",
        "target": target,
        "time": time,
        "type": "pf",        
    })
    
retreive_tigge_data()