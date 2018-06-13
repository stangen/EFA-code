#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:53:36 2017

@author: stangen
"""

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

# Create the year/month/day string for use in code
y = '2016'
m = '12'
d = range(1,32)
date_string = []
for day in d:
    day = format(day, '02d')
    day = str(day)
    date_string.append('%s-%s-%s' % (y, m, day))
print(date_string)

def retreive_tigge_data():
    dates = date_string
    times = ['00', '12']
    for date in dates:
         for time in times:
             target = '/home/disk/hot/stangen/Documents/ensembles/ecmwf/dec2016/%s_%s_ecmwf_sfc.nc' % (date, time)
             tigge_pf_sfc_request(date, time, target)
#             target = 'ecmwf_pl_%s_%s.grb' % (date, time)
#             tigge_pf_pl_request(date, time, target)
             
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
        "number": "1/TO/50",
        "origin": "ecmf",
        "format": "netcdf",
        "param": "136/167",
        "step": "0/TO/240/BY/6",
        "target": target,
        "time": time,
        "type": "pf",        
    })
    
retreive_tigge_data()