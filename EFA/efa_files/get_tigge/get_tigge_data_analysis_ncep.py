#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:53:36 2017

@author: stangen
"""

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

vardict = {'10' : 'oct',
           '11' : 'nov',
           '12' : 'dec',
           '01' : 'jan',
           '02' : 'feb',
           '03' : 'mar'
           }

# Create the year/month/day string for use in code- just change these for different months!
y = '2016'
m = '10'
d = range(1,32)
#------------------------------------------------------------------------------------------

date_string = []
for day in d:
    day = format(day, '02d')
    day = str(day)
    date_string.append('%s-%s-%s' % (y, m, day))
print(date_string)
date_str = ""
for i in date_string:
    date_str += str(i) + "/"
date_str = date_str[:-1]
print(date_str)

def retreive_tigge_data():
    target = '/home/disk/hot/stangen/Documents/ensembles/analysis/rawmonths/%s%s/%s%s_ncep_sfc.nc' % (vardict[m], y, vardict[m], y)
    print(target)
    tigge_cf_sfc_request(date_str, target)
    target = '/home/disk/hot/stangen/Documents/ensembles/analysis/rawmonths/%s%s/%s%s_ncep_pl.nc' % (vardict[m], y, vardict[m], y)
    print(target)
    tigge_cf_pl_request(date_str, target)
             
#full ecmwf from 0-240 hr should be about 1 Gb for 2 variables (t2m and pwat)
#ECMWF sfc: 167 is t2m, 136 is total column water
def tigge_cf_sfc_request(date_str, target):
    server.retrieve({
        "class": "ti",
        "dataset": "tigge",
        "date": date_str,
        "expver": "prod",
        "grid": "0.5/0.5",
        "area": "90/-180/0/179.5",
        "levtype": "sfc",
        "origin": "kwbc",
        "format": "netcdf",
        "param": "136/167",
        "step": "0",
        "target": target,
        "time": "00/06/12/18",
        "type": "cf",        
    })
    
# 156 is 500 mb height
def tigge_cf_pl_request(date_str, target):
    server.retrieve({
        "class": "ti",
        "dataset": "tigge",
        "date": date_str,
        "expver": "prod",
        "grid": "0.5/0.5",
        "area": "90/-180/0/179.5",
        "levelist": "500",
        "levtype": "pl",
        "origin": "kwbc",
        "format": "netcdf",
        "param": "156",
        "step": "0",
        "target": target,
        "time": "00/06/12/18",
        "type": "cf",        
    })
    
retreive_tigge_data()