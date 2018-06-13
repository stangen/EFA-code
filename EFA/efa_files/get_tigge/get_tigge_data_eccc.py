#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:51:05 2017

@author: stangen
"""
import os
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

# Create the year/month/day string for use in code
#10/1/16-10/3/16 missing for ECCC

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

y = '2013'
m = '04'
d = range(1,31)
date_string = []
for day in d:
    day = format(day, '02d')
    day = str(day)
    date_string.append('%s-%s-%s' % (y, m, day))
print(date_string)

save_dir = '/home/disk/hot/stangen/Documents/ensembles/eccc/%s%s/' % (vardict[m], y)

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
             target = save_dir+date+'_'+time+'_eccc_T_SP.nc' 
             tigge_pf_sfc_request(date, time, target)
             #target = save_dir+date+'_'+time+'_eccc_pl.nc' #%s_%s_ecmcf_pl.nc' % (vardict[m], y, date, time)
             #tigge_pf_pl_request(date, time, target)
             
             
#full ecmwf from 0-240 hr should be about 1 Gb for 2 variables (t2m and pwat)
#ECMWF sfc: 167 is t2m, 136 is total column water, 134 is surface pressure, 228228 is total precipitation
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
        "origin": "cwao",
        "format": "netcdf",
        "param": "134/167",#/228228",
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
        "origin": "cwao",
        "format": "netcdf",
        "param": "156",
        "step": "0/TO/240/BY/6",
        "target": target,
        "time": time,
        "type": "pf",        
    })
    
retreive_tigge_data()