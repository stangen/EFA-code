#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 13:19:47 2018

@author: stangen
"""

import os
from ecmwfapi import ECMWFDataServer
from datetime import datetime
server = ECMWFDataServer()

#----------------Change these accordingly--------------------------------------
start_date = datetime(2013,4,1,0)
ens = 'ncep'
#------------------------------------------------------------------------------

ens_dict = {
            'ecmwf' : 'ecmf',
            'ncep' : 'kwbc',
            'eccc' : 'cwao',
            'jma' : 'rjtd'       
            }

save_dir = '/home/disk/hot/stangen/Documents/tigge_ensembles/orography/'

#Create directories if they don't yet exit
if (os.path.isdir(save_dir)):
    pass
else:
    os.mkdir(save_dir)

def retreive_tigge_data(date):
    
    y = date.strftime('%Y')
    m = date.strftime('%m')
    d = date.strftime('%d')
    h = date.strftime('%H')
    date = y+'-'+m+'-'+d               
    target = save_dir+date+'_'+h+'_'+ens+'.nc'
    tigge_pf_sfc_request(date, h, target)
            
def tigge_pf_sfc_request(date, h, target):
    server.retrieve({
        "class": "ti",
        "dataset": "tigge",
        "date": date,
        "expver": "prod",
        "grid": "0.5/0.5",
        "area": "90/-180/0/179.5",
        "levtype": "sfc",
        "origin": ens_dict[ens],
        "format": "netcdf",
        "param": "228002",
        "step": "0",
        "target": target,
        "time": h,
        "type": "cf",        
    })
    
retreive_tigge_data(start_date)