#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:53:36 2017

@author: stangen
"""
import os
from ecmwfapi import ECMWFDataServer
from datetime import datetime
import surface_obs.madis_example.madis_utilities as mt
server = ECMWFDataServer()

###---------Change these variables accordingly---------------------------------
#start and end date to get ensembles. 
start_date = datetime(2015,11,10,0) #YYYY,m,d,h
end_date = datetime(2015,11,17,12)
start_fh = '0' #first forecast hour of each forecast
end_fh = '54' #last forecast hour of each forecast
fh_step = '6' #increment of forecast hour to retrieve
ens = 'ncep' #ensemble
param = ['Q','U','V'] #parameters to get
surface = False #True if surface variables, false if aloft (500mb?)
levels = '850'
hourstep = 12 #how often you want a new forecast initialization, usually 12 hr
#------------------------------------------------------------------------------


ens_dict = {
            'ecmwf' : 'ecmf',
            'ncep' : 'kwbc',
            'eccc' : 'cwao',
            'jma' : 'rjtd'       
            }

ens_num_dict = {
                'ecmwf' : '1/TO/50',
                'jma' : '1/TO/50',
                'ncep' : '1/TO/20',
                'eccc' : '1/TO/20'
                }

param_dict = {
            'T2M' : '167', #2m temperature
            'SP' : '134', #surface pressure
            'TCW' : '136', #total column water
            'PCP' : '228228', #total precipitation
            'MSLP' : '151', #mean sea level pressure
            'Z500' : '156', #500 mb height
            'Q' : '133', #specific humidity
            'U' : '131', #u-component of wind
            'V' :' 132' #v-component of wind
            }   

#construct the parameter string to be fed into tigge
#construct the variable string for saving the ensembles
param_string = ''
var_string = ''
for p in param[:-1]:  
    param_string = param_string+param_dict[p]+'/' 
    var_string = var_string+p+'_'
param_string = param_string+param_dict[param[-1]]
var_string = var_string+param[-1]

#construct the timestep string to be fed into tigge
timestep = start_fh+'/TO/'+end_fh+'/BY/'+fh_step

#a list of dates to loop through to load each forecast initialized on these dates
dates = mt.make_datetimelist(start_date,end_date,hourstep)    
    
def retreive_tigge_data(date):
    y = date.strftime('%Y')
    m = date.strftime('%m')
    d = date.strftime('%d')
    h = date.strftime('%H')
    save_dir = '/home/disk/hot/stangen/Documents/tigge_ensembles/'+ens+'/'+y+m+'/'
    
    #Create directories if they don't yet exit
    if (os.path.isdir(save_dir)):
        pass
    else:
        os.makedirs(save_dir)
        
    date = y+'-'+m+'-'+d      
    if surface == True:
        target = save_dir+date+'_'+h+'_'+ens+'_'+var_string+'_sfc.nc'
        tigge_pf_sfc_request(date, h, target)
    elif surface == False:
        target = save_dir+date+'_'+h+'_'+ens+'_'+var_string+'_pl.nc' 
        tigge_pf_pl_request(date, h, target)
             
#full ecmwf from 0-240 hr should be about 1 Gb for 2 variables (t2m and pwat)
#ECMWF sfc: 167 is t2m, 136 is total column water, 134 is surface pressure, 228228 is total precipitation
def tigge_pf_sfc_request(date, h, target):
    server.retrieve({
        "class": "ti",
        "dataset": "tigge",
        "date": date,
        "expver": "prod",
        "grid": "0.5/0.5",
        "area": "90/-180/0/179.5",
        "levtype": "sfc",
        "number": ens_num_dict[ens],
        "origin": ens_dict[ens],
        "format": "netcdf",
        "param": param_string,
        "step": timestep,
        "target": target,
        "time": h,
        "type": "pf",        
    })

# 156 is 500 mb height
def tigge_pf_pl_request(date, h, target):
    server.retrieve({
        "class": "ti",
        "dataset": "tigge",
        "date": date,
        "expver": "prod",
        "grid": "0.5/0.5",
        "area": "90/-180/0/179.5",
        "levelist": levels,
        "levtype": "pl",
        "number": ens_num_dict[ens],
        "origin": ens_dict[ens],
        "format": "netcdf",
        "param": param_string,
        "step": timestep,
        "target": target,
        "time": h,
        "type": "pf",        
    })
    
for date in dates:
    retreive_tigge_data(date)