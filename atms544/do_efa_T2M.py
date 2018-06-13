#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 12:10:56 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
import numpy as np
from nicks_files.operational_cfsv2 import get_cfsv2_ensemble
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import efa_files.cfs_utilities_st as ut
import efa_files.madis_example.madis_utilities as mt
import time
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from efa_xray.state.ensemble import EnsembleState
from efa_xray.observation.observation import Observation
from efa_xray.assimilation.ensrf import EnSRF
import efa_xray.postprocess.postprocess as pp
from do_efa import Run_efa



ensemble_type = 'ecmwf'
#change this later
datestr = '20130401_0000'
year = datestr[0:4]
month = datestr[4:6]
day = datestr[6:8]
hour = datestr[9:11]
variables = ['T2M', 'ALT', 'P6HR', 'TCW']
ob_type = 'T2M'

var_dict = {
        'ALT' : 'alts',
        'T2M' : 'temp'
        }

m_dict = {
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

efa = Run_efa(datestr,ensemble_type,variables,var_dict[ob_type])
statecls, lats, lons, elevs = efa.load_data()

obs = efa.load_obs(forecast_hour=6)

observations = []
for ob in obs:
    ob_split = ob.split(',')
    #get the lat/lon of the station
    ob_lat = float(ob_split[1])
    ob_lon = float(ob_split[2])
    #get longitude positive-definite- ie -130 lon is 230 E
    if ob_lon < 0:
        ob_lon = ob_lon + 360
    ob_elev = float(ob_split[3])
    ob_time = float(ob_split[4])
    time2 = datetime.utcfromtimestamp(ob_time)
    ob_value = float(ob_split[5])
    
    #call function to check elevation to see whether or not to assimilate ob
    TorF = efa.check_elevation(lats,lons,elevs,ob,ob_lat,ob_lon,ob_elev)
    
    obser = Observation(value=ob_value, time=time2,
                    lat=ob_lat,lon=ob_lon, obtype=ob_type, localize_radius=1000,
                    assimilate_this=TorF, error=1)
    observations.append(obser)

# Put the state class object and observation objects into EnSRF object
assimilator = EnSRF(statecls, observations, loc='GC')

# Update the prior with EFA- post_state is an EnsembleState object
post_state, post_obs = assimilator.update()

state=post_state
outdir = '/home/disk/hot/stangen/Documents/atms544/posterior/'+ensemble_type+'/'+m_dict[month]+year+'/'

#Create output directories if they don't yet exit
if (os.path.isdir(outdir)):
    pass
else:
    os.makedirs(outdir)
    
outfile = outdir+year+'-'+month+'-'+day+'_'+hour+'_'+ensemble_type+'T2M.nc'

tunit='seconds since 1970-01-01'
# Write ensemble forecast to netcdf
with Dataset(outfile,'w') as dset:
        dset.createDimension('time',None)
        dset.createDimension('lat',state.ny())
        dset.createDimension('lon',state.nx())
        dset.createDimension('ens',state.nmems())
        dset.createVariable('time','i4',('time',))
        dset.createVariable('lat','f8',('lat',))
        dset.createVariable('lon','f8',('lon'))
        dset.createVariable('ens','i4',('ens',))
        dset.variables['time'].units = tunit
        dset.variables['lat'].units = 'degrees_north'
        dset.variables['lon'].units = 'degrees_east'
        dset.variables['ens'].units = 'member_number'
        dset.variables['time'][:] = date2num(state.ensemble_times(),tunit)
        dset.variables['lat'][:] = state['lat'].values[:,0]
        dset.variables['lon'][:] = state['lon'].values[0,:]
        dset.variables['ens'][:] = state['mem'].values
        for var in state.vars():
            print('Writing variable {}'.format(var))
            dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
            dset.variables[var].units = ut.get_units(var)
            dset.variables[var][:] = state[var].values

#    
#    
#    #create observations
#    ob1 = Observation(value=5900, time=datetime(2017,9,6,12),lat=24.55,lon=278.21,
#                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
#                  error=10)