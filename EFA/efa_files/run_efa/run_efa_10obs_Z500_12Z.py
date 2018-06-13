#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 09:56:29 2017

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
import numpy as np
from nicks_files.operational_cfsv2 import get_cfsv2_ensemble
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import efa_files.cfs_utilities_st as ut
import time
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from efa_xray.state.ensemble import EnsembleState
from efa_xray.observation.observation import Observation
from efa_xray.assimilation.ensrf import EnSRF

# directory where the ensemble of all times is
infile = '/home/disk/hot/stangen/Documents/GEFS/ensembles' + \
            '/2017090600/2017090600_21mem_10days.nc'

# only variable I am using is 500mb height            
vrbls=['Z500','T500','RH500','U500','V500','Z700','T700','RH700','U700','V700', \
       'Z850','T850','RH850','U850','V850','Z925','T925','RH925','U925','V925', \
       'Z1000','T1000','RH1000','U1000','V1000','T2M','RH2M','U10M','V10M', \
       'PWAT','MSLP','P6HR']

# loading/accessing the netcdf data            
with Dataset(infile,'r') as ncdata:
    times = ncdata.variables['time']
    ftimes = num2date(times[:],
                      times.units)
    lats = ncdata.variables['lat'][:]
    lons = ncdata.variables['lon'][:]
    mems = ncdata.variables['ens'][:]
    #print(ncdata.variables)
    
    # storing the variable data in a dict (state?)
    allvars = {}
    for var in vrbls:
        allvars[var] = (['validtime','y','x','mem'],
                        ncdata.variables[var][:])
lonarr, latarr = np.meshgrid(lons, lats)

# Package into an EnsembleState object knowing the state and metadata
statecls = EnsembleState.from_vardict(allvars,
                                      {'validtime' : ftimes,
                                       'lat' : (['y','x'], latarr),
                                       'lon' : (['y','x'], lonarr),
                                       'mem' : mems,
                                       })

# Creating 2 dummy obs to test EFA- if exact coordinates in the lat/lon,
# interpolate will fail.   
#key west  
ob1 = Observation(value=5900, time=datetime(2017,9,6,12),lat=24.55,lon=278.21,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Miami
ob2 = Observation(value=5900, time=datetime(2017,9,6,12),lat=25.75,lon=279.62,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Tampa Bay
ob3 = Observation(value=5870, time=datetime(2017,9,6,12),lat=27.70,lon=277.6,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Tallahassee
ob4 = Observation(value=5850, time=datetime(2017,9,6,12),lat=30.45,lon=275.7,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Jacksonville
ob5 = Observation(value=5860, time=datetime(2017,9,6,12),lat=30.50,lon=278.3,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#BMX
ob6 = Observation(value=5780, time=datetime(2017,9,6,12),lat=33.16,lon=273.24,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Charleston
ob7 = Observation(value=5840, time=datetime(2017,9,6,12),lat=32.90,lon=279.97,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#LIX Slidell
ob8 = Observation(value=5850, time=datetime(2017,9,6,12),lat=30.34,lon=270.17,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Jackson
ob9 = Observation(value=5820, time=datetime(2017,9,6,12),lat=32.32,lon=269.92,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Nashville
ob10 = Observation(value=5690, time=datetime(2017,9,6,12),lat=36.25,lon=273.43,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)

# Put the observations into a list for EnSRF
observations = []
observations.append(ob1)
observations.append(ob2)
observations.append(ob3)
observations.append(ob4)
observations.append(ob5)
observations.append(ob6)
observations.append(ob7)
observations.append(ob8)
observations.append(ob9)
observations.append(ob10)

# Put the state class object and observation objects into EnSRF object
assimilator = EnSRF(statecls, observations, loc='GC')
print(assimilator)

# Update the prior with EFA- post_state is an EnsembleState object
post_state, post_obs = assimilator.update()

state=post_state
outfile = '/home/disk/hot/stangen/Documents/GEFS/posterior/' + \
            '2017090600/2017090600_21mem_10days_Z500_12Z_efa.nc'
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
            dset.createVariable(var, 'f8', ('time','lat','lon','ens',))
            dset.variables[var].units = ut.get_units(var)
            dset.variables[var][:] = state[var].values









#Get the required packages
#from netCDF4 import Dataset
#import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#
##Import the ncfile, assign a file handle, indicate read-only
#my_example_nc_file = '/home/disk/hot/stangen/Documents/GEFS/ensembles/2017081400_21mem_1days.nc'
#fh = Dataset(my_example_nc_file, mode='r')
##Print the variables to see what we have available
#print(fh.variables)