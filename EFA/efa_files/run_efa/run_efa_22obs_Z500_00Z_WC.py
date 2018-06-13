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
#Annette Island, AK  
ob1 = Observation(value=5820, time=datetime(2017,9,6,0),lat=55.05,lon=228.41,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Prince George, BC
ob2 = Observation(value=5890, time=datetime(2017,9,6,0),lat=53.9,lon=237.21,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Edmonton Stony Plain, AB
ob3 = Observation(value=5870, time=datetime(2017,9,6,0),lat=53.53,lon=245.9,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Port Hardy, BC
ob4 = Observation(value=5870, time=datetime(2017,9,6,0),lat=50.68,lon=232.64,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Quillayute
ob5 = Observation(value=5910, time=datetime(2017,9,6,0),lat=47.95,lon=235.45,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Spokane
ob6 = Observation(value=5940, time=datetime(2017,9,6,0),lat=47.68,lon=242.37,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Great Falls
ob7 = Observation(value=5920, time=datetime(2017,9,6,0),lat=47.46,lon=248.61,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Salem
ob8 = Observation(value=5930, time=datetime(2017,9,6,0),lat=44.91,lon=237.0001,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Medford
ob9 = Observation(value=5940, time=datetime(2017,9,6,0),lat=42.36,lon=237.14,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Boise
ob10 = Observation(value=5950, time=datetime(2017,9,6,0),lat=43.56,lon=243.79,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Riverton
ob11 = Observation(value=5920, time=datetime(2017,9,6,0),lat=43.06,lon=251.52,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Reno
ob12 = Observation(value=5920, time=datetime(2017,9,6,0),lat=39.56,lon=240.2,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Elko
ob13 = Observation(value=5950, time=datetime(2017,9,6,0),lat=40.86,lon=244.27,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Salt Lake City
ob14 = Observation(value=5940, time=datetime(2017,9,6,0),lat=40.77,lon=248.05,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Grand Junction
ob15 = Observation(value=5940, time=datetime(2017,9,6,0),lat=39.11,lon=251.47,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Oakland
ob16 = Observation(value=5880, time=datetime(2017,9,6,0),lat=37.73,lon=237.79,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Vandenburg AFB, CA
ob17 = Observation(value=5860, time=datetime(2017,9,6,0),lat=34.75,lon=239.44,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Las Vegas
ob18 = Observation(value=5930, time=datetime(2017,9,6,0),lat=36.05,lon=244.82,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Flagstaff
ob19 = Observation(value=5940, time=datetime(2017,9,6,0),lat=35.23,lon=248.18,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#San Diego
ob20 = Observation(value=5900, time=datetime(2017,9,6,0),lat=32.85,lon=242.88,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Phoenix
ob21 = Observation(value=5910, time=datetime(2017,9,6,0),lat=33.45,lon=248.05,
                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
                  error=10)
#Tucson
ob22 = Observation(value=5920, time=datetime(2017,9,6,0),lat=32.23,lon=249.04,
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
observations.append(ob11)
observations.append(ob12)
observations.append(ob13)
observations.append(ob14)
observations.append(ob15)
observations.append(ob16)
observations.append(ob17)
observations.append(ob18)
observations.append(ob19)
observations.append(ob20)
observations.append(ob21)
observations.append(ob22)


# Put the state class object and observation objects into EnSRF object
assimilator = EnSRF(statecls, observations, loc='GC')
print(assimilator)

# Update the prior with EFA- post_state is an EnsembleState object
post_state, post_obs = assimilator.update()

state=post_state
outfile = '/home/disk/hot/stangen/Documents/GEFS/posterior/' + \
            '2017090600/2017090600_21mem_10days_22obs_Z500_00Z_WC_efa.nc'
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