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
   
# 10 sea level pressure observations, 3 hours after model is run. Accuracy is
# approx 0.01 in mercury, 33.8 Pa, rounded to 34
#Ketchikan  
ob1 = Observation(value=101320, time=datetime(2017,9,6,6),lat=55.356,lon=228.286, #101210
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Seatac
ob2 = Observation(value=101420, time=datetime(2017,9,6,6),lat=47.445,lon=237.686, #101280
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Spokane
ob3 = Observation(value=101620, time=datetime(2017,9,6,6),lat=47.621,lon=242.472, #101630
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Great Falls
ob4 = Observation(value=102000, time=datetime(2017,9,6,6),lat=47.473,lon=248.618, #102140
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Portland
ob5 = Observation(value=101340, time=datetime(2017,9,6,6),lat=45.595,lon=237.391, #101180
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Medford/Rogue Valley
ob6 = Observation(value=101400, time=datetime(2017,9,6,6),lat=42.375,lon=237.123, #101220
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Boise
ob7 = Observation(value=101370, time=datetime(2017,9,6,6),lat=43.567,lon=243.759, #101360
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Riverton
ob8 = Observation(value=102220, time=datetime(2017,9,6,6),lat=43.062,lon=251.553, #102350
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Prince George, BC
ob9 = Observation(value=101693, time=datetime(2017,9,6,6),lat=53.929,lon=237.198, #101625
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Slave Lake, Alta, AB
ob10 = Observation(value=101880, time=datetime(2017,9,6,6),lat=55.3,lon=245.217, #101920
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Sacramento
ob11 = Observation(value=101250, time=datetime(2017,9,6,6),lat=38.700,lon=238.405, #101160
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Elko
ob12 = Observation(value=101390, time=datetime(2017,9,6,6),lat=40.824,lon=244.214, #101390
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Salt Lake City
ob13 = Observation(value=101270, time=datetime(2017,9,6,6),lat=40.770,lon=248.035, #101260
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Grand Junction
ob14 = Observation(value=101280, time=datetime(2017,9,6,6),lat=39.133,lon=251.461, #101170
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#San Francisco
ob15 = Observation(value=101480, time=datetime(2017,9,6,6),lat=37.619,lon=237.634, #101430
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#LAX
ob16 = Observation(value=101400, time=datetime(2017,9,6,6),lat=33.938,lon=241.611, #101270
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Las Vegas
ob17 = Observation(value=101400, time=datetime(2017,9,6,6),lat=36.071,lon=244.837, #100890
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Flagstaff
ob18 = Observation(value=101710, time=datetime(2017,9,6,6),lat=35.144,lon=248.334, #101810
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Tucson
ob19 = Observation(value=101080, time=datetime(2017,9,6,6),lat=32.131,lon=249.044, #100910
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)
#Yuma
ob20 = Observation(value=100890, time=datetime(2017,9,6,6),lat=32.659,lon=245.407, #100770
                  obtype = 'MSLP', localize_radius=2000, assimilate_this=True,
                  error=34)



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

# Put the state class object and observation objects into EnSRF object
assimilator = EnSRF(statecls, observations, loc='GC')
print(assimilator)

# Update the prior with EFA- post_state is an EnsembleState object
post_state, post_obs = assimilator.update()

state=post_state
outfile = '/home/disk/hot/stangen/Documents/GEFS/posterior/' + \
            '2017090600/2017090600_21mem_10days_20MSLP_WC_efa.nc'
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