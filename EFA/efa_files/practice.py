#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 13:04:42 2017

@author: stangen
"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
from nicks_files.operational_cfsv2 import get_cfsv2_ensemble
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import nicks_files.cfs_utilities as ut
import time
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from efa_xray.state.ensemble import EnsembleState
from efa_xray.observation.observation import Observation
from efa_xray.assimilation.ensrf import EnSRF

# Filepath of the precip analysis
#precip_analysis_path = '/home/disk/hot/stangen/Documents/ensembles/ecmwf/nov2016/2016-11-13_00_ecmwf_sfc.nc'
#field = []
#with Dataset(precip_analysis_path,'r') as ecmwf:
#    #print(ecmwf.variables)
#    field = ecmwf.variables['t2m'][:,:,:,:]
#ncep_path = '/home/disk/hot/stangen/Documents/ensembles/ncep/nov2016/2016-11-13_00_ncep_sfc.nc'
#with Dataset(ncep_path, 'r') as ncep:
#    field2 = ncep.variables['t2m'][:,:,:,:]
    
#field3 = np.stack((field, field2), axis = 1 )
#field3 = field + field2

#    field = ncdata.variables['t2m'][:,:,:,:]
#field = np.swapaxes(field, 1, 3)
#field = np.swapaxes(field, 1,2)
#field2 = np.swapaxes(field2, 1, 3)
#field2 = np.swapaxes(field2, 1,2)
#     print(field)      
#    print(ncdata.variables['t2m'][:,:,:,:])
#    tunit = ncdata.variables['time'].units
#    ftimes = num2date(ncdata.variables['time'][:],tunit)
#    print(tunit)
#    print(ftimes)
#    print(ftimes[40])

#path = '/home/disk/hot/stangen/Documents/GEFS/ensembles/2017090600/2017090600_21mem_10days.nc'
#with Dataset(path,'r') as ncdata:
#    #print(ncdata.variables)
#
#field3 = np.concatenate((field, field2), axis=3)

indir = '/home/disk/hot/stangen/Documents/ensembles/ecmwf/oct2016/2016-10-01_00_ecmwf_sfc_ivt.nc'
with Dataset(indir, 'r') as ncdata:
    print(ncdata.variables)