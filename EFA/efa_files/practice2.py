#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 08:00:39 2017

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
from datetime import datetime, timedelta
from subprocess import check_output
import numpy as np
import efa_files.cfs_utilities_st as ut
from efa_xray.state.ensemble import EnsembleState

#indir = '                                                                                 '
#indir = '/home/disk/hot/stangen/Documents/ensembles/eccc/dec2016/2016-12-31_12_eccc_sfc.nc'
#indirpf = '/home/disk/hot/stangen/Documents/ensembles/eccc/dec2016/2016-12-31_12_eccc_pf.nc'
#indir2 = '/home/disk/hot/stangen/Documents/ensembles/ecmwf/dec2016/2016-12-31_12_sfc.nc'
#indir2pf = '/home/disk/hot/stangen/Documents/ensembles/ecmwf/dec2016/2016-12-31_12_pf.nc'
#
#
#if indir[0:50] != indirpf[0:50]:
#    print('it will run the try loop')
#else:
#    print('skip the try loop')
#    #continue
#    

#indir = '/home/disk/hot/stangen/Documents/ensembles/analysis/oct2016-mar2017/oct2016_ncep_anl_pl.nc'
#indir2 = '/home/disk/hot/stangen/Documents/ensembles/analysis/oct2016-mar2017/oct2016_ncep_anl_sfc.nc' 

#with Dataset(indir, 'r') as ncdata:
#    print(len(ncdata.variables['time']))
#    print(ncdata.variables)
#    print(ncdata.variables['time'][2]) #in hours since 1900-01-01 00Z
    
#with Dataset(indir2, 'r') as ncdata:
#    print(ncdata.variables)

    # Here we have a dictionary translating some generic variable names to
    # what they correspond to in the CFSv2 netcdf files
vardict = {'Z500' : 'HGT_500mb',
           'gh' : 'Z500',
           't2m' : 'T2M',
           'tcw' : 'TCW',
           'time' : 'time',
           'lat' : 'latitude',
           'lon' : 'longitude',
           '10' : 'oct',
           '11' : 'nov',
           '12' : 'dec'
           }
outfile = '/home/disk/hot/stangen/Documents/ensembles/analysis/combined/oct-mar.nc' 

# directories for the ensembles
months = ['oct2016', 'nov2016', 'dec2016', 'jan2016', 'feb2016', 'mar2016']
memfiles = []
each_time_length = list()
state = []

indir = '/home/disk/hot/stangen/Documents/ensembles/analysis/rawmonths/oct2016-mar2017/*' 
# Get a list of filenames (each file is a different member)
command = 'ls -1a {}'.format(indir)
#print(command)
memfiles.extend(list(reversed(check_output([command],shell=True).split())))

print(memfiles)

ntimes = 0
priorens = '                                                                                  '
t2 = -1
# Count the number of ensemble members for allocation of state later
# try blocks are to filter out dates when TIGGE didn't have any data
#print('\nCreating superensemble for {}/{}/{} {}Z\n'.format(month,day,year,hour))
print('Counting ensemble members for allocation')
for t, times in enumerate(memfiles):
    times = times.decode('utf-8')
    # to not double-count ens members if from both sfc and pl
    print(times[0:90])
    if times[0:90] != priorens[0:90]:
        t2 += 1
        try:
            with Dataset(times, 'r') as ncdata:                   
                #keep track of number of members of each ensemble
                each_time_length.append(len(ncdata.variables['time']))
                #keep track of total length of all ensembles
                ntimes = ntimes + len(ncdata.variables['time'])
                print('times from {}: {}'.format(months[t2],each_time_length[t2]))
                print('total times: {}'.format(ntimes))
                #print(each_mem_length)
        except:
            print('Bad {} ensembles- not counting ensemble members'.format(months[t2]))

            pass
    priorens = times