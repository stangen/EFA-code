#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 12:10:56 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
import numpy as np
from datetime import datetime
import EFA.efa_files.cfs_utilities_st as ut
import surface_obs.madis_example.madis_utilities as mt
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from EFA.efa_xray.observation.observation import Observation
from EFA.efa_xray.assimilation.ensrf import EnSRF
from EFA.duplicate_madaus.load_data import Load_data
from EFA.duplicate_madaus.efa_functions import check_elev



ensemble_type = 'ecmwf'
#change this later
start_date = datetime(2013,4,1,0)
end_date = datetime(2013,4,1,0)
hourstep=12

#All variables in the netCDF
variables = ['T2M','ALT']#['T2M', 'ALT', 'P6HR', 'TCW']
#the ob type of the observations we are assimilating
ob_type = 'ALT'
#the variables in the netCDF we want to update
update_var= ob_type

#a list of dates to loop through to load each forecast initialized on these dates
dates = mt.make_datetimelist(start_date,end_date,timestep = 3600*hourstep) 

for date in dates:
    y = date.strftime('%Y')
    m = date.strftime('%m')
    d = date.strftime('%d')
    h = date.strftime('%H')

    #initialize an instance of the Load_data class (load in data)
    efa = Load_data(date,ensemble_type,variables,ob_type)
    statecls, lats, lons, elevs = efa.load_netcdfs()
    #load in the obs file
    obs = efa.load_obs()
    
    observations = []
    #loop through each line in the text file (loop through each observation)
    for ob in obs:
        #this gets the observation information from the text file
        ob_dict = mt.get_ob_info(ob)
        #ob_split = ob.split(',')
        #get the lat/lon of the station
        #ob_lat = float(ob_split[1])
        #ob_lon = float(ob_split[2])
        #get longitude positive-definite- ie -130 lon is 230 E
        if ob_dict['lon'] < 0:
            ob_dict['lon'] = ob_dict['lon'] + 360
        #ob_elev = float(ob_split[3])
        #ob_time = float(ob_split[4])
        utctime = datetime.utcfromtimestamp(ob_dict['time'])
        #ob_value = float(ob_split[5])
        
        #call function to check elevation to see whether or not to assimilate ob
        TorF = check_elev(lats,lons,elevs,ob_dict['lat'],ob_dict['lon'],ob_dict['elev'])
        
#        obser = Observation(value=ob_value, time=utctime,
#                        lat=ob_lat,lon=ob_lon, obtype=ob_type, localize_radius=1000,
#                        assimilate_this=TorF, error=1)
#        observations.append(obser)
#    
#    # Put the state class object and observation objects into EnSRF object
#    assimilator = EnSRF(statecls, observations, loc='GC')
#    
#    # Update the prior with EFA- post_state is an EnsembleState object
#    post_state, post_obs = assimilator.update()
#    
#    state=post_state
#    outdir = '/home/disk/hot/stangen/Documents/posterior_ensembles/'+ensemble_type+'/'+y+m+'/'
#    
#    #Create output directories if they don't yet exit
#    if (os.path.isdir(outdir)):
#        pass
#    else:
#        os.makedirs(outdir)
#        
#    outfile = outdir+y+'-'+m+'-'+d+'_'+h+'_'+ensemble_type+'_'+ob_type+'.nc'
#    
#    tunit='seconds since 1970-01-01'
#    # Write ensemble forecast to netcdf
#    with Dataset(outfile,'w') as dset:
#            dset.createDimension('time',None)
#            dset.createDimension('lat',state.ny())
#            dset.createDimension('lon',state.nx())
#            dset.createDimension('ens',state.nmems())
#            dset.createVariable('time','i4',('time',))
#            dset.createVariable('lat',np.float32,('lat',))
#            dset.createVariable('lon',np.float32,('lon'))
#            dset.createVariable('ens','i4',('ens',))
#            dset.variables['time'].units = tunit
#            dset.variables['lat'].units = 'degrees_north'
#            dset.variables['lon'].units = 'degrees_east'
#            dset.variables['ens'].units = 'member_number'
#            dset.variables['time'][:] = date2num(state.ensemble_times(),tunit)
#            dset.variables['lat'][:] = state['lat'].values[:,0]
#            dset.variables['lon'][:] = state['lon'].values[0,:]
#            dset.variables['ens'][:] = state['mem'].values
#            for var in state.vars():
#                print('Writing variable {}'.format(var))
#                dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
#                dset.variables[var].units = ut.get_units(var)
#                dset.variables[var][:] = state[var].values

#    
#    
#    #create observations
#    ob1 = Observation(value=5900, time=datetime(2017,9,6,12),lat=24.55,lon=278.21,
#                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
#                  error=10)