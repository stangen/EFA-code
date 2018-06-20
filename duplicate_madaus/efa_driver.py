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
import EFA.duplicate_madaus.efa_functions as ef



ensemble_type = 'ecmwf'
#change this later
date = datetime(2013,4,1,0)
#end_date = datetime(2013,4,1,0)
#hourstep=12
#self_update = True
#All variables in the netCDF
variables = ['T2M','ALT']#['T2M', 'ALT', 'P6HR', 'TCW']
##do we only want to update each variable using its own obs?
#if self_update == True:
#    for var in variables:
#        ob_type = [var]
#        update_var = [var]
#the ob type of the observations we are assimilating
ob_type = ['ALT']
#the variables in the netCDF we want to update
update_var= ['ALT']
#are the observations only updating their corresponding variable, or
#are they updating all variables? 
last_ob == False
if update_var == ob_type and len(ob_type ==1):
    self_update = True #obs only update their corresponding vars
else:
    self_update = False #each ob updates all ensemble vars
#localization type
loc_type = 'GC'
#localization radius (for Gaspari-Cohn)
localize_radius = 1000

#assign the radius to the localization_radius variable in Observation
#create a string for saving to specific directories
if loc_type == 'GC':
    loc_rad = localize_radius
    loc_str = str(loc_rad)
else:
    loc_rad = None
    loc_str = '_stat_sig'
    
time_lag = False

#a list of dates to loop through to load each forecast initialized on these dates

y = date.strftime('%Y')
m = date.strftime('%m')
d = date.strftime('%d')
h = date.strftime('%H')

for o, o_type in enumerate(ob_type):
    efa = Load_data(date,ensemble_type,variables,o_type,update_var)
    if o == 0:
        #initialize an instance of the Load_data class (load in data)
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
        TorF = ef.check_elev(lats,lons,elevs,ob_dict['lat'],ob_dict['lon'],ob_dict['elev'])
        
#        obser = Observation(value=ob_dict['ob'], time=utctime,
#                        lat=ob_dict['lat'],lon=ob_dict['lon'], obtype=o_type, localize_radius=loc_rad,
#                        assimilate_this=TorF, error=1)
#        observations.append(obser)
#    
## Put the state class object and observation objects into EnSRF object
#assimilator = EnSRF(statecls, observations, loc=loc_type)
#
## Update the prior with EFA- post_state is an EnsembleState object
#post_state, post_obs = assimilator.update()
#
#state=post_state
#
##build the string of where to save the files
#outdir = '/home/disk/hot/stangen/Documents/posterior_ensembles/
##are the observations only updating their corresponding variable, or
##are they updating all variables? 
#if self_update == True:
#    outdir += 'ob_update_self/' 
#elif self_update == False:
#    outdir += 'ob_update_all/' 
##directory for the type of localization used       
#outdir += 'loc_'+loc_str+'/'
##if dealing with time-lagged ensembles, add another directory
#if time_lag == True:
#    outdir += 'time_lag/'
##add directories for the ens type and year & month
#outdir += ensemble_type+'/'+y+m+'/'
#
#    
##Create output directories if they don't yet exist
#if (os.path.isdir(outdir)):
#    pass
#else:
#    os.makedirs(outdir)
#
#outdir_date_ens = outdir+y+'-'+m+'-'+d+'_'+h+'_'+ensemble_type
        
# If we're only processing one variable
if self_update == True:
    for var in variables:
        checkfile = outdir_date_ens+'*'
        if os.path.isfile(checkfile):
            print('Appending to existing file!')
            # append to the existing file
            with Dataset(checkfile,'a') as dset:
                dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
                dset.variables[var].units = ut.get_units(var)
                dset.variables[var][:] = state[var].values
            # Rename the checkfile so the filename no longer specifies a 
            # single variable type
            newfile = checkfile+'_'+ob_type
            if last_ob == True:
                newfile += '.nc'

            os.system('mv {} {}'.format(checkfile,newfile))
            # ALL DONE!!
        else:
            
            # If the checkfile does not exist, make a new file
            #outfile = outdir_date_ens+'_'+ef.var_string(ob_type)

            ##tunit='seconds since 1970-01-01'
            #tunit='hours since 1900-01-01'
            ## Write ensemble forecast to netcdf
            #with Dataset(outfile,'w') as dset:
            #        dset.createDimension('time',None)
            #        dset.createDimension('lat',state.ny())
            #        dset.createDimension('lon',state.nx())
            #        dset.createDimension('ens',state.nmems())
            #        dset.createVariable('time','i4',('time',))
            #        dset.createVariable('lat',np.float32,('lat',))
            #        dset.createVariable('lon',np.float32,('lon'))
            #        dset.createVariable('ens','i4',('ens',))
            #        dset.variables['time'].units = tunit
            #        dset.variables['lat'].units = 'degrees_north'
            #        dset.variables['lon'].units = 'degrees_east'
            #        dset.variables['ens'].units = 'member_number'
            #        dset.variables['time'][:] = date2num(state.ensemble_times(),tunit)
            #        dset.variables['lat'][:] = state['lat'].values[:,0]
            #        dset.variables['lon'][:] = state['lon'].values[0,:]
            #        dset.variables['ens'][:] = state['mem'].values
            #        for var in state.vars():
            #            print('Writing variable {}'.format(var))
            #            dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
            #            dset.variables[var].units = ut.get_units(var)
            #            dset.variables[var][:] = state[var].values

#    
#    
#    #create observations
#    ob1 = Observation(value=5900, time=datetime(2017,9,6,12),lat=24.55,lon=278.21,
#                  obtype = 'Z500', localize_radius=2000, assimilate_this=True,
#                  error=10)