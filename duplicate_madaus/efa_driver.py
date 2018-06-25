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
import glob
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from EFA.efa_xray.observation.observation import Observation
from EFA.efa_xray.assimilation.ensrf import EnSRF
from EFA.duplicate_madaus.load_data import Load_Data
import EFA.duplicate_madaus.efa_functions as ef



ensemble_type = 'ecmwf'
#change this later
date = datetime(2013,4,1,0)
#end_date = datetime(2013,4,1,0)
#hourstep=12
#All variables in the netCDF
variables = ['T2M','ALT']#['T2M', 'ALT', 'P6HR', 'TCW']
#the ob type of the observations we are assimilating
obs_type = ['T2M','ALT']
#the variables in the netCDF we want to update
update_vars= ['T2M','ALT']
#are the observations only updating their corresponding variable, or
#are they updating all variables? -ie t2m only updates t2m, alt only updates alt
self_update=True #true if you want the above updates, otherwise false
#localization type
loc_type = 'GC'
#localization radius (for Gaspari-Cohn)
localize_radius = 1000


#if update_var == ob_type and len(ob_type)==1:
#    self_update = True #obs only update their corresponding vars
#    #if it's the last observation type to update, make this true to add a .nc
#    #to end of file name.
#    last_ob = False
#else:
#    self_update = False #each ob updates all ensemble vars

def run_efa(ob_type,update_var):
    
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

    observations = []   
    for o, o_type in enumerate(ob_type):
        efa = Load_Data(date,ensemble_type,variables,o_type,update_var)
        #only need to load the netCDF once (first time through the obtype loop)
        if o == 0:
            #initialize an instance of the Load_data class (load in data)
            statecls, lats, lons, elevs = efa.load_netcdfs()   
        
        #load in the obs file
        obs = efa.load_obs()
        
        #loop through each line in the text file (loop through each observation)
        for ob in obs:
            #this gets the observation information from the text file
            ob_dict = mt.get_ob_info(ob)
            #get longitude positive-definite- ie -130 lon is 230 E
            if ob_dict['lon'] < 0:
                ob_dict['lon'] = ob_dict['lon'] + 360
            utctime = datetime.utcfromtimestamp(ob_dict['time'])
            #check elevation of 4 nearest gridpoints to see whether or not to assimilate ob
            TorF = ef.closest_points(ob_dict['lat'],ob_dict['lon'],lats,
                                          lons,ob_dict['elev'],elevs)            
            #fill the observation class object with information for assimilation
            obser = Observation(value=ob_dict['ob'], time=utctime,
                            lat=ob_dict['lat'],lon=ob_dict['lon'], obtype=o_type, localize_radius=loc_rad,
                            assimilate_this=TorF, error=1)
            observations.append(obser)
    print('loaded '+str(len(observations))+' obs for assimilation')
            
#    ob1 = Observation(value=10000.25, time=datetime(2013,4,1,6),lat=24.55,lon=278.21,
#                  obtype = ob_type[0], localize_radius=1000, assimilate_this=True,
#                  error=1)
    
    
    #observations.append(ob1)
    #    
    # Put the state class object and observation objects into EnSRF object
    assimilator = EnSRF(statecls, observations, loc=loc_type)
    
    # Update the prior with EFA- post_state is an EnsembleState object
    post_state, post_obs = assimilator.update()
    
    state=post_state
    
    #build the string of where to save the files
    outdir = '/home/disk/hot/stangen/Documents/posterior_ensembles/'
    #are the observations only updating their corresponding variable, or
    #are they updating all variables? 
    if self_update == True:
        outdir += 'ob_update_self/' 
    elif self_update == False:
        outdir += 'ob_update_all/' 
    #directory for the type of localization used       
    outdir += 'loc_'+loc_str+'/'
    #if dealing with time-lagged ensembles, add another directory
    if time_lag == True:
        outdir += 'time_lag/'
    #add directories for the ens type and year & month
    outdir += ensemble_type+'/'+y+m+'/'
            
    #Create output directories if they don't yet exist
    if (os.path.isdir(outdir)):
        pass
    else:
        os.makedirs(outdir)
    
    outdir_date_ens = outdir+y+'goodbye-'+m+'-'+d+'_'+h+'_'+ensemble_type
            
    #if we are only updating the variable type that matches the observation type
    if self_update == True:
        for var in state.vars():
            checkfile = outdir_date_ens+'*'
            # If we're writing a single variable and the other variable already exists,
            # append to that file
            existing_file = glob.glob(checkfile)
            #if the other variable already exists
            if existing_file != []:
                #when only working with 2 variables, this should hold true
                last_ob = True
                print('Appending to existing file!')
                #should only be one file in there, convert from list to string
                existing_file = existing_file[0]
                # append to the existing file
                with Dataset(existing_file,'a') as dset:
                    print('Writing variable {}'.format(var))
                    dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
                    dset.variables[var].units = ut.get_units(var)
                    dset.variables[var][:] = state[var].values
                # Rename the checkfile so the filename no longer specifies a 
                # single variable type
                newfile = existing_file+'_'+ef.var_string(ob_type)
                if last_ob == True:
                    newfile += '.nc'
                    print('Done!')
    
                os.system('mv {} {}'.format(existing_file,newfile))
                # ALL DONE!!
            else:
                
                # If the checkfile does not exist, make a new file
                outfile = outdir_date_ens+'_'+ef.var_string(ob_type)
                ef.make_netcdf(state,outfile)
    
    #if we are updating all variable types with the observation (regardless of its type)
    elif self_update == False:
        outfile = outdir_date_ens+'_'+ef.var_string(ob_type)+'.nc'
        ef.make_netcdf(state,outfile)

#----------Call the run_efa function. Use a for loop if self-updating variables
#----------This is to update each variable one at a time.----------------------
if self_update == True:
    for obs in obs_type:
        #put back into list format or won't run in the function properly (due to ob_type loop)
        #make the ob_type and update_var match, so order of entry shouldn't matter
        run_efa([obs],[update_vars[update_vars.index(obs)]])
elif self_update == False:
    run_efa(obs_type,update_vars)
