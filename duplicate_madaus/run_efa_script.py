#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 12:10:56 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
import numpy as np
from datetime import datetime, timedelta
import EFA.efa_files.cfs_utilities_st as ut
import surface_obs.madis_example.madis_utilities as mt
import os
import glob
import sys
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from EFA.efa_xray.observation.observation import Observation
from EFA.efa_xray.assimilation.ensrf import EnSRF
from EFA.duplicate_madaus.load_data import Load_Data
import EFA.duplicate_madaus.efa_functions as ef


start_time = datetime.now()
print('start time: ',start_time)

shell_script = True

if shell_script == False:
    ensemble_type = 'eccc' #'eccc', 'ecmwf', or 'ncep'
    
    #All variables in the prior netCDF
    variables = ['IWV','IVT','D-IVT']#['TCW']#['QF850','D-QF850']#['T2M','ALT']#['T2M', 'ALT', 'P6HR', 'TCW']
    
    #the ob type of the observations we are assimilating
    obs_type = ['IVT']#['IVT','IWV']#['TCW']#['QF850']#['ALT','ALT']
    
    #observation error variance associated with the variables defined in obs_type.
    #list must be same length as obs_type. i.e. ob error variance for IVT is 10000,
    #ob error variance for IWV is 20.     
    ob_err_var = ['1000']#['10000','20'] #or ['ensvar'] for using ensemble variance as ob error variance
    
    #the variables in the netCDF we want to update- if self_update=True, this should match obs_type
    update_vars= ['IVT']#['IVT','IWV']#['IWV','IVT','D-IVT']#['TCW']#['QF850','D-QF850']#['ALT'] #['T2M','ALT']
    
    #is each observation type only updating its corresponding variable, or
    #is it updating all variables? -ie t2m only updates t2m, alt only updates alt
    #true if you want obs to only update their own variable type, otherwise false.
    self_update=True 
    
    #localization type- 'GC' is Gaspari-Cohn
    #'hybrid' uses a different localization radius if an
    #observation is within the AR being studied from 11/10/15 - 11/15/15
    #'statsig' and 'statsig2' use a Spearman rank correlation to determine which points
    #covary with the ob estimate at a confidence threshold, defined in localize_radius
    #'statsig' uses the covariances from the ensemble which is being updated as observations are assimilated,
    #'statsig2' uses the prior ensemble to determine statistical significance
    loc_type = 'GC'
    
    #localization halfwidth in km (for loc_type = 'GC') or, confidence threshold for
    #statistical significance ('statsig'/'statsig2') in percent- i.e. 99 = 99% confidence threshold.
    #for 'hybrid', the localization radius for obs within the AR, other obs use 1000 km.
    localize_radius = 10000 
    
    #date to run efa
    date = datetime(2015,11,13,12)#2013,4,1,0)
    
    #inflation?
    inflation = 'none' #scalar value is all I have set up to deal with. 
    
    #what kind of observations are we using? MADIS or gridded future 0-hour
    #forecast "obs", sampled at some interval?
    ob_category = 'gridded' #'madis' or 'gridded'
    
    #if true, include observation error variance in the names of the posterior variables
    #in the netCDF files, and in the names of the netCDF files.
    use_oberrvar = True   
    
    new_format = True #old or new naming conventions?
    
    efh = 54 #end forecast hour of each forecast, used in file name
    
    #grid of observations (if gridded. Shouldn't matter what is here if 
    #ob_category is 'madis', since all code will only do stuff
    #with grid if ob_category is 'gridded', and if new_format = True)
    #format is west edge, east edge, north edge, south edge, sampling interval
    grid = [-180,180,90,0,3] 
    

elif shell_script == True:
    ensemble_type = sys.argv[1]
    variables = sys.argv[2].split(',')
    obs_type = sys.argv[3].split(',')
    ob_err_var = sys.argv[4].split(',')
    update_vars = sys.argv[5].split(',')
    self_update = sys.argv[6]
    if self_update == 'true':
        self_update=True #true if you want the above updates, otherwise false
    elif self_update =='false':
        self_update=False
    loc_type = sys.argv[7] #loc_type = 'GC' 
    localize_radius = int(sys.argv[8])
    datestr = sys.argv[9]
    #first date that EFA is run in this batch
    firstdate = datetime.strptime(datestr,'%Y%m%d_%H%M')
    #get index of job number (make index start at 0)
    dateind = int(sys.argv[10]) - 1
    #use the index to access files corresponding with 12*ind hours after first 
    #forecast datetime
    date = firstdate+timedelta(hours=dateind*12)    
    #if a number is passed, convert to a float, otherwise, inflation = None
    inflation = sys.argv[11]
    ob_category = sys.argv[12]
    new_format = sys.argv[13]
    if new_format == 'true':
        new_format=True
    elif self_update == 'false':
        new_format=False
    efh = int(sys.argv[14])
    grid = sys.argv[15].split(',')
    use_oberrvar = sys.argv[16]
    if use_oberrvar == 'true':
        use_oberrvar=True
    elif use_oberrvar == 'false':
        use_oberrvar=False
    
#deal with the inflation based on its type, and generate a string for saving    
try:
    inflation = float(inflation)
    inflation_str = str(inflation)
    inflation_str = inflation_str.replace('.','-') 
except:
    if inflation == 'none':
        inflation = None
        inflation_str = 'none'
    

def run_efa(ob_type,update_var,ob_err_var):
    """
    All of the inputs are lists.
    """
    
    #print values used 
    print('localization radius: '+str(localize_radius))
    print('observation error variance: '+str(ob_err_var))
       
    #assign the radius to the localization_radius variable in Observation
    #create a string for saving to specific directories
    if loc_type == 'GC':
        localize_type = 'GC'
        loc_rad = localize_radius
        loc_str = str(loc_rad)
    elif loc_type == 'hybrid':
        localize_type = 'GC'
        #set the localization string to be the special loc radius, plus hybrid (2000hybrid)        
        loc_str = str(localize_radius)+loc_type
        
    elif loc_type.startswith('statsig'):
        localize_type = loc_type
        loc_rad = localize_radius
        loc_str = str(localize_radius)+loc_type

    else:
        loc_rad = None
        loc_str = '_stat_sig'
        
    y = date.strftime('%Y')
    m = date.strftime('%m')
    d = date.strftime('%d')
    h = date.strftime('%H')
    
    observations = []
    
    if loc_type == 'hybrid':
        #need to access 12-hour forecast IVT to determine presence of AR
        efaIVT = Load_Data(date,ensemble_type,variables,'IVT',['IVT'],
            grid=grid,new_format=new_format,efh=efh)
        
        #initialize an instance of the Load_data class (load in data)
        IVT_statecls, lats, lons, elevs = efaIVT.load_netcdfs()  
        #get IVT 12 hours into the forecast, at time of observations 
        #(time_ind = 2 is forecast hour 12)
        IVT = IVT_statecls.variables['IVT'].values[2,:,:,:]
        
        #obtain the mean of the ensemble (nlats x nlons)
        ens_mean = IVT.mean(axis=-1)
        
        #obtain variance of the ensemble (nlats x nlons)
        variance = np.var(IVT,axis=-1,ddof=1)
    
    #loop through each observation type    
    for o, o_type in enumerate(ob_type):
        
        #if we are wanting to use ensemble variance as ob error variance
        if ob_err_var[o].startswith('ensvar'):
            use_ens_var = True
        else:
            use_ens_var = False
        
        efa = Load_Data(date,ensemble_type,variables,o_type,update_var,
                        grid=grid,new_format=new_format,efh=efh)
        #only need to load the netCDF once (first time through the obtype loop)
        if o == 0:
            #initialize an instance of the Load_data class (load in data)
            statecls, lats, lons, elevs = efa.load_netcdfs()
        
        #if we are assimilating MADIS observations 
        if ob_category == 'madis':
            #load in the obs file
            obs = efa.load_obs()
        #if we are assimilating next-cycle 0-hour forecast gridded "observations"
        elif ob_category == 'gridded':
            #load in the obs file
            obs = efa.load_obs(forecast_hour=12,madis=False,variance=use_ens_var)
            
        #loop through each line in the text file (loop through each observation)
        for ob in obs:
            #this gets the observation information from the text file
            ob_dict = mt.get_ob_info(ob,use_ens_var)
            #get longitude positive-definite- ie -130 lon is 230 E
            if ob_dict['lon'] < 0:
                ob_dict['lon'] = ob_dict['lon'] + 360
            #get ob time in datetime format
            utctime = datetime.utcfromtimestamp(ob_dict['time'])
            if ob_category == 'madis':
                #check elevation of 4 nearest gridpoints to see whether or not to assimilate ob-
                #returns true or false
                TorF = ef.closest_points(ob_dict['lat'],ob_dict['lon'],lats,
                                              lons,ob_dict['elev'],elevs)
            elif ob_category == 'gridded':
                #no need to check elevation, since using gridded obs
                TorF = True
                
            #if we are using ens variance as ob error variance
            if use_ens_var == True:
                #multiply by factor specified at end of 'ensvar'
                mult_factor = ob_err_var[o].replace('ensvar','')
                if mult_factor == '':
                    mult_factor = 1
                else:
                    mult_factor = float(mult_factor)
                print('multiplication factor: ',mult_factor)
                ob_var = ob_dict['variance']*mult_factor
            elif use_ens_var == False:
                ob_var = float(ob_err_var[o])
                
            #if we are using a hybrid localization radius- longer localization 
            #radius within the AR, 1000 km outside of it
            #if the lat/lon lie within some AR box, check if the ob point is in an AR-
            #if IVT > 250
            if loc_type == 'hybrid':    

                if ob_dict['lon'] >= 140 and ob_dict['lon'] <= 245 and ob_dict['lat'] >= 33 and ob_dict['lat'] <= 51:
                    
                    #get the ensemble value at the lat/lon pair (ob_variance isn't used)
                    ob_value, ob_variance = ef.closest_points(ob_dict['lat'],ob_dict['lon'],lats,lons,variable=ens_mean,
                                                 need_interp=True,gen_obs=True,variance=variance)
                    #if ob value is AR and in the grid area, set its localization radius to the input (hybrid) loc_rad
                    if ob_value >= 250:
                        loc_rad = localize_radius
                        #print(str(ob_value))
                    else:
                        loc_rad = 1000
                else:    
                    #set the default loc_rad to be 1000
                    loc_rad = 1000
            
            #check if it's working
#            print(str(ob_dict['ob']))
#            print(str(ob_dict['lat'])+' '+str(ob_dict['lon'])+' '+str(loc_rad))
            #fill the observation class object with information for assimilation
            obser = Observation(value=ob_dict['ob'], time=utctime,
                            lat=ob_dict['lat'],lon=ob_dict['lon'], obtype=o_type, localize_radius=loc_rad,
                            assimilate_this=TorF, error=ob_var)
            observations.append(obser)                

    print('loaded '+str(len(observations))+' obs for assimilation')
               
    # Put the state class object and observation objects into EnSRF object
    assimilator = EnSRF(statecls, observations, inflation=inflation, loc=localize_type)
    
    # Update the prior with EFA- state is an EnsembleState object and is the posterior, post_obs isn't used
    state, post_obs = assimilator.update()
        
#---build the string of which directory to save the file--------------------------------
    outdir = '/home/disk/hot/stangen/Documents/posterior_ensembles/'
    #what kind of observations are we assimilating?
    outdir += ob_category + '/'
    
    #are the observations only updating their corresponding variable, or
    #are they updating all variables? 
    if self_update == True:
        outdir += 'ob_update_self/' 
    elif self_update == False:
        outdir += 'ob_update_all/' 
    #directory for inflation used
    outdir += 'inf_'+inflation_str+'/'
    #directory for the type of localization used       
    outdir += 'loc_'+loc_str+'/'
    #add directories for the ens type and year & month
    outdir += ensemble_type+'/'+y+m+'/'
#----------------------------------------------------------------------------------------
            
    #Create output directories if they don't yet exist
    if (os.path.isdir(outdir)):
        pass
    else:
        os.makedirs(outdir)
    
    #add initialization date to filename
    outdir_date_ens = outdir+y+'-'+m+'-'+d+'_'+h
    
    #new format string: add nhrs and if gridded, grid dimensions to filename
    if new_format == True:
        outdir_date_ens += '_'+str(efh)+'hrs'
        if ob_category == 'gridded':
            outdir_date_ens += '_'+ef.var_string(grid)
            
    #convert 1-length ob err var list to a string with no decimals
    #this makes it so that multiple ob err vars that acted on one data type
    #can all be saved to one netCDF- makes managing files easier where I was testing
    #what observation error variance to use.
    #more specifically, this adds ob err var to the variable name in the netCDF.
    #ob_err_var_str is only called if self_update == True, so the ob err var
    #is never in the variable name within the netCDF if self_update == False.
    ob_err_var_str = ''
    if use_oberrvar == True:
        ob_err_var_str = str(ob_err_var[0]).replace('.','-')
    #if we don't want the observation error variance used in the name of the file
    #and/or the names of the variables in the file
    elif use_oberrvar == False:
        ob_err_var = ''
            
    #if we are only updating the variable type that matches the observation type
    if self_update == True:
        for var in state.vars():
            checkfile = outdir_date_ens+'*'
            # If we're self-updating variables and we are on the 2nd or later variable in
            # the loop, this will return with something
            existing_file = glob.glob(checkfile)

            #if other variables already exists, append to the netCDF
            if existing_file != []:
                print('Appending to existing file!')
                #should only be one file in there, convert from list to string
                existing_file = existing_file[0]
                
                # append to the existing file
                with Dataset(existing_file,'a') as dset:
                    print('Writing variable {}'.format(var))
                    dset.createVariable(var+ob_err_var_str, np.float32, ('time','lat','lon','ens',))
                    dset.variables[var+ob_err_var_str].units = ef.get_units(var)
                    dset.variables[var+ob_err_var_str][:] = state[var].values
                #if the filename already contains a .nc because we created the
                #file we are appending to separately from this run, delete
                #the .nc
                existing_file = existing_file.replace('.nc', '')
                # Rename the checkfile so the filename no longer specifies a 
                # single variable type- add new variable to filename
                newfile = existing_file+'_'+ef.var_num_string(ob_type,ob_err_var)+'.nc'
                os.system('mv {} {}'.format(existing_file+'.nc',newfile))
            else:                
                # If the checkfile does not exist, make a new file
                outfile = outdir_date_ens+'_'+ef.var_num_string(ob_type,ob_err_var)+'.nc'
                ef.make_netcdf(state,outfile,ob_err_var_str)
    
    #if we are updating all variable types with the observation (regardless of its type)
    #and use_oberrvar == True, observation error variance is not saved in the 
    #variable name in the netCDF, but it is saved in the filename.
    elif self_update == False:
        outfile = outdir_date_ens+'_'+ef.var_num_string(ob_type,ob_err_var)+'.nc'
        ef.make_netcdf(state,outfile)    

#----------Call the run_efa function. Use a for loop if self-updating variables
#----------This is to update each variable one at a time.----------------------
if self_update == True:
    for i, obs in enumerate(obs_type):
        #put back into list format or won't run in the function properly (due to ob_type loop)
        #make the ob_type and update_var match, so order of entry shouldn't matter
        #this runs EFA on one ob_type using its associated ob error var
        #at a time, and only updates the variable matching the ob type. The updates 
        #will be appended together into one netCDF if updating more than 1 variable.
        run_efa([obs],[update_vars[update_vars.index(obs)]],[ob_err_var[i]])
#this will update all variable types specified using all observations, of any type.
elif self_update == False:
    run_efa(obs_type,update_vars,ob_err_var)

print('Done!')
end_time = datetime.now()
print('end time: ',end_time)
print('total time elapsed: ',end_time-start_time)