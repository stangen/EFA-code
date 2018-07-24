#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 11:48:09 2018

@author: stangen
"""
import numpy as np
from datetime import datetime
import surface_obs.madis_example.madis_utilities as mt
import sys
from EFA.duplicate_madaus.load_data import Load_Data
import EFA.duplicate_madaus.efa_functions as ef

start_time = datetime.now()
print('start time: ',start_time)
#to run this in spyder, change this to False, to run in shell script, change to True
shell_script=True

if shell_script==False:
    ensemble_type = 'ecmwf'
    #used for loading the prior files-order matters for loading the file
    #used for prior ens stats comparison and for loading grid for use as observation
    #all the vars-look at the filename if unsure what they all are.
    prior_vrbls = ['QF850','D-QF850']#['T2M','ALT']  #['T2M','ALT']
    #used for loading the posterior files-put in ob err var as part of the string, order matters
    #all the vars-look at the filename if unsure what they all are.
    #if doing stats on prior, don't worry about these, they won't be loaded.
    post_vrbls = ['QF8501']#['ALT1','ALT0','ALT0-1']#['T2M','ALT']#
    #which ob type we're getting stats for, can be more than 1 ['T2M','ALT']
    ob_types = ['QF850']#['ALT']
    #observation error variance associated with the observation type
    #If prior, don't worry about value, will be set to [''] later
    ob_err_var = ['1']#[0] 
    #all obs we will have in the end- for saving the file
    allobs = ['QF850']#['ALT']
    #date range of ensembles used
    start_date = datetime(2015,11,12,0)#2013,4,4,0)
    end_date = datetime(2015,11,12,12)#2013,4,4,0)
    #hour increment between ensemble forecasts
    incr = 12
    
    #change for how far into the forecast to get observations, 1 is 6 hours in,
    #end_index is the last forecast hour to get obs for
    start_index = 3
    end_index = 4
    
    #if True, will load/do stats on posterior ensembles, if False, will load prior
    post=False
     
    #localization radius
    loc_rad = '1000'
    
    #inflation factor (string for loading files)
    inflation = 'none'
    
    #what kind of observations did we assimilate? MADIS or gridded future 0-hour
    #forecast "obs", sampled at some interval?
    ob_category = 'gridded'
    
    new_format = True #old or new naming conventions?
    efh = 54 #end forecast hour of each forecast
    grid = [-180,180,90,0,3] #grid of observations (if gridded)
    #if true, include observation error variance in the names of the posterior variables
    #in the netCDF files, and in the names of the netCDF files, to load them.
    use_oberrvar = True
    #did we self-update each variable
    self_update = False
    
    #create strings for saving file at the end
    sy = start_date.strftime('%Y')
    sm = start_date.strftime('%m')
    sd = start_date.strftime('%d')
    sh = start_date.strftime('%H')
    
    ey = end_date.strftime('%Y')
    em = end_date.strftime('%m')
    ed = end_date.strftime('%d')
    eh = end_date.strftime('%H')
    
    datestr=sy+sm+sd+sh+'-'+ey+em+ed+eh
    
elif shell_script==True:
    ensemble_type = sys.argv[1]
    prior_vrbls = sys.argv[2].split(',')
    post_vrbls = sys.argv[3].split(',')
    ob_types = sys.argv[4].split(',')
    ob_err_var = sys.argv[5].split(',')
    allobs = sys.argv[6].split(',')
    startstr = sys.argv[7]
    start_date = datetime.strptime(startstr,'%Y%m%d%H')
    endstr = sys.argv[8]
    end_date = datetime.strptime(endstr,'%Y%m%d%H')
    incr=int(sys.argv[9])
    start_index = int(sys.argv[10])
    end_index=int(sys.argv[11])
    boolstr=sys.argv[12]
    #change string to boolean for loading prior or posterior ensembles
    if boolstr == 'true':
        post=True
    elif boolstr =='false':
        post=False     
    loc_rad = sys.argv[13]  
    inflation = sys.argv[14]
    ob_category = sys.argv[15]
    new_format = sys.argv[16]
    if new_format == 'true':
        new_format=True
    elif self_update == 'false':
        new_format=False
    efh = int(sys.argv[17])
    grid = sys.argv[18].split(',')
    use_oberrvar = sys.argv[19]
    if use_oberrvar == 'true':
        use_oberrvar=True
    elif use_oberrvar == 'false':
        use_oberrvar=False
    #are the observations only updating their corresponding variable, or
    #are they updating all variables? -ie t2m only updates t2m, alt only updates alt
    self_update = sys.argv[20]
    if self_update == 'true':
        self_update=True #true if you want the above updates, otherwise false
    elif self_update =='false':
        self_update=False
    
    datestr = startstr+'-'+endstr



save_dir = '/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/mse_var_output/'    
#variable string (for saving the .txt file)
varstr = ef.var_string(allobs)
#make a list of netCDF variable names/strings from obs+obs error var
#post_vrbls contains the variable names in the filename in the posterior
netcdf_varnames = []
#make a list of variable names we want to save to the .txt file (this may include ob err variance,
#even when not doing self-update)
dict_varnames = []

for i, ob in enumerate(ob_types):
    #If we are loading prior to do stats on, we need to set ob error variance to 
    #[''] to make variable in netCDF load correctly.
    #also, if we did not previously use the obs err var in creating the posterior
    #ens name/variable names, or did not only self-update each variable,
    #set ob_err_var to [''], since ob err var was not included in names.    
    if post == False or use_oberrvar == False or self_update == False:
        netcdf_varnames.append(ob)
        #ob_err_var = ['']*len(ob_types)
    else:
        netcdf_varnames.append(ef.var_num_string([ob],[ob_err_var[i]]))
    #this defines the variable names in the dictionary/.txt file we are creating.
    if post == False or use_oberrvar == False:
        dict_varnames.append(ob)
    else:
        dict_varnames.append(ef.var_num_string([ob],[ob_err_var[i]]))
#prior/post string
if post==True:
    prior_or_post='loc'+loc_rad
if post==False:
    prior_or_post='prior'
    
#ob_update_self or ob_update_all
if self_update == True:
    ob_upd = 'ob_update_self'
elif self_update == False:
    ob_upd = 'ob_update_all'
    
#last forecast hour I want to get observations for
end_hour = 6*end_index

#the dict where all the data is
ob_dict = {}
#dict for the means from the data
data_dict = {}
#dicts to store SE and variance from all obs for the time period, not just stationary ones.
SE_dict = {}
variance_dict = {}
#dicts to store SE and variance from gridded obs, and from comparison with every gridpoint.
gridob_SE_dict = {}
gridob_variance_dict = {}

#creates a datelist which increments in 12 or 24 hour chunks
dates = mt.make_datetimelist(start_date, end_date, incr)
for date in dates:
    hour = date.strftime('%H')

    #create dict for variable type
    #ob_type is purely used in defining the variable name in the dictionary,
    #which in turn is used in saving the .txt file. It is NOT used as a parameter
    #in any of the functions or as a netCDF variable name. These roles are
    #filled by ob_types[i], netcdf_varnames[i], and post_vrbls
    for i, ob_type in enumerate(dict_varnames):
        ob_dict[ob_type] = ob_dict.get(ob_type,{})
        data_dict[ob_type] = data_dict.get(ob_type,{})
        SE_dict[ob_type] = SE_dict.get(ob_type,{})
        variance_dict[ob_type] = variance_dict.get(ob_type,{})
        
        gridob_SE_dict[ob_type] = gridob_SE_dict.get(ob_type,{})
        gridob_variance_dict[ob_type] = gridob_variance_dict.get(ob_type,{})
        
        #Load in the ensemble data for a given initialization
        efa = Load_Data(date,ensemble_type,prior_vrbls,ob_types[i],[netcdf_varnames[i]],post_vrbls=post_vrbls,
                        grid=grid,new_format=new_format,efh=efh)
        statecls, lats, lons, elevs = efa.load_netcdfs(post=post,ob_cat=ob_category,inf=inflation,lr=loc_rad,ob_upd=ob_upd)
        
        #Want to interpolate ALL observations valid during a given ensemble forecast.
        #Or not, to minimize runtime of script. 1 ens, var, and forecast hour takes
        #about 1 minute to run. 
        hour_step = 6
        j = start_index
        #fh = hour_step*j
        #sfh = str(fh)
        # j = the current forecast hour index. i.e. 1 = 6 hr, 2 = 12 hr, 3 = 18 hr, etc.
        while j <= end_index:
        #while fh <= end_hour:
            fh = j*hour_step
            sfh = str(fh)
            #print(fh)
            
            #only compare against MADIS obs if using T2M or ALT 
            if ob_types[i] == 'T2M' or ob_types[i] == 'ALT':
                obs = efa.load_obs(fh)
                    
                #create dictionary for forecast hour
                ob_dict[ob_type][sfh] = ob_dict[ob_type].get(sfh,{})
                data_dict[ob_type][sfh] = data_dict[ob_type].get(sfh,{})
                SE_dict[ob_type][sfh] = SE_dict[ob_type].get(sfh,[])            
                variance_dict[ob_type][sfh] = variance_dict[ob_type].get(sfh,[])
                
                #create dictionary for 00Z or 12Z
                ob_dict[ob_type][sfh][hour] = ob_dict[ob_type][sfh].get(hour,{})
                
                #Create dictionaries to hold the mean values to be calculated
                data_dict[ob_type][sfh]['station_all_mse'] = data_dict[ob_type][sfh].get('station_all_mse', np.array([]))
                data_dict[ob_type][sfh]['station_all_mean_variance'] = data_dict[ob_type][sfh].get('station_all_mean_variance', np.array([]))
                data_dict[ob_type][sfh]['station_all_error_variance'] = data_dict[ob_type][sfh].get('station_all_error_variance', np.array([]))
                data_dict[ob_type][sfh]['weight'] = data_dict[ob_type][sfh].get('weight',np.array([]))
                
                #can probably delete this
                obs_pass = []
                obs_allmaritime = []
                for ob_counter, ob in enumerate(obs):
                    #this gets the observation information from the text file
                    ob_info = mt.get_ob_info(ob)
                    #get longitude positive-definite- ie -130 lon is 230 E
                    if ob_info['lon'] < 0:
                        ob_info['lon'] = ob_info['lon'] + 360
                    utctime = datetime.utcfromtimestamp(ob_info['time'])
                    #find interpolated ob estimate, if it passes the terrain check.
                    #the terrain check is done within the closest_points function
                    interp, TorF = ef.closest_points(ob_info['lat'],ob_info['lon'],lats,lons,ob_info['elev'],
                                               elevs,utctime,statecls['validtime'].values,
                                               statecls.variables[netcdf_varnames[i]].values,need_interp=True)
    
                    ob_id = ob_info['name']
                    ob_value = ob_info['ob']
                    
                    #for stations which pass the terrain check (and were assimilated):
                    if TorF == True:
                    #if len(interp) > 0:
                        #calculate MSE
                        se = (np.mean(interp)-ob_value)**2
                        #calculate variance
                        hx_variance_unbiased = np.var(interp, ddof=1)
    
                        #only fill the dictionaries with stationary stations, so bias removal can work.
                        if ob_info['stationary'] == 0:
                        
                            #-----Added just to save obs used, to plot later-----------
                            obs_pass.append(ob)
                            
                            #create dictionary for station ID if it doesn't exist
                            ob_dict[ob_type][sfh][hour][ob_id] = ob_dict[ob_type][sfh][hour].get(ob_id,{})
                            ob_dict[ob_type][sfh][hour][ob_id]['variance_unbiased'] = ob_dict[ob_type][sfh][hour][ob_id].get('variance_unbiased',[])
                            ob_dict[ob_type][sfh][hour][ob_id]['se'] = ob_dict[ob_type][sfh][hour][ob_id].get('se',[])
                            ob_dict[ob_type][sfh][hour][ob_id]['error'] = ob_dict[ob_type][sfh][hour][ob_id].get('error',[])
                            
                            #add the data to the dictionaries
                            ob_dict[ob_type][sfh][hour][ob_id]['variance_unbiased'].append(hx_variance_unbiased)                
                            ob_dict[ob_type][sfh][hour][ob_id]['error'].append(np.mean(interp)-ob_value)
                            ob_dict[ob_type][sfh][hour][ob_id]['se'].append(se)
                            
                            #ob_counter +=1
                            #print("on observation "+str(ob_counter)+" out of "+str(len(obs)))  
                        
                        #add squared error and ensemble variance from all observations which pass terrain check
                        SE_dict[ob_type][sfh].append(se)
                        variance_dict[ob_type][sfh].append(hx_variance_unbiased)
                        obs_allmaritime.append(ob)
                                      
                print('Added stats from '+str(len(obs_pass))+' obs to stats dict')
                print('Added stats from '+str(len(obs_allmaritime))+'- including all maritime observations')
            #-------End of MADIS obs load------------------------------------------------------------------------
            
            #-------Start of gridded obs load--------------------------------------------------------------------
            #gridded observations are only available every 12 hours
            if ob_category =='gridded' and j % 2 == 0:
                #load in the gridded observations for comparison
                grid_obs = efa.load_obs(fh,madis=False)
                
                #create dictionary for forecast hour                
                gridob_SE_dict[ob_type][sfh] = gridob_SE_dict[ob_type].get(sfh,{})
                gridob_variance_dict[ob_type][sfh] = gridob_variance_dict[ob_type].get(sfh,{})
                
                #create empty lists to append SE and variance into
                gridob_SE_dict[ob_type][sfh]['obs'] = gridob_SE_dict[ob_type][sfh].get('obs',[])
                gridob_SE_dict[ob_type][sfh]['all_gridpoints'] = gridob_SE_dict[ob_type][sfh].get('all_gridpoints',[])
                gridob_SE_dict[ob_type][sfh]['gridpoints_region'] = gridob_SE_dict[ob_type][sfh].get('gridpoints_region',[])
                gridob_SE_dict[ob_type][sfh]['gridpoints_region2'] = gridob_SE_dict[ob_type][sfh].get('gridpoints_region2',[])
                gridob_variance_dict[ob_type][sfh]['obs'] = gridob_variance_dict[ob_type][sfh].get('obs',[])
                gridob_variance_dict[ob_type][sfh]['all_gridpoints'] = gridob_variance_dict[ob_type][sfh].get('all_gridpoints',[])
                gridob_variance_dict[ob_type][sfh]['gridpoints_region'] = gridob_variance_dict[ob_type][sfh].get('gridpoints_region',[])
                gridob_variance_dict[ob_type][sfh]['gridpoints_region2'] = gridob_variance_dict[ob_type][sfh].get('gridpoints_region2',[])
               
                for ob_counter, ob in enumerate(grid_obs):
                    #this gets the observation information from the text file
                    ob_info = mt.get_ob_info(ob)
                    #get longitude positive-definite- ie -130 lon is 230 E
                    if ob_info['lon'] < 0:
                        ob_info['lon'] = ob_info['lon'] + 360
                    utctime = datetime.utcfromtimestamp(ob_info['time'])
                    #find interpolated ob estimate, if it passes the terrain check.
                    #the terrain check is done within the closest_points function
                    interp, TorF = ef.closest_points(ob_info['lat'],ob_info['lon'],lats,lons,ob_info['elev'],
                                               elevs,utctime,statecls['validtime'].values,
                                               statecls.variables[netcdf_varnames[i]].values,need_interp=True)
                    
                    ob_value = ob_info['ob']
                    #calculate MSE
                    
                    se = (np.mean(interp)-ob_value)**2
                    #calculate variance
                    hx_variance_unbiased = np.var(interp, ddof=1)
                    #append SE/variance to list
                    gridob_SE_dict[ob_type][sfh]['obs'].append(se)
                    gridob_variance_dict[ob_type][sfh]['obs'].append(hx_variance_unbiased)
                    
                print('Added stats from '+str(len(grid_obs))+' gridded obs')        
                
                #comparison with all grid points
                #access the 0 hour forecast initialized fh hours after EFA'd forecast (the analysis/observation grid)
                anl_ens_mean, lats, lons = efa.load_ens_netcdf(fh)
                #find the ensemble mean at each grid point for the specified forecast hour (the updated forecast grid)
                fcst_ens_mean = (statecls.variables[netcdf_varnames[i]].values[j]).mean(axis=-1)
                #find the squared error of each grid point                
                se = (fcst_ens_mean-anl_ens_mean)**2
                #find variance of each forecast grid point for the specified forecast hour
                fcst_var = np.var(statecls.variables[netcdf_varnames[i]].values[j], axis=-1, ddof=1)
                
                #find average of northern hemisphere- make it so that high
                #latitudes aren't weighted unfairly high
                weights = np.cos(np.radians(lats))
                #np.average is the weighted average of all latitudes for each longitude, 
                #np.mean is the average of all the longitudes
                avg_se = np.mean(np.average(se, axis=0, weights=weights))
                avg_var = np.mean(np.average(fcst_var, axis=0, weights=weights))
                #now append each date's ensemble-average values to the corresponding list
                gridob_SE_dict[ob_type][sfh]['all_gridpoints'].append(avg_se)
                gridob_variance_dict[ob_type][sfh]['all_gridpoints'].append(avg_var)
                
                #comparison within a defined region of gridpoints- along west coast
                w = -135
                e = -115
                n = 55
                s = 30
                
                #convert to grid indices, add 1 to the right endpoints so python includes the last index
                l = int(abs(-180-w)*2)
                r = int(abs(-180-e)*2)+1
                t = int((90-n)*2)
                b = int((90-s)*2)+1
                
                se_region = se[t:b,l:r]
                var_region = fcst_var[t:b,l:r]
                weights_region = np.cos(np.radians(lats[t:b]))
                avg_se_region = np.mean(np.average(se_region, axis=0, weights=weights_region))
                avg_var_region = np.mean(np.average(var_region, axis=0, weights=weights_region))
                gridob_SE_dict[ob_type][sfh]['gridpoints_region'].append(avg_se_region)
                gridob_variance_dict[ob_type][sfh]['gridpoints_region'].append(avg_var_region)
                
                #comparison within a defined region specific to this AR
                w2 = -180
                e2 = -115
                n2 = 50
                s2 = 35
                
                #convert to grid indices, add 1 to the right endpoints so python includes the last index
                l2 = int(abs(-180-w2)*2)
                r2 = int(abs(-180-e2)*2)+1
                t2 = int((90-n2)*2)
                b2 = int((90-s2)*2)+1
                
                se_region2 = se[t2:b2,l2:r2]
                var_region2 = fcst_var[t2:b2,l2:r2]
                weights_region2 = np.cos(np.radians(lats[t2:b2]))
                avg_se_region2 = np.mean(np.average(se_region2, axis=0, weights=weights_region2))
                avg_var_region2 = np.mean(np.average(var_region2, axis=0, weights=weights_region2))
                gridob_SE_dict[ob_type][sfh]['gridpoints_region2'].append(avg_se_region2)
                gridob_variance_dict[ob_type][sfh]['gridpoints_region2'].append(avg_var_region2)
                        
            #update the forecast hour to load the next forecast, 6 hours later
            j = j + 1
            fh = hour_step*j
            #print(fh)
            sfh = str(fh)

#now go through all the data, convert to numpy arrays, and calculate the bias, 
#mse, and mean variance for each station and forecast hour

#create list to save for creating plots
stats_list = []

#if we can compare with MADIS obs
if ob_types[i] == 'T2M' or ob_types[i] == 'ALT': 
    for var in ob_dict:
        for f_h in ob_dict[var]:
            for h in ob_dict[var][f_h]:
                for sID in list(ob_dict[var][f_h][h].keys()):
                    ob_dict[var][f_h][h][sID]['variance_unbiased'] = np.array(ob_dict[var][f_h][h][sID]['variance_unbiased'])
                    ob_dict[var][f_h][h][sID]['se'] = np.array(ob_dict[var][f_h][h][sID]['se'])
                    ob_dict[var][f_h][h][sID]['error'] = np.array(ob_dict[var][f_h][h][sID]['error'])
                    
                    #Greg's method of finding variance of error being equal to mse - bias^2
                    error_var = np.var(ob_dict[var][f_h][h][sID]['error'])
                    #calculate mse
                    ob_dict[var][f_h][h][sID]['mse'] = np.mean(ob_dict[var][f_h][h][sID]['se'])
                    #calculate mean variance (using ddof=1)
                    ob_dict[var][f_h][h][sID]['mean_variance_unbiased'] = np.mean(ob_dict[var][f_h][h][sID]['variance_unbiased'])
    
                    ob_dict[var][f_h][h][sID]['error_variance'] = error_var
                                 
    
                    data_dict[var][f_h]['station_all_mse'] = np.append(data_dict[var][f_h]['station_all_mse'],ob_dict[var][f_h][h][sID]['mse'])
                    data_dict[var][f_h]['station_all_mean_variance'] = np.append(data_dict[var][f_h]['station_all_mean_variance'],ob_dict[var][f_h][h][sID]['mean_variance_unbiased'])
                    data_dict[var][f_h]['station_all_error_variance'] = np.append(data_dict[var][f_h]['station_all_error_variance'],ob_dict[var][f_h][h][sID]['error_variance'])
                    data_dict[var][f_h]['weight'] = np.append(data_dict[var][f_h]['weight'],len(ob_dict[var][f_h][h][sID]['se']))
                    del ob_dict[var][f_h][h][sID]
                    
            data_dict[var][f_h]['average_mse'] = np.average(data_dict[var][f_h]['station_all_mse'],weights=data_dict[var][f_h]['weight'])
            data_dict[var][f_h]['average_variance'] = np.average(data_dict[var][f_h]['station_all_mean_variance'],weights=data_dict[var][f_h]['weight'])
            data_dict[var][f_h]['average_error_variance'] = np.average(data_dict[var][f_h]['station_all_error_variance'],weights=data_dict[var][f_h]['weight'])
     
            #find MSE and variance of all observations, not just stationary obs, of time period.
            SE_list = np.array(SE_dict[var][f_h])
            #SE_list_sorted = np.sort(SE_list)
            MSE = np.mean(SE_list)
            variance_list = np.array(variance_dict[var][f_h])
            #variance_list_sorted = np.sort(variance_list)
            variance = np.mean(variance_list)
            
    #        SE_list_stationary = data_dict[var][f_h]['station_all_mse']
    #        SE_list_stationary_sorted = np.sort(SE_list_stationary)
            
            #gridded observations are only available every 12 hours
            #all of the forecast hours will go into the ob_dict/data_dicts. However,
            #only fh 12, 24, 36, etc will go into the grid dicts. This is why the try
            #block works- for these forecast hours, the grid dicts will have info,
            #but for fh 6, 18, 30, etc they will not, but it will loop through each 
            #forecast hour that is in ob_dict and data_dict, which may contain these
            #forecast hours. Therefore, the try block will execute for fh 12, 24, 36,
            #etc and the except block will execute for 6, 18, 30, etc, and the stats,
            #whether the full comparison with gridded obs or just stationary obs,
            #will be appended and saved in the .txt file.
            try:
                #MSE and variance using ensemble mean points as observations, and MSE and
                #variance of the entire grid
                MSE_gridob_obs_list = np.array(gridob_SE_dict[var][f_h]['obs'])
                MSE_gridob_obs = np.mean(MSE_gridob_obs_list)
                MSE_gridob_fullgrid_list = np.array(gridob_SE_dict[var][f_h]['all_gridpoints'])
                MSE_gridob_fullgrid = np.mean(MSE_gridob_fullgrid_list)
                
                variance_gridob_obs_list = np.array(gridob_variance_dict[var][f_h]['obs'])
                variance_gridob_obs = np.mean(variance_gridob_obs_list)
                variance_gridob_fullgrid_list = np.array(gridob_variance_dict[var][f_h]['all_gridpoints'])
                variance_gridob_fullgrid = np.mean(variance_gridob_fullgrid_list)
                
                #append to the stats list each statistic and info for one inflation, localization radius, ens type, variable+ob_error_var, and forecast hour.
                #order of stats is comparison with: stationary MADIS obs MSE/variance, all MADIS obs MSE/variance,
                #MSE/var at gridded observation locations, MSE/var for entire grid (averaged proportional to surface area of NH.)
                stats_list.append(inflation+','+prior_or_post+','+ensemble_type+','+var+','+f_h+','+str(data_dict[var][f_h]['average_mse'])+','+str(data_dict[var][f_h]['average_variance'])+','+str(MSE)+','+str(variance)+','+str(MSE_gridob_obs)+','+str(variance_gridob_obs)+','+str(MSE_gridob_fullgrid)+','+str(variance_gridob_fullgrid)+'\n')
            except:            
                #append to the stats list each statistic and info for one inflation, localization radius, ens type, variable+ob_error_var, and forecast hour.
                #Don't try to append the gridded obs stuff, since there are no observations to compare with at forecast hours 6, 18, 30, 42, etc.
                stats_list.append(inflation+','+prior_or_post+','+ensemble_type+','+var+','+f_h+','+str(data_dict[var][f_h]['average_mse'])+','+str(data_dict[var][f_h]['average_variance'])+','+str(MSE)+','+str(variance)+'\n')

#even if the dictionaries are empty (don't contain forecast hours), this will
#still work. This is because it just won't execute the block of code, since
#there are no f_h in var
else:
    for var in gridob_SE_dict:
        for f_h in gridob_SE_dict[var]:
            MSE_gridob_obs_list = np.array(gridob_SE_dict[var][f_h]['obs'])
            MSE_gridob_obs = np.mean(MSE_gridob_obs_list)
            MSE_gridob_fullgrid_list = np.array(gridob_SE_dict[var][f_h]['all_gridpoints'])
            MSE_gridob_fullgrid = np.mean(MSE_gridob_fullgrid_list)
            MSE_gridob_region_list = np.array(gridob_SE_dict[var][f_h]['gridpoints_region'])
            MSE_gridob_region = np.mean(MSE_gridob_region_list)
            MSE_gridob_region2_list = np.array(gridob_SE_dict[var][f_h]['gridpoints_region2'])
            MSE_gridob_region2 = np.mean(MSE_gridob_region2_list)
            
            variance_gridob_obs_list = np.array(gridob_variance_dict[var][f_h]['obs'])
            variance_gridob_obs = np.mean(variance_gridob_obs_list)
            variance_gridob_fullgrid_list = np.array(gridob_variance_dict[var][f_h]['all_gridpoints'])
            variance_gridob_fullgrid = np.mean(variance_gridob_fullgrid_list)
            variance_gridob_region_list = np.array(gridob_variance_dict[var][f_h]['gridpoints_region'])
            variance_gridob_region = np.mean(variance_gridob_region_list)
            variance_gridob_region2_list = np.array(gridob_variance_dict[var][f_h]['gridpoints_region2'])
            variance_gridob_region2 = np.mean(variance_gridob_region2_list)            
            
            stats_list.append(inflation+','+prior_or_post+','+ensemble_type+','+var+','+f_h+','+str(MSE_gridob_obs)+','+str(variance_gridob_obs)+','+str(MSE_gridob_fullgrid)+','+str(variance_gridob_fullgrid)+','+str(MSE_gridob_region)+','+str(variance_gridob_region)+','+str(MSE_gridob_region2)+','+str(variance_gridob_region2)+'\n')
        
#save the stats list after all forecast hours have been appended for one variable
print('Writing statistics to .txt file')
savedir_str = save_dir+datestr+'_'+varstr
if new_format == True:
    savedir_str += '_' + ef.var_string(grid)
if ob_category == 'gridded':
    savedir_str += '_gridobs'
savedir_str += '.txt'
f = open(savedir_str, 'a')
for s in stats_list:
    f.write(s)
f.close()

end_time = datetime.now()
print('end time: ',end_time)
print('total time elapsed: ',end_time-start_time)
print('Done!')

#want to save ens type, ob, forecast hour in string identifier, followed by mse, variance, and error variance

##this was to save all the stations that pass the elevation check so I could plot them. 
#f = open('/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/plots/2013040106obs_allmaritime.txt','w')
#for obser in obs_allmaritime:
##for obser in obs_pass:
#    f.write(obser)
#f.close()

        
            
