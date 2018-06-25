#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 11:48:09 2018

@author: stangen
"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import EFA.efa_files.cfs_utilities_st as ut
import surface_obs.madis_example.madis_utilities as mt
import time
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from EFA.efa_xray.observation.observation import Observation
from load_data import Load_Data
import efa_functions as ef



ensemble_type = 'ecmwf'
variables = ['T2M', 'ALT']
obs = ['T2M']
#change this later
start_date = datetime(2013,4,1,0)#'20130401_0000'
end_date = datetime(2013,4,2,12)#'20130401_0000'

#change for how far into the forecast to get observations, 1 is 6 hours
start_index = 1
end_index = 1

save_dir = '/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/mse_var_output/'

#last forecast hour I want to get observations for
end_hour = 6*end_index

#the dict where all the data is
ob_dict = {}
#dict for the means from the data
data_dict = {}
#create dict for variable type
for ob_type in obs:
    ob_dict[ob_type] = ob_dict.get(ob_type,{})
    data_dict[ob_type] = data_dict.get(ob_type,{})
    #creates a datelist which increments in 12 or 24 hour chunks
    dates = mt.make_datetimelist(start_date, end_date, 12)
    for date in dates:
        hour = date.strftime('%H')
        
        
        #Load in the ensemble data for a given initilization
        efa = Load_Data(date,ensemble_type,variables,ob_type,[ob_type])
        statecls, lats, lons, elevs = efa.load_netcdfs()
        
        #Want to interpolate ALL observations valid during a given ensemble forecast
        hour_step = 6
        j = start_index
        fh = hour_step*j
        sfh = str(fh)
        while fh <= end_hour:
            #print(fh)
            obs = efa.load_obs(fh)
    
    
    
            
            #create dictionary for forecast hour
            ob_dict[ob_type][sfh] = ob_dict[ob_type].get(sfh,{})
            data_dict[ob_type][sfh] = data_dict[ob_type].get(sfh,{})
            
            #create dictionary for 00Z or 12Z
            ob_dict[ob_type][sfh][hour] = ob_dict[ob_type][sfh].get(hour,{})
            
            #Create dictionaries to hold the mean values to be calculated
            data_dict[ob_type][sfh]['station_all_mse'] = data_dict[ob_type][sfh].get('station_all_mse', np.array([]))
            data_dict[ob_type][sfh]['station_all_mean_variance'] = data_dict[ob_type][sfh].get('station_all_mean_variance', np.array([]))
            #data_dict[ob_type][sfh]['station_all_mse_unbiased'] = data_dict[ob_type][sfh].get('station_all_mse_unbiased', np.array([]))
            #data_dict[ob_type][sfh]['station_all_mean_variance_unbiased'] = data_dict[ob_type][sfh].get('station_all_mean_variance_unbiased', np.array([]))
            #data_dict[ob_type][sfh]['station_all_mse_no_bias'] = data_dict[ob_type][sfh].get('station_all_mse_no_bias', np.array([]))
            data_dict[ob_type][sfh]['station_all_error_variance'] = data_dict[ob_type][sfh].get('station_all_error_variance', np.array([]))
            data_dict[ob_type][sfh]['weight'] = data_dict[ob_type][sfh].get('weight',np.array([]))
            #data_dict[ob_type][sfh]['total_weight'] = data_dict[ob_type][sfh].get('total_weight',0)
            
            #can probably delete this
            obs_pass = []
            #i = 0
            #ob_counter = 0
            #while i < 10:
            for ob_counter, ob in enumerate(obs):
            #ob = obs[0]
            #print(ob)
                #ob=obs[i]
                #this gets the observation information from the text file
                ob_info = mt.get_ob_info(ob)
                #get longitude positive-definite- ie -130 lon is 230 E
                if ob_info['lon'] < 0:
                    ob_info['lon'] = ob_info['lon'] + 360
                utctime = datetime.utcfromtimestamp(ob_info['time'])
                #find interpolated ob estimate, if it passes the terrain check.
                #the terrain check is done within the closest_points function
                interp = ef.closest_points(ob_info['lat'],ob_info['lon'],lats,lons,ob_info['elev'],
                                           elevs,utctime,statecls['validtime'].values,
                                           statecls.variables[ob_type].values,need_interp=True)            

                ob_id = ob_info['name']
                ob_value = ob_info['ob']
                
                #if isinstance(interp,bool) == True and ob_id != 'SHIP':
                if len(interp) > 0 and ob_id != 'SHIP':

                    
                    #-----Added just to save obs used, to plot later-----------
                    obs_pass.append(ob)
                    
                    #create dictionary for station ID if it doesn't exist
                    ob_dict[ob_type][sfh][hour][ob_id] = ob_dict[ob_type][sfh][hour].get(ob_id,{})
                    #create empty lists for the station ID if they don't exist yet
                    #ob_dict[ob_type][sfh][hour][ob_id]['hx'] = ob_dict[ob_type][sfh][hour][ob_id].get('hx',[])
                    #ob_dict[ob_type][sfh][hour][ob_id]['ob_all'] = ob_dict[ob_type][sfh][hour][ob_id].get('ob_all',[])
                    #ob_dict[ob_type][sfh][hour][ob_id]['obs'] = ob_dict[ob_type][sfh][hour][ob_id].get('obs',[])
                    #ob_dict[ob_type][sfh][hour][ob_id]['hx_error'] = ob_dict[ob_type][sfh][hour][ob_id].get('hx_error',[])
                    #ob_dict[ob_type][sfh][hour][ob_id]['variance'] = ob_dict[ob_type][sfh][hour][ob_id].get('variance',[])
                    ob_dict[ob_type][sfh][hour][ob_id]['variance_unbiased'] = ob_dict[ob_type][sfh][hour][ob_id].get('variance_unbiased',[])
                    #ob_dict[ob_type][sfh][hour][ob_id]['hx_mean'] = ob_dict[ob_type][sfh][hour][ob_id].get('hx_mean',[])
                    ob_dict[ob_type][sfh][hour][ob_id]['se'] = ob_dict[ob_type][sfh][hour][ob_id].get('se',[])
                    ob_dict[ob_type][sfh][hour][ob_id]['error'] = ob_dict[ob_type][sfh][hour][ob_id].get('error',[])
                    
                    #add the data to the dictionaries
                    #ob_dict[ob_type][sfh][hour][ob_id]['ob_all'].append(ob)
                    #ob_dict[ob_type][sfh][hour][ob_id]['obs'].append(ob_value)
                    hx_oneob = interp
                    #ob_dict[ob_type][sfh][hour][ob_id]['hx'].append(hx_oneob)
                    #ob_dict[ob_type][sfh][hour][ob_id]['hx_error'].append(hx_oneob-ob_value)
                    hx_variance_unbiased = np.var(hx_oneob, ddof=1)
                    #hx_variance = np.var(hx_oneob)
                    #ob_dict[ob_type][sfh][hour][ob_id]['variance'].append(hx_variance)
                    ob_dict[ob_type][sfh][hour][ob_id]['variance_unbiased'].append(hx_variance_unbiased)                
                    #ob_dict[ob_type][sfh][hour][ob_id]['hx_mean'].append(np.mean(hx_oneob))
                    ob_dict[ob_type][sfh][hour][ob_id]['error'].append(np.mean(hx_oneob)-ob_value)
                    ob_dict[ob_type][sfh][hour][ob_id]['se'].append((np.mean(hx_oneob)-ob_value)**2)
                    
                    #ob_counter +=1
                    print("on observation "+str(ob_counter)+" out of "+str(len(obs)))                
                    
                        
                #print(len(variance))
                #i = i +1
            
            
            #update the forecast hour to load the next forecast, 6 hours later
            j = j + 1
            fh = hour_step*j
            #print(fh)
            sfh = str(fh)

#now go through all the data, convert to numpy arrays, and calculate the bias, 
#mse, and mean variance for each station and forecast hour

for var in ob_dict:
    #create list to save for creating plots
    stats_list = []
    for f_h in ob_dict[var]:
        for h in ob_dict[var][f_h]:
            for sID in ob_dict[var][f_h][h]:
                #ob_dict[var][f_h][h][sID]['hx'] = np.array(ob_dict[var][f_h][h][sID]['hx'])
                #ob_dict[var][f_h][h][sID]['hx_error'] = np.array(ob_dict[var][f_h][h][sID]['hx_error'])
                #ob_dict[var][f_h][h][sID]['obs'] = np.array(ob_dict[var][f_h][h][sID]['obs'])
                #ob_dict[var][f_h][h][sID]['variance'] = np.array(ob_dict[var][f_h][h][sID]['variance'])
                ob_dict[var][f_h][h][sID]['variance_unbiased'] = np.array(ob_dict[var][f_h][h][sID]['variance_unbiased'])
                #ob_dict[var][f_h][h][sID]['hx_mean'] = np.array(ob_dict[var][f_h][h][sID]['hx_mean'])
                ob_dict[var][f_h][h][sID]['se'] = np.array(ob_dict[var][f_h][h][sID]['se'])
                ob_dict[var][f_h][h][sID]['error'] = np.array(ob_dict[var][f_h][h][sID]['error'])
                
                #vertical calculate bias
                #hx_vertmean = np.mean(ob_dict[var][f_h][h][sID]['hx'],axis=0)
                #ob_vertmean = np.mean(ob_dict[var][f_h][h][sID]['obs'])
                #each member's bias
                #hx_bias = hx_vertmean-ob_vertmean
                #ob_dict[var][f_h][h][sID]['hx_bias'] = hx_bias
                #remove each member's bias from hx
                #hx_no_bias = ob_dict[var][f_h][h][sID]['hx'] - hx_bias
                #calculate variance of ensemble at each time
                #hx_no_bias_var = np.var(hx_no_bias, axis=1)
                #calculate squared error of ensemble mean for each time
                #hx_no_bias_mean = np.mean(hx_no_bias, axis=1)
                #hx_no_bias_se = (hx_no_bias_mean-ob_dict[var][f_h][h][sID]['obs'])**2
                #Greg's method of finding variance of error being equal to mse - bias^2
                error_var = np.var(ob_dict[var][f_h][h][sID]['error'])
                
                

                #calculate mse
                ob_dict[var][f_h][h][sID]['mse'] = np.mean(ob_dict[var][f_h][h][sID]['se'])
                #calculate mean variance (using ddof=1)
                ob_dict[var][f_h][h][sID]['mean_variance_unbiased'] = np.mean(ob_dict[var][f_h][h][sID]['variance_unbiased'])

                #PRACTICE STUFF
                #ob_dict[var][f_h][h][sID]['mean_variance_ne+1/ne'] = np.mean(ob_dict[var][f_h][h][sID]['variance'])*(len(hx_bias)+1)/len(hx_bias)
                #ob_dict[var][f_h][h][sID]['mean_variance'] = np.mean(ob_dict[var][f_h][h][sID]['variance'])


                
                #calculate bias removal from mse (mse-bias^2)
                #calculate the bias
                #bias1 = np.mean(ob_dict[var][f_h][h][sID]['hx_mean']-ob_dict[var][f_h][h][sID]['obs'])
                #ob_dict[var][f_h][h][sID]['bias1'] = bias1
                #ob_dict[var][f_h][h][sID]['mse_no_bias'] = ob_dict[var][f_h][h][sID]['mse'] - (ob_dict[var][f_h][h][sID]['bias1'])**2
                
                #calculate variance and mse of unbiased hx
                #ob_dict[var][f_h][h][sID]['mse_unbiased'] = np.mean(hx_no_bias_se)
                #ob_dict[var][f_h][h][sID]['mean_variance_unbiased'] = np.mean(hx_no_bias_var) CHANGE the name if you ever want to use this
                ob_dict[var][f_h][h][sID]['error_variance'] = error_var
                             

                data_dict[var][f_h]['station_all_mse'] = np.append(data_dict[var][f_h]['station_all_mse'],ob_dict[var][f_h][h][sID]['mse'])
                data_dict[var][f_h]['station_all_mean_variance'] = np.append(data_dict[var][f_h]['station_all_mean_variance'],ob_dict[var][f_h][h][sID]['mean_variance_unbiased'])
                #data_dict[var][f_h]['station_all_mse_unbiased'] = np.append(data_dict[var][f_h]['station_all_mse_unbiased'],ob_dict[var][f_h][h][sID]['mse_unbiased'])
                #data_dict[var][f_h]['station_all_mean_variance_unbiased'] = np.append(data_dict[var][f_h]['station_all_mean_variance_unbiased'],ob_dict[var][f_h][h][sID]['mean_variance_unbiased'])
                #data_dict[var][f_h]['station_all_mse_no_bias'] = np.append(data_dict[var][f_h]['station_all_mse_no_bias'],ob_dict[var][f_h][h][sID]['mse_no_bias'])
                data_dict[var][f_h]['station_all_error_variance'] = np.append(data_dict[var][f_h]['station_all_error_variance'],ob_dict[var][f_h][h][sID]['error_variance'])
                data_dict[var][f_h]['weight'] = np.append(data_dict[var][f_h]['weight'],len(ob_dict[var][f_h][h][sID]['se']))
                #data_dict[var][f_h]['total_weight'] = data_dict[var][f_h]['total_weight'] + len(ob_dict[var][f_h][h][sID]['se'])    
#        #data_dict[var][f_h]['station_all_mse_no_bias'] = np.sort(data_dict[var][f_h]['station_all_mse_no_bias'])
#        #data_dict[var][f_h]['station_all_mse_unbiased'] = np.sort(data_dict[var][f_h]['station_all_mse_unbiased'])
#        data_dict[var][f_h]['station_all_mse'] = np.sort(data_dict[var][f_h]['station_all_mse'])
#        data_dict[var][f_h]['station_all_mean_variance'] = np.sort(data_dict[var][f_h]['station_all_mean_variance'])
#        #data_dict[var][f_h]['station_all_mean_variance_unbiased'] = np.sort(data_dict[var][f_h]['station_all_mean_variance_unbiased'])
#        data_dict[var][f_h]['station_all_error_variance'] = np.sort(data_dict[var][f_h]['station_all_error_variance'])
        
        data_dict[var][f_h]['average_mse'] = np.average(data_dict[var][f_h]['station_all_mse'],weights=data_dict[var][f_h]['weight'])
        data_dict[var][f_h]['average_variance'] = np.average(data_dict[var][f_h]['station_all_mean_variance'],weights=data_dict[var][f_h]['weight'])
        data_dict[var][f_h]['average_error_variance'] = np.average(data_dict[var][f_h]['station_all_error_variance'],weights=data_dict[var][f_h]['weight'])
        #data_dict[var][f_h]['average_variance_hx_each_bias_removed'] = np.mean(data_dict[var][f_h]['station_all_mean_variance_unbiased'])
 
        #append to the stats list each stuff from forecast hour, for one variable
        stats_list.append(ensemble_type+','+var+','+f_h+','+str(data_dict[var][f_h]['average_mse'])+','+str(data_dict[var][f_h]['average_variance'])+','+str(data_dict[var][f_h]['average_error_variance'])+'\n')#','+str(data_dict[var][f_h]['average_variance_hx_each_bias_removed'])+'\n')
        
        #save the stats list after all forecast hours have been appended for one variable
        f = open(save_dir+ensemble_type+'_'+var+'_'+str(start_index*6)+'-'+str(end_index*6)+'.txt', 'a')
        for s in stats_list:
            f.write(s)
        f.close()
#want to save ens type, ob, forecast hour in string identifier, followed by mse, variance, and error variance

##this was to save all the stations that pass the elevation check so I could plot them. 
#f = open('/home/disk/hot/stangen/Documents/atms544/obsdoog.txt','w')
#for obser in obs_pass:
#    f.write(obser)
#f.close()

     

    
    
#    f = open(base_dir+"/"+dates[0:6]+"/"+ob+"_"+var_short+"/"+ob+"_"+var_short+"_"+str(dates)+".txt","w")
#    for s in sstr_one_station_1hr:
#        f.write(s)
#    f.close()

            
    
    ###Strategy: make a dictionary containing all obs, everything. Make a subdictionary
    # for each station ID. Make sure to check if the key (stationID) already exists. If it
    # doesn't, initialize all the other keys, use lists, I think this is necessary. Then add
    # the relevant values (hx_oneob, hx_variance, hx_mean maybe, and se) 
    # Can convert to numpy arrays after all the values are in, with this structure,
    # I can find the bias for each ensemble member for each station to find the 
    # ens mean bias, can average all the other quantities across all stations and times
        
            
