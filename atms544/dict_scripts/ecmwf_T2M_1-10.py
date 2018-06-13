#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 11:48:09 2018

@author: stangen
"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
from nicks_files.operational_cfsv2 import get_cfsv2_ensemble
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import efa_files.cfs_utilities_st as ut
import efa_files.madis_example.madis_utilities as mt
import time
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from efa_xray.state.ensemble import EnsembleState
from efa_xray.observation.observation import Observation
from efa_xray.assimilation.ensrf import EnSRF
import efa_xray.postprocess.postprocess as pp
from efa_files.do_efa import Run_efa



ensemble_type = 'ecmwf'
variables = ['T2M', 'ALT']
obs = ['T2M']
#change this later
start_date = '20130401_0000'
end_date = '20130430_0000'

start_index = 1
end_index = 10

save_dir = '/home/disk/hot/stangen/Documents/atms544/stats/'

#last forecast hour I want to get observations for
end_hour = 6*end_index

var_dict = {
        'ALT' : 'alts',
        'T2M' : 'temp'
        }

m_dict = {
           '01' : 'jan',
           '02' : 'feb',
           '03' : 'mar',
           '04' : 'apr',
           '05' : 'may',
           '06' : 'jun',
           '07' : 'jul',
           '08' : 'aug',
           '09' : 'sep',
           '10' : 'oct',
           '11' : 'nov',
           '12' : 'dec'           
           }

##the dict where all the data is
#ob_dict = {}
##dict for the means from the data
#data_dict = {}
##create dict for variable type
#for ob_type in obs:
#    ob_dict[ob_type] = ob_dict.get(ob_type,{})
#    data_dict[ob_type] = data_dict.get(ob_type,{})
#    
#    #creates a datelist which increments in 12 or 24 hour chunks
#    dates = mt.make_datelist(start_date, end_date, torf=True, timestep=3600*24)
#    for datestr in dates:
#        year = datestr[0:4]
#        month = datestr[4:6]
#        day = datestr[6:8]
#        hour = datestr[9:11]
#        
#        
#        #Load in the ensemble data for a given initilization
#        efa = Run_efa(datestr,ensemble_type,variables,m_dict,var_dict[ob_type])
#        statecls, lats, lons, elevs = efa.load_data()
#        
#        #Want to interpolate ALL observations valid during a given ensemble forecast
#        hour_step = 6
#        j = start_index
#        fh = hour_step*j
#        sfh = str(fh)
#        while fh <= end_hour:
#            #print(fh)
#            obs = efa.load_obs(fh)
#    
#    
#    
#            
#            #create dictionary for forecast hour
#            ob_dict[ob_type][sfh] = ob_dict[ob_type].get(sfh,{})
#            data_dict[ob_type][sfh] = data_dict[ob_type].get(sfh,{})
#            
#            #create dictionary for 00Z or 12Z
#            ob_dict[ob_type][sfh][hour] = ob_dict[ob_type][sfh].get(hour,{})
#            
#            #Create dictionaries to hold the mean values to be calculated
#            data_dict[ob_type][sfh]['station_all_mse'] = data_dict[ob_type][sfh].get('station_all_mse', np.array([]))
#            data_dict[ob_type][sfh]['station_all_mean_variance'] = data_dict[ob_type][sfh].get('station_all_mean_variance', np.array([]))
#            data_dict[ob_type][sfh]['station_all_mse_unbiased'] = data_dict[ob_type][sfh].get('station_all_mse_unbiased', np.array([]))
#            data_dict[ob_type][sfh]['station_all_mean_variance_unbiased'] = data_dict[ob_type][sfh].get('station_all_mean_variance_unbiased', np.array([]))
#            data_dict[ob_type][sfh]['station_all_mse_no_bias'] = data_dict[ob_type][sfh].get('station_all_mse_no_bias', np.array([]))
#            data_dict[ob_type][sfh]['station_all_error_variance'] = data_dict[ob_type][sfh].get('station_all_error_variance', np.array([]))
#            
#    
#            #i = 0
#            #ob_counter = 0
#            #while i < 10:
#            for ob_counter, ob in enumerate(obs):
#            #ob = obs[0]
#            #print(ob)
#                #ob=obs[i]
#                ob_split = ob.split(',')
#                #get the lat/lon of the station
#                ob_id = ob_split[0]
#                ob_lat = float(ob_split[1])
#                ob_lon = float(ob_split[2])
#                #get longitude positive-definite- ie -130 lon is 230 E
#            #    if ob_lon < 0:
#            #        ob_lon = ob_lon + 360
#                ob_elev = float(ob_split[3])
#                ob_time = mt.timestamp2utc(float(ob_split[4]))
#                ob_value = float(ob_split[5])
#                
#                #call function to check elevation to see whether or not to assimilate ob
#                TorF = efa.check_elevation(lats,lons,elevs,ob,ob_lat,ob_lon,ob_elev)
#                
#                obser = Observation(value=ob_value, time=ob_time,
#                                lat=ob_lat,lon=ob_lon, obtype=ob_type, localize_radius=1000,
#                                assimilate_this=TorF, error=1)
#                
#                if TorF == True and ob_id != 'SHIP':
#                    
#                    #create dictionary for station ID if it doesn't exist
#                    ob_dict[ob_type][sfh][hour][ob_id] = ob_dict[ob_type][sfh][hour].get(ob_id,{})
#                    #create empty lists for the station ID if they don't exist yet
#                    ob_dict[ob_type][sfh][hour][ob_id]['hx'] = ob_dict[ob_type][sfh][hour][ob_id].get('hx',[])
#                    ob_dict[ob_type][sfh][hour][ob_id]['ob_all'] = ob_dict[ob_type][sfh][hour][ob_id].get('ob_all',[])
#                    ob_dict[ob_type][sfh][hour][ob_id]['obs'] = ob_dict[ob_type][sfh][hour][ob_id].get('obs',[])
#                    #ob_dict[ob_type][sfh][hour][ob_id]['hx_error'] = ob_dict[ob_type][sfh][hour][ob_id].get('hx_error',[])
#                    ob_dict[ob_type][sfh][hour][ob_id]['variance'] = ob_dict[ob_type][sfh][hour][ob_id].get('variance',[])
#                    ob_dict[ob_type][sfh][hour][ob_id]['variance_unbiased'] = ob_dict[ob_type][sfh][hour][ob_id].get('variance_unbiased',[])
#                    ob_dict[ob_type][sfh][hour][ob_id]['hx_mean'] = ob_dict[ob_type][sfh][hour][ob_id].get('hx_mean',[])
#                    ob_dict[ob_type][sfh][hour][ob_id]['se'] = ob_dict[ob_type][sfh][hour][ob_id].get('se',[])
#                    ob_dict[ob_type][sfh][hour][ob_id]['error'] = ob_dict[ob_type][sfh][hour][ob_id].get('error',[])
#                    
#                    #add the data to the dictionaries
#                    ob_dict[ob_type][sfh][hour][ob_id]['ob_all'].append(ob)
#                    ob_dict[ob_type][sfh][hour][ob_id]['obs'].append(ob_value)
#                    hx_oneob = obser.estimate(statecls)
#                    ob_dict[ob_type][sfh][hour][ob_id]['hx'].append(hx_oneob)
#                    #ob_dict[ob_type][sfh][hour][ob_id]['hx_error'].append(hx_oneob-ob_value)
#                    hx_variance = np.var(hx_oneob)
#                    hx_variance_unbiased = np.var(hx_oneob, ddof=1)
#                    ob_dict[ob_type][sfh][hour][ob_id]['variance'].append(hx_variance)
#                    ob_dict[ob_type][sfh][hour][ob_id]['variance_unbiased'].append(hx_variance_unbiased)                
#                    ob_dict[ob_type][sfh][hour][ob_id]['hx_mean'].append(np.mean(hx_oneob))
#                    ob_dict[ob_type][sfh][hour][ob_id]['error'].append(np.mean(hx_oneob)-ob_value)
#                    ob_dict[ob_type][sfh][hour][ob_id]['se'].append((np.mean(hx_oneob)-ob_value)**2)
#                    
#                    #ob_counter +=1
#                    print("on observation "+str(ob_counter)+" out of "+str(len(obs)))                
#                    
#    
#                    #I think bias will be calculated at the end, all these can be calculated as obs come in
#    
#    
#    
#    
#    
#    #                
#    #                ob_all.append(ob)
#    #                observations.append(ob_value)
#    #                hx_oneob = obser.estimate(statecls)
#    #                hx.append(hx_oneob)
#    #                hx_variance = np.var(hx_oneob)
#    #                #check for outlier data, 4 standard deviations from mean
#    #                std_dev = np.sqrt(hx_variance)
#    #                #if abs(ob_value-np.mean(hx_oneob)) < 10*std_dev:
#    #                variance.append(hx_variance)
#    #                hx_mean.append(np.mean(hx_oneob))
#    #                se.append((np.mean(hx_oneob)-ob_value)**2)
#                    
#                #print(len(variance))
#                #i = i +1
#            
#            
#            #update the forecast hour to load the next forecast, 6 hours later
#            j = j + 1
#            fh = hour_step*j
#            #print(fh)
#            sfh = str(fh)

#now go through all the data, convert to numpy arrays, and calculate the bias, 
#mse, and mean variance for each station and forecast hour
#station_all_mse = np.array([])
#station_all_mean_variance = np.array([])
#station_all_mse_unbiased = np.array([])
#station_all_mean_variance_unbiased = np.array([])
#station_all_mse_no_bias = np.array([])
#station_all_error_variance = np.array([])

for var in ob_dict:
    #create list to save for creating plots
    stats_list = []
    for f_h in ob_dict[var]:
        for h in ob_dict[var][f_h]:
            for sID in ob_dict[var][f_h][h]:
                ob_dict[var][f_h][h][sID]['hx'] = np.array(ob_dict[var][f_h][h][sID]['hx'])
                #ob_dict[var][f_h][h][sID]['hx_error'] = np.array(ob_dict[var][f_h][h][sID]['hx_error'])
                ob_dict[var][f_h][h][sID]['obs'] = np.array(ob_dict[var][f_h][h][sID]['obs'])
                ob_dict[var][f_h][h][sID]['variance'] = np.array(ob_dict[var][f_h][h][sID]['variance'])
                ob_dict[var][f_h][h][sID]['variance_unbiased'] = np.array(ob_dict[var][f_h][h][sID]['variance_unbiased'])
                ob_dict[var][f_h][h][sID]['hx_mean'] = np.array(ob_dict[var][f_h][h][sID]['hx_mean'])
                ob_dict[var][f_h][h][sID]['se'] = np.array(ob_dict[var][f_h][h][sID]['se'])
                ob_dict[var][f_h][h][sID]['error'] = np.array(ob_dict[var][f_h][h][sID]['error'])
                
                
                #calculate the bias
                bias1 = np.mean(ob_dict[var][f_h][h][sID]['hx_mean']-ob_dict[var][f_h][h][sID]['obs'])
                ob_dict[var][f_h][h][sID]['bias1'] = bias1
                #vertical calculate bias
                hx_vertmean = np.mean(ob_dict[var][f_h][h][sID]['hx'],axis=0)
                ob_vertmean = np.mean(ob_dict[var][f_h][h][sID]['obs'])
                #each member's bias
                hx_bias = hx_vertmean-ob_vertmean
                ob_dict[var][f_h][h][sID]['hx_bias'] = hx_bias
                #remove each member's bias from hx
                hx_no_bias = ob_dict[var][f_h][h][sID]['hx'] - hx_bias
                #calculate variance of ensemble at each time
                hx_no_bias_var = np.var(hx_no_bias, axis=1)
                #calculate squared error of ensemble mean for each time
                hx_no_bias_mean = np.mean(hx_no_bias, axis=1)
                obs = ob_dict[var][f_h][h][sID]['obs']
                hx_no_bias_se = (hx_no_bias_mean-ob_dict[var][f_h][h][sID]['obs'])**2
                #Greg's method of finding variance of error being equal to mse - bias^2
                error_var = np.var(ob_dict[var][f_h][h][sID]['error'])
                
                
#                bias2 = np.mean(hx_vertmean-ob_vertmean)
#                ob_dict[var][f_h][h][sID]['bias2'] = bias2
                #calculate mse
                ob_dict[var][f_h][h][sID]['mse'] = np.mean(ob_dict[var][f_h][h][sID]['se'])
                #calculate mean variance
                ob_dict[var][f_h][h][sID]['mean_variance'] = np.mean(ob_dict[var][f_h][h][sID]['variance'])
                
                #PRACTICE STUFF
                ob_dict[var][f_h][h][sID]['mean_variance_ub'] = np.mean(ob_dict[var][f_h][h][sID]['variance_unbiased'])
                ob_dict[var][f_h][h][sID]['mean_variance_ne+1/ne'] = np.mean(ob_dict[var][f_h][h][sID]['variance'])*(len(hx_bias)+1)/len(hx_bias)

                
                #calculate bias removal from mse
                ob_dict[var][f_h][h][sID]['mse_no_bias'] = ob_dict[var][f_h][h][sID]['mse'] - (ob_dict[var][f_h][h][sID]['bias1'])**2
                
                #calculate variance and mse of unbiased hx
                ob_dict[var][f_h][h][sID]['mse_unbiased'] = np.mean(hx_no_bias_se)
                ob_dict[var][f_h][h][sID]['mean_variance_unbiased'] = np.mean(hx_no_bias_var)
                ob_dict[var][f_h][h][sID]['error_variance'] = error_var
                
                #create lists to put the means from each station
#                station_all_mse = np.append(station_all_mse, ob_dict[var][f_h][h][sID]['mse'])
#                station_all_mean_variance = np.append(station_all_mean_variance, ob_dict[var][f_h][h][sID]['mean_variance'])
#                station_all_mse_unbiased = np.append(station_all_mse_unbiased, ob_dict[var][f_h][h][sID]['mse_unbiased'])
#                station_all_mean_variance_unbiased = np.append(station_all_mean_variance_unbiased, ob_dict[var][f_h][h][sID]['mean_variance_unbiased'])
#                station_all_mse_no_bias = np.append(station_all_mse_no_bias, ob_dict[var][f_h][h][sID]['mse_no_bias'])
#                station_all_error_variance = np.append(station_all_error_variance, ob_dict[var][f_h][h][sID]['error_variance'])
#                

                data_dict[var][f_h]['station_all_mse'] = np.append(data_dict[var][f_h]['station_all_mse'],ob_dict[var][f_h][h][sID]['mse'])
                data_dict[var][f_h]['station_all_mean_variance'] = np.append(data_dict[var][f_h]['station_all_mean_variance'],ob_dict[var][f_h][h][sID]['mean_variance_ub'])
                #data_dict[var][f_h]['station_all_mse_unbiased'] = np.append(data_dict[var][f_h]['station_all_mse_unbiased'],ob_dict[var][f_h][h][sID]['mse_unbiased'])
                data_dict[var][f_h]['station_all_mean_variance_unbiased'] = np.append(data_dict[var][f_h]['station_all_mean_variance_unbiased'],ob_dict[var][f_h][h][sID]['mean_variance_unbiased'])
                #data_dict[var][f_h]['station_all_mse_no_bias'] = np.append(data_dict[var][f_h]['station_all_mse_no_bias'],ob_dict[var][f_h][h][sID]['mse_no_bias'])
                data_dict[var][f_h]['station_all_error_variance'] = np.append(data_dict[var][f_h]['station_all_error_variance'],ob_dict[var][f_h][h][sID]['error_variance'])
                       
        
        data_dict[var][f_h]['station_all_mse'] = np.sort(data_dict[var][f_h]['station_all_mse'])
        data_dict[var][f_h]['station_all_mean_variance'] = np.sort(data_dict[var][f_h]['station_all_mean_variance'])
        data_dict[var][f_h]['station_all_mean_variance_unbiased'] = np.sort(data_dict[var][f_h]['station_all_mean_variance_unbiased'])
        data_dict[var][f_h]['station_all_error_variance'] = np.sort(data_dict[var][f_h]['station_all_error_variance'])
        
        data_dict[var][f_h]['average_mse'] = np.mean(data_dict[var][f_h]['station_all_mse'])
        data_dict[var][f_h]['average_variance'] = np.mean(data_dict[var][f_h]['station_all_mean_variance'])
        data_dict[var][f_h]['average_error_variance'] = np.mean(data_dict[var][f_h]['station_all_error_variance'])
        data_dict[var][f_h]['average_variance_hx_each_bias_removed'] = np.mean(data_dict[var][f_h]['station_all_mean_variance_unbiased'])
        
        mse_median = np.median(data_dict[var][f_h]['station_all_mse'])
        mean_variance_median = np.median(data_dict[var][f_h]['station_all_mean_variance'])
        mean_variance_unbiased_median = np.median(data_dict[var][f_h]['station_all_mean_variance_unbiased'])
        error_variance_median = np.median(data_dict[var][f_h]['station_all_error_variance'])
    
    
    
        #append to the stats list each stuff from forecast hour, for one variable
        stats_list.append(f_h+','+str(data_dict[var][f_h]['average_mse'])+','+str(data_dict[var][f_h]['average_variance'])+','+str(data_dict[var][f_h]['average_error_variance'])+','+str(data_dict[var][f_h]['average_variance_hx_each_bias_removed'])+','+str(mse_median)+','+str(mean_variance_median)+','+str(error_variance_median)+','+str(mean_variance_unbiased_median)+'\n')
        
    #save the stats list after all forecast hours have been appended for one variable
    f = open(save_dir+ensemble_type+'_'+var+'_'+str(start_index*6)+'-'+str(end_index*6)+'.txt', 'w')
    for s in stats_list:
        f.write(s)
    f.close()
    
    
#    f = open(base_dir+"/"+dates[0:6]+"/"+ob+"_"+var_short+"/"+ob+"_"+var_short+"_"+str(dates)+".txt","w")
#    for s in sstr_one_station_1hr:
#        f.write(s)
#    f.close()

####Next up: calculate averages across all stations for a given forecast hour, should be pretty easy!
#for ob in ob_type:
#    for fh in sfh:
        

#        variance = np.array(variance)
#        hx_mean = np.array(hx_mean)
#        se = np.array(se)
#        
#        ##calculate the mean squared error
#        mse = np.mean(se)
#        print(mse)
#        ##calculate the mean variance
#        mean_var = np.mean(variance)
#        print(mean_var)
            
    
    ###Strategy: make a dictionary containing all obs, everything. Make a subdictionary
    # for each station ID. Make sure to check if the key (stationID) already exists. If it
    # doesn't, initialize all the other keys, use lists, I think this is necessary. Then add
    # the relevant values (hx_oneob, hx_variance, hx_mean maybe, and se) 
    # Can convert to numpy arrays after all the values are in, with this structure,
    # I can find the bias for each ensemble member for each station to find the 
    # ens mean bias, can average all the other quantities across all stations and times
        
        
            
        #thing = pp.obs_assimilation_statistics(statecls, observations)
        
            
