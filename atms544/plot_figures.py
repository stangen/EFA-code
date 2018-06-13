#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 26 14:42:14 2018

@author: stangen
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap


obs = ['T2M','ALT']
ens = ['ecmwf','eccc','ncep']
times = ['6-60','66-120','126-180','186-234']

var_dict = {
            'ALT' : 'Alt',
            'T2M' : '2m temp'
            }
            
var_units = {
            'ALT' : 'hPa$^{2}$',
            'T2M' : 'K$^{2}$'       
            }

save_dir = '/home/disk/hot/stangen/Documents/atms544/plots/'

#read in the data from the 24 files and put into a dict
stats_dict = {}

for e in ens:
    stats_dict[e] = stats_dict.get(e,{})
    for o in obs:
        stats_dict[e][o] = stats_dict[e].get(o,{})
        for t in times:
            #read in the corresponding stats file
            stats_dict[e][o]['raw_stats'] = stats_dict[e][o].get('raw_stats',[])
            stats_dict[e][o]['forecast_hour'] = stats_dict[e][o].get('forecast_hour',[])
            stats_dict[e][o]['MSE'] = stats_dict[e][o].get('MSE',[])
            stats_dict[e][o]['Ensemble_Variance'] = stats_dict[e][o].get('Ensemble_Variance',[])
            stats_dict[e][o]['Error_Variance'] = stats_dict[e][o].get('Error_Variance',[])
            stats_dict[e][o]['Total_Variance'] = stats_dict[e][o].get('Total_Variance',[])
            
            
            #other things: median, dividing them by each other?
            stats_path = '/home/disk/hot/stangen/Documents/atms544/stats/'+e+'_'+o+'_'+t+'.txt'
            f1 = open(stats_path, 'r')
            stats = f1.readlines()
            for s in stats:
                stats_dict[e][o]['raw_stats'].append(s)
        #now, put the stats in their own variable names
        for fh in stats_dict[e][o]['raw_stats']:
            broken_stats = fh.split(',')
            stats_dict[e][o]['forecast_hour'].append(int(broken_stats[0]))
            stats_dict[e][o]['MSE'].append(float(broken_stats[1]))
            stats_dict[e][o]['Ensemble_Variance'].append(float(broken_stats[2]))
            stats_dict[e][o]['Total_Variance'].append(float(broken_stats[2])+1)
            stats_dict[e][o]['Error_Variance'].append(float(broken_stats[3]))
        
        for key in stats_dict[e][o]:
            stats_dict[e][o][key] = np.array(stats_dict[e][o][key])
        stats_dict[e][o]['Var_MSE_Ratio'] = stats_dict[e][o]['Total_Variance']/stats_dict[e][o]['MSE']
        stats_dict[e][o]['Var_Error_Var_Ratio'] = stats_dict[e][o]['Total_Variance']/stats_dict[e][o]['Error_Variance']
#        stats_dict[e][o]['forecast_hour'] = np.array(stats_dict[e][o]['forecast_hour'])  
#        stats_dict[e][o]['MSE'] = np.array(stats_dict[e][o]['MSE']) 
#        stats_dict[e][o]['Ensemble_Variance'] = np.array(stats_dict[e][o]['Ensemble_Variance'])
#        stats_dict[e][o]['Error_Variance'] = np.array(stats_dict[e][o]['Error_Variance'])
         

f_h = stats_dict['ecmwf']['T2M']['forecast_hour']

for o in obs:

    #MSE vs Variance, T2M  
    fig = plt.figure(figsize=(14,8))
    
    plt.plot(f_h, stats_dict['ecmwf'][o]['MSE'], color='blue',linestyle='solid',label='ECMWF MSE')
    plt.plot(f_h, stats_dict['eccc'][o]['MSE'], color='red',linestyle='solid',label='CMC MSE')
    plt.plot(f_h, stats_dict['ncep'][o]['MSE'], color='green',linestyle='solid',label='GEFS MSE')
    
    plt.plot(f_h, stats_dict['ecmwf'][o]['Total_Variance'], color='blue',linestyle='dashed',label='ECMWF Var')
    plt.plot(f_h, stats_dict['eccc'][o]['Total_Variance'], color='red',linestyle='dashed',label='CMC Var')
    plt.plot(f_h, stats_dict['ncep'][o]['Total_Variance'], color='green',linestyle='dashed',label='GEFS Var')
    
    plt.xticks(np.arange(min(f_h)-6, max(f_h), 24.0))
    plt.ylim(ymin=0)
    plt.grid()
    plt.title(var_dict[o]+' MSE and Variance',fontsize=20)
    plt.xlabel('Forecast Hour',fontsize=14)
    plt.ylabel(var_units[o],fontsize=14)
    plt.legend(loc = 'upper left')
    
    fig.savefig(save_dir+'MSE_Var_'+o+'.png')
    
    #MSE vs Variance zoomed in
    fig = plt.figure(figsize=(14,8))

    plt.plot(f_h[0:12], stats_dict['ecmwf'][o]['MSE'][0:12], color='blue',linestyle='solid',label='ECMWF MSE')
    plt.plot(f_h[0:12], stats_dict['eccc'][o]['MSE'][0:12], color='red',linestyle='solid',label='CMC MSE')
    plt.plot(f_h[0:12], stats_dict['ncep'][o]['MSE'][0:12], color='green',linestyle='solid',label='GEFS MSE')
    
    plt.plot(f_h[0:12], stats_dict['ecmwf'][o]['Total_Variance'][0:12], color='blue',linestyle='dashed',label='ECMWF Var')
    plt.plot(f_h[0:12], stats_dict['eccc'][o]['Total_Variance'][0:12], color='red',linestyle='dashed',label='CMC Var')
    plt.plot(f_h[0:12], stats_dict['ncep'][o]['Total_Variance'][0:12], color='green',linestyle='dashed',label='GEFS Var')
    
    plt.xticks(np.arange(min(f_h)-6, f_h[11]+1, 24.0))
    plt.ylim(ymin=0)
    plt.grid()
    plt.title(var_dict[o]+' MSE and Variance',fontsize=20)
    plt.xlabel('Forecast Hour',fontsize=14)
    plt.ylabel(var_units[o],fontsize=14)
    plt.legend(loc = 'upper left')
    
    fig.savefig(save_dir+'MSE_Var_'+o+'_zoomed.png')
    
        
    #Variance/MSE
    fig = plt.figure(figsize=(14,8))
    
    plt.plot(f_h, stats_dict['ecmwf'][o]['Var_MSE_Ratio'], color='blue',linestyle='solid',label='ECMWF')
    plt.plot(f_h, stats_dict['eccc'][o]['Var_MSE_Ratio'], color='red',linestyle='solid',label='CMC')
    plt.plot(f_h, stats_dict['ncep'][o]['Var_MSE_Ratio'], color='green',linestyle='solid',label='GEFS')
    
    plt.xticks(np.arange(min(f_h)-6, max(f_h), 24.0))
    plt.grid()
    plt.title(var_dict[o]+' Variance/MSE',fontsize=20)
    plt.xlabel('Forecast Hour',fontsize=14)
    if o == 'T2M':
        plt.legend(loc = 'upper left')
    elif o =='ALT':
        plt.legend(loc = 'upper right')
        
    fig.savefig(save_dir+'Var_div_MSE_'+o+'.png')
    
    
    #Error Variance vs Variance
    fig = plt.figure(figsize=(14,8))
                
    plt.plot(f_h, stats_dict['ecmwf'][o]['Error_Variance'], color='blue',linestyle='solid',label='ECMWF Error Var')
    plt.plot(f_h, stats_dict['eccc'][o]['Error_Variance'], color='red',linestyle='solid',label='CMC Error Var')
    plt.plot(f_h, stats_dict['ncep'][o]['Error_Variance'], color='green',linestyle='solid',label='GEFS Error Var')
    
    plt.plot(f_h, stats_dict['ecmwf'][o]['Total_Variance'], color='blue',linestyle='dashed',label='ECMWF Var')
    plt.plot(f_h, stats_dict['eccc'][o]['Total_Variance'], color='red',linestyle='dashed',label='CMC Var')
    plt.plot(f_h, stats_dict['ncep'][o]['Total_Variance'], color='green',linestyle='dashed',label='GEFS Var')
    
    plt.xticks(np.arange(min(f_h)-6, max(f_h), 24.0))
    plt.ylim(ymin=0)
    plt.grid()
    plt.title(var_dict[o]+' Error Variance and Variance',fontsize=20)
    plt.xlabel('Forecast Hour',fontsize=14)
    plt.ylabel(var_units[o],fontsize=14)
    plt.legend(loc = 'upper left')
    
    fig.savefig(save_dir+'Var_Error_Var_'+o+'.png')
    
    #Error Variance vs Variance zoomed in
    fig = plt.figure(figsize=(14,8))
    
    plt.plot(f_h[0:12], stats_dict['ecmwf'][o]['Error_Variance'][0:12], color='blue',linestyle='solid',label='ECMWF MSE')
    plt.plot(f_h[0:12], stats_dict['eccc'][o]['Error_Variance'][0:12], color='red',linestyle='solid',label='CMC MSE')
    plt.plot(f_h[0:12], stats_dict['ncep'][o]['Error_Variance'][0:12], color='green',linestyle='solid',label='GEFS MSE')
    
    plt.plot(f_h[0:12], stats_dict['ecmwf'][o]['Total_Variance'][0:12], color='blue',linestyle='dashed',label='ECMWF Var')
    plt.plot(f_h[0:12], stats_dict['eccc'][o]['Total_Variance'][0:12], color='red',linestyle='dashed',label='CMC Var')
    plt.plot(f_h[0:12], stats_dict['ncep'][o]['Total_Variance'][0:12], color='green',linestyle='dashed',label='GEFS Var')
    
    plt.xticks(np.arange(min(f_h)-6, f_h[11]+1, 24.0))
    plt.ylim(ymin=0)
    plt.grid()
    plt.title(var_dict[o]+' MSE and Variance',fontsize=20)
    plt.xlabel('Forecast Hour',fontsize=14)
    plt.ylabel(var_units[o],fontsize=14)
    plt.legend(loc = 'upper left')
    
    fig.savefig(save_dir+'Var_Error_Var_'+o+'_zoomed.png')

    
    #Variance/Error Variance
    fig = plt.figure(figsize=(14,8))
    
    plt.plot(f_h, stats_dict['ecmwf'][o]['Var_Error_Var_Ratio'], color='blue',linestyle='solid',label='ECMWF')
    plt.plot(f_h, stats_dict['eccc'][o]['Var_Error_Var_Ratio'], color='red',linestyle='solid',label='CMC')
    plt.plot(f_h, stats_dict['ncep'][o]['Var_Error_Var_Ratio'], color='green',linestyle='solid',label='GEFS')
    
    plt.xticks(np.arange(min(f_h)-6, max(f_h), 24.0))
    plt.grid()
    plt.title(var_dict[o]+' Variance/Error Variance',fontsize=20)
    plt.xlabel('Forecast Hour',fontsize=14)
    if o == 'T2M':
        plt.legend(loc = 'upper left')
    elif o =='ALT':
        plt.legend(loc = 'upper right')
        
    fig.savefig(save_dir+'Var_div_Error_Var_'+o+'.png')
        
        
    #Variance/MSE and Variance/Error Variance
    fig = plt.figure(figsize=(14,8))
    
    plt.plot(f_h, stats_dict['ecmwf'][o]['Var_MSE_Ratio'], color='blue',linestyle='solid',label='ECMWF Var/MSE')
    plt.plot(f_h, stats_dict['eccc'][o]['Var_MSE_Ratio'], color='red',linestyle='solid',label='CMC Var/MSE')
    plt.plot(f_h, stats_dict['ncep'][o]['Var_MSE_Ratio'], color='green',linestyle='solid',label='GEFS Var/MSE')
    
    plt.plot(f_h, stats_dict['ecmwf'][o]['Var_Error_Var_Ratio'], color='blue',linestyle='dashed',label='ECMWF Var/EV')
    plt.plot(f_h, stats_dict['eccc'][o]['Var_Error_Var_Ratio'], color='red',linestyle='dashed',label='CMC Var/EV')
    plt.plot(f_h, stats_dict['ncep'][o]['Var_Error_Var_Ratio'], color='green',linestyle='dashed',label='GEFS Var/EV')
    
    plt.xticks(np.arange(min(f_h)-6, max(f_h), 24.0))
    plt.grid()
    plt.title(var_dict[o]+' Variance/MSE and Variance/Error Variance',fontsize=20)
    plt.xlabel('Forecast Hour',fontsize=14)
    if o == 'T2M':
        plt.legend(loc = 'upper left')
    elif o =='ALT':
        plt.legend(loc = 'upper right')
        
    fig.savefig(save_dir+'Var_div_MSE_Var_div_Error_Var_'+o+'.png')


