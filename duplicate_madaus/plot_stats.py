#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 11:49:37 2018

@author: stangen
"""

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import efa_functions as ef


start_date = datetime(2013,4,1,0)
end_date = datetime(2013,6,30,12)

variables = ['ALT']

var_units = {
            'ALT' : 'hPa$^{2}$',
            'T2M' : 'K$^{2}$'       
            }

ls = {
      'prior' : 'solid',
      'loc1000' : 'dashed',
      'loc500' : 'dotted'
      }

clr = {
        'ecmwf' : 'r',
        'eccc' : 'k'
        }

filedir = '/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/mse_var_output/'

savedir = '/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/plots/'

#create strings for loading txt file containing the stats
sy = start_date.strftime('%Y')
sm = start_date.strftime('%m')
sd = start_date.strftime('%d')
sh = start_date.strftime('%H')

ey = end_date.strftime('%Y')
em = end_date.strftime('%m')
ed = end_date.strftime('%d')
eh = end_date.strftime('%H')

datestr=sy+sm+sd+sh+'-'+ey+em+ed+eh

varstr = ef.var_string(variables)

filepath = filedir+datestr+'_'+varstr+.txt'

f1 = open(filepath, 'r')
stats = f1.readlines()

stats_dict = {}

for line in stats:
    line_split = line.split(',')
    efa = line_split[0] #posterior or prior
    ens = line_split[1] #ensemble type
    var = line_split[2] #observation variable type
    fh = line_split[3] #forecast hour
    mse = line_split[4] #mean squared error
    variance = line_split[5] #variance
    errvar = line_split[6] #error variance
    #this block is to set up the dictionary structure
    stats_dict[var] = stats_dict.get(var,{})
    stats_dict[var][ens] = stats_dict[var].get(ens, {})
    stats_dict[var][ens][efa] = stats_dict[var][ens].get(efa, {})
    stats_dict[var][ens][efa]['Forecast_Hour'] = stats_dict[var][ens][efa].get('Forecast_Hour',[])
    stats_dict[var][ens][efa]['MSE'] = stats_dict[var][ens][efa].get('MSE',[])
    stats_dict[var][ens][efa]['Ensemble_Variance'] = stats_dict[var][ens][efa].get('Ensemble_Variance',[])
    stats_dict[var][ens][efa]['Error_Variance'] = stats_dict[var][ens][efa].get('Error_Variance',[])
    
    #add the data to the dictionary
    stats_dict[var][ens][efa]['Forecast_Hour'].append(int(fh))
    stats_dict[var][ens][efa]['MSE'].append(float(mse))
    stats_dict[var][ens][efa]['Ensemble_Variance'].append(float(variance))
    stats_dict[var][ens][efa]['Error_Variance'].append(float(errvar))

for v in stats_dict:
    for e in stats_dict[v]:
        for E in stats_dict[v][e]:
            #convert each array to a numpy array so it can be sorted
            stats_dict[v][e][E]['Forecast_Hour'] = np.array(stats_dict[v][e][E]['Forecast_Hour'])
            stats_dict[v][e][E]['MSE'] = np.array(stats_dict[v][e][E]['MSE'])
            stats_dict[v][e][E]['Ensemble_Variance'] = np.array(stats_dict[v][e][E]['Ensemble_Variance'])
            stats_dict[v][e][E]['Error_Variance'] = np.array(stats_dict[v][e][E]['Error_Variance'])
            
            #sort each array in order of ascending forecast hour
            sorted_ind = stats_dict[v][e][E]['Forecast_Hour'].argsort()
            stats_dict[v][e][E]['Forecast_Hour'] = stats_dict[v][e][E]['Forecast_Hour'][sorted_ind]
            stats_dict[v][e][E]['MSE'] = stats_dict[v][e][E]['MSE'][sorted_ind]
            stats_dict[v][e][E]['Ensemble_Variance'] = stats_dict[v][e][E]['Ensemble_Variance'][sorted_ind]
            stats_dict[v][e][E]['Error_Variance'] = stats_dict[v][e][E]['Error_Variance'][sorted_ind]
            
            
#plotting variables
stats_list = ['MSE','Ensemble_Variance']
#now plot 
#separate plot for each variable type
for v in variables:
    #separate plot for MSE and variance
    for s in stats_list:
        fig = plt.figure(figsize=(14,8))
        #line for each ensemble type
        for m in stats_dict[v]:
            #line for prior/localization radius
            for l in stats_dict[v][m]:
                plt.plot(stats_dict[v][m][l]['Forecast_Hour'],stats_dict[v][m][l][s],
                         linestyle=ls[l],marker='o',color=clr[m],label=l+' '+m)
        
        
#        plt.plot(stats_dict[v]['ecmwf']['prior']['Forecast_Hour'],stats_dict[v]['ecmwf']['prior'][s],
#                 linestyle='dashed',marker='o',color='r',label='Original ECMWF')
#        plt.plot(stats_dict[v]['ecmwf']['posterior']['Forecast_Hour'],stats_dict[v]['ecmwf']['posterior'][s],
#                 marker='o',color='r',label='Updated ECMWF')
#        plt.plot(stats_dict[v]['eccc']['prior']['Forecast_Hour'],stats_dict[v]['eccc']['prior'][s],
#                 linestyle='dashed',marker='o',color='k',label='Original CMC')
#        plt.plot(stats_dict[v]['eccc']['posterior']['Forecast_Hour'],stats_dict[v]['eccc']['posterior'][s],
#                 marker='o',color='k',label='Updated CMC')   
        
        plt.xticks(np.arange(min(stats_dict[var]['ecmwf']['prior']['Forecast_Hour'])-6, 
                             max(stats_dict[var]['ecmwf']['prior']['Forecast_Hour'])+6, 6))
        plt.grid()
        plt.legend(loc = 'upper left')
        plt.title(v+' '+s,fontsize=20)
        plt.xlabel('Forecast Hour',fontsize=14)
        plt.ylabel(var_units[v],fontsize=14)
        
        fig.savefig(savedir+v+'_'+s+'_'+datestr+'.png')
    