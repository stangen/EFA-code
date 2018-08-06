#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 11:49:37 2018

@author: stangen
"""

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import EFA.duplicate_madaus.efa_functions as ef


start_date = datetime(2015,11,10,0)
end_date = datetime(2015,11,15,12)

#which ob type(s) were assimilated (8/2)? used to load the .txt file
#prior to 8/2, naming convention was all types of obs the .txt file contained.
assim_obs = ['IVT']#['TCW']#

grid = [-180,180,90,0,3]

#want to plot each variable or control which are plotted?
control_vars = True

#plot_vars = ['IWV','IWV100','IWV1000','IWV10000']#['QF850','QF850100','QF850750','QF8501000']#,'QF850250','QF850500']
plot_vars = ['IVT','IVT100','IVT1000','IVT10000']

AR_specific = True

separate_plots = False

#prior_var_key = 'ALT'

var_units = {
            'ALT' : 'hPa$^{2}$',
            'T2M' : 'K$^{2}$',
            'QF850' : '(g/kg*m/s)$^{2}$',
            'TCW' : 'mm',
            'IWV' : 'mm',
            'IVT' : 'kg/m/s'
            }

ls = {
      'prior' : 'solid',
      'posterior' : 'dashed',
      'loc3000' : ' ',
      'loc2000' : '-',
      'loc1000' : 'dashed',
      'loc500' : 'dotted',
      'loc100' : '-.'
      
      }

ls2 = {
       'QF850' : '-', #solid
       'QF850ensvar' : '--', #dashed
       'QF8501000' : ':', #dotted
       'QF850100' : '-.', #dashdotted
       'QF85010' : ' ',
       'QF8501' : '' ,
       'QF850250' : '--',
       'QF850500' : '--',
       'QF850750' : '--',
       
       'TCW' : '-',
       'TCW0-1' : ' ',
       'TCW1' : '--',
       'TCW10' : '-.',
       'TCW100' : ':',
       
       'IWV' : '-',
       'IWV1' : '',
       'IWV5' : ':',
       'IWV10' : '--',
       'IWV20' : '-.',
       'IWV100' : ':',
       'IWV1000' : '',
       'IWV5000' : '-.',
       'IWV10000' : '--',
       'IWV20000' : ':',
       
       'IVT' : '-',
       'IVT1' : '',
       'IVT5' : ':',
       'IVT10' : '--',
       'IVT20' : '-.',
       'IVT100' : ':',
       'IVT1000' : '',
       'IVT5000' : '-.',
       'IVT10000' : '--',
       'IVT20000' : ':'
       
       }

md = {
       'QF850' : 'o', #solid
       'QF850ensvar' : 'o', #dashed
       'QF8501000' : 'o', #dotted
       'QF850100' : 'o', #dashdotted
       'QF850250' : 'x',
       'QF850500' : '*',
       'QF850750' : 's',
       'QF85010' : 'x',
       'QF8501' : '*' ,
       
       'TCW' : 'o',
       'TCW0-1' : 'o',
       'TCW1' : 'o',
       'TCW10' : 'o',
       'TCW100' : 'o',
              
       'IWV' : 'o',
       'IWV1' : 'o',
       'IWV5' : 'o',
       'IWV10' : 'o',
       'IWV20' : 'o',
       'IWV100' : 'x',
       'IWV1000' : 'o',
       'IWV5000' : 'o',
       'IWV10000' : 'o',
       'IWV20000' : 'o',
       
       'IVT' : 'o',
       'IVT1' : 'o',
       'IVT5' : 'o',
       'IVT10' : 'o',
       'IVT20' : 'o',
       'IVT100' : 'x',
       'IVT1000' : 'o',
       'IVT5000' : 'o',
       'IVT10000' : 'o',
       'IVT20000' : 'o'
       
      }

clr = {
        'ecmwf' : 'r',
        'eccc' : 'k',
        'ncep' : 'b'
        }
#color for separate plot for each ensemble type
clr_sp = {
        'QF850' : 'k', #solid
        'QF850ensvar' : 'b', #dashed
        'QF8501000' : 'r', #dotted
        'QF850100' : 'g', #dashdotted
        'QF85010' : 'c',
        'QF8501' : 'y', 
        'QF850250' : 'm',
        'QF850500' : 'y',
        'QF850750' : 'c',
        
        'TCW' : 'k',
        'TCW0-1' : 'y',
        'TCW1' : 'g',
        'TCW10' : 'b',
        'TCW100' : 'r',
        
        'IWV' : 'k',
        'IWV1' : 'y',
        'IWV5' : 'r',
        'IWV10' : 'g',
        'IWV20' : 'c',
        'IWV100' : 'b',
        'IWV1000' : 'r',
        'IWV5000' : 'y',
        'IWV10000' : 'g',
        'IWV20000' : 'c',
        
       
        'IVT' : 'k',
        'IVT1' : 'y',
        'IVT5' : 'r',
        'IVT10' : 'g',
        'IVT20' : 'c',
        'IVT100' : 'b',
        'IVT1000' : 'r',
        'IVT5000' : 'y',
        'IVT10000' : 'g',
        'IVT20000' : 'c'
               
        }

ens_dict = {
        'eccc': 'CMC',
        'ecmwf' : 'ECMWF',
        'ncep' : 'GEFS'
        
        }

title_dict = { 
        'QF850' : '850mb Moisture Flux ',
        'TCW' : 'Total Column Water ',
        'IWV' : 'Integrated Water Vapor ',
        'IVT' : 'Integrated Vapor Transport '
        }

filedir = '/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/mse_var_output/'

savedir = '/home/disk/hot/stangen/Documents/EFA/AR/plots/'

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

varstr = ef.var_string(assim_obs)

gridstr = ef.var_string(grid)

#if AR_specific == True:
#    filepath = filedir+datestr+'_'+varstr+'_'+gridstr+'_gridobs_ARonly.txt'
#else:
#    filepath = filedir+datestr+'_'+varstr+'_'+gridstr+'_gridobs.txt'

filepath = filedir+datestr+'_'+varstr+'_'+gridstr+'_gridobs.txt'

f1 = open(filepath, 'r')
stats = f1.readlines()

stats_dict = {}

#-----------Changed it so that var is now ensemble type, ens is now observation variable type------

for line in stats:
    line_split = line.split(',')
    inf = line_split[0] #inflation
    efa = line_split[1] #localization radius or prior
    var = line_split[2] #ensemble type
    ens = line_split[3] #observation variable type
    fh = line_split[4] #forecast hour
    mse_gridobs = line_split[5] #mean squared error of gridded obs
    variance_gridobs = line_split[6] #variance of gridded obs
    #if we are plotting the AR-specific region
    #if AR_specific == False:
    mse_gridall = line_split[7] #mean squared error of entire grid
    variance_gridall = line_split[8] #variance of entire grid
    mse_grid_reg1 = line_split[9] #region 1 = -135:-115W, 55:30 N
    variance_grid_reg1 = line_split[10] 
    mse_grid_reg2 = line_split[11] #region2 = -180:-115W, 50:35 N
    variance_grid_reg2 = line_split[12] 

#----To put loc radius back in, add [efa] after [ens]
    #this block is to set up the dictionary structure
    stats_dict[var] = stats_dict.get(var,{})
    stats_dict[var][ens] = stats_dict[var].get(ens, {})
    srt = stats_dict[var][ens]
    
    srt['Forecast_Hour_Gridded'] = srt.get('Forecast_Hour_Gridded',[])
    #add the data to the dictionary
    srt['Forecast_Hour_Gridded'].append(int(fh)) #forecast hours of just these ones
    
    if AR_specific == False:
        srt['MSE_Grid_Obs'] = srt.get('MSE_Grid_Obs',[])
        srt['Variance_Grid_Obs'] = srt.get('Variance_Grid_Obs',[])
        srt['MSE_Grid_Obs'].append(float(mse_gridobs))
        srt['Variance_Grid_Obs'].append(float(variance_gridobs))
    
    elif AR_specific == True:
        srt['MSE_AR_Gridpoints'] = srt.get('MSE_AR_Gridpoints',[])
        srt['Variance_AR_Gridpoints'] = srt.get('Variance_AR_Gridpoints',[])    
        srt['MSE_AR_Gridpoints'].append(float(mse_gridobs))
        srt['Variance_AR_Gridpoints'].append(float(variance_gridobs))
    
    srt['MSE_Grid_All'] = srt.get('MSE_Grid_All',[])
    srt['Variance_Grid_All'] = srt.get('Variance_Grid_All',[])
    srt['MSE_Region_1'] = srt.get('MSE_Region_1',[])
    srt['Variance_Region_1'] = srt.get('Variance_Region_1',[])
    srt['MSE_Region_2'] = srt.get('MSE_Region_2',[])
    srt['Variance_Region_2'] = srt.get('Variance_Region_2',[])

    srt['MSE_Grid_All'].append(float(mse_gridall))
    srt['Variance_Grid_All'].append(float(variance_gridall))
    srt['MSE_Region_1'].append(float(mse_grid_reg1))
    srt['Variance_Region_1'].append(float(variance_grid_reg1))
    srt['MSE_Region_2'].append(float(mse_grid_reg2))
    srt['Variance_Region_2'].append(float(variance_grid_reg2)) 
    
    

        
   

#each variable type
for v in stats_dict:
    #each ensemble type
    for e in stats_dict[v]:
        #each stat type (MSE, var, etc)
        for S in stats_dict[v][e]:
            #convert all lists to np arrays
            stats_dict[v][e][S] = np.array(stats_dict[v][e][S])
        #find indices of increasing forecast hour for Gridded obs
        sorted_ind_gridded = stats_dict[v][e]['Forecast_Hour_Gridded'].argsort()
        #go through each stat and sort in order of increasing forecast hour
        for S in stats_dict[v][e]:
            stats_dict[v][e][S] = stats_dict[v][e][S][sorted_ind_gridded]

            
#plotting variables
#stats_list = ['MSE_Stationary','Variance_Stationary']
#get all stats types in a list, remove forecast hour
stats_list = list(stats_dict[v][e].keys())
remove_list = ['Forecast_Hour_Gridded']
for i in remove_list:
    stats_list.remove(i)

#now plot 
#separate plot for each variable type- remove the prior variable type
stats_dict_vars = list(stats_dict.keys())
#for v in stats_dict_vars:

def plot_stats():
    if control_vars == True:
        variables = plot_vars
    else:
        variables = stats_dict[m]
    #each variable
    for v in variables:
        #if we are plotting gridded obs stats
        plt.plot(stats_dict[m][v]['Forecast_Hour_Gridded'],stats_dict[m][v][s],
                 linestyle=ls2[v],marker=md[v],color=clr[m],label=v+' '+ens_dict[m])
        
        plt.xticks(np.arange(min(stats_dict[m][v]['Forecast_Hour_Gridded']), 
        max(stats_dict[m][v]['Forecast_Hour_Gridded'])+12, 12))
#            ax.set_xticks(numpy.arange(0, 1, 0.1))
#            ax.set_yticks(numpy.arange(0, 1., 0.1))
        plt.grid(True)
        plt.legend(loc = 'upper left')
        if v.startswith('QF850'):
            vstr = v[0:5]
        else:
            vstr = v[0:3]
        plt.title(title_dict[vstr]+s,fontsize=20)
        plt.xlabel('Forecast Hour',fontsize=14)
        plt.ylabel(var_units[vstr],fontsize=14)

#if AR_specific == False:
    #separate plot for MSE and variance
for s in stats_list:
    if separate_plots == False:
        fig = plt.figure(figsize=(14,8))  
        #each ensemble type      
        for m in stats_dict_vars:
            plot_stats()
        #plt.show()
        #fig.savefig(savedir+'850mb_Moisture_Flux_'+s+'_'+datestr+'.png',frameon=False,bbox_inches='tight')
       
    elif separate_plots == True:
        for m in stats_dict_vars:
            fig = plt.figure(figsize=(14,8))  
            plot_stats()
            #plt.show()
            #fig.savefig(savedir+'850mb_Moisture_Flux_'+s+'_'+datestr+'.png',frameon=False,bbox_inches='tight')
 
                
#elif AR_specific == True:
#    #separate plot for MSE and variance
#    for s in stats_list:
#        #each ensemble type      
#        for m in stats_dict_vars:
#            fig = plt.figure(figsize=(14,8))  
#            #each variable
#            for v in stats_dict[m]:
#                #if we are plotting gridded obs stats
#                plt.plot(stats_dict[m][v]['Forecast_Hour_Gridded'],stats_dict[m][v][s],
#                         linestyle='-',color=clr_sp[v],marker='o',label=v)
#                
#                plt.xticks(np.arange(min(stats_dict[m][v]['Forecast_Hour_Gridded']), 
#                max(stats_dict[m][v]['Forecast_Hour_Gridded'])+12, 12))
#    #            ax.set_xticks(numpy.arange(0, 1, 0.1))
#    #            ax.set_yticks(numpy.arange(0, 1., 0.1))
#                plt.grid(True)
#                plt.legend(loc = 'upper left')
#                plt.title(ens_dict[m]+' 850mb Moisture Flux '+s,fontsize=20)
#                plt.xlabel('Forecast Hour',fontsize=14)
#                plt.ylabel(var_units['QF850'],fontsize=14)
    
    #        fig.savefig(savedir+ens_dict[m]+'_850mb_Moisture_Flux_'+s+'_'+datestr+'.png',frameon=False,bbox_inches='tight')
    #        
        
##separate plot for each variable type
#for v in variables:
#    #separate plot for MSE and variance
#    for s in stats_list:
#        #separate plot for ensemble type
#        for m in stats_dict[v]:
#            fig = plt.figure(figsize=(14,8))
#            for l in stats_dict[v][m]:
#                plt.plot(stats_dict[v][m][l]['Forecast_Hour'],stats_dict[v][m][l][s],
#                         linestyle=ls[l],marker='o',color=clr[m],label=l+' '+m)
#                
#            plt.xticks(np.arange(min(stats_dict[var]['ecmwf']['prior']['Forecast_Hour'])-6, 
#                     max(stats_dict[var]['ecmwf']['prior']['Forecast_Hour'])+6, 6))
#            plt.grid()
#            plt.legend(loc = 'upper left')
#            plt.title(v+' '+s,fontsize=20)
#            plt.xlabel('Forecast Hour',fontsize=14)
#            plt.ylabel(var_units[v],fontsize=14)
#            
#            fig.savefig(savedir+v+'_'+s+'_'+m+'_'+datestr+'.png')
    
#    