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

#--------Change these----------------------------------------------------------

start_date = datetime(2013,4,1,0)
end_date = datetime(2013,4,30,12)

#which variables are in the .txt file?
variables = ['ALT']

#did we assimilate MADIS or gridded observations? 'madis' or 'gridded'
ob_category = 'gridded'

#-----------------------------------------------------------------------------

#prior_var_key = 'ALT'

var_units = {
            'ALT' : 'hPa$^{2}$',
            'T2M' : 'K$^{2}$'       
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


filepath = filedir+datestr+'_'+varstr
if ob_category == 'gridded':
    filepath += '_gridobs.txt'
elif ob_category == 'madis':
    filepath += '.txt'

f1 = open(filepath, 'r')
stats = f1.readlines()

stats_dict = {}

for line in stats:
    line_split = line.split(',')
    inf = line_split[0] #inflation
    efa = line_split[1] #posterior or prior
    ens = line_split[2] #ensemble type
    var = line_split[3] #observation variable type
    fh = line_split[4] #forecast hour
    mse_stationary = line_split[5] #mean squared error of stationary MADIS obs
    variance_stationary = line_split[6] #variance of stationary MADIS obs
    mse_allmaritime = line_split[7] #mean squared error of all MADIS obs
    variance_allmaritime = line_split[8] #variance of stationary MADIS obs

    
    #this block is to set up the dictionary structure
    stats_dict[var] = stats_dict.get(var,{})
    stats_dict[var][ens] = stats_dict[var].get(ens, {})
    stats_dict[var][ens][efa] = stats_dict[var][ens].get(efa, {})
    srt = stats_dict[var][ens][efa]
    srt['Forecast_Hour_Madis'] = srt.get('Forecast_Hour_Madis',[])
    srt['MSE_Stationary'] = srt.get('MSE_Stationary',[])
    srt['Variance_Stationary'] = srt.get('Variance_Stationary',[])
    #new stuff
    srt['MSE_All_Maritime'] = srt.get('MSE_All_Maritime',[])
    srt['Variance_All_Maritime'] = srt.get('Variance_All_Maritime',[])
    
    #add the data to the dictionary
    srt['Forecast_Hour_Madis'].append(int(fh))
    srt['MSE_Stationary'].append(float(mse_stationary))
    srt['Variance_Stationary'].append(float(variance_stationary))
    srt['MSE_All_Maritime'].append(float(mse_allmaritime))
    srt['Variance_All_Maritime'].append(float(variance_allmaritime))

    
    try:
        mse_gridobs = line_split[9] #mean squared error of gridded obs
        variance_gridobs = line_split[10] #variance of gridded obs
        mse_gridall = line_split[11] #mean squared error of entire grid
        variance_gridall = line_split[12] #variance of entire grid
        srt['Forecast_Hour_Gridded'] = srt.get('Forecast_Hour_Gridded',[])
        srt['MSE_Grid_Obs'] = srt.get('MSE_Grid_Obs',[])
        srt['Variance_Grid_Obs'] = srt.get('Variance_Grid_Obs',[])
        srt['MSE_Grid_All'] = srt.get('MSE_Grid_All',[])
        srt['Variance_Grid_All'] = srt.get('Variance_Grid_All',[])
        srt['Forecast_Hour_Gridded'].append(int(fh)) #forecast hours of just these ones
        srt['MSE_Grid_Obs'].append(float(mse_gridobs))
        srt['Variance_Grid_Obs'].append(float(variance_gridobs))
        srt['MSE_Grid_All'].append(float(mse_gridall))
        srt['Variance_Grid_All'].append(float(variance_gridall))
    except:
        pass

#each variable type
for v in stats_dict:
    #each ensemble type
    for e in stats_dict[v]:
        #each localization radius
        for E in stats_dict[v][e]:
            #each stat type (MSE, var, etc)
            for S in stats_dict[v][e][E]:
                #convert all lists to np arrays
                stats_dict[v][e][E][S] = np.array(stats_dict[v][e][E][S])
            #find indices of increasing forecast hour for MADIS obs
            sorted_ind_madis = stats_dict[v][e][E]['Forecast_Hour_Madis'].argsort()
            #find indices of increasing forecast hour for Gridded obs
            sorted_ind_gridded = stats_dict[v][e][E]['Forecast_Hour_Gridded'].argsort()
            #go through each stat and sort in order of increasing forecast hour
            for S in stats_dict[v][e][E]:
                #if this is a MADIS stat
                if len(stats_dict[v][e][E][S]) == len(stats_dict[v][e][E]['Forecast_Hour_Madis']):
                    stats_dict[v][e][E][S] = stats_dict[v][e][E][S][sorted_ind_madis]      
                elif len(stats_dict[v][e][E][S]) == len(stats_dict[v][e][E]['Forecast_Hour_Gridded']):
                    stats_dict[v][e][E][S] = stats_dict[v][e][E][S][sorted_ind_gridded]

            
#plotting variables
#stats_list = ['MSE_Stationary','Variance_Stationary']
#get all stats types in a list, remove forecast hour
stats_list = list(stats_dict[v][e][E].keys())
remove_list = ['Forecast_Hour_Madis','Forecast_Hour_Gridded']
for i in remove_list:
    stats_list.remove(i)

#now plot 
#separate plot for each variable type- remove the prior variable type, if
#the prior variable type name differs from other variable type names (i.e.,
#if we added observation error variance to the variable name in the .txt file)
stats_dict_vars = list(stats_dict.keys())
for j in variables:
    count=0
    for k in stats_dict_vars:
        if k.startswith(j):
            count +=1
    if count > 1:
        stats_dict_vars.remove(j)

for v in stats_dict_vars:
#for v in variables:
    #separate plot for MSE and variance
    for s in stats_list:
        fig = plt.figure(figsize=(14,8))
        #line for each ensemble type
        for m in stats_dict[v]:
            #line for prior/localization radius
            for l in stats_dict[v][m]:
                #if we are plotting MADIS obs stats
                if len(stats_dict[v][m][l][s]) == len(stats_dict[v][m][l]['Forecast_Hour_Madis']):
                    plt.plot(stats_dict[v][m][l]['Forecast_Hour_Madis'],stats_dict[v][m][l][s],
                             linestyle=ls[l],marker='o',color=clr[m],label=l+' '+m)
                #if we are plotting gridded obs stats
                elif len(stats_dict[v][m][l][s]) == len(stats_dict[v][m][l]['Forecast_Hour_Gridded']):
                    plt.plot(stats_dict[v][m][l]['Forecast_Hour_Gridded'],stats_dict[v][m][l][s],
                             linestyle=ls[l],marker='o',color=clr[m],label=l+' '+m)
                
            #add the prior for each ens type to each variable and stat_type plot- only
            #necessary if variable name is different for prior (if we added ob err var to variable name in .txt file)
            if count > 1:
                for k in variables:
                    #match the prior variable name (ALT) with the assimilated variable name (ALT1)
                    if v[:len(k)] == k:
                        z = v[:len(k)]
                        if len(stats_dict[z][m]['prior'][s]) == len(stats_dict[z][m]['prior']['Forecast_Hour_Madis']):
                            plt.plot(stats_dict[z][m]['prior']['Forecast_Hour_Madis'],stats_dict[z][m]['prior'][s],
                                     linestyle=ls['prior'],marker='o',color=clr[m],label='prior '+m)
                            
                        elif len(stats_dict[z][m]['prior'][s]) == len(stats_dict[z][m]['prior']['Forecast_Hour_Gridded']):
                            plt.plot(stats_dict[z][m]['prior']['Forecast_Hour_Gridded'],stats_dict[z][m]['prior'][s],
                                     linestyle=ls['prior'],marker='o',color=clr[m],label='prior '+m)
            
#        plt.plot(stats_dict[v]['ecmwf']['prior']['Forecast_Hour'],stats_dict[v]['ecmwf']['prior'][s],
#                 linestyle='dashed',marker='o',color='r',label='Original ECMWF')
#        plt.plot(stats_dict[v]['ecmwf']['posterior']['Forecast_Hour'],stats_dict[v]['ecmwf']['posterior'][s],
#                 marker='o',color='r',label='Updated ECMWF')
#        plt.plot(stats_dict[v]['eccc']['prior']['Forecast_Hour'],stats_dict[v]['eccc']['prior'][s],
#                 linestyle='dashed',marker='o',color='k',label='Original CMC')
#        plt.plot(stats_dict[v]['eccc']['posterior']['Forecast_Hour'],stats_dict[v]['eccc']['posterior'][s],
#                 marker='o',color='k',label='Updated CMC')   
        
#        plt.xticks(np.arange(min(stats_dict[var]['eccc']['prior']['Forecast_Hour'])-6, 
#                             max(stats_dict[var]['eccc']['prior']['Forecast_Hour'])+6, 6))
        plt.grid()
        plt.legend(loc = 'upper left')
        plt.title(v+' '+s,fontsize=20)
        plt.xlabel('Forecast Hour',fontsize=14)
        #plt.ylabel(var_units[v],fontsize=14)
        
        #fig.savefig(savedir+v+'_allobs_'+s+'_'+datestr+'.png',frameon=False,bbox_inches='tight')
        
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