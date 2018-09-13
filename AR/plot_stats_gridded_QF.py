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

#--------_Change these--------------------------------------------------------

#start and end dates of statistical averaging are used to load the correct file
start_date = datetime(2015,11,10,0)
end_date = datetime(2015,11,15,12)

#which ob type(s) were assimilated? used to load the correct .txt filename
#or, if desired, naming convention could be all types of obs the .txt file contains.
assim_obs = ['IVT']#['TCW']#

#what ob types are we looking at? could differ from assim_obs if
#we used assim_obs to update more than one variable (i.e. we didn't self-update)
#used if plotting change in statistic (if plot_change_stats = True)
plot_ob = 'IVT'

#grid used to generate gridded obs, used for loading correct .txt filename
grid = [-180,180,90,0,3]

#want to plot each variable or control which are plotted?
control_vars = True


#if control_vars = True, plot_vars controls which variables we are going to plot. This can be useful 
#if there are lots of different variables in the txt file and we don't want to plot them all.
#plot_vars = ['QF850_prior','QF8501_loc1000','QF85010_loc1000','QF850100_loc1000','QF850250_loc1000','QF850500_loc1000','QF850750_loc1000','QF8501000_loc1000']
#plot_vars = ['TCW_prior','TCW0-1_loc1000','TCW1_loc1000','TCW10_loc1000','TCW100_loc1000']
plot_vars = ['IVT_prior','IVT100_loc1000','IVT1000_loc1000','IVT5000_loc1000','IVT10000_loc1000','IVT20000_loc1000']
#plot_vars = ['IWV_prior','IWV1_loc1000','IWV5_loc1000','IWV10_loc1000','IWV20_loc1000','IWV100_loc1000']
#plot_vars = ['IVT_prior','IVT10000_loc1000','IVT10000_loc2000hybrid','IVT10000_loc5000hybrid','IVT10000_loc10000hybrid','IVT10000_loc10000']
#plot_vars = ['IWV_prior','IWV20_loc1000','IWV20_loc2000hybrid','IWV20_loc5000hybrid','IWV20_loc10000hybrid','IWV20_loc10000']
#plot_vars = ['IVT_prior','IVT10000_loc1000','IVT10000_loc10000']
#plot_vars = ['IWV_prior', 'IWV20_loc1000','IWV20_loc10000']
#plot_vars = ['IVT_prior','IVT10000_loc99statsig','IVT10000_loc98statsig','IVT10000_loc95statsig','IVT10000_loc90statsig']
#plot_vars = ['IWV_prior','IWV20_loc99statsig','IWV20_loc98statsig','IWV20_loc95statsig','IWV20_loc90statsig']
#plot_vars = ['IVT_prior','IVT1000_loc1000','IVT1000_loc10000']

#are we wanting to look at statistics for MSE/variance within a specific AR?
AR_specific = True

#do we want a separate plot for each ensemble type?
separate_plots = False

#plot actual statistics, or change in statistics compared to prior?
plot_change_stats = False

#gridded or madis obs?

#-------------------------------------------------------------------------------

var_units = {
            'ALT' : 'hPa$^{2}$',
            'T2M' : 'K$^{2}$',
            'QF850' : '(g/kg*m/s)$^{2}$',
            'TCW' : 'mm$^{2}$',
            'IWV' : 'mm$^{2}$',
            'IVT' : '(kg/m/s)$^{2}$'
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
        'eccc' : 'k',
        'ncep' : 'b'
        }

plot_dict = {
       'QF850' : {'prior' : {'ls' : '-', 'mkr' : 'o', 'clr' : 'k'}},#solid
       'QF850ensvar' : {'loc1000' : {'ls' : '--', 'mkr' : 'o', 'clr' : 'b'}}, #dashed
       'QF8501000' : {'loc1000' : {'ls' : ':', 'mkr' : 'o', 'clr' : 'r'}}, #dotted
       'QF850100' : {'loc1000' : {'ls' : '-.', 'mkr' : 'o', 'clr' : 'g'}}, #dashdotted
       'QF85010' : {'loc1000' : {'ls' : ' ', 'mkr' : 'x', 'clr' : 'c'}},
       'QF8501' : {'loc1000' : {'ls' : '', 'mkr' : '*', 'clr' : 'y'}},
       'QF850250' : {'loc1000' : {'ls' : '--', 'mkr' : 'x', 'clr' : 'm'}},
       'QF850500' : {'loc1000' : {'ls' : '--', 'mkr' : '*', 'clr' : 'y'}},
       'QF850750' : {'loc1000' : {'ls' : '', 'mkr' : 's', 'clr' : 'c'}},
       
       'TCW' : {'prior' : {'ls' : '-', 'mkr' : 'o', 'clr' : 'k'}},
       'TCW0-1' : {'loc1000' : {'ls' : ' ', 'mkr' : 'o', 'clr' : 'y'}},
       'TCW1' : {'loc1000' : {'ls' : '--', 'mkr' : 'o', 'clr' : 'g'}},
       'TCW10' : {'loc1000' : {'ls' : '-.', 'mkr' : 'o', 'clr' : 'b'}},
       'TCW100' : {'loc1000' : {'ls' : ':', 'mkr' : 'o', 'clr' : 'r'}},
       
       'IWV' : {'prior' : {'ls' : '-', 'mkr' : 'o', 'clr' : 'k'}},
       'IWV1' : {'loc1000' : {'ls' : '', 'mkr' : 'o', 'clr' : 'y'}},
       'IWV5' : {'loc1000' : {'ls' : ':', 'mkr' : 'o', 'clr' : 'r'}},
       'IWV10' : {'loc1000' : {'ls' : '-.', 'mkr' : 'o', 'clr' : 'g'}},
       'IWV20' : {'loc1000' : {'ls' : '--', 'mkr' : 'o', 'clr' : 'c'},
                  'loc2000hybrid': {'ls' : ':', 'mkr' : 'o', 'clr' : 'b'},
                  'loc5000hybrid': {'ls' : '-.', 'mkr' : 'o', 'clr' : 'c'},
                  'loc10000hybrid': {'ls' : ':', 'mkr' : 'x', 'clr' : 'r'},
                  'loc99statsig' : {'ls' : ':', 'mkr' : 'o', 'clr' : 'b'},
                  'loc98statsig' : {'ls' : '--', 'mkr' : 'x', 'clr' : 'c'},
                  'loc95statsig' : {'ls' : '-.', 'mkr' : 'o', 'clr' : 'r'},
                  'loc90statsig' : {'ls' : ':', 'mkr' : 'x', 'clr' : 'g'},
                  'loc10000' : {'ls' : '-.', 'mkr' : 'x', 'clr' : 'y'}},
       'IWV100' : {'loc1000' : {'ls' : ':', 'mkr' : 'x', 'clr' : 'b'}},
       'IWV1000' : {'loc1000' : {'ls' : '', 'mkr' : 'o', 'clr' : 'r'}},
       'IWV5000' : {'loc1000' : {'ls' : '-.', 'mkr' : 'o', 'clr' : 'y'}},
       'IWV10000' : {'loc1000' : {'ls' : '--', 'mkr' : 'o', 'clr' : 'g'}},
       'IWV20000' : {'loc1000' : {'ls' : ':', 'mkr' : 'o', 'clr' : 'c'}},
       
       'IVT' : {'prior' : {'ls' : '-', 'mkr' : 'o', 'clr' : 'k', 'lw' : 2.5}},
       'IVT1' : {'loc1000' : {'ls' : '', 'mkr' : 'o', 'clr' : 'y'}},
       'IVT5' : {'loc1000' : {'ls' : ':', 'mkr' : 'o', 'clr' : 'r'}},
       'IVT10' : {'loc1000' : {'ls' : '--', 'mkr' : 'o', 'clr' : 'g'}},
       'IVT20' : {'loc1000' : {'ls' : '-.', 'mkr' : 'o', 'clr' : 'c'}},
       'IVT100' : {'loc1000' : {'ls' : '-', 'mkr' : '*', 'clr' : 'b', 'lw' : .75}},
       'IVT1000' : {'loc1000' : {'ls' : '-', 'mkr' : 'd', 'clr' : 'r', 'lw' : 1.5},
                    'loc10000hybrid' : {'ls' : ':', 'mkr' : 'o', 'clr' : 'r', 'lw' : 1.5},
                    'loc10000' : {'ls' : '-.', 'mkr' : 'x', 'clr' : 'y', 'lw' : 1.5}},
       'IVT5000' : {'loc1000' : {'ls' : '-.', 'mkr' : 'x', 'clr' : 'y', 'lw' : 1.5}},
       'IVT10000' : {'loc1000' : {'ls' : '--', 'mkr' : 'v', 'clr' : 'g', 'lw' : 1.5},
                  'loc2000hybrid': {'ls' : ':', 'mkr' : 'o', 'clr' : 'b'},
                  'loc5000hybrid': {'ls' : '-.', 'mkr' : 'o', 'clr' : 'c'},
                  'loc10000hybrid': {'ls' : ':', 'mkr' : 'x', 'clr' : 'r'},
                  'loc99statsig' : {'ls' : ':', 'mkr' : 'o', 'clr' : 'b'},
                  'loc98statsig' : {'ls' : '--', 'mkr' : 'x', 'clr' : 'c'},
                  'loc95statsig' : {'ls' : '-.', 'mkr' : 'o', 'clr' : 'r'},
                  'loc90statsig' : {'ls' : ':', 'mkr' : 'x', 'clr' : 'g'},
                  'loc10000' : {'ls' : '-.', 'mkr' : 'x', 'clr' : 'y'}},
       'IVT20000' : {'loc1000' : {'ls' : ':', 'mkr' : '+', 'clr' : 'c', 'lw' : 1.5}}
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

if AR_specific == True:
    filepath = filedir+datestr+'_'+varstr+'_'+gridstr+'_gridobs_ARspecific.txt'
else:
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
    
    #if we are testing a method of localization other than gaspari-cohn, make the ob variable type
    #be the localization radius.
    #if efa.endswith('hybrid') or efa.endswith('statsig') or efa.endswith('statsig2'):
    
    #combine the variable type and the localization radius/type together to make a descriptive dictionary key
    ens = ens+'_'+efa
    
    fh = line_split[4] #forecast hour
    mse_gridobs = line_split[5] #MSE of gridded obs if AR_specific == False, or MSE of AR gridpoints if AR_specific == True
    variance_gridobs = line_split[6] #variance of gridded obs if AR_specific == False, or variance of AR gridpoints if if AR_specific == True
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
#get all stats types in a list, remove forecast hour
stats_list = list(stats_dict[v][e].keys())
remove_list = ['Forecast_Hour_Gridded']
for i in remove_list:
    stats_list.remove(i)

#the different ensemble types in the .txt file (eccc, ecmwf, ncep) 
stats_dict_vars = list(stats_dict.keys())

#now plot
def plot_stats(separate):
    """
    takes in a boolean argument which affects which color table to use
    """
    if control_vars == True:
        variables = plot_vars
    else:
        variables = stats_dict[m]
    #each variable
    for v in variables:
        #for plotting
        varoev, loc = v.split('_')
        #for plotting title and ylabel
        if v.startswith('QF850'):
            vstr = v[0:5]
        else:
            vstr = v[0:3]
        #if we want to plot the change in statistics instead of actual statistics
        if plot_change_stats == True:
            compare_var = plot_ob + '_prior' 
            plot_stats = (stats_dict[m][v][s] - stats_dict[m][compare_var][s])#/stats_dict[m][v][s]
            title_str = 'Change vs Prior in '+title_dict[vstr]+s
        elif plot_change_stats == False:
            plot_stats = stats_dict[m][v][s]
            title_str = title_dict[vstr]+s
        #if we are plotting gridded obs stats
        if separate == False:
            plt.plot(stats_dict[m][v]['Forecast_Hour_Gridded'],plot_stats,
                     linestyle=plot_dict[varoev][loc]['ls'],marker=plot_dict[varoev][loc]['mkr'],color=clr[m],linewidth=plot_dict[varoev][loc]['lw'],label=v+' '+ens_dict[m])
        elif separate == True:
            plt.plot(stats_dict[m][v]['Forecast_Hour_Gridded'],plot_stats,
                     linestyle=plot_dict[varoev][loc]['ls'],marker=plot_dict[varoev][loc]['mkr'],color=plot_dict[varoev][loc]['clr'],linewidth=plot_dict[varoev][loc]['lw'],label=v+' '+ens_dict[m])            
        
        plt.xticks(np.arange(min(stats_dict[m][v]['Forecast_Hour_Gridded']), 
        max(stats_dict[m][v]['Forecast_Hour_Gridded'])+12, 12))
#            ax.set_xticks(numpy.arange(0, 1, 0.1))
#            ax.set_yticks(numpy.arange(0, 1., 0.1))
        plt.grid(True)
        plt.legend(loc = 'upper left')
        plt.title(title_str,fontsize=20)
        plt.xlabel('Forecast Hour',fontsize=14)
        plt.ylabel(var_units[vstr],fontsize=14)


for s in stats_list:
    if separate_plots == False:
        fig = plt.figure(figsize=(14,8))  
        #each ensemble type      
        for m in stats_dict_vars: #['eccc','ecmwf']:['eccc']:
            plot_stats(False)
            #plt.show()
            #fig.savefig(savedir+'850mb_Moisture_Flux_'+s+'_'+datestr+'.png',frameon=False,bbox_inches='tight')
       
    elif separate_plots == True:
        for m in stats_dict_vars:
            fig = plt.figure(figsize=(14,8))  
            plot_stats(True)
            #plt.show()
            #fig.savefig(savedir+'850mb_Moisture_Flux_'+s+'_'+datestr+'.png',frameon=False,bbox_inches='tight')
 
                
#elif AR_specific == True:
    ##separate plot for MSE and variance
    #for s in stats_list:
    #    #each ensemble type      
    #    for m in stats_dict_vars:
    #        fig = plt.figure(figsize=(14,8))  
    #        #each variable
    #        for v in stats_dict[m]:
    #            #if we are plotting gridded obs stats
    #            plt.plot(stats_dict[m][v]['Forecast_Hour_Gridded'],stats_dict[m][v][s],
    #                     linestyle='-',color=clr_sp[v],marker='o',label=v)
    #            
    #            plt.xticks(np.arange(min(stats_dict[m][v]['Forecast_Hour_Gridded']), 
    #            max(stats_dict[m][v]['Forecast_Hour_Gridded'])+12, 12))
    ##            ax.set_xticks(numpy.arange(0, 1, 0.1))
    ##            ax.set_yticks(numpy.arange(0, 1., 0.1))
    #            plt.grid(True)
    #            plt.legend(loc = 'upper left')
    #            plt.title(ens_dict[m]+' 850mb Moisture Flux '+s,fontsize=20)
    #            plt.xlabel('Forecast Hour',fontsize=14)
    #            plt.ylabel(var_units['QF850'],fontsize=14)
    
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