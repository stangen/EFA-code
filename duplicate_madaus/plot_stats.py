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
end_date = datetime(2013,4,1,12)

variables = ['T2M','ALT']

var_units = {
            'ALT' : 'hPa$^{2}$',
            'T2M' : 'K$^{2}$'       
            }

filedir = '/home/disk/hot/stangen/Documents/EFA/duplicate_madaus/mse_var_output/'

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

filepath = filedir+datestr+'_'+varstr+'.txt'

f1 = open(filepath, 'r')
stats = f1.readlines()

stats_dict = {}
#build the dictionary from each entry in the statistics file
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
    stats_dict[var][ens] = stats_dict.get(ens, {})
    stats_dict[var][ens][efa] = stats_dict.get(efa, {})
    stats_dict[var][ens][efa]['Forecast_Hour'] = stats_dict[var][ens][efa].get('Forecast_Hour',[])
    stats_dict[var][ens][efa]['MSE'] = stats_dict[var][ens][efa].get('MSE',[])
    stats_dict[var][ens][efa]['Ensemble_Variance'] = stats_dict[var][ens][efa].get('Ensemble_Variance',[])
    stats_dict[var][ens][efa]['Error_Variance'] = stats_dict[var][ens][efa].get('Error_Variance',[])
#now fill the dictionary with the statistics
for line in stats:
    line_split = line.split(',')
    efa = line_split[0] #posterior or prior
    ens = line_split[1] #ensemble type
    var = line_split[2] #observation variable type
    fh = line_split[3] #forecast hour
    mse = line_split[4] #mean squared error
    variance = line_split[5] #variance
    errvar = line_split[6] #error variance
    
    stats_dict[var][ens][efa]['Forecast_Hour'].append(int(fh))
    stats_dict[var][ens][efa]['MSE'].append(float(mse))
    stats_dict[var][ens][efa]['Ensemble_Variance'].append(float(variance))
    stats_dict[var][ens][efa]['Error_Variance'].append(float(errvar))
    
    
    