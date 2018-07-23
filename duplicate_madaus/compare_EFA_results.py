#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:24:41 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date, date2num
import numpy as np
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time
import efa_functions as ef

ens = 'ncep'
loc_rad = '1000'
forecast_time = datetime(2015,11,12,0) # when the forecast was initialized
analysis_time = datetime(2015,11,13,0) # when to compare with analysis
oberrvar = [1,10,100,1000]
efh = '54hrs'
grid = [-180,180,90,0,3]
prior_var = ['QF850','D-QF850']

post_var = ['QF850']
vrbl= 'QF850'

time_ind = 1
point_lon = -99#163
point_lat = 32#116

ay = analysis_time.strftime('%Y')
am = analysis_time.strftime('%m')
ad = analysis_time.strftime('%d')
ah = analysis_time.strftime('%H')

fy = forecast_time.strftime('%Y')
fm = forecast_time.strftime('%m')
fd = forecast_time.strftime('%d')
fh = forecast_time.strftime('%H')

prior_var_str = ef.var_string(prior_var)
grid_str = ef.var_string(grid)



# Filepath of the analysis at the desired time 
analysis_path = '/home/disk/hot/stangen/Documents/prior_ensembles/'+ens+'/'+ay+am+'/'+ay+'-'+am+'-'+ad+'_'+ah+'_'+efh+'_'+prior_var_str+'.nc'             
# Filepath of the prior forecast
prior_path = '/home/disk/hot/stangen/Documents/prior_ensembles/'+ens+'/'+fy+fm+'/'+fy+'-'+fm+'-'+fd+'_'+fh+'_'+efh+'_'+prior_var_str+'.nc'

for oev in oberrvar:
    post_varstr = vrbl+str(oev)
    # Filepath of the posterior forecast
    post_path = '/home/disk/hot/stangen/Documents/posterior_ensembles/gridded/ob_update_all/inf_none/loc_'+loc_rad+'/'+ens+'/'+fy+fm+'/'+fy+'-'+fm+'-'+fd+'_'+fh+'_'+efh+'_'+grid_str+'_'+post_varstr+'.nc'
