#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 14:31:58 2018

@author: stangen
"""
from datetime import datetime

from EFA.duplicate_madaus.load_data import Load_Data
import surface_obs.madis_example.madis_utilities as mt 
#file to call the functions to generate/save gridded obs
ens = 'ecmwf'
vrbls = ['T2M','ALT']
ob_type = 'ALT'
update_vars = ['ALT']
start_date = datetime(2013,4,1,0)
end_date = datetime(2013,4,1,0)

dates = mt.make_datetimelist(start_date,end_date,12)

for date in dates:
    efa = Load_Data(date,ens,vrbls,ob_type,update_vars)
    
    efa.save_gridded_obs()


#remember to add sum to the multiplication in the spaceweights to get 1 number