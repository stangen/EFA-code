#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 14:31:58 2018

@author: stangen

Reading the .txt file this script generates:
    First column: observation number
    Second column: latitude
    Third column: longitude
    Fourth column: Elevation, set arbitrarily to 0, since we aren't doing a 
        terrain check on the gridded obs, and elevation isn't important.
    Fifth column: Time of observation in epoch time (seconds since Jan 1, 1970)
    Sixth column: Observation value
    Seventh column: Observation type (GRIDDED in this case)
    Eighth column: Marker for stationary ob, 0 used arbitrarily since we 
    aren't looking at whether or not an observation is stationary
    Ninth column: if get_variance == True, then there will be a 9th column 
    containing variance.
    
"""
from datetime import datetime

from EFA.duplicate_madaus.load_data import Load_Data
import surface_obs.madis_example.madis_utilities as mt 

#---------------Change these---------------------------------------------------
ens = ['ncep']#,'eccc','ecmwf']
vrbls = ['IWV','IVT','D-IVT']#['TCW']#['QF850','D-QF850']#['T2M','ALT']
ob_type = 'IVT'#'QF850'#'ALT'
start_date = datetime(2015,11,10,0)
end_date = datetime(2015,11,17,12)
get_variance=False #if true, will save ensemble variance as well
#------------------------------------------------------------------------------

dates = mt.make_datetimelist(start_date,end_date,12)

for e in ens:
    for date in dates:
        efa = Load_Data(date,e,vrbls,ob_type,grid=[-180,180,90,0,3],
                        new_format=True,efh=54)
        
        efa.save_gridded_obs(get_variance=get_variance)

