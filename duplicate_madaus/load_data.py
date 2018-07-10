#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 16:19:04 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date
from datetime import timedelta
import numpy as np
import efa_functions as ef
#from old_ensemble_verification import error_vs_spread
from EFA.efa_xray.state.ensemble import EnsembleState

class Load_Data():
    """
    This class exists mainly because the two functions share some variables,
    and this helps reduce the amount of redundant inputs into each function. 
    #ens_type, ob_type are str; vrbls, update_var are lists of str, date is datetime object.
    date = date of the model initialization
    ens_type = eccc, ecmwf
    vrbls = all the variables in the netCDF
    ob_type = the observation variable type we will be assimilating (only 1)
    update_var = the variable type(s) we will be updating
    """
    
    def __init__(self,date,ens_type,vrbls,ob_type,update_var):
        self.date = date  
        self.y = self.date.strftime('%Y')
        self.m = self.date.strftime('%m')
        self.d = self.date.strftime('%d')
        self.h = self.date.strftime('%H')
        self.ens_type = ens_type
        self.var_string = ef.var_string(vrbls) #convert list of vrbls to string
        #self.vrbls = vrbls
        self.ob_type = ob_type
        self.update_var = update_var


    
    
    
    def load_netcdfs(self,post=False,ob_cat='madis',ob_upd='ob_update_self',inf='none',lr='1000'):
        """
        Loads the ensemble netCDF and the elevation netCDF. 
        Allows an option to package the ensemble data into an EnsembleState 
        object for use in Luke's code. If posterior ensemble, there are more 
        options for exactly what EFA took place for finding the right file.      
        """
        
        # directory where the ensemble of all times is
        if post==False:
            infile = '/home/disk/hot/stangen/Documents/prior_ensembles/'+self.ens_type+'/'+self.y+self.m+'/'+self.y+'-'+self.m+'-'+self.d+'_'+self.h+'_'+self.ens_type+'_'+self.var_string+'.nc' 
            prior_or_post='prior'
        elif post==True:
            infile = '/home/disk/hot/stangen/Documents/posterior_ensembles/'+ob_cat+'/'+ob_upd+('/inf_'+inf).replace('.','-')+'/loc_'+str(lr)+'/'+self.ens_type+'/'+self.y+self.m+'/'+self.y+'-'+self.m+'-'+self.d+'_'+self.h+'_'+self.ens_type+'_'+self.var_string+'.nc' 
            prior_or_post='posterior'
        print('loading netcdf file: '+prior_or_post+': '+self.ens_type+' '+self.y+self.m+self.d+'_'+self.h+'00')
        # loading/accessing the netcdf data            
        ncdata = Dataset(infile,'r')
        #print(ncdata.variables.keys())
        times = ncdata.variables['time']
        ftimes = num2date(times[:],
                          times.units)
        lats = ncdata.variables['lat'][:]
        lons = ncdata.variables['lon'][:]
        mems = ncdata.variables['ens'][:]
        #print(ncdata.variables)
        
        # directory where the orography file is
        orography = '/home/disk/hot/stangen/Documents/tigge_ensembles/orography/2013-04-01_00_'+self.ens_type+'.nc'
        orog_data = Dataset(orography,"r")
        #print(of.variables)
        elevs = orog_data.variables['orog'][0,:] # Elevation of ecmwf or eccc  
        
        # storing the variable data in a dict (state?)
        allvars = {}
        for var in self.update_var:
            allvars[var] = (['validtime','y','x','mem'],
                            ncdata.variables[var][:])
        lonarr, latarr = np.meshgrid(lons, lats)
        
        pack_str = ef.var_string(self.update_var)
        print('packaging '+pack_str+' into EnsembleState object')
        # Package into an EnsembleState object knowing the state and metadata
        statecls = EnsembleState.from_vardict(allvars,
                                      {'validtime' : ftimes,
                                       'lat' : (['y','x'], latarr),
                                       'lon' : (['y','x'], lonarr),
                                       'mem' : mems,
                                       })
        
        return statecls, lats, lons, elevs
             
    def load_obs(self, forecast_hour=6):
        """
        Loads the observations corresponding with n hours after ensemble was
        initialized (default 6 hours). Returns the observations from the text file.
        """
        
        #get observations from n hours after the model was initialized
        dt0 = self.date
        dt = dt0.replace(minute = 0, second=0, microsecond=0)
        dt = dt + timedelta(hours=forecast_hour)
        
        #convert back to strings
        dty = dt.strftime('%Y')
        dtm = dt.strftime('%m')
        dtd = dt.strftime('%d')
        dth = dt.strftime('%H')
        
        ob_str = ef.var_string([self.ob_type])        
        # directory where the observations are
        obs_file = '/home/disk/hot/stangen/Documents/surface_obs/MADIS/'+dty+dtm+'/combined_'+self.ob_type+'/'+self.ob_type+'_'+dty+dtm+dtd+'_'+dth+'00.txt'
        print('loading '+ob_str+' obs from '+dty+dtm+dtd+'_'+dth+'00')
        #print(obs_file)
        f1 = open(obs_file, 'r')
        obs = f1.readlines()
        
        return obs
    
    def load_gridded_obs(self, forecast_hour=12):
        """
        Loads the ensemble from 12 hours after the forecast was initialized
        and uses its 0-hour forecast as the "observation" grid. 
        """
        
        #get the list of points
        points = ef.get_ob_points()
        #get the 0-hour ensemble forecast, n hours after the model was initialized
        dt0 = self.date
        dt = dt0.replace(minute = 0, second=0, microsecond=0)
        dt = dt + timedelta(hours=forecast_hour)
        
        #convert back to strings
        dty = dt.strftime('%Y')
        dtm = dt.strftime('%m')
        dtd = dt.strftime('%d')
        dth = dt.strftime('%H')
        
        infile = '/home/disk/hot/stangen/Documents/prior_ensembles/'+self.ens_type+'/'+dty+dtm+'/'+dty+'-'+dtm+'-'+dtd+'_'+dth+'_'+self.ens_type+'_'+self.var_string+'.nc'
        print('loading netcdf gridded '+self.ob_type+ ' "obs": '+self.ens_type+' '+dty+dtm+dtd+'_'+dth+'00')
        
        ncdata = Dataset(infile,'r')
        #print(ncdata.variables.keys())
        #gridded ensemble-mean 0-hour forecasts n (12) hours after the ensemble
        #we are updating was initialized
        var = ncdata.variables[self.ob_type][0,:,:,:].mean(axis=-1)

        lats = ncdata.variables['lat'][:]
        lons = ncdata.variables['lon'][:]
        
        #change efa functions so that it takes an argument for need_elev, 
        #then create a separate interpolate function for this scenario: no
        #interpolation in time, no ensemble members (just ens mean), just use
        #4 closest grid points.
        
        
    
