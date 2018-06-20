#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 16:19:04 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date
from datetime import timedelta
import numpy as np
#from old_ensemble_verification import error_vs_spread
from EFA.efa_xray.state.ensemble import EnsembleState

class Load_data():
    """
    This class exists mainly because the two functions share some variables,
    and this helps reduce the amount of redundant inputs into each function. 
    """
    
    def __init__(self, date,ens_type,vrbls,ob_type):
        self.date = date  
        self.y = self.date.strftime('%Y')
        self.m = self.date.strftime('%m')
        self.d = self.date.strftime('%d')
        self.h = self.date.strftime('%H')
        self.ens_type = ens_type
        self.vrbls = vrbls
        self.ob_type = ob_type
        


    
    
    
    def load_netcdfs(self):
        """
        Loads the ensemble netCDF and the elevation netCDF. Packages the 
        ensemble data into an EnsembleState object for use in Luke's code.
        Returns this, lats, lons, and elevations
        """
        #block of code to access netCDF files, according to naming convention.
        #deals with names of variables, since names of variables are included
        #in the file name (for now).
        var_string = ''
        for v in self.vrbls[:-1]:  
            var_string = var_string+v+'_'
        var_string = var_string+self.vrbls[-1]
        
        # directory where the ensemble of all times is
        infile = '/home/disk/hot/stangen/Documents/prior_ensembles/'+self.ens_type+'/'+self.y+self.m+'/'+self.y+'-'+self.m+'-'+self.d+'_'+self.h+'_'+self.ens_type+'_'+var_string+'.nc' 
        
        print('loading netcdf file: '+self.y+self.m+self.d+'_'+self.h+'00')
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
        
        # storing the variable data in a dict (state?)
        allvars = {}
        for var in self.vrbls:
            allvars[var] = (['validtime','y','x','mem'],
                            ncdata.variables[var][:])
        lonarr, latarr = np.meshgrid(lons, lats)
        
        # Package into an EnsembleState object knowing the state and metadata
        statecls = EnsembleState.from_vardict(allvars,
                                      {'validtime' : ftimes,
                                       'lat' : (['y','x'], latarr),
                                       'lon' : (['y','x'], lonarr),
                                       'mem' : mems,
                                       })
        
        
        # directory where the orography file is
        orography = '/home/disk/hot/stangen/Documents/tigge_ensembles/orography/2013-04-01_00_'+self.ens_type+'.nc'
        orog_data = Dataset(orography,"r")
        #print(of.variables)
        elevs = orog_data.variables['orog'][0,:] # Elevation of ecmwf or eccc         
        
        return statecls, lats, lons, elevs

        
    def load_obs(self, forecast_hour=6):
        """
        Loads the observations corresponding with n hours after ensemble was
        initialized (default 6 hours). Returns the observations from the text file.
        """
        
        #get observations from n hours after the model was initialized
        dt0 = self.date
        dt = dt0.replace(minute = 0, second=0, microsecond=0)
        dt = dt + timedelta(0,3600*forecast_hour)
        
        #convert back to strings
        dty = dt.strftime('%Y')
        dtm = dt.strftime('%m')
        dtd = dt.strftime('%d')
        dth = dt.strftime('%H')
        
        # directory where the observations are
        obs_file = '/home/disk/hot/stangen/Documents/surface_obs/MADIS/'+dty+dtm+'/combined_'+self.ob_type+'/'+self.ob_type+'_'+dty+dtm+dtd+'_'+dth+'00.txt'
        print('loading '+self.ob_type+ ' from '+dty+dtm+dtd+'_'+dth+'00')
        #print(obs_file)
        f1 = open(obs_file, 'r')
        obs = f1.readlines()
        
        return obs
    
