#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 16:19:04 2018

@author: stangen
"""

from netCDF4 import Dataset, num2date
from datetime import timedelta
import pytz
import numpy as np
import EFA.duplicate_madaus.efa_functions as ef
import os
#from old_ensemble_verification import error_vs_spread
from EFA.efa_xray.state.ensemble import EnsembleState

class Load_Data():
    """
    This class exists mainly because the two functions share some variables,
    and this helps reduce the amount of redundant inputs into each function. 
    #ens_type, ob_type are str; vrbls, update_var, post_vrbls are lists of str, 
    #grid is list of ints, date is datetime object, new_format is boolean.
    date = date of the model initialization
    ens_type = eccc, ecmwf, ncep
    prior_vrbls = all the variables in the prior netCDF (T2M, ALT, IVT, etc)
    ob_type = the observation variable type we will be assimilating (only 1- T2M, ALT, IVT, etc)
    update_var = the variable type(s) we will be updating (T2M, ALT, IVT, etc)
    post_vrbls = all the variables in the posterior netCDF - only used when 
    running statistics on the posterior ensemble in mse_variance_gridded
        for use in mse_variance, ob_type is the variable type of the observation,
        without observation error variance attached with it. This is to load 
        observations, which do not have observation error variance as part of 
        their filename.
        update_var is used purely to load the netCDF we are doing stats on, 
        which may have ob error variance as part of the variable name. This can
        include the posterior, or the prior, if we are interested in the prior stats.
        ob_type and update_var should match if ob error variance was not used 
        in creating the posterior variable names, otherwise update_var contains ob err variance.
        (just 1, since it runs 1 at a time)
        
        for use in generate_gridded_obs, update_var and post_vrbls are not needed.
        
        for use in run_efa_script, post_vrbls is not needed.
        
    new stuff:
        grid: list of ints containing:
        l = left side longitude of the grid region
        r = right side longitude of the grid region
        t = top side latitude of the grid region
        b = bottom side latitude of the grid region
        s = number of degrees longitude apart obs are at equator, and number
        of degrees apart each row of obs is spaced latitudinally.
        grid is used to load filenames, and is also used in save_gridded_obs
        new_format = if True, file names will be called using the "new format"
        efh = end forecast hour (i.e. what's the last forecast hour in the forecast?)
    """
    
    def __init__(self,date,ens_type,prior_vrbls,ob_type,update_var=[],post_vrbls=[],
                 grid=[], new_format=False, efh=54):
        self.date = date  
        self.y = self.date.strftime('%Y')
        self.m = self.date.strftime('%m')
        self.d = self.date.strftime('%d')
        self.h = self.date.strftime('%H')
        self.ens_type = ens_type
        self.var_string = ef.var_string(prior_vrbls) #convert list of vrbls to string
        self.post_vrbls = post_vrbls
        #self.vrbls = vrbls
        self.ob_type = ob_type
        self.update_var = update_var
        #this allows an empty grid to be input
        try:
            self.l = grid[0]
            self.r = grid[1]
            self.t = grid[2]
            self.b = grid[3]
            self.s = grid[4]
        except:
            pass
        self.new_format = new_format
        self.efh = str(efh)+'hrs'


    
    
    
    def load_netcdfs(self,post=False,ob_cat='madis',ob_upd='ob_update_self',inf='none',lr='1000'):
        """
        Loads the ensemble netCDF and the elevation netCDF. 
        Packages and returns ensemble data into an EnsembleState (xarray)
        object for use in efa_xray code. Also returns latitudes, longitudes, 
        and elevations of the ensemble type. 
        If posterior ensemble, there are more options for exactly what EFA 
        took place for finding the right file. 
        
        post = boolean, if true, we are loading the posterior, if false, loading the prior.
        posterior options:
            ob_cat = observation category, either 'madis' or 'gridded' observations
            ob_upd = observation update, is either 'ob_update_self' or 'ob_update_all'-
            did we use the observations to update just their corresponding variables,
            or did we allow them to update all variable types in the ensemble?
            inf = inflation
            lr = localization radius we used
            
        """
        
        # directory where the ensemble of all times is
        if post==False:
            if self.new_format == False:
                infile = '/home/disk/hot/stangen/Documents/prior_ensembles/'+self.ens_type+'/'+self.y+self.m+'/'+self.y+'-'+self.m+'-'+self.d+'_'+self.h+'_'+self.ens_type+'_'+self.var_string+'.nc' 
            elif self.new_format == True:
                infile = '/home/disk/hot/stangen/Documents/prior_ensembles/'+self.ens_type+'/'+self.y+self.m+'/'+self.y+'-'+self.m+'-'+self.d+'_'+self.h+'_'+self.efh+'_'+self.var_string+'.nc'
            prior_or_post='prior'
        elif post==True:
            post_varstring = ef.var_string(self.post_vrbls)
            if self.new_format == False:
                infile = '/home/disk/hot/stangen/Documents/posterior_ensembles/'+ob_cat+'/'+ob_upd+('/inf_'+inf).replace('.','-')+'/loc_'+str(lr)+'/'+self.ens_type+'/'+self.y+self.m+'/'+self.y+'-'+self.m+'-'+self.d+'_'+self.h+'_'+self.ens_type+'_'+post_varstring+'.nc' 
            elif self.new_format == True:
                infile = '/home/disk/hot/stangen/Documents/posterior_ensembles/'+ob_cat+'/'+ob_upd+('/inf_'+inf).replace('.','-')+'/loc_'+str(lr)+'/'+self.ens_type+'/'+self.y+self.m+'/'+self.y+'-'+self.m+'-'+self.d+'_'+self.h+'_'+self.efh+'_'+str(self.l)+'_'+str(self.r)+'_'+str(self.t)+'_'+str(self.b)+'_'+str(self.s)+'_'+post_varstring+'.nc'                 
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
        elevs = orog_data.variables['orog'][0,:] # Elevation of ecmwf, eccc, ncep  
        
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
    
    def load_ens_netcdf(self, forecast_hour=0):
        """
        Loads the ensemble from n forecast hours after the forecast was 
        initialized and uses its 0-hour forecast as the "observation" grid. 
        Returns the variable of interest and the lats/lons. 
        The filepath loads non-EFA'd netCDF files (it loads the prior).
        var is  nlats x nlons x nmems. (at time t=0)
        lats, lons are 1-d float arrays (masked).
        """  
        dty, dtm, dtd, dth = ef.dt_str_timedelta(self.date,forecast_hour)
        
        print('loading analysis grid: '+self.ens_type+' '+dty+dtm+dtd+'_'+dth+'00')
        
        if self.new_format == False:
            infile = '/home/disk/hot/stangen/Documents/prior_ensembles/'+self.ens_type+'/'+dty+dtm+'/'+dty+'-'+dtm+'-'+dtd+'_'+dth+'_'+self.ens_type+'_'+self.var_string+'.nc'
        elif self.new_format == True:
            infile = '/home/disk/hot/stangen/Documents/prior_ensembles/'+self.ens_type+'/'+dty+dtm+'/'+dty+'-'+dtm+'-'+dtd+'_'+dth+'_'+self.efh+'_'+self.var_string+'.nc'
        ncdata = Dataset(infile,'r')
        #gridded ensemble 0-hour forecasts n forecast hours after the ensemble
        #we are updating was initialized
        var = ncdata.variables[self.ob_type][0,:,:,:]

        lats = ncdata.variables['lat'][:]
        lons = ncdata.variables['lon'][:]
        
        return var, lats, lons

        
    
        
    def load_obs(self, forecast_hour=6, madis=True, variance=False):
        """
        Loads the observations corresponding with n hours after ensemble was
        initialized (default 6 hours). Returns the observations from the text file.
        -Loads MADIS observations if madis=True, loads gridded observations if False.
        -If variance == True, loads gridded observations which also contain 
        ensemble variance at the ob location.
        -Gridded and MADIS observations must be generated/obtained before calling this function.
        Returns a list of observations. 
        """
        
        dty, dtm, dtd, dth = ef.dt_str_timedelta(self.date,forecast_hour)
        
        ob_str = ef.var_string([self.ob_type])        
        # directory where the observations are
        if madis==True:
            obs_file = '/home/disk/hot/stangen/Documents/surface_obs/MADIS/'+dty+dtm+'/combined_'+self.ob_type+'/'+self.ob_type+'_'+dty+dtm+dtd+'_'+dth+'00.txt'
        elif madis==False:
            if self.new_format == False:
                obs_file = '/home/disk/hot/stangen/Documents/gridded_obs/'+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/'+self.ob_type+'_'+dty+dtm+dtd+'_'+dth+'00'
            elif self.new_format == True:
                obs_file = '/home/disk/hot/stangen/Documents/gridded_obs/'+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/'+dty+dtm+dtd+'_'+dth+'00_'+str(self.l)+'_'+str(self.r)+'_'+str(self.t)+'_'+str(self.b)+'_'+str(self.s)
            if variance==True:
                obs_file += '_variance'
            obs_file += '.txt'
        print('loading '+ob_str+' obs from '+dty+dtm+dtd+'_'+dth+'00')
        f1 = open(obs_file, 'r')
        obs = f1.readlines()
        
        return obs
    
    def save_gridded_obs(self, forecast_hour=0, get_variance=False):
        """
        Loads the ensemble from n hours after the forecast was initialized
        and uses its 0-hour forecast as the "observation" grid. 
        Saves information about the generated gridded obs to a .txt file,
        with the same format as MADIS observations.
        
        If get_variance is true, will also save the variance of the ensemble
        at the observation locations at the end of each line.
        """
        basedir = '/home/disk/hot/stangen/Documents/gridded_obs/'
        
        #get the list of points
        points = ef.get_ob_points(self.l,self.r,self.t,self.b,self.s)
        #get the 0-hour ensemble forecast, n hours after the model was initialized
        dt0 = self.date
        dt = dt0.replace(minute = 0, second=0, microsecond=0)
        dt = dt + timedelta(hours=forecast_hour)
        
        #convert to epoch time for saving observation
        epoch = str(dt.replace(tzinfo=pytz.utc).timestamp())
        
        dty, dtm, dtd, dth = ef.dt_str_timedelta(self.date,forecast_hour)
        
        print('starting saving of '+self.ens_type+' gridded '+self.ob_type+ ' "obs" at: '+dty+dtm+dtd+'_'+dth+'00')
        
        var, lats, lons = self.load_ens_netcdf(forecast_hour)
            
        #obtain the mean of the ensemble (nlats x nlons)
        ens_mean = var.mean(axis=-1)
        
        #obtain variance of the ensemble (nlats x nlons)
        variance = np.var(var,axis=-1,ddof=1)
        
        #initialize the obs list to append to
        obs = []
        #loop through each point (lat/lon pair)
        for i, p in enumerate(points):
            ob_lat = p[0]
            ob_lon = p[1]
            #get the ensemble value at the lat/lon pair
            ob_value, ob_variance = ef.closest_points(ob_lat,ob_lon,lats,lons,variable=ens_mean,
                                         need_interp=True,gen_obs=True,variance=variance)
            obs.append(str(i)+','+str(ob_lat)+','+str(ob_lon)+','+str(0)+','+epoch+','+str(ob_value)+',GRIDDED,'+str(0))
            if get_variance == True:
                obs.append(','+str(ob_variance))
            obs.append('\n')
            
        #save directory for the observations
        if (os.path.isdir(basedir+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/')):
            pass
        else:
            os.makedirs(basedir+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/')
            
        #save the list of observations
        if self.new_format == False:
            f = open(basedir+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/'+self.ob_type+'_'+dty+dtm+dtd+'_'+dth+'00.txt',"w")
        elif self.new_format == True:
            savestr = basedir+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/'+dty+dtm+dtd+'_'+dth+'00_'+str(self.l)+'_'+str(self.r)+'_'+str(self.t)+'_'+str(self.b)+'_'+str(self.s)
            if get_variance == True:
                savestr += '_variance'
            savestr += '.txt'
            f = open(savestr,"w")
        for s in obs:
            f.write(s)
        f.close()
        
        
    
            

        
        
    
