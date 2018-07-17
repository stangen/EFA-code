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
    #ens_type, ob_type are str; vrbls, update_var are lists of str, date is datetime object.
    date = date of the model initialization
    ens_type = eccc, ecmwf
    prior_vrbls = all the variables in the prior netCDF
    ob_type = the observation variable type we will be assimilating (only 1)
    update_var = the variable type(s) we will be updating
    post_vrbls = all the variables in the posterior netCDF
        for use in mse_variance, ob_type is the variable type of the observation,
        without observation error variance attached with it. This is to load 
        observations, which do not have ob err var as part of the name.
        update_var is used purely to load the netCDF we are doing stats on, 
        which may have ob error variance as part of the variable name. This can
        include the posterior, or the prior, if we are interested in the prior stats.
        These should match if ob error variance was not used in creating the posterior
        variable names, otherwise update_var contains ob err variance.
        (just 1, since it runs 1 at a time)
    """
    
    def __init__(self,date,ens_type,prior_vrbls,ob_type,update_var,post_vrbls=[],
                 l=-180, r=180, t=90, b=0, s=3):
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
        self.l = l
        self.r = r
        self.t = t
        self.b = b
        self.s = s


    
    
    
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
            post_varstring = ef.var_string(self.post_vrbls)
            infile = '/home/disk/hot/stangen/Documents/posterior_ensembles/'+ob_cat+'/'+ob_upd+('/inf_'+inf).replace('.','-')+'/loc_'+str(lr)+'/'+self.ens_type+'/'+self.y+self.m+'/'+self.y+'-'+self.m+'-'+self.d+'_'+self.h+'_'+post_varstring+'.nc' 
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
    
    def load_ens_netcdf(self, forecast_hour=0):
        """
        Loads the ensemble from n forecast hours after the forecast was 
        initialized and uses its 0-hour forecast as the "observation" grid. 
        Returns the variable of interest and the lats/lons. 
        var is  nlats x nlons. (at time t=0, ens mean)
        lats, lons are 1-d float arrays (masked).
        """  
        dty, dtm, dtd, dth = ef.dt_str_timedelta(self.date,forecast_hour)
        
        infile = '/home/disk/hot/stangen/Documents/prior_ensembles/'+self.ens_type+'/'+dty+dtm+'/'+dty+'-'+dtm+'-'+dtd+'_'+dth+'_'+self.ens_type+'_'+self.var_string+'.nc'
        
        ncdata = Dataset(infile,'r')
        #gridded ensemble 0-hour forecasts n forecast hours after the ensemble
        #we are updating was initialized
        var = ncdata.variables[self.ob_type][0,:,:,:]

        lats = ncdata.variables['lat'][:]
        lons = ncdata.variables['lon'][:]
        
        #obtain the mean of the ensemble
        var = var.mean(axis=-1)
        
        return var, lats, lons
    
        
    def load_obs(self, forecast_hour=6, madis=True):
        """
        Loads the observations corresponding with n hours after ensemble was
        initialized (default 6 hours). Returns the observations from the text file.
        Loads MADIS observations if True, loads gridded observations if False.
        """
        
        dty, dtm, dtd, dth = ef.dt_str_timedelta(self.date,forecast_hour)
        
        ob_str = ef.var_string([self.ob_type])        
        # directory where the observations are
        if madis==True:
            obs_file = '/home/disk/hot/stangen/Documents/surface_obs/MADIS/'+dty+dtm+'/combined_'+self.ob_type+'/'+self.ob_type+'_'+dty+dtm+dtd+'_'+dth+'00.txt'
        elif madis==False:
            obs_file = '/home/disk/hot/stangen/Documents/gridded_obs/'+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/'+self.ob_type+'_'+dty+dtm+dtd+'_'+dth+'00.txt'
        print('loading '+ob_str+' obs from '+dty+dtm+dtd+'_'+dth+'00')
        #print(obs_file)
        f1 = open(obs_file, 'r')
        obs = f1.readlines()
        
        return obs
    
    def save_gridded_obs(self, forecast_hour=0):
        """
        Loads the ensemble from n hours after the forecast was initialized
        and uses its 0-hour forecast as the "observation" grid. 
        Saves information about the generated gridded obs to a .txt file,
        with the same format as MADIS observations.
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
        
        #initialize the obs list to append to
        obs = []
        #loop through each point (lat/lon pair)
        for i, p in enumerate(points):
            ob_lat = p[0]
            ob_lon = p[1]
            #get the ensemble value at the lat/lon pair
            ob_value = ef.closest_points(ob_lat,ob_lon,lats,lons,variable=var,
                                         need_interp=True,gen_obs=True)
            obs.append(str(i)+','+str(ob_lat)+','+str(ob_lon)+','+str(0)+','+epoch+','+str(ob_value)+',GRIDDED,'+str(0)+'\n')
            
        #save directory for the observations
        if (os.path.isdir(basedir+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/')):
            pass
        else:
            os.makedirs(basedir+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/')
            
        #save the list of observations
        #f = open(basedir+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/'+self.ob_type+'_'+dty+dtm+dtd+'_'+dth+'00.txt',"w")
        f = open(basedir+self.ens_type+'/'+dty+dtm+'/'+self.ob_type+'/'+dty+dtm+dtd+'_'+dth+'00_'+str(self.l)+'_'+str(self.r)+'_'+str(self.t)+'_'+str(self.b)+'_'+str(self.s)+'.txt',"w")
        for s in obs:
            f.write(s)
        f.close()
        
        
    
            

        
        
    
