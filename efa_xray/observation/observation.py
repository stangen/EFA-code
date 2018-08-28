#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import scipy.stats.mstats as ss
import pandas as pd
import xarray
import _pickle as cPickle
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, time
import pytz
from EFA.efa_xray.state.ensemble import EnsembleState
import sys, os
import EFA.duplicate_madaus.efa_functions as ef
sys.path.append('../data_parsers')


class Observation:
    def __init__(self,value=None,obtype=None,time=None,error=None,lat=None,
                 lon=None, vert=None,
                 prior_mean=None, post_mean=None, prior_var=None, post_var=None,
                 assimilate_this=False,description=None,localize_radius=None):
        self.value = value
        self.obtype = obtype
        self.time = time
        self.error = error
        self.lat = lat
        self.lon = lon
        self.vert = vert
        self.prior_mean = prior_mean
        self.post_mean = post_mean
        self.prior_var = prior_var
        self.post_var = post_var
        self.assimilate_this = assimilate_this
        self.assimilated = False
        self.description = description
        self.localize_radius = localize_radius
        


    def estimate(self,state):
        """ Given an EnsembleState, compute the ensemble
        estimate of this observation
        
        In the future, this could have the option to call some
        external function for more complex forward operators.
        For now, we assume that the state has a field that
        identically matches this observation type and simply
        interpolate to that point.
        """
        #ST I think this actually interpolates a value to the exact location
        #of the observation lat/lon, maybe from other nearby points?
        return state.interpolate(self.obtype, self.time, self.lat, self.lon)
        

    def distance_to_state(self,state):
        """ Return the distance from this ob to all locations in the state
        vector (in km) """
        return state.distance_to_point(self.lat, self.lon) 


    def localize(self,state,type='GC',full_state=False):
        """ Given a state vector object, assume location is in lat/lon and compute a
        Gaspari-Cohn weighting function with the specified halfwidth (in km) """
        # Get the localization halfwidth from 
        halfwidth = self.localize_radius

        if isinstance(state, EnsembleState):
            # Get distance to all points in the state
            distances = state.distance_to_point(self.lat, self.lon)
        else:
            # Must be a list of observations
            ourloc = (self.lat, self.lon)
            other_lats = [ob.lat for ob in state]
            other_lons = [ob.lon for ob in state]
            distances = np.array([haversine(ourloc,s) for s in zip(other_lats,
                                                                   other_lons)])
        
        # If halfwidth is None, return an array of ones
        if halfwidth is None:
            localization = np.ones(distances.shape)
       
        # Can have options here for other localization types
        # Gaspari-Cohn is most common
        if type == 'GC':
            localization = gaspari_cohn(distances, halfwidth)

            
        elif type.startswith('statsig'):
            #if not within 4000 km of the ob, don't let ob update the point.
            localization = np.ones(distances.shape)
            boolarray = distances > 4000
            localization[boolarray] = 0

            

        # Reshaping will be handled during the assimilation
        # For now, just return the array
        return localization



    def stat_sig(self,kcov,state,ye):
        #want to check if correlation is statistically significant- 
        #don't get the first value (this is a correlation of ye with itself,
        #we just want the correlation of ye with the state- must be same length as state)
        r_crit = ef.get_r_crit()
        #get ranks/rank perts of ob estimate and state
        ye_rank = np.argsort(ye)
        ye_rank = np.argsort(ye_rank)
        #sort once to obtain the indices of a sorted array
        state_1sort = np.argsort(state,axis=1)
        #this block obtains rank from sorted array
        n = np.arange(state_1sort.shape[0])[:, np.newaxis]
        state_rank = np.zeros(state_1sort.shape)
        state_rank[n,state_1sort] = np.arange(len(state_1sort[0]))
        
#        print(state[0])
#        print(state_1sort[0])
#        print(state_rank[0])
        
        ye_rank_pert = ye_rank - np.mean(ye_rank)
        state_rank_pert = state_rank - np.mean(state_rank,axis=1)[:,None]
        #calculate covariance of ranks, then spearman correlation
        cov_rank = np.dot(state_rank_pert,ye_rank_pert)/len(ye)
        spearmanr = cov_rank/(np.std(ye_rank_pert)*np.std(state_rank_pert,axis=1))
        #find which values are below the critical r value (and therefore fail)
        r_fail = abs(spearmanr) < r_crit[len(ye)][self.localize_radius]
#        print(r_crit[len(ye)][self.localize_radius])
        
        kcov[r_fail] = 0
        
        return kcov



    def map_localization(self, state, m, type='GC'):
        """ Function to map localization radius 
        Requires:
            state --> The state vector we are using
            m     --> A basemap instance for projecting the map

            """
        # Get the localization weights and reshape to ny x nx
        localization = self.localize(state, type=type)
        # Get map projected coordinates from the state
        gx, gy = state.project_coordinates(m)
        # Make the plot
        F = plt.figure(figsize=(10,8))
        plt.pcolormesh(gx,gy,localization,vmin=0.0,vmax=1.0)
        plt.colorbar()
        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()
        plt.title('Localization Weights for {:s} ({:5.3f},{:5.3f})'.format(self.description,
                                                                           self.lat,
                                                                           self.lon))
        plt.show()

def gaspari_cohn(distances, halfwidth):
    """ Compute Gaspari-Cohn weights from a distance array and
    a given halfwidth """
    r = np.divide(distances, abs(halfwidth))
    

    # Do the Gaspari Cohn weighting
    weights = np.zeros(r.shape)
    # Less than halfwidth
    weights[r <= 1.0] = ((((-0.25*r+0.5)*r+0.625)*r-5.0/3.0) * r**2 + 1.0)[r <= 1.0]
    # Between halfwidth and fullwidth
    weights[(r > 1.0) & (r < 2.0)] = (((((r/12.0 - 0.5)*r + 0.625) *r+\
                                    5.0/3.0)*r-5.0)*r + 4.0 -\
                                    2.0 / (3.0 * r))[(r > 1.0) & (r < 2.0)]
    return weights

    


def haversine(loc1,loc2):
    """ Use Haversine formula to compute the distance between two lat-lon
    coordinate pairs """
    R = 6371. # Radius of earth in kilometers
    lat1 = np.radians(loc1[0])
    lat2 = np.radians(loc2[0])
    dlat = lat2 - lat1
    dlon = np.radians(loc2[1] - loc1[1])

    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    return R * c

