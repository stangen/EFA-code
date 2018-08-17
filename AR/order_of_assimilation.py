#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:52:17 2018

@author: stangen
"""
import numpy as np
#---------------Change these---------------------------------------------------
#ens = np.array([18,19.5,20,18.5,21,19.5,21.5,19.5,20.5])
ens = np.array([70,72,69,71,73]) #ensemble forecast at point A
#ens_2 = np.array([24,23,22,24,19,20,19.5,21.5,20])
ens_2 = np.array([65,67,64,66,68]) #ensemble forecast at point B
ens_mean = np.mean(ens)
ens_2_mean = np.mean(ens_2)
#obs = np.array([10, 20, 20,20,20,20,20])
obs = np.array([73]) #observation(s)
R = 0 #observation error variance
loc_scalar = 1 #localization factor between 2 points
#------------------------------------------------------------------------------

pert = ens - ens_mean
pert2 = ens_2-ens_2_mean

cov = np.dot(pert,pert2) /(len(ens)-1)

print('Point A ensemble mean: ',ens_mean)
print('Point B ensemble mean: ',ens_2_mean)
print('Point A ensemble variance :',np.var(pert,ddof=1))
print('Point B ensemble variance :',np.var(pert2,ddof=1))
print('Covariance between Points A and B: ',cov)

for ob in obs:
    print('Point A perturbations before assimilation: ',pert)
    print('Point B perturbations before assimilation: ',pert2)
    print('Observation assimilating into point A: ',ob)
    ye_var = np.var(pert,ddof=1)
    
    cov = np.dot(pert,pert2) /(len(ens)-1)

    K = ye_var/(ye_var+R)
    K2 = loc_scalar*cov/(ye_var+R)
    
    innov = ob-ens_mean
    ens_mean = ens_mean+K*(innov)
    ens_2_mean = ens_2_mean+K2*innov
       
    B = 1/(1+np.sqrt(R/(ye_var+R)))
    
    pert2 = pert2 - B*K2*pert
    pert = pert - B*K*pert
    
    print('Point A ensemble mean after assimilation: ',ens_mean)
    print('Point A perturbations after assimilation :',pert)
    print('Point A ensemble variance after assimilation: ',np.var(pert,ddof=1))
    print('Point B ensemble mean after assimilation: ',ens_2_mean)
    print('Point B perturbations after assimilation: ',pert2)
    print('Point B ensemble variance after assimilation: ',np.var(pert2,ddof=1))

