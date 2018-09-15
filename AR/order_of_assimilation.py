#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:52:17 2018

@author: stangen
"""
import numpy as np
#---------------Change these---------------------------------------------------
#ens = np.array([18,19.5,20,18.5,21,19.5,21.5,19.5,20.5])
##ens = np.array([70,72,69,71,73]) #ensemble forecast at point A
#ens_2 = np.array([24,23,22,24,19,20,19.5,21.5,20])
##ens_2 = np.array([65,67,64,66,68]) #ensemble forecast at point B
#ens_3 = np.array([20,22,25,20.5,23,22,24,21.5,22.5])
ens = np.array([65,68,70,66,67,64,67,70,68])
ens_2 = np.array([61,63,62,64,61,61,61,63,64])
ens_3 = np.array([56,54,57,59,56,56,56,58,59])


ens_mean = np.mean(ens)
ens_2_mean = np.mean(ens_2)
ens_3_mean = np.mean(ens_3)
#obs = np.array([20, 20, 20,20,20,20,10])
obs = np.array([64])
obs2 = np.array([62.5])
#obs = np.array([73]) #observation(s)
R = .10#observation error variance
loc_scalar = 1 #localization factor between 2 points
#------------------------------------------------------------------------------

pert = ens - ens_mean
pert2 = ens_2-ens_2_mean
pert3 = ens_3-ens_3_mean

cov = np.dot(pert,pert2) /(len(ens)-1)
cov3 = np.dot(pert,pert3) /(len(ens)-1)
cov2 = np.dot(pert2,pert3) /(len(ens)-1)

print('Point A ensemble mean: ',ens_mean)
print('Point B ensemble mean: ',ens_2_mean)
print('Point C ensemble mean: ',ens_3_mean)
print('Point A ensemble variance :',np.var(pert,ddof=1))
print('Point B ensemble variance :',np.var(pert2,ddof=1))
print('Point C ensemble variance :',np.var(pert3,ddof=1))
print('Covariance between Points A and B: ',cov)
print('Correlation between points A and B: ',cov/(np.std(pert,ddof=1)*np.std(pert2,ddof=1)))
print('Correlation between points A and C: ',cov3/(np.std(pert,ddof=1)*np.std(pert3,ddof=1)))
print('Correlation between points B and C: ',cov2/(np.std(pert2,ddof=1)*np.std(pert3,ddof=1)))


for ob in obs:
    print('Point A perturbations before assimilation: ',pert)
    print('Point B perturbations before assimilation: ',pert2)
    print('Point C perturbations before assimilation: ',pert3)
    print('Observation assimilating into point A: ',ob)
    ye_var = np.var(pert,ddof=1)
    
    cov = np.dot(pert,pert2) /(len(ens)-1)
    cov3 = np.dot(pert,pert3) /(len(ens)-1)

    K = ye_var/(ye_var+R)
    K2 = loc_scalar*cov/(ye_var+R)
    K3 = loc_scalar*cov3/(ye_var+R)
    
    print('K: ',K)
    print('K2: ',K2)
    print('K3: ',K3)
    
    innov = ob-ens_mean
    print('innov: ',innov)
    ens_mean = ens_mean+K*(innov)
    ens_2_mean = ens_2_mean+K2*innov
    ens_3_mean = ens_3_mean+K3*innov
       
    B = 1/(1+np.sqrt(R/(ye_var+R)))
    print('beta :',B)
    
    pert2 = pert2 - B*K2*pert
    pert3 = pert3 - B*K3*pert
    pert = pert - B*K*pert
    
    print('Point A ensemble mean after assimilation: ',ens_mean)
    print('Point A perturbations after assimilation :',pert)
    print('Point A ensemble variance after assimilation: ',np.var(pert,ddof=1))
    print('Point B ensemble mean after assimilation: ',ens_2_mean)
    print('Point B perturbations after assimilation: ',pert2)
    print('Point B ensemble variance after assimilation: ',np.var(pert2,ddof=1))
    print('Point C ensemble mean after assimilation: ',ens_3_mean)
    print('Point C perturbations after assimilation: ',pert3)
    print('Point C ensemble variance after assimilation: ',np.var(pert3,ddof=1))
    
    print('Covariance between points after assimilation: ',np.dot(pert,pert2) /(len(ens)-1))
    print('Correlation between points A and B after assimilation: ',(np.dot(pert,pert2) /(len(ens)-1))/(np.std(pert,ddof=1)*np.std(pert2,ddof=1)))

    print('Correlation between points B and C after assimilation: ',(np.dot(pert2,pert3) /(len(ens)-1))/(np.std(pert2,ddof=1)*np.std(pert3,ddof=1)))



for ob in obs2:
    print('Point A perturbations before assimilation: ',pert)
    print('Point B perturbations before assimilation: ',pert2)
    print('Point C perturbations before assimilation: ',pert3)
    print('Observation assimilating into point B: ',ob)
    ye_var = np.var(pert2,ddof=1)
    
    cov = np.dot(pert2,pert) /(len(ens)-1)
    cov3 = np.dot(pert2,pert3) /(len(ens)-1)

    K = loc_scalar*cov/(ye_var+R)
    K2 = ye_var/(ye_var+R)   
    K3 = loc_scalar*cov3/(ye_var+R)
    
    print('K: ',K)
    print('K2: ',K2)
    print('K3: ',K3)
    
    innov = ob-ens_2_mean
    print('innov: ',innov)
    ens_mean = ens_mean+K*(innov)
    ens_2_mean = ens_2_mean+K2*innov
    ens_3_mean = ens_3_mean+K3*innov
       
    B = 1/(1+np.sqrt(R/(ye_var+R)))
    print('beta :',B)
    
    pert3 = pert3 - B*K3*pert2
    pert = pert - B*K*pert2
    pert2 = pert2 - B*K2*pert2
    
    print('Point A ensemble mean after assimilation: ',ens_mean)
    print('Point A perturbations after assimilation :',pert)
    print('Point A ensemble variance after assimilation: ',np.var(pert,ddof=1))
    print('Point B ensemble mean after assimilation: ',ens_2_mean)
    print('Point B perturbations after assimilation: ',pert2)
    print('Point B ensemble variance after assimilation: ',np.var(pert2,ddof=1))
    print('Point C ensemble mean after assimilation: ',ens_3_mean)
    print('Point C perturbations after assimilation: ',pert3)
    print('Point C ensemble variance after assimilation: ',np.var(pert3,ddof=1))
    
    print('Covariance between points after assimilation: ',np.dot(pert,pert2) /(len(ens)-1))
    print('Correlation between points A and B after assimilation: ',(np.dot(pert,pert2) /(len(ens)-1))/(np.std(pert,ddof=1)*np.std(pert2,ddof=1)))

    print('Correlation between points B and C after assimilation: ',(np.dot(pert2,pert3) /(len(ens)-1))/(np.std(pert2,ddof=1)*np.std(pert3,ddof=1)))


