#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:52:17 2018

@author: stangen
"""
import numpy as np
ens = np.array([18,19.5,20,18.5,21,19.5,21.5,19.5,20.5])
ens_2 = np.array([22,23,23.5,22,24,22.5,25,22,23])
ens_mean = np.mean(ens)
ens_2_mean = np.mean(ens_2)
obs = np.array([10, 20, 20,20,20,20,20])
R = 1
loc_scalar = 1

pert = ens - ens_mean
pert2 = ens_2-ens_2_mean

cov = np.dot(pert,pert2) /(len(ens)-1)

print(ens_mean)
print(ens_2_mean)
print(np.var(pert,ddof=1))
print(np.var(pert2,ddof=1))
print(cov)


for ob in obs:
    ye_var = np.var(pert,ddof=1)
    
    cov = np.dot(pert,pert2) /(len(obs)-1)

    K = ye_var/(ye_var+R)
    K2 = loc_scalar*cov/(ye_var+R)
    
    innov = ob-ens_mean
    ens_mean = ens_mean+K*(innov)
    ens_2_mean = ens_2_mean+K2*innov
    print('ens1 mean: ',ens_mean)
    print('ens2 mean: ',ens_2_mean)
    
    B = 1/(1+np.sqrt(R/(ye_var+R)))
    print(K2)
    #print(B*K2*pert)
    
    pert2 = pert2 - B*K2*pert
    pert = pert - B*K*pert
    #print(pert2)

    #print(pert2)
    
    print('ens1 variance: ',np.var(pert,ddof=1))
    print('ens2 variance: ',np.var(pert2,ddof=1))


