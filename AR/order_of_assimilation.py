#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:52:17 2018

@author: stangen
"""
import numpy as np
ens = np.array([18,19.5,20,18.5,21,19.5,21.5,19.5,20.5])
ens_mean = np.mean(ens)
obs = np.array([10, 20, 20,20,20,20,20])
R = .01

pert = ens - ens_mean
print(ens_mean)
print(np.var(pert,ddof=1))


for ob in obs:
    ye_var = np.var(pert,ddof=1)

    K = ye_var/(ye_var+R)
    
    innov = ob-ens_mean
    ens_mean = ens_mean+K*(innov)
    print(ens_mean)
    
    B = 1/(1+np.sqrt(R/(ye_var+R)))
    
    pert = pert - B*K*pert
    
    print(np.var(pert,ddof=1))


