#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from copy import deepcopy
from EFA.efa_xray.assimilation.assimilation import Assimilation


class EnSRF(Assimilation):
    """
    Class to do the EnSRF update
    Originator: G. J. Hakim

     Input variables are
     Xbp => Array of ensemble estimates of state perturbations from mean (num_state_vars x num_ens_mems)
     xbm => ensemble mean (vector of length num_state_vars)
     Y => Observations (num_observations)
     H => Translation matrix from state vector to obs (num_observations x num_state_vars)
     R => Observation error covariance matrix (num_observations x num_observations)
     loc => Localization
    #print "Y shape:", np.shape(Y)
    #print "H shape:", np.shape(H)
    #print "R shape:", np.shape(R)
    #print "loc shape:", np.shape(loc)

    # Modified to work with XRAY enkf --> L. Madaus 2/20/2015
    """
    
    def __init__(self, state, obs, nproc=1, inflation=None, verbose=True, loc=False):
        # Initialize the Assimilation inheritance
        Assimilation.__init__(self, state, obs, nproc, inflation, verbose)
        self.loc = loc

    def update(self):
        if self.verbose: print("Beginning update sequence")
        # Make a dummy localization here to allocate thisobs_err
        # array only once
        state_shape = self.prior.shape()[:-1] # Don't include mems
        #shape of state_shape is nvars x ntimes x nlats x nlons
        dum_localize = np.ones(state_shape)
        Nstate = self.prior.nstate()
        Nens = self.prior.nmems()

        # Do pre-processing to estimate obs and format
        # as state vector
        xam, Xap = self.format_prior_state()
        
        if self.loc == 'statsig2':
            Xap_start = deepcopy(Xap)
        
        #make a deep copy of original if using initial correlation between ob and state


        numobs_assim = 0
        # Now loop over all observations
        if self.verbose: print("Beginning observation loop")
        for obnum,ob in enumerate(self.obs):
#            if (obnum % 100==0) and self.verbose: print("    On ob:", obnum)
            print('on ob '+str(obnum))
            # Reset the mean and perturbations
            
            #ST
            #xbm a row vector, xbp column-state matrix with columns of ensemble members?
            xbm = xam
            Xbp = Xap
            # Make a vector of all of the ensemble members
            # estimates of state translated into observation
            # space (by H)
            #print "Mean", np.tile(xbm,(Nens,1))
            #print np.transpose(np.tile(xbm,(Nens,1))) + Xbp
            #print H
            H = np.zeros(xam.shape)
            H[Nstate+obnum] = 1.0
            #ST This gets the mean/pertubation of just the interpolated ob value.
            mye = np.dot(H, xbm)
            ye = np.dot(H, Xbp)

            ob.prior_mean = mye
            #print "ye", ye
            # Find the variance among the ensemble members
            varye = np.var(ye, ddof=1)
            ob.prior_var = varye

            # IMPORTANT --- here we check to see if we should actually
            # assimilate this ob
            if not ob.assimilate_this:
                ob.assimilated = False
                continue


            # And find the observation error variance from the R matrix
            # (Assumes ob errors are uncorrelated)
            obs_err = ob.error
            #print(obs_err)
            # Find the innovation --the difference between the ob value
            # and the ensemble mean estimate of the ob
            # This is y-HXb
            innov = ob.value - mye

            # Now find the innovation variance -- the sum of the variance of the ob
            # and the varaiance of the ensemble estimate
            # This goes into the denominator of the Kalman gain
            kdenom = (varye + obs_err)

            # The numerator of the Kalman gain is the covariance between
            # the ensemble members and the obs-transformed ensemble members
            #ST I think what happens here is that the ye is only a 1xnmems vector,
            #where the ob is interpolated to be, and then the covariance of the
            #entire grid of all the values of a certain variable with that point.
            kcov = np.dot(Xbp,np.transpose(ye)) / (Nens-1)


            # Option to localize the gain
            if self.loc not in [None, False]:
                # Project the localization
#                if self.loc == 'statsig':
#                    state_localize = ob.localize(self.prior, type=self.loc,Xbp,ye)
#                elif self.loc == 'statsig2':
#                    state_localize = ob.localize(self.prior, type=self.loc,Xap_start,ye)
#                
#                else:
                state_localize = ob.localize(self.prior, type=self.loc)
                #print(state_localize.shape)
                #print(dum_localize.shape)
                # LEM---CHECK TO BE SURE THIS LOGIC WORKS FOR 1-D LATLON!!!
                # This needs to be projected into the full state vector
                # Check for (ny,nx) or just 1-d
                # Project this across all vars, all times
                if len(state_localize.shape) == 2:
                    state_localize = (state_localize[None,None,:,:] * dum_localize).flatten()
                else:
                    state_localize = (state_localize[None,None,None,:] * dum_localize).flatten()
                # Now need to localize for obs
                obs_localize = ob.localize(self.obs, type=self.loc)
                state_localize = np.hstack((state_localize, obs_localize))
                kcov = np.multiply(state_localize,kcov)
                #kcov = np.dot(kcov,np.transpose(loc[ob,:]))
            
            
            
            if self.loc.startswith('statsig'):
                nvars = self.prior.nvars()
                ntimes = self.prior.ntimes()
                nlats = self.prior.ny()
                nlons = self.prior.nx()
                #count the number of gridpoints within 4000 km for each forecast hour
                npoints_less4000 = np.count_nonzero(kcov)/(nvars*ntimes)
#                print(npoints_less4000/(nlats*nlons))
                if self.loc == 'statsig':
                    kcov = ob.stat_sig(kcov,Xbp,ye)
                elif self.loc == 'statsig2':
                    kcov = ob.stat_sig(kcov,Xap_start,ye)
                #find percentage of points that were updated by ob for each forecast hour
                npoints_12hr = 0
                npoints_24hr = 0
                npoints_36hr = 0
                npoints_48hr = 0
                percent = {}
                percent['12'] = percent.get('12',[])
                percent['24'] = percent.get('24',[])
                percent['36'] = percent.get('36',[])
                percent['48'] = percent.get('48',[])
                for i in range(0,self.prior.nvars()):
                    npoints_12hr = npoints_12hr + np.count_nonzero(kcov[i*nvars*ntimes*nlats*nlons+2*nlats*nlons:i*nvars*ntimes*nlats*nlons+3*nlats*nlons])
#                    print(i*nvars*ntimes*nlats*nlons+2*nlats*nlons)
#                    print(i*nvars*ntimes*nlats*nlons+3*nlats*nlons)
                    npoints_24hr = npoints_24hr + np.count_nonzero(kcov[i*nvars*ntimes*nlats*nlons+4*nlats*nlons:i*nvars*ntimes*nlats*nlons+5*nlats*nlons])
                    npoints_36hr = npoints_36hr + np.count_nonzero(kcov[i*nvars*ntimes*nlats*nlons+6*nlats*nlons:i*nvars*ntimes*nlats*nlons+7*nlats*nlons])
                    npoints_48hr = npoints_48hr + np.count_nonzero(kcov[i*nvars*ntimes*nlats*nlons+8*nlats*nlons:i*nvars*ntimes*nlats*nlons+9*nlats*nlons])
                percent['12'].append(npoints_12hr/npoints_less4000)
                percent['24'].append(npoints_24hr/npoints_less4000)
                percent['36'].append(npoints_36hr/npoints_less4000)
                percent['48'].append(npoints_48hr/npoints_less4000)
#                print(npoints_less4000)
#                print(npoints_12hr)
#                print(npoints_24hr)
#                print(npoints_36hr)
#                print(npoints_48hr)
#                print(npoints_12hr/npoints_less4000)
#                print(npoints_24hr/npoints_less4000)
#                print(npoints_36hr/npoints_less4000)
#                print(npoints_48hr/npoints_less4000)


            
            
            # Compute the Kalman gain
            kmat = np.divide(kcov, kdenom)
            #kmat = np.divide(kcov,kdenom)
            #print "kmat", kmat.shape
            #print "innov", innov.shape

            # Now do the updates
            # First update the mean
            #xam = xbm + np.dot(np.dot(H,kmat),innov)
            #print "kmat", np.shape(kmat)
            #print "innov", np.shape(innov)
            #xam = xbm + np.dot(kmat,innov)
            xam = xbm + np.multiply(kmat,innov)
            

            # And each ensemble member perturbation
            # This is the "Square Root" 
            # step in the Kalman filter equations
            beta = 1./(1. + np.sqrt(obs_err/(varye+obs_err)))
            
            #ST artificially inflating the posterior pertubations by reducing beta.
            #this will reduce the subtraction from the prior pertubations
#            beta_const = 1.1
#            beta = beta/(beta_const)
            
            kmat = np.multiply(beta,kmat)

            ye = np.array(ye)[np.newaxis]
            kmat = np.array(kmat)[np.newaxis]

            Xap = Xbp - np.dot(kmat.T, ye)
            
            check_ob_estimate = False
            if check_ob_estimate == True:
                print('posterior ob estimate mean: ',np.dot(H,xam))
                print('posterior ob estimate perts:\n',np.dot(H,Xap))

            # For reference, grab the post mean and variance
            post_ye = np.dot(H,xam)
            post_var = np.var(np.dot(H,Xap), ddof=1)
            ob.post_mean = post_ye
            ob.post_var = post_var
            # Record that this ob was assimilated
            ob.assimilated = True
            numobs_assim = numobs_assim+1
#            print('total number of obs assimilated so far: ',numobs_assim)
        print('total number of obs assimilated: ',numobs_assim)
        
        
        if self.loc.startswith('statsig'):
            for i in percent:
                avg = np.mean(np.array(percent[i]))
                print('average number of points updated for forecast hour '+i+': '+str(avg))
            
        #print('beta constant: ',beta_const)
        # After having assimilated everything, rebuild the state
        return self.format_posterior_state(xam, Xap)
        
