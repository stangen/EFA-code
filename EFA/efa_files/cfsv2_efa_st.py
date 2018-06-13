# -*- coding: utf-8 -*-
"""
Created on Thu May 19 11:50:35 2016

@author: njweber2
"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
from nicks_files.operational_cfsv2 import get_cfsv2_ensemble
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import nicks_files.cfs_utilities as ut
import time
import os
#from old_ensemble_verification import error_vs_spread
# Luke's (super useful) assimilation tools:
from efa_xray.state.ensemble import EnsembleState
from efa_xray.observation.observation import Observation
from efa_xray.assimilation.ensrf import EnSRF

def checkdir(directory):
    """
    Checks to see if the directory exists. If not, the directory is created.
    """
    if not os.path.isdir(directory):
        os.system('mkdir {}'.format(directory))
        
def build_state(nlag, start, end, vrbls=['Z500'], fromfile=True, prior=True,
                loc='NO', coarse_obs=False, recent_mean=False, obtimes=1,
                inflated=False, historic_cov=False, postfile=None, oberr=None):
    """
    Acquires and returns a CFSv2 operational ensemble state.
    
    Requires:
    nlag -----> The number of lag days in the ensemble construction.
    start ----> A datetime object of the ensemble forecast initialization time.
    end ------> A datetime object of the ensemble forecast end time.
    vrbls ----> A list of variables we want to retrieve.
                Default: 500 mb heights
    fromfile -> A boolean object. If True, the ensemble state is loaded from
                a preprocessed ensemble netcdf file. If False, the ensemble is
                created with the get_cfsv2_ensemble function.
    prior ----> A boolean object. If True, the ensemble state of the prior is 
                loaded (i.e., the un-adjusted forecast). If false, the 
                posterior (adjust forecast) ensemble forecast is loaded.
    inflated -> Are we loading a posterior from an inflated prior?
                
    Returns:
    statecls -> An ensemble state object (see EnsembleState class)
    """
    ndays = (end-start).days
    if fromfile:
        # Load the ensemble state data from an ensemble netcdf file
        if prior:
            indir = '/home/disk/hot/stangen/Documents/GEFS/ensembles' + \
            '/2017081400_21mem_1days.nc'
#            if inflated:
#                indir += '/inflated_prior'
#            infile = '{}/{:%Y%m%d%H}_{:02d}mem_{}days.nc'.format(indir, start,
#                                                                 nlag*16, ndays)
        else:
            if postfile is None:
                indir = '/home/disk/vader2/njweber2/research/subseasonal/efa/' + \
                        'adjusted_ensembles'
                if historic_cov:
                    indir += '/historic_cov'
                if inflated:
                    indir += '/inflated_prior'
                if coarse_obs: 
                    indir += '/coarse_obs'
                if recent_mean: 
                    indir += '/most_recent_mean'
                if obtimes > 1: 
                    indir += '/{}obtimes'.format(obtimes)
                if len(vrbls)==1:
                    infile = '{}/{}_{:%Y%m%d%H}_{:02d}mem_{}days_{}loc'.format(indir,
                                                  vrbls[0], start, nlag*16, ndays, loc)
                else:
                    infile = '{}/{:%Y%m%d%H}_{:02d}mem_{}days_{}loc'.format(indir,
                                                            start, nlag*16, ndays, loc)
                if oberr is None:
                    infile += '.nc'
                else:
                    infile += '_{}oberror.nc'.format(int(oberr))
            else:
                infile = postfile
        with Dataset(infile,'r') as ncdata:
            times = ncdata.variables['time']
            ftimes = num2date(times[:],
                              times.units)
            lats = ncdata.variables['lat'][:]
            lons = ncdata.variables['lon'][:]
            mems = ncdata.variables['ens'][:]
            allvars = {}
            for var in vrbls:
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
    
    else:
        # Create the ensemble state with get_cfsv2_ensemble
        statecls =  get_cfsv2_ensemble(nlag, start, end, vrbls=vrbls,
                                       writenc=False)
    return statecls
    
def get_gdas_obs(atime, vrbls=['Z500'], sampling_interval=5, cutofflat=None,
                 loc='GC2000', err=100., errvar_file=None):
    """
    Gets a list of 'observations' from the CFSv2 f00 forecast (i.e., GDAS
    analysis) grid.
    
    Requires:
    atime ------------> A datetime object indicating the time of the desired
                        observations.
    vrbls ------------> A list of variables to get observations of.
                        Default: 500 mb heights.
    sampling_interval-> An integer specifying the interval at which to sample
                        the gridded analysis for observations.
                        Default: Sample every 5 degrees lat/lon. (1 deg grid)
    cutofflat --------> Grid points poleward of this latitude are ignored.
                        
    Returns:
    observations -----> A list of Observation objects.
    """
    # Set the localization radius
    if loc[:2]=='GC':
        locrad = int(loc[2:])
    else:
        locrad = 2000.0 # kilometers
    
    # Corresponding parameters in the analysis file
    apars = {'Z500'   : 'HGT_500mb',
             'CHI200' : 'VPOT_200mb'}
    # Set I/O variable(s)
    anl_file = '/home/disk/vader2/njweber2/research/subseasonal/efa/' + \
               'verification/analyses_01Dec2015-01Mar2016.nc'
               
    # Load the data
    vardata = {}
    with Dataset(anl_file,'r') as ncdata:
        # Load dimensions
        timenums = ncdata.variables['time']
        times = num2date(timenums[:],timenums.units)
        lats  = ncdata.variables['latitude'][:]
        lons  = ncdata.variables['longitude'][:]
        # Cut off polar grid points:
        if cutofflat is not None:
            yi = ut.nearest_ind(lats,-cutofflat)
            yf = ut.nearest_ind(lats,cutofflat) + 1
        else:
            yi = 0
            yf = len(lats)
        lats = lats[yi:yf]
        # Get a lat/lon value for every location
        lonarr, latarr = np.meshgrid(lons, lats)
        # Sample only a fraction of the points in the analysis
        lonarr = lonarr[::sampling_interval,::sampling_interval].flatten()
        latarr = latarr[::sampling_interval,::sampling_interval].flatten()
        # Find the desired time
        t_ind = np.where(times==atime)[0][0]
        # Load the variable
        for var in vrbls:
            vardata[var] = ncdata.variables[apars[var]][t_ind,yi:yf,:]
            # Sample only a fraction of the points in the analysis
            vardata[var] = vardata[var][::sampling_interval,
                                        ::sampling_interval].flatten()
            if var=='CHI200': vardata[var] /= 10**6     # to 10^6 m^2 s^-1
            
    # Load the ob error variance file, if necessary
    if errvar_file is not None:
        with Dataset(errvar_file, 'r') as ncdata:
            errvar_la = ncdata.variables['lat'][:]
            errvar_lo = ncdata.variables['lon'][:]
            errvar_xy = ncdata.variables['errvar'][:]

    # Create Observation objects for each of the sampled locations
    observations = []
    for x in range(len(lonarr)):
        # Assign the ob error variance
        if errvar_file is None:
            oberr = err
        else:
            la_i = ut.nearest_ind(errvar_la, latarr[x])
            lo_i = ut.nearest_ind(errvar_lo, lonarr[x])
            oberr = errvar_xy[la_i, lo_i]
        # Create an Observation object for each variable
        for var in vrbls:
            obsject = Observation(value=vardata[var][x],obtype=var,time=atime,
                                  error=oberr,lat=latarr[x],lon=lonarr[x],
                                  localize_radius=locrad,description='GDAS')
            obsject.assimilate_this = True
            observations.append(obsject)
    return observations
    
def map_ob_density(lats, lons, title='ob_density', save_fig=True):
    """
    Plots the locations of the observations on a map.
    
    Requires:
    lats -----> A list of latitudes at each location. 
    lons -----> A list of longitudes at each location.
    title ----> A string containing the figure title.
                * NO white space!! This is used in the png filename too.
    save_fig -> A boolean indicating whether or not save the figure to a
                png file.
                
    Returns:
    Nothing! Either saves the figure or displays it onscreen.
    """
    if save_fig: plt.ioff()
    fig = plt.figure()
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9]) #l, b, w, h
    m = Basemap(projection='kav7',lon_0=180,resolution='l',ax=ax)
    x,y = m(lons,lats)
    m.scatter(x,y,3,marker='o',color='k')
    m.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    plt.title(title,fontsize=12)
    if save_fig:
        outdir = '/home/disk/p/njweber2/research/subseasonal/figures/efa/' + \
                 'ob_density'
        outfile = '{}/{}.png'.format(outdir, title)
        plt.savefig(outfile)
        plt.close()
    else:
        plt.show()
    
def write_state(state, idate, numdays, nlag, loc='GC', coarse_obs=False,
                    tunit='seconds since 1970-01-01', most_recent_mean=False,
                    obtimes=[1], prior=False, inflation=None, cov_file=None,
                    postfile=None, errvar=None, monte=False, mc=1,
                    mclat=0.):
    """
    Writes an observation-adjusted CFSv2 ensemble forecast to a netcdf file.
    
    Requires:
    state -> An Ensemble State object for the ensemble state.
    idate -> A datetime object specifying the ensemble initialization time.
    tunit -> The units for the time variable to be stored in the netcdf file.
             Default: seconds since 1970-01-01
    most_recent_mean -> A boolean object. If True, indicates that the
                        assimilation was done with the ensemble mean from only
                        the most recent 16 members (1 day lag).
    prior ------------> True = this is the prior state.
                        False = this is the post state
    inflation --------> Specifies whether any inflation was used on the prior.
    
    Returns:
    Nothing! Ensemble is written to a netcdf file.
    """
    # I/O information
    if prior:
        outdir = '/home/disk/vader2/njweber2/research/subseasonal/efa/' + \
                 'ensembles'
    else:
        outdir = '/home/disk/vader2/njweber2/research/subseasonal/efa/' + \
                 'adjusted_ensembles'
    if monte:
        outdir += '/monte_carlo'
    if cov_file is not None:
        outdir += '/historic_cov'
        checkdir(outdir)
    if inflation is not None:
        outdir += '/inflated_prior'
        checkdir(outdir)
    if most_recent_mean:
        outdir += '/most_recent_mean'
        checkdir(outdir)
    if coarse_obs:
        outdir += '/coarse_obs'
        checkdir(outdir)
    if len(obtimes) > 1:
        outdir += '/{}obtimes'.format(len(obtimes))
        checkdir(outdir)
    # Format the output file
    if loc in [None,False]:
        loc = 'NO'
    # If we're writing a single variable and the other variable already exists,
    # append to that file
    varlist = np.array(['Z500','CHI200'])
    # If we're only processing one variable
    if state.nvars() == 1:
        for var in varlist:
            checkfile = '{}/{}_{:%Y%m%d%H}_{}mem_{}days'.format(outdir,
                          varlist[varlist != var][0],idate,nlag*16,numdays)
            if prior: checkfile += '.nc'
            else: checkfile += '_{}loc.nc'.format(loc)
            if state.vars()[0] == var and os.path.isfile(checkfile):
                print('Appending to existing file!')
                # append to the existing file
                with Dataset(checkfile,'a') as dset:
                    dset.createVariable(var, 'f8', ('time','lat','lon','ens',))
                    dset.variables[var].units = ut.get_units(var)
                    dset.variables[var][:] = state[var].values
                # Rename the checkfile so the filename no longer specifies a 
                # single variable type
                newfile = '{}/{:%Y%m%d%H}_{}mem_{}days'.format(outdir,
                                                  idate, nlag*16, numdays)
                if prior: newfile += '.nc'
                else: newfile += '_{}loc.nc'.format(loc)
                os.system('mv {} {}'.format(checkfile,newfile))
                # ALL DONE!!
                return
        # If the checkfile does not exist, make a new file
        outfile = '{}/{}_{:%Y%m%d%H}_{}mem_{}days'.format(outdir,
                             state.vars()[0],idate,nlag*16,numdays)
    # If we're processing both Z500 and CHI200
    else:
        outfile = '{}/{:%Y%m%d%H}_{}mem_{}days'.format(outdir,idate,
                                                  nlag*16,numdays)
    if not prior: outfile += '_{}loc'.format(loc)
    if monte:
        outfile += '_mc{}_lat{}'.format(mc, int(mclat))
    if errvar is not None:
        outfile += '_{}oberror.nc'.format(int(errvar[varlist[0]]))
    else:
        outfile += '.nc'
    
    # Write ensemble forecast to netcdf
    with Dataset(outfile,'w') as dset:
        dset.createDimension('time',None)
        dset.createDimension('lat',state.ny())
        dset.createDimension('lon',state.nx())
        dset.createDimension('ens',state.nmems())
        dset.createVariable('time','i4',('time',))
        dset.createVariable('lat','f8',('lat',))
        dset.createVariable('lon','f8',('lon'))
        dset.createVariable('ens','i4',('ens',))
        dset.variables['time'].units = tunit
        dset.variables['lat'].units = 'degrees_north'
        dset.variables['lon'].units = 'degrees_east'
        dset.variables['ens'].units = 'member_number'
        dts = [datetime.utcfromtimestamp(t.astype(int)*1e-9) for \
               t in state.ensemble_times()]
        dset.variables['time'][:] = date2num(dts,tunit)
        dset.variables['lat'][:] = state['lat'].values[:,0]
        dset.variables['lon'][:] = state['lon'].values[0,:]
        dset.variables['ens'][:] = state['mem'].values
        for var in state.vars():
            dset.createVariable(var, 'f8', ('time','lat','lon','ens',))
            dset.variables[var].units = ut.get_units(var)
            dset.variables[var][:] = state[var].values
    return

def plot_covariances(lons, lats, covs, kal, mask, varlist, olat, olon, ndays,
                     cov_file=None):
    import matplotlib.pyplot as plt
    
    outdir1 = '/home/disk/p/njweber2/research/subseasonal/figures/efa/covs'
    outdir1 += '/{}N_{}E'.format(olat, olon)
    checkdir(outdir1)
    outdir2 = '/home/disk/p/njweber2/research/subseasonal/figures/efa/kal'
    outdir2 += '/{}N_{}E'.format(olat, olon)
    checkdir(outdir2)
    par = varlist[0]
    assert par=='Z500'
    
    #plot dem covariances
    covs = np.reshape(covs[:-1], (ndays, 181, 360))
    kal = np.reshape(kal[:-1], (ndays, 181, 360))
    mask = np.reshape(mask[:-1], (ndays, 181, 360))
    clevs = np.arange(-400., 401., 50.)
    klevs = np.arange(-3., 3.1, .5)
    if cov_file is not None: clevs *= 10.
    
    from mpl_toolkits.basemap import Basemap
    print('Making plots...')
    for lead in range(ndays):
        print('Lead: day-{}'.format(lead+1))
        print(' Covariances')
        # Make the figure and basemap
        plt.figure(figsize=(10,5))
        m = Basemap(projection='kav7',lon_0=180,resolution='l')
        m.drawcoastlines()
        x,y = np.meshgrid(lons, lats)
        X,Y = m(x,y)
        # Plot the data and make a point for the desired lat/lon
        cs = m.contourf(X, Y, covs[lead,:,:], cmap=plt.cm.RdBu_r, 
                        levels=clevs,extend='both')
        m.contourf(X, Y, mask[lead,:,:], colors='none', extend='both',
                   hatches=[None,'.'], levels=[0,1])
        X,Y = m(olon,olat) 
        m.scatter(X, Y, 20, marker='o', color='g')
        plt.colorbar(cs,label='covariance')
        plt.title('{} day-{} covariances'.format(par,lead+1))
        # How many points were updated?
        txt = '{} pts updated ({}%)'.format(int(np.sum(mask[lead,:,:])),
                                  int(100.*np.sum(mask[lead,:,:]/(181.*360.))))
        plt.figtext(0.95, 0.95, txt, horizontalalignment='right')
        # Save the figure
        if cov_file is not None:
            savefile = '{}/historic_{}_{}N_{}E_day{}.png'.format(outdir1,par,olat,
                                                          olon,lead+1)
        else:
            savefile = '{}/{}_{}N_{}E_day{}.png'.format(outdir1,par,olat,
                                                          olon,lead+1)
        plt.savefig(savefile)
        plt.close()    
        
        #===
        
        print(' Kalman Gain')
        # Make the figure and basemap
        plt.figure(figsize=(10,5))
        m = Basemap(projection='kav7',lon_0=180,resolution='l')
        m.drawcoastlines()
        x,y = np.meshgrid(lons, lats)
        X,Y = m(x,y)
        # Plot the data and make a point for the desired lat/lon
        cs = m.contourf(X, Y, kal[lead,:,:], cmap=plt.cm.RdBu_r, 
                        levels=klevs,extend='both')
        m.contourf(X, Y, mask[lead,:,:], colors='none', extend='both',
                   hatches=[None,'.'], levels=[0,1])
        X,Y = m(olon,olat) 
        m.scatter(X, Y, 20, marker='o', color='g')
        plt.colorbar(cs,label='Kalman gain')
        plt.title('{} day-{} Kalman gain'.format(par,lead+1))
        # How many points were updated?
        txt = '{} pts updated ({}%)'.format(int(np.sum(mask[lead,:,:])),
                                  int(100.*np.sum(mask[lead,:,:]/(181.*360.))))
        plt.figtext(0.95, 0.95, txt, horizontalalignment='right')
        # Save the figure
        if cov_file is not None:
            savefile = '{}/historic_{}_{}N_{}E_day{}.png'.format(outdir2,par,olat,
                                                          olon,lead+1)
        else:
            savefile = '{}/{}_{}N_{}E_day{}.png'.format(outdir2,par,olat,
                                                          olon,lead+1)
        plt.savefig(savefile)
        plt.close()    
    print('TOTAL PTS UPDATED: {} ({}%)'.format(int(np.sum(mask)),
                                               int(100.*np.sum(mask)/np.size(mask))))
    
###############################################################################

if __name__ == '__main__':
    from random import shuffle
    #=== Edit these ===========================================================
    nlagdays = np.array([4])     # num of days to lag ensemble
    idates = [datetime(2015,12,1,0)]      # ensemble initialization date
    ndays = 21                         # num of days in forecast  
    ens_ncdf = True                    # load pre-processed ensemble from netcdf?
    varlist = ['Z500']                 # what variables to process?
    ob_interval = 10                   # select every [?] analysis points for obs
    cutofflat = 70.                    # don't sample analysis above this lat
    localization = 'GC2000'               # what kind of ob localization?
    use_most_recent_mean_state = False # use only the 16-mem ens mean?
    historic_covariances = False       # use historic covariances rather than ens.
    plot_covs = False                   # plot the covariances instead of 
    #                                    writing an adjusted forecast?      
    monteCarlo = 50                    # number of Monte-Carlo iterations
    #                                    (randomly shuffle obs list)                             
    one_ob = False                      # only assimilate one ob?
    if one_ob:
        olat = 60
        olon = 0
        
    # dates (forecasts) used for inflating the prior
    caldates = [datetime(2015,12,1,0),datetime(2015,12,16,0),datetime(2016,1,1,0),\
                        datetime(2016,1,16,0),datetime(2016,2,1,0),datetime(2015,12,8,0),\
                        datetime(2015,12,24,0),datetime(2016,1,8,0),datetime(2016,1,24,0),\
                        datetime(2016,2,8,0)]
    caldates = []
    #==========================================================================
    for mc in np.arange(monteCarlo)+1:
        print('!!!!!!! M-C {} !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'.format(mc))
        
        shuff_inds = np.arange(540)
        shuffle(shuff_inds)
        for errvar in [{'Z500':1.},{'Z500':10.},{'Z500':100.}]:
            print('-----ERRVAR: {} ----------'.format(errvar['Z500']))
            if historic_covariances:
                # find the proper file if we're using historic covariances
                cov_file = '/home/disk/vader2/njweber2/research/subseasonal/efa/'+\
                           'historic_covariances/{}_{}_covariances_{}days.nc'
                cov_file = cov_file.format(varlist[0], idates[0].strftime('%b'), ndays)
                cov_file = '/home/disk/vader2/njweber2/research/subseasonal/efa/'+\
                           'historic_covariances/{}_covariances_{}days.nc'
                cov_file = cov_file.format(varlist[0], ndays)
            else:
                cov_file = None
            
            ### MAIN BODY:
            if len(caldates)>0:
                inflate = True
            else:
                inflate = False
            
            if inflate:    
                print('\n(0) Fetching calibration/inflation data...')
                cal_dict = error_vs_spread(caldates, ndays, ens_sizes=nlagdays*16,
                                           vrbls=varlist, loc=localization,  
                                           recent_mean=use_most_recent_mean_state,
                                           coarse_obs=False, return_weights = True)
                    
            for idate in idates:
                print('\n##################### {:%Y%m%d} ########################'.format(idate))
                fdate = idate + timedelta(days=ndays)
                assim_times = [idate]
                for nlagday in nlagdays:
                    print('\n======= {}-mem ensemble =================='.format(nlagday*16))
            
                #==== BUILD THE STATE =========================================================
                    print('\n(1) Building state...')
                    state = build_state(nlagday, idate, fdate, vrbls=varlist,
                                        fromfile=ens_ncdf, oberr=errvar)
                    if use_most_recent_mean_state:
                        most_recent_state = build_state(1,idate,fdate,vrbls=varlist,
                                                        fromfile=ens_ncdf)
                    print('Dimensions/variables:')
                    for key in state.variables.keys():
                        print('{}{}'.format(key.ljust(10),np.shape(state.variables[key])))
                      
                    if inflate:  
                        print('Selecting appropriate inflation...')
                        cal_id = '{}{}'.format(int(nlagday*16), varlist[0])
                        inflations = {'validtime' : cal_dict[cal_id]}
                    else: 
                        inflations = None
                        
                #==== LOAD THE OBSERVATIONS ===================================================
                    print('\n(2) Getting observations from GDAS analyses...')
                    observations = []
                    for atime in assim_times:
                        print('{:%Y-%m-%d %H:00}'.format(atime))
                        newobs = get_gdas_obs(atime, vrbls=varlist, loc=localization,
                                        sampling_interval=ob_interval, 
                                        cutofflat=cutofflat, err=errvar)
                        observations += newobs
                    # shuffle the obs!
                    if monteCarlo > 1:
                        observations = list(np.array(observations)[shuff_inds])
                    mclat = observations[0].lat
                    print('Number of obs: {}'.format(len(observations)))
                    
                    # If we're only interested in a single ob, find it!
                    if one_ob:
                        olats = np.array([o.lat for o in observations])
                        olons = np.array([o.lon for o in observations])
                        assert olat in olats
                        assert olon in olons
                        lons = state['lon'][0,:]
                        lats = state['lat'][:,0]
                        o_ind = np.where((olats==olat)*(olons==olon))[0][0]
                        observations = [observations[o_ind]]
                        print(' Only assimilating the ob at {}N {}E'.format(olat,olon))
                 
                #==== ASSIMILATE THE OBS ======================================================
                    print('\n(3) Assimilating the observations...')
                    assimilator = EnSRF(state, observations, nproc=1, verbose=True, 
                                        loc=localization, inflation=inflations)
                    start = time.time()
                    if plot_covs:
                        assert one_ob
                        covs, kal, mask = assimilator.update(return_covs=True,
                                                             cov_file=cov_file)
                        
                        plot_covariances(lons, lats, covs, kal, mask, varlist, olat, 
                                         olon, ndays, cov_file=cov_file)    
                    else:
                        if use_most_recent_mean_state:
                            assimilator2 = EnSRF(most_recent_state, observations, nproc=1, 
                                             verbose=True, loc=localization)
                            xam, Xap = assimilator2.format_prior_state()
                            del assimilator2
                            post_state, post_obs, updated_pts = assimilator.update(mean_state=xam)
                        else:
                            post_state, post_obs, updated_pts = assimilator.update(cov_file=cov_file)
                        end = time.time()
                        print('Elapsed assimilation time: {:.2f} min'.format((end-start)/60.))
                        
                        # Save the updated points array to a text file
                        tfile = '/home/disk/vader2/njweber2/research/subseasonal/efa/updated_pts'
                        tfile += '/{}_{:%Y%m%d}_{}mem_{}loc.txt'.format(varlist[0],idate,
                                                                    nlagday*16, localization)
                        np.savetxt(tfile, updated_pts)             
                        
                        if monteCarlo > 1:
                            monte = True
                        else:
                            monte = False
                            
                        # Let's write the posterior state to a netcdf
                        print('(3.1) Writing posterior to netcdf...')
                        write_state(post_state, idate, ndays, nlagday, loc=localization, 
                                        most_recent_mean=use_most_recent_mean_state,
                                        coarse_obs=(ob_interval==30), obtimes=assim_times,
                                        inflation=inflations, cov_file=cov_file,
                                        errvar=errvar,
                                        monte=monte, mc=mc, mclat=mclat)   
        #            print('(3.1) Writing inflated prior to netcdf...')
        #            write_state(prior, idate, ndays, nlagday, loc=localization, 
        #                            most_recent_mean=use_most_recent_mean_state,
        #                            coarse_obs=(ob_interval==30), obtimes=assim_times,
        #                            prior=True, inflation=inflations)  
        
    print('DONE!')
