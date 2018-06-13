# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 11:39:33 2015

@author: njweber2
"""

import numpy as np
import subprocess
from subprocess import check_output, CalledProcessError
from netCDF4 import Dataset, num2date, date2num
import datetime
import matplotlib.pyplot as plt
import re

# Global variables/data
ninofile = '/home/disk/vader2/njweber2/research/subseasonal/data/'+\
           'enso/nino3.4_1982-2008.txt'
ndata = np.loadtxt(ninofile,usecols=(0,1,-1)).T
mjofile = '/home/disk/vader2/njweber2/research/subseasonal/data/'+\
           'mjo/RMM1RMM2_1982-2008.txt'
mdata = np.loadtxt(mjofile,usecols=(0,1,2,6)).T
sstfile = '/home/disk/vader2/njweber2/research/subseasonal/data/'+\
          'sst/cfsv2_aveSSTanom_timeseries.nc'
sdata = Dataset(sstfile)

##################### FUNCTIONS ###############################################
def checkdir(directory):
    """
    Checks to see if the directory exists. If not, the directory is created.
    """
    import os
    if not os.path.isdir(directory):
        os.system('mkdir {}'.format(directory))

def nearest_ind(array,value):
    """
    Finds the nearest index of the given value in the given array.
    """
    return int((np.abs(array-value)).argmin())
    
def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta
    
def check_input_var(var):
    """
    checks to see if the input string is a valid meteorological parameter,
    then returns the parameter string in the format listed in the CFSv2
    grib filenames.
    """
    if var=='T2M': return 'tmp2m'
    elif var=='Z500': return 'z500'
    elif var=='Z200': return 'z200'
    elif var=='Z850': return 'z850'
    elif var=='T850': return 't850'
    elif var=='MSLP': return 'prmsl'
    elif var=='PRATE': return 'prate'
    elif var=='SST': return 'tmpsfc'
    elif var=='CHI200': return 'chi200'
    elif var=='PWAT': return 'pwat'
    elif var in ['U200', 'V200']: return 'wnd200'
    elif var=='U850': return 'wnd850'
    elif var=='OLR': return 'ulwtoa'
    else: raise RuntimeError('Invalid atmospheric parameter')
    
def get_units(par):
    """
    Returns the units of the given geophysical parameter
    """
    if par in ['Z500','Z700','Z850','Z925','Z1000']: return 'meters'
    elif par in ['T500','T700','T850','T925','T1000','T2M','t2m']: return 'Kelvin'
    elif par in ['RH500','RH700','RH850','RH925','RH1000','RH2M']: return 'percent'
    elif par=='MSLP': return 'pascals',
    elif par =='ALT': return 'hPa',
    elif par=='CHI200': return '10$^{6}$ m$^{2}$ s$^{-1}$'
    elif par=='PRATE': return 'mm/day'
    elif par in ['PWAT','P6HR','tcw', 'TCW']: return 'mm'
    elif par in ['U500','V500','U700','V700','U850','V850','U925','V925', \
                 'U1000','V1000','U10M','V10M']: return 'm/s'
    elif par=='OLR': return 'W m$^{-2}$'
    else: raise IOError('Invalid field: {}'.format(par))
    
def check_fcst(fcsts,times,des_time,count=False):
    """
    If the desired forecast lead time is in the list of times, then the 
    corresponding forecast grid is returned. Otherwise, an array of nan's.
    """
    if (times == des_time).any():
        if count:
            return 1
        else:
            arr = fcsts[np.where(times==des_time)[0][0],:,:]
            return arr
    else:
        if count:
            return 0
        else:
            arr = np.empty(np.shape(fcsts[0,:,:]))
            arr[:] = np.nan
            return arr
    
def load_cfs_forecast(cfs_file,numdays):
    """
    Returns the first [numdays] days of a CFSv2 forecast initialized on the 
    specified date.
    """
    import pygrib as pyg
    try:
        grbs = pyg.open(cfs_file)
        grb = grbs.select(forecastTime=range(numdays*24+1))
        grbs.close()
        return np.array([g.values for g in grb])
    except subprocess.CalledProcessError:
        raise IOError()
        
def load_cfs_ensemble(fcst_dir,init_date,par,numdays,bli,tri,
                      intervals,des_inds,returnpar='mean',daily=True):
    """
    Returns the the ensemble mean averaged forecast for the first [numdays] 
    days of the forecast. [init_date] is the forecast initialization date and 
    the ensemble (4 members) includes the 00Z, 06Z, 12Z, and 18Z runs. 
    [returnpar] specifies whether to return the ensemble mean, ensemble spread,
    ensemble stdv, or the individual members. The averaging intervals are 
    specified by [intervals], which is a list of integers specifiying the 
    number of weeks to average the forecast (0weeks = 1day).
    """
    import pygrib as pyg
    datestr = init_date.strftime('%Y%m%d')
    ls_fils = '{}/{}.{}*.grb2'.format(fcst_dir,check_input_var(par),datestr)
    fils = check_output(['ls -1a '+ls_fils],shell=True).split()
#    if not len(fils)==4:
#        raise IOError('There are only {} files for this day!!'.format(len(fils)))
    if returnpar=='mem1':
        fils = [fils[0]]
    mems = len(fils)
    ensemble = np.empty((mems,len(intervals),len(des_inds),bli[0]-tri[0],tri[1]-bli[1]))
    m = 0
    for cfs_file in fils:
        try:
            hh = int(cfs_file[-12:-10])
            print('   member{}: {}Z'.format(m+1,cfs_file[-12:-10]))
            grbs = pyg.open(cfs_file)
            if par in ['U200', 'U850']:
                if hh == 18:
                    grb = grbs.select(name='U component of wind',forecastTime=range(numdays*24+1))
                else:
                    grb = grbs.select(name='U component of wind',forecastTime=range((numdays+1)*24+1))
            elif par in ['V200', 'V850']:
                if hh == 18:
                    grb = grbs.select(name='V component of wind',forecastTime=range(numdays*24+1))
                else:
                    grb = grbs.select(name='V component of wind',forecastTime=range((numdays+1)*24+1))
            else:
                if hh == 18:
                    grb = grbs.select(forecastTime=range(numdays*24+1))
                else:
                    grb = grbs.select(forecastTime=range((numdays+1)*24+1))
            grbs.close()
            if daily:
                # Take daily average
                fcst_daily, ftims, skp = daily_ave(np.array([g.values for g in grb]), hh)
            if par == 'PRATE':
                    fcst_daily *= 60*60*24
            # Loop over all the different averaging intervals (0-6 weeks)
            for n in range(len(intervals)):
                numints = intervals[n]
                # Compute the mean for the interval at each fcst time
                if numints == 0:
                    fcst_mean = fcst_daily[des_inds,:,:]
                else:
                    ave_interval = numints*7
                    #average over the time axis
                    fcst_mean = np.array([np.nanmean(fcst_daily[j:(j+ave_interval),\
                               :,:],axis=0) for j in des_inds])
                
                ensemble[m,n,:,:,:] =  fcst_mean[:,tri[0]:bli[0],bli[1]:tri[1]]
        except ValueError:
            print('Houston, we have a ValueError')
            raise ValueError('There was a ValueError')
            ensemble[m,n,:,:,:] = np.nan
        except subprocess.CalledProcessError:
            raise IOError('There was a CalledProcessError!!')
        m += 1
    if returnpar=='members':
        return ensemble[0,:,:],ensemble[1,:,:],ensemble[2,:,:],ensemble[3,:,:]
    elif returnpar in ['mean','mem1']:
        return np.nanmean(ensemble,axis=0)
    elif returnpar=='stdv':
        return np.nanstd(ensemble,axis=0)
    elif returnpar=='spread':
        return np.amax(ensemble,axis=0) - np.amin(ensemble,axis=0)
    else:
        raise IOError('Invalid return parameter specified.')
        
def load_cfs_member(fcst_dir,init_date,par,numdays,bli,tri,
                      intervals,des_inds,mem):
    import pygrib as pyg
    """
    Returns the member [mem] forecast for the first [numdays] 
    days of the forecast. [init_date] is the forecast initialization date and 
    the ensemble (4 members) includes the 00Z, 06Z, 12Z, and 18Z runs. 
    [returnpar] specifies whether to return the ensemble mean, ensemble spread,
    ensemble stdv, or the individual members. The averaging intervals are 
    specified by [intervals], which is a list of integers specifiying the 
    number of weeks to average the forecast (0weeks = 1day).
    """
    ave_fcsts = np.empty((len(intervals),len(des_inds),bli[0]-tri[0],tri[1]-bli[1]))
    datestr = init_date.strftime('%Y%m%d')
    ls_fils = '{}/{}.{}*.grb2'.format(fcst_dir,check_input_var(par),datestr)
    fils = check_output(['ls -1a '+ls_fils],shell=True).split()
    if not len(fils)==4:
        raise IOError('There are only {} files for this day!!'.format(len(fils)))
    try:
        # Load the desired member
        cfs_file = fils[mem-1]
        hh = int(cfs_file[-12:-10])
        grbs = pyg.open(cfs_file)
        if par in ['U200', 'U850']:
            if hh == 18:
                grb = grbs.select(name='U component of wind',forecastTime=range(numdays*24+1))
            else:
                grb = grbs.select(name='U component of wind',forecastTime=range((numdays+1)*24+1))
        elif par in ['V200', 'V850']:
                if hh == 18:
                    grb = grbs.select(name='V component of wind',forecastTime=range(numdays*24+1))
                else:
                    grb = grbs.select(name='V component of wind',forecastTime=range((numdays+1)*24+1))
        else:
            if hh == 18:
                grb = grbs.select(forecastTime=range(numdays*24+1))
            else:
                grb = grbs.select(forecastTime=range((numdays+1)*24+1))
        grbs.close()
        # Take daily average
        fcst_daily, ftims, skp = daily_ave(np.array([g.values for g in grb]),hh)
        if par == 'PRATE':
                fcst_daily *= 60*60*24
        # Loop over all the different averaging intervals (0-6 weeks)
        for n in range(len(intervals)):
            numints = intervals[n]
            # Compute the mean for the interval at each fcst time
            if numints == 0:
                fcst_mean = fcst_daily[des_inds,:,:]
            else:
                ave_interval = numints*7
                #average over the time axis
                fcst_mean = np.array([np.nanmean(fcst_daily[j:(j+ave_interval),\
                           :,:],axis=0) for j in des_inds])
            ave_fcsts[n,:,:,:] =  fcst_mean[:,tri[0]:bli[0],bli[1]:tri[1]]
        
    except subprocess.CalledProcessError:
        raise IOError('There was a CalledProcessError!!')
    return ave_fcsts
        
def load_cfs_forecasts(direc,year,datestr,par,fc_len):
    """
    Returns the most recent [fc_len]/4 days of CFSv2 forecasts valid on the 
    specified date.
    """
    variable = check_input_var(par)
    nc_var = check_input_var(par,ncvar=True)
    # List the forecasts valid at desired date
    lscommand = 'ls -1a '+direc+'/'+year+'*/'+variable+'*v.'+datestr+'*'
    try:
        cfs_files = np.array(check_output([lscommand],shell=True).split())
    except CalledProcessError:
        cfs_files = np.array([])
    # If within the first 45 days of the year, load previous year's fcsts
    dt = datetime.datetime(2000,int(datestr[4:6]),int(datestr[6:8]),0)
    if dt < datetime.datetime(2000,2,15,0) and year != '1999':
        lscommand = 'ls -1a '+direc+'/'+str(int(year)-1)+'*/'+variable+\
                    '*v.'+datestr+'*'
        prevyear = check_output([lscommand],shell=True).split()
        cfs_files = np.append(prevyear,cfs_files)
    # Load the data
    cfs_data = np.array([Dataset(f).variables[nc_var][0,::-1,:] \
                         for f in cfs_files[::-1]])     # ^reverse lat
    cfs_times = np.array([int(f[-7:-3]) for f in cfs_files[::-1]])
    # Use only the most recent [fc_len] forecasts; if any of those are
    # missing, fill with NaN array
    desired_times = (np.arange(fc_len)+1)*6
    cfs_fulldata = np.array([check_fcst(cfs_data,cfs_times,t) \
                    for t in desired_times])
    counts = [check_fcst(cfs_data,cfs_times,t,count=True) \
              for t in desired_times]
    return cfs_fulldata, counts


def load_analysis(direc,year,par,bl=None,tr=None):
    """
    Returns an entire year's worth of CFS analysis data for the desired
    parameter over a subdomain specified by [bl] and [tr].
    """
    cdf_file = direc+'/'+par+'.'+year+'.nc'
    if par=='CHI200':
        par = 'chi200'
    elif par=='PRATE':
        par = 'prate'
    elif par=='PWAT':
        par = 'pwat'
    ncdata = Dataset(cdf_file,mode='r')
    lats = ncdata.variables['latitude'][:]
    lons = ncdata.variables['longitude'][:]
    times = ncdata.variables['time'][:]
    data = ncdata.variables[par][:]
    ncdata.close()
    if bl is not None and tr is not None:
        return times, lats[tr[0]:bl[0]], lons[bl[1]:tr[1]],\
               data[:,tr[0]:bl[0],bl[1]:tr[1]]
    else:
        return times, lats, lons, data
        
def yearly_analysis(direc,year,var,bli,tri,nday,only00z=True,prevyear=False):
    """
    Returns the valid dates (datetime), lats, lons, and analyses for a year's
    worth (plus part of the following year) of analyses of a given paramer.
    """
    print('Loading yearly analysis data...')
    if var=='totallynotPRATE':
        perday = 1
    else:
        perday = 4
    # Load this year's analyses in N. Amer. subdomain
    anl_t,anl_lat,anl_lon,anl_data = load_analysis(direc,year,var,
                                                   bli,tri)
    if prevyear:
        # Append the last [ndays] of previous year's analyses
        t2,anl_lat,anl_lon,data2 = load_analysis(direc,str(int(year)-1),
                                                 var,bli,tri)
        anl_data = np.concatenate((data2[-(nday*perday):,:,:],anl_data),axis=0)
        anl_t = np.concatenate((t2[-(nday*perday):],anl_t))
    # Append the first [ndays] of next year's analyses
    t2,anl_lat,anl_lon,data2 = load_analysis(direc,str(int(year)+1),
                                             var,bli,tri)
    anl_data = np.concatenate((anl_data,data2[:(nday*perday),:,:]),axis=0)
    anl_t = np.concatenate((anl_t,t2[:(nday*perday)]))
    # Put the analysis times in datetime format
    if var=='totallynotPRATE':
        anl_t = np.array([num2date(d,'days since 1990-01-01').\
                          replace(hour=0) for d in anl_t])
    else:
        tims = np.array([num2date(n,'hours since 1800-01-01') for \
                          n in anl_t])
        # Take only the 00Z times, if appropriate
        if only00z:
            is00z = np.array(np.where([d.hour==0 for d in tims]))[0,:]
            anl_data = anl_data[is00z,:,:]
            anl_t = tims[is00z]
        else:
            anl_t = tims
    return anl_t, anl_lat, anl_lon, anl_data
    
def daily_ave(fcst,init_h):
    """
    Takes in a six-hourly forecast and its corresponding forecast times (in
    hours) and returns a daily-average forecast with forecast times
    """
    # Strip off the first, incomplete forecast day
    l = len(fcst[:,0,0])
    if init_h==0: 
        #remove first 06Z, 12Z, 18Z validation times
        fcst = np.delete(fcst,[0,1,2,l-1],axis=0)
        first_fh = 24
        skipped_first = True
    elif init_h==6: 
        #remove first 12Z, 18Z validation times
        fcst = np.delete(fcst,[0,1,l-1,l-2],axis=0)
        first_fh = 18
        skipped_first = True
    elif init_h==12: 
        #remove first 18Z validation time
        fcst = np.delete(fcst,[0,l-1,l-2,l-3],axis=0)
        first_fh = 12
        skipped_first = True
    elif init_h==18: 
        #first valid. time is 00Z; no removal necessary
        first_fh = 6
        skipped_first = False
    # Average the forecast daily
    nlat = len(fcst[0,:,0]); nlon = len(fcst[0,0,:]) 
    daily_fcst = np.nanmean(fcst.reshape(-1,4,nlat,nlon),axis=1)
    ftimes = np.arange(len(daily_fcst[:,0,0]))*24 + first_fh
    return daily_fcst, ftimes, skipped_first
    
def get_fcst_inds(all_t,daily_t):
    """
    Takes in an array of all the 6-hourly forecast times in a CFS forecast and
    returns the indices within that array corresponding to the given array
    of daily forecast times
    """
    ind = [np.where(all_t==t)[0][0] for t in daily_t]
    return np.array(ind)
    
def isnino(dt,highthresh=False):
    """
    Returns a string indicating the phase of the Nino3.4 index for the
    given date.
    """
    yr = int(dt.strftime('%Y'))
    mon = int(dt.strftime('%m'))
    ninoindex = ndata[2,:][(ndata[0,:]==yr)*(ndata[1,:]==mon)]
    thresh = 0.5
    if highthresh: thresh *= 2.
    if ninoindex[0] > thresh:
        return 'ElNino'
    elif ninoindex[0] < -thresh:
        return 'LaNina'
    else:
        return 'Neither'
    
def ismjo(dt,highthresh=False):
    """
    Returns True if the given date is during an MJO event (amplitude > 1).
    """
    yr = int(dt.strftime('%Y'))
    mon = int(dt.strftime('%m'))
    day = int(dt.strftime('%d'))
    mjoamp = mdata[3,:][(mdata[0,:]==yr)*(mdata[1,:]==mon)*(mdata[2,:]==day)]
    thresh = 1.0
    if highthresh: thresh *= 1.5
    return mjoamp[0] > thresh
    
def sstanom(dt,highthresh=False):
    """
    Returns a list of strings indicating the "types" of SST anomalies, if any,
    that are present on the given date. (SST1w, SST1c, SST2w, and SST2c)
    """
    # These will be timeseries of the averaged SST anomalies, loaded from a
    # file in the Global Variables section of this script
    sst1 = sdata.variables['SST1'][:]
    sst2 = sdata.variables['SST2'][:]
    sst3 = sdata.variables['SST3'][:]
#    # Choose a threshold value: 90th percentile?
#    perc = 90
#    lo1 = np.percentile(sst1,100-perc)
#    hi1 = np.percentile(sst1,perc)
#    lo2 = np.percentile(sst2,100-perc)
#    hi2 = np.percentile(sst2,perc)
    hi1 = 0.5; lo1 = -0.5
    hi2 = 0.5; lo2 = -0.5
    hi3 = 0.5; lo3 = -0.5
    if highthresh:
        hi1 *= 2.; lo1 *= 2.
        hi2 *= 2.; lo2 *= 2.
        hi3 *= 2.; lo3 *= 2.
    # Check if the SSTs on the current date exceed these thresholds
    t_ind = np.where(sdata.variables['validtime'][:]==date2num(dt,
                                         'hours since 1800-01-01'))[0][0]
    t1 = sst1[t_ind]
    t2 = sst2[t_ind]
    t3 = sst3[t_ind]
    return_list = []
    if t1 >= hi1:
        return_list.append('SST1w')
    if t1 <= lo1:
        return_list.append('SST1c')
    if t2 >= hi2:
        return_list.append('SST2w')
    if t2 <= lo2:
        return_list.append('SST2c')
    if t3 >= hi3:
        return_list.append('SST3w')
    if t3 <= lo3:
        return_list.append('SST3c')
    return np.array(return_list)
  
def draw_mercator():
    from mpl_toolkits.basemap import Basemap
    m = Basemap(projection='merc',llcrnrlat=16,urcrnrlat=55,\
            llcrnrlon=220,urcrnrlon=298,lat_ts=20,resolution='l')
    m.drawmeridians(np.arange(0,360,10),labels=[0,0,0,1])
    m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0])
    m.drawstates(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    return m
    
def map_pcolor(m,lons,lats,var,minlev,maxlev,title,cmap=plt.cm.jet):
    X,Y = np.meshgrid(lons,lats)
    m.pcolor(X,Y,var,cmap=cmap,latlon=True,vmin=minlev,vmax=maxlev)
    plt.title(title)
    
def map_contourf(m,lons,lats,var,minlev,maxlev,title,cmap=plt.cm.jet):
    X,Y = np.meshgrid(lons,lats)
    x,y = m(X,Y)
    levs = np.arange(minlev,maxlev,np.round((maxlev-minlev)/10))
    m.contourf(x,y,var,cmap=cmap,levels=levs,extend='both')
    plt.title(title)
    
def bandpass_filter(par, data, init_dates, leads, continuous_data=False, 
                    hi=20., lo=100., returnmeans=False):
    from scipy.fftpack import rfft, irfft, fftfreq
    assert(len(init_dates))==np.shape(data)[0]
    if continuous_data:
        w_min = 1./lo
        w_max = 1./hi
        dt = (init_dates[1]-init_dates[0]).days
        
        w = fftfreq(np.shape(data)[0], d=dt)
        
        data_hat = rfft(data, axis=0)
        data_hat[(w<w_min),:] = 0   # low pass
        data_hat[(w>w_max),:] = 0   # high pass
        data = irfft(data_hat, axis=0)    
    else:
        # ASSUME data is already weekly-averaged, so simply subtract the anomaly
        # from the past 120d to filter out the low frequencies
        meanfile = '/home/disk/vader2/njweber2/research/subseasonal/data/' + \
                   'all_forecasts/{}_120dMEANS_1982-2008.nc'.format(par)
        with Dataset(meanfile, 'r') as ncdata:
            mdates = ncdata.variables['init_dates']
            mdates = num2date(mdates[:], mdates.units)
            mleads = ncdata.variables['ftime'][:]
            l_inds = np.array([np.where(mleads==l)[0][0] for l in leads])
            means = ncdata.variables['forecasts'][:,l_inds, :, :]
            if len(leads)>1: assert np.shape(means)[1]==len(leads)
        if not returnmeans: 
            assert(len(leads))==np.shape(data)[1]
            assert np.shape(data)==np.shape(means)
            data -= means
    if returnmeans:
        return means.squeeze()
    else:
        return data
        
def plot_eq_land(ax, y, plus=False, no_SA=False):
    """
    Plots cute little lines on the bottom of your hovmoller diagram (or any
    plot with longitude as the x-axis) that indicate where equatorial land
    masses are.
    
    Requires:
    ax --> The axis object for your plot.
    y ---> The lower y-limit for your plot.
    """
    # Add 2 to y, so the lines are plotted *just* off the bottom
    if plus: yt = y + 2
    else: yt = y - 2
    # Make labels, locations, and colors, for the three land areas
    #         Africa,     Maritime Cont., South America
    labels = ['AF',      'MC',            'SA']
    lons =   [[5.,45.],  [95.,150.],     [280.,320.]]
    cols =   ['dodgerblue', 'firebrick', 'darkorchid']
    if no_SA:
        labels = labels[:2]; lons = lons[:2]; cols=cols[:2]
    # Loop through each one
    for lab, lon, col in zip(labels, lons, cols):
        # Plot the line
        ax.plot(lon, [y,y], color=col, linestyle='-', linewidth=10)
        # Plot the label just above the line
        ax.text(np.mean(lon), yt, lab, ha='center', va='bottom', color=col,
                fontsize=18, weight='bold')
    return

def interp_2d_grid(lon1, lat1, data, lon2, lat2):
    """
    Interpolates CFS analysis to CFS forecast resolution
    """
    from scipy import interpolate
    f = interpolate.interp2d(lon1, -1*lat1, data, kind='linear')
    interp_data = f(lon2, -1*lat2)
    return interp_data
    
def keep_like_elements(lists):
    """
    Sorts through lists of items and returns the indices of the items that are
    in all of the lists.
    Requires:
    lists ----> a list of lists or numpy arrays
    Returns:
    indices --> a list of indices corresponding to each list in [lists]
    """
    indices = []
    for list_i in lists:
        inAllLists = [all([item in eachlist for eachlist in lists]) for item in list_i]
        indices.append(np.where(inAllLists)[0])
    return indices

#--------------ADDED BY SPENCER------------------------------------------------
# These are used to sort numbers correctly, according to humans.
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
