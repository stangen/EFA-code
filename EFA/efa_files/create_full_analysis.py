 # -*- coding: utf-8 -*-
"""
Created on Wed May 11 14:19:56 2016

@author: njweber2/stangen
"""

from netCDF4 import Dataset, num2date, date2num
from datetime import datetime, timedelta
from subprocess import check_output
import numpy as np
import efa_files.cfs_utilities_st as ut
from efa_xray.state.ensemble import EnsembleState

def create_full_analysis(vrbls=['Z500'], writenc=True):
    """
    Populates and returns a CFSv2 global ensemble forecast using the
    archived operational forecasts on vader.atmos.washington.edu.
    
    Requires:    
    
    month -> string of month we want to retrieve.
    year  -> string of year we want to retrieve.
    day   -> string of day we want to retrieve.
    hour  -> string (00 or 12) we want to retrieve.
    
    vrbls -> A list of variables we want to retrieve.

    writenc -> A boolean object. If True, the ensemble data will be written
               out to a netcdf. If False, the ensemble data will simply be
               returned.
               
    Returns:
    statecls -> An ensemble state object (see EnsembleState class)
                *only if writenc==False
    """
    
    
    
    # Here we have a dictionary translating some generic variable names to
    # what they correspond to in the CFSv2 netcdf files
    vardict = {'Z500' : 'HGT_500mb',
               'gh' : 'Z500',
               't2m' : 'T2M',
               'tcw' : 'TCW',
               'time' : 'time',
               'lat' : 'latitude',
               'lon' : 'longitude',
               '10' : 'oct',
               '11' : 'nov',
               '12' : 'dec',
               '01' : 'jan',
               '02' : 'feb',
               '03' : 'mar'
               }
    outfile = '/home/disk/hot/stangen/Documents/ensembles/analysis/combined/oct-mar.nc' 
    
    # directories for the months, initialize some lists
    months = ['oct2016', 'nov2016', 'dec2016', 'jan2017', 'feb2017', 'mar2017']
    memfiles = []
    each_time_length = list()
    state = []
    
    # try is because may not have all months yet
    try:
        for mon in months:
            indir = '/home/disk/hot/stangen/Documents/ensembles/analysis/rawmonths/%s/*' % (mon) 
            # Get a list of filenames (each file is a different member)
            command = 'ls -1a {}'.format(indir)
            memfiles.extend(list(reversed(check_output([command],shell=True).split())))
    except:
        pass
    
    # more initializations
    ntimes = 0
    priorens = '                                                                                  '
    t2 = -1
    
    # Count the number of times for allocation of state later
    # try blocks are to filter out dates when TIGGE didn't have any data, probably not used in analysis data.
    print('Counting times for allocation')
    for t, times in enumerate(memfiles):
        times = times.decode('utf-8')
        # to not double-count times if from both sfc and pl
        if times[0:80] != priorens[0:80]:
            t2 += 1
            try:
                with Dataset(times, 'r') as ncdata:                   
                    #keep track of number of members of each ensemble
                    each_time_length.append(len(ncdata.variables['time']))
                    #keep track of total length of all ensembles
                    ntimes = ntimes + len(ncdata.variables['time'])
                    print('times from {}: {}'.format(months[t2],each_time_length[t2]))
                    print('total times: {}'.format(ntimes))
                    #print(each_mem_length)
            except:
                print('Bad {} ensembles- not counting ensemble members'.format(months[t2]))

                pass
        priorens = times
    
    #even more initializations
    timerange_sum = 0
    bad_times_range = 0    
    priorens = '                                                                                  '
    del_mems = []
    mnum2 = -1
    ftimes = []
    valid_times = np.zeros(ntimes)
    
    #mem naming is leftover from ensemble script, this loop adds data from each file
    for mnum, mem in enumerate(memfiles):
        mem = mem.decode('utf-8')
        
        try:
            # Read the netcdf
            with Dataset(mem,'r') as ncdata:
                # Find the indices corresponding to the start and end times
                tunit = ncdata.variables[vardict['time']].units
                ftimes = num2date(ncdata.variables[vardict['time']][:],tunit)
    
                # If this is the first month, calculate how large the state array
                # needs to be and allocate. Also set up the metadata.               
                if mnum == 0:
                    nvars = len(vrbls)
                    nlats = len(ncdata.dimensions[vardict['lat']])
                    nlons = len(ncdata.dimensions[vardict['lon']])
                    
                    # Allocate the state array
                    print('\nAllocating the state vector array...')
                    state = np.zeros((nvars,ntimes,nlats,nlons))
                    print('state contains {} variables, {} times, {} lats, {} lons'.format(nvars, ntimes, nlats, nlons))
                    # For the metadata, need a list of locations
                    lats = ncdata.variables[vardict['lat']][:][:,None]
                    lons = ncdata.variables[vardict['lon']][:][None,:]
                    # Do a 2d mesh of lat and lon
                    lonarr, latarr = np.meshgrid(lons, lats)
                    
                    #And an array of ensemble members
#                    memarr = np.arange(1,nmems+1)
                
                
                # only increase the ensemble range if running through a new center
                if mem[0:80] != priorens[0:80]:
                    mnum2 += 1
                    # get the lower and upper ranges of the ensembles
                    timerange_lower = timerange_sum
                    timerange_sum = timerange_sum + each_time_length[mnum2]
                    timerange_upper = timerange_sum
                    print(timerange_lower)
                    print(timerange_upper)
                    print(timerange_sum)
                    print('Adding {} to state'.format(months[mnum2]))
                    
                    # Convert times back to integers
                    valid_times[timerange_lower:timerange_upper] = date2num(ftimes, tunit)
                priorens = mem
                
                # cycle through each variable (can be multiple per file)
                for v, var in enumerate(vrbls):
                    #if from sfc and pl, not each variable will be in each file
                    try:
                        field = ncdata.variables[var][:,:,:]
                        
                        # see if the data will fit into state (filters out bad data)
                        try:
                            state[v,timerange_lower:timerange_upper,:,:] = field
                            print('Adding {} to {}'.format(var, months[mnum2]))
                        #if the ensembles are a bad shape(missing times, etc)    
                        except ValueError:
                            state[v,timerange_lower:timerange_upper,:,:] = np.nan
                            if v==0:
                                print('{}: bad forecast array shape- not adding to combined analysis'.format(months[mnum2]))
                                del_mems.append(range(timerange_lower,timerange_upper))
                                bad_times_range += timerange_upper-timerange_lower
                    except:
                        pass
        # this runs if there is no data for the month
        except:
            # only increase the time range if running through a month
            if mem[0:80] != priorens[0:80]:
                mnum2 += 1
                print('Bad {} month- not adding to combined analysis :('.format(months[mnum2]))
            priorens = mem
            pass
                    

#------Haven't touched this from ensemble script, just commented out b/c don't need it for now.-------                
#    # Remove the nans (the incomplete members)
#    if len(del_mems) > 0:
#        state = np.delete(state,del_mems,axis=-1)
#        nmems -= bad_mems_range
#        #nmems -= len(del_mems)
#        memarr = np.arange(1,nmems+1)
#        print(del_mems)
#        print(len(del_mems))
#        print(state.shape)
#        print(nmems)
#        print(memarr)
       
    # If we are writing this out...
    if writenc:
        print('\nWriting to netcdf...')

        
        # Write ensemble forecast to netcdf
        with Dataset(outfile,'w') as dset:
            dset.createDimension('time',None)
            dset.createDimension('lat',nlats)
            dset.createDimension('lon',nlons)
            dset.createVariable('time','i4',('time',))
            dset.createVariable('lat',np.float64,('lat',))
            dset.createVariable('lon',np.float64,('lon'))
            dset.variables['time'].units = tunit
            dset.variables['lat'].units = 'degrees_north'
            dset.variables['lon'].units = 'degrees_east'
            dset.variables['time'][:] = np.array(valid_times)
            dset.variables['lat'][:] = lats
            dset.variables['lon'][:] = lons
            for v,var in enumerate(vrbls):
                var = vardict[var]
                print('Writing variable {}'.format(var))
                dset.createVariable(var, np.float32, ('time','lat','lon'))
                dset.variables[var].units = ut.get_units(var)
                dset.variables[var][:] = state[v,:,:,:]
        # Free up memory held by the state array
        del state
        
    # If we are NOT writing this out...
    else:
        # Reshape 5D state into a dictionary of 4D arrays
        allvars = {}
        for v,var in enumerate(vrbls):
            allvars[var] = (['validtime','y','x','mem'], state[v,:,:,:,:])
        # Package into an EnsembleState object knowing the state and metadata
        statecls = EnsembleState.from_vardict(allvars,
                                              {'validtime' : ftimes,
                                               'lat' : (['y','x'], latarr),
                                               'lon' : (['y','x'], lonarr),
                                               'mem' : memarr,
                                               })
        
        # Free up memory held by the state array
        del state
        
        return statecls
        
###############################################################################

#if __name__ == '__main__':
#    # Specify the initialization date and the length of the forecast
#    idate = datetime(2016,2,25,0)
#    dates = [idate + timedelta(days=d) for d in range(4)]
#    for idate in dates:
#        print('\n======= {:%Y%m%d} =========='.format(idate))
#        fdate = idate + timedelta(days=21)
#        
#        # Create the ensembles for the desired day-lags
#        for i in [4]: #[2,3,4]:
#            print('\nCREATING {}-DAY LAGGED ENSEMBLE'.format(i))
#            get_cfsv2_ensemble(i, idate, fdate, writenc=True)

#idate = datetime(2017,9,6,0)
#fdate = datetime(2017,9,16,0)
#
#get_cfsv2_ensemble(0, idate, fdate, vrbls=['Z500','T500','RH500','U500','V500', \
#                   'Z700','T700','RH700','U700','V700', 'Z850','T850','RH850', \
#                   'U850','V850','Z925','T925','RH925','U925','V925','Z1000', \
#                   'T1000','RH1000','U1000','V1000','T2M','RH2M','U10M','V10M', \
#                   'PWAT','MSLP','P6HR'])
    
#year = '2016'
#mon = ['10', '11']
#h = ['00', '12']
#for month in mon:
#    if month == '10':
#        d = range(1,32)
#    elif month == '11':
#        d = range(1,31)
#    for day in d:
#        day = format(day, '02d')
#        day = str(day)
#        for hour in h:
#            create_full_ensemble(year, month, day, hour)
#            #print('%s%s%s%s' %(year, month, day, hour))
#        #day.append(da)    

create_full_analysis(vrbls = ['t2m','tcw', 'gh'])
