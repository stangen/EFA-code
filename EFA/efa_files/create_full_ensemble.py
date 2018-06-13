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

def create_full_ensemble(year, month, day, hour, vrbls=['Z500'], writenc=True):
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
# for each date I want to do this to:
    # for each of the ensembles (ecmwf, jma, ncep, eccc)
        # if the ensemble file for a date exists
            # put it into the big ensemble
    
    
    
    
    
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
    outfile = '/home/disk/hot/stangen/Documents/ensembles/all/%s%s/%s-%s-%s_%s.nc' % (vardict[month], year, year, month, day, hour)
    
    # directories for the ensembles, initialize lists used later
    ensembles = ['ecmwf', 'jma', 'ncep', 'eccc']
    memfiles = []
    each_mem_length = list()
    state = []
    
    # load the ensemble files into list called memfiles
    for ens in ensembles:    
        indir = '/home/disk/hot/stangen/Documents/ensembles/%s/%s%s/%s-%s-%s_%s_%s*' % (ens, vardict[month], year, year, month, day, hour, ens)
        # Get a list of filenames (each file is a different member)
        command = 'ls -1a {}'.format(indir)
        memfiles.extend(list(reversed(check_output([command],shell=True).split())))
    
    # more initializations....
    nmems = 0
    priorens = '                                                                                  '
    m2 = -1
    
    # Count the number of ensemble members for allocation of state later
    # try blocks are to filter out dates when TIGGE didn't have any data
    print('\nCreating superensemble for {}/{}/{} {}Z\n'.format(month,day,year,hour))
    print('Counting ensemble members for allocation')
    for m, mem in enumerate(memfiles):
        mem = mem.decode('utf-8')
        # to not double-count ens members if from both sfc and pl
        if mem[0:50] != priorens[0:50]:
            m2 += 1
            try:
                with Dataset(mem, 'r') as ncdata:                   
                    #keep track of number of members of each ensemble
                    each_mem_length.append(ncdata.variables[vrbls[0]][:,:,:,:].shape[1])
                    #keep track of total length of all ensembles
                    nmems = nmems + ncdata.variables[vrbls[0]][:,:,:,:].shape[1]
                    print('ensembles from {}: {}'.format(ensembles[m2],each_mem_length[m2]))
                    print('total members: {}'.format(nmems))
            # if there is no data from one of the ensembles, add a 0 to list of number of ensembles
            except:
                print('No data from {} ensembles'.format(ensembles[m2]))
                each_mem_length.append(0)
                #pass
        # keep track of file name to compare with next file
        priorens = mem
        
    # even more initialization...
    memrange_sum = 0
    bad_mems_range = 0    
    priorens = '                                                                                  '
    del_mems = []
    mnum2 = -1
    
    # loop to add data from each file
    for mnum, mem in enumerate(memfiles):
        mem = mem.decode('utf-8')
        
        try:
            # Read the netcdf
            with Dataset(mem,'r') as ncdata:
                # Find the indices corresponding to the start and end times
                tunit = ncdata.variables[vardict['time']].units
                ftimes = num2date(ncdata.variables[vardict['time']][:],tunit)
    
                # If this is the first member, calculate how large the state array
                # needs to be and allocate. Also set up the metadata.               
                if mnum == 0:
                    ntimes = len(ftimes)
                    nvars = len(vrbls)
                    nlats = len(ncdata.dimensions[vardict['lat']])
                    nlons = len(ncdata.dimensions[vardict['lon']])
                    
                    # Allocate the state array
                    print('\nAllocating the state vector array...')
                    state = np.zeros((nvars,ntimes,nlats,nlons,nmems))
                    print('state contains {} variables, {} times, {} lats, {} lons, {} ensembles'.format(nvars, ntimes, nlats, nlons, nmems))
                    # For the metadata, need a list of locations
                    lats = ncdata.variables[vardict['lat']][:][:,None]
                    lons = ncdata.variables[vardict['lon']][:][None,:]
                    # Do a 2d mesh of lat and lon
                    lonarr, latarr = np.meshgrid(lons, lats)
                    
                    #And an array of ensemble members
                    memarr = np.arange(1,nmems+1)
                
                
                # only increase the ensemble range if running through a new center
                if mem[0:50] != priorens[0:50]:
                    mnum2 += 1
                    # get the lower and upper ranges of the ensembles
                    memrange_lower = memrange_sum
                    memrange_sum = memrange_sum + each_mem_length[mnum2]
                    memrange_upper = memrange_sum
                    print('Adding {} to state'.format(ensembles[mnum2]))
                priorens = mem
                
                # cycle through each variable (can be multiple per file)
                for v, var in enumerate(vrbls):
                    #if from sfc and pl, not each variable will be in each file
                    try:
                        field = ncdata.variables[var][:,:,:,:]
                        
                        # see if the data will fit into state (filters out bad data)
                        try:
                            # make the ensembles at the end of state
                            field = np.swapaxes(field, 1, 3)
                            field = np.swapaxes(field, 1, 2)              
                            state[v,:,:,:,memrange_lower:memrange_upper] = field
                            print('Adding {} to {}'.format(var, ensembles[mnum2]))
                        #if the ensembles are a bad shape(missing times, etc)    
                        except ValueError:
                            state[v,:,:,:,memrange_lower:memrange_upper] = np.nan
                            if v==0:
                                print('{}: bad forecast array shape- not adding to superensemble'.format(ensembles[mnum2]))
                                del_mems.append(range(memrange_lower,memrange_upper))
                                bad_mems_range += memrange_upper-memrange_lower
                    except:
                        pass
        # this runs if there is no data for the ensemble, or it can't read the data
        except:
            # only increase the ensemble range if running through a new center
            if mem[0:50] != priorens[0:50]:
                mnum2 += 1
                print('Bad {} ensembles- not adding to superensemble'.format(ensembles[mnum2]))
            priorens = mem
            pass
                    
                
    # Remove the nans (the incomplete members)
    if len(del_mems) > 0:
        state = np.delete(state,del_mems,axis=-1)
        nmems -= bad_mems_range
        #nmems -= len(del_mems)
        memarr = np.arange(1,nmems+1)

            
    # If we are writing this out...
    if writenc:
        print('\nWriting to netcdf...')
        # Convert times back to integers
        valid_times = date2num(ftimes,tunit)
        #outfile = '{}/{:%Y%m%d%H}_{}mem_{}days.nc'.format(outdir,start,nmems4name,
        #                                                  (end-start).days)
        
        # Write ensemble forecast to netcdf
        with Dataset(outfile,'w') as dset:
            dset.createDimension('time',None)
            dset.createDimension('lat',nlats)
            dset.createDimension('lon',nlons)
            dset.createDimension('ens',nmems)
            dset.createVariable('time','i4',('time',))
            dset.createVariable('lat',np.float64,('lat',))
            dset.createVariable('lon',np.float64,('lon'))
            dset.createVariable('ens','i4',('ens',))
            dset.variables['time'].units = tunit
            dset.variables['lat'].units = 'degrees_north'
            dset.variables['lon'].units = 'degrees_east'
            dset.variables['ens'].units = 'member_number'
            dset.variables['time'][:] = np.array(valid_times)
            dset.variables['lat'][:] = lats
            dset.variables['lon'][:] = lons
            dset.variables['ens'][:] = memarr
            for v,var in enumerate(vrbls):
                var = vardict[var]
                print('Writing variable {}'.format(var))
                dset.createVariable(var, np.float32, ('time','lat','lon','ens',))
                dset.variables[var].units = ut.get_units(var)
                dset.variables[var][:] = state[v,:,:,:,:]
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

create_full_ensemble('2016', '12', '06', '12', vrbls = ['t2m','tcw', 'gh'])
