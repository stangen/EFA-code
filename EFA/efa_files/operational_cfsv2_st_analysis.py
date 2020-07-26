 # -*- coding: utf-8 -*-
"""
Created on Wed May 11 14:19:56 2016

@author: njweber2
"""

from netCDF4 import Dataset, num2date, date2num
from datetime import datetime, timedelta
from subprocess import check_output
import numpy as np
import efa_files.cfs_utilities_st as ut
from efa_xray.state.ensemble import EnsembleState

def get_cfsv2_ensemble(ndays, start, end, vrbls=['Z500'],
                       only00z=True, writenc=True):
    """
    Populates and returns a CFSv2 global ensemble forecast using the
    archived operational forecasts on vader.atmos.washington.edu.
    
    Requires:    
    ndays -> The number of init days we want to use to populate the ensemble.
             16 forecasts are initialized per day. That is, an ensemble with
             16*ndays members will be built.
    start -> A datetime object of the ensemble forecast initialization time.
    end   -> A datetime object of the ensemble forecast end time.
    vrbls -> A list of variables we want to retrieve.
    only00z -> A boolean object. If True, only the 00z validation times are
               saved. If False, all 6-hourly times are saved.
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
               'T2M' : 'TMP_2maboveground',
               'PWAT' : 'PWAT_entireatmosphere_consideredasasinglelayer_',
               'MSLP' : 'PRMSL_meansealevel',
               'P6HR' : 'APCP_surface',
               'time' : 'time',
               'lat' : 'latitude',
               'lon' : 'longitude',
               }
    
    # This is the input directory for the CFSv2 forecast netcdfs
    indir  = '/home/disk/hot/stangen/Documents/GEFS/analysis/2017090600_2017091600/netcdf/precip'
    # This is the output directory for the ncfile if writenc==True
    outdir = '/home/disk/hot/stangen/Documents/GEFS/analysis/2017090600_2017091600/ensembles/precip'
  
    # Get a list of filenames (each file is a different member)
    mem = []
# if you want to include time-lag to increase ensemble, uncomment the for loop
# and the timedelta part and indent after the for loop.
#    for i in range(4*ndays):
    idate = start #- timedelta(hours=i*6)
    # List the four members at this time and append to main list
    datestr = idate.strftime('%Y%m%d%H')
    command = 'ls -1a {}/*.nc'.format(indir,datestr)
    print(command)
    mem.extend(list(reversed(check_output([command],shell=True).split())))
        
    # Loop through the individual member files to load the forecasts
    print('Loading {} ensemble...'.format(start.strftime('%Y-%m-%d %H:00')))
    del_mems = []
    mem = mem[0].decode('utf-8')
    print(mem)
    # Read the netcdf
    with Dataset(mem,'r') as ncdata:
        # Find the indices corresponding to the start and end times
        tunit = ncdata.variables[vardict['time']].units
        ftimes = num2date(ncdata.variables[vardict['time']][:],tunit)
        tbeg = ut.nearest_ind(ftimes,start)
        tend = ut.nearest_ind(ftimes,end)
        ftimes = ftimes[tbeg:tend+1] #for metadata
        
        # If this is the first member, calculate how large the state array
        # needs to be and allocate. Also set up the metadata.
#                nmems = len(memfiles)
#                nmems4name = len(memfiles)
        ntimes = len(ftimes)
        nvars = len(vrbls)
        nlats = len(ncdata.dimensions[vardict['lat']])
        nlons = len(ncdata.dimensions[vardict['lon']])
        
        # Allocate the state array
        print('Allocating the state vector array...')
        state = np.zeros((nvars,ntimes,nlats,nlons))
        
        # For the metadata, need a list of locations
        lats = ncdata.variables[vardict['lat']][:][:,None]
        lons = ncdata.variables[vardict['lon']][:][None,:]
        # Do a 2d mesh of lat and lon
        lonarr, latarr = np.meshgrid(lons, lats)
        
        #And an array of ensemble members
        
        #field = ncdata.variables[vardict['Z500']][tbeg:tend,:,:]
        # Now to populate the state array
        for v, var in enumerate(vrbls):
            field = ncdata.variables[vardict[var]][:,:,:]#[tbeg:tend,:,:]
            #print(field)
            print('Adding variable {}'.format(var))
             # Populate its component of the state array
            try:
                state[v,:,:,:] = field
            except ValueError:
                state[v,:,:,:] = np.nan
                if v==0:
                    print('  member {}: bad forecast array shape'.format(mnum))
                    del_mems.append(mnum)
        # END of ncdata load
    
    # Grab only the 00z times, if appropriate
    if only00z:
        ftimes = ftimes[::4]
        state = state[:,::4,:,:,:]
    
    # Remove the nans (the incomplete members)
#    if len(del_mems) > 0:
#        state = np.delete(state,del_mems,axis=-1)
#        nmems -= len(del_mems)
#        memarr = np.arange(1,nmems+1)
            
    # If we are writing this out...
    if writenc:
        print('Writing to netcdf...')
        # Convert times back to integers
        valid_times = date2num(ftimes,tunit)
        outfile = '{}/{:%Y%m%d%H}_{}days.nc'.format(outdir,start,
                                                          (end-start).days)
        
        # Write ensemble forecast to netcdf
        with Dataset(outfile,'w') as dset:
            dset.createDimension('time',None)
            dset.createDimension('lat',nlats)
            dset.createDimension('lon',nlons)
#            dset.createDimension('ens',nmems)
            dset.createVariable('time','i4',('time',))
            dset.createVariable('lat','f8',('lat',))
            dset.createVariable('lon','f8',('lon'))
#            dset.createVariable('ens','i4',('ens',))
            dset.variables['time'].units = tunit
            dset.variables['lat'].units = 'degrees_north'
            dset.variables['lon'].units = 'degrees_east'
#            dset.variables['ens'].units = 'member_number'
            dset.variables['time'][:] = np.array(valid_times)
            dset.variables['lat'][:] = lats
            dset.variables['lon'][:] = lons
#            dset.variables['ens'][:] = memarr
            for v,var in enumerate(vrbls):
                print('Writing variable {}'.format(var))
                dset.createVariable(var, 'f8', ('time','lat','lon'))
                dset.variables[var].units = ut.get_units(var)
                dset.variables[var][:] = state[v,:,:,:]
        # Free up memory held by the state array
        del state
        
    # If we are NOT writing this out...
    else:
        # Reshape 5D state into a dictionary of 4D arrays
        allvars = {}
        for v,var in enumerate(vrbls):
            allvars[var] = (['validtime','y','x',], state[v,:,:,:])
        # Package into an EnsembleState object knowing the state and metadata
        statecls = EnsembleState.from_vardict(allvars,
                                              {'validtime' : ftimes,
                                               'lat' : (['y','x'], latarr),
                                               'lon' : (['y','x'], lonarr),
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

idate = datetime(2017,9,6,6)
fdate = datetime(2017,9,16,6)

get_cfsv2_ensemble(0, idate, fdate, vrbls=['Z500','T2M','PWAT','MSLP','P6HR'], only00z=False)