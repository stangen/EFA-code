# -*- coding: utf-8 -*-
"""
Created on Tue May 17 15:39:31 2016

@author: njweber2/spencer tangen
"""
from subprocess import check_output
import os
from datetime import datetime, timedelta
from multiprocessing import Pool
import nicks_files.cfs_utilities as ut

# 0-20 for 21 mems of gefs
def create_fcst_netcdfs(date, mems=range(0,21)):
    """
    Creates/saves a netcdf for each of the four CFSv2 forecasts (members)
    initialized on the given date. The forecasts are compiled from the
    archived grib files on Rick Steed's account on vader. Only the desired 
    parameters are extracted from the gribs and saved to the netcdfs.
    
    Requires:
    dates -> a datetime object specifying the initialization time of the 
             forecastst
    mems  -> the desired members at this initialization time
    
    Returns:
    Nothing! Will simply save the forecasts to netcdf files
    
    IMPORTANT: This will not run correctly in Spyder- run it instead in
    terminal.
    """
    
        # Designate the parameters and levels desired
    # *formatting is important here for the wgrib2 command!
    pars = "':(HGT|TMP|RH|UGRD|VGRD|PRES|PWAT|APCP):'"
    levs = "':(500 mb|700 mb|850 mb|925 mb|1000 mb|2 m above ground|10 m above" + \
           " ground|surface):'"
    allpars = "':(HGT:500 mb|TMP:500 mb|RH:500 mb|UGRD:500 mb|VGRD:500 mb|" + \
              "HGT:700 mb|TMP:700 mb|RH:700 mb|UGRD:700 mb|VGRD:700 mb|" + \
              "HGT:850 mb|TMP:850 mb|RH:850 mb|UGRD:850 mb|VGRD:850 mb|" + \
              "HGT:925 mb|TMP:925 mb|RH:925 mb|UGRD:925 mb|VGRD:925 mb|" + \
              "HGT:1000 mb|TMP:1000 mb|RH:1000 mb|UGRD:1000 mb|VGRD:1000 mb|" + \
              "TMP:2 m above ground|RH:2 m above ground|PRES:mean sea level|" + \
              "UGRD:10 m above ground|VGRD:10 m above ground|APCP|PWAT):'"
    # Only want the PWAT of the entire atmosphere
    exclude = "':(PWAT:30-0 mb above ground):'"
    
    # I/O directories
    datestr = date.strftime('%Y%m%d%H') # yyyymmddhh
    print(datestr)
    # ST changed indir/outdir to my directory
    indir = '/home/disk/hot/stangen/Documents/GEFS/data/{}'.format(datestr)
    outdir = '/home/disk/hot/stangen/Documents/GEFS/netcdf/{}/'.format(datestr)
    
    # Loop through the desired members
    for mem in mems:
        print(' Member {}'.format(mem))
        outfile = '{}pgbf_{}_{:02d}.nc'.format(outdir,datestr,mem)
    #     Remove the file if it exists (otherwise -append won't work)
        if os.path.isfile(outfile):
            os.system('rm -f {}'.format(outfile))
        
        # List all of the grb2 files for this member
        command = 'ls -1a {}/ge*{:02d}.t{:02d}z.pgrb2f*'.format(indir,mem,date.hour)
        files = list(check_output([command],shell=True).split())
        # Get rid of the utf-8 filetype
        files_no_utf8 = []
        for grb in files:
            grb = grb.decode('utf-8')
            files_no_utf8.append(grb)
        # Sort the files correctly in ascending time to append correctly
        files_no_utf8.sort(key=ut.natural_keys)
        
        # We only want the first 45 days of forecast
        #files = files[:45*4]
        # ST 1 day of forecast (0,6,12,18,24)
        #files = files[:5*21]
        
        
        # Convert all files to netcdfs, saving only the desired variables
        # The -append keyword allows us to combine all times into 1 netcdf
        for grb in files_no_utf8:
            command = 'wgrib2 {} -match {} -not {} -append -netcdf {} '
            #command = 'wgrib2 {} -match {} -match {} -append -netcdf {} '
            command = command.format(grb,allpars,exclude,outfile)
            #command = command.format(grb.decode('utf-8'),pars,levs,outfile)
            print(command)
            os.system(command)
            #check_output([command], shell=True)
        
###############################################################################
if True:            
#if __name__ == '__main__':
    #======== EDIT THESE ======================================================
    # Designate date range for forecast processing
    idate = datetime(2016,1,25,0)
    fdate = datetime(2016,2,1,0)
    # Number of processors for multiprocessing
    nproc = 16
    mproc = False
    #==========================================================================
    
    # Create list of dates with 6-hourly intervals
    #dates = list(ut.perdelta(idate,fdate,timedelta(hours=6)))
    dates = datetime(2017,9,6,0)
    
    if mproc:
        # Use multiprocessing to create all the netcdfs
        print('USING {} PROCESSORS'.format(nproc))
        pool =  Pool(processes=nproc)
        results = pool.map(create_fcst_netcdfs,dates)
        pool.close()
        pool.join()
    else:
        # Just loop through the dates with one processor
#        for date in dates:
#            create_fcst_netcdfs(date)
        create_fcst_netcdfs(dates)
            