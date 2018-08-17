#!/usr/bin/env python
import os
import sys
import netCDF4
from datetime import datetime
from time import strftime
from datetime import timedelta
import madis_utilities as mt


#-------------Modify these-----------------------------------------------------
#if true, obtain the observation type from MADIS server
metar = True
mesonet = False
maritime = True
raob = False

#start and end time- will get observations at 6-hour intervals between these times.
start_time = '20180614_0000'
end_time = '20180614_0000'
#------------------------------------------------------------------------------

#Set working directory
base_dir="/home/disk/hot/stangen/Documents/surface_obs"

#Get current date/time
dt_now = datetime.utcnow()
dt_hour = dt_now.replace(minute = 0, second=0, microsecond=0) # Truncate current datetime to hour precision
dt_hour = dt_hour - timedelta(0,3600) # Get previous hour (MADIS data is usually near-realtime, 1h behind)


#make a list of dates
#see make_datelist documentation for more information on this
dt_list = mt.make_datelist(start_time,end_time,torf=False,hour_before_after=True)

#Function to retrieve MADIS data, by observation type, from NOAA ftp server
def retrieve_madis_type(ob_type,ob_selected):
	if (ob_type == "mesonet"):
		#Mesonet data is stored in the LDAD folder on MADIS ftp server
		loc = "LDAD"; nc = "netCDF"
	else:
		#METAR, MARITIME, etc.. are in point folder on MADIS ftp server
		loc = "point"; nc = "netcdf"

	#Retrieve MADIS observations and write data to text (CSV) file
	if ((os.path.isfile(save_dir+"/"+ob_type+"."+fname_gz) == 0) and (ob_selected)):
		os.system(retrieval_string+" "+url_short+"/"+loc+"/"+ob_type+"/"+nc+"/"+fname_gz+" -o "+save_dir+"/"+ob_type+"."+fname_gz)
		os.system("gunzip -f "+save_dir+"/"+ob_type+"."+fname_gz)
		#Read MADIS data and extract QC'ed data to file
		os.system(save_dir+"/read_"+ob_type+".py "+ob_type+"."+fname)
    

#Loop through date list to retrieve observations hour by hour	
for d in range(0,len(dt_list)):

    #Parse date
    d2 = dt_list[d]
    yyyy = d2[0:4]
    mons = d2[4:6]
    dys = d2[6:8]
    hrs = d2[8:10] 

    #Define the name of the MADIS (gzipped) file
    fname_gz = yyyy+mons+dys+"_"+hrs+"00.gz"
    fname = yyyy+mons+dys+"_"+hrs+"00"

    #Create directories if they don't yet exit
    if (os.path.isdir(base_dir+"/MADIS/"+yyyy+mons+"/raw/")):
        pass
    else:
        os.makedirs(base_dir+"/MADIS/"+yyyy+mons+"/raw/")
    
    save_dir= base_dir+"/MADIS/"+yyyy+mons+"/raw/"

    #Retrieve data using curl account for connection timeout (allow retry) 
    retrieval_string = "curl --connect-timeout 20 --http1.0 --max-time 60 --retry 10 --retry-delay 10 --retry-max-time 120 "

    #Determine if MADIS observation data is archived or near-realtime
    #MADIS data more than 5 days old is stored on a different server
    if ((dt_hour - datetime.strptime(d2,'%Y%m%d%H')).total_seconds() >= 86400*5):
        url_short = "ftp://madis-data.cprk.ncep.noaa.gov/archive/"+str(yyyy)+"/"+str(mons)+"/"+str(dys)
    else:
        url_short = "ftp://madis-data.cprk.ncep.noaa.gov"

    #Retrieve MADIS data from surface observing networks and write QC'ed data to CSV files
    print("Retrieving MADIS data for "+yyyy+mons+dys+hrs)
    retrieve_madis_type("metar",metar)
    retrieve_madis_type("mesonet",mesonet)
    retrieve_madis_type("maritime",maritime)
    retrieve_madis_type("raob",raob)
	
