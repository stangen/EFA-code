#!/usr/bin/env python
import sys
import numpy as np
from datetime import datetime
from datetime import timedelta
import sys
import os
import netCDF4 
import madis_utilities as mt

base_dir = "/home/disk/hot/stangen/Documents/EFA/surface_obs/MADIS/"

#Set extent of MADIS data to extract from file
minLat = 0.01 
maxLat = 90
minLng = -180
maxLng = 180

# Change these to read different data and times
start_date = '20130401_0000'
end_date = '20130512_1800'
ob_type = ['metar','maritime'] #Can add 'raob'
ob_vars = ['altimeter','temperature'] #or 'altimeter'

#Call a function which creates a list of dates in 6 hour intervals
date_list = mt.make_datelist(start_date, end_date, torf = True, hour_before_after=False)


def read_madis(ob_type='metar', ob_var='altimeter', dtstr='20130401_0000'):
    """
    Reads and filters observations based on metar, maritime, or raob, and temperature 
    or altimeter as desired. Also takes in the date/hour of the observations
    """
    
    #Open the file
    fname = yyyy+mons+'/raw/'+ob_type+'.'+dtstr
    f1 = netCDF4.Dataset(base_dir+"/"+fname,"r")
    
    #Read in netcdf file 
    variableNames = f1.variables.keys() # print to show the names of all variables in the netCDF file
    if (ob_type =='metar' or ob_type =='maritime'):        
        lats = f1.variables['latitude'][:] # latitude of MADIS observation
        lngs = f1.variables['longitude'][:] # longitude of MADIS observation
        elevs = f1.variables['elevation'][:] # Elevation of MADIS observation 
        epoch = f1.variables['timeObs'][:] # Time of MADIS Observation
        stn_name = f1.variables['stationName'][:] # Name of MADIS observation station
    elif ob_type =='raob':
        lats = f1.variables['staLat'][:] # latitude of MADIS observation
        lngs = f1.variables['staLon'][:] # longitude of MADIS observation
        elevs = f1.variables['staElev'][:] # Elevation of MADIS observation
        epoch = f1.variables['synTime'][:] # Time of MADIS Observation
        stn_name = f1.variables['staName'][:] # Name of MADIS observation station
    if ob_type =='metar':
        station_type = np.zeros(len(lats))    #all metars are stationary        
        if ob_var == 'altimeter':
            var = f1.variables['altimeter'][:] # MADIS Altimeter Setting observation
            qc_var0 = f1.variables['altimeterDD'][:] # QC string for MADIS Altimeter Setting Observation
            qcr_var = f1.variables['altimeterQCR'][:] # QC boolean for MADIS Altimeter Setting Observation
        elif ob_var == 'temperature':
            var = f1.variables['temperature'][:] # MADIS Temperature observation
            qc_var0 = f1.variables['temperatureDD'][:] # QC string for MADIS Temperature Observation
            qcr_var = f1.variables['temperatureQCR'][:] # QC boolean for MADIS Temperature Observation
    elif ob_type =='maritime':
        station_type = f1.variables['dataPlatformType'][:] #see whether it's stationary or moving
        if ob_var =='altimeter':
            var = f1.variables['seaLevelPress'][:] # MADIS sea level pressure observation
            qc_var0 = f1.variables['seaLevelPressDD'][:] # QC string for MADIS sea level pressure Observation
            qcr_var = f1.variables['seaLevelPressQCR'][:] # QC boolean for MADIS sea level pressure Observation
        elif ob_var == 'temperature':
            var = f1.variables['temperature'][:] # MADIS Temperature observation
            qc_var0 = f1.variables['temperatureDD'][:] # QC string for MADIS Temperature Observation
            qcr_var = f1.variables['temperatureQCR'][:] # QC boolean for MADIS Temperature Observation
    elif ob_type =='raob':
        if ob_var =='altimeter':
            pres = f1.variables['prMan'][:,0] # MADIS Surface pressure observation (retrieve only first level (0) - i.e. surface altimter for RAOB)
            #Convert surface pressure to altimeter setting in mb
            presinHg = pres*.029528744
            var = (presinHg/((288-0.0065*elevs)/288)**5.2561)/.029528744
            qc_var0 = f1.variables['prManDD'][:,0] # QC string for MADIS Altimeter Setting Observation
            qcr_var = f1.variables['prManQCR'][:,0] # QC boolean for MADIS Altimeter Setting Observation
        elif ob_var =='temperature':
            var = f1.variables['tpMan'][:,0]
            qc_var0 = f1.variables['tpManDD'][:,0]
            qcr_var = f1.variables['prManQCR'][:,0]

    f1.close() # Close netCDF file
    
    #join station name (convert from list of bytes to string)
    stns = []
    for s in stn_name:
    	#decode list of bytes and join them to get the station identifier string
    	stns.append(b''.join(s).decode('utf-8'))
    
    #decode other character arrays (decoding is only necessary in Python 3)
    qc_var = [b.decode() for b in qc_var0]
    
    # Use QC Boolean to determine whether obs should be kept. If False interpret QC string
    QCR = False
    MADIS_QC = True
    
    # If QC Boolean == 0 than the given observation passed the QC checks applied by MADIS.
    
    # Its important to note that MADIS doesn't always apply rigorous QC checks. Some MADIS observations
    # may have only gone through a basic validity check (and not a more rigorous temporal or spatial consistency check)
    
    # If observation quality is paramount, it may be best to interpret the MADIS QC string and select only observations
    # undergoing at least 2 or 3 QC checks. MADIS QC check information for a given observation type (MESONET, METAR, etc.)
    # can be retrieved by running the ncdump -h command on the MADIS netCDF file.
    # The text below was copied from the contents of a ncdump -h from a METAR netCDF file:
    
    #:DD_long_name = "QC data descriptor model:  QC summary values" ;
    #:DD_reference = "AWIPS Technique Specification Package (TSP) 88-21-R2" ;
    #:DD_values = "Z,C,S,V,X,Q,K,k,G, or B" ;
    #:DD_value_Z = "No QC applied" ;
    #:DD_value_C = "Passed QC stage 1" ;
    #:DD_value_S = "Passed QC stages 1 and 2" ;
    #:DD_value_V = "Passed QC stages 1, 2 and 3" ;
    #:DD_value_X = "Failed QC stage 1" ;
    #:DD_value_Q = "Passed QC stage 1, but failed stages 2 or 3 " ;
    #:DD_value_K = "Passed QC stages 1, 2, 3, and 4" ;
    #:DD_value_k = "Passed QC stage 1,2, and 3, failed stage 4 " ;
    #:DD_value_G = "Included in accept list" ;
    #:DD_value_B = "Included in reject list" ;
    
    #:QCStage_long_name = "automated QC checks contained in each stage" ;
    #:QCStage_values = "1, 2, 3, or 4" ;
    #:QCStage_value_1 = "Validity and Position Consistency Check" ;
    #:QCStage_value_2 = "Internal, Temporal, and Model Consistency Checks" ;
    #:QCStage_value_3 = "Spatial and Statistical Spatial Consistency Checks" ;
    #:QCStage_value_4 = "Kalman Filter" ;
    
    #For more info about the stages of MADIS QC see: 
    #https://madis.ncep.noaa.gov/madis_sfc_qc_notes.shtml
    
    #In this script, when QCR is set to False, only METAR observations passing at least 2 stages of QC checks are retained.
    
    #Convert list to numpy array
    lats = np.array(lats)
    lngs = np.array(lngs)
    
    #Initialize string array (to store variables and write to ascii file)
#    sstr = []
    
    #Checking if stations are in China (psedoglobal coverage check)
#    z = 0
#    l = 0
#    for i in stns:
#        if i[0] == 'Z':
#            z = z+1
#    
#    for j in lngs:
#        if j<120 and j >100:
#            l = l+1
    
    #Loop through all MADIS observations. 
    for t in range(0,len(lats)):
        #if a maritime ob is missing elevation data, replace with 0.0
        if elevs[t] !='--':  
            pass
        else:
            if ob_type == 'maritime':
                elevs[t] = 0.0
        if MADIS_QC == True:
        	#If the MADIS ob is in the lat/lng bounding box and has passed MADIS QC checks add the altimeter observation to the string array
            #If "S" is not included, the QC will only accept passing at least QC 3 - for global coverage, must contain S
            if ((not QCR) and (float(lats[t]) >= minLat) and (float(lats[t]) <= maxLat) and (float(lngs[t]) >= minLng) and (float(lngs[t]) <= maxLng) and station_type[t] == 0 and ((qc_var[t] == "K") or (qc_var[t] == "V") or (qc_var[t] == "S"))):
        		#Save station identifier, latitude, longitude, elevation, epoch time, and observation
            #if (alts[t] != '--'):
                if ob_var == 'altimeter':
                    if var[t] > 94000 and var[t] < 106000:
                        sstr.append(str(stns[t])+","+str(lats[t])+","+str(lngs[t])+","+str(elevs[t])+","+str(epoch[t])+","+str(var[t]/100)+","+ob_type.upper()+"\n")
                else:
                    sstr.append(str(stns[t])+","+str(lats[t])+","+str(lngs[t])+","+str(elevs[t])+","+str(epoch[t])+","+str(var[t])+","+ob_type.upper()+"\n")
                #sstr.append(str(stns[t])+","+str(lats[t])+","+str(lngs[t])+","+str(elevs[t])+","+str(epoch[t])+","+str(var[t])+","+ob_type.upper()+"\n")
            #If QCR is set to "True", use the MADIS QC check
            elif ((QCR) and (var[t] != '--') and (float(lats[t]) >= minLat) and (float(lats[t]) <= maxLat) and (float(lngs[t]) >= minLng) and (float(lngs[t]) <= maxLng) and station_type[t] == 0 and (qcr_var[t] == 0)):
        		#If the MADIS ob exists, is in the lat/lng bounding box and has passed all MADIS QC checks applied
        		#add the altimeter observation to the string array
                if ob_var == 'altimeter':
                    if var[t] > 94000 and var[t] < 106000:
                        sstr.append(str(stns[t])+","+str(lats[t])+","+str(lngs[t])+","+str(elevs[t])+","+str(epoch[t])+","+str(var[t]/100)+","+ob_type.upper()+"\n")
                else:
                    sstr.append(str(stns[t])+","+str(lats[t])+","+str(lngs[t])+","+str(elevs[t])+","+str(epoch[t])+","+str(var[t])+","+ob_type.upper()+"\n")
        #Use Luke's method of QCing data: altimeter setting 880-1100 hPa, temp -40 to 40 C (233.15-313.15)
        elif MADIS_QC == False:
            #for altimeter setting---- CONVERTED TO hPa:
            if ob_var == 'altimeter':
                if ((var[t] != '--') and (float(lats[t]) >= minLat) and (float(lats[t]) <= maxLat) and (float(lngs[t]) >= minLng) and (float(lngs[t]) <= maxLng) and var[t] > 88000 and var[t] < 110000 and station_type[t] == 0):
                    sstr.append(str(stns[t])+","+str(lats[t])+","+str(lngs[t])+","+str(elevs[t])+","+str(epoch[t])+","+str(var[t]/100)+","+ob_type.upper()+"\n")
            #for temperature:
            elif ob_var == 'temperature':
                if ((var[t] != '--') and (float(lats[t]) >= minLat) and (float(lats[t]) <= maxLat) and (float(lngs[t]) >= minLng) and (float(lngs[t]) <= maxLng) and var[t] > 233.15 and var[t] < 313.15 and station_type[t] == 0):
                    sstr.append(str(stns[t])+","+str(lats[t])+","+str(lngs[t])+","+str(elevs[t])+","+str(epoch[t])+","+str(var[t])+","+ob_type.upper()+"\n")            
    
    #Sort the observations by station identifier 
    #list.sort(sstr)
    
    print("Non-QC'd number of "+ob_var+" observations for "+dtstr+" : "+str(len(lats)))
    print("Number of "+ob_type+" "+ob_var+" observations retrieved: "+str(len(sstr)))
    
    return(sstr)

#--------------End of Function--------------------------------------------------

# Loop through all the 6-hourly dates in the list      
for dates in date_list:
    for ob_var in ob_vars:

        #Shorten input ob variable to match naming convention to read file
        if ob_var == 'altimeter':
            var_short = 'alts'
        elif ob_var =='temperature':
            var_short = 'temp'
        
        #Convert 6-hourly date string to datetime object, useful for later but 
        #doesn't need to be redefined in the loop where it is accessed, so defined here.
        desired_time = datetime.strptime(dates, '%Y%m%d_%H%M')
    
        #Get the previous and following hour to also load. 
        #Convert dates string to datetime to add/subtract one hour
        dates_hour_before = datetime.strptime(dates,'%Y%m%d_%H%M') - timedelta(seconds = 3600)
        dates_hour_after = datetime.strptime(dates,'%Y%m%d_%H%M') + timedelta(seconds = 3600)
        #Convert back to string to have same style as dates
        dates_hour_before= dates_hour_before.strftime('%Y%m%d_%H00')
        dates_hour_after = dates_hour_after.strftime('%Y%m%d_%H00')
        #Combine the hour before and dates strings for use in looping to load data later
        hourlist = [dates_hour_before,dates,dates_hour_after]
        
        #Get the year and month of the 6-hourly dates
        yyyy = dates[0:4]
        mons = dates[4:6]
        
        
        #Loop through each ob_type (metar, maritime, raob)
        for ob in ob_type:
            #Create metar/maritime/raob (separate for each ob) directories if they don't yet exist
            if (os.path.isdir(base_dir+"/"+yyyy+mons+"/"+ob+"_"+var_short+"/")):
                pass
            else:
                os.makedirs(base_dir+"/"+yyyy+mons+"/"+ob+"_"+var_short+"/") 
            
            #Create combined directory (with all ob_types in it) if it doesn't yet exist
            if (os.path.isdir(base_dir+"/"+yyyy+mons+"/"+"combined_"+var_short+"/")):
                pass
            else:
                os.makedirs(base_dir+"/"+yyyy+mons+"/"+"combined_"+var_short+"/")
            
            #Create an empty list to append observations into. Must be done outside of hourlist
            #for loop to append hour before, hour after, and dates observations to the list.
            sstr = []
            #Create an empty list to append observations into, without duplicate stations. 
            #Must be done outside of hourlist for loop to append the files
            sstr_one_station = [] 
            
            #Loop through each hour (hour before, hour after, and dates)
            for hours in hourlist:
                #try is for if a file didn't gunzip properly
                try:
                
                #change yyyy and mons so that the proper file is loaded in the read_madis function.
                #this is because read_madis uses yyyy and mons and it needs to correlate with the correct hour.
                    yyyy = hours[0:4]
                    mons = hours[4:6]
                
#-------------------Call the function to read a MADIS netCDF file and filter the data.
                    #This function appends QC'd observations to sstr. 
                    sstr = read_madis(ob,ob_var,hours)
                except:
                    pass
                                   
            #sort the observation list by station name.
            list.sort(sstr)
      
            # Remove duplicate stations so that one station per observation is kept. 
            # This is done because some stations have up to ~10 duplicate observations,
            # while some have multiple observations per hour. This also only keeps
            # the observations closest to the dates time, if there are multiple obs
            # per hour, and only keeps the observations within an hour.
            
            #Loop through each observation (this is still within ob loop, so this
            #happens one ob_type at a time. Includes all qc'd obs from the hour before, 
            #hour of, and hour after.)
            for i, s in enumerate(sstr):
                #Convert each observation line to a string
                sstr_str = str(s)
                # always append the first station to the non-redundant list to get things started
                if i == 0:
                    sstr_one_station.append(sstr_str)
                else: 
                    #Split the current station string and last station string in the one station list,
                    #and get their names, lats, lons, and time to compare them. 
                    sstr_split = mt.get_ob_info(sstr_str)
                    sstr_one_station_prev = str(sstr_one_station[-1])
                    sstr_one_station_prev_split = mt.get_ob_info(sstr_one_station_prev)
                    
                    # check if this station is already in the one ob per station list.                
                    if sstr_one_station_prev_split['name'] == sstr_split['name']:
                        #Convert string times into datetime objects for comparison later
                        ob_time = mt.timestamp2utc(sstr_split['time'])
                        prev_ob_time = mt.timestamp2utc(sstr_one_station_prev_split['time'])
                        #if the name of the station is 'SHIP' it isn't a unique name- check lat/lon as well
                        if sstr_split['name'] == 'SHIP':
                            #if the lat/lon differ from the previous entry and are not less than 1 degree 
                            #from the previous entry (to eliminate duplicates of moving ships), append ob to one per station list.
                            sprev_lat = sstr_one_station_prev_split['lat']
                            sprev_lon = sstr_one_station_prev_split['lon']
                            if sprev_lat != sstr_split['lat'] or sprev_lon != sstr_split['lon']:
                                #check if the lat/lon are both off by less than 1 degree- if so, assume same ship
                                if abs(float(sprev_lat)-float(sstr_split['lat']))<1 and abs(float(sprev_lon)-float(sstr_split['lon']))<1:
                                    #find observation of ship that occurs closest to desired hour
                                    if (abs(ob_time-desired_time) < abs(prev_ob_time - desired_time)):
                                        sstr_one_station[-1] = sstr_str
                                else:
                                    sstr_one_station.append(sstr_str)
                            #If the name, lat, and lon are the same, need to do a time check:
                            #if current ob has a time closer to the desired hour, replace the 
                            #previous ob with the current ob    
                            else:
                                if (abs(ob_time-desired_time) < abs(prev_ob_time - desired_time)):
                                    sstr_one_station[-1] = sstr_str
                 
                        #if the names are the same, only keep the observation closest to 00,06,12,18Z
                        #by doing a time check
                        else:
                            #if current ob has a time closer to the desired hour, replace the 
                            #previous ob with the current ob. This means of observations with the 
                            #same time distance away from the desired hour, the first one only
                            #will be kept.
                            if (abs(ob_time-desired_time) < abs(prev_ob_time - desired_time)):
                                sstr_one_station[-1] = sstr_str
                                
                    #Append a station if the name isn't in the list yet. This works
                    #because the list was previously sorted by name. 
                    else:
                        sstr_one_station.append(sstr_str)
                        
            print("Removed duplicate stations: Number of obs: "+str(len(sstr_one_station)))
            
            #This removes observations not within a window of plus or minus 1 hour
            #of the 6-hourly time.
            sstr_one_station_1hr = []
            for i, s in enumerate(sstr_one_station):
                s_split = mt.get_ob_info(s)
                #Get timestamp of observation in epoch time and convert to utc datetime object
                s_time = mt.timestamp2utc(s_split['time']) 
                #if time is less than 1 hour away from 6-hour time, add observation
                #to final list of observations
                if abs(s_time-desired_time) <= timedelta(seconds=3600):                       
                    sstr_one_station_1hr.append(s)
                    #del sstr_one_station[i]
            print("Removed times outside of 1 hour window: New number of obs: " +str(len(sstr_one_station_1hr)))                                                             
            #sstr = sstr_one_station
            #Write selected MADIS obtype variable data into CSV file
            f = open(base_dir+"/"+dates[0:6]+"/"+ob+"_"+var_short+"/"+ob+"_"+var_short+"_"+str(dates)+".txt","w")
            for s in sstr_one_station_1hr:
                f.write(s)
            f.close()
            
            #Append MADIS variable setting to combined file
            if ob == 'metar':
                f = open(base_dir+"/"+dates[0:6]+"/combined_"+var_short+"/"+var_short+"_"+str(dates)+".txt","w")
                for s in sstr_one_station_1hr:
                        f.write(s)
                f.close()
            #This is to append to an already created file instead of writing another one
            elif ob == 'maritime' or ob =='raob':
                f = open(base_dir+"/"+dates[0:6]+"/combined_"+var_short+"/"+var_short+"_"+str(dates)+".txt","a")
                for s in sstr_one_station_1hr:
                        f.write(s)
                f.close()
