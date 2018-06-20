#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 12:29:55 2018

@author: stangen
"""

from datetime import datetime
from datetime import timedelta


def make_datelist(start_date = '20130401_0000', end_date = '20130406_0000', torf=False, hour_before_after=False, timestep=21600):
    """
    Creates a list of dates. Adjust value from 21600 to create dates spaced
    differently than 6 hours. Input values in the format 'yyyymmdd_hhhh'
    Includes the last value. If selected, also gets hour before and after 
    00,06,12,18. If torf == True, formats in %Y%m%d_%H00 format, otherwise
    formats in %Y%m%d%H format.
    """  
    
    #Convert string to datetime object
    dt0 = datetime.strptime(start_date,'%Y%m%d_%H%M')
    dt1 = datetime.strptime(end_date,'%Y%m%d_%H%M')
    
    #Truncate datetime to hour precision
    dt_start = dt0.replace(minute = 0, second=0, microsecond=0)
    dt_end = dt1.replace(minute = 0, second=0, microsecond=0)
    
    #Use while loop to generate list of date/times between start and end date every 6 hours
    dt = dt_start; dt_list = []
    while (dt <= dt_end):
        if torf == True:
            if hour_before_after == True:
                dt_list.append((dt - timedelta(seconds = 3600)).strftime('%Y%m%d_%H00'))
                dt_list.append(dt.strftime('%Y%m%d_%H00'))
                dt_list.append((dt + timedelta(seconds = 3600)).strftime('%Y%m%d_%H00'))
            elif hour_before_after == False:    
                dt_list.append(dt.strftime('%Y%m%d_%H00'))
    
        #print(dt)
        elif torf == False:
            if hour_before_after == True:
                dt_list.append((dt - timedelta(seconds = 3600)).strftime('%Y%m%d%H'))
                dt_list.append(dt.strftime('%Y%m%d%H'))
                dt_list.append((dt + timedelta(seconds = 3600)).strftime('%Y%m%d%H'))
            elif hour_before_after == False:
                dt_list.append(dt.strftime('%Y%m%d%H'))
    	#dt = dt + timedelta(0,3600)
        dt = dt + timedelta(0,timestep)
        
    return dt_list

def make_datetimelist(start_date = datetime(2013,4,1,0), end_date = datetime(2013,4,6,0), timestep=21600):
    """
    Similar to make_datelist, but this accepts and returns a list of datetime objects instead of 
    a list of strings.
    """
    
    #Truncate datetime to hour precision
    dt_start = start_date.replace(minute = 0, second=0, microsecond=0)
    dt_end = end_date.replace(minute = 0, second=0, microsecond=0)
    
    #Use while loop to generate list of date/times between start and end date every 6 hours
    dt = dt_start; dt_list = []
    while (dt <= dt_end):
        dt_list.append(dt)
        dt = dt + timedelta(0,timestep)
    
    return dt_list


def get_ob_info(sstr_str):
    """    
    Split the station string and return a dictionary containing their names, 
    lats, lons, time, elevation, and ob value. 
    **********Changed these to floats 6/19, not 100% sure how this will affect other code*****
    """
    
    sstr_split = sstr_str.split(',')
    sstr_name = sstr_split[0]
    sstr_lat = float(sstr_split[1])
    sstr_lon = float(sstr_split[2])
    sstr_elev = float(sstr_split[3]) 
    sstr_time = float(sstr_split[4])
    sstr_ob = float(sstr_split[5])
    
    info = {'name':sstr_name, 'lat':sstr_lat, 'lon':sstr_lon, 'time':sstr_time, 'elev':sstr_elev, 'ob':sstr_ob}
    return info

def timestamp2utc(timestamp):
    """
    Takes a timestamp in epoch time and returns a datetime object in utc.
    """
    return datetime.utcfromtimestamp(int(round(float(timestamp))))

def month_dict():
    """
    Returns a dictionary for month naming conventions
    """
    m_dict = {
           '01' : 'jan',
           '02' : 'feb',
           '03' : 'mar',
           '04' : 'apr',
           '05' : 'may',
           '06' : 'jun',
           '07' : 'jul',
           '08' : 'aug',
           '09' : 'sep',
           '10' : 'oct',
           '11' : 'nov',
           '12' : 'dec'           
           }
    return m_dict
    
    
