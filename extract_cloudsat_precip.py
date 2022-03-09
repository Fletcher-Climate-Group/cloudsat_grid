#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 23:19:00 2017

@author: Fraser King, March 2018
"""

import sys,os,csv,glob
import numpy as np
import pandas as pd
import h5py
import warnings
import math
warnings.filterwarnings("ignore")

# Check for arguments
if len(sys.argv) < 3:
    print ("ERROR: Please include an inpath, filename & outpath in this program call as arguments")
    sys.exit()

# Read inpath, filename & outpath
inpath = sys.argv[1]
filename = sys.argv[2]
outpath = sys.argv[3]

# Configuration variables
grid_size = 100
use_box = False
station_name = "cambridge"
station_lat = 69.10805556
station_lon = -105.13833333

# Helper functions
def haversine_np(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km
    
def get_box_points(lat, lon, offset):
    lat_offset  = (offset / 6371) * (180 / math.pi)
    lon_offset = (offset / 6371) * (180 / math.pi) / math.cos(lat * math.pi/180)
    top_right = (lat + lat_offset, lon + lon_offset)
    bottom_right = (lat - lat_offset, lon + lon_offset)
    top_left = (lat + lat_offset, lon - lon_offset)
    bottom_left = (lat - lat_offset, lon - lon_offset)
    return (top_right, bottom_right, bottom_left, top_left)

# Get box info
if use_box:
    print("Grid box: ", get_box_points(station_lat, station_lon, grid_size))

with h5py.File(inpath,'r') as hdf:
    #print(filepath)
    #extract year and day of year
    year = str(filename)[0:4]
    day_of_year = str(filename)[4:7]
    #extract hdf5 data using keys
    ls = list(hdf.keys()) #identify individual keys
    base_layer = list(hdf.items()) #top branch of hdf
    #print(base_layer)
    group1 = hdf.get('2C-RAIN-PROFILE')
    group1_items = list(group1.items())
    #print(group1_items)
    #extract lat, long, profile time and utc start time from geolocation fields
    geolocation_fields = group1.get('Geolocation Fields')
    #geolocation_fields = np.array(geolocation_fields)
    geolocation_fields_items = list(geolocation_fields.items())
    longitude = geolocation_fields.get('Longitude')
    #print(longitude)
    longitude = np.array(longitude)
    longitude = longitude.astype(np.float64)
    #total number of rows
    count = longitude.size
    
    latitude = geolocation_fields.get('Latitude')
    latitude = np.array(latitude)
    latitude = latitude.astype(np.float64)
    utc_time = geolocation_fields.get('UTC_start')
    utc_time = np.array(utc_time)
    profile_time = geolocation_fields.get('Profile_time')
    profile_time = np.array(profile_time) 

    data_fields = group1.get('Data Fields')
    data_fields_items = list(data_fields.items())
    
    precip_rate = data_fields.get('rain_rate')
    precip_rate = np.array(precip_rate)
    precip_rate = precip_rate.astype(np.float64)
    
    precip_rate_max = data_fields.get('rain_rate_uncertainity')
    precip_rate_max = np.array(precip_rate_max)
    precip_rate_max = precip_rate_max.astype(np.float64)

    precip_rate_min = data_fields.get('rain_status_flag')
    precip_rate_min = np.array(precip_rate_min)
    precip_rate_min = precip_rate_min.astype(np.float64)
    
    land_sea_flag = data_fields.get('Navigation_land_sea_flag')
    land_sea_flag = np.array(land_sea_flag)
    land_sea_flag = land_sea_flag.astype(np.float64)
    
    retrieval_info = data_fields.get('precip_flag')
    retrieval_info = np.array(retrieval_info)
    retrieval_info = retrieval_info.astype(np.float64)
    
    data_quality = data_fields.get('Data_quality')
    data_quality = np.array(data_quality)
    data_quality = data_quality.astype(np.float64)
    
    data_status = data_fields.get('Data_status')
    data_status = np.array(data_status)
    data_status = data_status.astype(np.float64)
    
    data_target_id = data_fields.get('Data_targetID')
    data_target_id = np.array(data_target_id)
    data_target_id = data_target_id.astype(np.float64)
    
    #assign negative values as nan
    precip_rate[precip_rate < 0] = np.nan

    precip_rate = list(precip_rate)
    precip_rate_min = list(precip_rate_min)
    #precip_rate_max = list(precip_rate_max)
    data_quality = list(data_quality)
    data_status = list(data_status)
    data_target_id = list(data_target_id)
    land_sea_flag = list(land_sea_flag)
    retrieval_info = list(retrieval_info)

    latitude = list(latitude)
    longitude = list(longitude)
    profile_time = [i for sub in profile_time for i in sub]
    utc_time = [i for sub in utc_time for i in sub]
    #utc time = profile time + utc start time
    new_list = [x+utc_time for x in profile_time]
    utc_time = [i for sub in new_list for i in sub]

    df_station = pd.DataFrame({"Year" : year,
                              "Day_of_Year" : day_of_year,
                              "Cloudsat_Lat" : latitude,
                              "Cloudsat_Lon" : longitude,
                              "UTC_Time" : utc_time,
                              "Rain_Rate" : precip_rate,
                              "Rain_Rate_Uncertainty" : precip_rate_min,
                              "Rain_Status_Flag" : precip_rate_max,
                              "Land_Sea_Flag" : land_sea_flag,
                              "Data_Quality" : data_quality,
                              "Data_Status" : data_status,
                              "Data_Target_ID": data_target_id,
                              "Retrieval_Info" : retrieval_info})
                              
    rect = get_box_points(station_lat, station_lon, grid_size)

    df_station['Ground_Lat'] = station_lat
    df_station['Ground_Lon'] = station_lon

    df_station = df_station[df_station.Rain_Rate > 0]
    df_station5 = []

    if not(use_box):
        df_station['Distance'] = haversine_np(df_station['Cloudsat_Lon'], df_station['Cloudsat_Lat'], df_station['Ground_Lon'], df_station['Ground_Lat'])
        df_station5 = df_station[df_station.Distance < grid_size]
    else:
        df_station2 = df_station[df_station.Cloudsat_Lat < rect[3][0]]
        df_station3 = df_station2[df_station2.Cloudsat_Lat > rect[2][0]]
        df_station4 = df_station3[df_station3.Cloudsat_Lon > rect[3][1]]
        df_station5 = df_station4[df_station4.Cloudsat_Lon < rect[0][1]]

    # #add file name to keep in track of date
    df_station5['File_Name'] = filename
    # #specify the order of columns                
    df_station5 = df_station5[["Year", "Day_of_Year", "Cloudsat_Lat", "Cloudsat_Lon", "UTC_Time", "Rain_Rate", "Rain_Rate_Uncertainty", "Rain_Status_Flag", "Land_Sea_Flag", "Data_Quality", "Data_Status", "Data_Target_ID", "Retrieval_Info", "File_Name"]]

    # check if empty result
    if not df_station5.empty:
        print ("Overpass hit")
        df_station5.to_csv((outpath + "/" + station_name + "/" + filename + ".csv"), index=True)



