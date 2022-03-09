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
grid_size = 200
use_box = True
#station_name1 = "cambridge"
#station_lat1 = 69.10805556
#station_lon1 = -105.13833333
#
#station_name2 = "resolute"
#station_lat2 = 74.71694444
#station_lon2 = -94.96944444
#
#station_name3 = "eureka"
#station_lat3 = 79.99444444
#station_lon3 = -85.81194444
#
#station_name4 = "iqaluit"
#station_lat4 = 63.74722222
#station_lon4 = -68.54444444
#
#station_name5 = "summit"
#station_lat5 = 72.579583
#station_lon5 = -38.459186

station_name1 = "owen"
station_lat1 = 44.568934
station_lon1 = -80.936059

station_name2 = "cochrane"
station_lat2 = 49.083099
station_lon2 = -81.029470

station_name3 = "thunder"
station_lat3 = 48.406576
station_lon3 = -89.276237

station_name4 = "sudbury"
station_lat4 = 46.483216
station_lon4 = -80.993322

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
#if use_box:
#    print("Grid box: ", get_box_points(station_lat, station_lon, grid_size))

with h5py.File(inpath,'r') as hdf:
    year = str(filename)[0:4]
    day_of_year = str(filename)[4:7]
    ls = list(hdf.keys())
    base_layer = list(hdf.items())
    group1 = hdf.get('2C-SNOW-PROFILE')
    group1_items = list(group1.items())
    geolocation_fields = group1.get('Geolocation Fields')
    geolocation_fields_items = list(geolocation_fields.items())
    longitude = geolocation_fields.get('Longitude')
    longitude = np.array(longitude)
    longitude = longitude.astype(np.float64)
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
    surface_snowfall_rate = data_fields.get('snowfall_rate_sfc')
    surface_snowfall_rate = np.array(surface_snowfall_rate)
    surface_snowfall_rate = surface_snowfall_rate.astype(np.float64)

    surface_snowfall_rate_uncertainty = data_fields.get('snowfall_rate_sfc_uncert')
    surface_snowfall_rate_uncertainty = np.array(surface_snowfall_rate_uncertainty)
    surface_snowfall_rate_uncertainty = surface_snowfall_rate_uncertainty.astype(np.float64)

    surface_snowfall_rate_confidence = data_fields.get('snowfall_rate_sfc_confidence')
    surface_snowfall_rate_confidence = np.array(surface_snowfall_rate_confidence)
    surface_snowfall_rate_confidence = surface_snowfall_rate_confidence.astype(np.float64)

    data_quality = data_fields.get('Data_quality')
    data_quality = np.array(data_quality)
    data_quality = data_quality.astype(np.float64)

    data_status = data_fields.get('Data_status')
    data_status = np.array(data_status)
    data_status = data_status.astype(np.float64)

    data_target_id = data_fields.get('Data_targetID')
    data_target_id = np.array(data_target_id)
    data_target_id = data_target_id.astype(np.float64)

    snow_retrieval_status = data_fields.get('snow_retrieval_status')
    snow_retrieval_status = np.array(snow_retrieval_status)
    snow_retrieval_status = snow_retrieval_status.astype(np.float64)

    dem_elevation = geolocation_fields.get('DEM_elevation')
    dem_elevation = np.array(dem_elevation)

    snowfall_rate = data_fields.get('snowfall_rate')
    snowfall_rate = np.array(snowfall_rate)
    snowfall_rate = snowfall_rate.astype(np.float64)
    elevation = geolocation_fields.get('Height')
    elevation = np.array(elevation)

    near_surface_snowlayer = np.zeros(shape = (count,125))

    snowfall_rate[snowfall_rate <0] = np.nan
    near_surface_snowlayer[near_surface_snowlayer<0] = np.nan
    near_surface_snowlayer[near_surface_snowlayer==0] = np.nan

    snowfall_rate_avg = np.nanmean(snowfall_rate, axis = 1)
    snowfall_rate_std = np.nanstd(snowfall_rate, axis = 1)
    near_surface_snow_layer_avg = np.nanmean(near_surface_snowlayer, axis = 1)
    near_surface_snowlayer_std = np.nanstd(near_surface_snowlayer, axis = 1)

    snowfall_rate_avg = list(snowfall_rate_avg)
    near_surface_snow_layer_avg = list(near_surface_snow_layer_avg)
    snowfall_rate_std = list(snowfall_rate_std)
    near_surface_snowlayer_std = list(near_surface_snowlayer_std)

    surface_snowfall_rate = list(surface_snowfall_rate)
    data_quality = list(data_quality)
    data_status = list(data_status)
    data_target_id = list(data_target_id)
    snow_retrieval_status = list(snow_retrieval_status)
    latitude = list(latitude)
    longitude = list(longitude)
    dem_elevation = list(dem_elevation)
    surface_snowfall_rate_uncertainty = list(surface_snowfall_rate_uncertainty)
    surface_snowfall_rate_confidence = list(surface_snowfall_rate_confidence)
    profile_time = [i for sub in profile_time for i in sub]
    utc_time = [i for sub in utc_time for i in sub]
    new_list = [x+utc_time for x in profile_time]
    utc_time = [i for sub in new_list for i in sub]

    for x in range(0, 4):
        station_name = station_name1
        station_lat = station_lat1
        station_lon = station_lon1

        if x == 1:
            station_name = station_name2
            station_lat = station_lat2
            station_lon = station_lon2
        elif x == 2:
            station_name = station_name3
            station_lat = station_lat3
            station_lon = station_lon3
        elif x == 3:
            station_name = station_name4
            station_lat = station_lat4
            station_lon = station_lon4


        df_station = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std})

        rect = get_box_points(station_lat, station_lon, grid_size)
        df_station['Ground_Lat'] = station_lat
        df_station['Ground_Lon'] = station_lon
        #df_station = df_station[df_station.Surface_Snowfall_Rate > 0]
        df_station5 = []

        if not(use_box):
            df_station['Distance'] = haversine_np(df_station['Cloudsat_Lon'], df_station['Cloudsat_Lat'], df_station['Ground_Lon'], df_station['Ground_Lat'])
            df_station5 = df_station[df_station.Distance < grid_size]
        else:
            df_station2 = df_station[df_station.Cloudsat_Lat < rect[3][0]]
            df_station3 = df_station2[df_station2.Cloudsat_Lat > rect[2][0]]
            df_station4 = df_station3[df_station3.Cloudsat_Lon > rect[3][1]]
            df_station5 = df_station4[df_station4.Cloudsat_Lon < rect[0][1]]

        df_station5['File_Name'] = filename
        df_station5 = df_station5[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD", "File_Name"]]

        # check if empty result
        if not df_station5.empty:
            print ("Overpass hit")
            df_station5.to_csv((outpath + "/" + station_name + "/" + filename + ".csv"), index=True)
