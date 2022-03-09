#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 23:19:00 2017

@author: Rithwik
"""

#code to extract snowfall rate from hdf5 cloudsat 2c-snow profile

import sys,os,csv,glob
import numpy as np
import pandas as pd
import h5py

# Check for arguments
if len(sys.argv) < 3:
    print ("ERROR: Please include an inpath, filename & outpath in this program call as arguments")
    sys.exit()

# Read inpath, filename & outpath
inpath = sys.argv[1]
filename = sys.argv[2]
outpath = sys.argv[3]

with h5py.File(inpath,'r') as hdf:
    #print(filepath)
    #extract year and day of year
    year = str(filename)[0:4]
    day_of_year = str(filename)[4:7]
    #extract hdf5 data using keys
    ls = list(hdf.keys()) #identify individual keys
    #print(ls)
    #2c_snow_profile = hdf.get('2C-SNOW-PROFILE')
    #data1_a = np.array(data1)
    base_layer = list(hdf.items()) #top branch of hdf
    #print(base_layer)
    group1 = hdf.get('2C-SNOW-PROFILE')
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

    #extract surfacesnowfall rate from datafields
    data_fields = group1.get('Data Fields')
    data_fields_items = list(data_fields.items())
    surface_snowfall_rate = data_fields.get('snowfall_rate_sfc')
    surface_snowfall_rate = np.array(surface_snowfall_rate)
    surface_snowfall_rate = surface_snowfall_rate.astype(np.float64)
    #surface snowfall uncertainty
    surface_snowfall_rate_uncertainty = data_fields.get('snowfall_rate_sfc_uncert')
    surface_snowfall_rate_uncertainty = np.array(surface_snowfall_rate_uncertainty)
    surface_snowfall_rate_uncertainty = surface_snowfall_rate_uncertainty.astype(np.float64)
    #surface snowfall confidence
    surface_snowfall_rate_confidence = data_fields.get('snowfall_rate_sfc_confidence')
    surface_snowfall_rate_confidence = np.array(surface_snowfall_rate_confidence)
    surface_snowfall_rate_confidence = surface_snowfall_rate_confidence.astype(np.float64)
    #data quality
    data_quality = data_fields.get('Data_quality')
    data_quality = np.array(data_quality)
    data_quality = data_quality.astype(np.float64)
    #data stats
    data_status = data_fields.get('Data_status')
    data_status = np.array(data_status)
    data_status = data_status.astype(np.float64)
    #data targetid
    data_target_id = data_fields.get('Data_targetID')
    data_target_id = np.array(data_target_id)
    data_target_id = data_target_id.astype(np.float64)
    #snow retrieval status
    snow_retrieval_status = data_fields.get('snow_retrieval_status')
    snow_retrieval_status = np.array(snow_retrieval_status)
    snow_retrieval_status = snow_retrieval_status.astype(np.float64)
    #dem elevation from geolocation
    dem_elevation = geolocation_fields.get('DEM_elevation')
    dem_elevation = np.array(dem_elevation)
    
    
    #extract snowfall rate and elevation for vertical bins
    snowfall_rate = data_fields.get('snowfall_rate')
    #snowfall_rate_a = list(snowfall_rate)
    snowfall_rate = np.array(snowfall_rate)
    snowfall_rate = snowfall_rate.astype(np.float64)
    #elevation
    elevation = geolocation_fields.get('Height')
    elevation = np.array(elevation)
    
    #extracting snowbin parameters
    snowfall_rate_near_surface_bin = np.zeros(shape=(count))    
    snowfall_rate_near_surface_endbin = np.zeros(shape=(count))
    elevation_near_surface_bin = np.zeros(shape= (count))
    elevation_near_surface_endbin = np.zeros(shape=(count))
    near_surface_bin = np.zeros(shape=(count))
    near_surface_endbin = np.zeros(shape=(count))
    #new array with near surface snow layer 
    near_surface_snowlayer = np.zeros(shape = (count,125))
    
    #loop to extract near surface snowfall layer parameters
    for x in range (0, count): #to extract near surface bin data
        for y in range (124, 0, -1):
            if (snowfall_rate[x, y] != -999 ):
                snowfall_rate_near_surface_bin[x] = snowfall_rate[x,y]
                elevation_near_surface_bin[x] = elevation[x,y]
                near_surface_bin[x] = y
                for l in range (y, 0, -1): #to extract near surface end bin data
                    near_surface_snowlayer[x,l] = snowfall_rate[x,l]
                    if (snowfall_rate[x,l] == -999):
                        snowfall_rate_near_surface_endbin[x] = snowfall_rate[x,l+1]
                        elevation_near_surface_endbin [x] = elevation[x,l+1]
                        near_surface_endbin[x] = l+1
                        break
                break
    
    #to extract mid snow layer parameters
    midlayer_snowfall_rate = np.zeros(shape = (count))
    difference_bin = near_surface_bin-near_surface_endbin
    midlayer_elevation = np.zeros(shape = count)
    for b in range (0, count):
        if ((near_surface_endbin[b] !=0) and difference_bin[b]>1):
            if (difference_bin[b] % 2 == 0):
                midlayer_snowfall_rate[b] = snowfall_rate[b, int((near_surface_bin[b]+near_surface_endbin[b])/2)]
                midlayer_elevation[b] = elevation[b, int((near_surface_bin[b]+near_surface_endbin[b])/2)]
            else:
                midlayer_snowfall_rate[b] = (float(snowfall_rate[b, int((near_surface_bin[b]-int(difference_bin[b]/2)))] + snowfall_rate[b, int((near_surface_endbin[b]+int(difference_bin[b]/2)))]))/2
                midlayer_elevation[b] = (float(elevation[b, int((near_surface_bin[b]-int(difference_bin[b]/2)))] + elevation[b, int((near_surface_endbin[b]+int(difference_bin[b]/2)))]))/2
       
    #assign negative values as nan
    snowfall_rate[snowfall_rate <0] = np.nan
    near_surface_snowlayer[near_surface_snowlayer<0] = np.nan
    near_surface_snowlayer[near_surface_snowlayer==0] = np.nan

    #calculate snowfall rate averaged over 125 vertical bins
    snowfall_rate_avg = np.nanmean(snowfall_rate, axis = 1) #axis = 1 => horizontal averaging
    snowfall_rate_std = np.nanstd(snowfall_rate, axis = 1)
    near_surface_snow_layer_avg = np.nanmean(near_surface_snowlayer, axis = 1)
    near_surface_snowlayer_std = np.nanstd(near_surface_snowlayer, axis = 1)

    #assign negative values as nan
    #snowfall_rate[snowfall_rate <0] = np.nan
    #create vertical bins and store snowfall rate and elevation for each bin against bin number
    #bin_number=range(0,125)
    #snowfall rate
    #for k in range(0,125):
    #    vars()["bin"+str(bin_number[k])]=list(snowfall_rate[:,k])
    #    vars()["ele"+str(bin_number[k])]=list(elevation[:,k])
    #calculate snowfall rate averaged over 125 vertical bins
    #snowfall_rate_avg = np.nanmean(snowfall_rate, axis = 1) #axis = 1 => horizontal averaging
    #snowfall_rate_ground = list(snowfall_rate_avg) #snowfall_rate in float32 format
    #snowfall_rate_ground = list(snowfall_rate[:,97]) #snowfall_rate in float32 format

    #merge tuples
    #snowfall_rate_ground = [i for sub in snowfall_rate_ground for i in sub]
    midlayer_snowfall_rate = list(midlayer_snowfall_rate)
    midlayer_elevation = list(midlayer_elevation)
    snowfall_rate_near_surface_bin = list(snowfall_rate_near_surface_bin)
    snowfall_rate_near_surface_endbin = list(snowfall_rate_near_surface_endbin)
    elevation_near_surface_bin = list(elevation_near_surface_bin)
    elevation_near_surface_endbin = list(elevation_near_surface_endbin)
    near_surface_bin = list(near_surface_bin)
    near_surface_endbin = list(near_surface_endbin)
    snowfall_rate_avg = list(snowfall_rate_avg)
    near_surface_snow_layer_avg = list(near_surface_snow_layer_avg)
    snowfall_rate_std = list(snowfall_rate_std)
    near_surface_snowlayer_std = list(near_surface_snowlayer_std)

    #surface_snowfall_rate = [i for sub in surface_snowfall_rate for i in sub]
    surface_snowfall_rate = list(surface_snowfall_rate)
    data_quality = list(data_quality)
    data_status = list(data_status)
    data_target_id = list(data_target_id)
    snow_retrieval_status = list(snow_retrieval_status)
    #latitude = [i for sub in latitude for i in sub]
    latitude = list(latitude)
    #longitude = [i for sub in longitude for i in sub]
    longitude = list(longitude)
    dem_elevation = list(dem_elevation)
    surface_snowfall_rate_uncertainty = list(surface_snowfall_rate_uncertainty)
    surface_snowfall_rate_confidence = list(surface_snowfall_rate_confidence)
    profile_time = [i for sub in profile_time for i in sub]
    utc_time = [i for sub in utc_time for i in sub]
    #utc time = profile time + utc start time
    new_list = [x+utc_time for x in profile_time]
    utc_time = [i for sub in new_list for i in sub]

    #create pandas data frame
    #dataframe for surface snowfall 
    #df = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Bin1" : bin0, "Bin5" : bin4, "Bin10" : bin9, "Bin15" : bin14, "Bin20" : bin19, "Bin25" : bin24, "Bin30" : bin29, "Bin35" : bin34, "Bin40" : bin39, "Bin45" : bin44, "Bin50" : bin49, "Bin55" : bin54, "Bin60" : bin59, "Bin65" : bin64, "Bin70" : bin69, "Bin75" : bin74, "Bin80" : bin79, "Bin85" : bin84, "Bin90" : bin89, "Bin95" : bin94, "Bin100" : bin99, "Bin105" : bin104, "Bin110" : bin109, "Bin115" : bin114, "Bin120" : bin119, "Bin125" : bin124})
    
    #function to calculate distance between cloudsat and ground
    def haversine_np(lon1, lat1, lon2, lat2):  
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
        c = 2 * np.arcsin(np.sqrt(a))
        km = 6367 * c
        return km
    
    
    #eureka
    df_eureka = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std,
                                   "Snowfall_Rate_Baselayer" : snowfall_rate_near_surface_bin, "Elevation_Baselayer" : elevation_near_surface_bin, "Bin_Number_Baselayer" : near_surface_bin, "Snowfall_Rate_Midlayer" : midlayer_snowfall_rate, "Elevation_Midlayer" : midlayer_elevation, "Snowfall_Rate_TopLayer": snowfall_rate_near_surface_endbin, "Elevation_TopLayer" : elevation_near_surface_endbin, "Bin_Number_TopLayer" : near_surface_endbin})
    df_eureka['Ground_Lat'] = 79.98916667
    df_eureka['Ground_Lon'] = -85.93388889
    #remove observations with -ve snowfall rates
    df_eureka = df_eureka[df_eureka.Surface_Snowfall_Rate > 0]
    #dataframe distance calculation
    df_eureka['Distance'] = haversine_np(df_eureka['Cloudsat_Lon'], df_eureka['Cloudsat_Lat'], df_eureka['Ground_Lon'], df_eureka['Ground_Lat'])
    #specifying distance threshold of 50km
    df_eureka = df_eureka[df_eureka.Distance < 100]
    #add file name to keep in track of date
    df_eureka['File_Name'] = filename
    #specify the order of columns                
    #df = df[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Bin1", "Bin5", "Bin10", "Bin15", "Bin20", "Bin25", "Bin30", "Bin35", "Bin40", "Bin45", "Bin50", "Bin55", "Bin60", "Bin65", "Bin70", "Bin75", "Bin80", "Bin85", "Bin90", "Bin95", "Bin100", "Bin105", "Bin110", "Bin115", "Bin120", "Bin125", "File_Name"]]
    df_eureka = df_eureka[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD","Snowfall_Rate_Baselayer", "Elevation_Baselayer", "Bin_Number_Baselayer", "Snowfall_Rate_Midlayer", "Elevation_Midlayer", "Snowfall_Rate_TopLayer", "Elevation_TopLayer", "Bin_Number_TopLayer", "File_Name"]]
    
    
    #eureka-a
    df_eureka_a = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std,
                                   "Snowfall_Rate_Baselayer" : snowfall_rate_near_surface_bin, "Elevation_Baselayer" : elevation_near_surface_bin, "Bin_Number_Baselayer" : near_surface_bin, "Snowfall_Rate_Midlayer" : midlayer_snowfall_rate, "Elevation_Midlayer" : midlayer_elevation, "Snowfall_Rate_TopLayer": snowfall_rate_near_surface_endbin, "Elevation_TopLayer" : elevation_near_surface_endbin, "Bin_Number_TopLayer" : near_surface_endbin})
    df_eureka_a['Ground_Lat'] = 79.99444444
    df_eureka_a['Ground_Lon'] = -85.81194444
    #remove observations with -ve snowfall rates
    df_eureka_a = df_eureka_a[df_eureka_a.Surface_Snowfall_Rate > 0]
    #dataframe distance calculation
    df_eureka_a['Distance'] = haversine_np(df_eureka_a['Cloudsat_Lon'], df_eureka_a['Cloudsat_Lat'], df_eureka_a['Ground_Lon'], df_eureka_a['Ground_Lat'])
    #specifying distance threshold of 50km
    df_eureka_a = df_eureka_a[df_eureka_a.Distance < 100]
    #add file name to keep in track of date
    df_eureka_a['File_Name'] = filename
    #specify the order of columns                
    #df = df[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Bin1", "Bin5", "Bin10", "Bin15", "Bin20", "Bin25", "Bin30", "Bin35", "Bin40", "Bin45", "Bin50", "Bin55", "Bin60", "Bin65", "Bin70", "Bin75", "Bin80", "Bin85", "Bin90", "Bin95", "Bin100", "Bin105", "Bin110", "Bin115", "Bin120", "Bin125", "File_Name"]]
    df_eureka_a = df_eureka_a[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD","Snowfall_Rate_Baselayer", "Elevation_Baselayer", "Bin_Number_Baselayer", "Snowfall_Rate_Midlayer", "Elevation_Midlayer", "Snowfall_Rate_TopLayer", "Elevation_TopLayer", "Bin_Number_TopLayer", "File_Name"]]
    
    
    #resolute
    df_resolute = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std,
                                   "Snowfall_Rate_Baselayer" : snowfall_rate_near_surface_bin, "Elevation_Baselayer" : elevation_near_surface_bin, "Bin_Number_Baselayer" : near_surface_bin, "Snowfall_Rate_Midlayer" : midlayer_snowfall_rate, "Elevation_Midlayer" : midlayer_elevation, "Snowfall_Rate_TopLayer": snowfall_rate_near_surface_endbin, "Elevation_TopLayer" : elevation_near_surface_endbin, "Bin_Number_TopLayer" : near_surface_endbin})
    df_resolute['Ground_Lat'] = 74.71694444
    df_resolute['Ground_Lon'] = -94.96944444
    #remove observations with -ve snowfall rates
    df_resolute = df_resolute[df_resolute.Surface_Snowfall_Rate > 0]
    #dataframe distance calculation
    df_resolute['Distance'] = haversine_np(df_resolute['Cloudsat_Lon'], df_resolute['Cloudsat_Lat'], df_resolute['Ground_Lon'], df_resolute['Ground_Lat'])
    #specifying distance threshold of 50km
    df_resolute = df_resolute[df_resolute.Distance < 100]
    #add file name to keep in track of date
    df_resolute['File_Name'] = filename
    #specify the order of columns                
    #df = df[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Bin1", "Bin5", "Bin10", "Bin15", "Bin20", "Bin25", "Bin30", "Bin35", "Bin40", "Bin45", "Bin50", "Bin55", "Bin60", "Bin65", "Bin70", "Bin75", "Bin80", "Bin85", "Bin90", "Bin95", "Bin100", "Bin105", "Bin110", "Bin115", "Bin120", "Bin125", "File_Name"]]
    df_resolute = df_resolute[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD","Snowfall_Rate_Baselayer", "Elevation_Baselayer", "Bin_Number_Baselayer", "Snowfall_Rate_Midlayer", "Elevation_Midlayer", "Snowfall_Rate_TopLayer", "Elevation_TopLayer", "Bin_Number_TopLayer", "File_Name"]]
    
    
    #egbert
    df_egbert = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std,
                                   "Snowfall_Rate_Baselayer" : snowfall_rate_near_surface_bin, "Elevation_Baselayer" : elevation_near_surface_bin, "Bin_Number_Baselayer" : near_surface_bin, "Snowfall_Rate_Midlayer" : midlayer_snowfall_rate, "Elevation_Midlayer" : midlayer_elevation, "Snowfall_Rate_TopLayer": snowfall_rate_near_surface_endbin, "Elevation_TopLayer" : elevation_near_surface_endbin, "Bin_Number_TopLayer" : near_surface_endbin})
    df_egbert['Ground_Lat'] = 44.23333333
    df_egbert['Ground_Lon'] = -79.78333333
    #remove observations with -ve snowfall rates
    df_egbert = df_egbert[df_egbert.Surface_Snowfall_Rate > 0]
    #dataframe distance calculation
    df_egbert['Distance'] = haversine_np(df_egbert['Cloudsat_Lon'], df_egbert['Cloudsat_Lat'], df_egbert['Ground_Lon'], df_egbert['Ground_Lat'])
    #specifying distance threshold of 50km
    df_egbert = df_egbert[df_egbert.Distance < 100]
    #add file name to keep in track of date
    df_egbert['File_Name'] = filename
    #specify the order of columns                
    #df = df[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Bin1", "Bin5", "Bin10", "Bin15", "Bin20", "Bin25", "Bin30", "Bin35", "Bin40", "Bin45", "Bin50", "Bin55", "Bin60", "Bin65", "Bin70", "Bin75", "Bin80", "Bin85", "Bin90", "Bin95", "Bin100", "Bin105", "Bin110", "Bin115", "Bin120", "Bin125", "File_Name"]]
    df_egbert = df_egbert[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD","Snowfall_Rate_Baselayer", "Elevation_Baselayer", "Bin_Number_Baselayer", "Snowfall_Rate_Midlayer", "Elevation_Midlayer", "Snowfall_Rate_TopLayer", "Elevation_TopLayer", "Bin_Number_TopLayer", "File_Name"]]
    
    
    #cambridge
    df_cambridge = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std,
                                   "Snowfall_Rate_Baselayer" : snowfall_rate_near_surface_bin, "Elevation_Baselayer" : elevation_near_surface_bin, "Bin_Number_Baselayer" : near_surface_bin, "Snowfall_Rate_Midlayer" : midlayer_snowfall_rate, "Elevation_Midlayer" : midlayer_elevation, "Snowfall_Rate_TopLayer": snowfall_rate_near_surface_endbin, "Elevation_TopLayer" : elevation_near_surface_endbin, "Bin_Number_TopLayer" : near_surface_endbin})
    df_cambridge['Ground_Lat'] = 69.10805556
    df_cambridge['Ground_Lon'] = -105.13833333
    #remove observations with -ve snowfall rates
    df_cambridge = df_cambridge[df_cambridge.Surface_Snowfall_Rate > 0]
    #dataframe distance calculation
    df_cambridge['Distance'] = haversine_np(df_cambridge['Cloudsat_Lon'], df_cambridge['Cloudsat_Lat'], df_cambridge['Ground_Lon'], df_cambridge['Ground_Lat'])
    #specifying distance threshold of 50km
    df_cambridge = df_cambridge[df_cambridge.Distance < 100]
    #add file name to keep in track of date
    df_cambridge['File_Name'] = filename
    #specify the order of columns                
    #df = df[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Bin1", "Bin5", "Bin10", "Bin15", "Bin20", "Bin25", "Bin30", "Bin35", "Bin40", "Bin45", "Bin50", "Bin55", "Bin60", "Bin65", "Bin70", "Bin75", "Bin80", "Bin85", "Bin90", "Bin95", "Bin100", "Bin105", "Bin110", "Bin115", "Bin120", "Bin125", "File_Name"]]
    df_cambridge = df_cambridge[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD","Snowfall_Rate_Baselayer", "Elevation_Baselayer", "Bin_Number_Baselayer", "Snowfall_Rate_Midlayer", "Elevation_Midlayer", "Snowfall_Rate_TopLayer", "Elevation_TopLayer", "Bin_Number_TopLayer", "File_Name"]]
    
    
    #Iqaluit
    df_iqaluit = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std,
                                   "Snowfall_Rate_Baselayer" : snowfall_rate_near_surface_bin, "Elevation_Baselayer" : elevation_near_surface_bin, "Bin_Number_Baselayer" : near_surface_bin, "Snowfall_Rate_Midlayer" : midlayer_snowfall_rate, "Elevation_Midlayer" : midlayer_elevation, "Snowfall_Rate_TopLayer": snowfall_rate_near_surface_endbin, "Elevation_TopLayer" : elevation_near_surface_endbin, "Bin_Number_TopLayer" : near_surface_endbin})
    df_iqaluit['Ground_Lat'] = 63.74722222
    df_iqaluit['Ground_Lon'] = -68.54444444
    #remove observations with -ve snowfall rates
    df_iqaluit = df_iqaluit[df_iqaluit.Surface_Snowfall_Rate > 0]
    #dataframe distance calculation
    df_iqaluit['Distance'] = haversine_np(df_iqaluit['Cloudsat_Lon'], df_iqaluit['Cloudsat_Lat'], df_iqaluit['Ground_Lon'], df_iqaluit['Ground_Lat'])
    #specifying distance threshold of 50km
    df_iqaluit = df_iqaluit[df_iqaluit.Distance < 100]
    #add file name to keep in track of date
    df_iqaluit['File_Name'] = filename
    #specify the order of columns                
    #df = df[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Bin1", "Bin5", "Bin10", "Bin15", "Bin20", "Bin25", "Bin30", "Bin35", "Bin40", "Bin45", "Bin50", "Bin55", "Bin60", "Bin65", "Bin70", "Bin75", "Bin80", "Bin85", "Bin90", "Bin95", "Bin100", "Bin105", "Bin110", "Bin115", "Bin120", "Bin125", "File_Name"]]
    df_iqaluit = df_iqaluit[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD","Snowfall_Rate_Baselayer", "Elevation_Baselayer", "Bin_Number_Baselayer", "Snowfall_Rate_Midlayer", "Elevation_Midlayer", "Snowfall_Rate_TopLayer", "Elevation_TopLayer", "Bin_Number_TopLayer", "File_Name"]]
    
    
    #alert
    df_alert = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std,
                                   "Snowfall_Rate_Baselayer" : snowfall_rate_near_surface_bin, "Elevation_Baselayer" : elevation_near_surface_bin, "Bin_Number_Baselayer" : near_surface_bin, "Snowfall_Rate_Midlayer" : midlayer_snowfall_rate, "Elevation_Midlayer" : midlayer_elevation, "Snowfall_Rate_TopLayer": snowfall_rate_near_surface_endbin, "Elevation_TopLayer" : elevation_near_surface_endbin, "Bin_Number_TopLayer" : near_surface_endbin})
    df_alert['Ground_Lat'] = 82.5
    df_alert['Ground_Lon'] = -62.33333333
    #remove observations with -ve snowfall rates
    df_alert = df_alert[df_alert.Surface_Snowfall_Rate > 0]
    #dataframe distance calculation
    df_alert['Distance'] = haversine_np(df_alert['Cloudsat_Lon'], df_alert['Cloudsat_Lat'], df_alert['Ground_Lon'], df_alert['Ground_Lat'])
    #specifying distance threshold of 50km
    df_alert = df_alert[df_alert.Distance < 100]
    #add file name to keep in track of date
    df_alert['File_Name'] = filename
    #specify the order of columns                
    #df = df[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Bin1", "Bin5", "Bin10", "Bin15", "Bin20", "Bin25", "Bin30", "Bin35", "Bin40", "Bin45", "Bin50", "Bin55", "Bin60", "Bin65", "Bin70", "Bin75", "Bin80", "Bin85", "Bin90", "Bin95", "Bin100", "Bin105", "Bin110", "Bin115", "Bin120", "Bin125", "File_Name"]]
    df_alert = df_alert[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD","Snowfall_Rate_Baselayer", "Elevation_Baselayer", "Bin_Number_Baselayer", "Snowfall_Rate_Midlayer", "Elevation_Midlayer", "Snowfall_Rate_TopLayer", "Elevation_TopLayer", "Bin_Number_TopLayer", "File_Name"]]
    
    
    #Churchill Climate Station
    df_churchill_climate = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std,
                                   "Snowfall_Rate_Baselayer" : snowfall_rate_near_surface_bin, "Elevation_Baselayer" : elevation_near_surface_bin, "Bin_Number_Baselayer" : near_surface_bin, "Snowfall_Rate_Midlayer" : midlayer_snowfall_rate, "Elevation_Midlayer" : midlayer_elevation, "Snowfall_Rate_TopLayer": snowfall_rate_near_surface_endbin, "Elevation_TopLayer" : elevation_near_surface_endbin, "Bin_Number_TopLayer" : near_surface_endbin})
    df_churchill_climate['Ground_Lat'] = 58.73333333
    df_churchill_climate['Ground_Lon'] = -94.06666667
    #remove observations with -ve snowfall rates
    df_churchill_climate = df_churchill_climate[df_churchill_climate.Surface_Snowfall_Rate > 0]
    #dataframe distance calculation
    df_churchill_climate['Distance'] = haversine_np(df_churchill_climate['Cloudsat_Lon'], df_churchill_climate['Cloudsat_Lat'], df_churchill_climate['Ground_Lon'], df_churchill_climate['Ground_Lat'])
    #specifying distance threshold of 50km
    df_churchill_climate = df_churchill_climate[df_churchill_climate.Distance < 100]
    #add file name to keep in track of date
    df_churchill_climate['File_Name'] = filename
    #specify the order of columns                
    #df = df[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Bin1", "Bin5", "Bin10", "Bin15", "Bin20", "Bin25", "Bin30", "Bin35", "Bin40", "Bin45", "Bin50", "Bin55", "Bin60", "Bin65", "Bin70", "Bin75", "Bin80", "Bin85", "Bin90", "Bin95", "Bin100", "Bin105", "Bin110", "Bin115", "Bin120", "Bin125", "File_Name"]]
    df_churchill_climate = df_churchill_climate[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD","Snowfall_Rate_Baselayer", "Elevation_Baselayer", "Bin_Number_Baselayer", "Snowfall_Rate_Midlayer", "Elevation_Midlayer", "Snowfall_Rate_TopLayer", "Elevation_TopLayer", "Bin_Number_TopLayer", "File_Name"]]
    
    
    #Churchill Climate Station
    df_churchill = pd.DataFrame({"Year" : year, "Day_of_Year" : day_of_year,"Cloudsat_Lat" : latitude, "Cloudsat_Lon" : longitude, "UTC_Time" : utc_time,"Surface_Snowfall_Rate" : surface_snowfall_rate, "Surface_Snowfall_Rate_Confidence" : surface_snowfall_rate_confidence, "Surface_Snowfall_Rate_Uncertainty" : surface_snowfall_rate_uncertainty, "DEM_Elevation" : dem_elevation, "Data_Quality" : data_quality, "Data_Status" : data_status, "Data_Target_ID": data_target_id, "Snow_Retrieval_Status" : snow_retrieval_status, "Snowfall_Rate_Mean" : near_surface_snow_layer_avg, "Snowfall_Rate_STD" : near_surface_snowlayer_std,"Overall_Snowfall_Rate_Mean" : snowfall_rate_avg, "Overall_Snowfall_Rate_STD" : snowfall_rate_std,
                                   "Snowfall_Rate_Baselayer" : snowfall_rate_near_surface_bin, "Elevation_Baselayer" : elevation_near_surface_bin, "Bin_Number_Baselayer" : near_surface_bin, "Snowfall_Rate_Midlayer" : midlayer_snowfall_rate, "Elevation_Midlayer" : midlayer_elevation, "Snowfall_Rate_TopLayer": snowfall_rate_near_surface_endbin, "Elevation_TopLayer" : elevation_near_surface_endbin, "Bin_Number_TopLayer" : near_surface_endbin})
    df_churchill['Ground_Lat'] = 58.73916667
    df_churchill['Ground_Lon'] = -94.06638889
    #remove observations with -ve snowfall rates
    df_churchill = df_churchill[df_churchill.Surface_Snowfall_Rate > 0]
    #dataframe distance calculation
    df_churchill['Distance'] = haversine_np(df_churchill['Cloudsat_Lon'], df_churchill['Cloudsat_Lat'], df_churchill['Ground_Lon'], df_churchill['Ground_Lat'])
    #specifying distance threshold of 50km
    df_churchill = df_churchill[df_churchill.Distance < 100]
    #add file name to keep in track of date
    df_churchill['File_Name'] = filename
    #specify the order of columns                
    #df = df[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Bin1", "Bin5", "Bin10", "Bin15", "Bin20", "Bin25", "Bin30", "Bin35", "Bin40", "Bin45", "Bin50", "Bin55", "Bin60", "Bin65", "Bin70", "Bin75", "Bin80", "Bin85", "Bin90", "Bin95", "Bin100", "Bin105", "Bin110", "Bin115", "Bin120", "Bin125", "File_Name"]]
    df_churchill = df_churchill[["Year", "Day_of_Year", "UTC_Time", "Cloudsat_Lat", "Cloudsat_Lon", "Ground_Lat", "Ground_Lon", "Distance", "Surface_Snowfall_Rate", "Surface_Snowfall_Rate_Confidence", "Surface_Snowfall_Rate_Uncertainty", "Data_Quality", "Data_Status", "Data_Target_ID", "Snow_Retrieval_Status","DEM_Elevation","Snowfall_Rate_Mean", "Snowfall_Rate_STD", "Overall_Snowfall_Rate_Mean", "Overall_Snowfall_Rate_STD","Snowfall_Rate_Baselayer", "Elevation_Baselayer", "Bin_Number_Baselayer", "Snowfall_Rate_Midlayer", "Elevation_Midlayer", "Snowfall_Rate_TopLayer", "Elevation_TopLayer", "Bin_Number_TopLayer", "File_Name"]]
    

    #write data frame to csv
    #df.to_csv("/Users/Rithwik/Desktop/p/2006/filepath.csv", index=True)
    #eureka
    if not df_eureka.empty:
        # TODO: include some indication of which superstation was hit here
        print ("Overpass hit_eureka")
        df_eureka.to_csv((outpath + "/eureka/" + filename + ".csv"), index=True)
        
    #eureka-a
    if not df_eureka_a.empty:
        # TODO: include some indication of which superstation was hit here
        print ("Overpass hit_eureka_a")
        df_eureka_a.to_csv((outpath + "/eureka_a/" + filename + ".csv"), index=True)
        
    #resolute
    if not df_resolute.empty:
        # TODO: include some indication of which superstation was hit here
        print ("Overpass hit_resolute")
        df_resolute.to_csv((outpath + "/resolute/" + filename + ".csv"), index=True)
        
    #egbert
    if not df_egbert.empty:
        # TODO: include some indication of which superstation was hit here
        print ("Overpass hit_egbert")
        df_egbert.to_csv((outpath + "/egbert/" + filename + ".csv"), index=True)
    
    #cambridge
    if not df_cambridge.empty:
        # TODO: include some indication of which superstation was hit here
        print ("Overpass hit_cambridge")
        df_cambridge.to_csv((outpath + "/cambridge/" + filename + ".csv"), index=True)
        
    #iqaluit
    if not df_iqaluit.empty:
        # TODO: include some indication of which superstation was hit here
        print ("Overpass hit_iqaluit")
        df_iqaluit.to_csv((outpath + "/iqaluit/" + filename + ".csv"), index=True)
    
    #alert
    if not df_alert.empty:
        # TODO: include some indication of which superstation was hit here
        print ("Overpass hit_alert")
        df_alert.to_csv((outpath + "/alert/" + filename + ".csv"), index=True)
        
    #churchill climate
    if not df_churchill_climate.empty:
        # TODO: include some indication of which superstation was hit here
        print ("Overpass hit_churchill_climate")
        df_churchill_climate.to_csv((outpath + "/churchill_climate/" + filename + ".csv"), index=True)
    
    #churchill
    if not df_churchill.empty:
        # TODO: include some indication of which superstation was hit here
        print ("Overpass hit_churchill")
        df_churchill.to_csv((outpath + "/churchill/" + filename + ".csv"), index=True)
