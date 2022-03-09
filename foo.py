import glob
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
import sys, os, datetime
import math
import calendar
from netCDF4 import Dataset
from datetime import timedelta
from datetime import datetime

run_combiner = False
run_extractor = False
run_extractor_one_degree = False
generate_netCDF = False
generate_netCDF_one_degree = True
gen_blended_accums = False
run_climatology_extractor = True

class BlendedObject(object):
    name = ""
    blended_swe_values = [[]]
    yearly_accumulations = [[]]
    accumulations = []
    days_in_period = 0
    event_count = 0
    total_accumulation = 0
    sd = 0
    se = 0
    CI_upper = 0
    CI_lower = 0
    accumulation_upper = 0
    accumulation_lower = 0

    # The class "constructor" - It's actually an initializer 
    def __init__(self, name, days_in_period):
        self.name = name
        self.blended_swe_values = [[] for x in range(9)]
        self.yearly_accumulations = np.zeros(9)
        self.accumulations = []
        self.days_in_period = days_in_period
        self.event_count = 0
        self.total_accumulation = 0
        self.sd = 0
        self.se = 0
        self.CI_upper = 0
        self.CI_lower = 0
        self.accumulation_upper = 0
        self.accumulation_lower = 0

##### RUN THIS TO EXTRACT SNOWFALL ACCUMULATIONS FOR EACH OF THE CLOUDSAT MASTER GRID CELL FILES
#####
##### ONE DEGREE VARIANT
#####

if run_climatology_extractor:
    # Range of cells to concat
    lat_range = (60, 82)
    lon_range = (-180, 180)

    lats = np.arange(lat_range[0], lat_range[1], 1)
    lons = np.arange(lon_range[0], lon_range[1], 1)
    climatology = list(range(12))
    cols = ['lat', 'lon', 'month', 'ovrps_count', 'sacc', 'me', 'uncert', 'qual', 'conf']
    master_df = pd.DataFrame(columns=cols)
    days_in_period = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    print("Sorting through lats/lons...")
    count = 1
    for lat in lats:
        for lon in lons:
            master_df = pd.DataFrame(columns=cols)
            print("Working on: (" + str(lat) + ", " + str(lon) + ") : " + str(count))
            filepaths = []
            filepaths.append('../master/half_degree_cells/' + str(lat+0.25) + '_' + str(lon+0.25) + '_cell.csv')
            filepaths.append('../master/half_degree_cells/' + str(lat+0.25) + '_' + str(lon+0.75) + '_cell.csv')
            filepaths.append('../master/half_degree_cells/' + str(lat+0.75) + '_' + str(lon+0.25) + '_cell.csv')
            filepaths.append('../master/half_degree_cells/' + str(lat+0.75) + '_' + str(lon+0.75) + '_cell.csv')
            count += 1
                             
            month = []
            for filepath in filepaths:
                if os.path.isfile(filepath):
                    df = pd.read_csv(filepath)
                    for index, row in df.iterrows():
                        if row['year'] == 2016 or row['year'] == 2006 or row['rate'] <= -999:
                            continue
                        month.append((row['utc_start'], row['rate'], row['rate_uncert'], row['confidence'], row['quality']))
                        
            month_count = 0
            
            if True:
                month.sort(key=lambda tup: tup[0])
                prev_start = -1
                overpass_median_rates = []
                overpass_median_uncertainties = []
                overpass_sounding_rates = []
                overpass_uncertainty_rates = []
                overpass_confidences = []
                overpass_qualities = []
                
                for sounding in month:
                    if prev_start == -1:
                        prev_start = sounding[0]

                    overpass_confidences.append(sounding[3])
                    overpass_qualities.append(sounding[4])

                    if prev_start == sounding[0]:
                        overpass_sounding_rates.append(sounding[1])
                        overpass_uncertainty_rates.append(sounding[2])
                    else:
                        overpass_median_rates.append(np.median(overpass_sounding_rates))
                        overpass_median_uncertainties.append(np.median(overpass_uncertainty_rates))
                        overpass_sounding_rates = []
                        overpass_uncertainty_rates = []
                        overpass_sounding_rates.append(sounding[1])
                        overpass_uncertainty_rates.append(sounding[2])
                        prev_start = sounding[0]
                        
                if len(overpass_sounding_rates) > 0:
                    overpass_median_rates.append(np.median(overpass_sounding_rates))
                if len(overpass_uncertainty_rates) > 0:
                    overpass_median_uncertainties.append(np.median(overpass_uncertainty_rates))

                monthly_accumulation = -9999
                monthly_uncertainty = -9999
                monthly_h = -9999
                monthly_confidence = -9999
                monthly_quality = -9999
                hours_in_month = 365 * 24

                if len(overpass_median_rates) > 1:
                    sd = np.std(overpass_median_rates, ddof=1) 
                    se = sd / math.sqrt(len(overpass_median_rates))
                    monthly_h = se * sp.stats.t._ppf((1+0.95)/2., len(overpass_median_rates)-1) * hours_in_month
                    monthly_accumulation = np.mean(overpass_median_rates) * hours_in_month
                    monthly_uncertainty = np.mean(overpass_median_uncertainties) * hours_in_month
                elif len(overpass_median_rates) == 1:
                    monthly_accumulation = np.mean(overpass_median_rates) * hours_in_month
                    monthly_uncertainty = np.mean(overpass_median_uncertainties) * hours_in_month

                if len(overpass_confidences) > 0:
                    monthly_confidence = int(sp.stats.mode(overpass_confidences)[0][0])

                if len(overpass_qualities) > 0:
                    monthly_quality = int(sp.stats.mode(overpass_qualities)[0][0])

                master_df = master_df.append({'lat':lat,
                                              'lon':lon,
                                              'month':month_count,
                                              'ovrps_count':len(overpass_median_rates),
                                              'sacc':monthly_accumulation,
                                              'me':monthly_h,
                                              'uncert':monthly_uncertainty,
                                              'qual':monthly_quality,
                                              'conf':monthly_confidence}, ignore_index=True)
                month_count += 1
#             break
#         break
                master_df.to_csv("../master/one_degree_clim/" + str(lat+0.5) + "_" + str(lon+0.5) + "_stats.csv", index=None, header=['lat', 'lon', 'month', 'ovrps_count', 'sacc', 'me', 'uncert', 'qual', 'conf'])

print("Done")
