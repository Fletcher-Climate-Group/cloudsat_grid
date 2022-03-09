import sys, os, datetime, json, math
import pandas as pd
import fnmatch
import calendar
import numpy as np
import load_hdf as lh

if len(sys.argv) < 2:
    print("ERROR: Please include a year to extract from...")
    sys.exit()

# Read filename
glob_year = sys.argv[1]

class Sounding(object):
    lat = 0
    lon = 0
    snowfall_rate = 0
    snowfall_uncert = 0
    snowfall_confidence = 0
    data_quality = 0
    utc_start = 0
    profile_time = 0
    day = 0
    month = 0
    year = 0

    # The class "constructor" - It's actually an initializer 
    def __init__(self, lat, lon, snowfall_rate, snowfall_uncert, snowfall_confidence, data_quality, utc_start, profile_time, day, month, year):
        self.lat = lat
        self.lon = lon
        self.snowfall_rate = snowfall_rate
        self.snowfall_uncert = snowfall_uncert
        self.snowfall_confidence = snowfall_confidence
        self.data_quality = data_quality
        self.utc_start = utc_start
        self.profile_time = profile_time
        self.day = day
        self.month = month
        self.year = year

#### HELPERS
def convert_DOY(y,jd,type):
    month = 1
    day = 0
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month = month + 1
    if type == "day":
        return jd
    return month

#### GLOBAL VARS
data_root = '../cloudsat/2C-SNOW-PROFILE.P1_R05/' + glob_year
min_lat = 60
max_lon = 180
prev_doy = -1
modified_lat_lon = []
day_grid = [[[] for j in range(721)] for i in range(60)]

# Walk through all the HDF data
print("Beginning extraction...")
for dir_, _, files in os.walk(data_root):
    for file in files:
        if fnmatch.fnmatch(file, '*.hdf'):
            rel_dir = os.path.relpath(dir_, data_root)
            rel_file = os.path.join(rel_dir, file)
            full_path = os.path.join(data_root, rel_file)

            vd_vars = ['Latitude', 'Longitude', 'Profile_time', 'UTC_start', 'snowfall_rate_sfc', 'snowfall_rate_sfc_uncert', 'snowfall_rate_sfc_confidence', 'Data_quality']
            vd = lh.load_vd(full_path, vd_vars)

            lats = np.transpose(vd[0])[0]
            lons = np.transpose(vd[1])[0]
            profile_times = np.transpose(np.array(vd[2]))[0]
            utc_start = np.transpose(np.array(vd[3]))[0][0]
            snowfall_rates_sfc = np.transpose(np.array(vd[4]))[0]
            snowfall_rates_sfc_uncert = np.transpose(np.array(vd[5]))[0]
            snowfall_rates_sfc_confidence = np.transpose(np.array(vd[6]))[0]
            data_quality = np.transpose(np.array(vd[7]))[0]

            year = int(file[:4])
            doy = int(file[4:7])
            month = convert_DOY(year, doy, "month")
            day = convert_DOY(year, doy, "day")

            print("\nReading", rel_file)
            for i,lat in enumerate(lats):
                if lat > min_lat:
                    grid_lat_pos = math.floor((lat - min_lat)*2)
                    grid_lon_pos = math.floor((lons[i] + max_lon)*2)
                    modified_lat_lon.append((grid_lat_pos, grid_lon_pos))
                    day_grid[grid_lat_pos][grid_lon_pos].append(Sounding(lat, lons[i], snowfall_rates_sfc[i], snowfall_rates_sfc_uncert[i], snowfall_rates_sfc_confidence[i], data_quality[i], utc_start, profile_times[i], day, month, year))

print("Saving results...")
modified_set = list(set(modified_lat_lon))
for i,grid_cell in enumerate(modified_set):
    soundings = day_grid[modified_set[i][0]][modified_set[i][1]]
    out_file = str((modified_set[i][0] / 2) + 60.25) + "_" + str((modified_set[i][1]) / 2 - 179.75) + "_cloudsat_cell.csv"
    cell_df = pd.DataFrame()

    print("Getting sounding data...")
    for j,sounding in enumerate(soundings):
        df = pd.DataFrame({ "Lat": sounding.lat, "Lon": sounding.lon, "Rate": sounding.snowfall_rate, "Uncert": sounding.snowfall_uncert, "Confidence": sounding.snowfall_confidence, "Data_quality": sounding.data_quality, "UTC_start": sounding.utc_start, "Profile_time": sounding.profile_time, "Day": sounding.day, "Month": sounding.month, "Year": sounding.year}, index=[0])
        cell_df = cell_df.append(df)

    print("Saving", i, len(modified_set))
    cell_df.to_csv(os.path.join(glob_year, out_file), mode='a', index=False, header=False, columns=["Lat", "Lon", "Rate", "Uncert", "Confidence", "Data_quality", "UTC_start", "Profile_time", "Day", "Month", "Year"])

print("Done!")

