import glob
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
from netCDF4 import Dataset
import os
import sys
import math

# Range of cells to concat
lat_range = (72.75, 90.25)
lon_range = (-179.75, 180.25)

lats = np.arange(lat_range[0], lat_range[1], 0.5)
lons = np.arange(lon_range[0], lon_range[1], 0.5)

print("Sorting through lats/lons...")
count = 1
for lat in lats:
    for lon in lons:
        print("Working on: (" + str(lat) + ", " + str(lon) + ") : " + str(count) + " / " + str(43200))
        # Find files across years
        files = []
        years = ["2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015"]

        for year in years:
                filepath = '/users/jk/15/fdmking/data/' + year + '/' + str(lat) + '_' + str(lon) + '_cloudsat_cell.csv'
                if os.path.isfile(filepath):
                        files.append(filepath)

        #for filename in glob.iglob('/users/jk/15/fdmking/data/**/' + str(lat) + '_' + str(lon) + '_cloudsat_cell.csv', recursive=True):
        #        print(filename)
        #        files.append(filename)

        frame = pd.DataFrame()
        # Combine files into single file
        if len(files) > 1:
            list_ = []
            for file_ in files:
                df = pd.read_csv(file_, index_col=None, header=None)
                list_.append(df)
            frame = pd.concat(list_)
        elif len(files) == 1:
            frame = pd.read_csv(files[0], index_col=None, header=None)

        # Save combined output to combined master cell csv
        if not(frame.empty):
            frame.to_csv('/users/jk/15/fdmking/master/' + str(lat) + '_' + str(lon) + '_cell.csv', index=None, header=["lat", "lon", "rate", "rate_uncert", "confidence", "quality", "utc_start", "profile_time", "day", "month", "year"])
        count += 1

print("Done!")
