#
# This does multiple calls to hFWI function from CFS over a list of points that time series have been downloaded
# and pre-processed for. It runs for ERA5 Landhourly time series, but could be modified for other weather data.
#
# This script was built using pieces of the NG-FWI tutorial
#
# By robert.field@columbia.edu
#


### Load packages ###
# Run `pip install` to install any you are missing.
import os
import csv
import logging
import sys

import pandas as pd
from datetime import datetime

### Load functions and data ###
# If the working directory is different from where you saved the FWI2025 scripts,
# you can add the path to the scripts with the sys package and sys.path.append().
import sys
#sys.path.append("CHANGE/PATH/TO/cffdrs-ng/FWI/Python")
sys.path.append("/Users/rfield1/My Drive/projects/EIS-Fire/EuropeFires2025/cffdrs-ng-main/FWI/Python/")

# Load the files containing the variables and functions to calculate FWI2025.
from NG_FWI import hFWI
from daily_summaries import generate_daily_summaries

# main local computer location. Change this according to where you'd like your output on your machine
machineRoot = "/Users/rfield1/data/observations/"

# name of dataset
dataset = "reanalysis-era5-land-timeseries"

# csv file with list of centroids to fetch. This could be anything really, with an ID, latitude and longitude variable specified in the loop below
# coordFN = machineRoot + "FEDS/LOCAL/NRT_Europe_West_Siberia_20201125_PM_lf_perimeter.fgb_20250401_20251231.nRegions.2.csv"

#regionName = "WestIberianPeninsula"
#fileDateString = "19900101.20251130"

regionName = "IberianPeninsula";
fileDateString = "20250101.20251130";

outputRoot = machineRoot + "ERA5Land/" + dataset + "/" + regionName + "/"

coordFN = outputRoot + regionName + ".gridSpacingDegrees.0.25.DownloadCoords.csv"

#
# All definitions and settings etc above here please
#

# specific location to put output. There will be subdirectories made in this, one for each centroid location
outputRoot = machineRoot + "/ERA5Land/" + dataset + "/" + regionName + "/"

# make the output directories if they doesn't exist. Maybe this can be done once, but I seem to have to do each directory and subdirectory
if not os.path.exists(machineRoot):
   os.mkdir(machineRoot)
if not os.path.exists(outputRoot):
   os.mkdir(outputRoot)
   
# open the csv and read it
with open(coordFN, newline='') as csvfile:
    coordinateTable = csv.DictReader(csvfile)
    for row in coordinateTable:
        # get current coordinates and an ID to name output files with. the currFireMaxArea is not being used, but is nice to see
        # variable names could be modified. It would be a good idea someday to replace these hard-coded variable names with variables whose values are set above 
 #       currID = row['fireID']
 #       currFireMaxArea = row['fireMaxArea']
 #       currLatCentroid = row['lat_centroid']
 #       currLonCentroid = row['lon_centroid']
       	currID = row['ID']
        currFireMaxArea = row['Area']
        currLatCentroid = row['Lat']
        currLonCentroid = row['Lon']  
  
        
        # print out info for current centroid
        print(currID,currLatCentroid,currLonCentroid)
        
        # make location to store output for 'currID'. There could be a lot of coordinates in coordinateTable, so it's nice to split the files up like this
        outputFN =machineRoot + "/ERA5Land/" + dataset + "/"  + regionName + "/" + currID + "/" + currID + "." + "HourlyFWI" + "."  + fileDateString + ".ERA5LAND.csv"
        inputFN = machineRoot + "/ERA5Land/" + dataset + "/" + regionName + "/" + currID + "/" + currID + "." + "HourlyWxInput" + "."  + fileDateString + ".ERA5LAND.csv"
        if not os.path.exists(outputFN):
             if os.path.exists(inputFN):
                 inputData = pd.read_csv(inputFN)
                 try:                 
                     data_fwi = hFWI(inputData)
                     data_fwi.to_csv(outputFN, index = False)
                 except ValueError:
                    print("Problem calculating " + inputFN + ". Sometimes it's because point was over water and there's no data. There are better ways of handling this than throwing an exception")
             else:
                 print(inputFN + " not found" )
        else:
             print(outputFN + " already exists")
