#
# This fetches ERA5Land single point time series for coordinates in a csv files for a list of variables.

# More info about the data here: 
# https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries?tab=overview
#
# There are also instructions for API downloading using the 'cdsapi' (Climate Data Store API)
# https://cds.climate.copernicus.eu/how-to-api
#
# Do a manual data selection and fetch to get the underlying data fetch python code. 
#
# This script was built by putting for loops around that kind of code
#
# By robert.field@columbia.edu
#

import os
import csv
import cdsapi
import logging
import sys


# main local computer location. Change this according to where you'd like your output on your machine
machineRoot = "/Users/rfield1/data/observations/"

# name of dataset
dataset = "reanalysis-era5-land-timeseries"

# csv file with list of centroids to fetch. This could be anything really, with an ID, latitude and longitude variable specified in the loop below
coordFN = machineRoot + "FEDS/LOCAL/NRT_Europe_West_Siberia_20201125_PM_lf_perimeter.fgb_20250401_20251231.nRegions.2.csv"

# this would be better if the string was constructed from two date variables
dateRangeString = "1950-01-02/2025-11-30"

# list of variables to fetch data for. These ones are for NG-FWI input
varNameList = ["2m_temperature","2m_dewpoint_temperature","surface_pressure","total_precipitation","surface_solar_radiation_downwards","snow_cover","10m_u_component_of_wind","10m_v_component_of_wind"]


#
# All definitions and settings etc above here please
#

# specific location to put output. There will be subdirectories made in this, one for each centroid location
outputRoot = machineRoot + "/ERA5Land/" + dataset + "/"

# By default the API displays a LOT of messages
# Tried to log INFO messages rather than display them on each call
# This did not work. Try again sometime...
#logging.basicConfig(filename = outputRoot + "logfile.txt")
#logger = logging.getLogger("cdsapi")
# only log really bad events
#logger.setLevel(logging.ERROR)

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
        currID = row['fireID']
        currFireMaxArea = row['fireMaxArea']
        currLatCentroid = row['lat_centroid']
        currLonCentroid = row['lon_centroid']
        
        # print out info for current centroid
        print(currID,currFireMaxArea,currLatCentroid,currLonCentroid)
        
        # make location to store output for 'currID'. There could be a lot of coordinates in coordinateTable, so it's nice to split the files up like this
        outputDir = outputRoot + "/" + currID + "/"
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)
        
        # download each variable as its own file. 
        # This is probably not efficient data-fetch wise. 
        # But I think the cdsapi outputs hard-to-understand groups of files named according to their internal variable organization, not according to the variable name specified
        # The csv files inside the 'target' zip files are oddly named FYI, like according to the internal variable organization. But the data column inside the csv has the correct name, which is their acronym for the names in the varNameList
        # So for now, this is easier
        # That date string should not be hard-coded
        # Do a manual data fetch and look a the API code it makes to understand where this comes from
        #
        # The returned file has lat/lon for each row. This is very redundant and the data are elswehere. It probably doubles the file size. It would be nice to omit this.
        #
        
        
        for currVar in varNameList:
            request = {
            "variable": currVar,
            "location": {"longitude": str(currLonCentroid), "latitude": str(currLatCentroid)},
            "date": [dateRangeString],
            "data_format": "csv"
            }
            
            # output file name
            target = outputDir + currID + "." + currVar + ".zip"
            print(target)
            
            # fetch the data if it doesn't already exist
            if not os.path.exists(target):
                client = cdsapi.Client()
                client.retrieve(dataset, request,target)
                
                # check that it was downloaded. For now, want to see the reason, maybe dea
                if not os.path.exists(target):
                    print(target + " not properly downloaded")
                    sys.exit()
                
            else:
                print("already exists")
	
	

	


