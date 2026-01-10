import cdsapi
import calendar
import os

startYear = 1990
endYear = 2026

#machineRoot = "/Users/rfield1/data/observations/"
machineRoot = "/autofs/brewer/rfield1/storage/observations/ERA5/"
dataset = "reanalysis-era5-land"


#regionName = "Canada"
#regionBox = [70, -180, 41, -50]
regionName = "Iberia"
regionBox = [44, -10, 35, -4]


outputRoot = machineRoot + dataset + "/"
if not os.path.exists(outputRoot):
   os.mkdir(outputRoot)

#,"2m_temperature","surface_pressure","total_precipitation","surface_solar_radiation_downwards","snow_cover","10m_u_component_of_wind","10m_v_component_of_wind"]


varNameList = ["2m_temperature","2m_dewpoint_temperature","surface_pressure","total_precipitation","surface_solar_radiation_downwards","snow_cover","10m_u_component_of_wind","10m_v_component_of_wind"]

for currVariable in varNameList:


a
    outputDir = outputRoot + regionName + "/" 
if not os.path.exists(outputDir):
   os.mkdir(outputDir)
outputDir = outputRoot + regionName + "/" + "inputData" + "/" 
if not os.path.exists(outputDir):
   os.mkdir(outputDir)
outputDir = outputRoot + regionName + "/" + "inputData" + "/"  + currVariable + "/"
if not os.path.exists(outputDir):
   os.mkdir(outputDir)



for currYear in range(startYear,endYear):
     for currMonth in range(1,13):
        currMonthInfo = calendar.monthrange(currYear,currMonth)
        currNumberOfDaysInMonth = currMonthInfo[1]
        currDateRangeStr = str(currYear) + "-" + str(currMonth).zfill(2) + "-" "01" + "/" + str(currYear) + "-" + str(currMonth).zfill(2) + "-" + str(currNumberOfDaysInMonth).zfill(2)
        request = {
                "variable": currVariable,
                "date": currDateRangeStr,
                "time": ["00:00", "01:00", "02:00","03:00", "04:00", "05:00","06:00", "07:00", "08:00","09:00", "10:00", "11:00","12:00", "13:00", "14:00","15:00", "16:00", "17:00","18:00", "19:00", "20:00","21:00", "22:00", "23:00"],
                "data_format": "netcdf",
                "download_format": "unarchived",
                "area": regionBox
        }
        target = outputDir + str(currYear) + str(currMonth).zfill(2) + "." + currVariable +  ".nc"
        print(target)
        client = cdsapi.Client()
        client.retrieve(dataset, request,target)
		