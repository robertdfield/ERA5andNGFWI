function main()
%{
    This takes a list of points, and for each,
        - reads in the ERA5Land raw time series fetched with the cdsapi
        with a python script
        - does some unit conversions
        - writes out new hourly csv weather file that can be read directly
        by the python scripts from the CFS (which is done in another
        script) to do 'next generation' fire weather index calculations
    By robert.field@columbia.edu
%}
    clear all;

    nProc = 4;
    oldPoolObj = gcp('nocreate');
    if (~isempty(oldPoolObj))
	    delete(oldPoolObj);
    end
    parpool(nProc)
    pctRunOnAll warning('off','all');
    
    % main location to process data
    machineRoot = '/Users/rfield1/data/observations/';

    % name of dataset...not using this?
    dataset = 'reanalysis-era5-land-timeseries';

    outputRoot = [machineRoot '/ERA5Land/reanalysis-era5-land-timeseries/'];
    % regionName = 'WestIberianPeninsula';
    regionName = 'IberianPeninsula';
    
    outputDir = [outputRoot regionName '/'];
    if (~exist(outputDir,'dir'))
        mkdir(outputDir);
    end

    % using the same table of coordinates that the python download script
    % used
    inputRootDir = [machineRoot '/ERA5Land/' dataset '/'];
    coordFN = [outputDir  '/' regionName '.gridSpacingDegrees.0.25.DownloadCoords.csv'];

    coordTable = readtable(coordFN);
    nCoords = size(coordTable,1);

    % the downloaded ERA5Land files might have more data than is wanted
    startDate = datenum(1980,1,1);
    endDate = datenum(2025,11,30);

    parfor currCoordI = 1:nCoords % parforable
        currID = char(coordTable{currCoordI,'ID'});
        stnLat = coordTable{currCoordI,'Lat'}; 	
        stnLon = coordTable{currCoordI,'Lon'};
        currInputDir = [inputRootDir regionName '/' currID '/'];
        stnTimeZone = -timezone(stnLon);

        % name of file in new format
        outputFN = [currInputDir currID '.HourlyWxInput.'  datestr(startDate,'yyyymmdd') '.' datestr(endDate,'yyyymmdd') '.ERA5LAND.csv'];
        zipFileName = [currInputDir currID '.zip'];
        
        if (~exist(outputFN,'file') & exist(zipFileName,'file'))
            try
                % unzip, which know in advance has 5 files
                cmdStr = ['unzip -o '  zipFileName ' -d ' currInputDir];
                system(cmdStr);

                % the variables of interest are bundled in different ways
                % in different files by the cdsapi. All of this and variable naming is
                % from inspection.
                currVarGroup = 'sfc-2m-temperature';
                dirContents = dir([currInputDir '*' currVarGroup '*.csv']);
                currCSVPath = [dirContents.folder '/' dirContents.name];
                currTable = readtable(currCSVPath); 
                
                % arbitrarily using temperature, dew point file for date
                % stuff
                dateVec = datenum(currTable{:,'valid_time'});
                goodDateI = find(dateVec >= startDate & dateVec <= endDate);
                dateVec = dateVec(goodDateI);
                nObs = size(dateVec,1);

                % adjust to local time
                for obI = 1:nObs
                    dateVec(obI,1) = addtodate(dateVec(obI),stnTimeZone,'hour');
                end

                tempTS = currTable{goodDateI,'t2m'};                
                dewpointTS = currTable{goodDateI,'d2m'};
                delete(currCSVPath);

                currVarGroup = 'sfc-pressure-precipitation';
                dirContents = dir([currInputDir '*' currVarGroup '*.csv']);
                currCSVPath = [dirContents.folder '/' dirContents.name];
                currTable = readtable(currCSVPath);
                precTS = currTable{goodDateI,'tp'};
                %surfPresTS = currTable{goodDateI,'sp'};
                delete(currCSVPath);

                currVarGroup = 'sfc-wind';
                dirContents = dir([currInputDir '*' currVarGroup '*.csv']);
                currCSVPath = [dirContents.folder '/' dirContents.name];
                currTable = readtable(currCSVPath);
                uwndTS = currTable{goodDateI,'u10'};
                vwndTS = currTable{goodDateI,'v10'};        
                delete(currCSVPath);

                %
                % make use of this someday. For now processing for the sake
                % of deleting the source csv files
                %
                currVarGroup = 'sfc-snow';
                dirContents = dir([currInputDir '*' currVarGroup '*.csv']);
                currCSVPath = [dirContents.folder '/' dirContents.name];
                currTable = readtable(currCSVPath);
                snowTS = currTable{goodDateI,'snowc'};
                delete(currCSVPath);

                currVarGroup = 'sfc-soil-water';
                dirContents = dir([currInputDir '*' currVarGroup '*.csv']);
                currCSVPath = [dirContents.folder '/' dirContents.name];
                currTable = readtable(currCSVPath);
                soilWater1 = currTable{goodDateI,'swvl1'};
                soilWater2 = currTable{goodDateI,'swvl2'};
                soilWater3 = currTable{goodDateI,'swvl2'};
                soilWater4 = currTable{goodDateI,'swvl2'};
                delete(currCSVPath);

                currVarGroup = 'sfc-radiation-heat';
                dirContents = dir([currInputDir '*' currVarGroup '*.csv']);
                currCSVPath = [dirContents.folder '/' dirContents.name];
                currTable = readtable(currCSVPath);
                downwardSolarTS = currTable{goodDateI,'ssrd'};
                delete(currCSVPath);
                
                tVec = NaN*ones(nObs,1);
                rhSimpleVec = NaN*ones(nObs,1);
                wdSpdVec = NaN*ones(nObs,1);
                precVec = NaN*ones(nObs,1);
            
                % do various conversions on each record
                for currObI = 1:nObs
                    currT = tempTS(currObI);
                    currDewP = dewpointTS(currObI);
                    currUWnd = uwndTS(currObI);
                    currVWnd = vwndTS(currObI);
                    currWdSpd = 3.6*sqrt(currUWnd^2+currVWnd^2); % m/s to kph
                    currPrec = 1000*precTS(currObI);; % m to mm
                    currT = currT -273.15;
                    currDewP = currDewP - 273.15;
                    currRH = RHCalc(currT,currDewP);
                    tVec(currObI,1) = currT;
                    rhSimpleVec(currObI,1) = currRH;
                    wdSpdVec(currObI,1) = currWdSpd;
                    precVec(currObI,1) = currPrec;
                end
                
                % corrections
                rhSimpleVec(find(rhSimpleVec>100))=100;
                precVec(find(precVec<0))=0;
                
                % convert date times for output
                dateVec = datevec(dateVec);
            
                % the format required for CFS NG-FWI scripts
                headerStr = {'id','lat','long','timezone','yr','mon','day','hr','temp','rh','ws','prec'};
                latVec = stnLat*ones(nObs,1);
                lonVec = stnLon*ones(nObs,1);
                timeZoneVec = stnTimeZone*ones(nObs,1);

                % NG-FWI scripts I'm pretty sure require character vector
                % for the 'station id', hence the 'FEDSID_' prefix
                stnIDVec = repmat(cellstr(['FEDSID_' currID]),nObs,1);

                outTable = [latVec lonVec timeZoneVec dateVec(:,1) dateVec(:,2) dateVec(:,3) dateVec(:,4) ...
                            tVec rhSimpleVec wdSpdVec precVec];
                outTable = num2cell(outTable);
                outTable = [stnIDVec outTable];
                outTable = [headerStr; outTable];
                outTable = table(outTable);
                writetable(outTable,outputFN,'WriteVariableNames',0);  
            catch
                disp([currID ' problem...file not there? partial zip file if fetch still in progress?']);
                pause;
            end
        else
            disp([currID ' already processed or zip file not downloaded']) 
        end
    end

    %{
    %
    % This takes the average over the gridpoints, but it's now done
    elsewhere
    %
    takeRegionalAveragesAfterCalculationsHaveBeenDoneInPython = 1;
    fwiVar = 'fwi';
    if (takeRegionalAveragesAfterCalculationsHaveBeenDoneInPython)
        [fireZones,allFireCentres, fireZoneShapes,fireZoneAreas]  = DefineWestIberianPeninsulas(machineRoot);

        nFireZones = size(fireZoneShapes,1);
        for currFireZoneI = 1:nFireZones
            currName = strrep(char(allFireCentres(currFireZoneI )),' ','');

            currPointsI = find(inpolygon(coordTable{:,'Lon'},coordTable{:,'Lat'},fireZoneShapes(currFireZoneI).Lon,fireZoneShapes(currFireZoneI).Lat));
            currNPoints = size(currPointsI,1);
            currFWITableAcrossPoints = [];
            for currPointInPolygonI = 1:currNPoints
                currPointInPolygonI = currPointsI(currPointInPolygonI);
                currID = char(coordTable{currPointInPolygonI,'ID'});

                currInputDir = [inputRootDir regionName '/' currID '/'];
                outputFN = [currInputDir currID '.HourlyFWI.'  datestr(startDate,'yyyymmdd') '.' datestr(endDate,'yyyymmdd') '.ERA5LAND.csv'];
                currFWITable = readtable(outputFN);
               
                currFWITableAcrossPoints = [currFWITableAcrossPoints currFWITable{:,fwiVar} ];
            end
            currFWIAverage = nanmean(currFWITableAcrossPoints,2);

            headerStr = {'id','yr','mon','day','hr','fwi'};            
            currOutputTable = [ currFWITable{:,'yr'} currFWITable{:,'mon'} currFWITable{:,'day'} currFWITable{:,'hr'} currFWIAverage];
            currOutputTable = [repmat(cellstr(currName),size(currFWITable,1),1) num2cell(currOutputTable)];
            currOutputTable = [headerStr;currOutputTable];

            currOutputTable = table(currOutputTable);
            currOutputDir = [inputRootDir regionName '/' currName '/'];
            if (~exist(currOutputDir,'dir'))
                mkdir(currOutputDir);
            end
            outputFN = [currOutputDir currName '.HourlyFWI.'  datestr(startDate,'yyyymmdd') '.' datestr(endDate,'yyyymmdd') '.ERA5LAND.csv'];
            writetable(currOutputTable,outputFN,'WriteVariableNames',0);
        end
    end
    %}
end



function RH = RHCalc(temp, dewp)

    RH = 100*exp(-5423*(temp- dewp)/((temp+ 273.16)^2));

end

%{
    currFN = [ERA5StnDir 'reanalysis-era5-land-timeseries-sfc-radiation-heatetxlwb9_.csv'];
    currFN = [ERA5StnDir 'reanalysis-era5-land-timeseries-sfc-snown399npsf.csv'];
%}

% run 'NG-FWI...might have to pause and do it manually
% pressure in Pa (not hPa nor mbar)
% temperature in K (not in degree Celsius)
% humidity:
% - partial pressure of water vapor in Pa (not hPa nor mbar)
% - specific humidity in kg/kg (not g/kg)
% - mixing ratio in kg/kg (not g/kg)
% - relative humidity in percent
% - dew point temperature in K (not degree Celsius)
% - virtual temperature in K (not degree Celsius)

function newCSVPath = unzipFileInComplicatedWay(currInputDir,currID, currVarName)
% Each file is special in its own way. The zip files are sensibly
% named, but the csv files in them are not really, plus they include a unique API
% request ID. So, it's handled in this complicated way of unzipping to a directory,
% renaming the file sensibly, moving it, then deleting the subdirectory. 
% Its return value is the path to the newly-named csv file
% Probably a better way. 

    newCSVPath = [currInputDir currID '.' currVarName '.csv'];
    currZipFile = [currID '.' currVarName ];
    currOutputDir = [currInputDir currZipFile '/'];

    if (~exist(newCSVPath,'file'))
        if (~exist("currOutputDir"))
            mkdir(currOutputDir)
        end
        cmdStr = ['unzip -o ' currInputDir currZipFile '.zip' ' -d ' currOutputDir];
        system(cmdStr);
    
        dirContents = dir([currOutputDir '*.csv']);
        oldCSVPath = [dirContents.folder '/' dirContents.name];
        cmdStr = ['mv ' oldCSVPath ' ' newCSVPath ];
        system(cmdStr);
    end
    if (exist(currOutputDir,'dir'))
        rmdir(currOutputDir);
    end
end



function [fireZones,allFireCentres, fireZoneShapes,fireZoneAreas]  =  DefineWestIberianPeninsulas(machineRoot)
% from Sanchez-Hernandez et al. 2025 Global Change Bio.
%
% use custom polygons made in QGIS
%
    % low res
    latRange = [33 44];
    lonRange = [-10 0];
    nGoodShapes = 0;

    nGoodShapes = nGoodShapes+1;
    currFN = [machineRoot '/shapefiles/NewIberia/NW-IP.shp'];    
    currShape = shaperead(currFN,"UseGeoCoords",true);
    currShapeInfo = readgeotable(currFN);
    currShapeInfo = currShapeInfo.Shape;
    fireZones(nGoodShapes) = cellstr('');
    allFireCentres(nGoodShapes) = cellstr('Northwest Iberian Peninsula'); %cellstr('Spain-CastLeonNorthwest');
    currLats = currShape.Lat';
    currLons = currShape.Lon';
    fireZoneShapes(nGoodShapes,1).Lon = currLons;
    fireZoneShapes(nGoodShapes,1).Lat = currLats;
    currArea = areaint(currLats,currLons,currShapeInfo.GeographicCRS.Spheroid); %m^2 ?
    fireZoneAreas(nGoodShapes,1) =  nansum(currArea);  
        

    nGoodShapes = nGoodShapes+1;
    currFN = [machineRoot '/shapefiles/NewIberia/SW-IP.shp'];    
    currShape = shaperead(currFN,"UseGeoCoords",true);
    currShapeInfo = readgeotable(currFN);
    currShapeInfo = currShapeInfo.Shape;
    fireZones(nGoodShapes) = cellstr('');
    allFireCentres(nGoodShapes) = cellstr('Southwest Iberian Peninsula'); %cellstr('Spain-CastLeonNorthwest');
    currLats = currShape.Lat';
    currLons = currShape.Lon';
    fireZoneShapes(nGoodShapes,1).Lon = currLons;
    fireZoneShapes(nGoodShapes,1).Lat = currLats;
    currArea = areaint(currLats,currLons,currShapeInfo.GeographicCRS.Spheroid); %m^2 ?
    fireZoneAreas(nGoodShapes,1) =  nansum(currArea);  
        

    goodI = [1:nGoodShapes];    
    fireZones = fireZones(goodI)';
    allFireCentres = allFireCentres(goodI)';
    fireZoneShapes = fireZoneShapes(goodI);   
    fireZoneAreas = fireZoneAreas/(1000*1000); % m2 to km2

end
