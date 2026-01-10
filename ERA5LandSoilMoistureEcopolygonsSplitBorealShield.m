
function main()
    clear all;


% data documentation
% https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation


% single-point time series download
%https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries?tab=overview


    nProc = 6;
    oldPoolObj = gcp('nocreate');
    if (~isempty(oldPoolObj))
        delete(oldPoolObj);
    end
    parpool(nProc)
    pctRunOnAll warning('off','all');

    machineDir = ['/Users/rfield1/data/'];
    rootDir = [machineDir '/observations/ERA5Land/'];
    regionName = 'Canada'; % descriptive name of downloaded data

    
    varLongName = 'volumetric_soil_water_layer';
    varPrefix = 'swvl';
    varAddAmount = 0;
    nLayers = 4;
    yLabel = 'Volumetric soil water (m^3 m^{-3})';
    yLim = [0 .8];%nanmax(totalSoilWater)];
    varLayerBottomDepths =[0,7,28,100,289]'; % their edges, in cm %swvl1 (0-7 cm)	swvl2 (7-28 cm)	swvl3 (28-100 cm)	swvl4 (100-289 cm)
    varLayersToActuallyTotal = 2; % number of layers (rather than edges) to include in downward total

    %{
    varLongName = '2m_temperature';
    varPrefix = 't2m';
    varAddAmount = -273.15; % K to C
    yLim = [-40 31];%nanmax(totalSoilWater)];
    nLayers = 1;
    yLabel = '2m temperature (C)';
    varLayerBottomDepths =[0 1]'; % single layer for unit-weighting
    varLayersToActuallyTotal = 1; % number of layers (rather than edges) to include in downward total
    %}

    inputDir = [rootDir regionName  '/inputData/' varLongName '/'];
    outputDir = [rootDir regionName  '/outputData/' varLongName '/'];
    mkdir(outputDir);

    startYear = 2002;
    endYear = 2025;
    colorStartYear= 2022;

    overlayStartMonth = 1;
    dateIncrement = NaN;
    dateIncrementUnit = NaN;
    plotSmoothFactor=1;
    fontSize = 30;                
    
    % no longer used, remove at some point
    lockDownStartMonth = NaN;
    lockDownStartDay = NaN;
    lockDownEndMonth = NaN;
    lockDownEndDay = NaN;

    doCumulative = 0;
    nCols = 2;

    %
    % 
    %
    % should rework to convert any  specific polygon into a generic
    % polygon, including rectangular boxes
    %
    currRegionShapes =  shaperead([[machineDir '/observations/shapefiles/ecozones/' 'ecozones.shp']],'UseGeoCoords',true);
    nRegionShapes = size(currRegionShapes,1);

    goodShapesI = [1:nRegionShapes];
    maskShapesI = 0*[1:nRegionShapes];
    for currShapeI = 1:nRegionShapes
        currShape = currRegionShapes(currShapeI);
        if (strcmp(currShape.ZONE_NAME,'Arctic Cordillera') | strcmp(currShape.ZONE_NAME,'Northern Arctic') | strcmp(currShape.ZONE_NAME,'Southern Arctic') | ...
            strcmp(currShape.ZONE_NAME,'Prairie') | strcmp(currShape.ZONE_NAME,'MixedWood Plain')  | strcmp(currShape.ZONE_NAME,'Pacific Maritime') | strcmp(currShape.ZONE_NAME,'Atlantic Maritime') | ...
            (currShape.ZONE_==23) | ...
            strcmp(currShape.ZONE_NAME,'Boreal Shield')); % excluding this big one, but adding on smaller ecoprovince-based ones that comprise it, from an adhoc selection
            goodShapesI(currShapeI) = 0;
            maskShapesI(currShapeI) = 1;
        end
    end

    currMaskShapes = currRegionShapes(find(maskShapesI));
    nMaskShapes = size(currMaskShapes,1);

    currRegionShapes = currRegionShapes(find(goodShapesI));
    nRegionShapes = size(currRegionShapes,1);

    % now add the extra QGIS-made ones for Boreal Shield
    nRegionShapes = nRegionShapes + 1;
    subregionName = 'BorealShieldWest';
    currExtraShape =  shaperead([machineDir '/observations/shapefiles/CanadaCustomEcoPolygons/' subregionName 'Dissolved.shp'],'UseGeoCoords',true);
    currRegionShapes(nRegionShapes).ZONE_ = '17';
    currRegionShapes(nRegionShapes).Geometry = currExtraShape.Geometry;
    currRegionShapes(nRegionShapes).Lat = currExtraShape.Lat;
    currRegionShapes(nRegionShapes).Lon = currExtraShape.Lon;
    currRegionShapes(nRegionShapes).BoundingBox = currExtraShape.BoundingBox;
    currRegionShapes(nRegionShapes).ZONE_NAME = subregionName;

    nRegionShapes = nRegionShapes + 1;
    subregionName = 'BorealShieldCentral';
    currExtraShape =  shaperead([machineDir '/observations/shapefiles/CanadaCustomEcoPolygons/' subregionName 'Dissolved.shp'],'UseGeoCoords',true);
    currRegionShapes(nRegionShapes).ZONE_ = '17';
    currRegionShapes(nRegionShapes).Geometry = currExtraShape.Geometry;
    currRegionShapes(nRegionShapes).Lat = currExtraShape.Lat;
    currRegionShapes(nRegionShapes).Lon = currExtraShape.Lon;
    currRegionShapes(nRegionShapes).BoundingBox = currExtraShape.BoundingBox;
    currRegionShapes(nRegionShapes).ZONE_NAME = subregionName;

    nRegionShapes = nRegionShapes + 1;
    subregionName = 'BorealShieldEast';
    currExtraShape =  shaperead([machineDir '/observations/shapefiles/CanadaCustomEcoPolygons/' subregionName 'Dissolved.shp'],'UseGeoCoords',true);
    currRegionShapes(nRegionShapes).ZONE_ = '17';
    currRegionShapes(nRegionShapes).Geometry = currExtraShape.Geometry;
    currRegionShapes(nRegionShapes).Lat = currExtraShape.Lat;
    currRegionShapes(nRegionShapes).Lon = currExtraShape.Lon;
    currRegionShapes(nRegionShapes).BoundingBox = currExtraShape.BoundingBox;
    currRegionShapes(nRegionShapes).ZONE_NAME = subregionName;

    nRegionShapes = nRegionShapes + 1;
    subregionName = 'BorealShieldAtlantic';
    currExtraShape =  shaperead([machineDir '/observations/shapefiles/CanadaCustomEcoPolygons/' subregionName 'Dissolved.shp'],'UseGeoCoords',true);
    currRegionShapes(nRegionShapes).ZONE_ = '17';
    currRegionShapes(nRegionShapes).Geometry = currExtraShape.Geometry;
    currRegionShapes(nRegionShapes).Lat = currExtraShape.Lat;
    currRegionShapes(nRegionShapes).Lon = currExtraShape.Lon;
    currRegionShapes(nRegionShapes).BoundingBox = currExtraShape.BoundingBox;
    currRegionShapes(nRegionShapes).ZONE_NAME = subregionName;
    
%{
    At this this point need structure with

    Mandatory
    - ZONE_NAME
    - ZONE_
    - Lat
    - Lon

    As of yet unused
    - Geometry ("Polygon")
    - Bounding Box
%}
    % read first file for lat/lon only
    currPath = [inputDir num2str(startYear) num2str(1,'%.02d') '.nc'];
    currLat = ncread(currPath,'latitude');
    currLon = ncread(currPath,'longitude');
    [YI,XI] = meshgrid(currLat,currLon);

    % assign wx grid cells to zone
    parfor currZoneI = 1:nRegionShapes        % parforable
        currPoly = [currRegionShapes(currZoneI).Lon' currRegionShapes(currZoneI).Lat'] ;        
        [inI] = find(inpolygon(XI(:),YI(:),currPoly(:,1),currPoly(:,2)));
        fireZoneIndexList(currZoneI).IndexList = inI;
    end
    fireZoneWxArray = NaN*ones(nRegionShapes,size(XI,1),size(XI,2));
    for currZoneI = 1:nRegionShapes
        inI =  fireZoneIndexList(currZoneI).IndexList;
        tempArray = 0*XI;
        tempArray(inI) = 1;
        fireZoneWxArray(currZoneI,:,:) = tempArray;
    end  


    %
    % really need to 
    % - parfor this.
    % - make generic to any dataset, with a bit of different logic for 2D
    % vs. 3D variables
    %
    dailyDates = [];
    dailyTimeSeriesZones = [];
    timerH = tic;
    for currYear = startYear:endYear
        for currMonth = 1:12
            currPath = [inputDir num2str(currYear) num2str(currMonth,'%.02d') '.nc']
            totalSoilWaterForAMonthZones = [];
            if (exist(currPath,'file'))
                currLat = ncread(currPath,'latitude');
                currLon = ncread(currPath,'longitude');
                currTime = ncread(currPath,'valid_time');
                nDays = size(currTime,1);

                era5SoilMoistureDatesForAMonth = [];
                dailyTimeSeriesForAMonth= [];
                dailyTimeSeriesForAMonthZones=[];

                %
                % daily spatial averages for each layer
                %
                for currLayerI = 1:nLayers
                    if (nLayers>1)
                        currGrid = ncread(currPath,[varPrefix num2str(currLayerI)]);
                    else
                        currGrid = ncread(currPath,varPrefix);
                    end
                
                    currGrid = currGrid  + varAddAmount;
                    for currDayI = 1:nDays
                        
                        era5SoilMoistureDatesForAMonth(currDayI,1) = datenum(currYear,currMonth,currDayI);

                        for currZone = 1:nRegionShapes
                            % find all grid cells for current region
                            % extract average over region
                            currFireZoneWxArray = squeeze(fireZoneWxArray(currZone,:,:));
                            currGridI = find(currFireZoneWxArray == 1);
                            if (~isempty(currGridI))
                                currSmallGrid = squeeze(currGrid(:,:,currDayI));
                                currSmallGrid = squeeze(currSmallGrid(currGridI));
                                currWeighting = squeeze(cos(deg2rad(YI(currGridI))));
                                goodI = find(~isnan(currSmallGrid));
                                if (~isempty(goodI))
                                    currSmallGrid = wmean(currSmallGrid(goodI),currWeighting(goodI));
                                    dailyTimeSeriesForAMonthZones(currZone,currDayI,currLayerI) = currSmallGrid;
                                end
                            else
                                dailyTimeSeriesForAMonthZones(currZone,currDayI,currLayerI) = NaN;
                            end
                        end
                    end
                end

                %
                % total soil moisture over layers
                %
                for currZone = 1:nRegionShapes
                    totalSoilWater = [];
                    for currvarLayerI = 1:varLayersToActuallyTotal
                        currWeightedSoilWater = squeeze(dailyTimeSeriesForAMonthZones(currZone,:,currvarLayerI));
                        currWeighting = (varLayerBottomDepths(currvarLayerI+1) - varLayerBottomDepths(currvarLayerI))/varLayerBottomDepths(varLayersToActuallyTotal+1);
                        totalSoilWater(currvarLayerI,:) = currWeighting*currWeightedSoilWater;
                    end
                    totalSoilWaterForAMonthZones(currZone,:) = nansum(totalSoilWater,1)';
                end

                dailyDates = [dailyDates era5SoilMoistureDatesForAMonth'];
                dailyTimeSeriesZones = [dailyTimeSeriesZones totalSoilWaterForAMonthZones];
            end
        end      
    end
    timerH = toc;
    disp ([num2str(timerH) ' s to process ERA5']);

    close all;

    nCols = 3;
    nRows = ceil(nRegionShapes/nCols);
    scaleFactor = 2;

    figH = figure('Visible','off');
    figPos = get(figH,'Position');
    set(figH,"Position",[figPos(1) figPos(2)  nCols*scaleFactor*figPos(3) nRows*scaleFactor*figPos(4)])    
    tileH = tiledlayout(nRows,nCols,'TileSpacing','compact','Padding','compact');
    clear columnLabels;

    for currSiteI = 1:nRegionShapes
        siteName = [char(currRegionShapes(currSiteI).ZONE_NAME) ' (ZONE ' num2str(currRegionShapes(currSiteI).ZONE_) ')' ]; 
        siteName = strrep(siteName,'_','');
        columnLabels(currSiteI) = cellstr([char(currRegionShapes(currSiteI).ZONE_NAME) '_' 'ZONE' num2str(currRegionShapes(currSiteI).ZONE_) ]);

        era5SoilMoisture = squeeze(dailyTimeSeriesZones(currSiteI,:,:));
        nexttile;
        titleStr = [cellstr([siteName])];% ' (' siteCoordString ')'])];
        PlotAnnualOverLay(dailyDates',era5SoilMoisture',titleStr,yLabel,startYear,endYear,colorStartYear,overlayStartMonth,dateIncrement,dateIncrementUnit,plotSmoothFactor,yLim,fontSize,...
                                figH,lockDownStartMonth,lockDownStartDay,lockDownEndMonth,lockDownEndDay,doCumulative);
    end

    if (varLayersToActuallyTotal>1)
        title(tileH,(['ERA5-Land, ' yLabel ' 0-' num2str(varLayerBottomDepths(varLayersToActuallyTotal+1)) ' cm, ' num2str(startYear) ' to ' num2str(endYear)]),'FontSize',fontSize,'FontWeight','bold');
    else
        title(tileH,['ERA5-Land, ' yLabel ', ' num2str(startYear) ' to ' num2str(endYear)],'FontSize',fontSize,'FontWeight','bold');        
    end

    outFN = [outputDir 'DailyERA5Land.' varLongName 'Overlay' '.nRegionShapes.' num2str(nRegionShapes) '.png'];
    saveas(gcf,outFN);

    dailyTimeSeriesZones = dailyTimeSeriesZones';
    outputDateVec = datevec(dailyDates);
    outputDateVec = outputDateVec(:,1:3);
    columnLabelsTable = [cellstr('YYYY') cellstr('MM') cellstr('DD')  columnLabels];

    outputTable = [outputDateVec dailyTimeSeriesZones];
    outputTable = num2cell(outputTable);
    outputTable = table([columnLabelsTable; outputTable]);
    outFN = [outputDir 'DailyERA5Land.' varLongName  '.nRegionShapes.' num2str(nRegionShapes)  '.varLayersToActuallyTotal.' num2str(varLayersToActuallyTotal) '.csv'];
    writetable(outputTable,outFN,'WriteVariableNames',false);

    %
    % calculate annual seasonal soil moisture...can't yet span calendar
    % years
    %
    seasonStartMonth = 5;
    seasonEndMonth = 9;

    clear annualMeans;
    for currYear = startYear:endYear
        currSeasonI = find( (currYear == outputDateVec(:,1) & outputDateVec(:,2)>=seasonStartMonth & outputDateVec(:,2)<=seasonEndMonth));
        currYearI = (currYear-startYear+1);
        annualMeans(currYearI,:) = nanmean(dailyTimeSeriesZones(currSeasonI,:));
    end

    plot([startYear:endYear],annualMeans,'LineWidth',2);
    legend(columnLabelsTable(2:end),'Location','EastOutside');

    columnLabelsTable = [cellstr('YYYY')  columnLabels];
    outputTable = [[startYear:endYear]' annualMeans];
    outputTable = num2cell(outputTable);
    outputTable = table([columnLabelsTable; outputTable]);
    outFN = [outputDir 'SeasonalERA5.'  varLongName  '.nRegionShapes.' num2str(nRegionShapes)  '.varLayersToActuallyTotal.' num2str(varLayersToActuallyTotal) ...
            '.seasonStartMonth.' num2str(seasonStartMonth) '.seasonStartMonth.' num2str(seasonEndMonth)  '.csv'];
    writetable(outputTable,outFN,'WriteVariableNames',false);

end


function PlotAnnualOverLay(allTimeStepDates,allTimeStepAvg,titleStr,yLabel,startYear,endYear,colorStartYear,overlayStartMonth,dateIncrement,dateIncrementUnit,plotSmoothFactor,yLim,fontSize,...
                            figH,lockDownStartMonth,lockDownStartDay,lockDownEndMonth,lockDownEndDay,doCumulative)


    dateTickFormat = 'mmm';
	%set(0, 'currentfigure', figH);
    currSubAx =gca;
    hold on;
    if (overlayStartMonth>1)
        plotYearList = [startYear:endYear-1];
    else
        plotYearList = [startYear:endYear];
    end
  
	nPlotYears = size(plotYearList,2);

    for currYearI = 1:nPlotYears
        if (mod(currYearI,2)==0)
            allLineStyle(currYearI) = cellstr('-');
            allLineWidth(currYearI) = 3;
        else
            allLineStyle(currYearI) = cellstr('-'); % dashed -- wasn't showing up 
            allLineWidth(currYearI) = 3;
        end

        % but thin grey line before colorStartYear
        if (plotYearList(currYearI)<=colorStartYear)
            allLineStyle(currYearI) = cellstr('-');
            allLineWidth(currYearI) = .5;
        end
    end
    
    annualLineColors = repmat([.5 .5 .5],nPlotYears, 1);
    nColorYears = endYear - colorStartYear+1;
    %currCMap = buildcmap('br');
    %currCMap = currCMap(1:3:end,:); 
    %currCMap = currCMap(1:nColorYears,:);
    annualLineColors(nPlotYears-nColorYears+1:end,:) = hsv(nColorYears); %
    %annualLineColors(nPlotYears,:) = [0 0 0];
    
    for currYearI = 1:nPlotYears
        currYear = plotYearList(currYearI);

        clear fakeDateNums;

        %{

        fakeDateStart = datenum(startYear,overlayStartMonth,1);
        fakeDateEnd = addtodate(fakeDateStart,12,'month');
        fakeDateEnd = addtodate(fakeDateEnd,-1,'day');

% this assumed evenly-spaced dates
        currFakeDateNum = fakeDateStart;
        currFakeStepI = 1;
        while(currFakeDateNum < fakeDateEnd)
            fakeDateNums(currFakeStepI) = currFakeDateNum;

            currFakeStepI = currFakeStepI + 1;
            currFakeDateNum = addtodate(currFakeDateNum,dateIncrement,dateIncrementUnit);
        end
        %}

        %fakeDateVec = datevec(fakeDateNums);
        %currI = find(~(fakeDateVec(:,2) == 2 & fakeDateVec(:,3) == 29));
        %fakeDateNums = fakeDateNums(currI);

        currStartDate = datenum(currYear,overlayStartMonth,1);
        currEndDate = addtodate(currStartDate,1,'year');

        currI = find(allTimeStepDates >= currStartDate & allTimeStepDates < currEndDate );
        currPlotDateNums = allTimeStepDates(currI);
        currYearNDates = size(currPlotDateNums,1);
        for currDateI = 1:currYearNDates
            currDateVec = datevec(currPlotDateNums(currDateI));
            fakeDateNums(currDateI,1)= datenum(startYear,currDateVec(2),currDateVec(3),currDateVec(4),currDateVec(5),currDateVec(6));
        end


        if (doCumulative)
            currPlotTS = cumsum(allTimeStepAvg(currI));
        else
            currPlotTS = allTimeStepAvg(currI);
        end

        %currNSteps = size(currPlotTS,1);
        %fakeDateNums = fakeDateNums(1:currNSteps);    
        % smoothed plot
        smoothDailyPlot = smooth(currPlotTS,plotSmoothFactor,'moving');

        goodDayI = find(~isnan(smoothDailyPlot));

        plotDates = fakeDateNums; % + (fakeDateNums(2)-fakeDateNums(1))/2; can't remember why this bump was here
        endI = min(size(plotDates,1),size(smoothDailyPlot,1));
        if (overlayStartMonth ==1)
            legendYearList(currYearI) = cellstr(num2str(currYear));
        else
            legendYearList(currYearI) = cellstr([num2str(currYear) '/' datestr(currEndDate,'yy') ]);
        end        
        if (~isempty(goodDayI))
            if (currYear < endYear)
                tsAx = plot(plotDates(1:endI),smoothDailyPlot(1:endI),'Color',annualLineColors(currYearI,:),'LineWidth',allLineWidth(currYearI),'LineStyle',char(allLineStyle(currYearI)));
            else
                tsAx = plot(plotDates(1:endI),smoothDailyPlot(1:endI),'Color',annualLineColors(currYearI,:),'LineWidth',allLineWidth(currYearI),'LineStyle',char(allLineStyle(currYearI)));
            end
            figHList(currYearI) = tsAx;
        end
    end


    if (isempty(find(isnan(yLim))))
        ylim(yLim);
    else
        currYLim = get(gca,'YLim');
        set(gca,'YLim',[0 currYLim(2)]);
    end
	newYLim = get(gca,'YLim');

    ylabel(yLabel);

 %   if (panelIndexCol == 2)
        legH = legend(figHList(nPlotYears-nColorYears+1:end),legendYearList(nPlotYears-nColorYears+1:end),'Location','Southeast');
        set(legH,'box','off');
%    end
    set(gca,'FontSize',fontSize);
    title(titleStr,'FontSize',fontSize);            

     set(currSubAx,'Color',[1 1 1]);% box off doesn't work,'box','off')
 %   set(currSubAx,'XColor',get(figH,'color'));
 %   set(currSubAx,'YColor',get(figH,'color'));
    set(currSubAx,'Units','inches');    
    datetick('x',dateTickFormat);
    xTickLabels = cellstr(get(gca,'XTickLabel'));
    set(gca,'XTickLabel',[xTickLabels(1:end-1); cellstr('')]);

   
    %currXPos = (xMargin +  (xGap+axWidth)*(panelIndexCol - 1)) ;
    %currYPos = (yMargin + (yGap+axHeight)*(panelIndexRow -1 )) ;
   % newAxPos = [currXPos currYPos axWidth axHeight];
   % set(currSubAx,'Position',newAxPos);    
    
  
end



function [cmap]=buildcmap(colors)
% [cmap]=buildcmap(colors)
%
% This function can be used to build your own custom colormaps. Imagine if
% you want to display rainfall distribution map. You want a colormap which
% ideally brings rainfall in mind, which is not achiveved by colormaps such
% as winter, cool or jet and such. A gradient of white to blue will do the
% task, but you might also use a more complex gradient (such as
% white+blue+red or colors='wbr'). This function can be use to build any
% colormap using main colors rgbcmyk. In image processing, w (white) can be
% used as the first color so that in the output, the background (usually
% with 0 values) appears white. In the example of rainfall map, 'wb' will
% produce a rainfall density map where the background (if its DN values are
% 0) will appear as white.
%
% Inputs:
%  colors: string (char) of color codes, any sequence of rgbcmywk
%  representing different colors (such as 'b' for blue) is acceptable. If a
%  gradient of white to blue is needed, colors would be 'wb'; a rainbow of
%  white+blue+red+green would be 'wbrg'.
%
% Example:
%  [cmap]=buildcmap('wygbr');
% %try the output cmap:
% im=imread('cameraman.tif');
% imshow(im), colorbar
% colormap(cmap) %will use the output colormap
%
% First version: 14 Feb. 2013
% sohrabinia.m@gmail.com
%--------------------------------------------------------------------------

if nargin<1
    colors='wrgbcmyk';
end

if ~ischar(colors)
    error(['Error! colors must be a variable of type char with '...
        'color-names, such as ''r'', ''g'', etc., '...
        'type ''help buildcmap'' for more info']);
end

ncolors=length(colors)-1;


bins=round(255/ncolors);
% diff1=255-bins*ncolors;

vec=zeros(300,3);

switch colors(1)
    case 'w'
        vec(1,:)=1;
    case 'r'
        vec(1,:)=[1 0 0];
    case 'g'
        vec(1,:)=[0 1 0];
    case 'b'
        vec(1,:)=[0 0 1];
    case 'c'
        vec(1,:)=[0 1 1];
    case 'm'
        vec(1,:)=[1 0 1];
    case 'y'
        vec(1,:)=[1 1 0];
    case 'k'
        vec(1,:)=[0 0 0];
end


for i=1:ncolors
 beG=(i-1)*bins+1;
 enD=i*bins+1; %beG,enD
 switch colors(i+1)
     case 'w'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD,
     case 'r'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';%colors(i+1),beG,enD
     case 'g'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';%colors(i+1),beG,enD
     case 'b'         
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD
     case 'c'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD
     case 'm'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';
     case 'y'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
     case 'k'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
 end
end
cmap=vec(1:bins*ncolors,:);
end %end of buildcmap

function y = wmean(x,w,dim)
%WMEAN   Weighted Average or mean value.
%   For vectors, WMEAN(X,W) is the weighted mean value of the elements in X
%   using non-negative weights W. For matrices, WMEAN(X,W) is a row vector 
%   containing the weighted mean value of each column.  For N-D arrays, 
%   WMEAN(X,W) is the weighted mean value of the elements along the first 
%   non-singleton dimension of X.
%
%   Each element of X requires a corresponding weight, and hence the size 
%   of W must match that of X.
%
%   WMEAN(X,W,DIM) takes the weighted mean along the dimension DIM of X. 
%
%   Class support for inputs X and W:
%      float: double, single
%
%   Example:
%       x = rand(5,2);
%       w = rand(5,2);
%       wmean(x,w)

if nargin<2
    error('Not enough input arguments.');
end

% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    error('Inputs x and w must be the same size.');
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
if (all(w(:)==0))
    error('At least one weight must be non-zero.');
end

if nargin==2, 
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end

y = sum(w.*x,dim)./sum(w,dim);

end



function boxString = GetCoordBoxString(latRange,lonRange)
    if (latRange(1) < 0)
        latMinStr = [num2str(abs(latRange(1))) 'S'];
    else
        latMinStr = [num2str(latRange(1)) 'N'];
    end    

    if (latRange(2) < 0)
        latMaxStr = [num2str(abs(latRange(2))) 'S'];
    else
        latMaxStr = [num2str(latRange(2)) 'N'];
    end

    if (lonRange(1) < 0)
        lonMinStr = [num2str(abs(lonRange(1))) 'W'];
    else
        lonMinStr = [num2str(lonRange(1)) 'E'];
    end

    if (lonRange(2) < 0)
        lonMaxStr = [num2str(abs(lonRange(2))) 'W'];
    else
        lonMaxStr = [num2str(lonRange(2)) 'E'];
    end    
    boxString = [latMinStr ' to ' latMaxStr ', ' lonMinStr ' to ' lonMaxStr ];
end
