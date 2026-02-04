function main()
%{
    This makes a table with coordinates on a grid over a lat/lon box at a  specified resolution
    in degrees and writes it to a csv file. Each grid point is labeled with a number
    from left to write and then down. It skips over points if not enough of
    an underlying land mask is classified as land.
    
    This is then used by a python script to download data from the 'Climate
    Data Store' using their single-point download tool for that table of
    coordinates. So far, it's been for ERA5-Land.

    More info ERA5-Land single point is here:
    https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries?tab=overview
%}

    % for running on different machines
    machineRoot = '/Users/rfield1/data/observations/';
   % machineRoot = ['/autofs/brewer/rfield1/storage/observations/'];

    % where to write the table
    outputRoot = [machineRoot '/ERA5Land/reanalysis-era5-land-timeseries/'];

    % This is the 'Land-sea' ERA5 mask under the 'Invariant' full gridded ERA5Land fields here:
    %https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=download
    landSeaMaskFN = [machineRoot 'ERA5Land/reanalysis-era5-land/LandSeaMask.nc'];

    % a descriptive name and the grid particulars
    regionName = 'IberianPeninsula';
    mapLatLim = [36 44];
    mapLonLim = [-10 4];    
    gridSpacingDegrees = .25;
    minLandFrac = 0.5;

%{    
    regionName = 'WestIberianPeninsula';
    mapLatLim = [36 44];
    mapLonLim = [-9 -4];
    gridSpacingDegrees = .25;
%}

    % make the outut directory if it doesn't exist
    outputDir = [outputRoot regionName '/'];
    if (~exist(outputDir,'dir'))
        mkdir(outputDir);
    end

    % read in the land-sea mask
    landSeaMask = ncread(landSeaMaskFN,'lsm');
    landSeaLat = ncread(landSeaMaskFN,'latitude');
    landSeaLon = ncread(landSeaMaskFN,'longitude');
    
    % this is the accumulated number of 'good points', might be less than
    % the rowsxcolumn dimensions because some points don't have enough land
    nGridPoints = 0;
    for currLat = mapLatLim(1):gridSpacingDegrees:mapLatLim(2)
        for currLon = mapLonLim(1):gridSpacingDegrees:mapLonLim(2)

            currLonInTheirCoords = currLon;
            if (currLon<0) 
                currLonInTheirCoords = 360+currLon; % barf
            end

            [latDist,closestLatI] = min((abs(landSeaLat-currLat)));
            [lonDist,closestLonI] = min((abs(landSeaLon-currLonInTheirCoords)));
            currMaskValue = landSeaMask(closestLonI,closestLatI);

            % add point to table if there's enough land
            if (currMaskValue > minLandFrac)
                nGridPoints = nGridPoints + 1;
                gridPointIDs(nGridPoints,1) = cellstr([regionName num2str(nGridPoints,'%.03d')]);
                gridPointLats(nGridPoints,1) = currLat;
                gridPointLons(nGridPoints,1) = currLon;
                gridPointAreas(nGridPoints,1) = NaN;
            else
                disp('not enough land')
            end
        end
    end
    close all;
    figure('Visible','on')
    scatter(gridPointLons,gridPointLats,'filled');
    % write the table to a csv file
    outFN = [outputDir regionName '.gridSpacingDegrees.' num2str(gridSpacingDegrees) '.DownloadCoords.csv'];
    outTable = [gridPointIDs num2cell(gridPointLats) num2cell(gridPointLons) num2cell(gridPointAreas)];
    columnHeaders = [cellstr('ID') cellstr('Lat') cellstr('Lon') cellstr('Area')];
    outTable = table([columnHeaders; outTable]);
    writetable(outTable,outFN,'WriteVariableNames',0);
end