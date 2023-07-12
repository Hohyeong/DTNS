% !!!!!!!!!!!!!!!!!!!!!!!!!!!GetElevationFromIndex 사용 버전임
% function [nominalWTP, valid, validity] = ScanMapGetWTP(predTraj)
%ScanMapGetWTP 지형지물 탐색 및 WTP 생성
% 입력받은 예측 비행경로를 이용해 지형지물 정보를 탐색하고 WCP를 생성한다
%
%   [nominalWTP, valid] = ScanMapGetWTP(predTraj)
%
% INPUT
% - predTraj(ALONG_TRACK_STEPS, 3, 3)
%  *ALONG_TRACK_STEPS = 81
%
% OUTPUT
% - nominalWTP(41, 1) : WTP array (0:0.5:20)
% - valid(41, 1)
% - validity

% -------------------------------------------------------------------------
ON_DEBUG = true;
ON_PLOT = true;
% -------------------------------------------------------------------------

% DEFINE CONSTANTS
RESOLUTION = 3;
ARCSEC2DEG = 0.000277777777777777;
RESDEG = RESOLUTION * ARCSEC2DEG;
% -------------------------------------------------------------------------
switch ON_DEBUG
    case true
        if exist('predTraj', 'var')
        else
            predTraj = PredictTrajectory();
        end
        if exist('DEMData', 'var')
        else
            load('C:\Users\ASCL\Documents\GitHub\DTNS\MAP\MATLAB\MAPDATA\DEM_N35_E128.mat');
        end

        scanned = [];

    otherwise

end
% -------------------------------------------------------------------------

% Initialize
nominalWTP = zeros(41, 1, 'single');
valid = zeros(41, 1, 'logical');
minLatIdx = int32(0);
minLonIdx = int32(0);
maxLatIdx = int32(0);
maxLonIdx = int32(0);

% For every 0.5 sec (from 0 to 20 sec)
for ii = 1:2:size(predTraj, 1)

    % Configure High and Low Indices
    if ii == 1
        lo = ii;
        hi = ii+1;
    elseif ii == size(predTraj, 1)
        lo = ii-1;
        hi = ii;
    else
        lo = ii-1;
        hi = ii+1;
    end
    
    % Position of Trapezoid. [Lat, Lon, Height]
    posLL = [predTraj(lo, 2, 1), predTraj(lo, 2, 2)];
    posUL = [predTraj(hi, 2, 1), predTraj(hi, 2, 2)];
    posUR = [predTraj(hi, 3, 1), predTraj(hi, 3, 2)];
    posLR = [predTraj(lo, 3, 1), predTraj(lo, 3, 2)];
    
    % Bounding Box of Vertices
    bBox = zeros(4,2);

    if posLL(1) >= posUL(1)
        bBox(1,1) = ceil(posLL(1)/RESDEG)*RESDEG;
        bBox(4,1) = floor(posUL(1)/RESDEG)*RESDEG;
    else
        bBox(1,1) = floor(posLL(1)/RESDEG)*RESDEG;
        bBox(4,1) = ceil(posUL(1)/RESDEG)*RESDEG;
    end

    if posLL(2) >= posLR(2)
        bBox(1,2) = ceil(posLL(2)/RESDEG)*RESDEG;
        bBox(2,2) = floor(posLR(2)/RESDEG)*RESDEG;
    else
        bBox(1,2) = floor(posLL(2)/RESDEG)*RESDEG;
        bBox(2,2) = ceil(posLR(2)/RESDEG)*RESDEG;
    end

    if posLR(1) >= posUR(1)
        bBox(2,1) = ceil(posLR(1)/RESDEG)*RESDEG;
        bBox(3,1) = floor(posUR(1)/RESDEG)*RESDEG;
    else
        bBox(2,1) = floor(posLR(1)/RESDEG)*RESDEG;
        bBox(3,1) = ceil(posUR(1)/RESDEG)*RESDEG;
    end

    if posUR(2) >= posUL(2)
        bBox(3,2) = ceil(posUR(2)/RESDEG)*RESDEG;
        bBox(4,2) = floor(posUL(2)/RESDEG)*RESDEG;
    else
        bBox(3,2) = floor(posUR(2)/RESDEG)*RESDEG;
        bBox(4,2) = ceil(posUL(2)/RESDEG)*RESDEG;
    end

    % Min/Max of Bounding Box of Vertices
    bBoxMinMax = [min(bBox(:,1)), min(bBox(:,2)); ...
                  max(bBox(:,1)), max(bBox(:,2))];

    bBoxMinMaxIdx = zeros(2, 2, 'int32');
%     bBoxMinMaxIdx(1,1) = coder.ceval('IndexOfLat', bBoxMinMax(1,1));
%     bBoxMinMaxIdx(2,1) = coder.ceval('IndexOfLat', bBoxMinMax(2,1));
%     bBoxMinMaxIdx(1,2) = coder.ceval('IndexOfLon', bBoxMinMax(1,2));
%     bBoxMinMaxIdx(2,2) = coder.ceval('IndexOfLon', bBoxMinMax(2,2));

    bBoxIdx = zeros(4, 2, 'int32');
%     bBoxIdx(1,1) = coder.ceval('IndexOfLat', bBox(1,1));
%     bBoxIdx(2,1) = coder.ceval('IndexOfLat', bBox(2,1));
%     bBoxIdx(3,1) = coder.ceval('IndexOfLat', bBox(3,1));
%     bBoxIdx(4,1) = coder.ceval('IndexOfLat', bBox(4,1));
%     bBoxIdx(1,2) = coder.ceval('IndexOfLon', bBox(1,2));
%     bBoxIdx(2,2) = coder.ceval('IndexOfLon', bBox(2,2));
%     bBoxIdx(3,2) = coder.ceval('IndexOfLon', bBox(3,2));
%     bBoxIdx(4,2) = coder.ceval('IndexOfLon', bBox(4,2));
    % -------------------------------------------------------------------------
    bBoxMinMaxIdx = [IndexOfLat(bBoxMinMax(1,1)), IndexOfLon(bBoxMinMax(1,2)); ...
                     IndexOfLat(bBoxMinMax(2,1)), IndexOfLon(bBoxMinMax(2,2))];
    bBoxIdx = [IndexOfLat(bBox(1,1)), IndexOfLon(bBox(1,2)); ...
               IndexOfLat(bBox(2,1)), IndexOfLon(bBox(2,2)); ...
               IndexOfLat(bBox(3,1)), IndexOfLon(bBox(3,2)); ...
               IndexOfLat(bBox(4,1)), IndexOfLon(bBox(4,2))];
    % -------------------------------------------------------------------------

    highestElev = 0.0;
    curElev = 0.0;

    for jj = bBoxMinMaxIdx(1,1):bBoxMinMaxIdx(2,1)
        for kk = bBoxMinMaxIdx(1,2):bBoxMinMaxIdx(2,2)
            
            testIdxLat = jj;
            testIdxLon = kk;

            if ((testIdxLat < 0) || (testIdxLon < 0))
                isIn = false;
            else
                isIn = inpolygon(jj, kk, bBoxIdx(:,1), bBoxIdx(:,2));
            end
            
            if isIn == true
                
%                 curElev = coder.ceval('GetElevationFromIndex', testIdxLat, testIdxLon);
                % -------------------------------------------------------------------------
                if ON_DEBUG
                    curElev = DEMData(testIdxLat, testIdxLon);
                    scanned = [scanned; DBinfo.minLat + RESDEG * (testIdxLat), DBinfo.minLon + RESDEG * (testIdxLon)];
                end
                % -------------------------------------------------------------------------
                if curElev > highestElev
                    highestElev = curElev;
                else
                end
            else
            end
        end
    end

    nominalWTP((ii-1)/2+1, 1) = single(highestElev);
    valid((ii-1)/2+1, 1) = true;

end

validity = all(valid);

% -------------------------------------------------------------------------

if ON_PLOT

    figure('WindowStyle', 'docked');
    hold on;
    plot(nominalWTP);   
    hold off;

    figure('WindowStyle', 'docked');
    hold on;
    plot(predTraj(:,1,2), predTraj(:,1,1), 'Color', '#0072BD', 'LineWidth', 1);
    plot(predTraj(:,2,2), predTraj(:,2,1), 'Color', '#77AC30', 'LineWidth', 2);
    plot(predTraj(:,3,2), predTraj(:,3,1), 'Color', '#D95319', 'LineWidth', 2);
    plot(scanned(:,2), scanned(:,1), 'ro');
    grid minor;
    xmin = min(predTraj(:,1,2), [], 'all') - 0.005;
    xmax = max(predTraj(:,1,2), [], 'all') + 0.005;
    ymin = min(predTraj(:,1,1), [], 'all') - 0.005;
    ymax = max(predTraj(:,1,1), [], 'all') + 0.01;
    axis([xmin xmax ymin ymax]);
    for ii = 1:81
    h = plot(predTraj(ii, [2, 3], 2), predTraj(ii, [2, 3], 1), '-k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    set(gca,'ytick', floor(predTraj(:,1,1)):1/240:ceil(predTraj(:,1,1)));
    set(gca,'xtick', floor(predTraj(:,1,2)):1/240:ceil(predTraj(:,1,2)));
    xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
    axis equal;
end
% -------------------------------------------------------------------------

% end


% -------------------------------------------------------------------------
function idxLat = IndexOfLat(lat)
%INDEXOFLAT
persistent RESOLUTION ARCSEC2DEG DEG2RAD MAPCACHE_SIZE leftBottomLat
if isempty(RESOLUTION)
    RESOLUTION = 3;
end
if isempty(ARCSEC2DEG)
    ARCSEC2DEG = 0.000277777777777777;
end
if isempty(DEG2RAD)
    DEG2RAD = 0.0174532925199433;
end
if isempty(MAPCACHE_SIZE)
    MAPCACHE_SIZE = 1201;
end
if isempty(leftBottomLat)
    leftBottomLat = 35;
end

d = RESOLUTION * ARCSEC2DEG * DEG2RAD;

idxLat = floor(((lat - leftBottomLat ) * DEG2RAD) / d);

if ( (idxLat < 0) || (idxLat >= MAPCACHE_SIZE) )
    idxLat = -1;
else
    idxLat = idxLat + 1; % One-based Index
end

end

% -------------------------------------------------------------------------
function idxLon = IndexOfLon(lon)
%INDEXOFLON
persistent RESOLUTION ARCSEC2DEG DEG2RAD MAPCACHE_SIZE leftBottomLon
if isempty(RESOLUTION)
    RESOLUTION = 3;
end
if isempty(ARCSEC2DEG)
    ARCSEC2DEG = 0.000277777777777777;
end
if isempty(DEG2RAD)
    DEG2RAD = 0.0174532925199433;
end
if isempty(MAPCACHE_SIZE)
    MAPCACHE_SIZE = 1201;
end
if isempty(leftBottomLon)
    leftBottomLon = 128;
end

d = RESOLUTION * ARCSEC2DEG * DEG2RAD;

idxLon = floor(((lon - leftBottomLon ) * DEG2RAD) / d);

if ( (idxLon < 0) || (idxLon >= MAPCACHE_SIZE) )
    idxLon = -1;
else
    idxLon = idxLon + 1; % One-based Index
end

end