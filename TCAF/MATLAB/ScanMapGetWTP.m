function [nominalWTP, valid, validity, offTerMap] = ScanMapGetWTP(predTraj, MAP_LEVEL)
%ScanMapGetWTP 지형지물 탐색 및 WTP 생성
% 입력받은 예측 비행경로를 이용해 지형지물 정보를 탐색하고 WCP를 생성한다
%
%   [nominalWTP, valid] = ScanMapGetWTP(predTraj)
%
% INPUT
% - predTraj(ALONG_TRACK_STEPS, 3, 3)
%  *ALONG_TRACK_STEPS = 81
% - MAP_LEVEL
%
% OUTPUT
% - nominalWTP(41, 1) : WTP array (0:0.5:20)
% - valid(41, 1)
% - validity
% - offTerMap

% -------------------------------------------------------------------------
ON_DEBUG = true;
ON_PLOT = true;
ON_PLOT2 = false;
% -------------------------------------------------------------------------
MAP_LEVEL = 2;
% DEFINE CONSTANTS
if MAP_LEVEL == 1
    RESOLUTION = 3;
elseif MAP_LEVEL == 2
    RESOLUTION = 1;
else
    RESOLUTION = 3;
end
ARCSEC2DEG = 0.000277777777777777;
RESDEG = RESOLUTION * ARCSEC2DEG;

% -------------------------------------------------------------------------
if ON_DEBUG
    % MAP data
    load('C:\Users\ASCL\Documents\GitHub\DTNS\MAP\MATLAB\MAPDATA\DEM_N35_E128.mat');
    % predicted Trajectory
    predTraj = PredictTrajectory();
    scanned = [];
end
% -------------------------------------------------------------------------

% Initialize
nominalWTP = zeros(41, 1, 'single');
valid = zeros(41, 1, 'logical');
validity = false;
offTerMap = false;

% For every 0.5 sec (from 0 to 20 sec)
for ii = 1:2:size(predTraj, 1)

    % Initialize current & highest elevation
    highestElev = 0.0;
    curElev = 0.0;

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
    trapezoid = ...
        [ predTraj(lo, 2, 1), predTraj(lo, 2, 2) ...
        ; predTraj(hi, 2, 1), predTraj(hi, 2, 2) ...
        ; predTraj(hi, 3, 1), predTraj(hi, 3, 2) ...
        ; predTraj(lo, 3, 1), predTraj(lo, 3, 2) ...
        ; predTraj(lo, 2, 1), predTraj(lo, 2, 2)];
    
    % Bounding box coordinated at map post
    bBox = [min(trapezoid); max(trapezoid)];
    bBox = [...
        floor(bBox(1,:)/RESDEG)*RESDEG; ...
        ceil( bBox(2,:)/RESDEG)*RESDEG];

    % number of map post
    bBoxLatNum = round(abs((bBox(2,1)-bBox(1,1)) / RESDEG));
    bBoxLonNum = round(abs((bBox(2,2)-bBox(1,2)) / RESDEG));
% -------------------------------------------------------------------------
    if ON_PLOT2
        figure;
        hold on;
        plot(trapezoid(:,2), trapezoid(:,1), '-ko', 'MarkerFaceColor', 'black');
        plot([bBox(1,2), bBox(2,2), bBox(2,2), bBox(1,2), bBox(1,2)], [bBox(1,1), bBox(1,1), bBox(2,1), bBox(2,1), bBox(1,1)], '-b');
        h = animatedline('Marker','o', 'LineStyle','none', 'MarkerEdgeColor', 'red');
        grid minor;
        set(gca,'ytick', floor( min(predTraj(:,1,1))):1/240:ceil(max(predTraj(:,1,1))) );
        set(gca,'xtick', floor( min(predTraj(:,1,2))):1/240:ceil(max(predTraj(:,1,2))) );
        xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
        axis equal;
        drawnow
    end
% -------------------------------------------------------------------------
    for jj = 0:bBoxLatNum
        for kk = 0:bBoxLonNum

            lat = bBox(1,1) + (bBox(2,1)-bBox(1,1))/bBoxLatNum * (jj);
            lon = bBox(1,2) + (bBox(2,2)-bBox(1,2))/bBoxLonNum * (kk);

            isIn = inpolygon(lat, lon, trapezoid(:,1), trapezoid(:,2));

            if ~isIn
                [dist, dx, dy] = distFromPoly(lat, lon, trapezoid(:,1), trapezoid(:,2));

                isAround = ...
                    ( dist - RESDEG < 0) ...
                    | ( (dx - RESDEG < 0) & (dy - RESDEG <0) );

                if isAround
                else
                    continue
                end
            else
            end
%                 curElev = coder.ceval('GetElevation', lat, lon);
% -------------------------------------------------------------------------
                if ON_PLOT2
                    addpoints(h, lon, lat);
                end
                if true
                    curElev = fix(100 * rand);
                    scanned = [scanned; lat, lon];
                end
% -------------------------------------------------------------------------
            if curElev == -1.0
                offTerMap = true;
            elseif curElev > highestElev
                highestElev = curElev;
            else
            end
        end
    end

    % Update
    nominalWTP((ii-1)/2+1, 1) = single(highestElev);
    valid((ii-1)/2+1, 1) = ~offTerMap;

end

validity = all(valid(1:11, 1));

% -------------------------------------------------------------------------
if ON_PLOT
    % Plot
    figure('WindowStyle', 'docked');
    hold on;
    plot(predTraj(:,1,2), predTraj(:,1,1), 'Color', '#0072BD', 'LineWidth', 1);
    plot(predTraj(:,2,2), predTraj(:,2,1), 'Color', '#77AC30', 'LineWidth', 2);
    plot(predTraj(:,3,2), predTraj(:,3,1), 'Color', '#D95319', 'LineWidth', 2);
    plot(scanned(:,2), scanned(:,1), 'ro');
    xmin = min(predTraj(:,1,2), [], 'all') - 0.005;
    xmax = max(predTraj(:,1,2), [], 'all') + 0.005;
    ymin = min(predTraj(:,1,1), [], 'all') - 0.005;
    ymax = max(predTraj(:,1,1), [], 'all') + 0.01;
    grid minor;
    set(gca,'ytick', floor( min(predTraj(:,1,1))):1/240:ceil(max(predTraj(:,1,1))) );
    set(gca,'xtick', floor( min(predTraj(:,1,2))):1/240:ceil(max(predTraj(:,1,2))) );
    xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
    axis equal;
end
% -------------------------------------------------------------------------
end
% Internal functions
function [dist, dx, dy] = distFromPoly(x, y, xv, yv)
%distFromPoly
% INPUT
% -  x      : point's x coordinate
% -  y      : point's y coordinate
% -  xv     : vector of polygon vertices x coordinates
% -  yv     : vector of polygon vertices y coordinates
% * (xv, yv) should be closed!
%
% OUTPUT
% - dist    : distance from point to polygon (defined as a minimal distance from 
%       point to any of polygon's ribs, positive if the point is outside the
%       polygon and negative otherwise)
% - dx      : 
% - dy      : 
% -------------------------------------------------------------------------

% linear parameters of segments that connect the vertices
% Ax + By + C = 0
A = -diff(yv);
B =  diff(xv);
C = yv(2:end).*xv(1:end-1) - xv(2:end).*yv(1:end-1);

% find the projection of point (x,y) on each rib
AB = 1./(A.^2 + B.^2);
vv = (A*x+B*y+C);
xp = x - (A.*AB).*vv;
yp = y - (B.*AB).*vv;

% Test for the case where a polygon rib is either horizontal or vertical
id_xv = find(diff(xv)==0);
xp(id_xv)=xv(id_xv);
id_yv = find(diff(yv)==0);
yp(id_yv)=yv(id_yv);

% find all cases where projected point is inside the segment
idx_x = (((xp>=xv(1:end-1)) & (xp<=xv(2:end))) | ((xp>=xv(2:end)) & (xp<=xv(1:end-1))));
idx_y = (((yp>=yv(1:end-1)) & (yp<=yv(2:end))) | ((yp>=yv(2:end)) & (yp<=yv(1:end-1))));
idx = idx_x & idx_y;

% distance from point (x,y) to the vertices
dv = sqrt((xv(1:end-1)-x).^2 + (yv(1:end-1)-y).^2);
if(~any(idx)) % all projections are outside of polygon ribs
    [dist,I] = min(dv);
    dx = abs(xv(I)-x);
    dy = abs(yv(I)-y);
else
    % distance from point (x,y) to the projection on ribs
    dp = sqrt((xp(idx)-x).^2 + (yp(idx)-y).^2);
    [min_dv,I] = min(dv);
    [min_dp,I2] = min(dp);
    [dist,I3] = min([min_dv min_dp]);
    if I3 == 1
        dx = abs(xv(I)-x);
        dy = abs(yv(I)-y);
    else
        dx = abs(xp(I2)-x);
        dy = abs(yp(I2)-y);
    end
end

if(inpolygon(x, y, xv, yv))
    dist = -dist;
end

end