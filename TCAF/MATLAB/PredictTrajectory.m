function predTraj = PredictTrajectory(eulerDeg, velXYZ, velNED, accXYZ, lla, uncHorPos, uncHdg)
%PREDICTTRAJECTORY 예측 비행경로를 생성
% 예측 비행경로 배열을 생성 (ALONG_TRACK_STEPS, 3, 3)
%                          (TimeStep, [-3σ,0,+3σ], LLA)
% 설계사양
% - 종방향 비행경로 예측시간: 20 sec + (alpha)
% - 종방향 비행경로 예측시간 간격: 0.25 sec
% - 횡방향 탐색 너비 :
%          3 x min( [σ ∝ TRN 수평위치 불확실도], [σ = 포화(35.56, 228.6)] )
% - 탐색 패턴: 부등변 사각형(trapezoid)
% 
% 운용조건
% - 운용속도: 관성속도 100 ~ 700 knots (51.4 ~ 360.1 m/s)

% Define constants
GRAV_ACC = 9.8;
ALONG_TRACK_STEPS = 81;
CROSS_TRACK_STEPS = 3;
DELTA_T = 0.25;

ON_DEBUG = true;
ON_PLOT = true;
% -------------------------------------------------------------------------
if ON_DEBUG
%         ALONG_TRACK_STEPS = 628;
        ALONG_TRACK_STEPS = 81;
        CROSS_TRACK_STEPS = 3;
        eulerDeg = [55, 0, 10];
        velXYZ = [300, 0, 0];
        velNED = [295.4423, 52.0945, 0];
        accXYZ = [0, 0, 0];
        lla = [35.4, 128.6, 1000.0];
        uncHorPos = 15;
        uncHdg = 1.5;
        
        eulerDeg = [-0.9393, -1.4172, -143.3661];
        eulerDeg = [-40, 0, 0];
%         eulerDeg = [-0.9393, -1.4172, 143.3661];
%         eulerDeg = [0.9393, -1.4172, -43.3661];

        
        velXYZ = [199.8955, 0.0127, -3.8525];
        velNED = convertBody2NED(eulerDeg, velXYZ);
        accXYZ = [0, 0.1429, 1.6764];

end
% -------------------------------------------------------------------------


% Initialize Predicted Trajectory
predTraj = zeros(ALONG_TRACK_STEPS, CROSS_TRACK_STEPS, 3);

% Sature Uncertainties 
uncHorPos = saturate(uncHorPos, 17.78, 228.6);
uncHdg = saturate(uncHdg, 0.1, 5);

% [1] Nominal Projected Trajcetory 명목 궤적 예측
predTraj(:, 1, :) = ...
    ProjectTrajcetory(eulerDeg, velXYZ, velNED, accXYZ, lla);

% [2] Inner Boundary (선회 안쪽)
predTraj(:, 2, :) = ...
    PredictInnerTrajectory(eulerDeg, velXYZ, velNED, accXYZ, lla, uncHorPos, uncHdg);

% [3] Outer Boundary (선회 바깥쪽)
predTraj(:, 3, :) = ...
    PredictOuterTrajectory(eulerDeg, velXYZ, velNED, accXYZ, lla, uncHorPos, uncHdg);

% Plot
%--------------------------------------------------------------------------
if ON_PLOT
    figure('WindowStyle', 'docked');
    hold on;
    plot(predTraj(1,1,2), predTraj(1,1,1), 'ro');
    plot(predTraj(:,1,2), predTraj(:,1,1), 'Color', '#0072BD', 'LineWidth', 1);
    plot(predTraj(:,2,2), predTraj(:,2,1), 'Color', '#77AC30', 'LineWidth', 2);
    plot(predTraj(:,3,2), predTraj(:,3,1), 'Color', '#D95319', 'LineWidth', 2);
    grid minor;
    axis equal;
    xmin = min(predTraj(:,1:3,2), [], 'all') - 0.005;
    xmax = max(predTraj(:,1:3,2), [], 'all') + 0.005;
    ymin = min(predTraj(:,1:3,1), [], 'all') - 0.005;
    ymax = max(predTraj(:,1:3,1), [], 'all') + 0.01;
    % disp([xmin xmax ymin ymax]);
    axis([xmin xmax ymin ymax]);
    for ii = 1:ALONG_TRACK_STEPS
        h = plot(predTraj(ii, [2, 3], 2), predTraj(ii, [2, 3], 1), '-k');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    set(gca,'ytick', floor(lla(1)):1/240:ceil(lla(1)));
    set(gca,'xtick', floor(lla(2)):1/240:ceil(lla(2)));
    xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
    legend('Start point', 'Projected flight path', 'Inner(Roll-maintained) path', 'Outer(Roll-recovered) path');

%     figure('WindowStyle', 'docked');
%     tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
%     nexttile;
%     hold on;
%     plot(predTraj(:,1,1), '-g', 'LineWidth', 1);
%     plot(predTraj(:,idx2,1), '-b', 'LineWidth', 1);
%     plot(predTraj(:,idx3,1), '-r', 'LineWidth', 1);
%     legend('Projected', 'Roll-Recoverd', 'Roll-Maintained');
%     ylabel('Latitude (deg)');
%     grid on;
%     nexttile;
%     hold on;
%     plot(predTraj(:,1,2), '-g', 'LineWidth', 1);
%     plot(predTraj(:,idx2,2), '-b', 'LineWidth', 1);
%     plot(predTraj(:,idx3,2), '-r', 'LineWidth', 1);
%     legend('Projected', 'Roll-Recoverd', 'Roll-Maintained');
%     ylabel('Longitude (deg)');
%     grid on;
%     nexttile;
%     hold on;
%     plot(predTraj(:,1,3), '-g', 'LineWidth', 1);
%     plot(predTraj(:,idx2,3), '-b', 'LineWidth', 1);
%     plot(predTraj(:,idx3,3), '-r', 'LineWidth', 1);
%     legend('Projected', 'Roll-Recoverd', 'Roll-Maintained');
%     grid on;
%     ylabel('Height (m)');
else
end

%--------------------------------------------------------------------------
end

%% Internal Functions

function [latLength, lonLength] = getLengthOfADegree(lat, h)
%GETLENGTHOFADEGREE
% Calculate the length of one degree of both latitude and longitude
% for a specified latitude and height
%
% == Inputs =========
% lat               - Latitude [deg]
% h                 - Height [m]
% 
% == Outputs ========
% latLength         - Length of 1 deg of latitude on the ref. spheroid [m]
% lonLength         - Length of 1 deg of longitude on the ref. spheroid [m]
% ===================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    h   = 0.0;
end

% WGS84 Parameters
eRadius = 6378137.0;
eEccen = 0.08181919;

% Radius of curvature in prime vertical
nE = eRadius / sqrt( 1 - (eEccen^2*(sind(lat))^2) );
% Meridian Radius of Curavature
mE = eRadius * (1 - eEccen^2) / (1 - eEccen^2*(sind(lat))^2)^(1.5);

lonLength = (nE + h) * cosd(lat) * (pi/180);
latLength = (mE + h) * (pi/180);

end

function nominalTraj = ProjectTrajcetory(eulerDeg, velXYZ, velNED, accXYZ, lla)
%ProjectTrajectory 예측 비행경로를 생성
% 예측 비행경로 배열을 생성 (81, 1, 3)

persistent DELTA_T

% Define constants
if isempty(DELTA_T)
    DELTA_T = 0.25;
end

% Define constants
GRAV_ACC = 9.8;
LIMIT_TAN_ROL = 4;
LIMIT_TURN_RATE = 0.222; % [rad/s]
RAD2DEG = 180/pi;
ALONG_TRACK_STEPS = 81;

% Initialize Predicted Trajectory
nominalTraj = zeros(ALONG_TRACK_STEPS, 1, 3);

% Input Handling
rol = eulerDeg(1);
pit = eulerDeg(2);
hdg = eulerDeg(3);

velX = velXYZ(1);
velY = velXYZ(2);
velZ = velXYZ(3);

velN = velNED(1);
velE = velNED(2);
velD = velNED(3);

lat = lla(1);
lon = lla(2);
alt = lla(3);

% Flight path angle (deg)
gamma = atan2d(-velD, sqrt(velN^2+velE^2));

% Inertial speed
iSpd = sqrt( velX^2 + velY^2 + velZ^2 );
    
% Saturate inertial speed (since TAS is unavailable, use inertial spd) 
iSpd = saturate(iSpd, 75.1089, 500);

% Calculate the length of a degree
[latLength, lonLength] = getLengthOfADegree(lat, alt);
velLat = velN / latLength;
velLon = velE / lonLength;

% Start point
nominalTraj(1, 1, 1) = lat;
nominalTraj(1, 1, 2) = lon;
nominalTraj(1, 1, 3) = alt;

for ii = 2:ALONG_TRACK_STEPS
    
    % [1] Update Attitude
    
    % Saturate tangent of roll (max bank angle = 75.9630 deg)
    tanRol = saturate(tand(rol), -LIMIT_TAN_ROL, LIMIT_TAN_ROL);

    % Saturate inertial speed (since TAS is unavailable, use inertial spd) 
    iSpd = saturate(iSpd, 75.1089, 500);
    
    % Calculate turn rate [rad/s]
    turnRate = GRAV_ACC * tanRol / iSpd;
    
    % Saturate turn rate 
    turnRate = saturate(turnRate, -LIMIT_TURN_RATE, LIMIT_TURN_RATE);
        
    % Update heading [deg]
    hdg = hdg + turnRate * DELTA_T * RAD2DEG;
    
    % Update Flight Path Angle [deg]
    %        reaction time: 1.5 sec
    %        rate of change: 10 deg/s
    if ii < 4
        gamma = gamma;
    elseif gamma < -(10*DELTA_T)
        gamma = gamma + 10 * DELTA_T;
    elseif gamma < 0.0
        gamma = 0.0;
    else
        gamma = gamma;
    end

    % [2] Update Velocity
    velNTmp = sqrt(velN^2+velE^2+velD^2)*cosd(hdg)*cosd(gamma);
    velE    = sqrt(velN^2+velE^2+velD^2)*sind(hdg)*cosd(gamma);
    velN    = velNTmp;
    velDOld = velD;
    velD    = - sqrt(velN^2+velE^2+velD^2)*sind(gamma);

    % [3] Update Position
    % Rate of change in geodetic coordinates [deg/s]
    [latLength, lonLength] = getLengthOfADegree(nominalTraj(ii, 1, 1), nominalTraj(ii, 1, 3));

    velLatOld = velLat;
    velLonOld = velLon;
    velLat  = velN / latLength;
    velLon  = velE / lonLength;
    
    nominalTraj(ii, 1, 1) = nominalTraj(ii-1, 1, 1) + (velLat + velLatOld)/2 * DELTA_T;
    nominalTraj(ii, 1, 2) = nominalTraj(ii-1, 1, 2) + (velLon + velLonOld)/2 * DELTA_T;
    nominalTraj(ii, 1, 3) = nominalTraj(ii-1, 1, 3) - (velD + velDOld)/2 * DELTA_T;

end

end

function outerTraj = PredictOuterTrajectory(eulerDeg, velXYZ, velNED, accXYZ, lla, uncHorPos, uncHdg)
%PredictOuterTrajectory 예측 비행경로를 생성

persistent DELTA_T

% Define constants
if isempty(DELTA_T)
    DELTA_T = 0.25;
end
GRAV_ACC = 9.8;
ALONG_TRACK_STEPS = 81;

% Input Handling
rol = eulerDeg(1);
pit = eulerDeg(2);
hdg = eulerDeg(3);

velX = velXYZ(1);
velY = velXYZ(2);
velZ = velXYZ(3);

velN = velNED(1);
velE = velNED(2);
velD = velNED(3);

lat = lla(1);
lon = lla(2);
alt = lla(3);

% Flight path angle (deg)
gamma = atan2d(-velD, sqrt(velN^2+velE^2));

% Saturate Uncertainties 
uncHorPos = saturate(uncHorPos, 35.56, 228.6);
uncHdg = saturate(uncHdg, 0.1, 5);

% Direction of uncertaitny addition
isTurnRight = (0 <= eulerDeg(1)) && (eulerDeg(1) < 180);
if isTurnRight
    uncHdg = (-1) * uncHdg;
    uncPosRot = -90;
else
    uncPosRot = 90;
end

% 3-sigma uncertainty for Heading
unitDirection = [velN; velE; 0] / norm([velN; velE; 0]);
unitDirection = rotateZ(3*uncHdg) * unitDirection;
hdg = hdg + 3*uncHdg;

% Update Velocity
velRef = sqrt(velN^2+velE^2+velD^2);
velN   = velRef*cosd(hdg)*cosd(gamma);
%         velNTmp = sqrt(velN^2+velE^2+velD^2)*cosd(hdg)*cosd(gamma);
velE    = velRef*sind(hdg)*cosd(gamma);
%         velN    = velNTmp;
velDOld = velD;
%         velD    = - sqrt(velN^2+velE^2+velD^2)*sind(gamma);
velD    = - velRef*sind(gamma);

% 3-sigma uncertainty for position
delPos = 3 * uncHorPos * rotateZ(uncPosRot) * unitDirection;

% Calculate the length of a degree
[latLength, lonLength] = getLengthOfADegree(lat, alt);
velLat  = velN / latLength;
velLon  = velE / lonLength;

% Initialize Predicted Trajectory
outerTraj = zeros(ALONG_TRACK_STEPS, 1, 3);

% Start point
outerTraj(1, 1, 1) = lat + delPos(1)/latLength;
outerTraj(1, 1, 2) = lon + delPos(2)/lonLength;
outerTraj(1, 1, 3) = alt;

rolRate = 120;
keepT = 2.5;

stage = [ floor(keepT/DELTA_T) + 1, ...
          floor(keepT/DELTA_T) + 1 + floor( (abs(rol)/rolRate)/DELTA_T ) + 1];

velXRef = saturate(velX, 75.1089, 360.111);

for ii = 2:ALONG_TRACK_STEPS
    
    if (stage(1) < ii) && (ii <= stage(2))
        % [2] roll recovery
        
        % Update Attitude
        % (roll)
        rol_old = rol;
        if abs(rol-0) <= abs(rolRate*DELTA_T)
            rol_new = 0;
        else
            rol_new = sign(rol) * (abs(rol) - rolRate*DELTA_T);
        end
        
        % (heading)
        if ( abs( tand(rol_old) ) > 4 )
            hdg = hdg + ...
                0.5 * ( GRAV_ACC * 4 * sign( tand(rol_old) ) / velXRef * DELTA_T) * 180/pi +...
                0.5 * ( GRAV_ACC * 4 * sign( tand(rol_new) ) / velXRef * DELTA_T) * 180/pi; 
        else
            hdg = hdg + ...
                0.5 * ( GRAV_ACC * tand(rol_old) / velXRef * DELTA_T) * 180/pi + ...
                0.5 * ( GRAV_ACC * tand(rol_new) / velXRef * DELTA_T) * 180/pi;
        end
        rol = rol_new;

        % Update Flight Path Angle [deg]
        %        reaction time: 1.5 sec
        %        rate of change: 10 deg/s
        if ii < 4
            gamma = gamma;
        elseif gamma < -(10*DELTA_T)
            gamma = gamma + 10 * DELTA_T;
        elseif gamma < 0.0
            gamma = 0.0;
        else
            gamma = gamma;
        end
        
        % [2] Update Velocity
        velRef = sqrt(velN^2+velE^2+velD^2);
        velN   = velRef*cosd(hdg)*cosd(gamma);
        %         velNTmp = sqrt(velN^2+velE^2+velD^2)*cosd(hdg)*cosd(gamma);
        velE    = velRef*sind(hdg)*cosd(gamma);
        %         velN    = velNTmp;
        velDOld = velD;
        %         velD    = - sqrt(velN^2+velE^2+velD^2)*sind(gamma);
        velD    = - velRef*sind(gamma);
        
        % Rate of change in geodetic coordinates [deg/s]
        [latLength, lonLength] = getLengthOfADegree(outerTraj(ii, 1, 1), outerTraj(ii, 1, 3));
        velLatOld = velLat;
        velLonOld = velLon;
        velLat  = velN / latLength;
        velLon  = velE / lonLength;
    
        outerTraj(ii, 1, 1) = outerTraj(ii-1, 1, 1) + (velLat + velLatOld)/2 * DELTA_T;
        outerTraj(ii, 1, 2) = outerTraj(ii-1, 1, 2) + (velLon + velLonOld)/2 * DELTA_T;
        outerTraj(ii, 1, 3) = outerTraj(ii-1, 1, 3) - (velD + velDOld)/2 * DELTA_T;
        
    else
        % [1] in Pilot reaction time
        % [3] wing-level flight

        % Update Attitude
        if ( abs( tand(rol) ) > 4 ) && ( abs(velX) > 0.1 )
            hdg = hdg + ( GRAV_ACC * 4 * sign( tand(rol) ) / velX * DELTA_T) * 180/pi; 
        elseif ( abs( tand(rol) ) <= 4 ) && ( abs(velX) > 0.1 )
            hdg = hdg + ( GRAV_ACC * tand(rol) / velX * DELTA_T) * 180/pi;
        else
            hdg = hdg + ( GRAV_ACC * tand(rol) / 0.1 * DELTA_T) * 180/pi;
        end
        
        % Update Flight Path Angle [deg]
        %        reaction time: 1.5 sec
        %        rate of change: 10 deg/s
        if ii < 4
            gamma = gamma;
        elseif gamma < -(10*DELTA_T)
            gamma = gamma + 10 * DELTA_T;
        elseif gamma < 0.0
            gamma = 0.0;
        else
            gamma = gamma;
        end
        
        % [2] Update Velocity
        velRef = sqrt(velN^2+velE^2+velD^2);
        velN   = velRef*cosd(hdg)*cosd(gamma);
        %         velNTmp = sqrt(velN^2+velE^2+velD^2)*cosd(hdg)*cosd(gamma);
        velE    = velRef*sind(hdg)*cosd(gamma);
        %         velN    = velNTmp;
        velDOld = velD;
        %         velD    = - sqrt(velN^2+velE^2+velD^2)*sind(gamma);
        velD    = - velRef*sind(gamma);
        % Rate of change in geodetic coordinates [deg/s]
        [latLength, lonLength] = getLengthOfADegree(outerTraj(ii, 1, 1), outerTraj(ii, 1, 3));

        velLatOld = velLat;
        velLonOld = velLon;
        velLat  = velN / latLength;
        velLon  = velE / lonLength;

        outerTraj(ii, 1, 1) = outerTraj(ii-1, 1, 1) + (velLat + velLatOld)/2 * DELTA_T;
        outerTraj(ii, 1, 2) = outerTraj(ii-1, 1, 2) + (velLon + velLonOld)/2 * DELTA_T;
        outerTraj(ii, 1, 3) = outerTraj(ii-1, 1, 3) - (velD + velDOld)/2 * DELTA_T;

    end

end


end


function innerTraj = PredictInnerTrajectory(eulerDeg, velXYZ, velNED, accXYZ, lla, uncHorPos, uncHdg)
%PredictInnerTrajectory 예측 비행경로를 생성
% 예측 비행경로 배열을 생성 (ALONG_TRACK_STEPS, 3, 3)
%       - Assumption 1: current aircraft state has uncertainty (3-sigma)
persistent DELTA_T
% Define constants
if isempty(DELTA_T)
    DELTA_T = 0.25;
end
GRAV_ACC = 9.8;
LIMIT_TAN_ROL = 4;
LIMIT_TURN_RATE = 0.222; % [rad/s]
RAD2DEG = 180/pi;
ALONG_TRACK_STEPS = 81;

% Constants for compare floating-point numbers
TOL = eps(0.5);

% Input Handling
rol = eulerDeg(1);
pit = eulerDeg(2);
hdg = eulerDeg(3);

velX = velXYZ(1);
velY = velXYZ(2);
velZ = velXYZ(3);

velN = velNED(1);
velE = velNED(2);
velD = velNED(3);

lat = lla(1);
lon = lla(2);
alt = lla(3);

% Flight path angle (deg)
gamma = atan2d(-velD, sqrt(velN^2+velE^2));

% Inertial speed
velXRef = sqrt( velX^2 + velY^2 + velZ^2 );

% Saturate inertial speed (since TAS is unavailable, use inertial spd)
velXRef = saturate(velXRef, 75.1089, 500);

% Saturate tangent of roll (max bank angle = 75.9630 deg)
tanRol = saturate(tand(rol), -LIMIT_TAN_ROL, LIMIT_TAN_ROL);

% Calculate turn rate [rad/s]
turnRate = GRAV_ACC * tanRol / velXRef;

% Saturate turn rate
turnRate = saturate(turnRate, -LIMIT_TURN_RATE, LIMIT_TURN_RATE);

% hdg = hdg + turnRate * RAD2DEG * 1;

% Saturate Uncertainties 
uncHorPos = saturate(uncHorPos, 35.56, 228.6);
uncHdg = saturate(uncHdg, 0.1, 5);

% Direction of uncertaitny addition
isTurnRight = (0 <= eulerDeg(1)) && (eulerDeg(1) < 180);
if isTurnRight
    uncPosRot = 90;
else
    uncHdg = (-1) * uncHdg;
    uncPosRot = -90;
end

% 3-sigma uncertainty for Heading
unitDirection = [velN; velE; 0] / norm([velN; velE; 0]);
unitDirection = rotateZ(3*uncHdg) * unitDirection;
hdg = hdg + 3*uncHdg;

% Update Velocity
velRef = sqrt(velN^2+velE^2+velD^2);
velN   = velRef*cosd(hdg)*cosd(gamma);
%         velNTmp = sqrt(velN^2+velE^2+velD^2)*cosd(hdg)*cosd(gamma);
velE    = velRef*sind(hdg)*cosd(gamma);
%         velN    = velNTmp;
velDOld = velD;
%         velD    = - sqrt(velN^2+velE^2+velD^2)*sind(gamma);
velD    = - velRef*sind(gamma);

% 3-sigma uncertainty for position
delPos = 3 * uncHorPos * rotateZ(uncPosRot + turnRate * RAD2DEG * 1) * unitDirection;

% Calculate the length of a degree
[latLength, lonLength] = getLengthOfADegree(lat, alt);
velLat  = velN / latLength;
velLon  = velE / lonLength;

% Initialize Predicted Trajectory
innerTraj = zeros(ALONG_TRACK_STEPS, 1, 3);

% Start point
innerTraj(1, 1, 1) = lat + delPos(1)/latLength;
innerTraj(1, 1, 2) = lon + delPos(2)/lonLength;
innerTraj(1, 1, 3) = alt;

for ii = 2:ALONG_TRACK_STEPS

    curT = (ii -1 ) * DELTA_T;

    if curT <= 20
        % [1] Update Attitude
        % Saturate tangent of roll (max bank angle = 75.9630 deg)
        tanRol = saturate(tand(rol), -LIMIT_TAN_ROL, LIMIT_TAN_ROL);

        % Saturate inertial speed (since TAS is unavailable, use inertial spd)
        velXRef = saturate(velXRef, 75.1089, 360.111);

        % Calculate turn rate [rad/s]
        turnRate = GRAV_ACC * tanRol / velXRef;

        % Saturate turn rate
        turnRate = saturate(turnRate, -LIMIT_TURN_RATE, LIMIT_TURN_RATE);

        % Update heading [deg]
        hdg = hdg + turnRate * DELTA_T * RAD2DEG;

        % Update Flight Path Angle [deg]
        %        reaction time: 1.5 sec
        %        rate of change: 10 deg/s
        if ii < 4
            gamma = gamma;
        elseif gamma < -(10*DELTA_T)
            gamma = gamma + 10 * DELTA_T;
        elseif gamma < 0.0
            gamma = 0.0;
        else
            gamma = gamma;
        end

        % [2] Update Velocity
        velRef = sqrt(velN^2+velE^2+velD^2);
        velN   = velRef*cosd(hdg)*cosd(gamma);
%         velNTmp = sqrt(velN^2+velE^2+velD^2)*cosd(hdg)*cosd(gamma);
        velE    = velRef*sind(hdg)*cosd(gamma);
%         velN    = velNTmp;
        velDOld = velD;
%         velD    = - sqrt(velN^2+velE^2+velD^2)*sind(gamma);
        velD    = - velRef*sind(gamma);
    else
        % [1] Update Attitude
        % No update

        % Update Flight Path Angle [deg]
        %        reaction time: 1.5 sec
        %        rate of change: 10 deg/s
        if ii < 4
            gamma = gamma;
        elseif gamma < -(10*DELTA_T)
            gamma = gamma + 10 * DELTA_T;
        elseif gamma < 0.0
            gamma = 0.0;
        else
            gamma = gamma;
        end

        % [2] Update Velocity
        velRef = sqrt(velN^2+velE^2+velD^2);
        velN   = velRef*cosd(hdg)*cosd(gamma);
%         velNTmp = sqrt(velN^2+velE^2+velD^2)*cosd(hdg)*cosd(gamma);
        velE    = velRef*sind(hdg)*cosd(gamma);
%         velN    = velNTmp;
        velDOld = velD;
%         velD    = - sqrt(velN^2+velE^2+velD^2)*sind(gamma);
        velD    = - velRef*sind(gamma);
    end

    % [3] Update Position
    % Rate of change in geodetic coordinates [deg/s]
    if abs( mod(ii, 41) - 0) < TOL % mod(ii, 41) == 0
        [latLength, lonLength] = getLengthOfADegree(innerTraj(ii, 1, 1), innerTraj(ii, 1, 3));
    end
    velLatOld = velLat;
    velLonOld = velLon;
    velLat  = velN / latLength;
    velLon  = velE / lonLength;

    innerTraj(ii, 1, 1) = innerTraj(ii-1, 1, 1) + (velLat + velLatOld)/2 * DELTA_T;
    innerTraj(ii, 1, 2) = innerTraj(ii-1, 1, 2) + (velLon + velLonOld)/2 * DELTA_T;
    innerTraj(ii, 1, 3) = innerTraj(ii-1, 1, 3) - (velD + velDOld)/2 * DELTA_T;

end

end

%% Internal Functions
function R = rotateZ(ang)
%ROTATEZ
% Rotation matrix for rotations around z-axis
% R = rotateZ(ang) creates a 3-by-3 matrix used to rotate a 3-by-1 vector
% or 3-by-N matrix of vectors around the z-axis by ang degrees.

%     R = zeros(3,3);
assert( isscalar(ang));
    
    R = [cosd(ang), -sind(ang), 0; ...
         sind(ang), cosd(ang),  0; ...
         0,         0,          1];
    
end

function out = saturate(in, min, max)
%saturate
% Saturates input value to [min, max]

if in < min
    out = min;
elseif in > max
    out = max;
else
    out = in;
end

end

function [vecNED, C_b_n] = convertBody2NED(eulerDeg, vecXYZ)
%CONVERTBODY2NED 동체 좌표계 벡터를 NED 좌표계 벡터로 변환
% Body2NED 변환 행렬 출력

% Initialize
vecNED = zeros(3,1);
C_n_b = zeros(3,3);
C_b_n = zeros(3,3);

% Trigonometirc functions
cRol = cosd(eulerDeg(1));
sRol = sind(eulerDeg(1));
cPit = cosd(eulerDeg(2));
sPit = sind(eulerDeg(2));
cYaw = cosd(eulerDeg(3));
sYaw = sind(eulerDeg(3));

% Transformation matrix
C_n_b   = [              cPit*cYaw,                cPit*sYaw,     -sPit;...
          sRol*sPit*cYaw-cRol*sYaw, sRol*sPit*sYaw+cRol*cYaw, sRol*cPit;...
          cRol*sPit*cYaw+sRol*sYaw, cRol*sPit*sYaw-sRol*cYaw, cRol*cPit];

C_b_n   = C_n_b';

% Transform
if size(vecXYZ, 1) == 3
    vecNED = C_b_n * vecXYZ;
elseif size(vecXYZ, 2) == 3
    vecNED = C_b_n * vecXYZ';
else
    vecNED = [0; 0; 0];
    warning('Invalid input vectors for <convertBody2NED>');
end

end