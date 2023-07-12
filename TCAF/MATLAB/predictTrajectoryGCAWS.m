function predTraj = predictTrajectoryGCAWS(GCAWSParam, predTraj, eulerDeg, velXYZ, velNED, accXYZ, lla, DELTA_T)
%PREDICTTRAJECTORYGCAWS GCAWS 기준/회복궤적 생성
%
% == Input =========
% GCAWSParam = [reactT; gMax; rollRate; gRate; gammaMax]
% predTraj - size of [41, 12, 3]
% eulerDeg = [roll; pitch; heading] (deg)
% velXYZ = [celX; velY; velZ] (mps)
% velNED = [velN; velE; velD] (mps)
% accXYZ = [accX; accY; accZ] (mps2)
% lla = [latitude_deg; longitude_deg; altitude_msl_m]
% DELTA_T = 0.5
% == Output ========
% predTraj = [41, 12, 3]
% 

% Define Constants
GRAV_ACC = 9.8;
% ALONG_TRACK_STEPS = size(predTraj, 1);
ON_DEBUG = 1;

switch ON_DEBUG
    case 1
        ALONG_TRACK_STEPS = 41;
        predTraj = zeros(41, 12, 3);
        GCAWSParam = [1.5, 4.0, 120, 4.0, 75];
        
        DELTA_T = 0.5;
        eulerDeg = [00, 0, 10];
        velXYZ = [300, 0, 0];
%         [velNED, ~] = convertBody2NED(eulerDeg, velXYZ);
        velNED = [295.4675, 50.9427, 10.1885];
        accXYZ = [0, 0, 0];
        lla = [35.4, 128.6, 1000.0];

        accXYZ = [0, 3, 0];
        
    otherwise
        
end

% % If aircraft sink, 
% % PGCAS initial altitude should be lower for 0.5 sec x down velocity
% if velNED(3) >= 2
%     lla(3) = lla(3) - velNED(3) * 0.5;
% end

%% Trajectory Prediction

% <1> Non-Pull-Up Scenario
predTraj(:, 1, :) = projectTrajcetory(eulerDeg, velXYZ, velNED, accXYZ, lla, ALONG_TRACK_STEPS, DELTA_T);

% <2> Pull-Up Scenario
for jj = 2:12
    % [jj=2: 바로 Pull-up 기동; jj=3: 0.5초 후 Pull-up 기동; ... ; jj=k: (k-2)*0.5초 후 Pull-up 기동] 
    predTraj(:, jj, :) = pullUpTrajectory(GCAWSParam, eulerDeg, velXYZ, velNED, accXYZ, lla, ALONG_TRACK_STEPS, DELTA_T, jj-2);
end

if ON_DEBUG
    
    figure('WindowStyle', 'docked');
    hold on;
    plot(1:41, predTraj(:, 1, 3), '-r');
    for jj = 2:12
        plot(1:41, predTraj(:, jj, 3), '-b');
    end
    legend('projected', 'pull-up');
    grid on;
    
    figure('WindowStyle', 'docked');
    hold on;
    plot3(predTraj(:, 1, 2), predTraj(:, 1, 1), predTraj(:, 1, 3), '-r');
    for jj = 2:12
        plot3(predTraj(:, jj, 2), predTraj(:, jj, 1), predTraj(:, jj, 3), '-b');
    end
    legend('projected', 'pull-up');
    grid on;
    xlabel('Latitude [deg]');
    ylabel('Longitude [deg]');
    zlabel('Height [m]');

    figure('WindowStyle', 'docked');
    tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');
    nexttile;
    plot(predTraj(:, 2, 2) - predTraj(:, 1, 2));
    ylabel('Longitude (deg)'); xlabel('Index');
    grid minor;
    nexttile;
    plot(predTraj(:, 2, 1) - predTraj(:, 1, 1));
    ylabel('Latitude (deg)'); xlabel('Index');
    grid minor;
    
end

end

function pullUp = pullUpTrajectory(gcawsParam, eulerDeg, velXYZ, velNED, accXYZ, lla, ALONG_TRACK_STEPS, DELTA_T, pullUpTIdx)
%PULLUPTRAJECTORY
% Pull-up Trajectry 풀업 비행경로 예측

% Define constants
GRAV_ACC = 9.8;
RAD2DEG = 180/pi;
LIMIT_TAN_ROL = 4;
LIMIT_TURN_RATE = 0.222; % [rad/s]

assert (ALONG_TRACK_STEPS <= 81);

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

accNED = zeros(3,1);

lat = lla(1);
lon = lla(2);
alt = lla(3);

reactT = gcawsParam(1);     % [sec]
if velD > 0
    % 항공기 하강 시에Reaction Time을 0.5초 * Sin(비행경로 각도)만큼 증가
    gamma = atan2(velD, sqrt(velN^2 + velE^2));
    reactT = gcawsParam(1) + 0.5 * sin(gamma);
else
    reactT = gcawsParam(1); % [sec]
end

upGMax = gcawsParam(2);     % [g]
rollRate = gcawsParam(3);   % [deg/s]
gRate = gcawsParam(4);      % [g/s]
gammaMax = gcawsParam(5);   % [deg]

% Inertial speed
iSpd = sqrt( velX^2 + velY^2 + velZ^2 );
    
% Saturate inertial speed (since TAS is unavailable, use inertial spd) 
iSpd = saturate(iSpd, 0.1, 500);

% Calculate the length of a degree
[latLength, lonLength] = getLengthOfADegree(lat, alt);

velLat = velN / latLength;
velLon = velE / lonLength;

% Initialize Predicted Trajectory
pullUp = zeros(ALONG_TRACK_STEPS, 1, 3);
pullUp(1, 1, :) = lla;

% Time slices
% - Pull-up Timing (t1Idx = 1)
t1Idx = 1 + pullUpTIdx;
% - Reaction Tiem (t2Idx = 4.2)
t2Idx = t1Idx + (reactT/DELTA_T);
% - Roll Recovery (t3Idx = 4.3667)
t3Idx = t2Idx + abs( rol/rollRate )/DELTA_T;

fprintf("pullUpTIdx: %d, t1Idx: %f, t2Idx: %f, t3Idx: %f \r\n", pullUpTIdx, t1Idx, t2Idx, t3Idx);

for ii = 2:ALONG_TRACK_STEPS
    
    if ii <= t2Idx
        % [1] Projected & Pilot Reaction 
        % [1.1] Update Attitude

        % Saturate tangent of roll (max bank angle = 75.9630 deg)
        tanRol = saturate(tand(rol), -LIMIT_TAN_ROL, LIMIT_TAN_ROL);

        % Saturate inertial speed (since TAS is unavailable, use inertial spd)
        iSpd = saturate(iSpd, 0.1, 500);

        % Calculate turn rate [rad/s]
        turnRate = GRAV_ACC * tanRol / iSpd;

        % Saturate turn rate
        turnRate = saturate(turnRate, -LIMIT_TURN_RATE, LIMIT_TURN_RATE);

        % Update heading [deg]
        hdg = hdg + turnRate * DELTA_T * RAD2DEG;

        % [1.2] Update Velocity
        velNTmp = sqrt(velN^2+velE^2)*cosd(hdg);
        velE    = sqrt(velN^2+velE^2)*sind(hdg);
        velN    = velNTmp;
        velDOld = velD;

        if ii == 2
            % Update Velocity with Acceleration only when first time
            % Body2NED to get accNED
            [accNED, ~] = convertBody2NED(eulerDeg, accXYZ);
            
            % Update velocity (Acceleration is assumed to be diminishing)
            velN = velN + (accNED(1) + 0)/2 * DELTA_T;
            velE = velE + (accNED(2) + 0)/2 * DELTA_T;
            velD = velD + (accNED(3) + 0)/2 * DELTA_T;
        end

        % [1.3] Update Position
        % Rate of change in geodetic coordinates [deg/s]
        if mod(ii, 41) == 0
            [latLength, lonLength] = getLengthOfADegree(pullUp(ii, 1, 1), pullUp(ii, 1, 3));
        end
        velLatOld = velLat;
        velLonOld = velLon;
        velLat  = velN / latLength;
        velLon  = velE / lonLength;
        
        pullUp(ii, 1, 1) = pullUp(ii-1, 1, 1) + (velLat + velLatOld)/2 * DELTA_T;
        pullUp(ii, 1, 2) = pullUp(ii-1, 1, 2) + (velLon + velLonOld)/2 * DELTA_T;
        pullUp(ii, 1, 3) = pullUp(ii-1, 1, 3) - (velD + velDOld)/2 * DELTA_T;
        
    elseif (ii > t2Idx) && ( ii <= t3Idx )
        % [2] Roll recovery
        if ii == ceil(t2Idx)
            DELTA_T_rolRec = DELTA_T * (ii - t2Idx);
        else
            DELTA_T_rolRec = DELTA_T;
        end

        % [2.1] Update Attitude
        rolOld = rol;
        rol = rol - sign(rol) * rollRate * DELTA_T_rolRec;
        
        if sign(rol) ~= sign(rolOld) 
            rol = 0;
        end
        
        % Saturate tangent of roll (max bank angle = 75.9630 deg)
        tanRol = saturate(tand(rol), -LIMIT_TAN_ROL, LIMIT_TAN_ROL);
        tanRolOld = saturate(tand(rolOld), -LIMIT_TAN_ROL, LIMIT_TAN_ROL);

        % Saturate inertial speed (since TAS is unavailable, use inertial spd) 
        iSpd = saturate(iSpd, 0.1, 500);
        
        % Calculate turn rate [rad/s]
        turnRate = GRAV_ACC * tanRol / iSpd;
        turnRateOld = GRAV_ACC * tanRolOld / iSpd;

        % Saturate turn rate 
        turnRate = saturate(turnRate, -LIMIT_TURN_RATE, LIMIT_TURN_RATE);
        turnRateOld = saturate(turnRateOld, -LIMIT_TURN_RATE, LIMIT_TURN_RATE);
        
        % Update heading [deg]
        hdg = hdg + (turnRate + turnRateOld)/2 * DELTA_T_rolRec * RAD2DEG;
        
        % [2.2] Update Velocity
        velNTmp = sqrt(velN^2+velE^2)*cosd(hdg);
        velE    = sqrt(velN^2+velE^2)*sind(hdg);
        velN    = velNTmp;
        velDOld = velD;
        
        u = sqrt(velN^2 + velE^2);
        pit = -velD/u * RAD2DEG;
        
        % [2.3] Update Position
        % Rate of change in geodetic coordinates [deg/s]
        if mod(ii, 41) == 0
            [latLength, lonLength] = getLengthOfADegree(pullUp(ii, 1, 1), pullUp(ii, 1, 3));
        end
        velLatOld = velLat;
        velLonOld = velLon;
        velLat  = velN / latLength;
        velLon  = velE / lonLength;
        
        pullUp(ii, 1, 1) = pullUp(ii-1, 1, 1) + (velLat + velLatOld)/2 * DELTA_T;
        pullUp(ii, 1, 2) = pullUp(ii-1, 1, 2) + (velLon + velLonOld)/2 * DELTA_T;
        pullUp(ii, 1, 3) = pullUp(ii-1, 1, 3) - (velD + velDOld)/2 * DELTA_T;

    elseif ii > t3Idx
        % [3&4] Pull-up Maneuver
        rolOld = rol;
        rol = 0;
        
        if ii == ceil(t3Idx)
            DELTA_T_pullUp = DELTA_T * (ii - t3Idx);
        else
            DELTA_T_pullUp = DELTA_T;
        end

        
        gNED = [0; 0; GRAV_ACC];
        [~, C_b_n] = convertBody2NED(eulerDeg, [1; 1; 1]);
        C_n_b = C_b_n';
        
        gXYZ = C_n_b * gNED;
        
        if abs( GRAV_ACC/cosd(pit) ) > upGMax * GRAV_ACC
            accXYZ	=  [0; 0; -upGMax*GRAV_ACC] + gXYZ;
        else
            accXYZ = [0; 0; -abs( GRAV_ACC/cosd(pit) )] + gXYZ;
        end
        
        az = accXYZ(3);
        pit = pit - ( (az + GRAV_ACC*cosd(pit) ) / velX) * RAD2DEG * DELTA_T_pullUp;
        
        velDOld = velD;
        if abs(pit) < gammaMax
        % third step - g on set rate with pull up G
            
            az = az - gRate*GRAV_ACC*DELTA_T_pullUp;
            az = saturate(az, -upGMax*GRAV_ACC, upGMax*GRAV_ACC);
            aTmp = [0, 0, az]';
            [accNED, ~] = convertBody2NED([rol, pit, hdg], aTmp);

            velN = velN + accNED(1)*DELTA_T_pullUp;
            velE = velE + accNED(2)*DELTA_T_pullUp;
            velD = velD + accNED(3)*DELTA_T_pullUp;
        
        else
        % fourth step - Maximum flight path angle rising
            pit = sign(pit) * gammaMax;
        end
        
        % [3.3] Update Position
        % Rate of change in geodetic coordinates [deg/s]
        if mod(ii, 41) == 0
            [latLength, lonLength] = getLengthOfADegree(pullUp(ii, 1, 1), pullUp(ii, 1, 3));
        end
        velLatOld = velLat;
        velLonOld = velLon;
        velLat  = velN / latLength;
        velLon  = velE / lonLength;
        
        pullUp(ii, 1, 1) = pullUp(ii-1, 1, 1) + (velLat + velLatOld)/2 * DELTA_T;
        pullUp(ii, 1, 2) = pullUp(ii-1, 1, 2) + (velLon + velLonOld)/2 * DELTA_T;
        pullUp(ii, 1, 3) = pullUp(ii-1, 1, 3) - (velD + velDOld)/2 * DELTA_T;
        
    end
    
end

end


function nominal = projectTrajcetory(eulerDeg, velXYZ, velNED, accXYZ, lla, ALONG_TRACK_STEPS, DELTA_T)
%PROJECTTRAJECTORY
% Nominal Trajcetory 명목 궤적 예측
%       - Assumption 1: current aircraft sate is 100% reliable
%       - Assumption 2: aircraft maneuver is maintained

% Define constants
GRAV_ACC = 9.8;
LIMIT_TAN_ROL = 4;
LIMIT_TURN_RATE = 0.222; % [rad/s]
RAD2DEG = 180/pi;

assert (ALONG_TRACK_STEPS <= 81);

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

% Inertial speed
iSpd = sqrt( velX^2 + velY^2 + velZ^2 );
    
% Saturate inertial speed (since TAS is unavailable, use inertial spd) 
iSpd = saturate(iSpd, 0.1, 500);

% Calculate the length of a degree
[latLength, lonLength] = getLengthOfADegree(lat, alt);

velLat = velN / latLength;
velLon = velE / lonLength;

% Initialize Predicted Trajectory
nominal = zeros(ALONG_TRACK_STEPS, 1, 3);

% Start point
nominal(1, 1, 1) = lat;
nominal(1, 1, 2) = lon;
nominal(1, 1, 3) = alt;

for ii = 2:ALONG_TRACK_STEPS
    
    % [1] Update Attitude
    
    % Saturate tangent of roll (max bank angle = 75.9630 deg)
    tanRol = saturate(tand(rol), -LIMIT_TAN_ROL, LIMIT_TAN_ROL);
    
    % Saturate inertial speed (since TAS is unavailable, use inertial spd) 
    iSpd = saturate(iSpd, 0.1, 500);
    
    % Calculate turn rate [rad/s]
    turnRate = GRAV_ACC * tanRol / iSpd;
    
    % Saturate turn rate 
    turnRate = saturate(turnRate, -LIMIT_TURN_RATE, LIMIT_TURN_RATE);
        
    % Update heading [deg]
    hdg = hdg + turnRate * DELTA_T * RAD2DEG;
        
    % [2] Update Velocity
    velNTmp = sqrt(velN^2+velE^2)*cosd(hdg);
    velE    = sqrt(velN^2+velE^2)*sind(hdg);
    velN    = velNTmp;
    velDOld = velD;

    % [3] Update Position
    % Rate of change in geodetic coordinates [deg/s]
    if mod(ii, 41) == 0
        [latLength, lonLength] = getLengthOfADegree(nominal(ii, 1, 1), nominal(ii, 1, 3));
    end
    velLatOld = velLat;
    velLonOld = velLon;
    velLat  = velN / latLength;
    velLon  = velE / lonLength;
    
    nominal(ii, 1, 1) = nominal(ii-1, 1, 1) + (velLat + velLatOld)/2 * DELTA_T;
    nominal(ii, 1, 2) = nominal(ii-1, 1, 2) + (velLon + velLonOld)/2 * DELTA_T;
    nominal(ii, 1, 3) = nominal(ii-1, 1, 3) - (velD + velDOld)/2 * DELTA_T;

end

end

function out = saturate(in, min, max)
%SATURATE
% Saturates input value to [min, max]
%
% == Inputs =========
% in                - input value
% min               - minimum value
% max               - maximum value
% == Outputs ========
% out               - saturated output value
% ===================

if in < min
    out = min;
elseif in > max
    out = max;
else
    out = in;
end


end

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