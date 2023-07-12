

%% [0]

% Define constants
GRAV_ACC = 9.8;

%% [1] INPUT/OUTPUIT INITIALZIATION & HANDLING 

% Initial conditions
lat = 35.3;                 % (deg)
lon = 128.3;                % (deg)
alt = 1000.0;               % (m)

rol = 00.0 * (pi/180);      % (rad)
pit = -850.0 * (pi/180);       % (rad)
hdg = 10.0 * (pi/180);      % (rad)

velU = 300.0;               % (m/s)
velV = 0.0;                 % (m/s)
velW = 0.0;                 % (m/s)

velN = 295.4423;            % (m/s)
velE = 52.0945;             % (m/s)
velD = 0.0;                 % (m/s)

% Get reference length of a geodetic degrees
[latLength, lonLength] = getLengthOfADegree(lat, alt);

% Define Vehicle struct
vehicle = struct();
vehicle.lat = lat;
vehicle.lon = lon;
vehicle.alt = alt;
vehicle.rol = rol;
vehicle.pit = pit;
vehicle.hdg = hdg;
vehicle.velU = velU;
vehicle.velV = velV;
vehicle.velW = velW;
vehicle.velN = velN;
vehicle.velE = velE;
vehicle.velD = velD;
vehicle.dt = 0.5;
vehicle.velLat = vehicle.velN / latLength; % [deg/s]
vehicle.velLon = vehicle.velE / lonLength; % [deg/s]
vehicle.latLength = latLength;
vehicle.lonLength = lonLength;
vehicle.accX = 0;
vehicle.accY = 0;
vehicle.accZ = 0;

% Simulation Time
dt = 0.04;
simT = 1;
time_arr = 0:dt:simT;
simLength = size(time_arr, 2);

hist_predTraj = zeros(81, 3, 3, simLength);
hist_nominalWTP = zeros(41, 1, simLength);

%% [2] Simulation
for ii = 1:simLength

    if ii == 1
    else
        vehicle = constantVelocity(vehicle);
    end

    eulerDeg = [vehicle.rol; vehicle.pit; vehicle.hdg] * 180/pi;
    velXYZ = [vehicle.velU; vehicle.velV; vehicle.velW];
    velNED = [vehicle.velN; vehicle.velE; vehicle.velD];
    accXYZ = [vehicle.accX; vehicle.accY; vehicle.accZ];
    lla = [vehicle.lat, vehicle.lon, vehicle.alt];

    predTraj = PredictTrajectory(eulerDeg, velXYZ, velNED, accXYZ, lla, 50, 1);
    [nominalWTP, valid, validity] = ScanMapGetWTP(predTraj, 1);

    hist_predTraj(:, :, :, ii) = predTraj;
    hist_nominalWTP(:, :, ii) = nominalWTP;

end

%% [3] Plot
figure;
wtp = hist_nominalWTP(:,1,1);
wtp_time = linspace(0,20,41);
h_wtp = plot(wtp_time, wtp, '-.b');
h_wtp.XDataSource = 'wtp_time';
h_wtp.YDataSource = 'wtp';
grid minor;
xlabel('Time (sec)'); ylabel('Height (m)');
axis([0, max(wtp_time), 0, max(hist_nominalWTP(:,1,:), [], 'all')])
for ii = 2:simLength
    wtp = hist_nominalWTP(:,1,ii);
    refreshdata
    pause(0.01);
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

function vehicle = constantVelocity(vehicle)
    % Model: CV

    % [1] Attitude update
    % → No update
    
    % [2] Velocity update
    % → No update
    
    % [3] Position update
    vehicle.lat = vehicle.lat + vehicle.velLat * vehicle.dt;
    vehicle.lon = vehicle.lon + vehicle.velLon * vehicle.dt;
    vehicle.alt = vehicle.alt + (-vehicle.velD) * vehicle.dt;    

end

function vehicle = constantTurn(vehicle)
    % Model: CT

    % [1] Attitude update
    if abs(GRAV_ACC * tan(vehicle.rol) ) > 4
        vehicle.hdg = vehicle.hdg + 4 * GRAV_ACC * sign( tan(vehicle.rol) ) / vehicle.velU * vehicle.dt;
    else
        vehicle.hdg = vehicle.hdg + 4 * GRAV_ACC * tan(vehicle.rol) / vehicle.velU * vehicle.dt;
    end
    
    % [2] Velocity update
    new_velN = sqrt(vehicle.velN^2 + vehicle.velE^2) * cos(vehicle.hdg);
    new_velE = sqrt(vehicle.velN^2 + vehicle.velE^2) * sin(vehicle.hdg);
    
    new_velLat = new_velN / vehicle.latLength; % [deg/s]
    new_velLon = new_velE / vehicle.lonLength; % [deg/s]
    
    % [3] Position update
    vehicle.lat = vehicle.lat + (vehicle.velLat + new_velLat)/2 * vehicle.dt;
    vehicle.lon = vehicle.lon + (vehicle.velLon + new_velLon)/2 * vehicle.dt;
    vehicle.alt = vehicle.alt + (-vehicle.velD) * vehicle.dt; 

end