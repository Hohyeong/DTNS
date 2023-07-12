%% ------------------------------------------------------------------
%  You can modify the values of the fields in BUS_FLCC_STATE_MATLABStruct
%  and evaluate this cell to create/update this structure
%  in the MATLAB base workspace.
% -------------------------------------------------------------------

initState = struct;
initState.latitude = 35.10;
initState.longitude = 127.70;
initState.altitude = 1000;
initState.roll = 0;
initState.pitch = 3*pi/180;
initState.yaw = 0;
initState.velX = 200;
initState.velY = 0;
initState.velZ = 0;
initState.omX = 0;
initState.omY = 0;
initState.omZ = 0;
initState.accX = 0;
initState.accY = 0;
initState.accZ = 0;
initState.omDotX = 0;
initState.omDotY = 0;
initState.omDotZ = 0;
initState.throttleN = 50000;
initState.elevatorDeg = 0;
initState.aileronDeg = 0;
initState.rudderDeg = 0;


version = -10;
% version = -3;
% version = -10;
switch version
% -------- TRN -----------------------------------------------------------%
    case 0
        % linear flight
        initState.latitude = 35.2;
        load('FTG_INPUT_SCENARIOS_linear.mat');
        initState.yaw = 0;
        initState.altitude = 2000;
    case 1
        % curve flight #2
        initState.latitude = 35.2;
        load('FTG_INPUT_SCENARIOS_curve.mat');
        initState.yaw = pi/2;
    case 2
        % curve flight #2
        initState.latitude = 35.8;
        load('FTG_INPUT_SCENARIOS_curve_v2.mat');
        initState.yaw = pi/2;
    case 101
        % Bank turn
        load('FTG_INPUT_SCENARIOS_bank.mat');
        initState.velX = 250;
% -------- TRN -----------------------------------------------------------%\
% -------- TRN_test ------------------------------------------------------%
    case 11
        initState.latitude = 37.235348;
        initState.longitude = 127.074891;
        load('FTG_INPUT_SCENARIOS_GPSinvalid_flat.mat');
        initState.yaw = -105*pi/180;
    case 12
        initState.latitude = 36.802863;
        initState.longitude = 127.07506;
        load('FTG_INPUT_SCENARIOS_GPSinvalid_smooth.mat');
        initState.yaw = 36*pi/180;
    case 13
        initState.latitude = 36.424351;
        initState.longitude = 127.142126;
        load('FTG_INPUT_SCENARIOS_GPSinvalid_moderate.mat');
        initState.yaw = 38*pi/180;
    case 14
        initState.latitude = 37.722089;
        initState.longitude = 127.07508;
        load('FTG_INPUT_SCENARIOS_GPSinvalid_rough.mat');
        initState.yaw = 45*pi/180;
    case 21
        initState.latitude = 36.965078;
        initState.longitude = 127.074929;
        load('FTG_INPUT_SCENARIOS_GPSvalid_flat.mat');
        initState.yaw = -38*pi/180;
    case 22
        initState.latitude = 36.532447;
        initState.longitude = 127.07489;
        load('FTG_INPUT_SCENARIOS_GPSvalid_smooth.mat');
        initState.yaw = -76*pi/180;
    case 23
        initState.latitude = 36.856776;
        initState.longitude = 127.209084;
        load('FTG_INPUT_SCENARIOS_GPSvalid_moderate.mat');
        initState.yaw = -165*pi/180;
    case 24
        initState.latitude = 35.581937;
        initState.longitude = 127.075048;
        load('FTG_INPUT_SCENARIOS_GPSvalid_rough.mat');
        initState.yaw = 154*pi/180;
        initState.altitude = 1700;
    case 31
        initState.latitude = 36.965073;
        initState.longitude = 127.074922;
        load('FTG_INPUT_SCENARIOS_GPSvalid_Pcode_flat.mat');
        initState.yaw = -43*pi/180;
    case 32
        initState.latitude = 37.85315;
        initState.longitude = 128.575053;
        load('FTG_INPUT_SCENARIOS_GPSvalid_Pcode_smooth.mat');
        initState.yaw = 27*pi/180;
    case 33
        initState.latitude = 35.311599;
        initState.longitude = 127.074912;
        load('FTG_INPUT_SCENARIOS_GPSvalid_Pcode_moderate.mat');
        initState.yaw = -127*pi/180;
    case 34
        initState.latitude = 35.527858;
        initState.longitude = 127.075026;
        load('FTG_INPUT_SCENARIOS_GPSvalid_Pcode_rough.mat');
        initState.yaw = 166*pi/180;
% -------- TRN_test ------------------------------------------------------%
% -------- TCAF ----------------------------------------------------------%
    case 1001
        % ATF
        load('FTG_INPUT_SCENARIOS_TCAFTest.mat');
        % ROUGH
%         load('TCAF_WHOLESIM_INPUT_SCENARIOS_ATF_N35_10E127_70.mat');
        initState.latitude = 37.39;
        initState.longitude = 128.5;
        initState.altitude = 1500;
        initState.velX = 250;
        initState.yaw = -145*pi/180;
    case 1101
        % AGCAS
%         load('TPF_WHOLESIM_INPUT_SCENARIOS_AGCAS_N35_10E127_70.mat');
        initState.yaw = -170*pi/180;
    case 1201
        % Feasibility Check
%         load('TPF_WHOLESIM_INPUT_SCENARIOS_FeasibilityCheck.mat');
        initState.altitude = 1000;
        initState.latitude = 35.16;
    case 1301
%         load('TCAF_WHOLESIM_INPUT_SCENARIOS_LevelFlight_N35_10E127_70.mat');
%         initState.latitude = 35.15;
        initState.latitude = 35.15;
        initState.altitude = 1000;
    case -22
        load('FTG_INPUT_SCENARIOS_TC022.mat');
        initState.latitude = 35.7;
        initState.longitude = 128;
        initState.altitude = 1300;
        initState.velX = 200;
        initState.yaw = 180*pi/180;
% -------- TPF -----------------------------------------------------------%
% -------- FTG -----------------------------------------------------------%
    case -1
        % FTG Test
        load('FTG_INPUT_SCENARIOS_FTGTest.mat');
        initState.latitude = 36.50;
        initState.longitude = 128.00;
        initState.altitude = 5000;
        initState.velX = 100;
        initState.yaw = 0.0;
        initState.pitch = 5 * pi/180;
    case -2
        % FTG Test on Bank-To-Turn
        load('FTG_INPUT_SCENARIOS_FTGTest2.mat');
        initState.latitude = 35.50;
        initState.longitude = 128.00;
        initState.altitude = 2500;
        initState.velX = 150;
        initState.yaw = 0.0;
        initState.pitch = 5 * pi/180;
    case -3
        % FTG Test on Bank-To-Turn Starting form ground
        load('FTG_INPUT_SCENARIOS_FTGTest3.mat');
        initState.latitude = 35.50;
        initState.longitude = 128.00;
        initState.altitude = 2500;
        initState.velX = 150;
        initState.yaw = 0.0;
        initState.pitch = 5 * pi/180;
    case 5
        % FTG Long
        load('FTG_INPUT_SCENARIOS_long.mat');
        initState.latitude = 35.30;
        initState.longitude = 127.70;
        initState.velX = 350;
        initState.altitude = 1000;
%         initState.velX = 320;
% -------- FTG -----------------------------------------------------------%
% -------- GROUND START --------------------------------------------------%
    case -10
        load('FTG_INPUT_SCENARIOS_GroundStart.mat');
        initState.latitude = 35.8;
        initState.longitude = 128.25;
        initState.velX = 0.0;
        initState.velY = 0.0;
        initState.velZ = 0.0;
        initState.velN = 0.0;
        initState.velE = 0.0;
        initState.velD = 0.0;
        initState.yaw = 0 * pi/180;
% -------- GROUND START --------------------------------------------------%
    otherwise
end
 

% Calculate (DCM) rotation matirx
% Define trigonometric functions
cRol = cos(initState.roll);
sRol = sin(initState.roll);
cPit = cos(initState.pitch);
sPit = sin(initState.pitch);
cYaw = cos(initState.yaw);
sYaw = sin(initState.yaw);

initState.C_b_n = ...
    [ cPit*cYaw, sRol*sPit*cYaw-cRol*sYaw, cRol*sPit*cYaw+sRol*sYaw;...
      cPit*sYaw, sRol*sPit*sYaw+cRol*cYaw, cRol*sPit*sYaw-sRol*cYaw;...
          -sPit,                sRol*cPit,                cRol*cPit];

velNED = initState.C_b_n ... 
    * [initState.velX; initState.velY; initState.velZ];
initState.velN = velNED(1);
initState.velE = velNED(2);
initState.velD = velNED(3);
windVelNED = zeros(3,1);
initState.windVelN = windVelNED(1);
initState.windVelE = windVelNED(2);
initState.windVelD = windVelNED(3);
[alpha, beta, TAS] = GetKinematicVariables(initState.C_b_n, velNED, windVelNED);
initState.alpha = alpha;
initState.beta = beta;
initState.Vt = TAS;
initState.f_ib_b = [0; 0; 9.8];


function [alpha, beta, TAS] = GetKinematicVariables(C_b_n, velNED, windVelNED)
%GetKinematicVariables

% Initialize
alpha   = 0.0;
beta    = 0.0;
TAS     = 0.0;

% Body to NED
% cRol = cos(attiEuler(1));
% sRol = sin(attiEuler(1));
% cPit = cos(attiEuler(2));
% sPit = sin(attiEuler(2));
% cYaw = cos(attiEuler(3));
% sYaw = sin(attiEuler(3));

% Coordinate Transformation
% NED to Body
C_n_b = C_b_n';
% C_n_b   = [              cPit*cYaw,                cPit*sYaw,     -sPit;...
%           sRol*sPit*cYaw-cRol*sYaw, sRol*sPit*sYaw+cRol*cYaw, sRol*cPit;...
%           cRol*sPit*cYaw+sRol*sYaw, cRol*sPit*sYaw-sRol*cYaw, cRol*cPit];

% Velocity relative to the wind
vRA_n = [velNED(1) - windVelNED(1); ...
    velNED(2) - windVelNED(2); ...
    velNED(3) - windVelNED(3)];
vRA_b = C_n_b * vRA_n;

TAS = norm(vRA_n);

% AOA
if (norm(vRA_b) < 0.5)
    alpha = 0.0;
else
    alpha = atan2(vRA_b(3), vRA_b(1));
end

% Sideslip Angle
if ( (abs(vRA_b(2)) < TAS) && (TAS > 0.5) )
    beta    = asin(vRA_b(2) / TAS);
else
    % Abnormal case
    beta    = 0.0;
end
end