%%TCAF_WHOLESIM_AFTERPLOT
% TCAF_WHOLESIM 시뮬레이션 종료 후 로그 데이터를 출력한다.
%
global ON_SAVE TIMESTR HEADER

ON_SAVE = true;

TIMESTR = datestr(now, 'yyyy-mm-dd_HHMMSS');
HEADER = ['logs\', TIMESTR];
mkdir(HEADER);

%% IDENTIFY

dataName = logsout.getElementNames;
numData = length(dataName);
for ii = 1:numData
    if strcmp(dataName{ii}, 'OUT_TRUE_STATE')
        idx_TRUE = ii;
    elseif strcmp(dataName{ii}, 'IN_TRNInput')
        idx_TRNInput = ii;
    elseif strcmp(dataName{ii}, 'OUT_TRNOutput')
        idx_TRNOutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_TRNDatapump')
        idx_TRNDatapump = ii;
    elseif strcmp(dataName{ii}, 'OUT_WTPGOutput')
        idx_WTPGOutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_WTPGDatapump')
        idx_WTPGDatapump = ii;
    elseif strcmp(dataName{ii}, 'IN_GCAWSInput')
        idx_GCAWSInput = ii;
    elseif strcmp(dataName{ii}, 'OUT_GCAWSOutput')
        idx_GCAWSOutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_GCAWSDatapump')
        idx_GCAWSDatapump = ii;
    elseif strcmp(dataName{ii}, 'IN_DTFInput')
        idx_DTFInput = ii;
    elseif strcmp(dataName{ii}, 'OUT_DTFOutput')
        idx_DTFOutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_DTFDatapump')
        idx_DTFDatapump = ii;
    elseif strcmp(dataName{ii}, 'Roughness')
        idx_R = ii;
    elseif strcmp(dataName{ii}, 'TRNEst_hDEM')
        idx_hDEM = ii;
    end
end


%% CALCULATION

%% PLOT

% % Plot TRN Data
plotTRNData(logsout, idx_TRUE, idx_TRNInput, idx_TRNOutput, idx_TRNDatapump, idx_R, idx_hDEM);

% % Plot DTF Data
plotDTFData(logsout, idx_WTPGOutput, idx_TRUE, idx_DTFDatapump, idx_DTFOutput, idx_TRNOutput, idx_WTPGDatapump, idx_DTFInput);

% % Plot GCAWS Data
plotGCAWSData(logsout, idx_GCAWSInput, idx_GCAWSOutput, idx_TRNOutput, idx_WTPGDatapump, idx_WTPGOutput, idx_GCAWSDatapump);

% -------------------------------------------------------------------------
%% INTERNAL FUNCTIONS -----------------------------------------------------
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

% -------------------------------------------------------------------------
%% PLOT FUINCTIONS --------------------------------------------------------
function plotTrueData(logsout, idx_TRUE)
global ON_SAVE TIMESTR HEADER

% /* True Position  */
posFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_TRUE}.Values.latitude);
ylabel('Latitude (deg)');
nexttile
plot(logsout{idx_TRUE}.Values.longitude);
ylabel('Longitude (deg)');
nexttile
plot(logsout{idx_TRUE}.Values.altitude);
ylabel('Altitude (m)');

% /* True Attitude  */
attFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_TRUE}.Values.roll * 180/pi);
ylabel('Roll (deg)');
nexttile
plot(logsout{idx_TRUE}.Values.pitch * 180/pi);
ylabel('Pitch (deg)');
nexttile
plot(logsout{idx_TRUE}.Values.yaw * 180/pi);
ylabel('Yaw (deg)');

% /* True Body Velocity  */
velFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_TRUE}.Values.velX);
ylabel('Body Velocity X (m/s)')
nexttile
plot(logsout{idx_TRUE}.Values.velY);
ylabel('Body Velocity Y (m/s)')
nexttile
plot(logsout{idx_TRUE}.Values.velZ);
ylabel('Body Velocity Z (m/s)')

% /* True Body Acceleration  */
accFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_TRUE}.Values.accX);
ylabel('Body Acceleration X (m/ss)');
nexttile
plot(logsout{idx_TRUE}.Values.accY);
ylabel('Body Acceleration Y (m/ss)');
nexttile
plot(logsout{idx_TRUE}.Values.accZ);
ylabel('Body Acceleration Z (m/ss)');

% /* Actuators  */
actFig = figure('WindowStyle', 'docked');
tiledlayout(4,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_TRUE}.Values.throttleN);
ylabel('Throttle (N)')
nexttile
plot(logsout{idx_TRUE}.Values.elevatorDeg);
ylabel('Elevator (deg)')
nexttile
plot(logsout{idx_TRUE}.Values.aileronDeg);
ylabel('Aileron (deg)')
nexttile
plot(logsout{idx_TRUE}.Values.rudderDeg);
ylabel('Rudder (deg)')


if ON_SAVE
    saveas(posFig, [HEADER, '\position(true)_', TIMESTR, '.fig']);
    saveas(attFig, [HEADER, '\attitude(true)_', TIMESTR, '.fig']);
    saveas(velFig, [HEADER, '\velocity(true)_', TIMESTR, '.fig']);
    saveas(accFig, [HEADER, '\acceleration(true)_', TIMESTR, '.fig']);
    saveas(actFig, [HEADER, '\actuator(true)_', TIMESTR, '.fig']);

    saveas(posFig, [HEADER, '\position(true)_', TIMESTR, '.png']);
    saveas(attFig, [HEADER, '\attitude(true)_', TIMESTR, '.png']);
    saveas(velFig, [HEADER, '\velocity(true)_', TIMESTR, '.png']);
    saveas(accFig, [HEADER, '\acceleration(true)_', TIMESTR, '.png']);
    saveas(actFig, [HEADER, '\actuator(true)_', TIMESTR, '.png']);
end

end

function plotTRNData(logsout, idx_TRUE, idx_TRNInput, idx_TRNOutput, idx_TRNDatapump, idx_R, idx_hDEM)
global ON_SAVE TIMESTR HEADER

% Time
trn_time = logsout{idx_TRNOutput}.Values.TRNTimeTag.Time;

% Errors
trn_LatEstErr = logsout{idx_TRNOutput}.Values.TRNEstOfLatitude.Data ...
    - logsout{idx_TRUE}.Values.latitude.Data;
trn_LonEstErr = logsout{idx_TRNOutput}.Values.TRNEstOfLongitude.Data ...
    - logsout{idx_TRUE}.Values.longitude.Data;
trn_AltEstErr = logsout{idx_TRNOutput}.Values.Longitude_INS.Data ...
    - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data ...
    - logsout{idx_TRUE}.Values.altitude.Data;
[latLength, lonLength] = getLengthOfADegree(...
    median(logsout{idx_TRNOutput}.Values.TRNEstOfLatitude.Data), ...
    median(logsout{idx_TRNOutput}.Values.Longitude_INS.Data));

trn_HorPosErr = sqrt( (trn_LatEstErr * latLength).^2 + ...
    (trn_LonEstErr * lonLength).^2 );

% RALTInvalidFlag
j = 1;

for i = 1 : length(trn_time)
    if logsout{idx_TRNInput}.Values.RALTValid.Data(i) == 0
        RALTinValid_Flag(j) = logsout{idx_TRNInput}.Values.RALTValid.Time(i);
        j = j + 1;
    end
end

if ~exist('RALTinValid_Flag', 'var')
    RALTinValid_Flag = 0;
end


% [1] TRN Position error + Roughness + Terrain height (4 by 1)
TRN_EstPosErrFig = figure('WindowStyle', 'docked');
tiledlayout(4,1, 'TileSpacing', 'tight', 'Padding', 'compact');
% [1 - (1,1)]
nexttile
plot(trn_time, trn_HorPosErr, '-b');
title('TRN Position Error [Horizon]')
xlabel('Time (sec)'); ylabel('Horizon (m)'); grid on;
yyaxis right
plot(trn_time, logsout{idx_TRNOutput}.Values.NavigationMode.Data);
ylabel('Navigation Mode')
yticks([0 1]); ylim([-0.2 1.2])
% [1 - (2,1)]
nexttile
plot(trn_time, trn_AltEstErr)
title('TRN Position Error [Vertical]')
xlabel('Time (sec)'), ylabel('Vertical (m)')
grid on
yyaxis right
plot(trn_time, logsout{idx_TRNOutput}.Values.NavigationMode.Data)
ylabel('Navigation Mode')
yticks([0 1]); ylim([-0.2 1.2])
% [1 - (3,1)]
nexttile
plot(trn_time, logsout{idx_R}.Values.Data);
title('Roughness')
xlabel('Time (sec)'), ylabel('[%]')
grid on
ylim([-5 70]);
% [1 - (4,1)]
nexttile
plot(trn_time, logsout{idx_hDEM}.Values.Data);
title('Terrain height (@TRN Est)')
xlabel('Time (sec)'), ylabel('(m)')
grid on
yyaxis right
plot(trn_time, logsout{idx_TRNInput}.Values.RALTValid.Data);
ylabel('RALT Valid');
yticks([0 1]); ylim([-0.2 1.2])


% 2 Sensor Flag (4 by 1)
TRN_SensorFig = figure('WindowStyle', 'docked');
tiledlayout(4,1, 'TileSpacing', 'tight', 'Padding', 'compact');
% [2 - (1,1)]
nexttile
hold on;
plot(trn_time, logsout{idx_TRNInput}.Values.GPSDataNotValid.Data, '-bo', 'DisplayName', 'GPSDataNotValid', 'MarkerIndices', 1:250:length(trn_time));
plot(trn_time, logsout{idx_TRNInput}.Values.INSDataNotValid.Data, '-g*', 'DisplayName', 'INSDataNotValid',  'MarkerIndices', 76:250:length(trn_time));
plot(trn_time, logsout{idx_TRNInput}.Values.RALTValid.Data, '-mx', 'DisplayName', 'RALTValid', 'MarkerIndices', 151:250:length(trn_time));
hold off;
title('Sensor Validity')
xlabel('Time (sec)'), ylabel('Flag')
legend('Location', 'best', 'Orientation', 'horizontal');
grid on;
ylim([-0.2, 2.0]);
% [2 - (2,1)]
nexttile
hold on;
plot(trn_time, logsout{idx_TRNInput}.Values.P_code.Data, '-bo', 'DisplayName', 'GPS P code', 'MarkerIndices', 1:250:length(trn_time));
plot(trn_time, logsout{idx_TRNInput}.Values.C_code.Data, '-g*', 'DisplayName', 'GPS CA code',  'MarkerIndices', 76:250:length(trn_time));
plot(trn_time, logsout{idx_TRNInput}.Values.FigureofMerit.Data, '-mx', 'DisplayName', 'GPS FOM', 'MarkerIndices', 151:250:length(trn_time));
hold off;
title('GPS Status')
xlabel('Time (sec)'), ylabel('Value')
legend('Location', 'best', 'Orientation', 'horizontal');
grid on;
ylim([-0.2, 9.0]);
% [2 - (3,1)]
nexttile
plot(trn_time, logsout{idx_TRNDatapump}.Values.MAD.Data, 'DisplayName', 'MAD');
title('TRN MAD');
xlabel('Time (sec)'), ylabel('MAD (m)')
legend('Location', 'best', 'Orientation', 'horizontal');
grid on;
% [2 - (4,1)]
nexttile
hold on;
plot(trn_time, logsout{idx_TRNDatapump}.Values.profileCnt.Data, 'DisplayName', 'Profile Count');
ylim([-1 25]);
yyaxis right
plot(trn_time, logsout{idx_TRNOutput}.Values.AircraftPositionCorrectionValid.Data, 'DisplayName', 'TRN Correction Valid');
title('TRN Profile Count');
xlabel('Time (sec)'), ylabel('count')
legend('Location', 'best', 'Orientation', 'horizontal');
grid on;
yticks([0 1]); ylim([-0.2 1.2])

% 3 Position (6 by 1)
TRN_PosFig = figure('WindowStyle', 'docked');
tiledlayout(6,1, 'TileSpacing', 'tight', 'Padding', 'compact');
% [3 - (1,1)]
nexttile
hold on;
plot(trn_time, logsout{idx_TRUE}.Values.latitude.Data, '-m', 'DisplayName', 'True');
plot(trn_time, logsout{idx_TRNInput}.Values.Latitude_INS.Data, '-g', 'DisplayName', 'INS');
plot(trn_time, logsout{idx_TRNInput}.Values.Latitude_GPS.Data, '-b', 'DisplayName', 'GPS');
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNEstOfLatitude.Data, '-r', 'DisplayName', 'TRN');
xlabel('Time (sec)'), ylabel('Latitude (deg)')
yyaxis right
xline(RALTinValid_Flag, 'r-', 'Alpha', 0.2, 'DisplayName', 'RALT Invalid')
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1]);
legend('Location', 'Best', 'Orientation', 'horizontal');
title('Latitude');
grid on;
% [3 - (2,1)]
nexttile
hold on;
plot(trn_time, logsout{idx_TRNInput}.Values.Latitude_INS.Data - logsout{idx_TRUE}.Values.latitude.Data, '-g', 'DisplayName', 'INS');
plot(trn_time, logsout{idx_TRNInput}.Values.Latitude_GPS.Data - logsout{idx_TRUE}.Values.latitude.Data, '-b', 'DisplayName', 'GPS');
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNEstOfLatitude.Data - logsout{idx_TRUE}.Values.latitude.Data, '-r', 'DisplayName', 'TRN');
xlabel('Time (sec)'), ylabel('Latitude (deg)')
title('Latitude Error');
legend('Location', 'Best', 'Orientation', 'horizontal');
grid on;
axisLim = axis; ylim([-1.2*max(abs(axisLim(3:4))) +1.2*max(abs(axisLim(3:4)))]);
yyaxis right
xline(RALTinValid_Flag, 'r-', 'Alpha', 0.2, 'DisplayName', 'RALT Invalid')
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1]);
% [3 - (3,1)]
nexttile
hold on;
plot(trn_time, logsout{idx_TRUE}.Values.longitude.Data, '-m', 'DisplayName', 'True');
plot(trn_time, logsout{idx_TRNInput}.Values.Longitude_INS.Data, '-g', 'DisplayName', 'INS');
plot(trn_time, logsout{idx_TRNInput}.Values.Longitude_GPS.Data, '-b', 'DisplayName', 'GPS');
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNEstOfLongitude.Data, '-r', 'DisplayName', 'TRN');
xlabel('Time (sec)'), ylabel('Longitude (deg)')
yyaxis right
xline(RALTinValid_Flag, 'r-', 'Alpha', 0.2, 'DisplayName', 'RALT Invalid')
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1]);
legend('Location', 'Best', 'Orientation', 'horizontal');
title('Longitude');
grid on;
% [3 - (4,1)]
nexttile
hold on;
plot(trn_time, logsout{idx_TRNInput}.Values.Longitude_INS.Data - logsout{idx_TRUE}.Values.longitude.Data, '-g', 'DisplayName', 'INS');
plot(trn_time, logsout{idx_TRNInput}.Values.Longitude_GPS.Data - logsout{idx_TRUE}.Values.longitude.Data, '-b', 'DisplayName', 'GPS');
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNEstOfLongitude.Data - logsout{idx_TRUE}.Values.longitude.Data, '-r', 'DisplayName', 'TRN');
xlabel('Time (sec)'), ylabel('Longitude (deg)')
title('Longitude Error');
legend('Location', 'Best', 'Orientation', 'horizontal');
grid on;
axisLim = axis; ylim([-1.2*max(abs(axisLim(3:4))) +1.2*max(abs(axisLim(3:4)))]);
yyaxis right
xline(RALTinValid_Flag, 'r-', 'Alpha', 0.2, 'DisplayName', 'RALT Invalid')
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1]);
% [3 - (5,1)]
nexttile
hold on;
plot(trn_time, logsout{idx_TRUE}.Values.altitude.Data, '-m', 'DisplayName', 'True');
plot(trn_time, logsout{idx_TRNInput}.Values.Altitude_INS.Data, '-g', 'DisplayName', 'INS');
plot(trn_time, logsout{idx_TRNInput}.Values.Altitude_GPS.Data, '-b', 'DisplayName', 'GPS');
plot(trn_time, logsout{idx_TRNInput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data, '-r', 'DisplayName', 'TRN');
xlabel('Time (sec)'), ylabel('Altitude (m)')
yyaxis right
xline(RALTinValid_Flag, 'r-', 'Alpha', 0.2, 'DisplayName', 'RALT Invalid')
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1]);
legend('Location', 'Best', 'Orientation', 'horizontal');
title('Altitude');
grid on;
% [3 - (6,1)]
nexttile
hold on;
plot(trn_time, logsout{idx_TRNInput}.Values.Altitude_INS.Data - logsout{idx_TRUE}.Values.altitude.Data, '-g', 'DisplayName', 'INS');
plot(trn_time, logsout{idx_TRNInput}.Values.Altitude_GPS.Data - logsout{idx_TRUE}.Values.altitude.Data, '-b', 'DisplayName', 'GPS');
plot(trn_time, logsout{idx_TRNInput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data - logsout{idx_TRUE}.Values.altitude.Data, '-r', 'DisplayName', 'TRN');
xlabel('Time (sec)'), ylabel('Altitude (m)')
axisLim = axis; ylim([-1.2*max(abs(axisLim(3:4))) +1.2*max(abs(axisLim(3:4)))]);
yyaxis right
xline(RALTinValid_Flag, 'r-', 'Alpha', 0.2, 'DisplayName', 'RALT Invalid')
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1]);
legend('Location', 'Best', 'Orientation', 'horizontal');
title('Altitude Error');
grid on;

% 4 Time Tags (4 by 1)
% TRN_TimeTagFig = figure('WindowStyle', 'docked');
% tiledlayout(4,1, 'TileSpacing', 'tight', 'Padding', 'compact');
% % [4 - (1,1)]
% nexttile
% hold on;
% plot(trn_time, logsout{idx_TRNInput}.Values.INSTimeTag.Data, '-g', 'DisplayName', 'INS');
% plot(trn_time, logsout{idx_TRNInput}.Values.GPStimeValid.Data, '-b', 'DisplayName', 'GPS');
% plot(trn_time, logsout{idx_TRNInput}.Values.RALTTimeTag.Data, '-r', 'DisplayName', 'RALT');
% plot(trn_time, logsout{idx_TRNOutput}.Values.TRNTimeTag.Data, '-c', 'DisplayName', 'TRN');
% xlabel('Time (sec)'), ylabel('TimeTag')
% legend('Location', 'Best', 'Orientation', 'horizontal');
% [4 - (2,1)]
% nexttile
% hold on;
% plot(trn_time, mod(logsout{idx_TRNInput}.Values.GPStimeValid.Data - logsout{idx_TRNInput}.Values.INSTimeTag.Data, 15626), '-b', 'DisplayName', 'GPS');
% plot(trn_time, mod(logsout{idx_TRNInput}.Values.RALTTimeTag.Data - logsout{idx_TRNInput}.Values.INSTimeTag.Data, 15626), '-r', 'DisplayName', 'RALT');
% plot(trn_time, mod(logsout{idx_TRNOutput}.Values.TRNTimeTag.Data - logsout{idx_TRNInput}.Values.INSTimeTag.Data, 15626), '-c', 'DisplayName', 'TRN');
% xlabel('Time (sec)'), ylabel('TimeTag Difference')
% legend('Location', 'Best', 'Orientation', 'horizontal');

if ON_SAVE
%     saveas(trnFig, [HEADER, '\TRNEstimate_', TIMESTR, '.fig']);
%     saveas(trnErrFig, [HEADER, '\TRNEstimateError_', TIMESTR, '.fig']);
% 
%     saveas(trnFig, [HEADER, '\TRNEstimate_', TIMESTR, '.png']);
%     saveas(trnErrFig, [HEADER, '\TRNEstimateError_', TIMESTR, '.png']);
end

end

function plotDTFData(logsout, idx_WTPGOutput, idx_TRUE, idx_DTFDatapump, idx_DTFOutput, idx_TRNOutput, idx_WTPGDatapump, idx_DTFInput)
global ON_SAVE TIMESTR HEADER

% Time
trn_time = logsout{idx_TRNOutput}.Values.TRNTimeTag.Time;
dtf_time = logsout{idx_DTFDatapump}.Values.referenceTrajectory.Time;

% MSD, SWH
dtf_msd = logsout{idx_DTFOutput}.Values.msdUsedWrap.Data;

% Worst-case Terrain profile
wtpg_wtp = logsout{idx_WTPGOutput}.Values.wtp.Data;

% violateTRN = (logsout{idx_WCPGOutput}.Values.WCPGWCP(:,1).Data + dbtf_msd) - (logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data + logsout{idx_TRNOutput}.Values.Altitude_INS.Data);
violateTrue = (wtpg_wtp(:,1) + dtf_msd) - logsout{idx_TRUE}.Values.altitude.Data;

% violateTRN = min(inf, max(0, violateTRN));
violateTrue = min(inf, max(0, violateTrue));

% violateTRNRatio = violateTRN ./ dbtf_msd * 100;
violateTrueRatio = violateTrue ./ dtf_msd * 100;

dtf_traj = zeros(size(dtf_time));

dtf_trajOneSec = zeros(size(dtf_time));
dtf_trajTwoSec = zeros(size(dtf_time));
dtf_trajThrSec = zeros(size(dtf_time));
dtf_trajFouSec = zeros(size(dtf_time));

for ii = 1:size(dtf_traj, 1)
    
    % 1초 전
    if  (ii > 25)
        idxOneSec = ii - 25;
        cursor = 3;
    else
        idxOneSec = ii;
        cursor = 1;
    end
    dtf_trajOneSec(ii) = logsout{idx_DTFDatapump}.Values.referenceTrajectory.Data(idxOneSec, cursor);
    
    % 2초 전
    if (ii > 50)
        idxTwoSec = ii - 50;
        cursor = 5;
    elseif (ii > 25)
        idxTwoSec = ii - 25;
        cursor = 3;
    else
        idxTwoSec = ii;
        cursor = 1;
    end
    dtf_trajTwoSec(ii) = logsout{idx_DTFDatapump}.Values.referenceTrajectory.Data(idxTwoSec, cursor);
    
    % 3초 전
    if (ii > 75)
        idxThrSec = ii - 75;
        cursor = 7;
    elseif (ii > 50)
        idxThrSec = ii - 50;
        cursor = 5;
    elseif (ii > 25)
        idxThrSec = ii - 25;
        cursor = 3;
    else
        idxThrSec = ii;
        cursor = 1;
    end
    dtf_trajThrSec(ii) = logsout{idx_DTFDatapump}.Values.referenceTrajectory.Data(idxThrSec, cursor);
    
    % 4초전
    if (ii > 100)
        idxFouSec = ii - 100;
        cursor = 9;
    elseif (ii > 75)
        idxFouSec = ii - 75;
        cursor = 7;
    elseif (ii > 50)
        idxFouSec = ii - 50;
        cursor = 5;
    elseif (ii > 25)
        idxFouSec = ii - 25;
        cursor = 3;
    else
        idxFouSec = ii;
        cursor = 1;
    end
    dtf_trajFouSec(ii) = logsout{idx_DTFDatapump}.Values.referenceTrajectory.Data(idxFouSec, cursor);
    
end

dtfTrajFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile
hold on;
plot(logsout{idx_TRUE}.Values.altitude, '-r');
plot(logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Time, logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data ,...
    '-b');
plot(dtf_time, wtpg_wtp(:,1)+dtf_msd, 'Color', [0.1 0.8 0.1], 'LineWidth', 3);

% plot(dbtf_time, dbtf_traj);
plot(dtf_time, dtf_trajOneSec, 'Color', '#960096');
plot(dtf_time, dtf_trajTwoSec, 'Color', '#A500A5');
plot(dtf_time, dtf_trajThrSec, 'Color', '#B400B4');
plot(dtf_time, dtf_trajFouSec, 'Color', '#CC00CC');
ylabel('Height (m)');
grid minor;
axis([dtf_time(1), dtf_time(end), min(wtpg_wtp(:,1)+dtf_msd)-100, max(dtf_trajOneSec)+100]);
legend('Flight Trajectory', 'TRN Alt Est.', 'Terrain + MSD', ...
    'DTF Ref. Traj.(before 1 sec)', 'DTF Ref. Traj.(before 2 sec)', 'DTF Ref. Traj.(before 3 sec)', 'DTF Ref. Traj.(before 4 sec)', ...
    'Location', 'best', 'FontSize', 7, 'orientation', 'horizontal');
nexttile
accX = logsout{idx_TRUE}.Values.accX.Data;
accY = logsout{idx_TRUE}.Values.accY.Data;
accZ = logsout{idx_TRUE}.Values.accZ.Data;
cRol = cos(logsout{idx_TRUE}.Values.roll.Data);
sRol = sin(logsout{idx_TRUE}.Values.roll.Data);
cPit = cos(logsout{idx_TRUE}.Values.pitch.Data);
sPit = sin(logsout{idx_TRUE}.Values.pitch.Data);
accVert = -sPit .* accX + sRol.*cPit .* accY + cRol.*cPit .* accZ;
hold on;
plot(dtf_time, logsout{idx_DTFOutput}.Values.verticalGCommand.Data*9.8);
plot(dtf_time, -1 * logsout{idx_TRUE}.Values.accZ.Data);
plot(dtf_time, -accVert);
legend('Normal Acceleration Command', 'Body Acceleration (-Z)', 'Vertical Acceleration (+Up)', 'FontSize', 7);
ylabel('Acceleration (m/s^2)');
grid minor;
axis([dtf_time(1), dtf_time(end), -15, 25]);
nexttile
hold on;
TRN_SD = logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data - wtpg_wtp(:,1);
TRUE_SD = logsout{idx_TRUE}.Values.altitude.Data - wtpg_wtp(:,1);
plot(trn_time, TRN_SD);
plot(trn_time, TRUE_SD);
plot(dtf_time, logsout{idx_DTFOutput}.Values.msdUsedWrap.Data);
axis([dtf_time(1), dtf_time(end), -50, max(TRN_SD)+100]);
ylabel('Seperation Distance (m)'); xlabel('Time (seconds)');
grid minor;
legend('Seperation distacne w.r.t. True', 'Seperation distacne w.r.t. TRN', 'MSD', 'FontSize', 7);

fprintf('DTF Statistics\r\n');
fprintf('Seperation Distance: mean = %f(m), standard deviation = %f(m)\r\n', mean(TRUE_SD(501:end)), std(TRUE_SD(501:end)));
fprintf('Acceleration Command: mean = %f(m/s^2), standard deviation = %f(m/s^2)\r\n', mean(abs(logsout{idx_DTFOutput}.Values.verticalGCommand.Data(501:end)*9.8)), std((logsout{idx_DTFOutput}.Values.verticalGCommand.Data(501:end)*9.8)));

% /* MSD Violation */
msdFig = figure('WindowStyle', 'docked');
tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
hold on;
title('MSD 침하율');
% plot(dbtf_time, violateTRNRatio, 'Marker', '.', 'Color', [72 155 123]/255);
plot(dtf_time, violateTrueRatio, 'Marker', '.','Color', [203 101 39]/255);
plot([dtf_time(1), dtf_time(end)], [25, 25], 'LineStyle', '--', 'Color', [0.8353, 0.2784, 0.5373]);
ylabel('Violation Ratio (%)');
legend('Violation w.r.t True Altitude', '25%');
hold off;
grid on;
axis([dtf_time(1), dtf_time(end), -1, 30]);
nexttile
hold on;
title('MSD 침하 고도')
% plot(dbtf_time, violateTRN, 'Marker', '.', 'Color', [72 155 123]/255);
plot(dtf_time, violateTrue, 'Marker', '.','Color', [203 101 39]/255);
plot(dtf_time, dtf_msd*0.25, 'LineStyle', '--', 'Color', [0.8353, 0.2784, 0.5373]);
ylabel('Violation Height (m)');
legend('Violation w.r.t True Altitude', '25%');
hold off;
grid on;
axis([dtf_time(1), dtf_time(end), -1, 60]);

% /* DBTF Command */
dtfFig = figure('WindowStyle', 'docked');
% Calculate Vertical Acceleration
accX = logsout{idx_TRUE}.Values.accX.Data;
accY = logsout{idx_TRUE}.Values.accY.Data;
accZ = logsout{idx_TRUE}.Values.accZ.Data;
cRol = cos(logsout{idx_TRUE}.Values.roll.Data);
sRol = sin(logsout{idx_TRUE}.Values.roll.Data);
cPit = cos(logsout{idx_TRUE}.Values.pitch.Data);
sPit = sin(logsout{idx_TRUE}.Values.pitch.Data);
accVert = -sPit .* accX + sRol.*cPit .* accY + cRol.*cPit .* accZ;
hold on;
plot(dtf_time, logsout{idx_DTFOutput}.Values.verticalGCommand.Data*9.8);
plot(dtf_time, -1 * logsout{idx_TRUE}.Values.accZ.Data);
plot(dtf_time, -accVert);
legend('Normal Acceleration Command', 'Body Acceleration (-Z)', 'Vertical Acceleration (+Up)');
ylabel('Acceleration (m/s^2)');
grid on;
axis([dtf_time(1), dtf_time(end), -30, 30]);

dtfOutFig = figure('WindowStyle', 'docked');
tiledlayout(7,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(dtf_time, logsout{idx_DTFOutput}.Values.valid.Data);
ylabel('DTF Valid Wraparound'); grid minor;
nexttile
plot(dtf_time, logsout{idx_DTFOutput}.Values.selectedWrap.Data);
ylabel('DTF Selected Wraparound'); grid minor;
nexttile
plot(dtf_time, logsout{idx_DTFOutput}.Values.rideQualityWrap.Data);
ylabel('Ride Quality Wraparound'); grid minor;
nexttile
plot(dtf_time, logsout{idx_DTFOutput}.Values.msdUsedWrap.Data);
ylabel('MSD Used Wraparound'); grid minor;
nexttile
plot(dtf_time, logsout{idx_DTFOutput}.Values.gLimitWarning.Data);
ylabel('G Limit Warning'); grid minor;
nexttile
plot(dtf_time, logsout{idx_DTFOutput}.Values.lowTFWarning.Data);
ylabel('Low TF Warning'); grid minor;
nexttile
plot(dtf_time, logsout{idx_DTFOutput}.Values.verticalGCommand.Data);
ylabel('Vertical G Command'); grid minor;
xlabel('Time (sec)');


% DTF Feasibility
dtfFeasibilityFig = figure('WindowStyle', 'docked');
tiledlayout(7,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.TRN_valid.Data, '-b');
ylabel('TRN Valid'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.AircraftPositionCorrectionValid.Data, '-b');
ylabel('TRN Correction Valid'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.AircraftOffTerrainMap.Data, '-b');
ylabel('TRN Off Terrain Map'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.NavigationMode.Data, '-b');
ylabel('TRN Navigation Mode'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
nexttile
plot(trn_time, logsout{idx_WTPGDatapump}.Values.isFeasible.Data, '-b');
ylabel('WTPG Feasibility'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
xlabel('Time (sec)');
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.TRN_VelocityX.Data, '-b');
ylabel('Speed (m/s)'); grid on; axis([trn_time(1), trn_time(end), 0, 450]);
xlabel('Time (sec)');
nexttile
plot(dtf_time, logsout{idx_DTFDatapump}.Values.isFeasible.Data, '-r');
ylabel('DTF Feasibility'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
xlabel('Time (sec)');


% DTF Basic Function Test
DTFBasicFig = figure('WindowStyle', 'docked');
tiledlayout(5,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
hold on;
plot(dtf_time, logsout{idx_DTFOutput}.Values.valid.Data, 'DisplayName', 'Validity');
plot(dtf_time, logsout{idx_DTFOutput}.Values.selectedWrap.Data, 'DisplayName', 'Selected Wrap');
ylabel('DTF Valid'); xlabel('Time (sec)');
grid on; axis([dtf_time(1), dtf_time(end), -0.1, 1.1]);
legend
nexttile
hold on;
plot(dtf_time, logsout{idx_DTFOutput}.Values.gLimitWarning.Data, 'DisplayName', 'G Limit');
plot(dtf_time, logsout{idx_DTFOutput}.Values.lowTFWarning.Data, 'DisplayName', 'Low TF');
ylabel('Warning'); xlabel('Time (sec)');
legend;
grid on; axis([dtf_time(1), dtf_time(end), -0.1, 1.1]);
nexttile
hold on;
plot(dtf_time, logsout{idx_DTFOutput}.Values.msdUsedWrap.Data, 'DisplayName', 'MSD Used Wrap');
plot(dtf_time, logsout{idx_DTFOutput}.Values.msdUsedWrap.Data*0.75, '--r', 'DisplayName', '75% MSD Used Wrap');
plot(dtf_time, logsout{idx_DTFInput}.Values.raltAlt.Data, 'DisplayName', 'RALT AGL');
ylabel('Height (m)'); xlabel('Time (sec)'); grid on;
axis([dtf_time(1), dtf_time(end), 0, max(logsout{idx_DTFInput}.Values.raltAlt.Data)+100]);
legend
nexttile
plot(dtf_time, logsout{idx_DTFOutput}.Values.rideQualityWrap.Data, 'DisplayName', 'Ride Quality');
ylabel('Ride Quality'); xlabel('Time (sec)'); grid on;
axis([dtf_time(1), dtf_time(end), 0, 4]);
legend
nexttile
plot(dtf_time, logsout{idx_DTFOutput}.Values.verticalGCommand.Data, 'DisplayName', 'G Commnad');
ylabel('G Command (g)'); xlabel('Time (sec)'); grid on;
axis([dtf_time(1), dtf_time(end), -4.0, 2.5]);
legend

if ON_SAVE
    saveas(dtfTrajFig, [HEADER, '\DTFTrajectories_', TIMESTR, '.fig']);
    saveas(msdFig, [HEADER, '\DTFMSD_', TIMESTR, '.fig']);
    saveas(dtfFig, [HEADER, '\DTFCommand_', TIMESTR, '.fig']);
    saveas(dtfOutFig, [HEADER, '\DTFOutput_', TIMESTR, '.fig']);
    saveas(dtfFeasibilityFig, [HEADER, '\DFTFeasibility_', TIMESTR, '.fig']);
    saveas(DTFBasicFig, [HEADER, '\DTFBasicTest_', TIMESTR, '.fig']);

    saveas(dtfTrajFig, [HEADER, '\DTFTrajectories_', TIMESTR, '.png']);
    saveas(msdFig, [HEADER, '\DTFMSD_', TIMESTR, '.png']);
    saveas(dtfFig, [HEADER, '\DTFCommand_', TIMESTR, '.png']);
    saveas(dtfOutFig, [HEADER, '\DTFOutput_', TIMESTR, '.png']);
    saveas(dtfFeasibilityFig, [HEADER, '\DFTFeasibility_', TIMESTR, '.png']);
    saveas(DTFBasicFig, [HEADER, '\DTFBasicTest_', TIMESTR, '.png']);
end

end


function plotGCAWSData(logsout, idx_GCAWSInput, idx_GCAWSOutput, idx_TRNOutput, idx_WTPGDatapump, idx_WTPGOutput, idx_GCAWSDatapump)
global ON_SAVE TIMESTR HEADER

% Time
trn_time = logsout{idx_TRNOutput}.Values.TRNTimeTag.Time;
gcaws_time = logsout{idx_GCAWSOutput}.Values.selectedWrap.Time;
wtpg_time = logsout{idx_WTPGOutput}.Values.wtp.Time;

% Worst-case Terrain profile
wtpg_wtp = logsout{idx_WTPGOutput}.Values.wtp.Data;

% SWH
gcaws_swh = logsout{idx_GCAWSOutput}.Values.swhUsedWrap.Data;

% GCAWS_1_ProfileFig
gcawsFig = figure('WindowStyle', 'docked');
tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
hold on;
plot(trn_time, (logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data), '-r');
% plot(trn_time, (logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data), '-b');
plot(wtpg_time, (wtpg_wtp(:,1)+gcaws_swh), 'Color', [0.1 0.8 0.1], 'LineWidth', 3);
plot(wtpg_time, wtpg_wtp(:,1), 'Color', [0.6 0.8 0.4], 'LineWidth', 3);
legend('Reference Flight Trajectory', 'Terrain(WTPG) + SWH', 'Terrain(WTPG)', 'Location', 'Best');
ylabel('Height (m)');
grid on;
axis([gcaws_time(1), gcaws_time(end), 0, max(logsout{idx_TRNOutput}.Values.Altitude_INS.Data)+100]);
nexttile
hold on;
plot(logsout{idx_GCAWSOutput}.Values.timeToGoToPullUp, 'r.');
% plot(dbtf_time, min(inf, max(0, 38.28 - dbtf_time)));
xlabel('Flight Time (sec)');
ylabel('Time Left (sec)');
axis([wtpg_time(1), wtpg_time(end), -1, 10]);
grid on;
title('Time To Pull-Up Left')
% legend('Time to Go to Pull-Up', 'Collision Left');

% GCAWS_2_Feasibility
gcawsFeasibilityFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
title('TRN Validity');
hold on;
plot(trn_time, logsout{idx_TRNOutput}.Values.TRN_valid.Data, '-bo', 'MarkerIndices',1:200:length(trn_time));
plot(trn_time, logsout{idx_TRNOutput}.Values.AircraftPositionCorrectionValid.Data, '-g*', 'MarkerIndices',51:200:length(trn_time));
plot(trn_time, logsout{idx_TRNOutput}.Values.AircraftOffTerrainMap.Data, '-mx', 'MarkerIndices',101:200:length(trn_time));
plot(trn_time, logsout{idx_TRNOutput}.Values.NavigationMode.Data, '-r+', 'MarkerIndices',151:200:length(trn_time));
ylabel('True/False'); grid on; axis([trn_time(1), trn_time(end), -0.2, 2.0]);
legend('TRN Valid', 'TRN Correction Valid', 'TRN Off Map', 'TRN Mode', 'Location', 'Best', 'Orientation', 'horizontal');
hold off;
nexttile
title('WTPG Validity');
hold on;
plot(trn_time, logsout{idx_WTPGDatapump}.Values.isFeasible.Data, '-bo', 'MarkerIndices',1:200:length(trn_time));
plot(trn_time, logsout{idx_WTPGOutput}.Values.wtpgValid.Data, '-g*', 'MarkerIndices',51:200:length(trn_time));
plot(trn_time, logsout{idx_WTPGOutput}.Values.wtpgOffTerMap.Data, '-mx', 'MarkerIndices',101:200:length(trn_time));
ylabel('True/False'); grid on; axis([trn_time(1), trn_time(end), -0.2, 2.0]);
xlabel('Time (sec)');
legend('WTPG Feasibility', 'WTPG Valid', 'WTPG Off Terrain Map', 'Location', 'Best', 'Orientation', 'horizontal');
nexttile
title('GCAWS Slected/Feasible')
hold on;
plot(gcaws_time, logsout{idx_GCAWSInput}.Values.selected.Data, '-bo', 'MarkerIndices',1:200:length(gcaws_time))
plot(gcaws_time, logsout{idx_GCAWSDatapump}.Values.isFeasible.Data, '-r*', 'MarkerIndices',51:200:length(gcaws_time));
ylabel('True/False'); grid on; axis([trn_time(1), trn_time(end), -0.2, 2.0]);
xlabel('Time (sec)');
legend('GCAWS Selected', 'GCAWS Feasible');

% GCAWS_3_BasicFunction
gcawsBasicFig = figure('WindowStyle', 'docked');
tiledlayout(5,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
hold on;
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.selectedWrap.Data, 'DisplayName', 'Selected Wrap');
grid on; axis([gcaws_time(1), gcaws_time(end), -0.1, 1.1]);
legend;
nexttile
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.collisionWarningValid.Data, 'DisplayName', 'Warning Valid');
legend;
grid on; axis([gcaws_time(1), gcaws_time(end), -0.1, 1.1]);
nexttile
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.collisionWarning.Data, 'DisplayName', 'Collision Warning');
ylabel('True/False');
legend;
grid on; axis([gcaws_time(1), gcaws_time(end), -0.1, 1.1]);
nexttile
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.timeToGoToPullUp.Data, '-r', 'DisplayName', 'TTGTPU');
ylabel('Time Left (sec)');
grid on; axis([gcaws_time(1), gcaws_time(end), -0.5, 5]);
legend;
nexttile
hold on;
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.swhUsedWrap.Data, 'DisplayName', 'SWH Used Wrap');
% plot(wtpg_time, logsout{idx_TRNInput}.Values.RALTAltitude.Data, 'DisplayName', 'RALT AGL');
grid on; axis([gcaws_time(1), gcaws_time(end), -50, max(logsout{idx_GCAWSOutput}.Values.swhUsedWrap.Data) + 200]);
legend;
xlabel('Time (sec)');

% GCAWS_4_Output
gcawsOutputFig = figure('WindowStyle', 'docked');
tiledlayout(5,1, 'TileSpacing', 'tight', 'Padding', 'compact');

nexttile
title('Output: Time-to-Go-to-Pull-Up')
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.timeToGoToPullUp.Data, '-r', 'DisplayName', 'TTGTPU');
ylabel('Time Left (sec)');
grid on; axis([gcaws_time(1), gcaws_time(end), -0.5, 5]);
legend;

nexttile
title('SWH Used Wrap')
hold on;
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.swhUsedWrap.Data, 'DisplayName', 'SWH Used Wrap');
grid on; axis([gcaws_time(1), gcaws_time(end), -50, max(logsout{idx_GCAWSOutput}.Values.swhUsedWrap.Data) + 200]);
legend;
xlabel('Time (sec)');

nexttile
title('Output: Collision Warning & Validity')
hold on;
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.collisionWarning.Data, 'DisplayName', 'Collision Warning');
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.collisionWarningValid.Data, 'DisplayName', 'Warning Valid');
ylabel('True/False');
grid on; axis([gcaws_time(1), gcaws_time(end), -0.2, 2.0]);
legend('Location', 'Best');

nexttile
title('Output: Wraparound');
hold on;
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.selectedWrap.Data, '-ro', 'DisplayName', 'Selected Wrap', 'MarkerIndices',1:350:length(gcaws_time));
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.paramSelectedValid.Data, '-b*', 'DisplayName', 'Parameter Selected Valid', 'MarkerIndices',51:350:length(gcaws_time));
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.reactionTimeWrap.Data, '-gx', 'DisplayName', 'Reaction Time wrap', 'MarkerIndices',101:350:length(gcaws_time));
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.rollRateWrap.Data, '-csquare', 'DisplayName', 'Roll rate wrap', 'MarkerIndices',151:350:length(gcaws_time));
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.gOnsetRateWrap.Data, '-mdiamond', 'DisplayName', 'g-onset rate wrap', 'MarkerIndices',201:350:length(gcaws_time));
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.pullUpGWrap.Data, '-y^', 'DisplayName', 'Pull-up g Wrap', 'MarkerIndices',251:350:length(gcaws_time));
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.maxFlightPathAngleWrap.Data, '-kv', 'DisplayName', 'Max Flight Path Angle Wrap', 'MarkerIndices',301:350:length(gcaws_time));
plot(gcaws_time, logsout{idx_GCAWSInput}.Values.gsAvailable.Data, 'Color', "#7E2F8E", 'LineStyle', ':', 'Marker', 'pentagram', 'DisplayName', 'Gs Available', 'MarkerIndices', 326:350:length(gcaws_time))
hold off;
legend('Location', 'Best', 'Orientation', 'horizontal', 'NumColumns', 4); grid on; axis([gcaws_time(1), gcaws_time(end), -0.5, 10.0]);

nexttile
title('TRN Estimate Uncertainties');
hold on;
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNHorUncertainty.Data, '-r', 'DisplayName', 'TRN Horizontal Unc.');
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNVerticalUncertainty.Data, '-b', 'DisplayName', 'TRN Vertical Unc.');
xlabel('Time (sec)');
ylabel('Uncertainty (m)');
grid minor;
axis([gcaws_time(1), gcaws_time(end), -10, max(logsout{idx_TRNOutput}.Values.TRNHorUncertainty.Data) + 50]);
legend('Location', 'best');


collisionIndices = (1 == diff(logsout{idx_GCAWSOutput}.Values.collisionWarning.Data));
hlIndices = find(collisionIndices);

if isempty(hlIndices)
    % pass
else

    for ff = 1:length(hlIndices)
        figure('WindowStyle', 'docked');
        tiledlayout(7,1, 'TileSpacing', 'tight', 'Padding', 'compact');
        
        % Index at Collision Warning Detected
        curIdx  = hlIndices(ff);

        startIdx = curIdx - 250;
        if startIdx <= 0
            startIdx  = 1;
        end
        endIdx = curIdx + 250;
        if endIdx > length(logsout{idx_GCAWSOutput}.Values.collisionWarning.Data)
            endIdx = length(logsout{idx_GCAWSOutput}.Values.collisionWarning.Data);
        end

        plot_time = logsout{idx_GCAWSOutput}.Values.collisionWarning.Time(startIdx:endIdx);
        
        % [1] TIME-TO-GO-TO-PULL-UP
        nexttile(1)
        hold on;
        plot(plot_time, logsout{idx_GCAWSOutput}.Values.timeToGoToPullUp.Data(startIdx:endIdx), '-r', 'DisplayName', 'TTGTPU');
        plot(logsout{idx_GCAWSOutput}.Values.collisionWarning.Time(curIdx), logsout{idx_GCAWSOutput}.Values.collisionWarning.Data(curIdx), 'bo', 'DisplayName', 'Current TTGTPU');
        ylabel('Time Left (sec)');
        grid on; axis([plot_time(1), plot_time(end), -0.5, 5.5]);
        legend('Location', 'best');

        % [2] GCAWS Trajetories at Collision Warning
        nexttile(2, [2 1])
        hold on;
        plot(logsout{idx_TRNOutput}.Values.Altitude_INS.Time(curIdx), logsout{idx_TRNOutput}.Values.Altitude_INS.Data(curIdx) - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data(curIdx), 'ro', 'DisplayName', 'TRN Alt. at Collision Warning');
        plot(plot_time, (logsout{idx_TRNOutput}.Values.Altitude_INS.Data(startIdx:endIdx) - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data(startIdx:endIdx)), '-b', 'DisplayName', 'TRN Trajectory');
        plot(plot_time, (wtpg_wtp(startIdx:endIdx,1)+gcaws_swh(startIdx:endIdx)), 'Color', [0.1 0.8 0.1], 'LineWidth', 3, 'DisplayName', 'Terrain Under TRN Est. + SWH');

        recIdx = curIdx;
        recTime = linspace(logsout{idx_GCAWSOutput}.Values.collisionWarning.Time(recIdx), logsout{idx_GCAWSOutput}.Values.collisionWarning.Time(recIdx)+20, 41);
%         plot(recTime(1), logsout{idx_GCAWSDatapump}.Values.predictedTrajectory.Data(1, 2, 3, recIdx), 'ro');
        for rr = 2:12
            temp = plot(recTime, logsout{idx_GCAWSDatapump}.Values.predictedTrajectory.Data(:, rr, 3, recIdx), ':k','LineWidth', 2, 'DisplayName', 'Recovery Trajectory');
            if rr ~= 2
                temp.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        end

        plot(recTime, logsout{idx_WTPGOutput}.Values.wtp.Data(recIdx, :) + gcaws_swh(recIdx), '-r', 'LineWidth', 3, 'DisplayName', 'WTP at Col. Warn. + SWH');
        title('GCAWS Trajectories (at Collision Warning)')
        legend('Location', 'Best', 'Orientation', 'horizontal', 'NumColumns', 4);
        ylabel('Height (m)');
        grid on;
        axis([plot_time(1), plot_time(end), min(wtpg_wtp(startIdx:endIdx,1))-200, logsout{idx_GCAWSDatapump}.Values.predictedTrajectory.Data(1, 1, 3, recIdx)+200]);
        hold off;
        
        nexttile(4, [2 1])
        hold on;
        if (recIdx-10) > 1
            plot(logsout{idx_TRNOutput}.Values.Altitude_INS.Time(curIdx-10), logsout{idx_TRNOutput}.Values.Altitude_INS.Data(curIdx-10) - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data(curIdx-10), 'ro', 'DisplayName', 'TRN Alt. at 0.4 sec before Collision Warning');
            plot(plot_time, (logsout{idx_TRNOutput}.Values.Altitude_INS.Data(startIdx:endIdx) - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data(startIdx:endIdx)), '-b', 'DisplayName', 'TRN Trajectory');
            plot(plot_time, (wtpg_wtp(startIdx:endIdx,1)+gcaws_swh(startIdx:endIdx)), 'Color', [0.1 0.8 0.1], 'LineWidth', 3, 'DisplayName', 'Terrain Under TRN Est. + SWH');

            plot(recTime-0.4, logsout{idx_WTPGOutput}.Values.wtp.Data(recIdx-10, :) + gcaws_swh(recIdx-10), '-r', 'LineWidth', 3, 'DisplayName', 'WTP at 0.4sec Before Col. Warn. + SWH');
        end
        if recIdx > 10
            for rr = 2:12
                temp = plot(recTime-0.4, logsout{idx_GCAWSDatapump}.Values.predictedTrajectory.Data(:, rr, 3, recIdx - 10), ':k', 'LineWidth', 2,'DisplayName', 'Recovery Trajectory Before 0.4sec');
                if rr ~= 2
                    temp.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
            end
        end
        title('GCAWS Trajectories (0.4 sec Before)')
        legend('Location', 'Best', 'Orientation', 'horizontal', 'NumColumns', 4);
        ylabel('Height (m)');
        grid on;
        axis([plot_time(1), plot_time(end), min(wtpg_wtp(startIdx:endIdx,1))-100, logsout{idx_GCAWSDatapump}.Values.predictedTrajectory.Data(1, 1, 3, recIdx)+200]);
        hold off;

        nexttile(6)
        title('Output: Collision Warning & Validity')
        hold on;
        plot(plot_time, logsout{idx_GCAWSOutput}.Values.collisionWarning.Data(startIdx:endIdx), 'DisplayName', 'Collision Warning');
        plot(plot_time, logsout{idx_GCAWSOutput}.Values.collisionWarningValid.Data(startIdx:endIdx), 'DisplayName', 'Warning Valid');
        ylabel('True/False');
        grid on; axis([plot_time(1), plot_time(end), -0.2, 2.0]);
        legend('Location', 'Best');
        %
        %         nexttile(5)
        %         title('Output: Wraparound');
        %         hold on;
        %         plot(cur_time, logsout{idx_GCAWSOutput}.Values.selectedWrap.Data(startIdx:endIdx), '-ro', 'DisplayName', 'Selected Wrap', 'MarkerIndices',1:30:length(cur_time));
        %         plot(cur_time, logsout{idx_GCAWSOutput}.Values.paramSelectedValid.Data(startIdx:endIdx), '-b*', 'DisplayName', 'Parameter Selected Valid', 'MarkerIndices',5:30:length(cur_time));
        %         plot(cur_time, logsout{idx_GCAWSOutput}.Values.reactionTimeWrap.Data(startIdx:endIdx), '-gx', 'DisplayName', 'Reaction Time wrap', 'MarkerIndices',9:30:length(cur_time));
        %         plot(cur_time, logsout{idx_GCAWSOutput}.Values.rollRateWrap.Data(startIdx:endIdx), '-csquare', 'DisplayName', 'Roll rate wrap', 'MarkerIndices',14:30:length(cur_time));
        %         plot(cur_time, logsout{idx_GCAWSOutput}.Values.gOnsetRateWrap.Data(startIdx:endIdx), '-mdiamond', 'DisplayName', 'g-onset rate wrap', 'MarkerIndices',18:30:length(cur_time));
        %         plot(cur_time, logsout{idx_GCAWSOutput}.Values.pullUpGWrap.Data(startIdx:endIdx), '-y^', 'DisplayName', 'Pull-up g Wrap', 'MarkerIndices',22:30:length(cur_time));
        %         plot(cur_time, logsout{idx_GCAWSOutput}.Values.maxFlightPathAngleWrap.Data(startIdx:endIdx), '-kv', 'DisplayName', 'Max Flight Path Angle Wrap', 'MarkerIndices',26:30:length(cur_time));
        %         hold off;
        %         legend('Location', 'Best', 'Orientation', 'horizontal'); grid on; axis([cur_time(1), cur_time(end), -0.5, 10.0]);

        nexttile(7)
        title('TRN Estimate Uncertainties');
        hold on;
        plot(plot_time, logsout{idx_TRNOutput}.Values.TRNHorUncertainty.Data(startIdx:endIdx), '-r', 'DisplayName', 'TRN Horizontal Unc.');
        plot(plot_time, logsout{idx_TRNOutput}.Values.TRNVerticalUncertainty.Data(startIdx:endIdx), '-b', 'DisplayName', 'TRN Vertical Unc.');
        xlabel('Time (sec)');
        ylabel('Uncertainty (m)');
        grid minor;
        axis([plot_time(1), plot_time(end), -10, max(logsout{idx_TRNOutput}.Values.TRNHorUncertainty.Data) + 50]);
        legend('Location', 'best');

        if ON_SAVE
            saveas(gcf, [HEADER, sprintf('\\GCAWSOutput_Highlighted_%03d', ff), TIMESTR, '.fig']);
            saveas(gcf, [HEADER, sprintf('\\GCAWSOutput_Highlighted_%03d', ff), TIMESTR, '.png']);
        end
        if ff > 100
            break;
        end
    end
end


if ON_SAVE
    saveas(gcawsFeasibilityFig, [HEADER, '\GCAWSFeasibility_', TIMESTR, '.fig']);
    saveas(gcawsOutputFig, [HEADER, '\GCAWSOutput_', TIMESTR, '.fig']);
    saveas(gcawsFig, [HEADER, '\GCAWS_', TIMESTR, '.fig']);

    saveas(gcawsFeasibilityFig, [HEADER, '\GCAWSFeasibility_', TIMESTR, '.png']);
    saveas(gcawsOutputFig, [HEADER, '\GCAWSOutput_', TIMESTR, '.png']);
    saveas(gcawsFig, [HEADER, '\GCAWS_', TIMESTR, '.png']);
  
end

end