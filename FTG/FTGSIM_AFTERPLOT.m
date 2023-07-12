%%FTGSIM_AFTERPLOT
% FTGSIM_AFTERPLOT 시뮬레이션 종료 후 로그 데이터를 출력한다.
%
global ON_SAVE TIMESTR HEADER

ON_SAVE = true;

TIMESTR = datestr(now, 'yyyy-mm-dd_HHMMSS');
HEADER = ['logs\', TIMESTR];
mkdir(HEADER);

%% IDENTIFY logging index

dataName = logsout.getElementNames;
numData = length(dataName);
for ii = 1:numData
    if strcmp(dataName{ii}, 'OUT_TRUE_STATE')
        idx_true = ii;
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
    elseif strcmp(dataName{ii}, 'IN_CRInput')
        idx_CRInput = ii;
    elseif strcmp(dataName{ii}, 'OUT_CROutput')
        idx_CROutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_CRDatapump')
        idx_CRDatapump = ii;
    elseif strcmp(dataName{ii}, 'IN_LOSRInput')
        idx_LOSRInput = ii;
    elseif strcmp(dataName{ii}, 'OUT_LOSROutput')
        idx_LOSROutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_LOSRDatapump')
        idx_LOSRDatapump = ii;
    end
end


%% CALCULATION

% Worst-case Terrain profile
wtpg_wtp = logsout{idx_WTPGOutput}.Values.wtp.Data;

% MSD, SWH
dtf_msd = logsout{idx_DTFOutput}.Values.msdUsedWrap.Data;
gcaws_swh = logsout{idx_GCAWSOutput}.Values.swhUsedWrap.Data;


%% PLOT
% Plot true flight data
% plotTrueData(logsout, idx_true);

% Plot TRN Data
% plotTRNData(trn_time, logsout, idx_true, idx_TRNInput, idx_TRNOutput);

% Plot DTF Data
% plotDTFData(logsout, idx_WTPGOutput, idx_true, idx_DTFDatapump, idx_DTFOutput, idx_TRNOutput, idx_WTPGDatapump, idx_DTFInput);

% Plot GCAWS Data
plotGCAWSData(logsout, idx_true, idx_GCAWSOutput, idx_TRNOutput, idx_WTPGDatapump, idx_WTPGOutput, idx_GCAWSDatapump);

% Plot AGR Data
% plotAGRData(logsout, idx_LOSROutput, idx_CROutput);

% /* */
% figure;
% hold on;
% pHull = plot(0:0.5:20, logsout{idx_GCAWSDatapump}.Values.hulledProfile.Data(1,:));
% pWTP = plot(0:0.5:20, wtpg_wtp(1,:));
% for ii = 2:size(wtpg_wtp, 1)
%     pHull.YData = logsout{idx_GCAWSDatapump}.Values.hulledProfile.Data(ii,:);
%     pWTP.YData = wtpg_wtp(ii,:);
%     drawnow
% end
%% Internal Functions -----------------------------------------------------
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
function plotTrueData(logsout, idx_true)
global ON_SAVE TIMESTR HEADER

% /* True Position  */
posFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_true}.Values.latitude);
ylabel('Latitude (deg)');
nexttile
plot(logsout{idx_true}.Values.longitude);
ylabel('Longitude (deg)');
nexttile
plot(logsout{idx_true}.Values.altitude);
ylabel('Altitude (m)');

% /* True Attitude  */
attFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_true}.Values.roll * 180/pi);
ylabel('Roll (deg)');
nexttile
plot(logsout{idx_true}.Values.pitch * 180/pi);
ylabel('Pitch (deg)');
nexttile
plot(logsout{idx_true}.Values.yaw * 180/pi);
ylabel('Yaw (deg)');

% /* True Body Velocity  */
velFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_true}.Values.velX);
ylabel('Body Velocity X (m/s)')
nexttile
plot(logsout{idx_true}.Values.velY);
ylabel('Body Velocity Y (m/s)')
nexttile
plot(logsout{idx_true}.Values.velZ);
ylabel('Body Velocity Z (m/s)')

% /* True Body Acceleration  */
accFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_true}.Values.accX);
ylabel('Body Acceleration X (m/ss)');
nexttile
plot(logsout{idx_true}.Values.accY);
ylabel('Body Acceleration Y (m/ss)');
nexttile
plot(logsout{idx_true}.Values.accZ);
ylabel('Body Acceleration Z (m/ss)');

% /* Actuators  */
actFig = figure('WindowStyle', 'docked');
tiledlayout(4,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(logsout{idx_true}.Values.throttleN);
ylabel('Throttle (N)')
nexttile
plot(logsout{idx_true}.Values.elevatorDeg);
ylabel('Elevator (deg)')
nexttile
plot(logsout{idx_true}.Values.aileronDeg);
ylabel('Aileron (deg)')
nexttile
plot(logsout{idx_true}.Values.rudderDeg);
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

function plotTRNData(logsout, idx_true, idx_TRNInput, idx_TRNOutput)
global ON_SAVE TIMESTR HEADER

% Time
trn_time = logsout{idx_TRNOutput}.Values.TRNTimeTag.Time;

% /* True, INS, TRN Position  */
trnFig = figure('WindowStyle', 'docked');
sgtitle('Estimate of Position');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
hold on;
plot(trn_time, logsout{idx_true}.Values.latitude.Data, '-b', 'Linewidth', 1);
plot(trn_time, logsout{idx_TRNOutput}.Values.Latitude_INS.Data, '-g', 'Linewidth', 1);
plot(trn_time, logsout{idx_TRNInput}.Values.Latitude_GPS.Data, '-m', 'Linewidth', 1);
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNEstOfLatitude.Data, '-r', 'Linewidth', 2);
ylabel('Latitude (deg)');
legend('True', 'INS', 'GPS', 'TRN Estimate');
grid on;
hold off;
nexttile
hold on;
plot(trn_time, logsout{idx_true}.Values.longitude.Data, '-b', 'Linewidth', 1);
plot(trn_time, logsout{idx_TRNOutput}.Values.Longitude_INS.Data, '-g', 'Linewidth', 1);
plot(trn_time, logsout{idx_TRNInput}.Values.Longitude_GPS.Data, '-m', 'Linewidth', 1);
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNEstOfLongitude.Data, '-r', 'Linewidth', 2);
ylabel('Longitude (deg)');
legend('True', 'INS', 'GPS', 'TRN Estimate');
grid on;
hold off;
nexttile
hold on;
plot(trn_time, logsout{idx_true}.Values.altitude.Data, '-b', 'Linewidth', 1);
plot(trn_time, logsout{idx_TRNOutput}.Values.Altitude_INS.Data, '-g', 'Linewidth', 1);
plot(trn_time, logsout{idx_TRNInput}.Values.Altitude_GPS.Data, '-m', 'Linewidth', 1);
plot(trn_time, logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data, '-r', 'Linewidth', 2);
ylabel('Altitude (m)');
legend('True', 'INS', 'GPS', 'TRN Estimate');
grid on;
hold off;

[latLength, lonLength] = getLengthOfADegree(logsout{idx_true}.Values.latitude.Data(1), logsout{idx_true}.Values.altitude.Data(1));

trnErrFig = figure('WindowStyle', 'docked');
sgtitle('TRN Errors');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
hold on;
plot(trn_time, (logsout{idx_true}.Values.latitude.Data - logsout{idx_TRNOutput}.Values.TRNEstOfLatitude.Data) * latLength);
ylabel('Latitude error (m)');
grid on;
hold off;
nexttile
hold on;
plot(trn_time, (logsout{idx_true}.Values.longitude.Data - logsout{idx_TRNOutput}.Values.TRNEstOfLongitude.Data) * lonLength);
ylabel('Longitude error (m)');
grid on;
hold off;
nexttile
hold on;
plot(trn_time, logsout{idx_true}.Values.altitude.Data - logsout{idx_TRNOutput}.Values.Altitude_INS.Data + logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data);
ylabel('Altitude error (m)');
grid on;
hold off;

if ON_SAVE
    saveas(trnFig, [HEADER, '\TRNEstimate_', TIMESTR, '.fig']);
    saveas(trnErrFig, [HEADER, '\TRNEstimateError_', TIMESTR, '.fig']);

    saveas(trnFig, [HEADER, '\TRNEstimate_', TIMESTR, '.png']);
    saveas(trnErrFig, [HEADER, '\TRNEstimateError_', TIMESTR, '.png']);
end

end

function plotDTFData(logsout, idx_WTPGOutput, idx_true, idx_DTFDatapump, idx_DTFOutput, idx_TRNOutput, idx_WTPGDatapump, idx_DTFInput)
global ON_SAVE TIMESTR HEADER

% Time
trn_time = logsout{idx_TRNOutput}.Values.TRNTimeTag.Time;
dtf_time = logsout{idx_DTFDatapump}.Values.referenceTrajectory.Time;

% MSD, SWH
dtf_msd = logsout{idx_DTFOutput}.Values.msdUsedWrap.Data;

% Worst-case Terrain profile
wtpg_wtp = logsout{idx_WTPGOutput}.Values.wtp.Data;

% violateTRN = (logsout{idx_WCPGOutput}.Values.WCPGWCP(:,1).Data + dbtf_msd) - (logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data + logsout{idx_TRNOutput}.Values.Altitude_INS.Data);
violateTrue = (wtpg_wtp + dtf_msd) - logsout{idx_true}.Values.altitude.Data;

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
plot(logsout{idx_true}.Values.altitude, '-r');
% plot(logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Time, logsout{idx_TRNOutput}.Values.Altitude_INS.Data + logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data ,...
%     '-b');
plot(dtf_time, wtpg_wtp(:,1)+dtf_msd, 'Color', [0.1 0.8 0.1], 'LineWidth', 3);

% plot(dbtf_time, dbtf_traj);
plot(dtf_time, dtf_trajOneSec, 'Color', '#960096');
plot(dtf_time, dtf_trajTwoSec, 'Color', '#A500A5');
plot(dtf_time, dtf_trajThrSec, 'Color', '#B400B4');
plot(dtf_time, dtf_trajFouSec, 'Color', '#CC00CC');
ylabel('Height (m)');
grid minor;
axis([dtf_time(1), dtf_time(end), 0, max(dtf_trajOneSec)+100]);
legend('Flight Trajectory', 'Terrain + MSD', ...
    'DTF Ref. Traj.(before 1 sec)', 'DTF Ref. Traj.(before 2 sec)', 'DTF Ref. Traj.(before 3 sec)', 'DTF Ref. Traj.(before 4 sec)', ...
    'Location', 'best', 'FontSize', 7);
nexttile
accX = logsout{idx_true}.Values.accX.Data;
accY = logsout{idx_true}.Values.accY.Data;
accZ = logsout{idx_true}.Values.accZ.Data;
cRol = cos(logsout{idx_true}.Values.roll.Data);
sRol = sin(logsout{idx_true}.Values.roll.Data);
cPit = cos(logsout{idx_true}.Values.pitch.Data);
sPit = sin(logsout{idx_true}.Values.pitch.Data);
accVert = -sPit .* accX + sRol.*cPit .* accY + cRol.*cPit .* accZ;
hold on;
plot(dtf_time, logsout{idx_DTFOutput}.Values.verticalGCommand.Data*9.8);
plot(dtf_time, -1 * logsout{idx_true}.Values.accZ.Data);
plot(dtf_time, -accVert);
legend('Normal Acceleration Command', 'Body Acceleration (-Z)', 'Vertical Acceleration (+Up)', 'FontSize', 7);
ylabel('Acceleration (m/s^2)');
grid minor;
axis([dtf_time(1), dtf_time(end), -15, 25]);
nexttile
hold on;
TRN_SD = logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data - wtpg_wtp(:,1);
TRUE_SD = logsout{idx_true}.Values.altitude.Data - wtpg_wtp(:,1);
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
accX = logsout{idx_true}.Values.accX.Data;
accY = logsout{idx_true}.Values.accY.Data;
accZ = logsout{idx_true}.Values.accZ.Data;
cRol = cos(logsout{idx_true}.Values.roll.Data);
sRol = sin(logsout{idx_true}.Values.roll.Data);
cPit = cos(logsout{idx_true}.Values.pitch.Data);
sPit = sin(logsout{idx_true}.Values.pitch.Data);
accVert = -sPit .* accX + sRol.*cPit .* accY + cRol.*cPit .* accZ;
hold on;
plot(dtf_time, logsout{idx_DTFOutput}.Values.verticalGCommand.Data*9.8);
plot(dtf_time, -1 * logsout{idx_true}.Values.accZ.Data);
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

function plotGCAWSData(logsout, idx_true, idx_GCAWSOutput, idx_TRNOutput, idx_WTPGDatapump, idx_WTPGOutput, idx_GCAWSDatapump)
global ON_SAVE TIMESTR HEADER

% Time
trn_time = logsout{idx_TRNOutput}.Values.TRNTimeTag.Time;
gcaws_time = logsout{idx_GCAWSOutput}.Values.selectedWrap.Time;
wtpg_time = logsout{idx_WTPGOutput}.Values.wtp.Time;

% Worst-case Terrain profile
wtpg_wtp = logsout{idx_WTPGOutput}.Values.wtp.Data;

% SWH
gcaws_swh = logsout{idx_GCAWSOutput}.Values.swhUsedWrap.Data;

% GCAWS Profile
gcawsFig = figure('WindowStyle', 'docked');
tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
hold on;
plot(trn_time, logsout{idx_true}.Values.altitude.Data./0.3048, '-r');
% plot(trn_time, (logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data), '-b');
plot(wtpg_time, (wtpg_wtp(:,1)+gcaws_swh)./0.3048, 'Color', [0.1 0.8 0.1], 'LineWidth', 3);
plot(wtpg_time, wtpg_wtp(:,1)./0.3048, 'Color', [0.6 0.8 0.4], 'LineWidth', 3);
legend('True Flight Trajectory', 'Terrain + SWH', 'Terrain', 'Location', 'Best');
ylabel('Height (ft)');
grid on;
axis([gcaws_time(1), gcaws_time(end), 0, max(logsout{idx_true}.Values.altitude.Data./0.3048)+500]);
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

% GCAWS Feasibility
gcawsFeasibilityFig = figure('WindowStyle', 'docked');
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
hold on;
plot(trn_time, logsout{idx_WTPGDatapump}.Values.isFeasible.Data, '-b');
plot(trn_time, logsout{idx_WTPGOutput}.Values.wtpgValid.Data, '-r');
plot(trn_time, logsout{idx_WTPGOutput}.Values.wtpgOffTerMap.Data, '-m');
ylabel('Validity'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
xlabel('Time (sec)');
legend('WTPG Feasibility', 'WTPG Valid', 'WTPG Off Terrain Map');
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.TRN_VelocityX.Data, '-b');
ylabel('Speed (m/s)'); grid on; axis([trn_time(1), trn_time(end), 0, 450]);
xlabel('Time (sec)');
nexttile
plot(gcaws_time, logsout{idx_GCAWSDatapump}.Values.isFeasible.Data, '-r');
ylabel('GCAWS Feasibility'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
xlabel('Time (sec)');


% GCAWS Basic Function Test
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

% GCAWS Parameters
gcawsParamFig = figure('WindowStyle', 'docked');
tiledlayout(6,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.paramSelectedValid.Data, 'DisplayName', 'Parameter Selected Valid');
legend; grid on; axis([gcaws_time(1), gcaws_time(end), -0.1, 5.0]);
nexttile
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.reactionTimeWrap.Data, 'DisplayName', 'Reaction Time wrap');
legend; grid on; axis([gcaws_time(1), gcaws_time(end), -0.1, 5.0]);
nexttile
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.rollRateWrap.Data, 'DisplayName', 'Roll rate wrap');
legend; grid on; axis([gcaws_time(1), gcaws_time(end), -0.1, 5.0]);
nexttile
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.gOnsetRateWrap.Data, 'DisplayName', 'g-onset rate wrap');
legend; grid on; axis([gcaws_time(1), gcaws_time(end), -0.1, 5.0]);
nexttile
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.pullUpGWrap.Data, 'DisplayName', 'Pull-up g Wrap');
legend; grid on; axis([gcaws_time(1), gcaws_time(end), -0.1, 5.0]);
nexttile
plot(gcaws_time, logsout{idx_GCAWSOutput}.Values.maxFlightPathAngleWrap.Data, 'DisplayName', 'Max Flight Path Angle Wrap');
legend; grid on; axis([gcaws_time(1), gcaws_time(end), -0.1, 5.0]);
xlabel('Time (sec)');

% TRN Uncertainty
GCAWSTRNUncFig = figure('WindowStyle', 'docked');
tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNHorUncertainty.Data);
grid on; axis([trn_time(1), trn_time(end), -5, max(logsout{idx_TRNOutput}.Values.TRNHorUncertainty.Data)+5]);
ylabel('Hor. Uncertainty (m)');
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.TRNVerticalUncertainty.Data);
grid on; axis([trn_time(1), trn_time(end), -5, max(logsout{idx_TRNOutput}.Values.TRNVerticalUncertainty.Data)+5]);
ylabel('Ver. Uncertainty (m)');
xlabel('Time (sec)');

if ON_SAVE
    saveas(gcawsFeasibilityFig, [HEADER, '\GCAWSFeasibility_', TIMESTR, '.fig']);
    saveas(gcawsBasicFig, [HEADER, '\GCAWSBasicTest_', TIMESTR, '.fig']);
    saveas(gcawsParamFig, [HEADER, '\GCAWSParameters_', TIMESTR, '.fig']);
    saveas(gcawsFig, [HEADER, '\GCAWSOutput_', TIMESTR, '.fig']);

    saveas(gcawsFeasibilityFig, [HEADER, '\GCAWSFeasibility_', TIMESTR, '.png']);
    saveas(gcawsBasicFig, [HEADER, '\GCAWSBasicTest_', TIMESTR, '.png']);
    saveas(gcawsParamFig, [HEADER, '\GCAWSParameters_', TIMESTR, '.png']);
    saveas(gcawsFig, [HEADER, '\GCAWSOutput_', TIMESTR, '.png']);
  
end

end

function plotAGRData(logsout, idx_LOSROutput, idx_CROutput)
% Time
cr_time = logsout{idx_CROutput}.Values.CRDataTagWrap.Time;
losr_time = logsout{idx_LOSROutput}.Values.LOSRDataTagWrap.Time;

losrFig = figure('WindowStyle', 'docked');
tiledlayout(6,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRSelectedWraparound.Data, 'DisplayName', 'LOSR selected wrap');
legend; grid on; 
axis([losr_time(1), losr_time(end), -1.0, 2.5]);
nexttile
hold on;
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRDataValid.Data, 'DisplayName', 'LOSR valid');
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetOffMap.Data, 'DisplayName', 'Target off map');
legend; grid on; axis([losr_time(1), losr_time(end), -1.0, 2.5]);
nexttile
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRDataTagWrap.Data, 'DisplayName', 'Data tag wrap');
legend; grid on; axis([losr_time(1), losr_time(end), -0.1, 5.0]);
nexttile
hold on;
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetElevation.Data, 'DisplayName', 'Target elevation');
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetRangeEast.Data, 'DisplayName', 'Target range north');
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetRangeNorth.Data, 'DisplayName', 'Target range east');
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetRangeVertical.Data, 'DisplayName', 'Target range vertical');
legend; grid on; 
axis([losr_time(1), losr_time(end), -20000, 20000]);
ylabel('Range (m)');
nexttile
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetLatitude.Data, 'DisplayName', 'Target latitude');
axis([losr_time(1), losr_time(end), min(nonzeros(logsout{idx_LOSROutput}.Values.LOSRTargetLatitude.Data))-0.05, max(nonzeros(logsout{idx_LOSROutput}.Values.LOSRTargetLatitude.Data))+0.05 ])
legend; grid on; 
ylabel('Latitude (deg)');
nexttile
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetLongitude.Data, 'DisplayName', 'Target longitude');
legend; grid on; 
axis([losr_time(1), losr_time(end), min(nonzeros(logsout{idx_LOSROutput}.Values.LOSRTargetLongitude.Data))-0.05, max(nonzeros(logsout{idx_LOSROutput}.Values.LOSRTargetLongitude.Data))+0.05 ])
ylabel('Longitude (deg)');

crFig = figure('WindowStyle', 'docked');
tiledlayout(5,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(cr_time, logsout{idx_CROutput}.Values.CRSelectedWraparound.Data, 'DisplayName', 'CR selected wrap');
legend; grid on; axis([cr_time(1), cr_time(end), -1.0, 2.5]);
nexttile
hold on;
plot(cr_time, logsout{idx_CROutput}.Values.CRValid.Data, 'DisplayName', 'CR valid');
plot(cr_time, logsout{idx_CROutput}.Values.CRElevationValid.Data, 'DisplayName', 'CR elevation valid');
legend; grid on; axis([cr_time(1), cr_time(end), -1.0, 2.5]);
nexttile
plot(cr_time, logsout{idx_CROutput}.Values.CRTargetElevation.Data, 'DisplayName', 'CR target elevation');
legend; grid on; axis([cr_time(1), cr_time(end), -100, max(logsout{idx_CROutput}.Values.CRTargetElevation.Data)+200]);
nexttile
plot(cr_time, logsout{idx_CROutput}.Values.CRTargetLOSAvailable.Data, 'DisplayName', 'CR target LOS Available');
legend; grid on; axis([cr_time(1), cr_time(end), -1.0, 2.5]);
nexttile
plot(cr_time, logsout{idx_CROutput}.Values.CRTargetLOSAvailableValid.Data, 'DisplayName', 'CR target LOS Available Valid');
legend; grid on; axis([cr_time(1), cr_time(end), -1.0, 2.5]);
end