%%KAISIM_GCAWS_AFTERPLOT
% KAISIM 시뮬레이션 종료 후 로그 데이터를 출력한다.
%
global ON_SAVE TIMESTR HEADER

ON_SAVE = true;

TIMESTR = datestr(now, 'yyyy-mm-dd_HHMMSS');
HEADER = ['GCAWS_Logs\', TIMESTR];
mkdir(HEADER);

%% IDENTIFY
dataName = logsout.getElementNames;
numData = length(dataName);
for ii = 1:numData
    if strcmp(dataName{ii}, 'OUT_TRUE_STATE')
        idx_true = ii;
    elseif strcmp(dataName{ii}, 'IN_TRNRef')
        idx_TRNRef = ii;
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
%     elseif strcmp(dataName{ii}, 'IN_DTFInput')
%         idx_DTFInput = ii;
%     elseif strcmp(dataName{ii}, 'OUT_DTFOutput')
%         idx_DTFOutput = ii;
%     elseif strcmp(dataName{ii}, 'OUT_DTFDatapump')
%         idx_DTFDatapump = ii;
    end
end

%% PLOT
% % Plot true flight data
% plotTrueData(logsout, idx_true);

% % Plot TRN Data
% plotTRNData(logsout, idx_true, idx_TRNInput, idx_TRNOutput);

% % Plot GCAWS Data
plotGCAWSData(logsout, idx_TRNRef, idx_GCAWSInput, idx_GCAWSOutput, idx_TRNOutput, idx_WTPGDatapump, idx_WTPGOutput, idx_GCAWSDatapump);

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

%% PLOT FUNCTIONS ---------------------------------------------------------
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

function plotGCAWSData(logsout, idx_TRNRef, idx_GCAWSInput, idx_GCAWSOutput, idx_TRNOutput, idx_WTPGDatapump, idx_WTPGOutput, idx_GCAWSDatapump)
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
plot(trn_time, logsout{idx_TRNRef}.Values.Altitude.Data, '-r');
% plot(trn_time, (logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data), '-b');
plot(wtpg_time, (wtpg_wtp(:,1)+gcaws_swh), 'Color', [0.1 0.8 0.1], 'LineWidth', 3);
plot(wtpg_time, wtpg_wtp(:,1), 'Color', [0.6 0.8 0.4], 'LineWidth', 3);
legend('Reference Flight Trajectory', 'Terrain(WTPG) + SWH', 'Terrain(WTPG)', 'Location', 'Best');
ylabel('Height (m)');
grid on;
axis([gcaws_time(1), gcaws_time(end), 0, max(logsout{idx_TRNRef}.Values.Altitude.Data)]);
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