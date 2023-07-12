%WTPG_SIM_AFTERPLOT

%% IDENTIFY
dataName = logsout.getElementNames;
numData = length(dataName);
for ii = 1:numData
    if strcmp(dataName{ii}, 'TRNOutput')
        idx_TRNOutput = ii;
    elseif strcmp(dataName{ii}, 'WTPGOutput')
        idx_WTPGOutput = ii;
    elseif strcmp(dataName{ii}, 'WTPGDatapump')
        idx_WTPGDatapump = ii;
    end
end

timeStr = datestr(now, 'yyyy-mm-dd_HHMMSS');
header = ['logs\', timeStr];
mkdir(header);

%%
% Worst-case Terrain profile
wtpg_wtp = logsout{idx_WTPGOutput}.Values.wtp.Data;
wtpg_time = logsout{idx_WTPGOutput}.Values.wtp.Time;


%% PLOT
wtpgFig = figure('WindowStyle', 'docked');
tiledlayout(4, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
tfValues = categorical({'TRN Valid', 'Correction Valid', 'TRN Off Map', 'TRN Mode', 'WTPG Valid', 'WTPG Off Map', 'Feasibility'});
tfValues = reordercats(tfValues, {'TRN Valid', 'Correction Valid', 'TRN Off Map', 'TRN Mode', 'WTPG Valid', 'WTPG Off Map', 'Feasibility'});

for ii=1:26

    projTraj = reshape(logsout{idx_WTPGDatapump}.Values.predictedTrajectory.Data(:, 1, :, ii), size(logsout{idx_WTPGDatapump}.Values.predictedTrajectory.Data(:, 1, :, ii), 1, 3));
    innerTraj = reshape(logsout{idx_WTPGDatapump}.Values.predictedTrajectory.Data(:, 2, :, ii), size(logsout{idx_WTPGDatapump}.Values.predictedTrajectory.Data(:, 2, :, ii), 1, 3));
    outerTraj = reshape(logsout{idx_WTPGDatapump}.Values.predictedTrajectory.Data(:, 3, :, ii), size(logsout{idx_WTPGDatapump}.Values.predictedTrajectory.Data(:, 3, :, ii), 1, 3));

    sgtitle(sprintf("TC-%03d",ii));
    if all(projTraj(:,1)<90)
        nexttile(1, [4, 3])
        geoplot(projTraj(:,1), projTraj(:,2), '--k', innerTraj(:,1), innerTraj(:,2), '-r', outerTraj(:,1), outerTraj(:,2), '-r');
        legend('Projected path', 'Predicted Boundary');
    else
        geoplot(37.5, 127.5, 'ro');
    end
    gx = gca;
    gx.LatitudeAxis.TickValuesMode = 'manual';
    gx.LatitudeAxis.TickValues = min(gx.LatitudeAxis.TickValues):5*0.000833333:max(gx.LatitudeAxis.TickValues);
    gx.LatitudeAxis.TickLabelFormat = 'dd';
    gx.LongitudeAxis.Visible = 'on';
    gx.LongitudeAxis.TickValuesMode = 'manual';
    gx.LongitudeAxis.TickValues = min(gx.LongitudeAxis.TickValues):5*0.000833333:max(gx.LongitudeAxis.TickValues);
    gx.LongitudeAxis.TickLabelFormat = 'dd';
    drawnow;

    nexttile(4, [2 3])
    plot(0:0.5:20, wtpg_wtp(ii,:));
    title('WTP');
    xlabel('Time(sec)'); ylabel('Height(m)');
    grid on;
    drawnow;

    nexttile(16, [2 3])

    TRN_valid = logsout{idx_TRNOutput}.Values.TRN_valid.Data(ii);
    AircraftPositionCorrectionValid = logsout{idx_TRNOutput}.Values.AircraftPositionCorrectionValid.Data(ii);
    AircraftOffTerrainMap = logsout{idx_TRNOutput}.Values.AircraftOffTerrainMap.Data(ii);
    NavigationMode = logsout{idx_TRNOutput}.Values.NavigationMode.Data(ii);

    wtpgValid = logsout{idx_WTPGOutput}.Values.wtpgValid.Data(ii);
    wtpgOffTerMap = logsout{idx_WTPGOutput}.Values.wtpgOffTerMap.Data(ii,:);
    isFeasible = logsout{idx_WTPGDatapump}.Values.isFeasible.Data(ii,:);
   

    bar(tfValues, [TRN_valid, AircraftPositionCorrectionValid, AircraftOffTerrainMap, NavigationMode, wtpgValid, wtpgOffTerMap, isFeasible]);
    drawnow


    input(sprintf("Current TC:%03d, next?", ii));
    saveas(wtpgFig, sprintf("%s\\WTPG_TC-%03d%s%s", header, ii, '.png'));

    
end

