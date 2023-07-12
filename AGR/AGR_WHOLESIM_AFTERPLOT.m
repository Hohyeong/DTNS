%%AGR_WHOLESIM_AFTERPLOT
% TCAF_WHOLESIM 시뮬레이션 종료 후 로그 데이터를 출력한다.
ON_SAVE = false;

%% IDENTIFY

dataName = logsout.getElementNames;
numData = length(dataName);
for ii = 1:numData
    if strcmp(dataName{ii}, 'OUT_TRUE_STATE')
        idx_true = ii;
    elseif strcmp(dataName{ii}, 'OUT_TRNOutput')
        idx_TRNOutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_TRNDatapump')
        idx_TRNDatapump = ii;

    elseif strcmp(dataName{ii}, 'IN_LOSRInput')
        idx_LOSRInput = ii;
    elseif strcmp(dataName{ii}, 'OUT_LOSROutput')
        idx_LOSROutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_LOSRDatapump')
        idx_LOSRDatapump = ii;
    elseif strcmp(dataName{ii}, 'IN_HRInput')
        idx_HRInput = ii;
    elseif strcmp(dataName{ii}, 'OUT_HROutput')
        idx_HROutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_HRDatapump')
        idx_HRDatapump = ii;
    elseif strcmp(dataName{ii}, 'IN_CRInput')
        idx_CRInput = ii;
    elseif strcmp(dataName{ii}, 'OUT_CROutput')
        idx_CROutput = ii;
    elseif strcmp(dataName{ii}, 'OUT_CRDatapump')
        idx_CRDatapump = ii;
    end
end

%% CALCULATE
trn_time    = logsout{idx_TRNOutput}.Values.TRNTimeTag.Time;
losr_time   = logsout{idx_LOSRInput}.Values.LOSRSelected.Time;
hr_time     = logsout{idx_HRInput}.Values.HRSelected.Time;
cr_time     = logsout{idx_CRInput}.Values.CRSelected.Time;

%% PLOT

% Feasibiltiy
feasibilityFig = figure('WindowStyle', 'docked');
tiledlayout(8,1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.TRN_valid.Data, '-b', 'LineWidth', 2);
ylabel('TRN Valid'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
nexttile
plot(trn_time, logsout{idx_TRNOutput}.Values.NavigationMode.Data, '-b', 'LineWidth', 2);
ylabel('TRN Navigation Mode'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
nexttile
plot(losr_time, logsout{idx_LOSRInput}.Values.LOSRSelected.Data, '-b', 'LineWidth', 2);
ylabel('LOSR Selected'); grid on; axis([losr_time(1), losr_time(end), -0.1, 1.1]);
nexttile
plot(losr_time, logsout{idx_LOSRDatapump}.Values.isFeasible.Data, '-r', 'LineWidth', 2);
ylabel('LOSR Feasibility'); grid on; axis([losr_time(1), losr_time(end), -0.1, 1.1]);
nexttile
plot(hr_time, logsout{idx_HRInput}.Values.HRSelected.Data, '-b', 'LineWidth', 2);
ylabel('HR Selected'); grid on; axis([hr_time(1), hr_time(end), -0.1, 1.1]);
nexttile
plot(hr_time, logsout{idx_HRDatapump}.Values.isFeasible.Data, '-r', 'LineWidth', 2);
ylabel('HR Feasibility'); grid on; axis([hr_time(1), hr_time(end), -0.1, 1.1]);
nexttile
plot(cr_time, logsout{idx_CRInput}.Values.CRSelected.Data, '-b', 'LineWidth', 2);
ylabel('CR Selected'); grid on; axis([cr_time(1), cr_time(end), -0.1, 1.1]);
nexttile
plot(cr_time, logsout{idx_CRDatapump}.Values.isFeasible.Data, '-r', 'LineWidth', 2);
ylabel('CR Feasibility'); grid on; axis([cr_time(1), cr_time(end), -0.1, 1.1]);

% LOSR-CR
temp = nonzeros(logsout{idx_LOSROutput}.Values.LOSRTargetLatitude.Data);
minLat = min(temp); maxLat = max(temp);
temp = nonzeros(logsout{idx_LOSROutput}.Values.LOSRTargetLongitude.Data);
minLon = min(temp); maxLon = max(temp);

losrCrFig = figure('WindowStyle', 'docked');
tiledlayout(6,1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetLatitude.Data, '-k', 'LineWidth', 1);
ylabel('LOSR Target Latitude (deg)'); grid on; axis([losr_time(1), losr_time(end), minLat, maxLat]);
nexttile
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetLongitude.Data, '-k', 'LineWidth', 1);
ylabel('LOSR Target Latitude (deg)'); grid on; axis([losr_time(1), losr_time(end), minLon, maxLon]);
nexttile
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRDataValid.Data, '-k', 'LineWidth', 2);
ylabel('LOSR Data Valid'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
nexttile
plot(cr_time, logsout{idx_CROutput}.Values.CRTargetLOSAvailable.Data, '-b', 'LineWidth', 2);
ylabel('CR LOS Available'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
nexttile
plot(cr_time, logsout{idx_CROutput}.Values.CRTargetLOSAvailableValid.Data, '-b', 'LineWidth', 2);
ylabel('CR LOS Available Valid'); grid on; axis([trn_time(1), trn_time(end), -0.1, 1.1]);
nexttile
temp = nonzeros(logsout{idx_CROutput}.Values.CRTargetElevation.Data);
maxElev = max(temp);
plot(cr_time, logsout{idx_CROutput}.Values.CRTargetElevation.Data./ 0.3048, '-b', 'LineWidth', 2);
ylabel('CR Target Elevation (ft)'); grid on; axis([cr_time(1), cr_time(end), 0, (maxElev+200)./ 0.3048]);
% nexttile
% hold on;
% plot(cr_time, terrainH, '-g', 'LineWidth', 2);
% ylabel('Terrain Height (m)'); grid on; axis([cr_time(1), cr_time(end), 0, maxElev+200]);
% plot(trn_time, logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data, '--r');

alt = logsout{idx_TRNOutput}.Values.Altitude_INS.Data - logsout{idx_TRNOutput}.Values.TRNEstOfINSAltitudeErr.Data;
hrFig = figure('WindowStyle', 'docked');
tiledlayout(1,1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile
hold on;
plot(trn_time, alt./ 0.3048, '-b', 'LineWidth', 2);
plot(losr_time, logsout{idx_LOSROutput}.Values.LOSRTargetElevation.Data ./ 0.3048, '-g', 'LineWidth', 2);
plot(hr_time, logsout{idx_HROutput}.Values.HRImpactPointElevation.Data ./ 0.3048, '-r', 'LineWidth', 2);
axis([losr_time(1), losr_time(end), 300./ 0.3048, 2200./ 0.3048]);
ylabel('Height (ft)');
legend('Aircraft', 'Terrain', 'Impact Point Height', 'Location', 'best');
grid on;


hrAllFig = figure('WindowStyle', 'docked');
tiledlayout(5,1, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile
plot(hr_time, logsout{idx_HROutput}.Values.HRElevationValid.Data, '-b');
axis([hr_time(1), hr_time(end), -0.1, 1.1]);
ylabel('HRElevationValid');
nexttile
plot(hr_time, logsout{idx_HROutput}.Values.HRTargetOffMap.Data, '-b');
axis([hr_time(1), hr_time(end), -0.1, 1.1]);
ylabel('HRTargetOffMap');
nexttile
plot(hr_time, logsout{idx_HROutput}.Values.HRSelectedWraparound.Data, '-b');
axis([hr_time(1), hr_time(end), -0.1, 1.1]);
ylabel('HRSelectedWraparound');
nexttile
plot(hr_time, logsout{idx_HROutput}.Values.HRDataTagWrap.Data, '-b');
axis([hr_time(1), hr_time(end), -0.1, 1.1]);
ylabel('HRDataTagWrap');
nexttile
plot(hr_time, logsout{idx_HROutput}.Values.HRImpactPointElevation.Data, '-r');
ylabel('Height (m)');
grid on;

