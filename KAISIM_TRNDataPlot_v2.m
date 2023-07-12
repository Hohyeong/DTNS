%   ******************************************************************************
%   *                                 ACUS Lab.
%   * File Name          : KAISIM_TRNData_plot_v2.m
%   * Description        : DTNS-TRN data plot
%   * Version            : v2
%   * Created on         : 2023. 04. 30
%   * Revised on         : 2023. 05. 10
%   * Author             : HSU
%   ******************************************************************************

%% Data load

global ON_SAVE TIMESTR HEADER

ON_SAVE = true;

TIMESTR = datestr(now, 'yyyy-mm-dd_HHMMSS');
HEADER = ['TRN_Logs\', TIMESTR];
mkdir(HEADER);

load("TRNDataOut.mat");
load("TRNWeightOut.mat");

%%

TRN_Time = Out_TRNData.TRNPosition_error.TRNPositionError_Horizon.Time;

j = 1;

for i = 1 : 1 : length(Out_TRNData.GPS_INSData_AcqCSU_Flag.RALTValid.Time)
    if Out_TRNData.GPS_INSData_AcqCSU_Flag.RALTValid.Data(i) == 0
        RALTinValid_Flag(j) = Out_TRNData.GPS_INSData_AcqCSU_Flag.RALTValid.Time(i);

        j = j + 1;

    end
end

if ~exist('RALTinValid_Flag', 'var')
    RALTinValid_Flag = 0;
end

%% Plot

% 1 TRN Position error + Roughness + Terrain height (4 by 1)
TRN_EstPosErrFig = figure('WindowStyle', 'docked');
tiledlayout(4,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.TRNPosition_error.TRNPositionError_Horizon.Data, '-b');
title('TRN Position Error [Horizon]')
xlabel('Time'); ylabel('Horizon [m]'); grid on;
yyaxis right
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.NavigationMode.Data)
ylabel('Navigation Mode')
yticks([0 1])

nexttile
plot(TRN_Time, Out_TRNData.TRNPosition_error.TRNPositionError_Vertical.Data)
title('TRN Position Error [Vertical]')
xlabel('Time'), ylabel('Vertical [m]')
grid on
yyaxis right
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.NavigationMode.Data)
ylabel('Navigation Mode')
yticks([0 1])

nexttile
plot(TRN_Time, Out_TRNData.TRNPosition_error.Roughness.Data)
title('Roughness')
xlabel('Time'), ylabel('[%]')
grid on

nexttile
plot(TRN_Time, Out_TRNData.TRNPosition_error.hDEM.Data)
title('Terrain height (TRN)')
xlabel('Time'), ylabel('[m]')
grid on

% 2 GPS Flag (4 by 1)
TRN_GPSFlagFig = figure('WindowStyle', 'docked');
tiledlayout(4,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.GPSDataNotValid.Data)
title('GPS Data Not Valid')
xlabel('Time'), ylabel('Flag')
grid on

nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.P_code.Data)
title('GPS P-code Flag')
xlabel('Time'), ylabel('Flag')
grid on

nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.C_code.Data)
title('GPS C-code Flag')
xlabel('Time'), ylabel('Flag')
grid on

nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.FigureofMerit.Data)
title('GPS Figure of Merit')
xlabel('Time'), ylabel('FoM')
grid on

% 3 INS-Acq Flag (4 by 1)
TRN_INSAcqFlagFig = figure('WindowStyle', 'docked');
tiledlayout(4,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.INSDataNotValid.Data)
title('INS Data Not Valid')
xlabel('Time'), ylabel('Flag')
grid on

nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.NavigationMode.Data)
title('Navigation Mode')
xlabel('Time'), ylabel('Flag')
grid on

nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.profileCnt.Data)
title('Acq CSU - Profile count')
xlabel('Time'), ylabel('Count')
grid on

nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.MAD.Data)
title('Acq CSU - MAD')
xlabel('Time'), ylabel('MAD')
grid on

% 4 etc Flag (4 by 1)
TRN_etcFlagFig = figure('WindowStyle', 'docked');
tiledlayout(4,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.TRNOperatingEnv_Flag.Data)
title('TRN Operating Enviroment Valid')
xlabel('Time'), ylabel('Flag')
grid on

nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.AircraftPositionCorrectionValid.Data)
title('Aircraft Position Correction Valid')
xlabel('Time'), ylabel('Flag')
grid on

nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.RALTValid.Data)
title('RALT Valid')
xlabel('Time'), ylabel('Flag')
grid on

nexttile
plot(TRN_Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.WeightOnWheels.Data)
title('Weight On Wheels')
xlabel('Time'), ylabel('Flag')
grid on

% 5 Attitude (3 by 1)
TRN_AttitudeFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Roll.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Roll.Data(:,2))
title('Roll')
xlabel('Time'), ylabel('[deg]')
grid on
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha', 0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('Ref', 'TRN')

nexttile
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Pitch.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Pitch.Data(:,2))
title('Pitch')
xlabel('Time'), ylabel('[deg]')
legend('Ref', 'TRN')
grid on
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha',0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('Ref', 'TRN')

nexttile
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Yaw.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Yaw.Data(:,2))
title('Yaw')
xlabel('Time'), ylabel('[deg]')
legend('Ref', 'TRN')
grid on
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha',0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('Ref', 'TRN')

% 6 Velocity (3 by 1)
TRN_VelocityFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_X.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_X.Data(:,2))
title('Velocity_X')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

nexttile
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Y.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Y.Data(:,2))
title('Velocity_Y')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

nexttile
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Z.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Z.Data(:,2))
title('Velocity_Z')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

% 7 Acceleration (3 by 1)
TRN_AccelerationFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_X.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_X.Data(:,2))
title('Acceleration_X')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

nexttile
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Y.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Y.Data(:,2))
title('Acceleration_Y')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

nexttile
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Z.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Z.Data(:,2))
title('Acceleration_Z')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

% 8 INS Position (3 by 1)
TRN_INSPosFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lat.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lat.Data(:,2))
title('INS Position Error [Latitude]')
xlabel('Time'), ylabel('[m]')
legend('TRN correcation', 'True/INS error')
grid on

nexttile
plot(TRN_Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lon.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lon.Data(:,2))
title('INS Position Error [Longitude]')
xlabel('Time'), ylabel('[m]')
legend('TRN correcation', 'True/INS error')
grid on

nexttile
plot(TRN_Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_alt.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_alt.Data(:,2))
title('INS Position Error [Altitude]')
xlabel('Time'), ylabel('[m]')
legend('TRN correcation', 'True/INS error')
grid on

% 9 Position (3 by 1)
TRN_PosFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.LLA.Latitude.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.LLA.Latitude.Data(:,2))
plot(TRN_Time, Out_TRNData.LLA.Latitude.Data(:,3))
plot(TRN_Time, Out_TRNData.LLA.Latitude.Data(:,4))
title('Latitude')
xlabel('Time'), ylabel('[Deg]')
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha',0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('True', 'INS', 'GPS', 'TRN')
grid on

nexttile
plot(TRN_Time, Out_TRNData.LLA.Longitude.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.LLA.Longitude.Data(:,2))
plot(TRN_Time, Out_TRNData.LLA.Longitude.Data(:,3))
plot(TRN_Time, Out_TRNData.LLA.Longitude.Data(:,4))
title('Longitude')
xlabel('Time'), ylabel('[Deg]')
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha',0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('True', 'INS', 'GPS', 'TRN')
grid on

nexttile
plot(TRN_Time, Out_TRNData.LLA.Altitude.Data(:,1))
hold on
plot(TRN_Time, Out_TRNData.LLA.Altitude.Data(:,2))
plot(TRN_Time, Out_TRNData.LLA.Altitude.Data(:,3))
plot(TRN_Time, Out_TRNData.LLA.Altitude.Data(:,4))
title('Altitude')
xlabel('Time'), ylabel('[m]')
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha',0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('True', 'INS', 'GPS', 'TRN')
grid on

% 10 RALT Altitude - TRNEstAGL (1 by 1)
TRN_EstAGLFig = figure('WindowStyle', 'docked');
tiledlayout(1,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.LLA.RALT_TRN_True_AGL.Data(:,1))
hold on
grid on
plot(TRN_Time, Out_TRNData.LLA.RALT_TRN_True_AGL.Data(:,2))
plot(TRN_Time, Out_TRNData.LLA.RALT_TRN_True_AGL.Data(:,3))
title('RALT - TRN AGL')
xlabel('Time'), ylabel('[m]')
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha',0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('RALT', 'TRN AGL', 'True AGL', 'Location','southeast')

% 11 Altitude comparison (2 by 1)
TRN_AltComparisonFig = figure('WindowStyle', 'docked');
tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.Altitude_Comparison.Ref_Alt.Data(:,1))
hold on
grid on
plot(TRN_Time, Out_TRNData.Altitude_Comparison.Altitude_INS.Data(:,1))
plot(TRN_Time, Out_TRNData.Altitude_Comparison.TRN_MSL.Data(:,1))
plot(TRN_Time, Out_TRNData.Altitude_Comparison.RALT_hDEM.Data(:,1))
title('Altitude comparison (MSL)')
xlabel('Time'), ylabel('[m]')
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha',0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('True', 'INS', 'TRN', 'RALT+hDEM','Location','southeast')

nexttile
plot(TRN_Time, Out_TRNData.Altitude_Comparison.Ref_hDEM.Data(:,1))
title('Terrain height (True)')
xlabel('Time'), ylabel('[m]')
grid on

% 12 INS-GPS Position Error (3 by 1)
for i = 1 : 1 : length(Out_TRNData.INSGPS_Err.INSGPS_Err.Data)
    if abs(Out_TRNData.INSGPS_Err.INSGPS_Err.Data(i,1)) > 1000
        Out_TRNData.INSGPS_Err.INSGPS_Err.Data(i,1) = 0;
    end

    if abs(Out_TRNData.INSGPS_Err.INSGPS_Err.Data(i,2)) > 1000
        Out_TRNData.INSGPS_Err.INSGPS_Err.Data(i,2) = 0;
    end

end

TRN_LLAErrFig = figure('WindowStyle', 'docked');
tiledlayout(3,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.INSGPS_Err.INSGPS_Err.Data(:,1))
title('INS-GPS Latitude Error')
xlabel('Time'), ylabel('[Deg]')
legend('INS-GPS Latitude Error')
grid on

nexttile
plot(TRN_Time, Out_TRNData.INSGPS_Err.INSGPS_Err.Data(:,2))
title('INS-GPS Longitude Error')
xlabel('Time'), ylabel('[Deg]')
legend('INS-GPS Longitude Error')
grid on

nexttile
plot(TRN_Time, Out_TRNData.INSGPS_Err.INSGPS_Err.Data(:,3))
title('INS-GPS Altitude Error')
xlabel('Time'), ylabel('[m]')
legend('INS-GPS Altitude Error')
grid on

% 13 GPS Uncertainty (2 by 1)
TRN_GPSUncertatintyFig = figure('WindowStyle', 'docked');
tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(TRN_Time, Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,1))
title('GPS Horizon Uncertainty')
xlabel('Time'), ylabel('Horizon [m]')
ylim([min(Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,1))-2 30])
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha',0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('GPS Horizon Uncertainty')
grid on

nexttile
plot(TRN_Time, Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,2))
title('GPS Vertical Uncertainty')
xlabel('Time'), ylabel('Vertical [m]')
ylim([min(Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,2))-2 30])
yyaxis right
xline(RALTinValid_Flag(:), 'r-','Alpha',0.2)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('GPS Vertical Uncertainty')
grid on

% 14 TRN Tracking Weight (1 by 1)

TRNWeight.GPS(1,:) = Out_TRNWeight.weight_GPS.Data(1,1,:);
TRNWeight.TRN(1,:) = Out_TRNWeight.weight_TRN.Data(1,1,:);

TRN_TrackingWeightFig = figure('WindowStyle', 'docked');
tiledlayout(1,1, 'TileSpacing', 'tight', 'Padding', 'compact');
nexttile
plot(Out_TRNWeight.weight_GPS.Time, TRNWeight.GPS(1,:))
hold on
plot(Out_TRNWeight.weight_TRN.Time, TRNWeight.TRN(1,:))
title('TRN Tracking Weight')
xlabel('Time'), ylabel('Weight')
ylim([-0.5 1.5])
legend('GPS weight', 'TRN weight')
grid on

%% .fig .png save
if ON_SAVE
    %     saveas(TRN_EstPosErrFig, [HEADER, '\TRN_EstPosErr_', TIMESTR, '.fig']);
    %     saveas(TRN_GPSFlagFig, [HEADER, '\TRN_GPSFlag_', TIMESTR, '.fig']);
    %     saveas(TRN_INSAcqFlagFig, [HEADER, '\TRN_INSAcqFlag_', TIMESTR, '.fig']);
    %     saveas(TRN_etcFlagFig, [HEADER, '\TRN_etcFlag_', TIMESTR, '.fig']);
    %     saveas(TRN_AttitudeFig, [HEADER, '\TRN_Attitude_', TIMESTR, '.fig']);
    %     saveas(TRN_VelocityFig, [HEADER, '\TRN_Velocity_', TIMESTR, '.fig']);
    %     saveas(TRN_AccelerationFig, [HEADER, '\TRN_Acceleration_', TIMESTR, '.fig']);
    %     saveas(TRN_INSPosFig, [HEADER, '\TRN_INSPos_', TIMESTR, '.fig']);
    %     saveas(TRN_PosFig, [HEADER, '\TRN_Pos_', TIMESTR, '.fig']);
    %     saveas(TRN_EstAGLFig, [HEADER, '\TRN_EstAGL_', TIMESTR, '.fig']);
    %     saveas(TRN_LLAErrFig, [HEADER, '\TRN_LLAErr_', TIMESTR, '.fig']);
    %     saveas(TRN_GPSUncertatintyFig, [HEADER, '\TRN_GPSUncertatinty_', TIMESTR, '.fig']);

    saveas(TRN_EstPosErrFig, [HEADER, '\TRN_EstPosErr_', TIMESTR, '.png']);
    saveas(TRN_GPSFlagFig, [HEADER, '\TRN_GPSFlag_', TIMESTR, '.png']);
    saveas(TRN_INSAcqFlagFig, [HEADER, '\TRN_INSAcqFlag_', TIMESTR, '.png']);
    saveas(TRN_etcFlagFig, [HEADER, '\TRN_etcFlag_', TIMESTR, '.png']);
    saveas(TRN_AttitudeFig, [HEADER, '\TRN_Attitude_', TIMESTR, '.png']);
    saveas(TRN_VelocityFig, [HEADER, '\TRN_Velocity_', TIMESTR, '.png']);
    saveas(TRN_AccelerationFig, [HEADER, '\TRN_Acceleration_', TIMESTR, '.png']);
    saveas(TRN_INSPosFig, [HEADER, '\TRN_INSPos_', TIMESTR, '.png']);
    saveas(TRN_PosFig, [HEADER, '\TRN_Pos_', TIMESTR, '.png']);
    saveas(TRN_EstAGLFig, [HEADER, '\TRN_EstAGL_', TIMESTR, '.png']);
    saveas(TRN_AltComparisonFig, [HEADER, '\TRN_AltComparison(MSL)_', TIMESTR, '.png']);
    saveas(TRN_LLAErrFig, [HEADER, '\TRN_LLAErr_', TIMESTR, '.png']);
    saveas(TRN_GPSUncertatintyFig, [HEADER, '\TRN_GPSUncertatinty_', TIMESTR, '.png']);
    saveas(TRN_TrackingWeightFig, [HEADER, '\TRN_TrackingWeight_', TIMESTR, '.png']);

end
