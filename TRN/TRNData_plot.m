%   ******************************************************************************
%   *                                 ACUS Lab.
%   * File Name          : TRNData_plot.m
%   * Description        : DTNS-TRN data plot
%   * Version            : 1.0.0
%   * Created on         : 2023. 04. 30
%   * Revised on         : 2023. 04. 30
%   * Author             : HSU
%   ******************************************************************************

%% Data clear && load
% clear; clc; close all;
	
cwd = pwd;
addpath(genpath(cwd));

load("TRNDataOut.mat");

j = 1;

for i = 1 : 1 : length(Out_TRNData.GPS_INSData_AcqCSU_Flag.RALTValid.Time)
    if Out_TRNData.GPS_INSData_AcqCSU_Flag.RALTValid.Data(i) == 0
        RALTValid_Flag(j) = Out_TRNData.GPS_INSData_AcqCSU_Flag.RALTValid.Time(i);
        
        j = j + 1;

    end
end

% BEGIN: Modification before making CD (230501)
if ~exist('RALTValid_Flag', 'var')
    RALTValid_Flag = 0;
end
% END: Modification before making CD (230501)

%% TRN Position error + Roughness + Terrain height (4 by 1)
figure(1)
subplot(4,1,1)
plot(Out_TRNData.TRNPosition_error.TRNPositionError_Horizon.Time, Out_TRNData.TRNPosition_error.TRNPositionError_Horizon.Data)
title('TRN Position Error [Horizon]')
xlabel('Time'), ylabel('Horizon [m]')
grid on
yyaxis right
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.NavigationMode.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.NavigationMode.Data)
ylabel('Navigation Mode')
yticks([0 1])

subplot(4,1,2)
plot(Out_TRNData.TRNPosition_error.TRNPositionError_Vertical.Time, Out_TRNData.TRNPosition_error.TRNPositionError_Vertical.Data)
title('TRN Position Error [Vertical]')
xlabel('Time'), ylabel('Vertical [m]')
grid on
yyaxis right
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.NavigationMode.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.NavigationMode.Data)
ylabel('Navigation Mode')
yticks([0 1])

subplot(4,1,3)
plot(Out_TRNData.TRNPosition_error.Roughness.Time, Out_TRNData.TRNPosition_error.Roughness.Data)
title('Roughness')
xlabel('Time'), ylabel('[%]')
grid on


subplot(4,1,4)
plot(Out_TRNData.TRNPosition_error.hDEM.Time, Out_TRNData.TRNPosition_error.hDEM.Data)
title('Terrain height')
xlabel('Time'), ylabel('[m]')
grid on

%% GPS, INS, Acq CSU, Flag (4 by 3)
figure(2)
subplot(4,3,1)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.GPSDataNotValid.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.GPSDataNotValid.Data)
title('GPS Data Not Valid')
xlabel('Time'), ylabel('Flag')
grid on

subplot(4,3,2)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.INSDataNotValid.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.INSDataNotValid.Data)
title('INS Data Not Valid')
xlabel('Time'), ylabel('Flag')
grid on

subplot(4,3,3)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.TRNOperatingEnv_Flag.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.TRNOperatingEnv_Flag.Data)
title('TRN Operating Enviroment Valid')
xlabel('Time'), ylabel('Flag')
grid on

subplot(4,3,4)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.P_code.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.P_code.Data)
title('GPS P-code Flag')
xlabel('Time'), ylabel('Flag')
grid on

subplot(4,3,5)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.NavigationMode.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.NavigationMode.Data)
title('Navigation Mode')
xlabel('Time'), ylabel('Flag')
grid on

subplot(4,3,6)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.AircraftPositionCorrectionValid.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.AircraftPositionCorrectionValid.Data)
title('Aircraft Position Correction Valid')
xlabel('Time'), ylabel('Flag')
grid on

subplot(4,3,7)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.C_code.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.C_code.Data)
title('GPS C-code Flag')
xlabel('Time'), ylabel('Flag')
grid on

subplot(4,3,8)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.profileCnt.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.profileCnt.Data)
title('Acq CSU - Profile count')
xlabel('Time'), ylabel('Count')
grid on

subplot(4,3,9)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.RALTValid.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.RALTValid.Data)
title('RALT Valid')
xlabel('Time'), ylabel('Flag')
grid on

subplot(4,3,10)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.FigureofMerit.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.FigureofMerit.Data)
title('GPS Figure of Merit')
xlabel('Time'), ylabel('FoM')
grid on

subplot(4,3,11)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.MAD.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.MAD.Data)
title('Acq CSU - MAD')
xlabel('Time'), ylabel('MAD')
grid on

subplot(4,3,12)
plot(Out_TRNData.GPS_INSData_AcqCSU_Flag.WeightOnWheels.Time, Out_TRNData.GPS_INSData_AcqCSU_Flag.WeightOnWheels.Data)
title('Weight On Wheels')
xlabel('Time'), ylabel('Flag')
grid on

%% TRN Uncertainty (2 by 1)
figure(3)
subplot(2,1,1)
plot(Out_TRNData.TRNUncertainty.TRNHorUncertainty.Time, Out_TRNData.TRNUncertainty.TRNHorUncertainty.Data)
title('TRN Horizon Uncertainty')
xlabel('Time'), ylabel('Horizon [m]')
grid on

subplot(2,1,2)
plot(Out_TRNData.TRNUncertainty.TRNVerticalUncertainty.Time, Out_TRNData.TRNUncertainty.TRNVerticalUncertainty.Data)
title('TRN Vertical Uncertainty')
xlabel('Time'), ylabel('Vertical [m]')
grid on

%% Attitude (3 by 1)
figure(4)
subplot(3,1,1)
plot(Out_TRNData.Attitude_Velocity_Acceleration.Roll.Time, Out_TRNData.Attitude_Velocity_Acceleration.Roll.Data(:,1))
hold on
plot(Out_TRNData.Attitude_Velocity_Acceleration.Roll.Time, Out_TRNData.Attitude_Velocity_Acceleration.Roll.Data(:,2))
title('Roll')
xlabel('Time'), ylabel('[deg]')
grid on
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('Ref', 'TRN')

subplot(3,1,2)
plot(Out_TRNData.Attitude_Velocity_Acceleration.Pitch.Time, Out_TRNData.Attitude_Velocity_Acceleration.Pitch.Data(:,1))
hold on
plot(Out_TRNData.Attitude_Velocity_Acceleration.Pitch.Time, Out_TRNData.Attitude_Velocity_Acceleration.Pitch.Data(:,2))
title('Pitch')
xlabel('Time'), ylabel('[deg]')
legend('Ref', 'TRN')
grid on
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('Ref', 'TRN')

subplot(3,1,3)
plot(Out_TRNData.Attitude_Velocity_Acceleration.Yaw.Time, Out_TRNData.Attitude_Velocity_Acceleration.Yaw.Data(:,1))
hold on
plot(Out_TRNData.Attitude_Velocity_Acceleration.Yaw.Time, Out_TRNData.Attitude_Velocity_Acceleration.Yaw.Data(:,2))
title('Yaw')
xlabel('Time'), ylabel('[deg]')
legend('Ref', 'TRN')
grid on
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('Ref', 'TRN')

%% Velocity + Acceleration (3 by 2)
figure(5)

subplot(3,2,1)
plot(Out_TRNData.Attitude_Velocity_Acceleration.Velocity_X.Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_X.Data(:,1))
hold on
plot(Out_TRNData.Attitude_Velocity_Acceleration.Velocity_X.Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_X.Data(:,2))
title('Velocity_X')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

subplot(3,2,3)
plot(Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Y.Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Y.Data(:,1))
hold on
plot(Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Y.Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Y.Data(:,2))
title('Velocity_Y')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

subplot(3,2,5)
plot(Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Z.Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Z.Data(:,1))
hold on
plot(Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Z.Time, Out_TRNData.Attitude_Velocity_Acceleration.Velocity_Z.Data(:,2))
title('Velocity_Z')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

subplot(3,2,2)
plot(Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_X.Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_X.Data(:,1))
hold on
plot(Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_X.Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_X.Data(:,2))
title('Acceleration_X')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

subplot(3,2,4)
plot(Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Y.Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Y.Data(:,1))
hold on
plot(Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Y.Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Y.Data(:,2))
title('Acceleration_Y')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

subplot(3,2,6)
plot(Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Z.Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Z.Data(:,1))
hold on
plot(Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Z.Time, Out_TRNData.Attitude_Velocity_Acceleration.Acceleration_Z.Data(:,2))
title('Acceleration_Z')
xlabel('Time'), ylabel('[m/s]')
legend('Ref', 'TRN')
grid on

%% INS Position - TRN Correcation Error (3 by 1)
figure(6)
subplot(3,1,1)
plot(Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lat.Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lat.Data(:,1))
hold on
plot(Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lat.Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lat.Data(:,2))
title('INS Position Error [Latitude]')
xlabel('Time'), ylabel('[m]')
legend('TRN correcation', 'INS error')
grid on

subplot(3,1,2)
plot(Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lon.Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lon.Data(:,1))
hold on
plot(Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lon.Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_lon.Data(:,2))
title('INS Position Error [Longitude]')
xlabel('Time'), ylabel('[m]')
legend('TRN correcation', 'INS error')
grid on

subplot(3,1,3)
plot(Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_alt.Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_alt.Data(:,1))
hold on
plot(Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_alt.Time, Out_TRNData.TRNCorr_INSErr.TRNCorr_INSErr_alt.Data(:,2))
title('INS Position Error [Altitude]')
xlabel('Time'), ylabel('[m]')
legend('TRN correcation', 'INS error')
grid on

%% Position (3 by 1)
figure(7)
subplot(3,1,1)
plot(Out_TRNData.LLA.Latitude.Time, Out_TRNData.LLA.Latitude.Data(:,1))
hold on
plot(Out_TRNData.LLA.Latitude.Time, Out_TRNData.LLA.Latitude.Data(:,2))
plot(Out_TRNData.LLA.Latitude.Time, Out_TRNData.LLA.Latitude.Data(:,3))
plot(Out_TRNData.LLA.Latitude.Time, Out_TRNData.LLA.Latitude.Data(:,4))
title('Latitude')
xlabel('Time'), ylabel('[Deg]')
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('True', 'INS', 'GPS', 'TRN')
grid on

subplot(3,1,2)
plot(Out_TRNData.LLA.Longitude.Time, Out_TRNData.LLA.Longitude.Data(:,1))
hold on
plot(Out_TRNData.LLA.Longitude.Time, Out_TRNData.LLA.Longitude.Data(:,2))
plot(Out_TRNData.LLA.Longitude.Time, Out_TRNData.LLA.Longitude.Data(:,3))
plot(Out_TRNData.LLA.Longitude.Time, Out_TRNData.LLA.Longitude.Data(:,4))
title('Longitude')
xlabel('Time'), ylabel('[Deg]')
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('True', 'INS', 'GPS', 'TRN')
grid on

subplot(3,1,3)
plot(Out_TRNData.LLA.Altitude.Time, Out_TRNData.LLA.Altitude.Data(:,1))
hold on
plot(Out_TRNData.LLA.Altitude.Time, Out_TRNData.LLA.Altitude.Data(:,2))
plot(Out_TRNData.LLA.Altitude.Time, Out_TRNData.LLA.Altitude.Data(:,3))
plot(Out_TRNData.LLA.Altitude.Time, Out_TRNData.LLA.Altitude.Data(:,4))
title('Altitude')
xlabel('Time'), ylabel('[m]')
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('True', 'INS', 'GPS', 'TRN')
grid on

%% RALT Altitude - TRNEstAGL (1 by 1)
figure(8)
plot(Out_TRNData.LLA.RALT_TRN_AGL.Time, Out_TRNData.LLA.RALT_TRN_AGL.Data(:,1))
hold on
grid on
plot(Out_TRNData.LLA.RALT_TRN_AGL.Time, Out_TRNData.LLA.RALT_TRN_AGL.Data(:,2))
title('RALT - TRN AGL')
xlabel('Time'), ylabel('[m]')
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('RALT', 'TRN AGL','Location','southeast')

%% INS-GPS Position Error (3 by 1)
figure(9)
subplot(3,1,1)
plot(Out_TRNData.INSGPS_Err.INSGPS_Err.Time, Out_TRNData.INSGPS_Err.INSGPS_Err.Data(:,1))
title('Latitude Error')
xlabel('Time'), ylabel('[Deg]')
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('INS-GPS')
grid on

subplot(3,1,2)
plot(Out_TRNData.INSGPS_Err.INSGPS_Err.Time, Out_TRNData.INSGPS_Err.INSGPS_Err.Data(:,2))
title('Longitude Error')
xlabel('Time'), ylabel('[Deg]')
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('INS-GPS')
grid on

subplot(3,1,3)
plot(Out_TRNData.INSGPS_Err.INSGPS_Err.Time, Out_TRNData.INSGPS_Err.INSGPS_Err.Data(:,3))
title('Altitude Error')
xlabel('Time'), ylabel('[m]')
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('INS-GPS')
grid on

%% GPS Uncertainty (2 by 1)
figure(10)
subplot(2,1,1)
plot(Out_TRNData.INSGPS_Err.GPSEstUncer.Time, Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,1))
title('GPS Horizon Uncertainty')
xlabel('Time'), ylabel('Horizon [m]')
ylim([min(Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,1))-2 max(Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,1))+2])
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('INS-GPS')
grid on

subplot(2,1,2)
plot(Out_TRNData.INSGPS_Err.GPSEstUncer.Time, Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,2))
title('GPS Vertical Uncertainty')
xlabel('Time'), ylabel('Vertical [m]')
ylim([min(Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,2))-2 max(Out_TRNData.INSGPS_Err.GPSEstUncer.Data(:,2))+2])
yyaxis right
xline(RALTValid_Flag(:), 'r-','Alpha',0.03)
ylabel('RALT Invalid', 'Color', 'r')
yticks([0 1])
legend('INS-GPS')
grid on
