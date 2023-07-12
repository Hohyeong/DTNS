clc;
load('KAIinput.mat');
format long

time = In_KAIdata.time.signal1.Data;

EGI = [In_KAIdata.In_True.latitude.Data In_KAIdata.In_True.longitude.Data In_KAIdata.In_True.velX.Data In_KAIdata.In_True.velY.Data In_KAIdata.In_True.velZ.Data ...
       In_KAIdata.In_True.accX.Data In_KAIdata.In_True.accY.Data In_KAIdata.In_True.accZ.Data In_KAIdata.In_True.altitude.Data ...
       In_KAIdata.In_True.roll.Data In_KAIdata.In_True.pitch.Data In_KAIdata.In_True.yaw.Data];

INS = [In_KAIdata.In_INS.Latitude_INS.Data In_KAIdata.In_INS.Longitude_INS.Data In_KAIdata.In_INS.VelocityX.Data ...
       In_KAIdata.In_INS.VelocityY.Data In_KAIdata.In_INS.VelocityZ.Data In_KAIdata.In_INS.AccelerationX.Data In_KAIdata.In_INS.AccelerationY.Data ...
       In_KAIdata.In_INS.AccelerationZ.Data In_KAIdata.In_INS.INSDataNotValid.Data In_KAIdata.In_INS.INS_Alignment_State.Data ...
       In_KAIdata.In_INS.Roll.Data In_KAIdata.In_INS.Pitch.Data In_KAIdata.In_INS.TrueHeading.Data double(In_KAIdata.In_INS.INSTimeTag.Data) In_KAIdata.In_INS.INS_Mech_Flag.Data];

GPS = [In_KAIdata.In_GPS.Latitude_GPS.Data In_KAIdata.In_GPS.Longitude_GPS.Data In_KAIdata.In_GPS.VelocityEast.Data, In_KAIdata.In_GPS.VelocityNorth.Data In_KAIdata.In_GPS.VelocityUp.Data ...
       In_KAIdata.In_GPS.Altitude_GPS.Data double(In_KAIdata.In_GPS.GPSEstHorErr.Data) double(In_KAIdata.In_GPS.GPSEstVerticalErr.Data) double(In_KAIdata.In_GPS.FigureofMerit.Data) ...
       In_KAIdata.In_GPS.GPSDataNotValid.Data double(In_KAIdata.In_GPS.GPS_State.Data) double(In_KAIdata.In_GPS.C_code.Data) double(In_KAIdata.In_GPS.P_code.Data) double(In_KAIdata.In_GPS.GPStimeValid.Data)];

RALT = [In_KAIdata.In_RALT.RALTAltitude.Data In_KAIdata.In_RALT.RALTValid.Data double(In_KAIdata.In_RALT.RALTTimeTag.Data)];

MC = In_KAIdata.In_MC.In_MC.Data;

Baro = In_KAIdata.In_Baro.Altitude_INS.Data;

dlmwrite('KAI_input.csv',[time EGI INS GPS RALT MC Baro],'precision','%.8f');