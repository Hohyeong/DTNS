% lat = logsout{21}.Values.latitude.Data;
% lon = logsout{21}.Values.longitude.Data;
% alt = logsout{21}.Values.altitude.Data;
% 
% tarLat = logsout{15}.Values.CRTargetLatitude.Data;
% tarLon = logsout{15}.Values.CRTargetLongitude.Data;
% tarAlt = double(logsout{2}.Values.CRTargetElevation.Data);
% 
% inputDCX = logsout{19}.Values.LOSRDirectionCosineX.Data;
% inputDCY = logsout{19}.Values.LOSRDirectionCosineY.Data;
% inputDCZ = logsout{19}.Values.LOSRDirectionCosineZ.Data;
%
% rol = logsout{21}.Values.roll.Data;
% pit = logsout{21}.Values.pitch.Data;
% yaw = logsout{21}.Values.yaw.Data;

load('TestFTGSWLOS_Data.mat');


acECEF = lla2ecef([lat, lon, alt]);
tarECEF = lla2ecef([tarLat, tarLon, tarAlt]);

losECEF = tarECEF - acECEF;
losECEF = losECEF ./ (vecnorm(losECEF, 2, 2));

% losECEF를 NED 좌표로 회전하여 표현
% ecef2ned 행렬
losNED = zeros(size(losECEF));
losBody = zeros(size(losECEF));
for ii = 1:size(losECEF, 1)
    R_en = [ -sind(lat(ii))*cosd(lon(ii)), -sind(lat(ii))*sind(lon(ii)), cosd(lat(ii)); ...
        -sind(lat(ii)), cosd(lon(ii)), 0; ...
        -cosd(lat(ii))*cosd(lon(ii)), -cosd(lat(ii))*sind(lon(ii)),-sind(lat(ii))];
    losNED(ii, :) = (R_en * losECEF(ii, :)' )';


%     [latLength, lonLength] = getLengthOfADegree(lat(ii), alt(ii));
%     losNED(ii, 1) = (tarLat(ii) - lat(ii)) * latLength;
%     losNED(ii, 2) = (tarLon(ii) - lon(ii)) * lonLength;
%     losNED(ii, 3) = -(tarAlt(ii) - alt(ii));
% 
%     losNED(ii, :) = losNED(ii, :) / norm(losNED(ii, :));

    cPit = cosd(pit(ii));
    sPit = sind(pit(ii));
    cRol = cosd(rol(ii));
    sRol = sind(rol(ii));
    cYaw = cosd(yaw(ii));
    sYaw = sind(yaw(ii));
    R_nb    = [              cPit*cYaw,                cPit*sYaw,     -sPit;...
              sRol*sPit*cYaw-cRol*sYaw, sRol*sPit*sYaw+cRol*cYaw, sRol*cPit;...
              cRol*sPit*cYaw+sRol*sYaw, cRol*sPit*sYaw-sRol*cYaw, cRol*cPit];

    losBody(ii, :) = (R_nb * losNED(ii, :)' )';


    

end



time = (1:size(losECEF, 1))*0.04;

figure;
tiledlayout(3,3, 'TileSpacing', 'tight', 'Padding', 'tight');
nexttile
plot(time, losBody(:,1) - inputDCX);
title('Differences');
ylabel('Direction Cosine X'); grid on;
axis([time(1), time(end), -1, 1]);
nexttile
plot(time, losBody(:,1));
grid on;
axis([time(1), time(end), -1, 1]);
title('KAIST 자체 연산');
nexttile
plot(time, inputDCX);
grid on;
axis([time(1), time(end), -1, 1]);
axis([time(1), time(end), -1, 1]);
title('FTG SW');
nexttile
plot(time, losBody(:,2) - inputDCY);
ylabel('Direction Cosine Y'); grid on;
nexttile
plot(time, losBody(:,2));
grid on;
nexttile
plot(time, inputDCY);
grid on;
axis([time(1), time(end), -1, 1]);
axis([time(1), time(end), -1, 1]);
nexttile
plot(time, losBody(:,3) - inputDCZ);
grid on;
axis([time(1), time(end), -1, 1]);
ylabel('Direction Cosine Z'); grid on;
xlabel('Time (sec)');
nexttile
plot(time, losBody(:,3));
grid on;
xlabel('Time (sec)');
axis([time(1), time(end), -1, 1]);
nexttile
plot(time, inputDCZ);
grid on;
xlabel('Time (sec)');
axis([time(1), time(end), -1, 1]);

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