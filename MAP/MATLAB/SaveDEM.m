function [DEMData, DBinfo] = SaveDEM(fileName)
%SAVEDEM
% Read DTED0, DTED1 and DTED2 and save it as .mat file
% Input :
% Name of a dt0, dt1, or dt2 file
% Output :
% DEMData   | (N+1) x (N+1) double. Elevation data. N=3600 for dt2
% DBinfo    | struct which contains DB information
%
% ex) DEM Matrix arrangement for DEM_N36_E128
%     | [ N36, E128 ] ... [ N36, E129]  |
%     |     ...               ...       |
%     | [ N37, E128 ] ... [ N37, E129]  |
%
% * Maintainer : Sungjoong Kim
% =========================================================================

% Call .dt* parser
[data, lat, lon] = readDT(fileName);

% Create DEM Data
DEMData     = data;

% Remove anomalies in DEM Data
for ii = 1:length(DEMData(:,1))
    for jj = 1:length(DEMData(1,:))
        if DEMData(ii,jj) > 9000
            DEMData(ii,jj) = 00;
        end
    end 
end

% Create DB Info data
DBinfo.GridDegDistMap       = (lat(end) - lat(1)) / (length(lat) -1);

DBinfo.minLat   = min(lat);
DBinfo.maxLat   = max(lat);
DBinfo.minLon   = min(lon);
DBinfo.maxLon   = max(lon);

DBinfo.initGridDegLonMap    = DBinfo.minLon;
DBinfo.initGridDegLatMap    = DBinfo.maxLat;

% Save as .mat file

name = sprintf('MAPDATA/DEM_N%02.f_E%03.f.mat', DBinfo.minLat, DBinfo.minLon);
save(name, 'DEMData', 'DBinfo');


end


function [data, latGrid, lonGrid] = readDT(varargin)

% Check Inpput Arguments
if isempty(varargin)
    [file, path] = ...
        uigetfile({'*.dt0;*.dt1;*.dt2', 'DTED(*.dt0,*.dt1,.*dt2)'},'File Selector');
    if file == 0
        warning('No DT files in the current folder');
        return
    end
    filnavn = strcat(path,file);
    dtedLevel = str2double(filnavn(end));
else
    filnavn = cell2mat(varargin);
    dtedLevel = str2double(filnavn(end));
end

minLatIdx = strfind(filnavn, ['n' + digitsPattern(2) + '_e'] );
minLat = str2double(filnavn(minLatIdx+1:minLatIdx+2));

switch dtedLevel
    case 0  %DTED 0
        if minLat >= 80
            disp('Out of range')
            return
            % N = [];
        elseif (minLat >= 70) && (minLat < 80)
            N = 31;
        elseif (minLat >= 50) && (minLat < 70)
            N = 61;   %                       [50 - 70>
        else
            N = 121;  %                       <-50 - 50>
        end
        M = 121;
        K = 254;
    case 1  %DTED 1
        if minLat >= 80
            N = 201;
        elseif (minLat >= 70) && (minLat < 80)
            N = 301;
        elseif (minLat >= 50) && (minLat < 70)
            N = 601;
        else
            N = 1201;
        end
        M = 1201;
        K = 2414;
    case 2 %DTED 2
        if minLat >= 80
            disp('Out of range')
            return
            %N = [];
        elseif (minLat >= 70) && (minLat < 80)
            N=901;
        elseif (minLat >= 50) && (minLat < 70)
            N = 1801;
        else
            N = 3601;
        end
        M = 3601;
        K = 7214;
end

fil1 = fopen(filnavn, 'rb');

try
    c = fread(fil1, [3428 1], 'uint8=>char');
catch
    disp('No data avaible')
    return
end

lonOrigin = str2double(c(5:7));
latOrigin = str2double(c(13:15));
latGrid = latOrigin:1/(M-1):latOrigin+1;
lonGrid = lonOrigin:1/(N-1):lonOrigin+1;

mtrix = fread(fil1, [K inf], 'uint8');
fclose(fil1);

el = mtrix(9:2*M+8,:);

map = byte_2_usign16(el);

data = map;

end

function data = byte_2_usign16(matrix)

[y x] = size(matrix);
matrix1 = matrix(1:2:end,:);
matrix1_temp =matrix1(:);
matrix2 = matrix(2:2:end,:);
matrix2_temp = matrix2(:);
f =  matrix1_temp * 2^8 + matrix2_temp;
data = reshape( f, y/2,x);

end

