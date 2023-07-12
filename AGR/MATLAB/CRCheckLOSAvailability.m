function [validity, losAvail, losAvailValid, testLLH, probe] = CRCheckLOSAvailability(targetECEF, deltaECEF, nSearch)

if (exist("targetECEF") && exist(deltaECEF) && exist(nSearch))

else
    targetECEF = [-3.19e+06, 4.13e+06, 3.65e+06];
    deltaECEF = [-15.4, 20, -16.1];
    nSearch = 105;
end

% ipath = ['-I' fullfile(pwd, 'MAP')];
% mapdLib = fullfile(pwd, 'MAP', 'MAPd.lib');
% ws2_32Lib = fullfile(pwd, 'MAP', 'WS2_32.lib');
% mapC = fullfile(pwd, 'MAP', 'MAP.c');
% dvofC = fullfile(pwd, 'MAP', 'DVOF.c');
% userDefinedC = fullfile(pwd, 'MAP', 'userDefined.c');
% 
% ipath = ['-I' fullfile(pwd)];
% mapdLib = fullfile(pwd, 'MAPd.lib');
% ws2_32Lib = fullfile(pwd, 'WS2_32.lib');
% mapC = fullfile(pwd, 'MAP.c');
% dvofC = fullfile(pwd, 'DVOF.c');
% userDefinedC = fullfile(pwd, 'userDefined.c');
% 
% mex('-v', ipath, mapC);

% Initialization
validity        = false;
losAvail        = false;
losAvailValid   = false;
isMapCached_int = int32(0);
testLLH         = zeros(3, 1);
testTerElev     = 0.0;
testObsElev     = 0.0;
testElev        = 0.0;
index           = 0;
latOut          = 0.0;
lonOut          = 0.0;
AMSLI           = 0;
AMSL            = 0;
probe           = 0;

for ii=0:nSearch
    % 1. Propagate search vector
    testECEF = targetECEF + (ii * deltaECEF); 

    % 2. Convert test ECEF coordinate to LLH
    testLLH = ecef2llh(testECEF);

    % 3-1. Check Map data is cached
    isMapCached_int = coder.ceval('IsScanBoundary', testLLH(1), testLLH(2));
    if ( isMapCached_int == 0 )
        losAvailValid = false;
        probe = ii;
        return
    end

    % 3-2. Get Terrain Elevation
    testTerElev = coder.ceval('GetElevation', testLLH(1), testLLH(2));

    % 3-3. Get Obstacle Info
    index = coder.ceval('GetObstTotalCount');
    AMSL  = 0;
    if (index >= 0)
        for jj=0:index-1
            latOut = coder.ceval('GetObstacleInfo', testLLH(1), testLLH(2), jj, 1);
            lonOut = coder.ceval('GetObstacleInfo', testLLH(1), testLLH(2), jj, 2);
            dist=sqrt(((testLLH(1)-latOut)*pi/180*6378137)^2+((testLLH(2)-lonOut)*pi/180*6378137)^2);
            if (dist<50)
                AMSLI = coder.ceval('GetObstacleInfo', testLLH(1), testLLH(2), jj, 4);
                if (AMSLI>AMSL)
                    AMSL=AMSLI;
                end
            end
        end
    end
    testObsElev   = AMSL;

    % 2-3. Get Highest Profile
    if ( testTerElev > testObsElev )
        testElev = testTerElev;
    else
        testElev = testObsElev;
    end

    % 3. Test LOS Availability
    if (testElev > testLLH(3))
        losAvail = false;
        probe = ii;
        return
    else
        
    end

end

validity = true;
losAvail = true;
losAvailValid = true;

end

function llh = ecef2llh(xyzECEF)
%ECEF2LLH
% Compute the {LLH} coordinates given {ECEF}
%
% == Inputs =========
% ecef              - {ECEF} coordinates [m, m, m]
%
% == Outputs ========
% llh               - Corresponding Latitude, Longitude, Height [deg,deg,m]
%

% Initialization
llh = zeros(3,1);

% WGS84 Parameters
r_e     = 6378137.0;
f       = 1/298.257223563;
r_b     = r_e * (1-f);
e       = 0.08181919;

% For convenience
x       = xyzECEF(1); 
y       = xyzECEF(2); 
z       = xyzECEF(3);

% Longitude
lon = atan2(y, x);

% Iteration for Latitude and Height
eps     = 1;
tol     = 1e-10;
p       = sqrt(x^2 + y^2);
lat     = atan(z/(p*(1-e^2)));
h       = 0.0;

while (eps > tol)
    N   = r_e^2 / sqrt(r_e^2*cos(lat)^2 + r_b^2*sin(lat)^2);
    h   = p / cos(lat) - N;
    lat0 = lat;
    lat = atan(z / ( p * ( 1 - e^2*N/(N+h) ) ) );
    eps = abs(lat - lat0);
end

% Parse outputs
llh = [lat*(180/pi); lon*(180/pi); h];


end