function llh = ecef2llh(xyzECEF)
%ECEF2LLH
% Compute the {LLH} coordinates given {ECEF}
%
% == Inputs =========
% ecef              -  {ECEF} coordinates [m, m, m]
%
% == Outputs ========
% llh               - Corresponding Latitude, Longitude, Height [rad,rad,m]
%

x       = xyzECEF(1); 
y       = xyzECEF(2); 
z       = xyzECEF(3);

% WGS84 Parameters
r_e     = 6378137.0;
f       = 1/298.257223563;
r_b     = r_e * (1-f);
e       = 0.08181919;
e2      = e*e;

% Longitude
lon = atan2(y, x);

% Iteration for Latitude and Height
eps     = 1;
tol     = 1e-10;
p       = sqrt(x^2 + y^2);
lat     = atan(z/(p*(1-e^2)));

while (eps > tol)
    N   = r_e^2 / sqrt(r_e^2*cos(lat)^2 + r_b^2*sin(lat)^2);
    h   = p / cos(lat) - N;
    lat0 = lat;
    lat = atan(z / ( p * ( 1 - e^2*N/(N+h) ) ) );
    eps = abs(lat - lat0);
end

llh = [lat; lon; h];

end