
acLLH  = [36.0, 127, 5000];
tarLLH = [36.2, 127, 0];

acECEF = llh2ecef(acLLH);
tarECEF = llh2ecef(tarLLH);

losLLH  = acLLH - tarLLH;
losECEF = acECEF - tarECEF;


ptLLH = zeros(100, 3);
ptLLH(:,1) = tarLLH(1) + linspace(0,1, 100) .* losLLH(1);
ptLLH(:,2) = tarLLH(2) + linspace(0,1, 100) .* losLLH(2);
ptLLH(:,3) = tarLLH(3) + linspace(0,1, 100) .* losLLH(3);

ptLLH_ECEF = zeros(100,3);
for ii = 1:100
    ptLLH_ECEF(ii, :) = llh2ecef(ptLLH(ii,:));
end

ptECEF = zeros(100, 3);
ptECEF(:,1) = tarECEF(1) + linspace(0,1, 100) .* losECEF(1);
ptECEF(:,2) = tarECEF(2) + linspace(0,1, 100) .* losECEF(2);
ptECEF(:,3) = tarECEF(3) + linspace(0,1, 100) .* losECEF(3);


e = referenceEllipsoid('wgs84', 'm');

figure;
hold on;
ellipsoid(0,0,0,e.SemimajorAxis, e.SemimajorAxis, e.SemiminorAxis);
plot3(ptECEF(:,1), ptECEF(:,2), ptECEF(:,3), 'b');
plot3(ptLLH_ECEF(:,1), ptLLH_ECEF(:,2), ptLLH_ECEF(:,3), '.r')
hold off;
axis equal;

figure;
for ii = 1:3
    subplot(3, 1, ii);
    hold on;
    plot(ptECEF(:,ii), '-r');
    plot(ptLLH_ECEF(:,ii), '-b');
    hold off;
    legend('ECEF', 'Geodetic');
end


figure;
for ii = 1:3
    sgtitle('Error');
    subplot(3, 1, ii);
    hold on;
    plot(ptECEF(:,ii) - ptLLH_ECEF(:,ii), '-r');
    hold off;
end


function llh = ecef2llh(xyzECEF)
%ECEF2LLH
% Compute the {LLH} coordinates given {ECEF}
%
% == Inputs =========
% ecef              -  {ECEF} coordinates [m, m, m]
%
% == Outputs ========
% llh               - Corresponding Latitude, Longitude, Height [deg,deg,m]
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

lat = lat * 180/pi;
lon = lon * 180/pi;

llh = [lat; lon; h];
  

end

function ecef = llh2ecef(llh)
%LLH2ECEF
% Compute the {ECEF} coordinates given {Geodetic=LLH}
%
% == Inputs =========
% llh               - Latitude, Longitude, Height [deg, deg, m]
%
% == Outputs ========
% ecef              - Corresponding {ECEF} coordinates [m]
%

ecef = zeros(3,1);

lat = llh(1);
lon = llh(2);
h   = llh(3);

% WGS84 Parameters
r_e = 6378137.0;
e   = 0.08181919;

% Prime vertical radius of curvature
ne = r_e / sqrt( 1 - e^2 * (sind(lat))^2 );

xECEF = (ne+h) * cosd(lat) * cosd(lon);
yECEF = (ne+h) * cosd(lat) * sind(lon);
zECEF = (ne * (1-e^2) + h ) * sind(lat);

% Output
ecef = [xECEF; yECEF; zECEF];

end
