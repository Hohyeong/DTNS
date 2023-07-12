function R = Terrain_level(Lat, Lon)

lat2m = 110961.7060516;
lon2m = 89476.51124;

d_Lat = 90/lat2m; d_Lon = 90/lon2m;
d = 90;

A=0.0; B=0.0; C=0.0; D=0.0; E=0.0; F=0.0; G=0.0; H=0.0; I=0.0; 

A = coder.ceval('MAPGetTerrainHeight', Lat+d_Lat, Lon-d_Lon);
B = coder.ceval('MAPGetTerrainHeight', Lat+d_Lat, Lon);
C = coder.ceval('MAPGetTerrainHeight', Lat+d_Lat, Lon+d_Lon);
D = coder.ceval('MAPGetTerrainHeight', Lat, Lon-d_Lon);
E = coder.ceval('MAPGetTerrainHeight', Lat, Lon);
F = coder.ceval('MAPGetTerrainHeight', Lat, Lon+d_Lon);
G = coder.ceval('MAPGetTerrainHeight', Lat-d_Lat, Lon-d_Lon);
H = coder.ceval('MAPGetTerrainHeight', Lat-d_Lat, Lon);
I = coder.ceval('MAPGetTerrainHeight', Lat-d_Lat, Lon+d_Lon);

SLP_CR=((C+F+I)-(A+D+G))/(6*d);
SLP_DR=((A+B+C)-(G+H+I))/(6*d);

mean_R=((mean(SLP_CR)-mean(SLP_DR))/2)^2;
std_R=((std(SLP_CR)^2)+(std(SLP_DR)^2))/2;

R=sqrt(mean_R+std_R)*100;

end