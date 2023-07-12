function [C_pb, vec_b] = getPLTF2Body(vec_p, alpha, rpy)
%getPLTF2Body Convert a vector in INS Platform axis to body frame
%
% Inputs
%   vec_p   : a vector in INS platform axis (3x1)
%   alpha   : wander angle "platform azimuth - true heading" (deg)
%   rpy     : roll, pitch, yaw(true heading) angles (deg)
%
% Outputs
%   C_pb    : computed tranformation matrix from INS Platform axis to body 
% axis
%   vec_b   : converted vec_p in in body frame

% Initialize
C_pb = zeros(3,3);
vec_b = zeros(3,1);

% Arrange
cRol = cosd(rpy(1));
sRol = sind(rpy(1));
cPit = cosd(rpy(2));
sPit = sind(rpy(2));
cYaw = cosd(rpy(3));
sYaw = sind(rpy(3));

% Transformation matrix form INS Platform axis to local NED navigation
% frame
C_pn = [ cosd(alpha), -sind(alpha),  0; ...
        -sind(alpha), -cosd(alpha),  0; ...
                   0,            0, -1];

% Transformation matrix from local NED to body
C_nb = [               cPit*cYaw,                cPit*sYaw,     -sPit;...
        sRol*sPit*cYaw-cRol*sYaw, sRol*sPit*sYaw+cRol*cYaw, sRol*cPit;...
        cRol*sPit*cYaw+sRol*sYaw, cRol*sPit*sYaw-sRol*cYaw, cRol*cPit];

C_pb = C_nb * C_pn;

% Convert a given vector to body
vec_b = C_pb * vec_p;

end