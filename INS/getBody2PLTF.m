function [C_bp, vec_p] = getBody2PLTF(vec_b, alpha, rpy)
%getBody2PLTF Convert a vector in body axis to INS Platform axis
%
% Inputs
%   vec_b   : a vector in body axis (3x1)
%   alpha   : wander angle "platform azimuth - true heading" (deg)
%   rpy     : roll, pitch, yaw(true heading) angles (deg)
%
% Outputs
%   C_bp    : computed tranformation matrix from body to INS Platform axis 
%   vec_p   : converted vec_b in INS platform axis

% Initialize
C_bp = zeros(3,3);
vec_p = zeros(3,1);

% Arrange
cRol = cosd(rpy(1));
sRol = sind(rpy(1));
cPit = cosd(rpy(2));
sPit = sind(rpy(2));
cYaw = cosd(rpy(3));
sYaw = sind(rpy(3));

% Transformation matrix form INS Platform axis to local NED navigation
% frame
C_np = [ cosd(alpha), -sind(alpha),  0; ...
        -sind(alpha), -cosd(alpha),  0; ...
                   0,            0, -1];

% Transformation matrix from local NED to body
C_bn   = [ cPit*cYaw, sRol*sPit*cYaw-cRol*sYaw, cRol*sPit*cYaw+sRol*sYaw;...
           cPit*sYaw, sRol*sPit*sYaw+cRol*cYaw, cRol*sPit*sYaw-sRol*cYaw;...
               -sPit,                sRol*cPit,                cRol*cPit];

C_bp = C_np * C_bn;

% Convert a given vector to platform axis
vec_p = C_bp * vec_b;

end

