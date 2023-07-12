function [C_np, vec_p] = getNED2PLTF(vec_n, alpha)
%getNED2PLTF Convert a vector in local NED navigation frame to INS
%platform axis
%
% Inputs
%   vec_n   : a vector in local NED navigation frame (3x1)
%   alpha   : wander angle "platform azimuth - true heading" (deg)
%
% Outputs
%   C_np    : computed tranformation matrix from INS Platform axis to local
% NED navigation frame
%   vec_p   : converted vec_P in local NED navigation frame

% Initialize
C_np = zeros(3,3);
vec_p = zeros(3,1);

% Transformation matrix from local NED navigation frame to INS Platform
% axis
C_np = [ cosd(alpha), -sind(alpha),  0; ...
        -sind(alpha), -cosd(alpha),  0; ...
                   0,            0, -1];

% C_pn = [ cosd(alpha), -sind(alpha),  0; ...
%         -sind(alpha), -cosd(alpha),  0; ...
%                    0,            0, -1];

% Convert a given vector to INS platform axis
vec_p = C_np * vec_n;

end

% Testcases
% vec_n = [1; 0; 0]; alpha = 0;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));
% 
% vec_n = [1; 0; 0]; alpha = 30;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));
% 
% vec_n = [1; 0; 0]; alpha = -45;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));
% 
% vec_n = [1; 0; 0]; alpha = 90;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));
% 
% vec_n = [1; 0; 0]; alpha = 180;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));
% 
% vec_n = [0; 1; 0]; alpha = 0;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));
% 
% vec_n = [0; 1; 0]; alpha = 30;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));
% 
% vec_n = [0; 1; 0]; alpha = -45;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));
% 
% vec_n = [0; 1; 0]; alpha = 90;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));
% 
% vec_n = [0; 1; 0]; alpha = 180;
% [~, vec_p] = getNED2PLTF(vec_n, alpha);
% fprintf("vec_n: [%f, %f, %f], alpha: %f (deg) \n", vec_n(1), vec_n(2), vec_n(3), alpha);
% fprintf(">> vec_p: [%f, %f, %f] \n", vec_p(1), vec_p(2), vec_p(3));