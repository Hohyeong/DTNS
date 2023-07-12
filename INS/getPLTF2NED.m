function [C_pn, vec_n] = getPLTF2NED(vec_p, alpha)
%getPLTF2NED Convert a vector in INS Platform axis to local NED
%navigation frame
%
% Inputs
%   vec_p   : a vector in INS platform axis (3x1)
%   alpha   : wander angle "platform azimuth - true heading" (deg)
%
% Outputs
%   C_pn    : computed tranformation matrix from INS Platform axis to local
% NED navigation frame
%   vec_n   : converted vec_P in local NED navigation frame

% Initialize
C_pn = zeros(3,3);
vec_n = zeros(3,1);

% Transformation matrix form INS Platform axis to local NED navigation
% frame
C_pn = [ cosd(alpha), -sind(alpha),  0; ...
        -sind(alpha), -cosd(alpha),  0; ...
                   0,            0, -1];

% Convert a given vector to NED
vec_n = C_pn * vec_p;

end

% % Testcases
% vec_p = [1; 0; 0]; alpha = 0;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));
% 
% vec_p = [1; 0; 0]; alpha = 30;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));
% 
% vec_p = [1; 0; 0]; alpha = -45;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));
% 
% vec_p = [1; 0; 0]; alpha = 90;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));
% 
% vec_p = [1; 0; 0]; alpha = 180;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));
%
% vec_p = [0; 1; 0]; alpha = 0;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));
% 
% vec_p = [0; 1; 0]; alpha = 30;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));
% 
% vec_p = [0; 1; 0]; alpha = -45;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));
% 
% vec_p = [0; 1; 0]; alpha = 90;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));
% 
% vec_p = [0; 1; 0]; alpha = 180;
% [~, vec_n] = getPLTF2NED(vec_p, alpha);
% fprintf("vec_p: [%f, %f, %f], alpha: %f (deg) \n", vec_p(1), vec_p(2), vec_p(3), alpha);
% fprintf(">> vec_n: [%f, %f, %f] \n", vec_n(1), vec_n(2), vec_n(3));


