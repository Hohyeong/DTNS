dt = 0.04;
old_llh_b = [0.61435727017776; 2.2383847656853; 900];
old_C_b_n = [0.998647, -1.68136940e-06, 0.05200122; ...
             1.55515e-06, 0.9999999, 2.4676967e-06; ...
             -0.0520012250, -2.383488256e-06, 0.9986470];

 old_v_eb_n = [220; -3.56e-05; 11.4];
 
 f_ib_b = [4.09; 0.00122; -7.62];
 
 w_ib_b = [9.5e-06; -0.0202; -6.11e-06];
 
  
 %% [0] Initialize Outputs
llh_b   = zeros(3,1);
C_b_n   = zeros(3,3);
euler   = zeros(3,1);
v_eb_n  = zeros(3,1);
v_eb_b  = zeros(3,1);
a_b     = zeros(3,1);

%% [1] Define Constant Parameters
w_ie    = 7.2921159E-5;     % Earth rotation rate [rad/s]
R_0     = 6378137;          % WGS84 Equatorial radius [m]
R_P     = 6356752.31425;    %WGS84 Polar radius [m]
e       = 0.0818191908425;  % WGS84 eccentricity
f       = 1 / 298.257223563;%WGS84 flattening
mu      = 3.986004418E14;   %WGS84 Earth gravitational constant [m^3 s^-2]

% Ib      = [inertia(1), -inertia(4), -inertia(5);...
%           -inertia(4),  inertia(2), -inertia(6);...
%           -inertia(5), -inertia(6),  inertia(3)];
%                           % Moment of inertia [kg m^2]

%% [2] ATTITUDE UPDATE

% Calculate w_ib_b: angular rate of body frame w.r.t. ECI frame,
%           resolved about body axes, averaged over time interval [rad/s]
% w_dot_ib_b  = Ib\(M_ib_b - cross(old_w_ib_b, Ib*old_w_ib_b));
% w_ib_b      = old_w_ib_b + 0.5 * dt * (w_dot_ib_b + old_w_dot_ib_b); 

% lumpedS     = [         1,        0;      -sin(euler(2));...
%                         0, cos(euler(
% eulerRate   = 

% First, Calculate attitude increment, magnitude, and skew-symmetric matrix
alpha_ib_b  = w_ib_b * dt;
mag_alpha   = sqrt(alpha_ib_b' * alpha_ib_b);
Alpha_ib_b  = [             0, -alpha_ib_b(3),  alpha_ib_b(2); ...
                alpha_ib_b(3),              0, -alpha_ib_b(1); ...
               -alpha_ib_b(2),  alpha_ib_b(1),             0];

w_ie_n      = w_ie * [cos(old_llh_b(1)); 0; -sin(old_llh_b(1))];

% From (5.44), determine the angular rate of the NED frame
% w.r.t the ECEF frame, resolved about NED
old_N_E     = R_0 * (1-e^2) / sqrt( 1 - (e * sin(old_llh_b(1)))^2 );
old_M_E     = R_0 / ( 1 - (e * sin(old_llh_b(1)))^2 )^1.5;

old_w_en_n = [old_v_eb_n(2) / (old_N_E + old_llh_b(3));...
             -old_v_eb_n(1) / (old_M_E + old_llh_b(3));...
             -old_v_eb_n(2) * tan(old_llh_b(1)) / (old_N_E + old_llh_b(3))];


%% [3] SPECIFIC FORCE FRAME TRANSFORMATION
% Calculate ave_C_b_e: the average body-to-ECEF-frame transformation matrix
%           over the update interval using (5.84) and (5.85)
temp_w_enie_n   = old_w_en_n + w_ie_n;
skew_w_enie_n   = [               0,-temp_w_enie_n(3), temp_w_enie_n(2);...
                   temp_w_enie_n(3),                0,-temp_w_enie_n(1);...
                  -temp_w_enie_n(2), temp_w_enie_n(1),                0];

if (mag_alpha > 1.E-8)
    ave_C_b_n = old_C_b_n * (eye(3) + (1 - cos(mag_alpha)) / mag_alpha^2 ...
        * Alpha_ib_b + (1 - sin(mag_alpha) / mag_alpha) / mag_alpha^2 ...
        * Alpha_ib_b * Alpha_ib_b) ...
        - 0.5 * skew_w_enie_n * old_C_b_n;
else
     ave_C_b_n = old_C_b_n ...
        - 0.5 * skew_w_enie_n * old_C_b_n;
end %if mag_alpha     

% Transform specific force to ECEF-frame resolving axes using (5.86)
f_ib_n = ave_C_b_n * f_ib_b;

%% [4] UPDATE VELOCITY
% Calculate acceleration due to gravity resolved about NED
grav_n = zeros(3,1);

% Calculate surface gravity using the Somigliana model, (2.134)
sinsqL = sin(old_llh_b(1))^2;
g_0 = 9.7803253359 * (1 + 0.001931853 * sinsqL) / sqrt(1 - e^2 * sinsqL);

% Calculate north gravity using (2.140)
grav_n(1) = -8.08E-9 * old_llh_b(3) * sin(2 * old_llh_b(1));

% East gravity is zero
grav_n(2) = 0.0;

% Calculate down gravity using (2.139)
grav_n(3) = g_0 ...
    * (1 - (2 / R_0) * (1 + f * (1 - 2 * sinsqL) + (w_ie^2 * R_0^2 * R_P / mu)) ...
    * old_llh_b(3) + (3 * old_llh_b(3)^2 / R_0^2));

% Calculate Skew-symmetric matrix for Coriolis and transport-rate term
w_CnT   = old_w_en_n + 2 * w_ie_n;
skew_w_CnT  = [       0, -w_CnT(3),  w_CnT(2); ...
               w_CnT(3),         0, -w_CnT(1); ...
              -w_CnT(2),  w_CnT(1),         0];

% From (5.54),
v_eb_n = old_v_eb_n ...
    + dt * f_ib_n ...
    + dt * grav_n ...
    - dt * skew_w_CnT * old_v_eb_n;

%% [5] UPDATE CURVILINEAR POSITION
% Update height using (5.56)
llh_b(3)    = old_llh_b(3) ...
    - dt/2 * (old_v_eb_n(3) + v_eb_n(3));

% Update latitude using (5.56)
llh_b(1)    = old_llh_b(1) ...
    + dt/2 * old_v_eb_n(1) / (old_M_E + old_llh_b(3)) ...
    + dt/2 * v_eb_n(1) / (old_M_E + llh_b(3));

% Calculate meridian and transverse radii of curvature
N_E         = R_0 * (1-e^2) / sqrt( 1 - (e * sin(llh_b(1)))^2 );
M_E         = R_0 / ( 1 - (e * sin(llh_b(1)))^2 )^1.5;

% Update longitude using (5.56)
llh_b(2)    = old_llh_b(2) ...
    + dt/2 * old_v_eb_n(2) / ((old_N_E + old_llh_b(3)) * cos(old_llh_b(1))) ...
    + dt/2 * v_eb_n(2) / ((N_E + llh_b(3)) * cos(llh_b(1)));

%% [6] UPDATE ATTITUDE
% From (5.44), determine the angular rate of the NED frame
% w.r.t the ECEF frame, resolved about NED frame
w_en_n  = [v_eb_n(2) / (N_E + llh_b(3));...
          -v_eb_n(1) / (M_E + llh_b(3));...
          -v_eb_n(2) * tan(llh_b(1)) / (N_E + llh_b(3))];

% Obtain coordinate transformation matrix from the new attitude w.r.t. an
% inertial frame to the old using Rodrigues' formula, (5.73)
if mag_alpha>1.E-8
    C_new_old = eye(3) + sin(mag_alpha) / mag_alpha * Alpha_ib_b +...
        (1 - cos(mag_alpha)) / mag_alpha^2 * Alpha_ib_b * Alpha_ib_b;
else
    C_new_old = eye(3) + Alpha_ib_b;
end %if mag_alpha

% Update attitude using (5.77)
ave_w_ie_n  = w_ie_n + 0.5 * w_en_n + 0.5 * old_w_en_n;
skew_w_ie_n = [            0, -ave_w_ie_n(3),  ave_w_ie_n(2); ...
               ave_w_ie_n(3),              0, -ave_w_ie_n(1); ...
              -ave_w_ie_n(2),  ave_w_ie_n(1),              0];

C_b_n   = (eye(3) - skew_w_ie_n * dt) * old_C_b_n * C_new_old;


%% [7] CALCULATE EULER ANGLES
C_n_b = C_b_n';

euler(1,1) = atan2(C_n_b(2,3),C_n_b(3,3));  % roll
euler(2,1) = - asin(C_n_b(1,3));        % pitch
euler(3,1) = atan2(C_n_b(1,2),C_n_b (1,1));  % yaw

%% [8] CALCULATE BODY VELOCITY
v_eb_b  = C_n_b * v_eb_n;
a_b = f_ib_b + ave_C_b_n' * grav_n;
