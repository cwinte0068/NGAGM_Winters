function [X] = EoM1(X0,F_thrust, F_aero, M_aero, dt, obj)
% EoM1 Is a trajectory solver for one step at a time
% This intial version is for a 3DoF pitch system
%   X is the state 

% These are in the inertial frame
x = X0(1); % x position
% y = X0(2); % y position = 0
z = X0(2); % z position
theta = X0(3); % angle of attack in the inertial frame
u = X0(4); % x axis velocity
w = X0(5); % z axis velocity
q = X0(6); % pitch axis angular velocity

% Body to inertial rotation matrix
R_BI = [cos(theta),sin(theta);   % * [F_x, F_z ]
        -sin(theta),cos(theta)]; 

% linear acceleration due to thrust, drag, and gravity
% Summing Forces in body frame
F_x_B = F_thrust - F_aero(1);
F_z_B = F_aero(2);
% a = F/m


% angular accelerations due to aero moments
% Find Moment of Inertia
g_c = 32.2;
% W_body = l*structure_density + W_inert; % there is a referece in the book for this
W_body = 2214; % lbm %%%%%%%%%%% TEMP VALUE
I_y = I_y_body + (W_body(x_body-x_cg)^2)/g_c +...
    I_y_motor + (W_motor(x_motor-x_cg)^2)/g_c;

M_y_cg = M_aero;
omega = q + (M_y_cg/I_y)*dt;

end

