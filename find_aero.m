function [F, M] = find_aero(M, alpha, x_CG, l, l_N, d, altitude, powered)
% find_aero calculates the aerodynamic lift, drag, and induced moment on
% the missile
%     l = length of missile in ft
%     d = diameter of missile in ft
%     M = Mach Number
%     q = dynamic pressure
%     l_N = length of nosecone
%     a = speed of sound
%     alpha = angle of attack
%     x_CG = center of mass location from nose
%     powered = logical 1 or 0 indicating powered or coasting flight

% Atmospheric Conditions
[Pa, Ta, Rhoa] = Atm(altitude); % psi, R, plb/ft^3
% Mu = 28.97; % [g/mol] air
Mu = 28.97; % [lbm/lbmol] air
R = 1545.3; % [ft*lbf / lbmol R]
R_bar = R / Mu;
gamma = 1.4; 
% a = sqrt(gamma * R_bar * Ta * 32.2) % 32.2 is to account for units
a = sqrt(gamma * R_bar * Ta);
% M = V / a;
q = 0.5*Rhoa*(M*a)^2; % lbf/ft^2

%% Drag Coeficients
C_D0_Friction = 0.053*(l/d)*(M/(q*l))^0.2 % Based on Jerger reference,...
%                        turbulent boundary layer, q in psf, l in ft
A_e = d; % for testing
S_ref = pi/4 * d^2; 
if M < 1
    if powered == true
        C_D0_Base = (1 - A_e/S_ref)*(0.12+0.13*M^2); % Powered
    else
        C_D0_Base = (0.12 + 0.13*M^2); % Coast
    end
    C_D0_Wave = 0; % no wave drag if subsonic
elseif M >= 1
    if powered == true
        C_D0_Base = (1 - A_e/S_ref)*(0.25/M); % Powered
    else
        C_D0_Base = 0.25/M; % Coast
    end
    C_D0_Wave = (1.59 + 1.83/M^2)*(atan(0.5/(l_N/d)))^1.69; % Based on...
    %                                   Bonney Ref, tan^-1 in rad
end
C_D0_Body = C_D0_Friction + C_D0_Base + C_D0_Wave;

%% No Boattail Considered

%% Body Normal Force Prediction
%   alpha = pitch about y axis
%   phi = roll about x axis
%   a = body major axis radius
%   b = body semi-major axis radius

% Coeficient of Normal Force (Z-X Plane, 0 Y component)
a = d/2;
b = a;
phi = 0;
C_N = ((a/b)*cos(phi)^2 + (b/a)*sin(phi)^2)*(abs(sin(2*alpha)*cos(alpha/2)) + 1.3*(l/d)*sin(alpha)^2);
% C_N = ((a/b)*cos(phi)^2 + (b/a)*sin(phi)^2)*((sin(2*alpha)*cos(alpha/2)) + 1.3*(l/d)*sin(alpha)^2);

% Lift/Drag = C_L/C_D 
LoD = @(alpha) (C_N*cos(alpha) - C_D0_Body*sin(alpha))/(C_N*sin(alpha) + C_D0_Body*cos(alpha));
%% Aerodynamic Center Prediction
% For conceptual design purposes, the effect of Mach number on the location
% of x_AC may be ignored
x_AC = l_N * (0.63*(1-sin(alpha)^2) + 0.5*(l/l_N)*sin(alpha)^2); % for wingless body

% NEED TO ADD IN TAIL FINS HERE
% will find x_AC/c_MAC 
% find hinge moment

% Planar surface drag coefficient
% Etc

%% Forces and Moments
% Assuming Normal force acts on the aerodynamic center or center of
% pressure (used synonamounsly in book at this point)
% Normal Force
F_N = q*C_N*(pi/4*d^2);
F_D = q*C_D0_Body*(pi/4*d^2);
M = F_N*(x_AC - x_CG); % Might be wrong direction 
if alpha > 0
    L = F_N*cos(alpha);
    M = [M];
else
    L = -F_N*cos(alpha);
    M = [-M];
end
D = -(F_N*sin(alpha) + F_D);
F = [D, L];
end

