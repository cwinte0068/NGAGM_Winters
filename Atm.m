function [Pa,Ta, Rhoa] = Atm(altitude)
% Atm_Conds This gives you atmostpheric conditions at specified altitude
%   Inputs: 
%     altitude in feet
%   Outputs:
%       Pa, pressure in psi
%       Ta, temperature in R
%       Rhoa, density in lbs/ft^3

% Pressure 
if altitude < 83000
    Pa = (-4.272981*10^(-14))*altitude.^3 + 0.000000008060081*altitude.^2 - ...
        0.0005482655*altitude+14.69241;
else 
    Pa = 0;
end

% Temperature 
if altitude < 32809
    Ta = -0.0036*altitude+518.000;
else
    Ta = 399;
end

% Density
if altitude < 82000
    Rhoa = 0.00000000001255*altitude.^2 - 0.0000019453*altitude + 0.07579;
else
    Rhoa = 0;
end
end