classdef motor
    % MOTOR class is a class for a solid rocket motor 
    %   This motor has a circular grain geometry and these calculations are
    %   in the inertial frame
    %
    % PROPERTIES 
    % a_0 - burn rate constant
    % n - propellant burning rate exponent
    % sigma_p - propellant temperature sensitivity
    % gamma_p - propellant specific heat ratio
    % gamma_a - air specific heat ratio
    % rho_p - propellant density
    % mw_air - molecular weight of air
    % R_u - Universsal gas constant
    % g_e - gravitational acceleration constant
    % Mu - molar mass of air
    %% MOTOR PROPERTIES
    properties
        % These are inital conditions
        a_0 = 0.03; % in/s [lbf/in^2]^-n
        n = 0.35;
        sigma_p = 0.001;% -/F
        gamma_p = 1.25;
        gamma_a = 1.4;
        rho_p = 0.065; % lbm/ft^3
        mw_air = 28.97; % lbm/lbmole
        R_u = 1545.3; % ftlbf lbmolR
        g_e = 32.2; %ft/s^2
        
        Mu = 28.97; % [g/mol] air
        R_u = 1545.3; % [ft*lbf / lbmol R]
        
        % need to move these
        R_bar = R_u / Mu;
        R_air = R_bar*g_e;
    end
    %% MOTOR FUNCTIONS
    methods
        function obj = motor(inputArg1,inputArg2)
            %MOTOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

