classdef motor
    % MOTOR class is a class for a solid rocket motor 
    %   This motor has a circular grain geometry and these calculations are
    %   in the inertial frame. Used to calculate specific burn
    %   charachteristics
    %
    % FIXED PROPERTIES 
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
    % 
    % DESIGN PROPERTIES
    % N - number of gains
    % ep - area ratio
    % ep_0 - initial area ratio
    % c_star - charachteristic velocity
    % T_b - inital propellant temperature
    % 
    % R1 - Grain inner diameter
    % R0 - Grain outer diameter
    % L0 - Grain length
    %% MOTOR PROPERTIES
    properties
        % These are fixed properties of propellant
        a_0 = 0.03; % in/s [lbf/in^2]^-n
        n = 0.35;
        sigma_p = 0.001;% -/F
        gamma_p = 1.25;
        gamma_a = 1.4;
        rho_p = 0.065; % lbm/ft^3
        mw_air = 28.97; % lbm/lbmole
        g_e = 32.2; %ft/s^2
        Mu = 28.97; % [g/mol] air
        R_u = 1545.3; % [ft*lbf / lbmol R]
        
        % need to move these
        R_bar = R_u / Mu;
        R_air = R_bar*g_e;

        % Design parameters
        N = 1; 
%         ep = 4;
        ep_0 = 4;
        c_star = 5210; % ft/s
        T_b = 70; % degrees F
        % R1, R0, L0, numGrains, ep_0, A_t0, mBal, initPropTemp
        R1 = 1; % in, inital grain radius
        R0 = 2.375; % in, final grain radius
        L0 = 8; % in, inital grain length

    end
    %% MOTOR FUNCTIONS
    methods
        function obj = motor(inputArg1,inputArg2)
            % MOTOR Construct an instance of this class
            %   Able to modify properties of motor with this function
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function B = MotorPerformace(obj,inputArg)
            % MotorPerformance Calculates motor performance 
            %   This function takes in design properties to calculate a
            %   burn profile matrix
            outputArg = obj.Property1 + inputArg;
        end
    end
end

