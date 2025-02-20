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
    %% Constant Properties
    properties(Constant)
        gamma_a = 1.4; % Specific heat constant
        mw_air = 28.97; % lbm/lbmole
        g_e = 32.2; %ft/s^2
        Mu = 28.97; % [g/mol] air
        R_u = 1545.3; % [ft*lbf / lbmol R]
    end

    %% MOTOR PROPERTIES
    properties
%         R_bar = R_u / Mu;
%         R_air = R_bar*g_e;
        
        % These are fixed properties of propellant
        a_0 = 0.03; % in/s [lbf/in^2]^-n
        n = 0.35;
        sigma_p = 0.001;% -/F
        gamma_p = 1.25; % Propellant specific heat ratio
        rho_p = 0.065; % Propellant Density lbm/ft^3

        % Design parameters
        N = 1; % Number of grains
        ep_0 = 4; % Nozzle expansion ration, A_e/A_t
        c_star = 5210; % ft/s
        T_b = 70; % degrees F
        R1 = 1; % in, inital grain radius
        R0 = 2.375; % in, final grain radius
        L0 = 8; % in, inital grain length

    end
    %% MOTOR FUNCTIONS
    methods
        function obj = motor(a_0, sigma_p, gamma_p, rho_p, N, ep_0,...
                c_star, T_b, T_b0, R1, R0, L0)
                obj.a_0 = a_0;
                obj.n = n;
                obj.sigma_p = sigma_p;
                obj.gamma_p = gamma_p;
                obj.rho_p = rho_p;
                obj.N = N;
                obj.ep_0 = ep_0;
                obj.c_star = c_star;
                obj.T_b0 = T_b0;
                obj.T_b = T_b;
                obj.R1 = R1;
                obj.R0 = R0;
                obj.L0 = L0;
        end

        function cf = CF(AK, AreaRatio, P1oP3)
            % '*******************************************************
            %     THRUST COEFFICIENT CODE 
            %     AK is Specific Heat Ratio
            %     AreaRatio is the Exit Area over the Nozzle Throat Area
            %     P1oP3 is the Chamber Pressure over the Ambient Pressure
            %     *********************************************************
                StopCriteria = 0.000001;                         %'Percent Rel Error to Stop
                EA = StopCriteria * 1.1;                         %'Doctoring the Stopping Criteria
                AM2 = 1.5;                                       %'Exit Mach Number Guess
                IterNo = 0;
                P3oP1 = 1 / P1oP3;                                %  'Setting the Iteration Counter
            %     Do While EA > StopCriteria And IterNo < 100     %'Loop.
                while (EA > StopCriteria && IterNo < 100)
                    IterNo = IterNo + 1;                         %'Increment Counter.
                    AFUN = (2 + (AK - 1) * AM2 ^ 2) / (AK + 1);
                    BFUN = (AK + 1) / (2 * (AK - 1));
                    CFUN = 1 / AFUN;
                    DFUN = 1 / AM2 ^ 2;
                    DERFUN = ((AFUN) ^ BFUN) * (CFUN - DFUN);       % 'Derivative of Mach # Function
                    FUNFUN = ((1 / AM2) * AFUN ^ BFUN) - AreaRatio;  %'Mach Number Root Function
                    AMOLD = AM2;                                     %'Old Solution
                    AM2 = AM2 - FUNFUN / DERFUN;                     %'New Solution via Newton Rhapson
                    EA = abs((AM2 - AMOLD) / AM2) * 100;             %'Percent Relative Error
                end
            
                P2oP1 = (1 + 0.5 * (AK - 1) * AM2 ^ 2) ^ (-AK / (AK - 1));
                TERM1 = 2 * AK * AK / (AK - 1);
                TERM2 = 2 / (AK + 1);
                TERM3 = (AK + 1) / (AK - 1);
                TERM4 = (AK - 1) / AK;
                cf = (TERM1 * (TERM2 ^ TERM3) * (1 - (P2oP1 ^ TERM4))) ^ 0.5 ...
                    + (P2oP1 - P3oP1) * AreaRatio;
            end
        
%         function [B] = MotorPerformace(obj,inputArg)
%             % MotorPerformance Calculates motor performance 
%             %   This function takes in design properties to calculate a
%             %   burn profile matrix
%             outputArg = obj.Property1 + inputArg;
%         end
        
        function a = temp_sensitivity(obj)
        % calculates the temperature sensitivity_coeficient
            a = obj.a_0*exp(obj.sigma_p*(obj.T_b - obj.T_b0));
        end

        function [B, T, H] = MotorPerformace(obj)
            % Calculates motor performance 
            R_bar = motor.R_u / motor.Mu;
            R_air = R_bar*motor.g_e;
            
            % Define Anonymous Functions
            % Port Area
            A_p = @(w) (pi*(obj.R1+w).^2);
            
            % Burn Area
            A_b = @(w) obj.N*(2*pi*(obj.R1+w).*(obj.L0-2*w)+2*pi*(obj.R0^2-(obj.R1+w).^2)); % in^2
            
            % Propellant Mass
            m_p = @(w) obj.rho_p*obj.N*pi*((obj.R0^2)-((obj.R1+w).^2)).*(obj.L0-2*w); % lb_m
            
            % Chamber Pressure
            a = temp_sensitivity(obj);
            P_c = @(w, A_t) ((a*obj.rho_p*A_b(w)*obj.c_star)./(A_t*32.2)).^(1/(1-obj.n)); % psi

            % Setting Up Analysis (Still working on this)
                % How to model a changing 

            %%%%%%%%%%%%% Web Distance
            dw = 0.01;
            w = 0:dw:(obj.R0 - obj.R1);
            w_last_step = obj.R0 - obj.R1 - w(end);
            w(end+1) = w(end) + w_last_step;
            
            %%%%%%%%%%%%% Burn Area
            A_b_1 = A_b(w);
            
            %%%%%%%%%%%%% Propellant Mass
            m_p_1 = m_p(w);
        end

    end
end

