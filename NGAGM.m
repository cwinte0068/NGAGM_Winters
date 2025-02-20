classdef NGAGM
    % This is the main file where I want to access the subclasses that I
    % thought I would have in the motor class
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        obj = motor(a_0, sigma_p, gamma_p, rho_p, N, ep_0,...
                c_star, T_b, T_b0, R1, R0, L0)
        
        function obj = untitled(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
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

