classdef Fluid
    % Fluid definition
    % Define th fluid properties  
    
    properties
        nu = 1e-6; % kinematic viscosity of water (m2*s-1*10-6)
        rho = 1027; % density of water (kg*m3) 
    end
    
    methods
        function f = Fluid(T, S)
            % S = salinity    [psu      (PSS-78)]
            % T = temperature [degree C (IPTS-68)]
            if nargin == 2
                %----------------------
                % DEFINE CONSTANTS
                %----------------------
                a0 = 999.842594;
                a1 =   6.793952e-2;
                a2 =  -9.095290e-3;
                a3 =   1.001685e-4;
                a4 =  -1.120083e-6;
                a5 =   6.536332e-9;
                
                rhoPure = a0 + (a1 + (a2 + (a3 + (a4 + a5*T).*T).*T).*T).*T; %density of pure water at atmospherical pressure
                
                %----------------------
                % DEFINE CONSTANTS
                %----------------------
                %  UNESCO 1983 eqn(13) p17.
                
                b0 =  8.24493e-1;
                b1 = -4.0899e-3;
                b2 =  7.6438e-5;
                b3 = -8.2467e-7;
                b4 =  5.3875e-9;
                
                c0 = -5.72466e-3;
                c1 = +1.0227e-4;
                c2 = -1.6546e-6;
                
                d0 = 4.8314e-4;
                
                
                f.rho = rhoPure + (b0 + (b1 + (b2 + (b3 + b4*T).*T).*T).*T).*S  ...
                    + (c0 + (c1 + c2*T).*T).*S.*sqrt(S) + d0*S.^2;
                
                f.nu = 1e-4*(17.91-0.5381*T+0.00694*T.^2+0.02305*S)./f.rho; %kinematic viscosity
            end
            
        end
        
    end
end
