classdef WaveCurrent
    %WaveCurrent calculates the non-linear interaction between waves and
    %currents for estimating the maximum bed shear stress
    properties
        wave;
        current;
%       particle;
%       roughness;
             
    end
    
    properties (Constant) 
        %Grant and Madsen 1979 model (table 9, pg. 91, Soulby)
        a = [0.11; 1.95; -0.49; -0.28];
        m = [0.65; -0.22; 0.15; 0.06];
        n = [0.71; -0.19; 0.17; -0.15];
        b = [0.73; 0.4; -0.23; -0.24];
        p = [-0.68; 0.13; 0.24; -0.07];
        q = [1.04; -0.56; 0.34; -0.27];
        I = 0.67;
        J = 0.5;
        %DATA13
%         b = [0.47; 0.69; -0.09; -0.08];
%         p = [-0.53; 0.47; 0.07; -0.02];
%         q = [2.34; -2.41; 0.45; -0.61];
%         J = 8.8;
        
    end
    
    methods
        function obj = WaveCurrent(dim, wave, current)% particle,roughness)
            if nargin == 1
                if isscalar(dim)
                    obj(dim, 1) = WaveCurrent; % Preallocate object array
                else
                    n = dim(1);
                    m = dim(2);
                    obj(n, m) = WaveCurrent;
                end
            elseif nargin >= 2
                [n , m] = size(wave);
                obj(n, m) = WaveCurrent;
                
                for i = 1:n
                    for j = 1:m
                        obj(i,j).wave = wave(i,j);
                        obj(i,j).current = current(i,j);
%                         obj(i,j).particle = particle;
%                         obj(i,j).roughness = roughness;
                    end
                end
            end
        end
        
        
        function x = phi(obj)%%Ângulo entre ondas e correntes (pg. 88 Soulby)
            [n , m] = size(obj);
            x = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    x(i,j) = mod((obj(i,j).wave.Dir)- (obj(i,j).current.uDir),90);
                    
                end
            end
        end
        
        %Maximum bed shear stress
        function x = tau_max(obj)
            [n , m] = size(obj);
            x = zeros(n, m);           
            for i = 1:n
                for j = 1:m
                    [tau_w, fw] = obj(i,j).wave.tau;
                    %[tau_c, Cd] = obj(i,j).current.tau;
                    tau_c = obj(i,j).current.tauFc;
                    Cd = obj(i,j).current.dragCoef;
                    X = tau_c/(tau_c+tau_w);                 
                    A = (obj(i,j).a(1) + obj(i,j).a(2)*abs(cosd(obj(i,j).phi)).^obj(i,j).I) + (obj(i,j).a(3)+obj(i,j).a(4)*abs(cosd(obj(i,j).phi)).^obj(i,j).I)*log10(fw/Cd);
                    M = (obj(i,j).m(1) + obj(i,j).m(2)*abs(cosd(obj(i,j).phi)).^obj(i,j).I) + (obj(i,j).m(3)+obj(i,j).m(4)*abs(cosd(obj(i,j).phi)).^obj(i,j).I)*log10(fw/Cd);
                    N = (obj(i,j).n(1) + obj(i,j).n(2)*abs(cosd(obj(i,j).phi)).^obj(i,j).I) + (obj(i,j).n(3)+obj(i,j).n(4)*abs(cosd(obj(i,j).phi)).^obj(i,j).I)*log10(fw/Cd);
                    
                    Z = 1 + A*(X^M)*(1-X)^N;
                    %x(i,j) = Z*(tau_c+tau_w);
                    x(i,j) = sqrt((obj(i,j).tau_mean+tau_w*abs(cosd(obj(i,j).phi)))^2+(tau_w*abs(sind(obj(i,j).phi)))^2);
                    %Data 13 tau max
                   
                    
                    
                end
            end
        end
        
        %Mean bed shear stress
        function x = tau_mean(obj)
            [n , m] = size(obj);
            x = zeros(n, m);           
            for i = 1:n
                for j = 1:m
                    [tau_w, fw] = obj(i,j).wave.tau;
                    %[tau_c, Cd] = obj(i,j).current.tau;
                    tau_c = obj(i,j).current.tauFc;
                    Cd = obj(i,j).current.dragCoef;
                    X = tau_c/(tau_c+tau_w);                 
                    
                    B = (obj(i,j).b(1) + obj(i,j).b(2)*abs(cosd(obj(i,j).phi)).^obj(i,j).J) + (obj(i,j).b(3)+obj(i,j).b(4)*abs(cosd(obj(i,j).phi)).^obj(i,j).J)*log10(fw./Cd);
                    P = (obj(i,j).p(1) + obj(i,j).p(2)*abs(cosd(obj(i,j).phi)).^obj(i,j).J) + (obj(i,j).p(3)+obj(i,j).p(4)*abs(cosd(obj(i,j).phi)).^obj(i,j).J)*log10(fw./Cd);
                    Q = (obj(i,j).q(1) + obj(i,j).q(2)*abs(cosd(obj(i,j).phi)).^obj(i,j).J) + (obj(i,j).q(3)+obj(i,j).q(4)*abs(cosd(obj(i,j).phi)).^obj(i,j).J)*log10(fw./Cd);
                    Y = X*(1+B*(X^P)*(1-X)^Q);
                    x(i,j) = Y*(tau_c+tau_w);
                end
            end
        end
        
        %Shields maximum
          function shields_max = shields_max(obj)
            [n, m] = size(obj);
            shields_max = zeros(size(obj));
            for i = 1:n
                for j = 1:m 
                   s = obj(i,j).current.roughness.particle.rho/obj(i,j).current.fluid.rho; 
                   shields_max(i,j) = obj(i,j).tau_max/(g*obj(i,j).current.fluid.rho*(s-1)*obj(i,j).current.roughness.particle.dn); % eq (36) Soulsby 1997
                end
            end
          end
          
           %Shields mean
           function shields_mean = shields_mean(obj)
               [n, m] = size(obj);
               shields_mean = zeros(size(obj));
               
               for i = 1:n
                   for j = 1:m
                       s = obj(i,j).current.roughness.particle.rho./obj(i,j).current.fluid.rho;
                       shields_mean(i,j) = obj(i,j).tau_mean/(g*obj(i,j).current.fluid.rho*(s-1)*obj(i,j).current.roughness.particle.dn); % eq (36) Soulsby 1997                      
                   end
               end
           end
           %uStar max (maximum shear velocity)
          function uStar_max = uStar_max(obj)
            [n, m] = size(obj);
            uStar_max = zeros(size(obj));
            
            for i = 1:n
                for j = 1:m
				  uStar_max(i,j) = (obj(i,j).tau_max/obj(i,j).current.fluid.rho).^(1/2);
                end
            end
          end
          
          %uStar mean (mean shear velocity)
          
          function uStar_mean = uStar_mean(obj)
            [n, m] = size(obj);
            uStar_mean = zeros(size(obj));
            
            for i = 1:n
                for j = 1:m
				  uStar_mean(i,j) = (obj(i,j).tau_mean/obj(i,j).current.fluid.rho).^(1/2);
                end
            end
          end
        
          %wave boundary layer thickness (Soulsby 1997; ph 148 eq. 115e)
          function zw = zw(obj)
            [n, m] = size(obj);
            zw = zeros(size(obj));
            
            for i = 1:n
                for j = 1:m
				  zw(i,j) = obj(i,j).uStar_max*obj(i,j).wave.T/(2*pi);
                end
            end
          end
                
    end
end

