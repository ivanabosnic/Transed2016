classdef Current
    %Definition of a current
    %   Defines a current with a logaritmic velocity profile
    %   estimates the depth-averaged current speed from a velocity at a
    %   given depth
    
    properties
        fluid;      % Fluid with the defined properties
        uz;         % Current speed at height z (m*s)
        uAvg;       % Depth-averaged current speed (m*s)
        z;          % Heigth above sea bed (m)
        h;          % Water depth (m)
        roughness;  % Roughness with its defined properties
        uDir;       % Depth-averaged current direction (deg)
        
        
        
        %k=0.4;     % von Karman's constant
    end
    
    methods
        function obj = Current(dim, fluid, uz, uDir, z, h, roughness)
            % dim - dimension
            % fluid - fluid object defining the properties of the fluid
            % uz - current speed at height z (m*s)
            % uDir - depth-averaged current direction (deg)
            % z - heigth above sea bed (m)
            % h - Water depth (m)
            
            if nargin == 1
                if isscalar(dim)
                    obj(dim, 1) = Current; % Preallocate object array
                else
                    n = dim(1);
                    m = dim(2);
                    obj(n, m) = Current;
                end
            elseif nargin >= 6
                obj.fluid = fluid;
                obj.uz = uz;
                obj.z = z;
                obj.roughness = roughness;
                obj.h = h;
                obj.uAvg = obj.uz*(log(obj.h/obj(1,1).roughness(1,1).z0)-1)/log(obj.z/obj(1,1).roughness(1,1).z0); % deduced from eq. 22, pg. 46 [Soulsby, 1997] 
                obj.uDir = uDir;
                
            end
            
        end
        
        function [x, fc] = tau_gm79(obj)
            %compute tau using fc from Grant and Madsen 1979 (pg. 95 Soulsby)
            [n , m] = size(obj);
            x = zeros(n, m);
            fc = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    fc(i,j) = (1/(4*log10(obj(i,j).z./obj(i,j).roughness.z0))).^2;
                    x(i,j) = 1/2*(obj(i,j).fluid.rho * obj(i,j).dragCoef *obj(i,j).uAvg);
                end
                
            end
        end
        
        function Cd = dragCoef(obj)
            % Drag coefficient (adimensional)
            [n , m] = size(obj);
            Cd = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    Cd(i,j) = (0.4./(1+log(obj(i,j).roughness.z0/obj(i,j).h)))^2;% eq.37, pg. 55[Soulsby, 1997] 
                end
            end
        end
        
        function tau = tau(obj)
            % Bed shear-stress (N/m2) 
            [n , m] = size(obj);
            tau = zeros(n,m);
            for i = 1:n
                for j = 1:m    
                    tau(i,j) = obj(i,j).fluid.rho * obj(i,j).dragCoef * obj(i,j).uAvg^2;% eq.30, pg. 53[Soulsby, 1997] 
                end
            end
        end
        
         function tauFc = tauFc(obj)
            % Bed shear-stress (N/m2)
            [n , m] = size(obj);
            tauFc = zeros(n,m);
            for i = 1:n
                for j = 1:m
                    tauFc(i,j) = obj(i,j).fluid.rho * (obj(i,j).uStarFc)^2;
                end
            end
         end
        
              
         function uStarFc = uStarFc(obj)
            % uStar - friction velocity (m/s)
            [n,m] = size(obj);
            uStarFc = zeros(n,m);
            for i=1:n
                for j =1:m
                    uStarFc(i,j) = obj(i,j).uz*sqrt(obj(i,j).fc/2);% eq.10, pg7, [Taborda 2015].
                end
            end
         end
        
        function fc = fc(obj)
            [n,m] = size(obj);
            fc = zeros(n,m);
            for i=1:n
                for j=1:m
                    fc(i,j) = 2 * (vk / log((obj(i,j).z)/obj(i,j).roughness.z0))^2;
                end
            end
        end
                 
         
         function uStar = uStar(obj)
            % uStar - friction velocity (m/s)
            [n,m] = size(obj);
            uStar = zeros(n,m);
            for i=1:n
                for j =1:m
                uStar(i,j) = sqrt(obj(i,j).tau./obj(i,j).fluid.rho);% eq.32, pg. 53[Soulsby, 1997]
                end
            end
         end
         
                
        function uAtZ=uAtZ(obj, Z)
            % Current velocity at Z heigth above sea bed (m)
            [n , m] = size(obj);
            uAtZ = zeros(n,m);
            for i = 1:n
                for j=1:m
                    uAtZ(i,j) = obj(i,j).uStar./0.4*log(Z/obj(i,j).roughness.z0);% eq.22, pg. 46 [Soulsby, 1997]
                end
            end
        end  
        
        function chezyCoeficient = chezyCoef(obj)
            % Chezy coeficient, see pg. 53
            [n , m] = size(obj);
            chezyCoeficient = zeros(n, m);
            for i=1:n
                for j= 1:m
                    chezyCoeficient(i,j) = sqrt(g/obj(i,j).dragCoef); % eq.31, pg. 53[Soulsby, 1997] 
                end
            end
        end
        
        function manningCoeficient = manningCoef(obj)
            % Manning coeficient, see pg. 53
            [n , m] = size(obj);
            manningCoeficient = zeros(n, m);
            for i=1:n
                for j= 1:m
                    manningCoeficient(i,j) = sqrt(obj(i,j).h.^(1/3).*obj(i,j).dragCoef/g); % eq.31, pg. 53[Soulsby, 1997] 
                end
            end
        end
        
        function shields = shields(obj)
            % shields - Shields parameter, dimensionless form of bed
            % shear-stress
            % (adimensional)
            [n , m] = size(obj);
            shields = zeros(n, m);
            for i=1:n
                for j= 1:m
                    shields(i,j) = obj(i,j).tau/(g*(obj(i,j).roughness.particle.rho - obj(i,j).fluid.rho)*obj(i,j).roughness.particle.dn); % eq.2a, pg. 9[Soulsby, 1997] 
                end
            end
        end
    end
end

