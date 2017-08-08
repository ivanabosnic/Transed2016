classdef Particle
    %Particle definition
    %   Define a particle using the representative diameters and the
    %   density.
    
    properties
        rho = 2650; % density of the particle (kg*m3)
        di          % intermediate diameter (m) 
        ds          % smaller diameter (m)
        dl          % large diameter (m)
        shape = 'elipsoidal';
    end
        
    
    methods
        function obj = Particle(di, ds, dl, rho, shape)
            % di - intermediate diameter (m)
            % ds - smaller diameter (m)
            % dl - large diameter (m)
            % rho - density of the particle (kg*m3)
            % shape - main shape of the particle: 'elipsoidal' [Default] or
            % 'paralelepipedic' 
            
            if nargin == 1
                obj.ds = di;
                obj.dl = di;
                obj.di = di;
            elseif nargin >= 3
                obj.ds = ds;
                obj.dl = dl;
                obj.di = di;
                if nargin >= 4
                    obj.rho = rho;
                    if nargin == 5
                        obj.shape = shape;
                    end
                end
            end
        end
        
        function dn = dn(obj)
            % dn - Nominal diameter (m); diameter of an equivalent sphere (i.e. a sphere with the same particle volume)
            switch obj.shape
                case 'elipsoidal'
                    dn = (obj.ds*obj.dl*obj.di)^(1/3);
                case 'paralelepipedic'
                    dn = (6*obj.ds*obj.dl*obj.di/pi)^(1/3);
            end
        end
        
        function da = da(obj,fluid)
            %adimensional diameter
            if nargin ==1
                fluid = Fluid;
            end
            s = obj.rho/fluid.rho;
            da = ((g*(s-1))/fluid.nu^2)^(1/3)*obj.dn; % eq. 75, pg. 104[Soulsby, 1997]
        end
        
        function ws = ws(obj,fluid)
            % ws - Settling velocity of isolated sediment grains (m/s)
            if nargin == 1
                fluid = Fluid;
            end
            ws = (fluid.nu/obj.dn)*(((10.36^2)+1.049*obj.da(fluid)^3)^(1/2)-10.36);% eq.102, pg. 134[Soulsby, 1997]
        end
        
              
        function tauCr = tauCr(obj,fluid)
            %tauCr - Threshold bed shear-stress for motion of sediment (N/m2)
            
            if nargin == 1
                fluid = Fluid;
            end
            tauCr = obj.shieldsCr(fluid)*g*(obj.rho-fluid.rho)*obj.dn; % eq.74, pg. 104[Soulsby, 1997] 
        end
        
        function shieldsCr = shieldsCr(obj,fluid)
            %shieldsCr - Threshold shields parameter (adimentional) 
            
            if nargin == 1
                fluid = Fluid;
            end
            shieldsCr = 0.3/(1+1.2*obj.da(fluid))+0.055*(1-exp(-0.020*obj.da(fluid))); % eq.77, pg. 106[Soulsby, 1997] 
        end
        
        function D_star = D_star(obj,fluid)
            %D_star - Dimensionless grain size 
            
            if nargin == 1
                fluid = Fluid;
            end
            s = obj.rho/fluid.rho; % grain density / fluid density 
            D_star = ((g*(s-1)/fluid.nu^2)^(1/3)) * obj.dn; % eq.75, pg. 104[Soulsby, 1997] 
        end
        
        function uCr = uCr(obj,fluid)
             %uCr critical shear velocity   
             
            if nargin == 1
                fluid = Fluid;
            end
            uCr = (obj.tauCr/fluid.rho)^(1/2); % eq.1b, pg. 9[Soulsby, 1997] 
        end
        
    end % methods
    
    
end % class

