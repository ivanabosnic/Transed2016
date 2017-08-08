classdef Wave
    %WAVE class definition 
    %   Includes several methods to compute wave parameters
    %   (e.g. wave length, energy and power)
    
    properties
        H = 0;
        T = 0;
        Dir = 0;
        DirBat = NaN;
        waveAngle = 0;
        depth = Inf;
        index = 1;
        x;
        y;
        breakPoint = false;
        gama = 0.78;
        Uw; %Ubm passou a se chamar Uw - notação Soulsby% Ivana%
        k;
        Date;
        nu = 1.1585e-06;
        rho = 1027;
        Hpar = 'Hs';
        Tpar = 'Tp';
        Dirpar = 'DirM';
        roughness;
        
        
    end
    
    methods
        function obj = Wave(dim, H, T, Dir, depth, DirBat, Uw, roughness)
            if nargin == 1
                if isscalar(dim)
                    obj(dim, 1) = Wave; % Preallocate object array
                else          
                    n = dim(1);
                    m = dim(2);
                    obj(n, m) = Wave;
                end
            elseif nargin >= 4
                obj.H = H;
                obj.T = T;
                obj.Dir = Dir;
                if nargin >= 5
                    obj.depth = depth;
                    if nargin >= 6
                        obj.DirBat = DirBat;
                        
                        if nargin >= 8
                            obj.roughness = roughness;
                            
                        end
                    end
                end
                obj.waveAngle = obj.DirBat - obj.Dir;                          
                if exist('Uw', 'var')
                    obj.Uw = Uw;
                else 
                    obj.Uw = NaN;
                end
                
            end
        end
        
        function x = L(obj)                        %wavelength
            [n , m] = size(obj);
            x = zeros(n, m);
            for i = 1:n
               for j = 1:m
                  x(i,j) = wl(obj(i,j).T, obj(i,j).depth);
               end
            end
        end
        
        function x = E(obj)                        %wave energy
            [n , m] = size(obj);
            x = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    if strcmp(obj(i,j).Hpar, 'Hs')
                        x(i,j) = 1/8 * obj(i,j).rho * 9.81 * obj(i,j).H.^2 / 2;  %
                    elseif strcmp(obj(i,j).Hpar, 'Hrms')
                        x(i,j) = 1/8 * obj(i,j).rho * 9.81 * obj(i,j).H.^2;  %
                    end
                end
            end
        end
        
        function x = c(obj)
            [n , m] = size(obj);
            x = zeros(n, m);
            for i = 1:n
               for j = 1:m
                  x(i,j) = obj(i,j).L ./ obj(i,j).T;
               end
            end
        end

        function x = cg(obj)
            [n , m] = size(obj);
            x = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    x(i,j) = obj(i,j).c .* obj(i,j).n;
                end
            end
            
        end
        
        function x = Q(obj, k, option) % Fernando, introduzco que sea función de K
            [n , m] = size(obj);
            x = zeros(n, m);
            
            if nargin == 1
                k = 0.39;
            end
            
            for i = 1:n
                for j = 1:m
                    x(i,j) = 0.233*k*obj(i,j).H^(5/2)*sind(obj(i,j).waveAngle * 2); % CEM Eq III-2-7b m3/s, % Fernando, introduzco K en vez de 0.39
                end
            end
            if nargin == 3
                switch option
                    case 'NaN to zero'
                        x(isnan(x)) = 0;
                end
            end
        end
        function x = Q2(obj,k2,vx,option) %Function that estimates the LST in function of wind longshore component
            [n , m] = size(obj);
            vl = obj.current(vx); %Fernando
            E = obj.E; %Fernando
            Cn =obj.cg; % Fernando
            ubm = obj.ubm; % Fernando
            Uw = obj.Uw;
            x = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    x(i,j) = k2*E(i,j)*Cn(i,j)*1*(vl(i,j)/Uw(i,j))/((2650-1027)*9.81*0.6); %m3/s Fernando; Komar, P.D., 1998. Beach processes and sedimentation
                end
            end
            if nargin == 4
                switch option
                    case 'NaN to zero'
                        x(isnan(x)) = 0;
                end
            end
        end
       
        function wts = toWaveTimeSeries(obj, j)
            wts = WaveTimeSeries;
            if nargin == 1
                j = 1;
            end
            [n , ~] = size(obj);
            wts.n = n;
            wts.Index = (1:n)';
            for i = 1:n
                 wts.H(i, 1) = obj(i, j).H;
                 wts.T(i, 1) = obj(i, j).T;
                 wts.Dir(i, 1) = obj(i, j).Dir;
                 wts.Date(i, 1) = obj(i, j).Date;
                 wts.depth(i, 1) = obj(i, j).depth;
                 wts.Index(i, 1) = obj(i,j).index;
            end
           
                    
        end
     
        function x = n(obj)
            [n , m] = size(obj);
            x = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    x(i,j) = ene(obj(i,j).L, obj(i,j).depth);
                end
            end
        end
        
        function x = P(obj)
            x = obj.E.*obj.n.*obj.c;
        end
        
        function obj = updateWaveAngle(obj)
            [n , m] = size(obj);
            for i = 1:n
                for j = 1:m
%                     obj(i,j).waveAngle = obj(i,j).DirBat - obj(i,j).Dir ;
%                     Fernando: I believe that these doesn't work for
%                     shorelines oriented to the east. Should be:
                      obj(i,j).waveAngle = obj(i,j).Dir -obj(i,j).DirBat;
                      
                end
            end
        end

        function waveAtDepth = atDepth(obj, d)
            waveAtDepth = obj;
            [n , m] = size(obj);
            for i = 1:n
                for j = 1:m
                    waveAtDepth(i,j).depth = d; 
  
                    if abs(obj(i,j).waveAngle) >= 90
                        waveAtDepth(i,j).H = NaN;
                        waveAtDepth(i,j).waveAngle = NaN;
                        waveAtDepth(i,j).Dir = NaN;
                    else
                        waveAtDepth(i,j).waveAngle = asind(sind(obj(i,j).waveAngle) *  waveAtDepth(i,j).c /obj(i,j).c);
                        waveAtDepth(i,j).Dir = obj(i,j).DirBat + waveAtDepth(i,j).waveAngle;
                        kr = sqrt(abs(cosd(obj(i,j).waveAngle) / cosd(waveAtDepth(i,j).waveAngle)));
                        ks = sqrt(obj(i,j).cg / waveAtDepth(i,j).cg);
                        waveAtDepth(i,j).H =  obj(i,j).H * kr * ks;
                    end
                end
            end
        end
                   
        function obj = computeUw(obj,spectra)
             if nargin==1
                spectra='mono';
            end
            
            [n , m] = size(obj);
            %obj.Uw = zeros(n, m);
            
            for i = 1:n
                for j = 1:m
                    
                    switch spectra
                        case 'poli' %Based on Soulsby e Smallman (1986)
                            Tn = sqrt(obj(i,j).depth/9.81);
                            t = Tn/obj(i,j).T;
                            A = (6500+(0.56+15.54*t)^6)^(1/6);
                            Urms = 0.25*obj(i,j).H/(Tn*(1+A*t^2)^3);
                            obj(i,j).Uw = Urms*sqrt(2);
                            
                        otherwise
                           obj(i,j).Uw = pi * obj(i,j).H / (obj(i,j).T * sinh(2*pi / obj(i,j).L * obj(i,j).depth));
                    end
                    
                    %obj(i,j).Uw= Uw(i,j);
                    
                end
            end
        end
        
        function [x, fw] = tau_gm79(obj)%compute tau using fw from Grant and Madsen 1979 (pg. 95 Soulsby)
            
            [n , m] = size(obj);
            x = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    
                    A = obj(i,j).Uw * obj(i,j).T/(2*pi);
                    %fw calculation using fixed point iteration
                    y = 0.4; %initial guess
                    es = 0.001; %tolerance
                    ea = 10; %randomly large relative approximate error
                    yold = y;
                    kn = 30*obj(i,j).roughness;
                    
                    while ea > es
                        y = 1/((log10(A/kn)-0.17)-log10(1/yold)+0.24*yold);
                        ea = abs((y-yold)/y)*100;
                        yold = y;
                    end
                    fw =(y/4)^2 ;
                    x(i,j) = 0.5 * obj(i,j).rho*fw*obj(i,j).Uw^2;
                end
            end
        end
        
        function [x, fw] = tau(obj)
            % x - (= tauW) 
            % A -(m) Amplitude of oscillatory bed shear-stress due
            % to waves (N/m2)
            % fw - Wave friction factor for rough bed (adimensional)
            [n , m] = size(obj);
            x = zeros(n, m);
            for i = 1:n
                for j = 1:m
                     
                    A(i,j) = (obj(i,j).Uw * obj(i,j).T)/(2*pi); %Semi-orbital excursion (pg. 77 Soulsby, 1997)
                    Rw(i,j) = (obj(i,j).Uw*A(i,j))/obj(i,j).nu; %Reynolds number Equation 58a from Soulsby 1997 pg. 77
                    B = 0.0521; %value from Soulby 1997 pg 79;
                    N= 0.187; %value from Soulby 1997 pg 79;
                    
                    fwr(i,j) = 1.39*(A(i,j)/obj(i,j).roughness.z0)^-0.52;  % eq.62a, pg. 78 [Soulsby, 1997] 
                    fws(i,j) = B*Rw(i,j)^(-N);
                    fw(i,j) = max(fwr(i,j),fws(i,j));
                end
                    x(i,j) = 0.5 * obj(i,j).rho *fw(i,j)*obj(i,j).Uw^2; % eq. 57, pg. 76 [Soulsby, 1997]
            end
            
        end
        
        function uStar = uStar(obj)
            % uStar - friction velocity (m/s)
            [n,m] = size(obj);
            uStar = zeros(n,m);
            for i=1:n
                for j =1:m
                uStar(i,j) = sqrt(obj(i,j).tau./obj(i,j).rho);% eq.32, pg. 53[Soulsby, 1997]
                end
            end
         end
        
        function A = A(obj) %Semi-orbital excursion (pg. 77 Soulsby, 1997)
            [n , m] = size(obj);
            A = zeros (n,m);
            for i =1:n
                for j=1:m
                    A(i,j) = (obj(i,j).Uw * obj(i,j).T)/(2*pi);
                end
            end
        end
        
        function shields = shields(obj)
            [n,m]=size(obj);
            shields = zeros(n,m);
            for i=1:n
                for j=1:m
                    shields(i,j)=obj(i,j).tau/(g*(obj(i,j).roughness.particle.rho - 1026)*obj(i,j).roughness.particle.dn); % eq.2a, pg. 9[Soulsby, 1997]
                end
            end
        end
	
		
        
        function waveAtBreaking = atBreaking(obj)
            waveAtBreaking = obj;
            [n , m] = size(obj);
            for i = 1:n
                for j = 1:m
                    
                    if abs(obj(i,j).waveAngle) > 90
                        waveAtBreaking(i,j).depth = NaN;
                        waveAtBreaking(i,j).waveAngle = NaN;
                        waveAtBreaking(i,j).breakPoint = NaN;
                        waveAtBreaking(i,j).Dir = NaN;
                        waveAtBreaking(i,j).H = NaN;
                    else
                        % Follows Larson et al, 2010 aproximation
                        % DOI: 10.1061/ASCEWW.1943-5460.0000030
                        alfa = (obj(i,j).c / sqrt(9.81 * obj(i,j).H)) ^ 4 / obj(i,j).n * obj(i,j).gama ^2; % eq (10)
                        lambda = (cosd(obj(i,j).waveAngle) / alfa) ^ (2./5); % eq (12)
                        epsilon = sind(obj(i,j).waveAngle) ^ 2 * lambda; % eq (14)
                        if epsilon > 0 && epsilon < 0.5
                            p = [2.8573 -1.6787 0.5948 0.1649 1];
                            delta = polyval(p, epsilon);
                        else
                            delta = 1;
                        end
                        lambda = lambda * delta;
                        waveAtBreaking(i,j).depth = lambda / 9.81 * obj(i,j).c ^ 2;
                        waveAtBreaking(i,j).waveAngle = asind(sind(obj(i,j).waveAngle) * sqrt(lambda));
                        waveAtBreaking(i,j).breakPoint = true;
                        waveAtBreaking(i,j).Dir = obj(i,j).DirBat + waveAtBreaking(i,j).waveAngle;
                        waveAtBreaking(i,j).H = waveAtBreaking(i,j).depth * waveAtBreaking(i,j).gama;
                    end
                end
            end
            
        end
        
      
       
        function plotWave(obj, auxVar)
            if nargin == 1
                auxVar = 'H';
            end
            plot([obj.x], [obj.(auxVar)])
        end
        
       
    function vl = current(obj,vx) %Fernando, the vx vector is the wind longshore component vector
			[n , m] = size(obj);
            vl = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    vl(i,j)=1.17*sqrt(9.81*obj(i,j).H)*sind(obj(i,j).waveAngle)*cosd(obj(i,j).waveAngle); % (m/s) Komar, P.D., 1998. Beach processes and sedimentation
                    %vl(i,j)=-18.734649+23.3472511*vl(i,j)+0.0948317*vx(i,j)+0.0062995*vx(i,j)^2+(-0.0015610)*vx(i,j)^3; %Empirical correlation between wind longshore component and predicted longshore current
                    vl(i,j)=1.7410457*vl(i,j)+0.1101198*vx(i,j)+0.0043746*vx(i,j)^2+(-0.0015508)*vx(i,j)^3; %Empirical correlation between wind longshore component and predicted longshore current
                end
            end
    end 
        
    function x = Cn(obj)                        %Fernando: Wave phase velocity
            [n , m] = size(obj);
            x = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    x(i,j) = (9.81*obj(i,j).H)^(1/2);  %
                end
            end
        end
        
    end
    
end

function n = ene(L, d)
%ENE Summary of this function goes here
%

if d == Inf
    n = 0.5;
else
    n = 0.5 * (1 + 4 * pi * d ./ L ./ sinh(4 * pi * d ./ L));
end
end

function L = wl(T, d)
%WL compute the wave length 
%   Guan & Hongmei(2005)formulation 
%   Chinese Journal of Oceanology and Limnology

g = 9.81;

L0 = g / 2 / pi * T .^2;

if d == Inf || d > L0/2
    L = L0;
    
else
    
    x = 2 * pi ./ T .* sqrt(d ./ g);
    y = x .* exp(-1.115 * x.^2) + x.^2 .* tanh(1.325 * x);
    L = 2 * pi * d ./ y;
    
end

        
end