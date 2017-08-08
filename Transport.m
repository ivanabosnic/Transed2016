classdef Transport
    %UNTITLED2 Summary of this class goes here
        
    properties
        wavecurrent_sf;
        wavecurrent_tot;
        particle;
        bedshearstress;
        current;
        date;
        beta=0; %slope of bed in streamwise direction (positive flow truns uphill)
        
    end
                 
    methods
        function obj = Transport(dim, wavecurrent_sf, wavecurrent_tot, particle, bedshearstress)
            if nargin == 1
                if isscalar(dim)
                    obj(dim, 1) = Transport; % Preallocate object array
                else
                    n = dim(1);
                    m = dim(2);
                    obj(n, m) = Transport;
                end
            elseif nargin >= 4
                [n , m] = size(wavecurrent_sf);
                obj(n, m) = Transport;
                for i = 1:n
                    for j = 1:m
                        obj(i,j).particle = particle;
                        obj(i,j).wavecurrent_sf = wavecurrent_sf(i,j);
                        obj(i,j).wavecurrent_tot = wavecurrent_tot(i,j);
                        obj(i,j).bedshearstress = bedshearstress(i,j);                        
                        
                    end
                end
            end             
        end
                 
         %Threshold current speed       
         function Ucr = Ucr(obj)
            [n , m] = size(obj);
            Ucr = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    Ucr(i,j) = 0.19*(obj(i,j).particle.dn^0.1)*log10(4*obj(i,j).wavecurrent_sf.current.h/obj(i,j).particle.dn); %for 0.1<=d50<=0.5mm - Van Rijn (1984)
                end                
            end
         end
         
         %%   
         %%%Suspended load%%
         
         %Soulsby-Van Rijn
         
         function qs_svAss = qs_svAss(obj) %Soulsby 1997, pg.183;
            [n , m] = size(obj);
            qs_svAss = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    um = obj(i,j).wavecurrent_sf.current.uAvg;
                    Cd = obj(i,j).wavecurrent_sf.current.dragCoef;
                    Uwc = sqrt(um.^2+0.018/Cd*obj(i,j).wavecurrent.wave.Uw^2);
                    T =Uwc - obj(i,j).Ucr;
                    if T > 0
                        qs_svAss(i,j) = obj(i,j).Ass*um*(T)^2.4*(1-1.6*tand(obj(i,j).beta));
                    else
                        qs_svAss(i,j) = 0;
                    end
                    
                end
                
            end
         end
         
         %General formulation 
         %Reference concentration from Smith and Mclean 1977 OR Van Rijn 1984.                  
       
         function [qss_CS, qss_AS, qss] = qss(obj)
             [n , m] = size(obj);
             qss = zeros(m, n);             
             uz = zeros(100, n);
             cz = zeros(100, n);
             b = zeros(1, n);
             b_max = zeros(1, n);
             T = zeros (1, n);
             zr = zeros (1, n);
             cr = zeros (1, n);
             zw = zeros (1, n);
             zmin = 1e-3; %minimum distance from the bottom
             zmax = 3.6; %maximum distance from the bed; i.e., top of the logarithmic layer;first measurement of the ADCP
             z = logspace(log10(zmin),log10(zmax),100); %velocity profile vector
             dz=[0 diff(z)]';
             
             for i = 1:n
                 for j = 1:m
                     z0s = obj(i,j).wavecurrent_sf.current.roughness.z0;
                     wc_sf = obj(i,j).wavecurrent_sf;
                     wc_tot = obj(i,j).wavecurrent_tot;
                     p = obj(i,j).particle;                     
                     uz(:,i)=(wc_sf.current.uStarFc./0.4)*log10(z./z0s)';%velocity profile vector
                                          
                     if wc_sf.shields_max>0.8 %the bed is expected to be flat with sheet flow, thus total bed shear stress is equal to skin friction in this situation.
                         wc_tot = wc_sf;
                     else
                         wc_tot = wc_tot;
                     end
                     
                     b(1,i) = p.ws(wc_tot.current.fluid)/(vk.*wc_tot.uStar_mean); %Rouse number
                     b_max(1,i) = p.ws(wc_tot.current.fluid)/(vk.*wc_tot.uStar_max); %Rouse number
                      
                     cb = 0.65;  % bed concentration(cb = 1-e; e  = porosity);
                     rc = 2.4e-3; %resuspension coefficient
                     T(1,i) = (wc_sf.tau_max-p.tauCr(wc_tot.current.fluid))/p.tauCr(wc_tot.current.fluid);
                     
                     if T(1,i)>0 && wc_sf.uStar_max>p.ws(wc_tot.current.fluid)
                         %ref. height examples 8.4 soulsby
                         %zr(1,i) = 2*p.dn;                          
                         %Ref. concencentration Van Rijn 1984, Soulsby 1997 pg 140
                         %cr(1,i) = (0.01*p.dn*T(1,i).^(3/2))/zr(1,i)*p.D_star^(0.3);
                         %Ref. concencentration Eq. 111, pg 140 Soulsby 1997(Zyserman and Fredsoe (1994))
                         %cr (1,i) = (0.331*(wc_sf.shields_max-0.045)^(1.75))/(1+0.72*(wc_sf.shields_max-0.045)^(1.75));
                         
                         %reference height by Smith and McLean 1977 (Soulsby, pg.140)
                         zr(1,i) = 26.3*p.tauCr(wc_tot.current.fluid)*T(1,i)/(wc_sf.current.fluid.rho*g*(2.58-1))+(p.dn/12);
                         %reference concentration by Smith and McLean 1977 (Soulsby, pg.140)
                         cr(1,i) = (cb*rc*T(1,i))/(1+0.0024*T(1,i));
                         % sediment concentration at the top of zw:                         
                         cw(1,i) = cr(1,i)*(wc_sf.zw/zr(1,i)).^-b_max(1,i);
                         
                         %sediment concentration profile
                         cz(wc_sf.zw>z,i) = cr(1,i)*(z(wc_sf.zw>z)/zr(1,i)).^-b_max(1,i); %Equation 115a pg 148 Soulsby 1997 for zr<z<zw;
                         cz(z>wc_sf.zw,i) = cw(1,i)*(z(z>wc_sf.zw)/wc_sf.zw).^-b(1,i); %Equation 115b pg 148 Soulsby 1997 for zw<z<h;                            
                                                      
                     else
                         cz(:,i) = 0;                         
                     end
                     
                     IntSedSuspP50(1,i) = trapz(z(z>z0s),cz(z>z0s,i).*uz(z>z0s,i));
                     qss(1,i) = IntSedSuspP50(1,i);
                     
                       DirBat=240;
                       
                       qss_CS(1,i) = sind(mod(DirBat-wc_sf.current.ADCProfile.dirCur(1),360))*qss(1,i);
                       qss_AS(1,i) = -cosd(mod(DirBat-wc_sf.current.ADCProfile.dirCur(1),360))*qss(1,i);
                          %qss(i,j) = sqrt((qb_CSx(i,j))^2+ (qb_ASx(i,j))^2);
                     
                 end
             end
             
         end
         
         %%
         %%%Bed load%%         
          %%Soulsby & Van Rijn
                  
         function qt_svAsb = qt_svAsb(obj) %formula for bed load transport rate
             [n , m] = size(obj);
             qt_svAsb = zeros(n, m);
             for i = 1:n
                 for j = 1:m
                     um = obj(i,j).wavecurrent.current.uAvg;
                     Cd = obj(i,j).wavecurrent.current.dragCoef;
                     Uwc = sqrt(um.^2+0.018/Cd*obj(i,j).wavecurrent.wave.Uw^2);
                     T =Uwc - obj(i,j).Ucr;
                     if T > 0
                         qt_svAsb(i,j) = obj(i,j).Asb*um*(T)^2.4*(1-1.6*tand(obj(i,j).beta));
                     else
                         qt_svAsb(i,j) = 0;
                     end
                     
                 end
                 
             end
         end
                               
          %CEM partIII ch06 (Madsen 2002 based in Madsen 1993)
          
          function qb_CEM = qb_CEM(obj)
              [n , m] = size(obj);
              qb_CEM = zeros(n, m);
              
              for i = 1:n
                  for j = 1:m
                      s = obj(i,j).particle.rho/obj(i,j).wavecurrent.current.fluid.rho;
                      
                      Shieldscr = obj(i,j).particle.shieldsCr;
                      Shields_max(i,j) = obj(i,j).wavecurrent.shields_max;
                      Shields_mean(i,j) = obj(i,j).wavecurrent.shields_mean;
                      
                      tauCr = obj(i,j).particle.tauCr;
                      tau_max(i,j)=obj(i,j).wavecurrent.tau_max;
                      
                     
                      if Shields_max(i,j) > Shieldscr
                          
                          phi_u(i,j) = 6*((obj(i,j).wavecurrent.current.shields)^(3/2))*(obj(i,j).wavecurrent.wave.uStar/obj(i,j).wavecurrent.current.uStarFc)*(1.5*cosd(obj(i,j).wavecurrent.phi));
                          phi_v(i,j) = 6*((obj(i,j).wavecurrent.current.shields)^(3/2))*(obj(i,j).wavecurrent.wave.uStar/obj(i,j).wavecurrent.current.uStarFc)*(sind(obj(i,j).wavecurrent.phi));
                          
                          qb_CEM_u(i,j) = phi_u(i,j)*sqrt(g*(s-1)*obj(i,j).particle.dn^3); %sediment transport in the direction of wave propagation
                          qb_CEM_v(i,j) = phi_v(i,j)*sqrt(g*(s-1)*obj(i,j).particle.dn^3); %sediment transport at 90 deg counterclockwise from the wave direction
                          qb_CEM(i,j) = max(qb_CEM_u(i,j),qb_CEM_v(i,j));
                      else
                          qb_CEM(i,j) = 0;
                      end
                  end
              end
          end
           
                     
          %Soulsby - Meyer Peter & Muller - with waves: Soulsby(1997),pg.
          %167.
          function [qb_CS, qb_AS, qb_SMP] = qb_SMP(obj) %[qb_CSx, qb_ASx, qb_CSy, qb_ASy, qb_SMPx, qb_SMPy,qb_CS, qb_AS, qb_SMP] = qb_SMP(obj)
              [n , m] = size(obj);
              %qb_SMP = zeros(n, m);
              for i = 1:n
                  for j = 1:m 
                      wc_sf = obj(i,j).wavecurrent_sf;
                      wc_tot = obj(i,j).wavecurrent_tot;
                     
                      Shieldscr = obj(i,j).particle.shieldsCr(wc_sf.current.fluid);
                      s = obj(i,j).particle.rho/wc_sf.current.fluid.rho;
                      
                      Shields_max(i,j) = wc_sf.shields_max;
                      Shields_mean(i,j) = wc_sf.shields_mean;
                      
                      
                      if Shields_max(i,j) >= Shieldscr
                          
                          o=obj(i,j);
                          phi_x1(i,j) = 12*(Shields_mean(i,j)^(1/2))*(Shields_mean(i,j)-Shieldscr);
                          phi_x2(i,j) = 12*(0.95+0.19*cosd(2*wc_sf.phi))*sqrt(obj(i,j).bedshearstress.shieldsW)*Shields_mean(i,j);
                          phi_x(i,j) = max(phi_x1(i,j),phi_x2(i,j));
                          qb_SMPx(i,j) = phi_x(i,j)*sqrt(g*(s-1)*obj(i,j).particle.dn^3); %bed load transport rate in the direction of the current m2/s
                          phi_y( i,j) = 12*(0.19*Shields_mean(i,j)*o.bedshearstress.shieldsW^2*sind(2*wc_sf.phi))/(o.bedshearstress.shieldsW^(3/2)+1.5*Shields_mean(i,j)^(3/2));
                          qb_SMPy(i,j) = phi_y(i,j)*sqrt(g*(s-1)*obj(i,j).particle.dn^3); %bed load transport rate at right angles to current m2/s
                          
                          DirBat=240;
                          qb_CSx(i,j) = sind(mod(DirBat-wc_sf.current.ADCProfile.dirCur(1),360))*qb_SMPx(i,j);
                          qb_ASx(i,j) = -cosd(mod(DirBat-wc_sf.current.ADCProfile.dirCur(1),360))*qb_SMPx(i,j);
                          qbx(i,j) = sqrt((qb_CSx(i,j))^2+ (qb_ASx(i,j))^2);
                          
                          qb_CSy(i,j) = -cosd(mod(DirBat-wc_sf.current.ADCProfile.dirCur(1),360))*qb_SMPy(i,j);
                          qb_ASy(i,j) = sind(mod(DirBat-wc_sf.current.ADCProfile.dirCur(1),360))*qb_SMPy(i,j);
                          qby(i,j) = sqrt((qb_CSy(i,j))^2+ (qb_ASy(i,j))^2);
                          
                          qb_CS(i,j) = qb_CSx(i,j) + qb_CSy(i,j);
                          qb_AS(i,j) = qb_ASx(i,j) + qb_ASy(i,j);
                          
                          qb_SMP(i,j) = sqrt((qb_CS(i,j))^2+(qb_AS(i,j))^2);
                          
                      elseif Shields_max(i,j) < Shieldscr
                          
                          qb_SMPx(i,j) = 0;
                          qb_SMPy(i,j) = 0;
                          qb_SMP(i,j) = 0;
                          qb_CSx(i,j) = 0;
                          qb_ASx(i,j) = 0;
                          qb_CSy(i,j) = 0;
                          qb_ASy(i,j) = 0;
                          qb_CS(i,j) = 0;
                          qb_AS(i,j) = 0;
                          
                      else
                          
                          qb_SMPx(i,j) = NaN;
                          qb_SMPy(i,j) = NaN;
                          qb_SMP(i,j) = NaN;
                          qb_CSx(i,j) = NaN;
                          qb_ASx(i,j) = NaN;
                          qb_CSy(i,j) = NaN;
                          qb_ASy(i,j) = NaN;
                          qb_CS(i,j) = NaN;
                          qb_AS(i,j) = NaN;
                          
                      end
                      
                  end
                  
              end
              
          end
          
                                     
         %%%Total load%%
         %Soulsby-VanRijn
         function qt_sv = qt_sv(obj) %formula for total load transport - Soulsby-VanRijn
            [n , m] = size(obj);
            qt_sv = zeros(n, m);
            for i = 1:n
                for j = 1:m
                    um = obj(i,j).wavecurrent.current.uAvg;
                    [~, Cd] = obj(i,j).wavecurrent.current.tau;
                    Uwc = sqrt(um.^2+0.018/Cd*obj(i,j).wavecurrent.wave.Uw^2);
                    T =Uwc - obj(i,j).Ucr;
                    if T > 0
                        qt_sv(i,j) = obj(i,j).As*um*(T)^2.4*(1-1.6*tand(obj(i,j).beta));
                    else
                        qt_sv(i,j) = 0;
                    end
                    
                end
                
            end
         end        
                 
    end %methods
end

