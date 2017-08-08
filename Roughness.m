classdef Roughness
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        particle;
        fluid;
        z0; %total
        z0s; %skin friction/"grain-related"
        z0s2;
        z0f; %form drag
        z0t; %sediment transport (at very high flow speeds)
        alpha_r; %
        bedform;
        
        
    end
    
    methods
        function obj = Roughness(particle, bedform, alpha_r) %,bedshearstress)
            
            %delta_r = ripple height
            %lambda_r = ripple wavelength
            
            if nargin >= 1
                n =1; %3273; %3153; %657;; %2433; %8537
                m = 1;
                obj(n,m)=Roughness;
                for i = 1:n
                    for j = 1:m
                        obj(i,j).particle = particle;
                        obj(i,j).z0s2 = 5.25e-5;% roughness of d50=6.3e-4 mm (in situ sediment SHORE)%0.3e-3;
                        obj(i,j).z0s = obj(i,j).particle.dn/12; %Equation 25 in Soulsby (1997)
                        obj(i,j).z0 = obj(i,j).z0s;
                    end
                end
                
                if nargin >= 2
                    
                    [n , m] = size(bedform);
                    obj(n, m) = Roughness;
                    
                    for i = 1:n
                        for j = 1:m
                            
                            obj(i,j).bedform = bedform(i,j);
                            eta = obj(i,j).bedform.eta;
                            lambda = obj(i,j).bedform.lambda;
                            obj(i,j).particle = obj(i,j).bedform.particle;
                            obj(i,j).z0s = obj(i,j).particle.dn/12; %Equation 25 in Soulsby (1997)
                            obj(i,j).z0f = alpha_r*eta^2/lambda; %Equation 90 in Soulsby(1997)
                            
                            if isnan(obj(i,j).z0f)
                                obj(i,j).z0 = obj(i,j).z0s;
                            else
                                obj(i,j).z0 = obj(i,j).z0s+obj(i,j).z0f;
                            end
                            
                        end
                        
                        %             elseif nargin >= 4
                        %                  obj.z0t(i,j) = (5*obj(i,j).bedshearstress.tau0s)/30*g*(obj.particle.rho-obj.fluid.rho) ; %Wilson (1989a)|Equation 42 in Soulsby(1997)
                        %                  obj.z0(i,j) = obj.z0s + obj.z0f + obj.z0t(i,j);
                        
                    end
                    
                    
                end
            end
            
        end
    end
end

