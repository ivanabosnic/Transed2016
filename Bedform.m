classdef Bedform
    %Calculates the bedform geometry (height and wave lenght)
    %
    
    properties
        particle;
        wave;
        current;
       
       
        
    end
    
    
    methods
        
        function obj = Bedform(dim,particle,wave,current)
            if nargin == 1
                if isscalar(dim)
                    obj(dim,1) = Bedform;
                else
                    n = dim(1);
                    m = dim(2);
                    obj(n,m) = Bedform;
                end
            elseif nargin >= 2
                [n, m] = size(wave);
                obj(n, m) = Bedform;
                for i = 1:n
                    for j = 1:m
                        obj(i,j).wave = wave(i,j);
                        obj(i,j).current = current(i,j);
                        obj(i,j).particle = particle;
                       
                    end
                end
            end
            
        end
        
    
        
        function lambda = lambda(obj) %Formulation of Soulsby & Whitehouse(2005)
            [n , m] = size(obj);
            lambda = zeros(n, m); %bedform wave lenght
            for i = 1:n
                for j = 1:m
                    if obj(i,j).wave.shields>obj(i,j).particle.shieldsCr
                        lambda(i,j) = obj(i,j).wave.A/(1+1.87*1e-3*obj(i,j).wave.A/obj(i,j).particle.di*(1-exp(-(2e-4*(obj(i,j).wave.A/obj(i,j).particle.di))^(1.5))));
                    else
                        if i==1;
                        lambda(i,j) = lambda(i,j);
                        else
                            lambda(i,j) = lambda(i-1,j);
                        end
                    end
                end
            end
        end
        
        function eta = eta(obj)
            [n , m] = size(obj);
            eta = zeros(n, m); %bedform height
            for i = 1:n
                for j = 1:m
                    if obj(i,j).wave.shields>obj(i,j).particle.shieldsCr
                        eta(i,j) = obj(i,j).lambda*0.15*(1-exp(-(5000/(obj(i,j).wave.A/(obj(i,j).particle.di)))^3.5));
                    else
                        if i==1;
                        eta(i,j)=eta(i,j);
                        else
                        eta(i,j)=eta(i-1,j);       
                        end
                    end
                end
            end
        end
        
        
        
        
    end
end

