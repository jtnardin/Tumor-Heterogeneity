 classdef RectSupport 
    properties 
       L = [0,1;0,1]  
    end
    
    properties (Dependent)
        origin;
        dim; 
    end
    
    methods 
        function obj = RectSupport(Li)
            obj.L = Li;
        end
        
        function val = get.dim(obj)
            val = size(obj.L,1);
        end
        function val = get.origin(obj)
            val = zeros(1,size(obj.L,1));
            for i=1:size(obj.L,1)
                val(i) = obj.L(i,1);
            end
        end
              
        function plot(obj,alpha)
            % alpha is the opacity 
            if(nargin<2)
                alpha = 0.5;
            end
            l = obj.L;
            if(obj.dim==2)
               x = [l(1,1),l(1,2),l(1,2),l(1,1),l(1,1)];
               y = [l(2,1),l(2,1),l(2,2),l(2,2),l(2,1)];
               line(x,y,'LineWidth',2);
               xlabel('$x$ (m)');ylabel('$z$ (m)');
            elseif(obj.dim==3)
               lx = l(1,2)-l(1,1);
               ly = l(2,2)-l(2,1);
               lz = l(3,2)-l(3,1);
               plotcube([lx,ly,lz],obj.origin,alpha,[0.1,0.1,0.2]);
               xlabel('$x$ (m)');ylabel('$y$ (m)');zlabel('$z$ (m)');  
            end
        end
        
        function r = bdyFun(obj,theta,z)
            % aziumuth angle theta, elevation angle phi
            r = obj.origin + [obj.radius*cos(theta),obj.radius*sin(theta),z];
        end
        
        function a = isBdy(obj,r)
            a = (norm(r-obj.origin)==obj.radius);
        end
            
        function a = isBdyPlus(obj,r,s)
            a = dot(r,s)>0;
        end
        
        function a = isBdyMinus(obj,r,s)
            a = dot(r,s)<0;
        end
        
        
        function X = unifGrid(obj,N)
            l = obj.L;
            xtmp = linspace(l(1,1),l(1,2),N(1));
            ytmp = linspace(l(2,1),l(2,2),N(2));
            if(obj.dim ==2)
                [xx,yy] = meshgrid(xtmp,ytmp);
                X = cat(3,xx,yy);
            elseif(obj.dim == 3)
                ztmp = linspace(l(3,1),l(3,2),N(3));
                [xx,yy,zz] = meshgrid(xtmp,ytmp,ztmp);
                X = cat(4,xx,yy,zz);
            end
        end
        
        function X = unifSample(obj,N)
            X = rand(N,obj.dim);
            l = obj.L; 
            dx = l(1,2)-l(1,1);
            dy = l(2,2)-l(2,1);
            X(:,1) = obj.L(1,1) + dx*X(:,1);
            X(:,2) = obj.L(2,1) + dy*X(:,2);
            if(obj.dim==3)
                dz = l(3,2)-l(3,1);
                X(:,3) = obj.L(3,1) + dz*X(:,3);
            end
        end    
    end
end
