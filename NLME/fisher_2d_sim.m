%fisher_2d_sim written 12-6-2019 to run the 2d fisher-KPP Equation for
%"A tutorial Review of Mathematical Techniques for Quantifying 
%Tumor Heterogeneity" by Everett et al. 

function yd = fisher_2d_sim(phi,X)

    
    if any(phi<0)
        
        %return nonsense, don't simulate with negative params
        yd = 10*ones(size(X,1),1);
        
    else
       
        
        %retrieve time and space domains from X
        t = unique(X(:,1));
        x = unique(X(:,2));
        y = unique(X(:,3));
        [X_d,Y_d] = meshgrid(x,y);

        
        
        %finer simulation domains
        N = 256;
        x_sim = linspace(x(1),x(end),N);
        y_sim = linspace(y(1),y(end),N);
        t_sim = t;
        [X_sim,Y_sim] = meshgrid(x_sim,y_sim);
        
        % %  Set up initial condition
        s  = 0.015;    % Std. dev of initial condition
        I0 = 1;        % Integral of intial cell density ("number of initial cells")
        initcond = @(x,y) I0*exp(-((x-0.5).^2+(y-0.5).^2)/(2*s^2))/(2*pi*s^2);
        u0 = initcond(X_sim(:,:,1),Y_sim(:,:,1));
            
        dx = x_sim(2) - x_sim(1);
        
        
        
        if length(phi) == 1
            %Estimating rho
            odefun = @(T,Y)myODEFun2D(Y,1e-6,phi,10.001,dx,length(x_sim));
            [T,Y] = ode45(odefun,t,u0(:));
        elseif length(phi) == 2
            %estimating (D,rho)
            odefun = @(T,Y)myODEFun2D(Y,phi(1),phi(2),10.001,dx,length(x_sim));
            [T,Y] = ode45(odefun,t,u0(:));
        end
        
        %interpolate computational grid back to same grid as data
        yd = zeros(length(x),length(y),length(t));
        for i = 1:length(t)
            yd(:,:,i) = interp2(X_sim,Y_sim,reshape(Y(i,:),N,N),X_d,Y_d,'cubic'); 
        end
        
        yd = yd(:);

    end

end
