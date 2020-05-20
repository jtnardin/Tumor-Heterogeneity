% -------------------------------------------------------------------------
% ****************** Tumor Heterogeneity Package **************************
% File:     fisher_2d_sim.m 
% Paper:    Everett et al. 'A tutorial review of mathematical techniques 
%           for quantifying tumor heterogeneity'. Math. Biosci. Eng, 2020
%           doi: 10.3934/mbe.2020207
% Date:     12-19, 5-2020
% Info:     Solves the 2D reaction-diffusion equation with scalar growth 
%           and diffusion values.          
% Inputs:   fisher_2d_sim(X,phi)
%           X = [T(:),X(:),Y(:)], assuming that T,X,Y were created using
%           meshgrid(t,x,y). 
%           phi: either phi = rho (the growth
%           constant), or phi = [D,rho] (the diffusion constant and the
%           growth constant).       
%               
% Contact:  nph@email.arizona.edu, jtnardin@ncsu.edu 
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------
function yd = fisher_2d_sim(phi,X)
    if any(phi<0)
        %return nonsense, don't simulate with negative params
        yd = 10*ones(size(X,1),1);
    else
        % Define the RDE solver object
        R = RDE; 
        
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
        
        % Set the grid in the RDE solver 
        R.grid = {x_sim,y_sim}; 
        
        % %  Set up initial condition
        s  = 0.015;    % Std. dev of initial condition
        I0 = 1;        % Integral of intial cell density ("number of initial cells")
        initcond = @(x,y) I0*exp(-((x-0.5).^2+(y-0.5).^2)/(2*s^2))/(2*pi*s^2);
        u0 = initcond(X_sim(:,:,1),Y_sim(:,:,1));
        
        % Set the initial condition in the RDE solver
        R.u0 = u0; 
        
        % Set the carrying capacity in the RDE solver 
        kappa_0 = 10.001; 
        R.kappa = kappa_0*ones(N); 
        
        if length(phi) == 1
            %Estimating rho
            R.D   = 1e-6*ones(N); 
            R.rho = phi*ones(N); 
            % Solve the RDE and extract the tumor cell density
            n     = R.Solve(t_sim).cell_density; 
        elseif length(phi) == 2
            %estimating (D,rho)
            R.D   = phi(1)*ones(N); 
            R.rho = phi(2)*ones(N); 
            n     = R.Solve(t_sim).cell_density; 
        end
        
        %interpolate computational grid back to same grid as data
        yd = zeros(length(x),length(y),length(t));
        for i = 1:length(t)
            yd(:,:,i) = interp2(X_sim,Y_sim,n(:,:,i),X_d,Y_d,'cubic'); 
        end
        
        yd = yd(:);
    end
end
