% -------------------------------------------------------------------------
% ****************** Tumor Heterogeneity Package **************************
% File:     RDE.m
% Author:   Nick Henscheid
% Paper:    Everett et al. 'A tutorial review of mathematical techniques 
%           for quantifying tumor heterogeneity'. Math. Biosci. Eng, 2020
%           doi: 10.3934/mbe.2020207
% Date:     10-2019, 2-2020, 5-2020
% Info:     This is the class definition file for the Matlab RDE solver           
% Inputs:         
%               
% Contact:  nph@email.arizona.edu, jtnardin@ncsu.edu 
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------
classdef RDE < handle
    properties 
        rho;   % Growth factor     
        kappa; % Carrying capacity 
        u0;    % Initial condition 
        D;     % Diffusion coefficient 
        %      NOTE: rho, kappa, 
        grid;  % Grid {x,y} or {x,y,z} (assume rectalinear grid)
    end
    
    properties (Dependent)
        dim;       % Ambient dimension (2 or 3 - inherited from grid)
        n_grid;    % Number of grid points in each direction e.g. [128,128]
        h_grid;    % Mesh width in each direction e.g. [0.01,0.01].  
    end
   
    methods
        function obj = RDE
            % Nothing is set by default - user must define all class
            % properties manually. See example file RDE_Virtual_Pop_2D. 
        end
        
        % Externally defined functions
        n = Solve(obj,T);   % Solve the RDE for the time vector T
        [n,rho,kappa,D] = Sample(obj,N,T); % Create a virtual population 
                                           % of N subjects via
                                           % randomization
        function Randomize(obj)
            % Randomize the coefficients.   Note: this assumes all of these
            % have a method called "randomize".  
            obj.rho.Randomize; 
            obj.kappa.Randomize; 
            obj.D.Randomize;  
        end
        
        % Get methods for dependent properties
        function y = get.dim(obj)
            y = length(obj.grid);  % Ambient dimension set by number of grid vectors
        end
        function y = get.n_grid(obj)
            if(obj.dim==2)
                y = [length(obj.grid{1}),length(obj.grid{2})];
            elseif(obj.dim==3)
                y = [length(obj.grid{1}),length(obj.grid{2}),length(obj.grid{3})];
            end
        end
        function y = get.h_grid(obj)
            % NOTE: this assumes that the grid is uniform. 
            if(obj.dim==2)
                y = [obj.grid{1}(2)-obj.grid{1}(1),...
                     obj.grid{2}(2)-obj.grid{2}(1)];
            elseif(obj.dim==3)
                y = [obj.grid{1}(2)-obj.grid{1}(1),...
                     obj.grid{2}(2)-obj.grid{2}(1),...
                     obj.grid{3}(2)-obj.grid{3}(1)];
            end
        end
    end
    
    methods (Static) 
       Y_o = ODEFun2D(Y_i,D,Aa,Kk,hx,hy,nx,ny); 
       Y_0 = ODEFun3D(Y_i,D,Aa,Kk,hx,hy,hz,nx,ny,nz); 
    end
end