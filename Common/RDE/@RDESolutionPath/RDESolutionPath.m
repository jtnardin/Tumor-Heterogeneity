% -------------------------------------------------------------------------
% ************************* RDEVPop Package *******************************
% File:     RDESolutionPath.m
% Author:   Nick Henscheid
% Date:     10-2019, 2-2020, 5-2020
% Info:     This is the class definition file for the RDESolutionPath
%           object.  An instance of RDESolutionPath contains the 
%           computational grid and the tumor cell density on that grid
%           The object also implements several useful methods for plotting
%           and analyzing a tumor cell density path. 
%           
% Inputs:         
%               
% Contact: nph@email.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

classdef RDESolutionPath
    
    properties
        grid;          % Grid definition {x_g,y_g} or {x_g,y_g,z_g} 
        cell_density;  % Cell density array (size [nx,ny,nt] or [nx,ny,nz,nt])
        times;         % Vector of solution times
    end
    
    properties (Dependent)
        dim;         % Ambient dimension (inherited from the grid)
        num_times;   % Number of solution times
        grid_size;   % [nx,ny] or [nx,ny,nz]
    end
    
    methods
        function obj = RDESolutionPath
            % Nothing is set by default - user must define all class
            % properties manually.  Typically, an instance of
            % RDESolutionPath is produced by the RDE.solve() method, which 
            % sets the properties automatically. 
        end
        
        function y = get.dim(obj)
            y = length(obj.grid);
        end
        function y = get.num_times(obj)
            y = length(obj.times);
        end
        function y = get.grid_size(obj)
            y = zeros(obj.dim,1);
            for i=1:obj.dim
                y(i) = length(obj.grid{i});
            end
        end
        
        function plot(obj,T,varargin)
            % Plots the solution path 
            % T is an integer index set e.g. 1:2:10
            
            if(nargin<2)
                nt = obj.num_times;
            else 
                nt = length(T);
            end
            
            if(obj.dim == 2)
                figure; 
                num_cols = 5;
                num_rows = ceil(nt/num_cols);
                ax = tight_subplot(num_rows,num_cols,[.03 .03],[.01 .01],[.03 .03]);
                N = obj.TumorBurden; 
                for i=1:nt
                    axes(ax(i)); 
                    imagesc(obj.grid{1},obj.grid{2},obj.cell_density(:,:,T(i)));
                    set(gca,'YDir','normal');
                    axis image;
                    title(sprintf('$N(%1.1f) = %1.2e$',obj.times(T(i)),N(T(i))));
                end
            end
        end  % plot
        
        function vols = TumorBurden(obj)
            % Estimates the total tumor burden time series in number of 
            % cells.  This assumes that the units of cell_density are
            % number per unit volume. 
            
            dx = obj.grid{1}(2) - obj.grid{1}(1); % Assume uniform
            dy = obj.grid{2}(2) - obj.grid{2}(1); % Assume uniform
            if(obj.dim==3)
                dz = obj.grid{3}(2) - obj.grid{3}(1);
            end
            
            vols = zeros(obj.num_times,1);
            if(obj.dim==2)
                for i=1:obj.num_times
                    vols(i) = sum(sum(obj.cell_density(:,:,i)))*dx*dy;
                end
            end
            if(obj.dim==3)
                for i=1:obj.num_times
                    vols(i) = sum(sum(sum(obj.cell_density(:,:,:,i))))*dx*dy*dz;
                end
            end
        end
       
    end % Methods
    
end % RDESolutionPath