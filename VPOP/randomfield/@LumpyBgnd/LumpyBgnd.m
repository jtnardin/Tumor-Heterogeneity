% -------------------------------------------------------------------------
% ******************* Tumor Heterogeneity Package *************************
% File:     LumpyBgnd.m
% Author:   Nick Henscheid
% Date:     9-2016, 5-2020 
% Info:     u = LumpyBgnd(varargin). The lumpy background object.
%           The name-value pairs in varargin determine the random field
%
%           Note: this is a simplified version of LumpyBgnd, using only
%           cpu-based Matlab code and no compiled mex files. 
%           See https://github.com/nhenscheid/RDEVPop/ for a more
%           sophisticated version that uses GPU-based mex files for much
%           faster computation. 
%      
%           A realization of a LumpyBgnd random field takes the form 
%                      f(x) = B0 + \sum b(j)*l(x-c(j);th(j)) 
%           where B0, b, c and th are the field parameters, some or all of
%           which may be randomized. See below. 
%
% Inputs:   The following name-value pair inputs are allowed.  All have 
%           default values.  
%
%           'b' (a scalar or vector)   The lump amplitude, defined so that  
%                              l_j(x) = b(j)*l(x;th_j)
%           'B0' (a scalar)   The DC offset 
%           'cov' The Gaussian lump covariance specification.  If cov is a
%           scalar, this specifies the isotropic lump variance.  If cov is
%           a dim-by-1 vector, it specifies a diagonal covariance matrix
%           (both elements must be positive).  If cov is a dim-by-dim
%           matrix, it specifies an anisotropic covariance (must be
%           positive-definite).
%           'Kbar' (a positive scalar, the mean number of lumps)
%           'N' (an integer or dim-by-1 vector of integers) the default
%           field evaluation grid size e.g. N=[256,256].
%               
% Contact: nph@email.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

classdef LumpyBgnd < handle
    properties (SetObservable = true)
        b       % Lump "Amplitude" vector (usually a constant
        B0      % DC offset i.e. 
        cov     % Lump covariance matrix 
        Kbar    % Mean number of lumps
        centers % Lump centers
        N       % Default number of eval. pts in ea. direction e.g. 256
        israndnumlumps = 1; % Randomize lumps?
    end
    
    properties (Dependent)
        L;% Bounding box for evaluation (lump centers can extend slightly beyond to avoid edge effect issues)
        dim; 
        K; 
    end

    properties (SetAccess=private)
        Support    % Support set for evaluation.  Can only set when LumpyBgnd is first initialized!
        padfactor = 3;      % Padding factor (so that boundary effects don't occur)
        showwarnings = 1;
    end
    
    % Standard methods
    methods
        function obj = LumpyBgnd(varargin)
            p = obj.ParseInputs(varargin{:});
            obj.Kbar    = p.Results.Kbar;
            obj.b       = p.Results.b;
            obj.B0      = p.Results.B0;
            obj.Support = p.Results.Support;
            if(numel(p.Results.cov)==1)
                obj.cov = p.Results.cov*eye(obj.Support.dim);
            else
                obj.cov = p.Results.cov;
            end
            obj.centers = p.Results.centers;
            obj.centers;
            
            obj.N       = p.Results.N;
            if(numel(obj.centers)==0)
                % Need to generate random lump centers
                obj.Randomize();
            end
            addlistener(obj,'Kbar','PostSet',@LumpyBgnd.HandlePropertyEvents);
            addlistener(obj,'cov','PostSet',@LumpyBgnd.HandlePropertyEvents);
        end
        
        % Get methods for dependent properties
        function val = get.L(obj)
            val = obj.Support.L;
        end
        function val = get.dim(obj)
            val = obj.Support.dim;
        end
        function val = get.K(obj)
            val = size(obj.centers,1);
        end
        function set.K(obj,value)
            if(isnumeric(value)&&value>=0&&(mod(value,1)==0))
                obj.K = value;
            else
                error('Invalid K value');
            end
        end
        
        % Misc utilities
        function TurnOffWarnings(obj)
            obj.showwarnings = 0;
        end
        
        function TurnOnWarnings(obj)
            obj.showwarnings = 1;
        end
        
        function SetPadFactor(obj,x)
            obj.padfactor = x;
        end
        
        % Externally defined functions
        p = ParseInputs(varargin);
        u = Eval(obj,X,XSize);
        obj = Randomize(obj)
        varargout = plot(obj,alpha);
        z = minus(x,y);  % Method to subtract two lumpy backgrounds
        z = plus(x,y);   % Method to add two lumpy backgrounds
        
        function U = Sample(obj,Ns)
            % Generates Ns samples of the lumpy background, returning it in
            % a cell array U (U{i} is a sample for each i)
            U = cell(Ns,1);
            for i=1:Ns
                obj.Randomize;
                U{i} = obj.Eval;
            end
        end
        
        function v = Copy(obj)
           % Makes a copy of the current object 
            v = LumpyBgnd('Support',obj.Support); v.TurnOffWarnings;
            v.centers = obj.centers;
            v.b       = obj.b;
            v.cov     = obj.cov;
        end  
    end
    
    % Static methods
    methods (Static)
        HandlePropertyEvents(src,evnt);   % Defined externally
    end
    
end