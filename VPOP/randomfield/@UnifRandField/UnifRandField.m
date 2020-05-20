% -------------------------------------------------------------------------
% ************************* RDEVPop Package *******************************
% File:     UnifRandField.m
% Author:   Nick Henscheid
% Date:     10-2019, 2-2020
% Info:     u = UnifRandField(varargin). The uniform random field object.
%           A uniform random field takes a constant value c, sampled from 
%           a distribution P_c.  Note that "uniform" applies to the spatial
%           variable!  The probability distribution P_c need not be
%           uniform!
%           The name-value pairs in varargin determine the random field
%           
% Inputs: 
%           'dist' (function handle) random number generator for the field
%                  value.  Must be of the form @()rng(params) i.e. a
%                  "no-input function handle" with all parameters specified 
%           'c' (scalar)   Realized field value
%           'N' (an integer or dim-by-1 vector of integers)
%           'S' a support object
%               
% Contact: nph@email.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

classdef UnifRandField < handle
    properties (SetObservable = true)
        c    % Realized field value f(x)\equiv c (can be set manually)
        N    % Default number of eval. pts in ea. direction
        dist % Distribution (function handle of the form @()gen_rand_num(params))
    end
    
    properties (Dependent)
        L;% Bounding box for evaluation 
        dim; 
    end

    properties (SetAccess=private)
        Support    % Support set for evaluation.  Can only set when first initialized!
    end
    
    % Standard methods
    methods
        function obj = UnifRandField(varargin)
            p            = obj.ParseInputs(varargin{:});
            obj.Support  = p.Results.Support;
            obj.N        = p.Results.N;
            obj.dist     = p.Results.dist;
            obj.Randomize();
        end
        
        function val = get.L(obj)
            val = obj.Support.L;
        end
        
        function val = get.dim(obj)
            val = obj.Support.dim;
        end
        
        % Externally defined functions
        p = ParseInputs(varargin);
        u = Eval(obj,X,XSize);
        obj = Randomize(obj,shuff);
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
           % Make a copy of the field
           v = UnifRandField('dist',obj.dist,'N',obj.N,'Support',obj.Support); 
           v.c = obj.c; % need to set manually so the value is the same
        end
    end
end