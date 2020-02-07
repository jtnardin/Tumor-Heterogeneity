function driver(is_discrete)
% Driver for the non-parametric estimation fram
%
%    INPUTS:
%        is_discrete: indicator whether you want to estimate using discrete
%            approximations (1) or continuous/spline approximations (0)
%
%    OUTPUTS:
%        precomputed_sol: a mat-file with the full forward solutions over  
%             the parameter space of interest (named 'fullsol')
%
% Written by Erica Rutter (July 2019)

%% declare variables and data
real_data ='DataGeneration/bimodal_rho_const_VP.mat';

Nr = 500; % quadrature nodes (if using splines)
nrs=10:10:500;    % number of nodes, M
param.D=1e-6;   
rmin=0; % min for r parameter
rmax=0.06; % max for r parameter

%% Precompute solns for all nr values (if not already exist)
load(real_data,'X_m','Y_m','T_m')
param.xs = squeeze(X_m(1,:,1))';
param.ys = Y_m(:,1,1);
ts=zeros(size(param.ys));
param.t = squeeze(T_m(1,1,:));
ts(1:8)=param.t;
X=[ts,param.xs,param.ys];

% precompute solutions if they don't already exist
if is_discrete==1
    for nr_val=1:length(nrs)
        Nr=nrs(nr_val);
        precomputed_sol = strcat('PrecomputedSolutions/fullsol_2par_longtime_with_nr_',int2str(Nr),'.mat');
        if ~exist(precomputed_sol, 'file')
            disp('Precomputing Solutions')
            param.bsr=linspace(rmin,rmax,Nr);
            forwardsolvers(param,precomputed_sol,X)
        end
    end
else
    precomputed_sol = strcat('PrecomputedSolutions/fullsol_2par_longtime_with_nr_',int2str(Nr),'.mat');
    if ~exist(precomputed_sol,'file')
        disp('Precomputing solutions')
        param.bsr=linspace(rmin,rmax,Nr);
        forwardsolvers(param,precomputed_sol,X)
    end
end

%% Use Akaike Information Criteria to select optimal node number
load(real_data,'Y')
num_pats=size(Y);
num_pats = num_pats(1);
% load and reshape the full solutions into the 'aggregate' total tumor
agg_sol = reshape(Y,num_pats,length(param.xs),length(param.ys),length(param.t));
agg_sol = squeeze(sum(agg_sol))/num_pats;
param.rmin=rmin;
param.rmax=rmax;
if is_discrete==0
    param.Nr=Nr;
end

% Find the optimal number of nodes for discrete and continuous
% approximations
BestM_value=strcat('Outputs/BestM_discrete_',int2str(is_discrete),'.mat');

% do the optimization routine if it hasn't already been done
if ~exist(BestM_value,'file')
    disp('Finding optimal M number of nodes')
    best_nr=akaike_finder(param,agg_sol,nrs,is_discrete,BestM_value);
end
load(BestM_value,'best_nr','optweights','best_nr_num');
optweight=optweights{best_nr_num};

%% Plot the solutions if wanted
param.nodesr=linspace(rmin,rmax,best_nr);
param.bsr=linspace(rmin,rmax,Nr);
load(precomputed_sol,'fullsol')

disp('plotting!')
if is_discrete==1
    plotter_disc(param,real_data,optweight);
else
    plotter_spline(param,fullsol,agg_sol,real_data,optweight);
end

