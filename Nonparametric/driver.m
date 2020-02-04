function driver(is_discrete)
% Driver for Optimization Stuff %



%% declare variables and data
real_data ='DataGeneration/bimodal_rho_const_VP.mat';

Nr = 500; % quadrature nodes (if using spline)
nrs=10:10:500;    % number of nodes, M
param.D=1e-6;   
rmin=0; % min for parameter 2
rmax=0.06; % max for parameter 2

%% Precompute solns for all nr values (if not already exist)
load(real_data,'X_m','Y_m','T_m')
param.xs = squeeze(X_m(1,:,1))';
param.ys = Y_m(:,1,1);
ts=zeros(size(param.ys));
param.t = squeeze(T_m(1,1,:));
ts(1:8)=param.t;
X=[ts,param.xs,param.ys];

% calculate solutions and reshape if they don't already exist
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

if ~exist(BestM_value,'file')
    disp('Finding optimal M number of nodes')
    best_nr=akaike_finder(param,agg_sol,nrs,is_discrete,BestM_value);
else
    load(BestM_value,'best_nr','optweights','best_nr_num');
    optweight=optweights{best_nr_num};
end

%% Plot the solutions if wanted
param.nodesr=linspace(rmin,rmax,best_nr);
param.bsr=linspace(rmin,rmax,Nr);
load(precomputed_sol,'fullsol')

% %titler_spline=strcat('Outputs/D_',distributions{dist1},'_nD',int2str(best_nD),'_ND_',int2str(ND),'_r_',distributions{dist2},'_nr',int2str(best_nr),'_Nr_',int2str(Nr),'_longtime_error_',int2str(100*err_level),'_xscale_',int2str(x_scale),'_tscale_',int2str(t_scale),'_randinits.mat');
% %titler_disc=strcat('Outputs/D_',distributions{dist1},'_nD',int2str(best_nD),'_r_',distributions{dist2},'_nr',int2str(best_nr),'_longtime_error_',int2str(100*err_level),'_xscale_',int2str(x_scale),'_tscale_',int2str(t_scale),'_randinits.mat');
% %precomputed_sol = strcat('PrecomputedSolutions/fullsol_2par_longtime_with_nr_',int2str(best_nr),'_nd_',int2str(best_nD),'.mat');
% %load(precomputed_sol,'fullsol');
% 
disp('plotting!')
keyboard
if is_discrete==1
    plotter_disc(param,fullsol,agg_sol,real_data,optweight);
else
    plotter_spline(param,fullsol,agg_sol,real_data,optweight);
end
% 
% elseif is_discrete==1
%     if ~exist(titler_disc,'file')
%         for curr_init =1:num_inits
%             disp(strcat('Initialization #', int2str(curr_init), ' out of #', int2str(num_inits)))
%             titles=strcat('Outputs/D_',distributions{dist1},'_nD',int2str(best_nD),'_r_',distributions{dist2},'_nr',int2str(best_nr),'_longtime_error_',int2str(100*err_level),'_xscale_',int2str(x_scale),'_tscale_',int2str(t_scale),'_init',int2str(curr_init),'.mat');
%             if ~exist(titles,'file')
%                 [curr_optweight,~,finish_flag] = prohorovit_2par_GLS(titles,param,fullsol,persol,x_scale,t_scale,err_level,1);
%                 if finish_flag~=0
%                     optweights(curr_init,:)=curr_optweight;
%                 end
%             end
%         end
%         save(titler_disc,'optweights')
%     end
% end

%% Plot the solutions if wanted

%if do_plot==1
%    plot_fits(persol,fullsol,optweights,param,t_scale,x_scale)
%end
%% perform UQ 