function [best_nr]=akaike_finder(param,agg_sol,nrs,is_discrete,BestM_value)
% AKAIKE_FINDER finds the optimal number of nodes for the parameter mesh of
% r by calculating the akaike information criteria (AIC) score for the
% various numbers of nodes
% 
%    INPUTS:
%        param: a variable that contains the following information:
%            -param.xs: the x vector for solutions
%            -param.ys: the y vector for solutions
%            -param.t:  the t vector for solutions
%            -param.nodesr: the vector of nodes for parameter 'r' for the
%                optimization 
%            -param.bsr: vector of nodes for parameter 'r' on which to solve
%                the pde (used for the spline approximations only)
%        agg_sol: the aggregate 'data' we are comparing with
%        nrs: a vector of the potential number of nodes to use in the
%            approximation
%        is_discrete: determining whether to do the discrete approximations
%            (1) or spline approximations (0)
%        BestM_value: the savename for the mat file containing the AIC
%            scores, optimal number of nodes, etc
%
%    OUTPUTS:
%        BestM_value: a mat-file with the full AIC scores, the optimal
%            number of nodes to choose, the optimal weights for those nodes
%
% Written by Erica Rutter (July 2019)

%% Optimize the solution for all meshes of nD and nr and calculate AICs
[datpoints,obs]=size(agg_sol);
AICs=10^20*ones(1,length(nrs));
optweights=cell(1,length(nrs));

if is_discrete==0
    % in spline case precomputed solution are always the same
    precomputed_sol = strcat('PrecomputedSolutions/fullsol_2par_longtime_with_nr_',int2str(param.Nr),'.mat');
    load(precomputed_sol,'fullsol');
    param.bsr = linspace(param.rmin,param.rmax,param.Nr);
end
for nr_val=1:length(nrs)
    nr=nrs(nr_val);
    disp(strcat('nr=',int2str(nr)))
    param.nodesr=linspace(param.rmin,param.rmax,nr);
    if is_discrete==1
        % in discrete case, we must load in precomputed solutions every
        % cycle and we call the discrete optimizer, prohorovit_r_GLS
        precomputed_sol = strcat('PrecomputedSolutions/fullsol_2par_longtime_with_nr_',int2str(nr),'.mat');
        load(precomputed_sol,'fullsol');
        titles='';
        [optweight,curr_error,finish_flag] = prohorovit_r_GLS(titles,param,fullsol,agg_sol,0);
    else
        % in the spline case, we call the spline optimizer,
        % prohorovit_spline_r_indy_GLS
        titles='';
        [optweight,curr_error,finish_flag] = prohorovit_spline_r_indy_GLS(titles,param,fullsol,agg_sol,0);
    end
    if finish_flag>0
        % if the opitimization procedure finished correctly, calculate AIC
        % scores in the context of least squares
        AICs(nr_val)=datpoints*obs*log(curr_error/(datpoints*obs))+datpoints*obs*(1+log(2*pi))+2*(length(nrs)+1);
        optweights{nr_val}=optweight;
    else
        disp(strcat('Did not converge for nr=',int2str(nr)))
    end
end

%% calculate the best Akaike Information Criteria Scores

mins=min(min(AICs));
[best_nr_num]=find(AICs==mins);

best_nr=nrs(best_nr_num);
save(BestM_value,'best_nr','best_nr_num','AICs','optweights')
