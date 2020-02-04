function [best_nr]=akaike_finder(param,agg_sol,nrs,is_discrete,BestM_value)

% Finds the optimal nr and nD value for a specific set of parameters and
% saves them.

[datpoints,obs]=size(agg_sol);

%% Optimize the solution for all meshes of nD and nr and calculate AICs
AICs=10^20*ones(1,length(nrs));
optweights=cell(1,length(nrs));

if is_discrete==0
    precomputed_sol = strcat('PrecomputedSolutions/fullsol_2par_longtime_with_nr_',int2str(param.Nr),'.mat');
    load(precomputed_sol,'fullsol');
    param.bsr = linspace(param.rmin,param.rmax,param.Nr);
end
for nr_val=1:length(nrs)
    nr=nrs(nr_val);
    disp(strcat('nr=',int2str(nr)))
    param.nodesr=linspace(param.rmin,param.rmax,nr);
    if is_discrete==1
        precomputed_sol = strcat('PrecomputedSolutions/fullsol_2par_longtime_with_nr_',int2str(nr),'.mat');
        load(precomputed_sol,'fullsol');
        titles='a.mat';
        [optweight,curr_error,finish_flag] = prohorovit_r_GLS(titles,param,fullsol,agg_sol,0);
    else
        titles='a.mat';
        [optweight,curr_error,finish_flag] = prohorovit_spline_r_indy_GLS(titles,param,fullsol,agg_sol,0);
    end
    if finish_flag>0
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
