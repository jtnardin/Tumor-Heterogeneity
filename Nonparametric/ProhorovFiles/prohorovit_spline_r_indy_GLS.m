function [optweight,current_err,converge_flag]=prohorovit_spline_r_indy_GLS(titles,param,fullsol,agg_sol,init_guess)
%
% PROHOROVIT_SPLINE_RK_INDY_GLS performs the inverse problem for 2
% paramters for the spline case of the Prohorov Metric. This assumes the
% existence of the following items:
%
%    The solution you are trying to match (in this case it is called
%        K_dist1_dist2_joint_longtime_smooth_error_5.mat)
%
%    Precomputed solutions to your differential equation at all the
%    quadrature nodes (in this case called
%    'fullsol_2par_longtime_with100.mat')
%
% This Function creates the following files:
%    outname: For each iteration of the GLS algorithm, I save the current
%        weights, the previous weights, the soltuion to the differential
%        equation, and the weights for the GLS algorithm
%
%
% Written by Erica Rutter (July 2017)

%% initialization of the algorithm and data

% Define model parameters
nr = length(param.nodesr);   % nodes
nD = 1;
Nr=length(param.bsr);   % quadrature nodes
maxdim=max(nr,nD);
numpar=1;

% Determine if scaling the solutions in x and t.
[x,y,t]=size(agg_sol);
% Calculate the splines for parameter 1
intoverr = zeros(length(param.nodesr),x,y,t);
myintr = zeros(1,length(param.nodesr));
fullljsr=zeros(length(param.nodesr),Nr);
for j=1:length(param.nodesr)
    if j==1
        nodesj=[param.nodesr(j),param.nodesr(j),param.nodesr(j+1)];
    elseif j==length(param.nodesr)
        nodesj=[param.nodesr(j-1),param.nodesr(j),param.nodesr(j)];
    else
        nodesj=param.nodesr(j-1:j+1);
    end
    ljs=make_splines(nodesj,param.bsr);
    diags=diag(ljs);
    for ll=1:t
        for k=1:x
            intoverr(j,k,:,ll)=trapz(param.bsr,squeeze(fullsol(k,:,ll,:))*diags,2);
        end
    end
    myintr(j)=trapz(param.bsr,ljs);
    fullljsr(j,:)=ljs;
end

%% Constrained Optimization

% These are the constrains for the optimization problem to ensure that the
% sum of the weights/splines add to 1 and are nonnegative
Aeq=zeros(numpar,numpar*maxdim);
Aeq(1,1:nr)=myintr;
beq = ones(numpar,1);
lb = zeros(numpar*maxdim,1);
ub = zeros(numpar*maxdim,1);
ub(1:nr) = 10^9*length(param.nodesr);
r0=zeros(numpar*maxdim,1);

if init_guess==1
    % initial guess is random
    r0(1:nr)=rand(nr,1);
    dividers=Aeq*r0;
    r0(1:nr)=r0(1:nr)/dividers(1);                         
elseif init_guess==0
    % initial guess is uniform
    r0(1:nr)=1./(param.nodesr(end)-param.nodesr(1));
end


% initial solve to pass through 
intoverr=sparse(reshape(intoverr,[length(param.nodesr),numel(intoverr)/length(param.nodesr)]));
weights=ones(size(agg_sol));
optweight=r0;


% For the regular algorithm
errsy = @(r)errorfunc_discrete_sparse(r,intoverr,agg_sol,weights);
options = optimset('MaxIter',25000,'MaxFunEvals',500000,'Display', 'None');
[optweight,~,converge_flag,~]=fmincon(errsy, optweight, [],[],Aeq,beq,lb,ub,[],options);
if converge_flag==0
    disp('Did not converge')
end
% finds the error and best approximation for the optimized weights 
[current_err,best_approx] = errorfunc_discrete_sparse(optweight,intoverr,agg_sol,weights);


%% Saving Information
% If you did not converge, return

if (ii>=maxits)&&(parchange>partol)
    disp('Did not converge')
    converge_flag=0;
end
save(titles,'optweight','current_err','best_approx','fullljsr','converge_flag')

end

