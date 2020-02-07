function [optweight,current_err,converge_flag]=prohorovit_spline_r_indy_GLS(titles,param,fullsol,agg_sol,init_guess)
%
% PROHOROVIT_SPLINE_RK_INDY_GLS finds the optimal weights for the inverse 
% problem for the parameter 'r' in the reaction-diffusion equation under 
% the Prohorov Metric Framework using continuous/spline approximations.
%
%    INPUTS:
%        title: the title of the matlab file to save information to
%        param: a variable that contains the following information:
%            -param.nodesr: the vector of nodes for parameter 'r' for the
%                optimization 
%            -param.bsr: the computational/numerical mesh of the parameter 
%                'r' (for the precomputed solutions)
%        fullsol: the precomputed solution 
%        agg_sol: the aggregate 'data' we are comparing with
%        init_guess: determining whether use a uniform distribution as an
%            initial guess at parameter distribution (1) or use random 
%            numbers (0)
%
%    OUTPUTS:
%        optweight: the optimal weights for the parameter mesh 'r'
%        current_err: the value of the error (OLS/GLS) for the optimal
%            weighted parameter 'r'
%        converge_flag: the converge flag given out by fmincon (0 signifies
%            minima found)
%        title: a mat-file with the optimal weight (optweight), the best
%            error (current_err), the spatiotemporal solution using
%            the optweights, and the converge flag for the optimization.
%

% Written by Erica Rutter (July 2017)

%% initialization of the algorithm and data

% Define model parameters
nr = length(param.nodesr);   % nodes
Nr=length(param.bsr);   % quadrature nodes

% Determine if scaling the solutions in x and t.
[x,y,t]=size(agg_sol);

% Calculate the splines for the parameter 'r' over the mesh
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
Aeq=zeros(1,nr);
Aeq(1,:)=myintr;
beq = ones(1,1);

% specify bounds. must be nonnegative 
lb = zeros(nr,1);
ub = zeros(nr,1);
ub(1:nr) = 10^9;
r0=zeros(nr,1);

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

weights=ones(size(agg_sol)); %Assuming an ordinary least squars approach

% Fit the parameters
errsy = @(r)errorfunc_discrete_sparse(r,intoverr,agg_sol,weights);
options = optimset('MaxIter',25000,'MaxFunEvals',500000,'Display', 'None');
[optweight,~,converge_flag,~]=fmincon(errsy, r0, [],[],Aeq,beq,lb,ub,[],options);
if converge_flag==0
    disp('Did not converge')
end
% finds the error and best approximation for the optimized weights 
[current_err,best_approx] = errorfunc_discrete_sparse(optweight,intoverr,agg_sol,weights);

%% Saving Information
%save(titles,'optweight','current_err','best_approx','fullljsr','converge_flag')

end

