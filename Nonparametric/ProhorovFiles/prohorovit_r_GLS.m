function [optweight,current_err,converge_flag]= prohorovit_r_GLS(titles,param,fullsol,agg_sol,init_guess)
%
% PROHOROVIT_R_GLS finds the optimal weights for the inverse problem for 
% the parameter 'r' in the reaction-diffuion equation under the Prohorov
% Metric Framework using discrete approximations.
%
%    INPUTS:
%        title: the title of the matlab file to save information to
%        param: a variable that contains the following information:
%            -param.nodesr: the vector of nodes for parameter 'r' for the
%                optimization 
%        fullsol: the precomputed solution for the current number of
%            parameter nodes
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

%% Initialize the inverse problem

% define model parameters such as the number of nodes estimated 
nr = length(param.nodesr); 

%% Constrained Optimization

% Constraints -- ensure that sum(ri)=1. Note: fmincon requires A
% to have same number of r0 as a column vector
Aeq=zeros(1,nr);
Aeq(1,:)=1;
beq = ones(1,1);

% specify bounds. must be nonnegative and cant be >1
lb = zeros(nr,1); 
ub = zeros(nr,1); 
ub(1:nr) = 1;
r0=zeros(nr,1);

if init_guess==1
    % initial guess is random.
    r0(1:nr)=rand(nr,1);
    r0(1:nr)=r0(1:nr)/sum(r0(1:nr));                         
elseif init_guess==0 
    % uniform initial guess
    r0(1:nr)=1/nr;
end

% permute the solution to be in correct order for matrix multiplication
fullsol=permute(fullsol,[4,1,2,3]);
fullsol=reshape(fullsol,[length(param.nodesr),numel(fullsol)/length(param.nodesr)]);

weights=ones(size(agg_sol)); % Assuming ordinary least squares (OLS)

errsy = @(r)errorfunc_discrete_sparse(r,fullsol,agg_sol,weights);
options = optimset('MaxIter',25000,'MaxFunEvals',50000,'Display', 'Iter');
[optweight,~,converge_flag,~]=fmincon(errsy, r0, [],[],Aeq,beq,lb,ub,[],options);
if converge_flag==0
    disp('Did not converge')
end
[current_err,~] = errorfunc_discrete_sparse(optweight,fullsol,agg_sol,weights);
%% Save output 
%save(titles,'optweight','current_err','best_approx','converge_flag')

end


