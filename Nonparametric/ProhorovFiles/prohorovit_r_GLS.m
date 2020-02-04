function [optweight,current_err,converge_flag]= prohorovit_r_GLS(titles,param,fullsol,persol,init_guess)
%
% PROHOROVIT_2PAR_GLS performs the inverse problem for 2
% paramters for the discrete case of the Prohorov Metric. This assumes the
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

%% Initialize the inverse problem

% define model parameters such as the nodes 
numpar=1;
nr = length(param.nodesr);   % nodes
nD = 1;
maxdim=max(nr,nD);

%% Constrained Optimization

% Constraints -- ensure that sum(ri)=1, sum(Ki)=1. Note: fmincon requires A
% to have same nubmer of r0 as a columna vector
Aeq=zeros(numpar,numpar*maxdim);
Aeq(1,1:nr)=1;
beq = ones(numpar,1);
lb = zeros(numpar*maxdim,1); % weights cannot be nonnegative
ub = zeros(numpar*maxdim,1);
ub(1:nr) = 1;
r0=zeros(numpar*maxdim,1);
if init_guess==1
    r0(1:nr)=rand(nr,1);
    r0(1:nr)=r0(1:nr)/sum(r0(1:nr));                         % initial guess is random.
elseif init_guess==0 % uniform initial guess
    r0(1:nr)=1/nr;
end

% permute the solution to be in correct order for matrix multiplication
fullsol=permute(fullsol,[4,1,2,3]);
fullsol=reshape(fullsol,[length(param.nodesr),numel(fullsol)/length(param.nodesr)]);

weights=ones(size(persol)); % for the statistical error model.
optweight=r0;

errsy = @(r)errorfunc_discrete_sparse(r,fullsol,persol,weights);
options = optimset('MaxIter',25000,'MaxFunEvals',50000,'Display', 'Off');
[optweight,~,converge_flag,~]=fmincon(errsy, optweight, [],[],Aeq,beq,lb,ub,[],options);
if converge_flag==0
    disp('Did not converge')
end
[current_err,best_approx] = errorfunc_discrete_sparse(optweight,fullsol,persol,weights);

%% Save output 
if (ii>=maxits)&&(parchange>partol)
    disp('Did not converge')
    converge_flag=0;
end
save(titles,'optweight','current_err','best_approx','converge_flag')

end


