function forwardsolvers(param,precomputed_sol,X)
%
% FORWARDSOLVERS performs the forward solutions over the parameter vector 
% of parameter 1 and parameter 2 for comparing with our actual simulated
% solutions
% 
%    INPUTS:
%        param: a variable that contains the following information:
%            -param.xs: the x vector for solutions
%            -param.ys: the y vector for solutions
%            -param.t:  the t vector for solutions
%            -param.bsr: vector of nodes for parameter 'r' on which to solve
%                the pde
%        X: a matrix of times, xs, and ys for the solution
%        precomputed_sol: the savename of the solutions
%
%    OUTPUTS:
%        precomputed_sol: a mat-file with the full forward solutions over  
%             the parameter space of interest (named 'fullsol')
%
% Written by Erica Rutter (July 2019)

fullsol=zeros(length(param.xs),length(param.ys),length(param.t),length(param.bsr));
for k =1:length(param.bsr)
    yd=fisher_2d_sim([param.D,param.bsr(k)],X);
    Y = reshape(yd,length(param.xs),length(param.ys),length(param.t));
    fullsol(:,:,:,k)=Y;
    k
end
save(precomputed_sol,'fullsol')


