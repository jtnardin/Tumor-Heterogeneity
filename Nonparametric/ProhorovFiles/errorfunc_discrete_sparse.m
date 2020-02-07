function [err,best_approx]=errorfunc_discrete_sparse(r,totsol,meansol,weights)
% ERROR_DISCRETE_SPARSE computes the error between the spatiotemporal
% dataset and the proposed solution for generalized least squares.
% 
%    INPUTS:
%        r: the weights for the parameter mesh 'r'
%        totsol: the weighted solution (different for spline and discrete)
%        meansol: the dataset we are trying to fit
%        weights: these are weights for the error fucntion using least
%            squares. For ordinarly least squares, weights=1.
%
%    OUTPUTS:
%        err: the OLS/GLS error comparing the computed solution and the
%            dataset
%        best_approx: the aggregate solution for the input 'r' 
%
% Written by Erica Rutter (July 2019)

[a,b,c]=size(meansol);

best_approx=reshape(r'*totsol,[a,b,c]);
err=sum(sum(sum(weights.*((meansol-best_approx).^2))));
