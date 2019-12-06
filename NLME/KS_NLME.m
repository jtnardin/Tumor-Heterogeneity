%KS_NLME written 12-6-2019 to compute the KS test for "A tutorial 
%"Review of Mathematical Techniques for Quantifying Tumor Heterogeneity" 
%by Everett et al. 

%load in patient data
load('bimodal_rho_const_VP.mat')

%load in results from NLME
load('VP_results_rho_const_D_known.mat')

% number benign, maligant
N = 15;
M = 10;

%inferred distn
NLME_dist = [beta(1) + sqrt(PSI(1,1))*randn(N,1);beta(1) + beta(2) + (sqrt(PSI(1,1)) + sqrt(PSI(2,2)))*randn(M,1)];

%true distn
true_dist = [r0 + r0sigma*randn(N,1);r0 + rM + rMsigma*randn(M,1)];

%KS Test
[h,p,ks2stat] = kstest2(NLME_dist,true_dist)

%plot empirical CDFs
ecdf(NLME_dist)
hold on
ecdf(true_dist)