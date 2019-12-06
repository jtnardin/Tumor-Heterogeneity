%fisher_NLME written 12-6-2019 to peformed NLME on a dataset of
%several fisher KPP simulaions with 2 subgroups, benign and malignant, as
%part of "A tutorial Review of Mathematical Techniques for Quantifying 
%Tumor Heterogeneity" by Everett et al. 

clear all; clc

disp('Simulation started ')
disp(datetime('now'))

%load in VP data with rho varying
load('bimodal_rho_const_VP.mat')


%%% # patients
N = size(Y,1);
%%% # measurements per patient
M = size(Y,2);

%spatiotemporal grid
X = repmat([T_m(:) X_m(:) Y_m(:)],N,1);


%1d time and space grids
t_uniqe = unique(T_m);
x_uniqe = unique(X_m);
y_uniqe = unique(Y_m);

%reconfigure y to y_vec so that each row of
%[X y_vec] corresponds to : [t_i x_j y(t_i,x_j)] for some patient
y_vec = Y';
y_vec = y_vec(:);


%label each row of y based on which patient it belongs to
for i = 1:N
    if i == 1
        NUMS = ones(numel(X_m),1);
    else
        NUMS = [NUMS; i*ones(numel(X_m),1)];
    end
end



%fixed design matrix
for i = 1:N
    A(:,:,i) = [1 metast_ind(i) 0 0 ;0 0 1 metast_ind(i) ];
end
%random design matrix is the same
B = A;

% initial guesses for parameter estimates
beta0 = [1e-6, 0, 0.02,.02]';

% If the tolerance is set too low, then the computation may take a long
% time. the current setting (1e0) takes about 4 days on one CPU.
opt = statset('TolFun',1e0);

%the model is the 2d Fisher-KPP
model = @(phi,t) fisher_2d_sim(phi,t);

tic
[beta,PSI,stats,b] = nlmefit(X,y_vec,NUMS,[],model,beta0,'FEGroupDesign',A,...
    'REGroupDesign',B,'Options',opt);
toc

%save 
save('VP_results_rho_const_D_known.mat','beta','PSI','stats','b')