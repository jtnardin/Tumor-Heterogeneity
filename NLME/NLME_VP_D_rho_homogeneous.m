%fisher_NLME written 3-20-18 by JTN to peformed NLME on a dataset of
%several fisher KPP simulaions with 2 subgroups

clear all; clc

disp('Simulation started ')
disp(datetime('now'))

%load in data
%
%fisher with rho and D varying
load('bimodal_rho_D_const_VP.mat')
%
%virtual population data
%load('bimodal_RDE_VP.mat')


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


% % % % % uncomment this to make sure X_m,Y_m,T_m and y_vec are formatted correctly
% % % % % if the plots make sense (each subplotcorresponds to a patient data)
% % % % % then it should be correct
% % % % %make sure dimensions correct
% % % figure
% % % j=1;
% % % y_vec1 = y_vec(NUMS==j);
% % % Tm_vec1 = T_m(NUMS==j);
% % % for i = 1:length(t_uniqe)
% % %      subplot(length(t_uniqe),1,i)
% % %      t_ind = Tm_vec1 == t_uniqe(i);
% % %      yplot = reshape(y_vec1(t_ind),[length(x_uniqe),length(y_uniqe)]);
% % %      surf(x_uniqe,y_uniqe,yplot,'edgecolor','none')
% % %      view(2)
% % % end
% % % title(['patient ' num2str(j) ', M = ' num2str(metast_ind(j))])


%fixed design matrix
for i = 1:N
    A(:,:,i) = [1 metast_ind(i) 0 0 ;0 0 1 metast_ind(i) ];
end
%random design matrix is the same
B = A;


% % initial guesses for parameter estimates
beta0 = [1e-6, 0, 0.02,.02]';
opt = statset('TolFun',1e0);
model = @(phi,t) fisher_2d_sim(phi,t);
tic
[beta,PSI,stats,b] = nlmefit(X,y_vec,NUMS,[],model,beta0,'FEGroupDesign',A,...
    'REGroupDesign',B,'Options',opt);
toc

save('VP_results_phi_const_D_const.mat','beta','PSI','stats','b')