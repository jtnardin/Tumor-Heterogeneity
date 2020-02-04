function plotter_disc(dist1,dist2,err_level,x_scale,t_scale)
% Driver for Optimization Stuff %

%% declare variables
load('DataGeneration/distributions.mat','distributions');
global distributions

Nr = 50; % quadrature nodes (if using spline)
ND = 50; % quadratrue nodes (if using spline)
dmin=0; % min for parameter 1
dmax=2e-4; %max for parameter 1
rmin=0; % min for parameter 2
rmax=2; % max for parameter 2
num_inits=100; % number of initializations to perform.
do_plot=0; % this is to plot results

%% Ensure all Runs have been completed before 
simulated_data = strcat('DataGeneration/D_',distributions{dist1},'_r_',distributions{dist2},'_joint_error_',int2str(100*err_level),'.mat');
if ~exist(simulated_data,'file')
    disp('Simulated Data Does Not Exist')
    %return;
end
BestM_value_spline=strcat('Outputs/D_',distributions{dist1},'_r_',distributions{dist2},'_joint_error_',int2str(100*err_level),'_discrete_0_xscale_',int2str(x_scale),'_tscale_',int2str(t_scale),'_bestMs.mat');
BestM_value_disc=strcat('Outputs/D_',distributions{dist1},'_r_',distributions{dist2},'_joint_error_',int2str(100*err_level),'_discrete_1_xscale_',int2str(x_scale),'_tscale_',int2str(t_scale),'_bestMs.mat');

if ~exist(BestM_value_disc,'file')
    disp('Best M vlaues for discrete does not exist')
    %return;
elseif ~exist(BestM_value_spline,'file')
    disp('Best M vlaues for splines does not exist')
    %return;
end

load(BestM_value_disc,'best_nr','best_nD')
best_nr_disc=best_nr;
best_nD_disc=best_nD;

load(BestM_value_spline,'best_nr','best_nD')
best_nr_spline=best_nr;
best_nD_spline=best_nD;

titler_spline=strcat('Outputs/D_',distributions{dist1},'_nD',int2str(best_nD_spline),'_ND_',int2str(ND),'_r_',distributions{dist2},'_nr',int2str(best_nr_spline),'_Nr_',int2str(Nr),'_longtime_error_',int2str(100*err_level),'_xscale_',int2str(x_scale),'_tscale_',int2str(t_scale),'_randinits.mat');
titler_disc=strcat('Outputs/D_',distributions{dist1},'_nD',int2str(best_nD_disc),'_r_',distributions{dist2},'_nr',int2str(best_nr_disc),'_longtime_error_',int2str(100*err_level),'_xscale_',int2str(x_scale),'_tscale_',int2str(t_scale),'_randinits.mat');
if ~exist(titler_spline,'file')
    disp('Random Initializations Do Not Exist for Spline Case')
    %return;
elseif ~exist(titler_disc,'file')
    disp('Random Initializations Do Not Exist for Discrete Case')
    %return;
end
load(simulated_data,'rs','Ds','weightsr','weightsD','persol')
load(titler_disc)
nodesr_disc = linspace(rmin,rmax,best_nr_disc);
nodesD_disc = linspace(dmin,dmax,best_nD_disc);
optweights_disc = optweights;
maxdim=max(best_nr_disc,best_nD_disc);
optweights_disc = reshape(optweights_disc,100,maxdim,2);
%% Plot the pdfs
figure;
subplot(2,1,1)
[hAx,hLine1,hLine2] =plotyy(rs,weightsr,nodesr_disc,optweights_disc(:,1:best_nr_disc,1));
xlim([0,2])
xlim(hAx(2), [0,2])
set(hAx(2),'XTickLabel',[])
set(hAx(2),'FontSize',24)
%set(hAx(2),'xlim',([
legend('Actual','Estimated','Location','NorthEast')
set(gca,'Fontsize',24,'linewidth',2)
%title([distributions{dist2},' with ',int2str(nr),' nodes'])
xlabel('\rho','Fontsize',24)
ylabel('Probability','Fontsize',24)
set(gca,'Fontsize',24,'linewidth',2)
xlim([nodesr_disc(1),nodesr_disc(end)])
set(hLine2,'LineStyle','--','Marker','+','LineWidth',1,'color',[0.5,0.1,0.1])
ylim(hAx(2),[0,max(max(optweights_disc(:,1:best_nr_disc,1)))])
if dist2==4 || dist2==5
    set(hLine1,'LineStyle','None','Marker','*')
    ylim(hAx(1),[0,1])
else
    set(hLine1,'LineWidth',2)
end
subplot(2,1,2)
[hAx,hLine1,hLine2] =plotyy(Ds,weightsD,nodesD_disc,optweights_disc(:,1:best_nD_disc,2));
xlim([0,2e-4])
xlim(hAx(2), [0,2e-4])
set(hAx(2),'XTickLabel',[])
set(hAx(2),'FontSize',24)
set(gca,'Fontsize',24,'linewidth',2)
xlabel('D','Fontsize',24)
set(gca,'Fontsize',24,'linewidth',2)
xlim([nodesD_disc(1),nodesD_disc(end)])
set(hLine2,'LineStyle','--','Marker','+','LineWidth',1,'color',[0.5,0.1,0.1])
ylim(hAx(2),[0,max(max(optweights_disc(:,1:best_nD_disc,2)))])
if dist1==4 || dist1==5
    set(hLine1,'LineStyle','None','Marker','*')
    ylim(hAx(1),[0,1])

else
    set(hLine1,'LineWidth',2)
end
keyboard;
%title(strcat('D ',distributions{dist1},' r ',distributions{dist2},' with ',int2str(100*err_level),'% error, discrete 1 xscale ',int2str(x_scale),' tscale ',int2str(t_scale)))
title(strcat('D ',distributions{dist1},' r ',distributions{dist2},' with ',int2str(100*err_level),'% error, discrete'))

titles=strcat('Figures/D_',distributions{dist1},'_r_',distributions{dist2},'_joint_error_',int2str(100*err_level),'_discrete_',int2str(1),'_xscale_',int2str(x_scale),'_tscale_',int2str(t_scale),'_pds');
print(titles,'-dpdf')
savefig(strcat(titles,'.fig'))
%% Plot the solutions if wanted

%if do_plot==1
%    plot_fits(persol,fullsol,optweights,param,t_scale,x_scale)
%end
%% perform UQ 