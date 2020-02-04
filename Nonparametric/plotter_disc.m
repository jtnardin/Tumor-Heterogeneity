function plotter_disc(param,fullsol,agg_sol,real_data,optweight)
% Driver for Optimization Stuff %

%% load the acutla rho values and bin according to 'best fit'
load(real_data,'rho_vals');
hist_vals1=0;
hist_vals2=linspace(0.5*10^-3,5.5*10^-3,50);
hist_vals3=linspace(.02,.06,10);
real_rho2=histc(rho_vals,hist_vals2);
real_rho3=histc(rho_vals,hist_vals3);
histy=histc(optweight,hist_vals2)

%real_rho=real_rho/sum(real_rho);
figure;
histogram(rho_vals,hist_vals2)
%plot(param.nodesr,real_rho,'k+-')
%hold on
%plot(param.nodesr,optweight,'r*-')
keyboard
%% Plot the pdfs
figure;
%subplot(2,1,1)
[hAx,hLine1,hLine2] =plotyy(histy,real_rho,param.nodesr,optweight);
xlim([min(param.nodesr),max(param.nodesr)])
xlim(hAx(2), [min(param.nodesr),max(param.nodesr)])
set(hAx(2),'XTickLabel',[])
set(hAx(2),'FontSize',24)
%set(hAx(2),'xlim',([
legend('Actual','Estimated','Location','NorthEast')
set(gca,'Fontsize',24,'linewidth',2)
%title([distributions{dist2},' with ',int2str(nr),' nodes'])
xlabel('\rho','Fontsize',24)
ylabel('Probability','Fontsize',24)
set(gca,'Fontsize',24,'linewidth',2)
set(hLine2,'LineStyle','--','Marker','+','LineWidth',2,'color',[0.5,0.1,0.1])
%ylim(hAx(2),[0,max(max(optweights_disc(:,1:best_nr_disc,1)))])
set(hLine1,'LineStyle','--','Marker','*','LineWidth',2)


%title(strcat('D ',distributions{dist1},' r ',distributions{dist2},' with ',int2str(100*err_level),'% error, discrete 1 xscale ',int2str(x_scale),' tscale ',int2str(t_scale)))
keyboard;
titles=strcat('Figures/Best_discrete_',int2str(1),'_pdf');
print(titles,'-dpdf')
savefig(strcat(titles,'.fig'))
%% Plot the solutions if wanted

%if do_plot==1
%    plot_fits(persol,fullsol,optweights,param,t_scale,x_scale)
%end
%% perform UQ 