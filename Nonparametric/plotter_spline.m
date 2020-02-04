function plotter_spline(param,fullsol,agg_sol,real_data,optweight)
% Driver for Optimization Stuff %

% Define model parameters
nr = length(param.nodesr);   % nodes
nD = 1;
Nr=length(param.bsr);   % quadrature nodes
maxdim=max(nr,nD);
numpar=1;
err_level=0;
load(real_data,'rho_vals');
real_rho=histc(rho_vals,param.nodesr);
real_rho=real_rho/sum(real_rho);

% Determine if scaling the solutions in x and t.
[x,y,t]=size(agg_sol);
% Calculate the splines for parameter 1
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

%% Plot the pdfs

figure;
thepdfr=optweight'*fullljsr;
 % if underlying distribution was discrete
hist_vals1=0;
hist_vals2=linspace(0.5*10^-3,1.5*10^-3,10);
hist_vals3=linspace(.02,.06,10);
histogram(rho_vals,hist_vals2)
keyboard;
 
[hAx,hLine1,hLine2]=plotyy(param.nodesr,real_rho,param.bsr,thepdfr);
xlim([min(param.nodesr),max(param.nodesr)])
xlim(hAx(2), [min(param.bsr),max(param.bsr)])
ylim([0,max(real_rho)])
ylim(hAx(2),[0,max(thepdfr)])
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

titles=strcat('Figures/Best_discrete_',int2str(0),'_pdf');
print(titles,'-dpdf')
savefig(strcat(titles,'.fig'))
%% Plot the solutions if wanted

%if do_plot==1
%    plot_fits(persol,fullsol,optweights,param,t_scale,x_scale)
%end
%% perform UQ 