addpath('~/Documents/matlab_custom/')
set(0,'defaultfigurecolor',[1 1 1])

%load in patient data
load('bimodal_rho_const_VP.mat')

%load in results from NLME
load('VP_results_phi_const_D_known_2.mat')


benign_data = rho_vals(metast_ind==0);
metast_data = rho_vals(metast_ind==1);

benign_num = numel(benign_data);
metast_num = numel(metast_data);

benign_domain = linspace(.9*min(benign_data),1.1*max(benign_data),100);
metast_domain = linspace(.9*min(metast_data),1.1*max(metast_data),100);

benign_pdf = normpdf(benign_domain,beta(1),sqrt(PSI(1,1)));
metast_pdf = normpdf(metast_domain,beta(1),sqrt(PSI(1,1))) ...
           + normpdf(metast_domain,beta(2),sqrt(PSI(2,2)));
       

benign_pdf_true = normpdf(benign_domain,r0,r0sigma);
metast_pdf_true = normpdf(metast_domain,r0,r0sigma) ...
           + normpdf(metast_domain,rM,rMsigma);
       
figure();
f=subplot(1,2,1);

left_color = [0 0 0];
right_color = [0 0 0];
set(f,'defaultAxesColorOrder',[[1 0 0]; right_color])

yyaxis left
hold on
plot(benign_domain,benign_pdf/500,'r--','linewidth',2)
plot(benign_domain,benign_pdf_true/500,'-','color',[0 .5 0],'linewidth',2)
ylabel('Frequency','interpreter','latex')

yyaxis right
hold on
hist(benign_data);
h = findobj(gca,'Type','patch');
h.FaceAlpha = 0.125;
yticks([])
line([0.5e-3,1.5e-3],[4 4],'color','k')
axis([0.5e-3 1.5e-3 0 4])
set(gca,'FontSize',20)
xlabel('$\rho_0$','interpreter','latex')
title(['Benign Patients, N = ' num2str(benign_num)],'interpreter','latex')

g=legend('Inferred Distribution','True Distribution','{\boldmath$\rho_k$}','location','southwest');
set(g,'interpreter','latex')

ax = gca;
ax.YAxis(1).Color = [0,0,0];
ax.YAxis(2).Color = [0,0,0];


f = subplot(1,2,2);
set(f,'defaultAxesColorOrder',[[1 0 0]; right_color])

yyaxis left
hist(metast_data)
h = findobj(gca,'Type','patch');
h.FaceAlpha = 0.125;
axis([0.02 .055 0 3])
yyaxis right
hold on
plot(metast_domain,metast_pdf,'r--','linewidth',2)
plot(metast_domain,metast_pdf_true,'-','color',[0 .5 0],'linewidth',2)
yticks([])
set(gca,'defaultAxesColorOrder',[left_color; right_color])

xlabel('$\rho_M$','interpreter','latex')
yyaxis left
ylabel('Frequency','interpreter','latex')
title(['Malignant Patients, N = ' num2str(metast_num)],'interpreter','latex')
set(gca,'FontSize',20)
ax = gca;
ax.YAxis(1).Color = [0,0,0];
ax.YAxis(2).Color = [0,0,0];

%saveas(gcf,'NLME_plot.png','png')
export_fig 'NLME_plot_ef.png' 
savefig('NLME_plot.fig')
