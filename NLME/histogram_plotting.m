%histogram_plotting written 12-6-2019 to plot figure 5 from "A tutorial 
%"Review of Mathematical Techniques for Quantifying Tumor Heterogeneity" 
%by Everett et al. 

%background color
set(0,'defaultfigurecolor',[1 1 1])

%load in patient data
load('bimodal_rho_const_VP.mat')

%load in results from NLME
load('VP_results_rho_const_D_known.mat')


% benign and malignant datasets
benign_data = rho_vals(metast_ind==0);
metast_data = rho_vals(metast_ind==1);

% number of benign, malignant tumors
benign_num = numel(benign_data);
metast_num = numel(metast_data);

% domains for plotting
benign_domain = linspace(.9*min(benign_data),1.1*max(benign_data),100);
metast_domain = linspace(.9*min(metast_data),1.1*max(metast_data),100);

%inferred distributions
benign_pdf = normpdf(benign_domain,beta(1),sqrt(PSI(1,1)));
metast_pdf = normpdf(metast_domain,beta(1),sqrt(PSI(1,1))) ...
           + normpdf(metast_domain,beta(2),sqrt(PSI(2,2)));
       
%true distributions
benign_pdf_true = normpdf(benign_domain,r0,r0sigma);
metast_pdf_true = normpdf(metast_domain,r0,r0sigma) ...
           + normpdf(metast_domain,rM,rMsigma);
       
figure();

%plot benign data
f=subplot(1,2,1);

%distributions
yyaxis left
hold on
plot(benign_domain,benign_pdf/500,'r--','linewidth',2)
plot(benign_domain,benign_pdf_true/500,'-','color',[0 .5 0],'linewidth',2)
ylabel('Frequency','interpreter','latex')
yticks([0,1,2,3,4])

%histogram
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
title(['Benign Tumors, N = ' num2str(benign_num)],'interpreter','latex')

%change axes to black
ax = gca;
ax.YAxis(1).Color = [0,0,0];
ax.YAxis(2).Color = [0,0,0];

%malignant tumors
f = subplot(1,2,2);

%distributions
yyaxis left
hold on
plot(metast_domain,metast_pdf,'r--','linewidth',2)
plot(metast_domain,metast_pdf_true,'-','color',[0 .5 0],'linewidth',2)
yticks(linspace(0,4/3*max(metast_pdf),5))
yticklabels([0,1,2,3,4])
axis([0.02 .055 0 4/3*max(metast_pdf)])

%histogram
yyaxis right
hist(metast_data)
h = findobj(gca,'Type','patch');
h.FaceAlpha = 0.125;
axis([0.02 .055 0 4])
yticks([])
xlabel('$\rho_M$','interpreter','latex')
yyaxis left
ylabel('Frequency','interpreter','latex')
title(['Malignant Tumors, N = ' num2str(metast_num)],'interpreter','latex')
set(gca,'FontSize',20)

%add black horizontal line on top
plot([.02 .055],4*ones(2,1),'k-')

%legend
g=legend('Inferred Distribution','True Distribution','{\boldmath$\rho_k$}');
set(g,'interpreter','latex')

%make axes black
ax = gca;
ax.YAxis(1).Color = [0,0,0];
ax.YAxis(2).Color = [0,0,0];

%save
export_fig 'NLME_plot.png' 
savefig('NLME_plot.fig')
