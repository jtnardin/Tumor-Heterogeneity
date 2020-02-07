function plotter_disc(param,real_data,optweight)
% PLOTTER_DISC plots the best discrete approximation versus the correct
% underlying distribution (Fig 6 in the associated paper)
% 
%    INPUTS:
%        param: a variable that contains the following information:
%            -param.nodesr: the vector of nodes for parameter 'r' for the
%                optimization 
%        real_data: the filename where the 'ground truth' underlying 
%            distribution is 
%        optweight: the optimum
%
%    OUTPUTS:
%        A figure displaying the predicted best discrete approximation
%            versus the actual underlying distribution, as seen in Figure 6
%            in the associated tumor heterogeneity paper.
%
% Written by Erica Rutter (January 2020)


%% load the actual rho values and bin according to 'best fit'
load(real_data,'rho_vals');
hist_vals2=linspace(0.5*10^-3,1.5*10^-3,10);
hist_vals3=linspace(.02,.06,10);

% Find locations that match to the histogram locations
[~,a] = min(abs(hist_vals2(1)-param.nodesr));
[~,b] = min(abs(hist_vals2(end)-param.nodesr));
[~,c] = min(abs(hist_vals3(1)-param.nodesr));
[~,d] = min(abs(hist_vals3(end)-param.nodesr));

% find out how much proportion is given in the histogram intervals. This is
% to rescale solutions after interpolation
prop_2=sum(optweight(a:b));
prop_3=sum(optweight(c:d));

figure;
subplot(1,2,1)
histogram(rho_vals,hist_vals2,'facecolor',[.78,.74,.86])
% interpolate the optimal weights into the histogram bin intervals, then
% scale
zz=interp1(param.nodesr,optweight,hist_vals2);
zz=zz/sum(zz)*prop_2*length(rho_vals);
set(gca,'Fontsize',24,'linewidth',2)
ylabel('Frequency')
hold on
plot(hist_vals2,zz,'r*--','linewidth',2)
xlabel('\rho_0')

subplot(1,2,2)
histogram(rho_vals,hist_vals3,'facecolor',[.78,.74,.86])
% interpolate the optimal weights into the 2nd histogram bin intervals, then
% scale
zz=interp1(param.nodesr,optweight,hist_vals3);
zz=zz/sum(zz)*prop_3*length(rho_vals);
set(gca,'Fontsize',24,'linewidth',2)
hold on
plot(hist_vals3,zz,'r*--','linewidth',2)
legend('\rho_k',['Discrete' newline 'Approximation'])
xlabel('\rho_M')
