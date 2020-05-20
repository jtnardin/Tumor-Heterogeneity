% File: UnifRandField_example.m
% Author: Nick Henscheid 
% Date: 5-2020
% Purpose: Demonstrates the UnifRandField object functionality (d = 2)

close all

%% Test the basic functionality
fprintf('Plotting a single realization of the default UnifRandField\n');
U = UnifRandField('N',64);  
plot(U); set(gca,'CLim',[0,1]); 
title(sprintf('Realization of a uniform random field; $U(x)\\equiv$ %2.4e',U.c),'Interpreter','latex','FontSize',16); 
uans = input('Continue? (y/n) or enter for yes ','s'); 
if(strcmp(uans,'n'))
    close all;return;
else 
    close all;
end
%% Plot some realizations
fprintf('Plotting multiple realization of the default UnifRandField\n');
fig1 = figure; 
for i=1:8
    U.Randomize; 
    subplot(2,4,i); 
    plot(U); set(gca,'CLim',[0,1]); 
    title(sprintf('$U(x)\\equiv$ %2.4e',U.c),'Interpreter','latex','FontSize',16); 
end
uans = input('Continue? (y/n) or enter for yes ','s'); 
if(strcmp(uans,'n'))
    close all;return;
else
    close all;
end
%% Estimate the mean field, compute a histogram estimate of the one-point PDF
fprintf('Estimating the mean field and computing a histogram estimate of the one-point PDF\n');
nsamp = 256; 
Y  = U.Sample(nsamp); 
mu = zeros(U.N); 
y  = zeros(nsamp,1); 
for i=1:nsamp
    mu = mu + Y{i}; 
    y(i) = Y{i}(1,1); % For histogram estimate
end
mu = mu/nsamp; 
fig2 = figure; fig2.Position(3) = 1000; fig2.Position(4) = 500; 
x = linspace(0,1,U.N); 
subplot(1,2,1); 
imagesc(x,x,mu); set(gca,'YDir','normal'); axis image; 
title(sprintf('$\\mu(x)\\equiv$ %2.4e',mu(1,1)),'Interpreter','latex','FontSize',16); 
subplot(1,2,2); 
histogram(y,'Normalization','pdf'); 
title('Histogram estimate of the one-point PDF','FontSize',16); 
uans = input('Continue? (y/n) or enter for yes ','s'); 
if(strcmp(uans,'n'))
    close all;return;
else 
    close all;
end
%% Change the distribution and repeat the last 2 blocks 
fprintf('Changing the PDF and repeating the sample, mean field and histogram estimate\n');
mu0 = 1;  sig0 = 0.25; 
U.dist = @()mu0+sig0*randn(); 

fig3 = figure; 
for i=1:8
    U.Randomize; 
    subplot(2,4,i); 
    plot(U); set(gca,'CLim',[0,2]); 
    title(sprintf('$U(x)\\equiv$ %2.4e',U.c),'Interpreter','latex','FontSize',16); 
end

nsamp = 256; 
Y2  = U.Sample(nsamp); 
mu2 = zeros(U.N); 
y2  = zeros(nsamp,1); 
for i=1:nsamp
    mu2 = mu2 + Y2{i}; 
    y2(i) = Y2{i}(1,1); % For histogram estimate
end
mu2 = mu2/nsamp; 
fig4 = figure; fig4.Position(3) = 1200; fig4.Position(4) = 600; 
subplot(1,2,1); 
x = linspace(0,1,U.N); 
imagesc(x,x,mu2); set(gca,'YDir','normal'); axis image; 
title(sprintf('$\\mu(x)\\equiv$ %2.4e',mu2(1,1)),'Interpreter','latex','FontSize',16); 

subplot(1,2,2); 
histfit(y2); 
title('Histogram estimate of the one-point PDF','FontSize',16); 
str = sprintf('$\\hat{\\mu}=$%2.2e\n$\\hat{\\sigma}=$%2.2e',mean(y2),std(y2)); 
annotation('textbox',[0.8,0.3,0.4,0.4],'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',16);
fprintf('Done!\n'); 