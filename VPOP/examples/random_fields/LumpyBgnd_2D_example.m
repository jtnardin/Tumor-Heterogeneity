% File: LumpyBgnd_2D_example.m
% Author: Nick Henscheid 
% Date: 4-2019, 5-2020
% Purpose: Demonstrates the lumpy background object functionality (d = 2)

% Initialize a default lumpy background object & display its properties
close all;
fprintf('Constructing a default LumpyBgnd and plotting it.\n'); 
L = LumpyBgnd('Support',RectSupport([0,1;0,1]));  
L.TurnOffWarnings;  % Don't show warning messages (comment out if you want to see them).
L.Kbar = 200;
L.N = 256;

% Display a realization 
plot(L);title('Realization of LumpyBgnd');

uans = input('Continue? (n for no, enter for yes) ','s'); 
if(strcmp(uans,'n'))
    close all;return;
end
%%
% Randomize the field and re-plot
fprintf('Randomizing the field with the same parameters\n');
L.Randomize;

plot(L);title('Realization of LumpyBgnd'); figure(gcf); 
uans = input('Continue? (n for no, enter for yes) ','s'); 
if(strcmp(uans,'n'))
    close all;return;
else 
    close all;
end
%%
% Change the lump width, but keep it isotropic and keep the same centers
% Note the warning displayed in the command window! The original lump
% centers were selected based on the original lump width 
fprintf('Changing the lump size and re-plotting (maintaining lump centers)\n'); 
fig = figure; fig.Position(3) = 1200; fig.Position(4) = 600; 
L.cov = 0.001;
subplot(1,3,1);
plot(L);title('L.cov = 1e-3','FontSize',16);

L.cov = 0.0001;
subplot(1,3,2);
plot(L);title('L.cov = 1e-4','FontSize',16);

L.cov = 0.1;
subplot(1,3,3);
plot(L);title('L.cov = 0.1','FontSize',16);
uans = input('Continue? (n for no, enter for yes) ','s'); 
if(strcmp(uans,'n'))
    close all;return;
else 
    close all;
end
%%
% Change the lump to be anisotropic
fprintf('Making the lumps anisotropic and re-randomizing\n'); 
fig = figure; fig.Position(3) = 1200; fig.Position(4) = 600; 
L.cov = diag([0.001,0.01]);
L.Randomize; 
subplot(1,2,1);
plot(L);title('L.cov = diag([1e-3,1e-2])','FontSize',16);
L.cov = diag([0.05,0.0005]);
L.Randomize; 
subplot(1,2,2);
plot(L);title('L.cov = diag([5e-2,5e-4])','FontSize',16);
uans = input('Continue? (n for no, enter for yes) ','s'); 
if(strcmp(uans,'n'))
    close all;return;
else 
    close all;
end
%%
fprintf('Rotating the lumps by a fixed rotation\n'); 
% Rotate lumps by a fixed rotation matrix 
theta = pi/4;  % Rotation angle
var1  = 0.01;   % var (rotated x-hat)
var2  = 0.001; % var (rotated y-hat)

R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
K = R*diag([var1,var2])*R';

L.cov = K;
L.Randomize; 
plot(L);title(sprintf('L.cov = [%1.1e,%1.1e;%1.1e,%1.1e]',L.cov(1,1),L.cov(1,2),L.cov(2,1),L.cov(2,2)),'FontSize',16);
uans = input('Continue? (n for no, enter for yes) ','s'); 
if(strcmp(uans,'n'))
    close all;return;
else 
    close all;
end
%% Make the lumps non-uniform in size (randomly generated isotropic covariance)
fprintf('Making the lumps non-uniform in size (random isotropic lump covariance)\n'); 
lump_mean = log(0.07);
lump_std  = log(0.12);
L.cov = exp(lump_mean+lump_std*randn())*eye(2);  % Lognormally distributed
L.b   = 1/sqrt((2*pi)^2*det(L.cov));  %PDF/Gaussian mixture model scaling
for i=1:L.K-1
    covtemp = exp(lump_mean+lump_std*randn())*eye(2);
    L.cov = [L.cov,covtemp];
    L.b  = [L.b,1/sqrt((2*pi)^2*det(covtemp))];
end

plot(L);title('Indepdent cov for ea. lump','FontSize',16);
fprintf('Done!\n'); 

