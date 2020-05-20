% File: RDE_Virtual_Pop_2D.m
% Author: Nick Henscheid 
% Date: 4-2019, 5-2020
%
% Purpose: Solves the 2D reaction diffusion (Fisher-Kolmogorov) problem
%          u_t = div(D*grad(u)) + rho*u*(1-u/k) using the method of lines 
%          and random field parameters (D,rho and k).
% Notes:   The Laplacian is discretized using the standard 2nd order 
%          5-point stencil. The resulting ODE system is solved using ode45.
%          Both of these discreizations are sub-optimal, but tend to work
%          in reasonable parameter regimes.
%          Units of u(x,t) are cells/cm^2 with t in days (2D)
%          Units of D are thus cm^2/day
%          Units of alpha are cells/day
%          Units of k are cells/cm^2 

% See documentation for full details. 

%%
close all
plt = input('Would you like to display the figures? (y/n) ','s');
if(strcmp(plt,'y'))
    plotting = 1;  % Display plots during simulation?
else
    plotting = 0; 
end
%% Set up sample storage
% Number of samples to generate
Nsamp = input('How many samples would you like to generate?\n(Hit Enter for default (4 samples)) '); 
if(isempty(Nsamp))
    Nsamp = 4;   % Number of samples to generate
end

%% Set up spatial grid
Ngrid = input('How many spatial grid points (in ea. direction) would you like to use?\n(Hit Enter for default (256)) ');
if(isempty(Ngrid))
    Ngrid = 256;       % Number of grid points in ea. direction
end
h = 1/(Ngrid+1);
xx = linspace(0,1,Ngrid+2);  
[x,y] = meshgrid(xx);  % Uniform grid on [0,1]^2  (Spatial units are cm)

%% Set up time grid
Nt = input('How many time points would you like to use?\n(Hit Enter for default (50)) ');
if(isempty(Nt))
    Nt = 50;        % Number of time points to keep 
end
tmax = input('How many days would you like to simulate?\n(Hit Enter for default (3*365 days)) '); 
if(isempty(tmax))
    tmax = 3*365;   % Time (days)
end
t = linspace(0,tmax,Nt);
%%
disp('Generating initial condition and random coefficient fields');
%% Set up initial condition
s  = 0.01;    % Std. dev of initial condition
I0 = 5;        % Integral of intial cell density ("number of initial cells")
initcond = @(x,y,xc,yc) I0*exp(-((x-xc).^2+(y-yc).^2)/(2*s^2))/(2*pi*s^2);
u0 = initcond(x,y,0.5,0.5);
if(plotting)
    thetafig = figure; set(thetafig,'Position',[thetafig.Position(1),thetafig.Position(2),1200,800]);
    subplot(2,2,1);
    imagesc(xx,xx,u0);axis image; set(gca,'YDir','normal'); colorbar;
    set(gca,'Position',[0.05,0.55,0.4,0.4])
    title(sprintf('Initial Condition (I0 = %f)',I0),'FontSize',14);
end
%% Lumpy isotropic diffusion coefficient D(x,y)
L_D = LumpyBgnd('N',Ngrid+2,'Kbar',20,'b',1e-7,'cov',0.04);

if(plotting)
    figure(thetafig);subplot(2,2,2);
    imagesc(xx,xx,L_D.Eval);axis image; set(gca,'YDir','normal'); colorbar;
    set(gca,'Position',[0.55,0.55,0.4,0.4])
    title('Sample diffusion coefficient $D(x,y)$ (units cm$^2/$day)','FontSize',14);
end
%% Lumpy growth function rho(x,y)
% Construct lumpy background object.  Can set parameters here or later w/
% e.g. L.Kbar, L.cov, etc.
L_rho = LumpyBgnd('N',Ngrid+2,'Kbar',200,'b',0.25,'cov',0.002);
% Plot a sample from the growth function
if(plotting)
    figure(thetafig); subplot(2,2,3); 
    imagesc(xx,xx,L_rho.Eval);
    axis image;set(gca,'YDir','normal');colorbar;
    set(gca,'Position',[0.05,0.05,0.4,0.4])
    title('Sample growth rate $\rho(x,y)$ (units cells/day)','FontSize',14); 
end
%% Lumpy carrying cap function
L_kappa = LumpyBgnd('N',Ngrid+2,'Kbar',100,'b',5e7,'cov',0.1); 
K0 = 10;  %NOTE: This isn't the carrying cap, just a scale factor.
if(plotting)
    figure(thetafig); subplot(2,2,4);
    imagesc(xx,xx,L_kappa.Eval);axis image;set(gca,'YDir','normal');colorbar;
    kappaaxis = gca;
    set(gca,'Position',[0.55,0.05,0.4,0.4]);hold on;
    title('Sample carrying capacity function $\kappa(x,y)$','FontSize',14);
    xlabel('$x$'); ylabel('$y$'); zlabel('$\kappa(x,y)$ (Units cells/cm$^2$)','FontSize',14);
end
%% Set up RDE Solver
R = RDE; 

R.rho = L_rho; 
R.kappa = L_kappa;
R.D = L_D;
R.u0 = u0;
R.grid = {xx,xx};
%% Solve
[n,rho,kappa,D] = R.Sample(Nsamp,t);

N = zeros(Nt,Nsamp);
for i=1:Nsamp
    N(:,i) = n{i}.TumorBurden;
end
%% Plot one solution path and the estimated tumor burden over time

n{1}.plot(1:10); 
figure; 
plot(t,N); title('Total tumor burden versus time');
ylabel('\# cells'); 
xlabel('Time (days)'); 
%% Save
fname = input('Enter a file name to save the results, or hit Enter to cancel.\nAccepted formats = .mat, .dat and .npy ','s');
if(~isempty(fname))
    format = fname(end-3:end); 
    if(~strcmp(format,'.mat')&&~strcmp(format,'.dat')&&~strcmp(format,'.npy'))
        error('Incorrect file format!'); 
    end
    if(strcmp(format,'.mat')||strcmp(format,'.dat'))
        save(fname,'n','R','rho','kappa','D','N');
    elseif(strcmp(format,'.npy'))
        % Export to .npy format 
        disp('Saving to .npy file (tumor cell density data only!)')
    end
end

fprintf('Done!\n');