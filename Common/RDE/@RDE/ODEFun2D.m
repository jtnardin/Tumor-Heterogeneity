% -------------------------------------------------------------------------
% ************************* RDEVPop Package *******************************
% File:     ODEFun2D.m
% Author:   Nick Henscheid
% Date:     10-2019, 2-2020
% Info:     Computes the ODE RHS for the Fisher-KPP equation, assuming
%           method-of-lines semidiscrete scheme with 
%           conservative centered spatial differencing w/ a 5-point stencil 
%           (see Hundsdorfer section 6.2).
%           Assumptions: 
%           1. Y_i and Y_o have size (Nx+2)*(Ny+2) where Nx and Ny are the 
%              number of *interior* solution nodes.
%           2. D is a 2-by-1 cell array; the first cell is the 
%              x-shifted coefficient, the second cell is the y-shifted 
%              coefficient. D{1} is (Nx+1)-by-(Ny).  D{2} is (Nx)-by-(Ny+1) 
%           
% Inputs:         
%               
% Contact: nph@email.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

function Y_o = ODEFun2D(Y_i,D,Aa,Kk,hx,hy,nx,ny)

% First put the input into a matrix
Y_i = reshape(Y_i,[nx,ny]);
Y_o = Y_i;

% Set Dirichlet bdy conditions 
Y_o(1,:) = 0;
Y_o(nx,:) = 0;
Y_o(:,1) = 0;
Y_o(:,ny) = 0;

Dx = D{1}; 
Dy = D{2}; 
% Compute diffusion operator on the interior nodes
Y_o(2:nx-1,2:ny-1) = (1/hx^2)*(Dx(2:end,:).*(Y_i(3:end,2:end-1) - Y_i(2:end-1,2:end-1)) -...
                            Dx(1:end-1,:).*(Y_i(2:end-1,2:end-1) - Y_i(1:end-2,2:end-1))) +...
                     (1/hy^2)*(Dy(:,2:end).*(Y_i(2:end-1,3:end) - Y_i(2:end-1,2:end-1)) -...
                            Dy(:,1:end-1).*(Y_i(2:end-1,2:end-1) - Y_i(2:end-1,1:end-2)));    
 
% Add reaction term 
Y_o(2:nx-1,2:ny-1) = Y_o(2:nx-1,2:ny-1) + Aa(2:end-1,2:end-1).*Y_i(2:end-1,2:end-1).*(1-(Y_i(2:end-1,2:end-1)./Kk(2:end-1,2:end-1)));

% Unroll
Y_o = Y_o(:);

end