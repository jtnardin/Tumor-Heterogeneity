%myODEFun2D.m written 12-6-2019 to provide the RHS for the 2d Fisher-KPP 
%equation for "A tutorial Review of Mathematical Techniques for Quantifying 
%Tumor Heterogeneity" by Everett et al. 

function Y_o = myODEFun2D(Y_i,D,Aa,Kk,h,n)

% First put into matrix
Y_i = reshape(Y_i,[n,n]);
Y_o = Y_i;

% Set bdy conditions
Y_o(1,:) = 0;
Y_o(n,:) = 0;
Y_o(:,1) = 0;
Y_o(:,n) = 0;


Y_o(2:n-1,2:n-1) = (D/h^2)*(Y_i(3:end,2:end-1) + Y_i(1:end-2,2:end-1) + ...
                            Y_i(2:end-1,3:end) + Y_i(2:end-1,1:end-2) -...
                            4*Y_i(2:n-1,2:n-1)) +...
                            Aa*Y_i(2:end-1,2:end-1).*(1-(Y_i(2:end-1,2:end-1)/Kk));

% Unroll
Y_o = Y_o(:);

end