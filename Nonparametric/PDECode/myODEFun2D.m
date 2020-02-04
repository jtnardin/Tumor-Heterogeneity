function Y_o = myODEFun2D(Y_i,D,Aa,Kk,h,n)
% First put into matrix
Y_i = reshape(Y_i,[n,n]);
Y_o = Y_i;

% Set bdy conditions
Y_o(1,:) = 0;
Y_o(n,:) = 0;
Y_o(:,1) = 0;
Y_o(:,n) = 0;

%Aa is a scalaer
if numel(Aa) == 1

    Y_o(2:n-1,2:n-1) = (D/h^2)*(Y_i(3:end,2:end-1) + Y_i(1:end-2,2:end-1) + ...
                                Y_i(2:end-1,3:end) + Y_i(2:end-1,1:end-2) -...
                                4*Y_i(2:n-1,2:n-1)) +...
                                Aa*Y_i(2:end-1,2:end-1).*(1-(Y_i(2:end-1,2:end-1)/Kk));

else %Aa is a matrix

    Y_o(2:n-1,2:n-1) = (D/h^2)*(Y_i(3:end,2:end-1) + Y_i(1:end-2,2:end-1) + ...
                                Y_i(2:end-1,3:end) + Y_i(2:end-1,1:end-2) -...
                                4*Y_i(2:n-1,2:n-1)) +...
                                Aa(2:end-1,2:end-1).*Y_i(2:end-1,2:end-1).*...
                                (1-(Y_i(2:end-1,2:end-1)./Kk(2:end-1,2:end-1)));
end

    
% Unroll
Y_o = Y_o(:);

end