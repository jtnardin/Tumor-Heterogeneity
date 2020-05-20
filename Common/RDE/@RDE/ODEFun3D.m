function Y_o = ODEFun3D(Y_i,D,alpha,K,h,n)

    % First put into matrix
    Y_i = reshape(Y_i,[n,n,n]);
    Y_o = Y_i;

    % Set Dirichlet bdy conditions
    Y_o(1,:,:) = 0;
    Y_o(n,:,:) = 0;
    Y_o(:,1,:) = 0;
    Y_o(:,n,:) = 0;
    Y_o(:,:,1) = 0;
    Y_o(:,:,n) = 0;

    Y_o(2:n-1,2:n-1,2:n-1) = (D/h^2)*(Y_i(3:end,2:end-1,2:end-1) + Y_i(1:end-2,2:end-1,2:end-1) + ...
                                      Y_i(2:end-1,3:end,2:end-1) + Y_i(2:end-1,1:end-2,2:end-1) + ...
                                      Y_i(2:end-1,2:end-1,3:end) + Y_i(2:end-1,2:end-1,1:end-2) - ...                       
                                      6*Y_i(2:n-1,2:n-1)) +...
                                alpha(2:end-1,2:end-1,2:end-1).*Y_i(2:end-1,2:end-1,2:end-1).*(1-(Y_i(2:end-1,2:end-1,2:end-1)./K(2:end-1,2:end-1,2:end-1)));
    % Unroll
    Y_o = Y_o(:);
end