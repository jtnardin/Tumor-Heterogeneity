function [n,rho,kappa,D] = Sample(obj,N,T)
    % Generate N sample paths
    % Each element of n is an RDESolutionPath object
    n = cell(N,1);
    rho = cell(N,1);
    kappa = cell(N,1);
    D = cell(N,1);

    for i=1:N
        fprintf('Computing sample %i/%i\n',i,N);
        obj.Randomize;
        % Need to use the "copy" method since rand field models are
        % handle class
        rho{i} = obj.rho.Copy;
        kappa{i} = obj.kappa.Copy; 
        D{i} = obj.D.Copy;
        n{i} = obj.Solve(T);
    end
end