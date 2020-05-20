function u = Eval(obj,X,XSize)
    % This function will evaluate the texture for the sample points
    % defined in the array X.
    % This should be updated once the GridData object is ready to
    % go!
    L = obj.L;   % Bounding box
    N = obj.N;
    if(numel(N)==1)
        N = N*ones(1,obj.dim);
    end
    if(nargin<2) % No array provided
        if(obj.dim==2)
            xtemp = linspace(L(1,1),L(1,2),N(1));
            ytemp = linspace(L(2,1),L(2,2),N(2));
            [xx,yy] = meshgrid(xtemp,ytemp);
            X = [xx(:),yy(:)];
            XSize = N;
        elseif(obj.dim==3)
            xtemp = linspace(L(1,1),L(1,2),N(1));
            ytemp = linspace(L(2,1),L(2,2),N(2));
            ztemp = linspace(L(3,1),L(3,2),N(3));
            [xx,yy,zz] = meshgrid(xtemp,ytemp,ztemp);
            X = [xx(:),yy(:),zz(:)];
            XSize = N;
        end
    end
    %X = X';
    if(numel(obj.b)==1)
        b_vec = obj.b*ones(1,obj.K);
    elseif(numel(obj.b)==obj.K)
        b_vec = reshape(obj.b,[1,obj.K]);
    else
        error('Incorrectly formatted b! Must be 1-by-K or K-by-1');
    end
    
    % Process cov into the vectorized form 
    if(all(size(obj.cov)==[obj.dim,obj.dim]))
        % Uniform covariance
        cov = repmat(obj.cov,[1,obj.K]); 
    elseif(size(obj.cov,2)==obj.dim*obj.K)
        % Non-uniform covariance
        cov = obj.cov; 
    else
        error('Incorrectly formatted cov! Must be dim-by-dim or dim-by-K*dim.  Size(cov) = (%i,%i)\n dim = %i\n K = %i',size(obj.cov,1),size(obj.cov,2),obj.dim,obj.K);
    end
    
    % Note: the lump function is an un-normalized Gaussian PDF, i.e. 
    %            l(x) = A*exp(-0.5(x-x0)*inv(S)*(x-x0)')
    lumpfun = @(xx,A,mu,S) reshape(A*((2*pi)^(obj.dim/2))*sqrt(det(S))*mvnpdf(xx,mu,S),XSize);
    
    u = zeros(XSize); 
    
    % Compute the random field.  Note that this is the slowest step; a mex
    % implementation is available.  
    for i=1:obj.K
        S_i  = cov(:,(obj.dim*(i-1)+1):obj.dim*i); 
        A_i  = b_vec(i); 
        mu_i = obj.centers(i,:); 
        
        u = u + lumpfun(X,A_i,mu_i,S_i); 
    end  
end
