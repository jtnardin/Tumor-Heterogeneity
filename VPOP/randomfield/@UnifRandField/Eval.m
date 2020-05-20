function u = Eval(obj,X,XSize)
    % 
    N = obj.N;
    if(numel(N)==1)
        N = N*ones(1,obj.dim);
    end
    
    if(nargin==2) % Array size provided
        N = XSize;
    end
    
    u = obj.c*ones(N);
end


 