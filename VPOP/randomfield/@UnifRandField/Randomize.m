

function Randomize(obj,shuff)
    % #mosttrivialmethodever
    if(nargin>0)
        rng('shuffle');
    end
    obj.c = obj.dist();   
end