function p = ParseInputs(~,varargin)
    p = inputParser;
    p.CaseSensitive = false;
    dist_default    = @()rand(); 
    N_default       = 128;
    S_default       = RectSupport([0,1;0,1]);
    ispos           = @(x) (x>0);
    addParameter(p,'dist',dist_default,@(x)isa(x,'function_handle')); 
    addParameter(p,'N',N_default,ispos);
    addParameter(p,'Support',S_default,@(x) 1);
    parse(p,varargin{:});
end