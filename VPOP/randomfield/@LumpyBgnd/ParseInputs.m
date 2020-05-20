function p = ParseInputs(~,varargin)
    p = inputParser;
    p.CaseSensitive = true;
    Kbar_default    = 100;
    b_default      = 1;
    B0_default      = 0;
    cov_default     = 0.005;
    centers_default = [];
    N_default       = 128;
    L_default       = [0,1;0,1];
    S_default       = RectSupport([0,1;0,1]);
    isnonneg        = @(x) (x>=0);
    ispos           = @(x) (x>0);
    isbinary        = @(x) (x==0||x==1);
    addParameter(p,'b',b_default,isnonneg);
    addParameter(p,'B0',B0_default,isnonneg);
    addParameter(p,'cov',cov_default,ispos);
    addParameter(p,'Kbar',Kbar_default,ispos);
    addParameter(p,'centers',centers_default,@isnumeric);
    addParameter(p,'N',N_default,ispos);
    addParameter(p,'L',L_default,@isnumeric);
    addParameter(p,'Support',S_default,@(x) 1);
    addParameter(p,'pdf',false,@islogical);
    parse(p,varargin{:});
end