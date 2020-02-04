function lj=make_splines(nodesj,ts)
%
% MAKE_SPLINES calculates the spline functions for given node set of nodes
% (nodesj) and quadrature nodes (ts).
%
%    ARGUMENTS:
%        nodesj: A set of 3 numbers which are the nodes you are trying to
%            spline over
%        ts: The full quadrature node space 
%
%    RETURNS:
%        lj: The spline approximation for the points nodesj on the space ts
%
% Written by Erica Rutter (July 2017)

lj=zeros(length(ts),1);

for k=1:length(ts)
    t=ts(k);
    if (t>=nodesj(1) && t<nodesj(2))
        lj(k) = (t-nodesj(1))./(nodesj(2)-nodesj(1));
    elseif(t>=nodesj(2) && t<nodesj(3))
        lj(k) = (nodesj(3)-t)./(nodesj(3)-nodesj(2));
    else
        lj(k)=0;
    end
end




