function [err,best_approx]=errorfunc_discrete_sparse(r,totsol,meansol,weights)
% solve differential equation
% calculate error
% maxdim=length(r)/2;                                  % dimensions of what we estimate
% %keyboard;
% 
% weightsol=zeros(size(totsol));
% %keyboard;
% for k=1:length(nodesr)
%     for j=maxdim+1:maxdim+length(nodesK)
%         weightsol(k,j-maxdim,:,:)=r(k)*r(j)*squeeze(totsol(k,j-maxdim,:,:));
%     end
% end
% best_approx=squeeze(sum(sum(weightsol)));                  % this is our weighted average
[a,b,c]=size(meansol);

best_approx=reshape(r'*totsol,[a,b,c]);
err=sum(sum(sum(weights.*((meansol-best_approx).^2))));
