function varargout = plot(obj,alpha)
    plot(obj.Support,0.25)
    n = obj.N;

    if(obj.dim == 2)
        X = obj.Support.unifGrid([n,n]);
        xplot = X(1,:,1);
        imagesc(xplot,xplot,obj.Eval);set(gca,'YDir','normal');axis image;
    elseif(obj.dim == 3)
        fprintf('Plotting 3D Isosurface...this may take a moment...\n');
        X = obj.Support.unifGrid([n,n,n]);
        Xeval = [reshape(X(:,:,:,1),[n^3,1]),reshape(X(:,:,:,2),[n^3,1]),reshape(X(:,:,:,3),[n^3,1])];
        u  = obj.Eval(Xeval,[n,n,n]);
        if(nargin<2)
            alpha = mean(u(:));
        end
        s = isosurface(X(:,:,:,1),X(:,:,:,2),X(:,:,:,3),u,alpha);
        p = patch(s);
        p.FaceColor = [1,0,0];
        p.EdgeAlpha = 0.2;
    end
    if(nargout==1)
        varargout{1} = gcf;
    elseif(nargout==2)
        varargout{1} = gcf;
        varargout{2} = p;
    end
    
end