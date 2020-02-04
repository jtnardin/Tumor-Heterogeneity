x = linspace(0,1,10);
y = linspace(0,1,10);
t = linspace(0,1000,10);

[X,Y,T] = meshgrid(x,y,t);

X_all = [T(:) X(:) Y(:)];

Y = fisher_2d_sim(.02,X_all);

Y = reshape(Y,length(t),length(x),length(y));

for i = 1:10
    subplot(2,5,i)
    surf(Y(:,:,i),'edgecolor','none'); view(2)
end


