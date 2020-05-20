function z = minus(x,y)
% Computes the difference of two lumpy backgrounds

suppx = x.S;
suppy = y.S;
Kx    = x.K;
Ky    = y.K;

if(~all(suppx.L(:)==suppy.L(:)))
    error('Incompatible supports!')
end
    
    
z = LumpyBgnd('S',suppx);
z.israndnumlumps = 0;

z.centers = [x.centers;y.centers];
z.K = Kx + Ky;

if(numel(x.b0)==1)
    b0x = x.b0*ones(Kx,1);
elseif(numel(x.b0)==Kx)
    b0x = x.b0(:);
else
    error('Incorrectly formatted b0!');
end

if(numel(y.b0)==1)
    b0y = y.b0*ones(Ky,1);
elseif(numel(y.b0)==Ky)
    b0y = y.b0(:);
else
    error('Incorrectly formatted b0!');
end


z.b0 = [b0x;-b0y];

end