function n = Solve(obj,T)
    % Solves the PDE on the grid for the times specified
    % Returns an object of type RDESolutionPath

    % Set up the mesh grid
    x = obj.grid{1}; y = obj.grid{2};
    [xx,yy] = meshgrid(x,y);
    nx = obj.n_grid(1); ny = obj.n_grid(2); 
    hx = obj.h_grid(1); hy = obj.h_grid(2);
    n_time = length(T);

    % Check to make sure the initial condition is correctly specified 
    if(ismethod(obj.u0,'Eval'))
        u0_grid = obj.u0.Eval([xx(:),yy(:)],[nx,ny]);
    elseif(size(obj.u0,1)==nx&&size(obj.u0,2)==ny)
        u0_grid = obj.u0; 
    else 
        error('The intial condition u0 is incorrectly specified; must be either a class with an Eval method, or a matrix of dimension [nx,ny] = [%i,%i]\n',nx,ny); 
    end 
        
    % Compute fields on the grid OR use fixed parameter matrices if user
    % has specified them. 
    if(ismethod(obj.rho,'Eval'))
        rho_grid = obj.rho.Eval([xx(:),yy(:)],[nx,ny]);
    elseif(size(obj.rho,1)==nx&&size(obj.rho,2)==ny)
        rho_grid = obj.rho; 
    else 
        error('Rho is incorrectly specified; must be either a class with an Eval method, or a matrix of dimension [nx,ny] = [%i,%i]\n',nx,ny); 
    end
    if(ismethod(obj.kappa,'Eval'))
        kappa_grid = obj.kappa.Eval([xx(:),yy(:)],[nx,ny]);
    elseif(size(obj.kappa,1)==nx&&size(obj.kappa,2)==ny)
        kappa_grid = obj.kappa; 
    else 
        error('Kappa is incorrectly specified; must be either a class with an Eval method, or a matrix of dimension [nx,ny] = [%i,%i]\n',nx,ny);
    end
    % NOTE: Dx and Dy must be evaluated on a shifted grid.  See
    % documentation. 
    if(ismethod(obj.D,'Eval'))
        Dx = obj.D.Eval([xx(:)+hx/2,yy(:)],[nx,ny]);  Dx = Dx(1:end-1,2:end-1);
        Dy = obj.D.Eval([xx(:),yy(:)+hy/2],[nx,ny]);  Dy = Dy(2:end-1,1:end-1);
        Dcell{1} = Dx;
        Dcell{2} = Dy;
        disp(size(Dx))
        disp(size(Dy))
    elseif(size(obj.D,1)==nx&&size(obj.D,2)==ny)
        Dcell{1} = obj.D(1:end-1,2:end-1); 
        Dcell{2} = obj.D(2:end-1,1:end-1); 
        if(numel(unique(obj.D))~=1)
            warning('Using user-specified diffusion coefficient matrix D is non-constant. Be aware that the finite difference code assumes a shifted grid for D. See documentation.');
        end
    else 
        error('D is incorrectly specified; must be either a class with an Eval method, or a matrix of dimension [nx,ny] = [%i,%i]\n',nx,ny);
    end
        
    % Set up solution

    n = RDESolutionPath;
    n.grid = obj.grid;
    n.times = T;
    n.cell_density = zeros(nx,ny,length(T));

    % Set up the ODE right-hand-side
    options = odeset; 
    odeset(options,'AbsTol',1e-10,'RelTol',1e-6); 

    if(obj.dim==2)
        odefun = @(T,Y)obj.ODEFun2D(Y,Dcell,rho_grid,kappa_grid,hx,hy,nx,ny);
    elseif(obj.dim==3)
        odefun = @(T,Y)obj.ODEFun3D(Y,Dcell,rho_grid,kappa_grid,hx,hy,nx,ny);
    else
        error('Something is wrong; RDE.dim should only be 2 or 3!'); 
    end

    % Solve the ODE system
    fprintf('Solving RDE...'); tic;
    [~,Y] = ode23(odefun,T,u0_grid(:),options);  
    fprintf('done! (%1.2fs)\n',toc);
    if(n_time==2)
            n.cell_density(:,:,1) = reshape(Y(1,:),[nx,ny]);
            n.cell_density(:,:,2) = reshape(Y(end,:),[nx,ny]);
    else
        for i=1:n_time
            n.cell_density(:,:,i) = reshape(Y(i,:),[nx,ny]); 
        end
    end
end