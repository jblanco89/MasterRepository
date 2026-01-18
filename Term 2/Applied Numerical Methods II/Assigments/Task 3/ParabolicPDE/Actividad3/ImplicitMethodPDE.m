function [U, x, t] = ImplicitMethodPDE(alpha, ci, a, b, nx, Tmax, nt)
    % ImplicitMethodPDE Solves a parabolic PDE using the implicit method
    %
    % Inputs:
    %   alpha - Coefficient in the PDE
    %   ci - Initial condition function handle
    %   a - Left boundary of the spatial domain
    %   b - Right boundary of the spatial domain
    %   nx - Number of spatial points
    %   Tmax - Maximum time
    %   nt - Number of time steps
    %
    % Outputs:
    %   U - Solution matrix (spatial points x time steps)
    %   x - Spatial grid points
    %   t - Time grid points

    % Compute spatial and temporal step sizes
    h=(b-a)/nx;     x=a:h:b;
    k=Tmax/nt;      t=0:k:Tmax;
    
    % Stability criterion
    lambda=k*alpha^2/h^2;


    % Initialize solution matrix and set initial condition
    U = zeros(nx+1, nt);
    cix= feval(ci ,x);
    U(:,1)=cix';
    
    
    % Initialize tridiagonal matrix components
    ds = -lambda*ones(nx-1,1);
    di = ds;
    ds(1) = -2*lambda;
    dp = ones(nx,1);
    d = ones(nx,1);
    
    % Time-stepping loop to solve the PDE
    for j = 2:nt+1
        dp(1) = 1+(2*lambda) - k*(t(j)^2);
        dp(nx) = 1+(2*lambda) - k*(t(j)^2);
        dp(2:nx-1) = 1+(2*lambda) - k*(t(j)^2);
        
        d(1) =  k*(x(1) * cos(x(1) * t(j))) - (2 * h * lambda * t(j));
        d(nx) = lambda * sin(t(j)) + k*(x(nx) * cos(x(nx) * t(j)));
        d(2:nx-1) = k.*(x(2:nx-1) .* cos(x(2:nx-1) * t(j)));
        dj = U(1:nx,j-1) + d;
        
        % Solve the tridiagonal system using Crout's method
        z = Crout(dp, ds, di, dj);
        U(1:nx, j) = z;
        U(nx+1, j) = sin(t(j));
    end
    
end
