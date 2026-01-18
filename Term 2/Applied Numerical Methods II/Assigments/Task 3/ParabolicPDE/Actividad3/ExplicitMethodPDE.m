function [u, x, t] = ExplicitMethodPDE(alpha, ci, a, b, nx, Tmax, nt)
    % ExplicitMethodPDE Solves a parabolic PDE using the explicit method with non-Dirichlet conditions
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
    %   u - Solution matrix (spatial points x time steps)
    %   x - Spatial grid points
    %   t - Time grid points


    % Compute spatial and temporal step sizes
    h=(b-a)/nx;     
    x=a:h:b;
    k = Tmax / nt;
    t = 0:k:Tmax;

    % Evaluate the initial condition at spatial grid points
    cix = feval(ci, x);
    u = zeros(nx+1, nt);
    u(:, 1) = cix';
    
    % Stability criterion
    lambda = k * alpha^2 / h^2;
    if lambda > 1/2
        disp('No se cumple el criterio de convergencia')
    else
        disp('sin problema')
    end
    
    for j = 2:nt
        u(1, j+1) = (1 - 2*lambda + k*(t(j)^2)) * u(1, j) + lambda * (2*u(2, j) - (2*h*t(j))) + k * x(1) * cos(x(1) * t(j));
        u(2:nx, j+1) = (1 - 2*lambda + k * (t(j)^2)) .* u(2:nx, j) + lambda * (u(3:nx+1, j) + u(1:nx-1, j)) + k * x(2:nx)' .* cos(x(2:nx)' * t(j));
        u(nx+1, j+1) = sin(t(j+1));
    end
end
