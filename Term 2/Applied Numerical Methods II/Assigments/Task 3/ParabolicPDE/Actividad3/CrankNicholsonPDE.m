function [U,x, t]=CrankNicholsonPDE(alpha,ci,a,b,nx,Tmax,nt)
    % CrankNicholsonPDE Solves a parabolic PDE using the Crank-Nicholson method
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
    h=(b-a)/nx;
    x=a:h:b;
    k=Tmax/nt;
    t=0:k:Tmax;
    
    
    % Initialize solution matrix and set initial condition
    U=zeros(nx+1,nt);
    cix=feval(ci,x);
    U(:,1)=cix';
    
    % Stability criterion
    lambda=k*alpha^2/h^2;
    
    % Initialize tridiagonal matrix components for the Crank-Nicholson method
    ds=(-lambda/2)*ones(nx-1,1);
    di=ds;
    ds(1)=-lambda;
    
    bs=(lambda/2)*ones(nx-1,1);
    bi=bs;
    bs(1)=lambda;
    d = ones(nx,1);


    % Time-stepping loop to solve the PDE
    for j=1:nt

        dp=(1+lambda)-((k/2)*t(j)^2)*ones(nx,1);
        bp=(1-lambda)+((k/2)*t(j)^2)*ones(nx,1);

        % Construct the B matrix
        B=diag(bp)+diag(bs,1)+diag(bi,-1);

        % Construct the dj vector
        d(1)=(-lambda*h)*(t(j)+(t(j+1))) + (k/2)*(x(1)*(cos(x(1)*t(j)))+cos(x(1)*t(j+1)));
        d(nx)=(lambda/2)*(sin(t(j))+sin((t(j+1)))) + (k/2)*(x(nx)*(cos(x(nx)*t(j)))+cos(x(nx)*t(j+1)));
        d(2:nx-1)=(k/2)*(x(2:nx-1).*(cos(x(2:nx-1)*t(j)))+cos(x(2:nx-1)*t(j+1)));
        dj=(B*U(1:nx,j)) + d;


        % Solve the tridiagonal system using Crout's method
        z=Crout(dp,ds,di,dj);
    
        U(1:nx,j+1)=z;
        U(nx+1,j+1)=sin(t(j));
    end
end
