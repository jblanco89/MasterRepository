function x = Jacobi(A,b,x0,Tol,maxiter)
% JACOBI Solves a linear system Ax = b using the Jacobi iterative method.
% x = Jacobi(A, b, x0, Tol, maxiter)
%   Inputs:
%       A:      Coefficient matrix of the system.
%       b:      Right-hand side vector of the system.
%       x0:     Initial guess for the solution.
%       Tol:    Tolerance for convergence (stop criterion).
%       maxiter: Maximum number of iterations allowed.
%
%   Output:
%       x:      Approximate solution vector.
%
%   Example:
%       A = [4 -1 0; -1 4 -1; 0 -1 4];
%       b = [15; 10; 10];
%       x0 = [0; 0; 0];
%       Tol = 1e-6;
%       maxiter = 1000;
%       x = Jacobi(A, b, x0, Tol, maxiter);
%
%   Note:
%       The function prints the iteration number, the updated solution
%       vector, and the norm of the difference between consecutive
%       solutions in each iteration.

K=0;
n = length(x0);
Norm = 1;
x = zeros(1, n);
while Norm > Tol
    K = K + 1;
    fprintf('%2d', K);
    for i=1:n
        sum = 0;
        for j=1:n
            if i ~= j
                sum = sum + A(i,j)*x0(j);
            end
        end
        x(i) = (b(i) - sum) / A(i, i);
    fprintf('%10.6f', x(i));
    end
    Norm = norm(x0-x, inf);
    fprintf('%10.6f\n', Norm);
    x0 = x;
    if K > maxiter
        disp('Convergencia no alcanzada con el máximo número de iteraciones')
        break
    end
end
end
