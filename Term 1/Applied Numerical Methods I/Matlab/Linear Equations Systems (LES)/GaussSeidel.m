function x =GaussSeidel(A,b,x0,Tol,maxiter)
% GAUSSSEIDEL Solves a linear system Ax = b using the Gauss-Seidel iterative method.
%   x = GaussSeidel(A, b, x0, Tol, maxiter)
%
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
%       x = GaussSeidel(A, b, x0, Tol, maxiter);
%
%   Note:
%       The function checks the convergence condition for the Gauss-Seidel
%       method by examining the eigenvalues of the iterative matrix. If the
%       largest eigenvalue modulus is not less than 1, the method may not
%       converge. The iteration process prints the iteration number, the
%       updated solution vector, and the norm of the difference between
%       consecutive solutions in each iteration.
%
k=0;
Norm = 1;
D = diag(diag(A));
L = tril(A)- D;
U = triu(A)- D;
% check for convergence condition for Gauss-Seidel method
e=max(eig(-inv(D+L)*(U)));
if abs(e) >= 1
    disp ('Dado que el módulo del valor propio más grande de la matriz iterativa no es menor que 1.') 
    disp ('El proceso no será convergente')
    return
end
x(:,1) = x0;
while Norm > Tol && k < maxiter
    k = k + 1;
    x(:,k+1) = -inv(D+L)*(U)*x(:,k) + (D+L)\b';% Gauss-Seidel formula
    Norm = norm(x(:,k+1) - x(:,k), inf);
    fprintf('%2d', k);
    fprintf('%10.6f', x(:,k)');
    fprintf('%10.6f\n', Norm);
    if k > maxiter
        disp('Convergencia no alcanzada con el máximo número de iteraciones')
       break
    end
end
x = x(:,k)';
end