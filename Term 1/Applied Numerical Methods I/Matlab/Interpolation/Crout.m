function x=Crout(dP,dS,dI,b)
% Crout - Solves a linear system Ax = b using the Crout algorithm.
%
% Usage:
%   x = Crout(dP, dS, dI, b)
%
% Inputs:
%   dP: Main diagonal of matrix A
%   dS: Upper diagonal of matrix A
%   dI: Lower diagonal of matrix A
%   b: Right-hand side vector of the system Ax = b
%
% Output:
%   x: Solution of the linear system Ax = b
%
% Description:
%   This function uses the Crout algorithm to solve a linear system Ax = b,
%   where A is a tridiagonal matrix defined by its main (dP), upper (dS), and
%   lower (dI) diagonals. The vector b contains the coefficients of the
%   right-hand side of the system. The function returns the solution x of the
%   system.
%
% Example:
%   % Define a tridiagonal matrix
%   dP = [2, 3, 4];
%   dS = [1, 2, 3];
%   dI = [1, 2, 3];
%   
%   % Define the right-hand side vector
%   b = [1; 2; 3];
%   
%   % Solve the linear system using Crout
%   x = Crout(dP, dS, dI, b);
%
% References:
%   - Golub, G. H., and Van Loan, C. F. (1996). Matrix Computations.
%     The Johns Hopkins University Press.
%   - Saad, Y. (2003). Iterative Methods for Sparse Linear Systems.
%     Society for Industrial and Applied Mathematics.
%
% Author:
%   Dr. Paular Triguero
%   Msc. Javier Blanco
%
% Creation Date:
%   Feb-24-2024

n=length(dP);
% 1. Obtención de las matrices L y U tales que A = LU

l(1)=dP(1);
u(1)=dS(1)/l(1);

for i=2:n-1
    l(i)=dP(i)-dI(i-1)*u(i-1);
    u(i)=dS(i)/l(i);
end

l(n)=dP(n)-dI(n-1)*u(n-1);
% 2. Solución del sistema Lz = b
z(1) =b(1) /l(1) ;
for i =2:n
    z(i)=(1/l(i))*(b(i)-dI(i-1)*z(i-1));
end
% 3. Solución del sistema Ux = z
x(n)=z(n);
for i=n-1:-1:1
    x(i)=z(i)-u(i)*x(i+1);
end
x=x(:);
