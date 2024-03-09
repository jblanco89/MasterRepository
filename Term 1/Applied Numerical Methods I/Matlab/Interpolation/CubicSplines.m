function [ai, bi, ci, di, p] = CubicSplines(xi, fi)
% CubicSplines - Computes natural cubic spline coefficients and interpolates points.
%
% Usage:
%   [ai, bi, ci, di, p] = CubicSplines(xi, fi)
%
% Inputs:
%   xi: Vector of x-coordinates of the data points.
%   fi: Vector of corresponding function values at the data points.
%
% Outputs:
%   ai: Coefficients representing the function values at data points.
%   bi: Coefficients for linear terms in the cubic spline equations.
%   ci: Coefficients for quadratic terms in the cubic spline equations.
%   di: Coefficients for cubic terms in the cubic spline equations.
%   p: Array of symbolic cubic spline polynomials for each segment.
%
% Description:
%   This function takes vectors of x-coordinates (xi) and corresponding
%   function values (fi) and computes the coefficients for a natural cubic
%   spline interpolation. It returns the coefficients ai, bi, ci, and di,
%   which define the cubic spline equations. Additionally, it provides an
%   array of symbolic cubic spline polynomials (p) for each segment.
%
%   The cubic spline is represented by the equation:
%   S_i(x) = ai + bi * (x - xi) + ci * (x - xi)^2 + di * (x - xi)^3
%   for each segment i, where xi is the starting point of the segment.
%
% Example:
%   % Define data points
%   xi = [1971, 1981, 1991, 2001, 2011];
%   fi = [33.956, 37.743, 39.434, 40.847, 46.816];
%
%   % Compute cubic spline coefficients
%   [ai, bi, ci, di, p] = CubicSplines(xi, fi);
%
%   % Evaluate the cubic spline at a specific point (e.g., x = 2005)
%   x_value = 2005;
%   interpolated_value = subs(p(4), x, x_value);
%
% References:
%   - De Boor, Carl. (1978). A Practical Guide to Splines.
%     Springer-Verlag.
%   - Burden, R. L., & Faires, J. D. (2010). Numerical Analysis.
%     Brooks/Cole.
%
% Author:
%   Msc. Javier Blanco
%
% Creation Date:
%   Feb-24-2024

n = length(fi);
ai = fi';
bi = zeros(n,1);
di = zeros(n,1);
b = zeros(n,1);
A = zeros(n,n);
h = xi(2) - xi(1);
%check this in case points are not equispaced. 
% for i = 2:n-1
%     h(i) = xi(i) - xi(i-1);
% end
A(1,1) = 1; A(n,n) = 1;
b(1) = 0; b(n) = 0;
for i=2:n-1
    A(i, i-1) = h;
    A(i, i) = 2*(h + h);
    A(i, i+1) = h;
    b(i) = 3*(((ai(i+1) - ai(i)) / h) - ((ai(i) - ai(i-1))/ h));
end
 
dP = diag(A);
[L, U] = lu(A);
dI = L;
dS = U;
ci = Crout(dP, dS, dI, b);
% ci = A\b;
ci(1) = 0; ci(n) = 0;
for i = 1:n-1
    bi(i) = ((ai(i+1) - ai(i))/h) - (h*(2*ci(i)+ci(i+1))/3);
    di(i) = (ci(i+1) - ci(i)) / (3*h);
    syms x;
    p(i) = ai(i) + bi(i)*(x-xi(i)) + ci(i)*(x - xi(i))^2 + di(i)*(x-xi(i))^3;

end
pretty(p)
end