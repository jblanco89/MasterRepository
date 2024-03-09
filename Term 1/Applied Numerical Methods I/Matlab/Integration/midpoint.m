function [I, err] = midpoint(f, y, a, b, n)
% Middle point method for numerical integration
% Inputs:

    % f: function handle; example: f=@(x) sin (x).* exp (-x);
    % y: function image from data table, example: y = [f(x1), f(x2),
    % ...f(xn)]
    % a, b: integration limits
    % n: number of subintervals or values of table

% Return:
%   I: Result of numerical integration using the middle point method
% err: Error result between analytical and numeric calculation. If function 
% is unknown, err = ((b-a)^3 / (24 * n^2)) * f''(xi) (see Theorem 3.5 https://openstax.org/books/c%C3%A1lculo-volumen-2/pages/3-6-integracion-numerica#fs-id1165042234318)

%   Example:
%       f = @(x) sin(x) .* exp(-x);
%       a = 0;
%       b = pi;
%       n = 100;
%       result = middlepoint(f, a, b, n);

h = (b-a)/ n;
x = a:h:b;
w = ones(1, n+1);
if  y == 0
    y = f(x);
    I = h * sum(w .* y);
    syms x
    I_actual = double(int(f(x), a, b));
    % err = abs((I - I_actual) / I_actual)*100; % porcentual relative error
    err = abs(I - I_actual); % absolute error
else 
    I = h .* sum(w .* y);
    second_derivative = diff(y, 2);
    err = ((b-a)^3 / (24 * n^2)) * max(second_derivative);
end
end