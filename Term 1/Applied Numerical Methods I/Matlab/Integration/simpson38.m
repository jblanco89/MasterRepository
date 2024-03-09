function [I, err] = simpson38(f, y, a, b, n)
% Simpson's method for numerical integration
% Inputs:

    % f: function handle; example: f=@(x) sin (x).* exp (-x);
    % y: function image from data table, example: y = [f(x1), f(x2),
    % ...f(xn)]
    % a, b: integration limits
    % n: number of subintervals or values of table

% Notice:

%   n sub-intervals in Simpson's method must be multiple of 3
%   set f = 0 if function handle is unknown and you only have a data table

% Return:
%   I: Result of numerical integration using the Simpson's method
% err: Error result between analytical and numeric calculation. If function 
% is unknown, err = (n * h^5 / 80) * f''''(xi) (see: https://es.wikipedia.org/wiki/Regla_de_Simpson)

%   Example:
%       f = @(x) sin(x) .* exp(-x);
%       a = 0;
%       b = pi;
%       n = 100;
%       result = simpson38(f, a, b, n);

% Check if the number of subintervals is even
if mod(n, 3) ~= 0
    error('Number of subintervals (n) must be multiple of 3.');
end
h = (b-a)/ n;
x = a:h:b;
w = ones(1, n-1) * 3;
w([1, end+2]) = 1;
w(4:3:end+1) = 2;
if  y == 0
    y = f(x);
    I = (3/8) * h * (sum(w .* y));
    syms x
    I_actual = double(int(f(x), a, b));
    % err = abs((I - I_actual) / I_actual)*100; % porcentual relative error
    err = abs(I - I_actual); % absolute error
else
    I = (3/8) * h * (sum(w .* y));
    fourth_derivative = diff(y, 4);
    err = (n * h^5 * max(fourth_derivative)) / 80;
end
end