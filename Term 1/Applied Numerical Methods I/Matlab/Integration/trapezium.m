function [I, err] = trapezium(f, y, a, b, n)
% Trapezoidal method for numerical integration
% Inputs:

    % f: function handle; example: f=@(x) sin (x).* exp (-x);
    % y: function image from data table, example: y = [f(x1), f(x2),
    % ...f(xn)]
    % a, b: integration limits
    % n: number of subintervals or values of table

% Return:
%   I: Result of numerical integration using the trapezoidal method
% err: Error result between analytical and numeric calculation. If function 
% is unknown, err = - (h^2 / 12) * (b-a) * f''(xi) (see Lesson 5 PDF, pg. 9)

%   Example:
%       f = @(x) sin(x) .* exp(-x);
%       a = 0;
%       b = pi;
%       n = 100;
%       result = trapezium(f, a, b, n);

h = (b-a)/ n;
x = a:h:b;
if  y == 0
    y = f(x);
    w = [1 2*ones(1, n-1) 1];
    I = h/2 * sum(w .* y);
    syms x
    I_actual = double(int(f(x), a, b));
    % err = abs((I - I_actual) / I_actual)*100; % porcentual relative error
    err = abs(I - I_actual); % absolute error
else 
    I = h/2 * ((y(1) + y(end)) + 2*(sum(y(2:end-1))));
    second_derivative = diff(y, 2);
    err = - (h^2 * (b-a) * max(second_derivative)) / 12;
end
end