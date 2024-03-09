function [I, err] = simpson13(f, y, a, b, n)
% Simpson's method for numerical integration
% Inputs:

    % f: function handle; example: f=@(x) sin (x).* exp (-x);
    % y: function image from data table, example: y = [f(x1), f(x2),
    % ...f(xn)]
    % a, b: integration limits
    % n: number of subintervals or values of table

% Notice:

%   n sub-intervals in Simpson's method must be even
%   set f = 0 if function handle is unknown and you only have a data table

% Return:
%   I: Result of numerical integration using the Simpson's method
% err: Error result between analytical and numeric calculation. If function 
% is unknown, err = - (h^4 / 180) * (b-a) * f''''(xi) (see Lesson 5 PPT, pg. 24)

%   Example:
%       f = @(x) sin(x) .* exp(-x);
%       a = 0;
%       b = pi;
%       n = 100;
%       result = simpson13(f, a, b, n);

% Check if the number of subintervals is even
if mod(n, 2) ~= 0
    error('Number of subintervals (n) must be even.');
end
h = (b-a)/ (n-1);
x = a:h:b;
w = ones(1, n);
w(2:2:end-1) = 4;
w(3:2:end-2) = 2;
if  y == 0
    y = f(x);
    I = h/3 * sum(w .* y);
    syms x
    I_actual = double(int(f(x), a, b));
    % err = abs((I - I_actual) / I_actual)*100; % porcentual relative error
    err = abs(I - I_actual); % absolute error
else 
%     I = h/3 * ((y(1) + y(end)) + 2*(sum(y(2:end-1))));
    I = h/3 * sum(w .* y);
    disp(w')
    fourth_derivative = diff(y, 4);
    err = - (h^4 * (b-a) * max(fourth_derivative)) / 180;
end
end