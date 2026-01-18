function [t,y] = Euler(f, a, b, N, ya)
%EULER Solves an ordinary differential equation (ODE) using the Euler method.
%   [t, y] = Euler(f, a, b, N, ya) solves the ODE y' = f(t, y) over the interval t=[a, b]
%   with N steps using the Euler method. The initial condition is y(a) = ya.
%
%   Input:
%       f: Function defining the ODE y' = f(t, y).
%       a: Lower limit of the interval.
%       b: Upper limit of the interval.
%       N: Number of discretization steps.
%       ya: Initial value y(a).
%
%   Output:
%       t: Vector of discrete time points.
%       y: Vector of approximate solutions corresponding to the times in t.
%
%   Example:
%       % Solve the ODE y' = -y with y(0) = 1 over the interval [0, 1] with 10 steps.
%       f = @(t, y) -y;
%       [t, y] = Euler(f, 0, 1, 10, 1);
%       plot(t, y);
%       xlabel('Time (t)');
%       ylabel('Solution (y)');
%       title('Solution of the ODE y'' = -y using the Euler method');
%
h = (b - a)/ N;
t = a:h:b;
t = t(:);
y = zeros(N+1, 1);
y(1) = ya;
for k = 1:N
    y(k+1) = y(k) + h*feval(f, t(k), y(k));
end
end