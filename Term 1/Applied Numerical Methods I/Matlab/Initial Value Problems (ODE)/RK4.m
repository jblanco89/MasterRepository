function [t, y] = RK4(f, a, b, N, ya)
% 4th Order Runge-Kutta for Ordinary Differential Equation.
% Inputs:
%   - f: Represents the ordinary differential equation.
%   - a: Initial time.
%   - b: Final time.
%   - N: Number of time steps.
%   - ya: Initial conditions of the system.
% Outputs:
%   - t: Time vector ranging from a to b with n+1 points.
%   - y: Vector containing the values of the state variables at each time step.
% Example:
    % % Call the RK4systems function to solve the system
    % [t, y] = RK4(f, a, b, N, Ya);

h = (b - a) / N; 
t = a:h:b; 
t = t(:); 
y = zeros(N+1, 1); 
y(1) = ya;
for i  = 1:N
    k1 = feval(f,t(i), y(i));
    k2 = feval(f,t(i) + h/2, y(i) + h/2*k1);
    k3 = feval(f,t(i) + h/2, y(i) + h/2*k1);
    k4 = feval(f,t(i) + h, y(i) + h*k3);
    y(i+1) = y(i) + h / 6 * (k1 + 2*k2 + 2*k3 + k4);
end
end