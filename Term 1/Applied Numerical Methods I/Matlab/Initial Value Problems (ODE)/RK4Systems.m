function [t, Y] = RK4Systems(f, a, b, N, Ya)
% 4th Order Runge-Kutta for Differential Equation Systems.
% Inputs:
%   - function_systems: Represents the system of differential equations. It should have the form:
%                       dydt = function_systems(t, y)
%                       where t is time and y is the vector of state variables.
%   - a: Initial time.
%   - b: Final time.
%   - n: Number of time steps.
%   - Ya: Initial conditions of the system.
% Outputs:
%   - t: Time vector ranging from a to b with n+1 points.
%   - Y_values: Matrix containing the values of the state variables at each time step.
% Example:
    % Define the system of differential equations
    % function dydt = example_system(t, y)
    %     dydt = zeros(2, 1); dydt(1) = -0.1*y(1) + 0.2*y(2); dydt(2) = -0.2*y(1) - 0.1*y(2);
    % end
    % Define the initial conditions and the time interval
    % a = 0; b = 10; n = 100; Ya = [1; 0];
    % % Call the RK4systems function to solve the system
    % [t, Y_values] = RK4systems(@example_system, a, b, n, Ya);

h = (b - a) / N; 
t = a:h:b; 
t = t(:); 
Y = zeros(N+1, length(Ya)); 
Y(1,:) = Ya;
for i  = 1:N
    k1 = feval(f,t(i), Y(i, :))';
    k2 = feval(f,t(i) + h/2, Y(i,:) + h/2*k1)';
    k3 = feval(f,t(i) + h/2, Y(i,:) + h/2*k1)';
    k4 = feval(f,t(i) + h, Y(i,:) + h*k3)';
    Y(i+1,:) = Y(i, :) + h / 6 * (k1 + 2*k2 + 2*k3 + k4);
end
end