clear;clc
format shortG
% Numerical Methods 1
% Computer Science and Mathematical Engineering Master
% Authors: Javier Blanco, Pablo Castaneda y Jesús Alberto Gallinal
% Matlab Version: R2021b 
% Date: Jan-17-2024
% UNIR - 2024
% Space Ship Data
G = 6.672e-11; % m^3/Kg*s^2 Gravitational Constant
Me = 5.9742e24; % Kg Earth mass
Re = 6378.14; % Km Earth's radius
Vo = 6700; % m/s Velocity
H = 772; % Km Heigh
% Initial Conditions
r_ini = (Re + H) * 1000; % converting to meters
r_dot = 0;
theta_ini = 0;
theta_dot = Vo / r_ini;
% IVP Solution
Y0 = [r_ini; r_dot; theta_ini; theta_dot];
[time, Y_sol] = RK4systems(@space_ship_system, 0, 1200, 100, Y0);
t = table(time, Y_sol(:,1),Y_sol(:,2),Y_sol(:,3),Y_sol(:,4), ...
    'VariableNames',{'Tiempo (s)','r (m)', 'r_dot (m/s)', ...
    'theta (rad)', 'theta_dot (rad/s)'});
disp(t)
figure;
% r vs time
subplot(2, 2, [1, 2]);
plot(time, Y_sol(:,1), 'LineWidth', 2);
title('Posición en función del tiempo');
xlabel('Tiempo (s)');
ylabel('r (m)');
grid on;
% Phase space r vs r_dot
subplot(2, 2, 3);
plot(Y_sol(:,1), Y_sol(:,2), 'LineWidth', 2);
title('Espacio de fase: r vs r\_dot');
xlabel('r (m)');
ylabel('r\_dot (m/s)');
grid on;
% Phase space theta vs theta_dot
subplot(2, 2, 4);
plot(Y_sol(:,3), Y_sol(:,4), 'LineWidth', 2);
title('Espacio de fase: \theta vs \theta\_dot');
xlabel('\theta (rad)');
ylabel('\theta\_dot (rad/s)');
grid on;
sgtitle('Nave Espacial Runge Kutta 4to Orden');
% tightfig;
saveas(gcf, 'plot_output.png')
% Velocity Impact
r_dot_impact = interp1(Y_sol(:,1),Y_sol(:,2),Re*1000,'spline');
fprintf('La velocidad de impacto cuando r = Re por el método de RK-4to Orden es: %.3f m/s\n', r_dot_impact);
% defining Equation System function
function dy = space_ship_system(~, y)
dy = zeros(1,4);
dy(1) = y(2);
dy(2) = y(1) * y(4)^2 - 3.986e14 / y(1)^2;
dy(3) = y(4);
dy(4) = -2 * y(2) * y(4) / y(1);
end
