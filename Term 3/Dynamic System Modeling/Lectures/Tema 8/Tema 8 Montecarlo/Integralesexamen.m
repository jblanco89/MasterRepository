clear, clc
format shortG

seed = 142;
rng(seed);


iterations = 100000;
x_1 = 0;
x_2 = @(y) (6-y) / 2;
y_1 = 0;
y_2 = 2;

f = @(x,y) (x-y).^2;
y_int = (y_2 - y_1) .* rand(iterations, 1) + y_1;
x_int = ((x_2(y_int)) - x_1) .* rand(iterations, 1) + x_1;

fv = f(x_int, y_int);

% Volumen MC
Vhp = (mean(x_2(y_int)) - x_1)*(y_2 - y_1);
V_mc = Vhp.*(sum(fv)/iterations);
%V_num = integral2(f, r_lim(1), r_lim(2), theta_lim(1), theta_lim(2));
% err = (abs(V_mc - V_num) / V_num) * 100

%Plot
figure;
hold on;
scatter3(x_int, y_int, fv, 2, 'blue', 'filled')
[X,Y] = meshgrid(linspace(min(x_int), max(x_int),50), linspace(y_1, y_2, 50));
Z = f(X,Y);
surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
colormap(parula); % Profesional colormap
colorbar; % Colorbar
% pcolor(X, Y, Z);
shading interp;
xlabel('x')
ylabel('y')
zlabel('z')
view(45,20)
title('Volumen de la función con Montecarlo');
grid on;