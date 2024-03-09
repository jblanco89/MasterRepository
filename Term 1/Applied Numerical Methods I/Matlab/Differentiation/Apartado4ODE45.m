% Script que resuelve el Apartado 4 para ode45.
% nota: requiere haber cargado previamente el script "Apartado3ODE45.m"
% represento el espacio de fases de (r, dr/dt)...
subplot(1,2,1); p = plot(y(:, 1), y(:, 2), '-');
%set(gca, 'XDir', 'reverse');
p.Color = 'blue'; grid on;
xlabel('$r$','interpreter', 'latex'); ylabel('$\dot{r}$','interpreter', 'latex');
title(['Espacio de fases ', '$(r,\ \dot{r})$'], 'interpreter', 'latex');
% represento el espacio de fases de (theta, dtheta/dt)...
subplot(1,2,2); p = plot(y(:, 3), y(:, 4), '-');
p.Color = 'red'; grid on;
xlabel('$\theta$','interpreter', 'latex'); ylabel('$\dot{\theta}$','interpreter', 'latex');
title(['Espacio de fases ', '$(\theta,\ \dot{\theta})$'], 'interpreter', 'latex');
