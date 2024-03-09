% Script que resuelve el Apartado 6 para ode45.
% nota: requiere haber cargado previamente el script "Apartado3ODE45.m"
% obtiene la función "drdt" que interpola el espacio de fases (r, dr/dt)...
% nota: se ha aplicado "flip" porque los valores "y(:, 1)" son decrecientes y "griddedInterpolant" requiere que sean crecientes...
drdt = griddedInterpolant(flip(y(:, 1)), flip(y(:, 2)));
R_e = 6378140; % radio terrestre en metros...
% estimamos la velocidad de impacto cuando r=R_e (se alcanza el suelo)...
drdt_impacto = drdt(R_e); disp(vpa(drdt_impacto));
% represento el espacio de fases de (r, dr/dt) usando "drdt"...
r = linspace(min(y(:, 1)), max(y(:, 1)), 100);
p1 = plot(r, drdt(r), '-'); hold on;
p2 = plot(y(:, 1), y(:, 2), '.'); hold on;
p3 = plot(R_e, drdt(R_e), 'o');
legend('Puntos interpolados', 'Puntos calculados', 'Punto de impacto');
%set(gca, 'XDir', 'reverse');
p1.Color = '#0072BD'; p2.Color = 'blue'; p3.Color = 'red'; grid on;
xlabel('$r$','interpreter', 'latex'); ylabel('$\dot{r}$','interpreter', 'latex');
title(['Espacio de fases ', '$(r,\ \dot{r})$'], 'interpreter', 'latex');
