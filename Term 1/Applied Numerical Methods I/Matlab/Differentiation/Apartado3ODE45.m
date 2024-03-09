% Script que resuelve el Apartado 3 con ode45.
% establezco formato de salida...
format long;
% carga datos...
tspan = [0, 1200]; % intervalo de tiempo...
y0 = [7150140; 0; 0; 0.000937]; % condiciones iniciales...
% resuelve el sistema de ecuaciones diferenciales de primer orden con "ode45"...
[t, y] = ode45(@sistemaODE45, tspan, y0);
