function dydt = sistemaODE45(t, y)
    % Establece el sistema de ecuaciones diferenciales de primer orden para pasar por "ode45"...
    dydt = [y(2); y(1)*y(4)^2 - 3.986*10^14 / y(1)^2; y(4); -2*y(2)*y(4) / y(1)];
end