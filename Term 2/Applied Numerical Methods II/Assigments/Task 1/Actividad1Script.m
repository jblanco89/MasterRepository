clear, clc
digits(100)
format longG

%PARAMETERS
a = 1.0;
b = 2.0;
alfa = 1/2;
beta = 17/2;
n = 10;
tol = 1e-8;
maxiter = 10;
I = 0:1:n;

%METHODS

[xs,ys,ts,iter_s]=DisparoSecanteAct1('actividad1_s',a,b,alfa,beta,n,tol,maxiter);
[xn,yn,tn,iter_n]=DisparoNewtonAct1('actividad1_n',a,b,alfa,beta,n,tol,maxiter);


t = table(I', xs,ys(:,1),ys(:,2),yn(:,1),yn(:,2), ...
    'VariableNames',{'I','x','Y Secante', ...
    'dY Secante', 'Y Newton', 'dY Newton'});
disp(t)


% hold on
% plot(xs, ys(:,1),'*b')
% % plot(xs, ys(:,2), '-r')
% plot(xn, yn(:,1), 'go')
% % plot(xn, yn(:,2), '-b')
% grid on
% legend('Y Secante', 'Y Newton')

% Gráfico
hold on
plot(xs, ys(:,1), '*b', 'LineWidth', 1.5) % Estilo de línea mejorado
plot(xn, yn(:,1), 'go', 'MarkerSize', 9) % Aumenta el tamaño del marcador
grid on
legend('Secante', 'Newton', 'Location', 'best') % Mueva la leyenda a una posición óptima
xlabel('Nodos') % Etiqueta del eje X
ylabel('Solución de Y(x)') % Etiqueta del eje Y
title('Método del Disparo con Secante y Newton', 'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom') % Título del gráfico