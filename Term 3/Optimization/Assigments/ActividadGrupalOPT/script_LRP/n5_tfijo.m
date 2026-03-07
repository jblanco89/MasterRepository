clc; clear; close all;
format longg;

n = 5;
c = 1/(2*n);

f = @(u) ((u(2)+u(1)) * sqrt(1 + (n*(u(2)-u(1)))^2) + ...
          (u(3)+u(2)) * sqrt(1 + (n*(u(3)-u(2)))^2) + ...
          (u(4)+u(3)) * sqrt(1 + (n*(u(4)-u(3)))^2) + ...
          (u(5)+u(4)) * sqrt(1 + (n*(u(5)-u(4)))^2) + ...
          (u(6)+u(5)) * sqrt(1 + (n*(u(6)-u(5)))^2)) * c;

grad_f = @(u) [(c - u(1)*n*(u(2)-u(1))) / sqrt(1 + (n*(u(2)-u(1)))^2);
               (c - u(2)*n*(u(3)-u(2))) / sqrt(1 + (n*(u(3)-u(2)))^2) + ...
               (c + u(2)*n*(u(2)-u(1))) / sqrt(1 + (n*(u(2)-u(1)))^2);
               (c - u(3)*n*(u(4)-u(3))) / sqrt(1 + (n*(u(4)-u(3)))^2) + ...
               (c + u(3)*n*(u(3)-u(2))) / sqrt(1 + (n*(u(3)-u(2)))^2);
               (c - u(4)*n*(u(5)-u(4))) / sqrt(1 + (n*(u(5)-u(4)))^2) + ...
               (c + u(4)*n*(u(4)-u(3))) / sqrt(1 + (n*(u(4)-u(3)))^2);
               (c - u(5)*n*(u(6)-u(5))) / sqrt(1 + (n*(u(6)-u(5)))^2) + ...
               (c + u(5)*n*(u(5)-u(4))) / sqrt(1 + (n*(u(5)-u(4)))^2);
               (c + u(6)*n*(u(6)-u(5))) / sqrt(1 + (n*(u(6)-u(5)))^2)];

disp(['Función a minimizar: f(u)' newline])

u0 = ones(n+1, 1); % punto de inicio

t_k = 0.15; % paso
delta = 1e-6; % tolerancia de convergencia
max_iter = 1000; % número máximo de iteraciones
iter = 0;

error_1 = delta + 1; % error ||u_{k+1} - u_k||
error_2 = delta + 1; % error ||grad_f(u_{k+1})||
u_k = u0;
u_hist = u_k'; % guarda historial de iteraciones
x0 = 0;
xn = 1;

tabla_iteraciones = [];

% Preparar grabación de video
% Crear carpeta para guardar resultados
output_folder = 'resultados_JG';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

video_filename = fullfile(output_folder, sprintf('gradiente_n%d_t%.3f.mp4', n, t_k));
video_writer = VideoWriter(video_filename, 'MPEG-4');
open(video_writer);


%while (error_1 > delta || error_2 > delta) && iter < max_iter
while (error_1 > delta && error_2 > delta) && iter < max_iter
    iter = iter + 1;
    gradiente = grad_f(u_k);

    u_k_plus_1 = u_k - t_k .* gradiente;
    u_k_plus_1(1) = 1; u_k_plus_1(n+1) = 1; % exijo las condiciones de contorno!!!
    error_1 = norm(u_k_plus_1 - u_k)
    error_2 = norm(grad_f(u_k_plus_1))
    tabla_iteraciones = [tabla_iteraciones; iter, u_k(1), u_k(2), u_k(3), u_k(4), u_k(5), u_k(6), f(u_k), error_1, error_2];
    u_k = u_k_plus_1;
    u_hist = [u_hist; u_k'];

    x_values = linspace(x0, xn, n+1);
            
    % Graficar puntos discretos
    plot(x_values, u_k, '-o', 'LineWidth', 1.5);
    grid on;
    title(['Iteración: ', num2str(iter), ' (n=', num2str(n), ', t_k=', num2str(t_k), ')']);
    xlabel('Intervalo X');
    ylabel('Valores de U');
    xlim([x0 xn]);
    ylim([min(u_k) - 0.1, max(u_k) + 0.1]); % Ajustar los límites de Y dinámicamente
    pause(0.1); % Pausa para animación

        % Guardar frame en el video
    frame = getframe(gcf);
    writeVideo(video_writer, frame);
end

T = array2table(tabla_iteraciones, 'VariableNames', {'iter', 'u_1', 'u_2', 'u_3', 'u_4', 'u_5', 'u_6', 'f(u)', 'error_1', 'error_2'});
disp(T);
u_opt = u_k;
f_opt = f(u_opt);
fprintf('Solución óptima encontrada en la iteración %d:\n', iter);
fprintf('u_1 = %f\n', u_opt(1));
fprintf('u_2 = %f\n', u_opt(2));
fprintf('u_3 = %f\n', u_opt(3));
fprintf('u_4 = %f\n', u_opt(4));
fprintf('u_5 = %f\n', u_opt(5));
fprintf('u_6 = %f\n', u_opt(6));
fprintf('Valor mínimo de f(u) = %f\n', f_opt);
