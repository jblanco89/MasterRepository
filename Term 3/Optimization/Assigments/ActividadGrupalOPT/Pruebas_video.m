%% Script: Pruebas.m
clc; clear; close all;
format longg;

%% parámetros...
n_vector = [5, 10, 20, 50];
t_vector = [0.17, 0.06, 0.03, 0.01];
delta_1 = 1e-6;
delta_2 = 0.75;
max_k = 10000;

%% salidas...
l = length(n_vector);

for i = 2
    n = n_vector(i);
    t = t_vector(i);

    [tabla_iteraciones, u_min_aprox, u_hist, k, texec] = metodo_gradiente_paso_fijo(n, t, delta_1, delta_2, max_k);

    % calculo valores exactos...
    x = 0:1/n:1;
    a = 0.848338;
    u_min_exacta = a*cosh((x-1/2)/a); % catenaria exacta
    I_min_exacta = (a/2)*(1 + a*sinh(1/a)); % se ha obviado el 2*pi

    fprintf('* Solución óptima encontrada en la iteración %d, para n = %d, t = %f\n\n', k, n, t);
    
    % tabulo...
    disp(array2table(tabla_iteraciones(end-3:end,:), 'VariableNames', {'k', 'F(u)', '||u_{k+1} - u_k||', '||grad_F(u_{k+1})||'}));
    F_min_aprox = F(u_min_aprox);
    fprintf('Valor mínimo de F(u) = F_min_aprox = %f\n|F_min_aprox - I_min_exacta| = %f\n\n', F_min_aprox, abs(F_min_aprox - I_min_exacta));

    % grafico...
    hold on; grid on;
    s = size(u_hist);

    output_folder = 'videos';

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    video_filename = fullfile(output_folder, sprintf('demo_n%d.mp4', n));
    video_writer = VideoWriter(video_filename, 'MPEG-4');
    open(video_writer);

    for v = 1:round(n/2):s(1)
        p0 = plot(x, u_hist(v,:)', '.-', 'Color', 'r', 'MarkerSize', 10);
        p2 = plot(x, u_min_exacta, '.-', 'Color', 'b', 'MarkerSize', 10);
        xlabel('$$x$$', 'Interpreter', 'latex'); ylabel('$$u(x)$$', 'Interpreter', 'latex');
        legend([p0, p2], '$$u_{k}(x)$$', '$$u_{\min}(x)$$', 'Interpreter', 'latex');
        title(['iter = ', num2str(v), ',   n = ', num2str(n), ',   t = ', num2str(t), ',   |F(u) - I(u)| = ', num2str(abs(F(u_hist(v,:)) - I_min_exacta))]);
        pause(0.5);
        frame = getframe(gcf);
        writeVideo(video_writer, frame);
        delete(p0);
    end

    delete(p0);
    p1 = plot(x, u_min_aprox', '.-', 'Color', 'r', 'MarkerSize', 10);
    p2 = plot(x, u_min_exacta, '.-', 'Color', 'b', 'MarkerSize', 10);
    xlabel('$$x$$', 'Interpreter', 'latex'); ylabel('$$u(x)$$', 'Interpreter', 'latex');
    legend([p1, p2], '$$u_{k}(x)$$', '$$u_{\min}(x)$$', 'Interpreter', 'latex');
    title(['iter = ', num2str(k), ',   n = ', num2str(n), ',   t = ', num2str(t), ',   |F(u) - I(u)| = ', num2str(abs(F_min_aprox - I_min_exacta))]);
    pause(0.5);
    frame = getframe(gcf);
    writeVideo(video_writer, frame);

    close(video_writer);
end

%% método del gradiente...
function [tabla_iteraciones, u, u_hist, k, texec] = metodo_gradiente_paso_fijo(n, t, delta_1, delta_2, max_k)
    k = 0; % índice de iteración
    error_1 = delta_1 + 1; % error ||u_{k+1} - u_k||
    error_2 = delta_2 + 1; % error ||grad_f(u_{k+1})||
    u = ones(n+1, 1); % vector de inicio
    u_hist = u'; % guarda historial de iteraciones
    tabla_iteraciones = [];

    tic;
    while (error_1 > delta_1 || error_2 > delta_2) && k < max_k
        k = k + 1;
        u_nuevo = u - t * grad_F(u);
        u_nuevo(1) = 1; u_nuevo(n+1) = 1; % exijo las condiciones de contorno!!!
        error_1 = norm(u_nuevo - u);
        error_2 = norm(grad_F(u_nuevo));
        tabla_iteraciones = [tabla_iteraciones; k, F(u), error_1, error_2];
        u = u_nuevo;
        u_hist = [u_hist; u'];
    end
    texec = toc;
end

%% función integral...
function integral = F(u)
    n = length(u) - 1;
    integral = 0;

    for j = 1:n
        integral = integral + (u(j+1)+u(j)) * sqrt(1 + (n*(u(j+1)-u(j)))^2);
    end

    integral = integral/(2*n);
end

%% función gradiente...
function gradiente = grad_F(u)
    l = length(u);
    n = l - 1;
    c = 1/(2*n);

    gradiente(1) = (c - u(1)*n*(u(2)-u(1))) / sqrt(1 + (n*(u(2)-u(1)))^2);

    for j = 2:n
        gradiente(j) = (c + u(j)*n*(u(j)-u(j-1))) / sqrt(1 + (n*(u(j)-u(j-1)))^2) + ...
                    (c - u(j)*n*(u(j+1)-u(j))) / sqrt(1 + (n*(u(j+1)-u(j)))^2);
    end

    gradiente(l) = (c + u(l)*n*(u(l)-u(l-1))) / sqrt(1 + (n*(u(l)-u(l-1)))^2);

    gradiente = gradiente';
end
