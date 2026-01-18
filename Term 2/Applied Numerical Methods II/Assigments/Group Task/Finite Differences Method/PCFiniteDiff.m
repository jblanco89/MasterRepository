% Este script de MATLAB realiza un análisis comparativo entre dos métodos numéricos, Gauss y Crout,
% para resolver ecuaciones diferenciales no lineales. Utiliza diferentes tamaños de paso y mide el
% tiempo de ejecución de cada método para estos tamaños de paso. Los resultados se presentan en una
% tabla y un gráfico, mostrando los tiempos medios de ejecución de ambos métodos para cada tamaño de
% paso. Esto proporciona una evaluación de la eficiencia relativa de los métodos en función del tamaño
% de paso utilizado.

clear, clc
digits(100)
format shortG

f = @(x, y, z) -z.^2 + 4.*y + 2;
fy = @(x, y, z) 4;
fz = @(x, y, z) -2.*z;
a = 1;
b = 2;
alfa = 0;
beta = -3;
maxiter = 50;
tol = 1e-5;

Ns = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000];
num_iterations = 10;

tiempos_Gauss = zeros(num_iterations, length(Ns));
tiempos_Crout = zeros(num_iterations, length(Ns));

for j = 1:length(Ns)
    N = Ns(j);
    for k = 1:num_iterations
        tic;
        [~, ~, ~, ~] = DifnolinearGauss(f, fy, fz, a, b, alfa, beta, N, maxiter, tol);
        tiempos_Gauss(k, j) = toc;
        
        tic;
        [~, ~, ~, ~] = DifnolinearCrout(f, fy, fz, a, b, alfa, beta, N, maxiter, tol);
        tiempos_Crout(k, j) = toc;
    end
end

tiempos_medios_Gauss = mean(tiempos_Gauss, 1);
tiempos_medios_Crout = mean(tiempos_Crout, 1);

% Tabla de tiempos medios
tabla_tiempos_medios = table(Ns', tiempos_medios_Gauss', tiempos_medios_Crout', ...
    'VariableNames', {'N', 'Tiempo Medio Gauss', 'Tiempo Medio Crout'});
disp(tabla_tiempos_medios);

% Graficar tiempos medios
figure;
bar(Ns, [tiempos_medios_Gauss; tiempos_medios_Crout]');
grid on;
xlabel('Subintervalos (N)');
ylabel('Tiempo Medio (s)');
legend('Gauss', 'Crout', 'Location','north');
title('Tiempos Medios de Ejecución para Gauss y Crout');
