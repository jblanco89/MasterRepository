clear, clc
digits(16)
format shortG

% ********************************************************************%
% ******************* BENCHMARK **************************************%
%% Con tolerancia constante y 16 digitos de precisión

max_iter = 10000;   % máximo de iteraciones
tol = 1e-4; % tolerancia máxima
n_values = [15, 20, 25, 50]; % rectas inexactas no tolera n < 15
integral_methods = {'simpson', 'trapezium'};
t_opt_methods = {'rec', 'fibo', 'dico'};

results = table('Size', [numel(n_values)*numel(integral_methods)*numel(t_opt_methods), 5], ...
    'VariableTypes', {'double', 'double', 'double', 'string', 'string'}, ...
    'VariableNames', {'n', 'min_area', 't_opt', 'integral_method', 't_opt_method'});

row = 1;
for n = n_values
    for integral_method = integral_methods
        for t_opt_method = t_opt_methods
            [u_opt, min_area, t_opt] = min_rev_surface(n, max_iter, tol, integral_method{1}, t_opt_method{1});
            results.n(row) = n;
            results.min_area(row) = min_area;
            results.t_opt(row) = t_opt;
            results.integral_method(row) = string(integral_method{1});
            results.t_opt_method(row) = string(t_opt_method{1});
            row = row + 1;
        end
    end
end

disp(results)