clear, clc
digits(32)
format shortG

n = 50;           % número de puntos
max_iter = 10000;   % máximo de iteraciones
tol = 1e-8; % tolerancia máxima
integral_method = 'simpson';
% integral_method = 'trapezium';
% t_opt_metod = 'rec'; % rectas inexactas
t_opt_metod = 'dico';
% t_opt_metod = 'fibo';
[u_opt, min_area, alpha_opt] = min_rev_surface(n, max_iter, tol,integral_method, t_opt_metod);