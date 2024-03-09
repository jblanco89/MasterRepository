function [ci, xi] = Coeficientes_Nodos_Gauss_Chebyshev(n)
    % Esta función calcula los coeficientes (ci) y los nodos (xi)
    % para la cuadratura de Gauss-Chebyshev con n nodos.

    % Definimos la variable simbólica
    syms x;

    % Calculamos los nodos (raíces del polinomio de Chebyshev de primera especie)
    Tn = chebyshevT(n, x);
    xi = double(solve(Tn == 0));

    % Calculamos los coeficientes
    ci = pi / n;

end
