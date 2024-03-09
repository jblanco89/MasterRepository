function [ci, xi] = Coeficientes_Nodos_Gauss_Legendre(n)
    % Esta función calcula los coeficientes (ci) y los nodos (xi)
    % para la cuadratura de Gauss-Legendre con n nodos.

    % Definimos la variable simbólica
    syms x;

    % Calculamos los nodos (raíces del polinomio de Legendre)
    Pn = legendreP(n, x);
    xi = double(solve(Pn == 0));

    % Calculamos los coeficientes
    ci = zeros(n, 1);
    for i = 1:n
        % Calculamos la derivada del polinomio de Legendre en xi(i)
        Pn_prime = diff(legendreP(n, x), x);
        derivative_at_xi = subs(Pn_prime, x, xi(i));

        % Calculamos el peso correspondiente
        ci(i) = 2 / ((1 - xi(i)^2) * derivative_at_xi^2);
    end

    % Normalizamos los coeficientes
    ci = ci / sum(ci);

    % Ajuste para que coincidan con la función de cuadratura estándar
    ci = ci * sqrt(pi);

end
