function [P, err] = newton_interpolation(x, y, n, t)
    % Fits Newton interpolation polynomial of degree n
    %
    % Inputs:
    %
    %   x: x-coordinates vector
    %   y: y-coordinates vector
    %   n: Newton polynomial degree
    %   t: Points where the polynomial will be evaluated
    %
    % Returns:
    %
    %   P: Newton polynomial values at points t (Coeffs)

    if length(x) ~= length(y)
        error('Error, ensure length(x) == length(y)');
    end

    F = zeros(n, n);
    F(:, 1) = y(1:n)';

    for j = 2:n
        for i = j:n
            F(i, j) = (F(i, j - 1) - F(i - 1, j - 1)) / (x(i) - x(i - j + 1));
        end
    end

    P = zeros(size(t));
    product_term = 1;

    for j = 1:n
        P = P + F(j, j) * product_term;
        product_term = product_term .* (t - x(j));
    end
   poly_fit = polyfit(x, y, n);
   poly_val = polyval(poly_fit, t);
   err = abs((P - poly_val) ./ poly_val);
end
