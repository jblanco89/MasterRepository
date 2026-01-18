function [P, err, coefficients] = lagrange_interpolation(x, y, n, t)
    % Fits Lagrange interpolation polynomial of degree n
    %
    % Parameters:
    %
    %   x: x-coordinates vector
    %   y: y-coordinates vector
    %   n: Lagrange polynomial degree 
    %   t: Points where the polynomial will be evaluated
    %
    % Return:
    %
    %   P: Lagrange polynomial values at points t
    % err: relative error between polyfit values and P
    % coefficients: Coefficients of the Lagrange polynomial
    %
    % Numerical Method 1 - UNIR
    % First Quarter
    % Cohort 2023 - 2024

    format longG;

    if length(x) ~= length(y) || n >= length(x)
        error(['please check input dim because must be' ...
            'length(x) == length(y) and n < length(x)']);
    end

    coefficients = zeros(1, n + 1);
    
    P = zeros(size(t));

    for k = 1:length(x)
        Lk = ones(size(t));

        for j = 1:length(x)
            if j ~= k
                Lk = Lk .* (t - x(j)) / (x(k) - x(j));
            end
        end
        P = P + y(k) * Lk;
    end
   
   poly_fit = polyfit(x, y, n);
   coefficients = poly_fit;
   poly_val = polyval(poly_fit, t);
   err = abs((P - poly_val) ./ poly_val)*100;

end
