function dy = finite_diff(Y,n,h,order)
    % calculates the second order finite difference approximation of 
    % the derivative of a given vector 'Y' with respect to the index, 
    % using different finite difference schemes. 
    % 
    % Parameters:
    % 
    % Y: The input vector.
    % n: The size of the vector Y.
    % h: The step size for the finite difference.
    % order: Specifies the finite difference order approximation 
    % ('first', 'upper')
    
    % Numerical Method 1 - UNIR
    % First Quarter
    % Cohort 2023 - 2024

    if strcmp(order, 'upper') == 1
        dy(1) = ((4*Y(2) - 3*Y(1) - Y(3)) / (2*h));
        dy(n) = ((3*Y(n) - 4*Y(n-1) + Y(n-2)) / (2*h));
        dy(2) = ((4*Y(3) - 3*Y(2) - Y(4)) / (2*h));
        dy(n-1) = ((3*Y(n-1) - 4*Y(n-2) + Y(n-3)) / (2*h));
        for i=3:n-2
            dy(i) = ((-Y(i+2) + 8*Y(i+1) - 8*Y(i-1) + Y(i-2))/(12*h));
            % % dy(1,i) = ((Y(i+1) - Y(i)) / (2*h));
        end
    elseif strcmp(order, 'first') == 1
        dy(1) = (Y(2) - Y(1)) / h;
        dy(n) = (Y(n) - Y(n-1)) / h;
        for i=2:n-1
            dy(i) = (Y(i+1) - Y(i-1)) / 2*h;
        end
    end
end