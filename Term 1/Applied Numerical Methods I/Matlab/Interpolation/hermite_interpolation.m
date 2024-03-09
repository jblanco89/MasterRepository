function [P, err, M] = hermite_interpolation(x, y, dy, n, t)
%   Computes Hermite interpolation using divided differences method.
%   It interpolates the values of 't'
%   using the Hermite interpolation method.

%   Inputs:
%       x - Vector of x-coordinates of the data points
%       y - Vector of y-coordinates of the data points
%       dy - Vector of derivatives at the corresponding data points
%       n - Number of data points
%       t - Vector of points where interpolation is desired
%
%   Outputs:
%       P - Vector of interpolated y-values corresponding to the points in t
%       err - Vector of relative errors between interpolated and true values
%       M - Matrix containing the divided differences
%
%   Polynomial grade is 2n + 1.
%
%   Example:
%       x = [1, 2, 3, 4];
%       y = [3, 6, 10, 12];
%       dy = [2, 4, 8, 10];
%       t = [2.5, 3.5];
%       [P, err, M] = hermite_interpolation(x, y, dy, length(x), t);
%
% Polynomial grade is 2n + 1
% Divided differences method (See Newton's method)
if length(x) ~= length(y)
    error('Error, ensure length(x) == length(y)');
end

z = zeros(2*n,1);
Q = zeros(2*n, 2*n);

for i = 0:n-1
    z(2*i+1)=x(i+1);
    z(2*i+2)=x(i+1);
    Q(2*i+1,1)=y(i+1);
    Q(2*i+2,1)=y(i+1);
    Q(2*i+2,2)=dy(i+1);
    if i ~= 0
        Q(2*i+1,2)=(Q(2*i+1,1)-Q(2*i,1))/(z(2*i+1)-z(2*i));
    end
end
for i = 2 : 2*(n-1)+1
    for j = 2 : i
        Q(i+1,j+1)=(Q(i+1,j)-Q(i,j))/(z(i+1)-z(i-j+1));
    end
end
M = [z Q];
% disp(M)
k = 2 * n;
P = 0;
for i = 1:k
    prod = 1;
    for j = 1:i - 1
        prod = prod .* (t - M(j, 1));
    end
    P = P + M(i, i+1) * prod;
%     coefficients(i) = M(i, i+1);
%     P = P + coefficients(i) * prod;

end
poly_fit = polyfit(x, y, n);
poly_val = polyval(poly_fit, t);
err = (poly_val - P) ./ poly_val;
end