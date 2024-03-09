function [I, err] = Milne(f, y, a, b, n)
%Milne's method o also known as Boole's one, calculates area under the curve
%of a function by using numerical method. 
% In this case We're using composite Boole’s formula taken from:
% Sablonnière, P.; Sbibih, D.; Tahrichi, M. (2010). 
% "Error estimate and extrapolation of a quadrature formula derived from a 
% quartic spline quasi-interpolant". BIT Numerical Mathematics. 50:
% 843–862. doi:10.1007/s10543-010-0278-0

if n < 4
    disp('Número de puntos a evaluar debe ser mínimo 4 para este método')
end
h = (b - a) / n;
x = a:h:b;
w = ones(1, n+1);
w(1) = 7;
w(n+1) = 7;
w(2:2:n) = 32;
w(3:4:n-1) = 12;
w(5:4:n-1) = 14;
if y == 0
    I = 2*h*sum(w.*f(x))/45;
    syms x;
    solex = double(int(f(x), a, b));
    err = abs(solex - I);
else
    % y = f(x);
    y = y(1:length(w));
    I = 2*h*sum(w.*y)/45;
    der_6 = max(diff(y, 6));
    err = (2/945)*(h.^6)*(b-a) * der_6;  
end
end