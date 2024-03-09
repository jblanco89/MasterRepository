function [N1, N2] = richardson(f,x,h, order)
% Computes finite difference approximations using Richardson extrapolation
% Basically, it calculates finite difference approximations
% for the derivative of a given function f at a specified point x using Richardson
% extrapolation. It provides two approximations, N1 and N2, with different orders
% of accuracy.

%   Input:
%       f:      Function handle representing the function for which the derivative
%               is to be approximated.
%       x:      Point at which the derivative is approximated.
%       h:      Step size for finite difference calculations.
%       order:  Order of accuracy desired (4 or 6).
%
%   Output:
%       N1:     Finite difference approximation of the derivative using the
%               progressive order 2 method.
%       N2:     Finite difference approximation of the derivative with higher order
%               accuracy specified by the input 'order' (4 or 6).


% Diferencia finita progresiva de orden 2
N1 = ((4*f(x + h) - 3*f(x) - f(x + 2*h)) / (2*h));
if order == 4
    % Richardson O(h^4)
    N2 = (1/(6 * h)) * (f(x + 2 * h) - 12*f(x + h) + 32*f(x + 0.5 * h) - 21*f(x));
elseif order == 6
    % Richardson O(h^6)
    N2 = (1/(90 * h)) * (-f(x + (2*h)) + 44*f(x + h) - 416*f(x + (0.5*h)) + 1024*f(x + (0.25 * h)) - 651*f(x));
end
end