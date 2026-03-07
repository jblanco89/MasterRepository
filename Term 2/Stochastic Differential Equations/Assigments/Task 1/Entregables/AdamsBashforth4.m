function [t,y] = AdamsBashforth4(f,a,b, N, ya)
% Adams-Bashforth method for four steps
% y_k+1 = y_k + (h/24) * 55*f(t_k, y_k) - 59*f(t_k-1, y_k-1) + 37*f(t_k-2,
% y_k-2) - 9*f(t_k-3, y_k-3)
% Global error = forth order O(h^4)
% y_0 => known. 
% y_1, y_2, y_3 => calculated by Runge-Kutta 4th Order
h = (b - a) / N;
t = a:h:b;
t = t(:);
y = zeros(N+1,1);
y(1) = ya;
for k = 1:3
    % y_k (y(k)) are estimated by Runge-Kutta
    k1 = feval(f,t(k), y(k));
    k2 = feval(f, t(k) + h/2, y(k) + h*k1/2);
    k3 = feval(f,t(k)+h/2, y(k) + h*k2/2);
    k4 = feval(f,t(k+1), y(k)+ h*k3);
    y(k+1) = y(k) + h*(k1 + 2*k2 + 2*k3 + k4) / 6;
end
for k = 4:N
    y(k+1) = y(k) + (h/24) * (55 * feval(f,t(k), y(k)) - ...
        59 * feval(f, t(k-1), y(k-1)) + ...
        37 * feval(f, t(k-2), y(k-2)) - ...
        9 * feval(f, t(k-3), y(k-3)));
end
end