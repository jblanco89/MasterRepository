function [t, Y] = AdamsBashforth4Systems(f,a,b,N,Ya)
% Adams-Bashforth method for four steps
% y_k+1 = y_k + (h/24) * 55*f(t_k, y_k) - 59*f(t_k-1, y_k-1) + 37*f(t_k-2,
% y_k-2) - 9*f(t_k-3, y_k-3)
% Global error = second order O(h^4)
% y_0 => known. 
% y_1, y_2, y_3 => calculated by Runge-Kutta 4th Order
h = (b - a) / N;
t = a:h:b;
t = t(:);
Y = zeros(N+1,length(Ya));
Y(1,:) = Ya;
for k = 1:3
    % y_k (y(k)) are estimated by Runge-Kutta
    k1 = feval(f,t(k), Y(k,:))';
    k2 = feval(f, t(k) + h/2, Y(k,:) + h*k1/2)';
    k3 = feval(f,t(k)+h/2, Y(k,:) + h*k2/2)';
    k4 = feval(f,t(k+1), Y(k,:)+ h*k3)';
    Y(k+1,:) = Y(k,:) + h*(k1 + 2*k2 + 2*k3 + k4) / 6;
end
for k = 4:N
    Y(k+1,:) = Y(k,:)' + (h/24) * (55 * feval(f,t(k), Y(k,:)) - ...
        59 * feval(f, t(k-1), Y(k-1,:)) + ...
        37 * feval(f, t(k-2), Y(k-2,:)) - ...
        9 * feval(f, t(k-3), Y(k-3,:)));
end
end