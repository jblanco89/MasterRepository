function [t, y] = AdamsBashforth2(f,a,b, N, ya)
% Adams-Bashforth methos for two steps
% y_k+1 = y_k + 3* (h/2)* f(t_k, y_k) - (h/2) * f(t_k-1, y_k-1)
% Global error = second order O(h^2)
% y_0 => known. 
% y_1, y_2, y_3 => calculated by Heun method
h = (b - a) / N;
t = a:h:b;
t = t(:);
y = zeros(N+1,1);
y(1) = ya;
k1 = h*feval(f,t(1), y(1));
k2 = h*feval(f, t(2), y(1) + k1);
y(2) = y(1) + (k1 + k2) / 2;
for k = 2:N
    y(k+1) = y(k) + 3 * (h/2) * feval(f,t(k), y(k)) - (h/2) * feval(f, t(k-1), y(k-1));
end
end