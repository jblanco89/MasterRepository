function [t,y] = AdamsMoulton2(f,a, b, N, ya, Tol, maxiter)
% Adams-Moulton methos for 1 step
% y_k+1 = y_k + (h/2)* f(t_k+1, y_k+1) + f(t_k, y_k)
% Global error = second order O(h^2)
% y_0 => known. 
% y_k => Calculated by Heun method
% y_k +1 => Calculated by solving non-lineal equation (implicit method)
% ********************************************************************
% g(y_k+1) = y_k+1 - y_k - (h/2)*(f(t_k+1, y_k+1) + f(t_k, y_k)) = 0
% ********************************************************************
%
% This non-linear equation is usually solved by Newton-Raphson method
% ********************************************************************
% x_n+1 = x_n - (g(x_n) / g'(x_n))
% ********************************************************************
h = (b - a) / N;
t = a:h:b;
t = t(:);
y = zeros(N+1,1);
y(1) = ya;
tol = Tol;
% Heun Method
k1=h*feval(f,t(1),y(1));
k2=h*feval(f,t(2),y(1)+k1);
y(2)=y(1)+k1/2+k2/2;
for k = 2:N
    ff = feval(f,t(k),y(k));
    iter = 1;
    dif = tol + 1;
    x0 = y(k);
    while and(iter<maxiter, dif>tol)
        [fx0,dfx0] = feval(f,t(k+1),x0);
        g = x0 - y(k) - h/2*(fx0 + ff);
        dg = 1 - h/2*dfx0;
        x1 = x0 - g / dg;
        dif=abs(x1-x0);
        iter = iter + 1;
        x0 = x1;
    end
    y(k+1) = y(k) + h/2*(feval(f,t(k+1),x0) + ff);
end
end