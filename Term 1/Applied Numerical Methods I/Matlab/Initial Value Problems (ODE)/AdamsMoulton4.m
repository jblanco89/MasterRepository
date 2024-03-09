function [t,y] = AdamsMoulton4(f,a, b, N, ya, Tol, maxiter)
% Adams-Moulton methos for 3 step
% y_k+1 = y_k + (h/24) * f(t_k-2, y_k-2) - 5*f(t_k-1, y_k-1) + 19*f(t_k,
% y_k) + 9*(t_k+1, y_k+1)
% Global error = forth order O(h^4)
% y_0 => known. 
% y_2, y_1 => Calculated by Runge-Kutta method
% y_k+1 => Calculated by solving non-lineal equation (implicit method)
% ********************************************************************
% g(y_k+1) = y_k+1 - y_k - (h/24)*(f(t_k-2, y_k-2) + 5*f(t_k-1, y_k-1) + 19*f(t_k, y_k) + 9*f(t_k+1, y_k+1)) = 0
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
% Runge-Kutta Method
for k=1:2
    k1 = feval(f,t(k), y(k));
    k2 = feval(f, t(k) + h/2, y(k) + h*k1/2);
    k3 = feval(f,t(k)+h/2, y(k) + h*k2/2);
    k4 = feval(f,t(k+1), y(k)+ h*k3);
    y(k+1) = y(k) + h*(k1 + 2*k2 + 2*k3 + k4)/6; 
end
for k = 3:N
    fk = feval(f,t(k),y(k));
    fkm1 = feval(f,t(k-1),y(k-1));
    fkm2 = feval(f,t(k-2),y(k-2));
    iter = 1;
    dif = tol + 1;
    x0 = y(k);
    while and(iter<maxiter, dif>tol)
        [fx0,dfx0] = feval(f,t(k+1),x0);
        g = x0 - y(k) - h/24*(fkm2 -5*fkm1 + 19*fk + 9*fx0);
        dg = 1 - h/24*9*dfx0;
        x1 = x0 - g / dg;
        dif=abs(x1-x0);
        iter = iter + 1;
        x0 = x1;
    end
    y(k+1) = y(k) + h/24*(fkm2 -5*fkm1 + 19*fk + 9*feval(f, t(k+1), x0));
end
end