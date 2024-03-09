function [t,Y] = AdamsMoulton4Systems(f,a, b, N, Ya)
% Adams-Moulton methos for 3 step for equation systems
% y_k+1 = y_k + (h/24) * f(t_k-2, y_k-2) - 5*f(t_k-1, y_k-1) + 19*f(t_k,
% y_k) + 9*(t_k+1, y_k+1)
% Global error = forth order O(h^4)
% y_0 => known. 
% y_3, y_2, y_1 => Calculated by Runge-Kutta method
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
Y = zeros(N+1,length(Ya));
Y(1,:) = Ya;
maxiter = 10;
tol = 1e-6;
% Runge-Kutta Method
for k=1:3
    k1 = feval(f,t(k), Y(k,:))';
    k2 = feval(f, t(k) + h/2, Y(k,:) + h*k1/2)';
    k3 = feval(f,t(k)+h/2, Y(k,:) + h*k2/2)';
    k4 = feval(f,t(k+1), Y(k,:)+ h*k3)';
    Y(k+1,:) = Y(k,:) + h*(k1 + 2*k2 + 2*k3 + k4)/6; 
end
for k = 3:N
    fk = feval(f,t(k),Y(k,:));
    fkm1 = feval(f,t(k-1),Y(k-1,:));
    fkm2 = feval(f,t(k-2),Y(k-2,:));
    iter = 1;
    dif = tol + 1;
    x0 = Y(k,:)';
    while and(iter<maxiter, dif>tol)
        [fx0,dfx0] = feval(f,t(k+1),x0);
        g = x0 - Y(k,:) - h/24*(fkm2 -5*fkm1 + 19*fk + 9*fx0);
        dg = 1 - h/24*9*dfx0;
        x1 = x0 - dg \ g;
        dif=norm(x1-x0, inf);
        iter = iter + 1;
        x0 = x1;
    end
    % y(k+1) = y(k) + h/24*(fkm2 -5*fkm1 + 19*fk + 9*feval(f, t(k+1), x0));
    Y(k+1,:) = Y(k,:)' + h/24*(fkm2 -5*fkm1 + 19*fk + 9*feval(f, t(k+1), x0));
    % Y(k+1,:) = x1';
end
end