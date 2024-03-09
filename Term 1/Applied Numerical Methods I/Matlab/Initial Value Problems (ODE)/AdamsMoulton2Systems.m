function [t,Y] = AdamsMoulton2Systems(f,a, b, N, Ya)
% Adams-Moulton method for 1 step for systems
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
Y = zeros(N+1,length(Ya));
Y(1,:) = Ya;
maxiter = 10;
tol = 1e-6;
% Heun Method
k1=h*feval(f,t(1),Y(1,:))';
k2=h*feval(f,t(2),Y(1,:)+k1)';
Y(2,:)=Y(1,:)+k1/2+k2/2;
for k = 2:N
    ff = feval(f,t(k),Y(k,:));
    iter = 1;
    dif = tol + 1;
    x0 = Y(k,:)';
    while and(iter<maxiter, dif>tol)
        [fx0,dfx0] = feval(f,t(k+1),x0);
        g = x0 - Y(k,:) - h/2*(fx0 + ff);
        dg = 1 - h/2*dfx0;
        x1 = x0 - dg \ g;
        dif=norm(x1-x0, inf);
        iter = iter + 1;
        x0 = x1;
    end
    % Y(k+1,:) = x1';
    Y(k+1,:) = Y(k,:)' + h/2*(feval(f,t(k+1),x0) + ff);
end
end