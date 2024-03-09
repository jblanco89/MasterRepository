function [t, Y] = AdamsBashforth2Systems(f,a,b,N, Ya)
% Adams-Bashforth methos for two steps in Equation Systems
% y_k+1 = y_k + 3* (h/2)* f(t_k, y_k) - (h/2) * f(t_k-1, y_k-1)
% Global error = second order O(h^2)
% Y_0 => known. 
% Y_1 => calculated by Heun method
h = (b - a) / N;
t = a:h:b;
t = t(:);
Y = zeros(N+1,length(Ya));
Y(1,:) = Ya;
for k = 1:2
    % y_k (y(k)) are estimated by Heun
    k1 = h*feval(f,t(k), Y(k,:))';
    k2 = h*feval(f, t(k+1), Y(k,:) + k1)';
    Y(k+1,:) = Y(k,:) + (k1 + k2) / 2;
end
for k = 3:N
    Y(k+1,:) = Y(k,:)' + 3 * (h/2) * feval(f,t(k), Y(k,:)) - (h/2) * feval(f, t(k-1), Y(k-1,:));
end
end