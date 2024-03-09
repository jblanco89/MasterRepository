function [t,y] = ABM2(f,a,b,N,ya)
% Predictor-Corrector method for two steps
%
% *****************************************
% Predictor: Adams-Bashforth order 2
% y_k+1 = y_k + 3* (h/2)* f(t_k, y_k) - (h/2) * f(t_k-1, y_k-1)
% *****************************************
%
% *******************************************
% Corrector: Adams-Moulton order 2
% y_k+1 = y_k + (h/2)* f(t_k+1, y_k+1) + f(t_k, y_k)
% *******************************************
%

h=(b-a)/N;
t=a:h:b;
t=t(:);
y=zeros(N+1,1) ;
y(1)=ya;
% Heun Method
k1 = h*feval(f,t(1), y(1));
k2 = h*feval(f, t(2), y(1) + k1);
y(2) = y(1) + (k1 + k2) / 2;
for k=2:N
    k1 = feval(f,t(k),y(k));
    k2 = feval(f,t(k-1) ,y(k-1));
    yp = y(k) + (h/2) *(3*k1-k2);
    y(k+1) =y(k)+h /2*(feval(f,t(k+1),yp)+k1);
end

end