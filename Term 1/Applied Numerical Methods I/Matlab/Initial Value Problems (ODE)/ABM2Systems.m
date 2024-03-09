function [t,Y] = ABM2Systems(f,a,b,N,Ya)
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
Y=zeros(N+1,length(Ya));
Y(1,:)=Ya;
% Heun Method
k1 = h*feval(f,t(1) ,Y(1,:))';
k2 = h*feval(f,t(2) ,Y(1,:)+k1)';
Y(2,:) = Y(1,:) +(k1+k2)/2;
for k=2:N
    k1 = feval(f,t(k),Y(k,:));
    k2 = feval(f,t(k-1),Y(k-1,:));
    yp = Y(k,:) + h/2*(3*k1-k2);
    Y(k+1,:) =Y(k,:)'+ h/2*(feval(f,t(k+1),yp)+k1);
end

end