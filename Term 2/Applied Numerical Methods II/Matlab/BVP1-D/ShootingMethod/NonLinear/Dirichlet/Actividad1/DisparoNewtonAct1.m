function [nodos,solaprox,t,iter] = DisparoNewtonAct1(funcion,a,b,alfa,beta,n,tol,maxiter)

h=(b-a)/(n+1); 
x=a:h:b; 
% t0 =(beta-alfa)/(b-a);
t0=0;

% [x,Y]= ode45(funcion,x,[alfa,t0,0,1]');
% [x,Y]= Heun_Systems(funcion,a,b,n, [t0,((1/2) + t0),1,1/2]);
[x,Y]= AdamsBashforth4Systems(funcion,a,b,n, [t0,(1/2 + t0),1,1]);
% [x,Y]= ABM2Systems(funcion,a,b,n, [alfa,t0,0,1]);
yb1=Y(end,1); 
zb1=Y(end,3);
dzb1 = Y(end,4);
dy1 = Y(end,2);
iter=1; 
% incre=abs(Y(end,1)-beta);
incre = abs(dy1 + yb1 - beta);
% zb1 = dy1 + yb1;

while incre > tol && iter < maxiter
    t= t0 - (dy1+yb1-beta)/(zb1 + dzb1); %Metodo de Newton
%     [x,Y]= Heun_Systems(funcion,a,b,n,[t,((1/2) + t),1,1/2]);
    [x,Y]= AdamsBashforth4Systems(funcion,a,b,n, [t,(1/2 + t),1,1]);
%     [x,Y]= ode45(funcion,x,[alfa,t,0,1]');
%     [x,Y]= ABM2Systems(funcion,a,b,n, [alfa,t,0,1]);
%     incre=abs(Y(end,1)-beta);
    t0=t;
    iter= iter+1; 
    yb1=Y(end,1); 
    zb1=Y(end,3);
    dzb1 = Y(end,4);
    dy1 = Y(end,2);
    incre = abs(Y(end,2) + Y(end,1) - beta);
end
if incre <= tol
nodos =x; 
solaprox=Y;
else
disp('se necesitan mas iteraciones ')
end
end