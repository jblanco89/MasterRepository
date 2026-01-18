function [nodos,solaprox,t,iter] = DisparoSecante(funcion,a,b,alfa,beta,n,tol,maxiter)
h=(b-a)/(n+1);  
x=a:h:b;  
x=x(:);

t0=0;
% [x,Y]= ode45(funcion,x,[alfa,t0]');
%[x,Y]= Heun_Systems(funcion,a,b,n, [t0,((1/2)+t0)]);
[x,Y]= AdamsBashforth4Systems(funcion,a,b,n,[t0,(-t0/2)]);
disp(Y)
% [x,Y]= ABM2Systems(funcion,a,b,n,[alfa,t0]);
yb0=Y(end,1);
yb0prima = Y(end,2);

% t1=(beta - alfa)/(b-a);
t1=1;
% [x,Y]= ode45(funcion,x,[alfa ,t1]');
% [x,Y]= Heun_Systems(funcion,a,b,n, [t1,((1/2)+t1)]);
[x,Y]= AdamsBashforth4Systems(funcion,a,b,n,[t1,(-t1/2)]);
% [x,Y]= ABM2Systems(funcion,a,b,n,[alfa,t1]);
yb1=Y(end,1);
yb1prima = Y(end,2);
iter=1; 
% incre=abs(yb1-beta);
incre = abs(yb1 - 2*yb1prima - beta);
% disp(incre)

% while incre > tol && iter < maxiter
while incre > tol && iter < maxiter
%     t=t1 -(t1 -t0)*(yb1 -beta )/(yb1 -yb0);
    t = t1 - (t1-t0)*(yb1 - 2*yb1prima - beta)/((yb1 - 2*yb1prima - beta) - (yb0 - 2*yb0prima - beta));
%    [x,Y]= Heun_Systems(funcion,a,b,n,[t,((1/2)+t)]);
     [x,Y]= AdamsBashforth4Systems(funcion,a,b,n,[t,-(t/2)]);
%     [x,Y]= ABM2Systems(funcion,a,b,n,[alfa,t]);
%     [x,Y]= ode45(funcion ,x ,[alfa ,0, t]');
    incre =abs(Y(end,1) - 2*Y(end,2) - beta); 
    iter= iter+1;
    t0=t1; 
    yb0=yb1;
    yb0prima = yb1prima;
    t1=t; 
    yb1=Y(end,1);
    yb1prima=Y(end,2);
end
% disp(incre)
if incre <= tol
nodos=x;
solaprox=Y;
else
disp('se necesitan mas iteraciones ')
end
end
