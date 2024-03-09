function [sol,iter,ACOC, incre1, incre2] = GR_Sist(F,x0,tol,maxit)
digits(200);
% Inicializacion de las variables
iter=1;
incre1=tol+1;
incre2=tol+1;
x0=x0(:);
[Fx,dFx]=feval(F,x0);
a = (- 1 + sqrt(5)) / 2;
b = (3 + sqrt(5)) / 2;
% Criterio de parada
while incre1>tol && incre2>tol && iter<maxit
    % Expresion del metodo de Golden Ratio
    u=a*(dFx\Fx);
    y0=x0-u;
    [Fy0, dFy0] = feval(F, y0);
    v = b*(dFx\Fy0);
    x1 = x0 - v;
    
    %vpa(x,5);
    incre1=norm(x1-x0);
    %incX=vpa(incre1,5)
    I(iter)=incre1;
    
    % Actualizacion de la estimacion inicial
    x0=x1;
    
    [Fx,dFx]=feval(F,x0);
    incre2=norm(Fx);
    %incF=vpa(incre2,5)
    %pause
    % Incremento del contador de iteraciones
    iter=iter+1;
end

iter=iter-1;
if length(I)>2
    sol=x0;
    ACOC=log(I(3:end)./I(2:end-1))./log(I(2:end-1)./I(1:end-2));
else
    disp('necesito más iteraciones')
end
 sol=vpa(x0,6);
 incre2=vpa(incre2,6);
 incre1=vpa(incre1,6);
 ACOC=vpa(ACOC,6);
end

