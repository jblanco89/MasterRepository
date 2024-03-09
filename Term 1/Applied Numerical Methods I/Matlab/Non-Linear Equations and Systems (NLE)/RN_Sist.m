function [sol,iter,ACOC, incre1, incre2] = RN_Sist(F,x0,tol,maxit)
digits(200);
% Inicializacion de las variables
iter=1;
incre1=tol+1;
incre2=tol+1;
x0=x0(:);
[Fx,dFx]=feval(F,x0);

% Criterio de parada
while incre1>tol && incre2>tol && iter<maxit
    % EXPRESION NEWTON JARRAT
    u=dFx\Fx;
    z=x0-2/3*u;
    
    [Fz,dFz]=feval(F,z);
    
    w=(3*dFz-dFx)\(3*dFz+dFx);
    
    y=x0-1/2*(w)*u;
    
    [Fy,dFy]=feval(F,y);
    %x=y-dFy\Fy;
    a=-0.5; b=1.5;
    x=y-(a*dFx-b*dFz)\Fy;
    
    
    
    
    %vpa(x,5);
    incre1=norm(x-x0);
    %incX=vpa(incre1,5)
    I(iter)=incre1;
    
    % Actualizacion de la estimacion inicial
    x0=x;
    
    [Fx,dFx]=feval(F,x0);
    incre2=norm(Fx);
    %incF=vpa(incre2,5)
    %pause
    % Incremento del contador de iteraciones
    iter=iter+1;
end

iter=iter-1;
if length(I)>2
    sol=x;
    ACOC=log(I(3:end)./I(2:end-1))./log(I(2:end-1)./I(1:end-2));
else
    disp('necesito mas iteraciones')
end
 sol=vpa(x0,6);
 incre2=vpa(incre2,6);
 incre1=vpa(incre1,6);
 ACOC=vpa(ACOC,6);
end

