function [sol,iter,ACOC,inc] = Newton(fun,x0,tol,maxiter)
digits(200)
x0=x0(:);
iter=0;
[fx0,dfx0]=feval(fun,x0);
incre1=tol+1;
incre2=tol+1;
p=[];
% inc = tol+1;
% while incre2>tol && incre1>tol && iter<maxiter
while incre2 + incre1 > tol && iter<maxiter
    
    %Linea NEWTON
    x1=x0-fx0/dfx0;
    %
    
    %actualizo criterio de parada
    incre1=norm(x1-x0);
    p=[p incre1];
    x0=x1;
    [fx0,dfx0]=feval(fun,x0);
    incre2=norm(fx0);
    iter=iter+1;
end
% calculo de ACOC
ACOC=log(p(3:end)./p(2:end-1))./log(p(2:end-1)./p(1:end-2));

sol=x1;
incre1=vpa(incre1,6);
incre2=vpa(incre2,6);
inc = vpa((incre1 + incre2), 6);
ACOC=vpa(ACOC,6);
ACOC=ACOC(:);
sol=vpa(sol,10);
end
