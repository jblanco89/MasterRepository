function [sol,iter,ACOC,incre1,incre2] = Secante(fun,x0, x1,tol,maxiter)
digits(200)
x0=x0(:);
x1=x1(:);
iter=0;
[fx0,dfx0]=feval(fun,x0);
[fx1,dfx1]=feval(fun,x1);
incre1=tol+1;
incre2=tol+1;
p=[];
while incre2>tol && incre1>tol && iter<maxiter
    
    %Linea Secante
    x2 = x1 - (fx1 * (x1 - x0)) / (fx1 - fx0);
    %
    
    %actualizo criterio de parada
    incre1=norm(x2-x1);
    p=[p incre1];
    x0 = x1;
    fx0 = fx1;
    x1 = x2;
    fx1 = feval(fun, x2);
    incre2 = norm(fx1);
    iter = iter + 1;
end
% calculo de ACOC
ACOC=log(p(3:end)./p(2:end-1))./log(p(2:end-1)./p(1:end-2));

sol=x1;
incre1=vpa(incre1,6);
incre2=vpa(incre2,6);
ACOC=vpa(ACOC,6);
ACOC=ACOC(:);
sol=vpa(sol,10);
end