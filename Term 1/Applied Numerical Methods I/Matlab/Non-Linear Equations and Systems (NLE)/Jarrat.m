function [sol,iter,ACOC,incre1,incre2] = Jarrat(fun,x0,tol,maxiter)
digits(200)
x0=x0(:);
iter=0;
[fx0,dfx0, d2fx0]=feval(fun,x0);
incre1=tol+1;
incre2=tol+1;
% p=[];
p = incre1;
while incre2>tol && incre1>tol && iter<maxiter
    
    %Linea Método de Jarrat
    y0 = x0 - (2/3)*(fx0/dfx0);
    [fy0, dfy0, d2fy0] = feval(fun, y0);
    x1=x0-((3*dfy0 + dfx0)/(6*dfy0 - 2*dfx0))*(fx0/dfx0);
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
ACOC=vpa(ACOC,6);
ACOC=ACOC(:);
sol=vpa(sol,10);
end