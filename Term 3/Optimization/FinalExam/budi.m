function [x_min,f_min,k_opt,prec, lapse] = budi(fun,dfun, a, b, tol)
epsilon = 1e-6;
max_iter = 50;
k = 0;

tic;
x1 = (a+b)/2;
xa = x1-(epsilon/2);
xb = x1 + (epsilon/2);

while abs(b-a) >= tol && k < max_iter

if fun(xa) <= fun(xb)
    b=xb; a=a;
elseif fun(xa) > fun(xb)
    a=xa; b=b;
end

x1 = (a + b) /2;
xa = x1-(epsilon/2);
xb = x1 + (epsilon/2);
fl = fun(a);
fr=fun(b);

k = k+1;

end
x_min = (a+b)/2;
f_min = fun(x_min);
k_opt = k;
prec = abs(dfun(x_min));
lapse = toc;
lapse = vpa(lapse);
end