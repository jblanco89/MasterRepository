function [x_min,f_min,k, prec, lapse] = new(fun,dfun,hess,x0,tol)
    max_iter = 50;
    x = x0;
    k = 0;
    tic;
    while k < max_iter
        D = dfun(x);
        n1 = norm(D);
        h = hess(x)\D;
        n2 = norm(hess(x));
        x = x - h';
        if (n1 < tol) || (n2 < tol)
        break
        end
        k = k+1;
    end
    x_min = x;
    f_min = fun(x_min);
    k = k;
    prec = abs(dfun(x_min));
    lapse = toc;
    lapse = vpa(lapse);
end