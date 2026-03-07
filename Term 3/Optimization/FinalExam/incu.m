function [x_min,f_min,k_opt,prec, lapse] = incu(fun,dfun, a, b, tol)
    max_iter = 50;
    epsilon = 1e-6;
    k = 0;
    tic;
    while abs(b - a) > tol && k < max_iter
        c = (a + b)/2;
        x= [a c b];
        A = vander(x);
        y= [fun(a); fun(c); fun(b)];
        p = A\y;
        x_opt = -p(2)/(2*p(1));
    
        if x(2) < x_opt
            a = x(2);
            b = x(3);
            c = (a + b)/2;
            x= [a c b];
            A = vander(x);
            y= [fun(x(1)); fun(x(2)); fun(x(3))];
            p = A\y;
            x_opt = -p(2)/(2*p(1));
        elseif x(2) > x_opt
            a = x(1);
            b = x(2);
            c = (a + b)/2;
            x= [a c b];
            A = vander(x);
            y= [fun(x(1)); fun(x(2)); fun(x(3))];
            p = A\y;
            x_opt = -p(2)/(2*p(1));
        end
            k = k+1;
    end
    x_min = x_opt;
    f_min = fun(x_min);
    k_opt = k; 
    prec = abs(b - a);
    lapse = toc;
    lapse = vpa(lapse);
end