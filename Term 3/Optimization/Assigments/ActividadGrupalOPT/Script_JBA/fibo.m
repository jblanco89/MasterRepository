function [x_min,f_min,k, prec, lapse] = fibo(fun,dfun,a,b,tol)
    I0 = [a, b];
    % F0 = fun(I0);
    F0 = [fun(I0(1)), fun(I0(2))];
    w1 = I0(2) - I0(1);
    max_iter = 50;
    epsilon = 1e-6;
    
    tic;
    % Fibonacci serie implementation
    F = zeros(max_iter, 1);
    F(1) = 1; F(2) = 1; n = 2;

    while w1 > tol*F(n)
        n = n+1;
        F(n) = F(n-1) + F(n-2);
    
    end

    % Fibonacci search algorithm
    w = w1*(F(n-2)/F(n-1));
    xa = I0(2) - w; 
    xb = I0(1) + w;
    fa = fun(xa); 
    fb = fun(xb);
    k = 0;
    
    for n=2:n-2
        w = w*(F(n-k-1)/F(n-k));
        if fa > fb
            I0(1) = xa; F0(1) = fa; xa = xb; fa = fb; xb = xa + w; fb = fun(xb);
        else
            I0(2) = xb; F0(2) = fb; xb = xa; fb = fa; xa = xb - w; fa = fun(xa); 
        end
    
        width_interval = abs(I0(2) - I0(1));
    
        if width_interval < tol
            break;
        end
    
    end
    x_min = (I0(2) + I0(1))/2;
    f_min = fun(x_min);
    k = n;
    prec = abs(dfun(x_min));
    lapse = toc;
    lapse = vpa(lapse);
end


