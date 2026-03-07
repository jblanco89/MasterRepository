function [x_opt, f_opt, lambda_opt, k_opt, prec_opt, lapse] = rein(fun, grad, x0)
    % Line search with Wolfe conditions (Armijo & Goldstein conditions)
    
    % Constants
    max_iter = 50;
    alpha = 1e-4;  % Armijo parameter (first condition)
    beta = 1e-1;  % Goldstein parameter (second condition)
    epsilon = 1e-6;
    
    % Initial values
    x_k = x0;
    k = 0;
    tic;
    while k < max_iter
        %k = k+1;
        g_k = grad(x_k);
        if norm(grad(x_k)) > epsilon
            d_k = -g_k;
            lambda = 1;
            while true
                x_k2 = x_k + lambda * d_k;

                % Armijo
                lower_bound =  fun(x_k) + alpha * lambda * g_k' * d_k;
                % Goldstein
                upper_bound = grad(x_k2)'*d_k;
                
                if fun(x_k2) <= lower_bound && upper_bound <= beta * g_k' * d_k 
                    break
                end
                lambda = lambda / 2;  
            end
            x_k = x_k + lambda * d_k;
            k = k + 1;
        else
            break;
        end
    end
    x_opt = x_k;
    f_opt = fun(x_opt);
    lambda_opt = lambda;
    k_opt = k;
    prec_opt = abs(grad(x_opt));
    lapse = toc;
    lapse = vpa(lapse);
end