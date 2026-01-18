function [X_res, Y_res, iter, incre] = DifCroutPablo(f,fy,fz, a, b, alfa, beta, N, maxiter, tol)

    h = (b-a)/(N+1);
    k = (beta-alfa)/(N+1);
    X = a:h:b;
    Y = alfa:k:beta;


    incre =tol +1;
    iter = 0;
    
    while incre>tol && iter<maxiter
       Z = (Y(3:N+2)-Y(1:N))/(2* h); 
       Z=[(alfa-Y(1))/2 Z (Y(end)+beta)/2];


       f_eval = feval(f,X,Y,Z);
       fy_eval = feval(fy, X, Y,Z);
       fz_eval = feval(fz, X, Y, Z);
   
       dp=zeros(N+2,1);
       dp(2:end-1) = 2 + h^2*fy_eval;
       dp(1) = 2-h+h^2*fy_eval(1) - (1/2)*h^2*fz_eval(1);
       dp(end) = 2-h+h^2*fy_eval(1) + (1/2)*h^2*fz_eval(1);

       ds = -1+h/2*fz_eval(1:end-1);
       ds(1) = -2;

       di = -1-h/2*fz_eval(2:end);
       di(end) = -2;
 
       d(2:N+1) = -(-diff(Y,2) + h^2*f_eval(2:N+1));
       d(1) = -((2-h) * Y(1) + h^2*f_eval(1) -2*Y(2));
       d(N+2) = -((2-h) * Y(N+2) - 2*Y(N+1) + h^2*f_eval(N+2) - 3*h);
       
       v_c = Crout(dp',ds,di,d);

       Y = Y+v_c';
       incre = norm(v_c);
       iter = iter+1;
       if iter == maxiter
            disp('No se ha encontrado solución')
       end
       Y_res = Y';
       X_res = X';
       
    end
end
       

