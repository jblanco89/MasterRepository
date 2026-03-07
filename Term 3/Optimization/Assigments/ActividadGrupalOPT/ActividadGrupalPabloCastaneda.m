
%% Buscamos encontrar la curva que uno los puntos (0,1) y (1,1) que minimice el area de revolución
syms x;
f = @(x) x^2-x+1;

x_i = 0:0.01:1;
n=length(x_i);
u_i = ones(1,n);
du_i = zeros(1,n);
T_i = zeros(1,n);
g = zeros(1,n);
maxiter = 100000;
tol = 10^-8;
t=0.005;
iter = 0;
error = 100;

while error>tol && iter<maxiter
    disp(iter)
    for i=1:n-1
        %% Aproximamos la derivada
        du_i(i) = (u_i(i+1)-u_i(i))*n;
        %% Calculamos el valor de la integral
        T_i(i) = (u_i(i+1)+u_i(i))*sqrt(1+n^2*(u_i(i+1)-u_i(i))^2)/(2*n);
        
        %% Calculamos el gradiente
            
        if i>1 
            u_r1 = sqrt(1+n^2*(u_i(i)-u_i(i-1))^2);
            u_r2 = sqrt(1+n^2*(u_i(i+1)-u_i(i))^2);
            g(i) = (1/(2*n)) * ((u_r1 + n^2*(u_i(i)^2-u_i(i-1)^2)/u_r1) + (u_r2 - n^2*(u_i(i+1)^2-u_i(i)^2)/u_r2));
            
        end

    
    end
    u_r = sqrt(1+n^2*(u_i(2)-u_i(1))^2);
    g(1) = (1/(2*n)) * (u_r - n^2*(u_i(2)^2-u_i(1)^2)/u_r);
    u_rn = sqrt(1+n^2*(u_i(end)-u_i(end-1))^2);
    g(n) = (1/(2*n))* (u_rn + n^2*(u_i(end)^2 - u_i(end-1)^2)/u_rn);
    
    u_i_new = u_i - t*g;
    u_i_new(1)=1;
    u_i_new(end)=1;

    error = norm(u_i-u_i_new);
    res = sum(T_i);
    disp(res)
    disp(vpa(error))
    iter = iter+1;
    
    u_i = u_i_new;

    if norm(g) < tol
        break
    end

end
disp(iter)
u_analytic = @(x) 0.8484*cosh((x-1/2)/0.8484);
u_i_analytic = subs(u_analytic,x_i);
figure
plot(x_i, u_i, 'b')
hold on
plot(x_i, u_i_analytic, 'r')

