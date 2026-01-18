% clear, clc
digits(100)
format shortG

f=@(x,y,z)-z.^2 + 4.*y + 2;
fy=@(x,y,z)4;
fz=@(x,y,z)-2.*z;
a=1;
b=2;
alfa=0;
beta=-3;
N=100;
h=(b-a)/(N+1);
maxiter=50;
tol=1e-5;
i = 0:1:N+1;
[X_G,Y_G,iter_G,incr_G]=DifnolinearGauss(f,fy,fz,a,b,alfa,beta,N,maxiter,tol);
[X_C,Y_C,iter_C,incr_C]=DifCroutPablo(f,fy,fz, a, b, alfa, beta, N, maxiter, tol);

t = table(i', X_G, Y_G, Y_C);
t.Properties.VariableNames = {'I', 'X', 'Y Gauss', 'Y Crout'};
t.ERROR = abs((t.("Y Gauss") - t.("Y Crout")));
disp(t)

figure
hold on
plot(X_G,Y_G,'-g', X_C,Y_C, '-r')
xlabel('Eje X');
ylabel('Eje Y');
grid on
title('Problema de contorno');
legend('Gauss', 'Crout')