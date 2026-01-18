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
maxiter=50;
tol=1e-5;
i = 0:1:N+1;
[X_G,Y_G,iter_G,incr_G]=DifnolinearGauss(f,fy,fz,a,b,alfa,beta,N,maxiter,tol);
[X_C,Y_C,iter_C,incr_C]=DifnolinearCrout(f,fy,fz,a,b,alfa,beta,N,maxiter,tol);

t = table(i(end-9:end)', X_G(end-9:end), Y_G(end-9:end), Y_C(end-9:end));
t.Properties.VariableNames = {'I', 'X', 'Y Gauss', 'Y Crout'};
t.ERROR = abs((t.("Y Gauss") - t.("Y Crout")));
disp(t)

figure
hold on
% plot(X_G,Y_G,'-g', X_C,Y_C, '-r')
plot(X_G,Y_G,'-b', 'LineWidth', 1.0)
xlabel('Eje X intervalo [a, b]');
ylabel('Valores de Y');
grid on
title('Problema de contorno con Eliminación de Gauss');
% legend('Gauss', 'Crout')
