clear;
clc;

f=@(x,y,z) -z.*y;
fy=@(x,y,z) -z.*ones(size(x));
fz=@(x,y,z) -y.*ones(size(x));


% PARAMETROS
% ************************************************
a=1;
b=2;
alfa=0;
beta=-3/2;
N=19;
maxiter=500;
tol=1e-5;

% ************************************************

% LLAMADA A FUNCIÓN
[X,Y,iter,incr, Y_diff] = Difnolin3(f,fy,fz,a,b,alfa,beta,N,maxiter,tol);
% %%
% ex=@(x) sin(x);
% exacta=sin(X);
% Error=abs(exacta-Y);


% GRÁFICA
% ************************************************
hold on
plot(X,Y,'*-r')
% plot(X,exacta,'b')
% ************************************************

% X=vpa(X,6); Y=vpa(Y,6); Error=vpa(Error,6); exacta=vpa(exacta,6);


% RESULTADOS
% ************************************************
results = [X,Y];
t = table(X,Y, 'VariableNames',{'Nodos', 'Solución Y'});
disp(t)
disp('Iteraciones totales:')
disp(iter)
disp('Incremento:')
disp(incr)
% ************************************************
