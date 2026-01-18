function [X,Y,iter,incr]= DifnolinearGauss(f,fy ,fz ,a,b,alfa , beta ,N,maxiter ,tol)
% Resuelve una EDO de segundo orden no lineal con condiciones naturales.
% por diferencias finitas. 
% Esta función se basa en el método de Newton con eliminación de Gauss

% [X, Y, iter, incr] = DifnolinearGauss(f, fy, fz, a, b, alfa, beta, N, maxiter, tol)
% 
% Parámetros de entrada:
%   - f: función de la EDO, f(x, y, z)
%   - fy: derivada parcial de f con respecto a y
%   - fz: derivada parcial de f con respecto a z (donde y' = z)
%   - a, b: intervalos en x
%   - alfa, beta: condiciones de frontera en y
%   - N: número de intervalos
%   - maxiter: número máximo de iteraciones
%   - tol: tolerancia para la convergencia
% 
% Parámetros de salida:
%   - X: vector de puntos x
%   - Y: vector de puntos y
%   - iter: número de iteraciones realizadas
%   - incr: incremento máximo en la última iteración
% 
% La función implementa un método iterativo para resolver la EDO no lineal
% usando diferencias centrales para aproximar la derivada y el método de
% Gauss para resolver el sistema resultante.

% Inicialización de variables
% f = f(x,y,z) es la ecuacion diferencial
% fy = fy(x,y,z) es la parcial de f respecto a y
% fz = fz(x,y,z) es la parcial de f respecto a ... y'= z;
h=(b-a)/(N+1); 
k=(beta-alfa)/(N+1);
X=a:h:b; 
Y=alfa:k:beta;
x=X(1:N+2); 
y=Y(1:N+2);

incr=tol +1; % Inicializar parametros
iter =0;
while incr >tol && iter < maxiter
z=(Y(3:N+2)-Y(1:N))/(2* h); %Estimacion de la
%derivada por diferencias centrales ...
z=[((alfa-Y(1))/2) z ((Y(end)+beta)/2)];

fe= feval(f,x,y,z);
fye= feval(fy,x,y,z);
fze= feval(fz,x,y,z);
% fye = fye*ones(N+2,1);


dp=zeros(N+2,1);
dp(2:end-1)=2+(h^2)*fye;
% dp(1)=2*(1-h)+h^2*fye(1)-h^2*fze(1);
dp(1)=(2-h)+h^2*fye(1)-((1/2)*h^2*fze(1));
dp(end)=(2-h)+h^2*fye(end)+((1/2))*h^2*fze(end);
% dp(end)=(2-h)+h^2*fye(end)+h^2/2*fze(end);

ds =-1+h/2*fze(1:end-1);
ds(1)=-2;

di=-1-h/2*fze(2:end);
di(end)=-2;

d(2:N+1)=-(-diff(Y,2)+h^2*fe(2:N+1));
d(1)=-2*y(2)+(2-h)*y(1)+h^2*fe(1)+2*h*alfa;
d(1)=-d(1);
% d(1)=2*(1-h)*y(1)-2*y(2)+h^2*fe(1)+2*h*alfa;
% d(1)=-d(1);


d(N+2)=(2-h)*y(N+2)-2*y(N+1)+h^2*fe(N+2)-3*h;
d(N+2)=-d(N+2);
% d(N+2)= -2*y(N+1)+(2-h)*y(N+2)+h^2*fe(N+2)+beta*h;
% d(N+2)=-d(N+2);

D = diag(dp,0) + diag(ds, 1) + diag(di,-1);

v = D\d';
y=y+v';
Y=y;
incr=max(abs(v));
iter= iter +1;
end
X=X(:);
Y=Y(:);
end