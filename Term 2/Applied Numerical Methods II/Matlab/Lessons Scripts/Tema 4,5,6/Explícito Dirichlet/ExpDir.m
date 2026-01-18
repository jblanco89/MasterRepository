%Ejemplo 2 a)
alpha=1;
L=1;
ci=@(x)sin(pi.*x);
h1=@(t)0;
h2=@(t)0;
nx=L/0.1; %nx=(L-0)/h;
%% Estable
Tmax=0.5
nt=(0.5-0)/0.0005; %nt=(Tmax-0)/k=1000
[ua,xa,ta] = CalorExpl (alpha ,ci ,h1 ,h2 ,L,nx ,Tmax ,nt)
%% b) Inestable
Tmax=0.5
nt=(0.5-0)/0.01;
[ub,xb,tb] = CalorExpl (alpha ,ci ,h1 ,h2 ,L,nx ,Tmax ,nt)
%% EXACTA
exacta=@(x,t)exp(-pi^2.*t).*sin(pi.*x);
ex=exacta(xa,0.5)

hold on
plot(xa,ex,'r')
plot(xa,ua(:,end),'*--b')
%%
hold on
plot(xa,ex,'r')
plot(xb,ub(:,end),'*--g')
%%
format long e
[xa',ua(:,end),abs((exacta(xa,0.5))'- ua(:,end)),ub(:,end),abs((exacta(xb,0.5))'- ub(:,end))]