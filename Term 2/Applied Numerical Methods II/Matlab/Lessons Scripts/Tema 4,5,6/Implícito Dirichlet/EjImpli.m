%Ejemplo 2 a)
alfa=1;
a=0;
b=1;
ci=@(x)sin(pi*x);
cc1=@(t)0*ones(length(t),1);
cc2=@(t)0*ones(length(t),1);
nx=10;
nt=1000;
Tmax=0.5
[U,x] = CalorImpl(alfa ,ci ,cc1 ,cc2 ,a,b,nx ,Tmax ,nt)
%% EXACTA
exacta=@(x,t)exp(-pi^2.*t).*sin(pi.*x);
ex=exacta(x,0.5)

hold on
plot(x,ex','r')
plot(x,U(:,end),'*--b')
%%
format short e
[x',U(:,end),ex',abs((exacta(x,0.5))'- U(:,end))]