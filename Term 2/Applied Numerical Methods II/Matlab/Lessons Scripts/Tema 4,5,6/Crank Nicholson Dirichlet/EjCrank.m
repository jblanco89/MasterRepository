%Ejemplo 2 a)
alfa=1;
a=0;
b=1;
ci=@(x)sin(pi*x);
cc1=@(t)0*ones(length(t),1);
cc2=@(t)0*ones(length(t),1);
h=0.1;
k=0.01;
nx=(b-a)/h;
Tmax=0.5
nt=Tmax/k;
[xC,UC]=CrankNicholson(cc2,cc1,ci,a,b,nx,Tmax,nt,alfa)
%% EXACTA
exacta=@(x,t)exp(-pi^2.*t).*sin(pi.*x);
ex=exacta(xC,0.5)

hold on
plot(xC,ex,'r')
plot(xC,UC(:,end),'o-m')
%%
format short e
[xC',UC(:,end),ex',abs((exacta(xC,0.5))'- UC(:,end))]