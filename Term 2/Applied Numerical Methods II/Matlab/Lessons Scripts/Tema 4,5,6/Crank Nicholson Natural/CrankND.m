% Ejemplo
alpha=1;
ci=@(x)1-x;
L=1;
a=0;
b=1;
Tmax=1;
nx=10;
nt=2000;
[xC,UC]=CrankNicholsonND(ci,a,b,nx,Tmax,nt,alpha)

%%
[xC',UC(:,(end-1)/4),UC(:,(end-1)/2),UC(:,3*(end-1)/4) ,UC(:,end)]
hold on
plot(xC,UC(:,end),'*--b')