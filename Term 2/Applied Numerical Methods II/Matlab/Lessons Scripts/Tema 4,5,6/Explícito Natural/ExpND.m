% Ejemplo
alpha=1;
ci=@(x)1-x;
L=1;
Tmax=1;
nx=10;
nt=2000;
[u,x,t] = CalorExplND (alpha ,ci ,L,nx ,Tmax ,nt)

%%
[x',u(:,(end-1)/4),u(:,(end-1)/2),u(:,3*(end-1)/4) ,u(:,end)]

plot(x,u(:,end))
