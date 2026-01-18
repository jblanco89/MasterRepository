% Ejemplo
alpha=1;
ci=@(x)1-x;
L=1;
a=0;
b=1;
Tmax=1;
nx=10;
nt=2000;
[U,x] = CalorImplND(alpha ,ci,0,1,nx ,Tmax ,nt)

%%
[x',U(:,(end-1)/4),U(:,(end-1)/2),U(:,3*(end-1)/4) ,U(:,end)]

%%
hold on
plot(x,U(:,end),'r')