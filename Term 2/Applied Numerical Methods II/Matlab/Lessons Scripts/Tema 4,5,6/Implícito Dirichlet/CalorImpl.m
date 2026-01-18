function [U,x] = CalorImpl(alfa ,ci ,cc1 ,cc2 ,a,b,nx ,Tmax ,nt)
h=(b-a)/nx;     x=a:h:b;
k=Tmax/nt;      t=0:k:Tmax;

U=zeros(nx+1,nt+1);
U(1,:)=feval(cc1,t);
U(nx+1,:)=feval(cc2,t);

U(:,1)=feval(ci,x);

u(1,1) =(ci(1)+cc1(1))/2;  
u(end,1) =(cc2(1)+ci(end))/2;  

lambda=alfa^2*k/h^2;

dp=(1+2*lambda)*ones(nx-1,1);
ds=-lambda*ones(nx-2,1);
di=ds;  

for j=2:nt+1
    d=U(2:nx,j-1); 

    d(1)=d(1)+lambda*U(1,j-1);
    d(end)=d(end)+lambda*U(nx+1,j-1);   
    
    z=Crout(dp,ds,di,d);
    U(2:nx,j)=z;
end
end