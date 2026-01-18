function [U,x] = CalorImplND(alfa ,ci ,a,b,nx ,Tmax ,nt)
h=(b-a)/nx;     x=a:h:b;
k=Tmax/nt;      t=0:k:Tmax;

U=zeros(nx+1,nt+1);
cix= feval (ci ,x);
U(: ,1)=cix';

lambda=alfa^2*k/h^2;

dp=(1+2*lambda)*ones(nx+1,1);
ds=-lambda*ones(nx,1);
di=ds; 

dp(end)=(1+2*lambda+2*lambda*h);
ds(1)=-2*lambda;
di(end)=-2*lambda;
 
for j=2:nt+1
    d=U(1:nx+1,j-1)+k*exp(-t(j))*ones(nx+1,1);
    z=Crout(dp,ds,di,d);
    U(1:nx+1,j)=z;
end
end