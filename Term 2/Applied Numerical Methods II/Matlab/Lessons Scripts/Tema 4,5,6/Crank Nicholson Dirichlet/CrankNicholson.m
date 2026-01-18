function [x,U]=CrankNicholson(cc2,cc1,ci,a,b,nx,Tmax,nt,alfa)

h=(b-a)/nx;
x=a:h:b;
k=Tmax/nt;
t=0:k:Tmax;

U=zeros(nx+1,nt+1);
U(1,:)=feval(cc1,t);
U(nx+1,:)=feval(cc2,t);
U(:,1)=feval(ci,x);

lambda=k*alfa^2/h^2;

dp=(1+lambda)*ones(nx-1,1);
ds=(-lambda/2)*ones(nx-2,1);
di=ds;


bp=(1-lambda)*ones(nx-1,1);
bs=(lambda/2)*ones(nx-2,1);
bi=bs;

B=diag(bp)+diag(bs,1)+diag(bi,-1);

for j=1:nt
    d=B*U(2:nx,j);
    z=Crout(dp,ds,di,d);
    U(2:nx,j+1)=z;
end
end
