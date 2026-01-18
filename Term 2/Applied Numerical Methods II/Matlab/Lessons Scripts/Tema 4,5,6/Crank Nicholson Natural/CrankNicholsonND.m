function [x,U]=CrankNicholsonND(ci,a,b,nx,Tmax,nt,alfa)

h=(b-a)/nx;
x=a:h:b;
k=Tmax/nt;
t=0:k:Tmax;
U=zeros(nx+1,nt+1);
U(:,1)=feval(ci,x);

lambda=k*alfa^2/h^2;

dp=(1+lambda)*ones(nx+1,1);
dp(end)=1+lambda+lambda*h;

ds=(-lambda/2)*ones(nx,1);
di=ds;
ds(1)=-lambda;
di(end)=-lambda;

bp=(1-lambda)*ones(nx+1,1);
bp(end)=1-lambda-lambda*h;
bs=(lambda/2)*ones(nx,1);
bi=bs;
bs(1)=lambda;
bi(end)=lambda;

B=diag(bp)+diag(bs,1)+diag(bi,-1);

for j=1:nt
    d=B*U(1:nx+1,j);
    d=d+k/2*(exp(-t(j))+exp(-t(j+1)))*ones(nx+1,1);
    z=Crout(dp,ds,di,d);
    U(1:nx+1,j+1)=z;
end
end
