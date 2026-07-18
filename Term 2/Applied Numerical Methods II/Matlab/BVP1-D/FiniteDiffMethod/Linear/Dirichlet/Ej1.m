%Ejemplo 1
p=@(x)-2./x;
q=@(x)2./(x.^2);
r=@(x)(sin(log(x)))./(x.^2);
a=1;
b=2;
alfa=1;
beta=2;
N=9;
[xi,yi] = DiFiLinealEj1 (p,q,r,a,b,alfa ,beta ,N)
%%
c1=1/70*(69+4*cos(log(2))+12*sin(log(2)));
c2=4/35-2/35*cos(log(2))-6/35*sin(log(2));
ex=@(x)c1.*x+c2./(x.^2)-3/10.*sin(log(x))-1/10.*cos(log(x));

exacta=ex(xi);
Error=vpa(abs(exacta-yi),6);
[vpa(xi),vpa(yi),vpa(exacta),vpa(Error)]
%%
hold on
plot(xi,yi,'*-r')
plot(xi,exacta,'b')