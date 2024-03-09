function [ci, xi] = Coeficientes_Nodos_Gauss_Laguerre(n)
%Con esta función de Matlab vamos a calcular los coeficientes y los 
%nodos para la cuadratura de Gauss-Laguerre. 

%La variable de entrada es el número de nodos que queremos calcular: n

%Las variables de salida son los coeficientes (ci) y los nodos (xi)

%--------------------------------------------------
%Definimos la variable simbólica
syms x 

%Definimos los dos primeros polinomios: (como matlab no puede trabajar 
%                                        con la componente 0 empezamos con
%                                        la primera)
pk=1;
pk1=(1-x);
%(En los polinomios añadimos una unidad porque matlab no trabaja con la componente 0)
%Calculamos un polinomio más que el número de nodos para los coeficientes
for k=0:n-1
    pk2=simplify((2*k+3-x).*pk1-(k+1)^2.*pk);
    pk=pk1;
    pk1=pk2;
end
xi=double(solve(pk==0));
ci=(factorial(n)^2.*xi)./(double(subs(pk1,x,xi)).^2);
end

