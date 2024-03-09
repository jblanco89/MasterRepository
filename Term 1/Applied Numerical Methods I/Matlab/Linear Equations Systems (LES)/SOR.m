function x = SOR(A,b,x0,w,Tol,maxiter)
k=0;
Norm = 1;
D = diag(diag(A));
L = tril(A)- D;
U = triu(A)- D;
% check for convergence condition for Gauss-Seidel method
e= max(eig(inv(D+w*L)*(D*(1-w) -w*U)));
if abs(e) >= 1
    disp ('Dado que el módulo del valor propio más grande de la matriz iterativa no es menor que 1.') 
    disp ('El proceso no será convergente')
    return
end
x(:,1) = x0;
while Norm > Tol && k < maxiter
    k = k + 1;
    x(:,k+1) = (D+w*L)\((- w*U + D*(1-w))*x(:,k) + w*b');% SOR formula
    Norm = norm(x(:,k+1) - x(:,k), inf);
    fprintf('%2d', k);
    fprintf('%10.6f', x(:,k)');
    fprintf('%10.6f\n', Norm);
    if k > maxiter
        disp('Convergencia no alcanzada con el máximo número de iteraciones')
       break
    end
end
x = x(:,k)';

end