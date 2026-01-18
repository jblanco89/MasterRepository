clear, clc
digits(100)
format shortG

%PARAMETERS
a = 1;
b = 2;
alfa = 0;
beta = log(2);
n = 10;
tol = 1e-8;
maxiter = 100;

%METHODS

[xs,ys,ts,iter_s]=DisparoSecante('ejercicio1_a_s',a,b,alfa,beta,n,tol,maxiter);
[xn,yn,tn,iter_n]=DisparoNewton('ejercicio1_a_n',a,b,alfa,beta,n,tol,maxiter);

% EXACT SOLUTION
% yex=@(x) 1./(x+3);
yex=@(x) log(x);

% TABLE RESULTS
exacta=yex(xn);
Error_n=abs(yn(:,1)-exacta);
Error_s=abs(ys(:,1)-exacta);
% ['x','Newton','Secante','Exacta','Error Newton', 'Error Secante']
% result = [xn,yn(:,1),ys(:,1),exacta,Error_n, Error_s];

t = table(xn,yn(:,1),ys(:,1),exacta,Error_n, Error_s, ...
    'VariableNames',{'x','Y Newton','Y Secante', ...
    'Y Exacta','Error Newton','Error Secante'});
disp(t)

%% PLOT SOLUTION
hold on
plot(xs,exacta,'r')
plot(xs,ys(:,1),'*b')
plot(xn,yn(:,1),'go')
legend('Exacta', 'Secante', 'Newton')
grid on;



% %% RICHARDSON EXTRAPOLATION
% 
% h1 = 0.1;
% h2 = 0.05;
% h3 = 0.025;
% p = 2;
% 
% n1 = (b-a) / h1;
% n2 = (b-a) / h2;
% n3 = (b-a) / h3;
% 
% [xQ0,yQ0,tn,iter_n]=DisparoNewton('ejercicio5_n',a,b,alfa,beta,n1,tol,maxiter);
% [xQ1,yQ1,tn,iter_n]=DisparoNewton('ejercicio5_n',a,b,alfa,beta,n2,tol,maxiter);
% [xQ2,yQ2,tn,iter_n]=DisparoNewton('ejercicio5_n',a,b,alfa,beta,n3,tol,maxiter);
% 
% Q0 = yQ0(end,1);
% Q1 = yQ1(end,1);
% Q2 = yQ2(end,1);
% 
% H = (h1/h2).^p;
% H2 = (h1/h3).^p;
% 
% Qh1 = Q0 + (Q0 - Q1) / (H - 1);
% Qh2 = Q0 + (Q0 - Q2) / (H2 - 1);
% 
% 
% ErrQ0 = abs(exacta(end,1) - Qh1);
% ErrQ1 = abs(exacta(end,1) - Qh2);
% 
% vpa(Qh1, 10);
% vpa(Qh2, 10);
% vpa(ErrQ0, 10);
% vpa(ErrQ1, 10);
% 
% 
% Qtable = table(Error_n(end,1),ErrQ0,ErrQ1, 'VariableNames', ...
%     {'Error Newton 0.1', 'Error Richardson 0.05', 'Error Richardson 0.025'});
% 
% disp(Qtable)