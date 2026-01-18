clear, clc
digits(100)
format shortG

%PARAMETERS
a = 1;
b = 2;
alfa = 2;
beta = 1;
n = 19;
tol = 1e-5;
maxiter = 100;

%METHODS

% [xs,ys,ts,iter_s]=DisparoSecante('final_exam_s',a,b,alfa,beta,n,tol,maxiter);
[xn,yn,tn,iter_n]=DisparoNewton('final_exam_n',a,b,alfa,beta,n,tol,maxiter);

% EXACT SOLUTION
% yex=@(x) 1./(x+3);
% yex=@(x) log(x);

% TABLE RESULTS
% exacta=yex(xn);
% Error_n=abs(yn(:,1)-exacta);
% Error_s=abs(ys(:,1)-exacta);
% ['x','Newton','Secante','Exacta','Error Newton', 'Error Secante']
% result = [xn,yn(:,1),ys(:,1),exacta,Error_n, Error_s];

% t = table(xn,yn(:,1),ys(:,1),exacta,Error_n, Error_s, ...
%     'VariableNames',{'x','Y Newton','Y Secante', ...
%     'Y Exacta','Error Newton','Error Secante'});


% results = ['x','Newton','Secante'];
% result = [xn,yn(:,1),ys(:,1)];
% 
% t = table(xn,yn(:,1),ys(:,1),'VariableNames',{'x','Y Newton','Y Secante'});

results = ['x','Newton'];
result = [xn,yn(:,1)];

t = table(xn,yn(:,1), yn(:,2),'VariableNames',{'x','Y Newton', 'dY Newton'});
disp(t)
disp(tn)
disp(iter_n)

%% PLOT SOLUTION
hold on
plot(xn,yn(:,1),'-b', xn,yn(:,2),'-r')
legend('Y Newton', 'dY Newton')
grid on;


%% PLOT SOLUTION
% hold on
% plot(xs,exacta,'r')
% plot(xs,ys(:,1),'*b')
% plot(xs,yn(:,1),'go')
% legend('Exacta', 'Secante', 'Newton')
% grid on;

