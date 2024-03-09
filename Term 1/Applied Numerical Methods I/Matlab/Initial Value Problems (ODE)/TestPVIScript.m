format shortG
clear,clc

%%
% Lecture 7 Test
% *** Question 9***
% y'(t) = y(t) - t^2 + 1, t=[0,2], y(0) = 0.5, h = 1
% select iterative method whose solution is y(t) = [0.5, 2.25, 4.875]
fun = "testFunction";
a = 0;
b = 2;
N = 3;
ya = 0.5;
Tol = 1e-8;
maxiter = 30;

[t_AB2, y_AB2] = AdamsBashforth2(fun, a,b,N,ya);
[t_AB4, y_AB4] = AdamsBashforth4(fun, a,b,N,ya);
[t_AM2, y_AM2] = AdamsMoulton2(fun, a,b,N,ya, Tol, maxiter);
[t_AM4, y_AM4] = AdamsMoulton4(fun, a,b,N,ya, Tol, maxiter);
[t_ABM2, y_ABM2] = ABM2(fun,a,b,N,ya);
[t_RK, y_RK] = RK4(fun, a,b,N,ya);
[t_Eu, y_Eu] = Euler(fun, a,b,N,ya);
[t_EuI, y_EuI] = EulerImplicit(fun, a,b,N,ya,Tol, maxiter);

t1 = table(y_AB2, y_AB4, y_AM2, y_AM4,y_ABM2, y_RK, y_Eu, ...
    y_EuI, 'VariableNames',{'Adams-Bashforth 2', ...
    'Adams-Bashforth 4', 'Adams-Moulton 2', 'Adams-Moulton 4', ...
    'Predictor-Corrector', 'Runge-Kutta 4', 'Euler', 'Euler Implícito'});
disp(t1)
hold on;
plot(t_AB2, y_AB2, t_AB2, y_AM2, ...
    t_AB2, y_AM4, t_ABM2, y_ABM2, ...
    t_RK,y_RK, t_Eu, y_Eu, t_EuI, y_EuI)
legend('Adams-Bashforth 2', ...
    'Adams-Moulton 2', 'Adams-Moulton 4', 'Predictor-Corrector', ...
    'Runge-Kutta 4','Euler', 'Euler Implícito')
grid on;
% hold off;


%%
% Lecture 7; Example 01 - pg. 18 of slides
% we are going to use testFunctionSystem.m file
[t_AB2s, y_AB2s] = AdamsBashforth2Systems("testFunctionSystem", 0, 2, 10, [0 0.25]);
[t_AB4s, y_AB4s] = AdamsBashforth4Systems("testFunctionSystem", 0, 2, 10, [0 0.25]);
[t_AM2s, y_AM2s] = AdamsMoulton2Systems("testFunctionSystem", 0, 2, 10, [0 0.25]);
[t_AM4s, y_AM4s] = AdamsMoulton4Systems("testFunctionSystem", 0, 2, 10, [0 0.25]);
[t_Heuns, y_Heuns] = Heun_Systems("testFunctionSystem", 0, 2, 10, [0 0.25]);
[t_ABM2s, y_ABM2s] = ABM2Systems("testFunctionSystem", 0, 2, 10, [0 0.25]);
[t_RKs, y_RKs] = RK4Systems("testFunctionSystem", 0, 2, 10, [0 0.25]);

t2 = table(y_AB2s(:,1), y_AB4s(:,1), y_AM2s(:,1), y_AM4s(:,1), y_Heuns(:,1), ...
    y_ABM2s(:,1), y_RKs(:,1) ,'VariableNames',{'Adams-Bashforth 2', ...
    'Adams-Bashforth 4', 'Adams-Moulton 2', 'Adams-Moulton 4', 'Heun', ...
    'Predictor-Corrector', 'Runge-Kutta 4'});
disp(t2)
% plot(t_AB2s, y_AB2s(:,1), t_AB2s, y_AB4s(:,1), t_AB2s, y_AM2s(:,1), ...
%     t_AB2s, y_AM4s(:,1), t_Heuns, y_Heuns(:,1), t_ABM2s, y_ABM2s(:,1), ...
%     t_RKs,y_RKs(:,1))
% legend('Adams-Bashforth 2', 'Adams-Bashforth 4', ...
%     'Adams-Moulton 2', 'Adams-Moulton 4', 'Heun', 'Predictor-Corrector', ...
%     'Runge-Kutta 4')
% grid on;
% We test PVI method using flame propagation problem 
% this problem is well-known as a type of stiff ODE. 

