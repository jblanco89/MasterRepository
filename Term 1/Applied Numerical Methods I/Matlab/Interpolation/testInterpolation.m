format shortG;
clear, clc
% x = [0 0.5 1];
% y = [1 0.9385 .7652];
% dy = [0 .2423 .4401];

load('input_data.mat');
load('dy.mat');
x = data(1,:);
y = data(2,:);
t = [0, 3, 5, 8, 10, 12, 13, 15, 18, 20, 22, 25]; %value points
% t = [10, 12, 15];
% dy = dBdt(1:end);
xi = x(1:end-1);
fi = y(1:end-1);
% testing interpolation methods
% [P_Hermite, err_Hermite, M] = hermite_interpolation(x, y, dy, length(x), t);

% t = 2005;
% xi =1971:10:2011;
% fi =[33.956 37.743 39.434 40.847 46.816];


[P_Newton, err_Newton] = newton_interpolation(xi, fi, length(xi), t);
[P_Lagrange, err_Lagrange] = lagrange_interpolation(xi, fi, length(xi) - 1, t);
[ai, bi, ci, di, p] = CubicSplines(xi, fi);
syms x
P_Spline = double(subs(p(4),x,t));
P_SplinesMatlab = spline(xi, fi, t);


% Result = [P_Hermite', P_Newton', P_Lagrange'];
Result = table(t',P_SplinesMatlab', P_Spline', P_Newton', P_Lagrange', ...
    'VariableNames',{'Points','Matlab','Splines', 'Newton', 'Lagrange'});
disp(Result);
plot(xi, fi,'-+',t,P_SplinesMatlab,'-o',t, P_Spline, '-*', t, P_Lagrange, '-^');
legend('Data Points', 'Matlab','Cubic Splines', 'Lagrange Interpolation')
title('Comparison of Interpolation Methods')
xlabel('x')
ylabel('y')
grid on;