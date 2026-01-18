clc;
clear;
digits(6)
% PARAMETERS
% **************************
a=1;
b=2;
alfa=2;
beta=(log(2)+1)/2;
N = 19;
tol=10e-05;
maxiter=500;
% ***************************


[xi ,yi ,t,iter,incre] = DisparoSecanteP1('final_exam_s',a,b,alfa,beta,N, tol, maxiter);

% RESULTS
format long
r = [xi,yi(:,1)];
disp(r)
disp(iter)
disp(incre)


% PLOT
hold on
plot(xi,yi(:,1),'*-b')

