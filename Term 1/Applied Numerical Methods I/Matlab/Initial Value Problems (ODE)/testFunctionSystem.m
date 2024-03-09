function [fun, dfun] = testFunctionSystem(t,theta)
t1=theta(1); t2=theta(2);
fun=[t2; 9.81/0.5*t1];
dfun = [1; 19.62]; %9.81/0.5
end