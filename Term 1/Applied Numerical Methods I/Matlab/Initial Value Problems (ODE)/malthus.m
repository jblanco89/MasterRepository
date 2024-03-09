function [fun, dfun] = malthus(t, y)
K = -10;
%ya = 1
%a = 0;
%b = 3;
fun = K*y;
dfun = K;
end