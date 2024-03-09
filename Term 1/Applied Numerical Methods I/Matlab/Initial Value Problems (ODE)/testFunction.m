function [fun, dfun] = testFunction(t,y)
% y=(1-2*t).*y;
% dy = 1-2*t;
fun = y + 1 - (t^2);
dfun = 1;
end