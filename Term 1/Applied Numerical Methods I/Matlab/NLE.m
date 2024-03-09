function [fun,dfun] = NLE(x)
fun = sin(x) - x.^2 + 1;
dfun = cos(x) - 2*x;
end