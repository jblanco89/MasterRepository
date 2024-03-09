function [fun,dfun] = testFunction2(t,y)
fun = 1 + (1/t).*y + (1/t^2).*y^2;
dfun = (2*y)/t^2 + 1/t;
end