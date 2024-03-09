function [fun,dfun] = testFunctionNLSE(X)
    x=X(1); y=X(2);
    fun=[exp(x).*exp(y)+x.*cos(y); x+y-1];
    dfun=[exp(x).*exp(y)+cos(y), ...
    exp(x).*exp(y)-x.*sin(y); 1, 1];
end