format short
clear, clc
x = linspace(1,5,12);
f = @(x) log(x);
n = length(x);
h = (x(end) - x(1)) / n;
% h = 0.1;
N1 = zeros(size(n));
N2 = zeros(size(n));
for i=1:n
    [N1(i),N2(i)]  = richardson(f, x(i), h, 6);
end
der = 1 ./ x;
e1 = abs((der - N1));
e2 = abs((der - N2));
t = table(x', N1', N2', der', e1', e2', 'VariableNames', {'x','Diff Finita', 'Richardson', 'Real', 'Error N1', 'Error N2'});
disp(t)
plot(x,N1, '-', x, N2, '-', x, der, '-')
legend('Diff Finita', 'Richardson', 'Real');