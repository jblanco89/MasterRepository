function dy = ejercicio1_a_n(x,y)
dy = [y(2);-y(2).^2 - y(1) + log(x); y(4); -y(3) - 2.*y(2).*y(4)];
end