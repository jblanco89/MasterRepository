function dy = ejercicio1_b_n(x,y)
dy = [y(2); y(1).^3 - y(1).*y(2); y(4);3.*y(1).^2 - y(2).*y(3) - y(1).*y(3)];
end