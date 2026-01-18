function dy = ejercicio1_c_n(x,y)
dy = [y(2); 2.*y(1).^3 - 6.*y(1) - 2.*x.^3;y(4); 6.*y(3).*y(1).^2 - 6.*y(3)];
end