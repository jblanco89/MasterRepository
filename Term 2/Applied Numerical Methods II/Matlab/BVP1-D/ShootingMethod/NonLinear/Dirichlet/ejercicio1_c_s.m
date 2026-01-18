function dy = ejercicio1_c_s(x,y)
dy = [y(2); 2.*y(1).^3 - 6.*y(1) - 2.*x.^3];
end