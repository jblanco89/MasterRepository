function dy = actividad1_n(x,y)
dy = [y(2);2.*x.*y(1)-(y(2).*y(1))+2;y(4);(-y(1).*y(4)) + 2.*x.*y(3) - (y(3).*y(2))];
end