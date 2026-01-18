function dy = actividadGrupalFunction_s(x,y)
dy = [y(2);-y(2).^2 + 4.*y(1)+2];
end