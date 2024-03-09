% Condiciones Iniciales y Parámetros
G = 6.672*10^-11; M = 5.9742*10^24; R = 6378140; H = 772000; v_0 = 6700; h=1;
t = transpose(0:h:1200);
r_0 = R + H; theta_0 = 0; r_0_dot = 0; theta_0_dot = v_0/r_0;
% Posibles sistemas de Ecuaciones
% Descripcción original, variables ()
%r_dot2 = r*theta_dot^2 - G*M/r^2;
%theta_dot2 = - 2*r_dot*theta_dot/r;
y_1_0 = r_0; y_2_0 = r_0_dot; y_3_0 = theta_0; y_4_0 = theta_0_dot;
y_1 = zeros(1200,1); y_2 = zeros(1200,1); y_3 = zeros(1200,1); y_4 = zeros(1200,1);
y_1(1) = y_1_0; y_2(1) = y_2_0; y_3(1) = y_3_0; y_4(1) = y_4_0;
%y_1_dot  = y_2; y_2_dot = y_1*y_2^2 - G*M/y_1^2; y_3_dot = y_4; y_4_dot = -2*y_2*y_4/y_1;
solved = Euler_Method(y_1,y_2,y_3,y_4, 1, 1200);
%%%Representamos las trayectorias
figure(1)
plot(t,solved(:,1))
g = title('Espacio (m) en función del tiempo (s)');
set(g, 'Interpreter', 'latex')
%%%Representamos los espacios de fase
figure(2)
tiledlayout(1,2)
nexttile
plot(solved(:,1), solved(:,2))
t = title('Espacio de fases para r vs $\dot{r}$');
set(t, 'Interpreter', 'latex')
nexttile
plot(solved(:,3),solved(:,4))
h = title('Espacio de fase para $\theta$ vs $\dot{\theta}$');
set(h, 'Interpreter', 'latex')
velocidad_de_impacto = crash(solved(:,1), solved(:,2));
function result = Euler_Method(y_1,y_2,y_3,y_4, h, finish)
    G = 6.672*10^-11; M = 5.9742*10^24;
    for i = 2:h:finish+1
        y_1(i) = y_1(i-1) + h*y_2(i-1);
        y_2(i) = y_2(i-1) + h*(y_1(i-1)*(y_4(i-1))^2 - G*M/(y_1(i-1))^2);
        y_3(i) = y_3(i-1) + h*(y_4(i-1));
        y_4(i) = y_4(i-1) - h*(2*y_2(i-1)*y_4(i-1)/y_1(i-1));
      result = [y_1, y_2,y_3,y_4];
    end
end
function result = crash(r, r_dot)
    R = 6378140;
    for i = 2:1:length(r)
        if r(i-1) > R && r(i)<R
            result = r_dot(i-1) + (r_dot(i) - r_dot(i-1))/(r(i) - r(i-1))*(R - r(i-1));
        end
    end
end
