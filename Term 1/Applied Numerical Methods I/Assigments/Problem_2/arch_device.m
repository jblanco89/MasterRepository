clear, clc
format longG;
% Archery Device - Activity 2
%
%
% This is the matlab script with Activity 2 solution. 
%
% Numerical Methods 1, UNIR
% Computer Science and Mathematical Engineering Master
%
% Author / Dev: Javier Blanco
% Date: Dic-19-2023
% Cohort: 2023-2024
% 
% 
% Matlab Version: R2021b 

%Archery Device Data
archery_data = readtable("arch_data.dat");


x = archery_data.x_m_;
F = archery_data.F_N_;


% PART 1: INTEGRATION EVALUATION
% Hints: Trapezium and Simpson's rules are needed
% Calling trapezium and simpsons matlab functions
% Bonus: Simpson 3/8 and MidPoint methods have been developed and used too 

% INPUT PARAMETERS

a = x(1); % lower bound 
b = x(end); % upper bound
n = length(x); % number of intervals
F = F'; % function of x (f(x))
m = 0.075; % mass in kg of arrow

% CALLING NUMERIC METHODS:
% Trapezium
% Simpson 1/3
% Simpson 3/8
% MidPoint

[I_trapezium, err_trapezium] = trapezium(0, F, a, b, n);
[I_simpson, err_simpson] = simpson13(0, F, a, b, 10); %n even
[I_simpson38, err_simpson38] = simpson38(0, F, a, b, 9); %n multiple of 3
[I_mp, err_mp] = midpoint(0, F, a, b, n);

Ints = [I_trapezium, I_simpson, I_simpson38, I_mp];
Errs = [err_trapezium, err_simpson, err_simpson38, err_mp];

part_1_table = table(Ints', Errs','RowNames', ...
    {'Trapezium', 'Simpson 1/3', 'Simpson 3/8', 'MidPoint'}, ...
    'VariableNames',{'Integral', 'Error'});
disp(part_1_table);

% PART 2: DISCUSSION:
% (In Spanish because of asignation requirement)
% El método de Simpson 1/3 exhibe una mayor precisión en comparación 
% con el método del trapecio y el método de Simpson 3/8. 
% Además, el método del punto medio demuestra un error inferior 
% en comparación con el método del trapecio. 
% La jerarquía de precisión entre los métodos evaluados 
% para esta actividad se puede establecer de la siguiente manera: 
% Simpson 1/3 > Simpson 3/8 > Punto Medio > Trapecio.

% PART 3: VELOCITY ESTIMATION GIVEN ALL X VALUES
% Evaluation Integration according method's requirements individually

In = zeros([length(x) 1]);
for i=2:length(x)
    y_temp = F(1:i);
    b = x(i);
    if mod(i, 2) ~= 0 
        In(i) = simpson13(0, y_temp, a, b, i-1);
    else
        In(i) = trapezium(0, y_temp, a, b, i);
    end

end

% Velocoity at x = 0.5 m with m = 0.075 kg. 

% Ec = W
% Ec = 1/2 * (m * V^2)
% W = F * x  ==> I
% I = 1/2 (m * V^2)
% 2 * I / m = V^2
% V = sqrt((2*I) / m)

V_numeric = sqrt((2 .* In) / m);

part_2_table = table(x, F', In, V_numeric, 'VariableNames', ...
    {'x (m)', 'F (N)' ,'I', 'V calc. m/s'});
disp(part_2_table)

% V @ 0.5 m
V_05 = interp1(part_2_table.("x (m)"), part_2_table.("V calc. m/s"), 0.5, 'linear');

disp(['When x = ', num2str(0.5),' kg', ', V = ', num2str(V_05), ' m/s']);

figure;
plot(x, F, 'k-', 'LineWidth', 2);
hold on;
x_fill = [a, x(2:end)', b];
y_fill = [0, F(2:end), 0];
fill(x_fill, y_fill, 'b', 'FaceAlpha', 0.2);%area under the curve
xlabel('x (m)');
ylabel('F (N)');
title('Numeric Area Under the Curve (AUC)');
legend('f(x)', 'Numeric AUC', Location='best');
text(b, 25, ['Numeric AUC: ', num2str(In(end))], 'HorizontalAlignment', 'right');

%% BONUS PART %%

t = linspace(a, b, 100);
[p, err_lagrange, coeffs] = lagrange_interpolation(x, F, 2, t);

syms z; 
p_sym = poly2sym(coeffs, z);% symbolic polynomial
integral_sym = int(p_sym, z); % Symbolic polynomial integration
symb_auc = double(subs(integral_sym, z, b) - subs(integral_sym, z, a)); % Barrow's Rule

fprintf('Analytical AUC between %.2f y %.2f: %.4f\n', a, b, symb_auc);


figure;
plot(x, F, 'o', t, p, '-', 'LineWidth', 2);
legend('Actual', 'Fitted', Location='best');
xlabel('x (m)');
ylabel('F (N)');
title('Polynomial Data Fitting');
text(b, 25, ['Analytical AUC: ', num2str(symb_auc)], ...
    'HorizontalAlignment', 'right');

err_rel = (abs(symb_auc - In(end)) / symb_auc) * 100;

fprintf('Relative Error AUC: %.4f\n', err_rel);

figure;
stacked_fig_table = array2table([F', V_numeric], 'VariableNames', ...
    {'Force (N)', 'Velocity (m/s)'});
stackedplot(stacked_fig_table, 'cyan-', LineWidth=1, ...
    Marker='o', MarkerSize=4, MarkerFaceColor=[1.00,0.41,0.16], ...
    MarkerEdgeColor =[1.00,0.41,0.16], FontSize=11, GridVisible='on', ...
    Title='Velocity and Force Behaviour');

