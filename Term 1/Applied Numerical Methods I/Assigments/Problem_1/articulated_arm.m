% This script solves articulate arm problem by finite differences,
% This script corresponds to Activity 1 of Numerical Method 1
% it also uses polynomial interpolation methods to estimate intermediate
% values.
% UNIR
% Author: Javier Blanco
% Course: Numerical Method 1, UNIR
% Cohort: 2023 - 2024

% Date: November 19, 2023
% Version: 1.0
% MATLAB Version: R2021b or later

clear;clc
format long;
%Upload data

data = load('input_data.mat', 'data');
array_size = size(data.data);
n = array_size(2);

low_boundary = data.data(1,1); %alpha = 0
high_boundary = data.data(1,n);
h = (high_boundary - low_boundary) / (n - 1);
alphas = data.data(1,:);
bethas = data.data(2,:);
dadt = 25.0;


%transform bethas radians into degrees
bethas_deg = rad2deg(bethas);
%disp(bethas);
%disp(h)

% PART 1

y_1 = finite_diff(bethas_deg,n,h,'forward');
y_2 =  finite_diff(bethas_deg, n, h, 'backward');
y_n =  finite_diff(bethas_deg, n, h, 'central');

dB = y_n;
dB(1,1) = y_1;
dB(1,n) = y_2;

% According activity, we have this relation:
% dBdt = dBda*dadt

dBdt = dB*dadt;

T = array2table([alphas' bethas', bethas_deg' dBdt(1,:)'], ...
    "VariableNames", ...
    {'Alpha (degree)','Betha (radians)', 'Betha (degree)', 'dB/dt (rad/seg)'});

disp(T);

fig_table = array2table([alphas', bethas_deg' dBdt(1,:)'], ...
"VariableNames", ...
    {'Alpha (degree)','Betha (degree)', 'dB/dt (rad/seg)'});

% PART 2

betha_interpolation = lagrange_interpolation(alphas, bethas, 6, 12);
disp(['Betha value at Alpha = 12º: ', ...
    num2str(betha_interpolation,'%.8f'), ' radians']);


% PART 3
interpolation_point = 12;
[~, nearest_indexes] = mink(abs(alphas - interpolation_point), 3);
nearest_alphas = alphas(nearest_indexes);

[lagrange_coeffs, lagrange_err] = lagrange_interpolation(alphas, bethas, 2, nearest_alphas);
lagrange_second_derivative_point = polyder(polyder(lagrange_coeffs));
lagrange_velocity = (polyder(lagrange_coeffs))*dadt; % just because of the relation dBdt = dBda*dadt
lagrange_acceleration = polyder(lagrange_velocity);

[newton_coeffs, newton_err] = newton_interpolation(alphas, bethas, 2, nearest_alphas);
newton_second_derivative_point = polyder(polyder(newton_coeffs));
newton_velocity = (polyder(newton_coeffs))*dadt; % just because of the relation dBdt = dBda*dadt
newton_acceleration = polyder(newton_velocity);

% disp(["second derivative of Betha position by Lagrange's in Alpha = 12º: ", ...
%     num2str(lagrange_second_derivative_point,'%.8f'),' radians']);
% 
% disp(["second derivative of Betha position by Newton's in Alpha = 12º: ", ...
%     num2str(lagrange_second_derivative_point,'%.8f'),' radians']);
% 
% disp(["Angular acceleration by Lagrange's in Alpha = 12º: ", ...
%     num2str(lagrange_acceleration,'%.8f'),' rad/seg^2']);
% 
% disp(["Angular acceleration by Newton's in Alpha = 12º: ", ...
%     num2str(newton_acceleration, '%.8f'), ' rad/seg^2']);

% Crear una tabla con los resultados
results_table = array2table( ...
    [lagrange_second_derivative_point, newton_second_derivative_point, ...
    lagrange_acceleration, newton_acceleration], ...
    "VariableNames", {'f 2nd Lagrange', 'f 2nd Newton', ...
                      'Lagranges (rad/seg^2)', ...
                      'Newtons (rad/seg^2)'});

% Mostrar la tabla
disp('Results for Alpha = 12º:');
disp(results_table);


% Saving results
save('output_vars.mat')
writetable(T, 'output_results.dat', 'Delimiter',' ');

% stackedplot(fig_table);
