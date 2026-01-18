%%              APPLIED NUMERICAL METHODS 2 
%     Parabolic partial differential equation solution
%                       Activity 3
%                   Author: Javier Blanco

%****************************************************************%
%                           PDE                                  %
%                                                                %
%    u_t = u_{xx} + t^2 * u + x * cos(xt),  0 < x < 1,  t > 0    %
%                                                                %
%    u_x(0, t) = t,  u(1, t) = sin(t),  t > 0                    %
%                                                                %
%    u(x, 0) = 0,  0 < x < 1                                     %
%****************************************************************%

format shortG
alpha=1;
a=0;
b=1;
Tmax=0.5;
Tmax_cn = 1;
h=0.1;
k=0.0005;
nx=(b-a)/h;
nt=Tmax/k;
nt_cn=Tmax_cn/k;
ci=@(x) 0*ones(1,nx+1); %<--- Initial condition

% Methods (Explicit, Implicit, Crank-Nicholson)
[U_exp,x_exp,t_exp] = ExplicitMethodPDE(alpha,ci,a,b,nx,Tmax,nt);
[U_im,x_im,t_im] = ImplicitMethodPDE(alpha, ci, a,b, nx,Tmax,nt);
[U_cn, x_cn, t_cn]=CrankNicholsonPDE(alpha, ci,a,b,nx,Tmax_cn,nt);

%% RESULTS

% Results for slices of time t = 1/8*Tmax, 1/4*Tmax, 1/2*Tmax, 3/4*Tmax and Tmax
% ['x',      'U(x,Tmax/8)', 'U(x,Tmax/4)', 'U(x,Tmax/2)', 'U(x,3Tmax/4)', 'U(x,Tmax)']
% [x_exp',U_exp(:,(end-1)/8),U_exp(:,(end-1)/4),U_exp(:,(end-1)/2),U_exp(:,3*(end-1)/4),U_exp(:,end)]
% 
% [x_im',U_im(:,(end-1)/8),U_im(:,(end-1)/4),U_im(:,(end-1)/2),U_im(:,3*(end-1)/4),U_im(:,end)]
% 
% [x_cn',U_cn(:,(end-1)/8),U_cn(:,(end-1)/4),U_cn(:,(end-1)/2),U_cn(:,3*(end-1)/4),U_cn(:,end)]


headers = {'X', 'U_1/8', 'U_1/4', 'U_1/2', 'U_3/4', 'U_end'};

% Exp table
table_exp = array2table([x_exp', U_exp(:, (end-1)/8), U_exp(:, (end-1)/4), U_exp(:, (end-1)/2), U_exp(:, 3*(end-1)/4), U_exp(:, end)], ...
                        'VariableNames', headers);

% Im table
table_im = array2table([x_im', U_im(:, (end-1)/8), U_im(:, (end-1)/4), U_im(:, (end-1)/2), U_im(:, 3*(end-1)/4), U_im(:, end)], ...
                       'VariableNames', headers);

% CN table
table_cn = array2table([x_cn', U_cn(:, (end-1)/8), U_cn(:, (end-1)/4), U_cn(:, (end-1)/2), U_cn(:, 3*(end-1)/4), U_cn(:, end)], ...
                       'VariableNames', headers);

error_headers = {'X','Error_Tmax/8', 'Error_Tmax/4', 'Error_Tmax/2', 'Error_3Tmax/4', 'Error_Tmax'};

% Calcular el error absoluto entre las columnas de table_exp y table_im
error_abs = abs(table_exp{:, 2:end} - table_im{:, 2:end});

% Crear la tabla de errores
table_error = array2table([x_exp', error_abs], 'VariableNames', error_headers);


disp('Exp Table:');
disp(table_exp);
disp('Im Table:');
disp(table_im);

disp('Tabla de Errores Absolutos:')
disp(table_error);

disp('CN Table:');
disp(table_cn);


%% 2D Plot in Tmax
hold on
plot(x_exp, U_exp(:, end), '-r', 'DisplayName', 'Explicit', 'LineWidth', 1.5)
plot(x_im, U_im(:, end), '*b', 'DisplayName', 'Implicit', 'MarkerSize', 8)
plot(x_cn, U_cn(:, end), '--k', 'DisplayName', 'Crank-Nicholson', 'LineWidth', 1.5)
xlabel('Spatial Dimention (x)')
ylabel('U(x, t) Solution in t = Tmax')
title('Explicit, Implicit and Crank-Nicholson Methods')
legend('show', 'Location', 'best', 'FontSize', 11)
grid on
hold off


%% 3D Plot

% Meshgrid for plotting
[X_exp, T_exp] = meshgrid(x_exp, t_exp);
[X_im, T_im] = meshgrid(x_im, t_im);
[X_cn, T_cn] = meshgrid(x_cn, t_cn);

% Transposed
U_imp_T = U_im';
U_exp_T = U_exp';
U_cn_T = U_cn';

% % 3D surface plot
% figure;
% surf(X_exp, T_exp, U_exp_T);
% xlabel('Spatial Coordinate (x)');
% ylabel('Time (t)');
% zlabel('Explicit U(x,t)');
% title('Explicit Method 3D Surface Plot of U(x,t)');
% shading interp; % <--- For smoother shading
% colormap(jet)
% colorbar('eastoutside')
% 
% figure;
% surf(X_im, T_im, U_imp_T);
% xlabel('Spatial Coordinate (x)');
% ylabel('Time (t)');
% zlabel('Implicit U(x,t)');
% title('Implicit Method 3D Surface Plot of U(x,t)');
% shading interp; % <--- For smoother shading
% colormap(jet);
% colorbar('eastoutside');
% 
% 
% figure;
% surf(X_cn, T_cn, U_cn_T);
% xlabel('Spatial Coordinate (x)');
% ylabel('Time (t)');
% zlabel('Crank-Nicholson U(x,t)');
% title('Crank-Nicholson Method 3D Surface Plot of U(x,t)');
% shading interp; % <--- For smoother shading
% colormap(jet);
% colorbar('eastoutside');


% Create a figure for subplots
figure;

% Explicit Method subplot
subplot(1, 3, 1);
surf(X_exp, T_exp, U_exp_T);
xlabel('Spatial Coordinate (x)');
ylabel('Time (t)');
zlabel('Explicit U(x,t)');
title('Explicit Method 3D Surface Plot of U(x,t)');
shading interp; % For smoother shading
colormap(jet);

% Implicit Method subplot
subplot(1, 3, 2);
surf(X_im, T_im, U_imp_T);
xlabel('Spatial Coordinate (x)');
ylabel('Time (t)');
zlabel('Implicit U(x,t)');
title('Implicit Method 3D Surface Plot of U(x,t)');
shading interp; % For smoother shading
colormap(jet);

% Crank-Nicholson Method subplot
subplot(1, 3, 3);
surf(X_cn, T_cn, U_cn_T);
xlabel('Spatial Coordinate (x)');
ylabel('Time (t)');
zlabel('Crank-Nicholson U(x,t)');
title('Crank-Nicholson Method 3D Surface Plot of U(x,t)');
shading interp; % For smoother shading
colormap(jet);
