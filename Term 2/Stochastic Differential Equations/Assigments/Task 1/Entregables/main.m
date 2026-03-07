clear, clc
format shortG

%% INITIAL PARAMETERS
% Random Initial Value Problem (RIVP)
% X'(t) = a*(b-X(t)); X(0) = X0
% t € [0, 5]

a=0;
b = 5;
N = 20;

%% MONTECARLO IMPLEMENTATION
% Simulations
sims = 600;
rv = unifrnd(0,1,[sims 1]);
X_sol_sym = zeros(N+1, sims);
h = (b - a) / N;
syms X(t)
a_p=1;
b_p=1;
ode_sym = diff(X,t) == a_p.*(b_p-X);
figure;
hold on;
for i =1:sims
    tic;
    t = a:h:b;
    cond = X(0) == rv(i);
    X_sol_sym(:,i) = double(subs(dsolve(ode_sym, cond),t));
    plot(t, X_sol_sym(:,i))
    time_elapsed_for_sym(i) = toc;
end
grid on;
xlabel('Time (t)');
ylabel('X(t)');
title('MonteCarlo Solution for RIVP with different X0');

% Expected Value and Variance
mean_X = mean(X_sol_sym, 2);
var_X = var(X_sol_sym, 0, 2);



% Visual Results
figure;
t1 = a:h:b;
plot(t1, mean_X, 'b-', 'LineWidth', 2);
hold on;
plot(t1, mean_X + sqrt(var_X), 'g--', 'LineWidth', 1);
plot(t1, mean_X - sqrt(var_X), 'r--', 'LineWidth', 1);
grid on;
xlabel('Time (t)');
ylabel('X(t)');
title('Mean and Variance of RIVP Solution');
legend('Mean', 'Mean + Std', 'Mean - Std');


%% NUMERICAL SOLUTION (JUST FOR FUN)
X_sol = zeros(N+1,sims);
t_sol = zeros(N+1,sims);

figure;
hold on;
for i=1:sims
    tic;
    [t_sol, X_sol(:,i)] = AdamsBashforth4('logitEqn', a, b, N, rv(i));
    plot(t_sol, X_sol(:,i));
    time_elapsed_for_num(i) = toc;
end
grid on;
xlabel('Time (t)');
ylabel('X(t)');
title('MonteCarlo Solution for Numerical RIVP with different X0');

%% STATISTIC ANALYSIS
% Histograms for each time step
num_steps = 6; % Number of time steps
figure;
for step = 1:num_steps
    subplot(2,3,step); % 2x3 subplots
    histogram(X_sol_sym(step,:), 10, 'Normalization', 'probability', 'FaceColor', 'b');
    xlabel(sprintf('X(%d)', step-1));
    ylabel('Probability');
    title(sprintf('Frecuency for X(%d) with different X0', step-1));
end


% Cumulative Histogram of Solutions
figure;
[N, edges] = histcounts(X_sol_sym, 'Normalization', 'probability');
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
cumulative_probs = cumsum(N, 1);
plot(bin_centers, cumulative_probs, 'b-', 'LineWidth', 2);
xlabel('X(t)');
ylabel('Cumulative Probability');
title('Cumulative Histogram of Solutions for RIVP with different X0 (Symbolic)');
grid on;

%% COMPUTATIONAL EFICIENCY
% Time elapsed execution plots
figure;
bar([mean(time_elapsed_for_sym); mean(time_elapsed_for_num)]);
xlabel('Simulation');
ylabel('Time (s)');
title('Mean Time Elapsed for Each Simulation');
xticklabels({'Symbolic','Numerical'});
