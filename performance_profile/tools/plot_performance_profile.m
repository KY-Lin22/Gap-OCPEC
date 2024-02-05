function plot_performance_profile(performance_matrix, solver_name)
%UNTITLED33 Summary of this function goes here
%   Detailed explanation goes here
% performance_matrix: rows are the benchmark problems, n_p, 
%                     columns are the solvers, n_s

% number of problem and solver
n_p = size(performance_matrix, 1);
n_s = size(performance_matrix, 2);
% find the best performance on each problem
best_performance_vector = min(performance_matrix, [], 2);
best_performance_vector(best_performance_vector == inf) = 1; % specify the problem that all solver fails as one
% compute performance ratio matrix
performance_ratio_matrix = performance_matrix ./ best_performance_vector;
% replace inf in p_r_mat_sort by max time factor tau for the comparison with p_r_mat
performance_ratio_matrix_sorted = sort(performance_ratio_matrix);
tauMax = max(performance_ratio_matrix(performance_ratio_matrix < inf),[],'all');
performance_ratio_matrix_sorted(performance_ratio_matrix_sorted == inf) = tauMax;
% plot
lineStyles = {'-','--',':','-.'};
X_limit = 10*tauMax;
figure(1)
for i = 1 : n_s
    p_r_mat_i = performance_ratio_matrix(:, i);
    p_r_mat_sorted_i = performance_ratio_matrix_sorted(:, i);
    distribu_func_i = zeros(n_p, 1);
    for j = 1 : n_p
        distribu_func_i(j) = length(p_r_mat_i(p_r_mat_i<=p_r_mat_sorted_i(j)))/n_p;
    end
    if p_r_mat_sorted_i(1) == 1
        Y_Init = distribu_func_i(1); % means solver i dominates some problems
    else
        Y_Init = 0;
    end
    stairs([1; p_r_mat_sorted_i; X_limit], [Y_Init; distribu_func_i; distribu_func_i(end)],...
        'LineWidth', 2, 'LineStyle',lineStyles{mod(i,4) + 1})
    hold on
end
legend(solver_name,'Location','southeast');
set(gca, 'XScale', 'log');
xlim([1 X_limit]);
ylim([0 1]);
xlabel('\tau')
ylabel('fraction of problem solved')
end