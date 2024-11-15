
close all
clear all
clc

% run test (it needs about   hours)
% [Rec, solver_name] = run_test_KKT_IPOPT_vs_Gap_DynSys();

save('Data_performance_test.mat', 'Rec', 'solver_name')

%%
Data_performance_test = load('Data_performance_test.mat');

%
performance_matrix = Data_performance_test.Rec.time;
solver_name = Data_performance_test.solver_name;

%
[n_p, n_s] = size(performance_matrix);
% problem_ploted_index = reshape([9 : 8: n_p], [], 1);
problem_ploted_index = 1 : n_p;
solver_ploted_index = 1 : 7; % 1 : n_s

plot_performance_profile(performance_matrix(problem_ploted_index, solver_ploted_index), solver_name(solver_ploted_index))