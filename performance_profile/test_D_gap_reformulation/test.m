close all
clear all
clc

% run test (need about 3 hours)
% [Rec, NLP_reformulation_name] = run_test_D_gap_reformulation();

save('Data_performance_test.mat', 'Rec', 'NLP_reformulation_name')

%%
% Data_performance_test = load('Data_performance_test.mat');

%
performance_matrix = Rec.time;
solver_name = NLP_reformulation_name;

%
n_p = size(performance_matrix, 1);
% problem_ploted_index = reshape([1 : 4: 44; 2 : 4: 44], [], 1);
% solver_ploted_index = [1, 2, 3, 4, 5, 6, 7, 8];
solver_ploted_index = 1 : 16;
plot_performance_profile(performance_matrix(:, solver_ploted_index), solver_name(solver_ploted_index))