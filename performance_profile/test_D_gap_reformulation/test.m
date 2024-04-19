close all
clear all
clc

% run test (need about 4 hours)
% [Rec, NLP_reformulation_name] = run_test_D_gap_reformulation();

save('Data_performance_test.mat', 'Rec', 'NLP_reformulation_name')

%%
Data_performance_test = load('Data_performance_test.mat');

%
performance_matrix = Data_performance_test.Rec.time;
solver_name = Data_performance_test.NLP_reformulation_name;

%
n_p = size(performance_matrix, 1);
% problem_ploted_index = reshape([1 : 4: 44; 2 : 4: 44], [], 1);
solver_ploted_index = [7:9, 16:18, 25:27];
plot_performance_profile(performance_matrix(:, solver_ploted_index), solver_name(solver_ploted_index))
