
close all
clear all
clc

% run test (it needs about 7 hours)
% [Rec, NLP_reformulation_name] = run_test_gap_reformulation();

save('Data_performance_test.mat', 'Rec', 'NLP_reformulation_name')

%%
Data_performance_test = load('Data_performance_test.mat');

%
performance_matrix = Data_performance_test.Rec.time;
solver_name = Data_performance_test.NLP_reformulation_name;

%
[n_p, n_s] = size(performance_matrix);
% problem_ploted_index = reshape([1 : 4: 44; 2 : 4: 44], [], 1);
solver_ploted_index = 3:7;

plot_performance_profile(performance_matrix(:, solver_ploted_index), solver_name(solver_ploted_index))
