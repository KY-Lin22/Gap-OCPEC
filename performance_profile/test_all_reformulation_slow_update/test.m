
close all
clear all
clc

% run test (it needs about 30 - 40 min)
[Rec, NLP_reformulation_name] = run_test_all_reformulation_slow_update();

save('Data_performance_test.mat', 'Rec', 'NLP_reformulation_name')

%%
Data_performance_test = load('Data_performance_test.mat');

%
performance_matrix = Data_performance_test.Rec.time;
solver_name = Data_performance_test.NLP_reformulation_name;

%
n_p = size(performance_matrix, 1);
% problem_ploted_index = reshape([1 : 4: 44; 2 : 4: 44], [], 1);
solver_ploted_index = 1:7;

plot_performance_profile(performance_matrix(:, solver_ploted_index), solver_name(solver_ploted_index))


