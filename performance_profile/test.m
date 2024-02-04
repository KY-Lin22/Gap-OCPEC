
close all
clear all
clc

% create OCPEC problem
OCPEC_problem_set = create_OCPEC_problem_set();
% create relaxed NLP problem
[NLP_problem_set, NLP_name] = create_NLP_problem_set(OCPEC_problem_set);
% create solver
solver_set = create_solver_set(OCPEC_problem_set, NLP_problem_set);
% run test
Rec = run_test(solver_set);

save('Data_performance_test.mat', 'Rec', 'NLP_name')

%%
% Data_performance_test = load('Data_performance_test.mat');

%
performance_matrix = Rec.time;
solver_name = NLP_name;

%
n_p = size(performance_matrix, 1);
% problem_ploted_index = reshape([1 : 4: 44; 2 : 4: 44], [], 1);
solver_ploted_index = [1, 4, 5, 6, 7];

plot_performance_profile(performance_matrix(:, solver_ploted_index), solver_name(solver_ploted_index))


