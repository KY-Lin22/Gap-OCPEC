function [Rec, NLP_reformulation_name] = run_test_gap_reformulation()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% OCPEC example to be tested 
OCPEC_func_handle = {...
    @Vieira_LCS_Analytic_1_With_Penalty_1;...
    @Vieira_LCS_Analytic_1_With_Penalty_2;...
    @Vieira_LCS_Analytic_1_Without_Penalty;...
    @Vieira_LCS_Analytic_2_With_Penalty_1;...
    @Vieira_LCS_Analytic_2_With_Penalty_2;...
    @Vieira_LCS_Analytic_2_Without_Penalty;...
    @Vieira_LCS_Rel_Deg_One_With_Penalty_1;...
    @Vieira_LCS_Rel_Deg_One_With_Penalty_2;...
    @Vieira_LCS_Rel_Deg_One_Without_Penalty;...
    @Vieira_LCS_High_Dim_With_Penalty_1;...
    @Vieira_LCS_High_Dim_With_Penalty_2;...
    @Vieira_LCS_High_Dim_Without_Penalty;...
    @Vieira_LCS_Init_Peak_With_Penalty_1;...
    @Vieira_LCS_Init_Peak_With_Penalty_2;...
    @Vieira_LCS_Init_Peak_Without_Penalty;...
    @Vieira_LCS_Control_Jump_With_Penalty_1;...
    @Vieira_LCS_Control_Jump_With_Penalty_2;...
    @Vieira_LCS_Control_Jump_Without_Penalty;...
    @Vieira_LCS_State_Jump_With_Penalty_1;...
    @Vieira_LCS_State_Jump_With_Penalty_2;...
    @Vieira_LCS_State_Jump_Without_Penalty};

% nStages_sequ = {40, 50, 80, 100, 125, 200};
nStages_sequ = {100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200};

%% relaxed NLP reformulation to be tested
NLP_option_set = {};
NLP_reformulation_name = {};
% gap (primal, symbolic)
param_c = {1, 0.7, 0.5, 0.3, 0.1};
for i = 1 : numel(param_c)
    NLP_option_set{end + 1} = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
    'gap_func_implementation', 'symbolic',...
    'primal_gap_param_c', param_c{i});
    NLP_reformulation_name{end + 1} = ['Gap (primal, symbolic, c = ' num2str(param_c{i}) ')'];
end
% gap (primal, codegen_fd)
for i = 1 : numel(param_c)
    NLP_option_set{end + 1} = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
    'gap_func_implementation', 'codegen_fd',...
    'primal_gap_param_c', param_c{i});
    NLP_reformulation_name{end + 1} = ['Gap (primal, codegen-fd, c = ' num2str(param_c{i}) ')'];
end

% gap (D, symbolic)
param_a = {0.1, 0.3, 0.5, 0.7, 0.9};
param_b = {10,  3.3, 2,   1.4, 1.1};
for i = 1 : numel(param_a)
    NLP_option_set{end + 1} = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_D_gap',...
    'gap_func_implementation', 'symbolic',...
    'D_gap_param_a', param_a{i}, 'D_gap_param_b', param_b{i});
    NLP_reformulation_name{end + 1} = ['Gap (D, symbolic, a = ' num2str(param_a{i}) ', b = ' num2str(param_b{i}) ')'];
end
% gap (D, codegen_fd)
for i = 1 : numel(param_a)
    NLP_option_set{end + 1} = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_D_gap',...
    'gap_func_implementation', 'codegen_fd',...
    'D_gap_param_a', param_a{i}, 'D_gap_param_b', param_b{i});
    NLP_reformulation_name{end + 1} = ['Gap (D, codegen-fd, a = ' num2str(param_a{i}) ', b = ' num2str(param_b{i}) ')'];
end

%% solver option set
solver_option_set = {};
% option
ipopt_tol = 1e-6; % default 1e-8
max_iter = 2000; % default 3000
hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory';
base_point = 0.8;
VI_nat_res_tol = 1e-2;

% gap (primal, for symbolic and codegen_fd)
solver_option_primal_gap_symbolic_codegen_fd = IPOPT_Based_Solver.create_Option();
solver_option_primal_gap_symbolic_codegen_fd.NLP_Solver.ipopt.tol = ipopt_tol; 
solver_option_primal_gap_symbolic_codegen_fd.NLP_Solver.ipopt.max_iter = max_iter;
solver_option_primal_gap_symbolic_codegen_fd.NLP_Solver.ipopt.hessian_approximation = hessian_approximation; 
solver_option_primal_gap_symbolic_codegen_fd.Homotopy.kappa_s_times = base_point^2;
solver_option_primal_gap_symbolic_codegen_fd.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

for i = 1 : (2*numel(param_c))
    solver_option_set{end + 1} = solver_option_primal_gap_symbolic_codegen_fd;
end

% gap (D, for all param set, for symbolic and codegen_fd)
solver_option_D_gap_symbolic_codegen_fd = IPOPT_Based_Solver.create_Option();
solver_option_D_gap_symbolic_codegen_fd.NLP_Solver.ipopt.tol = ipopt_tol;
solver_option_D_gap_symbolic_codegen_fd.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_D_gap_symbolic_codegen_fd.NLP_Solver.ipopt.hessian_approximation = hessian_approximation; 
solver_option_D_gap_symbolic_codegen_fd.Homotopy.kappa_s_times = base_point^2;
solver_option_D_gap_symbolic_codegen_fd.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

for i = 1 : (2*numel(param_a))
    solver_option_set{end + 1} = solver_option_D_gap_symbolic_codegen_fd;
end

%% generate solver set
solver_set = create_solver_set(OCPEC_func_handle, nStages_sequ, NLP_option_set, solver_option_set);

%% run test
param_set = {};
% parameter
s_End = 1e-10;
% gap (primal, symbolic): s_Init = 1 - c/2
for i = 1 : numel(param_c)
    param_set{end + 1} = struct('p_Init', 1-(param_c{i})/2, 'p_End', s_End);
end
% gap (primal, codegen_fd): s_Init = 1 - c/2
for i = 1 : numel(param_c)
    param_set{end + 1} = struct('p_Init', 1-(param_c{i})/2, 'p_End', s_End);
end
% gap (D, symbolic): s_Init = 1 - a/2 - 1/(2b)
for i = 1 : numel(param_a)
    param_set{end + 1} = struct('p_Init', 1-(param_a{i})/2-1/(2*param_b{i}), 'p_End', s_End); 
end
% gap (D, codegen_fd): s_Init = 1 - a/2 - 1/(2b)
for i = 1 : numel(param_a)
    param_set{end + 1} = struct('p_Init', 1-(param_a{i})/2-1/(2*param_b{i}), 'p_End', s_End); 
end
% solve
Rec = run_solver_test(solver_set, param_set);

end