function [Rec, NLP_reformulation_name] = run_test_gap_reformulation()
%UNTITLED2 Summary of this function goes here

%% OCPEC example to be tested 
OCPEC_func_handle = {...
    @Vieira_LCS_Analytic_1;...
    @Vieira_LCS_Analytic_2;...
    @Vieira_LCS_Rel_Deg_One;...
    % @Vieira_LCS_High_Dim;...
    @Vieira_LCS_Without_Penalty;...
    @Vieira_LCS_With_Penalty_1;...
    @Vieira_LCS_With_Penalty_2;...
    @Vieira_LCS_With_Penalty_3;...
    @Vieira_LCS_Control_Jump};

nStages_sequ = {40, 50, 80, 100, 125, 200, 250, 400};

%% relaxed NLP reformulation to be tested
% gap (primal, symbolic)
NLP_option_primal_gap_symbolic = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
    'gap_func_implementation', 'symbolic',...
    'primal_gap_param_c', 1);
% gap (primal, codegen_fd)
NLP_option_primal_gap_codegen_fd = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
    'gap_func_implementation', 'codegen_fd',...
    'primal_gap_param_c', 1);

NLP_option_set = {...
    NLP_option_primal_gap_symbolic,...
    NLP_option_primal_gap_codegen_fd};
NLP_reformulation_name = {...
    'Gap (primal, symbolic, c = 1)', ...
    'Gap (primal, codegen-fd, c = 1)'};

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
ipopt_tol = 1e-6; % default 1e-8
max_iter = 2000; % default 3000
base_point = 0.8;
VI_nat_res_tol = 1e-2;

% gap (primal, for symbolic and codegen_fd)
solver_option_primal_gap_symbolic_codegen_fd = IPOPT_Based_Solver.create_Option();
solver_option_primal_gap_symbolic_codegen_fd.NLP_Solver.ipopt.tol = ipopt_tol; 
solver_option_primal_gap_symbolic_codegen_fd.NLP_Solver.ipopt.max_iter = max_iter;
solver_option_primal_gap_symbolic_codegen_fd.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
solver_option_primal_gap_symbolic_codegen_fd.Homotopy.kappa_s_times = base_point^2;
solver_option_primal_gap_symbolic_codegen_fd.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

solver_option_set = {...
    solver_option_primal_gap_symbolic_codegen_fd,...
    solver_option_primal_gap_symbolic_codegen_fd};

% gap (D, for all param set, for symbolic and codegen_fd)
solver_option_D_gap_symbolic_codegen_fd = IPOPT_Based_Solver.create_Option();
solver_option_D_gap_symbolic_codegen_fd.NLP_Solver.ipopt.tol = ipopt_tol;
solver_option_D_gap_symbolic_codegen_fd.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_D_gap_symbolic_codegen_fd.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
solver_option_D_gap_symbolic_codegen_fd.Homotopy.kappa_s_times = base_point^2;
solver_option_D_gap_symbolic_codegen_fd.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

for i = 1 : numel(param_a)
    solver_option_set{end + 1} = solver_option_D_gap_symbolic_codegen_fd;
end
for i = 1 : numel(param_a)
    solver_option_set{end + 1} = solver_option_D_gap_symbolic_codegen_fd;
end

%% generate solver set
solver_set = create_solver_set(OCPEC_func_handle, nStages_sequ, NLP_option_set, solver_option_set);

%% run test
% parameter
s_End = 1e-10;
% gap (primal)
param_set = {...
    struct('p_Init', 0.5, 'p_End', s_End),...% Gap (primal, symbolic)
    struct('p_Init', 0.5, 'p_End', s_End)...% Gap (primal, codegen_fd)
    };
% gap (D, symbolic)
for i = 1 : numel(param_a)
    % Gap (D), s_Init = 1 - a/2 - 1/(2b)
    param_set{end + 1} = struct('p_Init', 1-(param_a{i})/2-1/(2*param_b{i}), 'p_End', s_End); 
end
% gap (D, codegen_fd)
for i = 1 : numel(param_a)
    % Gap (D), s_Init = 1 - a/2 - 1/(2b)
    param_set{end + 1} = struct('p_Init', 1-(param_a{i})/2-1/(2*param_b{i}), 'p_End', s_End); 
end
% solve
Rec = run_solver_test(solver_set, param_set);

end