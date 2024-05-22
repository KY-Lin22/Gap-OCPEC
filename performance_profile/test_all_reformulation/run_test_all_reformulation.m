function [Rec, NLP_reformulation_name] = run_test_all_reformulation()
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%% OCPEC example to be tested 
OCPEC_func_handle = {...
    @Vieira_LCS_Analytic_1;...
    @Vieira_LCS_Analytic_2;...
    @Vieira_LCS_Rel_Deg_One;...
    @Vieira_LCS_High_Dim;...
    @Vieira_LCS_Without_Penalty;...
    @Vieira_LCS_With_Penalty_1;...
    @Vieira_LCS_With_Penalty_2;...
    @Vieira_LCS_With_Penalty_3;...
    @Vieira_LCS_Control_Jump};

nStages_sequ = {40, 50, 80, 100, 125, 200, 250, 400};

%% relaxed NLP reformulation to be tested
% KKT
NLP_option_KKT_Scholtes = struct('relaxation_problem', 'KKT_based',...
    'KKT_complementarity_relaxation_strategy', 'Scholtes');
NLP_option_KKT_Lin_Fukushima = struct('relaxation_problem', 'KKT_based',...
    'KKT_complementarity_relaxation_strategy', 'Lin_Fukushima');
NLP_option_KKT_Kadrani = struct('relaxation_problem', 'KKT_based',...
    'KKT_complementarity_relaxation_strategy', 'Kadrani');
NLP_option_KKT_Steffensen_Ulbrich = struct('relaxation_problem', 'KKT_based',...
    'KKT_complementarity_relaxation_strategy', 'Steffensen_Ulbrich');
NLP_option_KKT_Kanzow_Schwartz = struct('relaxation_problem', 'KKT_based',...
    'KKT_complementarity_relaxation_strategy', 'Kanzow_Schwartz');
% gap (primal)
NLP_option_primal_gap_symbolic = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
    'gap_func_implementation', 'symbolic',...
    'primal_gap_param_c', 1);
% NLP_option_primal_gap_codegen_jac = struct('relaxation_problem', 'gap_constraint_based',...
%     'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
%     'gap_func_implementation', 'codegen_jac',...
%     'primal_gap_param_c', 1);
% gap (D)
param_a = 0.1;
param_b = 10;
NLP_option_D_gap_symbolic = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_D_gap',...
    'gap_func_implementation', 'symbolic',...
    'D_gap_param_a', param_a, 'D_gap_param_b', param_b);
% NLP_option_D_gap_codegen_jac = struct('relaxation_problem', 'gap_constraint_based',...
%     'gap_constraint_relaxation_strategy', 'generalized_D_gap',...
%     'gap_func_implementation', 'codegen_jac',...
%     'D_gap_param_a', param_a, 'D_gap_param_b', param_b);
NLP_option_set = {...
    NLP_option_KKT_Scholtes,...
    NLP_option_KKT_Lin_Fukushima,...
    NLP_option_KKT_Kadrani,...
    NLP_option_KKT_Steffensen_Ulbrich,...
    NLP_option_KKT_Kanzow_Schwartz,...
    NLP_option_primal_gap_symbolic,...
    NLP_option_D_gap_symbolic};
NLP_reformulation_name = {'KKT (Scholtes)', 'KKT (Lin-Fuku)', 'KKT (Kadrani)', ...
    'KKT (Steffensen-Ulbrich)', 'KKT (Kanzow-Schwartz)', ...
    'Gap (primal, symbolic, c = 1)', ...
    ['Gap (D, symbolic, a = ' num2str(param_a) ', b = ' num2str(param_b) ')']};

%% solver option set
ipopt_tol = 1e-6; % default 1e-8
max_iter = 2000; % default 3000
base_point = 0.8;
VI_nat_res_tol = 1e-2;

% KKT (Scholtes)
solver_option_KKT_Scholtes = IPOPT_Based_Solver.create_Option();
solver_option_KKT_Scholtes.NLP_Solver.ipopt.tol = ipopt_tol;
solver_option_KKT_Scholtes.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_KKT_Scholtes.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
solver_option_KKT_Scholtes.Homotopy.kappa_s_times = base_point^2;
solver_option_KKT_Scholtes.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

% KKT (Lin-Fukushima)
solver_option_KKT_Lin_Fukushima = IPOPT_Based_Solver.create_Option();
solver_option_KKT_Lin_Fukushima.NLP_Solver.ipopt.tol = ipopt_tol;
solver_option_KKT_Lin_Fukushima.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_KKT_Lin_Fukushima.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
solver_option_KKT_Lin_Fukushima.Homotopy.kappa_s_times = base_point^1;
solver_option_KKT_Lin_Fukushima.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

% KKT (Kadrani)
solver_option_KKT_Kadrani = IPOPT_Based_Solver.create_Option();
solver_option_KKT_Kadrani.NLP_Solver.ipopt.tol = ipopt_tol; 
solver_option_KKT_Kadrani.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_KKT_Kadrani.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
solver_option_KKT_Kadrani.Homotopy.kappa_s_times = base_point^1;
solver_option_KKT_Kadrani.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

% KKT (Steffensen-Ulbrich)
solver_option_KKT_Steffensen_Ulbrich = IPOPT_Based_Solver.create_Option();
solver_option_KKT_Steffensen_Ulbrich.NLP_Solver.ipopt.tol = ipopt_tol; 
solver_option_KKT_Steffensen_Ulbrich.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_KKT_Steffensen_Ulbrich.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
solver_option_KKT_Steffensen_Ulbrich.Homotopy.kappa_s_times = base_point^1;
solver_option_KKT_Steffensen_Ulbrich.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

% KKT (Kanzow-Schwartz)
solver_option_KKT_Kanzow_Schwartz = IPOPT_Based_Solver.create_Option();
solver_option_KKT_Kanzow_Schwartz.NLP_Solver.ipopt.tol = ipopt_tol; 
solver_option_KKT_Kanzow_Schwartz.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_KKT_Kanzow_Schwartz.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
solver_option_KKT_Kanzow_Schwartz.Homotopy.kappa_s_times = base_point^1;
solver_option_KKT_Kanzow_Schwartz.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

% gap (primal, symbolic)
solver_option_primal_gap_symbolic = IPOPT_Based_Solver.create_Option();
solver_option_primal_gap_symbolic.NLP_Solver.ipopt.tol = ipopt_tol; 
solver_option_primal_gap_symbolic.NLP_Solver.ipopt.max_iter = max_iter;
solver_option_primal_gap_symbolic.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
solver_option_primal_gap_symbolic.Homotopy.kappa_s_times = base_point^2;
solver_option_primal_gap_symbolic.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

% gap (primal, codegen_jac)
% solver_option_primal_gap_codegen_jac = IPOPT_Based_Solver.create_Option();
% solver_option_primal_gap_codegen_jac.NLP_Solver.ipopt.tol = ipopt_tol; 
% solver_option_primal_gap_codegen_jac.NLP_Solver.ipopt.max_iter = max_iter; 
% solver_option_primal_gap_codegen_jac.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
% solver_option_primal_gap_codegen_jac.Homotopy.kappa_s_times = base_point^2;
% solver_option_primal_gap_codegen_jac.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

% gap (D, symbolic)
solver_option_D_gap_symbolic = IPOPT_Based_Solver.create_Option();
solver_option_D_gap_symbolic.NLP_Solver.ipopt.tol = ipopt_tol;
solver_option_D_gap_symbolic.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_D_gap_symbolic.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
solver_option_D_gap_symbolic.Homotopy.kappa_s_times = base_point^2;
solver_option_D_gap_symbolic.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

% gap (D, codegen_jac)
% solver_option_D_gap_codegen_jac = IPOPT_Based_Solver.create_Option();
% solver_option_D_gap_codegen_jac.NLP_Solver.ipopt.tol = ipopt_tol;
% solver_option_D_gap_codegen_jac.NLP_Solver.ipopt.max_iter = max_iter;
% solver_option_D_gap_codegen_jac.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
% solver_option_D_gap_codegen_jac.Homotopy.kappa_s_times = base_point^2;
% solver_option_D_gap_codegen_jac.Homotopy.VI_nat_res_tol = VI_nat_res_tol;

solver_option_set = {...
    solver_option_KKT_Scholtes,...
    solver_option_KKT_Lin_Fukushima,...
    solver_option_KKT_Kadrani,...
    solver_option_KKT_Steffensen_Ulbrich,...
    solver_option_KKT_Kanzow_Schwartz,...
    solver_option_primal_gap_symbolic,...
    solver_option_D_gap_symbolic};

%% generate solver set
solver_set = create_solver_set(OCPEC_func_handle, nStages_sequ, NLP_option_set, solver_option_set);

%% run test
% parameter
s_End = 1e-10;
param_set = {...
    struct('p_Init', 1, 'p_End', s_End),...% KKT (Scholtes)
    struct('p_Init', 1, 'p_End', s_End),...% KKT (Lin-Fuku)
    struct('p_Init', 1, 'p_End', s_End),...% KKT (Kadrani)
    struct('p_Init', 2*pi/(pi-2), 'p_End', s_End),...% KKT (Steffensen-Ulbrich)
    struct('p_Init', 1, 'p_End', s_End),...% KKT (Kanzow-Schwartz)
    struct('p_Init', 0.5, 'p_End', s_End),...% Gap (primal, symbolic)
    struct('p_Init', 1-(param_a)/2-1/(2*param_b), 'p_End', s_End)...% Gap (D, symbolic), s_Init = 1 - a/2 - 1/(2b)
    };
% solve
Rec = run_solver_test(solver_set, param_set);

end