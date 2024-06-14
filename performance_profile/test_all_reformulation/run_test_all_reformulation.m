function [Rec, NLP_reformulation_name] = run_test_all_reformulation()
%UNTITLED8 Summary of this function goes here
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

% KKT
KKT_relaxation_strategy = {'Scholtes', 'Lin_Fukushima', 'Kadrani', ...
    'Steffensen_Ulbrich', 'Kanzow_Schwartz'};
KKT_relaxation_name = {'KKT (Scholtes)', 'KKT (Lin-Fukushima)', 'KKT (Kadrani)', ...
    'KKT (Steffensen-Ulbrich)', 'KKT (Kanzow-Schwartz)'};
for i = 1 : numel(KKT_relaxation_strategy)
    NLP_option_set{end + 1} = struct('relaxation_problem', 'KKT_based',...
    'KKT_complementarity_relaxation_strategy', KKT_relaxation_strategy{i});
    NLP_reformulation_name{end + 1} = KKT_relaxation_name{i};
end

% gap (primal, symbolic)
param_c = {1, 0.1};
for i = 1 : numel(param_c)
    NLP_option_set{end + 1} = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
    'gap_func_implementation', 'symbolic',...
    'primal_gap_param_c', param_c{i});    
    NLP_reformulation_name{end + 1} =  ['Gap (primal, symbolic, c = ' num2str(param_c{i}) ')'];
end

% gap (D, symbolic)
param_a = {0.1};
param_b = {10};
for i = 1 : numel(param_a)
    NLP_option_set{end + 1} = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_D_gap',...
    'gap_func_implementation', 'symbolic',...
    'D_gap_param_a', param_a{i}, 'D_gap_param_b', param_b{i});
    NLP_reformulation_name{end + 1} = ['Gap (D, symbolic, a = ' num2str(param_a{i}) ', b = ' num2str(param_b{i}) ')'];
end

%% solver option set
solver_option_set = {};
% option
ipopt_tol = 1e-6; % default 1e-8
max_iter = 2000; % default 3000
hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory';
base_point = 0.5;
VI_nat_res_tol = 1e-2;

% KKT (Scholtes): base_point^2
solver_option_KKT_Scholtes = IPOPT_Based_Solver.create_Option();
solver_option_KKT_Scholtes.NLP_Solver.ipopt.tol = ipopt_tol;
solver_option_KKT_Scholtes.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_KKT_Scholtes.NLP_Solver.ipopt.hessian_approximation = hessian_approximation;
solver_option_KKT_Scholtes.Homotopy.kappa_s_times = base_point^2;
solver_option_KKT_Scholtes.Homotopy.VI_nat_res_tol = VI_nat_res_tol;
solver_option_set{end + 1} = solver_option_KKT_Scholtes;

% KKT (Lin-Fukushima, Kadrani, Steffensen-Ulbrich, Kanzow-Schwartz): base_point^1
solver_option_KKT_LF_Kad_SU_KS = IPOPT_Based_Solver.create_Option();
solver_option_KKT_LF_Kad_SU_KS.NLP_Solver.ipopt.tol = ipopt_tol;
solver_option_KKT_LF_Kad_SU_KS.NLP_Solver.ipopt.max_iter = max_iter; 
solver_option_KKT_LF_Kad_SU_KS.NLP_Solver.ipopt.hessian_approximation = hessian_approximation;
solver_option_KKT_LF_Kad_SU_KS.Homotopy.kappa_s_times = base_point^1;
solver_option_KKT_LF_Kad_SU_KS.Homotopy.VI_nat_res_tol = VI_nat_res_tol;
for i = 1 : (numel(KKT_relaxation_strategy) - 1)
    solver_option_set{end + 1} = solver_option_KKT_LF_Kad_SU_KS;
end

% gap (primal and D): base_point^2
solver_option_gap_symbolic = IPOPT_Based_Solver.create_Option();
solver_option_gap_symbolic.NLP_Solver.ipopt.tol = ipopt_tol; 
solver_option_gap_symbolic.NLP_Solver.ipopt.max_iter = max_iter;
solver_option_gap_symbolic.NLP_Solver.ipopt.hessian_approximation = hessian_approximation;
solver_option_gap_symbolic.Homotopy.kappa_s_times = base_point^2;
solver_option_gap_symbolic.Homotopy.VI_nat_res_tol = VI_nat_res_tol;
for i = 1 : (numel(param_c) + numel(param_a))
    solver_option_set{end + 1} = solver_option_gap_symbolic;
end

%% generate solver set
solver_set = create_solver_set(OCPEC_func_handle, nStages_sequ, NLP_option_set, solver_option_set);

%% run test
param_set = {};
% init weight counter
init_weight_counter = 3;
init_weight = (base_point)^init_weight_counter;
% parameter
s_End = 1e-8;
% KKT
KKT_param_set = {...
    struct('p_Init', 1*init_weight, 'p_End', s_End),...% KKT (Scholtes)
    struct('p_Init', 1*init_weight, 'p_End', s_End),...% KKT (Lin-Fuku)
    struct('p_Init', 1*init_weight, 'p_End', s_End),...% KKT (Kadrani)
    struct('p_Init', (2*pi/(pi-2))*init_weight, 'p_End', s_End),...% KKT (Steffensen-Ulbrich)
    struct('p_Init', 1*init_weight, 'p_End', s_End)...% KKT (Kanzow-Schwartz)
    };
for i = 1 : numel(KKT_relaxation_strategy)
    param_set{end + 1} = KKT_param_set{i};
end
% gap (primal): s_Init = 1 - c/2
for i = 1 : numel(param_c)
    param_set{end + 1} = struct('p_Init', (1-(param_c{i})/2)*init_weight, 'p_End', s_End);
end
% gap (D): s_Init = 1 - a/2 - 1/(2b)
for i = 1 : numel(param_a)
    param_set{end + 1} = struct('p_Init', (1-(param_a{i})/2-1/(2*param_b{i}))*init_weight, 'p_End', s_End); 
end
% solve
Rec = run_solver_test(solver_set, param_set);

end