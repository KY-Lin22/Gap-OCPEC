function [Rec, NLP_reformulation_name] = run_test_all_reformulation_fast_update()
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
    @Vieira_LCS_Control_Jump;...
    @Vieira_LCS_State_Jump_1;...
    @Vieira_LCS_State_Jump_2;...
    @Stewart_FDI_Init_Val;...
    @Affine_DVI_Control_Jump};

nStages_sequ = {50, 80, 100, 200, 250, 400};

%% relaxed NLP reformulation to be tested
NLP_option_set = {...
    struct('relaxation_problem', 'KKT_based',...
    'KKT_complementarity_relaxation_strategy', 'Scholtes'),...
    struct('relaxation_problem', 'KKT_based',...
    'KKT_complementarity_relaxation_strategy', 'Lin_Fukushima'),...
    struct('relaxation_problem', 'KKT_based',...
    'KKT_complementarity_relaxation_strategy', 'Kadrani'),...
    struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
    'penalty_gap_func_auxiliary_variable', 'none',...
    'gap_func_smooth_param', 0),...   
    struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
    'penalty_gap_func_auxiliary_variable', 'L2',...
    'gap_func_smooth_param', 0),...
    struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_D_gap',...
    'penalty_gap_func_auxiliary_variable', 'none',...
    'gap_func_smooth_param', 0.001,...
    'D_gap_param_a', 0.1, 'D_gap_param_b', 10),...
    struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_D_gap',...
    'penalty_gap_func_auxiliary_variable', 'L2',...
    'gap_func_smooth_param', 0.001,...
    'D_gap_param_a', 0.1, 'D_gap_param_b', 10),...
    };
NLP_reformulation_name = {'KKT (Scholtes)', 'KKT (Lin-Fuku)', 'KKT (Kadrani)', ...
    'Gap (primal)', 'Gap (primal, pen)', 'Gap (D)', 'Gap (D, pen)'};

%% solver option
solver_Option.kappa_s_times = 0.2;
solver_Option.kappa_s_exp = 1.5;
solver_Option.VI_nat_res_tol = 1e-2;

%% generate solver set
solver_set = create_solver_set(OCPEC_func_handle, nStages_sequ, NLP_option_set, solver_Option);

%% run test
% parameter
s_Init = 1e0;
s_End = 1e-8;
mu_Init = 1e0;
mu_End = 1e5;
p_Init = [s_Init; mu_Init];
p_End = [s_End; mu_End];
% solve
Rec = run_solver_test(solver_set, p_Init, p_End);

end