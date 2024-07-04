function [Rec, solver_name] = run_test_Gap_DynSys()
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

nStages_sequ = {100, 200, 300, 400, 500};

%% relaxed NLP reformulation to be tested
NLP_option_set = {};

% primal gap 
param_c = {1};
primal_gap_name = {};
for i = 1 : numel(param_c)
    NLP_option_set{end + 1} = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_primal_gap',...
    'gap_func_implementation', 'symbolic',...
    'primal_gap_param_c', param_c{i});    
    primal_gap_name{end + 1} =  ['Gap (primal, c = ' num2str(param_c{i}) ')'];
end

% D gap
param_a = {0.1, 0.3, 0.5, 0.7, 0.9};
param_b = {10,  3.3, 2,   1.4, 1.1};
D_gap_name = {};
for i = 1 : numel(param_a)
    NLP_option_set{end + 1} = struct('relaxation_problem', 'gap_constraint_based',...
    'gap_constraint_relaxation_strategy', 'generalized_D_gap',...
    'gap_func_implementation', 'symbolic',...
    'D_gap_param_a', param_a{i}, 'D_gap_param_b', param_b{i});
    D_gap_name{end + 1} = ['Gap (D, a = ' num2str(param_a{i}) ', b = ' num2str(param_b{i}) ')'];
end

%% solver option set
% DynSys-Based integration method
integration_method = 'explitic_Euler'; % 'explitic_Euler', 'RK4'
% dynamical system parameter
dtau = 0.001;
epsilon = 1000;
% tolerance
KKT_error_tol = 1e-4;
VI_nat_res_tol = 1e-2;
% relaxation parameter init and end value
s_Init = 1e0;
s_End = 1e-12;
% FB smoothing parameter init and end value
sigma_Init = 1e-2;
sigma_End = 1e-6;
% parameter update
kappa_s_times = 0.9;% slow
kappa_s_exp = 1.0;
kappa_sigma_times = 0.9;
kappa_sigma_exp = 1.0;

% init
solver_option_set = {};
solver_name = {};
solver_type = {};

% DynSys-based solver for primal gap reformulation
Solver_option_gap_DynSys = DynSys_Based_Solver.create_Option();
Solver_option_gap_DynSys.Continuation.dtau = dtau;
Solver_option_gap_DynSys.Continuation.epsilon = epsilon;
Solver_option_gap_DynSys.Continuation.integration_method = integration_method;
Solver_option_gap_DynSys.Continuation.tol.KKT_error = KKT_error_tol;
Solver_option_gap_DynSys.Continuation.tol.VI_nat_res = VI_nat_res_tol;
Solver_option_gap_DynSys.Continuation.kappa_s_times = kappa_s_times;
Solver_option_gap_DynSys.Continuation.kappa_s_exp = kappa_s_exp;
Solver_option_gap_DynSys.Continuation.sigma_Init = sigma_Init;
Solver_option_gap_DynSys.Continuation.sigma_End = sigma_End;
Solver_option_gap_DynSys.Continuation.kappa_sigma_times = kappa_sigma_times;
Solver_option_gap_DynSys.Continuation.kappa_sigma_exp = kappa_sigma_exp;

for i = 1 : numel(param_c)
    solver_option_set{end + 1} = Solver_option_gap_DynSys;
    solver_name{end + 1} = [primal_gap_name{i} ' DynSys-Based'];
    solver_type{end + 1} = 'DynSys-Based';
end

% DynSys-based solver for D gap reformulation
for i = 1 : numel(param_a)
    solver_option_set{end + 1} = Solver_option_gap_DynSys;
    solver_name{end + 1} = [D_gap_name{i} ' DynSys-Based'];
    solver_type{end + 1} = 'DynSys-Based';
end

%% generate solver set
solver_set = create_solver_set(OCPEC_func_handle, nStages_sequ, NLP_option_set, solver_option_set, solver_type);

%% run test
param_set = {};
for i = 1 : (numel(param_c) + numel(param_a))
    param_set{end + 1} = struct('s_Init', s_Init, 's_End', s_End);
end
% solve
Rec = run_solver_test(solver_set, param_set);
end