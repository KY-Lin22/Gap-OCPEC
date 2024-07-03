clear all
clc

%% create OCPEC (all simulation use the same OCPEC)
timeHorizon = 1;
nStages = 1000;
OCPEC = OCPEC_Vieira_LCS_analytic();
OCPEC.timeHorizon = timeHorizon;
OCPEC.nStages = nStages;
OCPEC.timeStep = OCPEC.timeHorizon ./ OCPEC.nStages;

%% specify NLP and solver parameter
% primal gap parameter
param_c = {1};
% D gap parameter
param_a = {0.1, 0.3, 0.5, 0.7, 0.9};
param_b = {10,  3.3, 2,   1.4, 1.1};
% parameter init and end value
s_Init = 1e0;
s_End = 1e-16;
sigma_Init = 1e-2;
sigma_End = 1e-6;
% parameter update (slow)
kappa_s_times = 0.9;
kappa_s_exp = 1.0;
kappa_sigma_times = 0.9;
kappa_sigma_exp = 1.0;
% tolerance
KKT_error_tol = 1e-16;
VI_nat_res_tol = 1e-16;
% dynamical system parameter
dtau = 0.001;
epsilon = 1000;
% init
solver_name = {};
solver_set = {};

%% create solver set (gap, DynSys_Based)
% solver option
Solver_option_gap_DynSys = DynSys_Based_Solver.create_Option();
Solver_option_gap_DynSys.Continuation.dtau = dtau;
Solver_option_gap_DynSys.Continuation.epsilon = epsilon;
Solver_option_gap_DynSys.Continuation.integration_method = 'RK4'; % 'explitic_Euler', 'RK4'
Solver_option_gap_DynSys.Continuation.tol.KKT_error = KKT_error_tol;
Solver_option_gap_DynSys.Continuation.tol.VI_nat_res = VI_nat_res_tol;
Solver_option_gap_DynSys.Continuation.kappa_s_times = kappa_s_times;
Solver_option_gap_DynSys.Continuation.kappa_s_exp = kappa_s_exp;
Solver_option_gap_DynSys.Continuation.sigma_Init = sigma_Init;
Solver_option_gap_DynSys.Continuation.sigma_End = sigma_End;
Solver_option_gap_DynSys.Continuation.kappa_sigma_times = kappa_sigma_times;
Solver_option_gap_DynSys.Continuation.kappa_sigma_exp = kappa_sigma_exp;
% solver set (primal gap)
for i = 1 : numel(param_c)
    % discretize OCPEC into a NLP problem
    NLP_option_i.relaxation_problem = 'gap_constraint_based'; 
    NLP_option_i.gap_constraint_relaxation_strategy = 'generalized_primal_gap'; 
    NLP_option_i.gap_func_implementation = 'symbolic';
    NLP_option_i.primal_gap_param_c = param_c{i};
    NLP_i = NLP_Formulation(OCPEC, NLP_option_i);
    % create DynSys_Based_Solver
    solver_i = DynSys_Based_Solver(OCPEC, NLP_i, Solver_option_gap_DynSys);
    % save
    solver_name{end + 1} = 'Gap (primal)';
    solver_set{end + 1} = solver_i;    
end
% solver set (D gap)
for i = 1 : numel(param_a)
    % discretize OCPEC into a NLP problem
    NLP_option_i.relaxation_problem = 'gap_constraint_based'; 
    NLP_option_i.gap_constraint_relaxation_strategy = 'generalized_D_gap'; 
    NLP_option_i.gap_func_implementation = 'symbolic';
    NLP_option_i.D_gap_param_a = param_a{i};
    NLP_option_i.D_gap_param_b = param_b{i};
    NLP_i = NLP_Formulation(OCPEC, NLP_option_i);
    % create DynSys_Based_Solver
    solver_i = DynSys_Based_Solver(OCPEC, NLP_i, Solver_option_gap_DynSys);
    % save
    solver_name{end + 1} =  ['Gap (D, a = ' num2str(param_a{i}) ', b = ' num2str(param_b{i}) ')'];
    solver_set{end + 1} = solver_i;    
end

%% run test
% init record
rec.name = solver_name;
rec.z_Init = {};
rec.z_Opt = {};
rec.Info = {};
% run
repeat_num = 3;
for i = 1 : numel(solver_set)
    solver_i = solver_set{i};
    z_Init_i = ones(solver_i.NLP.Dim.z, 1);
    % solve
    for j = 1 : repeat_num
        [z_Opt_i, Info_i] = solver_i.solve_NLP(z_Init_i, s_Init, s_End);
    end    
    rec.z_Init{end + 1} = z_Init_i;
    rec.z_Opt{end + 1} = z_Opt_i;
    rec.Info{end + 1} = Info_i;
end
save('Data_test_converge_Gap_DynSys.mat', 'rec')
