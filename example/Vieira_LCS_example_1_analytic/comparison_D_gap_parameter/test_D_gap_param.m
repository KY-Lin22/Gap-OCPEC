clear all
clc

%% create solver set based on D gap parameter sequence
% OCPEC time parameter
timeHorizon = 1;
nStages = 100;
% D gap parameter
param_a = {0.1, 0.3, 0.5, 0.7, 0.9};
param_b = {10,  3.3, 2,   1.4, 1.1};
% relaxation and smoothing parameter
s_Init = 1e0;
s_End = 1e-16; 
sigma_Init = 1e-1;
sigma_End = 1e-6;
kappa_s_times = 0.95;
kappa_s_exp = 1.0;
kappa_sigma_times = 0.8;
kappa_sigma_exp = 1.2;
% tolerance
KKT_error_tol = 1e-4;
VI_nat_res_tol = 1e-4;
% dynamical system parameter
dtau = 0.001;
epsilon = 1000;
% init
D_gap_reformulation_name = cell(1, numel(param_a));
solver_set = cell(1, numel(param_a));
for i = 1 : numel(param_a)
    % construct OCPEC problem
    OCPEC_i = OCPEC_Vieira_LCS_analytic();
    OCPEC_i.timeHorizon = timeHorizon;
    OCPEC_i.nStages = nStages;
    OCPEC_i.timeStep = OCPEC_i.timeHorizon ./ OCPEC_i.nStages;
    % discretize OCPEC into a NLP problem
    NLP_option_i.relaxation_problem = 'gap_constraint_based'; 
    NLP_option_i.gap_constraint_relaxation_strategy = 'generalized_D_gap'; 
    NLP_option_i.gap_func_implementation = 'symbolic';
    NLP_option_i.strongly_convex_func = 'quadratic';
    NLP_option_i.D_gap_param_a = param_a{i};
    NLP_option_i.D_gap_param_b = param_b{i};
    NLP_i = NLP_Formulation(OCPEC_i, NLP_option_i);
    D_gap_reformulation_name{i} = ['Gap (D, symbolic, a = ' num2str(param_a{i}) ', b = ' num2str(param_b{i}) ')'];
    % create DynSys_Based_Solver
    Solver_option_i = DynSys_Based_Solver.create_Option();
    Solver_option_i.Continuation.dtau = dtau;
    Solver_option_i.Continuation.epsilon = epsilon;
    Solver_option_i.Continuation.integration_method = 'RK4'; % 'explitic_Euler', 'RK4'
    Solver_option_i.Continuation.tol.KKT_error = KKT_error_tol;
    Solver_option_i.Continuation.tol.VI_nat_res = VI_nat_res_tol;
    Solver_option_i.Continuation.kappa_s_times = kappa_s_times;
    Solver_option_i.Continuation.kappa_s_exp = kappa_s_exp;
    Solver_option_i.Continuation.sigma_Init = sigma_Init;
    Solver_option_i.Continuation.sigma_End = sigma_End;
    Solver_option_i.Continuation.kappa_sigma_times = kappa_sigma_times;
    Solver_option_i.Continuation.kappa_sigma_exp = kappa_sigma_exp;
    solver_i = DynSys_Based_Solver(OCPEC_i, NLP_i, Solver_option_i);
    % save
    solver_set{i} = solver_i;
end

%% run test
% init record
rec.name = D_gap_reformulation_name;
rec.z_Init = cell(1, numel(param_a));
rec.z_Opt = cell(1, numel(param_a));
rec.Info = cell(1, numel(param_a));
% run
for i = 1 : numel(param_a)
    solver_i = solver_set{i};
    z_Init_i = ones(solver_i.NLP.Dim.z, 1);
    % solve
    [z_Opt_i, Info_i] = solver_i.solve_NLP(z_Init_i, s_Init, s_End);
    rec.z_Init{i} = z_Init_i;
    rec.z_Opt{i} = z_Opt_i;
    rec.Info{i} = Info_i;
end
save('Data_test_D_gap_param.mat', 'rec')

