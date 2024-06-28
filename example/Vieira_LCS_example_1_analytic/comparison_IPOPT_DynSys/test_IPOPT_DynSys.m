clear all
clc

%% create OCPEC (all simulation use the same OCPEC)
timeHorizon = 1;
nStages = 100;
OCPEC = OCPEC_Vieira_LCS_analytic();
OCPEC.timeHorizon = timeHorizon;
OCPEC.nStages = nStages;
OCPEC.timeStep = OCPEC.timeHorizon ./ OCPEC.nStages;

%% create solver set (specify parameter)
% KKT relaxation strategy for fast and slow update
KKT_relaxation_strategy = {'Scholtes', 'Lin_Fukushima', 'Kadrani', ...
    'Steffensen_Ulbrich', 'Kanzow_Schwartz'};
KKT_relaxation_name_fast = {...
    'KKT (Scholtes), fast, IPOPT-based',...
    'KKT (Lin-Fukushima), fast, IPOPT-based',...
    'KKT (Kadrani), fast, IPOPT-based', ...
    'KKT (Steffensen-Ulbrich), fast, IPOPT-based',...
    'KKT (Kanzow-Schwartz), fast, IPOPT-based'};
KKT_relaxation_name_slow = {...
    'KKT (Scholtes), slow, IPOPT-based',...
    'KKT (Lin-Fukushima), slow, IPOPT-based',...
    'KKT (Kadrani), slow, IPOPT-based',...
    'KKT (Steffensen-Ulbrich), slow, IPOPT-based',...
    'KKT (Kanzow-Schwartz), slow, IPOPT-based'};
% primal gap parameter
param_c = {1};
% D gap parameter
param_a = {0.1, 0.3, 0.5, 0.7, 0.9};
param_b = {10,  3.3, 2,   1.4, 1.1};
% relaxation and smoothing parameter
s_Init_KKT = {...
    1,...% KKT (Scholtes)
    1,...% KKT (Lin-Fuku)
    1,...% KKT (Kadrani)
    (2*pi/(pi-2)),...% KKT (Steffensen-Ulbrich)
    1 ... % KKT (Kanzow-Schwartz)
    };
s_Init_set = {};
for i = 1 : numel(s_Init_KKT)
    s_Init_set{end + 1} = s_Init_KKT{i}; % KKT fast
end
for i = 1 : numel(s_Init_KKT)
    s_Init_set{end + 1} = s_Init_KKT{i}; % KKT slow
end
for i = 1 : numel(param_c)
    s_Init_set{end + 1} = (1-(param_c{i})/2); % primal gap
end
for i = 1 : numel(param_a)
    s_Init_set{end + 1} = (1-(param_a{i})/2-1/(2*param_b{i})); % D gap
end
s_End = 1e-16;
sigma_Init = 1e-1;
sigma_End = 1e-6;
kappa_s_times_fast = 0.1; % update s fast
kappa_s_exp_fast = 1.0; % update s fast
kappa_s_times_slow = 0.9;% update s slow
kappa_s_exp_slow = 1.0; % update s slow
kappa_sigma_times = 0.8;
kappa_sigma_exp = 1.2;
% tolerance
KKT_error_tol = 1e-4;
VI_nat_res_tol = 1e-3;
% dynamical system parameter
dtau = 0.001;
epsilon = 1000;
% init
solver_name = {};
solver_set = {};

%% create solver set (IPOPT_Based, fast)
% solver option
Solver_option_IPOPT_fast = IPOPT_Based_Solver.create_Option();
Solver_option_IPOPT_fast.Continuation.kappa_s_times = kappa_s_times_fast;
Solver_option_IPOPT_fast.Continuation.kappa_s_exp = kappa_s_exp_fast;
Solver_option_IPOPT_fast.Continuation.tol.VI_nat_res = VI_nat_res_tol;
% solver set
for i = 1 : numel(s_Init_KKT)
    % discretize OCPEC into a NLP problem
    NLP_option_i.relaxation_problem = 'KKT_based'; 
    NLP_option_i.KKT_complementarity_relaxation_strategy = KKT_relaxation_strategy{i};
    NLP_i = NLP_Formulation(OCPEC, NLP_option_i);
    % create IPOPT_Based_Solver
    solver_i = IPOPT_Based_Solver(OCPEC, NLP_i, Solver_option_IPOPT_fast);
    % save
    solver_name{end + 1} = KKT_relaxation_name_fast{i};
    solver_set{end + 1} = solver_i;    
end

%% create solver set (IPOPT_Based, slow)
% solver option
Solver_option_IPOPT_slow = IPOPT_Based_Solver.create_Option();
Solver_option_IPOPT_slow.Continuation.kappa_s_times = kappa_s_times_slow;
Solver_option_IPOPT_slow.Continuation.kappa_s_exp = kappa_s_exp_slow;
Solver_option_IPOPT_slow.Continuation.tol.VI_nat_res = VI_nat_res_tol;
% solver set
for i = 1 : numel(s_Init_KKT)
    % discretize OCPEC into a NLP problem
    NLP_option_i.relaxation_problem = 'KKT_based'; 
    NLP_option_i.KKT_complementarity_relaxation_strategy = KKT_relaxation_strategy{i};
    NLP_i = NLP_Formulation(OCPEC, NLP_option_i);
    % create IPOPT_Based_Solver
    solver_i = IPOPT_Based_Solver(OCPEC, NLP_i, Solver_option_IPOPT_slow);
    % save
    solver_name{end + 1} = KKT_relaxation_name_slow{i};
    solver_set{end + 1} = solver_i;    
end

%% create solver set (DynSys_Based, slow)
% solver option
Solver_option_DynSys_slow = DynSys_Based_Solver.create_Option();
Solver_option_DynSys_slow.Continuation.dtau = dtau;
Solver_option_DynSys_slow.Continuation.epsilon = epsilon;
Solver_option_DynSys_slow.Continuation.integration_method = 'RK4'; % 'explitic_Euler', 'RK4'
Solver_option_DynSys_slow.Continuation.tol.KKT_error = KKT_error_tol;
Solver_option_DynSys_slow.Continuation.tol.VI_nat_res = VI_nat_res_tol;
Solver_option_DynSys_slow.Continuation.kappa_s_times = kappa_s_times_slow;
Solver_option_DynSys_slow.Continuation.kappa_s_exp = kappa_s_exp_slow;
Solver_option_DynSys_slow.Continuation.sigma_Init = sigma_Init;
Solver_option_DynSys_slow.Continuation.sigma_End = sigma_End;
Solver_option_DynSys_slow.Continuation.kappa_sigma_times = kappa_sigma_times;
Solver_option_DynSys_slow.Continuation.kappa_sigma_exp = kappa_sigma_exp;
% solver set (primal gap)
for i = 1 : numel(param_c)
    % discretize OCPEC into a NLP problem
    NLP_option_i.relaxation_problem = 'gap_constraint_based'; 
    NLP_option_i.gap_constraint_relaxation_strategy = 'generalized_primal_gap'; 
    NLP_option_i.gap_func_implementation = 'symbolic';
    NLP_option_i.primal_gap_param_c = param_c{i};
    NLP_i = NLP_Formulation(OCPEC, NLP_option_i);
    % create DynSys_Based_Solver
    solver_i = DynSys_Based_Solver(OCPEC, NLP_i, Solver_option_DynSys_slow);
    % save
    solver_name{end + 1} = ['Gap (primal, symbolic, c = ' num2str(param_c{i}), ' slow, DynSys-Based'];
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
    solver_i = DynSys_Based_Solver(OCPEC, NLP_i, Solver_option_DynSys_slow);
    % save
    solver_name{end + 1} =  ['Gap (D, symbolic, a = ' num2str(param_a{i}) ', b = ' num2str(param_b{i}) ')', ' slow, DynSys-Based'];
    solver_set{end + 1} = solver_i;    
end

%% run test
% init record
rec.name = solver_name;
rec.z_Init = {};
rec.z_Opt = {};
rec.Info = {};
% run
for i = 1 : numel(solver_set)
    solver_i = solver_set{i};
    z_Init_i = ones(solver_i.NLP.Dim.z, 1);
    s_Init_i = s_Init_set{i};
    % solve
    [z_Opt_i, Info_i] = solver_i.solve_NLP(z_Init_i, s_Init_i, s_End);
    rec.z_Init{end + 1} = z_Init_i;
    rec.z_Opt{end + 1} = z_Opt_i;
    rec.Info{end + 1} = Info_i;
end
save('Data_test_IPOPT_DynSys.mat', 'rec')
