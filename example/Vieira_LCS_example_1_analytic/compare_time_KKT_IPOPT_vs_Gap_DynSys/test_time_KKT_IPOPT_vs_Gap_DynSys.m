clear all
clc

%% create OCPEC (all simulation use the same OCPEC)
timeHorizon = 1;
nStages = 100;
OCPEC = OCPEC_Vieira_LCS_analytic();
OCPEC.timeHorizon = timeHorizon;
OCPEC.nStages = nStages;
OCPEC.timeStep = OCPEC.timeHorizon ./ OCPEC.nStages;

%% specify NLP and solver parameter
% KKT relaxation strategy for all update
KKT_relaxation_strategy = {'Scholtes', 'Lin_Fukushima', 'Kadrani', ...
    'Steffensen_Ulbrich', 'Kanzow_Schwartz'};
KKT_IPOPT_name = {...
    'KKT (Scholtes), IPOPT-based',...
    'KKT (Lin-Fukushima), IPOPT-based',...
    'KKT (Kadrani), IPOPT-based', ...
    'KKT (Steffensen-Ulbrich), IPOPT-based',...
    'KKT (Kanzow-Schwartz), IPOPT-based'};
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
% parameter update
kappa_s_times_KKT = 0.1;% fast
kappa_s_times_gap = 0.9;% slow
kappa_s_exp = 1.0;
kappa_sigma_times = 0.9;
kappa_sigma_exp = 1.0;
% tolerance
KKT_error_tol = 1e-6;
VI_nat_res_tol = 1e-2;
% dynamical system parameter
dtau = 0.001;
epsilon = 1000;
% init
solver_name = {};
solver_set = {};

%% create solver set (KKT, IPOPT_Based)
% solver option
Solver_option_KKT_IPOPT = IPOPT_Based_Solver.create_Option();
Solver_option_KKT_IPOPT.Continuation.kappa_s_times = kappa_s_times_KKT;
Solver_option_KKT_IPOPT.Continuation.kappa_s_exp = kappa_s_exp;
Solver_option_KKT_IPOPT.Continuation.tol.VI_nat_res = VI_nat_res_tol;
% solver set
for i = 1 : numel(KKT_relaxation_strategy)
    % discretize OCPEC into a NLP problem
    NLP_option_i.relaxation_problem = 'KKT_based'; 
    NLP_option_i.KKT_complementarity_relaxation_strategy = KKT_relaxation_strategy{i};
    NLP_i = NLP_Formulation(OCPEC, NLP_option_i);
    % create IPOPT_Based_Solver
    solver_i = IPOPT_Based_Solver(OCPEC, NLP_i, Solver_option_KKT_IPOPT);
    % save
    solver_name{end + 1} = KKT_IPOPT_name{i};
    solver_set{end + 1} = solver_i;    
end

%% create solver set (gap, DynSys_Based)
% solver option
Solver_option_gap_DynSys = DynSys_Based_Solver.create_Option();
Solver_option_gap_DynSys.Continuation.dtau = dtau;
Solver_option_gap_DynSys.Continuation.epsilon = epsilon;
Solver_option_gap_DynSys.Continuation.integration_method = 'RK4'; % 'explitic_Euler', 'RK4'
Solver_option_gap_DynSys.Continuation.tol.KKT_error = KKT_error_tol;
Solver_option_gap_DynSys.Continuation.tol.VI_nat_res = VI_nat_res_tol;
Solver_option_gap_DynSys.Continuation.kappa_s_times = kappa_s_times_gap;
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
    solver_name{end + 1} = ['Gap (primal, c = ' num2str(param_c{i}),')',  ' DynSys-Based'];
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
    solver_name{end + 1} =  ['Gap (D, a = ' num2str(param_a{i}) ', b = ' num2str(param_b{i}) ')', ' DynSys-Based'];
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
save('Data_test_time_KKT_IPOPT_vs_Gap_DynSys.mat', 'rec')

%% Print time
Data_IPOPT_DynSys = load('Data_test_time_KKT_IPOPT_vs_Gap_DynSys');
disp('-----------------------------------------------------------------------------------------------------------------------------------------')
disp('  ID  | time(1:end)[s] | time(1)[s] | time(2:end)[s] | step | time(ave, 2:end)[s] |    cost    | KKT error  | VI_nat_res |  method name ')
for i = 1 : numel(Data_IPOPT_DynSys.rec.name)
    time_i    = Data_IPOPT_DynSys.rec.Info{i}.time;
    timeLog_i = Data_IPOPT_DynSys.rec.Info{i}.Log.time;
    step_i    = Data_IPOPT_DynSys.rec.Info{i}.continuationStepNum;
    cost_i      = Data_IPOPT_DynSys.rec.Info{i}.cost;
    KKT_error_i = Data_IPOPT_DynSys.rec.Info{i}.KKT_error;
    VI_natres_i = Data_IPOPT_DynSys.rec.Info{i}.VI_natural_residual;
    disp([' ',...
        num2str(i, '%10.4d'),  ' |   ', ...
        num2str(time_i, '%10.4e'), '   | ',...
        num2str(timeLog_i(1), '%10.4e'), ' |   ',...
        num2str(time_i - timeLog_i(1), '%10.4e'), '   | ',...
        num2str(step_i, '%10.4d'), ' |    ', ...
        num2str((time_i - timeLog_i(1))/step_i, '%10.4e'), '       | ', ...
        num2str(cost_i, '%10.4d'), ' | ', ...
        num2str(KKT_error_i, '%10.4d'), ' | ', ...
        num2str(VI_natres_i, '%10.4d'), ' | ', ...
        Data_IPOPT_DynSys.rec.name{i}])
end