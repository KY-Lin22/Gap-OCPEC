clear all
clc

%% create OCPEC (all simulation use the same OCPEC)
timeHorizon = 1;
nStages = 100;
OCPEC = OCPEC_Vieira_LCS_analytic();
OCPEC.timeHorizon = timeHorizon;
OCPEC.nStages = nStages;
OCPEC.timeStep = OCPEC.timeHorizon ./ OCPEC.nStages;

%% specify NLP parameter
% primal gap parameter
param_c = {1};
% D gap parameter
param_a = {0.1,  0.5, 0.9};
param_b = {10,   2,   1.1};

%% create solver set (gap, DynSys_Based)
% solver option
Solver_option_gap_DynSys = DynSys_Based_Solver.create_Option();
Solver_option_gap_DynSys.Continuation.s_Init = 1e0;
Solver_option_gap_DynSys.Continuation.s_End = 1e-6;
Solver_option_gap_DynSys.Continuation.sigma_Init = 1e-2;
Solver_option_gap_DynSys.Continuation.sigma_End = 1e-6;
Solver_option_gap_DynSys.Continuation.epsilon_T = 100;
Solver_option_gap_DynSys.Continuation.epsilon_p = 10;
Solver_option_gap_DynSys.Continuation.dtau = 0.01;
Solver_option_gap_DynSys.Continuation.l_Max = 500;
Solver_option_gap_DynSys.Continuation.integration_method = 'RK4'; % 'explitic_Euler', 'RK4'
Solver_option_gap_DynSys.Continuation.tol.KKT_error = 1e-16;
Solver_option_gap_DynSys.Continuation.tol.VI_nat_res = 1e-16;
% init
solver_name = {};
solver_set = {};
% solver set (primal gap)
for i = 1 : numel(param_c)
    % discretize OCPEC into a NLP problem
    NLP_option_i.relaxation_problem = 'gap_constraint_based'; 
    NLP_option_i.gap_constraint_relaxation_strategy = 'generalized_primal_gap'; 
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
repeat_num = 1;
for i = 1 : numel(solver_set)
    solver_i = solver_set{i};
    z_Init_i = ones(solver_i.NLP.Dim.z, 1);
    % solve
    for j = 1 : repeat_num
        [z_Opt_i, Info_i] = solver_i.solve_NLP(z_Init_i);
    end    
    rec.z_Init{end + 1} = z_Init_i;
    rec.z_Opt{end + 1} = z_Opt_i;
    rec.Info{end + 1} = Info_i;
end
save('Data_test_converge_Gap_DynSys.mat', 'rec')
