clear all
clc

%% create OCPEC (all simulation use the same OCPEC)
timeHorizon = 1;
nStages = 2000;
OCPEC = OCPEC_Vieira_LCS_lambda_penalty();
OCPEC.timeHorizon = timeHorizon;
OCPEC.nStages = nStages;
OCPEC.timeStep = OCPEC.timeHorizon ./ OCPEC.nStages;

%% specify NLP parameter
% KKT relaxation strategy for all update
KKT_relaxation_strategy = {'Scholtes', 'Lin_Fukushima', 'Kadrani', ...
    'Steffensen_Ulbrich', 'Kanzow_Schwartz'};
KKT_IPOPT_name = {'Scholtes','Lin-Fukushima','Kadrani', ...
    'Steffensen-Ulbrich','Kanzow-Schwartz'};

%% create solver set (KKT, IPOPT_Based)
% solver option
Solver_option_KKT_IPOPT = IPOPT_Based_Solver.create_Option();
Solver_option_KKT_IPOPT.Continuation.s_Init = 1e0;
Solver_option_KKT_IPOPT.Continuation.s_End = 1e-3;
Solver_option_KKT_IPOPT.Continuation.update_rule = 'dynamics'; % 'times_exp', 'dynamics'
Solver_option_KKT_IPOPT.Continuation.epsilon_s = 10;
Solver_option_KKT_IPOPT.Continuation.dtau = 0.20;
Solver_option_KKT_IPOPT.Continuation.l_Max = 25;
Solver_option_KKT_IPOPT.Continuation.tol.VI_nat_res = 1e-16;
% init
solver_name = {};
solver_set = {};
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
save('Data_test_converge_KKT_IPOPT.mat', 'rec')

