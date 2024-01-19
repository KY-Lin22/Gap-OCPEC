clear all
clc

%% construct OCPEC problem
OCPEC = OCPEC_Vieira_LCS_state_jump();

%% discretize OCPEC into a NLP problem
Option.relaxProbType = 'generalized_D_gap_constraint_based'; % 'generalized_primal_gap_constraint_based', 
                                                         % 'generalized_D_gap_constraint_based'
                                                         % 'KKT_based'
Option.D_gap_param.a = 0.9;
Option.D_gap_param.b = 1.1;
Option.stronglyConvexFuncType = 'quadratic';

NLP = NLP_Formulation(OCPEC, Option);

%% create solver
solver = IPOPT_Based_Solver(OCPEC, NLP);
solver.Option.Homotopy.kappa_times = 0.9;
solver.Option.Homotopy.kappa_exp = 1.1;

%% problem solve
z_Init = zeros(NLP.Dim.z, 1);

s_Init = 1e-1;
s_End = 1e-3; 

p_Init = s_Init;
p_End = s_End;
[z_Opt, Info] = solver.solveNLP(z_Init, p_Init, p_End);

%%
plotResult_Vieira_LCS_state_jump(OCPEC, NLP, z_Opt)
