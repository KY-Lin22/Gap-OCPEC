clear all
clc

%% construct OCPEC problem
OCPEC = OCPEC_AffineDVI();

%% discretize OCPEC into a NLP problem
Option.relaxProbType = 'generalized_D_gap_constraint_based'; % 'generalized_primal_gap_constraint_based', 
                                                         % 'generalized_D_gap_constraint_based'
                                                         % 'KKT_based'
Option.D_gap_param.a = 0.9;
Option.D_gap_param.b = 1.3;
Option.stronglyConvexFuncType = 'quadratic';

NLP = NLP_Formulation(OCPEC, Option);

%% create solver
solver = IPOPT_Based_Solver(OCPEC, NLP);

%% problem solve
z_Init = ones(NLP.Dim.z, 1);

s_Init = 1e-1;
s_End = 1e-4; 

p_Init = s_Init;
p_End = s_End;
[z_Opt, Info] = solver.solveNLP(z_Init, p_Init, p_End);

%%
plotResult_AffineDVI(OCPEC, NLP, z_Opt)
