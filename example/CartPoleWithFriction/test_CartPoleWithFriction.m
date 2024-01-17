clear all
clc

%% construct OCPEC problem
OCPEC = OCPEC_CartPoleWithFriction();

%% discretize OCPEC into a NLP problem
Option.relaxProbType = 'generalized_primal_gap_constraint_based'; % 'generalized_primal_gap_constraint_based', 
                                                         % 'generalized_D_gap_constraint_based'
                                                         % 'KKT_based'
Option.D_gap_param.a = 0.9;
Option.D_gap_param.b = 1.3;
Option.stronglyConvexFuncType = 'quadratic';

NLP = NLP_Formulation(OCPEC, Option);

%% create solver
solver = IPOPT_Based_Solver(OCPEC, NLP);

%% problem solve
Z_Init = zeros(NLP.Dim.z_Node(4), OCPEC.nStages);
Z_Init(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :) = randn(OCPEC.Dim.u, OCPEC.nStages);
z_Init = reshape(Z_Init, [], 1);

s_Init = 5e-1;
s_End = 1e-3; 

p_Init = s_Init;
p_End = s_End;
[z_Opt, Info] = solver.solveNLP(z_Init, p_Init, p_End);

%%
plotResult_CartPoleWithFriction(OCPEC, NLP, z_Opt)
