clear all
clc

%% create OCPEC, NLP, and solver
% construct OCPEC problem
OCPEC = OCPEC_Vieira_LCS_analytic();

% discretize OCPEC into a NLP problem
Option.relaxation_problem = 'gap_constraint_based'; % 'gap_constraint_based', 'KKT_based'
Option.gap_constraint_relaxation_strategy = 'generalized_primal_gap'; % 'generalized_primal_gap', 'generalized_D_gap'
Option.KKT_complementarity_relaxation_strategy = 'Scholtes'; % 'Scholtes', 'Lin_Fukushima', 'Kadrani'
Option.strongly_convex_func = 'quadratic';
Option.gap_func_smooth_param = 0.001;
Option.D_gap_param_a = 0.8;
Option.D_gap_param_b = 1.2;
Option.gap_func_auxiliary_variable_penalty = 'L2'; % 'none', 'L1', 'L2'
NLP = NLP_Formulation(OCPEC, Option);

% create solver
solver = IPOPT_Based_Solver(OCPEC, NLP);

% create initial guess
% z_Init = solver.create_initial_guess();
z_Init = zeros(NLP.Dim.z, 1);

%% parameter and problem solve
% parameter
s_Init = 1e-1;
s_End = 1e-8; 
mu_Init = 1e0;
mu_End = 1e4;
p_Init = [s_Init; mu_Init];
p_End = [s_End; mu_End];
solver.Option.Homotopy.kappa_s_times = 0.8;
solver.Option.Homotopy.kappa_s_exp = 1.2;
solver.Option.Homotopy.VI_nat_res_tol = 1e-2;

% solve
[z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End);

Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCPEC.nStages);

%% show result
plotResult_Vieira_LCS_analytic(OCPEC, NLP, z_Opt)
