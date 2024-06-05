clear all
clc
delete ThreeCarts.mp4

%% create OCPEC, NLP, and solver
% construct OCPEC problem
OCPEC = OCPEC_ThreeCartsSystem();

% discretize OCPEC into a NLP problem
NLP_option.relaxation_problem = 'gap_constraint_based'; % 'gap_constraint_based', 'KKT_based'
NLP_option.gap_constraint_relaxation_strategy = 'generalized_primal_gap'; % 'generalized_primal_gap', 'generalized_D_gap'
NLP_option.KKT_complementarity_relaxation_strategy = 'Scholtes'; % 'Scholtes', 'Lin_Fukushima', 'Kadrani', 'Steffensen_Ulbrich', 'Kanzow_Schwartz'
NLP_option.gap_func_implementation = 'symbolic'; % 'symbolic', 'codegen_fd', 'codegen_jac'
NLP_option.strongly_convex_func = 'quadratic';
NLP_option.primal_gap_param_c = 0.1;
NLP_option.D_gap_param_a = 0.1;
NLP_option.D_gap_param_b = 10;
NLP = NLP_Formulation(OCPEC, NLP_option);

% create solver
Solver_option = IPOPT_Based_Solver.create_Option();
Solver_option.NLP_Solver.ipopt.tol = 1e-6;% default 1e-8
Solver_option.NLP_Solver.ipopt.max_iter = 2000; % default 3000
Solver_option.NLP_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
Solver_option.Homotopy.kappa_s_times = 0.8;
Solver_option.Homotopy.VI_nat_res_tol = 1e-2;
solver = IPOPT_Based_Solver(OCPEC, NLP, Solver_option);

% create initial guess
z_Init = randn(NLP.Dim.z, 1);

%% parameter and problem solve
% parameter
s_Init = 1e1;
s_End = 1e-8; 
p_Init = s_Init;
p_End = s_End;

% solve
[z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End);

Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCPEC.nStages);

%% show result
plotResult_ThreeCartsSystem(OCPEC, NLP, z_Opt)

animateTrajectory_ThreeCartsSystem(OCPEC, NLP, z_Opt)