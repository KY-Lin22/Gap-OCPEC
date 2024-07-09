clear all
clc

%% create OCPEC, NLP and solver
% construct OCPEC problem
OCPEC = OCPEC_Vieira_LCS_analytic();

% discretize OCPEC into a NLP problem
NLP_option.relaxation_problem = 'KKT_based'; % 'gap_constraint_based', 'KKT_based'
NLP_option.gap_constraint_relaxation_strategy = 'generalized_primal_gap'; % 'generalized_primal_gap', 'generalized_D_gap'
NLP_option.KKT_complementarity_relaxation_strategy = 'Scholtes'; % 'Scholtes', 'Lin_Fukushima', 'Kadrani', 'Steffensen_Ulbrich', 'Kanzow_Schwartz'
NLP_option.gap_func_implementation = 'symbolic'; % 'symbolic', 'codegen_fd', 'codegen_jac'
NLP_option.strongly_convex_func = 'quadratic';
NLP_option.primal_gap_param_c = 1;
NLP_option.D_gap_param_a = 0.1;
NLP_option.D_gap_param_b = 10;
NLP = NLP_Formulation(OCPEC, NLP_option);

% create solver
solver_type = 'DynSys_Based'; % 'IPOPT_Based', 'DynSys_Based'
switch solver_type
    case 'IPOPT_Based'
        % create IPOPT_Based_Solver
        Solver_option = IPOPT_Based_Solver.create_Option();
        Solver_option.Continuation.s_Init = 1e0;
        Solver_option.Continuation.s_End = 1e-12;
        Solver_option.Continuation.kappa_s_times = 0.1; % fast update
        Solver_option.Continuation.kappa_s_exp = 1;
        Solver_option.Continuation.tol.VI_nat_res = 1e-4;
        solver = IPOPT_Based_Solver(OCPEC, NLP, Solver_option);
    case 'DynSys_Based'
        % create DynSys_Based_Solver
        Solver_option = DynSys_Based_Solver.create_Option();
        Solver_option.Continuation.s_Init = 1e0;
        Solver_option.Continuation.s_End = 1e-12;
        Solver_option.Continuation.sigma_Init = 1e-2;
        Solver_option.Continuation.sigma_End = 1e-6;
        Solver_option.Continuation.epsilon_T = 100;
        Solver_option.Continuation.epsilon_p = 50;
        Solver_option.Continuation.dtau = 0.01;
        Solver_option.Continuation.l_Max = 500;
        Solver_option.Continuation.integration_method = 'RK4'; % 'explitic_Euler', 'RK4'
        Solver_option.Continuation.tol.KKT_error = 1e-6;
        Solver_option.Continuation.tol.VI_nat_res = 1e-4;
        solver = DynSys_Based_Solver(OCPEC, NLP, Solver_option);
end

%% parameter and problem solve
% create initial guess
z_Init = ones(NLP.Dim.z, 1);

% solve
[z_Opt, Info] = solver.solve_NLP(z_Init);

%% show result
plotResult_Vieira_LCS_analytic(OCPEC, NLP, z_Opt)
