clear all
clc

%% create OCPEC, NLP, and solver
% construct OCPEC problem
OCPEC = OCPEC_Vieira_LCS_higher_dim();

% discretize OCPEC into a NLP problem
NLP_option.relaxation_problem = 'gap_constraint_based'; % 'gap_constraint_based', 'KKT_based'
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
        Solver_option.Continuation.s_End = 1e-3;
        Solver_option.Continuation.kappa_s_times = 0.1; % fast update
        Solver_option.Continuation.kappa_s_exp = 1;
        Solver_option.Continuation.tol.VI_nat_res = 1e-2;
        solver = IPOPT_Based_Solver(OCPEC, NLP, Solver_option);
    case 'DynSys_Based'
        % create DynSys_Based_Solver
        Solver_option = DynSys_Based_Solver.create_Option();
        Solver_option.Continuation.s_Init = 1e-1;
        Solver_option.Continuation.s_End = 1e-3;
        Solver_option.Continuation.epsilon_T = 100;
        Solver_option.Continuation.epsilon_s = 10;
        Solver_option.Continuation.dtau = 0.01;
        Solver_option.Continuation.l_Max = 500;
        Solver_option.Continuation.integration_method = 'explitic_Euler'; % 'explitic_Euler', 'RK4'
        Solver_option.Continuation.tol.KKT_error = 1e-4;
        Solver_option.Continuation.tol.VI_nat_res = 1e-2;
        Solver_option.KKT.RegParam.nu_h = 1e-6; 
        Solver_option.KKT.RegParam.nu_c = 1e-6; 
        Solver_option.KKT.RegParam.nu_H = 1e-3; 
        Solver_option.KKT.Hessian_approximation = 'Gauss_Newton'; % 'Exact', 'Gauss_Newton', 'Quasi_Newton'
        solver = DynSys_Based_Solver(OCPEC, NLP, Solver_option);
end


%% parameter and problem solve
% create initial guess
z_Init = randn(NLP.Dim.z, 1);

% solve
[z_Opt, Info] = solver.solve_NLP(z_Init);

%% show result
plotResult_Vieira_LCS_higher_dim(OCPEC, NLP, z_Opt)
