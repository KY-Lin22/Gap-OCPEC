function Option = create_Option()
%UNTITLED19 Summary of this function goes here
%   Detailed explanation goes here

%% function and gradient evaluation
% singularity regularization parameter for KKT matrix 
Option.KKT.RegParam.nu_h = 1e-6; % rank-deficiency of equality constraint h
Option.KKT.RegParam.nu_c = 1e-6; % non-negative definite of diagonal matrix related to inequality constraint c
Option.KKT.RegParam.nu_H = 1e-6; % non-positive definite of Hessian matrix
% Lagrangian Hessian approximation
Option.KKT.Hessian_approximation = 'Exact'; % 'Exact', 'Gauss_Newton', 'Quasi_Newton'

%% Option for stage 1: solve first parameterized NLP by non-interior-point (NIP) method
Option.NIP.printLevel = 2; % 0: print nothing;  
                       % 1: print results
                       % 2: print results and iteration log (should specified recordLevel as 1)
% tolerance
Option.NIP.maxIterNum = 300;
Option.NIP.tol.KKT_error = 1e-2;
Option.NIP.tol.dYNorm = 1e-6;
% evaluate search direction

% merit line search
Option.NIP.LineSearch.beta_Init = 1; % initial penalty parameter
Option.NIP.LineSearch.rho = 0.1; % desired extend for the negativity of merit function directional derivative
Option.NIP.LineSearch.stepSize_Min = 1e-4;
Option.NIP.LineSearch.stepSize_DecayRate = 0.5;% choose in (0,1)
Option.NIP.LineSearch.nu_D = 1e-4;% desired merit function reduction, default 1e-4 

%% Option for stage 1: solve first parameterized NLP by IPOPT (ref: https://coin-or.github.io/Ipopt/OPTIONS.html)
% print
Option.IPOPT_Solver.print_time = true;
Option.IPOPT_Solver.record_time = true;
Option.IPOPT_Solver.ipopt.print_level = 3;
% tolerance
Option.IPOPT_Solver.ipopt.tol = 1e-2;% default 1e-8
% Option.IPOPT_Solver.ipopt.compl_inf_tol = 1e-4; % default 1e-4
Option.IPOPT_Solver.ipopt.max_iter = 3000; % default 3000
Option.IPOPT_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
% barrier parameter
Option.IPOPT_Solver.ipopt.mu_strategy = 'monotone'; % default: 'monotone', 'adaptive'
% Option.IPOPT_Solver.ipopt.mu_min = 5e-3; % default 0 (for 'adaptive' mu strategy, need to specified by 0.5*sigma_Init^2)

%% Option for stage 2: solve differential equation
% fictitious time 
Option.Continuation.dtau = 0.001;
% method to solve the first parameterized NLP and differential equation
Option.Continuation.first_NLP_solve = 'IPOPT'; % 'non_interior_point', 'IPOPT'
Option.Continuation.differential_equation_solve = 'direct'; % 'FDGMRES', 'direct'
% tolerance
Option.Continuation.tol.KKT_error = 1e-2;
Option.Continuation.tol.VI_nat_res = 1e-2;
% homotopy (relaxation parameter)
Option.Continuation.kappa_s_times = 0.8; % update
% homotopy (smoothing parameter)
Option.Continuation.sigma_Init = 1e-1;
Option.Continuation.sigma_End = 1e-3;
Option.Continuation.kappa_sigma_times = 0.8; % update
Option.Continuation.kappa_sigma_exp = 1.2; % update
% FDGMRES method
Option.Continuation.FDGMRES.h_FD = 1e-9; % stepsize of forward difference approximation for the product of Jacobians and vectors
Option.Continuation.FDGMRES.k_max = 10; % GMRES max iteration number

end