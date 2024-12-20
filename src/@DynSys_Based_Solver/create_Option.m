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

%% Option for stage 1: solve first parameterized NLP by IPOPT (ref: https://coin-or.github.io/Ipopt/OPTIONS.html)
% print
Option.IPOPT_Solver.print_time = true;
Option.IPOPT_Solver.record_time = true;
Option.IPOPT_Solver.ipopt.print_level = 3;
% tolerance
Option.IPOPT_Solver.ipopt.tol = 1e-8;% default 1e-8
% Option.IPOPT_Solver.ipopt.compl_inf_tol = 1e-4; % default 1e-4
Option.IPOPT_Solver.ipopt.max_iter = 3000; % default 3000
Option.IPOPT_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
% barrier parameter
Option.IPOPT_Solver.ipopt.mu_strategy = 'monotone'; % default: 'monotone', 'adaptive'
% Option.IPOPT_Solver.ipopt.mu_min = 5e-3; % default 0 (for 'adaptive' mu strategy, need to specified by 0.5*sigma_Init^2)

%% Option for stage 2: solve differential equation
% Continuation (relaxation and smoothing parameter)
Option.Continuation.s_Init = 1e0;
Option.Continuation.s_End = 1e-12;
Option.Continuation.sigma_Init = 1e-2;
Option.Continuation.sigma_End = 1e-6;
% stabilization parameter for KKT and parameter dynamics
Option.Continuation.epsilon_T = 100;
Option.Continuation.epsilon_p = 10; 
% fictitious integration timestep and number of integration step (continuation step)
Option.Continuation.dtau = 0.01;
Option.Continuation.l_Max = 500;
% explicit method to integrating the KKT differential equation
Option.Continuation.integration_method = 'RK4'; % 'explitic_Euler', 'RK4'
% tolerance
Option.Continuation.tol.KKT_error = 1e-4;
Option.Continuation.tol.VI_nat_res = 1e-2;

end