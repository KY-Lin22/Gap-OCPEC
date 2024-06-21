function Option = create_Option()
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%% Basic Option
Option.printLevel = 2; % 0: print nothing;  
                       % 1: print results
                       % 2: print results and iteration log (should specified recordLevel as 1)

%% Option for stage 1: non-interior-point method
% tolerance
Option.maxIterNum = 200;

Option.KKT_scaling_max = 1;
Option.tol.KKT_error_primal = 1e-2;
Option.tol.KKT_error_dual = 1e-2;
Option.tol.KKT_error_complementarity = 1e-2;
Option.tol.KKT_error_total = 1e-2;
Option.tol.dYNorm = 1e-6;

% singularity regularization parameter for KKT matrix 
Option.RegParam.nu_h = 1e-7; % rank-deficiency of equality constraint h
Option.RegParam.nu_c = 1e-8; % non-negative definite of diagonal matrix related to inequality constraint c
Option.RegParam.nu_H = 1e-8; % non-positive definite of Hessian matrix

% Lagrangian Hessian approximation
Option.HessianApproximation = 'Gauss_Newton';

% evaluate search direction

% merit line search
Option.LineSearch.betaInit = 1; % initial penalty parameter
Option.LineSearch.rho = 0.1; % desired extend for the negativity of merit function directional derivative
Option.LineSearch.stepSize_Min = 0.001;
Option.LineSearch.stepSize_DecayRate = 0.5;% choose in (0,1)
Option.LineSearch.nu_D = 1e-4;% desired merit function reduction, default 1e-4 

%% Option for stage 2: C/GMRES method
% homotopy (relaxation parameter)
Option.Homotopy.kappa_s_times = 0.8;
Option.Homotopy.VI_nat_res_tol = 1e-2;

% homotopy (FB smoothing parameter)
Option.Homotopy.sigma_Init = 1e-1;
Option.Homotopy.sigma_End = 1e-3;
Option.Homotopy.kappa_sigma_times = 0.8;
Option.Homotopy.kappa_sigma_exp = 1.2;

% CGMRES
Option.CGMRES.dtau = 0.001; % fictitious time step
Option.CGMRES.h_FD = 1e-9; % forward difference approximation step
Option.CGMRES.k_max = 10;
end