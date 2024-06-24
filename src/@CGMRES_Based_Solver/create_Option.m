function Option = create_Option()
%UNTITLED19 Summary of this function goes here
%   Detailed explanation goes here

%% Option for stage 1: non-interior-point (NIP) method
Option.NIP.printLevel = 2; % 0: print nothing;  
                       % 1: print results
                       % 2: print results and iteration log (should specified recordLevel as 1)
% tolerance
Option.NIP.maxIterNum = 200;

Option.NIP.KKT_scaling_max = 1;
Option.NIP.tol.KKT_error_primal = 1e-3;
Option.NIP.tol.KKT_error_dual = 1e-3;
Option.NIP.tol.KKT_error_complementarity = 1e-3;
Option.NIP.tol.KKT_error_total = 1e-3;
Option.NIP.tol.dYNorm = 1e-6;

% singularity regularization parameter for KKT matrix 
Option.NIP.RegParam.nu_h = 1e-7; % rank-deficiency of equality constraint h
Option.NIP.RegParam.nu_c = 1e-8; % non-negative definite of diagonal matrix related to inequality constraint c
Option.NIP.RegParam.nu_H = 1e-8; % non-positive definite of Hessian matrix

% Lagrangian Hessian approximation
Option.NIP.HessianApproximation = 'Gauss_Newton'; % 'Exact', 'Gauss_Newton', 'Quasi_Newton'

% evaluate search direction

% merit line search
Option.NIP.LineSearch.beta_Init = 1; % initial penalty parameter
Option.NIP.LineSearch.rho = 0.1; % desired extend for the negativity of merit function directional derivative
Option.NIP.LineSearch.stepSize_Min = 0.001;
Option.NIP.LineSearch.stepSize_DecayRate = 0.5;% choose in (0,1)
Option.NIP.LineSearch.nu_D = 1e-4;% desired merit function reduction, default 1e-4 

%% Option for stage 2: C/GMRES method



end