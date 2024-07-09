function Option = create_Option()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% option for NLP solver (IPOPT)
% refer to: https://coin-or.github.io/Ipopt/OPTIONS.html

% print
Option.IPOPT_Solver.print_time = false;
Option.IPOPT_Solver.record_time = true;

Option.IPOPT_Solver.ipopt.print_level = 0;

% tolerance
Option.IPOPT_Solver.ipopt.tol = 1e-8;% default 1e-8
Option.IPOPT_Solver.ipopt.max_iter = 3000; % default 3000
Option.IPOPT_Solver.ipopt.hessian_approximation = 'exact'; % 'exact' (default), 'limited-memory'
% Option.NLP_Solver.ipopt.acceptable_tol = 1e-6; % default 1e-6

% barrier parameter
% Option.NLP_Solver.ipopt.mu_strategy = 'adaptive';

% NLP bound 
% Option.NLP_Solver.ipopt.bound_relax_factor = 0;

%% option for continuation method
% Continuation (relaxation parameter)
Option.Continuation.s_Init = 1e0;
Option.Continuation.s_End = 1e-12;
% update relaxation parameter
Option.Continuation.kappa_s_times = 0.1; 
Option.Continuation.kappa_s_exp = 1; 
% tolerance
Option.Continuation.tol.VI_nat_res = 1e-4;

end

