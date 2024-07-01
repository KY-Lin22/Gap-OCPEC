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
% tolerance
Option.Continuation.tol.VI_nat_res = 1e-4;
% relaxation parameter
Option.Continuation.kappa_s_times = 0.8; % update
Option.Continuation.kappa_s_exp = 1; % update

end

