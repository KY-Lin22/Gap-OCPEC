function Option = create_Option(self)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% option for NLP solver (IPOPT)
% refer to: https://coin-or.github.io/Ipopt/OPTIONS.html

% print
Option.NLP_Solver.print_time = false;
Option.NLP_Solver.record_time = true;

Option.NLP_Solver.ipopt.print_level = 0;

% tolerance
Option.NLP_Solver.ipopt.tol = 1e-8;% default 1e-8
Option.NLP_Solver.ipopt.max_iter = 3000; % default 3000
% Option.NLP_Solver.ipopt.acceptable_tol = 1e-6; % default 1e-6

% barrier parameter
% Option.NLP_Solver.ipopt.mu_strategy = 'adaptive';

% NLP bound 
% Option.NLP_Solver.ipopt.bound_relax_factor = 0;

%% option for homotopy
Option.Homotopy.kappa_s_times = 0.01;
Option.Homotopy.VI_nat_res_tol = 1e-2;

end

