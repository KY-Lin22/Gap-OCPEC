function Option = createSolverOption(self)
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

%% option for homotopy
Option.Homotopy.kappa_times = 0.9;
Option.Homotopy.kappa_exp = 1.1;

end

