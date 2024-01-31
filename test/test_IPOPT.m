clear all
clc

import casadi.*
% variables
lambda = SX.sym('lambda', 1, 1);
eta = SX.sym('eta', 1, 1);
w = SX.sym('w', 1, 1);
% parameter
s = SX.sym('s', 1, 1);
% cost
J = (lambda - 1)^2 + (eta - 1)^2;
% constraint
relax_prob_type = 'primal_gap';
switch relax_prob_type
    case 'primal_gap'
        phi = 0.5 * (eta^2 - (max(0, eta - lambda))^2);
        g = [phi - w;...
            lambda;...
            s - w];      
        lbg = [0; 0; 0];        
        ubg = [0; inf; inf];        
    case 'D_gap'    

end
phi_func = Function('phi_func', {eta, lambda}, {phi}, {'eta', 'lambda'}, {'phi'});
% solver
Prob = struct('x', [lambda; eta; w], 'f', J, 'g', g, 'p', s);
Option = struct;
Option.ipopt.bound_relax_factor = 0;
% Option.ipopt.constr_viol_tol = 1e-12;
% Option.ipopt.print_level = 0;
Option.ipopt.tol = 1e-8; % default 1e-8
% Option.ipopt.max_iter = 3000; % default 3000
solver = casadi.nlpsol('solver', 'ipopt', Prob, Option);

%% solve
s_0 = 1e-12;
x_0 = [-1; -1; 1];
solution = solver('x0', x_0, 'p', s_0, 'lbg', lbg, 'ubg', ubg);
% extract solution
x_Opt = full(solution.x);
g_Opt = full(solution.g);
phi_Opt = full(phi_func(x_Opt(1), x_Opt(2)));