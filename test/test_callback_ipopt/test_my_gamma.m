clear all
clc
import casadi.*

% opti = casadi.Opti();
% x = opti.variable();
% opti.set_initial(x, 2);
% 
% opti.minimize(mygamma(x));
% opti.subject_to(mygamma(sin(x))>=0.5);
% 
% opti.solver('ipopt');
% sol = opti.solve();
mygamma = MyGamma('myGamma', struct('enable_fd',true));
x = MX.sym('x', 1, 1);
x_0 = 2;
f = mygamma(x);
g = mygamma(sin(x));
lbg = 0;
ubg = inf;
Prob = struct('x', x, 'f', f, 'g', g);
Option = struct;
Option.ipopt.bound_relax_factor = 0;
% Option.ipopt.constr_viol_tol = 1e-12;
% Option.ipopt.print_level = 0;
Option.ipopt.tol = 1e-8; % default 1e-8
% Option.ipopt.max_iter = 3000; % default 3000
solver = casadi.nlpsol('solver', 'ipopt', Prob, Option);

%%
solution = solver('x0', x_0, 'lbg', lbg, 'ubg', ubg);
% extract solution
x_Opt = full(solution.x);
g_Opt = full(solution.g);