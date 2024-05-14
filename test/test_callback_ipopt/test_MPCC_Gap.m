clear all
clc

import casadi.*
% variable
lambda = MX.sym('lambda', 1, 1);
eta = MX.sym('eta', 1, 1);
% symbolic expression for primal gap function
param_c = 1;
stationary_point_c = lambda - (1/param_c) * eta;
omega_c = max(0, stationary_point_c);
p_c = 0.5*lambda^2 - 0.5*omega_c^2 + lambda*(omega_c - lambda);
phi_c = eta' * (lambda - omega_c) + param_c * p_c;
phi_c_func = Function('phi_c_func', {lambda, eta}, {phi_c}, {'lambda', 'eta'}, {'phi_c'});
% relaxation parameter
s = 1e-8;
% cost
J = (lambda - 1)^2 + (eta - 1)^2;
% constraint
gap_func_implementation = 'symbolic';
switch gap_func_implementation
    case 'symbolic'
        phi = phi_c;
        g = [lambda;...
            s - phi];
    case 'callback'
        
        opts = struct('enable_fd',true);
        phi = Gap_func('Gap_func', phi_c_func, opts);
        disp(phi)
        g = [lambda;...
            s - phi(lambda, eta)];
end
lbg = [0; 0];
ubg = [inf; inf];

% solver
Prob = struct('x', [lambda; eta], 'f', J, 'g', g);
Option = struct();
% Option.ipopt.print_level = 0;
% Option.ipopt.tol = 1e-8; % default 1e-8
% Option.ipopt.max_iter = 3000; % default 3000
Option.ipopt.hessian_approximation = 'limited-memory'; % 'exact', 'limited-memory'
solver = casadi.nlpsol('solver', 'ipopt', Prob, Option);

%% solve
x_0 = [0; 0.5];
solution = solver('x0', x_0, 'lbg', lbg, 'ubg', ubg);
% extract solution
disp('optimal solution: ')
disp(full(solution.x))
disp('inequality constraint residual: ')
disp(full(solution.g))
disp('gap function residual: ')
disp(full(phi_c_func(solution.x(1), solution.x(2))));