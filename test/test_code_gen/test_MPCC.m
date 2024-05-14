clear all
clc

import casadi.*
% variable
lambda = MX.sym('lambda', 1, 1);
eta = MX.sym('eta', 1, 1);
% relaxation parameter
s = 1e-8;
% cost
J = (lambda - 1)^2 + (eta - 1)^2;
% symbolic expression for primal gap function
param_c = 1;
stationary_point_c = lambda - (1/param_c) * eta;
omega_c = max(0, stationary_point_c);
p_c = 0.5*lambda^2 - 0.5*omega_c^2 + lambda*(omega_c - lambda);
phi_c = eta' * (lambda - omega_c) + param_c * p_c;
% function object for primal gap function
phi_c_func = Function('phi_c_func', {lambda, eta}, {phi_c});
gap_func_implementation = 'codegen_fd';
switch gap_func_implementation
    case 'symbolic'
        % direct use symbolic expression
        phi = phi_c;
    case 'codegen_fd'
        % code gen
        codeGenerator = CodeGenerator('phi_c_func_codegen_fd', struct('mex',true));
        codeGenerator.add(phi_c_func);
        codeGenerator.generate;
        mex phi_c_func_codegen_fd.c -largeArrayDims
        % load back as external function
        phi_c_func_import = Importer('phi_c_func_codegen_fd.mexw64','dll');
        phi_func = external('phi_c_func',phi_c_func_import, struct('enable_fd', true));
        disp(phi_func)
        phi = phi_func(lambda, eta);        
    case 'codegen_jac'
        % NOTE: only available for 'limit_memory' hessian approximation !!!
        % custom jacobian for primal gap function (currently use automatic symbolic one, later use the explicit one)
        phi_c_jac_lambda = jacobian(phi_c, lambda);
        phi_c_jac_eta = jacobian(phi_c, eta);     
        % function object for custom jacobian of primal gap function 
        nominal_out = MX(1,1); % a empty symbolic variable for 1 x 1 phi_c and not make use of it 
        phi_c_jac_func = Function('jac_phi_c_func',...% naming convention: names as jac_fname
            {lambda, eta, nominal_out}, {phi_c_jac_lambda, phi_c_jac_eta});
        % code gen
        codeGenerator = CodeGenerator('phi_c_func_codegen_jac', struct('mex',true));
        codeGenerator.add(phi_c_func);
        codeGenerator.add(phi_c_jac_func);
        codeGenerator.generate;
        mex phi_c_func_codegen_jac.c -largeArrayDims
        phi_c_func_import = Importer('phi_c_func_codegen_jac.mexw64', 'dll');
        phi_func = external('phi_c_func', phi_c_func_import);
        disp(phi_func)
        phi = phi_func(lambda, eta);        
end
% constraint
g = [lambda;...
    s - phi];
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